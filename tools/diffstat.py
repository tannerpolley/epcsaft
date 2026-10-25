#!/usr/bin/env python3
"""Runtime lines per layer, as laid out by one pinned formatter, against ARCHITECTURE.yaml code budgets.

    tools/diffstat.py --budgets        # current full runtime tree; a ceiling above or below its band fails
"""
from __future__ import annotations
import ast, bisect, io, re, shutil, sys, tokenize
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
CPP_STYLE = TOOLS / "line-count.clang-format"
RUFF = TOOLS.parent / ".venv/bin/ruff"
# A disabled formatter keeps the packed layout, so a disable directive is refused, not counted.
FORMATTER_DISABLE = re.compile(r"clang-format\s+(off|on)\b|#\s*(fmt|yapf)\s*:\s*(off|on|skip|disable|enable)\b")


def band(count: int) -> int:
    return -(-count // 100)


def budget_state(count: int, ceiling: int, layer_sum: int | None = None) -> str:
    """`layer_sum` bounds the total, which may not exceed the sum of the layer ceilings."""
    if count > ceiling:
        return f"FAIL over the ceiling by {count - ceiling}"
    if ceiling - count > band(count):
        return f"FAIL headroom exceeds the band: lower the ceiling to {min(count + band(count), layer_sum or ceiling)}"
    return "ok"


def python_lines(source: str) -> int:
    """Mask actual docstring byte spans, then count lines occupied by executable tokens."""
    tree = ast.parse(source)
    lines = source.encode("utf-8").split(b"\n")
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)) and node.body:
            first = node.body[0]
            if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
                for number in range(first.lineno, first.end_lineno + 1):
                    line = lines[number - 1]
                    start = first.col_offset if number == first.lineno else 0
                    end = first.end_col_offset if number == first.end_lineno else len(line)
                    lines[number - 1] = line[:start] + b" " * (end - start) + line[end:]
    layout = {tokenize.COMMENT, tokenize.NL, tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT,
              tokenize.ENDMARKER, tokenize.ENCODING}
    occupied = set()
    for token in tokenize.generate_tokens(io.StringIO(b"\n".join(lines).decode("utf-8")).readline):
        if token.type == tokenize.ERRORTOKEN and not token.string.isspace():
            raise ValueError(f"unrecognised Python token at {token.start}")
        if token.type not in layout and not (token.type == tokenize.ERRORTOKEN and token.string.isspace()):
            occupied.update(range(token.start[0], token.end[0] + 1))
    return len(occupied)


def cpp_lines(source: str) -> int:
    """Scan preprocessing tokens after line splicing; raw-string data keeps its original bytes."""
    chars, positions, line_numbers = [], [], []
    number = 1
    for match in re.finditer(r"\\\r?\n|[\s\S]", source):
        value = match[0]
        if not value.startswith("\\") or not value.endswith("\n"):
            chars.append(value)
            positions.append(match.start())
            line_numbers.append(number)
        number += value.count("\n")
    text = "".join(chars)
    occupied = set()
    i = 0
    prefix = ""
    while i < len(text):
        start = i
        if text.startswith("//", i):
            end = text.find("\n", i + 2)
            i = len(text) if end < 0 else end
            continue
        if text.startswith("/*", i):
            end = text.find("*/", i + 2)
            if end < 0:
                raise ValueError("unterminated C++ comment")
            i = end + 2
            continue
        if text[i].isspace():
            prefix = "" if text[i] == "\n" else prefix + " "
            i += 1
            continue
        raw = re.match(r'(?:u8|u|U|L)?R"([^ ()\\\t\r\n]{0,16})\(', text[i:])
        if raw:
            # Translation-phase splicing is reverted within a raw-string payload.
            opening = positions[i + raw.end() - 1]
            closing = ")" + raw[1] + '"'
            end = source.find(closing, opening + 1)
            if end < 0:
                raise ValueError("unterminated C++ raw string")
            end += len(closing)
            occupied.update(range(line_numbers[i], line_numbers[i] + source[positions[i]:end].count("\n") + 1))
            i = bisect.bisect_left(positions, end)
            prefix = ""
            continue
        if re.match(r'(?:u8|u|U|L)?R"', text[i:]):
            raise ValueError("invalid C++ raw-string delimiter")
        # A header-name is one token even if it contains comment-like text.
        if text[i] == "<" and re.search(
            r"(?:^\s*(?:#|%:)\s*(?:include|include_next|import)\s*|\b__has_include(?:_next)?\s*\(\s*)$",
            prefix,
        ):
            end = text.find(">", i + 1)
            if end < 0 or "\n" in text[i:end]:
                raise ValueError("unterminated C++ header name")
            i = end + 1
        elif text[i] in "\"'":
            quote = text[i]
            i += 1
            while i < len(text) and text[i] != quote:
                if text[i] in "\r\n":
                    raise ValueError("newline in C++ literal")
                i += 2 if text[i] == "\\" else 1
            if i >= len(text):
                raise ValueError("unterminated C++ literal")
            i += 1
        elif text[i].isdigit() or (text[i] == "." and text[i + 1:i + 2].isdigit()):
            # Preprocessing numbers include exponent signs and C++ digit separators.
            token = re.match(r"(?:[eEpP][+-]|[\w.]|'[\w])+", text[i:])
            i += len(token[0])
        elif text[i].isalpha() or text[i] == "_":
            token = re.match(r"\w+", text[i:])
            i += len(token[0])
        elif text[i] in "{}[]();:?.~!+-*/%<>=&|^,#":
            i += 1
        else:
            raise ValueError(f"unrecognised C++ input on line {line_numbers[i]}")
        occupied.update(line_numbers[start:i])
        prefix += text[start:i]
    return len(occupied)


def formatted_lines(path: Path) -> int:
    """Executable formatted lines; source formatting and formatter-disable refusals remain unchanged."""
    import subprocess
    source = path.read_bytes()
    # Newline only: splitlines() also breaks at form feeds and would split a directive the formatter honors.
    for number, line in enumerate(source.decode("utf-8").split("\n"), 1):
        if FORMATTER_DISABLE.search(line):
            raise ValueError(f"{path}:{number}: formatter-disable directive bypasses line counting")
    if path.suffix in {".py", ".pyx", ".pxd"}:
        if not RUFF.is_file():
            raise ValueError(f"missing {RUFF}: run bash .codex/environments/setup.sh")
        # A magic trailing comma does not choose the Python layout (clang-format still honours one in a braced
        # list; code-budgets.md); 120 columns in both languages.
        command = [str(RUFF), "format", "--isolated", "--line-length=120", "--target-version=py313",
                   "--config=format.skip-magic-trailing-comma=true", "--config=lint.isort.split-on-trailing-comma=false",
                   f"--stdin-filename={path.name}", "-"]
    else:
        if not shutil.which("clang-format-18"):
            raise ValueError("missing clang-format-18: install it (README Prerequisites)")
        command = ["clang-format-18", f"--style=file:{CPP_STYLE}", f"--assume-filename={path.name}"]
    done = subprocess.run(command, input=source, capture_output=True, check=False)
    if done.returncode:
        raise ValueError(f"cannot format {path} for counting: {done.stderr.decode(errors='replace').strip()}")
    formatted = done.stdout.decode("utf-8")
    return python_lines(formatted) if path.suffix in {".py", ".pyx", ".pxd"} else cpp_lines(formatted)


def code_budgets(root=None) -> int:
    """Count the complete runtime tree; budgets have one owner in ARCHITECTURE."""
    import json
    import re
    from concurrent.futures import ThreadPoolExecutor
    root = TOOLS.parent if root is None else Path(root)
    architecture = (root / "ARCHITECTURE.yaml").read_text()
    sections = list(re.finditer(r"^code_budgets:", architecture, re.MULTILINE))
    if len(sections) != 1:
        raise ValueError("expected one code_budgets JSON flow mapping in ARCHITECTURE.yaml")
    budgets, _ = json.JSONDecoder().raw_decode(architecture[sections[0].end():].lstrip())
    expected = {"core", "contributions", "properties", "equilibrium", "regression", "surface", "bindings", "python", "total"}
    if set(budgets) != expected:
        raise ValueError("code budgets must classify every runtime layer and total")
    if any(type(n) is not int or n < 0 for n in budgets.values()):
        raise ValueError("each code budget is one non-negative integer ceiling")
    if budgets["total"] > sum(n for k, n in budgets.items() if k != "total"):
        raise ValueError("total budget must not exceed the sum of layer budgets")
    files = []
    suffixes = {".c", ".cc", ".cpp", ".cxx", ".h", ".hh", ".hpp", ".hxx", ".inc", ".inl", ".tpp", ".py", ".pyx", ".pxd"}
    for folder in (root / "engine/native", root / "engine/src/epcsaft"):
        if folder.is_symlink():
            raise ValueError(f"runtime root symlink bypasses source accounting: {folder}")
        if not folder.is_dir():
            raise ValueError(f"missing runtime tree: {folder}")
        for path in sorted(folder.rglob("*")):
            if path.is_symlink():
                raise ValueError(f"runtime symlink bypasses source accounting: {path}")
            if not path.is_file() or (path.suffix == ".pyc" and "__pycache__" in path.parts):
                continue
            relative = path.relative_to(folder)
            layer = "python" if folder.name == "epcsaft" else relative.parts[0]
            if layer not in expected - {"total", "python"} and folder.name != "epcsaft":
                raise ValueError(f"unclassified native source: {path}")
            if path.suffix not in suffixes:
                # The extension's type stub declares native types; the architecture check owns its names.
                if path.name not in {".gitkeep", "py.typed", "AGENTS.md", "README.md", "_core.pyi"}:
                    raise ValueError(f"unclassified runtime file: {path}")
                continue
            files.append((layer, path))
    counts = dict.fromkeys(expected, 0)
    with ThreadPoolExecutor() as pool:
        for (layer, _), lines in zip(files, pool.map(formatted_lines, [path for _, path in files])):
            counts[layer] += lines
    counts["total"] = sum(counts.values())
    failed = False
    print("runtime formatted lines: count / ceiling (band = 1% of count, rounded up)")
    for layer, ceiling in budgets.items():
        state = budget_state(counts[layer], ceiling, sum(budgets.values()) - budgets["total"] if layer == "total" else None)
        print(f"{counts[layer]:5} / {ceiling:5} (band {band(counts[layer]):3})  {layer:13} {state}")
        failed |= state != "ok"
    return int(failed)


def budget_self_check() -> None:
    """One executable check of both band boundaries, formatter-normalized and untracked-tree accounting,
    and the refusal of formatter-disable directives."""
    import contextlib
    import tempfile
    import json
    python_fixture = '''"""module doc"""; x = "# data"
# comment
class C:
    """class doc"""
    field = "kept"
    "attribute prose counts"
    def f(self):
        """sync doc"""
        return "kept"
async def f():
    """async doc"""
    return "kept"
x = """data

# literal data"""
'''
    if python_lines(python_fixture) != 11 or python_lines('x = f"""\n\n"""\n') != 3 or python_lines('"é"; x = 1\n') != 1:
        raise AssertionError("Python docstring spans or executable string data counted incorrectly")
    cpp_fixture = r'''// line comment
/* block
comment */ int x; // beside code
const auto s = u8R"tag(// literal

/* data */ \
)tag";
const char *q = "escaped \" // data";
auto c = '\''; auto n = 1'000;
#include /* note */ <path//header>
// spliced \
still a comment
/\
* spliced block */
#define X \
    42 /* note */
#if 0
int inactive;
#endif
'''
    if cpp_lines(cpp_fixture) != 13:
        raise AssertionError("C++ literals, headers, comments or continuations counted incorrectly")
    for source, expected in (
        ("#if __has_include(<missing/*header>)\nint a;\nint b;\n#endif\n// */\nint c;\n", 5),
        ("%:include <path/*header>\nint a;\nint b;\n// */\nint c;\n", 4),
    ):
        if cpp_lines(source) != expected:
            raise AssertionError("C++ header-name contents hid executable lines")
    for invalid in ('/* open', '"open', 'R"tag(open', 'R"bad delimiter(x)bad delimiter"', "`"):
        try:
            cpp_lines(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError(f"unrecognised C++ source escaped refusal: {invalid!r}")
    layers = ["contributions", "properties", "equilibrium", "regression", "surface", "bindings"]
    for count, ceiling, fails in ((1000, 1000, False), (1000, 1010, False), (1000, 1011, True), (1001, 1000, True)):
        if (budget_state(count, ceiling) != "ok") != fails:  # exactly one band of headroom passes
            raise AssertionError(f"band edge wrong at count {count}, ceiling {ceiling}")
    if not budget_state(1000, 1020, 1005).endswith("to 1005"):
        raise AssertionError("total suggestion exceeds the sum of the layer ceilings")
    with tempfile.TemporaryDirectory(prefix="epcsaft-budget-") as temporary:
        root = Path(temporary)
        (root / "engine/native/core").mkdir(parents=True)
        (root / "engine/src/epcsaft").mkdir(parents=True)
        (root / "engine/src/epcsaft/packed.py").write_text("x = 1; y = 2\n")
        limits = {"core": 102, "python": 2, "total": 104} | dict.fromkeys(layers, 0)
        (root / "ARCHITECTURE.yaml").write_text("code_budgets: " + json.dumps(limits))
        source = root / "engine/native/core/untracked.cpp"
        with contextlib.redirect_stdout(io.StringIO()):
            # Packed statements count as the formatter lays them out: two per line here, two lines each.
            for statements, fails in ((100, True), (101, False), (102, False), (103, True)):
                source.write_text("int a; int b;\n" * (statements // 2) + "int c;\n" * (statements % 2))
                if code_budgets(root) != int(fails):
                    raise AssertionError(f"budget band/accounting failure at {statements} statements")
            source.write_text("// clang-format off\nint a; int b;\n")
            try:
                code_budgets(root)
            except ValueError as error:
                if "untracked.cpp:1: formatter-disable" not in str(error):
                    raise
            else:
                raise AssertionError("a formatter-disable directive escaped the count")
            source.unlink()
            packed = root / "engine/src/epcsaft/packed.py"
            packed.write_text("x = 1; y = 2  # fmt:\fskip\n")
            try:
                code_budgets(root)
            except ValueError as error:
                if "packed.py:1: formatter-disable" not in str(error):
                    raise
            else:
                raise AssertionError("a form-feed formatter directive escaped the count")
            packed.write_text("x = 1; y = 2\n")
            source = root / "engine/native/unknown/new.cpp"
            source.parent.mkdir()
            source.write_text("int hidden;\n")
            try:
                code_budgets(root)
            except ValueError as error:
                if "unclassified native source" not in str(error):
                    raise
            else:
                raise AssertionError("unclassified source escaped the budget")


def main() -> int:
    if sys.argv[1:] != ["--budgets"]:
        raise SystemExit("usage: tools/diffstat.py --budgets")
    budget_self_check()
    return code_budgets()


if __name__ == "__main__":
    raise SystemExit(main())
