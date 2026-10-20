import re
import subprocess
from pathlib import Path
from urllib.parse import quote, urlsplit

project = "ePC-SAFT"
REPO = Path(__file__).resolve().parents[2]
COMMIT = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
TRACKED = set(subprocess.check_output(["git", "ls-files"], cwd=REPO, text=True).splitlines())
extensions = [
    "autoapi.extension",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "myst_nb",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
]
html_theme = "pydata_sphinx_theme"
autoapi_dirs = [str(REPO / "engine/src/epcsaft")]
autoapi_file_patterns = ["*.pyi", "*.py"]
autoapi_generate_api_docs = False
autoapi_options = ["members", "undoc-members", "private-members", "show-inheritance"]
napoleon_google_docstring = False
napoleon_numpy_docstring = True
nb_execution_mode = "force"
nb_execution_allow_errors = False
nb_execution_raise_on_error = True
nb_execution_timeout = 120
nb_execution_cache_path = str(REPO / "build/docs/jupyter-cache")
exclude_patterns = ["**/.ipynb_checkpoints"]
html_title = "ePC-SAFT"
html_theme_options = {"navigation_depth": 3, "show_toc_level": 2}
intersphinx_mapping = {
    "python": ("https://docs.python.org/3.13", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pyomo": ("https://pyomo.readthedocs.io/en/stable", None),
    "idaes": ("https://idaes-pse.readthedocs.io/en/stable", None),
}

napoleon_use_ivar = True

myst_enable_extensions = ["dollarmath", "amsmath"]
myst_heading_anchors = 4

THEORY = {
    REPO / "engine/docs/equations.md": "equations.md",
    REPO / "engine/docs/science/algorithms.md": "algorithms.md",
    REPO / "engine/docs/science/full-eos.md": "eos.md",
    REPO / "engine/docs/equilibrium/formulation.md": "equilibrium.md",
    REPO
    / "engine/docs/regression/science/pure-saturation-regression.md": "parameter-estimation.md",
}


def include_source(app, relative_path, parent_docname, content):
    """Adapt canonical theory for MyST without changing fenced code or source files."""
    origin = (Path(app.srcdir) / relative_path).resolve()
    if origin not in THEORY:
        return

    def link(match):
        label, target = match.groups()
        if urlsplit(target).scheme or target.startswith("#"):
            return match[0]
        path, separator, fragment = target.partition("#")
        resolved = (origin.parent / path).resolve()
        if resolved in THEORY:
            url = THEORY[resolved]
        else:
            relative = resolved.relative_to(REPO).as_posix()
            tracked = relative in TRACKED or any(p.startswith(relative + "/") for p in TRACKED)
            if not tracked:
                issue = re.search(r"issue-(\d+)", relative)
                number = issue[1] if issue else "695"
                return f"[{label} — retained local evidence: `{relative}`](https://github.com/tannerpolley/ePC-SAFT/issues/{number})"
            kind = "tree" if resolved.is_dir() else "blob"
            url = f"https://github.com/tannerpolley/ePC-SAFT/{kind}/{COMMIT}/" + quote(relative)
        return f"[{label}]({url}{separator}{fragment})"

    def prose(text):
        text = re.sub(r'<a id="([^"]+)"></a>', lambda m: "\n(" + m[1] + ")=\n", text)
        text = text.replace(r"\[", "\n$$\n").replace(r"\]", "\n$$\n")
        text = text.replace(r"\(", "$").replace(r"\)", "$")
        return re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link, text)

    output, pending, fence = [], [], None
    for line in content[0].splitlines(keepends=True):
        marker = re.match(r"^ {0,3}(`{3,}|~{3,})", line)
        if fence:
            output.append(line)
            if (
                marker
                and marker[1][0] == fence[0]
                and len(marker[1]) >= len(fence)
                and not line[marker.end() :].strip()
            ):
                fence = None
            continue
        if marker:
            output.append(prose("".join(pending)))
            pending = []
            fence = marker[1]
            output.append(line)
        else:
            pending.append(line)
    output.append(prose("".join(pending)))
    content[0] = "".join(output)


def setup(app):
    app.connect("include-read", include_source)
