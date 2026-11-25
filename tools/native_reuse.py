"""Index native result storage and operand readers from the production compile database."""

from __future__ import annotations

import ctypes as C
import ctypes.util
import json
import re
import shlex
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, order=True)
class Site:
    kind: str
    owner: str
    operand: str
    file: str
    line: int

    @property
    def key(self):
        return self.kind, self.owner, self.operand


class Cursor(C.Structure):
    _fields_ = [("kind", C.c_uint), ("data_index", C.c_int), ("data", C.c_void_p * 3)]


class Type(C.Structure):
    _fields_ = [("kind", C.c_uint), ("data", C.c_void_p * 2)]


class String(C.Structure):
    _fields_ = [("data", C.c_void_p), ("flags", C.c_uint)]


class Location(C.Structure):
    _fields_ = [("data", C.c_void_p * 2), ("integer", C.c_uint)]


class UnsavedFile(C.Structure):
    _fields_ = [("filename", C.c_char_p), ("contents", C.c_char_p), ("length", C.c_ulong)]


P = C.c_void_p
Visit = C.CFUNCTYPE(C.c_uint, Cursor, Cursor, P)


class Clang:
    def __init__(self):
        library = ctypes.util.find_library("clang-18")
        if library is None:
            raise RuntimeError("native reuse inventory requires libclang 18")
        self.library = C.CDLL(library)
        for name, result, args in (
            ("createIndex", P, [C.c_int, C.c_int]),
            ("disposeIndex", None, [P]),
            ("parseTranslationUnit", P, [P, C.c_char_p, C.POINTER(C.c_char_p), C.c_int, P, C.c_uint, C.c_uint]),
            ("disposeTranslationUnit", None, [P]),
            ("getNumDiagnostics", C.c_uint, [P]),
            ("getDiagnostic", P, [P, C.c_uint]),
            ("getDiagnosticSeverity", C.c_uint, [P]),
            ("formatDiagnostic", String, [P, C.c_uint]),
            ("disposeDiagnostic", None, [P]),
            ("getTranslationUnitCursor", Cursor, [P]),
            ("visitChildren", C.c_uint, [Cursor, Visit, P]),
            ("getCursorKindSpelling", String, [C.c_uint]),
            ("getCursorSpelling", String, [Cursor]),
            ("getCursorType", Type, [Cursor]),
            ("getCursorResultType", Type, [Cursor]),
            ("getCanonicalType", Type, [Type]),
            ("getTypeSpelling", String, [Type]),
            ("getTypeDeclaration", Cursor, [Type]),
            ("getPointeeType", Type, [Type]),
            ("getArrayElementType", Type, [Type]),
            ("Type_getNumTemplateArguments", C.c_int, [Type]),
            ("Type_getTemplateArgumentAsType", Type, [Type, C.c_uint]),
            ("getSpecializedCursorTemplate", Cursor, [Cursor]),
            ("getCursorSemanticParent", Cursor, [Cursor]),
            ("getCursorReferenced", Cursor, [Cursor]),
            ("getCursorLocation", Location, [Cursor]),
            ("getExpansionLocation", None, [Location, C.POINTER(P), C.POINTER(C.c_uint), P, P]),
            ("getFileName", String, [P]),
            ("Location_isInSystemHeader", C.c_int, [Location]),
            ("Cursor_isNull", C.c_int, [Cursor]),
            ("Cursor_getStorageClass", C.c_uint, [Cursor]),
            ("getCString", C.c_char_p, [String]),
            ("disposeString", None, [String]),
        ):
            fn = getattr(self.library, "clang_" + name)
            fn.restype = result
            fn.argtypes = args
            setattr(self, name, fn)
        self.kinds = {}

    def string(self, value):
        try:
            return (self.getCString(value) or b"").decode()
        finally:
            self.disposeString(value)

    def kind(self, cursor):
        if cursor.kind not in self.kinds:
            self.kinds[cursor.kind] = self.string(self.getCursorKindSpelling(cursor.kind))
        return self.kinds[cursor.kind]

    def spelling(self, cursor):
        return self.string(self.getCursorSpelling(cursor))

    def location(self, cursor):
        location = self.getCursorLocation(cursor)
        if self.Location_isInSystemHeader(location):
            return "", 0
        file, line = P(), C.c_uint()
        self.getExpansionLocation(location, C.byref(file), C.byref(line), None, None)
        return (self.string(self.getFileName(file)), line.value) if file else ("", 0)

    def children(self, cursor):
        children = []

        @Visit
        def collect(child, parent, data):
            children.append(child)
            return 1  # CXChildVisit_Continue

        self.visitChildren(cursor, collect, None)
        return children


FUNCTION = {"FunctionDecl", "CXXMethod", "Constructor", "Destructor", "ConversionFunction", "FunctionTemplate"}
RECORD = {"StructDecl", "ClassDecl", "UnionDecl", "ClassTemplate", "ClassTemplatePartialSpecialization"}


@dataclass(frozen=True)
class Operand:
    name: str
    type: str
    records: frozenset[str]


class Inventory:
    def __init__(self, clang, native):
        self.clang = clang
        self.native = str(native) + "/"
        self.directory = native.parent
        self.locations = {}
        self.results = set()
        self.record_operands = defaultdict(set)
        self.native_parameters = set()
        self.storage = set()
        self.parameters = set()
        self.uses = set()

    def location(self, cursor):
        path, line = self.clang.location(cursor)
        if not path:
            return "", line
        if path not in self.locations:
            self.locations[path] = str((self.directory / path).resolve()).removeprefix(self.native)
        relative = self.locations[path]
        return ("", line) if relative.startswith("/") else (relative, line)

    def identity(self, cursor):
        clang = self.clang
        parts = []
        while not clang.Cursor_isNull(cursor) and clang.kind(cursor) != "TranslationUnit":
            name = clang.spelling(cursor)
            if not name:
                name = "(" + self.location(cursor)[0] + ")"
            if clang.kind(cursor) in FUNCTION:
                signature = clang.string(clang.getTypeSpelling(clang.getCanonicalType(clang.getCursorType(cursor))))
                name += " [" + signature + "]"
            parts.append(name)
            cursor = clang.getCursorSemanticParent(cursor)
        return "::".join(reversed(parts))

    def record_types(self, ty, seen=None):
        """Unwrap aliases, arrays, pointers and type arguments using Clang's type graph."""
        if not ty.kind:
            return set()
        clang = self.clang
        ty = clang.getCanonicalType(ty)
        key = ty.kind, *ty.data
        if seen is None:
            seen = set()
        if key in seen:
            return set()
        seen.add(key)
        records = set()
        decl = clang.getTypeDeclaration(ty)
        if clang.kind(decl) in RECORD and self.location(decl)[0]:
            template = clang.getSpecializedCursorTemplate(decl)
            records.add(self.identity(decl if clang.Cursor_isNull(template) else template))
        records.update(self.record_types(clang.getPointeeType(ty), seen))
        records.update(self.record_types(clang.getArrayElementType(ty), seen))
        for i in range(clang.Type_getNumTemplateArguments(ty)):
            records.update(self.record_types(clang.Type_getTemplateArgumentAsType(ty, i), seen))
        return records

    def operand(self, cursor):
        clang = self.clang
        ty = clang.getCursorType(cursor)
        return Operand(clang.spelling(cursor), clang.string(clang.getTypeSpelling(clang.getCanonicalType(ty))),
                       frozenset(self.record_types(ty)))

    def visit(self, cursor, function="", call=""):
        clang = self.clang
        kind = clang.kind(cursor)
        file, line = self.location(cursor)
        if not file and kind != "TranslationUnit":
            return
        if kind in FUNCTION:
            function = self.identity(cursor)
            ty = clang.getCursorResultType(cursor)
            if not clang.getPointeeType(ty).kind:
                self.results.update(self.record_types(ty))
            call = ""
        if kind in {"FieldDecl", "C++ base class specifier", "VarDecl"}:
            parent = clang.getCursorSemanticParent(cursor)
            parent_kind = clang.kind(parent)
            persistent = parent_kind in {"Namespace", "TranslationUnit"} or clang.Cursor_getStorageClass(cursor) == 3
            operand = self.operand(cursor)
            if kind != "VarDecl" or persistent or operand.type.startswith("std::map<"):
                owner = self.identity(parent)
                if parent_kind in RECORD:
                    self.record_operands[owner].update(operand.records)
                form = "store" if parent_kind in RECORD else "static" if persistent else "local"
                self.storage.add((form, owner, operand, file, line))
        if kind == "ParmDecl" and function:
            operand = self.operand(cursor)
            self.parameters.add((function, operand, file, line))
            if operand.type.endswith((" *", " &")):
                for child in clang.children(cursor):
                    if clang.kind(child) in {"TemplateRef", "TypeRef"}:
                        declaration = clang.getCursorReferenced(child)
                        if clang.kind(declaration) in RECORD and self.location(declaration)[0]:
                            self.native_parameters.add((function, operand, self.identity(declaration)))
        if kind == "CallExpr" and function:
            called = clang.getCursorReferenced(cursor)
            if self.location(called)[0]:
                call = self.identity(called)
        if kind in {"DeclRefExpr", "MemberRefExpr"} and function:
            ref = clang.getCursorReferenced(cursor)
            if clang.kind(ref) in {"ParmDecl", "FieldDecl", "VarDecl"} and self.location(ref)[0] and (call or kind == "MemberRefExpr"):
                operand = self.operand(ref)
                if operand.records or kind == "MemberRefExpr":
                    origin = self.identity(clang.getCursorSemanticParent(ref))
                    self.uses.add((function, origin, operand, call, file, line))
        for child in clang.children(cursor):
            self.visit(child, function, call)

    def sites(self):
        # A native returned record remains a result when buried under an alias, a
        # container, a base or a preparation/owner that holds it by pointer.
        retained = set(self.results)
        while True:
            owners = {owner for owner, operands in self.record_operands.items() if operands & retained}
            if owners <= retained:
                break
            retained.update(owners)
        sites = set()
        stored_types = set().union(*self.record_operands.values())
        stored_parameters = {(reader, operand) for reader, operand, record in self.native_parameters if record in stored_types}
        for form, owner, operand, file, line in self.storage:
            if operand.records & retained or operand.type.startswith("std::function<") or (
                form == "local" and operand.type.startswith("std::map<std::vector<double>")
            ):
                sites.add(Site(form, owner, f"{operand.name}: {operand.type}", file, line))
        for owner, operand, file, line in self.parameters:
            numerical_input = operand.type.startswith(("const std::vector<", "const std::array<")) and operand.type.endswith(" *")
            if operand.records & retained or numerical_input or (owner, operand) in stored_parameters:
                sites.add(Site("reader", owner, f"{operand.name}: {operand.type}", file, line))
        for reader, owner, operand, call, file, line in self.uses:
            if operand.records & retained or (owner, operand) in stored_parameters:
                edge = f"{call} ← " if call else ""
                sites.add(Site("forward" if call else "read", reader, edge + f"{owner}::{operand.name}", file, line))
        by_key = {}
        for site in sorted(sites):
            by_key.setdefault(site.key, site)
        return list(by_key.values())


def production_commands(root, build):
    native = root / "engine/native"
    declared = set(re.findall(r"['\"](native/[^'\"\s]+\.cpp)['\"]", (root / "engine/meson.build").read_text()))
    sources = {str(path.relative_to(root / "engine")) for path in native.rglob("*.cpp")}
    if sources != declared:
        raise RuntimeError(f"native sources and Meson disagree: {sorted(sources ^ declared)}")
    commands = []
    seen = set()
    for entry in json.loads((build / "compile_commands.json").read_text()):
        directory = Path(entry["directory"])
        source = (directory / entry["file"]).resolve()
        if not source.is_relative_to(native):
            continue
        relative = str(source.relative_to(root / "engine"))
        seen.add(relative)
        args = list(entry["arguments"]) if "arguments" in entry else shlex.split(entry["command"])
        if Path(args[0]).name in {"ccache", "sccache"}:
            args.pop(0)
        args.pop(0)  # libclang supplies its own compiler driver
        filtered = iter(args)
        compile_args = []
        for arg in filtered:
            if arg in {"-o", "-MF", "-MQ", "-MT"}:
                next(filtered)
            elif arg not in {"-c", "-MD", "-MMD"} and arg != entry["file"]:
                compile_args.append(arg)
        commands.append((directory, source, tuple(compile_args)))
    if seen != declared:
        raise RuntimeError(f"production compile database coverage differs from Meson: {sorted(seen ^ declared)}")
    return sorted(set(commands))


def scan(root, commands, overlays=None):
    clang = Clang()
    resource = subprocess.check_output(["clang++-18", "-print-resource-dir"], text=True).strip()
    inventory = Inventory(clang, root / "engine/native")
    unsaved = (UnsavedFile * len(overlays or {}))(*(
        UnsavedFile(str(path).encode(), contents, len(contents)) for path, contents in (overlays or {}).items()
    ))
    index = clang.createIndex(0, 0)
    try:
        for directory, source, args in commands:
            inventory.directory = directory
            inventory.locations.clear()
            arguments = [f"-working-directory={directory}", f"-resource-dir={resource}", *args]
            argv = (C.c_char_p * len(arguments))(*(arg.encode() for arg in arguments))
            unit = clang.parseTranslationUnit(index, str(source).encode(), argv, len(argv), unsaved, len(unsaved), 0)
            if not unit:
                raise RuntimeError(f"Clang could not inspect {source.relative_to(root)}")
            try:
                errors = []
                for i in range(clang.getNumDiagnostics(unit)):
                    diagnostic = clang.getDiagnostic(unit, i)
                    if clang.getDiagnosticSeverity(diagnostic) >= 3:  # CXDiagnostic_Error
                        errors.append(clang.string(clang.formatDiagnostic(diagnostic, 0)))
                    clang.disposeDiagnostic(diagnostic)
                if errors:
                    raise RuntimeError(f"Clang could not inspect {source.relative_to(root)}: " + "\n".join(errors))
                inventory.visit(clang.getTranslationUnitCursor(unit))
            finally:
                clang.disposeTranslationUnit(unit)
    finally:
        clang.disposeIndex(index)
    return inventory.sites()


def inspect(root: Path, build: Path) -> tuple[int, list[Site]]:
    commands = production_commands(root, build)
    return len(commands), scan(root, commands)


def table(revision: str, units: int, sites: list[Site]) -> str:
    groups = defaultdict(list)
    for site in sites:
        groups[site.kind, site.owner].append(site)
    lines = [
        "# Native retained-result inventory",
        "",
        f"Generated from production revision `{revision}` with Clang 18.",
        f"Inspected {units} translation units and flagged {len(sites)} operand sites at {len(groups)} boundaries.",
        "",
        "| Kind | Owner or reader | Operand sites | Locations under engine/native |",
        "| --- | --- | ---: | --- |",
    ]
    for (kind, owner), group in groups.items():
        owner = owner.replace("|", r"\|")
        locations = ", ".join(f"`{file}:{line}`" for file, line in sorted({(site.file, site.line) for site in group}))
        lines.append(f"| {kind} | `{owner}` | {len(group)} | {locations} |")
    return "\n".join(lines) + "\n"


def self_test(root, build):
    commands = production_commands(root, build)

    entries = json.loads((build / "compile_commands.json").read_text())
    without_asl = [entry for entry in entries if Path(entry["file"]).name != "asl_adapter.cpp"]
    with tempfile.TemporaryDirectory() as temporary:
        incomplete = Path(temporary)
        (incomplete / "compile_commands.json").write_text(json.dumps(without_asl))
        try:
            production_commands(root, incomplete)
        except RuntimeError as error:
            assert "native/bindings/asl_adapter.cpp" in str(error), error
        else:
            raise AssertionError("omitted native ASL translation unit was accepted")
    yield "compile database missing the production ASL translation unit refused"

    def additions(file, mutate):
        selected = [command for command in commands if command[1].name == file]
        assert selected, f"missing mutation translation unit {file}"
        source = selected[0][1]
        before = scan(root, selected)
        after = scan(root, selected, {source: mutate(source.read_text()).encode()})
        baseline = {site.key for site in before}
        return [site for site in after if site.key not in baseline]

    for name in ("column_buffer", "trial_dimension"):
        def alias(text, name=name):
            text = text.replace("struct ObservationEvaluator::Impl final {",
                f"struct ObservationEvaluator::Impl final {{\n"
                f"    using Column = std::optional<Prediction>;\n"
                f"    Column {name};\n", 1)
            marker = '        core::Scope scope("fit_evaluation");'
            assert marker in text
            return text.replace(marker,
                f"        if ({name} && {name}->value) {{ out.predictions = {{*{name}->value}}; return true; }}\n{marker}", 1)

        sites = additions("fit.cpp", alias)
        assert any(site.kind == "store" and "ObservationEvaluator::Impl" in site.owner
                   and site.operand.startswith(name + ":") for site in sites), sites
        yield f"aliased std::optional<Prediction> in ObservationEvaluator::Impl::{name} flagged"

    def stale(text):
        current = "density_branches_at(compiled, prepared, anchors, initial, scope, limits, &roots)"
        assert current in text
        return text.replace(current,
            "density_branches_at(compiled, *result.prepared_endpoint, anchors, initial, scope, limits, &roots)", 1)

    sites = additions("phase_boundary.cpp", stale)
    assert any(site.kind == "forward" and "::check_density_branches " in site.owner
               and "::density_branches_at " in site.operand and "::prepared_endpoint" in site.operand
               for site in sites), sites
    assert not any(site.kind == "reader" for site in sites), sites
    yield "stale PhasePoint forwarding flagged in check_density_branches with its declaration unchanged"
