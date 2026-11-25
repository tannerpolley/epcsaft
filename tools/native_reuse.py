"""Index native result storage and operand readers from the production compile database."""

from __future__ import annotations

import ctypes as C
import ctypes.util
import json
import re
import shlex
import subprocess
import tempfile
from collections import Counter, defaultdict
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
    _fields_ = [
        ("filename", C.c_char_p),
        ("contents", C.c_char_p),
        ("length", C.c_ulong),
    ]


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
            (
                "parseTranslationUnit",
                P,
                [P, C.c_char_p, C.POINTER(C.c_char_p), C.c_int, P, C.c_uint, C.c_uint],
            ),
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
            (
                "getExpansionLocation",
                None,
                [Location, C.POINTER(P), C.POINTER(C.c_uint), P, P],
            ),
            ("getFileName", String, [P]),
            ("Location_isInSystemHeader", C.c_int, [Location]),
            ("Cursor_isNull", C.c_int, [Cursor]),
            ("Cursor_getStorageClass", C.c_uint, [Cursor]),
            ("getCursorTLSKind", C.c_uint, [Cursor]),
            ("isConstQualifiedType", C.c_uint, [Type]),
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


FUNCTION = {
    "FunctionDecl",
    "CXXMethod",
    "Constructor",
    "Destructor",
    "ConversionFunction",
    "FunctionTemplate",
}
RECORD = {
    "StructDecl",
    "ClassDecl",
    "UnionDecl",
    "ClassTemplate",
    "ClassTemplatePartialSpecialization",
}

# Numerical records that can be supplied separately from the calculation that
# produced them. C++ aliases and containers are resolved by record_types().
RESULTS = frozenset(
    [
        "AssociationAction",
        "AssociationJacobian",
        "CheckResult",
        "CompiledSystem",
        "DensityRoot",
        "DensityScanCache",
        "ElectrolyteActivity",
        "EquilibriumNumericalEstimates",
        "EquilibriumResult",
        "FitResult",
        "GibbsDifference",
        "GibbsDifferenceReference",
        "HelmholtzDerivativeValues",
        "InfiniteDilutionReference",
        "IsothermMoments",
        "NeutralReferenceReactionDiagnostic",
        "NeutralReferenceTrace",
        "ObservationEvaluation",
        "OptimizerState",
        "ParameterStandard",
        "PhasePointValues",
        "PointEvaluationValues",
        "Prediction",
        "PredictionGradient",
        "PreparedSuppliedPoints",
        "PreparedSuppliedTrials",
        "PureIsothermalValues",
        "ReferenceCache",
        "ReferenceCalorics",
        "ReferenceTaylor",
        "ReferenceTaylorBatch",
        "ReferenceTaylorPreparation",
        "SolvedStateActionWorkspace",
        "SolvedStateCaloricCenter",
        "SolventActivity",
        "SpinodalSearch",
        "StabilityMeasure",
        "StabilityReport",
        "StabilityTrial",
        "StandardTerm",
        "State",
    ]
)

# Objects reused across requests or callbacks. A field holding one of these
# owners, or a static holding a native record, extends that storage graph.
OWNERS = frozenset(
    [
        "CompiledSystem",
        "DensityScanCache",
        "EquilibriumContext",
        "GibbsMinimizer",
        "HelmholtzEvaluation",
        "HelmholtzModel",
        "HelmholtzTape",
        "ObservationEvaluator::Impl",
        "Prepared",
        "PropertyContext",
        "ReferenceCache",
        "ReferenceTaylorBatch",
        "ReferenceTaylorPreparation",
        "SolvedStateActionWorkspace",
        "SolvedStateCaloricCenter",
        "Trial",
    ]
)


def record_named(name, names):
    return name.rsplit("::", 1)[-1] in names or any(name.endswith("::" + n) for n in names if "::" in n)


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
        self.record_operands = defaultdict(set)
        self.storage = set()
        self.parameters = set()
        self.uses = set()
        self.callbacks = set()
        self.defaults = set()

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
        records = self.record_types(ty)
        if clang.kind(cursor) == "ParmDecl":
            # Dependent argument types expose their declaration through
            # TemplateRef even when getTypeDeclaration cannot resolve it.
            for child in clang.children(cursor):
                if clang.kind(child) in {"TemplateRef", "TypeRef"}:
                    declaration = clang.getCursorReferenced(child)
                    if clang.kind(declaration) in RECORD and self.location(declaration)[0]:
                        records.add(self.identity(declaration))
        return Operand(
            clang.spelling(cursor),
            clang.string(clang.getTypeSpelling(clang.getCanonicalType(ty))),
            frozenset(records),
        )

    def is_static(self, cursor):
        clang = self.clang
        parent = clang.getCursorSemanticParent(cursor)
        return (
            clang.kind(parent) in RECORD | {"Namespace", "TranslationUnit"}
            or clang.Cursor_getStorageClass(cursor) == 3
            or clang.getCursorTLSKind(cursor)
        )

    def callable_target(self, cursor, inherited):
        clang = self.clang
        kind = clang.kind(cursor)
        if kind == "ReturnStmt":
            return "return value"
        if kind in {"FieldDecl", "VarDecl"}:
            if kind == "FieldDecl" or self.is_static(cursor):
                return self.identity(cursor)
            return ""
        if kind in {"BinaryOperator", "CallExpr"}:
            children = clang.children(cursor)
            if (kind == "BinaryOperator" or clang.spelling(clang.getCursorReferenced(cursor)) == "operator=") and (
                len(children) >= 2 and clang.kind(children[0]) == "MemberRefExpr"
            ):
                return self.identity(clang.getCursorReferenced(children[0]))
        return inherited

    def raw_capture(self, children):
        clang = self.clang
        for child in children:
            if clang.kind(child) == "CompoundStmt":
                pending = [child]
                while pending:
                    expression = pending.pop()
                    if clang.kind(expression) == "CXXThisExpr":
                        return True
                    pending.extend(clang.children(expression))
                return False
            if clang.kind(child) == "CXXThisExpr":
                return True
            if clang.kind(child) in {"DeclRefExpr", "VariableRef"}:
                captured = clang.getCursorReferenced(child)
                if clang.getPointeeType(clang.getCursorType(captured)).kind:
                    return True
        return False

    def retained_default(self, condition, alternatives):
        clang = self.clang

        def descendants(cursor):
            yield cursor
            for child in clang.children(cursor):
                yield from descendants(child)

        def numeric_literal(cursor):
            while clang.kind(cursor) in {"UnexposedExpr", "ParenExpr", "UnaryOperator"}:
                children = clang.children(cursor)
                if len(children) != 1:
                    return False
                cursor = children[0]
            return clang.kind(cursor) in {"FloatingLiteral", "IntegerLiteral"}

        if not any(numeric_literal(alt) for alt in alternatives):
            return None
        for node in descendants(condition):
            if clang.kind(node) != "MemberRefExpr":
                continue
            ref = clang.getCursorReferenced(node)
            if not self.location(ref)[0] or clang.kind(ref) != "FieldDecl":
                continue
            owner = self.identity(clang.getCursorSemanticParent(ref))
            if record_named(owner, RESULTS | OWNERS) and any(record_named(name, RESULTS) for name in self.operand(ref).records):
                return self.identity(ref)
        return None

    def visit(self, cursor, function="", call="", target="", record=""):
        clang = self.clang
        kind = clang.kind(cursor)
        file, line = self.location(cursor)
        if not file and kind != "TranslationUnit":
            return
        if kind in RECORD:
            record = self.identity(cursor)
        if kind in FUNCTION:
            function = self.identity(cursor)
            call = target = ""
            parameters = [child for child in clang.children(cursor) if clang.kind(child) == "ParmDecl"]
            for position, child in enumerate(parameters, 1):
                arg_file, arg_line = self.location(child)
                self.parameters.add((function, position, self.operand(child), arg_file, arg_line))
        if kind in {"FieldDecl", "C++ base class specifier"} or (kind == "VarDecl" and self.is_static(cursor)):
            operand = self.operand(cursor)
            if kind == "VarDecl":
                if not clang.isConstQualifiedType(clang.getCursorType(cursor)):
                    self.storage.add(
                        (
                            "static",
                            self.identity(clang.getCursorSemanticParent(cursor)),
                            operand,
                            file,
                            line,
                        )
                    )
            else:
                self.record_operands[record].update(operand.records)
                self.storage.add(("member", record, operand, file, line))
        if kind == "CallExpr" and function:
            called = clang.getCursorReferenced(cursor)
            if self.location(called)[0]:
                call = self.identity(called)
        if kind == "MemberRefExpr" and function and call:
            ref = clang.getCursorReferenced(cursor)
            if clang.kind(ref) == "FieldDecl" and self.location(ref)[0]:
                operand = self.operand(ref)
                if operand.records:
                    source = self.identity(clang.getCursorSemanticParent(ref))
                    self.uses.add((function, source, operand, call, file, line))
        target = self.callable_target(cursor, target)
        children = clang.children(cursor)
        if kind == "ConditionalOperator" and function and len(children) == 3:
            field = self.retained_default(children[0], children[1:])
            if field:
                self.defaults.add(Site("default", function, f"{field} selects a numeric fallback", file, line))
        if kind == "LambdaExpr" and function and target:
            if self.raw_capture(children):
                self.callbacks.add(Site("callback", function, target, file, line))
            target = ""
        for child in children:
            self.visit(child, function, call, target, record)

    def sites(self):
        retained = {name for _, _, op, _, _ in self.storage for name in op.records if record_named(name, RESULTS)}
        owners = {name for name in self.record_operands if record_named(name, OWNERS)}
        owners.update(name for form, _, operand, _, _ in self.storage if form == "static" for name in operand.records)
        pending = list(owners)
        while pending:
            for record in self.record_operands.get(pending.pop(), ()):
                if record not in owners and record not in retained:
                    owners.add(record)
                    pending.append(record)

        def is_retained(operand):
            return any(record_named(name, RESULTS) for name in operand.records)

        sites = self.callbacks | self.defaults
        for form, owner, operand, file, line in self.storage:
            callable_member = operand.type.startswith("std::function<")
            if form == "static" or owner in owners or callable_member:
                sites.add(Site(form, owner, f"{operand.name}: {operand.type}", file, line))

        readers = set()
        for owner, position, operand, file, line in self.parameters:
            if is_retained(operand):
                readers.add(owner)
                sites.add(
                    Site(
                        "reader",
                        owner,
                        f"parameter {position}: {operand.type}",
                        file,
                        line,
                    )
                )
        for reader, owner, operand, call, file, line in self.uses:
            if call in readers and is_retained(operand) and (owner in owners or record_named(owner, RESULTS)):
                sites.add(
                    Site(
                        "forward",
                        reader,
                        f"{call} ← {owner}::{operand.name}",
                        file,
                        line,
                    )
                )
        by_key = {}
        for site in sorted(sites):
            by_key.setdefault(site.key, site)
        return list(by_key.values())


def production_commands(root, build):
    native = root / "engine/native"
    declared = set(
        re.findall(
            r"['\"](native/[^'\"\s]+\.cpp)['\"]",
            (root / "engine/meson.build").read_text(),
        )
    )
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
    unsaved = (UnsavedFile * len(overlays or {}))(
        *(UnsavedFile(str(path).encode(), contents, len(contents)) for path, contents in (overlays or {}).items())
    )
    index = clang.createIndex(0, 0)
    try:
        for directory, source, args in commands:
            inventory.directory = directory
            inventory.locations.clear()
            arguments = [
                f"-working-directory={directory}",
                f"-resource-dir={resource}",
                *args,
            ]
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


def inspect(root: Path, build: Path, revision: str | None = None) -> tuple[int, list[Site]]:
    commands = production_commands(root, build)
    overlays = {}
    if revision:
        revision = subprocess.check_output(["git", "rev-parse", "--verify", revision], cwd=root, text=True).strip()
        changed = subprocess.check_output(
            ["git", "diff", "--name-status", "--no-renames", "-z", revision, "--", "engine/native", "engine/meson.build"], cwd=root
        )
        fields = iter(changed.decode().split("\0"))
        for status, path in zip(fields, fields):
            if path == "engine/meson.build" or path.endswith(".cpp") and status in {"A", "D"}:
                raise RuntimeError("the requested revision needs a production compile database of its own")
            if status != "A":
                overlays[root / path] = subprocess.check_output(["git", "show", f"{revision}:{path}"], cwd=root)
    return len(commands), scan(root, commands, overlays)


def table(revision: str, units: int, sites: list[Site]) -> str:
    counts = Counter(site.kind for site in sites)
    lines = [
        "# Native retained-result inventory",
        "",
        f"Generated from production revision `{revision}` with Clang 18.",
        f"Inspected {units} translation units and flagged {len(sites)} entries.",
        "",
        "| Kind | Entries |",
        "| --- | ---: |",
        *(f"| {kind} | {count} |" for kind, count in sorted(counts.items())),
        "",
        "| Kind | Owner or reader | Operand | Location under engine/native |",
        "| --- | --- | --- | --- |",
    ]
    for site in sites:
        owner = site.owner.replace("|", r"\|")
        operand = site.operand.replace("|", r"\|")
        lines.append(f"| {site.kind} | `{owner}` | `{operand}` | `{site.file}:{site.line}` |")
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
            text = text.replace(
                "struct ObservationEvaluator::Impl final {",
                f"struct ObservationEvaluator::Impl final {{\n"
                f"    using Column = std::optional<Prediction>;\n"
                f"    Column {name};\n",
                1,
            )
            marker = '        core::Scope scope("fit_evaluation");'
            assert marker in text
            return text.replace(
                marker,
                f"        if ({name} && {name}->value) {{ out.predictions = {{*{name}->value}}; return true; }}\n{marker}",
                1,
            )

        sites = additions("fit.cpp", alias)
        assert any(
            site.kind == "member" and "ObservationEvaluator::Impl" in site.owner and site.operand.startswith(name + ":")
            for site in sites
        ), sites
        yield f"aliased std::optional<Prediction> in ObservationEvaluator::Impl::{name} flagged"

    def stale(text):
        current = "density_branches_at(compiled, prepared, anchors, initial, scope, limits, &roots)"
        assert current in text
        return text.replace(
            current,
            "density_branches_at(compiled, *result.prepared_endpoint, anchors, initial, scope, limits, &roots)",
            1,
        )

    sites = additions("phase_boundary.cpp", stale)
    assert any(
        site.kind == "forward"
        and "::check_density_branches " in site.owner
        and "::density_branches_at " in site.operand
        and "::prepared_endpoint" in site.operand
        for site in sites
    ), sites
    assert not any(site.kind == "reader" for site in sites), sites
    yield "stale PhasePoint forwarding flagged in check_density_branches with its declaration unchanged"

    def static_cache(text):
        marker = '        core::Scope scope("fit_evaluation");'
        assert marker in text
        return text.replace(
            marker,
            """        static thread_local std::optional<Prediction> saved_row;
        if (saved_row && saved_row->value) {
            out.predictions = {*saved_row->value};
            return true;
        }
"""
            + marker,
            1,
        )

    sites = additions("fit.cpp", static_cache)
    assert any(site.kind == "static" and site.operand.startswith("saved_row:") for site in sites), sites
    yield "thread-local static Prediction cache flagged"

    def raw_pointer(text):
        return (
            text
            + """
namespace epcsaft {
std::function<double()> retained_pressure(const State* state) {
    return [state] { return state->pressure; };
}
}
"""
        )

    sites = additions("phase_boundary.cpp", raw_pointer)
    assert any(site.kind == "callback" and "::retained_pressure " in site.owner for site in sites), sites
    yield "returned callable capturing a raw State pointer flagged"

    def refused_projection(text):
        return text + """
namespace epcsaft::regression {
double replay_activity(const Prediction& prediction) {
    return prediction.activity ? prediction.activity->mean_ionic_molality_scale : 1.0;
}
}
"""

    sites = additions("observations.cpp", refused_projection)
    assert any(site.kind == "default" and "::replay_activity " in site.owner and "::Prediction::activity" in site.operand
               for site in sites), sites
    yield "numeric fallback from an absent retained activity flagged"

    def shared_spinodal_output(text):
        marker = "    static auto& begin_diagnostics(const CompiledSystem& c)"
        assert marker in text
        plant = """    static void shared_output(const CompiledSystem& c, const SpinodalSearch& carried) {
        spinodal_search(c, nullptr, carried.seed, &carried, &*c.spinodal_diagnostics_);
    }
"""
        return text.replace(marker, plant + marker, 1)

    sites = additions("phase_boundary.cpp", shared_spinodal_output)
    assert any(site.kind == "forward" and "::shared_output " in site.owner and "::spinodal_diagnostics_" in site.operand
               for site in sites), sites
    yield "carried spinodal output forwarded into a compiled result slot flagged"
