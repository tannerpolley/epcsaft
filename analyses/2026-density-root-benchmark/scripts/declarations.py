"""Reuse the retained declarations; all thermodynamic solves stay in the Engine."""

import ast
import importlib.util
import json
import os
import sys
from importlib.machinery import SourceFileLoader
from pathlib import Path

# Never let analysis imports compile an editable Engine.
for finder in sys.meta_path:
    if getattr(finder, "_name", None) == "epcsaft" and finder._work_to_do(os.environ.copy()):
        raise RuntimeError("Engine import would rebuild; prepare the selected environment first")


ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]


def load(path, name):
    loader = SourceFileLoader(name, str(ROOT / path))
    spec = importlib.util.spec_from_loader(name, loader)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


def cross_check():
    module = load("analyses/2026-competitor-cross-check/scripts/engine.py", "density_cross_check")
    inputs = json.loads((ROOT / "analyses/2026-competitor-cross-check/inputs/cases.json").read_text())
    return module, inputs


def key_solves():
    """Call time-key-solves' declaration section without its response pipe or solves."""
    tool = load("tools/time-key-solves", "density_key_solves")
    tree = ast.parse((ROOT / "tools/time-key-solves").read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "worker")
    body = []
    for node in function.body[2:]:
        if isinstance(node, ast.For):
            break
        body.append(node)
    body += ast.parse(
        "return {**declarations, 'khudaida/sequential-8':[(model,problem) for _,model,problem in sequence], 'mea-nominally-unloaded/353.15-K':unloaded}"
    ).body
    function.name, function.body = "declaration", body
    exec(  # noqa: S102 -- reuse repository-owned input declarations
        compile(
            ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])),
            str(ROOT / "tools/time-key-solves"),
            "exec",
        ),
        tool.__dict__,
    )
    return tool, tool.declaration()


def fit_declaration():
    """Extract the replay's setup before its fit call, as time-key-solves does for inputs."""
    path = "analyses/2026-pure-saturation-regression/scripts/replay_propane_fit.py"
    module = load(path, "density_propane_fit")
    tree = ast.parse((ROOT / path).read_text())
    function = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "main")
    body = []
    for node in function.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "results" for t in node.targets):
            break
        body.append(node)
    body += ast.parse("return document, coordinates, observations(document, training), optimizer_starts, controls").body
    function.name, function.body, function.returns = "declaration", body, None
    exec(  # noqa: S102 -- execute only the repository-owned declaration AST
        compile(
            ast.fix_missing_locations(ast.Module(body=[function], type_ignores=[])),
            str(ROOT / path),
            "exec",
        ),
        module.__dict__,
    )
    return module.declaration()


def models():
    from epcsaft import Parameters, Pure

    module, inputs = cross_check()
    requests = {q["id"]: q for q in inputs["requests"]}
    result = {
        cid: module.prepare(inputs, requests[rid], declare_only=True)[0]
        for cid, rid in [
            ("C1", "C1-300"),
            ("C2", "C2-0.5"),
            ("C3", "C3-sigma-1"),
            ("C4", "C4-298.15"),
            ("C5", "C5-300"),
        ]
    }
    # Water is the same cross-check packet's neutral component, without salt.
    result["water"] = Pure(
        Parameters.from_mapping(
            module.document(inputs["models"]["C3"], requests["C3-sigma-1"]),
            components=["water"],
        )
    )
    bubble = load(
        "analyses/2026-binary-methane-ethane-bubble-point/scripts/collect_results.py",
        "density_bubble",
    )
    result["methane-ethane"] = bubble.collect_engine(declare_only=True)[0]
    _, key = key_solves()
    result["MEA"] = key["case-7/flash-5918.5286-Pa"][0]
    return result


def declarations():
    return models(), key_solves(), fit_declaration()
