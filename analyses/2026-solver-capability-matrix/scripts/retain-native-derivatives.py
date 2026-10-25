"""Capture existing E5 operands at the named release; no thermodynamic calculation."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

RELEASE = "8298c6029ebb8b467929badc4a223c580ec591ab"
source = Path(sys.argv[1]).resolve()
out = Path(__file__).resolve().parents[1] / "results/main-8298c602-higher-order"
out.mkdir(parents=True, exist_ok=True)
(out / "summary.json").unlink(missing_ok=True)
(out / "failure.json").unlink(missing_ok=True)
env = {k: v for k, v in os.environ.items() if k not in ("LD_LIBRARY_PATH", "LD_PRELOAD", "EPCSAFT_HSL_LIBRARY")}
def run(*args):
    return subprocess.run(args, cwd=source, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def identity():
    assert run("git", "rev-parse", "HEAD").stdout.strip() == RELEASE, "wrong release"
    dirty = run("git", "status", "--porcelain", "--untracked-files=all", "--", "engine", "tools", ".codex", "uv.lock", "pyproject.toml")
    assert dirty.returncode == 0 and not dirty.stdout, dirty.stdout
    return {name: run("git", "rev-parse", f"HEAD:{name}").stdout.strip() for name in ("engine/native", "engine/tests/native", "tools/native-run")}
try:
    before = identity()
    build = Path((source / "build/environment-wheel/meson-build-dir").read_text().strip()).resolve()
    info = json.loads((build / "meson-info/meson-info.json").read_text())
    assert Path(info["directories"]["source"]).resolve() == source / "engine" and Path(info["directories"]["build"]).resolve() == build
    command = ["nice", "-n", "10", "tools/native-run", "equilibrium_checks", "--list"]
    built = run(*command)
    (out / "build-list.log").write_text(built.stdout)
    assert built.returncode == 0, "native build/list failed"
    executable, library = build / "equilibrium_checks", build / "libepcsaft_native.so"
    linked = run("ldd", str(executable))
    (out / "dependencies.log").write_text(linked.stdout)
    assert linked.returncode == 0 and Path(re.search(r"libepcsaft_native\.so => (\S+)", linked.stdout)[1]).resolve() == library.resolve()
    products = {str(p): sha(p) for p in (executable, library)}
    captures = {}
    for case in ("e5-neutral", "e5-electrolyte", "e5-reactive"):
        result = run("nice", "-n", "10", str(executable), "native.equilibrium." + case)
        (out / f"{case}.log").write_text(result.stdout)
        captures[case] = {"returncode": result.returncode, "passed": f"test_id=native.equilibrium.{case} execution=Passed" in result.stdout}
    assert identity() == before and products == {p: sha(Path(p)) for p in products}, "source or products changed"
    assert all(r["returncode"] == 0 and r["passed"] for r in captures.values()), captures
    phase = ["amount", "volume", "molar_density", "packing_fraction", "pressure"] + [f"{q}[{i}]" for i in range(3) for q in ("amount", "mole_fraction", "residual_chemical_potential_over_rt", "log_activity")]
    neutral = ["temperature", "pressure"] + [f"phase[{p}].{q}" for p in range(2) for q in phase]
    caloric = ["phase[0].property/" + q for q in ("ResidualEnthalpy", "ResidualIsobaricHeatCapacity", "TotalEnthalpy", "TotalIsobaricHeatCapacity")]
    fixtures = {"e5-1/neutral-ternary": neutral, "e5-5/electrolyte": caloric + [f"phase[0].log_activity[{i}]" for i in range(3)], "e5-5/electrolyte-residual-only": caloric[:2]}
    expected = {f"{f}/{q}/A{k}" for f, quantities in fixtures.items() for q in quantities for k in (2, 3)}
    analytic = {f"e5-3/{f}{s}/A{k}" for f in ("isomerisation", "isomerisation-two-mol") for s in ("", "/rescaled") for k in (2, 3)}
    records, rows = {}, []
    for case in captures:
        for line in (out / f"{case}.log").read_text().splitlines():
            if not line.startswith("EPCSAFT_COMPARE_JSON "):
                continue
            r = json.loads(line.split(" ", 1)[1])
            assert r.get("verdict", True) is True, r
            key = r["comparison_id"]
            assert key not in records, "duplicate comparison: " + key
            records[key] = dict(r, test_id="native.equilibrium." + case)
    selected = {key for key in records if re.search(r"/A[23]$", key) and any(key.startswith(f + "/") for f in (*fixtures, "e5-3/isomerisation", "e5-3/isomerisation-two-mol"))}
    assert selected == expected | analytic, {"missing": sorted((expected | analytic) - selected), "extra": sorted(selected - expected - analytic)}
    for key in sorted(selected):
        r, order = records[key], int(key[-1])
        v = r["inputs"]
        assert r["quantity"] == "action" and v["status"] == "available", r
        companion = records[key + "/rate"] if key in expected else None
        c = v["allowance"] if companion else r["atol"]
        operands = [r["model"], r["reference"], c] + ([v["term_scale"], v["reference_uncertainty"], companion["model"], companion["reference"], *v["ladder"], *v["reference_steps"]] if companion else [v["amount_action"], v["total"], v["bar"]])
        assert all(math.isfinite(x) for x in operands) and c > 0, r
        assert not companion or (companion["verdict"] is True and v["reference_uncertainty"] <= c), r
        for observable in ([key.rsplit("/A", 1)[0].split("/", 2)[2]] if companion else ["x_B", "n_B"]):
            amount = observable == "n_B"
            action, reference, allowance = (v["amount_action"], v["total"] * r["reference"], v["total"] * c) if amount else (r["model"], r["reference"], c)
            error = abs(action - reference)
            rule = "absolute_inclusive" if amount else r["rule"]
            assert error < allowance if rule == "strict_absolute" else error <= allowance, (key, observable)
            rows.append(dict(case=r["test_id"].rsplit(".", 1)[1], order=order, comparison_id=key, observable=observable, action=action, reference=reference, allowance=allowance, absolute_error=error, ratio=error / allowance, unit=(f"mol/K^{order}" if amount else f"1/K^{order}") if not companion else r["unit"], rule=rule, method=r["method"], native_record=json.dumps(r), refinement_record=json.dumps(companion)))
    assert len(rows) == 106 and sum(key + "/rate" in records for key in expected) == 90
    with (out / "comparisons.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    groups = {f"{case}/A{k}": [r for r in rows if r["case"] == case and r["order"] == k] for case in captures for k in (2, 3)}
    groups.update({f"pooled/A{k}": [r for r in rows if r["order"] == k] for k in (2, 3)})
    summary = dict(release="main-" + RELEASE, revision=RELEASE, evidence_type="numerical verification", source=str(source), build=str(build), source_identity=before, products_sha256=products, build_command=command, build_jobs=int((source / "build/environment-wheel/build-jobs").read_text()), compilers=json.loads((build / "meson-info/intro-compilers.json").read_text()), build_options=json.loads((build / "meson-info/intro-buildoptions.json").read_text()), captures=captures, row_count=len(rows), refinement_count=90, maxima={g: dict(count=len(rs), **max(rs, key=lambda r: r["ratio"])) for g, rs in groups.items()}, hashes={p.name: sha(p) for p in out.glob("*.log")})
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
except Exception as error:
    (out / "failure.json").write_text(json.dumps({"error": repr(error), "captures": globals().get("captures", {})}, indent=2) + "\n")
    raise
