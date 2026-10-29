"""Admit the frozen hydrogen extension only after complete numerical replay."""

import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tomllib

import pandas as pd
from epcsaft import Parameters

ROOT = Path(__file__).resolve().parents[3]
LEAF = Path(__file__).resolve().parents[1]
R = LEAF / "results/packet-readback-r2"
FROZEN = LEAF / "results/issue-901-profile-r2"
BASE = ROOT / "data/packets/figiel-2025-chloride-mixture-refit/2"
EVIDENCE = ROOT.parents[1] / "build/issue-900/build-r2"
DECL = LEAF / "inputs/fit-declaration-issue-901.json"


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


def read(path):
    return pd.read_csv(path, float_precision="round_trip").fillna("")


def save(name, value):
    (R / name).write_text(json.dumps(value, indent=2) + "\n")
    shutil.copy2(R / name, EVIDENCE / name)


def equal(a, b, tolerance=False):
    assert type(a) is type(b), (type(a), type(b))
    if isinstance(a, dict):
        assert a.keys() == b.keys()
        for key in a:
            equal(a[key], b[key], tolerance)
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b)
        for x, y in zip(a, b, strict=True):
            equal(x, y, tolerance)
    elif isinstance(a, float):
        assert (
            abs(a - b) <= max(1e-12, 1e-10 * abs(a))
            if tolerance
            else struct.pack("d", a) == struct.pack("d", b)
        ), (a, b)
    else:
        assert a == b, (a, b)


def execute(destination, mode):
    destination.mkdir(parents=True)
    subprocess.run(
        [
            str(ROOT / "build/native-drivers/activity/driver"),
            str(LEAF),
            str(destination),
            mode,
        ],
        env={**os.environ, "LD_LIBRARY_PATH": str(ROOT / "build/engine-meson/shared")},
        check=True,
    )


def replay():
    subprocess.run(
        [
            sys.executable,
            str(LEAF / "scripts/run.py"),
            str(DECL),
            str(R / "transport"),
            "--packet",
            str(STAGED),
        ],
        check=True,
    )
    comparisons = []
    for name, mode, filename, expected, count in [
        (
            "training",
            "profile-replay",
            "training.csv",
            FROZEN / "replay/training.csv",
            84,
        ),
        *(
            (
                s,
                f"source-{s}:100000",
                "source.csv",
                FROZEN / f"point-{int(ENDPOINT.point)}" / s / "source.csv",
                n,
            )
            for s, n in [
                ("roy1980-hcl-mg", 58),
                ("khoo1977-hcl-mg", 30),
                ("jiang1996-hcl-li", 40),
                ("khoo1977-hcl-ca", 30),
            ]
        ),
    ]:
        execute(R / name, mode)
        a, b = read(expected), read(R / name / filename)
        assert len(a) == len(b) == count and list(a.columns) == list(b.columns)
        columns = [c for c in a if c not in ("point", "start")]
        equal(a[columns].to_dict("records"), b[columns].to_dict("records"), True)
        execute(R / "cold" / name, mode)
        equal(
            b.to_dict("records"),
            read(R / "cold" / name / filename).to_dict("records"),
            True,
        )
        assert b.verified.eq(1).all() if name == "training" else b.passes.eq(1).all()
        comparisons.append(
            dict(family=name, rows=count, matches=True, cold_matches=True)
        )
    from score import profile_scores
    import score
    score.D = json.loads(DECL.read_text())
    profile_scores(R / "training")
    metrics = read(R / "training/metrics.csv")
    guards = metrics[metrics.quantity != "Jiang-I5"]
    assert len(guards) == 6 and (guards["mean"] <= guards.mean_limit).all()
    assert (guards.worst <= guards.worst_limit).all() and guards.verified.all()
    save("replay.json", comparisons)
    save("binary-guards.json", guards.to_dict("records"))
    execute(R / "numerical", "numerical")
    execute(R / "cold/numerical", "numerical")
    a, b = read(R / "numerical/numerical.csv"), read(R / "cold/numerical/numerical.csv")
    assert len(a) == 14 and a.passes.eq(1).all()
    equal(a.to_dict("records"), b.to_dict("records"), True)
    save(
        "numerical.json",
        dict(cases=14, max_identity_error=float(a.error.max()), cold_matches=True),
    )


def preservation_arm(path, destination):
    destination.mkdir(parents=True, exist_ok=False)
    document = json.loads(path.read_text())
    folder = ROOT / "analyses/2026-chloride-brine-refit/scripts"
    sys.path.insert(0, str(folder))
    chloride = module("evaluate", folder / "evaluate.py")
    compare = chloride.fit.comparison()
    scenarios = compare.read(compare.INPUTS / "frozen-realizations.csv")
    assert len(scenarios) == 1467
    chloride.evaluate(
        compare,
        destination,
        destination / "mixed",
        document,
        scenarios,
        chloride.agreement_limits(compare, scenarios),
    )
    for name in ("fit", "evaluate", "consistency"):
        sys.modules.pop(name, None)
    folder = ROOT / "analyses/2026-lithium-chloride-refit/scripts"
    sys.path.insert(0, str(folder))
    lithium = module("cation_mixing", folder / "cation_mixing.py")
    lithium.OUT = destination / "lithium"
    lithium.OUT.mkdir(parents=True)
    p = Parameters.from_mapping(document, components=lithium.li.COMPONENTS)
    trace = Parameters.from_mapping(document, components=lithium.TRACE_COMPONENTS)
    lithium.parameters = lambda physical, components: trace
    scores = lithium.score(p, None)
    water = lithium.checks.Water(p)
    for source in ("Yao", "Robinson-Covington"):
        nominal, rounds, buckets = lithium.checks.held_out(water, source, "preserved")
        (destination / f"{source}.json").write_text(
            json.dumps(dict(nominal=nominal, rounds=rounds, buckets=buckets))
        )
        if source == "Yao":
            assert len({r["realization"] for r in rounds}) == 455
    scores["trace_divalent"] = []
    for salt, ion in [("MgCl2", "magnesium-cation"), ("CaCl2", "calcium-cation")]:
        activity = lithium.Activities(
            Parameters.from_mapping(
                document, components=("water", "lithium-cation", ion, "chloride-anion")
            ),
            salt,
        )
        scores["trace_divalent"] += [
            dict(
                salt=salt,
                I=I,
                trace=activity.trace(I / 3),
                binary=activity.binary_check(I / 3),
            )
            for I in (1.0, 3.0, 6.0, 8.3)
        ]
    (destination / "scores.json").write_text(json.dumps(scores))
    projections = []
    for components in [compare.COMPONENTS] if hasattr(compare, "COMPONENTS") else []:
        projections.append(list(components))
    projections += [list(p.component_ids), list(trace.component_ids)]
    projections += [
        ["water", "lithium-cation", ion, "chloride-anion"]
        for ion in ("magnesium-cation", "calcium-cation")
    ]
    projections += [
        [
            "water",
            "sodium-cation",
            "magnesium-cation",
            "calcium-cation",
            "chloride-anion",
        ]
    ]
    admitted = []
    for components in projections:
        d = Parameters.from_mapping(document, components=components)._document
        admitted.append(
            {
                k: list(getattr(d, k))
                for k in ("component_ids", "site_component", "site_multiplicity")
            }
        )
        admitted[-1]["entries"] = [
            {
                k: getattr(e, k)
                for k in ("identity", "field", "index", "magnitude", "unit")
            }
            for e in d.entries
        ]
        admitted[-1]["choices"] = [
            getattr(d, k)
            for k in (
                "base_choice",
                "association_choice",
                "electrolyte_choice",
                "ionic_dispersion_choice",
                "permittivity_choice",
                "polar_choice",
            )
        ]
    (destination / "projections.json").write_text(json.dumps(admitted))


def preserve(original, extended, loaded):
    arms = []
    for name, document in [
        ("control", original),
        ("extended", extended),
        ("readback", loaded),
    ]:
        path = R / f"{name}.json"
        path.write_text(json.dumps(document))
        destination = R / "preservation" / name
        subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "--preserve",
                str(path),
                str(destination),
            ],
            check=True,
        )
        arms.append(destination)
    paths = [
        {
            p.relative_to(arm)
            for p in arm.rglob("*")
            if p.suffix in (".csv", ".json") and "packet" not in p.relative_to(arm).parts
        }
        for arm in arms
    ]
    assert paths[0] == paths[1] == paths[2]
    comparisons = []
    for path in sorted(paths[0]):
        if path.name in ("timing.csv", "summary.json"):
            continue
        values = [
            read(arm / path).to_dict("records")
            if path.suffix == ".csv"
            else json.loads((arm / path).read_text())
            for arm in arms
        ]
        equal(values[0], values[1])
        equal(values[0], values[2])
        comparisons.append(dict(path=str(path), binary64_equal=True))
    save("preservation.json", comparisons)


def main():
    global STAGED, ENDPOINT
    R.mkdir(parents=True, exist_ok=False)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    ENDPOINT = read(FROZEN / "frozen.csv").iloc[0]
    assert ENDPOINT.point == 11 and ENDPOINT.start == 2
    original = json.loads((BASE / "parameters/parameters.json").read_text())
    document = copy.deepcopy(original)
    candidate = json.loads((LEAF / "inputs/candidate.json").read_text())
    declaration = json.loads(DECL.read_text())
    hydrogen = candidate["components"][7]
    for coefficient in hydrogen["coefficients"]:
        if coefficient["family"] in declaration["H_fixed"]:
            coefficient["value"]["magnitude"] = declaration["H_fixed"][
                coefficient["family"]
            ]
        coefficient["provenance"] = dict(
            source_id="hydrogen-chloride-calibration",
            locator="Held 2008 Table 2; declared diameter conventions and binary calibration",
        )
    hydrogen["fixed"]["charge_number"]["provenance"]["locator"] = (
        "Table 3, H+ ion label"
    )
    document["components"].append(hydrogen)
    for pair, value in zip(
        candidate["pairs"][-2:], [ENDPOINT.u, ENDPOINT.c], strict=True
    ):
        pair["coefficients"][0]["value"]["magnitude"] = float(value)
        pair["coefficients"][0]["provenance"] = dict(
            source_id="hydrogen-chloride-calibration",
            locator="Frozen two-coordinate calibration endpoint",
        )
        document["pairs"].append(pair)
    source = dict(
        source_id="hydrogen-chloride-calibration",
        citation="Hamer and Wu 1972; Jiang 1996; Held 2008; native two-coordinate HCl calibration",
        use_basis="70 fitted observations and 14 adjusted-water guards of zero fitting weight",
    )
    document["sources"].append(source)
    document["document_version"] = 3
    assert document["components"] == original["components"] + [hydrogen]
    assert document["pairs"][: len(original["pairs"])] == original["pairs"]
    assert len(document["pairs"]) == len(original["pairs"]) + 2
    assert document["sources"][: len(original["sources"])] == original["sources"]
    assert {
        k: v for k, v in document.items() if k not in ("components", "pairs", "sources", "document_version")
    } == {k: v for k, v in original.items() if k not in ("components", "pairs", "sources", "document_version")}
    packet = tomllib.loads((BASE / "packet.toml").read_text())
    summary = (
        "HCl-water through 7 mol/kg and HCl-LiCl at ionic strength 5 are calibration at 298.15 K and 100000 Pa. "
        "The two dispersion coordinates are one point on a strongly correlated ridge, not a uniquely identified pair. "
        "Comparison accuracy is report-only. MgCl2 mixtures are overpredicted on average; HCl-LiCl underprediction grows with ionic strength (about -10% at 7). "
        "No claim for brines above the measured mixture range, including 12.8-14.6 mol/kg; no propagated uncertainty or independent physical validation. "
    )
    accuracy = read(FROZEN / "transfer.csv")
    accuracy = accuracy[
        (accuracy.point == ENDPOINT.point) & (accuracy.start == ENDPOINT.start)
    ]
    for row in accuracy.itertuples():
        summary += f"{row.source}, {row.role}, line {row.group}: mean/worst/signed mean {row.mean:.9f}/{row.worst:.9f}/{row.signed_mean:.9f} percent. "
    domain = dict(
        domain_id="hydrogen-chloride-298-15-k",
        kind="reported-conditions",
        temperature_min=298.15,
        temperature_max=298.15,
        temperature_unit="kelvin",
        pressure_min=100000,
        pressure_max=100000,
        pressure_unit="pascal",
        description=summary,
    )
    manifest = {
        k: packet[k]
        for k in ("title", "description", "license_or_use_basis", "tolerance_rationale")
    }
    manifest.update(
        title="Chloride mixture parameters with calibrated hydrogen",
        description=summary,
        tolerance_rationale="Binary osmotic RMS <=0.011; numerical replay relative 1e-10, absolute 1e-12; inherited results binary64 identical. Comparison errors are stated, not gated.",
    )
    sys.path.insert(0, str(ROOT / "analyses/2026-chloride-brine-refit/scripts"))
    writer = module(
        "packets", ROOT / "analyses/2026-chloride-brine-refit/scripts/packets.py"
    )
    writer.DATA = R / "staging/data"
    writer.DATA.mkdir(parents=True)
    shutil.copy2(ROOT / "data/packet-locks.json", writer.DATA / "packet-locks.json")
    writer.write(
        BASE.parent.name,
        "3",
        manifest,
        packet["sources"] + [source],
        packet["domains"] + [domain],
        document,
        [
            str((BASE / "parameters/parameters.json").relative_to(ROOT)),
            str((FROZEN / "frozen.csv").relative_to(ROOT)),
            str(DECL.relative_to(ROOT)),
        ],
        "Unchanged chloride basis plus the frozen hydrogen calibration endpoint; no refit.",
    )
    STAGED = writer.DATA / "packets" / BASE.parent.name / "3/parameters/parameters.json"
    from validate_packets import validate_packet

    errors = validate_packet(STAGED.parents[1] / "packet.toml")
    assert not errors, errors
    Parameters.from_json(STAGED)
    replay()
    preserve(original, document, json.loads(STAGED.read_text()))
    final = BASE.parent / "3"
    assert not final.exists()
    shutil.copytree(STAGED.parents[1], final)
    assert (final / "parameters/parameters.json").read_bytes() == STAGED.read_bytes()
    shutil.copy2(writer.DATA / "packet-locks.json", ROOT / "data/packet-locks.json")
    save(
        "admission.json",
        dict(
            published=True,
            byte_identical=True,
            sha256=hashlib.sha256(STAGED.read_bytes()).hexdigest(),
        ),
    )


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--preserve":
        preservation_arm(Path(sys.argv[2]), Path(sys.argv[3]))
        sys.exit(0)
    try:
        main()
    except Exception as error:
        if R.exists():
            save(
                "failure.json",
                dict(type=type(error).__name__, cause=str(error), published=False),
            )
        raise
