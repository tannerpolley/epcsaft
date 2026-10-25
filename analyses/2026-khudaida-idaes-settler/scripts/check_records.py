"""Several immutable equilibrium records in one ASL process (#452).

`run` evaluates the F1-F4 call sequences against the installed non-editable wheel, one fresh
single-threaded process per mode, and retains every call: value and the requested real-control
derivatives as float.hex, Errmsg, and the stderr lines written between its marker and the next.
"""

import argparse
import csv
import ctypes
import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path

NAME = "epcsaft_equilibrium_phase_component_amount"
MARK = "#452 call "
PREPARED = "epcsaft_equilibrium prepared record "
DIAGNOSTIC = "epcsaft_equilibrium diagnostic: "
F0 = (0.7198946513021484, 0.03377833329223333, 0.21864518731791907, 0.02768182808769922)
F1 = (F0[0], F0[1], F0[2] * 0.95, F0[3] * 0.95)
G0, G1 = (0.7, 0.3), (0.7, 0.285)
B_TEMPERATURE_K, B_PRESSURE_PA = 298.15, 1.0e5
FIELDS = ("mode", "call", "tag", "record", "order", "p", "c", "controls", "value", "gradient", "hessian", "errmsg", "prepared", "diagnostic")


def calls(record, blocks, components, tag=""):
    """(tag, record, controls, order, p, c) for each (controls, order, rows) block."""
    return [(tag, record, list(F), order, r // components, r % components) for F, order, rows in blocks for r in rows]


S_A = calls("A", [(F0, 2, range(10)), (F1, 2, range(10)), (F0, 0, [0]), (F0, 2, range(9, -1, -1))], 5)
S_B = calls("B", [(G0, 2, range(4)), (G1, 2, range(4)), (G0, 0, [0]), (G0, 2, range(3, -1, -1))], 2)
FIRST_A, FIRST_B = calls("A", [(F0, 2, range(10))], 5, "first"), calls("B", [(G0, 2, range(4))], 2, "first")
ERRORS = {
    "a": ("a", "B", [0.7, 0.3, 0.7, 0.3], 2, 0, 0),
    "b": ("b", "A", list(F0[:2]), 2, 0, 0),
    "c": ("c", "B", list(G0), 2, 0, 2),
    "dA": ("dA", "A", list(F0), 2, 2, 0),
    "dB": ("dB", "B", list(G0), 2, 2, 0),
}
MODES = {
    "A": S_A,
    "B": S_B,
    "AB": [call for pair in zip(S_A, S_B) for call in pair] + S_A[len(S_B) :],
    "errors-A": FIRST_A + [ERRORS[k] for k in ("b", "dA")] + [("recovery", *c[1:]) for c in FIRST_A],
    "errors-B": FIRST_B + [ERRORS[k] for k in ("a", "c", "dB")] + [("recovery", *c[1:]) for c in FIRST_B],
    "errors-AB": FIRST_A + FIRST_B + list(ERRORS.values()) + [("recovery", *c[1:]) for c in FIRST_A + FIRST_B],
    "limit": [(tag, f"B{k}", list(G0), 0, 0, 0) for tag in ("", "recovery") for k in range(17)],
}


def records():
    """Record A is the live Khudaida record; B_k is the engine/README.md water/isobutanol record with max_steps 4 + k."""
    import run_live_settler as live
    from epcsaft import Mixture, Parameters
    from epcsaft import equilibrium as eq

    reference = live.Reference()
    model = Mixture(
        Parameters.from_json(
            live.fixed.REPO / "data/packets/khudaida-2026-figure-2-electrolyte-lle/7/parameters/parameters.json",
            components=["water", "isobutanol"],
        )
    )
    both = ("water", "isobutanol")
    problem = eq.Problem(
        phases=[
            eq.Phase("L1", kind="liquid", support=both, composition_guess=(0.46, 0.54)),
            eq.Phase("L2", kind="liquid", support=both, composition_guess=(0.99, 0.01)),
        ],
        T=B_TEMPERATURE_K,
        P=B_PRESSURE_PA,
        feed=eq.Amounts({"water": G0[0], "isobutanol": G0[1]}),
    )
    result = {"A": reference.record}
    for k in range(17):
        options = eq.Continuation(initial_step=0.25, minimum_step=0.25, max_steps=4 + k, max_retries_per_step=1)
        result[f"B{k}"] = eq.asl_record(model, problem, options, [0.0, 0.0], [[1.0, 0.0], [0.0, 1.0]], list(G0), 1, 0.1)
    result["B"] = result["B0"]
    wheel, wheel_sha256, library = live.neutral.installed_identity()
    return result, {"wheel": str(wheel), "wheel_sha256": wheel_sha256, "library": str(library), "library_sha256": live.fixed.sha256(library)}


def mode(spec_path, output_path):
    """Child process: raw _ARGLIST calls in order, no early stop."""
    from pyomo.core.base.external import _ARGLIST
    from pyomo.environ import ExternalFunction

    spec = json.loads(Path(spec_path).read_text())
    external = ExternalFunction(library=spec["library"], function=NAME)
    external.load_library()
    function = external._known_functions[NAME][0]
    rows = []
    for index, (tag, record, controls, order, p, c) in enumerate(spec["calls"]):
        print(f"{MARK}{index}", file=sys.stderr, flush=True)
        packed = _ARGLIST([*controls, spec["records"][record], str(p), str(c)], fgh=order)
        value = function(ctypes.byref(packed))
        n = len(controls)
        rows.append(
            {
                "call": index,
                "tag": tag,
                "record": record,
                "order": order,
                "p": p,
                "c": c,
                "controls": " ".join(map(float.hex, controls)),
                "value": float.hex(value),
                "gradient": " ".join(float.hex(packed.derivs[i]) for i in range(n)) if order >= 1 else "",
                "hessian": " ".join(float.hex(packed.hes[i]) for i in range(n * (n + 1) // 2)) if order == 2 else "",
                "errmsg": packed.Errmsg.decode() if packed.Errmsg else "",
            }
        )
    Path(output_path).write_text(json.dumps(rows))


def run(output):
    output.mkdir(parents=True, exist_ok=True)
    texts, identity = records()
    identity["record_sha256"] = {name: hashlib.sha256(texts[name].encode()).hexdigest() for name in ("A", "B", *(f"B{k}" for k in range(17)))}
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    table = []
    for name, sequence in MODES.items():
        spec, result, log = output / f".{name}.spec.json", output / f".{name}.calls.json", output / f"stderr-{name}.log"
        spec.write_text(json.dumps({"library": identity["library"], "records": {r: texts[r] for r in {c[1] for c in sequence}}, "calls": sequence}))
        with log.open("w") as stream:
            subprocess.run([sys.executable, __file__, "mode", str(spec), str(result)], env=env, stderr=stream, check=True)
        attributed, current = {}, None
        for line in log.read_text().splitlines():
            if line.startswith(MARK):
                current = int(line[len(MARK) :])
            elif current is not None:
                attributed.setdefault(current, []).append(line)
        for row in json.loads(result.read_text()):
            lines = attributed.get(row["call"], [])
            table.append(
                {
                    "mode": name,
                    **row,
                    "prepared": " ".join(line[len(PREPARED) :].split()[0] for line in lines if line.startswith(PREPARED)),
                    "diagnostic": " | ".join(line[len(DIAGNOSTIC) :] for line in lines if line.startswith(DIAGNOSTIC)),
                }
            )
        spec.unlink()
        result.unlink()
    with (output / "calls.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    (output / "identity.json").write_text(json.dumps(identity, indent=2, sort_keys=True) + "\n")


def entries(row):
    return [row["value"], *row["gradient"].split(), *row["hessian"].split()]


def finite(row):
    n = len(row["controls"].split())
    expected = 1 + (n if int(row["order"]) >= 1 else 0) + (n * (n + 1) // 2 if int(row["order"]) == 2 else 0)
    values = entries(row)
    return len(values) == expected and all(math.isfinite(float.fromhex(v)) for v in values)


def failed(row):
    return all(math.isnan(float.fromhex(v)) for v in entries(row))


def recovered(first, again):
    return all(
        abs(float.fromhex(a) - float.fromhex(b)) <= 1e-11 + 1e-9 * abs(float.fromhex(a)) for a, b in zip(entries(first), entries(again), strict=True)
    )


EXPECTED = {
    "a": ("defining_equations_not_satisfied:phase[0].amount[0]", "control/observable dimension mismatch"),
    "b": ("defining_equations_not_satisfied:phase[0].amount[0]", "control/observable dimension mismatch"),
    "c": ("defining_equations_not_satisfied:phase[0].amount[2]", "control/observable dimension mismatch"),
    "dA": ("defining_equations_not_satisfied:request", "invalid observable selector"),
    "dB": ("defining_equations_not_satisfied:request", "invalid observable selector"),
}


def expected_error(row):
    status, diagnostic = EXPECTED[row["tag"] or row["record"]]
    return {"2": '"', "1": "'"}.get(row["order"], "") + "epcsaft_equilibrium:" + status, diagnostic


def evaluate(wheel, by, comparisons, errors):
    """Criterion outcomes of one wheel; appends the F1 per-call rows and the F3/F4 error rows."""
    for x in "AB":
        for a, b in zip([r for r in by["AB"] if r["record"] == x], by[x]):
            comparisons.append(
                {"check": "F1", "wheel": wheel, "record": x, "a_mode": "AB", "a_call": a["call"], "b_mode": x, "b_call": b["call"],
                 "order": a["order"], "p": a["p"], "c": a["c"], "a": " ".join(entries(a)), "b": " ".join(entries(b)),
                 "bitwise": entries(a) == entries(b) and finite(a)}
            )
    prepared = {m: [int(k) for r in rows for k in r["prepared"].split()] for m, rows in by.items()}
    outcome = {"preflight": all(finite(r) for m in "AB" for r in by[m])}
    outcome["F1"] = len(by["AB"]) == len(by["A"]) + len(by["B"]) and all(r["bitwise"] for r in comparisons if r["check"] == "F1" and r["wheel"] == wheel)
    outcome["F2"] = prepared["A"] == [1] and prepared["B"] == [1] and prepared["AB"] == [1, 2]
    outcome["F3 errors-AB prepares no third record"] = prepared["errors-AB"] == [1, 2]
    for m in ("errors-A", "errors-B", "errors-AB", "limit"):
        ok = m != "limit" or prepared[m] == list(range(1, 18))
        for r in by[m]:
            if (r["tag"] or r["record"]) in EXPECTED:
                message, diagnostic = expected_error(r)
                match = failed(r) and r["errmsg"] == message and r["diagnostic"] == diagnostic
                errors.append({"wheel": wheel, "mode": m, "call": r["call"], "tag": r["tag"], "record": r["record"], "controls": r["controls"],
                               "order": r["order"], "p": r["p"], "c": r["c"], "errmsg": r["errmsg"], "diagnostic": r["diagnostic"],
                               "expected_errmsg": message, "expected_diagnostic": diagnostic, "nonfinite_outputs": failed(r), "match": match})
                ok = ok and match
        first, again = [r for r in by[m] if r["tag"] == "first"], [r for r in by[m] if r["tag"] == "recovery"]
        if m == "limit":
            first, ok = by[m][:17], ok and all(finite(r) for r in by[m])
        outcome["F4" if m == "limit" else f"F3 {m}"] = ok and all(finite(b) and recovered(a, b) for a, b in zip(first, again, strict=True))
    failures = [r for r in by["AB"] if r["record"] == "B" and failed(r) and r["diagnostic"] == "immutable equilibrium record changed"]
    return outcome, {"prepared": prepared, "AB_B_calls_rejected_as_record_changed": len(failures)}


def drop_timing(value):
    """Every *_s key is a time.monotonic() duration."""
    if isinstance(value, dict):
        return {k: drop_timing(v) for k, v in value.items() if not k.endswith("_s")}
    return [drop_timing(v) for v in value] if isinstance(value, list) else value


def khudaida(main, changed):
    rows = [{"campaign": "khudaida", "file": name, "comparison": "bytes", "equal": (main / name).read_bytes() == (changed / name).read_bytes()}
            for name in ("transport.csv", "ladders.csv", "endpoints.csv", "record.netstrings")]

    def summary(path):
        document = json.loads(path.read_text())
        document.pop("identity")
        for case in (document["baseline"], *document["optimization"]):
            for key in ("elapsed_s", "option_file_sha256"):
                case.pop(key)
            case["options"].pop("output_file")
        return json.dumps(document, sort_keys=True)

    rows.append({"campaign": "khudaida", "file": "summary.json", "comparison": "JSON without identity and per-case elapsed_s, options.output_file, option_file_sha256",
                 "equal": summary(main / "summary.json") == summary(changed / "summary.json")})
    sha = [json.loads((d / "identity.json").read_text())["record_sha256"] for d in (main, changed)]
    rows.append({"campaign": "khudaida", "file": "identity.json", "comparison": "record_sha256 only", "equal": sha[0] == sha[1]})
    for path in sorted(main.glob("failure-*.log")):
        before, after = path.read_text().splitlines(), (changed / path.name).read_text().splitlines()
        if path.name == "failure-record_change.log":
            before = [line.replace("immutable equilibrium record changed", "trailing equilibrium record data") for line in before]
        rows.append({"campaign": "khudaida", "file": path.name, "comparison": "lines; changed wheel adds exactly one prepared line" + (", record-changed diagnostic becomes the parse message" if "record_change" in path.name else ""),
                     "equal": [line for line in after if not line.startswith(PREPARED)] == before and sum(line.startswith(PREPARED) for line in after) == 1})
    return rows


def lithium(runs):
    """Every .json/.csv/.opt/.netstrings file of each run, compared across runs after the timing and library-identity exclusions."""
    names = {label: sorted(str(p.relative_to(d)) for p in d.rglob("*") if p.suffix in (".json", ".csv", ".opt", ".netstrings")) for label, d in runs.items()}
    rows, timings, libraries = [], {}, {}
    for label, d in runs.items():
        identity = json.loads((d / "identity.json").read_text())
        libraries[label] = identity["engine_library_sha256"]
        assert json.loads((d / "derivatives/summary.json").read_text())["library_sha256"] == libraries[label], label
        campaign = json.loads((d / "campaign.json").read_text())
        assert campaign["status"] == "accepted" and len(campaign["units"]) == 13, label
        assert all(len((d / n).read_text().splitlines()) > 1 for n in names[label] if n.endswith("stage_rows.csv")), label
    reference = names[next(iter(runs))]
    rows.append({"campaign": "lithium", "file": "(file set)", "comparison": f"{len(reference)} files", "equal": all(n == reference for n in names.values())})
    for name in reference:
        texts = []
        for label, d in runs.items():
            text = (d / name).read_text()
            if name.endswith(".json"):
                document = json.loads(text)
                for path, value in flatten(document):
                    if path.endswith("_s"):
                        timings.setdefault((name, path), {})[label] = value
                document = drop_timing(document)
                for key in {"identity.json": ("engine_wheel", "engine_library_sha256"), "derivatives/summary.json": ("library_sha256",)}.get(name, ()):
                    document.pop(key)
                text = json.dumps(document, sort_keys=True)
            texts.append(text)
        rows.append({"campaign": "lithium", "file": name, "comparison": "JSON without *_s keys" + {"identity.json": ", engine_wheel, engine_library_sha256", "derivatives/summary.json": ", library_sha256"}.get(name, "") if name.endswith(".json") else "bytes",
                     "equal": all(t == texts[0] for t in texts)})
    return rows, [{"file": f, "key": k, **v} for (f, k), v in timings.items()], libraries


def flatten(value, prefix=""):
    if isinstance(value, dict):
        for k, v in value.items():
            yield from flatten(v, f"{prefix}.{k}" if prefix else k) if isinstance(v, (dict, list)) else [(f"{prefix}.{k}" if prefix else k, v)]
    elif isinstance(value, list):
        for i, v in enumerate(value):
            yield from flatten(v, f"{prefix}[{i}]")


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def report(records, khudaida_runs, lithium_runs):
    """records holds main/ and slice/ run directories; writes the comparison tables and summary.json beside them."""
    by = {}
    for wheel in ("main", "slice"):
        for r in csv.DictReader((records / wheel / "calls.csv").open()):
            by.setdefault(wheel, {}).setdefault(r["mode"], []).append(r)
    comparisons, errors = [], []
    outcomes = {wheel: evaluate(wheel, rows, comparisons, errors) for wheel, rows in by.items()}
    for x in "AB":
        for a, b in zip(by["main"][x], by["slice"][x], strict=True):
            comparisons.append({"check": "single-record", "wheel": "main vs slice", "record": x, "a_mode": x, "a_call": a["call"], "b_mode": x, "b_call": b["call"],
                                "order": a["order"], "p": a["p"], "c": a["c"], "a": " ".join(entries(a)), "b": " ".join(entries(b)), "bitwise": entries(a) == entries(b)})
    files = khudaida(*khudaida_runs)
    lithium_rows, timings, libraries = lithium(lithium_runs)
    files += lithium_rows
    labels = list(lithium_runs)
    for row in timings:  # each slice run against the main runs immediately before and after it
        for i in range(1, len(labels) - 1, 2):
            spread = (row[labels[i - 1]], row[labels[i + 1]])
            row[f"{labels[i]}_within_main_spread"] = min(spread) <= row[labels[i]] <= max(spread)
    identities = {wheel: json.loads((records / wheel / "identity.json").read_text()) for wheel in by}
    summary = {
        "identity": {w: {k: v for k, v in i.items() if k != "library"} for w, i in identities.items()},
        "record_sha256_agree": identities["main"]["record_sha256"] == identities["slice"]["record_sha256"],
        "lithium_engine_library_sha256": libraries,
        "outcomes": {w: o for w, (o, _) in outcomes.items()},
        "observations": {w: obs for w, (_, obs) in outcomes.items()},
        "expected_main": "preflight passes; F1 fails because B is rejected as record changed; no prepared lines; F3 errors-A and errors-B pass, errors-AB setup unreachable; F4 unreachable",
        "single_record_callback_bitwise": all(r["bitwise"] for r in comparisons if r["check"] == "single-record"),
        "files_equal": all(r["equal"] for r in files),
        "lithium_slice_times_outside_main_spread": {
            k[: -len("_within_main_spread")]: sum(not r[k] for r in timings) for k in timings[0] if k.endswith("_within_main_spread")
        },
        "lithium_timing_fields": len(timings),
        "limits": "Numerical verification of the callback only. Bitwise identity in one machine and build; per-record memory unmeasured; the Khudaida campaign omits the boundary and unavailable probes (#461); two main timing samples give no statistical bound.",
    }
    write_csv(records / "comparisons.csv", comparisons)
    write_csv(records / "errors.csv", errors)
    write_csv(records / "files.csv", files)
    write_csv(records / "timings.csv", timings)
    (records / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: summary[k] for k in ("outcomes", "record_sha256_agree", "single_record_callback_bitwise", "files_equal", "lithium_slice_times_outside_main_spread")}, indent=1))


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "mode":
        mode(sys.argv[2], sys.argv[3])
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--output-dir", type=Path, help="run every mode with the installed wheel into this directory")
        parser.add_argument("--report", type=Path, help="records directory holding main/ and slice/ runs")
        parser.add_argument("--khudaida", type=Path, nargs=2, metavar=("MAIN", "SLICE"))
        parser.add_argument("--lithium", type=Path, nargs="+", help="runs in time order: main-1 slice-1 main-2 [slice-2 main-3]")
        args = parser.parse_args()
        if args.report:
            report(args.report.resolve(), args.khudaida, dict(zip(("main-1", "slice-1", "main-2", "slice-2", "main-3"), args.lithium, strict=False)))
        else:
            run(args.output_dir.resolve())
