#!/usr/bin/env python3
"""Compare warm FeOs and Clapeyron density calls on the retained state grid."""

import copy
import hashlib
import importlib.metadata
import importlib.util
import json
import math
import os
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import run as run_utils

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
CACHE = Path.home() / ".cache/epcsaft-comparators/222"
THREAD_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "JULIA_NUM_THREADS",
    "JULIA_NUM_PRECOMPILE_TASKS",
    "RAYON_NUM_THREADS",
)
ENV = os.environ | {key: "1" for key in THREAD_VARS}
ENV.update(PYTHONNOUSERSITE="1", PYTHONPATH="")
ENV["JULIA_DEPOT_PATH"] = str(CACHE / "depot")

LIBRARIES = ("FeOs", "Clapeyron")
PRIMARY_MODE = "unique states (one pass)"
REPEATED_MODE = "repeated identical input (rewards caching)"
ACCURACY_LIMIT = 1e-10
FEOS_SCREENING_NOTE = (
    "FeOs revised Ionic uses 0.88sigma packing diameters; no sigma screening option"
)
MEA_NOTE = "No cross-check parameter conversion for the nine-species MEA packet"


def validate_cases(document):
    cases = [q for q in document["cases"] if q["tier"] in ("A", "B")]
    assert cases, "no Tier A/B cases"
    for q in cases:
        assert int(q["N"]) == 32 and int(q["warmups"]) == 5, (
            f"{q['id']}: protocol fields differ"
        )
        states = q["states"]
        assert len(states) == 32, f"{q['id']}: expected 32 states"
        indices = [int(state["index"]) for state in states]
        assert len(set(indices)) == 32, f"{q['id']}: repeated state index"
        signatures = set()
        for state in states:
            for key in (
                "T_K",
                "P_Pa",
                "x",
                "branch",
                "oracle_rho_mol_m3",
                "oracle_branch",
                "oracle_stable_root_count",
            ):
                assert key in state, f"{q['id']} state {state['index']}: missing {key}"
            assert math.isfinite(float(state["T_K"])) and float(state["T_K"]) > 0
            assert math.isfinite(float(state["P_Pa"])) and float(state["P_Pa"]) > 0
            assert (
                math.isfinite(float(state["oracle_rho_mol_m3"]))
                and float(state["oracle_rho_mol_m3"]) > 0
            )
            signature = (
                float(state["T_K"]),
                float(state["P_Pa"]),
                tuple(float(value) for value in state["x"]),
                str(state["branch"]),
                float(state["sigma_angstrom"]) if q["model"] == "water" else None,
            )
            assert signature not in signatures, f"{q['id']}: states are not unique"
            signatures.add(signature)
            if q["model"] == "water":
                assert (
                    math.isfinite(float(state["sigma_angstrom"]))
                    and float(state["sigma_angstrom"]) > 0
                )
    return cases


def unavailable_reason(library, case):
    if case["model"] == "MEA":
        return MEA_NOTE
    if library == "FeOs" and case["model"] == "C3":
        return FEOS_SCREENING_NOTE
    return None


def request_for(case):
    return {
        "case": case["model"],
        "diameter": "sigma",
        "electrostatic": "matched",
        "kcl": -0.25,
        "dipole_convention": "matched",
    }


def density_sample(mode, sample_index, state, elapsed_us, rho, status, note=""):
    oracle = float(state["oracle_rho_mol_m3"])
    delta = abs(float(rho) / oracle - 1.0) if rho is not None else None
    return {
        "mode": mode,
        "sample_index": sample_index,
        "state_index": int(state["index"]),
        "elapsed_us": elapsed_us,
        "rho_mol_m3": float(rho) if rho is not None else None,
        "oracle_rho_mol_m3": oracle,
        "oracle_branch": state["oracle_branch"],
        "oracle_stable_root_count": int(state["oracle_stable_root_count"]),
        "relative_density_difference": delta,
        "accuracy_within_1e10": delta <= ACCURACY_LIMIT if delta is not None else None,
        "status": status,
        "note": note,
    }


def feos_kwargs(state, np, kelvin, pascal):
    return {
        "temperature": float(state["T_K"]) * kelvin,
        "pressure": float(state["P_Pa"]) * pascal,
        "composition": np.asarray(state["x"], dtype=float),
        "density_initialization": str(state["branch"]),
    }


def prepare_feos_model(conversion, document, case, state):
    inputs = copy.deepcopy(document["comparator_inputs"])
    if case["model"] == "water":
        inputs["models"]["water"]["components"][0]["sigma"] = float(
            state["sigma_angstrom"]
        )
    eos, _, _ = conversion.prepare(inputs, request_for(case))
    return eos


def feos_runtime(eos, state, State, np, kelvin, pascal):
    return eos, feos_kwargs(state, np, kelvin, pascal), State


def warm_feos(runtime):
    eos, kwargs, State = runtime
    State(eos, **kwargs)


def timed_feos(runtime, mode, sample_index, state, rho_units):
    eos, kwargs, State = runtime
    begun = time.perf_counter_ns()
    try:
        result = State(eos, **kwargs)
    except Exception as error:  # Retain the measured failed call.
        elapsed = (time.perf_counter_ns() - begun) / 1000
        return density_sample(
            mode,
            sample_index,
            state,
            elapsed,
            None,
            "solve_failed",
            f"{type(error).__name__}: {error}",
        )
    elapsed = (time.perf_counter_ns() - begun) / 1000
    try:
        rho = float(result.density / rho_units)
    except Exception as error:  # Density extraction is outside the timed call.
        return density_sample(
            mode,
            sample_index,
            state,
            elapsed,
            None,
            "solve_failed",
            f"{type(error).__name__}: {error}",
        )
    if not math.isfinite(rho) or rho <= 0:
        return density_sample(
            mode,
            sample_index,
            state,
            elapsed,
            None,
            "solve_failed",
            f"invalid density {rho}",
        )
    return density_sample(mode, sample_index, state, elapsed, rho, "represented")


def feos_case(document, case, conversion, State, np, kelvin, pascal, rho_units):
    states = case["states"]
    primary, preparation_errors, models = [], [], {}
    for state in states:
        try:
            sigma = float(state["sigma_angstrom"]) if case["model"] == "water" else None
            if sigma not in models:
                models[sigma] = prepare_feos_model(conversion, document, case, state)
            primary.append(
                feos_runtime(models[sigma], state, State, np, kelvin, pascal)
            )
            preparation_errors.append("")
        except Exception as error:  # Retain failed preparation against this state.
            primary.append(None)
            preparation_errors.append(f"{type(error).__name__}: {error}")
    try:
        repeated_model = prepare_feos_model(conversion, document, case, states[0])
        repeated = feos_runtime(repeated_model, states[0], State, np, kelvin, pascal)
        repeated_error = ""
    except Exception as error:
        repeated, repeated_error = None, f"{type(error).__name__}: {error}"

    warmup_errors = {PRIMARY_MODE: [], REPEATED_MODE: []}
    for _ in range(5):
        if primary[0] is None:
            warmup_errors[PRIMARY_MODE].append(preparation_errors[0])
        else:
            try:
                warm_feos(primary[0])
            except Exception as error:
                warmup_errors[PRIMARY_MODE].append(f"{type(error).__name__}: {error}")
        if repeated is None:
            warmup_errors[REPEATED_MODE].append(repeated_error)
        else:
            try:
                warm_feos(repeated)
            except Exception as error:
                warmup_errors[REPEATED_MODE].append(f"{type(error).__name__}: {error}")

    samples = []
    for state, runtime, prep_error in zip(
        states, primary, preparation_errors, strict=True
    ):
        if runtime is None:
            samples.append(
                density_sample(
                    PRIMARY_MODE, 0, state, None, None, "solve_failed", prep_error
                )
            )
        else:
            samples.append(timed_feos(runtime, PRIMARY_MODE, 0, state, rho_units))
    for sample_index in range(int(case["N"])):
        state = states[0]
        if repeated is None:
            samples.append(
                density_sample(
                    REPEATED_MODE,
                    sample_index,
                    state,
                    None,
                    None,
                    "solve_failed",
                    repeated_error,
                )
            )
        else:
            samples.append(
                timed_feos(repeated, REPEATED_MODE, sample_index, state, rho_units)
            )
    return {"id": case["id"], "samples": samples, "warmup_errors": warmup_errors}


def feos_worker(replica, interpreter):
    sys.path = [
        path
        for path in sys.path
        if Path(path or ".").resolve() != Path(__file__).parent
    ]
    import numpy as np
    from feos import State, set_num_threads
    from si_units import KELVIN, METER, MOL, PASCAL

    set_num_threads(1)
    source = ROOT / "analyses/2026-competitor-cross-check/scripts/feos.py"
    spec = importlib.util.spec_from_file_location("density_feos_conversion", source)
    conversion = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(conversion)
    identity = {
        "type": "ready",
        "replica": replica,
        "pid": os.getpid(),
        "version": importlib.metadata.version("feos"),
        "effective_nice": os.getpriority(os.PRIO_PROCESS, 0),
        "interpreter": interpreter,
    }
    print(json.dumps(identity), flush=True)
    for line in sys.stdin:
        command = {}
        try:
            command = json.loads(line)
            case = command["case"]
            result = feos_case(
                command["document"],
                case,
                conversion,
                State,
                np,
                KELVIN,
                PASCAL,
                MOL / METER**3,
            )
            print(json.dumps(result, allow_nan=False), flush=True)
        except Exception as error:
            print(
                json.dumps(
                    {
                        "id": command.get("case", {}).get("id"),
                        "error_type": type(error).__name__,
                        "error": str(error),
                    }
                ),
                flush=True,
            )


def clapeyron_worker_command(replica, interpreter):
    julia = str((CACHE / "julia-1.13.1/bin/julia").absolute())
    script = str((HERE / "scripts/competitors.jl").absolute())
    project = str((ROOT / "analyses/2026-competitor-cross-check/inputs").absolute())
    return [
        julia,
        "--threads=1",
        "--startup-file=no",
        "--project=" + project,
        script,
        "--worker",
        str(replica),
        interpreter,
        "10",
    ]


def launch_worker(library, replica, inherited_nice):
    print(f"starting {library} replica {replica}", flush=True)
    log = tempfile.TemporaryFile(mode="w+t")
    if library == "FeOs":
        interpreter = str((CACHE / "feos/bin/python").absolute())
        command = [
            interpreter,
            str(Path(__file__).absolute()),
            "--feos-worker",
            str(replica),
            interpreter,
        ]
    else:
        interpreter = str((CACHE / "julia-1.13.1/bin/julia").absolute())
        command = clapeyron_worker_command(replica, interpreter)
    process = subprocess.Popen(
        ["nice", "-n", str(10 - inherited_nice), *command],
        cwd=ROOT,
        env=ENV,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=log,
        text=True,
        bufsize=1,
    )
    line = process.stdout.readline()
    if not line:
        log.seek(0)
        detail = log.read()
        process.wait()
        raise RuntimeError(
            f"{library} replica {replica} failed during startup: {detail}"
        )
    identity = json.loads(line)
    assert identity["type"] == "ready" and int(identity["replica"]) == replica
    assert int(identity["pid"]) == process.pid, (
        f"{library} replica {replica}: process identity changed"
    )
    effective_nice = os.getpriority(os.PRIO_PROCESS, process.pid)
    assert effective_nice == 10, (
        f"{library} replica {replica}: expected nice 10, got {effective_nice}"
    )
    assert int(identity["effective_nice"]) == 10, (
        f"{library} replica {replica}: worker did not inherit nice 10"
    )
    identity.update(
        library=library, interpreter=interpreter, nice=effective_nice, threadpools=1
    )
    return {"process": process, "log": log, "identity": identity}


def request_case(worker, document, case):
    process = worker["process"]
    process.stdin.write(
        json.dumps({"document": document, "case": case}, allow_nan=False) + "\n"
    )
    process.stdin.flush()
    line = process.stdout.readline()
    if not line:
        worker["log"].seek(0)
        raise RuntimeError(
            f"worker {worker['identity']['library']} replica {worker['identity']['replica']} stopped: {worker['log'].read()}"
        )
    result = json.loads(line)
    if result.get("error"):
        raise RuntimeError(f"{case['id']}: {result['error_type']}: {result['error']}")
    assert result["id"] == case["id"]
    return result


def stop_workers(workers):
    for worker in workers.values():
        process = worker["process"]
        if process.stdin and not process.stdin.closed:
            process.stdin.close()
    for worker in workers.values():
        process = worker["process"]
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.terminate()
            process.wait()
        worker["log"].close()


def unavailable_case(case, library, reason, version):
    return {
        "id": case["id"],
        "tier": case["tier"],
        "model": case["model"],
        "library": library,
        "version": version,
        "N": int(case["N"]),
        "processes": 3,
        "warmups": int(case["warmups"]),
        "status": "unrepresentable",
        "states": case["states"],
        "unavailable_states": [
            {
                **state,
                "state_index": int(state["index"]),
                "status": "unrepresentable",
                "note": reason,
            }
            for state in case["states"]
        ],
        "warmup_errors": [],
        "samples": [],
    }


def percentile_text(samples, p):
    return run_utils.percentile(samples, p) if samples else ""


def summarize_case(case, library, version, samples, warmups, unavailable_states):
    timed = [sample for sample in samples if sample["elapsed_us"] is not None]
    primary_times = [float(s["elapsed_us"]) for s in timed if s["mode"] == PRIMARY_MODE]
    repeated_times = [
        float(s["elapsed_us"]) for s in timed if s["mode"] == REPEATED_MODE
    ]
    represented = [s for s in samples if s["rho_mol_m3"] is not None]
    differences = [float(s["relative_density_difference"]) for s in represented]
    if not samples:
        status = "unrepresentable" if unavailable_states else "solve_failed"
    elif unavailable_states or any(s["status"] == "unrepresentable" for s in samples):
        status = "unrepresentable"
    elif any(s["status"] != "represented" for s in samples):
        status = "partial" if represented else "solve_failed"
    else:
        status = "represented"
    initial_index = int(case["states"][0]["index"])
    initial_densities = [
        float(s["rho_mol_m3"])
        for s in samples
        if s["mode"] == PRIMARY_MODE
        and int(s["state_index"]) == initial_index
        and s["rho_mol_m3"] is not None
    ]
    accuracy = (
        ""
        if not represented
        else str(all(delta <= ACCURACY_LIMIT for delta in differences)).lower()
    )
    notes = [entry["note"] for entry in unavailable_states if entry.get("note")]
    notes.extend(s["note"] for s in samples if s.get("note"))
    notes.extend(
        error
        for item in warmups
        for errors in item["errors"].values()
        for error in errors
    )
    notes = list(dict.fromkeys(notes))
    notes.append(
        "rho_mol_m3 is the median primary density at the first state; accuracy covers every represented timed density in both modes."
    )
    return {
        "id": case["id"],
        "library": library,
        "version": version,
        "status": status,
        "rho_mol_m3": statistics.median(initial_densities) if initial_densities else "",
        "median_us": statistics.median(primary_times) if primary_times else "",
        "p10_us": percentile_text(primary_times, 0.1),
        "p90_us": percentile_text(primary_times, 0.9),
        "N": int(case["N"]),
        "processes": 3,
        "warmups": int(case["warmups"]),
        "max_relative_density_difference": max(differences) if differences else "",
        "accuracy_within_1e10": accuracy,
        "repeated_identical_median_us": statistics.median(repeated_times)
        if repeated_times
        else "",
        "note": "; ".join(notes),
    }


def attach_process(sample, identity):
    sample.update(
        library=identity["library"],
        version=identity["version"],
        replica=int(identity["replica"]),
        pid=int(identity["pid"]),
        interpreter=identity["interpreter"],
        nice=int(identity["nice"]),
    )


def main():
    cases_path = HERE / "inputs/cases.json"
    cases_bytes = cases_path.read_bytes()
    document = json.loads(cases_bytes)
    cases = validate_cases(document)
    inherited_nice = os.getpriority(os.PRIO_PROCESS, 0)
    if inherited_nice > 10:
        raise RuntimeError(
            f"inherited niceness {inherited_nice} cannot be lowered to 10"
        )
    workers = {}
    results = {
        (q["id"], library): {"samples": [], "warmups": [], "unavailable_states": []}
        for q in cases
        for library in LIBRARIES
    }
    dispatch = []
    try:
        for library in LIBRARIES:
            for replica in range(3):
                workers[(library, replica)] = launch_worker(
                    library, replica, inherited_nice
                )
        for case_index, case in enumerate(cases):
            for replica in range(3):
                order = (
                    LIBRARIES
                    if (case_index + replica) % 2 == 0
                    else tuple(reversed(LIBRARIES))
                )
                dispatch.append(
                    {"id": case["id"], "replica": replica, "library_order": list(order)}
                )
                for library in order:
                    result = results[(case["id"], library)]
                    reason = unavailable_reason(library, case)
                    if reason:
                        if not result["unavailable_states"]:
                            result["unavailable_states"] = [
                                {
                                    **state,
                                    "state_index": int(state["index"]),
                                    "status": "unrepresentable",
                                    "note": reason,
                                }
                                for state in case["states"]
                            ]
                        continue
                    print(f"{case['id']} {library} replica {replica}", flush=True)
                    worker = workers[(library, replica)]
                    response = request_case(worker, document, case)
                    identity = worker["identity"]
                    for sample in response["samples"]:
                        attach_process(sample, identity)
                    result["samples"].extend(response["samples"])
                    result["warmups"].append(
                        {
                            "replica": replica,
                            "pid": identity["pid"],
                            "errors": response["warmup_errors"],
                        }
                    )
    finally:
        stop_workers(workers)

    raw_cases, summary_rows = [], []
    versions = {
        identity["library"]: identity["version"]
        for identity in (w["identity"] for w in workers.values())
    }
    for case in cases:
        case_results = []
        for library in LIBRARIES:
            data = results[(case["id"], library)]
            if data["unavailable_states"]:
                version = versions[library]
                result = unavailable_case(
                    case, library, data["unavailable_states"][0]["note"], version
                )
            else:
                version = versions[library]
                result = {
                    "id": case["id"],
                    "tier": case["tier"],
                    "model": case["model"],
                    "library": library,
                    "version": version,
                    "N": int(case["N"]),
                    "processes": 3,
                    "warmups": int(case["warmups"]),
                    "status": "pending summary",
                    "states": case["states"],
                    "unavailable_states": [],
                    "warmup_errors": data["warmups"],
                    "samples": data["samples"],
                }
            summary = summarize_case(
                case,
                library,
                version,
                result["samples"],
                result["warmup_errors"],
                result["unavailable_states"],
            )
            result["status"] = summary["status"]
            summary_rows.append(summary)
            case_results.append(result)
        raw_cases.extend(case_results)

    raw = {
        "protocol": {
            "cases_sha256": hashlib.sha256(cases_bytes).hexdigest(),
            "N": 32,
            "warmups_per_mode": 5,
            "processes_per_library": 3,
            "modes": [PRIMARY_MODE, REPEATED_MODE],
            "threadpools": 1,
            "nice": 10,
            "dispatch_order": dispatch,
            "density_accuracy_limit_relative": ACCURACY_LIMIT,
        },
        "processes": [worker["identity"] for worker in workers.values()],
        "cases": raw_cases,
    }
    samples_path = HERE / "results/competitor-samples.json"
    samples_path.write_text(json.dumps(raw, indent=2, allow_nan=False) + "\n")
    run_utils.write_csv(HERE / "results/competitors.csv", summary_rows)
    print(f"wrote {samples_path}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--feos-worker":
        feos_worker(int(sys.argv[2]), sys.argv[3])
    else:
        main()
