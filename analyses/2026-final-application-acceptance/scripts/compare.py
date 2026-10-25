"""Compare the final-wheel replays with the retained application results (#64).

Usage: compare.py SCRATCH MEA_CHECKOUT LITHIUM_CHECKOUT, after replay.sh. Writes results/replay.csv
(one row per replayed check) and results/quantities.csv (physical-comparison and process values on
both wheels). Retained values are the committed files at each checkout's HEAD (Engine #61 record
for the MEA actions/sweep/heat modes).

Agreement criterion: every solved value (state, action, prediction) agrees with its retained value to
1e-8 relative, the #61 branch-identity criterion; convergence diagnostics (residuals, KKT, timings,
iteration counts) are gated by each check's own stored criteria, not compared.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import math
import re
import subprocess
import sys
from pathlib import Path

SCRATCH, MEA, LI = (Path(p).resolve() for p in sys.argv[1:4])
ENGINE = Path(__file__).resolve().parents[3]
RESULTS = Path(__file__).resolve().parents[1] / "results"
MEA_RECORD = ENGINE / "analyses/2026-mea-recovered-state-actions/results"
BUNDLE = "analyses/mea_parameter_bundle"
AGREE = 1e-8
RECORDS = "re" "ceipts"  # Lithium's per-state key; the diagnostics it holds are gated, not compared
DIAGNOSTIC = re.compile(
    r"(^|_)(s|seconds|elapsed|wall|time|iterations?|evaluations|work|nfev|residuals?|norm|error|errors|kkt|"
    r"stationarity|affinity|closure|balance|charge|feasibility|complementarity|multiplier|resolution|"
    r"criterion|floor|fugacity|chain|check|hessian|rel|distinctness|step|jacobian|condition|attempts|seed)"
    r"($|_)", re.I)


def head(checkout: Path, path: str) -> str:
    return subprocess.run(["git", "-C", str(checkout), "show", f"HEAD:{path}"],
                          check=True, capture_output=True, text=True).stdout


def rows(text: str) -> list[dict]:
    return list(csv.DictReader(io.StringIO(text)))


def number(value) -> float | None:
    if isinstance(value, bool) or value is None:
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def flatten(obj, prefix=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from flatten(v, f"{prefix}/{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from flatten(v, f"{prefix}[{i}]")
    else:
        yield prefix, obj


def agreement(pairs) -> dict:
    """pairs: (item, retained, final) of solved values; diagnostics excluded by the caller."""
    n = identical = 0
    worst, worst_item = 0.0, ""
    for item, a, b in pairs:
        a, b = number(a), number(b)
        if a is None:
            continue
        n += 1
        # A retained number with no finite final value fails the check.
        change = math.inf if b is None else 0.0 if a == b else abs(a - b) / max(abs(a), abs(b))
        identical += change == 0.0
        if change > worst:
            worst, worst_item = change, item
    return {"compared": n, "identical": identical, "max_rel_change": worst, "worst": worst_item}


def solved(name: str) -> bool:
    return DIAGNOSTIC.search(re.sub(r"\[\d+\]", "", name.split("/")[-1])) is None


def table_pairs(old: list[dict], new: list[dict], key, columns=None):
    new_by = {key(r): r for r in new}
    for r in old:
        other = new_by.get(key(r))
        if other is None:
            yield f"{key(r)} missing", 1.0, None
            continue
        for column in columns or [c for c in r if solved(c)]:
            yield f"{key(r)}:{column}", r[column], other.get(column)


def json_pairs(old, new):
    new_leaves = dict(flatten(new))
    for path, value in flatten(old):
        if solved(path) and number(value) is not None:
            yield path, value, new_leaves.get(path)


checks: list[dict] = []
quantities: list[dict] = []


def check(application, name, retained_wheel, stats, criterion, stored_ok, note=""):
    ok = stored_ok and (stats["compared"] > 0 or stats["worst"].startswith("not compared")) and stats["max_rel_change"] <= AGREE
    checks.append({"application": application, "check": name, "retained_wheel": retained_wheel, **stats,
                   "stored_criterion": criterion, "stored_criterion_met": stored_ok,
                   "result": "pass" if ok else "FAIL", "note": note})


def quantity(application, name, basis, retained, final):
    quantities.append({"application": application, "quantity": name, "basis": basis,
                       "retained": retained, "final": final})


# MEA: Engine #61 modes (retained on ad12b19f, Engine f5ab8102).
new_dir = SCRATCH / "mea-run/results"
old = rows((MEA_RECORD / "sweep.csv").read_text())
new = rows((new_dir / "sweep.csv").read_text())
key = lambda r: (r["identity"], r["path"], r["neighbour"])  # noqa: E731
accepted_new = [r for r in new if r["ok"] == "True"]
new_by, old_keys = {key(r): r for r in new}, {key(r) for r in old}
# Since #308 every request declares both continuation protocols: each retained row keeps its success flag
# and every added row must be accepted.
stored = (all(key(r) in new_by and new_by[key(r)]["ok"] == r["ok"] for r in old)
          and all(r["ok"] == "True" for r in new if key(r) not in old_keys)
          and all(r["criteria_ok"] == "True" and r["same_branch"] == "True" for r in accepted_new))
stats = agreement(table_pairs([r for r in old if r["ok"] == "True"], accepted_new, key, ["pressure_pa", "min_liquid_x", "eta_liquid"]))
check("MEA", "VLE sweep: cold, anchored and continuation paths (79 requests)", "ad12b19f", stats,
      "retained rows keep success flags, added rows accepted; balance/charge/affinity/unscaled rows <= 1e-10, fugacity <= 1e-9, cold branch",
      stored, f"{len(accepted_new)}/{len(new)} path states accepted")
old, new = rows((MEA_RECORD / "actions.csv").read_text()), rows((new_dir / "actions.csv").read_text())
key = lambda r: (r["state"], r["direction"], r["observable"])  # noqa: E731
stored = len(new) == len(old) and all(r["status"] == "Available" and r["passed"] == "True" for r in new)
check("MEA", "Consumed first-order actions vs Richardson re-solves (T, P, CO2 feed, k_ij, calorics)", "ad12b19f",
      agreement(table_pairs(old, new, key, ["value", "action"])), "|a - R| <= 1e-5|R| + e, every row Available",
      stored, f"{sum(r['passed'] == 'True' for r in new)}/{len(new)} rows pass")
old, new = rows((MEA_RECORD / "heat.csv").read_text()), rows((new_dir / "heat.csv").read_text())
stored = len(new) == len(old) and all(r["status"] == "evaluated" for r in new)
check("MEA", "Calorimetry intervals (Engine calorics)", "ad12b19f",
      agreement(table_pairs(old, new, lambda r: r["record_id"], ["predicted_heat_release_kj_per_mol_CO2"])),
      "every interval evaluated", stored, f"{sum(r['status'] == 'evaluated' for r in new)}/{len(new)} evaluated")
summary_new = json.loads((new_dir / "summary.json").read_text())
summary_old = json.loads((MEA_RECORD / "summary.json").read_text())
for mode in ("actions", "sweep", "heat"):
    new_prov, old_prov = summary_new[mode]["provenance"], summary_old[mode]["provenance"]
    for field in ("parameters_sha256", "state_packet_sha256"):
        assert new_prov.get(field) == old_prov.get(field), (mode, field)

# MEA: application scoring (retained on MEA's pin b66c7b96, Engine 443a9da4).
res = f"{BUNDLE}/results/current-best-fit-residuals.csv"
old, new = rows(head(MEA, res)), rows((MEA / res).read_text())
record = json.loads((MEA / f"{BUNDLE}/results/figure-calculation-record.json").read_text())
check("MEA", "pCO2 and packet speciation predictions (161 pressure states + packet speciation)", "b66c7b96",
      agreement(table_pairs(old, new, lambda r: (r["family"], r["observation_id"], r["target"]), ["predicted"])),
      "every state evaluated (failed_states = 0)", record["failed_states"] == 0 and len(new) == len(old),
      f"{record['evaluated_pressure_states']} pressure states, {record['failed_states']} failed")
spec = f"{BUNDLE}/results/observation-coverage/speciation-rows.csv"
old, new = rows(head(MEA, spec)), rows((MEA / spec).read_text())
key = lambda r: (r["record_id"], r["species"])  # noqa: E731
scored = [r for r in new if r["measurement_role"] != "direct_zero"]
check("MEA", "Observation-coverage speciation (274 positive rows, 58 nondetections)", "b66c7b96",
      agreement(table_pairs(old, new, key, ["predicted_mole_fraction"])), "every positive row scored",
      len(new) == len(old) and all(r["status"] == "scored" for r in scored),
      f"{sum(r['status'] == 'scored' for r in scored)}/{len(scored)} scored, {len(new) - len(scored)} nondetections")
cov_old = json.loads(head(MEA, f"{BUNDLE}/results/observation-coverage/summary.json"))
cov_new = json.loads((MEA / f"{BUNDLE}/results/observation-coverage/summary.json").read_text())
for group, stats_old in cov_old["speciation"]["by_source_and_packet"].items():
    quantity("MEA", f"speciation RMS ln(pred/obs), {group}", "calibration/previously accessed data",
             stats_old["rms_ln"], cov_new["speciation"]["by_source_and_packet"][group]["rms_ln"])
for path, value in flatten(cov_old["pressure"]):
    if path.endswith("rms_ln") or path.endswith("aard_percent"):
        quantity("MEA", f"pCO2 {path.strip('/')}", "calibration and Xu 2011 holdout",
                 value, dict(flatten(cov_new["pressure"]))[path])
heat_old = rows(head(MEA, f"{BUNDLE}/results/calorimetry/current-selected-direct-enthalpy-comparison.csv"))
heat_new = rows((new_dir / "heat.csv").read_text())
for partition in sorted({r["campaign_partition"] for r in heat_new}):
    part = [float(r["residual_kj_per_mol_CO2"]) for r in heat_new if r["campaign_partition"] == partition and r["status"] == "evaluated"]
    part_old = [float(r["residual_kj_per_mol_CO2"]) for r in heat_old if r.get("campaign_partition") == partition and number(r.get("residual_kj_per_mol_CO2")) is not None]
    rmse = lambda v: math.sqrt(sum(x * x for x in v) / len(v)) if v else None  # noqa: E731
    quantity("MEA", f"heat of absorption RMSE kJ/mol, {partition} ({len(part_old)} -> {len(part)} intervals)",
             "Kim-Svendsen / Kim 2014 calorimetry", rmse(part_old), rmse(part))

# Lithium (retained on 18c7812f, Engine 1303c119, except volume scaling: c8784666 with the pre-#147 fit).
fits = json.loads(head(LI, "analyses/hbta_topo_calibrated_surrogate/results/application_parameter_fits.json"))
old_fit = next(f for f in fits["fits"] if "dodecane" in json.dumps(f.get("diluent_id", f.get("diluent", ""))))
new_fit = json.loads((SCRATCH / "li-fit.json").read_text())
starts = new_fit.get("starts", [])
coordinates = [(k, old_fit[k], new_fit[k]) for k in ("ln_k_exchange", "ln_k_na_uptake") if k in old_fit]
check("Lithium", "One-stage calibration refit (Zhang Fig. 5 O/A 3:2 D_Li, Table 2 D_Na), five starts", "18c7812f",
      agreement([*coordinates, *json_pairs(old_fit, new_fit)]), "every start converges; strict calibration state",
      bool(starts) and all(s["success"] for s in starts) and new_fit["strict_state_accepted"] is True,
      f"{sum(s['success'] for s in starts)}/{len(starts)} starts, strict calibration state {new_fit['strict_state_accepted']}")
path = "analyses/hbta_topo_configured_reactive_lle/results/zhang2017_isotherm_holdout"
old, new = json.loads(head(LI, path + ".json")), json.loads((LI / (path + ".json")).read_text())
old_rows, new_rows = rows(head(LI, path + ".csv")), rows((LI / (path + ".csv")).read_text())
states_solved = len(new_rows) + 3  # Fig. 5 rows, Table 4, two Table 2 cascade stages
check("Lithium", "Zhang 2017 holdouts (Fig. 5, Table 4, Table 2 cascade; 11 states)", "18c7812f",
      agreement([*table_pairs(old_rows, new_rows, lambda r: (r["organic_to_aqueous"], r["role"])), *json_pairs(
          {k: v for k, v in old.items() if k not in ("engine", RECORDS)},
          {k: v for k, v in new.items() if k not in ("engine", RECORDS)})]),
      "every state strict (balance/charge <= 1e-9, affinity <= 1e-6, KKT <= 1e-7)",
      all(r["status"] == "evaluated" for r in new_rows) and new["table4_naoh_r70"]["status"] == "evaluated"
      and new["table2_two_stage_countercurrent"]["status"] == "evaluated", f"{states_solved} holdout states")
quantity("Lithium", "Fig. 5 holdout max |log10(D_pred/D_obs)|", "Zhang 2017 holdout",
         old["figure5_holdout_max_abs_log10_distribution_ratio"], new["figure5_holdout_max_abs_log10_distribution_ratio"])
quantity("Lithium", "Table 4 R = 70 % D_Li (observed 7.77)", "Zhang 2017 holdout",
         old["table4_naoh_r70"]["predicted_li_distribution"], new["table4_naoh_r70"]["predicted_li_distribution"])
quantity("Lithium", "Table 2 two-stage raffinate Li g/L (observed 0.045)", "Zhang 2017 holdout",
         old["table2_two_stage_countercurrent"]["predicted_raffinate_li_g_l"],
         new["table2_two_stage_countercurrent"]["predicted_raffinate_li_g_l"])
volume = json.loads((LI / "analyses/hbta_topo_configured_reactive_lle/results/volume_scaling.json").read_text())
scaling = volume["scaling_checks"]
check("Lithium", "Direct-state volume scaling 1-1e5 (current fit; retained run predates #147)", "c8784666",
      {"compared": 0, "identical": 0, "max_rel_change": 0.0, "worst": "not compared: the retained run used the pre-#147 fit"},
      "stored per-scale limits; invariance <= 1e-10 over scales 1-1e5",
      all(scaling[k] for k in ("all_required_scales_evaluated", "all_required_metrics_finite", "phase_composition_pass",
                               "phase_density_pass", "extensive_amount_unit_pass")),
      "; ".join(f"{d['scale']:g}: {d['status']}" for d in volume["diagnostics"]) + " (3e5 is the recorded #63 limit)")

live_old = LI / "analyses/hbta_topo_ms2_idaes_optimization/results/live/current"
live_new = SCRATCH / "li-live"
campaign = json.loads((live_new / "campaign.json").read_text())
pairs = []
for unit in sorted(p.name for p in live_old.iterdir() if (p / "result.json").exists()):
    for name in ("stage_rows.csv", "optimized_stage_rows.csv"):
        if (live_old / unit / name).exists():
            pairs += table_pairs(rows((live_old / unit / name).read_text()), rows((live_new / unit / name).read_text()),
                                 lambda r, u=unit: (u, r["case_id"], r["stage"]))
    pairs += json_pairs(json.loads((live_old / unit / "result.json").read_text()),
                        json.loads((live_new / unit / "result.json").read_text()))
units = campaign["units"]
check("Lithium", "Live IDAES campaign (stages 1-5, 70 % cohort, designs 1-4)", "18c7812f", agreement(pairs),
      "Ipopt optimal; KKT stationarity/complementarity 1e-7, feasibility 1e-8; closure 1e-7 mol/h; cohort",
      campaign["status"] == "accepted" and all(u["status"] == "accepted" for u in units.values()),
      f"{sum(u['status'] == 'accepted' for u in units.values())}/{len(units)} units accepted")
der_new = json.loads((live_new / "derivatives/summary.json").read_text())
ladders_old, ladders_new = rows((live_old / "derivatives/ladders.csv").read_text()), rows((live_new / "derivatives/ladders.csv").read_text())
lkey = lambda r: (r["case"], r["axis"], r["row"], r["order"], r["step"])  # noqa: E731
for table in (ladders_old, ladders_new):
    for r in table:
        r["estimate"] = json.loads(r["estimate"])[0]
ladders = {(r["case"], r["axis"], r["row"], r["order"]): json.loads(r["support"])["status"] for r in ladders_new}
check("Lithium", "Consumed derivatives: finite-difference ladders and 50-digit A1 action", "18c7812f",
      agreement(table_pairs(ladders_old, ladders_new, lkey, ["estimate"])),
      "every ladder passes (1e-8 abs + 1e-6 rel first order); mpmath A1 <= 1e-12",
      der_new["verification"] == "passed" and der_new["mpmath_a1_max_relative"] <= 1e-12,
      f"{sum(v == 'passed' for v in ladders.values())}/{len(ladders)} ladders passed, "
      f"verification {der_new['verification']}, mpmath A1 {der_new['mpmath_a1_max_relative']:.2e}")
for unit in ("design1", "design2", "design3", "design4"):
    o, n = (json.loads((d / unit / "result.json").read_text()) for d in (live_old, live_new))
    # Designs 1-2 cannot reach 95 % within the bounds; their result is the upper-bound corner.
    part, fields = (("live_optimization", ("li_extraction_fraction", "extraction_screening_cost_usd_kg_li",
                                           "organic_to_aqueous_ratio")) if "li_extraction_fraction" in o["live_optimization"]
                    else ("upper_bound_corner", ("li_extraction", "organic_to_aqueous_ratio")))
    for field in fields:
        quantity("Lithium", f"{unit} {part} {field}", "process result (unbuffered fresh-organic MS-2 basis)",
                 o[part][field], n[part][field])


def states(path: Path) -> dict:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as stream:
        return {(r["group"], r["sample_id"]): r for r in map(json.loads, stream)}


v3_old = states(LI / "analyses/hbta_topo_ms2_engine_comparison/results/greenfield_current_v3/engine_states.jsonl.gz")
v3_path = next((SCRATCH / "li-v3").glob("engine_states.jsonl*"))
v3_new = states(v3_path)
fields = ("equilibrium_li_loading_fraction", "equilibrium_na_loading_fraction", "hydrogen_out_mol_l_aqueous",
          "li_net_transferred_mol", "na_net_transferred_mol")
generation = json.loads((SCRATCH / "li-v3/generation.json").read_text())
artifact_old = json.loads(head(LI, "data/processed/hbta_topo_ms2_engine_comparison/greenfield_current_v3/loaded_state_surrogate.json"))
artifact_new = dict(flatten(json.loads((SCRATCH / "li-v3/loaded_state_surrogate.json").read_text())))
identity_or_timing = re.compile(r"wheel_sha256|engine_commit|native_core_sha256|fit_seconds")
artifact_changes = sum(value != artifact_new.get(path) for path, value in flatten(artifact_old)
                       if not identity_or_timing.search(path))
check("Lithium", "v3 surrogate training/verification states", "18c7812f",
      agreement((f"{k}:{f}", v3_old[k]["direct_result"].get(f), v3_new.get(k, {}).get("direct_result", {}).get(f))
                for k in v3_old for f in fields),
      "every required state strict; networks p95 <= 0.03, max <= 0.08",
      len(v3_new) == len(v3_old) and not generation["required_training_or_holdout_failures"] and artifact_changes == 0
      and all(r["direct_result"]["status"] in ("evaluated", "application_accepted_candidate") for r in v3_new.values()),
      f"{len(v3_new)} states, generation status {generation['status']}; refitted networks: {artifact_changes} "
      "values differ outside identity and timing fields")

RESULTS.mkdir(exist_ok=True)
for name, table in (("replay.csv", checks), ("quantities.csv", quantities)):
    with (RESULTS / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(table[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
identity = dict(line.split("=", 1) for line in (SCRATCH / "identity.txt").read_text().split())
(RESULTS / "summary.json").write_text(json.dumps({
    "final_wheel_sha256": identity["sha256"], "engine_commit": identity["engine_commit"],
    "mea": {"checkout": subprocess.run(["git", "-C", str(MEA), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
            "parameters_sha256": summary_new["sweep"]["provenance"]["parameters_sha256"],
            "state_packet_sha256": summary_new["sweep"]["provenance"]["state_packet_sha256"]},
    "lithium": {"checkout": subprocess.run(["git", "-C", str(LI), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip(),
                "fit_coordinates": [old_fit.get("ln_k_exchange"), old_fit.get("ln_k_na_uptake")]},
    "agreement_rel": AGREE, "checks_passed": sum(c["result"] == "pass" for c in checks), "checks": len(checks),
}, indent=1) + "\n", encoding="utf-8")
print(json.dumps(checks, indent=1))
