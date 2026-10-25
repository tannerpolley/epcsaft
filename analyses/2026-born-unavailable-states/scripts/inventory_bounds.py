"""Reproduce #234 a1 using the shared native conserved-inventory bounds."""

import csv
import hashlib
import json
from pathlib import Path

import epcsaft
from epcsaft import equilibrium as eq

from global_scan import inventory_limits
from reproduce import HERE, load_case, problem_from_request, table


def main():
    out = HERE / "results/inventory-bounds"
    out.mkdir(parents=True, exist_ok=True)
    comparison = out / "comparison.csv"
    reference_path = comparison if comparison.exists() else HERE / "results/global-scan/inventory-limits.csv"
    old = {(int(r["case"]), int(r["coordinate"])): r for r in csv.DictReader(reference_path.open())}
    floor = next(t.value for t in epcsaft.tolerances() if t.name == "EQUILIBRIUM_BOX_FLOOR")
    rows, limits, witnesses = [], [], {}
    for case in (1, 2):
        mapping, request, anchor = load_case(case)
        model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping))
        c = eq.compile_problem(model, problem_from_request(request, anchor))
        # Exercise the converted consumer, then retain its native witnesses separately.
        limits.extend((case, *r) for r in inventory_limits(c, request, floor))
        for r in eq.inventory_bounds(c, [7, 5], True, floor):
            reference = old[case, r.species]
            error = r.maximum - float(reference.get("old_fraction", reference.get("x_upper")))
            log_error = r.log_fraction - float(reference.get("old_log_fraction", reference.get("ln_x_upper")))
            assert abs(error) <= 1e-10 and abs(log_error) <= 1e-8
            assert max(r.balance_residual, r.normalization_residual, r.lower_bound_residual, r.finite_residual) <= 1e-10
            rows.append((case, r.species, reference["species"], float(reference.get("old_fraction", reference.get("x_upper"))), r.maximum,
                         error, float(reference.get("old_log_fraction", reference.get("ln_x_upper"))), r.log_fraction, log_error,
                         r.balance_residual, r.normalization_residual, r.lower_bound_residual,
                         r.finite_residual, r.lp_error, r.status, r.termination,
                         r.feasibility_termination, r.limiting))
            witnesses[f"{case}/{r.species}"] = {"maximizer_y_tau": r.witness,
                "maximizer_amounts_mol": r.amounts_mol, "finite_y_tau": r.finite_witness}
    table(out, "comparison.csv", ("case", "coordinate", "species", "old_fraction", "native_fraction",
          "fraction_error", "old_log_fraction", "native_log_fraction", "log_error",
          "balance_residual", "normalization_residual", "lower_bound_residual", "finite_residual",
          "lp_error", "status", "termination", "feasibility_termination", "limiting"), rows)
    table(HERE / "results/global-scan", "inventory-limits.csv",
          ("case", "coordinate", "species", "x_upper", "ln_x_upper", "LP_row_max"), sorted(limits))
    root = Path(__file__).resolve().parents[3]
    wheel, = (root / "build/environment-wheel").glob("epcsaft-*.whl")
    provenance = {"wheel_sha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
        "input_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                         for p in sorted((HERE / "inputs").glob("*.json"))},
        "floor": floor, "units": "aggregate mole fraction; witness tau=N_feed/N", "witnesses": witnesses}
    (out / "witnesses.json").write_text(json.dumps(provenance, separators=(",", ":")) + "\n")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
