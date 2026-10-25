"""Compute native conserved-inventory fraction bounds for the two Born cases."""

import argparse
import csv

import epcsaft
from epcsaft import equilibrium as eq

from reproduce import HERE, load_case, problem_from_request, table


def inventory_limits(compiled, request, floor):
    species = request["reaction_system"]["species_ids"]
    records = []
    for r in eq.inventory_bounds(compiled, [7, 5], True, floor):
        assert r.status == "optimal", (r.status, r.termination)
        error = max(r.balance_residual, r.normalization_residual, r.lower_bound_residual)
        assert error <= 1e-10
        records.append((r.species, species[r.species], r.maximum, r.log_fraction, error))
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=int, choices=(1, 2), required=True)
    args = parser.parse_args()
    case = args.case
    mapping, request, anchor = load_case(case)
    model = epcsaft.Mixture(epcsaft.Parameters.from_mapping(mapping))
    compiled = eq.compile_problem(model, problem_from_request(request, anchor))
    policy = {t.name: t.value for t in epcsaft.tolerances()}
    limits = inventory_limits(compiled, request, policy["EQUILIBRIUM_BOX_FLOOR"])
    out = HERE / "results/global-scan"
    out.mkdir(parents=True, exist_ok=True)
    limit_path = out / "inventory-limits.csv"
    old = list(csv.DictReader(limit_path.open())) if limit_path.exists() else []
    stored = [
        (
            r["case"],
            r["coordinate"],
            r["species"],
            r["x_upper"],
            r["ln_x_upper"],
            r["LP_row_max"],
        )
        for r in old
        if int(r["case"]) != case
    ]
    table(
        out,
        "inventory-limits.csv",
        ("case", "coordinate", "species", "x_upper", "ln_x_upper", "LP_row_max"),
        stored + [(case, *r) for r in limits],
    )


if __name__ == "__main__":
    main()
