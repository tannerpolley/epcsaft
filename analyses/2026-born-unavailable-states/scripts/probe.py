"""Historical fixed-coordinate row/Jacobian replay; see results/current-main/diagnosis.md.

The 026b3031/e9b1aed5 wheel comparison identifies the #385 reference change.
Current main rejects these old rows; this command does not correct chart roots.
"""

import argparse
import csv
import json
import math
import posixpath
from pathlib import Path

import epcsaft
from epcsaft import equilibrium as eq
import numpy as np

from reproduce import HERE, archived, changed, load_case, problem_from_request


def rows(name):
    # Only the direct-solve rows stay loose under results/ (the README links them).
    text = archived()[name] if name in archived() else (HERE / "results" / name).read_text()
    return list(csv.DictReader(text.splitlines()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", nargs="*", help="retired-chart archive coordinate member names")
    parser.add_argument(
        "--output-dir", type=Path, default=HERE / "results/diagnostics/native/points"
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    names = args.points or sorted(
        n
        for n in archived()
        if n.endswith("-coordinates.csv") and not n.startswith("newton/")
    )
    policy = {t.name: t.value for t in epcsaft.tolerances()}
    for name in names:
        folder, file = posixpath.split(name)
        tag = file.replace("-failed-coordinates.csv", "").replace(
            "-coordinates.csv", ""
        )
        stem = posixpath.join(folder, tag)
        summary = json.loads(archived()[stem + ".json"])
        if "convergence" not in summary:
            continue
        # Default replay targets the previously retained point diagnostics only.
        if not args.points and stem + "-point-diagnostics.json" not in archived():
            continue
        out = args.output_dir / folder
        out.mkdir(parents=True, exist_ok=True)
        case, value = summary["case"], summary["value"]
        base, request, anchor = load_case(case)
        model = epcsaft.Mixture(
            epcsaft.Parameters.from_mapping(changed(base, case, value))
        )
        compiled = eq.compile_problem(model, problem_from_request(request, anchor))
        x = np.array([float(r["value"]) for r in rows(name)])
        ev = compiled.evaluate(x)
        assert (
            max(
                abs(a - b["residual"])
                for a, b in zip(
                    ev.values,
                    (
                        {k: float(v) if k == "residual" else v for k, v in row.items()}
                        for row in rows(stem + "-rows.csv")
                    ),
                    strict=True,
                )
            )
            < 2e-10
        )
        J = np.array(compiled.jacobian(x)).reshape(len(ev.values), len(x))
        row_scale = np.maximum(np.max(abs(J), axis=1), np.finfo(float).tiny)
        scaled = J / row_scale[:, None]
        col_scale = np.maximum(np.max(abs(scaled), axis=0), np.finfo(float).tiny)
        scaled /= col_scale
        sv = np.linalg.svd(scaled, compute_uv=False)
        d = {
            "row_recheck_max": max(abs(v) for v in ev.values),
            "unscaled_condition": float(np.linalg.cond(J)),
            "equilibrated_condition": float(sv[0] / sv[-1]),
            "equilibrated_singular_min": float(sv[-1]),
            "equilibrated_singular_max": float(sv[0]),
            "equilibrated_rank": int(np.linalg.matrix_rank(scaled)),
            "bounds": "square pinned composition cutoff; no inventory box or inventory upper cap",
            "bound_coordinate_file": tag + "-bounds.csv",
        }
        (out / (tag + "-point-diagnostics.json")).write_text(
            json.dumps(d, indent=2) + "\n"
        )
        labels = (
            [f"liquid ln x({s})" for s in request["reaction_system"]["species_ids"]]
            + ["liquid amount mol", "liquid ln packing"]
            + [
                f"vapor ln x({s})"
                for s in request["phases"][1]["support"]["component_ids"]
            ]
            + ["vapor ln packing", "pressure Pa"]
            + [f"moiety multiplier {i}" for i in range(compiled.rank)]
        )
        lo_x = math.log(policy["EQUILIBRIUM_BOX_FLOOR"])
        bounds = (
            [(lo_x, 0)] * 9
            + [
                (0, None),
                (
                    math.log(policy["PACKING_FRACTION_FLOOR"]),
                    math.log(policy["PACKING_FRACTION_CEILING"]),
                ),
            ]
            + [(lo_x, 0)] * 3
            + [
                (
                    math.log(policy["PACKING_FRACTION_FLOOR"]),
                    math.log(policy["PACKING_FRACTION_CEILING"]),
                ),
                (np.finfo(float).eps, None),
            ]
            + [(None, None)] * compiled.rank
        )
        with (out / (tag + "-bounds.csv")).open("w") as f:
            w = csv.writer(f, lineterminator="\n")
            w.writerow(
                (
                    "coordinate",
                    "quantity",
                    "value",
                    "lower",
                    "upper",
                    "distance_to_lower",
                    "distance_to_upper",
                )
            )
            for i, (label, v, (lo, hi)) in enumerate(
                zip(labels, x, bounds, strict=True)
            ):
                w.writerow(
                    (
                        i,
                        label,
                        v,
                        lo,
                        hi,
                        None if lo is None else v - lo,
                        None if hi is None else hi - v,
                    )
                )
        if summary["success"]:
            state_rows = rows(stem + "-state.csv")
            with (out / (tag + "-terms.csv")).open("w") as f:
                w = csv.writer(f, lineterminator="\n")
                w.writerow(
                    (
                        "phase",
                        "born_over_RT",
                        "debye_huckel_over_RT",
                        "bulk_relative_permittivity",
                        "net_charge_mol",
                    )
                )
                for phase, (a, b, N) in enumerate(((0, 9, x[9]), (11, 14, 0.0))):
                    fractions = np.exp(x[a:b])
                    fractions /= fractions.sum()
                    full = list(fractions) + ([0.0] * 6 if phase else [])
                    rho = float(
                        next(
                            r["rho_mol_m3"]
                            for r in state_rows
                            if int(r["phase"]) == phase
                        )
                    )
                    state = model.state(
                        request["temperature"]["value"], rho=rho, x=full
                    )
                    charge = N * sum(
                        z * f
                        for z, f in zip(
                            request["reaction_system"]["charges"], full, strict=True
                        )
                    )
                    w.writerow(
                        (
                            phase,
                            state.born,
                            state.debye_huckel,
                            state.bulk_relative_permittivity,
                            charge,
                        )
                    )
        if tag.endswith("direct"):
            with (out / "balance-basis.csv").open("w") as f:
                w = csv.writer(f, lineterminator="\n")
                w.writerow(("moiety", "species", "coefficient"))
                for k in range(compiled.rank):
                    for i, s in enumerate(request["reaction_system"]["species_ids"]):
                        w.writerow((k, s, -J[i, compiled.primal_size + k]))
        print(tag, d, flush=True)


if __name__ == "__main__":
    main()
