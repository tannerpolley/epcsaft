"""NIST Shomate and ideal-gas reaction comparisons; all evaluations use the native Engine."""

import csv
import hashlib
import json
import math
from pathlib import Path
import epcsaft as ep
from epcsaft import equilibrium as eq

LEAF = Path(__file__).resolve().parents[1]
SOURCE = LEAF / "inputs/source"


def read(name):
    return list(csv.DictReader((SOURCE / name).open()))


def declarations():
    coefficients = read("shomate-coefficients.csv")
    masses = {r["species"]: float(r["molar_mass_g_per_mol"]) for r in read("molar-masses.csv")}
    intervals = {}
    for r in coefficients:
        values = [float(r[k]) for k in ("A", "B", "C", "D", "E", "F", "G", "H_coefficient_kj_per_mol")]
        interval = ep.IdealInterval(
            float(r["temperature_min_k"]),
            float(r["temperature_max_k"]),
            True,
            not (r["species"] == "H2" and r["temperature_max_k"] == "1000"),
            ep.IdealCorrelation(ep.IdealShomate(values, values[-1] * 1000), float(r["standard_pressure_pa"])),
        )
        intervals.setdefault(r["species"], []).append(interval)

    def model(names):
        document = {
            "schema": "epcsaft.parameters",
            "schema_version": 1,
            "components": [
                {"component_id": s, "fixed": {"molar_mass": {"value": {"magnitude": masses[s], "unit": "gram / mole"}}}}
                for s in names
            ],
            "model_families": [
                {"kind": k, "choice": "ideal-gas" if k == "base" else "none"}
                for k in ("base", "association", "electrolyte", "permittivity")
            ],
        }
        record = ep.ThermochemistryRecord(100000, [intervals[s] for s in names])
        parameters = ep.Parameters.from_mapping(document)
        return (ep.Pure(parameters, thermochemistry=record) if len(names) == 1 else ep.Mixture(parameters, thermochemistry=record)), record

    cases = {}
    for name in ("water-gas-shift", "n2o4-no2"):
        source = json.loads((LEAF / f"inputs/cases/{name}/source/source.json").read_text())
        names = [r["id"] for r in source["species"]]
        mixture, record = model(names)
        cases[name] = mixture, record, names, source
    return model, cases


def main():
    model, cases = declarations()
    rows, solutions = [], []

    def compare(case, T, P, quantity, value, reference, tolerance, label):
        error = abs(value - reference) if value is not None else None
        outcome = "pass" if error is not None and math.isfinite(error) and error <= tolerance else "fail"
        fields = "case T_K P_Pa quantity value reference absolute_error tolerance outcome evidence_label".split()
        rows.append(dict(zip(fields, (case, T, P, quantity, value, reference, error, tolerance, outcome, label))))

    for r in read("shomate-reported-properties.csv"):
        name, T, P = r["species"], float(r["temperature_k"]), float(r["standard_pressure_pa"])
        pure, record = model([name])
        state = pure.state(T, P=P, phase="vapor")
        formation = next(x for x in read("shomate-coefficients.csv") if x["species"] == name)
        for quantity, value, ref, limit in (
            ("Cp_J_mol_K", state.ideal_isobaric_heat_capacity, "cp_j_per_mol_k", "cp_rounding_half_unit_j_per_mol_k"),
            ("S_J_mol_K", state.ideal_entropy, "entropy_j_per_mol_k", "entropy_rounding_half_unit_j_per_mol_k"),
            (
                "H_increment_kJ_mol",
                state.ideal_enthalpy / 1000 - float(formation["H_coefficient_kj_per_mol"]),
                "enthalpy_increment_kj_per_mol",
                "enthalpy_increment_rounding_half_unit_kj_per_mol",
            ),
        ):
            compare(name, T, P, quantity, value, float(r[ref]), float(r[limit]), r["evidence_label"])
    for r in read("reaction-thermochemistry-reference.csv"):
        name, T = r["reaction"], float(r["temperature_k"])
        mixture, record, names, source = cases[name]
        result = ep.reaction_thermochemistry_density(record, source["reaction_stoichiometry"], T)
        assert result.status == ep.ThermochemistryStatus.Available and result.k_pressure > 0, result.cause
        expected_g = float(r["delta_g_pressure_standard_j_per_mol"])
        label = r["reference_kind"]
        compare(name, T, 100000, "reaction_Gp_J_mol", result.pressure_gibbs, expected_g, 1e-6, label)
        compare(
            name, T, 100000, "Kp_relative", result.k_pressure / float(r["K_pressure"]), 1, 1e-10, r["reference_kind"]
        )
        conversion = (
            result.ln_k_density
            - result.ln_k_pressure
            - sum(source["reaction_stoichiometry"]) * math.log(100000 / (ep.GAS_CONSTANT_J_PER_MOL_K * T))
        )
        compare(name, T, 100000, "lnK_standard_conversion", conversion, 0, 1e-10, "numerical verification")

    references = read("analytic-amount-reference.csv")
    for name, T, P in dict.fromkeys(
        (r["case_id"], float(r["temperature_k"]), float(r["pressure_pa"])) for r in references
    ):
        mixture, record, names, source = cases[name]
        feed = [1, 1, 0, 0] if name == "water-gas-shift" else [1, 0]
        problem = eq.Problem(
            phases=[eq.Phase("V", kind="vapor")],
            T=T,
            P=P,
            feed=eq.Amounts(dict(zip(names, feed))),
            reactions=[eq.Reaction(dict(zip(names, source["reaction_stoichiometry"])), name=name)],
            row_tolerance=1e-10,
        )
        result = eq.solve_equilibrium(mixture, problem)
        fields = "success message solver_status amounts mole_fractions molar_densities residuals classified_residuals iterations evaluations".split()
        solutions.append(dict(case=name, T_K=T, P_Pa=P, **{k: getattr(result, k) for k in fields}))
        amounts = [result.amounts[0] * x for x in result.mole_fractions] if result.amounts else [None] * len(names)
        for r in (
            r for r in references if (r["case_id"], float(r["temperature_k"]), float(r["pressure_pa"])) == (name, T, P)
        ):
            value, expected = amounts[names.index(r["species"])], float(r["expected_amount_mol"])
            compare(name, T, P, r["species"] + "_amount_mol", value, expected, 2e-9, r["reference_kind"])
        if result.amounts:
            state = mixture.state(T, rho=result.molar_densities[0], x=result.mole_fractions)
            solutions[-1].update(pressure_slope=state.pressure_density_derivative, Z=state.compressibility_factor)
            for quantity, value, tolerance in (
                ("raw_residual", max(map(abs, result.residuals)), 1e-10),
                ("classified_residual", max(map(abs, result.classified_residuals)), 1e-10),
                ("ideal_pressure_relative", state.pressure / P - 1, 1e-10),
                ("mu_res_RT", max(abs(x) for x in state.residual_chemical_potential_over_rt), 0),
            ):
                compare(name, T, P, quantity, value, 0, tolerance, "numerical solution verification")
            atoms = [[1, 0, 1, 0], [1, 1, 2, 0], [0, 2, 0, 2]] if name == "water-gas-shift" else [[2, 1], [4, 2]]
            for i, atom in enumerate(atoms):
                defect = sum(a * (n - f) for a, n, f in zip(atom, amounts, feed))
                compare(name, T, P, f"atom_{i}_balance_mol", defect, 0, 1e-10, "numerical conservation verification")
            assert state.pressure_density_derivative > 0 and all(n > 0 for n in amounts)
    assert len(solutions) == 6 and len({(r["case"], r["T_K"], r["P_Pa"], r["quantity"]) for r in rows}) == len(rows)
    output = LEAF / "results"
    with (output / "comparison.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (output / "solutions.json").write_text(json.dumps(solutions, indent=2) + "\n")
    library = Path(ep._core.__file__)
    identity = dict(
        engine=list(ep.provenance()),
        native_library=str(library),
        native_library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
        source_sha256={
            str(p.relative_to(LEAF)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((LEAF / "inputs").rglob("*"))
            if p.is_file()
        },
    )
    (output / "run-identity.json").write_text(json.dumps(identity, indent=2) + "\n")
    print(f"{len(rows)} comparisons; {sum(r['outcome'] == 'pass' for r in rows)} pass")


if __name__ == "__main__":
    main()
