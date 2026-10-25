"""Pure MDEA PC-SAFT S6/S9 comparisons using native states and observations."""

import csv
import hashlib
import json
import math
from pathlib import Path
import epcsaft as ep
from epcsaft import equilibrium as eq

LEAF = Path(__file__).resolve().parents[1]
SOURCE = LEAF / "inputs/source"


def declarations():
    source = json.loads((SOURCE / "si-reference.json").read_text())
    p = source["mdea_pcsaft"]
    document = ep._core.ParameterDocument()
    document.component_ids = ["MDEA"]
    document.site_component = [0, 0]
    document.site_multiplicity = [1, 2]
    document.base_choice = "pc-saft"
    document.association_choice = "general-site"
    document.electrolyte_choice = document.permittivity_choice = document.polar_choice = "none"
    entries = []
    for field, value, unit in [
        ("molar_mass", p["molar_mass_kg_mol"], "kilogram / mole"),
        ("charge_number", 0, "elementary-charge"),
        ("segment_count", p["m"], "dimensionless"),
        ("sigma_constant", p["sigma_angstrom"], "angstrom"),
        ("dispersion_energy_over_k", p["epsilon_over_k_K"], "kelvin"),
        ("association_energy_over_k", p["epsilon_assoc_over_k_K"], "kelvin"),
        ("association_volume", p["kappa_assoc"], "dimensionless"),
    ]:
        e = ep._core.ParameterEntry()
        e.identity, e.field, e.magnitude, e.unit = "MDEA/" + field, field, value, unit
        e.index = [0, 1] if field.startswith("association_") else [0]
        entries.append(e)
    document.entries = entries
    model = ep.Pure(
        ep.Parameters(ep._core.build_parameters(document), ("MDEA",), (("MDEA", "donor"), ("MDEA", "acceptor")))
    )
    return model, source, p


def main():
    model, source, document = declarations()
    rows, diagnostics = [], []

    def compare(context, quantity, value, reference, numerical, cause):
        limit = source["comparison_tolerances"][quantity]
        error = abs(value - reference) if value is not None else None
        passed = error is not None and math.isfinite(error) and error <= limit
        outcome = ('pass' if passed else 'fail') if numerical else 'undetermined'
        keys = (
            "table row T_K P_Pa quantity value reference absolute_error tolerance numerical_pass outcome cause".split()
        )
        values = (*context, quantity, value, reference, error, limit, numerical, outcome, cause)
        rows.append(dict(zip(keys, values)))

    for index, r in enumerate(csv.DictReader((SOURCE / "si-table-s6-mdea.csv").open())):
        T, P = float(r["T_K"]), float(r["P_MPa"]) * 1e6
        state = model.state(T, P=P, phase="liquid", derived_properties=True)
        density = state.molar_density * model.parameters.molar_masses[0]  # mol/m3 to kg/m3 reporting conversion
        closure = abs(state.pressure / P - 1)
        numerical = closure <= 1e-8 and state.pressure_density_derivative > 0 and math.isfinite(density)
        diagnostics.append(
            dict(
                zip(
                    "table row T_K P_Pa pressure_relative pressure_slope".split(),
                    ("S6", index, T, P, closure, state.pressure_density_derivative),
                )
            )
        )
        for quantity, field, factor in [
            ("rho_kg_m3", "molar_density", model.parameters.molar_masses[0]),
            ("alpha_kK_inv", "isobaric_expansivity", 1000),
            ("kappa_TPa_inv", "isothermal_compressibility", 1e12),
        ]:
            value = getattr(state, field)
            reference = float(r[quantity.replace("_", "_pcsaft_", 1)])
            compare(
                ("S6", index, T, P),
                quantity,
                None if value is None else value * factor,
                reference,
                numerical,
                "native liquid",
            )
    for index, r in enumerate(csv.DictReader((SOURCE / "si-table-s9-mdea.csv").open())):
        T = float(r["T_K"])
        problem = model.saturation(T)
        problem.row_tolerance = 1e-10
        result = eq.solve_equilibrium(model, problem)
        fields = "success message solver_status pressure amounts mole_fractions molar_densities residuals classified_residuals iterations evaluations requested_tolerance_met".split()
        diagnostic = dict(table="S9", row=index, T_K=T, **{k: getattr(result, k) for k in fields})
        numerical = False
        density = None
        if result.pressure is not None and result.pressure > 0 and len(result.molar_densities) == 2:
            states = [model.state(T, rho=rho) for rho in result.molar_densities]
            checks = []
            for state, rho in zip(states, result.molar_densities):
                pressure, slope = state.pressure, state.pressure_density_derivative
                ulp = math.ulp(rho)
                scale = abs(slope) * ulp
                allowance = max(1e-8 * result.pressure, 8 * scale)
                defect = abs(pressure - result.pressure)
                keys = "pressure_Pa rho_mol_m3 slope density_ulp ulp_pressure_Pa allowance_Pa scaled_defect relative_defect density_backward_error".split()
                checks.append(
                    dict(
                        zip(
                            keys,
                            (
                                pressure,
                                rho,
                                slope,
                                ulp,
                                scale,
                                allowance,
                                defect / allowance,
                                defect / result.pressure,
                                defect / (rho * abs(slope)),
                            ),
                        )
                    )
                )
            chemical = [
                math.log(rho) + s.residual_chemical_potential_over_rt[0]
                for s, rho in zip(states, result.molar_densities)
            ]
            coexistence = abs(chemical[0] - chemical[1])
            density = result.molar_densities[0] * model.parameters.molar_masses[0]
            numerical = (
                result.success
                and result.requested_tolerance_met
                and max(map(abs, result.classified_residuals)) <= 1e-10
                and coexistence <= 1e-8
                and result.molar_densities[0] > result.molar_densities[1]
                and all(
                    c["scaled_defect"] <= 1
                    and c["slope"] > 0
                    and c["pressure_Pa"] > 0
                    and all(map(math.isfinite, c.values()))
                    for c in checks
                )
            )
            rho = result.molar_densities[0]
            diagnostic.update(
                pressure_checks=checks,
                common_standard_mu_RT=chemical,
                coexistence_defect=coexistence,
                liquid_neighbor_pressures_Pa=[
                    model.state(T, rho=math.nextafter(rho, d)).pressure for d in (0, math.inf)
                ],
            )
        diagnostics.append(diagnostic)
        for quantity, value, reference in [
            ("P_sat_Pa", result.pressure, float(r["P_pcsaft_Pa"])),
            ("rho_kg_m3", density, float(r["rho_pcsaft_kg_m3"])),
        ]:
            compare(("S9", index, T, result.pressure), quantity, value, reference, numerical, result.message)
        print("S9", T, result.pressure, numerical, flush=True)
    assert len(rows) == 380 and len(diagnostics) == 130
    output = LEAF / "results"
    with (output / "comparison.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (output / "diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    library = Path(ep._core.__file__)
    identity = dict(
        engine=list(ep.provenance()),
        native_library=str(library),
        native_library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
        declared_model=document,
        evidence_label="numerical verification against printed PC-SAFT model",
        numerical_limits=dict(pressure_relative=1e-8, density_ulps=8, common_standard_mu_RT=1e-8, defining_rows=1e-10),
        source_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCE.iterdir() if p.is_file()},
    )
    (output / "run-identity.json").write_text(json.dumps(identity, indent=2) + "\n")


if __name__ == "__main__":
    main()
