"""One native imposed-IAPWS equilibrium and raw-row diagnostic contraction."""

import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import epcsaft as ep
from epcsaft import equilibrium as eq
from epcsaft.parameters import parameter_document

LEAF = Path(__file__).resolve().parents[1]
ROOT = LEAF.parents[2]


def declarations():
    source = json.loads((LEAF / "inputs/cases/water-self-ionization/source/source.json").read_text())
    names = source["provider_component_ids"]
    document = parameter_document(ROOT / "data/packets", [("held-cameretti-sadowski-2008/4", names)])
    water = document["components"][0]["fixed"]["molar_mass"]
    original = water["value"].copy()
    water["value"] = dict(magnitude=source["solvent_molar_mass_kg_per_mol"], unit="kilogram / mole")
    model = ep.Mixture(ep.Parameters.from_mapping(document))
    T, P = source["temperature_k"], source["pressure_pa"]
    law = eq.ReactionLogPolynomial(source["ln_kw"], 0, 0, 0, T, T, T, source["source_standard_state_id"])
    reference = eq.ReactionReference(eq.ReferenceSourceBasis.CommonMolalityInfiniteDilution, P)
    reaction = eq.Reaction(dict(zip(names, [-2, 1, 1])), correlation=law, reference=reference)
    neutral = eq.NeutralReference(
        solvent=names[0],
        solute_molality_ratios=dict.fromkeys(names[1:], 1),
        standard_molality_mol_per_kg=source["standard_molality_mol_per_kg"],
    )
    feed = eq.Amounts({names[0]: 1})
    problem = eq.Problem(
        T=T,
        P=P,
        phases=[eq.Phase("L", kind="liquid")],
        feed=feed,
        reactions=[reaction],
        neutral_reference=neutral,
        row_tolerance=1e-12,
    )
    return model, problem, source, original, document


def main():
    model, problem, source, original, document = declarations()
    build = ROOT / "build/issue-789/2019-iapws"
    build.mkdir(parents=True, exist_ok=True)
    driver = build / "reaction_rows"
    subprocess.run(
        ["nice", "-n", "10", "c++", "-std=c++17", LEAF / "scripts/reaction_rows.cpp", "-o", driver], check=True
    )
    result = eq.solve_equilibrium(model, problem)
    compiled = eq.compile_problem(model, problem)
    evaluation = compiled.evaluate(result.coordinates).values
    fields = "success message solver_status amounts mole_fractions molar_densities residuals classified_residuals coordinates iterations evaluations requested_tolerance_met".split()
    diagnostic = {k: getattr(result, k) for k in fields}
    diagnostic["temperature_K"] = result.temperature
    diagnostic["rows"] = [
        dict(kind=r.kind.name, phase=r.phase, index=r.index, eqid=r.eqid, raw=v, retained=raw, classified=c)
        for r, v, raw, c in zip(result.rows, evaluation, result.residuals, result.classified_residuals)
    ]

    def row_id(r):
        return r.kind, r.phase, r.index, r.eqid

    assert list(map(row_id, compiled.rows)) == list(map(row_id, result.rows))
    request = " ".join(map(str, [result.amounts[0], *result.mole_fractions, len(result.rows)])) + "\n"
    request += "".join(
        f"{int(r.kind == ep._core.RowKind.Chemical)} {r.phase} {r.index} {v!r} {raw!r}\n"
        for r, v, raw in zip(result.rows, evaluation, result.residuals)
    )
    diagnostic.update(
        json.loads(subprocess.run([driver], input=request, text=True, capture_output=True, check=True).stdout)
    )
    reference = result.neutral_reference_diagnostics[0]
    ref_fields = "reaction_index terminal_contraction observed_terminal_change reference_pressure_pa message".split()
    diagnostic["reference"] = {k: getattr(reference, k) for k in ref_fields}
    diagnostic["reference"]["status"] = reference.status.name
    state = model.state(source["temperature_k"], rho=result.molar_densities[0], x=result.mole_fractions)
    diagnostic.update(pressure_pa=state.pressure, pressure_slope=state.pressure_density_derivative)
    diagnostic.update(
        pressure_relative=state.pressure / source["pressure_pa"] - 1,
        reference_terminal_change=reference.observed_terminal_change,
    )
    limits = json.loads((LEAF / "inputs/live-reproduction.json").read_text())["numerical_limits"]
    checks = [(k, diagnostic[k], limit) for k, limit in limits.items()]
    positive_ions = all(math.isfinite(n) and n > 0 for n in diagnostic["species_amounts_mol"][1:])
    positive_slope = math.isfinite(state.pressure_density_derivative) and state.pressure_density_derivative > 0
    available = reference.status == eq.NeutralReferenceStatus.Available
    finite_reference = reference.terminal_contraction is not None and math.isfinite(reference.terminal_contraction)
    flags = dict(positive_ions=positive_ions, positive_pressure_slope=positive_slope)
    flags.update(
        reference_available=available and finite_reference,
        solver_available=result.success and result.requested_tolerance_met,
    )
    checks.extend((k, int(not ok), 0) for k, ok in flags.items())
    rows = []
    for k, v, limit in checks:
        error = None if v is None else abs(v)
        outcome = "pass" if error is not None and math.isfinite(error) and error <= limit else "fail"
        row = dict(zip("quantity value absolute_error tolerance outcome".split(), (k, v, error, limit, outcome)))
        row["evidence_label"] = "numerical solution verification with imposed IAPWS constant"
        rows.append(row)
    output = LEAF / "results"
    with (output / "comparison.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    (output / "diagnostics.json").write_text(json.dumps(diagnostic, indent=2) + "\n")
    library = Path(ep._core.__file__)
    inputs = [
        *sorted((LEAF / "inputs").rglob("*.json")),
        ROOT / "data/packets/held-cameretti-sadowski-2008/4/parameters/parameters.json",
    ]
    identity = dict(
        engine=list(ep.provenance()),
        native_library=str(library),
        native_library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),
    )
    identity.update(packet_water_mass=original, admitted_water_mass_kg_mol=model.parameters.molar_masses[0])
    identity["source_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    identity["unavailable"] = "source molalities, activities and independent ln Q unavailable pending #793"
    (output / "run-identity.json").write_text(json.dumps(identity, indent=2) + "\n")
    print(result.message, [(r["quantity"], r["value"], r["outcome"]) for r in rows])


if __name__ == "__main__":
    main()
