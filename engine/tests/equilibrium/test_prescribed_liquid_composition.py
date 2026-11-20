"""Verify the prescribed liquid through the native declaration and snapshot bindings."""

import json
from dataclasses import replace
from math import log

from epcsaft import GAS_CONSTANT_J_PER_MOL_K as R, Mixture, Parameters, equilibrium as eq
from tests.packets import parameters_json

SPECIES = ("water", "ethanol", "isobutanol", "sodium-cation", "chloride-anion")
T = 293.15
P = 100_000.0
FEED = (0.646925133827633, 0.07045345803743744, 0.22564721508764596, 0.028487096523641797, 0.028487096523641797)
ROW = 1e-10
STATE = 1e-12
F = 6.02214076e23 * 1.602176634e-19
CHARGES = (0, 0, 0, 1, -1)
PHASES = (
    eq.Phase("water-rich", kind="liquid", composition_guess=(0.8, 0.05, 0.05, 0.05, 0.05)),
    eq.Phase("organic-rich", kind="liquid", composition_guess=(0.2, 0.1, 0.6998, 0.0001, 0.0001)),
)


def flash(model, feed, pressure=P, phases=PHASES):
    problem = eq.Problem(T=T, P=pressure, feed=eq.Amounts(dict(zip(SPECIES, feed, strict=True))),
                         phases=phases, row_tolerance=ROW)
    compiled = eq.compile_problem(model, problem)
    result = compiled.solve()
    assert result.success and result.requested_tolerance_met, result.message
    return result, estimates(compiled, result)


def estimates(compiled, result, input_errors=()):
    return compiled.numerical_estimates(
        result.coordinates, result.temperature, result.pressure, result.amounts, result.mole_fractions,
        result.molar_densities, result.packing_fractions, input_errors,
    )


def endpoint(model, reference):
    b = reference.phase_mole_fractions(1)
    b[0] += 1e-3
    b[2] -= 1e-3
    phases = (PHASES[0], replace(PHASES[1], composition_guess=tuple(b)))
    problem = eq.Problem(T=T, P=eq.Free(P + 1_000), phases=phases, row_tolerance=ROW,
                         prescribed_liquid_composition=dict(zip(SPECIES, reference.phase_mole_fractions(0), strict=True)))
    compiled = eq.compile_problem(model, problem)
    result = eq.solve_equilibrium(model, problem)
    assert result.success and result.requested_tolerance_met, result.message
    assert result.amounts == []
    return problem, compiled, result


def fugacity_evidence(model, result, errors):
    x = [result.phase_mole_fractions(p) for p in range(2)]
    states = [model.state(T, rho=result.molar_densities[p], x=x[p], derived_properties=True) for p in range(2)]
    charge = max(abs(side[3] - side[4]) for side in x)
    norm = max(abs(sum(side) - 1.0) for side in x)
    pressure = max(abs(state.pressure - result.pressure) / (rho * R * T)
                   for state, rho in zip(states, result.molar_densities, strict=True))
    assert charge <= STATE and norm <= STATE and pressure <= ROW
    drop = [log(x[0][i] * states[0].pressure / (x[1][i] * states[1].pressure))
            + states[0].log_fugacity_coefficient[i] - states[1].log_fugacity_coefficient[i]
            for i in range(len(SPECIES))]
    row_floors = [[floor for row, floor in zip(result.rows, errors.row_floor, strict=True)
                   if row.kind.name == "Chemical" and row.phase == p] for p in range(2)]
    weights = result.amounts or (1.0, 1.0)
    allowance = [
        sum((ROW if not result.amounts else ROW / (weights[p] * x[p][i])) + row_floors[p][i] for p in range(2))
        for i in range(len(SPECIES))
    ]
    voltage = result.electric_potential_difference_v
    mismatch = [abs(d - z * F * voltage / (R * T)) for d, z in zip(drop, CHARGES, strict=True)]
    assert all(value <= limit for value, limit in zip(mismatch, allowance, strict=True))
    print(f"charge={charge:.17g} norm={norm:.17g} pressure_row={pressure:.17g} "
          f"match={max(mismatch):.17g} allowance={allowance} voltage={voltage:.17g}")
    return x, states, allowance


def test_prescribed_liquid_nonzero_voltage_and_midpoint():
    model = Mixture(Parameters.from_json(parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7), components=SPECIES))
    reference, reference_errors = flash(model, FEED)
    problem, compiled, result = endpoint(model, reference)
    assert all(e.resolved for e in reference_errors.composition[:len(SPECIES)])
    estimated = estimates(compiled, result, [e.absolute for e in reference_errors.composition[:len(SPECIES)]])
    radius = estimated.pressure.absolute
    x, states, allowance = fugacity_evidence(model, result, estimated)
    pin = max(abs(a - b) for a, b in zip(x[0], reference.phase_mole_fractions(0), strict=True))
    assert pin <= STATE
    assert abs(result.electric_potential_difference_v) > min(allowance[3:]) * R * T / F
    info = sum(y * v for y, v in zip(x[1], states[0].partial_molar_volume, strict=True)) - 1.0 / states[1].molar_density
    print(f"pin={pin:.17g} P={result.pressure:.17g} interval=[{result.pressure-radius:.17g}, "
          f"{result.pressure+radius:.17g}] radius={radius:.17g} D={info:.17g} rank={estimated.rank} condition={estimated.condition:.17g}")
    assert estimated.rank == estimated.dimension == 3 * len(SPECIES) + 4
    assert estimated.pressure.resolved and abs(result.pressure - P) <= radius <= 1_000, "pressure unresolved"

    midpoint = [(a + b) / 2 for a, b in zip(*x, strict=True)]
    split_phases = tuple(replace(phase, composition_guess=tuple(guess)) for phase, guess in zip(PHASES, x, strict=True))
    split, split_error = flash(model, midpoint, result.pressure, split_phases)
    back, _, split_allowance = fugacity_evidence(model, split, split_error)
    discrepancy = max(abs(a - b) for a, b in zip((*x[0], *x[1]), (*back[0], *back[1]), strict=True))
    voltage_error = abs(split.electric_potential_difference_v - result.electric_potential_difference_v)
    voltage_allowance = min(a + b for a, b in zip(allowance[3:], split_allowance[3:], strict=True)) * R * T / F
    print(f"midpoint amounts={split.amounts} radii={[e.absolute for e in split_error.amounts]} "
          f"composition_delta={discrepancy:.17g} voltage_delta={voltage_error:.17g} voltage_allowance={voltage_allowance:.17g}")
    assert discrepancy <= 1e-9 and voltage_error <= voltage_allowance
    # Both-phase composition and amount estimates are from the independent flash.
    assert all(abs(a - 0.5) <= e.absolute for a, e in zip(split.amounts, split_error.amounts, strict=True))

    warm = eq.solve_equilibrium(model, problem, start=result)
    assert warm.success and warm.amounts == [], warm.message
    assert max(abs(a - b) for a, b in zip(warm.mole_fractions, result.mole_fractions, strict=True)) <= STATE


def test_distinct_symmetric_ions_return_zero_voltage():
    mapping = json.loads(parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7).read_text())
    components = {c["component_id"]: c for c in mapping["components"]}
    cation, anion = SPECIES[-2:]
    ion = json.loads(json.dumps(components[cation]).replace(cation, anion))
    ion["name"] = components[anion]["name"]
    ion["formula"] = components[anion]["formula"]
    ion["fixed"]["charge_number"] = components[anion]["fixed"]["charge_number"]
    mapping["components"][mapping["components"].index(components[anion])] = ion
    pairs = {frozenset((p["component_id_a"], p["component_id_b"])): p for p in mapping["pairs"]}
    for solvent in SPECIES[:3]:
        positive = pairs[frozenset((cation, solvent))]
        negative = pairs[frozenset((anion, solvent))]
        for source, target in zip(positive["coefficients"], negative["coefficients"], strict=True):
            target["value"] = source["value"]
    model = Mixture(Parameters.from_mapping(mapping, components=SPECIES))
    feed = [0.998 * v / sum(FEED[:3]) for v in FEED[:3]] + [0.001, 0.001]
    reference, _ = flash(model, feed)
    _, compiled, result = endpoint(model, reference)
    error = estimates(compiled, result)
    x, _, allowance = fugacity_evidence(model, result, error)
    difference = abs(x[0][0] - x[1][0])
    separation_error = error.composition[0].absolute + error.composition[len(SPECIES)].absolute
    voltage_allowance = R * T * sum(allowance[3:]) / (2 * F)
    print(f"symmetric separation={difference:.17g} separation_error={separation_error:.17g} "
          f"voltage={result.electric_potential_difference_v:.17g} allowance={voltage_allowance:.17g}")
    assert difference > max(1e-6, separation_error)
    assert abs(result.electric_potential_difference_v) <= voltage_allowance
