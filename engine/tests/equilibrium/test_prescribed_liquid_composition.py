"""Verify the prescribed liquid through the native declaration and snapshot bindings."""

import json
from collections.abc import Sequence
from dataclasses import replace
from math import isnan, log, ulp

import pytest

from epcsaft import GAS_CONSTANT_J_PER_MOL_K as R, Mixture, Parameters, State, _core, equilibrium as eq, trace
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


@pytest.fixture(scope="module")
def khudaida_model() -> Mixture:
    return Mixture(Parameters.from_json(parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7), components=SPECIES))


def flash(
    model: Mixture, feed: Sequence[float], pressure: float = P, phases: Sequence[eq.Phase] = PHASES,
) -> tuple[_core.EquilibriumResult, _core.EquilibriumNumericalEstimates]:
    problem = eq.Problem(T=T, P=pressure, feed=eq.Amounts(dict(zip(SPECIES, feed, strict=True))),
                         phases=phases, row_tolerance=ROW)
    compiled = eq.compile_problem(model, problem)
    result = compiled.solve()
    assert result.success and result.requested_tolerance_met, result.message
    return result, estimates(compiled, result)


def estimates(
    compiled: _core.CompiledSystem, result: _core.EquilibriumResult, input_errors: Sequence[float] = (),
) -> _core.EquilibriumNumericalEstimates:
    assert result.temperature is not None and result.pressure is not None
    return compiled.numerical_estimates(
        result.coordinates, result.temperature, result.pressure, result.amounts, result.mole_fractions,
        result.molar_densities, result.packing_fractions, input_errors,
    )


def endpoint(
    model: Mixture, reference: _core.EquilibriumResult,
) -> tuple[eq.Problem, _core.CompiledSystem, _core.EquilibriumResult]:
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


def value(number: float | None) -> float:
    assert number is not None
    return number


def fugacity_evidence(
    model: Mixture, result: _core.EquilibriumResult, errors: _core.EquilibriumNumericalEstimates,
) -> tuple[list[list[float]], list[State], list[float]]:
    x = [result.phase_mole_fractions(p) for p in range(2)]
    states = [model.state(T, rho=result.molar_densities[p], x=x[p], derived_properties=True) for p in range(2)]
    charge = max(abs(side[3] - side[4]) for side in x)
    norm = max(abs(sum(side) - 1.0) for side in x)
    pressure = max(abs(state.pressure - value(result.pressure)) / (rho * R * T)
                   for state, rho in zip(states, result.molar_densities, strict=True))
    assert charge <= STATE and norm <= STATE and pressure <= ROW
    phi = [state.log_fugacity_coefficient for state in states]
    phi_a, phi_b = phi
    assert phi_a is not None and phi_b is not None
    drop = [log(x[0][i] * states[0].pressure / (x[1][i] * states[1].pressure))
            + value(phi_a[i]) - value(phi_b[i])
            for i in range(len(SPECIES))]
    row_floors = [[floor for row, floor in zip(result.rows, errors.row_floor, strict=True)
                   if row.kind == _core.RowKind.Chemical and row.phase == p] for p in range(2)]
    weights = result.amounts or (1.0, 1.0)
    allowance = [
        sum((ROW if not result.amounts else ROW / (weights[p] * x[p][i])) + row_floors[p][i] for p in range(2))
        for i in range(len(SPECIES))
    ]
    voltage = value(result.electric_potential_difference_v)
    mismatch = [abs(d - z * F * voltage / (R * T)) for d, z in zip(drop, CHARGES, strict=True)]
    assert all(value <= limit for value, limit in zip(mismatch, allowance, strict=True))
    print(f"charge={charge:.17g} norm={norm:.17g} pressure_row={pressure:.17g} "
          f"match={max(mismatch):.17g} allowance={allowance} voltage={voltage:.17g}")
    return x, states, allowance


def test_prescribed_liquid_nonzero_voltage_and_midpoint(khudaida_model: Mixture) -> None:
    model = khudaida_model
    reference, reference_errors = flash(model, FEED)
    problem, compiled, result = endpoint(model, reference)
    assert isnan(result.objective) and isnan(compiled.evaluate(result.coordinates).objective)
    assert all(e.resolved for e in reference_errors.composition[:len(SPECIES)])
    estimated = estimates(compiled, result, [e.absolute for e in reference_errors.composition[:len(SPECIES)]])
    radius = estimated.pressure.absolute
    x, states, allowance = fugacity_evidence(model, result, estimated)
    pin = max(abs(a - b) for a, b in zip(x[0], reference.phase_mole_fractions(0), strict=True))
    assert pin <= STATE
    voltage = value(result.electric_potential_difference_v)
    pressure = value(result.pressure)
    assert abs(voltage) > min(allowance[3:]) * R * T / F
    partial_volumes = states[0].partial_molar_volume
    assert partial_volumes is not None
    info = sum(y * value(v) for y, v in zip(x[1], partial_volumes, strict=True)) - 1.0 / states[1].molar_density
    print(f"pin={pin:.17g} P={pressure:.17g} interval=[{pressure-radius:.17g}, "
          f"{pressure+radius:.17g}] radius={radius:.17g} D={info:.17g} rank={estimated.rank} condition={estimated.condition:.17g}")
    assert estimated.rank == estimated.dimension == 3 * len(SPECIES) + 4
    assert estimated.pressure.resolved and abs(pressure - P) <= radius <= 1_000, "pressure unresolved"

    midpoint = [(a + b) / 2 for a, b in zip(*x, strict=True)]
    split_phases = tuple(replace(phase, composition_guess=tuple(guess)) for phase, guess in zip(PHASES, x, strict=True))
    split, split_error = flash(model, midpoint, pressure, split_phases)
    back, _, split_allowance = fugacity_evidence(model, split, split_error)
    discrepancy = max(abs(a - b) for a, b in zip((*x[0], *x[1]), (*back[0], *back[1]), strict=True))
    voltage_error = abs(value(split.electric_potential_difference_v) - voltage)
    voltage_allowance = min(a + b for a, b in zip(allowance[3:], split_allowance[3:], strict=True)) * R * T / F
    print(f"midpoint amounts={split.amounts} radii={[e.absolute for e in split_error.amounts]} "
          f"composition_delta={discrepancy:.17g} voltage_delta={voltage_error:.17g} voltage_allowance={voltage_allowance:.17g}")
    assert discrepancy <= 1e-9 and voltage_error <= voltage_allowance
    contrast = [abs(a - b) for a, b in zip(*x, strict=True)]
    i = max(range(len(contrast)), key=contrast.__getitem__)
    rounding = [ulp(n) / 2 for n in midpoint]
    # Propagate feed rounding to each amount using the most separated species' lever rule.
    for p, (a, e) in enumerate(zip(split.amounts, split_error.amounts, strict=True)):
        midpoint_error = (rounding[i] + x[1 - p][i] * sum(rounding)) / contrast[i]
        assert abs(a - 0.5) <= e.absolute + midpoint_error

    warm = eq.solve_equilibrium(model, problem, start=result)
    assert warm.success and warm.amounts == [], warm.message
    assert max(abs(a - b) for a, b in zip(warm.mole_fractions, result.mole_fractions, strict=True)) <= STATE


def test_previous_prescribed_liquid_in_declared_starts(khudaida_model: Mixture) -> None:
    reference, _ = flash(khudaida_model, FEED)
    problem, _, result = endpoint(khudaida_model, reference)
    starts = eq.solve_equilibrium_starts(
        khudaida_model, problem, [eq.EquilibriumStart(previous=result)], eq.EquilibriumWorkLimits(),
    )
    assert starts.selected == 0, starts.causes
    warm = starts.results[0]
    assert warm is not None and warm.amounts == []
    assert max(abs(a - b) for a, b in zip(warm.mole_fractions, result.mole_fractions, strict=True)) <= STATE


def test_prescribed_liquid_seeds_feed_split(khudaida_model: Mixture) -> None:
    reference, _ = flash(khudaida_model, FEED)
    _, _, result = endpoint(khudaida_model, reference)
    x = [result.phase_mole_fractions(p) for p in range(2)]
    midpoint = {species: (a + b) / 2 for species, a, b in zip(SPECIES, *x, strict=True)}
    problem = eq.Problem(T=T, P=value(result.pressure), feed=eq.Amounts(midpoint), phases=PHASES, row_tolerance=ROW)
    with trace(diagnostics=True) as captured:
        split = eq.solve_equilibrium(khudaida_model, problem, start=result)
    assert split.success and split.requested_tolerance_met, split.message
    events, spans = captured.events, captured.spans
    initial = next(e.values for e in events if e.name == "row_values"
                   and (parent := spans[e.span].parent) is not None and spans[parent].name == "ipopt_eval_g")
    for row, residual in zip(split.rows, initial, strict=True):
        if row.kind == _core.RowKind.Moiety:
            # Without previous amounts the default is the feed total per phase.
            i = row.index
            assert abs(residual - (x[0][i] + x[1][i] - midpoint[SPECIES[i]] / sum(midpoint.values()))) <= ROW
        if row.kind == _core.RowKind.Pressure:
            assert abs(residual) <= ROW
    assert split.amounts == pytest.approx([0.5, 0.5], rel=0, abs=STATE)


def test_prescribed_uncertainty_cannot_change_support(khudaida_model: Mixture) -> None:
    phases = (
        replace(PHASES[0], composition_guess=(0.85, 0.0, 0.05, 0.05, 0.05)),
        replace(PHASES[1], composition_guess=(0.2, 0.0, 0.7998, 0.0001, 0.0001)),
    )
    reference, _ = flash(khudaida_model, [*FEED[:1], 0.0, *FEED[2:]], phases=phases)
    _, compiled, result = endpoint(khudaida_model, reference)
    assert estimates(compiled, result, [0.0] * len(SPECIES)).pressure.resolved
    errors = estimates(compiled, result, [0.0, STATE, 0.0, 0.0, 0.0])
    assert not errors.pressure.resolved and "support" in errors.cause


def test_distinct_symmetric_ions_return_zero_voltage() -> None:
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
    assert abs(value(result.electric_potential_difference_v)) <= voltage_allowance
