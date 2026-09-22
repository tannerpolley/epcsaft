"""S1 (#51): positional native inputs against the Python name mapping. Manufactured check."""
from collections.abc import Sequence

import pytest

from epcsaft import Mixture, Parameters, _core
from epcsaft.equilibrium import Amounts, Phase, Problem

from tests.packets import parameters_json

ORDER = ("water", "sodium-cation", "chloride-anion")
PACKET = parameters_json("figiel-2025-reference-electrolytes", 3)
# charge, segment count, sigma / angstrom, epsilon / K, solvation factor, ion diameter, Born diameter
SPECIES = (
    (0, 1.2047, 2.7927, 353.95, 1.5, None, None),
    (1, 1.0, 2.8232, 230.0, 1.0, 2.484416, 3.445),
    (-1, 1.0, 2.756, 170.0, 1.0, 2.42528, 4.1),
)


def _entry(field: str, index: Sequence[int], magnitude: float, unit: str = "dimensionless") -> _core.ParameterEntry:
    item = _core.ParameterEntry()
    item.identity, item.field, item.index = f"{field}/{tuple(index)}", field, list(index)
    item.magnitude, item.unit = magnitude, unit
    return item


def _positional_model() -> _core.HelmholtzModel:
    """The packet's water / Na+ / Cl- as a C++ caller states it: positions, units and every z."""
    document = _core.ParameterDocument()
    document.component_ids = list(ORDER)
    document.site_component, document.site_multiplicity = [0, 0], [1, 1]
    document.base_choice = "pc-saft"
    document.association_choice = "general-site"
    document.electrolyte_choice = "born"
    document.ionic_dispersion_choice = "exclude-same-sign-ion-pairs"
    document.permittivity_choice = "ion-fraction-suppression"
    entries = [_entry("molar_mass", (0,), 18.0153, "gram / mole")]
    for i, (charge, segment, sigma, epsilon, solvation, diameter, born) in enumerate(SPECIES):
        entries += [
            _entry("charge_number", (i,), charge, "elementary-charge"),
            _entry("segment_count", (i,), segment),
            _entry("sigma_constant", (i,), sigma, "angstrom"),
            _entry("dispersion_energy_over_k", (i,), epsilon, "kelvin"),
            _entry("solvation_factor", (i,), solvation),
        ]
        if diameter is not None and born is not None:
            entries += [
                _entry("packing_diameter", (i,), diameter, "angstrom"),
                _entry("debye_huckel_diameter", (i,), diameter, "angstrom"),
                _entry("born_diameter", (i,), born, "angstrom"),
            ]
    entries += [
        _entry("relative_permittivity", (0,), 78.09),
        _entry("sigma_term_amplitude", (0, 0), 10.11, "angstrom"),
        _entry("sigma_term_exponent", (0, 0), -0.01775, "1 / kelvin"),
        _entry("sigma_term_amplitude", (1, 0), -1.417, "angstrom"),
        _entry("sigma_term_exponent", (1, 0), -0.01146, "1 / kelvin"),
        _entry("k_ij", (0, 1), -0.3),
        _entry("k_ij", (0, 2), -0.3),
        _entry("k_ij", (1, 2), 0.8),
        _entry("born_c_shell", (), 1.0),
        _entry("born_c_dielectric", (), 1.0),
        _entry("ionic_region_relative_permittivity", (), 8.0),
        _entry("ion_suppression_coefficient", (), 7.01),
        _entry("association_energy_over_k", (0, 1), 2425.7, "kelvin"),
        _entry("association_volume", (0, 1), 0.04509),
    ]
    document.entries = entries
    return _core.HelmholtzModel(_core.build_parameters(document))


def test_positional_document_equals_the_named_packet() -> None:
    x, rho = (0.96, 0.02, 0.02), 5.7e4  # compressed liquid: the 1 bar liquid root is 5.63e4 mol/m^3
    positional = _positional_model().state(298.15, rho, x)
    named = Mixture(Parameters.from_json(PACKET, components=tuple(reversed(ORDER)))).state(
        298.15, rho=rho, x=tuple(reversed(x))
    )
    assert positional.pressure == pytest.approx(named.pressure, rel=1e-12, abs=1e-12)
    assert positional.residual_helmholtz == pytest.approx(named.residual_helmholtz, rel=1e-12, abs=1e-12)
    assert positional.residual_chemical_potential_over_rt == pytest.approx(
        list(reversed(named.residual_chemical_potential_over_rt)), rel=1e-12, abs=1e-12
    )


def test_positional_problem_equals_scrambled_names() -> None:
    model = Mixture(Parameters.from_json(PACKET, components=ORDER))
    scrambled = ("sodium-cation", "chloride-anion", "water")
    liquids = [Phase(name, kind="liquid", support=scrambled) for name in ("L1", "L2")]
    named = Problem(phases=[*liquids, Phase("V", kind="vapor", support=("water",))], T=298.15, P=1.0e5,
                    feed=Amounts({"chloride-anion": 1.0, "water": 50.0, "sodium-cation": 1.0}))
    positional = _core.EquilibriumProblem()
    positional.temperature, positional.pressure = _core.FreeValue.fixed(298.15), _core.FreeValue.fixed(1.0e5)
    positional.feed.values = [50.0, 1.0, 1.0]
    phases = []
    for name, kind, support in (("L1", "Liquid", [0, 1, 2]), ("L2", "Liquid", [0, 1, 2]), ("V", "Vapor", [0])):
        phase = _core.PhaseSpec()
        phase.name, phase.kind, phase.support = name, getattr(_core.PhaseKind, kind), support
        phases.append(phase)
    positional.phases = phases
    # Each route uses its own model, so the charge row and z*phi terms carry each route's charges.
    systems = [_core.compile_problem(model._model, named._native(model)),
               _core.compile_problem(_positional_model(), positional)]
    assert [system.charge_row_count for system in systems] == [1, 1]
    point = systems[0].generic_point()
    for k in range(systems[0].rank + 1):
        point[systems[0].primal_size + k] = 0.3 + 0.2 * k
    rows = [system.evaluate(point).values for system in systems]
    assert rows[0] == pytest.approx(rows[1], rel=1e-12, abs=1e-12)
