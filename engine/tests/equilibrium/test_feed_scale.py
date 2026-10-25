"""#243 feed homogeneity: the water/isobutanol split is the same intensive state at any feed total."""
from epcsaft import Mixture, Parameters
from epcsaft import equilibrium as eq
from tests.packets import parameters_json


def test_two_liquid_split_is_feed_homogeneous() -> None:
    """At 100 mol the solve coalesced into one liquid 0.0216 RT per mol above the split when Ipopt
    measured amounts in mol; each free amount now enters the optimizer in units of its start amount."""
    model = Mixture(Parameters.from_json(
        parameters_json("khudaida-2026-figure-2-electrolyte-lle", 7), components=("water", "isobutanol")))

    def solve(total: float) -> eq.EquilibriumResult:
        result = eq.solve_equilibrium(model, eq.Problem(
            T=293.15, P=1e5, row_tolerance=1e-10, feed=eq.Amounts({"water": .6 * total, "isobutanol": .4 * total}),
            phases=[eq.Phase("Lw", kind="liquid", composition_guess=(.999, .001)),
                    eq.Phase("Lo", kind="liquid", composition_guess=(.3, .7))]))
        assert result.success and result.requested_tolerance_met, result.message
        return result

    unit, large = solve(1.), solve(100.)
    assert abs(unit.mole_fractions[0] - unit.mole_fractions[2]) > .5  # distinct liquids
    assert max(abs(a - b) for a, b in zip(large.mole_fractions, unit.mole_fractions, strict=True)) <= 1e-9
    assert max(abs(a / 100. - b) for a, b in zip(large.amounts, unit.amounts, strict=True)) <= 1e-12
    assert max(abs(a / b - 1.) for a, b in zip(large.molar_densities, unit.molar_densities, strict=True)) <= 1e-9
    assert abs(large.objective / 100. - unit.objective) <= 1e-10 * max(1., abs(unit.objective))
