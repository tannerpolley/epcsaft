"""Analytic conserved inventories in #237; these bounds do not impose equilibrium."""

import math
from itertools import combinations
from typing import Any

from epcsaft import Mixture, _core
from epcsaft import equilibrium as eq
from tests.native_inputs import ETHANE, parameters, same


def compiled(feed: dict[str, float], reaction: dict[str, float], support: list[str] | None = None
             ) -> tuple[Mixture, _core.CompiledSystem]:
    model = Mixture(parameters((ETHANE,) * len(feed),
                    [("k_ij", pair, 0.0) for pair in combinations(range(len(feed)), 2)]))
    problem = eq.Problem(phases=[eq.Phase("L", kind="liquid", support=support)],
                         T=300.0, P=1e5, feed=eq.Amounts(feed),
                         reactions=[eq.Reaction(reaction, equilibrium_constant=1.0)])
    return model, eq.compile_problem(model, problem)


def test_amount_fraction_and_feed_scaling(evidence: Any) -> None:
    # n=(2-z,1-z,z), 0<=z<=1; a fraction floor shrinks this interval.
    for scale in (1.0, 7.0):
        model, c = compiled({"s0": 2*scale, "s1": scale, "s2": 0.0},
                            {"s0": -1.0, "s1": -1.0, "s2": 1.0})
        for fraction, floor, expected in (
            (False, 0.0, (2*scale, scale, scale)),
            (True, 0.0, (2/3, 1/3, 1/2)),
            (True, 0.01, ((2-0.01)/3, (1-0.02)/3, (1-0.03)/2)),
        ):
            for r, value in zip(eq.inventory_bounds(c, [0, 1, 2], fraction, floor), expected, strict=True):
                assert r.status == "optimal", r.termination
                same(evidence, f"inventory/{scale}/{fraction}/{floor}/{r.species}", r.maximum,
                     value, atol=1e-10, rtol=0.0, unit="1" if fraction else "mol")
                assert max(r.balance_residual, r.normalization_residual, r.lower_bound_residual) <= 1e-10
                assert len(r.amounts_mol) == 3
                evidence.observe(f"inventory/{scale}/{fraction}/{floor}/{r.species}/diagnostics",
                                 balance=r.balance_residual, normalization=r.normalization_residual,
                                 lower_bound=r.lower_bound_residual, maximizing_witness=r.witness,
                                 amounts_mol=r.amounts_mol, finite_witness=r.finite_witness,
                                 cause=r.termination, lp_error=r.lp_error)
                if fraction:
                    assert r.finite_witness[-1] > 0.0
                    assert r.log_fraction is not None
                    assert abs(r.log_fraction - math.log(value)) <= 1e-8
                else:
                    assert r.feed_normalized_amount is not None
                    assert abs(r.feed_normalized_amount - value/(3*scale)) <= 1e-10


def test_empty_finite_set_and_limiting_supremum(evidence: Any) -> None:
    for feed in ({"s0": 1.0 + delta, "s1": 1.0} for delta in (0.0, 1.0, 1e-10, 1e-12)):
        model, c = compiled(feed, {"s0": 1.0, "s1": 1.0})
        amount = eq.inventory_bounds(c, [0, 1], False, 0.0)
        assert all(r.status == "unbounded_amount" and r.maximum is None for r in amount)
        evidence.observe(f"inventory/signed/{feed['s0']}/amounts", causes=[r.termination for r in amount])
        bounds = eq.inventory_bounds(c, [0, 1], True, 0.0)
        if feed["s0"] == feed["s1"]:
            same(evidence, "inventory/signed/zero/species0", bounds[0].maximum, 0.5,
                 atol=1e-10, rtol=0.0, unit="1")
        r = bounds[1]
        same(evidence, f"inventory/signed/{feed['s0']}/supremum", r.maximum, 0.5,
             atol=1e-10, rtol=0.0, unit="1")
        assert r.finite_witness[-1] > 0.0
        if feed["s0"] != feed["s1"]:
            assert r.limiting and not r.amounts_mol
        else:
            # tau=0 can also maximize a fraction that is attained at finite amount.
            same(evidence, "inventory/signed/zero/finite_fraction", r.finite_witness[1], 0.5,
                 atol=1e-10, rtol=0.0, unit="1")
        evidence.observe(f"inventory/signed/{feed['s0']}/witnesses", maximizing_y_tau=r.witness,
                         finite_y_tau=r.finite_witness, cause=r.termination, limiting=r.limiting,
                         balance=r.balance_residual, normalization=r.normalization_residual,
                         lower_bound=r.lower_bound_residual, finite_residual=r.finite_residual)
        assert max(r.balance_residual, r.normalization_residual, r.lower_bound_residual) <= 1e-10
        empty = eq.inventory_bounds(c, [1], True, 0.5)[0]
        if feed["s0"] != feed["s1"]:
            assert empty.status == "infeasible_inventory" and empty.maximum is None
        else:
            assert empty.maximum == 0.5 and empty.finite_witness[-1] > 0.0
        evidence.observe(f"inventory/signed/{feed['s0']}/half_floor", disposition=empty.status,
                         maximizing_y_tau=empty.witness, finite_y_tau=empty.finite_witness,
                         cause=empty.feasibility_termination, finite_residual=empty.finite_residual,
                         lp_error=empty.lp_error)
