"""Analytic conserved inventories in #237; these bounds do not impose equilibrium."""

import math

from epcsaft import Mixture
from epcsaft import equilibrium as eq
from tests.native_inputs import ETHANE, parameters, same


def compiled(feed, reaction, support=None):
    model = Mixture(parameters((ETHANE,) * len(feed)))
    problem = eq.Problem(phases=[eq.Phase("L", kind="liquid", support=support)],
                         T=300.0, P=1e5, feed=eq.Amounts(feed),
                         reactions=[eq.Reaction(reaction, equilibrium_constant=1.0)])
    return model, eq.compile_problem(model, problem)


def test_amount_fraction_and_feed_scaling(evidence):
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
                if fraction:
                    assert r.finite_witness[-1] > 0.0
                    assert abs(r.log_fraction - math.log(value)) <= 1e-8
                else:
                    assert abs(r.feed_normalized_amount - value/(3*scale)) <= 1e-10


def test_empty_finite_set_and_limiting_supremum(evidence):
    for feed in ({"s0": 1.0, "s1": 1.0}, {"s0": 2.0, "s1": 1.0}):
        model, c = compiled(feed, {"s0": 1.0, "s1": 1.0})
        amount = eq.inventory_bounds(c, [0, 1], False, 0.0)
        assert all(r.status == "unbounded_amount" and r.maximum is None for r in amount)
        r = eq.inventory_bounds(c, [1], True, 0.0)[0]
        same(evidence, f"inventory/signed/{feed['s0']}/supremum", r.maximum, 0.5,
             atol=1e-10, rtol=0.0, unit="1")
        assert r.finite_witness[-1] > 0.0
        assert r.limiting == (feed["s0"] != feed["s1"])
        assert max(r.balance_residual, r.normalization_residual, r.lower_bound_residual) <= 1e-10
        empty = eq.inventory_bounds(c, [1], True, 0.5)[0]
        if feed["s0"] != feed["s1"]:
            assert empty.status == "infeasible_inventory" and empty.maximum is None
        else:
            assert empty.maximum == 0.5 and empty.finite_witness[-1] > 0.0

