"""Issue #275: the retained cryogenic mass-action failure must resolve without a skip."""
from typing import Any

from epcsaft import Mixture, trace

from ..native_inputs import parameters, same


def test_ethanol_butane_cold_scan_resolves_pinned_node(evidence: Any) -> None:
    model = Mixture(parameters(
        [(46.06844, 2.3827, 3.1771, 198.24), (58.1222, 2.3316, 3.7086, 222.88)],
        [("k_ij", (0, 1), 0.028), ("association_energy_over_k", (0, 1), 2653.4),
         ("association_volume", (0, 1), 0.032384)], ((0, 0), (1, 1))))
    # Chapman1990 Eq. 15, equal complementary populations: X=2/(1+sqrt(1+4q)).
    with trace(diagnostics=True) as node:
        model.state(45.0, rho=2709.4906792105835, x=(0.5, 0.5))
    fractions = next(event.values for event in node.events if event.name == "site_fractions")
    for index, fraction in enumerate(fractions):
        same(evidence, f"cold-node/site/{index}", fraction, 4.8632143028761764e-12,
             atol=0.0, rtol=1e-12)
    assert all(span.convergence.residual <= 2e-13 for span in node.spans
               if span.name == "association_newton" and span.convergence is not None)
    with trace(diagnostics=True) as scan:
        state = model.state(45.0, P=1e-4, x=(0.5, 0.5), phase="liquid")
    assert scan.trace is not None and not scan.trace.incomplete
    assert not any(event.name == "density_scan_skip" for event in scan.events)
    assert sum(event.name == "density_scan_point" for event in scan.events) == 41
    same(evidence, "cold-scan/density", state.molar_density, 5.3454360146661902e-7,
         atol=0.0, rtol=1e-12)
