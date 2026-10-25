"""Issue #275: the liquid-branch density scan of equimolar ethanol/butane at 45 K and 1e-4 Pa.
One scan node, 2709.49 mol/m3, once failed the ethanol association mass action; its site
fractions must equal Chapman's closed form at the original fixed-density node."""
import math
from typing import Any

from epcsaft import Mixture, trace

from tests.native_inputs import parameters, same


def test_ethanol_butane_cold_scan_resolves_pinned_node(evidence: Any) -> None:
    model = Mixture(parameters(
        [(46.06844, 2.3827, 3.1771, 198.24), (58.1222, 2.3316, 3.7086, 222.88)],
        [("k_ij", (0, 1), 0.028), ("association_energy_over_k", (0, 1), 2653.4),
         ("association_volume", (0, 1), 0.032384)], ((0, 0), (1, 1))))
    # Chapman1990 Eqs. 15 and 18 for ethanol's two sites (sigma^3 convention,
    # engine/docs/science/full-eos.md): X = 2/(1+sqrt(1+4q)), q = rho x_0 Delta,
    # Delta = g_hs(d_0/2) sigma_0^3 kappa (exp(epsilon/T) - 1), d = sigma(1 - 0.12 exp(-3 eps_i/T)).
    temperature, density, x = 45.0, 2709.4906792105835, (0.5, 0.5)
    segments, sigma, energy = (2.3827, 2.3316), (3.1771, 3.7086), (198.24, 222.88)
    d = [s * (1 - 0.12 * math.exp(-3 * e / temperature)) for s, e in zip(sigma, energy)]
    number_density = density * 6.02214076e23 * 1e-30  # per cubic angstrom
    zeta2, zeta3 = (math.pi / 6 * number_density * sum(xi * mi * di**n for xi, mi, di in zip(x, segments, d))
                    for n in (2, 3))
    contact = (1 / (1 - zeta3) + d[0] / 2 * 3 * zeta2 / (1 - zeta3) ** 2
               + (d[0] / 2) ** 2 * 2 * zeta2**2 / (1 - zeta3) ** 3)
    q = number_density * x[0] * contact * sigma[0] ** 3 * 0.032384 * math.expm1(2653.4 / temperature)
    chapman = 2 / (1 + math.sqrt(1 + 4 * q))
    with trace(diagnostics=True) as node:
        model.state(temperature, rho=density, x=x)
    fractions = next(event.values for event in node.events if event.name == "site_fractions")
    for index, fraction in enumerate(fractions):
        same(evidence, f"cold-node/site/{index}", fraction, chapman, atol=0.0, rtol=1e-12)
    assert all(span.convergence.residual <= 2e-13 for span in node.spans
               if span.name == "association_newton" and span.convergence is not None)
