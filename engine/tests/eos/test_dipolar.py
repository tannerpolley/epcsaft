"""#28: optional Gross-Vrabec (2006) dipole-dipole term (ADR 0011, full EOS 4.1).

a_dd is isolated as a^r(polar family) - a^r(same parameters without it) at fixed
(T, rho, x). The checks falsify the equations against teqp 0.23.2 rows, the
second-virial limit written independently here from GV2006 Eq. 8 and Table 1,
the derivative chains, and the absent-dipole limit.
"""

import copy
import csv
import json
import math
from pathlib import Path
from typing import Any

from epcsaft import DensityActionDirection, Mixture, Parameters, PropertyActionStatus, PropertyObservable, Pure
from tests.native_inputs import mixed_difference
from tests.packets import parameters_json

ROWS = Path(__file__).resolve().parents[3] / "analyses/2026-polar-dipole/results/teqp-dipolar-rows.csv"
AVOGADRO, BOLTZMANN, EPS0, DEBYE = 6.02214076e23, 1.380649e-23, 8.8541878128e-12, 1.0e-21 / 299792458.0


def _document(names: tuple[str, ...], polar: bool) -> dict[str, Any]:
    """GV2006 Table 2 records plus Gross-Sadowski 2001 n-pentane (mu = 0), all k_ij = 0."""
    gv = json.loads(parameters_json("gross-vrabec-2006", 1).read_text())
    gs = json.loads(parameters_json("gross-sadowski-2001", 3).read_text())
    pentane = copy.deepcopy(next(c for c in gs["components"] if c["component_id"] == "pentane"))
    pentane["component_id"] = "n-pentane"
    pentane["coefficients"] += [
        {"identity": f"n-pentane/{family}", "family": family, "value": {"magnitude": 0.0, "unit": unit}}
        for family, unit in (("dipole_moment", "debye"), ("dipole_count", "dimensionless"))]
    document: dict[str, Any] = copy.deepcopy(gv)
    document["components"] = [c for c in gv["components"] if c["component_id"] in names] + [pentane]
    document["pairs"] = [
        {"component_id_a": a, "component_id_b": b,
         "coefficients": [{"identity": f"{a}/{b}/k_ij", "family": "k_ij", "value": {"magnitude": 0.0, "unit": "dimensionless"}}]}
        for i, a in enumerate(names) for b in names[i + 1:]]
    if not polar:
        document["model_families"] = [f for f in document["model_families"] if f["kind"] != "polar"]
        for component in document["components"]:
            component["coefficients"] = [
                c for c in component["coefficients"] if c["family"] not in ("dipole_moment", "dipole_count")]
    return document


def _model(names: tuple[str, ...], polar: bool) -> Mixture | Pure:
    parameters = Parameters.from_mapping(_document(names, polar), components=names)
    return Pure(parameters) if len(names) == 1 else Mixture(parameters)


def _residual(model: Mixture | Pure, temperature: float, density: float, x: tuple[float, ...]) -> float:
    if isinstance(model, Pure):
        return model.state(temperature, rho=density).residual_helmholtz
    return model.state(temperature, rho=density, x=x).residual_helmholtz


def _dipolar(names: tuple[str, ...], temperature: float, density: float, x: tuple[float, ...]) -> float:
    return (_residual(_model(names, True), temperature, density, x)
            - _residual(_model(names, False), temperature, density, x))


def test_dipolar_helmholtz_matches_teqp() -> None:
    # teqp received mu*^2 from the Engine's SI constants; both codes then agree to ~1e-15.
    with ROWS.open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 8
    for row in rows:
        names = tuple(row["species"].split("+"))
        x = tuple(float(v) for v in row["x"].split())
        value = _dipolar(names, float(row["T_K"]), float(row["rho_mol_m3"]), x)
        reference = float(row["a_dipolar_teqp"])
        assert abs(value - reference) <= 1e-12 * abs(reference), row


def test_dipolar_second_virial_limit() -> None:
    # Eq. 8 at eta = 0: a_dd / rho_N -> -pi sigma^3 q^2 J2(0), with m capped at 2 (Eq. 15).
    record = {c["family"]: c["value"]["magnitude"] for c in next(
        c for c in json.loads(parameters_json("gross-vrabec-2006", 1).read_text())["components"]
        if c["component_id"] == "acetone")["coefficients"]}
    m, sigma, epsilon, mu = (record[k] for k in ("segment_count", "segment_diameter", "dispersion_energy_over_k", "dipole_moment"))
    temperature, chain = 300.0, min(record["segment_count"], 2.0)
    first, second = (chain - 1.0) / chain, (chain - 1.0) * (chain - 2.0) / chain**2
    j2 = (0.3043504 + first * 0.9534641 + second * -1.1610080
          + (0.2187939 + first * -0.5873164 + second * 3.4869576) * epsilon / temperature)
    q = mu**2 * DEBYE**2 / (4.0 * math.pi * EPS0 * BOLTZMANN * temperature * m * (sigma * 1e-10) ** 3)
    virial = -math.pi * sigma**3 * q * q * j2
    errors = []
    for density in (1.0e-1, 1.0, 10.0):
        number_density = density * AVOGADRO * 1.0e-30
        value = _dipolar(("acetone",), temperature, density, (1.0,)) / number_density
        errors.append(value / virial - 1.0)
    # The limit holds and the deviation is O(rho) (A3/A2 and J2(eta) - J2(0) are linear in rho).
    intercept = errors[0] - (errors[1] - errors[0]) / 9.0  # linear extrapolation to rho = 0
    assert abs(intercept) < 1e-8 and abs(errors[0]) > 1e-6
    for coarse, fine in zip(errors[1:], errors[:-1]):
        assert math.isclose(coarse / fine, 10.0, rel_tol=2e-2)


def test_dipolar_derivative_actions_and_identities() -> None:
    names, temperature, density, x = ("acetone", "n-pentane"), 372.7, 9000.0, (0.4, 0.6)
    model = _model(names, True)
    assert isinstance(model, Mixture)
    state = model.state(temperature, rho=density, x=x)
    raw = ((1.0, 0.0, (0.0, 0.0)), (0.0, 20.0, (0.0, 0.0)), (0.0, 0.0, (0.02, -0.02)))
    directions = tuple(DensityActionDirection(*item) for item in raw)
    nonpolar = _model(names, False)
    assert isinstance(nonpolar, Mixture)

    def enthalpy(t: float, rho: float, z: tuple[float, ...]) -> float:
        return model.state(t, rho=rho, x=z).residual_enthalpy

    # Richardson (h, h/2) mixed stencils; predeclared relative criteria 1e-8, 1e-6, 1e-4.
    for order, tolerance in ((1, 1e-8), (2, 1e-6), (3, 1e-4)):
        batch = model.density_actions(temperature, density, x=x,
                                      observables=(PropertyObservable.ResidualEnthalpy,), directions=directions[:order])
        result = batch.results[0]
        assert result.status == PropertyActionStatus.Available and result.action is not None
        step = 0.2
        coarse = mixed_difference(enthalpy, temperature, density, x, raw[:order], step)
        fine = mixed_difference(enthalpy, temperature, density, x, raw[:order], step / 2)
        reference = (4.0 * fine - coarse) / 3.0
        assert abs(result.action - reference) <= tolerance * abs(reference), (order, result.action, reference)
        without = nonpolar.density_actions(temperature, density, x=x,
                                           observables=(PropertyObservable.ResidualEnthalpy,), directions=directions[:order])
        # The dipolar share of each action exceeds its criterion by many orders (non-vacuous).
        assert without.results[0].action is not None
        assert abs(result.action - without.results[0].action) > 1e2 * tolerance * abs(reference)

    # Z - 1 = rho a_rho and u_res = -R T^2 a_T through the dipolar chain.
    def slope(f: Any, center: float, h: float) -> float:
        return float(4.0 * (f(center + h / 2) - f(center - h / 2)) / h - (f(center + h) - f(center - h)) / (2 * h)) / 3.0

    a_rho = slope(lambda r: model.state(temperature, rho=r, x=x).residual_helmholtz, density, 1e-3 * density)
    assert abs(state.compressibility_factor - 1.0 - density * a_rho) < 1e-9
    a_t = slope(lambda t: model.state(t, rho=density, x=x).residual_helmholtz, temperature, 1e-3 * temperature)
    u = -8.31446261815324 * temperature**2 * a_t
    assert abs(state.residual_internal_energy - u) < 1e-8 * abs(u)


def test_absent_dipolar_species_limit() -> None:
    # With acetone absent the declared family adds nothing, and a_dd ~ A2 is quadratic in x_acetone,
    # so the dipolar share of the acetone chemical potential vanishes at infinite dilution.
    names, temperature, density = ("acetone", "n-pentane"), 330.0, 8000.0
    polar, nonpolar = _model(names, True), _model(names, False)
    assert isinstance(polar, Mixture) and isinstance(nonpolar, Mixture)
    with_family = polar.state(temperature, rho=density, x=(0.0, 1.0))
    without = nonpolar.state(temperature, rho=density, x=(0.0, 1.0))
    assert with_family.residual_helmholtz == without.residual_helmholtz
    assert with_family.pressure == without.pressure
    assert with_family.residual_chemical_potential_over_rt == without.residual_chemical_potential_over_rt
    small = [_dipolar(names, temperature, density, (s, 1.0 - s)) for s in (1e-4, 2e-4)]
    assert math.isclose(small[1] / small[0], 4.0, rel_tol=1e-3)
