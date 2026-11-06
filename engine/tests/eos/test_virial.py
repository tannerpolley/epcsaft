"""Gross2001/2002 neutral virials, independently evaluated with teqp 0.23.2 and Clapeyron 0.6.29.

Hydrocarbons: teqp get_Bnvir(3) at zero density. Water 2B: ForwardDiff of
Clapeyron a_res in 256-bit arithmetic at 1e-20 and 1e-30 mol/m³ (identical
rounded coefficients), plus its public second_virial_coefficient. Parameters
are Gross2001 Table 2 and Gross2002 water, constant sigma, all k_ij=l_ij=0.
"""
import csv
from pathlib import Path

import pytest

from epcsaft import Mixture, Pure
from tests.native_inputs import aqueous_electrolyte, parameters

ROWS = list(csv.DictReader((Path(__file__).parents[1] / "oracle/neutral-virial.csv").open()))
SPECIES = {"Methane": (16.04, 1., 3.7039, 150.03), "Propane": (44.1, 2.002, 3.6184, 208.11),
           "Water2B": (18.01528, 1.0656, 3.0007, 366.51)}


def model(name: str) -> Mixture:
    names = name.split("+")
    pairs = [("k_ij", (0, 1), 0.)] if len(names) == 2 else []
    if name == "Water2B":
        pairs += [("association_energy_over_k", (0, 1), 2500.7), ("association_volume", (0, 1), .034868)]
    p = parameters([SPECIES[n] for n in names], pairs,
                   sites=([0, 0], [1, 1]) if name == "Water2B" else ([], []))
    return Pure(p) if len(names) == 1 else Mixture(p)


@pytest.mark.parametrize("row", ROWS, ids=lambda r: f"{r['model']}-{r['T']}")
def test_neutral_virials_match_independent_models(row: dict[str, str]) -> None:
    fluid = model(row["model"])
    v = fluid.virial_coefficients(float(row["T"]), None if isinstance(fluid, Pure) else [.5, .5])
    for value, ref, floor in [(v.second, float(row["B"]), 1e-12), (v.third, float(row["C"]), 1e-16),
                              (v.second, float(row["B_callable"]), 1e-12)]:
        assert abs(value - ref) <= max(1e-8 * abs(ref), floor)


def test_water_virial_remainder_and_admission() -> None:
    water = model("Water2B")
    v = water.virial_coefficients(300.)
    previous = float("inf")
    for rho in (1., .1, .01, .001):
        z = water.state(300., rho=rho).compressibility_factor
        remainder = abs((z - 1. - v.second * rho - v.third * rho**2) / rho**2)
        assert remainder < previous
        previous = remainder
    assert previous < 1e-10
    with pytest.raises(ValueError):
        water.state(300., rho=0.)
    for temperature in (0., -1., float("nan"), float("inf")):
        with pytest.raises(ValueError):
            water.virial_coefficients(temperature)
    mixture = model("Methane+Propane")
    face = mixture.virial_coefficients(300., [1., 0.])
    pure = model("Methane").virial_coefficients(300.)
    assert face.second == pure.second and face.third == pure.third
    for x in ([.2, .2], [1.1, -.1], [float("nan"), 0.]):
        with pytest.raises(ValueError):
            mixture.virial_coefficients(300., x)
    charged = Mixture(aqueous_electrolyte())
    for x in ([.98, .01, .01], [1., 0., 0.]):
        with pytest.raises(ValueError, match="charged declaration is nonanalytic"):
            charged.virial_coefficients(300., x)


@pytest.mark.parametrize("family", ["dipole", "quadrupole"])
def test_polar_third_virial_retains_three_body_term(family: str) -> None:
    # A finite-density verification of the second Taylor coefficient catches the
    # zero-center Padé branch that used to preserve only the two-body C term.
    if family == "dipole":
        from tests.eos.test_dipolar import _model
        names = ("acetone",)
    else:
        from tests.eos.test_quadrupolar import _model
        names = ("carbon-dioxide",)
    polar, unpolar = _model(names, True), _model(names, False)
    a, b = polar.virial_coefficients(300.), unpolar.virial_coefficients(300.)
    first, second = a.second - b.second, a.third - b.third
    errors = []
    for rho in (1., .1, .01):
        value = polar.state(300., rho=rho).dipolar
        estimate = 2. * (value / rho - first) / rho
        errors.append(abs(estimate / second - 1.))
    assert errors[-1] < 1e-4 and errors[-1] < errors[0] / 50.
    # An absent polar species must not activate an undefined quotient on a face.
    if family == "dipole":
        face = _model(("acetone", "n-pentane"), True).virial_coefficients(300., [0., 1.])
        plain = _model(("acetone", "n-pentane"), False).virial_coefficients(300., [0., 1.])
    else:
        face = _model(("carbon-dioxide", "methane"), True).virial_coefficients(300., [0., 1.])
        plain = _model(("carbon-dioxide", "methane"), False).virial_coefficients(300., [0., 1.])
    assert face.second == plain.second and face.third == plain.third
