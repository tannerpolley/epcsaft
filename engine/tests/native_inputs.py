"""Explicit inputs of the retained native fixtures, restated through the public package (#59).

`ROWS` holds the C++ values recorded by oracle/record_native_rows.py for these same inputs.
"""
import hashlib
import itertools
import math
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from epcsaft import (
    GAS_CONSTANT_J_PER_MOL_K as R,
    IdealCorrelation,
    IdealDIPPR107,
    IdealInterval,
    IdealNASA7,
    IdealNASA9,
    IdealPolynomial,
    IdealShomate,
    Mixture,
    Parameters,
    ThermochemistryRecord,
    _core,
)
from epcsaft.equilibrium import (
    Amounts, Free, NeutralReference, Phase, Pinned, Problem, Reaction, ReactionLogPolynomial, ReactionReference,
    ReferenceSourceBasis,
)

ROWS = {
    line.split("\t")[1]: float(line.split("\t")[3])
    for line in (Path(__file__).parent / "oracle" / "native-rows.tsv").read_text().splitlines()[1:]
}
def same(evidence: Any, name: str, model: float | None, reference: float | None = None,
         atol: float = 1e-12, rtol: float = 1e-12, **inputs: Any) -> None:
    """Record and assert |model - reference| <= max(atol, rtol max(|model|, |reference|)).

    The default is the retained same-owner transport rule against the recorded C++ row `name`,
    without the absolute floor when that row lies below 1e-9.
    """
    if reference is None:  # a recorded row below the absolute floor must still match relatively
        reference, atol = ROWS[name], 0.0 if abs(ROWS[name]) < 1e-9 else atol
    evidence.compare(name, model=model, reference=reference, **inputs,
                     criterion={"type": "scalar", "rule": "symmetric_close", "atol": atol, "rtol": rtol})
    assert model is not None and math.isclose(model, reference, rel_tol=rtol, abs_tol=atol), name


def packet_metadata(path: Path, components: Sequence[str]) -> dict[str, Any]:
    """Evidence identity of data/packets/<id>/<version>/parameters/parameters.json: id, version and hash."""
    return {
        "packet_id": path.parts[-4],
        "packet_version": int(path.parts[-3]),
        "parameter_file": {
            "path": str(path.relative_to(path.parents[5])),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "size_bytes": path.stat().st_size,
        },
        "selected_components": list(components),
    }


def parameters(
    species: Sequence[tuple[float, float, float, float]],  # g/mol, segments, sigma / angstrom, epsilon / K
    pairs: Sequence[tuple[str, tuple[int, ...], float]] = (),  # k_ij and association site-pair entries
    sites: tuple[Sequence[int], Sequence[int]] = ((), ()),
    charges: Sequence[int] | None = None,
    extra: Sequence[tuple[str, tuple[int, ...], float, str]] = (),
    choices: tuple[str, str, str] = ("none", "", "none"),  # electrolyte, ionic dispersion, permittivity
) -> Parameters:
    document = _core.ParameterDocument()
    document.component_ids = [f"s{i}" for i in range(len(species))]
    document.site_component, document.site_multiplicity = list(sites[0]), list(sites[1])
    document.base_choice = "pc-saft"
    document.association_choice = "general-site" if sites[0] else "none"
    document.electrolyte_choice, document.ionic_dispersion_choice, document.permittivity_choice = choices
    rows: list[tuple[str, tuple[int, ...], float, str]] = list(extra)
    for i, (mass, segments, sigma, epsilon) in enumerate(species):
        rows += [("molar_mass", (i,), mass, "gram / mole"), ("segment_count", (i,), segments, "dimensionless"),
                 ("sigma_constant", (i,), sigma, "angstrom"), ("dispersion_energy_over_k", (i,), epsilon, "kelvin")]
        rows.append(("charge_number", (i,), float(charges[i]) if charges else 0.0, "elementary-charge"))
    rows += [(field, index, value, "kelvin" if field == "association_energy_over_k" else "dimensionless")
             for field, index, value in pairs]
    entries = []
    for field, index, magnitude, unit in rows:
        item = _core.ParameterEntry()
        item.identity, item.field, item.index = f"{field}/{index}", field, list(index)
        item.magnitude, item.unit = magnitude, unit
        entries.append(item)
    document.entries = entries
    site_ids = tuple((document.component_ids[c], f"site{k}") for k, c in enumerate(sites[0]))
    return Parameters(_core.build_parameters(document), tuple(document.component_ids), site_ids)


ETHANE = (30.07, 1.6069, 3.5206, 191.42)


def budget() -> Any:
    """The former implicit fit budget, now explicit: 100 iterations and 600 s."""
    return _core.FitControls(maximum_iterations=100, maximum_elapsed_time_seconds=600.0)


def aqueous_electrolyte(*more: tuple[str, tuple[int, ...], float, str], c_shell: float = 0.0) -> Parameters:
    """properties_support::aqueous_electrolyte(): water, Na+ and Cl- with Debye-Huckel and Born terms."""
    extra: list[tuple[str, tuple[int, ...], float, str]] = [*more,
        ("relative_permittivity", (0,), 78.09, "dimensionless"), ("ionic_region_relative_permittivity", (), 8.0, "dimensionless"),
        ("born_c_shell", (), c_shell, "dimensionless"), ("born_c_dielectric", (), 0.0, "dimensionless"),
        ("ion_suppression_coefficient", (), 7.01, "dimensionless"),
        ("sigma_term_amplitude", (0, 0), 10.11, "angstrom"), ("sigma_term_exponent", (0, 0), -0.01775, "1 / kelvin"),
        ("sigma_term_amplitude", (1, 0), -1.417, "angstrom"), ("sigma_term_exponent", (1, 0), -0.01146, "1 / kelvin")]
    for i, (solvation, diameter, born) in enumerate(((1.5, None, None), (1.0, 2.484416, 3.445), (1.0, 2.42528, 4.1))):
        extra.append(("solvation_factor", (i,), solvation, "dimensionless"))
        if diameter is not None and born is not None:
            extra += [("packing_diameter", (i,), diameter, "angstrom"), ("debye_huckel_diameter", (i,), diameter, "angstrom"),
                      ("born_diameter", (i,), born, "angstrom")]
    pairs = [("k_ij", (0, 1), -0.3), ("k_ij", (0, 2), 0.0), ("k_ij", (1, 2), 0.8),
             ("association_energy_over_k", (0, 1), 2425.7), ("association_volume", (0, 1), 0.04509)]
    species = ((18.0153, 1.2047, 2.7927, 353.95), (22.98977, 1.0, 2.8232, 230.0), (35.4527, 1.0, 2.756, 170.0))
    return parameters(species, pairs, ((0, 0), (1, 1)), (0, 1, -1), extra,
                      ("born", "exclude-same-sign-ion-pairs", "ion-fraction-suppression"))


def associating_binary(shift: float = 0.0) -> Parameters:
    """properties_support::associating_binary(); `shift` / K moves the reciprocal site pair (0, 1) energy."""
    edges = (((0, 1), 2383.4744 + shift, 0.118488), ((0, 3), 2093.8222, 0.1540467874409226),
             ((1, 2), 2093.8222, 0.1540467874409226), ((2, 3), 1804.17, 0.2039))
    pairs: list[tuple[str, tuple[int, int], float]] = [("k_ij", (0, 1), -0.052)]
    for index, energy, volume in edges:
        pairs += [("association_energy_over_k", index, energy), ("association_volume", index, volume)]
    return parameters(((61.08, 4.5354, 2.6019, 204.0438), (18.01528, 2.1945, 2.229, 141.66)), pairs,
                      ((0, 0, 1, 1), (2, 1, 2, 2)))


def p3_forms() -> dict[str, Any]:
    """properties_support::p3_forms(): cp = 30 J/(mol K), h(300 K) = 1000 J/mol, s(300 K) = 50 J/(mol K)."""
    nasa = (1000.0 / R - 30.0 * 300.0 / R, 50.0 / R - 30.0 * math.log(300.0) / R)
    return {
        "polynomial": IdealPolynomial([30.0, 0.0], 300.0, 1000.0, 50.0),
        "shomate": IdealShomate([30.0, 0.0, 0.0, 0.0, 0.0, -9.0, 50.0 - 30.0 * math.log(0.3), 0.0], 1000.0),
        "nasa7": IdealNASA7([30.0 / R, 0.0, 0.0, 0.0, 0.0, *nasa]),
        "nasa9": IdealNASA9([0.0, 0.0, 30.0 / R, 0.0, 0.0, 0.0, 0.0, *nasa]),
        "dippr107": IdealDIPPR107([0.0, 30.0, 0.0, 0.0, 0.0], 300.0, 1000.0, 50.0),
    }


def record(*forms: Any, lower: float = 250.0, upper: float = 350.0) -> ThermochemistryRecord:
    """One closed interval per component at 1 bar, in component order."""
    return ThermochemistryRecord(1.0e5, [[IdealInterval(lower, upper, True, True, IdealCorrelation(form, 1.0e5))]
                                         for form in forms])


def surface_model(**inputs: Any) -> Mixture:
    """The C++ example's model: constant cp with anchors at 300 K and 1 bar on [250, 450] K."""
    forms = (IdealPolynomial([80.0, 0.0], 300.0, 1000.0, 50.0), IdealPolynomial([35.0, 0.0], 300.0, 2000.0, 70.0))
    return Mixture(associating_binary(), thermochemistry=record(*forms, upper=450.0), **inputs)


def bubble(T: float) -> Problem:
    """equilibrium_support::boundary(): L free, V pinned 0, P free from 1e5 Pa, equimolar amounts."""
    return Problem(phases=[Phase("L", kind="liquid"), Phase("V", kind="vapor", amount=Pinned(0.0))],
                   T=T, P=Free(1.0e5), feed=Amounts({"s0": 0.5, "s1": 0.5}))


def mixed_difference(
    evaluate: Callable[[float, float, tuple[float, ...]], float],
    temperature: float,
    coordinate: float,
    composition: tuple[float, ...],
    directions: Sequence[tuple[float, float, tuple[float, ...]]],
    step: float,
) -> float:
    total = 0.0
    for signs in itertools.product((-1.0, 1.0), repeat=len(directions)):
        varied_temperature = temperature + step * sum(
            sign * direction[0] for sign, direction in zip(signs, directions)
        )
        varied_coordinate = coordinate + step * sum(
            sign * direction[1] for sign, direction in zip(signs, directions)
        )
        varied_composition = tuple(
            value + step * sum(
                sign * direction[2][index]
                for sign, direction in zip(signs, directions)
            )
            for index, value in enumerate(composition)
        )
        total += math.prod(signs) * evaluate(
            varied_temperature, varied_coordinate, varied_composition
        )
    return total / (2.0 * step) ** len(directions)


# Reference-calorics fixture: identical ethane 2A <=> B; only A carries an ideal-gas record,
# cp = 52 + 0.12 (T - 300) J/(mol K) anchored at h = s = 0 at 300 K.
DIMER_LAW = ReactionLogPolynomial(math.log(0.3), 200.0, 3.0, -0.004, 300.0, 250.0, 350.0, "analytic-source")
CP_A = (52.0, 0.12)
A_RECORD = [IdealInterval(200.0, 400.0, True, True, IdealCorrelation(IdealPolynomial(CP_A, 300.0, 0.0, 0.0), 1.0e5))]


def records(*species: list[IdealInterval]) -> ThermochemistryRecord:
    return ThermochemistryRecord(1.0e5, list(species))


def twin_ethane(*species: list[IdealInterval]) -> Mixture:
    return Mixture(parameters((ETHANE, ETHANE), [("k_ij", (0, 1), 0.0)]), thermochemistry=records(*species))


def dimer_problem(temperature: float, reference_pressure: float | None,
                  feed: dict[str, float] | None = None, pressure: float = 3.0e6) -> Problem:
    return Problem(
        T=temperature, P=pressure, feed=Amounts(feed or {"s0": 2.0, "s1": 0.2}), phases=[Phase("L", kind="liquid")],
        neutral_reference=NeutralReference(solvent="s0", solute_molality_ratios={"s1": 1.0},
                                           standard_molality_mol_per_kg=1.0),
        reactions=[Reaction({"s0": -2.0, "s1": 1.0}, correlation=DIMER_LAW, reference=ReactionReference(
            ReferenceSourceBasis.RawMoleFractionInfiniteDilution, reference_pressure))],
    )
