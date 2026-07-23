from __future__ import annotations

from typing import assert_type

import pint

from epcsaft import EPCSAFT, EosResult, ParameterBundle
from epcsaft.parameters import BundlePurpose
from epcsaft.records import (
    AssociationParameterFamily,
    CorrelationFamily,
    DomainKind,
    ModelParameterFamily,
    PairParameterFamily,
    SingleParameterFamily,
)


def check_public_types(
    model: EPCSAFT,
    temperature: pint.Quantity[float],
    pressure: pint.Quantity[float],
    molar_density: pint.Quantity[float],
) -> None:
    assert_type(ParameterBundle.catalog(), tuple[ParameterBundle, ...])
    assert_type(
        model.evaluate(
            temperature=temperature,
            pressure=pressure,
            mole_fractions=(1.0,),
            phase="vapor",
        ),
        EosResult,
    )
    assert_type(
        model.evaluate(
            temperature=temperature,
            molar_density=molar_density,
            mole_fractions=(1.0,),
        ),
        EosResult,
    )
    result = model.evaluate(
        temperature=temperature,
        pressure=pressure,
        mole_fractions=(1.0,),
    )
    assert_type(result.component_ids, tuple[str, ...])
    assert_type(result.mole_fractions, tuple[float, ...])
    assert_type(result.parameter_fingerprint, str)

    model.evaluate(  # type: ignore[call-overload]
        temperature=temperature,
        mole_fractions=(1.0,),
    )
    model.evaluate(  # type: ignore[call-overload]
        temperature=temperature,
        pressure=pressure,
        molar_density=molar_density,
        mole_fractions=(1.0,),
    )
    model.evaluate(  # type: ignore[call-overload]
        temperature=temperature,
        molar_density=molar_density,
        mole_fractions=(1.0,),
        phase="liquid",
    )
    model.evaluate(  # type: ignore[call-overload]
        temperature=temperature,
        pressure=pressure,
        mole_fractions=(1.0,),
        phase="solid",
    )
    model.reference_state(
        temperature=temperature,
        pressure=pressure,
        solvent_mole_fractions=(1.0,),
        phase="liquid",
    )
    purpose: BundlePurpose = "user-provided"
    domain_kind: DomainKind = "fit-range"
    single: SingleParameterFamily = "segment_count"
    pair: PairParameterFamily = "k_ij"
    association: AssociationParameterFamily = "association_volume"
    model_family: ModelParameterFamily = "ionic_region_relative_permittivity"
    correlation: CorrelationFamily = "relative_permittivity"
    assert_type(purpose, BundlePurpose)
    assert_type(domain_kind, DomainKind)
    assert_type(single, SingleParameterFamily)
    assert_type(pair, PairParameterFamily)
    assert_type(association, AssociationParameterFamily)
    assert_type(model_family, ModelParameterFamily)
    assert_type(correlation, CorrelationFamily)
