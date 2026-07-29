from __future__ import annotations

from typing import assert_type

import pint

from epcsaft import Mixture, Parameters, State
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
    model: Mixture,
    parameters: Parameters,
    temperature: pint.Quantity[float],
    pressure: pint.Quantity[float],
    molar_density: pint.Quantity[float],
) -> None:
    assert_type(parameters.component_ids, tuple[str, ...])
    assert_type(
        model.state(
            T=temperature,
            P=pressure,
            x=(1.0,),
            phase="vapor",
        ),
        State,
    )
    assert_type(
        model.state(
            T=temperature,
            rho=molar_density,
            x=(1.0,),
        ),
        State,
    )
    result = model.state(
        T=temperature,
        P=pressure,
        x=(1.0,),
    )
    assert_type(result.component_ids, tuple[str, ...])
    assert_type(result.mole_fractions, tuple[float, ...])
    assert_type(result.parameter_fingerprint, str)

    model.state(  # type: ignore[call-overload]
        T=temperature,
        x=(1.0,),
    )
    model.state(  # type: ignore[call-overload]
        T=temperature,
        P=pressure,
        rho=molar_density,
        x=(1.0,),
    )
    model.state(  # type: ignore[call-overload]
        T=temperature,
        rho=molar_density,
        x=(1.0,),
        phase="liquid",
    )
    model.state(  # type: ignore[call-overload]
        T=temperature,
        P=pressure,
        x=(1.0,),
        phase="solid",
    )
    model.reference(
        T=temperature,
        P=pressure,
        solvent_x=(1.0,),
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
