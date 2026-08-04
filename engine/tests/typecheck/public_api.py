from __future__ import annotations

from typing import assert_type

import pint

from epcsaft import (
    Mixture,
    Parameters,
    SourceReferenceDeclaration,
    SourceReferenceTransfer,
    State,
    source_reference_transfer,
)
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
    assert_type(
        Parameters.from_dictionary(
            {
                "schema": "epcsaft.parameters",
                "schema_version": 1,
                "components": ("sample",),
                "parameters": {
                    "m": (1.0,),
                    "s": (3.0,),
                    "e": (200.0,),
                    "k_ij": ((0.0,),),
                },
                "options": {"permittivity_model": "none"},
                "validity": {"kind": "unknown"},
            }
        ),
        Parameters,
    )
    assert_type(
        Mixture.from_dictionary(
            {
                "schema": "epcsaft.parameters",
                "schema_version": 1,
                "components": ("sample",),
                "parameters": {
                    "m": (1.0,),
                    "s": (3.0,),
                    "e": (200.0,),
                    "k_ij": ((0.0,),),
                },
                "options": {"permittivity_model": "none"},
                "validity": {"kind": "unknown"},
            }
        ),
        Mixture,
    )
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
    declaration = SourceReferenceDeclaration(
        reference_state_id="source-v1",
        component_ids=model.component_ids,
        solvent_mole_fractions=tuple(1.0 for _ in model.component_ids),
        activity_convention_id="molality-infinite-dilution-v1",
        standard_molality_mol_per_kg=1.0,
        reference_pressure_pa=100000.0,
        required_derivatives=(),
    )
    assert_type(
        source_reference_transfer(
            model,
            declaration,
            T=temperature,
            P=pressure,
        ),
        SourceReferenceTransfer,
    )

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
    domain_kind: DomainKind = "fit-range"
    single: SingleParameterFamily = "segment_count"
    pair: PairParameterFamily = "k_ij"
    association: AssociationParameterFamily = "association_volume"
    model_family: ModelParameterFamily = "ionic_region_relative_permittivity"
    correlation: CorrelationFamily = "relative_permittivity"
    assert_type(domain_kind, DomainKind)
    assert_type(single, SingleParameterFamily)
    assert_type(pair, PairParameterFamily)
    assert_type(association, AssociationParameterFamily)
    assert_type(model_family, ModelParameterFamily)
    assert_type(correlation, CorrelationFamily)
