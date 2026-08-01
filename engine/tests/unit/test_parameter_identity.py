from __future__ import annotations

from dataclasses import replace

from epcsaft import Parameters
from epcsaft.parameters import _ParameterBundle
from epcsaft.records import (
    ConstantCorrelation,
    ModelParameterRecord,
    PairParameterRecord,
    SingleParameterRecord,
)

from tests._raw_cases import NEUTRAL_BINARY, neutral_binary


def test_parameter_fingerprint__depends_on_values_order_and_validity_not_provenance() -> None:
    original = Parameters.from_dictionary(NEUTRAL_BINARY)
    bundle = original._bundle
    renamed = _ParameterBundle.from_records(
        bundle_id="different-container",
        bundle_version=7,
        purpose="package-test-fixture",
        sources=(
            replace(
                bundle.sources[0],
                source_id="different-source",
                citation="Different source metadata",
            ),
        ),
        domains=(replace(bundle.domains[0], domain_id="different-domain"),),
        components=bundle.components,
        singles=(
            replace(
                record,
                record_id=f"renamed-{index}",
                source_id="different-source",
                locator=f"different locator {index}",
                domain_id="different-domain",
            )
            for index, record in enumerate(bundle.records)
            if record.__class__.__name__ == "SingleParameterRecord"
        ),
        pairs=(
            replace(
                record,
                record_id=f"renamed-pair-{index}",
                source_id="different-source",
                locator=f"different pair locator {index}",
                domain_id="different-domain",
            )
            for index, record in enumerate(bundle.records)
            if record.__class__.__name__ == "PairParameterRecord"
        ),
        models=(
            replace(
                record,
                record_id=f"renamed-model-{index}",
                source_id="different-source",
                locator=f"different model locator {index}",
                domain_id="different-domain",
            )
            for index, record in enumerate(bundle.records)
            if record.__class__.__name__ == "ModelParameterRecord"
        ),
    ).select(original.component_ids)

    assert renamed.fingerprint == original.fingerprint
    assert renamed._bundle.fingerprint != original._bundle.fingerprint
    assert (
        Parameters.from_dictionary(neutral_binary(("ethane", "methane"))).fingerprint
        != original.fingerprint
    )

    fit_range = neutral_binary(("methane", "ethane"))
    fit_range["validity"] = {"kind": "fit-range", "temperature_min_k": 200.0}
    reported = neutral_binary(("methane", "ethane"))
    reported["validity"] = {"kind": "reported-conditions", "temperature_min_k": 200.0}
    assert (
        Parameters.from_dictionary(fit_range).fingerprint
        != Parameters.from_dictionary(reported).fingerprint
    )


def test_parameter_containers__are_copied_to_immutable_tuples() -> None:
    mapping = neutral_binary(("methane", "ethane"))
    parameters = Parameters.from_dictionary(mapping)

    mapping["components"] = ("changed",)
    mapping["parameters"]["m"] = (999.0, 999.0)

    assert parameters.component_ids == ("methane", "ethane")
    assert isinstance(parameters.records, tuple)


def test_parameter_fingerprint__normalizes_constant_correlation_representation() -> None:
    direct = Parameters.from_dictionary(NEUTRAL_BINARY)
    bundle = direct._bundle
    diameter_records = tuple(
        record
        for record in bundle.records
        if isinstance(record, SingleParameterRecord) and record.family == "segment_diameter"
    )
    represented_as_correlations = _ParameterBundle.from_records(
        bundle_id="constant-correlation-representation",
        bundle_version=1,
        purpose="package-test-fixture",
        sources=bundle.sources,
        domains=bundle.domains,
        components=bundle.components,
        singles=(
            record
            for record in bundle.records
            if isinstance(record, SingleParameterRecord) and record.family != "segment_diameter"
        ),
        pairs=(record for record in bundle.records if isinstance(record, PairParameterRecord)),
        correlations=(
            ConstantCorrelation(
                record.record_id,
                record.component_id,
                "segment_diameter",
                record.value,
                record.source_id,
                record.locator,
                record.domain_id,
            )
            for record in diameter_records
        ),
        models=(record for record in bundle.records if isinstance(record, ModelParameterRecord)),
    ).select(direct.component_ids)

    assert represented_as_correlations.fingerprint == direct.fingerprint
