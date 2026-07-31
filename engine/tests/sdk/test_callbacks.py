from __future__ import annotations

from ._support import *


def test_resolved_model_artifact_is_native_loadable(tmp_path: Path) -> None:
    model = _mixture_model()
    artifact = tmp_path / "model.json"

    export_native_model(model, artifact)

    payload = artifact.read_text(encoding="utf-8")
    assert '"schema":"epcsaft-resolved-native-model-v1"' in payload
    assert model.parameter_fingerprint in payload


@pytest.mark.parametrize(
    "topology,component,temperature_k,density_mol_m3",
    (
        ("neutral", "methane", 190.0, 1000.0),
        ("neutral", "ethane", 250.0, 1000.0),
        ("neutral", "propane", 300.0, 1000.0),
    ),
    ids=("methane", "ethane", "propane"),
)
def test_native_sdk_phase_block_matches_public_thermodynamics(
    topology: str,
    component: str,
    temperature_k: float,
    density_mol_m3: float,
) -> None:
    model = _model(topology, component)
    capsule, table = _table(model)
    assert capsule is not None
    assert table.abi_version == ABI_VERSION
    assert table.table_size >= ACTIVE_KIJ_PREFIX_SIZE
    assert table.result_size == ctypes.sizeof(PhaseBlockResult)
    assert table.parameterized_result_size == ctypes.sizeof(ParameterizedPhaseBlockResult)
    assert table.model_context
    assert table.evaluate_pure_phase
    assert table.evaluate_pure_phase_parameters

    amount_mol = 1.0
    volume_m3 = amount_mol / density_mol_m3
    block = _evaluate(table, temperature_k, amount_mol, volume_m3)
    assert block.status == STATUS_OK
    assert block.error == b""
    assert block.parameter_fingerprint.decode() == model.parameter_fingerprint

    state = model.state(
        T=temperature_k * unit_registry.kelvin,
        rho=density_mol_m3 * unit_registry.mole / unit_registry.meter**3,
        x=(1.0,),
    )
    expected_helmholtz = amount_mol * (math.log(density_mol_m3) - 1.0 + state._residual_helmholtz)
    expected_chemical_potential = (
        math.log(density_mol_m3) + state.residual_chemical_potential_over_rt[0]
    )
    assert block.helmholtz_over_rt_reference_amount == pytest.approx(
        expected_helmholtz, rel=2.0e-12, abs=2.0e-12
    )
    assert block.pressure_pa == pytest.approx(
        state.pressure.to("pascal").magnitude, rel=2.0e-11, abs=1.0e-5
    )
    assert block.chemical_potential_over_rt == pytest.approx(
        expected_chemical_potential, rel=2.0e-11, abs=2.0e-11
    )


def test_native_sdk_derivative_tensors_match_directional_oracles() -> None:
    model = _model("neutral", "methane")
    capsule, table = _table(model)
    del model
    gc.collect()
    assert capsule is not None
    temperature_k = 190.0
    point = (1.2, 1.2e-3)
    direction = (0.05, -3.0e-5)

    def shifted(scale: float) -> PhaseBlockResult:
        result = _evaluate(
            table,
            temperature_k,
            point[0] + scale * direction[0],
            point[1] + scale * direction[1],
        )
        assert result.status == STATUS_OK
        return result

    center = shifted(0.0)
    step = 2.0e-4
    lower = shifted(-step)
    upper = shifted(step)

    gradient_contraction = sum(center.gradient[index] * direction[index] for index in range(2))
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert gradient_contraction == pytest.approx(value_difference, rel=2.0e-8)

    for row in range(2):
        hessian_vector = sum(
            center.hessian[row * 2 + column] * direction[column] for column in range(2)
        )
        gradient_difference = (upper.gradient[row] - lower.gradient[row]) / (2.0 * step)
        assert hessian_vector == pytest.approx(gradient_difference, rel=2.0e-7)

    def hessian_contraction(result: PhaseBlockResult) -> float:
        return sum(
            direction[row] * result.hessian[row * 2 + column] * direction[column]
            for row in range(2)
            for column in range(2)
        )

    third_contraction = sum(
        direction[i] * direction[j] * direction[k] * center.third[(i * 2 + j) * 2 + k]
        for i in range(2)
        for j in range(2)
        for k in range(2)
    )
    hessian_difference = (hessian_contraction(upper) - hessian_contraction(lower)) / (2.0 * step)
    assert third_contraction == pytest.approx(hessian_difference, rel=2.0e-6)


def test_native_sdk_describes_the_exact_neutral_binary_pair_derivative_contracts() -> None:
    model = _mixture_model()
    capsule, table = _capability_table(model)
    assert capsule is not None
    assert table.capability_count == 2
    assert table.capabilities

    descriptor = table.capabilities[0]
    assert descriptor.struct_size == ctypes.sizeof(NativeCapabilityDescriptor)
    assert (
        descriptor.schema_version,
        descriptor.capability,
        descriptor.parameter_family,
        descriptor.parameter_identity,
        descriptor.observation_contract,
        descriptor.model_domain,
        descriptor.tensor_layout,
        descriptor.derivative_order,
        descriptor.maturity,
        descriptor.authority_effect,
    ) == (1, 1, 1, 1, 1, 1, 1, 2, 1, 0)
    assert descriptor.unsupported_status == STATUS_UNSUPPORTED_MODEL
    assert descriptor.domain_status == STATUS_DOMAIN_ERROR
    assert (
        descriptor.state_coordinate_count,
        descriptor.active_parameter_count,
        descriptor.coordinate_count,
        descriptor.component_count,
    ) == (3, 1, 4, 2)
    assert tuple(
        descriptor.component_ids[index].decode() for index in range(descriptor.component_count)
    ) == ("methane", "ethane")
    assert descriptor.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert descriptor.topology_fingerprint.decode() == (
        "sha256:172e3c9dce280763cf785db5d6fc94deeb30a23359346b2222a836077930351a"
    )
    assert descriptor.topology_fingerprint != descriptor.parameter_fingerprint
    assert descriptor.helmholtz_basis_id.decode() == (
        "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3"
    )
    assert math.isnan(descriptor.temperature_min_k)
    assert math.isnan(descriptor.temperature_max_k)
    assert model._model.temperature_min_k is None
    assert model._model.temperature_max_k is None

    coordinates = tuple(
        descriptor.coordinates[index] for index in range(descriptor.coordinate_count)
    )
    assert all(
        coordinate.struct_size == ctypes.sizeof(NativeCapabilityCoordinate)
        for coordinate in coordinates
    )
    assert tuple(coordinate.kind for coordinate in coordinates) == (1, 1, 2, 3)
    assert tuple(coordinate.unit.decode() for coordinate in coordinates) == (
        "mol",
        "mol",
        "m3",
        "dimensionless",
    )
    assert tuple(coordinate.component_index for coordinate in coordinates) == (
        0,
        1,
        -1,
        -1,
    )
    assert tuple(
        (coordinate.pair_component_index_a, coordinate.pair_component_index_b)
        for coordinate in coordinates
    ) == ((-1, -1), (-1, -1), (-1, -1), (0, 1))

    _, reversed_table = _capability_table(_mixture_model(("ethane", "methane")))
    assert reversed_table.capability_count == 2
    assert reversed_table.capabilities[0].topology_fingerprint != descriptor.topology_fingerprint

    lij = table.capabilities[1]
    assert (
        lij.schema_version,
        lij.capability,
        lij.parameter_family,
        lij.parameter_identity,
        lij.observation_contract,
        lij.model_domain,
        lij.tensor_layout,
        lij.derivative_order,
        lij.maturity,
        lij.authority_effect,
    ) == (1, 2, 2, 1, 1, 1, 1, 2, 1, 0)
    assert lij.component_count == 2
    assert tuple(lij.component_ids[index].decode() for index in range(lij.component_count)) == (
        "methane",
        "ethane",
    )
    assert lij.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert lij.topology_fingerprint.decode() == descriptor.topology_fingerprint.decode()
    lij_coordinates = tuple(lij.coordinates[index] for index in range(lij.coordinate_count))
    assert tuple(coordinate.kind for coordinate in lij_coordinates) == (
        1,
        1,
        2,
        4,
    )
    assert tuple(coordinate.unit.decode() for coordinate in lij_coordinates) == (
        "mol",
        "mol",
        "m3",
        "dimensionless",
    )


def test_native_sdk_describes_exact_scalar_pure_parameter_contracts() -> None:
    model = _model("neutral", "methane")
    _, table = _capability_table(model)
    assert table.capability_count == 3
    expected = (
        (3, 3, 5, "dimensionless"),
        (4, 4, 6, "angstrom"),
        (5, 5, 7, "kelvin"),
    )
    for descriptor, (
        capability_code,
        family_code,
        coordinate_kind,
        unit,
    ) in zip(table.capabilities[:3], expected, strict=True):
        assert (
            descriptor.schema_version,
            descriptor.capability,
            descriptor.parameter_family,
            descriptor.parameter_identity,
            descriptor.observation_contract,
            descriptor.model_domain,
            descriptor.tensor_layout,
            descriptor.derivative_order,
            descriptor.maturity,
            descriptor.authority_effect,
        ) == (1, capability_code, family_code, 2, 1, 2, 1, 2, 1, 0)
        assert (
            descriptor.state_coordinate_count,
            descriptor.active_parameter_count,
            descriptor.coordinate_count,
            descriptor.component_count,
        ) == (2, 1, 3, 1)
        coordinates = tuple(
            descriptor.coordinates[index] for index in range(descriptor.coordinate_count)
        )
        assert tuple(item.kind for item in coordinates) == (
            1,
            2,
            coordinate_kind,
        )
        assert tuple(item.unit.decode() for item in coordinates) == (
            "mol",
            "m3",
            unit,
        )
        assert coordinates[-1].component_index == 0
        assert descriptor.parameter_fingerprint.decode() == model.parameter_fingerprint


def test_native_sdk_describes_exact_associating_pure_parameter_contracts() -> None:
    model = _model("associating", "ethanol")
    _, table = _capability_table(model)
    assert table.capability_count == 3
    expected = (
        (15, 9, 11, "kelvin"),
        (16, 10, 12, "dimensionless"),
    )
    for descriptor, (
        capability_code,
        family_code,
        coordinate_kind,
        unit,
    ) in zip(table.capabilities[:2], expected, strict=True):
        assert (
            descriptor.schema_version,
            descriptor.capability,
            descriptor.parameter_family,
            descriptor.parameter_identity,
            descriptor.observation_contract,
            descriptor.model_domain,
            descriptor.tensor_layout,
            descriptor.derivative_order,
            descriptor.maturity,
            descriptor.authority_effect,
        ) == (1, capability_code, family_code, 3, 1, 7, 1, 2, 1, 0)
        coordinates = tuple(
            descriptor.coordinates[index] for index in range(descriptor.coordinate_count)
        )
        assert tuple(item.kind for item in coordinates) == (
            1,
            2,
            coordinate_kind,
        )
        assert tuple(item.unit.decode() for item in coordinates) == (
            "mol",
            "m3",
            unit,
        )
        assert coordinates[-1].component_index == -1
        assert descriptor.topology_fingerprint

    joint = table.capabilities[2]
    assert (
        joint.schema_version,
        joint.capability,
        joint.parameter_family,
        joint.parameter_identity,
        joint.observation_contract,
        joint.model_domain,
        joint.tensor_layout,
        joint.derivative_order,
        joint.maturity,
        joint.authority_effect,
    ) == (1, 21, 13, 3, 1, 7, 1, 2, 1, 0)
    assert (
        joint.state_coordinate_count,
        joint.active_parameter_count,
        joint.coordinate_count,
        joint.component_count,
    ) == (2, 5, 7, 1)
    joint_coordinates = tuple(joint.coordinates[index] for index in range(joint.coordinate_count))
    assert tuple(item.kind for item in joint_coordinates) == (
        1,
        2,
        5,
        6,
        7,
        11,
        12,
    )
    assert tuple(item.unit.decode() for item in joint_coordinates) == (
        "mol",
        "m3",
        "dimensionless",
        "angstrom",
        "kelvin",
        "kelvin",
        "dimensionless",
    )
    assert joint.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert joint.topology_fingerprint


def test_native_sdk_evaluates_ion_fraction_suppression_exact_derivative() -> None:
    model = _electrolyte_model(("water", "sodium-cation", "bromide-anion"))
    _, table = _ion_fraction_suppression_table(model)
    assert table.ion_fraction_suppression_result_size == ctypes.sizeof(IonFractionSuppressionResult)
    callback = EvaluateIonFractionSuppression(table.evaluate_ion_fraction_suppression)

    def evaluate(
        coefficient: float,
        *,
        fingerprint: bytes | None = None,
        temperature: float = 298.15,
        pressure: float = 100000.0,
        mole_fraction: float = 0.1,
    ) -> IonFractionSuppressionResult:
        result = IonFractionSuppressionResult(
            struct_size=ctypes.sizeof(IonFractionSuppressionResult)
        )
        status = callback(
            table.model_context,
            fingerprint or model.parameter_fingerprint.encode(),
            temperature,
            pressure,
            mole_fraction,
            coefficient,
            ctypes.byref(result),
        )
        assert status == result.status
        return result

    coefficient = 7.01
    mole_fraction = 0.1
    salt_free = 78.09
    denominator = 1.0 + coefficient * mole_fraction
    result = evaluate(coefficient)
    assert result.status == STATUS_OK
    assert result.salt_free_relative_permittivity == salt_free
    assert result.bulk_relative_permittivity == pytest.approx(salt_free / denominator)
    assert result.derivative == pytest.approx(-salt_free * mole_fraction / denominator**2)
    assert result.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert not result.error

    step = 1.0e-5
    finite_difference = (
        evaluate(coefficient + step).bulk_relative_permittivity
        - evaluate(coefficient - step).bulk_relative_permittivity
    ) / (2.0 * step)
    assert result.derivative == pytest.approx(finite_difference, rel=2.0e-10)

    for failure, expected in (
        (evaluate(coefficient, fingerprint=b"sha256:wrong"), STATUS_UNSUPPORTED_MODEL),
        (evaluate(coefficient, temperature=300.0), STATUS_DOMAIN_ERROR),
        (evaluate(coefficient, pressure=101325.0), STATUS_DOMAIN_ERROR),
        (evaluate(coefficient, mole_fraction=-0.1), STATUS_DOMAIN_ERROR),
        (evaluate(-1.0), STATUS_DOMAIN_ERROR),
    ):
        assert failure.status == expected
        assert math.isnan(failure.bulk_relative_permittivity)
        assert math.isnan(failure.derivative)
        assert failure.parameter_fingerprint == b""
        assert failure.error

    methanol = _electrolyte_model(("methanol", "sodium-cation", "bromide-anion"))
    _, methanol_table = _ion_fraction_suppression_table(methanol)
    methanol_result = IonFractionSuppressionResult(
        struct_size=ctypes.sizeof(IonFractionSuppressionResult)
    )
    methanol_status = EvaluateIonFractionSuppression(
        methanol_table.evaluate_ion_fraction_suppression
    )(
        methanol_table.model_context,
        methanol.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        mole_fraction,
        coefficient,
        ctypes.byref(methanol_result),
    )
    assert methanol_status == methanol_result.status == STATUS_OK
    assert methanol_result.salt_free_relative_permittivity == 78.09
    assert (
        methanol_result.bulk_relative_permittivity / methanol_result.salt_free_relative_permittivity
    ) == pytest.approx(1.0 / denominator)


def test_native_sdk_mixture_phase_block_matches_public_thermodynamics_and_derivatives() -> None:
    model = _mixture_model()
    capsule, table = _table(model)
    assert capsule is not None
    assert table.component_count == 2
    assert table.mixture_result_size == ctypes.sizeof(MixturePhaseBlockResult)
    assert table.evaluate_mixture_phase

    temperature_k = 220.0
    point = (0.4, 0.6, 1.25e-4)
    direction = (0.02, -0.01, -2.0e-6)

    def shifted(
        scale: float,
    ) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
        return _evaluate_mixture(
            table,
            temperature_k,
            (point[0] + scale * direction[0], point[1] + scale * direction[1]),
            point[2] + scale * direction[2],
        )

    center, gradient, hessian = shifted(0.0)
    assert center.status == STATUS_OK
    assert center.error == b""
    assert center.coordinate_count == 3
    assert center.gradient_capacity == 3
    assert center.hessian_capacity == 9
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint

    total = point[0] + point[1]
    density = total / point[2]
    composition = (point[0] / total, point[1] / total)
    state = model.state(
        T=temperature_k * unit_registry.kelvin,
        rho=density * unit_registry.mole / unit_registry.meter**3,
        x=composition,
    )
    expected = sum(amount * (math.log(amount / point[2]) - 1.0) for amount in point[:2])
    expected += total * state._residual_helmholtz
    assert center.helmholtz_over_rt_reference_amount == pytest.approx(
        expected, rel=2.0e-12, abs=2.0e-12
    )
    assert center.pressure_pa == pytest.approx(
        state.pressure.to("pascal").magnitude, rel=2.0e-11, abs=1.0e-5
    )

    step = 2.0e-4
    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    gradient_contraction = sum(value * delta for value, delta in zip(gradient, direction))
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert gradient_contraction == pytest.approx(value_difference, rel=3.0e-8)
    for row in range(3):
        hessian_direction = sum(
            hessian[row * 3 + column] * direction[column] for column in range(3)
        )
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=3.0e-7)


def test_native_sdk_mixture_phase_block_is_component_permutation_invariant() -> None:
    forward_capsule, forward_table = _table(_mixture_model())
    reverse_capsule, reverse_table = _table(_mixture_model(("ethane", "methane")))
    assert forward_capsule is not None
    assert reverse_capsule is not None
    forward, forward_gradient, forward_hessian = _evaluate_mixture(
        forward_table, 220.0, (0.4, 0.6), 1.25e-4
    )
    reverse, reverse_gradient, reverse_hessian = _evaluate_mixture(
        reverse_table, 220.0, (0.6, 0.4), 1.25e-4
    )
    permutation = (1, 0, 2)
    assert reverse.helmholtz_over_rt_reference_amount == pytest.approx(
        forward.helmholtz_over_rt_reference_amount, rel=2.0e-13, abs=2.0e-13
    )
    assert reverse.pressure_pa == pytest.approx(forward.pressure_pa, rel=2.0e-13)
    for row in range(3):
        assert reverse_gradient[row] == pytest.approx(
            forward_gradient[permutation[row]], rel=2.0e-13, abs=2.0e-13
        )
        for column in range(3):
            assert reverse_hessian[row * 3 + column] == pytest.approx(
                forward_hessian[permutation[row] * 3 + permutation[column]],
                rel=2.0e-13,
                abs=2.0e-13,
            )


@pytest.mark.parametrize(
    ("temperature_k", "point", "direction"),
    (
        (220.0, (0.4, 0.6, 1.25e-4, 0.0), (0.02, -0.01, -2.0e-6, 0.01)),
        (280.0, (0.7, 0.3, 2.5e-4, -0.05), (-0.015, 0.02, 3.0e-6, -0.008)),
    ),
)
def test_native_sdk_active_kij_matches_fixed_and_directional_oracles(
    temperature_k: float,
    point: tuple[float, float, float, float],
    direction: tuple[float, float, float, float],
) -> None:
    model = _mixture_model()
    capsule, table = _table(model)
    assert capsule is not None
    assert table.evaluate_mixture_phase_kij

    def shifted(
        scale: float,
    ) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
        values = tuple(value + scale * delta for value, delta in zip(point, direction))
        return _evaluate_mixture_kij(
            table,
            temperature_k,
            (values[0], values[1]),
            values[2],
            values[3],
        )

    center, gradient, hessian = shifted(0.0)
    assert center.status == STATUS_OK
    assert center.error == b""
    assert center.coordinate_count == 4
    assert center.gradient_capacity == 4
    assert center.hessian_capacity == 16
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert center.pressure_pa == pytest.approx(
        -8.31446261815324 * temperature_k * gradient[2], rel=2.0e-13
    )

    if point[3] == 0.0:
        fixed, fixed_gradient, fixed_hessian = _evaluate_mixture(
            table, temperature_k, point[:2], point[2]
        )
        assert center.helmholtz_over_rt_reference_amount == pytest.approx(
            fixed.helmholtz_over_rt_reference_amount, rel=2.0e-13, abs=2.0e-13
        )
        assert center.pressure_pa == pytest.approx(fixed.pressure_pa, rel=2.0e-13)
        for row in range(3):
            assert gradient[row] == pytest.approx(fixed_gradient[row], rel=2.0e-13, abs=2.0e-13)
            for column in range(3):
                assert hessian[row * 4 + column] == pytest.approx(
                    fixed_hessian[row * 3 + column], rel=2.0e-13, abs=2.0e-13
                )

    step = 2.0e-4
    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    gradient_direction = sum(value * delta for value, delta in zip(gradient, direction))
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert gradient_direction == pytest.approx(value_difference, rel=3.0e-8)
    for row in range(4):
        hessian_direction = sum(
            hessian[row * 4 + column] * direction[column] for column in range(4)
        )
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=3.0e-7)


@pytest.mark.parametrize(
    ("temperature_k", "point", "direction"),
    (
        (220.0, (0.4, 0.6, 1.25e-4, 0.0), (0.02, -0.01, -2.0e-6, 0.01)),
        (280.0, (0.7, 0.3, 2.5e-4, -0.02), (-0.015, 0.02, 3.0e-6, -0.008)),
    ),
)
def test_native_sdk_active_lij_matches_fixed_and_directional_oracles(
    temperature_k: float,
    point: tuple[float, float, float, float],
    direction: tuple[float, float, float, float],
) -> None:
    model = _mixture_model()
    capsule, table = _lij_table(model)
    assert capsule is not None

    def shifted(
        scale: float,
    ) -> tuple[MixturePhaseBlockResult, tuple[float, ...], tuple[float, ...]]:
        values = tuple(value + scale * delta for value, delta in zip(point, direction, strict=True))
        return _evaluate_mixture_lij(
            table,
            temperature_k,
            (values[0], values[1]),
            values[2],
            values[3],
        )

    center, gradient, hessian = shifted(0.0)
    assert center.status == STATUS_OK
    assert center.error == b""
    assert center.coordinate_count == 4
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert center.pressure_pa == pytest.approx(
        -8.31446261815324 * temperature_k * gradient[2], rel=2.0e-13
    )
    if point[3] == 0.0:
        fixed, fixed_gradient, fixed_hessian = _evaluate_mixture(
            table, temperature_k, point[:2], point[2]
        )
        assert center.helmholtz_over_rt_reference_amount == pytest.approx(
            fixed.helmholtz_over_rt_reference_amount,
            rel=2.0e-13,
            abs=2.0e-13,
        )
        assert center.pressure_pa == pytest.approx(fixed.pressure_pa, rel=2.0e-13)
        for row in range(3):
            assert gradient[row] == pytest.approx(fixed_gradient[row], rel=2.0e-13, abs=2.0e-13)
            for column in range(3):
                assert hessian[row * 4 + column] == pytest.approx(
                    fixed_hessian[row * 3 + column],
                    rel=2.0e-13,
                    abs=2.0e-13,
                )

    step = 2.0e-4
    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    gradient_direction = sum(value * delta for value, delta in zip(gradient, direction))
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert gradient_direction == pytest.approx(value_difference, rel=3.0e-8)
    for row in range(4):
        hessian_direction = sum(
            hessian[row * 4 + column] * direction[column] for column in range(4)
        )
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=3.0e-7)


def test_native_sdk_electrolyte_phase_matches_public_state_and_directional_derivatives() -> None:
    model = _electrolyte_model()
    capsule, table = _electrolyte_table(model)
    assert capsule is not None
    assert table.table_size >= ctypes.sizeof(ElectrolyteNativeSdkTable) == 112
    assert table.component_count == 3
    assert tuple(table.component_ids[index].decode() for index in range(3)) == model.component_ids
    assert tuple(table.component_charges[index] for index in range(3)) == (0, 1, -1)
    assert table.evaluate_electrolyte_phase

    temperature_k = 298.15
    point = (0.98, 0.01, 0.01, 1.0e-3)
    center, gradient, hessian = _evaluate_electrolyte(table, temperature_k, point[:3], point[3])
    assert center.status == STATUS_OK
    assert center.parameter_fingerprint.decode() == model.parameter_fingerprint
    total_amount = math.fsum(point[:3])
    density = total_amount / point[3]
    fractions = tuple(amount / total_amount for amount in point[:3])
    public = model.state(
        T=temperature_k * unit_registry.kelvin,
        rho=density * unit_registry.mole / unit_registry.meter**3,
        x=fractions,
    )
    expected_value = (
        math.fsum(amount * (math.log(amount / point[3]) - 1.0) for amount in point[:3])
        + total_amount * public._residual_helmholtz
    )
    assert center.helmholtz_over_rt_reference_amount == pytest.approx(
        expected_value, rel=3.0e-12, abs=3.0e-12
    )
    assert center.pressure_pa == pytest.approx(
        public.pressure.to("pascal").magnitude, rel=3.0e-11, abs=2.0e-5
    )
    assert center.pressure_pa == pytest.approx(
        -8.31446261815324 * temperature_k * gradient[-1], rel=3.0e-13
    )

    direction = (0.02, -0.003, -0.003, 2.0e-5)
    step = 2.0e-4

    def shifted(scale: float):
        state = tuple(value + scale * delta for value, delta in zip(point, direction, strict=True))
        return _evaluate_electrolyte(table, temperature_k, state[:3], state[3])

    lower, lower_gradient, _ = shifted(-step)
    upper, upper_gradient, _ = shifted(step)
    value_difference = (
        upper.helmholtz_over_rt_reference_amount - lower.helmholtz_over_rt_reference_amount
    ) / (2.0 * step)
    assert math.fsum(
        left * right for left, right in zip(gradient, direction, strict=True)
    ) == pytest.approx(value_difference, rel=2.0e-7, abs=2.0e-8)
    for row in range(4):
        hessian_direction = math.fsum(
            hessian[row * 4 + column] * direction[column] for column in range(4)
        )
        gradient_difference = (upper_gradient[row] - lower_gradient[row]) / (2.0 * step)
        assert hessian_direction == pytest.approx(gradient_difference, rel=2.0e-6, abs=2.0e-6)

    permuted = _electrolyte_model(("water", "chloride-anion", "sodium-cation"))
    _, permuted_table = _electrolyte_table(permuted)
    permuted_result, permuted_gradient, permuted_hessian = _evaluate_electrolyte(
        permuted_table, temperature_k, (point[0], point[2], point[1]), point[3]
    )
    permutation = (0, 2, 1, 3)
    assert permuted_result.helmholtz_over_rt_reference_amount == pytest.approx(
        center.helmholtz_over_rt_reference_amount, rel=3.0e-13
    )
    assert permuted_result.pressure_pa == pytest.approx(center.pressure_pa, rel=3.0e-13)
    assert permuted_gradient == pytest.approx(tuple(gradient[index] for index in permutation))
    assert permuted_hessian == pytest.approx(
        tuple(hessian[row * 4 + column] for row in permutation for column in permutation)
    )


def test_native_sdk_packing_fraction_matches_independent_derivatives_and_permutation() -> None:
    canonical_ids = (
        "water",
        "sodium-cation",
        "chloride-anion",
    )
    temperature_k = 298.15
    amounts = (0.98, 0.01, 0.01)
    volume_m3 = 1.0e-3
    model = _electrolyte_model(canonical_ids)

    def evaluate(
        selected_model: Mixture,
        selected_amounts: tuple[float, ...],
    ) -> tuple[object, PackingFractionNativeSdkTable, float, tuple[float, ...], tuple[float, ...]]:
        capsule, table = _packing_fraction_table(selected_model)
        coordinates = len(selected_amounts) + 1
        amount_values = (ctypes.c_double * len(selected_amounts))(*selected_amounts)
        packing_fraction = ctypes.c_double()
        gradient = (ctypes.c_double * coordinates)()
        hessian = (ctypes.c_double * (coordinates * coordinates))()
        status = EvaluatePackingFraction(table.evaluate_packing_fraction)(
            table.model_context,
            selected_model.parameter_fingerprint.encode(),
            temperature_k,
            amount_values,
            len(selected_amounts),
            volume_m3,
            ctypes.byref(packing_fraction),
            gradient,
            coordinates,
            hessian,
            coordinates * coordinates,
        )
        assert status == STATUS_OK
        return capsule, table, packing_fraction.value, tuple(gradient), tuple(hessian)

    capsule, table, value, gradient, hessian = evaluate(model, amounts)
    assert capsule is not None
    assert table.table_size >= ctypes.sizeof(PackingFractionNativeSdkTable) == 128
    assert PackingFractionNativeSdkTable.evaluate_packing_fraction.offset == 120

    water_sigma = (
        2.7927
        + 10.11 * math.exp(-0.01775 * temperature_k)
        - 1.417 * math.exp(-0.01146 * temperature_k)
    )
    diameters = tuple(
        0.88 * sigma if charge else sigma * (1.0 - 0.12 * math.exp(-3.0 * energy / temperature_k))
        for sigma, energy, charge in zip(
            (water_sigma, 2.8232, 2.7560),
            (353.95, 230.0, 170.0),
            (0, 1, -1),
            strict=True,
        )
    )
    weights = tuple(
        math.pi / 6.0 * 6.02214076e23 / 1.0e30 * segments * diameter**3
        for segments, diameter in zip((1.2047, 1.0, 1.0), diameters, strict=True)
    )
    expected_value = (
        math.fsum(amount * weight for amount, weight in zip(amounts, weights, strict=True))
        / volume_m3
    )
    expected_gradient = tuple(weight / volume_m3 for weight in weights) + (
        -expected_value / volume_m3,
    )
    coordinates = len(weights) + 1
    expected_hessian = [0.0] * (coordinates * coordinates)
    volume_index = coordinates - 1
    for component, weight in enumerate(weights):
        expected_hessian[component * coordinates + volume_index] = -weight / volume_m3**2
        expected_hessian[volume_index * coordinates + component] = -weight / volume_m3**2
    expected_hessian[volume_index * coordinates + volume_index] = (
        2.0 * expected_value / volume_m3**2
    )
    assert value == pytest.approx(expected_value, rel=2.0e-15)
    assert gradient == pytest.approx(expected_gradient, rel=2.0e-15, abs=1.0e-15)
    assert hessian == pytest.approx(expected_hessian, rel=3.0e-15, abs=1.0e-14)

    direction = (0.01, -0.002, -0.004, 2.0e-6)
    directional_gradient = math.fsum(
        left * right for left, right in zip(gradient, direction, strict=True)
    )
    oracle_directional_gradient = math.fsum(
        left * right for left, right in zip(expected_gradient, direction, strict=True)
    )
    assert directional_gradient == pytest.approx(oracle_directional_gradient, rel=2.0e-15)
    for row in range(coordinates):
        actual = math.fsum(
            hessian[row * coordinates + column] * direction[column] for column in range(coordinates)
        )
        expected = math.fsum(
            expected_hessian[row * coordinates + column] * direction[column]
            for column in range(coordinates)
        )
        assert actual == pytest.approx(expected, rel=3.0e-15, abs=1.0e-14)

    permutation = (0, 2, 1, 3)
    permuted_model = _electrolyte_model(tuple(canonical_ids[index] for index in permutation[:-1]))
    permuted_capsule, _, permuted_value, permuted_gradient, permuted_hessian = evaluate(
        permuted_model, tuple(amounts[index] for index in permutation[:-1])
    )
    assert permuted_capsule is not None
    assert permuted_value == pytest.approx(value, rel=2.0e-15)
    assert permuted_gradient == pytest.approx(tuple(gradient[index] for index in permutation))
    assert permuted_hessian == pytest.approx(
        tuple(hessian[row * coordinates + column] for row in permutation for column in permutation)
    )


@pytest.mark.parametrize(
    ("active_ion", "counter_ion", "trial_diameter"),
    (
        ("sodium-cation", "chloride-anion", 3.2),
        ("bromide-anion", "sodium-cation", 4.8),
    ),
)
def test_native_sdk_active_born_derivative_matches_centered_difference(
    active_ion: str,
    counter_ion: str,
    trial_diameter: float,
) -> None:
    _, table = _active_born_table(_electrolyte_model(("water", active_ion, counter_ion)))
    center = _evaluate_active_born(table, trial_diameter)

    estimates = []
    for step in (2.0e-4, 1.0e-4):
        upper = _evaluate_active_born(table, trial_diameter + step)
        lower = _evaluate_active_born(table, trial_diameter - step)
        estimates.append(
            (upper.solvation_gibbs_j_per_mol - lower.solvation_gibbs_j_per_mol) / (2.0 * step)
        )
    tolerance = max(
        1.0e-3,
        20.0 * abs(estimates[0] - estimates[1]),
        2.0e-8 * abs(center.derivative_j_per_mol_per_angstrom),
    )
    assert abs(center.derivative_j_per_mol_per_angstrom - estimates[1]) <= tolerance


def test_native_sdk_aqueous_miac_active_parameters_match_public_value_and_exact_derivatives() -> (
    None
):
    model = _electrolyte_model(("water", "sodium-cation", "bromide-anion"))
    _, table = _aqueous_miac_batch_table(model)
    pair_values = model._model.k_ij
    k_ij = (pair_values[1], pair_values[2], pair_values[5])
    molalities = (1.0,)
    kij_center = _evaluate_aqueous_miac_kij_batch(
        table, model.parameter_fingerprint, molalities, k_ij
    )
    solvation_center = _evaluate_aqueous_miac_solvation_factor(
        table, model.parameter_fingerprint, 1.0, 1.5
    )
    reference = model.reference(
        T=298.15 * unit_registry.kelvin,
        P=1.0 * unit_registry.bar,
        solvent_x=(1.0, 0.0, 0.0),
        phase="liquid",
    )
    public = reference.state(molality=1.0 * unit_registry.mole / unit_registry.kilogram)
    expected = math.log(public.activity.mean_molality)

    assert table.table_size >= AQUEOUS_MIAC_KIJ_PREFIX_SIZE == 184
    assert table.aqueous_miac_kij_result_size == ctypes.sizeof(AqueousMiacKijResult) == 296
    assert table.aqueous_miac_solvation_factor_result_size == 280
    assert kij_center[0].status == solvation_center.status == STATUS_OK
    # The derivative callbacks project their values from the certified density
    # midpoint to the requested pressure; this bound covers only that correction.
    assert kij_center[0].log_mean_ionic_activity_coefficient_molality == pytest.approx(
        expected, abs=8.0e-11
    )
    assert solvation_center.log_mean_ionic_activity_coefficient_molality == pytest.approx(
        expected, abs=8.0e-11
    )
    assert solvation_center.reference_molality_mol_per_kg == pytest.approx(
        reference.reference_molality.to("mole / kilogram").magnitude,
        rel=0.0,
        abs=0.0,
    )
    assert solvation_center.reference_convergence_error == pytest.approx(
        reference.convergence_error, rel=0.0, abs=0.0
    )
    assert math.isfinite(solvation_center.reference_derivative_convergence_error)
    assert kij_center[0].parameter_fingerprint.decode() == model.parameter_fingerprint
    assert solvation_center.parameter_fingerprint.decode() == model.parameter_fingerprint

    kij_direction = (0.7, -0.4, 0.2)
    kij_step = 5.0e-5
    kij_high = _evaluate_aqueous_miac_kij_batch(
        table,
        model.parameter_fingerprint,
        molalities,
        tuple(value + kij_step * delta for value, delta in zip(k_ij, kij_direction, strict=True)),
    )
    kij_low = _evaluate_aqueous_miac_kij_batch(
        table,
        model.parameter_fingerprint,
        molalities,
        tuple(value - kij_step * delta for value, delta in zip(k_ij, kij_direction, strict=True)),
    )
    kij_estimate = (
        kij_high[0].log_mean_ionic_activity_coefficient_molality
        - kij_low[0].log_mean_ionic_activity_coefficient_molality
    ) / (2.0 * kij_step)
    kij_directional = math.fsum(
        kij_center[0].derivative[index] * kij_direction[index] for index in range(3)
    )
    assert abs(kij_directional - kij_estimate) <= max(1.0e-8, 2.0e-8 * abs(kij_directional))

    for solvation_factor in (1.35, 1.65):
        center = _evaluate_aqueous_miac_solvation_factor_batch(
            table,
            model.parameter_fingerprint,
            molalities,
            solvation_factor,
        )[0]
        estimates = []
        for step in (1.0e-4, 5.0e-5):
            high = _evaluate_aqueous_miac_solvation_factor_batch(
                table,
                model.parameter_fingerprint,
                molalities,
                solvation_factor + step,
            )[0]
            low = _evaluate_aqueous_miac_solvation_factor_batch(
                table,
                model.parameter_fingerprint,
                molalities,
                solvation_factor - step,
            )[0]
            estimates.append(
                (
                    high.log_mean_ionic_activity_coefficient_molality
                    - low.log_mean_ionic_activity_coefficient_molality
                )
                / (2.0 * step)
            )
        tolerance = max(
            1.0e-8,
            20.0 * abs(estimates[0] - estimates[1]),
            2.0e-8 * abs(center.derivative),
        )
        assert abs(center.derivative - estimates[1]) <= tolerance


def test_native_sdk_neutral_reference_matches_public_owner_and_is_permutation_safe() -> None:
    model = _electrolyte_model()
    capsule, table = _neutral_reference_table(model)
    assert capsule is not None
    assert table.table_size >= NEUTRAL_REFERENCE_PREFIX_SIZE == 208
    assert table.neutral_reference_basis_row_count == 2
    assert table.neutral_reference_result_size == ctypes.sizeof(NeutralReferenceResult) == 440

    result, basis, contractions, reference_composition = _evaluate_neutral_reference(
        table, model.parameter_fingerprint
    )
    public_reference = model.reference(
        T=298.15 * unit_registry.kelvin,
        P=1.0 * unit_registry.bar,
        solvent_x=(1.0, 0.0, 0.0),
        phase="liquid",
    )
    pure_water = model.state(
        T=298.15 * unit_registry.kelvin,
        P=1.0 * unit_registry.bar,
        x=(1.0, 0.0, 0.0),
        phase="liquid",
    )
    assert pure_water._log_fugacity_coefficient is not None

    assert result.status == STATUS_OK
    assert tuple(table.component_ids[index].decode() for index in range(3)) == model.component_ids
    assert tuple(table.component_charges[index] for index in range(3)) == (0, 1, -1)
    assert basis == (1.0, 0.0, 0.0, 0.0, 1.0, 1.0)
    assert reference_composition == (1.0, 0.0, 0.0)
    assert contractions[0] == pytest.approx(pure_water._log_fugacity_coefficient[0], abs=2e-12)
    assert contractions[1] == pytest.approx(
        math.fsum(public_reference.infinite_dilution_log_fugacity_coefficient),
        abs=2e-12,
    )
    assert result.temperature_k == 298.15
    assert result.pressure_pa == 100000.0
    assert model._model.pressure_min_pa == result.pressure_pa
    assert model._model.pressure_max_pa == result.pressure_pa
    assert result.solvent_molar_mass_kg_per_mol == model._model.molar_mass_kg_per_mol[0]
    assert result.reference_amount_mol == 1.0
    assert result.reference_number_density_mol_per_m3 == 1.0
    assert result.reference_molality_mol_per_kg == 1.0e-12
    assert result.reference_convergence_error == public_reference.convergence_error
    assert result.derivative_availability == 0
    assert result.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert result.helmholtz_basis_id.decode() == (
        "A_over_RT_reference_amount:n_ref=1mol:rho_ref=1mol_per_m3"
    )
    assert result.error == b""

    permuted = _electrolyte_model(("chloride-anion", "water", "sodium-cation"))
    _, permuted_table = _neutral_reference_table(permuted)
    permuted_result, permuted_basis, permuted_contractions, permuted_composition = (
        _evaluate_neutral_reference(permuted_table, permuted.parameter_fingerprint)
    )
    assert (
        tuple(permuted_table.component_ids[index].decode() for index in range(3))
        == permuted.component_ids
    )
    assert tuple(permuted_table.component_charges[index] for index in range(3)) == (
        -1,
        0,
        1,
    )
    assert permuted_basis == (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    assert permuted_composition == (0.0, 1.0, 0.0)
    assert permuted_contractions == pytest.approx(contractions, abs=3e-12)
    assert permuted_result.reference_convergence_error == pytest.approx(
        result.reference_convergence_error, abs=5e-12
    )
    assert permuted_result.parameter_fingerprint.decode() == permuted.parameter_fingerprint

    callback = EvaluateNeutralReference(permuted_table.evaluate_neutral_reference)
    basis_buffer = (ctypes.c_double * 6)(*([123.0] * 6))
    contraction_buffer = (ctypes.c_double * 2)(123.0, 123.0)
    composition_buffer = (ctypes.c_double * 3)(123.0, 123.0, 123.0)

    def call(
        *,
        fingerprint: bytes = permuted.parameter_fingerprint.encode(),
        temperature: float = 298.15,
        pressure: float = 100000.0,
        basis_capacity: int = 6,
    ) -> NeutralReferenceResult:
        failure = NeutralReferenceResult(
            struct_size=ctypes.sizeof(NeutralReferenceResult),
            component_count=3,
            neutral_basis_row_count=2,
            neutral_basis_capacity=basis_capacity,
            contraction_capacity=2,
            reference_composition_capacity=3,
            neutral_basis=basis_buffer,
            log_fugacity_contractions=contraction_buffer,
            reference_composition=composition_buffer,
        )
        status = callback(
            permuted_table.model_context,
            fingerprint,
            temperature,
            pressure,
            ctypes.byref(failure),
        )
        assert status == failure.status
        return failure

    for failure, expected_status in (
        (call(fingerprint=b"sha256:wrong"), STATUS_UNSUPPORTED_MODEL),
        (call(temperature=299.0), STATUS_DOMAIN_ERROR),
        (call(pressure=101325.0), STATUS_DOMAIN_ERROR),
        (call(basis_capacity=5), STATUS_INVALID_ARGUMENT),
    ):
        assert failure.status == expected_status
        assert failure.parameter_fingerprint == b""
        assert failure.helmholtz_basis_id == b""
        assert failure.error
        assert tuple(basis_buffer) == (123.0,) * 6
        assert tuple(contraction_buffer) == (123.0,) * 2
        assert tuple(composition_buffer) == (123.0,) * 3

    undersized = NeutralReferenceResult(struct_size=8, status=123)
    before = bytes(undersized)
    status = callback(
        permuted_table.model_context,
        permuted.parameter_fingerprint.encode(),
        298.15,
        100000.0,
        ctypes.byref(undersized),
    )
    assert status == STATUS_ABI_MISMATCH
    assert bytes(undersized) == before

    unsupported = _mixture_model()
    _, unsupported_table = _neutral_reference_table(unsupported)
    assert unsupported_table.neutral_reference_basis_row_count == 0
    unsupported_result = NeutralReferenceResult(struct_size=ctypes.sizeof(NeutralReferenceResult))
    unsupported_status = EvaluateNeutralReference(unsupported_table.evaluate_neutral_reference)(
        unsupported_table.model_context,
        unsupported.parameter_fingerprint.encode(),
        220.0,
        100000.0,
        ctypes.byref(unsupported_result),
    )
    assert unsupported_status == unsupported_result.status == STATUS_UNSUPPORTED_MODEL
    assert unsupported_result.error


def test_native_sdk_neutral_reference__uses_charge_role_not_solvent_name() -> None:
    model = _electrolyte_model(("generic-solvent", "generic-cation", "generic-anion"))
    capsule, table = _neutral_reference_table(model)

    result, _, _, composition = _evaluate_neutral_reference(
        table,
        model.parameter_fingerprint,
    )

    assert capsule is not None
    assert result.status == STATUS_OK
    assert composition == (1.0, 0.0, 0.0)


def test_native_sdk_internal_reference_and_density_failure_contracts(tmp_path: Path) -> None:
    compiler = shutil.which("c++")
    assert compiler is not None
    native_root = Path(__file__).parents[2] / "src" / "epcsaft" / "_native"
    source = r"""
#include "density.hpp"
#include "electrolyte_reference.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

struct State { double pressure; double pressure_density_derivative; };

int main() {
    bool rejected_nonfinite_reference = false;
    int call = 0;
    try {
        epcsaft::detail::electrolyte_reference_sequence(
            std::vector<double>{1.0, 0.0, 0.0},
            std::vector<double>{0.018, 0.007, 0.035},
            std::array<std::size_t, 2>{1, 2},
            std::array<int, 2>{1, 1},
            [&call](const std::vector<double>&) {
                ++call;
                std::vector<double> derivative(6, 1.0);
                if (call == 3) {
                    derivative[4] = std::numeric_limits<double>::quiet_NaN();
                }
                return epcsaft::detail::ReferenceObservation{{0.0, 0.0}, derivative};
            }
        );
    } catch (const std::runtime_error&) {
        rejected_nonfinite_reference = true;
    }
    if (!rejected_nonfinite_reference) return 1;

    bool domain_classified = false;
    try {
        epcsaft::detail::solve_density(
            1.0, 1.0, 10.0, "liquid",
            [](double density) { return State{density, 1.0}; }
        );
    } catch (const epcsaft::detail::DensityDomainError&) {
        domain_classified = true;
    }
    if (!domain_classified) return 2;

    bool evaluation_classified = false;
    try {
        epcsaft::detail::refine_density_root(
            1.0, 0.5, 1.5, -0.5,
            [](double) {
                return State{std::numeric_limits<double>::quiet_NaN(), 1.0};
            }
        );
    } catch (const epcsaft::detail::DensityDomainError&) {
        return 3;
    } catch (const std::runtime_error&) {
        evaluation_classified = true;
    }
    if (!evaluation_classified) return 4;

    const double adjacent_lower = 59456.796534681715;
    const double adjacent_upper = std::nextafter(
        adjacent_lower,
        std::numeric_limits<double>::infinity()
    );
    const auto adjacent_root = epcsaft::detail::refine_density_root(
        1.0, adjacent_lower, adjacent_upper, -2.0e-6,
        [adjacent_lower](double density) {
            return State{
                density <= adjacent_lower ? 1.0 - 2.0e-6 : 1.0 + 2.0e-6,
                1.0,
            };
        }
    );
    if (!(adjacent_root.molar_density >= adjacent_lower
          && adjacent_root.molar_density <= adjacent_upper)) return 5;
    return 0;
}
"""
    output = tmp_path / "owner-contract"
    compiled = subprocess.run(
        (compiler, "-std=c++17", f"-I{native_root}", "-x", "c++", "-", "-o", str(output)),
        input=source,
        text=True,
        capture_output=True,
        check=False,
    )
    assert compiled.returncode == 0, compiled.stderr
    executed = subprocess.run((str(output),), capture_output=True, text=True, check=False)
    assert executed.returncode == 0, executed.stderr


def test_native_sdk_inverse_packing_geometry_matches_forward_geometry_and_derivatives() -> None:
    model = _electrolyte_model()
    capsule, table = _inverse_packing_table(model)
    capability_capsule, capability_table = _capability_table(model)

    assert capsule is not None and capability_capsule is not None
    assert table.inverse_packing_geometry_result_size == ctypes.sizeof(InversePackingGeometryResult)
    assert table.evaluate_inverse_packing_geometry
    descriptors = tuple(
        capability_table.capabilities[index]
        for index in range(capability_table.capability_count)
        if capability_table.capabilities[index].capability in (19, 20)
    )
    topology_fingerprints = {descriptor.topology_fingerprint.decode() for descriptor in descriptors}
    assert len(topology_fingerprints) == 1
    topology = topology_fingerprints.pop()
    amounts = (0.98, 0.01, 0.01)
    target_eta = 0.25

    ordinary, ordinary_gradient, ordinary_hessian = _evaluate_inverse_packing_geometry(
        table,
        model.parameter_fingerprint,
        topology,
        amounts,
        math.log(target_eta),
        (),
    )
    assert ordinary.status == STATUS_OK
    assert ordinary.coordinate_count == 4
    assert ordinary.active_parameter_count == 0
    assert ordinary.volume_m3 > 0.0
    assert all(math.isfinite(value) for value in ordinary_gradient)
    assert all(math.isfinite(value) for value in ordinary_hessian)
    geometry_only, _, _ = _evaluate_inverse_packing_geometry(
        table,
        model.parameter_fingerprint,
        topology,
        (0.1, 0.45, 0.45),
        math.log(target_eta),
        (),
    )
    assert geometry_only.status == STATUS_OK
    _, packing_table = _packing_fraction_table(model)
    packing_value = ctypes.c_double()
    packing_gradient = (ctypes.c_double * 4)()
    packing_hessian = (ctypes.c_double * 16)()
    packing_status = EvaluatePackingFraction(packing_table.evaluate_packing_fraction)(
        packing_table.model_context,
        model.parameter_fingerprint.encode(),
        298.15,
        (ctypes.c_double * 3)(*amounts),
        3,
        ordinary.volume_m3,
        ctypes.byref(packing_value),
        packing_gradient,
        4,
        packing_hessian,
        16,
    )
    assert packing_status == STATUS_OK
    assert abs(packing_value.value - target_eta) <= math.ulp(target_eta)

    sigma = float(model._model.segment_diameter[1].value.to("angstrom").magnitude)
    epsilon = float(model._model.dispersion_energy_over_k[1])
    requests = (
        _active_request(4, 1, sigma + 0.01, b"angstrom"),
        _active_request(5, 1, epsilon + 1.0, b"kelvin"),
    )
    result, gradient, hessian = _evaluate_inverse_packing_geometry(
        table,
        model.parameter_fingerprint,
        topology,
        amounts,
        math.log(target_eta),
        requests,
    )
    assert result.status == STATUS_OK
    assert result.coordinate_count == 6
    assert result.active_parameter_count == 2
    assert result.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert result.topology_fingerprint.decode() == topology
    assert result.volume_m3 > 0.0
    assert gradient[3] == pytest.approx(-result.volume_m3, rel=2.0e-14)
    assert hessian[3 * 6 + 3] == pytest.approx(result.volume_m3, rel=2.0e-14)
    assert sum(amounts[index] * gradient[index] for index in range(3)) == pytest.approx(
        result.volume_m3, rel=3.0e-14
    )
    reacting_capsule, reacting_table = _reacting_sensitivity_table(model)
    assert reacting_capsule is not None
    forward = _evaluate_reacting_phase_parameters(
        reacting_table,
        model.parameter_fingerprint,
        topology,
        amounts,
        result.volume_m3,
        requests,
    )
    assert abs(forward[0].packing_fraction - target_eta) <= math.ulp(target_eta)

    direction = (0.2, -0.004, 0.003, 0.1, 0.02, 2.0)
    step = 1.0e-5
    shifted = []
    for sign in (-1.0, 1.0):
        shifted_amounts = tuple(
            amounts[index] + sign * step * direction[index] for index in range(3)
        )
        shifted_requests = (
            _active_request(
                4,
                1,
                requests[0].value + sign * step * direction[4],
                b"angstrom",
            ),
            _active_request(
                5,
                1,
                requests[1].value + sign * step * direction[5],
                b"kelvin",
            ),
        )
        shifted.append(
            _evaluate_inverse_packing_geometry(
                table,
                model.parameter_fingerprint,
                topology,
                shifted_amounts,
                math.log(target_eta) + sign * step * direction[3],
                shifted_requests,
            )
        )
    lower, upper = shifted
    numerical_value_direction = (upper[0].volume_m3 - lower[0].volume_m3) / (2.0 * step)
    assert numerical_value_direction == pytest.approx(
        sum(gradient[index] * direction[index] for index in range(6)),
        rel=2.0e-8,
        abs=2.0e-12,
    )
    numerical_gradient_direction = tuple(
        (upper[1][row] - lower[1][row]) / (2.0 * step) for row in range(6)
    )
    exact_gradient_direction = tuple(
        sum(hessian[row * 6 + column] * direction[column] for column in range(6))
        for row in range(6)
    )
    assert numerical_gradient_direction == pytest.approx(
        exact_gradient_direction, rel=3.0e-7, abs=3.0e-9
    )

    reversed_result, reversed_gradient, reversed_hessian = _evaluate_inverse_packing_geometry(
        table,
        model.parameter_fingerprint,
        topology,
        amounts,
        math.log(target_eta),
        tuple(reversed(requests)),
    )
    permutation = (0, 1, 2, 3, 5, 4)
    assert reversed_result.volume_m3 == result.volume_m3
    assert reversed_gradient == pytest.approx(tuple(gradient[index] for index in permutation))
    assert reversed_hessian == pytest.approx(
        tuple(hessian[row * 6 + column] for row in permutation for column in permutation)
    )
