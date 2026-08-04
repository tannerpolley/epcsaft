from __future__ import annotations

from ._support import *

_A0 = (
    0.9105631445,
    0.6361281449,
    2.6861347891,
    -26.547362491,
    97.759208784,
    -159.59154087,
    91.297774084,
)
_A1 = (
    -0.3084016918,
    0.1860531159,
    -2.5030047259,
    21.419793629,
    -65.255885330,
    83.318680481,
    -33.746922930,
)
_A2 = (
    -0.0906148351,
    0.4527842806,
    0.5962700728,
    -1.7241829131,
    -4.1302112531,
    13.776631870,
    -8.6728470368,
)
_B0 = (
    0.7240946941,
    2.2382791861,
    -4.0025849485,
    -21.003576815,
    26.855641363,
    206.55133841,
    -355.60235612,
)
_B1 = (
    -0.5755498075,
    0.6995095521,
    3.8925673390,
    -17.215471648,
    192.67226447,
    -161.82646165,
    -165.20769346,
)
_B2 = (
    0.0976883116,
    -0.2557574982,
    -9.1558561530,
    20.642075974,
    -38.804430052,
    93.626774077,
    -29.666905585,
)


def _independent_fixed_topology_phase(
    point: tuple[float, ...],
    *,
    temperature_k: float,
    site_count: int,
    multiplicities: tuple[int, ...],
    pairs: tuple[tuple[int, int], ...],
) -> float:
    """Equation-level oracle independent of the native and public EOS paths."""
    amount, volume, segment_count, sigma, epsilon = point[:5]
    density = amount / volume
    number_density = density * 6.02214076e-7
    diameter = sigma * (1.0 - 0.12 * math.exp(-3.0 * epsilon / temperature_k))
    zeta = tuple(
        math.pi / 6.0 * number_density * segment_count * diameter**order for order in range(4)
    )
    eta = zeta[3]
    one_minus_eta = 1.0 - eta
    hard_sphere = (
        3.0 * zeta[1] * zeta[2] / one_minus_eta
        + zeta[2] ** 3 / (zeta[3] * one_minus_eta**2)
        + (zeta[2] ** 3 / zeta[3] ** 2 - zeta[0]) * math.log(one_minus_eta)
    ) / zeta[0]
    contact = (
        1.0 / one_minus_eta
        + 1.5 * diameter * zeta[2] / one_minus_eta**2
        + 0.5 * diameter**2 * zeta[2] ** 2 / one_minus_eta**3
    )
    hard_chain = segment_count * hard_sphere - (segment_count - 1.0) * math.log(contact)

    chain_1 = (segment_count - 1.0) / segment_count
    chain_2 = chain_1 * (segment_count - 2.0) / segment_count
    a = tuple(x + chain_1 * y + chain_2 * z for x, y, z in zip(_A0, _A1, _A2, strict=True))
    b = tuple(x + chain_1 * y + chain_2 * z for x, y, z in zip(_B0, _B1, _B2, strict=True))
    integral_1 = sum(coefficient * eta**order for order, coefficient in enumerate(a))
    integral_2 = sum(coefficient * eta**order for order, coefficient in enumerate(b))
    c1 = 1.0 / (
        1.0
        + segment_count * (8.0 * eta - 2.0 * eta**2) / one_minus_eta**4
        + (1.0 - segment_count)
        * (20.0 * eta - 27.0 * eta**2 + 12.0 * eta**3 - 2.0 * eta**4)
        / ((one_minus_eta * (2.0 - eta)) ** 2)
    )
    first_moment = segment_count**2 * sigma**3 * epsilon / temperature_k
    second_moment = segment_count**2 * sigma**3 * (epsilon / temperature_k) ** 2
    dispersion = (
        -2.0 * math.pi * number_density * integral_1 * first_moment
        - math.pi * number_density * segment_count * c1 * integral_2 * second_moment
    )

    deltas = [[0.0] * site_count for _ in range(site_count)]
    for pair_index, (left, right) in enumerate(pairs):
        energy, association_volume = point[5 + 2 * pair_index : 7 + 2 * pair_index]
        delta = contact * math.expm1(energy / temperature_k) * sigma**3 * association_volume
        deltas[left][right] = delta
        deltas[right][left] = delta
    fractions = [1.0] * site_count
    for _ in range(2000):
        fixed = [
            1.0
            / (
                1.0
                + number_density
                * sum(
                    multiplicities[other] * fractions[other] * deltas[site][other]
                    for other in range(site_count)
                )
            )
            for site in range(site_count)
        ]
        next_fractions = [0.5 * (old + new) for old, new in zip(fractions, fixed, strict=True)]
        if (
            max(abs(new - old) for old, new in zip(fractions, next_fractions, strict=True))
            < 2.0e-14
        ):
            fractions = next_fractions
            break
        fractions = next_fractions
    else:
        raise AssertionError("independent association oracle did not converge")
    association = sum(
        multiplicity * (math.log(fraction) - 0.5 * fraction + 0.5)
        for multiplicity, fraction in zip(multiplicities, fractions, strict=True)
    )
    return amount * (math.log(density) - 1.0 + hard_chain + dispersion + association)


def _assert_derivatives_match_oracle(
    result: FixedTopologyPhaseResult,
    gradient: tuple[float, ...],
    hessian: tuple[float, ...],
    point: tuple[float, ...],
    steps: tuple[float, ...],
    oracle: Callable[[tuple[float, ...]], float],
) -> None:
    cache: dict[tuple[float, ...], float] = {}

    def value(candidate: tuple[float, ...]) -> float:
        if candidate not in cache:
            cache[candidate] = oracle(candidate)
        return cache[candidate]

    def displaced(*changes: tuple[int, float]) -> tuple[float, ...]:
        candidate = list(point)
        for coordinate, delta in changes:
            candidate[coordinate] += delta
        return tuple(candidate)

    center = value(point)
    assert result.helmholtz_over_rt_reference_amount == pytest.approx(
        center, rel=3.0e-12, abs=3.0e-12
    )
    for coordinate, step in enumerate(steps):
        lower = value(displaced((coordinate, -step)))
        upper = value(displaced((coordinate, step)))
        assert gradient[coordinate] == pytest.approx(
            (upper - lower) / (2.0 * step), rel=3.0e-6, abs=3.0e-8
        )
        for other in range(coordinate + 1):
            other_step = steps[other]
            if other == coordinate:
                expected = (upper - 2.0 * center + lower) / step**2
            else:
                expected = (
                    value(displaced((coordinate, step), (other, other_step)))
                    - value(displaced((coordinate, step), (other, -other_step)))
                    - value(displaced((coordinate, -step), (other, other_step)))
                    + value(displaced((coordinate, -step), (other, -other_step)))
                ) / (4.0 * step * other_step)
            assert hessian[coordinate * len(point) + other] == pytest.approx(
                expected, rel=3.0e-3, abs=3.0e-5
            )
            assert hessian[other * len(point) + coordinate] == pytest.approx(
                expected, rel=3.0e-3, abs=3.0e-5
            )


def _descriptor(
    model: Mixture,
) -> tuple[object, FixedTopologyNativeSdkTable, FixedTopologyAssociationDescriptor]:
    capsule, table = _fixed_topology_table(model)
    assert table.fixed_topology_association
    return capsule, table, table.fixed_topology_association.contents


def test_fixed_topology_associating_mixture_exposes_generic_active_slots() -> None:
    model = Mixture.from_dictionary(association_mixture_case(polar=True))
    capsule, table = _fixed_topology_mixture_table(model)
    assert capsule is not None
    descriptor = table.fixed_topology_association.contents
    assert table.fixed_topology_associating_mixture_result_size == ctypes.sizeof(
        FixedTopologyMixturePhaseResult
    )
    assert table.fixed_topology_association_result_size == 0
    assert not table.evaluate_fixed_topology_association
    assert descriptor.model_domain == 9
    assert descriptor.component_count == 2
    assert descriptor.state_coordinate_count == 3
    assert descriptor.temperature_min_k == 290.0
    assert descriptor.temperature_max_k == 320.0
    assert tuple(
        descriptor.component_ids[index].decode() for index in range(descriptor.component_count)
    ) == ("component-a", "component-b")
    slots = tuple(descriptor.slots[index] for index in range(descriptor.slot_count))
    assert tuple(slot.component_index for slot in slots) == (
        0,
        0,
        0,
        1,
        1,
        1,
        -1,
        -1,
    )
    temperature_k = 305.0
    point = (0.4, 0.6, 1.0 / 5500.0, 2.2, 0.027)

    def evaluate(
        candidate: tuple[float, ...],
    ) -> tuple[
        FixedTopologyMixturePhaseResult,
        tuple[float, ...],
        tuple[float, ...],
    ]:
        requests = _fixed_topology_requests(
            (slots[3], slots[7]),
            candidate[3:],
        )
        return _evaluate_fixed_topology_mixture(
            table,
            model.parameter_fingerprint,
            descriptor.topology_fingerprint.decode(),
            temperature_k,
            candidate[:2],
            candidate[2],
            requests,
        )

    result, gradient, hessian = evaluate(point)
    assert result.status == STATUS_OK
    assert len(gradient) == 5
    assert len(hessian) == 25
    assert result.pressure_pa == pytest.approx(-8.31446261815324 * temperature_k * gradient[2])
    assert result.parameter_fingerprint.decode() == model.parameter_fingerprint
    assert result.topology_fingerprint == descriptor.topology_fingerprint
    assert result.component_order_fingerprint.decode().startswith("sha256:")
    assert result.domain_fingerprint.decode().startswith("sha256:")
    assert result.artifact_fingerprint == descriptor.artifact_fingerprint

    def public_state(candidate: tuple[float, ...]):
        parameters = association_mixture_case(polar=True)
        parameters["parameters"]["m"] = (1.4, candidate[3])
        parameters["parameters"]["association"][0]["association_volume"] = candidate[4]
        public_model = Mixture.from_dictionary(parameters)
        total_amount = candidate[0] + candidate[1]
        return public_model.state(
            T=temperature_k * unit_registry.kelvin,
            rho=(total_amount / candidate[2]) * unit_registry.mole / unit_registry.meter**3,
            x=(candidate[0] / total_amount, candidate[1] / total_amount),
        )

    def public_value(candidate: tuple[float, ...]) -> float:
        total_amount = candidate[0] + candidate[1]
        return sum(
            amount * (math.log(amount / candidate[2]) - 1.0) for amount in candidate[:2]
        ) + total_amount * float(public_state(candidate).ares(basis="reduced"))

    state = public_state(point)
    expected = public_value(point)
    assert result.helmholtz_over_rt_reference_amount == pytest.approx(
        expected, rel=3.0e-12, abs=3.0e-12
    )
    assert result.pressure_pa == pytest.approx(
        state.pressure.to("pascal").magnitude, rel=3.0e-11, abs=1.0e-5
    )

    steps = (2.0e-5, 2.0e-5, 2.0e-9, 2.0e-4, 2.0e-6)
    for column, step in enumerate(steps):
        lower_point = list(point)
        upper_point = list(point)
        lower_point[column] -= step
        upper_point[column] += step
        finite_gradient = (public_value(tuple(upper_point)) - public_value(tuple(lower_point))) / (
            2.0 * step
        )
        assert gradient[column] == pytest.approx(finite_gradient, rel=5.0e-7, abs=5.0e-7)
        for row in range(5):
            if row == column:
                finite_hessian = (
                    public_value(tuple(upper_point))
                    - 2.0 * expected
                    + public_value(tuple(lower_point))
                ) / (step * step)
            else:
                row_step = steps[row]
                shifted_values = []
                for row_direction, column_direction in (
                    (1.0, 1.0),
                    (1.0, -1.0),
                    (-1.0, 1.0),
                    (-1.0, -1.0),
                ):
                    shifted = list(point)
                    shifted[row] += row_direction * row_step
                    shifted[column] += column_direction * step
                    shifted_values.append(public_value(tuple(shifted)))
                finite_hessian = (
                    shifted_values[0] - shifted_values[1] - shifted_values[2] + shifted_values[3]
                ) / (4.0 * row_step * step)
            assert hessian[row * 5 + column] == pytest.approx(
                finite_hessian, rel=2.0e-4, abs=2.0e-3
            )


def test_fixed_topology_associating_ternary_is_permutation_safe() -> None:
    base_order = ("component-a", "component-b", "component-c")
    permuted_order = ("component-c", "component-a", "component-b")
    amounts_by_id = {
        "component-a": 0.2,
        "component-b": 0.5,
        "component-c": 0.3,
    }

    def evaluate(order: tuple[str, ...]):
        model = Mixture.from_dictionary(association_mixture_case(order))
        _, table = _fixed_topology_mixture_table(model)
        descriptor = table.fixed_topology_association.contents
        assert (
            tuple(
                descriptor.component_ids[index].decode()
                for index in range(descriptor.component_count)
            )
            == order
        )
        slots = tuple(descriptor.slots[index] for index in range(descriptor.slot_count))
        selected = tuple(
            slot
            for slot in slots
            if slot.parameter_id.decode() == "component:component-b/segment_count"
            or slot.parameter_family == 10
        )
        assert len(selected) == 3
        result, gradient, hessian = _evaluate_fixed_topology_mixture(
            table,
            model.parameter_fingerprint,
            descriptor.topology_fingerprint.decode(),
            305.0,
            tuple(amounts_by_id[item] for item in order),
            1.0 / 5200.0,
            _fixed_topology_requests(selected, (2.2, 0.027, 0.027)),
        )
        assert result.status == STATUS_OK
        return result, gradient, hessian

    base, base_gradient, base_hessian = evaluate(base_order)
    permuted, permuted_gradient, permuted_hessian = evaluate(permuted_order)
    assert permuted.helmholtz_over_rt_reference_amount == pytest.approx(
        base.helmholtz_over_rt_reference_amount, rel=3.0e-13, abs=3.0e-13
    )
    base_to_permuted = (1, 2, 0, 3, 4, 5, 6)
    for base_index, permuted_index in enumerate(base_to_permuted):
        assert permuted_gradient[permuted_index] == pytest.approx(
            base_gradient[base_index], rel=2.0e-11, abs=2.0e-9
        )
        for base_column, permuted_column in enumerate(base_to_permuted):
            assert permuted_hessian[permuted_index * 7 + permuted_column] == pytest.approx(
                base_hessian[base_index * 7 + base_column],
                rel=3.0e-10,
                abs=3.0e-7,
            )


def test_fixed_topology_associating_mixture_fails_closed() -> None:
    model = Mixture.from_dictionary(association_mixture_case())
    capsule, table = _fixed_topology_mixture_table(model)
    assert capsule is not None
    descriptor = table.fixed_topology_association.contents
    slots = tuple(descriptor.slots[index] for index in range(descriptor.slot_count))
    requests = _fixed_topology_requests((slots[3], slots[7]), (2.2, 0.027))

    for boundary in (290.0, 320.0):
        result, gradient, hessian = _evaluate_fixed_topology_mixture(
            table,
            model.parameter_fingerprint,
            descriptor.topology_fingerprint.decode(),
            boundary,
            (0.4, 0.6),
            1.0 / 5500.0,
            requests,
        )
        assert result.status == STATUS_OK
        assert gradient and hessian

    sentinel = 12345.0
    rejected, gradient, hessian = _evaluate_fixed_topology_mixture(
        table,
        model.parameter_fingerprint,
        descriptor.topology_fingerprint.decode(),
        289.99,
        (0.4, 0.6),
        1.0 / 5500.0,
        requests,
        buffer_value=sentinel,
    )
    assert rejected.status == STATUS_DOMAIN_ERROR
    assert {*gradient, *hessian} == {sentinel}

    wrong_order, gradient, hessian = _evaluate_fixed_topology_mixture(
        table,
        model.parameter_fingerprint,
        descriptor.topology_fingerprint.decode(),
        305.0,
        (0.4, 0.6),
        1.0 / 5500.0,
        requests,
        component_ids=("component-b", "component-a"),
        buffer_value=sentinel,
    )
    assert wrong_order.status == STATUS_INVALID_ARGUMENT
    assert {*gradient, *hessian} == {sentinel}

    wrong_fingerprint, gradient, hessian = _evaluate_fixed_topology_mixture(
        table,
        f"sha256:{'0' * 64}",
        descriptor.topology_fingerprint.decode(),
        305.0,
        (0.4, 0.6),
        1.0 / 5500.0,
        requests,
        buffer_value=sentinel,
    )
    assert wrong_fingerprint.status == STATUS_UNSUPPORTED_MODEL
    assert {*gradient, *hessian} == {sentinel}
    wrong_artifact, gradient, hessian = _evaluate_fixed_topology_mixture(
        table,
        model.parameter_fingerprint,
        descriptor.topology_fingerprint.decode(),
        305.0,
        (0.4, 0.6),
        1.0 / 5500.0,
        requests,
        buffer_value=sentinel,
        artifact_fingerprint=f"sha256:{'f' * 64}".encode(),
    )
    assert wrong_artifact.status == STATUS_UNSUPPORTED_MODEL
    assert {*gradient, *hessian} == {sentinel}

    ordered, ordered_gradient, ordered_hessian = _evaluate_fixed_topology_mixture(
        table,
        model.parameter_fingerprint,
        descriptor.topology_fingerprint.decode(),
        305.0,
        (0.4, 0.6),
        1.0 / 5500.0,
        requests,
    )
    reversed_slots, reversed_gradient, reversed_hessian = _evaluate_fixed_topology_mixture(
        table,
        model.parameter_fingerprint,
        descriptor.topology_fingerprint.decode(),
        305.0,
        (0.4, 0.6),
        1.0 / 5500.0,
        _fixed_topology_requests((slots[7], slots[3]), (0.027, 2.2)),
    )
    assert reversed_slots.status == STATUS_OK
    assert reversed_slots.helmholtz_over_rt_reference_amount == pytest.approx(
        ordered.helmholtz_over_rt_reference_amount, rel=2.0e-13, abs=2.0e-13
    )
    order = (0, 1, 2, 4, 3)
    for row, reversed_row in enumerate(order):
        assert reversed_gradient[reversed_row] == pytest.approx(ordered_gradient[row])
        for column, reversed_column in enumerate(order):
            assert reversed_hessian[reversed_row * 5 + reversed_column] == pytest.approx(
                ordered_hessian[row * 5 + column]
            )
    outside_slot_domain, gradient, hessian = _evaluate_fixed_topology_mixture(
        table,
        model.parameter_fingerprint,
        descriptor.topology_fingerprint.decode(),
        305.0,
        (0.4, 0.6),
        1.0 / 5500.0,
        _fixed_topology_requests((slots[3],), (0.0,)),
        buffer_value=sentinel,
    )
    assert outside_slot_domain.status == STATUS_DOMAIN_ERROR
    assert {*gradient, *hessian} == {sentinel}
    nonfinite_amount, gradient, hessian = _evaluate_fixed_topology_mixture(
        table,
        model.parameter_fingerprint,
        descriptor.topology_fingerprint.decode(),
        305.0,
        (math.nan, 0.6),
        1.0 / 5500.0,
        requests,
        buffer_value=sentinel,
    )
    assert nonfinite_amount.status == STATUS_INVALID_ARGUMENT
    assert {*gradient, *hessian} == {sentinel}

    wrong_amount_count, gradient, hessian = _evaluate_fixed_topology_mixture(
        table,
        model.parameter_fingerprint,
        descriptor.topology_fingerprint.decode(),
        305.0,
        (0.4,),
        1.0 / 5500.0,
        requests,
        buffer_value=sentinel,
    )
    assert wrong_amount_count.status == STATUS_INVALID_ARGUMENT
    assert {*gradient, *hessian} == {sentinel}
    wrong_capacity, gradient, hessian = _evaluate_fixed_topology_mixture(
        table,
        model.parameter_fingerprint,
        descriptor.topology_fingerprint.decode(),
        305.0,
        (0.4, 0.6),
        1.0 / 5500.0,
        requests,
        buffer_value=sentinel,
        gradient_capacity=4,
    )
    assert wrong_capacity.status == STATUS_ABI_MISMATCH
    assert {*gradient, *hessian} == {sentinel}

    _, pure_table = _fixed_topology_mixture_table(Mixture.from_dictionary(ASSOCIATING_PURE))
    assert pure_table.fixed_topology_associating_mixture_result_size == 0
    assert not pure_table.evaluate_fixed_topology_associating_mixture
    _, ordinary_table = _fixed_topology_mixture_table(_mixture_model())
    assert ordinary_table.fixed_topology_associating_mixture_result_size == 0
    assert not ordinary_table.evaluate_fixed_topology_associating_mixture


def test_fixed_topology_descriptor_replaces_legacy_association_contracts() -> None:
    model = Mixture.from_dictionary(ASSOCIATING_PURE)
    capsule, table, descriptor = _descriptor(model)
    assert capsule is not None
    assert table.reserved_abi_size_0 == 0
    assert not table.reserved_abi_slot_1
    assert not table.reserved_abi_slot_2
    assert table.fixed_topology_association_result_size == ctypes.sizeof(FixedTopologyPhaseResult)
    assert descriptor.schema_version == 1
    assert descriptor.derivative_order == 2
    assert descriptor.site_count == 2
    assert tuple(descriptor.site_ids[index].decode() for index in range(2)) == (
        "a",
        "b",
    )
    assert descriptor.slot_count == 5
    assert descriptor.artifact_fingerprint.decode().startswith("sha256:")
    slots = tuple(descriptor.slots[index] for index in range(descriptor.slot_count))
    assert tuple(slot.parameter_family for slot in slots) == (3, 4, 5, 9, 10)
    assert tuple(slot.component_index for slot in slots) == (0, 0, 0, -1, -1)
    assert all(slot.lower_bound_exclusive == 0.0 for slot in slots)
    assert all(math.isinf(slot.upper_bound_exclusive) for slot in slots)
    assert (
        slots[3].parameter_id.decode().endswith("a--associating-fluid/b/association_energy_over_k")
    )
    _, capability_table = _capability_table(model)
    assert capability_table.capability_count == 0


def test_fixed_topology_2b_matches_independent_value_gradient_hessian_oracle() -> None:
    values = (1.0, 3.0, 200.0, 2400.0, 0.03)
    model = Mixture.from_dictionary(ASSOCIATING_PURE)
    _, table, descriptor = _descriptor(model)
    slots = tuple(descriptor.slots[index] for index in range(descriptor.slot_count))
    point = (1.0, 1.0 / 6000.0, *values)
    result, gradient, hessian = _evaluate_fixed_topology(
        table,
        model.parameter_fingerprint,
        descriptor.topology_fingerprint.decode(),
        300.0,
        point[0],
        point[1],
        _fixed_topology_requests(slots, values),
    )
    _assert_derivatives_match_oracle(
        result,
        gradient,
        hessian,
        point,
        (2.0e-4, 2.0e-8, 2.0e-4, 2.0e-4, 2.0e-2, 2.0e-1, 2.0e-6),
        lambda candidate: _independent_fixed_topology_phase(
            candidate,
            temperature_k=300.0,
            site_count=2,
            multiplicities=(1, 1),
            pairs=((0, 1),),
        ),
    )


def test_fixed_topology_non_2b_value_gradient_hessian_subsets_and_identity() -> None:
    values = (1.4, 3.2, 220.0, 2400.0, 0.030, 1800.0, 0.025)
    model = Mixture.from_dictionary(three_site_association_case(values))
    _, table, descriptor = _descriptor(model)
    assert tuple(descriptor.site_ids[index].decode() for index in range(3)) == (
        "a",
        "b",
        "c",
    )
    assert tuple(descriptor.site_multiplicities[index] for index in range(3)) == (
        1,
        1,
        2,
    )
    slots = tuple(descriptor.slots[index] for index in range(descriptor.slot_count))
    assert descriptor.slot_count == 7
    assert tuple((slot.site_index_a, slot.site_index_b) for slot in slots[3::2]) == (
        (0, 1),
        (0, 2),
    )
    requests = _fixed_topology_requests(slots, values)
    temperature_k = 300.0
    point = (1.0, 1.0 / 6000.0, *values)
    topology_fingerprint = descriptor.topology_fingerprint.decode()
    result, gradient, hessian = _evaluate_fixed_topology(
        table,
        model.parameter_fingerprint,
        topology_fingerprint,
        temperature_k,
        point[0],
        point[1],
        requests,
    )
    assert result.status == STATUS_OK
    assert result.artifact_fingerprint == descriptor.artifact_fingerprint
    assert result.pressure_pa == pytest.approx(-8.31446261815324 * temperature_k * gradient[1])
    assert result.chemical_potential_over_rt == pytest.approx(gradient[0])

    steps = (
        2.0e-4,
        2.0e-8,
        2.0e-4,
        2.0e-4,
        2.0e-2,
        2.0e-1,
        2.0e-6,
        2.0e-1,
        2.0e-6,
    )

    _assert_derivatives_match_oracle(
        result,
        gradient,
        hessian,
        point,
        steps,
        lambda candidate: _independent_fixed_topology_phase(
            candidate,
            temperature_k=temperature_k,
            site_count=3,
            multiplicities=(1, 1, 2),
            pairs=((0, 1), (0, 2)),
        ),
    )

    for indices in ((0,), (3, 4), (0, 4, 5), tuple(range(7))):
        subset = tuple(requests[index] for index in indices)
        subset_result, subset_gradient, subset_hessian = _evaluate_fixed_topology(
            table,
            model.parameter_fingerprint,
            topology_fingerprint,
            temperature_k,
            point[0],
            point[1],
            subset,
        )
        assert subset_result.status == STATUS_OK
        assert len(subset_gradient) == 2 + len(indices)
        assert len(subset_hessian) == (2 + len(indices)) ** 2
        full_coordinates = (0, 1, *(2 + index for index in indices))
        assert subset_result.helmholtz_over_rt_reference_amount == pytest.approx(
            result.helmholtz_over_rt_reference_amount,
            rel=2.0e-13,
            abs=2.0e-13,
        )
        assert subset_result.pressure_pa == pytest.approx(
            result.pressure_pa, rel=2.0e-13, abs=2.0e-8
        )
        assert subset_result.chemical_potential_over_rt == pytest.approx(
            result.chemical_potential_over_rt, rel=2.0e-13, abs=2.0e-13
        )
        assert subset_gradient == pytest.approx(
            tuple(gradient[index] for index in full_coordinates),
            rel=2.0e-11,
            abs=2.0e-10,
        )
        assert subset_hessian == pytest.approx(
            tuple(
                hessian[left * len(point) + right]
                for left in full_coordinates
                for right in full_coordinates
            ),
            rel=2.0e-10,
            abs=2.0e-8,
        )

    shared_energy_requests = tuple(
        FixedTopologyParameterRequest(
            struct_size=requests[index].struct_size,
            parameter_family=requests[index].parameter_family,
            component_index=requests[index].component_index,
            site_index_a=requests[index].site_index_a,
            site_index_b=requests[index].site_index_b,
            unit=requests[index].unit,
            value=2100.0,
        )
        for index in (3, 5)
    )
    shared_result, shared_gradient, shared_hessian = _evaluate_fixed_topology(
        table,
        model.parameter_fingerprint,
        topology_fingerprint,
        temperature_k,
        point[0],
        point[1],
        shared_energy_requests,
    )
    shared_point = (*point[:5], 2100.0, point[6], 2100.0, point[8])
    _assert_derivatives_match_oracle(
        shared_result,
        shared_gradient,
        shared_hessian,
        (shared_point[0], shared_point[1], shared_point[5], shared_point[7]),
        (2.0e-4, 2.0e-8, 2.0e-1, 2.0e-1),
        lambda candidate: _independent_fixed_topology_phase(
            (
                candidate[0],
                candidate[1],
                shared_point[2],
                shared_point[3],
                shared_point[4],
                candidate[2],
                shared_point[6],
                candidate[3],
                shared_point[8],
            ),
            temperature_k=temperature_k,
            site_count=3,
            multiplicities=(1, 1, 2),
            pairs=((0, 1), (0, 2)),
        ),
    )

    failures = (
        _evaluate_fixed_topology(
            table,
            f"sha256:{'0' * 64}",
            topology_fingerprint,
            temperature_k,
            point[0],
            point[1],
            requests,
        ),
        _evaluate_fixed_topology(
            table,
            model.parameter_fingerprint,
            f"sha256:{'1' * 64}",
            temperature_k,
            point[0],
            point[1],
            requests,
        ),
        _evaluate_fixed_topology(
            table,
            model.parameter_fingerprint,
            topology_fingerprint,
            temperature_k,
            point[0],
            point[1],
            requests,
            expected_artifact_fingerprint=f"sha256:{'2' * 64}",
        ),
        _evaluate_fixed_topology(
            table,
            model.parameter_fingerprint,
            topology_fingerprint,
            temperature_k,
            point[0],
            point[1],
            (requests[5], requests[4]),
        ),
        _evaluate_fixed_topology(
            table,
            model.parameter_fingerprint,
            topology_fingerprint,
            temperature_k,
            point[0],
            point[1],
            (requests[3], requests[3]),
        ),
    )
    for failed, failed_gradient, failed_hessian in failures:
        assert failed.status in (STATUS_INVALID_ARGUMENT, STATUS_UNSUPPORTED_MODEL)
        assert not any(failed_gradient)
        assert not any(failed_hessian)

    failed_inner_request = FixedTopologyParameterRequest(
        struct_size=requests[3].struct_size,
        parameter_family=requests[3].parameter_family,
        component_index=requests[3].component_index,
        site_index_a=requests[3].site_index_a,
        site_index_b=requests[3].site_index_b,
        unit=requests[3].unit,
        value=1.0e6,
    )
    failed_inner, failed_gradient, failed_hessian = _evaluate_fixed_topology(
        table,
        model.parameter_fingerprint,
        topology_fingerprint,
        temperature_k,
        point[0],
        point[1],
        (failed_inner_request,),
    )
    assert failed_inner.status == STATUS_EVALUATION_ERROR
    assert not any(failed_gradient)
    assert not any(failed_hessian)

    overflowing_request = FixedTopologyParameterRequest(
        struct_size=requests[2].struct_size,
        parameter_family=requests[2].parameter_family,
        component_index=requests[2].component_index,
        site_index_a=requests[2].site_index_a,
        site_index_b=requests[2].site_index_b,
        unit=requests[2].unit,
        value=1.0e308,
    )
    overflowing, overflowing_gradient, overflowing_hessian = _evaluate_fixed_topology(
        table,
        model.parameter_fingerprint,
        topology_fingerprint,
        temperature_k,
        point[0],
        point[1],
        (overflowing_request,),
        output_sentinel=17.0,
    )
    assert overflowing.status == STATUS_EVALUATION_ERROR
    assert overflowing_gradient == (17.0,) * 3
    assert overflowing_hessian == (17.0,) * 9

    guarded_gradient = (ctypes.c_double * 3)(11.0, 12.0, 13.0)
    guarded_hessian = (ctypes.c_double * 9)(*range(21, 30))
    malformed = FixedTopologyPhaseResult(
        struct_size=ctypes.sizeof(FixedTopologyPhaseResult),
        coordinate_count=3,
        active_parameter_count=1,
        gradient_capacity=2,
        hessian_capacity=9,
        gradient=guarded_gradient,
        hessian=guarded_hessian,
    )
    request = (FixedTopologyParameterRequest * 1)(requests[3])
    status = EvaluateFixedTopologyAssociation(table.evaluate_fixed_topology_association)(
        table.model_context,
        model.parameter_fingerprint.encode(),
        topology_fingerprint.encode(),
        descriptor.artifact_fingerprint,
        temperature_k,
        point[0],
        point[1],
        request,
        1,
        ctypes.byref(malformed),
    )
    assert status == malformed.status == STATUS_ABI_MISMATCH
    assert tuple(guarded_gradient) == (11.0, 12.0, 13.0)
    assert tuple(guarded_hessian) == tuple(range(21, 30))


def test_fixed_topology_canonicalizes_site_order_and_preserves_multiplicity() -> None:
    values = (1.4, 3.2, 220.0, 2400.0, 0.030, 1800.0, 0.025)
    base = Mixture.from_dictionary(three_site_association_case(values))
    _, base_table, base_descriptor = _descriptor(base)
    base_slots = tuple(base_descriptor.slots[index] for index in range(base_descriptor.slot_count))
    base_result, base_gradient, base_hessian = _evaluate_fixed_topology(
        base_table,
        base.parameter_fingerprint,
        base_descriptor.topology_fingerprint.decode(),
        300.0,
        1.0,
        1.0 / 6000.0,
        _fixed_topology_requests(base_slots, values),
    )

    permuted = Mixture.from_dictionary(
        three_site_association_case(values, site_order=("c", "a", "b"))
    )
    _, table, descriptor = _descriptor(permuted)
    assert tuple(descriptor.site_ids[index].decode() for index in range(3)) == (
        "a",
        "b",
        "c",
    )
    assert tuple(descriptor.site_multiplicities[index] for index in range(3)) == (
        1,
        1,
        2,
    )
    slots = tuple(descriptor.slots[index] for index in range(descriptor.slot_count))
    assert tuple((slot.site_index_a, slot.site_index_b) for slot in slots[3::2]) == (
        (0, 1),
        (0, 2),
    )
    result, gradient, hessian = _evaluate_fixed_topology(
        table,
        permuted.parameter_fingerprint,
        descriptor.topology_fingerprint.decode(),
        300.0,
        1.0,
        1.0 / 6000.0,
        _fixed_topology_requests(slots, values),
    )
    assert result.helmholtz_over_rt_reference_amount == pytest.approx(
        base_result.helmholtz_over_rt_reference_amount, rel=2.0e-13, abs=2.0e-13
    )
    assert gradient == pytest.approx(
        base_gradient,
        rel=2.0e-11,
        abs=2.0e-10,
    )
    assert hessian == pytest.approx(
        base_hessian,
        rel=2.0e-10,
        abs=2.0e-8,
    )
