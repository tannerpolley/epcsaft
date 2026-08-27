from __future__ import annotations

import hashlib
import json
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any, cast

import epcsaft

from .. import _core
from .._trials import NonEvaluableTrial, NonEvaluableTrialCode
from ._api import (
    ActiveParameterTrialReceipt,
    AllComponents,
    ChemicalEquilibriumProblem,
    EquilibriumContinuationState,
    EquilibriumPhaseState,
    _continuation_start,
    _contract_phase_output,
    _support_active_parameters,
)


@dataclass(frozen=True, slots=True)
class _ReactiveLiquidDiagnostics:
    provider_evaluation_count: int
    balance_inf_norm: float
    charge_inf_norm: float
    reaction_affinity_inf_norm: float
    pressure_relative_residual: float
    solver_status: str
    provider_domain_status: str
    local_minimum_status: str = "passed"


@dataclass(frozen=True, slots=True)
class _ReactiveLiquidState:
    species_ids: tuple[str, ...]
    temperature_k: float
    pressure_pa: float
    amounts_mol: tuple[float, ...]
    volume_m3: float
    packing_fraction: float
    chemical_potential_over_rt: tuple[float, ...]
    log_fugacities: tuple[float, ...]
    chemical_potential_directions: dict[str, tuple[float, ...]]
    log_fugacity_directions: dict[str, tuple[float, ...]]
    log_fugacity_parameter_partials: dict[str, tuple[float, ...]]
    resolved_policy: tuple[tuple[str, float], ...]
    diagnostics: _ReactiveLiquidDiagnostics

    @property
    def mole_fractions(self) -> tuple[float, ...]:
        total = math.fsum(self.amounts_mol)
        return tuple(amount / total for amount in self.amounts_mol)


def _canonical_fingerprint(schema: str, value: object) -> str:
    return (
        "sha256:"
        + hashlib.sha256(
            json.dumps(
                {"schema": schema, "value": value},
                allow_nan=False,
                separators=(",", ":"),
                sort_keys=True,
            ).encode()
        ).hexdigest()
    )


@dataclass(frozen=True, slots=True)
class _ProviderVapor:
    """The declared neutral subset evaluated by the installed nonideal EOS."""

    reference_id: str

    def __post_init__(self) -> None:
        if not self.reference_id:
            raise ValueError("nonideal vapor model requires a reference identity")


@dataclass(frozen=True, slots=True)
class _IdealVapor:
    """Explicit ideal-vapor comparison model; never a failure fallback."""

    reference_id: str

    def __post_init__(self) -> None:
        if not self.reference_id:
            raise ValueError("ideal-vapor comparison requires a reference identity")


ReactiveVaporModel = _ProviderVapor | _IdealVapor


@dataclass(frozen=True, slots=True)
class _FiniteIncipientKernel:
    """Numerical inputs for one finite-liquid/incipient-vapor solve."""

    liquid_phase_identity: str
    liquid_model: epcsaft.Mixture
    liquid_packing_fraction_interval: tuple[float, float]
    reaction_system: ChemicalEquilibriumProblem
    vapor_phase_identity: str
    vapor_component_ids: tuple[str, ...]
    vapor_model: ReactiveVaporModel
    pressure_interval_pa: tuple[float, float]
    pressure_starts_pa: tuple[float, ...]
    continuation_state: EquilibriumContinuationState | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "vapor_component_ids", tuple(self.vapor_component_ids))
        object.__setattr__(self, "pressure_interval_pa", tuple(self.pressure_interval_pa))
        object.__setattr__(self, "pressure_starts_pa", tuple(self.pressure_starts_pa))
        species = self.reaction_system.species_ids
        charges = self.reaction_system.charges
        charge_by_id = dict(zip(species, charges, strict=True))
        if (
            not self.liquid_phase_identity
            or not self.vapor_phase_identity
            or not isinstance(self.liquid_model, epcsaft.Mixture)
            or len(self.liquid_packing_fraction_interval) != 2
            or not isinstance(self.reaction_system, ChemicalEquilibriumProblem)
            or (
                self.continuation_state is not None
                and len(self.continuation_state.phase(self.liquid_phase_identity).mole_fractions)
                != len(species)
            )
            or not self.vapor_component_ids
            or len(set(self.vapor_component_ids)) != len(self.vapor_component_ids)
            or any(identity not in charge_by_id for identity in self.vapor_component_ids)
            or any(charge_by_id[identity] != 0 for identity in self.vapor_component_ids)
            or len(self.pressure_interval_pa) != 2
            or not 0.0 < self.pressure_interval_pa[0] < self.pressure_interval_pa[1]
            or not self.pressure_starts_pa
            or len(set(self.pressure_starts_pa)) != len(self.pressure_starts_pa)
            or any(
                not math.isfinite(value)
                or not self.pressure_interval_pa[0] <= value <= self.pressure_interval_pa[1]
                for value in self.pressure_starts_pa
            )
        ):
            raise ValueError("finite-incipient kernel inputs are invalid")


@dataclass(frozen=True)
class _Closure:
    pressure_pa: float
    residual: float
    residual_log_pressure_derivative: float
    equilibrium: _ReactiveLiquidState
    trial: ActiveParameterTrialReceipt
    liquid_state_identity: str
    phases: tuple[EquilibriumPhaseState, EquilibriumPhaseState]
    boundary_directional_derivatives: dict[str, float]
    boundary_direct_derivatives: dict[str, float]
    vapor_fraction_directional_derivatives: dict[str, tuple[float, ...]]
    vapor_fraction_direct_derivatives: dict[str, tuple[float, ...]]
    partner_parameter_fingerprint: str
    partner_component_ids: tuple[str, ...]
    chemical_amount_directions: dict[str, tuple[float, ...]]
    chemical_volume_directions: dict[str, float]
    balance_direction_ids: tuple[str, ...]


class _IncipientFailure(RuntimeError):
    def __init__(self, code: NonEvaluableTrialCode, diagnostic: str) -> None:
        super().__init__(diagnostic)
        self.code = code


def _liquid_trial_state(
    problem: _FiniteIncipientKernel,
    temperature: Any,
    descriptor: Any,
    active_parameters: epcsaft.ActiveParameterSet | None,
    pressure_pa: float,
) -> tuple[
    _ReactiveLiquidState,
    ActiveParameterTrialReceipt,
    str,
    dict[str, tuple[float, ...]],
    dict[str, float],
    tuple[str, ...],
]:
    reaction_system = problem.reaction_system
    if any(
        amount <= reaction_system.strict_interior_amount_floor_mol
        for amount in reaction_system.feed_amounts_mol
    ):
        raise _IncipientFailure(
            "active_bound_ambiguity",
            "the declared liquid feed is outside the strict positive derivative chart",
        )
    from ._grepe import (
        EquilibriumOutput,
        FinitePhaseStart,
        Fixed,
        GeneralReactiveEquilibriumProblem,
        ProviderModel,
        ReactivePhase,
        _compile_general_reactive_equilibrium,
        _evaluate_finite_phase,
        _FinitePhaseCompilation,
    )

    reference = problem.continuation_state
    continuation_start = _continuation_start(
        reaction_system,
        None if reference is None else reference.phase(problem.liquid_phase_identity),
    )
    start = (
        None
        if continuation_start is None
        else FinitePhaseStart(
            continuation_start[0],
            continuation_start[1] / math.fsum(continuation_start[0]),
        )
    )
    component_count = len(reaction_system.species_ids)
    outputs = tuple(
        EquilibriumOutput(
            f"incipient-liquid-log-fugacity/{identity}",
            "phase.log_fugacity",
            "dimensionless",
            "provider-helmholtz-coordinate-basis",
            problem.liquid_phase_identity,
            tuple(1.0 if index == species else 0.0 for index in range(component_count)),
        )
        for species, identity in enumerate(reaction_system.species_ids)
    )
    owner_problem = GeneralReactiveEquilibriumProblem(
        f"incipient-liquid/{problem.liquid_phase_identity}",
        Fixed(temperature),
        Fixed(pressure_pa * epcsaft.unit_registry.pascal),
        (
            ReactivePhase(
                problem.liquid_phase_identity,
                "liquid",
                "finite",
                AllComponents(),
                ProviderModel(
                    problem.liquid_packing_fraction_interval,
                    "installed-provider-eos",
                ),
                start,
            ),
        ),
        reaction_system,
        (problem.liquid_phase_identity,),
        outputs,
    )
    compiled = _compile_general_reactive_equilibrium(
        problem.liquid_model,
        owner_problem,
        active_parameters=active_parameters,
    )
    if not isinstance(compiled, _FinitePhaseCompilation):
        raise _IncipientFailure(
            "owner_contract_failure", "incipient liquid did not compile to a finite phase"
        )
    native: list[dict[str, object]] = []
    try:
        result = _evaluate_finite_phase(
            problem.liquid_model,
            owner_problem,
            compiled,
            active_parameters,
            native_sink=native,
        )
    except (ArithmeticError, RuntimeError, ValueError) as error:
        raise _IncipientFailure("numerical_convergence_failure", str(error)) from error
    active_ids = descriptor.active_parameter_identities
    if result.status != "evaluated" or result.trial is None or not native:
        raise _IncipientFailure(
            "numerical_convergence_failure",
            "reactive finite-phase liquid solve was non-evaluable",
        )
    block = native[0]
    amounts = tuple(next(iter(cast(Sequence[Sequence[float]], block["phase_amounts_mol"]))))
    if any(amount <= reaction_system.strict_interior_amount_floor_mol for amount in amounts):
        raise _IncipientFailure(
            "active_bound_ambiguity",
            "the reactive liquid contains a trace-bound species at the coupled state",
        )
    input_order = tuple(str(value) for value in cast(Sequence[object], block["input_order"]))
    width = len(input_order)
    amount_derivatives = tuple(cast(Sequence[float], block["input_amount_derivatives"]))
    volume_derivatives = tuple(cast(Sequence[float], block["input_volume_derivatives"]))
    logf_derivatives = tuple(cast(Sequence[float], block["input_log_fugacity_derivatives"]))
    balance_ids = tuple(
        f"compiled-balance-total/{index}"
        for index, identity in enumerate(input_order)
        if identity.startswith("balance_total[")
    )
    names = (*balance_ids, "pressure_pa", *active_ids)
    columns = {
        **{
            identity: index
            for index, identity in enumerate(input_order)
            if identity in {"pressure_pa", *active_ids}
        },
        **{
            identity: index
            for identity, index in zip(
                balance_ids,
                (i for i, value in enumerate(input_order) if value.startswith("balance_total[")),
                strict=True,
            )
        },
    }
    if any(identity not in columns for identity in names):
        raise _IncipientFailure(
            "parameter_topology_mismatch",
            "reactive finite-phase sensitivity order changed",
        )
    amount_rows = {
        identity: tuple(
            amount_derivatives[species * width + columns[identity]]
            for species in range(component_count)
        )
        for identity in names
    }
    volume_rows = {identity: volume_derivatives[columns[identity]] for identity in names}
    logf_rows = {
        identity: tuple(
            logf_derivatives[species * width + columns[identity]]
            for species in range(component_count)
        )
        for identity in names
    }
    parameter_partials = tuple(cast(Sequence[float], block["output_parameter_partials"]))
    logf_partials = {
        identity: tuple(
            parameter_partials[species * len(active_ids) + parameter]
            for species in range(component_count)
        )
        for parameter, identity in enumerate(active_ids)
    }
    temperature_k = float(temperature.to("kelvin").magnitude)
    equilibrium = _ReactiveLiquidState(
        reaction_system.species_ids,
        temperature_k,
        pressure_pa,
        amounts,
        float(next(iter(cast(Sequence[float], block["phase_volumes_m3"])))),
        float(next(iter(cast(Sequence[float], block["phase_packing_fractions"])))),
        tuple(
            value - math.log(8.31446261815324 * temperature_k)
            for value in next(
                iter(cast(Sequence[Sequence[float]], block["phase_log_fugacities"]))
            )
        ),
        tuple(next(iter(cast(Sequence[Sequence[float]], block["phase_log_fugacities"])))),
        logf_rows,
        logf_rows,
        logf_partials,
        result.resolved_policy,
        _ReactiveLiquidDiagnostics(
            int(cast(int, block["provider_evaluation_count"])),
            float(cast(float, block["balance_inf_norm"])),
            abs(
                math.fsum(
                    charge * amount
                    for charge, amount in zip(reaction_system.charges, amounts, strict=True)
                )
            ),
            float(cast(float, block["reaction_affinity_inf_norm"])),
            float(cast(float, block["pressure_relative_inf_norm"])),
            str(block["solver_status"]),
            str(block["provider_domain_status"]),
        ),
    )
    trial = result.trial
    state_identity = _canonical_fingerprint(
        "epcsaft.incipient-liquid-state",
        {
            "phase_identity": problem.liquid_phase_identity,
            "temperature_k": equilibrium.temperature_k,
            "pressure_pa": equilibrium.pressure_pa,
            "amounts_mol": equilibrium.amounts_mol,
            "volume_m3": equilibrium.volume_m3,
            "trial": trial.fingerprint,
        },
    )
    return equilibrium, trial, state_identity, amount_rows, volume_rows, balance_ids


@dataclass(frozen=True)
class _IncipientVaporPoint:
    mole_fractions: tuple[float, ...]
    volume_m3: float
    pressure_pa: float
    packing_fraction: float
    chemical_potential_over_rt: tuple[float, ...]
    mu_log_ratio_at_fixed_volume: tuple[tuple[float, ...], ...]
    mu_log_volume_derivative: tuple[float, ...]
    mu_log_ratio_jacobian: tuple[tuple[float, ...], ...]
    pressure_log_ratio_derivative: tuple[float, ...]
    pressure_log_volume_derivative: float
    mu_pressure_derivative: tuple[float, ...]
    mu_parameter_derivatives: dict[str, tuple[float, ...]]


def _softmax_reference_chart(log_ratios: tuple[float, ...]) -> tuple[float, ...]:
    logits = (*log_ratios, 0.0)
    maximum = max(logits)
    weights = tuple(math.exp(value - maximum) for value in logits)
    total = math.fsum(weights)
    fractions = tuple(value / total for value in weights)
    if any(not math.isfinite(value) or value <= 0.0 for value in fractions):
        raise _IncipientFailure(
            "active_bound_ambiguity",
            "declared vapor support left its strict positive composition chart",
        )
    return fractions


def _solve_incipient_linear_system(
    matrix: tuple[tuple[float, ...], ...],
    right_hand_side: tuple[float, ...],
) -> tuple[float, ...]:
    """Solve one small dense tangent system with fail-closed pivot checks."""

    dimension = len(right_hand_side)
    if dimension == 0:
        return ()
    work = [list(row) + [value] for row, value in zip(matrix, right_hand_side, strict=True)]
    scale = max(1.0, max(abs(value) for row in matrix for value in row))
    minimum_pivot = math.inf
    maximum_pivot = 0.0
    for column in range(dimension):
        pivot = max(range(column, dimension), key=lambda row: abs(work[row][column]))
        pivot_value = abs(work[pivot][column])
        if not math.isfinite(pivot_value) or pivot_value <= 1.0e-13 * scale:
            raise _IncipientFailure(
                "singular_sensitivity",
                "incipient-vapor composition Jacobian is singular",
            )
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
        minimum_pivot = min(minimum_pivot, pivot_value)
        maximum_pivot = max(maximum_pivot, pivot_value)
        for row in range(column + 1, dimension):
            multiplier = work[row][column] / work[column][column]
            for entry in range(column, dimension + 1):
                work[row][entry] -= multiplier * work[column][entry]
    if maximum_pivot / minimum_pivot > 1.0e10:
        raise _IncipientFailure(
            "ill_conditioned_sensitivity",
            "incipient-vapor composition Jacobian is ill-conditioned",
        )
    solution = [0.0] * dimension
    for row in range(dimension - 1, -1, -1):
        solution[row] = (
            work[row][dimension]
            - math.fsum(
                work[row][column] * solution[column] for column in range(row + 1, dimension)
            )
        ) / work[row][row]
    if any(not math.isfinite(value) for value in solution):
        raise _IncipientFailure(
            "singular_sensitivity",
            "incipient-vapor linear solve returned nonfinite values",
        )
    return tuple(solution)


def _evaluate_incipient_vapor_point(
    *,
    model: epcsaft.Pure | epcsaft.Mixture,
    provider_model: Callable[[], _core.NativeModel],
    active_parameters: epcsaft.ActiveParameterSet | None,
    all_active_identities: tuple[str, ...],
    temperature_k: float,
    log_ratios: tuple[float, ...],
    log_volume: float,
    packing_fraction_interval: tuple[float, float],
) -> _IncipientVaporPoint:
    """Evaluate one explicit composition/volume point and its exact tangent blocks."""

    fractions = _softmax_reference_chart(log_ratios)
    active_ids = () if active_parameters is None else active_parameters.identities
    active_position = {identity: index for index, identity in enumerate(active_ids)}
    active_native = (
        model._native
        if active_parameters is None
        else model._native.with_active_coordinates(model._active_coordinates(active_parameters))
    )
    component_count = len(fractions)
    state_count = component_count + 1
    volume_coordinate = component_count
    try:
        volume = math.exp(log_volume)
        lower_volume, upper_volume = (
            float(value)
            for value in active_native.molar_volume_bounds(
                temperature_k,
                fractions,
                *packing_fraction_interval,
                1.0,
            )
        )
        if not lower_volume < volume < upper_volume:
            raise _IncipientFailure(
                "active_bound_ambiguity",
                "declared vapor volume left its packing-fraction chart",
            )
        block = _core.evaluate_parameterized_phase(
            provider_model(),
            temperature_k,
            fractions,
            volume,
            *packing_fraction_interval,
        )
    except _IncipientFailure:
        raise
    except (RuntimeError, TypeError, ValueError) as error:
        raise _IncipientFailure("provider_domain_rejection", str(error)) from error

    gradient = tuple(
        float(cast(Any, value)) for value in cast(Sequence[object], block["state_gradient"])
    )
    hessian = tuple(
        float(cast(Any, value)) for value in cast(Sequence[object], block["state_hessian"])
    )
    direct_mu = tuple(
        float(cast(Any, value))
        for value in cast(
            Sequence[object], block["chemical_potential_parameter_derivatives_over_rt"]
        )
    )
    pressure_parameter = tuple(
        float(cast(Any, value))
        for value in cast(Sequence[object], block["pressure_parameter_derivatives_pa"])
    )
    if (
        len(gradient) != state_count
        or len(hessian) != state_count * state_count
        or len(direct_mu) != component_count * len(active_ids)
        or len(pressure_parameter) != len(active_ids)
    ):
        raise _IncipientFailure(
            "owner_contract_failure", "vapor Provider derivative layout changed"
        )

    h_vv = hessian[volume_coordinate * state_count + volume_coordinate]
    h_scale = max(
        1.0,
        max(abs(hessian[volume_coordinate * state_count + index]) for index in range(state_count)),
    )
    if not math.isfinite(h_vv) or h_vv <= 1.0e-14 * h_scale:
        raise _IncipientFailure(
            "second_order_failure",
            "declared vapor root lacks strict mechanical stability",
        )

    reference = component_count - 1
    chart = tuple(
        tuple(
            fractions[component]
            * ((1.0 if component == coordinate else 0.0) - fractions[coordinate])
            for coordinate in range(reference)
        )
        for component in range(component_count)
    )
    mu_log_ratio_at_fixed_volume = tuple(
        tuple(
            math.fsum(
                hessian[component * state_count + species] * chart[species][coordinate]
                for species in range(component_count)
            )
            for coordinate in range(reference)
        )
        for component in range(component_count)
    )
    mu_log_volume = tuple(
        hessian[component * state_count + volume_coordinate] * volume
        for component in range(component_count)
    )
    rt = 8.31446261815324 * temperature_k
    pressure_log_ratio = tuple(
        -rt
        * math.fsum(
            hessian[volume_coordinate * state_count + component] * chart[component][coordinate]
            for component in range(component_count)
        )
        for coordinate in range(reference)
    )
    pressure_log_volume = -rt * h_vv * volume
    volume_log_ratio = tuple(
        -math.fsum(
            hessian[volume_coordinate * state_count + component] * chart[component][coordinate]
            for component in range(component_count)
        )
        / h_vv
        for coordinate in range(reference)
    )
    mu_log_ratio = tuple(
        tuple(
            mu_log_ratio_at_fixed_volume[component][coordinate]
            + hessian[component * state_count + volume_coordinate] * volume_log_ratio[coordinate]
            for coordinate in range(reference)
        )
        for component in range(component_count)
    )
    volume_pressure = -1.0 / (rt * h_vv)
    mu_pressure = tuple(
        hessian[component * state_count + volume_coordinate] * volume_pressure
        for component in range(component_count)
    )
    mu_parameters: dict[str, tuple[float, ...]] = {}
    for identity in all_active_identities:
        parameter = active_position.get(identity)
        if parameter is None:
            mu_parameters[identity] = (0.0,) * component_count
            continue
        volume_parameter = pressure_parameter[parameter] / (rt * h_vv)
        mu_parameters[identity] = tuple(
            direct_mu[component * len(active_ids) + parameter]
            + hessian[component * state_count + volume_coordinate] * volume_parameter
            for component in range(component_count)
        )
    return _IncipientVaporPoint(
        mole_fractions=fractions,
        volume_m3=volume,
        pressure_pa=float(cast(Any, block["pressure_pa"])),
        packing_fraction=float(cast(Any, block["packing_fraction"])),
        chemical_potential_over_rt=gradient[:component_count],
        mu_log_ratio_at_fixed_volume=mu_log_ratio_at_fixed_volume,
        mu_log_volume_derivative=mu_log_volume,
        mu_log_ratio_jacobian=mu_log_ratio,
        pressure_log_ratio_derivative=pressure_log_ratio,
        pressure_log_volume_derivative=pressure_log_volume,
        mu_pressure_derivative=mu_pressure,
        mu_parameter_derivatives=mu_parameters,
    )


def _solve_declared_incipient_vapor(
    *,
    model: epcsaft.Pure | epcsaft.Mixture,
    active_parameters: epcsaft.ActiveParameterSet | None,
    all_active_identities: tuple[str, ...],
    temperature: Any,
    temperature_k: float,
    pressure_pa: float,
    liquid_potentials: tuple[float, ...],
    continuation_fractions: tuple[float, ...] | None,
    continuation_volume_m3: float | None,
    packing_fraction_interval: tuple[float, float],
) -> _IncipientVaporPoint:
    """Solve the declared composition and low-density volume coordinates together."""

    component_count = len(liquid_potentials)
    reference = component_count - 1
    native_model: _core.NativeModel | None = None

    def provider_model() -> _core.NativeModel:
        nonlocal native_model
        if native_model is None:
            try:
                native_model = model._native_model(active_parameters)
            except (RuntimeError, TypeError, ValueError) as error:
                raise _IncipientFailure("provider_domain_rejection", str(error)) from error
        return native_model

    if continuation_fractions is not None:
        if continuation_volume_m3 is None:
            raise _IncipientFailure("branch_loss", "continued vapor state lacks its volume")
        starts = [continuation_fractions]
    else:
        shifted = tuple(value - max(liquid_potentials) for value in liquid_potentials)
        weights = tuple(math.exp(max(-700.0, value)) for value in shifted)
        weight_total = math.fsum(weights)
        starts = [tuple(value / weight_total for value in weights)]
        starts.append((1.0 / component_count,) * component_count)
        minor_fraction = 1.0e-4 / max(1, component_count - 1)
        starts.extend(
            tuple(
                1.0 - minor_fraction * (component_count - 1)
                if component == dominant
                else minor_fraction
                for component in range(component_count)
            )
            for dominant in range(component_count)
        )

    active_native = (
        model._native
        if active_parameters is None
        else model._native.with_active_coordinates(model._active_coordinates(active_parameters))
    )
    cache: dict[tuple[float, ...], _IncipientVaporPoint] = {}

    def evaluate(coordinates: tuple[float, ...]) -> _IncipientVaporPoint:
        if coordinates not in cache:
            cache[coordinates] = _evaluate_incipient_vapor_point(
                model=model,
                provider_model=provider_model,
                active_parameters=active_parameters,
                all_active_identities=all_active_identities,
                temperature_k=temperature_k,
                log_ratios=coordinates[:-1],
                log_volume=coordinates[-1],
                packing_fraction_interval=packing_fraction_interval,
            )
        return cache[coordinates]

    def residual(point: _IncipientVaporPoint) -> tuple[float, ...]:
        if point.pressure_pa <= 0.0:
            raise _IncipientFailure("no_incipient_phase", "declared vapor pressure is not positive")
        return (
            *(
                (
                    point.chemical_potential_over_rt[index]
                    - point.chemical_potential_over_rt[reference]
                )
                - (liquid_potentials[index] - liquid_potentials[reference])
                for index in range(reference)
            ),
            math.log(point.pressure_pa / pressure_pa),
        )

    def certify_local_minimum(point: _IncipientVaporPoint) -> None:
        if reference == 0:
            return
        chart = tuple(
            tuple(
                point.mole_fractions[component]
                * ((1.0 if component == coordinate else 0.0) - point.mole_fractions[coordinate])
                for coordinate in range(reference)
            )
            for component in range(component_count)
        )
        raw = tuple(
            tuple(
                math.fsum(
                    chart[component][row] * point.mu_log_ratio_jacobian[component][column]
                    for component in range(component_count)
                )
                for column in range(reference)
            )
            for row in range(reference)
        )
        scale = max(1.0, max(abs(value) for row in raw for value in row))
        if any(
            abs(raw[row][column] - raw[column][row]) > 2.0e-9 * scale
            for row in range(reference)
            for column in range(row)
        ):
            raise _IncipientFailure(
                "first_order_failure",
                "incipient tangent-plane Hessian is not symmetric",
            )
        matrix = [
            [0.5 * (raw[row][column] + raw[column][row]) for column in range(reference)]
            for row in range(reference)
        ]
        lower = [[0.0] * reference for _ in range(reference)]
        tolerance = 1.0e-12 * scale
        for row in range(reference):
            for column in range(row + 1):
                value = matrix[row][column] - math.fsum(
                    lower[row][inner] * lower[column][inner] for inner in range(column)
                )
                if row == column:
                    if not math.isfinite(value) or value <= tolerance:
                        raise _IncipientFailure(
                            "second_order_failure",
                            "incipient composition is not a strict local tangent-plane minimum",
                        )
                    lower[row][column] = math.sqrt(value)
                else:
                    lower[row][column] = value / lower[column][column]

    accepted: list[tuple[_IncipientVaporPoint, float]] = []
    last_failure: _IncipientFailure | None = None
    for fractions in starts:
        log_ratios = tuple(
            math.log(fractions[index] / fractions[reference]) for index in range(reference)
        )
        lower_volume, upper_volume = (
            float(value)
            for value in active_native.molar_volume_bounds(
                temperature_k,
                fractions,
                *packing_fraction_interval,
                1.0,
            )
        )
        proposed_volume = (
            continuation_volume_m3
            if continuation_fractions is not None
            else 8.31446261815324 * temperature_k / pressure_pa
        )
        assert proposed_volume is not None
        volume = min(
            math.nextafter(upper_volume, lower_volume),
            max(math.nextafter(lower_volume, upper_volume), proposed_volume),
        )
        coordinates = (*log_ratios, math.log(volume))
        point: _IncipientVaporPoint | None = None
        final_norm = math.inf
        for _iteration in range(36):
            try:
                point = evaluate(coordinates)
                values = residual(point)
            except _IncipientFailure as error:
                last_failure = error
                point = None
                break
            final_norm = max(abs(value) for value in values)
            if final_norm <= 2.0e-10:
                break
            jacobian = (
                *tuple(
                    (
                        *(
                            point.mu_log_ratio_at_fixed_volume[row][column]
                            - point.mu_log_ratio_at_fixed_volume[reference][column]
                            for column in range(reference)
                        ),
                        point.mu_log_volume_derivative[row]
                        - point.mu_log_volume_derivative[reference],
                    )
                    for row in range(reference)
                ),
                tuple(
                    derivative / point.pressure_pa
                    for derivative in (
                        *point.pressure_log_ratio_derivative,
                        point.pressure_log_volume_derivative,
                    )
                ),
            )
            step = _solve_incipient_linear_system(jacobian, tuple(-value for value in values))
            composition_step = max((abs(value) for value in step[:-1]), default=0.0)
            step_scale = min(
                1.0,
                2.0 / composition_step if composition_step > 2.0 else 1.0,
                0.1 / abs(step[-1]) if abs(step[-1]) > 0.1 else 1.0,
            )
            accepted_step = False
            for reduction in range(10):
                factor = step_scale * 0.5**reduction
                candidate = tuple(
                    value + factor * delta for value, delta in zip(coordinates, step, strict=True)
                )
                try:
                    evaluate(
                        tuple(
                            0.5 * (left + right)
                            for left, right in zip(coordinates, candidate, strict=True)
                        )
                    )
                    trial = evaluate(candidate)
                    trial_values = residual(trial)
                except _IncipientFailure:
                    continue
                trial_norm = max(abs(value) for value in trial_values)
                if trial_norm < final_norm:
                    coordinates = candidate
                    accepted_step = True
                    break
            if not accepted_step:
                last_failure = _IncipientFailure(
                    "numerical_convergence_failure",
                    "incipient composition-volume Newton step did not reduce the residual",
                )
                point = None
                break
        if point is None or final_norm > 2.0e-10:
            continue
        try:
            certify_local_minimum(point)
        except _IncipientFailure as error:
            last_failure = error
            continue
        if not any(
            max(
                abs(math.log(left / right))
                for left, right in zip(point.mole_fractions, existing.mole_fractions, strict=True)
            )
            <= 2.0e-8
            for existing, _ in accepted
        ):
            accepted.append((point, final_norm))
    if not accepted:
        if continuation_fractions is not None:
            detail = "continued incipient composition-volume corrector failed"
            if last_failure is not None:
                detail = str(last_failure)
            raise _IncipientFailure("branch_loss", detail)
        if last_failure is not None and last_failure.code in {
            "second_order_failure",
            "singular_sensitivity",
            "ill_conditioned_sensitivity",
        }:
            raise last_failure
        raise _IncipientFailure(
            "no_incipient_phase",
            "no declared strict-stable vapor composition satisfies fugacity equality",
        )
    if continuation_fractions is not None:
        return accepted[0][0]

    ranked = sorted(
        accepted,
        key=lambda item: (
            item[0].chemical_potential_over_rt[reference] - liquid_potentials[reference],
            item[0].mole_fractions,
        ),
    )
    best, best_norm = ranked[0]
    if len(ranked) > 1:
        next_point, next_norm = ranked[1]
        best_lambda = best.chemical_potential_over_rt[reference] - liquid_potentials[reference]
        next_lambda = (
            next_point.chemical_potential_over_rt[reference] - liquid_potentials[reference]
        )
        if abs(next_lambda - best_lambda) <= 2.0e-10 + max(best_norm, next_norm):
            raise _IncipientFailure(
                "multiple_unresolved_branches",
                "multiple strict local incipient minima have unresolved tangent-plane distance",
            )
    return best


def _nonideal_closure(
    problem: _FiniteIncipientKernel,
    temperature: Any,
    descriptor: Any,
    active_parameters: epcsaft.ActiveParameterSet | None,
    pressure_pa: float,
) -> _Closure:
    """Close the declared vapor against the exact certified reactive liquid.

    Phase count and roles are immutable inputs. Chemical equilibrium owns the
    liquid ``(T,n,V)`` state; this operation solves only the declared neutral
    vapor composition and its strict-stable Provider density branch.
    """

    equilibrium, trial, state_identity, amount_rows, volume_rows, balance_ids = _liquid_trial_state(
        problem, temperature, descriptor, active_parameters, pressure_pa
    )
    active_ids = descriptor.active_parameter_identities
    partner_parameters = problem.liquid_model._parameters.select(problem.vapor_component_ids)
    partner_model: epcsaft.Pure | epcsaft.Mixture = (
        epcsaft.Pure(partner_parameters)
        if len(problem.vapor_component_ids) == 1
        else epcsaft.Mixture(partner_parameters)
    )
    partner_active = _support_active_parameters(partner_parameters, active_parameters)
    if partner_model.component_ids != problem.vapor_component_ids:
        raise _IncipientFailure(
            "parameter_topology_mismatch",
            "declared vapor Provider does not match the ordered neutral support",
        )
    component_ids = problem.reaction_system.species_ids
    component_position = {identity: index for index, identity in enumerate(component_ids)}
    vapor_positions = tuple(
        component_position[identity] for identity in problem.vapor_component_ids
    )
    liquid_potentials = tuple(
        equilibrium.chemical_potential_over_rt[index] for index in vapor_positions
    )
    direction_ids = tuple(amount_rows)
    liquid_mu_directions = {
        identity: tuple(
            equilibrium.chemical_potential_directions[identity][component]
            for component in vapor_positions
        )
        for identity in direction_ids
    }

    vapor = _solve_declared_incipient_vapor(
        model=partner_model,
        active_parameters=partner_active,
        all_active_identities=active_ids,
        temperature=temperature,
        temperature_k=equilibrium.temperature_k,
        pressure_pa=pressure_pa,
        liquid_potentials=liquid_potentials,
        continuation_fractions=(
            None
            if problem.continuation_state is None
            else problem.continuation_state.phase(problem.vapor_phase_identity).mole_fractions
        ),
        continuation_volume_m3=(
            None
            if problem.continuation_state is None
            else problem.continuation_state.phase(
                problem.vapor_phase_identity
            ).molar_volume_m3_per_mol
        ),
        packing_fraction_interval=problem.liquid_packing_fraction_interval,
    )

    vapor_count = len(problem.vapor_component_ids)
    reference = vapor_count - 1
    g_jacobian = tuple(
        tuple(
            vapor.mu_log_ratio_jacobian[row][column]
            - vapor.mu_log_ratio_jacobian[reference][column]
            for column in range(reference)
        )
        for row in range(reference)
    )
    h_log_ratio = vapor.mu_log_ratio_jacobian[reference]
    chart = tuple(
        tuple(
            vapor.mole_fractions[component]
            * ((1.0 if component == coordinate else 0.0) - vapor.mole_fractions[coordinate])
            for coordinate in range(reference)
        )
        for component in range(vapor_count)
    )
    closure_directions: dict[str, float] = {}
    vapor_fraction_directions: dict[str, tuple[float, ...]] = {}
    for identity in direction_ids:
        vapor_base = (
            vapor.mu_pressure_derivative
            if identity == "pressure_pa"
            else vapor.mu_parameter_derivatives[identity]
            if identity in active_ids
            else (0.0,) * vapor_count
        )
        potential_direction = tuple(
            vapor_value - liquid_value
            for vapor_value, liquid_value in zip(
                vapor_base, liquid_mu_directions[identity], strict=True
            )
        )
        g_direction = tuple(
            potential_direction[index] - potential_direction[reference]
            for index in range(reference)
        )
        log_ratio_direction = _solve_incipient_linear_system(
            g_jacobian, tuple(-value for value in g_direction)
        )
        closure_directions[identity] = potential_direction[reference] + math.fsum(
            derivative * direction
            for derivative, direction in zip(h_log_ratio, log_ratio_direction, strict=True)
        )
        vapor_fraction_directions[identity] = tuple(
            math.fsum(
                chart[component][coordinate] * log_ratio_direction[coordinate]
                for coordinate in range(reference)
            )
            for component in range(vapor_count)
        )

    residual = vapor.chemical_potential_over_rt[reference] - liquid_potentials[reference]
    residual_log_pressure_derivative = pressure_pa * closure_directions["pressure_pa"]
    if not math.isfinite(residual_log_pressure_derivative):
        raise _IncipientFailure("first_order_failure", "bubble closure derivative is non-finite")

    # _row_values consumes the generic B(P,x(theta,P))-P Schur layout.  For
    # the reduced chemical-potential closure h=0, encode B_P-1=P*h_P and
    # B_theta=P*h_theta; the resulting total dP/dtheta is identical.
    boundary_directional = {
        "pressure_pa": 1.0 + residual_log_pressure_derivative,
        **{
            identity: pressure_pa * closure_directions[identity]
            for identity in direction_ids
            if identity != "pressure_pa"
        },
    }
    boundary_direct = {identity: 0.0 for identity in active_ids}
    vapor_direct = {identity: (0.0,) * vapor_count for identity in active_ids}

    total = math.fsum(equilibrium.amounts_mol)
    liquid_molar_volume = equilibrium.volume_m3 / total
    liquid_phase = EquilibriumPhaseState(
        role="liquid",
        support="all_components",
        supported_component_ids=component_ids,
        mole_fractions=equilibrium.mole_fractions,
        molar_volume_m3_per_mol=liquid_molar_volume,
        molar_density_mol_m3=1.0 / liquid_molar_volume,
        packing_fraction=equilibrium.packing_fraction,
        pressure_pa=equilibrium.pressure_pa,
        chemical_potential_component_ids=component_ids,
        chemical_potential_over_rt=equilibrium.chemical_potential_over_rt,
        mechanical_class="strict_stable_local_minimum",
    )
    vapor_phase = EquilibriumPhaseState(
        role="vapor",
        support="declared_neutral_components",
        supported_component_ids=problem.vapor_component_ids,
        mole_fractions=vapor.mole_fractions,
        molar_volume_m3_per_mol=vapor.volume_m3,
        molar_density_mol_m3=1.0 / vapor.volume_m3,
        packing_fraction=vapor.packing_fraction,
        pressure_pa=vapor.pressure_pa,
        chemical_potential_component_ids=problem.vapor_component_ids,
        chemical_potential_over_rt=vapor.chemical_potential_over_rt,
        mechanical_class="strict_stable_vapor",
    )
    return _Closure(
        pressure_pa,
        residual,
        residual_log_pressure_derivative,
        equilibrium,
        trial,
        state_identity,
        (liquid_phase, vapor_phase),
        boundary_directional,
        boundary_direct,
        vapor_fraction_directions,
        vapor_direct,
        partner_model.parameter_fingerprint,
        partner_model.component_ids,
        amount_rows,
        volume_rows,
        balance_ids,
    )


def _ideal_closure(
    problem: _FiniteIncipientKernel,
    temperature: Any,
    descriptor: Any,
    active_parameters: epcsaft.ActiveParameterSet | None,
    pressure_pa: float,
) -> _Closure:
    equilibrium, trial, state_identity, amount_rows, volume_rows, balance_ids = _liquid_trial_state(
        problem, temperature, descriptor, active_parameters, pressure_pa
    )
    active_ids = descriptor.active_parameter_identities
    amounts = equilibrium.amounts_mol
    component_ids = problem.reaction_system.species_ids
    component_index = {identity: index for index, identity in enumerate(component_ids)}
    vapor_indices = tuple(component_index[identity] for identity in problem.vapor_component_ids)
    fugacities = tuple(math.exp(equilibrium.log_fugacities[index]) for index in vapor_indices)
    boundary_pressure = math.fsum(fugacities)
    if not math.isfinite(boundary_pressure) or boundary_pressure <= 0.0:
        raise _IncipientFailure(
            "no_incipient_phase", "ideal comparison has no positive incipient pressure"
        )
    direction_ids = tuple(amount_rows)
    boundary_directional: dict[str, float] = {}
    vapor_directional: dict[str, tuple[float, ...]] = {}
    for identity in direction_ids:
        log_fugacity_direction = tuple(
            equilibrium.log_fugacity_directions[identity][component] for component in vapor_indices
        )
        boundary_derivative = math.fsum(
            fugacity * derivative
            for fugacity, derivative in zip(fugacities, log_fugacity_direction, strict=True)
        )
        boundary_directional[identity] = boundary_derivative
        vapor_directional[identity] = tuple(
            fugacity
            * (derivative * boundary_pressure - boundary_derivative)
            / (boundary_pressure * boundary_pressure)
            for fugacity, derivative in zip(fugacities, log_fugacity_direction, strict=True)
        )
    boundary_direct: dict[str, float] = {}
    vapor_direct: dict[str, tuple[float, ...]] = {}
    for identity in active_ids:
        mu_partial = tuple(
            equilibrium.log_fugacity_parameter_partials[identity][component]
            for component in vapor_indices
        )
        boundary_derivative = math.fsum(
            fugacity * derivative
            for fugacity, derivative in zip(fugacities, mu_partial, strict=True)
        )
        boundary_direct[identity] = boundary_derivative
        vapor_direct[identity] = tuple(
            fugacity
            * (derivative * boundary_pressure - boundary_derivative)
            / (boundary_pressure * boundary_pressure)
            for fugacity, derivative in zip(fugacities, mu_partial, strict=True)
        )
    total = math.fsum(amounts)
    molar_volume = equilibrium.volume_m3 / total
    liquid_phase = EquilibriumPhaseState(
        role="liquid",
        support="all_components",
        supported_component_ids=component_ids,
        mole_fractions=equilibrium.mole_fractions,
        molar_volume_m3_per_mol=molar_volume,
        molar_density_mol_m3=1.0 / molar_volume,
        packing_fraction=equilibrium.packing_fraction,
        pressure_pa=equilibrium.pressure_pa,
        chemical_potential_component_ids=component_ids,
        chemical_potential_over_rt=equilibrium.chemical_potential_over_rt,
        mechanical_class="strict_stable_local_minimum",
    )
    vapor_fractions = tuple(value / boundary_pressure for value in fugacities)
    vapor_volume = 8.31446261815324 * equilibrium.temperature_k / boundary_pressure
    vapor_phase = EquilibriumPhaseState(
        role="vapor",
        support="declared_neutral_components",
        supported_component_ids=problem.vapor_component_ids,
        mole_fractions=vapor_fractions,
        molar_volume_m3_per_mol=vapor_volume,
        molar_density_mol_m3=1.0 / vapor_volume,
        packing_fraction=0.0,
        pressure_pa=boundary_pressure,
        chemical_potential_component_ids=problem.vapor_component_ids,
        chemical_potential_over_rt=tuple(
            equilibrium.chemical_potential_over_rt[index] for index in vapor_indices
        ),
        mechanical_class="strict_stable_ideal_gas",
    )
    residual = math.log(boundary_pressure / pressure_pa)
    residual_derivative = (
        pressure_pa / boundary_pressure * boundary_directional["pressure_pa"] - 1.0
    )
    if not math.isfinite(residual_derivative):
        raise _IncipientFailure(
            "first_order_failure", "ideal bubble closure derivative is non-finite"
        )
    return _Closure(
        pressure_pa,
        residual,
        residual_derivative,
        equilibrium,
        trial,
        state_identity,
        (liquid_phase, vapor_phase),
        boundary_directional,
        boundary_direct,
        vapor_directional,
        vapor_direct,
        _canonical_fingerprint(
            "epcsaft.ideal-vapor-comparison",
            {
                "reference_id": problem.vapor_model.reference_id,
                "component_ids": problem.vapor_component_ids,
            },
        ),
        problem.vapor_component_ids,
        amount_rows,
        volume_rows,
        balance_ids,
    )


def _closure(
    problem: _FiniteIncipientKernel,
    temperature: Any,
    descriptor: Any,
    active_parameters: epcsaft.ActiveParameterSet | None,
    pressure_pa: float,
) -> _Closure:
    if isinstance(problem.vapor_model, _IdealVapor):
        return _ideal_closure(problem, temperature, descriptor, active_parameters, pressure_pa)
    return _nonideal_closure(problem, temperature, descriptor, active_parameters, pressure_pa)


def _fixed_pressure_coordinate_derivatives(
    closure: _Closure,
    balance_direction: tuple[float, ...],
) -> tuple[float, dict[str, float]]:
    balance_ids = closure.balance_direction_ids
    if len(balance_ids) != len(balance_direction) or not balance_ids:
        raise _IncipientFailure(
            "parameter_topology_mismatch",
            "fixed-pressure incipient closure lacks the conserved-balance directions",
        )
    vapor = closure.phases[1]
    residual_scale = (
        vapor.pressure_pa
        if vapor.mechanical_class == "strict_stable_ideal_gas"
        else closure.pressure_pa
    )
    if not math.isfinite(residual_scale) or residual_scale <= 0.0:
        raise _IncipientFailure("first_order_failure", "incipient residual scale is invalid")

    def residual_direction(identity: str) -> float:
        return (
            closure.boundary_directional_derivatives[identity]
            + closure.boundary_direct_derivatives.get(identity, 0.0)
        ) / residual_scale

    coordinate_slope = math.fsum(
        residual_direction(identity) * direction
        for identity, direction in zip(balance_ids, balance_direction, strict=True)
    )
    slope_scale = max(
        1.0,
        *(abs(residual_direction(identity)) for identity in balance_ids),
    )
    if abs(coordinate_slope) <= 1.0e-12 * slope_scale:
        raise _IncipientFailure(
            "singular_sensitivity", "fixed-pressure affine-feed closure is singular"
        )
    if abs(coordinate_slope) <= 1.0e-8 * slope_scale:
        raise _IncipientFailure(
            "ill_conditioned_sensitivity",
            "fixed-pressure affine-feed closure is ill-conditioned",
        )
    return coordinate_slope, {
        identity: residual_direction(identity)
        for identity in closure.boundary_directional_derivatives
    }


def _candidate_distance(
    left: _Closure,
    right: _Closure,
) -> float:
    left_vapor = left.phases[1]
    right_vapor = right.phases[1]
    return max(
        abs(math.log(left.pressure_pa / right.pressure_pa)),
        max(
            abs(math.log(a / b))
            for a, b in zip(left_vapor.mole_fractions, right_vapor.mole_fractions, strict=True)
        ),
    )


def _validate_same_liquid_state(closure: _Closure) -> None:
    """Verify the phase receipt carries the exact chemical ``(T,n,V)`` liquid."""

    equilibrium = closure.equilibrium
    liquid_phase = closure.phases[0]
    if equilibrium is None:
        raise _IncipientFailure("owner_contract_failure", "chemical result is missing")
    chemical_molar_volume = equilibrium.volume_m3 / math.fsum(equilibrium.amounts_mol)
    if (
        liquid_phase.role != "liquid"
        or liquid_phase.support != "all_components"
        or liquid_phase.supported_component_ids != equilibrium.species_ids
        or any(
            abs(actual - expected) > 2.0e-13
            for actual, expected in zip(
                liquid_phase.mole_fractions, equilibrium.mole_fractions, strict=True
            )
        )
        or not math.isclose(
            liquid_phase.molar_volume_m3_per_mol,
            chemical_molar_volume,
            rel_tol=2.0e-9,
            abs_tol=2.0e-14,
        )
        or not math.isclose(
            liquid_phase.pressure_pa,
            equilibrium.pressure_pa,
            rel_tol=2.0e-9,
            abs_tol=2.0e-5,
        )
    ):
        raise _IncipientFailure(
            "branch_loss",
            "coupled result does not bind to the certified chemical liquid branch",
        )


def _solve_starts(
    problem: _FiniteIncipientKernel,
    temperature: Any,
    descriptor: Any,
    active_parameters: epcsaft.ActiveParameterSet | None,
) -> tuple[
    tuple[_Closure, ...],
    tuple[NonEvaluableTrial, ...],
    int,
    int,
    int,
]:
    lower, upper = problem.pressure_interval_pa
    starts = (
        (problem.continuation_state.pressure_pa,)
        if problem.continuation_state is not None
        else problem.pressure_starts_pa
    )
    cache: dict[float, _Closure] = {}

    def evaluate(pressure: float) -> _Closure:
        if pressure not in cache:
            cache[pressure] = _closure(
                problem, temperature, descriptor, active_parameters, pressure
            )
        return cache[pressure]

    failures: list[NonEvaluableTrial] = []
    accepted: list[_Closure] = []
    accepted_attempts = 0
    for start in starts:
        pressure = start
        final: _Closure | None = None
        failure: _IncipientFailure | None = None
        for _iteration in range(18):
            try:
                current = evaluate(pressure)
            except _IncipientFailure as error:
                failure = error
                break
            if abs(current.residual) <= 2.0e-10:
                final = current
                break
            slope = current.residual_log_pressure_derivative
            if abs(slope) < 1.0e-10:
                failure = _IncipientFailure(
                    "singular_sensitivity", "bubble closure slope is singular"
                )
                break
            step = max(-1.0, min(1.0, -current.residual / slope))
            candidate_log = min(
                math.log(upper),
                max(math.log(lower), math.log(pressure) + step),
            )
            candidate = math.exp(candidate_log)
            if candidate == pressure:
                failure = _IncipientFailure(
                    "no_incipient_phase", "declared pressure interval excludes the Newton step"
                )
                break
            pressure = candidate
        if final is not None:
            accepted_attempts += 1
            duplicate = next(
                (
                    index
                    for index, existing in enumerate(accepted)
                    if _candidate_distance(final, existing) <= 2.0e-8
                ),
                None,
            )
            if duplicate is None:
                accepted.append(final)
            else:
                existing = accepted[duplicate]
                final_key = (
                    final.pressure_pa,
                    final.phases[1].mole_fractions,
                    final.phases[0].molar_volume_m3_per_mol,
                    final.phases[1].molar_volume_m3_per_mol,
                )
                existing_key = (
                    existing.pressure_pa,
                    existing.phases[1].mole_fractions,
                    existing.phases[0].molar_volume_m3_per_mol,
                    existing.phases[1].molar_volume_m3_per_mol,
                )
                if final_key < existing_key:
                    accepted[duplicate] = final
        else:
            code = "numerical_convergence_failure" if failure is None else failure.code
            diagnostic = (
                "bubble Newton iteration limit reached" if failure is None else str(failure)
            )
            failures.append(NonEvaluableTrial(code, diagnostic))
    provider_evaluations = sum(
        candidate.equilibrium.diagnostics.provider_evaluation_count + 1
        for candidate in cache.values()
    )
    return (
        tuple(accepted),
        tuple(failures),
        len(cache),
        provider_evaluations,
        accepted_attempts,
    )


def _select_branch(
    problem: _FiniteIncipientKernel,
    descriptor: Any,
    candidates: tuple[_Closure, ...],
    failures: tuple[NonEvaluableTrial, ...],
) -> _Closure:
    if not candidates:
        priority: tuple[NonEvaluableTrialCode, ...] = (
            "multiple_unresolved_branches",
            "branch_loss",
            "active_bound_ambiguity",
            "component_order_mismatch",
            "parameter_topology_mismatch",
            "provider_domain_rejection",
            "source_reference_derivative_unavailable",
            "source_reference_limit_failure",
            "source_reference_failure",
            "ill_conditioned_sensitivity",
            "singular_sensitivity",
            "first_order_failure",
            "physical_rejection",
            "no_incipient_phase",
            "numerical_convergence_failure",
        )
        diagnostics = {failure.code: failure.diagnostic for failure in failures}
        for code in priority:
            if code in diagnostics:
                raise _IncipientFailure(code, diagnostics[code])
        raise _IncipientFailure("no_incipient_phase", "no declared pressure start found a branch")
    reference = problem.continuation_state
    if reference is None:
        if len(candidates) != 1:
            raise _IncipientFailure(
                "multiple_unresolved_branches", "declared starts found multiple bubble branches"
            )
        return candidates[0]
    first_trial = candidates[0].trial
    if (
        first_trial is None
        or reference.problem_fingerprint != descriptor.problem_fingerprint
        or reference.topology_fingerprint != descriptor.topology_fingerprint
        or reference.parameter_fingerprint != first_trial.parameter_fingerprint
    ):
        raise _IncipientFailure(
            "parameter_topology_mismatch",
            "finite-incipient continuation state does not match this trial topology",
        )
    if len(candidates) != 1:
        raise _IncipientFailure(
            "multiple_unresolved_branches",
            "continued incipient solve returned more than one anchor-seeded branch",
        )
    return candidates[0]


def _row_values(
    problem: _FiniteIncipientKernel,
    descriptor: Any,
    closure: _Closure,
    derivative_ids: tuple[str, ...] | None = None,
    fixed_balance_direction: tuple[float, ...] | None = None,
) -> tuple[tuple[float, tuple[float, ...]], ...]:
    direction_ids = (
        descriptor.active_parameter_identities if derivative_ids is None else derivative_ids
    )
    equilibrium = closure.equilibrium
    assert equilibrium is not None
    amounts = equilibrium.amounts_mol
    volume = equilibrium.volume_m3
    liquid_amount_total: dict[str, tuple[float, ...]] = {}
    liquid_volume_total: dict[str, float] = {}
    vapor = closure.phases[1]
    vapor_derivatives: dict[str, tuple[float, ...]] = {}
    if fixed_balance_direction is None:
        pressure_partial = closure.boundary_directional_derivatives["pressure_pa"]
        schur_pressure = pressure_partial - 1.0
        schur_scale = max(1.0, abs(pressure_partial))
        if not math.isfinite(schur_pressure) or abs(schur_pressure) < 1.0e-12 * schur_scale:
            raise _IncipientFailure(
                "singular_sensitivity", "coupled bubble Schur complement is singular"
            )
        if abs(schur_pressure) < 1.0e-8 * schur_scale:
            raise _IncipientFailure(
                "ill_conditioned_sensitivity",
                "coupled bubble Schur complement is ill-conditioned",
            )
        pressure_derivatives = {
            identity: -(
                closure.boundary_direct_derivatives.get(identity, 0.0)
                + closure.boundary_directional_derivatives[identity]
            )
            / schur_pressure
            for identity in direction_ids
        }
        for identity in direction_ids:
            pressure_derivative = pressure_derivatives[identity]
            liquid_amount_total[identity] = tuple(
                direct + pressure * pressure_derivative
                for direct, pressure in zip(
                    closure.chemical_amount_directions[identity],
                    closure.chemical_amount_directions["pressure_pa"],
                    strict=True,
                )
            )
            liquid_volume_total[identity] = (
                closure.chemical_volume_directions[identity]
                + closure.chemical_volume_directions["pressure_pa"] * pressure_derivative
            )
            vapor_derivatives[identity] = tuple(
                direct + directional + pressure * pressure_derivative
                for direct, directional, pressure in zip(
                    closure.vapor_fraction_direct_derivatives.get(
                        identity, (0.0,) * len(vapor.mole_fractions)
                    ),
                    closure.vapor_fraction_directional_derivatives[identity],
                    closure.vapor_fraction_directional_derivatives["pressure_pa"],
                    strict=True,
                )
            )
    else:
        if derivative_ids is not None:
            raise _IncipientFailure(
                "parameter_topology_mismatch",
                "fixed-pressure rows require active parameter directions",
            )
        slope, residual_directions = _fixed_pressure_coordinate_derivatives(
            closure, fixed_balance_direction
        )
        coordinate_derivatives = {
            identity: -residual_directions[identity] / slope for identity in direction_ids
        }
        pressure_derivatives = {identity: 0.0 for identity in direction_ids}
        balance_ids = closure.balance_direction_ids
        coordinate_amount = tuple(
            math.fsum(
                closure.chemical_amount_directions[identity][component] * direction
                for identity, direction in zip(balance_ids, fixed_balance_direction, strict=True)
            )
            for component in range(len(problem.reaction_system.species_ids))
        )
        coordinate_volume = math.fsum(
            closure.chemical_volume_directions[identity] * direction
            for identity, direction in zip(balance_ids, fixed_balance_direction, strict=True)
        )
        coordinate_vapor = tuple(
            math.fsum(
                closure.vapor_fraction_directional_derivatives[identity][component] * direction
                for identity, direction in zip(balance_ids, fixed_balance_direction, strict=True)
            )
            for component in range(len(vapor.mole_fractions))
        )
        for identity in direction_ids:
            dq = coordinate_derivatives[identity]
            liquid_amount_total[identity] = tuple(
                direct + coordinate * dq
                for direct, coordinate in zip(
                    closure.chemical_amount_directions[identity],
                    coordinate_amount,
                    strict=True,
                )
            )
            liquid_volume_total[identity] = (
                closure.chemical_volume_directions[identity] + coordinate_volume * dq
            )
            vapor_derivatives[identity] = tuple(
                direct + directional + coordinate * dq
                for direct, directional, coordinate in zip(
                    closure.vapor_fraction_direct_derivatives.get(
                        identity, (0.0,) * len(vapor.mole_fractions)
                    ),
                    closure.vapor_fraction_directional_derivatives[identity],
                    coordinate_vapor,
                    strict=True,
                )
            )
    results: list[tuple[float, tuple[float, ...]]] = []
    for row in descriptor.rows:
        if row.selector == "system.pressure":
            value = closure.pressure_pa
            jacobian = tuple(pressure_derivatives[identity] for identity in direction_ids)
        elif (
            row.phase_identity == problem.vapor_phase_identity
            and row.selector == "phase.mole_fraction"
        ):
            value = math.fsum(
                coefficient * fraction
                for coefficient, fraction in zip(
                    row.coefficients, vapor.mole_fractions, strict=True
                )
            )
            jacobian = tuple(
                math.fsum(
                    coefficient * derivative
                    for coefficient, derivative in zip(
                        row.coefficients, vapor_derivatives[identity], strict=True
                    )
                )
                for identity in direction_ids
            )
        elif (
            row.phase_identity == problem.vapor_phase_identity
            and row.selector == "phase.partial_pressure"
        ):
            vapor_fraction = math.fsum(
                coefficient * fraction
                for coefficient, fraction in zip(
                    row.coefficients, vapor.mole_fractions, strict=True
                )
            )
            value = closure.pressure_pa * vapor_fraction
            jacobian = tuple(
                pressure_derivatives[identity] * vapor_fraction
                + closure.pressure_pa
                * math.fsum(
                    coefficient * derivative
                    for coefficient, derivative in zip(
                        row.coefficients, vapor_derivatives[identity], strict=True
                    )
                )
                for identity in direction_ids
            )
        else:
            try:
                value, jacobian = _contract_phase_output(
                    row.selector,
                    row.coefficients,
                    row.solvent_mass_coefficients_kg_per_mol,
                    amounts,
                    volume,
                    tuple(liquid_amount_total[identity] for identity in direction_ids),
                    tuple(liquid_volume_total[identity] for identity in direction_ids),
                )
            except ValueError as error:
                raise _IncipientFailure("observation_out_of_support", str(error)) from error
        if not math.isfinite(value) or any(not math.isfinite(item) for item in jacobian):
            raise _IncipientFailure(
                "owner_contract_failure", "finite-incipient output returned non-finite data"
            )
        if row.support == "positive" and value <= 0.0:
            raise _IncipientFailure(
                "observation_out_of_support", "positive finite-incipient output is not positive"
            )
        results.append((value, jacobian))
    return tuple(results)


def _physical_evidence(closure: _Closure) -> tuple[tuple[str, float | int | str], ...]:
    diagnostics = closure.equilibrium.diagnostics
    liquid_phase, vapor_phase = closure.phases
    potential_residual = max(
        abs(left - right)
        for left, right in zip(
            (
                liquid_phase.chemical_potential_over_rt[
                    closure.equilibrium.species_ids.index(identity)
                ]
                for identity in vapor_phase.chemical_potential_component_ids
            ),
            vapor_phase.chemical_potential_over_rt,
            strict=True,
        )
    )
    phase_pressure_residual = abs(liquid_phase.pressure_pa - vapor_phase.pressure_pa)
    if diagnostics.local_minimum_status != "passed" or any(
        "strict_stable" not in phase.mechanical_class for phase in closure.phases
    ):
        raise _IncipientFailure(
            "second_order_failure",
            "coupled bubble state lacks strict local chemical/mechanical stability",
        )
    required = (
        diagnostics.balance_inf_norm,
        diagnostics.charge_inf_norm,
        diagnostics.reaction_affinity_inf_norm,
        diagnostics.pressure_relative_residual,
    )
    if any(value is None or not math.isfinite(value) for value in required):
        raise _IncipientFailure("physical_rejection", "coupled physical residuals are incomplete")
    if (
        abs(closure.residual) > 2.0e-10
        or abs(closure.residual_log_pressure_derivative) < 1.0e-10
        or not math.isfinite(potential_residual)
        or not math.isfinite(phase_pressure_residual)
    ):
        raise _IncipientFailure("first_order_failure", "coupled bubble closure is uncertified")
    return (
        ("bubble_log_closure_abs", abs(closure.residual)),
        ("bubble_log_pressure_derivative", closure.residual_log_pressure_derivative),
        ("balance_inf_norm", diagnostics.balance_inf_norm),
        ("electroneutrality_inf", diagnostics.charge_inf_norm),
        ("reaction_affinity_inf_norm", diagnostics.reaction_affinity_inf_norm),
        ("pressure_relative_inf_norm", diagnostics.pressure_relative_residual),
        ("phase_pressure_residual_inf_pa", phase_pressure_residual),
        ("phase_chemical_potential_residual_inf", potential_residual),
    )


def _fixed_pressure_coordinate_result(
    problem: _FiniteIncipientKernel,
    descriptor: Any,
    closure: _Closure,
    balance_direction: tuple[float, ...],
) -> tuple[
    tuple[tuple[float, tuple[float, ...]], ...],
    tuple[tuple[str, float | int | str], ...],
]:
    """Return generic output values and evidence for one fixed-pressure state."""

    _validate_same_liquid_state(closure)
    rows = _row_values(problem, descriptor, closure, fixed_balance_direction=balance_direction)
    return rows, _physical_evidence(closure)


def _evaluate_incipient(
    problem: _FiniteIncipientKernel,
    temperature: Any,
    descriptor: Any,
    *,
    active_parameters: epcsaft.ActiveParameterSet | None = None,
) -> tuple[
    _Closure | None,
    NonEvaluableTrial | None,
    tuple[tuple[float, tuple[float, ...]], ...],
    tuple[tuple[float, tuple[float, ...]], ...],
    tuple[tuple[str, float | int | str], ...],
]:
    """Solve one locally certified finite–incipient state for generic GREPE."""

    candidates, failures, closure_count, provider_count, accepted_count = _solve_starts(
        problem, temperature, descriptor, active_parameters
    )
    try:
        closure = _select_branch(problem, descriptor, candidates, failures)
        _validate_same_liquid_state(closure)
        rows = _row_values(problem, descriptor, closure)
        balance_rows = _row_values(problem, descriptor, closure, closure.balance_direction_ids)
        evidence = _physical_evidence(closure) + (
            ("closure_evaluation_count", closure_count),
            ("provider_evaluation_count", provider_count),
            ("accepted_candidate_count", accepted_count),
            ("distinct_branch_count", len(candidates)),
        )
    except _IncipientFailure as error:
        return None, NonEvaluableTrial(error.code, str(error)), (), (), ()
    return closure, None, rows, balance_rows, evidence
