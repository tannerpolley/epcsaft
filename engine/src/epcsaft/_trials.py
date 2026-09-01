from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, get_args

type NonEvaluableTrialCode = Literal[
    "component_order_mismatch",
    "parameter_topology_mismatch",
    "unsupported_topology",
    "unsupported_derivative",
    "eos_domain_rejection",
    "source_reference_failure",
    "source_reference_derivative_unavailable",
    "source_reference_limit_failure",
    "numerical_convergence_failure",
    "physical_rejection",
    "first_order_failure",
    "singular_sensitivity",
    "ill_conditioned_sensitivity",
    "active_bound_ambiguity",
    "second_order_failure",
    "multiple_unresolved_branches",
    "branch_loss",
    "no_incipient_phase",
    "unsupported_observable",
    "observation_out_of_support",
    "owner_contract_failure",
    "native_residual_failure",
    "optimizer_failure",
]

NON_EVALUABLE_TRIAL_CODES = frozenset(get_args(NonEvaluableTrialCode.__value__))


@dataclass(frozen=True, slots=True)
class NonEvaluableTrial:
    """Closed scientific failure code with non-normative diagnostic detail."""

    code: NonEvaluableTrialCode
    diagnostic: str

    def __post_init__(self) -> None:
        if self.code not in NON_EVALUABLE_TRIAL_CODES:
            raise ValueError("non-evaluable trial code is outside the public vocabulary")
        if not self.diagnostic:
            raise ValueError("non-evaluable trial diagnostic must be nonempty")


__all__ = (
    "NON_EVALUABLE_TRIAL_CODES",
    "NonEvaluableTrial",
    "NonEvaluableTrialCode",
)
