"""Declare native phase, reaction and solid equilibrium calculations.

Use ``Problem`` with explicit temperature (K), pressure (Pa), feed basis and
phase kinds, then call ``solve_equilibrium`` with a fluid model. Solved free phase amounts use mol; pins use dimensionless feed fractions. Total enthalpy and entropy constraints use J
and J/K, respectively. A returned stopping reason and equation diagnostics
must be inspected before using a solution. Finite phase discovery without a
descent witness is not a global stability certificate (ALG-STABILITY).

Examples
--------
>>> from epcsaft.equilibrium import Free, Pinned
>>> Free(300.0).guess
300.0
>>> Pinned(0.0).value
0.0
"""

from .. import _core
from .problem import (
    Amounts,
    IdealMoleFraction,
    effective_water_amounts,
    reaction_activities,
    Continuation,
    Discover,
    Feed,
    Free,
    MassFractions,
    Molalities,
    MoleFractions,
    NeutralReference,
    Pinned,
    Phase,
    Problem,
    Reaction,
    Solid,
    bubble_point,
    compile_problem,
    continue_equilibrium,
    dew_point,
    problem_from_mapping,
    solve_equilibrium,
    solve_equilibrium_starts,
    solved_state_actions,
)

EquilibriumStart = _core.EquilibriumStart
EquilibriumWorkLimits = _core.EquilibriumWorkLimits
EquilibriumStartsResult = _core.EquilibriumStartsResult
inventory_bounds = _core.inventory_bounds
InventoryBound = _core.InventoryBound

EquilibriumResult = _core.EquilibriumResult
ContinuationResult = _core.ContinuationResult
ContinuationStop = _core.ContinuationStop
CriticalDeclaration = _core.CriticalDeclaration
EnvelopeCondition = _core.EnvelopeCondition
EnvelopeCoordinate = _core.EnvelopeCoordinate
ReactionLogPolynomial = _core.ReactionLogPolynomial
ReactionReference = _core.ReactionReference
ReferenceSourceBasis = _core.ReferenceSourceBasis
SpeciesReference = _core.SpeciesReference
NeutralReferenceStatus = _core.NeutralReferenceStatus
SolvedStateActionBatch = _core.SolvedStateActionBatch
SolvedStateActionDiagnostics = _core.SolvedStateActionDiagnostics
SolvedStateActionDirection = _core.SolvedStateActionDirection
SolvedStateActionRequest = _core.SolvedStateActionRequest
SolvedStateActionResult = _core.SolvedStateActionResult
SolvedStateActionStatus = _core.SolvedStateActionStatus
SolvedStateObservable = _core.SolvedStateObservable
SolvedStateObservableKind = _core.SolvedStateObservableKind
EOS_STANDARD_STATE_ID = _core.EOS_STANDARD_STATE_ID
IDEAL_MOLE_FRACTION_STANDARD_STATE_ID = _core.IDEAL_MOLE_FRACTION_STANDARD_STATE_ID
IDEAL_MOLARITY_STANDARD_STATE_ID = _core.IDEAL_MOLARITY_STANDARD_STATE_ID
solved_state_first_actions = _core.solved_state_first_actions
asl_record = _core.equilibrium_asl_record
cooling_onset = _core.cooling_onset
stability_test = _core.stability_test
StabilityOptions = _core.StabilityOptions

__all__ = (
    "Amounts",
    "IdealMoleFraction",
    "effective_water_amounts",
    "reaction_activities",
    "IDEAL_MOLE_FRACTION_STANDARD_STATE_ID",
    "IDEAL_MOLARITY_STANDARD_STATE_ID",
    "inventory_bounds",
    "InventoryBound",
    "cooling_onset",
    "stability_test",
    "StabilityOptions",
    "asl_record",
    "Continuation",
    "ContinuationResult",
    "ContinuationStop",
    "CriticalDeclaration",
    "EnvelopeCondition",
    "EnvelopeCoordinate",
    "Discover",
    "EOS_STANDARD_STATE_ID",
    "EquilibriumResult",
    "EquilibriumStart",
    "EquilibriumWorkLimits",
    "EquilibriumStartsResult",
    "Feed",
    "Free",
    "MassFractions",
    "Molalities",
    "MoleFractions",
    "NeutralReference",
    "NeutralReferenceStatus",
    "SpeciesReference",
    "Pinned",
    "Phase",
    "Problem",
    "Reaction",
    "ReactionLogPolynomial",
    "ReactionReference",
    "ReferenceSourceBasis",
    "Solid",
    "SolvedStateActionBatch",
    "SolvedStateActionDiagnostics",
    "SolvedStateActionDirection",
    "SolvedStateActionRequest",
    "SolvedStateActionResult",
    "SolvedStateActionStatus",
    "SolvedStateObservable",
    "SolvedStateObservableKind",
    "bubble_point",
    "compile_problem",
    "continue_equilibrium",
    "dew_point",
    "problem_from_mapping",
    "solve_equilibrium",
    "solve_equilibrium_starts",
    "solved_state_actions",
    "solved_state_first_actions",
)
