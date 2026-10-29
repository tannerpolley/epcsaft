"""Native thermodynamic declarations and results.

Notes
-----
All component indices are zero based in declared component_ids order. Molar
quantities use mol as their amount basis; mass-based inputs name their solvent
or mixture basis explicitly. Arrays are ordinary Python sequences.

Every native enum provides name (the symbolic string) and value (the native
integer code). These introspection properties have no physical unit. Individual
enum declarations link here rather than duplicating this convention.
"""
from typing import Any, Mapping, Optional, Sequence, Union, overload

GAS_CONSTANT_J_PER_MOL_K: float
"""Molar gas constant in J/(mol K), fixed by the native physical_constants owner."""
IDEAL_MOLE_FRACTION_STANDARD_STATE_ID: str
IDEAL_MOLARITY_STANDARD_STATE_ID: str
EOS_STANDARD_STATE_ID: str
"""Native reaction activity identity: (rho*x/(1 mol/m^3))*exp(mu_res/(RT)), dimensionless."""

class ParameterError(ValueError):
    """A parameter value, unit, shape or declaration is inadmissible."""
    ...
class EosError(RuntimeError):
    """A native thermodynamic evaluation cannot supply the requested result."""
    ...
class IonicReferenceBackground(ParameterError):
    pass

class DensityStructureError(EosError):
    """The density search cannot establish the requested root structure.

    Notes
    -----
    A failed or incomplete scan is not evidence of a unique stable root. See ALG-DENSITY in engine/docs/science/algorithms.md.
    """
    ...
    record: Optional[Convergence]
    """Final native refusal, or None when unmeasured."""
    first_numerical_failure: Optional[Convergence]
    """Earlier native refusal, when retained separately from final convergence."""


class ConvergenceError(EosError):
    """A native iterative calculation stops without an admitted result.

    Notes
    -----
    The exception message retains the residual, stopping criterion and available convergence diagnostics.
    """
    ...
    record: Convergence
    """Final native refusal, or None when unmeasured."""
    first_numerical_failure: Optional[Convergence]
    """Earlier native refusal, when retained separately from final convergence."""


class DeclarationError(ParameterError):
    """A required physical declaration is absent or inadmissible."""
    ...
    record: Optional[Convergence]
    """Final native refusal, or None when unmeasured."""
    first_numerical_failure: Optional[Convergence]
    """Earlier native refusal, when retained separately from final convergence."""


class Stop:
    """Reason an iterative calculation stopped.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    converged: Stop
    """Residual meets the requested criterion."""
    converged_to_floor: Stop
    """Residual meets the estimated numerical floor criterion; the requested tolerance need not be
    met.
    """
    stalled: Stop
    """Residual contraction does not support further progress."""
    budget_exhausted: Stop
    """An iteration, evaluation or elapsed-time limit was reached."""
    left_domain: Stop
    """The calculation left its admitted physical or numerical domain."""
    singular: Stop
    """A required linear system is singular."""
    infeasible: Stop
    """The declared constraints could not be satisfied by this calculation; no global infeasibility
    proof is implied.
    """
    nonfinite: Stop
    """A required value or diagnostic is nonfinite."""

class Category:
    """Provenance category of a numerical criterion.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    ThermodynamicConstraint: Category
    """A constraint imposed by the thermodynamic formulation."""
    AlgorithmTolerance: Category
    """An algorithmic stopping or admission tolerance."""

class Tolerance:
    """Named immutable numerical criterion and its source."""
    name: str
    """Policy identifier, not a physical quantity."""
    value: float
    """Criterion magnitude in the residual's units or dimensionless normalization; no common unit
    applies to all policies.
    """
    category: Category
    """Whether the criterion expresses a thermodynamic constraint or an algorithm tolerance."""
    source: str
    """Source citation or derivation identifying the criterion."""
    factor: float
    """Declared amplification factor, dimensionless."""
    upstream: Optional[Tolerance]
    """Parent criterion from which this value derives, or None for an independent criterion."""

def tolerances() -> list[Tolerance]:
    """Read immutable source defaults for numerical policies.

    Returns
    -------
    list[Tolerance]
        List of Tolerance records in native policy order; effective caller settings are reported by each result.solver_options.
    """
    ...

class Convergence:
    """Residual and work evidence for one iterative calculation."""
    residual: float
    """Final residual magnitude in the producing calculation's units or normalization."""
    criterion: Optional[Tolerance]
    """Effective requested tolerance, or None when no criterion is attached."""
    floor: Optional[float]
    """Estimated attainable residual floor in residual units, or None when not estimated."""
    amplification: Optional[float]
    """Estimated dimensionless numerical error amplification, or None when not estimated."""
    contraction: Optional[float]
    """Dimensionless geometric residual contraction per iteration, or None without sufficient
    history.
    """
    iterations: int
    """Number of iterations performed, a count."""
    evaluations: int
    """Number of evaluations performed, a count; its scope is the producing calculation."""
    stop: Stop
    """Stopping reason; convergence to a floor is distinct from meeting the requested tolerance."""
    @property
    def residual_over_criterion(self) -> Optional[float]:
        """Read the final residual relative to the requested criterion.

        Returns
        -------
        Optional[float]
            Dimensionless residual/criterion.value, or None without an attached criterion.
        """
        ...
    @property
    def residual_over_floor(self) -> Optional[float]:
        """Read the final residual relative to the estimated floor.

        Returns
        -------
        Optional[float]
            Dimensionless residual/floor, or None without an estimated floor.
        """
        ...

    cause: str
    """Producing operation and actual refused operands; empty without an observation."""

    rank: Optional[NumericalRank]
    """Rank in the named matrix basis, or None without a factor observation."""

    backward_error: Optional[float]
    """Dimensionless measured linear backward error, or None when unavailable."""


class Event:
    """One named numerical probe within a captured span."""
    name: str
    """Probe identifier supplied by the native producer."""
    span: int
    """Zero-based index into Trace.spans identifying the parent span."""
    iteration: int
    """Producer's iteration counter; zero may denote an uniterated probe."""
    values: list[float]
    """Probe values in producer-defined order and units; there is no common physical basis for
    different event names.
    """
    reason: str
    """Producer's explanation, possibly empty."""

class Span:
    """Timing and convergence of one captured native operation."""
    name: str
    """Native operation identifier."""
    parent: Optional[int]
    """Zero-based parent index in Trace.spans, or None for a root span."""
    convergence: Optional[Convergence]
    """Attached convergence evidence, or None if this operation did not attach a record."""
    wall_us: float
    """Elapsed wall time in microseconds."""

class Trace:
    """Ordered native timing spans and numerical probes."""
    spans: list[Span]
    """Spans in creation order; parent indices refer to this list."""
    events: list[Event]
    """Events in observation order; span indices refer to spans."""
    current: Optional[int]
    """Currently active span index, or None when no span is open."""
    incomplete: bool
    """True if opening the requested Ipopt log failed; diagnostic allocation failure propagates normally."""

    ipopt_log_error: str
    """Opening error for the requested Ipopt destination; empty without an error."""

class Capture:
    """Capture native numerical traces in a context manager.

    Notes
    -----
    Capture observes calculations without changing their stopping criteria. Detailed probes are enabled by diagnostics; exceptions from the enclosed calculation are not suppressed.
    """
    def __init__(self, diagnostics: bool = False, ipopt_log: str = "") -> None:
        """Initialize Capture from the declared inputs.

        Parameters
        ----------
        diagnostics : bool
            Enable detailed numerical probes in addition to operation spans.
        ipopt_log : str
            Optional absolute destination for the original calculation's sequential Ipopt log.

        Returns
        -------
        None
            None.
        """
        ...
    trace: Optional[Trace]
    """Completed trace, or None until the context closes."""
    spans: list[Span]
    """Completed trace spans; access before closure raises RuntimeError."""
    events: list[Event]
    """Completed trace events; access before closure raises RuntimeError."""
    def __enter__(self) -> Capture:
        """Begin recording native operations.

        Returns
        -------
        Capture
            This capture instance; its completed trace is unavailable until closure.
        """
        ...
    def __exit__(self, exc_type: object, exc_value: object, traceback: object) -> bool:
        """Close capture and retain the recorded trace.

        Parameters
        ----------
        exc_type : object
            Exception type supplied by the context-manager protocol, or None.
        exc_value : object
            Exception instance supplied by the protocol, or None.
        traceback : object
            Exception traceback supplied by the protocol, or None.

        Returns
        -------
        bool
            False, so an enclosed exception is propagated.
        """
        ...

class ParameterEntry:
    """One source-identified parameter magnitude and declared unit."""
    identity: str
    """Source identity used to identify this entry in diagnostics."""
    field: str
    """Native parameter field name; it selects the required index arity and unit."""
    index: list[int]
    """Zero-based component or site indices in document order; scalar fields use an empty list."""
    magnitude: float
    """Declared value in unit; native construction performs the admitted conversion."""
    unit: str
    """Explicit unit string accepted for field, such as angstrom, kelvin or dimensionless."""
    def __init__(self) -> None:
        """Initialize ParameterEntry from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...

class ParameterDocument:
    """Ordered component, site, model-choice and coefficient declarations.

    Notes
    -----
    build_parameters owns admission; an omitted required coefficient is not inferred from an application name.
    """
    component_ids: list[str]
    """Component identities in the order used by all component arrays, length n."""
    site_component: list[int]
    """Owning component index for each association site, length s."""
    site_multiplicity: list[int]
    """Positive count of equivalent sites for each site type, length s."""
    base_choice: str
    """Declared base Helmholtz model identifier."""
    association_choice: str
    """Declared association model identifier."""
    electrolyte_choice: str
    """Declared electrolyte contribution identifier."""
    ionic_dispersion_choice: str
    """Declared ionic dispersion rule identifier."""
    permittivity_choice: str
    """Declared bulk relative-permittivity model identifier."""
    polar_choice: str
    """Declared polar contribution identifier."""
    entries: list[ParameterEntry]
    """Source-identified scalar entries with explicit indices and units."""
    def __init__(self) -> None:
        """Initialize ParameterDocument from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...

class Parameters:
    """Admitted native EOS parameters in declared component and site order.

    Notes
    -----
    Lengths used throughout this module are n = component_count, s = site_count, and a = the number of active parameter directions. Construction uses build_parameters, not direct mutation.
    """
    ideal_gas: bool
    @property
    def component_count(self) -> int:
        """Read the number of declared components.

        Returns
        -------
        int
            Dimensionless component count n.
        """
        ...
    @property
    def site_count(self) -> int:
        """Read the number of declared association-site types.

        Returns
        -------
        int
            Dimensionless site-type count s.
        """
        ...
    @property
    def molar_masses(self) -> list[float]:
        """Read component molar masses.

        Returns
        -------
        list[float]
            Length n in component_ids order, kg/mol; quantities requiring mass refuse missing or nonpositive occupied-species masses.
        """
        ...
    @property
    def temperature_law_coefficient_identities(self) -> list[str]:
        """Read identities of native temperature-law coefficient slots.

        Returns
        -------
        list[str]
            Coefficient identity strings in native slot order.
        """
        ...

def build_parameters(document: ParameterDocument) -> Parameters:
    """Admit a source-identified parameter document.

    Parameters
    ----------
    document : ParameterDocument
        Complete physical declarations, ordered identities and coefficients with explicit units.

    Returns
    -------
    Parameters
        Native Parameters in the document component and site order.

    Notes
    -----
    Missing or inconsistent fields, shapes, units and physical values raise ParameterError; required coefficients are not inferred.
    """
    ...
def prepare_thermochemistry(record: ThermochemistryRecord, source_gas_constant: float) -> ThermochemistryRecord:
    """Convert source NASA9 coefficients to the Engine gas-constant convention.

    Parameters
    ----------
    record : ThermochemistryRecord
        Source thermochemistry in component order.
    source_gas_constant : float
        Positive gas constant used by the source, J/(mol K).

    Returns
    -------
    ThermochemistryRecord
        Prepared thermal record; source intervals and pressure standards are retained.

    Notes
    -----
    Every interval must contain NASA9; other ideal forms are refused by this preparation operation. Each coefficient is multiplied by source_gas_constant/GAS_CONSTANT_J_PER_MOL_K. Intervals and source pressures remain unchanged; no entropy pressure transfer or parameter fit occurs here.
    """
    ...
def fitted_document(start: ParameterDocument, fitted: Parameters) -> ParameterDocument:
    """Write fitted native parameter magnitudes into a starting source document.

    Parameters
    ----------
    start : ParameterDocument
        Source document retaining ordered identities and unit declarations.
    fitted : Parameters
        Fitted native Parameters to reproduce.

    Returns
    -------
    ParameterDocument
        Document that reconstructs fitted bitwise through build_parameters.

    Notes
    -----
    Refuses a document that cannot represent the supplied parameters exactly; source identities are not changed into evidence of physical validation.
    """
    ...

class IdealPolynomial:
    """Ideal-gas centered polynomial molar heat-capacity correlation.

    Notes
    -----
    See cp_ideal_polynomial in engine/docs/equations.md and ALG-CALORICS in engine/docs/science/algorithms.md. Coefficients do not define their own validity interval.
    """
    coefficients: list[float]
    """Coefficients a[k] in Cp = sum(a[k]*(T-Tref)**k), ascending k; units J/(mol K**(k+1)),
    nonempty.
    """
    reference_temperature: float
    """Enthalpy and entropy anchor temperature in K, within the owning interval."""
    reference_enthalpy: float
    """Ideal-gas molar enthalpy at the anchor, J/mol."""
    reference_entropy: float
    """Ideal-gas molar entropy at the anchor and source pressure, J/(mol K)."""
    def __init__(self, coefficients: Sequence[float], reference_temperature: float,
                 reference_enthalpy: float, reference_entropy: float) -> None:
        """Initialize IdealPolynomial from the declared inputs.

        Parameters
        ----------
        coefficients : Sequence[float]
            Initial value of IdealPolynomial.coefficients; see its attribute description for units and ordering.
        reference_temperature : float
            Initial value of IdealPolynomial.reference_temperature; see its attribute description for units and ordering.
        reference_enthalpy : float
            Initial value of IdealPolynomial.reference_enthalpy; see its attribute description for units and ordering.
        reference_entropy : float
            Initial value of IdealPolynomial.reference_entropy; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class IdealShomate:
    """NIST Shomate ideal-gas molar thermochemistry.

    Notes
    -----
    See cp_ideal_shomate in engine/docs/equations.md; preserve the source pressure, integration constants and validity interval.
    """
    coefficients: list[float]
    """Eight unmodified NIST values [A, B, C, D, E, F, G, H] with t=T/1000 K: Cp and S use J/(mol
    K), the enthalpy increment uses kJ/mol before native conversion to J/mol.
    """
    formation_enthalpy: float
    """Source formation enthalpy in J/mol added to the Shomate enthalpy increment; not inferred
    from the coefficient H.
    """
    def __init__(self, coefficients: Sequence[float], formation_enthalpy: float) -> None:
        """Initialize IdealShomate from the declared inputs.

        Parameters
        ----------
        coefficients : Sequence[float]
            Initial value of IdealShomate.coefficients; see its attribute description for units and ordering.
        formation_enthalpy : float
            Initial value of IdealShomate.formation_enthalpy; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class IdealNASA7:
    """NASA7 ideal-gas molar thermochemistry.

    Notes
    -----
    See cp_ideal_nasa7 in engine/docs/equations.md. Coefficients must be enclosed in an explicit IdealInterval.
    """
    coefficients: list[float]
    """Seven coefficients [a0, a1, a2, a3, a4, a5, a6] in source order, with T
    in K. Cp/R uses powers T**0 through T**4; a5/T and a6 are the h/(RT)
    and s/R integration terms. Units are dimensionless for a0 and a6,
    K**(-k) for a[k] with k=1 through 4, and K for a5. Native evaluation
    uses GAS_CONSTANT_J_PER_MOL_K; prepare_thermochemistry accepts NASA9 only."""
    def __init__(self, coefficients: Sequence[float]) -> None:
        """Initialize IdealNASA7 from the declared inputs.

        Parameters
        ----------
        coefficients : Sequence[float]
            Initial value of IdealNASA7.coefficients; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class IdealNASA9:
    """NASA9 ideal-gas molar thermochemistry.

    Notes
    -----
    See cp_ideal_nasa9 in engine/docs/equations.md. Coefficients must be enclosed in an explicit IdealInterval.
    """
    coefficients: list[float]
    """9 values in source order a0 through a8; temperature is in K. Heat capacity divided by R uses
    powers T**-2 through T**4; a7/T and a8 are the enthalpy and entropy integration terms.
    Coefficient units follow these powers; retain the source gas constant until
    prepare_thermochemistry.
    """
    def __init__(self, coefficients: Sequence[float]) -> None:
        """Initialize IdealNASA9 from the declared inputs.

        Parameters
        ----------
        coefficients : Sequence[float]
            Initial value of IdealNASA9.coefficients; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class IdealDIPPR107:
    """Ideal-gas DIPPR equation 107 molar heat-capacity correlation.

    Notes
    -----
    See cp_ideal_dippr_107 in engine/docs/equations.md and ALG-CALORICS in engine/docs/science/algorithms.md. Coefficients do not define their own validity interval.
    """
    coefficients: list[float]
    """Five values [A, B, C, D, E] for Cp = A+B*((C/T)/sinh(C/T))**2+D*((E/T)/cosh(E/T))**2; A/B/D
    in J/(mol K), C/E in K.
    """
    reference_temperature: float
    """Enthalpy and entropy anchor temperature in K, within the owning interval."""
    reference_enthalpy: float
    """Ideal-gas molar enthalpy at the anchor, J/mol."""
    reference_entropy: float
    """Ideal-gas molar entropy at the anchor and source pressure, J/(mol K)."""
    def __init__(self, coefficients: Sequence[float], reference_temperature: float,
                 reference_enthalpy: float, reference_entropy: float) -> None:
        """Initialize IdealDIPPR107 from the declared inputs.

        Parameters
        ----------
        coefficients : Sequence[float]
            Initial value of IdealDIPPR107.coefficients; see its attribute description for units and ordering.
        reference_temperature : float
            Initial value of IdealDIPPR107.reference_temperature; see its attribute description for units and ordering.
        reference_enthalpy : float
            Initial value of IdealDIPPR107.reference_enthalpy; see its attribute description for units and ordering.
        reference_entropy : float
            Initial value of IdealDIPPR107.reference_entropy; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

IdealForm = Union[IdealPolynomial, IdealShomate, IdealNASA7, IdealNASA9, IdealDIPPR107]
"""Typing union of the admitted ideal-gas correlation forms; not a native runtime class."""

class IdealCorrelation:
    """One ideal-gas functional form with its source pressure.

    Notes
    -----
    Native thermochemistry transfers entropy to ThermochemistryRecord.standard_pressure; it does not change source enthalpy anchors.
    """
    form: IdealForm
    """Selected ideal-gas correlation; its coefficients retain their source convention."""
    source_pressure: float
    """Positive pressure of the source entropy standard in Pa."""
    def __init__(self, form: IdealForm, source_pressure: float) -> None:
        """Initialize IdealCorrelation from the declared inputs.

        Parameters
        ----------
        form : IdealForm
            Initial value of IdealCorrelation.form; see its attribute description for units and ordering.
        source_pressure : float
            Initial value of IdealCorrelation.source_pressure; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class IdealInterval:
    """Explicit temperature domain and endpoint ownership for one correlation.

    Notes
    -----
    Adjacent intervals must meet with exactly one owning endpoint. Gaps and overlaps are rejected. Endpoint values do not imply derivatives through an unqualified join.
    """
    lower: float
    """Positive lower temperature in K."""
    upper: float
    """Upper temperature in K, greater than lower."""
    lower_inclusive: bool
    """Whether the lower endpoint belongs to this interval."""
    upper_inclusive: bool
    """Whether the upper endpoint belongs to this interval."""
    correlation: IdealCorrelation
    """Ideal-gas correlation valid on this interval."""
    def __init__(self, lower: float, upper: float, lower_inclusive: bool,
                 upper_inclusive: bool, correlation: IdealCorrelation) -> None:
        """Initialize IdealInterval from the declared inputs.

        Parameters
        ----------
        lower : float
            Initial value of IdealInterval.lower; see its attribute description for units and ordering.
        upper : float
            Initial value of IdealInterval.upper; see its attribute description for units and ordering.
        lower_inclusive : bool
            Initial value of IdealInterval.lower_inclusive; see its attribute description for units and ordering.
        upper_inclusive : bool
            Initial value of IdealInterval.upper_inclusive; see its attribute description for units and ordering.
        correlation : IdealCorrelation
            Initial value of IdealInterval.correlation; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class ThermochemistryRecord:
    """Species ideal-gas caloric records on one pressure standard.

    Notes
    -----
    Ideal and total calorics consume this record. Missing or out-of-range records do not remove independently available residual quantities.
    """
    standard_pressure: float
    """Positive common ideal-gas standard pressure in Pa."""
    component_intervals: list[list[IdealInterval]]
    """Ragged array of n interval lists in component_ids order; each retained species list states
    its own temperature intervals and endpoint ownership.
    """
    def __init__(self, standard_pressure: float,
                 component_intervals: Sequence[Sequence[IdealInterval]]) -> None:
        """Initialize ThermochemistryRecord from the declared inputs.

        Parameters
        ----------
        standard_pressure : float
            Initial value of ThermochemistryRecord.standard_pressure; see its attribute description for units and ordering.
        component_intervals : Sequence[Sequence[IdealInterval]]
            Initial value of ThermochemistryRecord.component_intervals; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class ReactionLogPolynomial:
    """Source equilibrium law ln K = a + b/T + c ln(T/Tref) + d T + e/T² + f T².

    Notes
    -----
    See reaction_constant_ln_polynomial in engine/docs/equations.md. Without an explicit ReactionReference, equilibrium requires EOS_STANDARD_STATE_ID. With a source reference and the problem neutral_reference, native code converts the source standard. An isothermal K alone supplies no caloric derivative.
    """
    a: float
    """Natural-log constant coefficient, dimensionless."""
    b: float
    """Reciprocal-temperature coefficient, K."""
    c: float
    """Log-temperature coefficient, dimensionless."""
    d: float
    """Linear-temperature coefficient, 1/K."""
    reference_temperature: float
    """Positive logarithm reference temperature Tref in K."""
    temperature_min: float
    """Inclusive lower validity temperature in K."""
    temperature_max: float
    """Inclusive upper validity temperature in K."""
    standard_state_id: str
    """Identity of the source activity standard; K is dimensionless on this standard."""
    e: float
    """Inverse-square-temperature coefficient, K²; defaults to zero."""
    f: float
    """Square-temperature coefficient, 1/K²; defaults to zero."""
    def __init__(self, a: float, b: float, c: float, d: float, reference_temperature: float,
                 temperature_min: float, temperature_max: float, standard_state_id: str,
                 e: float = 0.0, f: float = 0.0) -> None:
        """Initialize ReactionLogPolynomial from the declared inputs.

        Parameters
        ----------
        a : float
            Initial value of ReactionLogPolynomial.a; see its attribute description for units and ordering.
        b : float
            Initial value of ReactionLogPolynomial.b; see its attribute description for units and ordering.
        c : float
            Initial value of ReactionLogPolynomial.c; see its attribute description for units and ordering.
        d : float
            Initial value of ReactionLogPolynomial.d; see its attribute description for units and ordering.
        reference_temperature : float
            Initial value of ReactionLogPolynomial.reference_temperature; see its attribute description for units and ordering.
        temperature_min : float
            Initial value of ReactionLogPolynomial.temperature_min; see its attribute description for units and ordering.
        temperature_max : float
            Initial value of ReactionLogPolynomial.temperature_max; see its attribute description for units and ordering.
        standard_state_id : str
            Initial value of ReactionLogPolynomial.standard_state_id; see its attribute description for units and ordering.
        e : float
            Inverse-square-temperature coefficient in K².
        f : float
            Square-temperature coefficient in 1/K².

        Returns
        -------
        None
            None.
        """
        ...

class ParameterFamily:
    """Native independent parameter family and coordinate units.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.

    Indices select native storage slots, not component tuples. A symmetric tied parameter may identify multiple flattened slots. Native parameter construction owns the selected model and admitted slots.
    """
    BornCDielectric: ParameterFamily
    """Born dielectric coefficient, dimensionless; scalar slot."""
    IonicRegionRelativePermittivity: ParameterFamily
    """Ionic-region relative permittivity, dimensionless; scalar slot."""
    BornDiameter: ParameterFamily
    """Born diameter in angstrom, component slots."""
    DebyeHuckelDiameter: ParameterFamily
    """Debye-Huckel closest-approach diameter in angstrom, component slots."""
    PackingDiameter: ParameterFamily
    """Ionic packing diameter in angstrom, component slots."""
    IonSuppressionCoefficient: ParameterFamily
    """Ion suppression coefficient, dimensionless, in charged-component slots."""
    SegmentCount: ParameterFamily
    """Segment count per molecule, dimensionless, component slots."""
    SigmaConstant: ParameterFamily
    """Constant segment diameter in angstrom, component slots."""
    DispersionEnergyOverK: ParameterFamily
    """Dispersion energy divided by Boltzmann's constant, K, component slots."""
    Kij: ParameterFamily
    """Unlike dispersion interaction coefficient, dimensionless, flattened component-pair slots."""
    Lij: ParameterFamily
    """Unlike size interaction coefficient, dimensionless, flattened component-pair slots."""
    AssociationEnergyOverK: ParameterFamily
    """Association energy divided by Boltzmann's constant, K, flattened site-pair slots."""
    AssociationVolume: ParameterFamily
    """Association bonding-volume coefficient, dimensionless, flattened site-pair slots."""
    KijReciprocalTemperatureSlope: ParameterFamily
    """b_ij in k_ij(T)=k_ij+b_ij*(1/T-1/Tref), K, flattened component-pair slots."""
    DielectricVolume: ParameterFamily
    """Solvent dielectric molar-volume coefficient in m^3/mol, component slots."""
    DielectricTemperature: ParameterFamily
    """Solvent dielectric temperature coefficient in K, component slots."""
    TemperatureLawCoefficient: ParameterFamily
    """Coefficient slot of a declared temperature law; units depend on the declared form and
    polynomial power.
    """
    CationPairMixing: ParameterFamily
    """Declared cation-pair mixing coefficient in kg/mol; scalar slot."""
    CationTripleMixing: ParameterFamily
    """Declared cation-triple mixing coefficient in kg^2/mol^2; scalar slot."""

class ActiveParameter:
    """One independent parameter value assigned to native storage slots."""
    family: ParameterFamily
    """Family that determines physical meaning and units; see ParameterFamily."""
    indices: list[int]
    """Native zero-based storage slots sharing this value, length at least one; pair slots are
    flattened in row-major order.
    """
    def __init__(self, family: ParameterFamily, indices: Sequence[int]) -> None:
        """Initialize ActiveParameter from the declared inputs.

        Parameters
        ----------
        family : ParameterFamily
            Initial value of ActiveParameter.family; see its attribute description for units and ordering.
        indices : Sequence[int]
            Initial value of ActiveParameter.indices; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class DerivativeOrder:
    """Requested derivative information from the native energy owner.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    Value: DerivativeOrder
    """Central value only; derivative arrays are empty."""
    Gradient: DerivativeOrder
    """Value and first partial derivatives."""
    PressureSlope: DerivativeOrder
    """Reduced energy information sufficient for a pressure slope; not a full Hessian and not
    admitted by neutral_property.
    """
    Hessian: DerivativeOrder
    """Value, gradient and second partial derivatives."""

class NeutralProperty:
    """Observable for the restricted single-neutral-component derivative reader.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    ResidualHelmholtz: NeutralProperty
    """Residual Helmholtz energy divided by nRT, dimensionless."""
    Pressure: NeutralProperty
    """Pressure in Pa."""

class NeutralPropertyDerivatives:
    """Single-component property partial derivatives on (T, rho)."""
    value: float
    """Requested observable: dimensionless residual Helmholtz energy or pressure in Pa."""
    gradient: list[float]
    """Length 2 in [T, rho] order when requested, otherwise empty; observable units divided by K or
    mol/m^3.
    """
    hessian: list[float]
    """Flattened row-major (2, 2) array in [T, rho] order when requested, otherwise empty; units
    divide the observable unit by both coordinate units.
    """

class PropertyObservable:
    """Homogeneous-state observable and its molar units.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.

    Heat capacities hold composition fixed. Ideal and total properties require admitted thermochemistry. Actions retain these observable units, because path variables are dimensionless.
    """
    IdealEnthalpy: PropertyObservable
    """Ideal molar enthalpy, J/mol."""
    IdealInternalEnergy: PropertyObservable
    """Ideal molar internal energy, J/mol."""
    IdealEntropy: PropertyObservable
    """Ideal molar entropy, J/(mol K)."""
    IdealIsobaricHeatCapacity: PropertyObservable
    """Ideal molar isobaric heat capacity, J/(mol K)."""
    IdealIsochoricHeatCapacity: PropertyObservable
    """Ideal molar isochoric heat capacity, J/(mol K)."""
    ResidualEnthalpy: PropertyObservable
    """Residual molar enthalpy, J/mol."""
    ResidualInternalEnergy: PropertyObservable
    """Residual molar internal energy, J/mol."""
    ResidualEntropy: PropertyObservable
    """Residual molar entropy, J/(mol K)."""
    ResidualIsobaricHeatCapacity: PropertyObservable
    """Residual molar isobaric heat capacity, J/(mol K)."""
    ResidualIsochoricHeatCapacity: PropertyObservable
    """Residual molar isochoric heat capacity, J/(mol K)."""
    TotalEnthalpy: PropertyObservable
    """Total molar enthalpy, J/mol."""
    TotalInternalEnergy: PropertyObservable
    """Total molar internal energy, J/mol."""
    TotalEntropy: PropertyObservable
    """Total molar entropy, J/(mol K)."""
    TotalIsobaricHeatCapacity: PropertyObservable
    """Total molar isobaric heat capacity, J/(mol K)."""
    TotalIsochoricHeatCapacity: PropertyObservable
    """Total molar isochoric heat capacity, J/(mol K)."""
    LogFugacityCoefficient: PropertyObservable
    """Natural logarithm of the component fugacity coefficient, dimensionless; pressure-chart only,
    one returned result per component.
    """

class PropertyActionStatus:
    """Availability or refusal cause for a homogeneous-state property action.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    Available: PropertyActionStatus
    """Requested quantity is available on the admitted branch."""
    MissingThermochemistry: PropertyActionStatus
    """Required species thermochemistry is absent."""
    OutsideInterval: PropertyActionStatus
    """Temperature lies outside a required source interval."""
    ThermalBoundary: PropertyActionStatus
    """Requested thermal derivative crosses an unqualified source-interval endpoint."""
    SupportChange: PropertyActionStatus
    """Requested derivative leaves the admitted composition support or meets a singular insertion."""
    IonicInsertion: PropertyActionStatus
    """Requested charged insertion is unavailable on an ion-free support."""
    AssociationSingular: PropertyActionStatus
    """The implicit association-site Jacobian is singular."""
    AssociationIllConditioned: PropertyActionStatus
    """Association-site conditioning fails the declared criterion."""
    AssociationResidual: PropertyActionStatus
    """The association residual does not meet its admission criterion."""
    DensityUnavailable: PropertyActionStatus
    """A required pressure-closed density is unavailable."""
    PressureSingular: PropertyActionStatus
    """The pressure-density slope does not admit the requested implicit action."""
    CentralPropertyUnavailable: PropertyActionStatus
    """The central observable is unavailable."""
    Nonfinite: PropertyActionStatus
    """A required value or derivative is nonfinite."""

class ParameterActionDirection:
    """Physical increment of one native parameter direction."""
    family: ParameterFamily
    """Parameter family and its physical unit; see ParameterFamily."""
    indices: list[int]
    """Native storage slots receiving the increment, in the same convention as
    ActiveParameter.indices.
    """
    value: float
    """Physical increment in the selected family's native units, not a dimensionless derivative
    normalization.
    """
    def __init__(self, family: ParameterFamily, indices: Sequence[int], value: float) -> None:
        """Initialize ParameterActionDirection from the declared inputs.

        Parameters
        ----------
        family : ParameterFamily
            Initial value of ParameterActionDirection.family; see its attribute description for units and ordering.
        indices : Sequence[int]
            Initial value of ParameterActionDirection.indices; see its attribute description for units and ordering.
        value : float
            Initial value of ParameterActionDirection.value; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class DensityActionDirection:
    """Direction on the homogeneous density chart.

    Notes
    -----
    The path variable is dimensionless. Mixed actions contract all supplied directions; divide by each physical increment to form a partial derivative.
    """
    temperature: float
    """Temperature increment in K."""
    molar_density: float
    """Molar density increment in mol/m^3."""
    composition: list[float]
    """Dimensionless mole-fraction increments, length n in component_ids order.
    Supply n zeros for no composition motion; an empty vector is rejected by
    the native action readers. Increments must preserve normalization and charge balance."""
    parameter: list[ParameterActionDirection]
    """Physical parameter increments, in declaration order; empty means no parameter motion."""
    def __init__(self, temperature: float = 0.0, molar_density: float = 0.0,
                 composition: Sequence[float] = (),
                 parameter: Sequence[ParameterActionDirection] = ()) -> None:
        """Initialize DensityActionDirection from the declared inputs.

        Parameters
        ----------
        temperature : float
            Initial value of DensityActionDirection.temperature; see its attribute description for units and ordering.
        molar_density : float
            Initial value of DensityActionDirection.molar_density; see its attribute description for units and ordering.
        composition : Sequence[float]
            Initial value of DensityActionDirection.composition; see its attribute description for units and ordering.
        parameter : Sequence[ParameterActionDirection]
            Initial value of DensityActionDirection.parameter; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class PressureActionDirection:
    """Direction on the homogeneous pressure chart.

    Notes
    -----
    The path variable is dimensionless. Mixed actions contract all supplied directions; divide by each physical increment to form a partial derivative.
    """
    temperature: float
    """Temperature increment in K."""
    pressure: float
    """Pressure increment in Pa."""
    composition: list[float]
    """Dimensionless mole-fraction increments, length n in component_ids order.
    Supply n zeros for no composition motion; an empty vector is rejected by
    the native action readers. Increments must preserve normalization and charge balance."""
    parameter: list[ParameterActionDirection]
    """Physical parameter increments, in declaration order; empty means no parameter motion."""
    def __init__(self, temperature: float = 0.0, pressure: float = 0.0,
                 composition: Sequence[float] = (),
                 parameter: Sequence[ParameterActionDirection] = ()) -> None:
        """Initialize PressureActionDirection from the declared inputs.

        Parameters
        ----------
        temperature : float
            Initial value of PressureActionDirection.temperature; see its attribute description for units and ordering.
        pressure : float
            Initial value of PressureActionDirection.pressure; see its attribute description for units and ordering.
        composition : Sequence[float]
            Initial value of PressureActionDirection.composition; see its attribute description for units and ordering.
        parameter : Sequence[ParameterActionDirection]
            Initial value of PressureActionDirection.parameter; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class PropertyActionResult:
    """Central property and one mixed directional contraction."""
    observable: PropertyObservable
    """Observable identifying the value and action units; see PropertyObservable."""
    component: Optional[int]
    """Component index for log fugacity, or None for a mixture caloric property."""
    value: Optional[float]
    """Central molar property, or None when unavailable; units follow observable."""
    action: Optional[float]
    """Mixed contraction over all supplied directions, or None with the status cause; retains the
    observable unit.
    """
    status: PropertyActionStatus
    """Availability or specific refusal cause; an available central value can survive an
    unavailable action.
    """
    energy_order: int
    """Highest Helmholtz derivative order required for this result, a count."""

class DensityBranch:
    """Requested or established side of the pressure-density structure.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.

    A mechanically stable root alone does not establish global phase stability. See ALG-DENSITY in engine/docs/science/algorithms.md.
    """
    Unspecified: DensityBranch
    """No liquid/vapor branch classification established; this is not a third physical phase."""
    Vapor: DensityBranch
    """Requested vapor-side root."""
    Liquid: DensityBranch
    """Requested liquid-side root."""
    def __int__(self) -> int:
        """Read the native enum code.

        Returns
        -------
        int
            Integer enum value; see the module enum-property convention.
        """
        ...

class ReactionThermochemistry:
    """Species-contracted reaction standards at the requested temperature.

    Values are defined only with Available status. Gibbs energy and enthalpy use
    J/mol of reaction, heat capacity J/(mol K), and Gibbs derivatives J/(mol K**k).
    Density standard is 1 mol/m³; pressure standard is the declared record pressure.
    ln K is dimensionless; K is absent if exponentiation is not positive and finite.
    Unavailable results retain the component, temperature and cause in ``cause``.
    """
    gibbs: float
    enthalpy: float
    heat_capacity: float
    pressure_gibbs: float
    ln_k_pressure: float
    ln_k_density: float
    k_pressure: float | None
    k_density: float | None
    gibbs_derivatives: list[float]
    status: ThermochemistryStatus
    cause: str

def reaction_thermochemistry_density(record: ThermochemistryRecord, stoichiometry: Sequence[float], temperature: float) -> ReactionThermochemistry:
    """Contract participant records and convert pressure to density standards natively.

    Parameters
    ----------
    record : ThermochemistryRecord
        Thermal records in native component order, with pressure standard in Pa.
    stoichiometry : Sequence[float]
        Products positive, in the same component order.
    temperature : float
        Requested temperature in K.

    Returns
    -------
    ReactionThermochemistry
        Available reaction values or the retained source-unavailability cause.
    """
    ...

class ThermochemistryStatus:
    """Availability of ideal-gas thermochemistry at the state temperature.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    MissingRecord: ThermochemistryStatus
    """No required thermochemistry record was supplied."""
    Available: ThermochemistryStatus
    """The required ideal-gas calorics are available."""
    OutsideInterval: ThermochemistryStatus
    """A required species temperature interval excludes the state."""
    Nonfinite: ThermochemistryStatus
    """A required caloric evaluation is nonfinite."""

class PropertyActionDiagnostics:
    """Numerical work and conditioning for one property-action batch."""
    site_rank: int
    """Rank of the association-site Jacobian, a count."""
    site_factorizations: int
    """Number of association Jacobian factorizations."""
    site_solves: int
    """Number of association linear solves."""
    site_pivot_ratio: float
    """Dimensionless largest/smallest absolute association pivot ratio when the
    smallest pivot is positive. Zero when not computed or when the smallest pivot
    is zero; inspect the action status before interpreting conditioning."""
    association_residual: float
    """Dimensionless association mass-action residual magnitude."""
    pressure_density_derivative: float
    """Pressure-chart central dP/drho at fixed T and composition, Pa m^3/mol.
    Zero when no pressure root was supplied; the density chart does not compute
    this diagnostic."""
    pressure_residual: float
    """Pressure-chart EOS pressure minus prescribed pressure, Pa. Zero when no
    pressure root was supplied; this default does not establish pressure closure."""
    stability_determinant: float
    """Pressure-chart dimensionless slope (dP/drho)/(RT); mechanical regularity
    evidence rather than global stability. Zero when no pressure root was
    supplied; this default does not identify a spinodal."""
    pressure_factorizations: int
    """Number of implicit pressure factorizations."""
    pressure_solves: int
    """Number of implicit pressure linear solves."""
    thermal_evaluations: int
    """Number of thermal-correlation evaluations."""
    recordings: int
    """Number of automatic-differentiation recordings."""
    directional_paths: int
    """Number of polarized directional paths."""
    forward_sweeps: int
    """Number of forward derivative sweeps."""
    branch: DensityBranch
    """Established density-branch label, or Unspecified."""
    stable_root_count: Optional[int]
    """Number of mechanically stable roots established by the search, or None when not established."""

    first_numerical_failure: Optional[Convergence]
    """Earliest observed numerical refusal, separate from final convergence."""

class PropertyActionBatch:
    """Property actions sharing one central state and implicit work."""
    solver_options: list[dict[str, Any]]
    """Effective numerical settings as named dictionaries with values, provenance and units where
    declared.
    """
    results: list[PropertyActionResult]
    """Results in observable order; log-fugacity requests expand into component_ids order."""
    maximum_energy_order: int
    """Reported Helmholtz derivative-order requirement for the batch. It can be
    populated when evaluation is refused and does not establish completed
    derivative work."""
    diagnostics: PropertyActionDiagnostics
    """Shared conditioning, work and pressure-closure diagnostics."""
    molar_density: Optional[float]
    """Pressure-closed central density in mol/m^3, or None when unavailable or not supplied by this
    chart.
    """
    density_action: list[float]
    """First density action per direction, mol/m^3; empty when density actions are unavailable."""
    density_mixed_action: Optional[float]
    """Mixed density action over two or three directions, mol/m^3, or None when not requested or
    unavailable.
    """
    density_status: PropertyActionStatus
    """Availability or refusal cause for the implicit density actions."""
    density_message: str
    """Native density refusal explanation, possibly empty."""
    density_convergence: Optional[Convergence]
    """Density-closure failure record when a ConvergenceFailure was caught. None
    for successful pressure closure and for failures without that record."""

class State:
    """One homogeneous state with SI molar quantities.

    Notes
    -----
    Component arrays use the declared component_ids order. Heat capacities and derived response functions hold composition fixed; they are not reacting equilibrium heat capacities. See full-eos.md, ALG-CALORICS and ALG-ACTIONS in engine/docs/science/.
    """
    solver_options: list[dict[str, Any]]
    """Effective named solver settings with provenance."""
    hard_chain: float
    """Hard-chain contribution to A_res/(nRT), dimensionless."""
    dispersion: float
    """Dispersion contribution to A_res/(nRT), dimensionless."""
    association: float
    """Association contribution to A_res/(nRT), dimensionless."""
    debye_huckel: float
    """Debye-Huckel contribution to A_res/(nRT), dimensionless."""
    born: float
    """Born contribution to A_res/(nRT), dimensionless."""
    cation_mixing: float
    """Cation mixing contribution to A_res/(nRT), dimensionless."""
    dipolar: float
    """Declared polar contribution to A_res/(nRT), dimensionless."""
    bulk_relative_permittivity: float
    """Bulk relative permittivity, dimensionless."""
    residual_helmholtz: float
    """Residual molar Helmholtz energy divided by RT, dimensionless."""
    compressibility_factor: float
    """Compressibility factor P/(rho RT), dimensionless."""
    pressure: float
    """Calculated EOS pressure in Pa; a pressure-closed state retains its closure residual."""
    residual_chemical_potential_over_rt: list[Optional[float]]
    """Length n in component_ids order: mu_res/(RT) against ideal gas at the same T and rho; an
    entry is None for unavailable charged insertion.
    """
    log_fugacity_coefficient: Optional[list[Optional[float]]]
    """Dimensionless component ln(phi_i), length n in component_ids order. Direct-density
    states subtract ln(compressibility_factor) from mu_res_i/(RT); pressure-closed
    states subtract ln(P_target/(rho RT)) using the prescribed pressure. The list
    is None unless the EOS compressibility factor is positive and finite; unavailable
    component chemical potentials give None entries."""
    component_property_status: list[PropertyActionStatus]
    """Length n, component_ids order; explains component-property availability, including
    IonicInsertion.
    """
    log_fugacity_parameter_derivative: list[Optional[float]]
    """Flattened row-major (n, a) first partials of ln(phi_i) at fixed T/P/x, in
    inverse native parameter units. Populated only when the fugacity list exists
    and dP/drho is finite and positive; otherwise empty. Within a populated array,
    unavailable component insertions have None entries."""
    residual_entropy: float
    """Residual molar entropy in J/(mol K). Residual means minus ideal gas at the same T and rho."""
    residual_internal_energy: float
    """Residual molar internal energy in J/mol. Residual means minus ideal gas at the same T and
    rho.
    """
    residual_enthalpy: float
    """Residual molar enthalpy in J/mol. Residual means minus ideal gas at the same T and rho."""
    residual_isochoric_heat_capacity: float
    """Residual molar isochoric heat capacity in J/(mol K). Residual means minus ideal gas at the
    same T and rho.
    """
    residual_isobaric_heat_capacity: Optional[float]
    """Residual molar isobaric heat capacity in J/(mol K). Residual means minus ideal gas at the
    same T and rho. None when the pressure slope does not admit the Cp calculation.
    """
    ideal_enthalpy: Optional[float]
    """Ideal molar enthalpy in J/mol. None when required ideal-gas thermochemistry is unavailable."""
    ideal_internal_energy: Optional[float]
    """Ideal molar internal energy in J/mol. None when required ideal-gas thermochemistry is
    unavailable.
    """
    ideal_entropy: Optional[float]
    """Ideal molar entropy in J/(mol K). None when required ideal-gas thermochemistry is
    unavailable.
    """
    ideal_isobaric_heat_capacity: Optional[float]
    """Ideal molar isobaric heat capacity in J/(mol K). None when required ideal-gas
    thermochemistry is unavailable.
    """
    ideal_isochoric_heat_capacity: Optional[float]
    """Ideal molar isochoric heat capacity in J/(mol K). None when required ideal-gas
    thermochemistry is unavailable.
    """
    total_enthalpy: Optional[float]
    """Total molar enthalpy in J/mol. None when required ideal-gas thermochemistry is unavailable."""
    total_internal_energy: Optional[float]
    """Total molar internal energy in J/mol. None when required ideal-gas thermochemistry is
    unavailable.
    """
    total_entropy: Optional[float]
    """Total molar entropy in J/(mol K). None when required ideal-gas thermochemistry is
    unavailable.
    """
    total_isobaric_heat_capacity: Optional[float]
    """Total molar isobaric heat capacity in J/(mol K). None when required ideal-gas
    thermochemistry is unavailable.
    """
    total_isochoric_heat_capacity: Optional[float]
    """Total molar isochoric heat capacity in J/(mol K). None when required ideal-gas
    thermochemistry is unavailable.
    """
    molar_density: float
    """Molar density in mol/m^3. Multiply by mixture molar mass in kg/mol to obtain kg/m^3."""
    pressure_density_derivative: float
    """dP/drho at fixed temperature and composition, Pa m^3/mol."""
    site_fractions: list[float]
    """Length s of dimensionless unbonded fractions in declared site order."""
    density_convergence: Optional[Convergence]
    """Pressure-root convergence evidence, or None for a direct density state."""
    thermochemistry_status: ThermochemistryStatus
    """Availability of the ideal-gas thermochemistry used by ideal and total quantities."""
    density_branch: DensityBranch
    """Established density-root branch or Unspecified for a state without a classified root."""
    requested_density_branch: DensityBranch
    """Branch requested by pressure closure; Unspecified for direct density states."""
    returned_density_branch: DensityBranch
    """Side from existing finite-search evidence or the locally admitted anchor label.
    Unspecified when no side is established. A wrongly labelled external anchor can conceal
    a mismatch; this attribute does not independently validate its physical phase.
    """
    density_branch_mismatch: bool
    """True only for explicit, opposite requested and returned branches. False means no
    established mismatch, including an unknown side or an Unspecified request.
    """
    stable_root_count: Optional[int]
    """Number of mechanically stable roots established by a pressure search, or None when not
    established.
    """
    isothermal_compressibility: Optional[float]
    """Fixed-composition isothermal compressibility in 1/Pa; None unless derived properties were
    requested and the pressure slope is positive.
    """
    isobaric_expansivity: Optional[float]
    """Fixed-composition isobaric expansivity in 1/K; None unless requested and the pressure slope
    is positive.
    """
    heat_capacity_ratio: Optional[float]
    """Total Cp/Cv, dimensionless; None unless requested and the required total heat capacities and
    positive pressure slope are available.
    """
    isentropic_compressibility: Optional[float]
    """Fixed-composition isentropic compressibility in 1/Pa; None unless requested and total heat
    capacities and positive pressure slope are available.
    """
    speed_of_sound: Optional[float]
    """Thermodynamic sound speed at frozen composition in m/s; None unless requested, the pressure
    slope and calorics are admitted, and occupied components have positive molar masses.
    """
    joule_thomson_coefficient: Optional[float]
    """Fixed-composition Joule-Thomson coefficient in K/Pa; None unless requested and the required
    total Cp and positive pressure slope are available.
    """
    partial_molar_volume: Optional[list[Optional[float]]]
    """Partial molar volumes in m^3/mol, in component_ids order. None when derived
    properties were not requested; otherwise a length-n list. Entries are None
    without a positive finite pressure slope, for unavailable insertions, or for nonfinite results."""
    partial_molar_enthalpy: Optional[list[Optional[float]]]
    """Partial molar enthalpies in J/mol, in component_ids order. None when derived
    properties were not requested; otherwise a length-n list. Entries are None
    without a positive finite pressure slope or admitted thermochemistry, for unavailable
    insertions, or for nonfinite results."""
    partial_molar_entropy: Optional[list[Optional[float]]]
    """Partial molar entropies in J/(mol K), in component_ids order. None when derived
    properties were not requested; otherwise a length-n list. Entries are None
    without a positive finite pressure slope or admitted thermochemistry, for zero-fraction
    species or unavailable insertions, or for nonfinite results."""

    first_numerical_failure: Optional[Convergence]
    """Earliest observed numerical refusal, separate from final convergence."""

class DensityDiagnosticPoint:
    """One sampled or refined point on a fixed-composition isotherm."""
    molar_density: float
    """Point density in mol/m^3."""
    packing_fraction: float
    """Dimensionless packing fraction; may be nonfinite for an unevaluated or failed point."""
    pressure: float
    """EOS pressure in Pa; may be nonfinite after evaluation failure."""
    pressure_density_derivative: float
    """dP/drho at fixed T/x, Pa m^3/mol; may be nonfinite after evaluation failure."""
    lower: float
    """Lower density bracket endpoint, mol/m^3."""
    upper: float
    """Upper density bracket endpoint, mol/m^3."""
    residual: float
    """Producer's root residual: pressure in Pa for pressure crossings or dimensionless
    (dP/drho)/(RT) for spinodals; may be nonfinite after refusal.
    """
    mechanical_sign: int
    """Sign of dP/drho for ordinary evaluated samples and pressure roots: +1 for
    positive slope and -1 for negative slope. Spinodal records force this value
    to zero even when their retained slope is not exactly zero; unevaluated
    records also retain zero."""
    linked: bool
    """Whether the diagnostic linked this point through its scan/refinement logic; not a global
    phase-stability statement.
    """
    declared: bool
    """False for extra samples introduced by the adaptive slope search; True for
    records retaining the declared-evidence default, including failed declared
    nodes and unresolved brackets. This marker does not indicate successful resolution."""
    cause: str
    """Native failure or refinement explanation, possibly empty."""
    stop: str
    """Native diagnostic termination name."""

class DensityDiagnosticBracket:
    """A declared pressure-crossing density bracket."""
    lower: float
    """Lower molar density in mol/m^3."""
    upper: float
    """Upper molar density in mol/m^3."""
    lower_residual: float
    """Pressure minus target pressure at lower, Pa."""
    upper_residual: float
    """Pressure minus target pressure at upper, Pa."""

class DensityDiagnostics:
    """Finite pressure-root and mechanical-loop search evidence."""
    solver_options: list[dict[str, Any]]
    """Effective named numerical settings and provenance."""
    samples: list[DensityDiagnosticPoint]
    """Sampled points in the retained evaluation order."""
    pressure_roots: list[DensityDiagnosticPoint]
    """Refined pressure crossings; each point retains its bracket and cause."""
    spinodals: list[DensityDiagnosticPoint]
    """Refined zeros of dP/drho, not global coexistence states."""
    declared_brackets: list[DensityDiagnosticBracket]
    """Pressure-crossing brackets declared by the finite search."""
    unresolved: list[DensityDiagnosticPoint]
    """Points or intervals not resolved within the declared work and domain."""
    evaluation_order: DerivativeOrder
    """Derivative order used to test finite evaluability."""
    finite_domain: Optional[DensityDiagnosticPoint]
    """Refined finite-domain boundary point, or None if no such search was requested or retained.
    Only its density bounds, stopping cause and availability cause are defined.
    """
    geometric_cap: float
    """Packing-derived upper density bound in mol/m^3; not proof of finite evaluability below it."""
    declared_crossings: int
    """Number of strict sign-changing pressure brackets between connected declared
    nodes; equals len(declared_brackets). Exact-zero samples and additional
    adaptive crossings do not increase this count."""
    evaluations: int
    """Number of diagnostic EOS evaluations."""
    stop: str
    """Diagnostic termination explanation; an incomplete search does not establish uniqueness."""

    first_numerical_failure: Optional[Convergence]
    """Earliest observed site or density refusal during this search."""

class SweepCounters:
    """Cumulative derivative work counts for a HelmholtzModel."""
    forward: int
    """Number of forward derivative sweeps."""
    reverse: int
    """Number of reverse derivative sweeps."""
    hessian_calls: int
    """Number of Hessian evaluations."""
    recordings: int
    """Number of automatic-differentiation recordings."""

class HelmholtzModel:
    """Native Helmholtz EOS, exact actions and pressure closure.

    Notes
    -----
    Parameters and optional thermochemistry remain fixed for this model. Active parameters declare exact derivative coordinates; no finite differences are used for calculation. See ALG-DENSITY and ALG-ACTIONS in engine/docs/science/algorithms.md.
    """
    solver_options: list[dict[str, Any]]
    """Effective named numerical settings with values and provenance."""
    def __init__(self, parameters: Parameters, active_parameters: Sequence[ActiveParameter] = (),
                 thermochemistry: Optional[ThermochemistryRecord] = None, solver_options: Mapping[str, float] | None = None) -> None:
        """Initialize HelmholtzModel from the declared inputs.

        Parameters
        ----------
        parameters : Parameters
            Admitted native EOS parameters.
        active_parameters : Sequence[ActiveParameter]
            Independent native parameter slots in requested derivative-column order.
        thermochemistry : Optional[ThermochemistryRecord]
            Optional species thermal records in component order; absence leaves ideal/total consumers unavailable.
        solver_options : Mapping[str, float] | None
            Initial value of HelmholtzModel.solver_options; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...
    def state(self, temperature: float, molar_density: float, mole_fractions: Sequence[float],
              derived_properties: bool = False) -> State:
        """Evaluate one homogeneous state on the density chart.

        Parameters
        ----------
        temperature : float
            Positive absolute temperature in K.
        molar_density : float
            Positive molar density in mol/m^3.
        mole_fractions : Sequence[float]
            Length n, dimensionless mole fractions in component_ids order; normalized and charge neutral.
        derived_properties : bool
            Request fixed-composition response functions and partial molar quantities in addition to the ordinary state.

        Returns
        -------
        State
            State at the supplied T, rho and composition.

        Notes
        -----
        No pressure root is solved. Invalid physical inputs raise ParameterError; native evaluation failures raise EosError. Derived quantities preserve their independent unavailability causes.
        """
        ...
    def neutral_property(self, temperature: float, molar_density: float, property: NeutralProperty, order: DerivativeOrder) -> NeutralPropertyDerivatives:
        """Evaluate a single-neutral-component property and exact partials.

        Parameters
        ----------
        temperature : float
            Positive absolute temperature in K.
        molar_density : float
            Positive molar density in mol/m^3.
        property : NeutralProperty
            Select dimensionless residual Helmholtz energy or pressure in Pa.
        order : DerivativeOrder
            Value, Gradient or Hessian; PressureSlope is refused here.

        Returns
        -------
        NeutralPropertyDerivatives
            Value, gradient and row-major Hessian on [T (K), rho (mol/m^3)].

        Notes
        -----
        Restricted to one neutral component without association sites, active parameters, electrolyte or permittivity declarations. Arrays not requested by order are empty.
        """
        ...
    def density_actions(self, temperature: float, molar_density: float,
                        mole_fractions: Sequence[float], observables: Sequence[PropertyObservable],
                        directions: Sequence[DensityActionDirection]) -> PropertyActionBatch:
        """Contract exact homogeneous-property derivatives on the density chart.

        Parameters
        ----------
        temperature : float
            Positive absolute temperature in K.
        molar_density : float
            Positive molar density in mol/m^3.
        mole_fractions : Sequence[float]
            Length n, dimensionless mole fractions in component_ids order; normalized and charge neutral.
        observables : Sequence[PropertyObservable]
            Ordered property selections; values and contractions use PropertyObservable units.
        directions : Sequence[DensityActionDirection]
            Ordered physical increments defining one mixed property contraction.

        Returns
        -------
        PropertyActionBatch
            PropertyActionBatch with central values, contractions and per-result causes.

        Notes
        -----
        One through three directions form the mixed contraction. All action units equal the observable unit. LogFugacityCoefficient is refused on this chart; see ALG-ACTIONS.
        """
        ...
    def pressure_actions(self, temperature: float, pressure: float,
                         mole_fractions: Sequence[float], branch: DensityBranch,
                         observables: Sequence[PropertyObservable],
                         directions: Sequence[PressureActionDirection],
                         anchor_density: Optional[float] = None) -> PropertyActionBatch:
        """Contract exact property derivatives through pressure closure.

        Parameters
        ----------
        temperature : float
            Positive absolute temperature in K.
        pressure : float
            Positive prescribed pressure in Pa.
        mole_fractions : Sequence[float]
            Length n, dimensionless mole fractions in component_ids order; normalized and charge neutral.
        branch : DensityBranch
            Requested density branch; a stable root need not be globally stable.
        observables : Sequence[PropertyObservable]
            Ordered property selections; values and contractions use PropertyObservable units.
        directions : Sequence[PressureActionDirection]
            Ordered physical increments; mixed actions contract all directions, whereas each first density action is reported separately.
        anchor_density : Optional[float]
            Optional previous density in mol/m^3 for the same branch; None requests an unanchored search.

        Returns
        -------
        PropertyActionBatch
            PropertyActionBatch with pressure-root evidence and density actions.

        Notes
        -----
        One through three directions form the mixed contraction on the admitted branch. Each log-fugacity observable expands into component order. A central value can remain available when an action is refused.
        """
        ...
    def packing_fraction_coefficient(self, temperature: float, composition: list[float]) -> float:
        """Evaluate the coefficient in eta = kappa(T,x) rho.

        Parameters
        ----------
        temperature : float
            Positive absolute temperature in K.
        composition : list[float]
            Length n, dimensionless mole fractions in component_ids order; normalized and charge neutral.

        Returns
        -------
        float
            Packing coefficient kappa in m^3/mol; eta is dimensionless.
        """
        ...
    def state_at_pressure(self, temperature: float, pressure: float, mole_fractions: Sequence[float], branch: DensityBranch,
                          anchor_density: Optional[float] = None, derived_properties: bool = False) -> State:
        """Close the EOS pressure on a declared density branch.

        Parameters
        ----------
        temperature : float
            Positive absolute temperature in K.
        pressure : float
            Positive prescribed pressure in Pa.
        mole_fractions : Sequence[float]
            Length n, dimensionless mole fractions in component_ids order; normalized and charge neutral.
        branch : DensityBranch
            Requested density branch; a stable root need not be globally stable.
        anchor_density : Optional[float]
            Optional previous density in mol/m^3 for the same branch; None requests an unanchored search.
        derived_properties : bool
            Request fixed-composition response functions and partial molar quantities in addition to the ordinary state.

        Returns
        -------
        State
            State at the resolved molar density, including density convergence and root-count evidence.

        Notes
        -----
        The search follows ALG-DENSITY. A structurally incomplete root search raises DensityStructureError; an iterative refusal may raise ConvergenceError. A returned root is not a global stability certificate.
        """
        ...
    def molar_density_upper_bound(self, temperature: float, composition: Sequence[float]) -> float:
        """Read the packing-derived density limit.

        Parameters
        ----------
        temperature : float
            Positive absolute temperature in K.
        composition : Sequence[float]
            Length n, dimensionless mole fractions in component_ids order; normalized and charge neutral.

        Returns
        -------
        float
            Geometric molar-density bound in mol/m^3; finite EOS evaluation throughout the interval is not established.
        """
        ...
    def density_diagnostics(self, temperature: float, pressure: float, composition: Sequence[float],
                            nodes: Sequence[float], density_tolerance: float, log_tolerance: float,
                            max_evaluations: int, max_seconds: float, adaptive: bool,
                            finite_anchor: Optional[float], finite_upper: float,
                            evaluation_order: DerivativeOrder) -> DensityDiagnostics:
        """Inspect a declared finite isotherm grid and refine crossings.

        Parameters
        ----------
        temperature : float
            Positive absolute temperature in K.
        pressure : float
            Positive prescribed pressure in Pa.
        composition : Sequence[float]
            Length n, dimensionless mole fractions in component_ids order; normalized and charge neutral.
        nodes : Sequence[float]
            At least two strictly increasing positive densities in mol/m^3 below the geometric cap.
        density_tolerance : float
            Positive absolute density bracket width in mol/m^3.
        log_tolerance : float
            Positive dimensionless logarithmic density bracket width.
        max_evaluations : int
            Positive count limiting diagnostic EOS evaluations.
        max_seconds : float
            Positive elapsed-time limit in seconds.
        adaptive : bool
            Whether to add adaptive density samples.
        finite_anchor : Optional[float]
            Optional finite density in mol/m^3 starting the evaluability-boundary search.
        finite_upper : float
            Upper density in mol/m^3 for the optional finite-domain search.
        evaluation_order : DerivativeOrder
            Gradient, PressureSlope or Hessian for finite-domain admission.

        Returns
        -------
        DensityDiagnostics
            DensityDiagnostics with pressure roots, slope zeros and unresolved intervals.

        Notes
        -----
        Both absolute and logarithmic bracket-width criteria apply. Unresolved finite searches retain their cause and do not prove a unique root.
        """
        ...
    def counters(self) -> SweepCounters:
        """Read cumulative differentiation work.

        Returns
        -------
        SweepCounters
            SweepCounters since model construction or reset_counters.
        """
        ...
    def reset_counters(self) -> None:
        """Reset cumulative derivative-work counts.

        Returns
        -------
        None
            None.
        """
        ...

class PhaseKind:
    """Declared fluid phase branch.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    Liquid: PhaseKind
    """Liquid-side density branch; not proof of global stability."""
    Vapor: PhaseKind
    """Vapor-side density branch; not proof of global stability."""

class AmountKind:
    """Whether a phase amount is solved or prescribed.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    Free: AmountKind
    """Phase amount is an unknown."""
    Pinned: AmountKind
    """Phase amount is fixed as a fraction of the resolved feed amount; zero declares an incipient
    phase.
    """

class FeedKind:
    """Physical basis used to interpret the feed values.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    MoleFractions: FeedKind
    """Dimensionless component mole fractions, giving a one-mole total feed basis."""
    Amounts: FeedKind
    """Absolute component amounts in mol."""
    MassFractions: FeedKind
    """Dimensionless component mass fractions on a one-kg mixture basis; divided by molar masses in
    kg/mol to obtain native molar amounts.
    """
    Molalities: FeedKind
    """Solute amounts in mol/kg of declared solvent; solvent fractions and solvent_mass_kg supply
    the absolute basis.
    """

class AmountSpec:
    """Free phase amount or prescribed feed-relative fraction."""
    kind: AmountKind
    """Whether to solve the amount or pin its feed-relative fraction."""
    fraction: float
    """Dimensionless pinned amount divided by resolved total feed mol; ignored for a free amount.
    Zero denotes an incipient phase.
    """
    def __init__(self) -> None:
        """Initialize AmountSpec from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...
    @staticmethod
    def free() -> AmountSpec:
        """Declare a solved phase amount.

        Returns
        -------
        AmountSpec
            Free AmountSpec; no fixed amount fraction.
        """
        ...
    @staticmethod
    def pinned(value: float) -> AmountSpec:
        """Declare a feed-relative pinned phase amount.

        Parameters
        ----------
        value : float
            Dimensionless phase amount divided by resolved total feed mol.

        Returns
        -------
        AmountSpec
            Pinned AmountSpec; zero declares an incipient phase.
        """
        ...

class FreeValue:
    """Fixed intensive value or solved intensive with an optional seed."""
    free: bool
    """True when the intensive is solved rather than held fixed."""
    value: float
    """Fixed value or declared initial guess; K for temperature and Pa for pressure. A missing
    variable guess is retained natively, not a physical zero.
    """
    def __init__(self) -> None:
        """Initialize FreeValue from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...
    @staticmethod
    def fixed(value: float) -> FreeValue:
        """Declare a fixed temperature or pressure.

        Parameters
        ----------
        value : float
            Positive temperature in K or pressure in Pa, according to its consuming field.

        Returns
        -------
        FreeValue
            Fixed FreeValue.
        """
        ...
    @staticmethod
    def variable(guess: Optional[float] = None) -> FreeValue:
        """Declare a solved intensive and optional seed.

        Parameters
        ----------
        guess : Optional[float]
            Optional positive temperature in K or pressure in Pa; None requests native initialization.

        Returns
        -------
        FreeValue
            Variable FreeValue retaining whether a guess was supplied.
        """
        ...

class PhaseSpec:
    """One declared fluid phase, support and initialization.

    Notes
    -----
    Support and charge admission are owned by compile_problem. A seed supplies no new thermodynamic constraints.
    """
    name: str
    """Caller phase label; it does not select a physical model."""
    kind: PhaseKind
    """Declared liquid or vapor density branch."""
    support: list[int]
    """Ordered component indices allowed in this phase; an explicit support cannot acquire extra
    components.
    """
    support_declared: bool
    """Whether support is explicit; an omitted support may be completed, while a declared empty
    support is rejected.
    """
    amount: AmountSpec
    """Free amount or dimensionless feed-relative pinned fraction."""
    composition_guess: list[float]
    """Dimensionless mole-fraction seed in support order, length len(support), or empty for native initialization. For direct supplied stability trials, use n entries in component order with zero outside support.
    """
    packing_guess: Optional[float]
    """Optional dimensionless packing-fraction seed, not a density in mol/m^3."""
    def __init__(self) -> None:
        """Initialize PhaseSpec from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...

class FeedSpec:
    """Feed declaration with explicit amount or composition basis.

    Notes
    -----
    Native conversion enforces normalization, required molar masses and electroneutrality. Amounts remain absolute; mole fractions establish a one-mole feed.
    """
    kind: FeedKind
    """Basis selecting the interpretation and native conversion of values."""
    values: list[float]
    """Length n in component_ids order: mol for Amounts, dimensionless for mole/mass fractions, or
    relative solvent molar amounts plus solute mol/kg solvent for Molalities.
    """
    solvent: list[int]
    """Solvent component indices for a molality feed; non-solvent entries of values are solute
    molalities.
    """
    solvent_mass_kg: float
    """Positive solvent mass in kg setting the absolute scale of a molality feed; not a total
    mixture mass.
    """
    def __init__(self) -> None:
        """Initialize FeedSpec from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...

class SpeciesReference:
    """Species standard selected for a neutral reaction reference.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    InfiniteDilution: SpeciesReference
    """Solute at infinite dilution in the declared solvent and dilution path."""
    PureLiquid: SpeciesReference
    """Pure-component liquid reference at the source pressure."""
    PureGas: SpeciesReference
    """Pure-component vapor reference at the source pressure."""

class NeutralReferenceDeclaration:
    """Solvent and charge-neutral dilution path for source-standard transfer.

    Notes
    -----
    Mixed pure-species overrides require the admitted raw mole-fraction source basis. See ALG-REFERENCE in engine/docs/science/algorithms.md.
    """
    solvent_index: Optional[int]
    """Index of the single pure neutral solvent, or None before declaration."""
    solute_molality_ratios: list[float]
    """Length n of dimensionless all-solute molality ratios in component_ids order; solvent entry
    is zero and the solute path is charge neutral.
    """
    standard_molality_mol_per_kg: float
    """Positive source standard molality in mol/kg solvent."""
    species_reference: list[SpeciesReference]
    """Optional length-n map in component_ids order; empty retains the default solvent liquid and
    solute infinite-dilution references.
    """
    def __init__(self) -> None:
        """Initialize NeutralReferenceDeclaration from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...

class ReferenceSourceBasis:
    """Declared activity convention of a source equilibrium constant.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    CommonMolalityInfiniteDilution: ReferenceSourceBasis
    """Common standard molality with solute infinite-dilution references; solvent-mass conversion
    is applied once.
    """
    RawMoleFractionInfiniteDilution: ReferenceSourceBasis
    """Raw mole-fraction reference convention; its basis terms cancel in the native
    density-standard conversion.
    """

class NeutralReferenceStatus:
    """Availability or refusal cause of a source-reference calculation.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    Available: NeutralReferenceStatus
    """Requested reference quantity is available."""
    InvalidDeclaration: NeutralReferenceStatus
    """The solvent, path or reference declaration is inadmissible."""
    NonNeutralDeclaration: NeutralReferenceStatus
    """The solvent or solute dilution path violates charge neutrality."""
    UnsupportedSourceBasis: NeutralReferenceStatus
    """The requested source/species reference combination is unsupported."""
    FreeTemperature: NeutralReferenceStatus
    """This reference consumer does not admit the requested free-temperature declaration."""
    FreePressure: NeutralReferenceStatus
    """This reference consumer does not admit the requested free-pressure declaration."""
    DensityUnavailable: NeutralReferenceStatus
    """A required reference density could not be resolved."""
    BranchChange: NeutralReferenceStatus
    """The reference path did not retain its declared density branch."""
    Nonfinite: NeutralReferenceStatus
    """A reference quantity is nonfinite."""
    """Consecutive terminal reference estimates differ beyond the admitted criterion."""

class NeutralReferenceReactionDiagnostic:
    """Reference transfer for one declared reaction."""
    solver_options: list[dict[str, Any]]
    """Effective named numerical settings and provenance."""
    reaction_index: int
    """Zero-based reaction index in the declaration."""
    terminal_contraction: Optional[float]
    """Dimensionless terminal reaction contraction used for standard-state transfer, or None when
    unavailable.
    """
    observed_terminal_change: Optional[float]
    """Dimensionless change between consecutive terminal extrapolations, or None when unavailable;
    a local numerical estimate, not a universal error bound.
    """
    terminal_contraction_pressure_derivative: Optional[float]
    """First partial of the contraction with respect to system pressure at fixed temperature, 1/Pa; exactly zero for a fixed source pressure, or None when the requested column is unavailable or was not evaluated.
    """
    terminal_contraction_pressure_second_derivative: Optional[float]
    """Second partial with respect to system pressure at fixed temperature, 1/Pa^2; exactly zero for a fixed source pressure, or None when the requested column is unavailable or was not evaluated."""
    reference_pressure_pa: float
    """Resolved source reference pressure in Pa."""
    status: NeutralReferenceStatus
    """Reference availability or specific refusal cause."""
    message: str
    """Native explanation retaining the failing reference calculation."""

    diagnostics: PropertyActionDiagnostics
    density_convergence: Convergence | None

class ReactionReference:
    """Source activity basis and pressure of a reaction constant."""
    source_basis: ReferenceSourceBasis
    """Declared source convention used for native standard-state transfer."""
    reference_pressure_pa: Optional[float]
    """Fixed source pressure in Pa, or None to use current system pressure; it does not silently
    become a fixed standard.
    """
    def __init__(self, source_basis: ReferenceSourceBasis,
                 reference_pressure_pa: Optional[float] = None) -> None:
        """Initialize ReactionReference from the declared inputs.

        Parameters
        ----------
        source_basis : ReferenceSourceBasis
            Initial value of ReactionReference.source_basis; see its attribute description for units and ordering.
        reference_pressure_pa : Optional[float]
            Initial value of ReactionReference.reference_pressure_pa; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class ReactionSpec:
    """One stoichiometric reaction and exactly one thermodynamic target.

    Notes
    -----
    No target silently overrides another. A scalar isothermal constant supplies no temperature derivative. See reaction_affinity and reaction_standard_state_conversion in engine/docs/equations.md.
    """
    name: str
    """Reaction identity used for source matching and diagnostics."""
    stoichiometry: list[float]
    """Length n in component_ids order, products positive and reactants negative; dimensionless
    stoichiometric numbers.
    """
    equilibrium_constant: Optional[float]
    """Optional positive dimensionless K on the native EOS standard without reference, or on the declared source basis with reference. Mutually exclusive with affinity and correlation.
    """
    affinity_j_per_mol: Optional[float]
    """Optional reaction affinity target in J/mol of reaction; mutually exclusive with K and
    correlation.
    """
    correlation: Optional[ReactionLogPolynomial]
    """Optional source ln K(T) law with interval and standard identity; mutually exclusive with
    scalar targets.
    """
    reference: Optional[ReactionReference]
    """Optional explicit source-reference declaration; requires the problem's neutral_reference."""
    def __init__(self) -> None:
        """Initialize ReactionSpec from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...

class SolidSpec:
    """Pure stoichiometric solid in equilibrium with declared fluid species.

    Notes
    -----
    Native classification retains solid amount and saturation-index complementarity. See solid_standard_gibbs, solid_saturation_index and solid_amount_complementarity in engine/docs/equations.md.
    """
    name: str
    """Solid identity used for declaration and diagnostics."""
    stoichiometry: list[float]
    """Length n of dissolution numbers for solid -> sum(nu_i A_i), in component_ids order;
    dimensionless and charge neutral.
    """
    equilibrium_constant: Optional[float]
    """Optional positive dimensionless dissolution K; exactly one of this and correlation is
    required.
    """
    correlation: Optional[ReactionLogPolynomial]
    """Optional source ln K(T) law on an explicit standard and interval."""
    reference: Optional[ReactionReference]
    """Source convention and pressure used for native solid-standard conversion."""
    molar_volume_m3_per_mol: float
    """Solid molar volume in m^3/mol for its pressure correction."""
    amount: AmountSpec
    """Free solid amount or Pinned(0) absent solid; other pinned fractions are not admitted."""
    def __init__(self) -> None:
        """Initialize SolidSpec from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...

class StabilityOptions:
    """Finite declared trial search for a tangent-plane descent witness.

    Notes
    -----
    A finite search with no resolved descent is not a global stability certificate. See ALG-STABILITY in engine/docs/science/algorithms.md.
    """
    evaluate_supplied_trials: bool
    """Evaluate supplied trial points directly rather than optimizing their compositions."""
    trials: Optional[list[PhaseSpec]]
    """Explicit trial phases, or None for the admitted default supports and seeds; an empty list
    remains an explicit declaration.
    """
    declaration_budget: Optional[int]
    """Optional count limiting candidate topology declarations."""
    trial_budget: Optional[int]
    """Optional count limiting trial searches."""
    standard_uncertainty_over_rt: list[float]
    """Optional length-n dimensionless species-standard uncertainty estimates in component order."""
    def __init__(self) -> None:
        """Initialize StabilityOptions from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...

class StabilityTrial:
    """One tangent-plane trial, numerical diagnostics and finite-split evidence."""
    phase: PhaseSpec
    """Declared trial phase and support."""
    seed: str
    """Native seed label, not an additional physical constraint."""
    evidence: str
    """Native classification of the observed trial evidence."""
    message: str
    """Trial explanation, including unresolved or refused work."""
    optimizer: EquilibriumResult
    """Retained trial solve and numerical evidence; a solve success alone does not prove descent."""
    raw_coordinates: list[float]
    """Native trial chart coordinates in its compiled order; dimensionless logarithms and scaled
    coordinates are not mole fractions.
    """
    raw_residuals: list[float]
    """Defining-row residuals in the trial's native row normalization."""
    composition: list[float]
    """Length n, normalized dimensionless trial mole fractions in component order when
    reconstructed; otherwise empty.
    """
    split_coordinates: list[float]
    """Finite-split seed blocks for existing phases in declaration order followed by the added trial phase: supported log mole fractions, phase amount in mol and log packing. Multipliers are omitted; empty without a retained finite split."""
    basis: list[list[float]]
    """Rectangular (b, n) conserved trial basis, with b basis rows and n component columns in component order; coefficients are dimensionless. Empty when no compiled trial basis is retained.
    """
    row_transformation: list[list[float]]
    """Rows mapping the native trial constraints into the retained basis; dimensionless matrix in
    producer order.
    """
    quantities: dict[str, float]
    """Named trial quantities; physical units are encoded by each native key and no shared unit
    applies. Error estimates are local, not enclosures.
    """
    split_checks: list[dict[str, float]]
    """Finite donor/split comparisons in trial order; each dictionary retains named values and
    units in its keys.
    """
    donor: Optional[int]
    """Candidate donor phase index, or None without an admitted donor."""
    resolved_descent: bool
    """True only when a finite descent witness exceeds its numerical uncertainty; not a global
    search result.
    """

class StabilityReport:
    """Finite tangent-plane trial evidence for one candidate."""
    completed: bool
    """Whether the declared finite search completed; does not certify global stability."""
    message: str
    """Search conclusion and any stopping limitation."""
    trials: list[StabilityTrial]
    """Trial records in attempted order, including failed or unresolved trials."""
    candidate_coordinates: list[float]
    """Repaired candidate phase blocks in declaration order: supported log mole fractions, phase amount in mol and log packing. Includes an amount entry for every phase and omits multipliers and free-intensive coordinates; this is not a full CompiledSystem point."""
    standard_over_rt: list[float]
    """Length n of species standard Gibbs terms divided by RT, dimensionless, in component order."""
    plane_over_rt: list[float]
    """Length n of the tangent-plane coefficients divided by RT, dimensionless, in component order."""
    present_phases: list[dict[str, float]]
    """Named present-phase consistency quantities per phase; units follow native dictionary keys."""
    work: dict[str, float]
    """Counts for this finite search; absent keys have no retained count. root_calls counts solve_density trace spans, root_evaluations counts energy-evaluation spans inside those roots, and energy_evaluations counts all energy-evaluation spans. trial_invocations counts rooted, distinct trials admitted to direct evaluation or search; reconstruction_evaluations counts completed composition reconstructions. driver_invocations and trial_driver_invocations count shared-driver trial calls; iterations, evaluations and hessian_calls sum their returned work. trial_tangent_plane_invocations, trial_tangent_plane_iterations and trial_tangent_plane_newton_steps report the separate reduced tangent-plane search work. association_diagnostic_evaluations counts extra association diagnostics; finite_split_evaluations counts finite donor/split energy comparisons. All are counts; no elapsed time is retained.
    """

    first_numerical_failure: Optional[Convergence]
    """Earliest observed numerical refusal, separate from final convergence."""

    last_rank: Optional[NumericalRank]
    """Rank of the candidate amount correction, when measured."""

class DiscoveryRecord:
    """Retained candidate topology and the finite discovery sequence.

    Notes
    -----
    Discovery changes declared topology only through admitted finite witnesses and work limits. It does not supply a global equilibrium certificate.
    """
    declaration: EquilibriumProblem
    """Final phase declaration; its phase order owns the returned candidate arrays."""
    candidate: EquilibriumResult
    """Retained candidate and its numerical evidence."""
    rounds: list[StabilityReport]
    """Stability reports in discovery-round order."""
    declarations: list[EquilibriumResult]
    """Attempted candidate solves in retained attempt order."""
    declaration_roles: list[str]
    """Role label per retained declaration attempt."""
    declaration_attempt_coordinates: list[list[float]]
    """Native attempted coordinates per declaration; preserved even when no state was admitted."""
    work: dict[str, float]
    """Accumulated search and shared-driver counts, using StabilityReport.work meanings. driver_invocations, iterations, evaluations and hessian_calls include declaration solves; declaration_invocations counts those driver calls and declaration_initialization_failures counts attempts that failed before invoking the driver. root_calls, root_evaluations and energy_evaluations are recounted over the whole discovery trace. When a phase removal is attempted, removal_component_error_mol and removal_moiety_error_mol retain the latest maximum absolute component and moiety conservation errors in mol. No elapsed time is retained."""


def cooling_onset(compiled: CompiledSystem) -> tuple[EquilibriumResult, list[StabilityReport]]:
    """Run the native cooling-onset search for an admitted compiled declaration.

    Parameters
    ----------
    compiled : CompiledSystem
        Admitted native compiled system; its HelmholtzModel must remain alive.

    Returns
    -------
    tuple[EquilibriumResult, list[StabilityReport]]
        Tuple of the retained EquilibriumResult and ordered StabilityReports.

    Notes
    -----
    The native ALG-EQUILIBRIUM cooling-onset consumer owns eligible declarations, step selection and refusal. Reports retain the finite stability-search limits; no kinetic onset or global certificate is implied.
    """
    ...

def stability_test(compiled: CompiledSystem, central: EquilibriumResult,
                   options: StabilityOptions = ...) -> StabilityReport:
    """Search declared tangent-plane trials around an admitted candidate.

    Parameters
    ----------
    compiled : CompiledSystem
        Admitted native compiled system; its HelmholtzModel must remain alive.
    central : EquilibriumResult
        Admitted candidate from the same compiled declaration.
    options : StabilityOptions
        Explicit finite trial supports, seeds and work limits.

    Returns
    -------
    StabilityReport
        StabilityReport with finite-search evidence, work and unresolved causes.

    Notes
    -----
    Resolved descent is a finite witness. No witness after a finite search is not a global stability proof; see ALG-STABILITY.
    """
    ...

class CriticalDeclaration:
    """Seed for the neutral critical-equation candidate.

    Notes
    -----
    The admitted critical cell has one free neutral phase with free T and P. A critical-equation root is not a stability certificate.
    """
    density_guess_mol_per_m3: float
    """Positive molar-density guess in mol/m^3."""
    mode_guess: list[float]
    """Optional unnormalized amount-perturbation vector in component order; empty selects the
    lowest seed mode.
    """
    def __init__(self) -> None:
        """Initialize CriticalDeclaration from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...

class EnvelopeCoordinate:
    """Dimensionless local coordinate selecting an envelope condition.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    LogTemperature: EnvelopeCoordinate
    """ln(T/(1 K))."""
    LogPressure: EnvelopeCoordinate
    """ln(P/(1 Pa))."""
    LogK: EnvelopeCoordinate
    """ln(x_incipient[k])-ln(x_feed[k]) for the selected component."""

class EnvelopeCondition:
    """One local equation used to trace an incipient phase envelope."""
    kind: EnvelopeCoordinate
    """Logarithmic local coordinate of the condition."""
    component: int
    """Component index used by LogK; ignored for LogTemperature and LogPressure."""
    target: float
    """Dimensionless target value of the selected logarithmic coordinate."""
    def __init__(self) -> None:
        """Initialize EnvelopeCondition from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...

class CriticalRecord:
    """Spectrum and conditioning of a critical-equation candidate.

    Notes
    -----
    These diagnostics support a local critical-equation candidate; they certify no global phase stability. See critical_equations in engine/docs/equations.md.
    """
    eigenvalues: list[float]
    """Dimensionless eigenvalues of the native scaled stability matrix M, in eigensolver order."""
    mode: list[float]
    """Critical mode vector in native component/support order, dimensionless after normalization."""
    cubic: float
    """Native dimensionless cubic-form contraction along the mode."""
    zero_band: float
    """Dimensionless numerical zero band used to interpret the spectrum."""
    condition_estimate: float
    """Dimensionless augmented defining-Jacobian condition estimate."""
    backward_error: float
    """Dimensionless numerical backward-error estimate."""
    rank: int
    """Numerical rank of the augmented defining Jacobian, a count."""

class NumericalRank:
    """Most recent existing factorization's rank decision and relative threshold."""
    basis: str
    dimension: int
    rank: int
    threshold: float

class Quantity:
    """Comparison chart: linear, positive physical quantity, or stored logarithm."""
    Linear: Quantity
    Positive: Quantity
    LogPositive: Quantity

class ComparisonMode:
    """Strict ordinary allowance or the accepted conditional numerical-error rule."""
    Strict: ComparisonMode
    EstimatedError: ComparisonMode

class NumericalEstimate:
    """Local numerical-error radius in the stored quantity and its comparison chart."""
    absolute: float
    chart: float
    available: bool
    resolved: bool
    cause: str
    def __init__(self) -> None: ...

class ComparisonEligibility:
    """Original declaration, branch, inventory, residual and admission outcomes."""
    declaration: bool
    branch: bool
    phase_matching: bool
    inventory: bool
    residual_a: bool
    residual_b: bool
    admission_a: bool
    admission_b: bool
    def __init__(self) -> None: ...

class EquilibriumComparison:
    """Native judgment and differences; the caller retains the comparison operands."""
    judgment: str
    cause: str
    difference: float
    chart_difference: float
    radius_sum: float

class EquilibriumNumericalEstimates:
    """Requested local error estimates at the returned equilibrium point, not enclosures."""
    loading: NumericalEstimate

    basis: str
    cause: str
    source_reference_method: str
    rank: int
    dimension: int
    precision_bits: int
    threshold: float
    condition: float
    inverse_norm: float
    backward_error: float
    residuals: list[float]
    row_terms: list[float]
    row_inner: list[float]
    row_floor: list[float]
    physical_scale: list[float]
    temperature: NumericalEstimate
    pressure: NumericalEstimate
    densities: list[NumericalEstimate]
    packing: list[NumericalEstimate]
    composition: list[NumericalEstimate]
    amounts: list[NumericalEstimate]

def compare_equilibrium_quantity(
    quantity: Quantity, a: float, b: float, ordinary_allowance: float,
    estimate_a: NumericalEstimate, estimate_b: NumericalEstimate,
    eligibility: ComparisonEligibility, mode: ComparisonMode,
) -> EquilibriumComparison:
    """Apply the accepted ordinary or conditional comparison in the native owner."""
    ...

def project_numerical_estimate(value: float, mapped: float, radius: float, positive: bool, constant: bool) -> NumericalEstimate:
    """Project a native log or linear radius, including output displacement and rounding."""
    ...

def equilibrium_numerical_estimates(model: HelmholtzModel | IdealMoleFraction, result: EquilibriumResult) -> EquilibriumNumericalEstimates:
    """Estimate numerical error at the frozen returned state of the declared model."""
    ...

class EquilibriumProblem:
    """Declared fluid/reaction/solid equilibrium and extensive constraints.

    Notes
    -----
    Native compilation owns support, conservation, charge, rank and specification admission. Extensive constraints require the consumed species calorics; missing calorics are not inferred from an isothermal constant. See ALG-EQUILIBRIUM.
    """
    @property
    def solver_options(self) -> list[dict[str, Any]]:
        """Read or set effective named numerical solver settings.

        Returns
        -------
        list[dict[str, Any]]
            Named settings with values and provenance when read; None on assignment.
        """
        ...
    @solver_options.setter
    def solver_options(self, value: Mapping[str, float] | None) -> None:
        """Read or set effective named numerical solver settings.

        Parameters
        ----------
        value : Mapping[str, float] | None
            Optional mapping of admitted numerical setting names to values; None selects defaults.

        Returns
        -------
        None
            Named settings with values and provenance when read; None on assignment.
        """
        ...
    discovery: Optional[StabilityOptions]
    """Optional finite stability/discovery request; None retains the declared topology solve."""
    phases: list[PhaseSpec]
    """Fluid phase declarations in the order used by every phase array."""
    temperature: FreeValue
    """Fixed or free absolute temperature, K."""
    pressure: FreeValue
    """Fixed or free system pressure, Pa."""
    feed: FeedSpec
    """Feed values with an explicit molar, mass or solvent-molality basis."""
    reactions: list[ReactionSpec]
    """Ordered reaction declarations and source targets."""
    standard_gibbs_j_per_mol: list[float]
    """Optional length-n supplied species Gibbs standards, J/mol in component_ids order; not
    independent of conflicting reaction targets.
    """
    topology_declared: bool
    """Whether the supplied phase topology is declared for the solve."""
    row_tolerance: float
    """Caller criterion applied to the maximum absolute classified row residual. Most rows are dimensionless; absent free-solid amount rows retain mol. Numerical floor admission is reported separately.
    """
    enthalpy_j: float | None
    """Optional extensive enthalpy specification in J on the absolute feed basis; not J/mol."""
    entropy_j_per_k: float | None
    """Optional extensive entropy specification in J/K on the absolute feed basis; not J/(mol K)."""
    critical: Optional[CriticalDeclaration]
    """Optional neutral critical-equation declaration."""
    envelope: Optional[EnvelopeCondition]
    """Optional incipient envelope condition."""
    azeotrope: bool
    """Request the admitted binary bubble/dew azeotrope condition; feed composition is its seed."""
    neutral_reference: Optional[NeutralReferenceDeclaration]
    """Optional explicit solvent/path declaration for source-standard transfers."""
    solids: list[SolidSpec]
    """Solid declarations in the order used by solid_amounts and solid_affinity_decades."""
    def __init__(self) -> None:
        """Initialize EquilibriumProblem from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...

class InventoryBound:
    """Conservation-only upper-bound evidence for one species.

    Notes
    -----
    These bounds use conservation and nonnegative amounts only. They do not prove phase equilibrium or stability.
    """
    species: int
    """Component index in component_ids order."""
    status: str
    """Native classification of the bound and its witness; not a thermodynamic equilibrium result."""
    termination: str
    """Linear-program termination explanation."""
    maximum: Optional[float]
    """Maximum species amount in mol, or dimensionless mole fraction when fraction=True; None when
    unresolved or unbounded.
    """
    feed_normalized_amount: Optional[float]
    """Dimensionless maximum amount divided by feed_total for amount bounds; None for a fraction
    bound or unavailable result.
    """
    log_fraction: Optional[float]
    """Natural log of a positive maximum mole fraction, dimensionless; None for zero, amount mode
    or unavailable result.
    """
    balance_residual: float
    """Conservation residual magnitude on the native feed-normalized linear-program chart,
    dimensionless.
    """
    normalization_residual: float
    """Dimensionless fractional-chart normalization residual."""
    lower_bound_residual: float
    """Dimensionless lower-bound violation in the normalized linear program."""
    witness: list[float]
    """Linear-program witness in native column order; fraction mode includes its normalization
    variable, not just component molar amounts.
    """
    amounts_mol: list[float]
    """Finite reconstructed component amounts in mol, component_ids order, when available;
    otherwise empty.
    """
    finite_witness: list[float]
    """Witness of the finite-feasibility check in native normalized column order, or empty without
    an optimal finite witness.
    """
    limiting: bool
    """Whether the result is a limiting bound rather than an attained finite witness."""
    supported: bool
    """Whether the finite inventory can support this species."""
    floor: float
    """Caller-declared dimensionless inventory support floor."""
    lp_error: float
    """Dimensionless native linear-program error estimate; may be nonfinite if unresolved."""
    finite_residual: float
    """Dimensionless maximum finite-witness feasibility residual; may be nonfinite if unresolved."""
    feasibility_termination: str
    """Termination explanation for the finite-witness calculation."""

def inventory_bounds(compiled: CompiledSystem, species: Sequence[int], fraction: bool,
                     floor: float) -> list[InventoryBound]:
    """Bound species inventory using conservation and nonnegativity.

    Parameters
    ----------
    compiled : CompiledSystem
        Admitted native compiled system; its HelmholtzModel must remain alive.
    species : Sequence[int]
        Ordered component indices whose inventory bounds are requested.
    fraction : bool
        True requests dimensionless mole-fraction bounds; False requests absolute mol bounds.
    floor : float
        Finite nonnegative mole-fraction floor; must be zero when fraction=False. In fraction mode, the number of supported species times floor must not exceed one.

    Returns
    -------
    list[InventoryBound]
        One InventoryBound per requested species, preserving input order.

    Notes
    -----
    This calculation has no thermodynamic objective and supplies no phase-stability result.
    """
    ...

class RowKind:
    """Physical equation represented by a compiled equilibrium row.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    ReactionMolarity: RowKind
    """Direct logarithmic molarity mass action, dimensionless."""
    ReactionMassAction: RowKind
    """Direct logarithmic mole-fraction mass action, dimensionless."""
    Chemical: RowKind
    """Chemical-potential stationarity, normalized by RT."""
    Pressure: RowKind
    """Pressure equality in the native dimensionless row normalization."""
    SaturationPressure: RowKind
    """Declared target pressure, normalized by that pressure."""
    Moiety: RowKind
    """Conserved moiety balance normalized on the feed basis."""
    Electro: RowKind
    """Per-phase electroneutrality, dimensionless charge-number contraction."""
    Normalization: RowKind
    """Phase mole-fraction normalization, dimensionless."""
    Enthalpy: RowKind
    """Extensive enthalpy specification in its native dimensionless normalization."""
    Entropy: RowKind
    """Extensive entropy specification in its native dimensionless normalization."""
    CriticalMode: RowKind
    """Critical-mode null-vector equation, native dimensionless normalization."""
    CriticalNormalization: RowKind
    """Critical-mode normalization, dimensionless."""
    CriticalCubic: RowKind
    """Critical cubic-form condition, native dimensionless normalization."""
    EnvelopeCondition: RowKind
    """Selected dimensionless logarithmic envelope condition."""
    AzeotropeCondition: RowKind
    """Binary azeotrope composition equality, dimensionless."""
    SolidAffinity: RowKind
    """Present free solid: -A_s/(RT), dimensionless. Absent free solid: n_s, in mol. Pinned(0) solids have no row."""

class Row:
    """Identity of a defining or classification equation."""
    kind: RowKind
    """Physical equation family."""
    phase: int
    """Zero-based declared phase index where applicable; not meaningful for a system-wide row."""
    index: int
    """Family-specific component, conserved-balance, mode or solid ordinal; not a common species
    index for all rows.
    """
    eqid: str
    """Canonical equation-registry identity owning the row's mathematical definition."""

class PointEvaluation:
    """Defining rows, objective and numerical floor at a compiled point."""
    values: list[float]
    """Length m of defining-row residuals in CompiledSystem.rows order. Most rows are dimensionless; an absent free-solid row is its amount in mol."""
    objective: float
    """Extensive Gibbs objective G/(RT), in mol on the absolute feed basis; a value alone does not establish equilibrium."""
    floor: float
    """Dimensionless estimated residual floor at this point."""

class IdealMoleFraction:
    activity_basis: str
    solution_volume_m3: float | None
    unit_activity_solvent: str | None
    component_ids: list[str]
    charge_number: list[int]
    solver_options: list[dict[str, Any]]
    def __init__(self, component_ids: Sequence[str], charge_number: Sequence[int], solver_options: Mapping[str, float] | None = ..., activity_basis: str | None = ..., solution_volume_m3: float | None = ..., unit_activity_solvent: str | None = ...) -> None: ...

def effective_water_amounts(model: IdealMoleFraction, result: EquilibriumResult, water: str, ions: Sequence[str], molecular_species: Sequence[str]) -> dict[str, float]: ...

class CompiledSystem:
    """Admitted equilibrium chart and exact defining-equation readers.

    Notes
    -----
    The model must outlive this compiled object. The full point has d = primal_size + rank + charge_row_count entries: each fluid phase's log compositions, free amount when present and log packing; then free-solid amounts in solid declaration order; then mode and free-intensive coordinates where present; then dimensionless conserved and electrical multipliers. Free amounts use mol, free temperature uses K and free pressure uses Pa. See ALG-EQUILIBRIUM and equilibrium/formulation.md.
    """
    rows: list[Row]
    """Ordered defining/classification rows; this order owns residual and Jacobian rows."""
    pinned_count: int
    """Number of pinned fluid phases."""
    free_intensive_count: int
    """Number of solved temperature/pressure coordinates."""
    rank: int
    """Rank of the independent conserved-moiety basis, a count."""
    charge_row_count: int
    """Number of per-phase charge constraints."""
    primal_size: int
    """Number of primal state coordinates before conserved/electrical multipliers."""
    generic_nullity: int
    """Nullity measured at the generic seed; not a guarantee of rank at a solution."""
    feed_total: float
    """Total resolved feed amount in mol; one mol for MoleFractions."""
    has_ions: bool
    """Whether the admitted model contains charged species."""
    has_reaction: bool
    """Whether reaction declarations are present."""
    conserves_total_moles: bool
    """Whether the admitted stoichiometry conserves total molar amount."""
    last_rank: NumericalRank | None
    jacobian_rank: int
    """Defining-Jacobian rank measured at the generic point."""
    removed_species: list[int]
    """Component indices eliminated by exact inventory/support restrictions."""
    def solve_saturation(self) -> EquilibriumResult: ...

    def solve(self) -> EquilibriumResult:
        """Solve this declaration from a fresh start through the native equilibrium driver."""
        ...
    def numerical_estimates(self, coordinates: Sequence[float], temperature: float, pressure: float,
                            amounts: Sequence[float], composition: Sequence[float], densities: Sequence[float],
                            packing: Sequence[float]) -> EquilibriumNumericalEstimates:
        """Estimate a frozen endpoint in this declaration; admission remains the reader's original outcome."""
        ...
    def generic_point(self) -> list[float]:
        """Construct the native generic chart seed.

        Returns
        -------
        list[float]
            Full point of length d in the compiled chart; a seed, not an equilibrium result.
        """
        ...
    def evaluate(self, point: Sequence[float]) -> PointEvaluation:
        """Evaluate native defining rows at a full chart point.

        Parameters
        ----------
        point : Sequence[float]
            Full native point in compiled chart order, length d including multipliers; not a flattened physical composition.

        Returns
        -------
        PointEvaluation
            PointEvaluation with m defining residuals, objective and floor.
        """
        ...
    def jacobian(self, point: Sequence[float]) -> list[float]:
        """Differentiate defining rows with respect to the full chart.

        Parameters
        ----------
        point : Sequence[float]
            Full native point in compiled chart order, length d including multipliers; not a flattened physical composition.

        Returns
        -------
        list[float]
            Flattened row-major (m, d) Jacobian; rows follow rows and columns follow the full point.
        """
        ...
    def lever_jacobian(self, point: Sequence[float]) -> list[float]:
        """Differentiate defining rows with respect to external physical controls.

        Parameters
        ----------
        point : Sequence[float]
            Full native point in compiled chart order, length d including multipliers; not a flattened physical composition.

        Returns
        -------
        list[float]
            Flattened row-major (m, 2+n+a) matrix in [T, P, feed amounts, active parameters] order; each column is row units per K, Pa, mol or native parameter unit.

        Notes
        -----
        A free T/P uses its solved state coordinate; the corresponding external control is unused. Source-reference lever requests retain native refusal limits. See ALG-ACTIONS.
        """
        ...
    def objective_gradient(self, point: Sequence[float]) -> list[float]:
        """Differentiate the native optimizer objective; feasibility cells return a zero vector.

        Parameters
        ----------
        point : Sequence[float]
            Full native point in compiled chart order, length d including multipliers; not a flattened physical composition.

        Returns
        -------
        list[float]
            Length d in full-point order. Gibbs cells differentiate G/(RT), in mol, with respect to each chart coordinate. Feasibility cells use a zero optimizer objective and do not differentiate PointEvaluation.objective.
        """
        ...
    def hessian(self, point: Sequence[float], constraint_multipliers: Sequence[float] = ...) -> list[float]:
        """Evaluate exact native Lagrangian curvature.

        Parameters
        ----------
        point : Sequence[float]
            Full native point in compiled chart order, length d including multipliers; not a flattened physical composition.
        constraint_multipliers : Sequence[float]
            Optional optimizer constraint multipliers in native constraint order; empty requests the default curvature.

        Returns
        -------
        list[float]
            Flattened row-major (d, d) optimizer Lagrangian Hessian in full-point order. Gibbs cells include G/(RT) curvature; feasibility cells have a zero objective and return constraint curvature only.
        """
        ...

    first_numerical_failure: Optional[Convergence]
    """Earliest observed numerical refusal, separate from final convergence."""

    rank_records: list[NumericalRank]
    """Latest observed rank per named operation; no extra admission rule."""


class LinearSolveDiagnostics:
    """Local first-order error estimates in the named system's unscaled coordinates.

    Attributes
    ----------
    system : str
        Defining equations and coordinate basis.
    reason : str
        Precision limitation or rejected-attempt cause, when present.
    residual_norm : float
        Achieved infinity norm of the defining rows.
    evaluation_floor : float
        Estimated row-evaluation error, including association.
    jacobian_inverse_norm : float
        Infinity norm of the inverse in the original row and coordinate scales.
    condition_estimate : float
        Jacobian infinity norm times its inverse norm.
    forward_error_estimate : float
        Estimated absolute coordinate error; for log densities, expm1 gives relative density error.
    precision_bits : int
        Working significand bits.
    refinement_iterations : int
        Extended corrections actually taken.
    """
    system: str
    reason: str
    residual_norm: float
    evaluation_floor: float
    jacobian_inverse_norm: float
    condition_estimate: float
    forward_error_estimate: float
    precision_bits: int
    refinement_iterations: int

class IpoptIteration:
    """Last callback tuple in Ipopt's internal optimization and scaling basis."""
    mode: int
    """Ipopt callback mode, a count or mode ordinal."""
    iteration: int
    """Ipopt callback iteration, a count or mode ordinal."""
    line_search_trials: int
    """Ipopt callback line search trials, a count or mode ordinal."""
    objective: float
    """Unmodified Ipopt callback objective; units follow Ipopt's internal basis."""
    inf_pr: float
    """Unmodified Ipopt callback inf_pr; units follow Ipopt's internal basis."""
    inf_du: float
    """Unmodified Ipopt callback inf_du; units follow Ipopt's internal basis."""
    mu: float
    """Unmodified Ipopt callback mu; units follow Ipopt's internal basis."""
    step_norm: float
    """Unmodified Ipopt callback step_norm; units follow Ipopt's internal basis."""
    regularization: float
    """Unmodified Ipopt callback regularization; units follow Ipopt's internal basis."""
    alpha_du: float
    """Unmodified Ipopt callback alpha_du; units follow Ipopt's internal basis."""
    alpha_pr: float
    """Unmodified Ipopt callback alpha_pr; units follow Ipopt's internal basis."""

class IpoptPoint:
    """Original TNLP final point before polish; absent without finalization."""
    tnlp_coordinates: list[float]
    """Original TNLP tnlp coordinates; see ALG-TERMINATION for the chart."""
    coordinate_scale: list[float]
    """Original TNLP coordinate scale; see ALG-TERMINATION for the chart."""
    tnlp_lower_multipliers: list[float]
    """Original TNLP tnlp lower multipliers; see ALG-TERMINATION for the chart."""
    tnlp_upper_multipliers: list[float]
    """Original TNLP tnlp upper multipliers; see ALG-TERMINATION for the chart."""
    tnlp_constraint_multipliers: list[float]
    """Original TNLP tnlp constraint multipliers; see ALG-TERMINATION for the chart."""
    declaration: EquilibriumProblem
    """Originating declaration; it may differ from the returned state's declaration."""
    rows: list[Row]
    """Originating defining rows."""
    constraint_rows: list[int]
    """Indices into this point's own rows, in TNLP multiplier order."""
    objective_factor: float
    """Requested objective normalization, 1/feed_total; not applied again to the arrays."""

class IpoptDiagnostics:
    """Bounded observations of one invoked Ipopt attempt."""
    last: Optional[IpoptIteration]
    """Last callback tuple, or None if no callback ran."""
    first_restoration_iteration: Optional[int]
    """Observed first restoration iteration, or None when unavailable."""
    initialization_return: Optional[int]
    """Observed initialization return, or None when unavailable."""
    application_return: Optional[int]
    """Observed application return, or None when unavailable."""
    restoration_entries: int
    """Attempt-local restoration entries, a count."""
    callbacks: int
    """Attempt-local callbacks, a count."""
    warm_initialized: bool
    """Whether this attempt was warm initialized."""
    predicted: bool
    """Whether this attempt was predicted."""
    physical_amounts: bool
    """Whether this attempt was physical amounts."""
    failed_point: Optional[IpoptPoint]
    """Original TNLP tuple on unsuccessful termination, or None without one."""

class EquilibriumResult:
    rejected_reduction: EquilibriumResult | None
    terminal_numerical_failure: Convergence | None
    last_full_evaluation: Convergence | None
    phase_pressures: list[float]

    def saturation_values(self) -> dict[str, Any]: ...

    """Solved phases and numerical evidence on the declared feed basis.

    Notes
    -----
    For discovery, discovery.declaration.phases owns the final phase ordering. Temperature/pressure are None and physical state lists are empty after a refused solve. Diagnostics can survive without a state. See ALG-EQUILIBRIUM; a successful local solve is not a global stability or physical-validation claim.
    """
    solver_settings: dict[str, float]
    """Effective named Ipopt settings; units and meanings follow the named solver option."""
    diagnostic_options_file: str
    """Native diagnostic options-file path, or an empty string without one."""
    linear_solver: str
    """Actual Ipopt linear solver name; empty when no Ipopt solve ran."""
    solver_options: list[dict[str, Any]]
    """Effective named Engine numerical settings and provenance."""
    discovery: Optional[DiscoveryRecord]
    """Finite discovery sequence, or None for an ordinary declared solve."""
    neutral_reference_diagnostics: list[NeutralReferenceReactionDiagnostic]
    """Reference-transfer evidence in declared reaction order where consumed."""
    success: bool
    """Whether the native solve admitted a state; does not establish global equilibrium or physical
    validation.
    """
    message: str
    """Native conclusion and refusal explanation."""
    solver_status: str
    """Optimizer termination name; separate from residual and topology admission."""
    temperature: Optional[float]
    """Resolved temperature in K, or None without an admitted state."""
    pressure: Optional[float]
    """Resolved pressure in Pa, or None without an admitted state."""
    amounts: list[float]
    """Length p of fluid phase amounts in mol on the feed basis, in declaration order; empty
    without a state.
    """
    mole_fractions: list[float]
    """Flattened row-major (p, n) dimensionless mole fractions; entry phase*n+component, empty
    without a state.
    """
    molar_densities: list[float]
    """Length p of phase molar densities in mol/m^3, empty without a state."""
    packing_fractions: list[float]
    """Length p of dimensionless packing fractions, empty without a state."""
    multipliers: list[float]
    """Dimensionless conserved-moiety and electrical multipliers in compiled order."""
    constraint_multipliers: list[float]
    """Gibbs cells: reconstructed extensive-objective multipliers for rows excluding Chemical, Pressure and SolidAffinity, preserving rows order. Material and charge values are -feed_total*lambda and -phase_amount*phi; other retained rows are zero. Feasibility cells: raw Ipopt constraint multipliers, or zeros for a Newton result.
    """
    residuals: list[float]
    """Unclassified defining rows in rows order. Most are dimensionless; an absent free-solid row is its amount in mol."""
    classified_residuals: list[float]
    """Rows used for admission, including free-amount weighting and trace-inventory normalization where applicable. Most are dimensionless; absent free-solid amount rows retain mol.
    """
    coordinates: list[float]
    """Full native state chart including multipliers; follows the compiled point order and is not a
    flat physical composition array.
    """
    phase_count: int
    """Number p of declared fluid phases."""
    iterations: int
    """Iterations including retained seed sub-solves where the producer accumulates them."""
    evaluations: int
    """Native equilibrium evaluation count."""
    hessian_calls: int
    """Exact Hessian evaluation count."""
    objective: float
    """Extensive Gibbs objective G/(RT), in mol on the absolute feed basis, at the retained candidate."""
    reduced_hessian_minimum: float
    """Minimum eigenvalue of the native reduced Hessian; local curvature evidence, not a global
    certificate.
    """
    convergence: Convergence
    """Residual, effective criterion, floor and stopping evidence."""
    last_rank: NumericalRank | None
    """Most recent evaluated rank decision; its basis identifies the matrix."""
    numerical_estimates: EquilibriumNumericalEstimates | None
    """Endpoint estimates on detailed capture; None when unrequested."""
    linear_solve_diagnostics: LinearSolveDiagnostics | None
    """Reduced-system accuracy at the returned state, or explicitly rejected-attempt evidence on failure."""
    rows: list[Row]
    """Equation identities in residual order."""
    unknown_count: int
    """Full state-chart dimension including multipliers."""
    jacobian_rank: int
    """Reported defining-Jacobian rank."""
    generic_nullity: int
    """Nullity of the generic compiled chart."""
    requested_tolerance_met: Optional[bool]
    """Whether the finite final residual meets the caller criterion, or None without a finite final
    residual; floor convergence need not be True.
    """
    topology_event: TopologyEvent
    """First vanished free phase, otherwise first coincident pair, when detected."""
    topology_phase: Optional[int]
    """Index of the vanished phase, or None without this event."""
    topology_pair: Optional[tuple[int, int]]
    """Declared phase-index pair for coalescence, or None without this event."""
    critical: Optional[CriticalRecord]
    """Critical-equation diagnostics, or None outside a retained critical candidate."""
    solid_amounts: list[float]
    """Amount in mol per declared solid, declaration order; empty without retained solid values."""
    solid_affinity_decades: list[float]
    """Saturation index in log10 decades per declared solid; interpreted with solid amount and
    complementarity.
    """
    def phase_mole_fractions(self, phase: int) -> list[float]:
        """Read one declared fluid phase composition.

        Parameters
        ----------
        phase : int
            Zero-based phase ordinal in the result declaration; invalid indices are refused.

        Returns
        -------
        list[float]
            Length n of dimensionless mole fractions in component_ids order.
        """
        ...

    first_numerical_failure: Optional[Convergence]
    """Earliest observed numerical refusal, separate from final convergence."""

    rank_records: list[NumericalRank]
    """Latest observed rank per named operation; no extra admission rule."""

    ipopt_diagnostics: Optional[IpoptDiagnostics]
    """Latest actually invoked Ipopt attempt, including initialization-only attempts."""

    recovered_ipopt_diagnostics: Optional[IpoptDiagnostics]
    """Earliest superseded restoration-bearing attempt."""


class TopologyEvent:
    """Candidate topology change detected by the native solve.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    None_: TopologyEvent
    """No vanished-phase or coincident-phase event detected."""
    PhaseDisappeared: TopologyEvent
    """First vanished free phase in declaration order."""
    PhasesCoalesced: TopologyEvent
    """First coincident pair in lexicographic phase-index order, if no free phase vanished."""

class ContinuationStop:
    """Terminal condition of a declared continuation or envelope trace.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    reached_target: ContinuationStop
    """Ordinary path reached its target, or an envelope trace completed its closing step."""
    initial_failure: ContinuationStop
    """No admitted starting state."""
    invalid_path: ContinuationStop
    """Path declarations or controls are inadmissible."""
    singular_jacobian: ContinuationStop
    """Required defining Jacobian is singular."""
    branch_lost: ContinuationStop
    """The declared density branch was lost."""
    phase_disappeared: ContinuationStop
    """A free phase vanished."""
    phases_coalesced: ContinuationStop
    """Two declared phases became coincident."""
    left_domain: ContinuationStop
    """The path left an admitted physical/source domain."""
    infeasible: ContinuationStop
    """A corrector could not satisfy the declared constraints; no global infeasibility proof."""
    nonfinite: ContinuationStop
    """A required path value is nonfinite."""
    stalled: ContinuationStop
    """The corrector or path stopped making admitted progress."""
    step_too_small: ContinuationStop
    """A halved retry step fell below minimum_step."""
    budget_exhausted: ContinuationStop
    """A step/work/retry budget was exhausted."""
    critical_endpoint: ContinuationStop
    """Envelope tracing reached its declared critical stopping condition."""

class ContinuationOptions:
    """Step and retry controls for continuation and phase envelopes.

    Notes
    -----
    Envelope mode ignores initial_step and the target problem after the shared rejection of solids in either endpoint. Explicit T/P guesses are optional; absent free-intensive guesses use native initialization. After an accepted envelope step the predictor density error sizes the next step; retryable ordinary rejection halves it. A failed critical crossing, closing step or terminal refusal ends the trace with accepted states retained. See ALG-EQUILIBRIUM.
    """
    initial_step: float
    """Largest dimensionless fraction of the ordinary start-to-target path; rejected steps halve
    and accepted ordinary steps never grow. Ignored for envelope tracing.
    """
    minimum_step: float
    """Positive minimum step in the active path coordinate; a smaller halved step stops the trace."""
    max_steps: int
    """Maximum number of continuation steps, a count."""
    max_retries_per_step: int
    """Maximum consecutive retryable rejections; exceeding it stops with budget_exhausted and
    preserves the cause.
    """
    envelope_max_log_step: float
    """Maximum dimensionless envelope log-coordinate increment; bounds ordinary envelope steps."""
    envelope_straddle: float
    """Positive near-critical absolute log-K straddle threshold; requires 2*envelope_straddle <=
    envelope_max_log_step.
    """
    @overload
    def __init__(self) -> None:
        """Initialize ContinuationOptions from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...
    @overload
    def __init__(self, *, initial_step: float, minimum_step: float, max_steps: int, max_retries_per_step: int,
                 envelope_max_log_step: float = 0.0, envelope_straddle: float = 0.0) -> None:
        """Initialize ContinuationOptions from the declared inputs.

        Parameters
        ----------
        initial_step : float
            Initial value of ContinuationOptions.initial_step; see its attribute description for units and ordering.
        minimum_step : float
            Initial value of ContinuationOptions.minimum_step; see its attribute description for units and ordering.
        max_steps : int
            Initial value of ContinuationOptions.max_steps; see its attribute description for units and ordering.
        max_retries_per_step : int
            Initial value of ContinuationOptions.max_retries_per_step; see its attribute description for units and ordering.
        envelope_max_log_step : float
            Initial value of ContinuationOptions.envelope_max_log_step; see its attribute description for units and ordering.
        envelope_straddle : float
            Initial value of ContinuationOptions.envelope_straddle; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class ContinuationRequest:
    """Start, target and physical parameter endpoints of a native trace."""
    start_problem: EquilibriumProblem
    """Starting physical declaration; an envelope condition selects envelope tracing."""
    target_problem: EquilibriumProblem
    """Ordinary-path endpoint. Envelope tracing ignores it after the shared rejection of solids in either endpoint."""
    start_result: Optional[EquilibriumResult]
    """Optional admitted starting result on the same declaration, or None to solve the start."""
    options: ContinuationOptions
    """Declared steps, retry limits and envelope controls."""
    parameter_targets: Sequence[float]
    """Physical endpoint values in model active-parameter order, with each family's native units;
    empty means no parameter path.
    """
    def __init__(self) -> None:
        """Initialize ContinuationRequest from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...

class ContinuationResult:
    """Accepted states, work and terminal cause of a continuation."""
    success: bool
    """Whether the trace completed its declared success condition; inspect stop and retained
    states.
    """
    message: str
    """Trace conclusion or terminal refusal explanation."""
    stop: ContinuationStop
    """Terminal continuation condition."""
    path_values: list[float]
    """One dimensionless path fraction per accepted ordinary state, or the envelope condition
    target per accepted envelope state.
    """
    parameter_values: list[list[float]]
    """Physical active-parameter vector per accepted state only for a declared parameter path, in model active order and family units; empty for ordinary feed/intensive paths and envelope traces.
    """
    states: list[EquilibriumResult]
    """Accepted EquilibriumResults in path order; retained on later failure."""
    reached_path: float
    """Path fraction or envelope target of the last accepted corrector step, dimensionless; remains zero before the first accepted step, even when an initialized envelope has a nonzero starting target."""
    attempted_path: Optional[float]
    """Path fraction or envelope target of the most recent rejected attempt, retained across later accepted steps; None if no attempt was rejected."""
    event_phase: Optional[int]
    """Vanished phase index, or None without this topology event."""
    event_pair: Optional[tuple[int, int]]
    """Coalesced phase-index pair, or None without this event."""
    attempts: int
    """Number of attempted corrector steps."""
    accepted_steps: int
    """Number of accepted continuation steps."""
    retries: int
    """Number of rejected retry attempts."""
    last_result: Optional[EquilibriumResult]
    """Result retained from the most recent rejected attempt, possibly refused; None if no rejected attempt retained a result. The latest accepted result is in states."""
    cause_stop: Optional[ContinuationStop]
    """Underlying retry/terminal stop, or None when no distinct cause is retained."""
    cause_convergence: Optional[Convergence]
    """Underlying corrector convergence evidence, or None without a record."""
    cause_message: Optional[str]
    """Underlying refusal explanation, or None without a distinct cause."""
    conditions: list[EnvelopeCondition]
    """Envelope condition per accepted envelope state; empty for ordinary continuation."""
    critical_points: list[EquilibriumResult]
    """Critical-equation candidate results in encounter order, not globally validated critical
    points.
    """
    critical_seeds: list[list[float]]
    """Per-candidate [ln(T/(1 K)), ln(P/(1 Pa)), ln(rho/(1 mol/m^3))] seeds, dimensionless."""

    first_numerical_failure: Optional[Convergence]
    """Earliest observed numerical refusal, separate from final convergence."""

    rank_records: list[NumericalRank]
    """Latest observed rank per named operation; no extra admission rule."""

    ipopt_diagnostics: Optional[IpoptDiagnostics]
    """Latest actually invoked Ipopt attempt, including initialization-only attempts."""

    recovered_ipopt_diagnostics: Optional[IpoptDiagnostics]
    """Earliest superseded restoration-bearing attempt."""


def continue_equilibrium(model: HelmholtzModel, request: ContinuationRequest) -> ContinuationResult:
    """Trace a declared equilibrium path or incipient phase envelope.

    Parameters
    ----------
    model : HelmholtzModel
        Native Helmholtz model with declared parameters, active coordinates and optional thermochemistry.
    request : ContinuationRequest
        Start/target declarations, physical parameter endpoints and explicit continuation controls.

    Returns
    -------
    ContinuationResult
        ContinuationResult retaining accepted states and the terminal cause.

    Notes
    -----
    Ordinary continuation follows the start-to-target path. Both endpoint declarations are first checked for unsupported solids; an envelope start then selects the native envelope controller and ignores the target declaration.
    """
    ...


class SolvedStateObservableKind:
    """Quantity extracted from the admitted equilibrium state.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    Temperature: SolvedStateObservableKind
    """System temperature, K."""
    Pressure: SolvedStateObservableKind
    """System pressure, Pa."""
    PhaseAmount: SolvedStateObservableKind
    """Fluid phase amount, mol on the absolute feed basis."""
    PhaseComponentAmount: SolvedStateObservableKind
    """Component amount in one phase, mol."""
    PhaseMoleFraction: SolvedStateObservableKind
    """Component mole fraction in one phase, dimensionless."""
    PhaseVolume: SolvedStateObservableKind
    """Fluid phase volume, m^3."""
    PhaseMolarDensity: SolvedStateObservableKind
    """Fluid molar density, mol/m^3."""
    PhasePackingFraction: SolvedStateObservableKind
    """Phase packing fraction, dimensionless."""
    PhasePressure: SolvedStateObservableKind
    """EOS phase pressure, Pa."""
    PhaseResidualChemicalPotentialOverRt: SolvedStateObservableKind
    """Component residual chemical potential divided by RT, dimensionless."""
    PhaseProperty: SolvedStateObservableKind
    """Selected homogeneous molar property in PropertyObservable units; log fugacity is not
    admitted by this solved-state consumer.
    """
    PhaseLogActivity: SolvedStateObservableKind
    """ln[(rho*x_i/(1 mol/m^3))*exp(mu_res_i/(RT))], dimensionless on EOS_STANDARD_STATE_ID;
    excludes an added source-standard term.
    """
    PhaseGibbsMolarVolume: SolvedStateObservableKind
    """Gibbs-conjugate phase molar volume, m^3/mol; includes admitted source-standard pressure
    dependence.
    """
    SystemEnthalpy: SolvedStateObservableKind
    """Extensive system enthalpy, J on the absolute feed basis."""
    EquilibriumIsobaricHeatCapacity: SolvedStateObservableKind
    """Reacting equilibrium isobaric heat capacity normalized by feed mol, J/(mol K); not
    frozen-composition Cp.
    """
    SolidAmount: SolvedStateObservableKind
    """Solid amount in mol; phase selects the solid ordinal."""
    SolidSaturationIndex: SolvedStateObservableKind
    """Solid saturation index in log10 decades; phase selects the solid ordinal."""

class SolvedStateObservable:
    """Observable with declared phase/component selection."""
    kind: SolvedStateObservableKind
    """Quantity and unit selected; see SolvedStateObservableKind."""
    phase: int
    """Fluid phase ordinal, or solid ordinal for SolidAmount/SolidSaturationIndex; ignored for a
    system quantity.
    """
    component: int
    """Component index for component observables; ignored otherwise."""
    property: PropertyObservable
    """Molar property selection when kind is PhaseProperty; ignored for other kinds."""
    def __init__(self, kind: SolvedStateObservableKind, phase: int = 0, component: int = 0,
                 property: PropertyObservable = ...) -> None:
        """Initialize SolvedStateObservable from the declared inputs.

        Parameters
        ----------
        kind : SolvedStateObservableKind
            Initial value of SolvedStateObservable.kind; see its attribute description for units and ordering.
        phase : int
            Initial value of SolvedStateObservable.phase; see its attribute description for units and ordering.
        component : int
            Initial value of SolvedStateObservable.component; see its attribute description for units and ordering.
        property : PropertyObservable
            Initial value of SolvedStateObservable.property; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class SolvedStateActionDirection:
    """Physical external increments for an implicit equilibrium action.

    Notes
    -----
    Every path variable is dimensionless, so first and mixed contractions retain observable units. Directions do not authorize differentiation through a support, topology or source-interval change.
    """
    temperature_k: float
    """External temperature increment in K; must be zero when temperature is solved."""
    pressure_pa: float
    """External pressure increment in Pa; must be zero when pressure is solved."""
    feed_amounts_mol: list[float]
    """Length n of feed-amount increments in mol, component order; not mole-fraction increments."""
    active_parameters: list[float]
    """Physical parameter increments in model active order, each family's native units."""
    reaction_coefficients: list[list[float]]
    """One [da, db, dc, dd] vector per declared reaction law, or empty for no source motion; units
    are 1, K, 1 and 1/K, evaluated along the temperature path.
    """
    enthalpy_j: float
    """Increment of the extensive enthalpy target in J."""
    entropy_j_per_k: float
    """Increment of the extensive entropy target in J/K."""
    def __init__(self, temperature_k: float, pressure_pa: float, feed_amounts_mol: Sequence[float],
                 active_parameters: Sequence[float], reaction_coefficients: Sequence[Sequence[float]] = ...,
                 enthalpy_j: float = ..., entropy_j_per_k: float = ...) -> None:
        """Initialize SolvedStateActionDirection from the declared inputs.

        Parameters
        ----------
        temperature_k : float
            Initial value of SolvedStateActionDirection.temperature_k; see its attribute description for units and ordering.
        pressure_pa : float
            Initial value of SolvedStateActionDirection.pressure_pa; see its attribute description for units and ordering.
        feed_amounts_mol : Sequence[float]
            Initial value of SolvedStateActionDirection.feed_amounts_mol; see its attribute description for units and ordering.
        active_parameters : Sequence[float]
            Initial value of SolvedStateActionDirection.active_parameters; see its attribute description for units and ordering.
        reaction_coefficients : Sequence[Sequence[float]]
            Initial value of SolvedStateActionDirection.reaction_coefficients; see its attribute description for units and ordering.
        enthalpy_j : float
            Initial value of SolvedStateActionDirection.enthalpy_j; see its attribute description for units and ordering.
        entropy_j_per_k : float
            Initial value of SolvedStateActionDirection.entropy_j_per_k; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class SolvedStateActionRequest:
    """Observables and ordered directions for one implicit mixed action.

    Notes
    -----
    solved_state_first_actions is the separate batching operation for independent first directions. Availability depends on the selected energy order, reference terms and branch regularity; see ALG-ACTIONS.
    """
    observables: list[SolvedStateObservable]
    """Ordered quantity selections; each result follows this order."""
    directions: list[SolvedStateActionDirection]
    """Ordered physical increments for the mixed contraction, not independent first-order batches."""
    def __init__(self, observables: Sequence[SolvedStateObservable],
                 directions: Sequence[SolvedStateActionDirection]) -> None:
        """Initialize SolvedStateActionRequest from the declared inputs.

        Parameters
        ----------
        observables : Sequence[SolvedStateObservable]
            Initial value of SolvedStateActionRequest.observables; see its attribute description for units and ordering.
        directions : Sequence[SolvedStateActionDirection]
            Initial value of SolvedStateActionRequest.directions; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class SolvedStateActionStatus:
    """Availability or refusal cause for an implicit solved-state action.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    Available: SolvedStateActionStatus
    """Requested quantity is available on the admitted branch."""
    CentralSolveUnavailable: SolvedStateActionStatus
    """No admitted central equilibrium state is available."""
    DefiningEquationsNotSatisfied: SolvedStateActionStatus
    """The central defining-row residual fails action admission."""
    PhaseBoundaryNotDifferentiable: SolvedStateActionStatus
    """A vanished-phase boundary does not admit the requested local derivative."""
    NonuniqueCoalescedPhases: SolvedStateActionStatus
    """Coincident phase labels make the requested derivative nonunique."""
    RankDeficient: SolvedStateActionStatus
    """The defining Jacobian lacks full required rank."""
    NonfiniteJacobian: SolvedStateActionStatus
    """A required defining-Jacobian entry is nonfinite."""
    LinearSolveResidual: SolvedStateActionStatus
    """The implicit linear solve fails its residual criterion."""
    SupportChange: SolvedStateActionStatus
    """Requested derivative leaves the admitted composition support or meets a singular insertion."""
    IonicInsertion: SolvedStateActionStatus
    """Requested charged insertion is unavailable on an ion-free support."""
    EnergyOrderUnavailable: SolvedStateActionStatus
    """The requested action needs an energy order not supplied by this consumer."""
    DensityUnavailable: SolvedStateActionStatus
    """A required pressure-closed density is unavailable."""
    PressureSingular: SolvedStateActionStatus
    """The pressure-density slope does not admit the requested implicit action."""
    AssociationSingular: SolvedStateActionStatus
    """The implicit association-site Jacobian is singular."""
    AssociationIllConditioned: SolvedStateActionStatus
    """Association-site conditioning fails the declared criterion."""
    AssociationResidual: SolvedStateActionStatus
    """The association residual does not meet its admission criterion."""
    MissingThermochemistry: SolvedStateActionStatus
    """Required species thermochemistry is absent."""
    ReferenceUnavailable: SolvedStateActionStatus
    """A consumed source-reference value or derivative is unavailable."""
    OutsideInterval: SolvedStateActionStatus
    """Temperature lies outside a required source interval."""
    ThermalBoundary: SolvedStateActionStatus
    """Requested thermal derivative crosses an unqualified source-interval endpoint."""
    CentralPropertyUnavailable: SolvedStateActionStatus
    """The central observable is unavailable."""
    Nonfinite: SolvedStateActionStatus
    """A required value or derivative is nonfinite."""
    IllConditioned: SolvedStateActionStatus
    """The defining system fails the declared conditioning criterion."""

class SolvedStateActionResult:
    """Central observable and admitted mixed implicit contraction."""
    observable: SolvedStateObservable
    """Requested quantity and phase/component selection."""
    unit: str
    """Physical unit shared by value and action; path variables are dimensionless."""
    order: int
    """Number of contracted directions, a derivative-order count."""
    value: Optional[float]
    """Central observable in unit, or None if unavailable; may survive an unavailable outer action."""
    action: Optional[float]
    """Mixed contraction in unit, or None with the status cause."""
    status: SolvedStateActionStatus
    """Availability or specific refusal cause."""
    property_diagnostics: Optional[PropertyActionDiagnostics]
    """Consumed homogeneous-property work/conditioning, or None when none is retained."""
    diagnostic_message: Optional[str]
    """Specific native action/reference refusal explanation, or None without one."""

class SolvedStateActionDiagnostics:
    """Defining-equation conditioning and shared implicit work."""
    equation_count: int
    """Number of defining equations."""
    state_dimension: int
    """Number of full compiled state coordinates."""
    maximum_energy_order: int
    """Maximum required Helmholtz derivative order."""
    state_factorizations: int
    """Number of shared defining-Jacobian factorizations."""
    linear_solves: int
    """Number of implicit linear solves."""
    rank: Optional[int]
    """Numerical defining-Jacobian rank, or None without an admitted factorization."""
    condition_estimate: Optional[float]
    """Dimensionless condition estimate, or None when unavailable."""
    maximum_linear_residual: Optional[float]
    """Maximum dimensionless normalized implicit linear residual, or None without a solve."""
    central_residual: Optional[float]
    """Dimensionless central defining-row residual, or None without a retained evaluation."""
    central_floor: Optional[float]
    """Dimensionless estimated central residual floor, or None without an estimate."""

class SolvedStateActionBatch:
    """Ordered solved-state actions sharing one prepared center."""
    results: list[SolvedStateActionResult]
    """Per-observable results in request order, each with its own cause and units."""
    diagnostics: SolvedStateActionDiagnostics
    """Shared defining-equation work and conditioning."""

def solved_state_actions(compiled: CompiledSystem, central: EquilibriumResult,
                         request: SolvedStateActionRequest) -> SolvedStateActionBatch:
    """Contract implicit equilibrium derivatives on one admitted branch.

    Parameters
    ----------
    compiled : CompiledSystem
        Admitted native compiled system; its HelmholtzModel must remain alive.
    central : EquilibriumResult
        Admitted result from the same compiled declaration.
    request : SolvedStateActionRequest
        Ordered observables and physical mixed directions.

    Returns
    -------
    SolvedStateActionBatch
        SolvedStateActionBatch in observable order with per-action availability.

    Notes
    -----
    Contractions use dimensionless path variables and retain observable units. This operation does not reinterpret mixed directions as independent first-action batches. See ALG-ACTIONS.
    """
    ...
def solved_state_first_actions(compiled: CompiledSystem, central: EquilibriumResult,
                               observables: Sequence[SolvedStateObservable],
                               directions: Sequence[SolvedStateActionDirection]) -> list[SolvedStateActionBatch]:
    """Evaluate independent first directions at one shared solved center.

    Parameters
    ----------
    compiled : CompiledSystem
        Admitted native compiled system; its HelmholtzModel must remain alive.
    central : EquilibriumResult
        Admitted result from the same compiled declaration.
    observables : Sequence[SolvedStateObservable]
        Ordered solved-state quantity selections.
    directions : Sequence[SolvedStateActionDirection]
        Independent physical first directions; these are separate columns, not one mixed contraction.

    Returns
    -------
    list[SolvedStateActionBatch]
        One SolvedStateActionBatch per direction in input order.
    """
    ...
def linear_solver() -> str:
    """Read the compiled native equilibrium linear-solver selection.

    Returns
    -------
    str
        Linear solver name; a result.linear_solver separately records an actual solve.
    """
    ...
def compile_problem(model: HelmholtzModel | IdealMoleFraction, problem: EquilibriumProblem) -> CompiledSystem:
    """Admit and compile the physical equilibrium declaration.

    Parameters
    ----------
    model : HelmholtzModel
        Native Helmholtz model with declared parameters, active coordinates and optional thermochemistry.
    problem : EquilibriumProblem
        Native physical declaration with ordered phases, feed, reactions, solids and constraints.

    Returns
    -------
    CompiledSystem
        CompiledSystem with rows, native chart and structural rank evidence.

    Notes
    -----
    Invalid support, units, charge, rank or constraint combinations are refused before solving. Keep model alive while the compiled object is used.
    """
    ...
def equilibrium_asl_record(model: Any, problem: Any, continuation: ContinuationOptions, feed_offset: Sequence[float],
                           feed_matrix: Sequence[Sequence[float]], anchor_controls: Sequence[float],
                           identity_component: int, minimum_gap: float) -> str:
    """Serialize the admitted two-liquid equilibrium and affine feed map for ASL.

    Parameters
    ----------
    model : Any
        Model declaration exposing its source parameter document.
    problem : Any
        Fixed-T/P two-free-liquid declaration with explicit supports and absolute molar feed.
    continuation : ContinuationOptions
        Declared continuation controls for the adapter.
    feed_offset : Sequence[float]
        Length n of absolute feed offsets in mol.
    feed_matrix : Sequence[Sequence[float]]
        Row-major (n, k) affine coefficients mapping controls into feed mol; units are mol per caller control unit.
    anchor_controls : Sequence[float]
        Length k of nonnegative anchor controls in the caller-declared units.
    identity_component : int
        Component index used to preserve two-phase identity.
    minimum_gap : float
        Positive dimensionless mole-fraction separation required for that identity.

    Returns
    -------
    str
        Native record string for the equilibrium ASL adapter.

    Notes
    -----
    The map n_feed = feed_offset + feed_matrix*controls must preserve charge and reproduce the declared anchor feed. Source-reference and unsupported reaction forms are refused. Serialization is not a solve.
    """
    ...
def property_asl_record(model: Any) -> str:
    """Serialize a declared model for native ASL property callbacks.

    Parameters
    ----------
    model : Any
        Python model/declaration exposing the source parameter document, thermochemistry and numerical options consumed by the native record writer.

    Returns
    -------
    str
        Native record string containing the parameter and thermochemistry declarations.

    Notes
    -----
    Serialization supplies no state prediction. The native ASL adapter owns subsequent values and exact derivatives.
    """
    ...
def sensitivity(compiled: CompiledSystem, point: Sequence[float], lever: Sequence[float]) -> list[float]:
    """Solve the native implicit first-sensitivity equations.

    Parameters
    ----------
    compiled : CompiledSystem
        Admitted native compiled system; its HelmholtzModel must remain alive.
    point : Sequence[float]
        Full native point in compiled chart order, length d including multipliers; not a flattened physical composition.
    lever : Sequence[float]
        Length 2+n+a of physical increments in [T (K), P (Pa), feed (mol), active parameter] order, optionally extended by [da, db, dc, dd] per reaction and extensive H (J), S (J/K) increments.

    Returns
    -------
    list[float]
        Length d of state-coordinate actions, or an empty list when the defining solve is singular or not admitted.

    Notes
    -----
    Solves J dy = -F_controls dcontrols using exact native columns. An empty result is unavailable sensitivity, not a zero derivative.
    """
    ...
def solve_equilibrium(
    model: HelmholtzModel | IdealMoleFraction,
    problem: EquilibriumProblem,
    start: Optional[EquilibriumResult] = ...,
) -> EquilibriumResult:
    """Solve the declared equilibrium problem natively.

    Parameters
    ----------
    model : HelmholtzModel
        Native Helmholtz model with declared parameters, active coordinates and optional thermochemistry.
    problem : EquilibriumProblem
        Native physical declaration; phase order owns result arrays.
    start : Optional[EquilibriumResult]
        Optional compatible admitted result used only as a seed.

    Returns
    -------
    EquilibriumResult
        EquilibriumResult with admitted physical state or retained refusal diagnostics.

    Notes
    -----
    Uses the native equilibrium driver, exact derivatives and declared numerical settings. A returned success is local numerical admission; global stability requires separate evidence. See ALG-EQUILIBRIUM.
    """
    ...

class ResidualEnergyAction:
    action: float | None
    status: PropertyActionStatus
    diagnostics: PropertyActionDiagnostics

class InfiniteDilutionReference:
    """Analytic ion reference, with solvent closure and association diagnostics (ALG-REFERENCE)."""
    solver_options: list[dict[str, Any]]
    parameter_derivative: list[tuple[ResidualEnergyAction, ResidualEnergyAction]]
    """Ion log-fugacity partials in inverse native parameter units, independently unavailable."""
    log_fugacity: tuple[float, float]
    molality: float
    """Zero mol/kg at the infinite-dilution limit."""
    convergence_error: float | None
    """None: no terminal sequence is measured; inner diagnostics are separate."""
    pure_solvent_log_fugacity: float | None
    pure_solvent_state: State | None
    solvent_state: State | None
    diagnostics: list[PropertyActionDiagnostics]
    solvent_status: PropertyActionStatus
    solvent_cause: str
    solvent_convergence: Convergence | None

class SolventActivity:
    """Same-temperature/pressure solvent activity and practical osmotic coefficient.

    Notes
    -----
    The osmotic coefficient uses mol/kg solvent and the pure-solvent reference at the same T/P. Preserve None and the precision-floor Stop value rather than substituting zero; see practical_osmotic_coefficient in engine/docs/equations.md.
    """
    solvent_activity: Optional[float]
    """Dimensionless solvent activity against the pure-solvent liquid reference, or None with
    status/cause.
    """
    osmotic_coefficient: Optional[Union[float, Stop]]
    """Dimensionless practical molality osmotic coefficient; None when unavailable, or
    Stop.converged_to_floor when precision is unresolved below the declared molality floor.
    """
    solvent_activity_parameter_derivative: list[Optional[float]]
    """Length a of first activity partials in inverse parameter units, with None per unavailable
    column.
    """
    osmotic_parameter_derivative: list[Optional[float]]
    """Length a of first osmotic-coefficient partials in inverse parameter units; columns also
    require an empty osmotic_cause.
    """
    parameter_status: list[PropertyActionStatus]
    """Length a of activity-column availability causes."""
    parameter_cause: list[str]
    """Length a of activity-column explanations."""
    status: PropertyActionStatus
    """Solvent-activity availability or refusal cause."""
    cause: str
    """Solvent-activity refusal explanation, possibly empty."""
    osmotic_cause: str
    """Independent osmotic-coefficient refusal explanation, possibly empty."""
    state: State
    """Retained mixture state, or the pure-solvent state for a pure-solvent input,
    when evaluated. After an early refusal this may be a default State or carry
    only density-convergence evidence; inspect status and cause before consuming
    its physical fields."""

class HenryConstant:
    """Infinite-dilution Henry pressure constant on the mole-fraction scale.

    Notes
    -----
    H_i^px = R*T*rho0*exp(raw_log_reference), with rho0 = 1 mol/m^3. See henry_constant_px in engine/docs/equations.md.
    """
    component_index: int
    """Solute component index in component_ids order."""
    raw_log_reference: Optional[float]
    """Dimensionless raw reference contraction on the native density standard, or None when
    unavailable.
    """
    henry_constant_px: Optional[float]
    """H_i^px = lim(f_i/x_i) in Pa, or None unless the reference is available."""
    observed_terminal_change: Optional[float]
    """Dimensionless change between the terminal reference extrapolations, or None without the
    comparison.
    """
    status: NeutralReferenceStatus
    """Reference availability or refusal cause."""
    message: str
    """Native reference explanation, possibly empty."""

def henry_constants(model: HelmholtzModel, temperature: float, pressure: float,
                    solvent_fractions: Sequence[float], solute_indices: Sequence[int]) -> list[HenryConstant]:
    """Evaluate dilute-solute Henry pressure constants.

    Parameters
    ----------
    model : HelmholtzModel
        Native Helmholtz model with declared parameters, active coordinates and optional thermochemistry.
    temperature : float
        Positive absolute temperature in K.
    pressure : float
        Positive prescribed pressure in Pa.
    solvent_fractions : Sequence[float]
        Length n in component_ids order, with exactly one pure neutral solvent fraction equal to 1 and all other entries zero.
    solute_indices : Sequence[int]
        Nonempty sequence of distinct known solute component indices excluding the solvent; native output is sorted by component index.

    Returns
    -------
    list[HenryConstant]
        One HenryConstant per requested solute, returned in ascending component-index order.

    Notes
    -----
    See henry_constant_px and ALG-REFERENCE. Constants are in Pa at the supplied T/P with one pure-solvent liquid reference. All model components must be neutral; a charged model returns NonNeutralDeclaration results.
    """
    ...

def solvent_reference(model: HelmholtzModel, temperature: float, pressure: float,
    solvent_fractions: Sequence[float]) -> InfiniteDilutionReference:
    """Build the same-T/P pure-solvent reference used for solvent activity.

    Parameters
    ----------
    model : HelmholtzModel
        Native Helmholtz model with declared parameters, active coordinates and optional thermochemistry.
    temperature : float
        Positive absolute temperature in K.
    pressure : float
        Positive prescribed pressure in Pa.
    solvent_fractions : Sequence[float]
        Length n of normalized, charge-neutral solvent mole fractions in component_ids order; no formula ions are inserted here.

    Returns
    -------
    InfiniteDilutionReference
        InfiniteDilutionReference carrying the pure-solvent state and its independent availability.

    Notes
    -----
    The native solvent-reference consumer requires one pure neutral solvent; it does not construct an ionic gauge or infer a mixed-solvent osmotic convention.
    """
    ...
def solvent_activity(model: HelmholtzModel, temperature: float, pressure: float,
    composition: Sequence[float], reference: InfiniteDilutionReference) -> SolventActivity:
    """Evaluate solvent activity and the practical osmotic coefficient.

    Parameters
    ----------
    model : HelmholtzModel
        Native Helmholtz model with declared parameters, active coordinates and optional thermochemistry.
    temperature : float
        Positive absolute temperature in K.
    pressure : float
        Positive prescribed pressure in Pa.
    composition : Sequence[float]
        Length n, dimensionless mole fractions in component_ids order; normalized and charge neutral.
    reference : InfiniteDilutionReference
        Reference carrying an admitted pure-solvent state for the same model, T/P and numerical settings; solvent_reference is sufficient and no ionic limit is required.

    Returns
    -------
    SolventActivity
        SolventActivity with a finite-composition state and independent quantity/column causes.

    Notes
    -----
    Requires an admitted pure-solvent reference at the same T/P. Precision-floor refusal of the osmotic coefficient does not remove other available properties.
    """
    ...

class ElectrolyteActivity:
    mean_ionic_parameter_derivative: list[ResidualEnergyAction]
    """One charge-neutral salt activity on an infinite-dilution convention.

    Notes
    -----
    The formula molality is mol of salt formula units per kg of declared solvent. Explicit ion amounts follow stoichiometry and charge balance. See ALG-REFERENCE in engine/docs/science/algorithms.md.
    """
    solver_options: list[dict[str, Any]]
    """Effective numerical settings for the calculation."""
    single_ion: tuple[float, float]
    """Two dimensionless mole-fraction-scale coefficients in (cation, anion) order; their
    convention is the supplied reference, not independently measurable single-ion activities.
    """
    mean_ionic_mole_fraction_scale: float
    """Stoichiometric mean ionic activity coefficient on the mole-fraction scale, dimensionless."""
    mean_ionic_log_molality_scale: float
    """Natural logarithm of the mean ionic molality-scale activity coefficient, dimensionless."""
    mean_ionic_molality_scale: float
    """Mean ionic activity coefficient on the mol/kg-solvent scale, dimensionless."""
    osmotic_coefficient: Optional[Union[float, Stop]]
    """Dimensionless practical molality osmotic coefficient; None without an admitted solvent
    reference, or Stop.converged_to_floor when unresolved below the declared precision floor.
    """

def infinite_dilution_reference(model: HelmholtzModel, temperature: float, pressure: float, solvent_fractions: Sequence[float], ion_indices: tuple[int, int], stoichiometry: tuple[int, int]) -> InfiniteDilutionReference:
    """Follow the declared charge-neutral salt dilution path.

    Parameters
    ----------
    model : HelmholtzModel
        Native Helmholtz model with declared parameters, active coordinates and optional thermochemistry.
    temperature : float
        Positive absolute temperature in K.
    pressure : float
        Positive prescribed pressure in Pa.
    solvent_fractions : Sequence[float]
        Length n of normalized, charge-neutral solvent mole fractions in component_ids order; no formula ions are inserted here.
    ion_indices : tuple[int, int]
        Two component indices in (cation, anion) order.
    stoichiometry : tuple[int, int]
        Two positive formula-unit ion counts in (cation, anion) order; the resulting formula must be charge neutral.

    Returns
    -------
    InfiniteDilutionReference
        InfiniteDilutionReference in (cation, anion) order, including terminal change and per-parameter availability.

    Notes
    -----
    The retained molality is finite and source-convention specific; terminal agreement is numerical evidence, not a universal limit-error bound. See ALG-REFERENCE.
    """
    ...
def electrolyte_activity(model: HelmholtzModel, temperature: float, pressure: float, solvent_fractions: Sequence[float], ion_indices: tuple[int, int], stoichiometry: tuple[int, int], molality: float, reference: InfiniteDilutionReference) -> ElectrolyteActivity:
    """Evaluate salt activity at a declared formula molality.

    Parameters
    ----------
    model : HelmholtzModel
        Native Helmholtz model with declared parameters, active coordinates and optional thermochemistry.
    temperature : float
        Positive absolute temperature in K.
    pressure : float
        Positive prescribed pressure in Pa.
    solvent_fractions : Sequence[float]
        Length n of normalized, charge-neutral solvent mole fractions in component_ids order; no formula ions are inserted here.
    ion_indices : tuple[int, int]
        Two component indices in (cation, anion) order.
    stoichiometry : tuple[int, int]
        Two positive formula-unit ion counts in (cation, anion) order; the resulting formula must be charge neutral.
    molality : float
        Finite positive formula-unit concentration in mol/kg declared solvent; zero is rejected.
    reference : InfiniteDilutionReference
        Infinite-dilution reference for the same model, T/P, solvent, ions and stoichiometry.

    Returns
    -------
    ElectrolyteActivity
        ElectrolyteActivity with dimensionless ionic coefficients and a separately available osmotic coefficient.

    Notes
    -----
    Ion mole fractions are formed natively from formula stoichiometry and solvent mass. This operation requires a matching admitted reference; it does not change source parameter validity.
    """
    ...

class ObservationKind:
    """Measured or source-reported quantity supplied to native parameter estimation.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.

    See ALG-FIT and engine/docs/regression/science/pure-saturation-regression.md for observation algebra, bases and admission. Parameter fitting is calibration, not independent physical validation.
    """
    SaturationCappedAdditionEnthalpy: ObservationKind

    MeanIonicActivityCoefficient: ObservationKind
    """Dimensionless mean ionic activity coefficient on the declared molality scale."""
    OsmoticCoefficient: ObservationKind
    """Dimensionless single-solvent practical osmotic coefficient of the formula salt."""
    MixtureOsmoticCoefficient: ObservationKind
    """Dimensionless practical osmotic coefficient -ln(a_w)/(M_w*sum(m_j)) for a full
    electroneutral liquid composition and one solvent.
    """
    SolvationGibbsEnergy: ObservationKind
    """Single-ion infinite-dilution solvation Gibbs energy, J/mol on the native reference
    convention.
    """
    TransferGibbsEnergy: ObservationKind
    """Single-ion transfer Gibbs energy, solvent B minus reference solvent A, J/mol on the matched
    convention.
    """
    BulkRelativePermittivity: ObservationKind
    """Dimensionless bulk relative permittivity, or a dimensionless two-state ratio when a
    reference composition is declared.
    """
    LogEquilibriumConstant: ObservationKind
    """Natural logarithm of the named source K at fixed problem temperature, dimensionless; no
    equilibrium solve or reference transfer.
    """
    IsobaricHeatCapacity: ObservationKind
    """Total homogeneous molar Cp at fixed composition, J/(mol K)."""
    EquilibriumIsobaricHeatCapacity: ObservationKind
    """Reacting equilibrium isobaric heat capacity per feed mole, J/(mol-feed K)."""
    Density: ObservationKind
    """Pressure-closed molar density in mol/m^3, or mass density in kg/m^3 with mass_basis."""
    PhaseDensity: ObservationKind
    """Selected equilibrium phase density in mol/m^3, or kg/m^3 with mass_basis."""
    EquilibriumPressure: ObservationKind
    """Solved bubble/dew/saturation pressure in Pa."""
    EquilibriumTemperature: ObservationKind
    """Solved bubble/dew/saturation temperature in K."""
    PhaseMoleFraction: ObservationKind
    """Sum of selected component mole fractions in one phase, dimensionless."""
    PhaseFraction: ObservationKind
    """Selected phase amount divided by the sum of solved fluid-phase amounts, dimensionless; not
    divided by feed_total.
    """
    PartialPressure: ObservationKind
    """Solved system pressure times the selected phase mole-fraction sum, Pa."""
    PhaseAmountRatio: ObservationKind
    """Selected numerator component amount divided by selected denominator amount, dimensionless;
    denominator must be nonzero.
    """
    AdditionEnthalpy: ObservationKind
    """Endpoint extensive-enthalpy difference less added-feed enthalpy, divided by added feed mol,
    J/mol of total added feed; negative is exothermic.
    """

class ResidualForm:
    """Dimensionless residual formed from a prediction and its observation.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    Difference: ResidualForm
    """(prediction-observed)/scale, with positive scale in the observation unit."""
    LogRatio: ResidualForm
    """ln(prediction/observed); both arguments must be positive, and scale is unused."""

class ObservationStatus:
    """Availability or refusal cause of an observation row.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    NotEvaluated: ObservationStatus
    """This row was skipped or evaluation stopped before reaching it."""
    ReferenceUnavailable: ObservationStatus
    """A required reference value or limiting column is unavailable."""
    ThermochemistryUnavailable: ObservationStatus
    """Required caloric data are absent, outside their interval or nonfinite."""
    Available: ObservationStatus
    """Requested prediction or derivative is available."""
    DensityUnavailable: ObservationStatus
    """The declared pressure-closed density was refused."""
    CentralSolveUnavailable: ObservationStatus
    """No admitted central equilibrium state is available."""
    TopologyEvent: ObservationStatus
    """A phase disappeared or phases coalesced, so the declared observation is unavailable."""
    MissingMolarMass: ObservationStatus
    """Mass conversion requires an occupied-species molar mass that is missing or nonpositive."""
    NonpositiveLogArgument: ObservationStatus
    """A logarithmic prediction/residual requires a positive argument."""
    ActionUnavailable: ObservationStatus
    """A required exact native action was refused."""
    Nonfinite: ObservationStatus
    """A required prediction, residual or column is nonfinite."""
    InadmissibleTrial: ObservationStatus
    """Native admission or evaluation raised at this trial; the original explanation is retained."""

class Observation:
    """Inputs, reported value and residual convention for one fitting row.

    Notes
    -----
    Only the fields consumed by kind define its observation. A failed prediction is unavailable with a cause, never a successful zero residual. Endpoint reaction/reference identities and component balances must match for AdditionEnthalpy.
    """
    basis_amount: float
    loading_cap: float
    loading_bounds: tuple[float, float]
    loading_grid: list[float]

    ion_indices: tuple[int, int]
    """Two component indices in (cation, anion) order for the formula-salt reference."""
    stoichiometry: tuple[int, int]
    """Two positive ion counts per formula unit in (cation, anion) order, with zero net formula
    charge.
    """
    molality: float
    finite_solute_molalities: list[float] | None
    """Formula-unit salt concentration in mol/kg solvent for formula-salt observations."""
    reference_composition: list[float]
    """Optional length-n reference solvent composition, dimensionless in component order; used for
    transfer A or a relative-permittivity denominator.
    """
    reaction: str
    """Named source reaction for LogEquilibriumConstant, matched within problem.reactions."""
    kind: ObservationKind
    """Observation quantity, unit and native consumer; see ObservationKind."""
    temperature: float
    """Pressure-closed or feed-state temperature in K; equilibrium observations use
    problem.temperature instead.
    """
    pressure: float
    """Pressure-closed or added-feed-state pressure in Pa; equilibrium observations use
    problem.pressure instead.
    """
    composition: list[float]
    """Length n, normalized dimensionless composition in component order; electrolyte kinds
    interpret it as the declared solvent or full mixture as specified by kind.
    """
    branch: DensityBranch
    """Density branch for a pressure-closed state or added-feed state."""
    problem: EquilibriumProblem
    """Declared equilibrium specification, or initial endpoint for AdditionEnthalpy; also owns the
    named source law for ln K observations.
    """
    final_problem: EquilibriumProblem
    """Final equilibrium endpoint for AdditionEnthalpy; unused by single-endpoint kinds."""
    feed_amount: float
    """Positive added-feed amount in mol for AdditionEnthalpy."""
    phase: int
    """Selected declared fluid-phase index."""
    species: list[int]
    """Nonempty selected component-index set within the phase support; mole-fraction and
    partial-pressure kinds sum the set.
    """
    denominator_phase: int
    """Denominator fluid-phase index for PhaseAmountRatio."""
    denominator_species: list[int]
    """Nonempty denominator component-index set for PhaseAmountRatio."""
    mass_basis: bool
    """For Density/PhaseDensity, request kg/m^3 rather than mol/m^3 through declared molar masses."""
    observed: float
    """Reported measured/source value in the selected prediction unit; not a residual."""
    form: ResidualForm
    """Dimensionless residual transformation."""
    scale: float
    """Positive Difference nondimensionalization scale in the observation unit; not an uncertainty
    unless explicitly declared as such by the measurement model.
    """

class Prediction:
    """One native observation prediction with physical units and failure evidence."""
    saturation: Optional[SaturationSearch]
    saturation_loading: Optional[float]
    endpoint_loading: Optional[float]
    added_amount: Optional[float]
    cap_branch: str
    saturation_loading_gradient: list[float]
    endpoint_loading_gradient: list[float]

    solver_options: list[dict[str, Any]]
    """Effective named numerical settings and provenance."""
    value: Optional[float]
    """Predicted physical quantity in unit, or None with status/message."""
    unit: str
    """Physical unit of the predicted observation."""
    status: ObservationStatus
    """Availability or specific refusal cause."""
    message: str
    """Native explanation, including the failing state/reference when relevant."""
    residual: Optional[float]
    """Unweighted dimensionless residual in the declared ResidualForm, or None when unavailable."""
    convergence: Optional[Convergence]
    """Consumed equilibrium/density convergence record, or None without one."""
    branch: Optional[DensityBranch]
    """Resolved density branch, or None without pressure-closure branch evidence; Unspecified is
    not a phase classification.
    """
    stable_root_count: Optional[int]
    """Mechanically stable root count if established, otherwise None."""
    requested_tolerance_met: Optional[bool]
    """Whether a retained central solve meets the caller criterion, or None without that evidence."""

    first_numerical_failure: Optional[Convergence]
    """First producer refusal, including recovered causes."""

class PredictionGradient:
    saturation_loading_gradient: list[float]
    endpoint_loading_gradient: list[float]
    parameter_actions: list[ResidualEnergyAction]
    """Exact first partials of one physical prediction.

    Notes
    -----
    These are prediction columns, not residual or weighted columns. A missing derivative is not a numerical zero.
    """
    value: list[float]
    """Physical prediction partials in model active-parameter order, with requested
    reaction-coordinate columns where the consumer supplies them; unit divided by each native
    coefficient unit. Empty on refusal.
    """
    status: ObservationStatus
    """Availability or blocking native action cause."""
    message: str
    """Derivative refusal explanation, possibly empty."""

    first_numerical_failure: Optional[Convergence]
    """First native action refusal, including a successful retry."""

class ReactionCoefficient:
    """Source reaction coordinate and its physical unit.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    A: ReactionCoefficient
    """Natural-log constant coefficient a, dimensionless."""
    B: ReactionCoefficient
    """Reciprocal-temperature coefficient b, K."""
    C: ReactionCoefficient
    """Log-temperature coefficient c, dimensionless."""
    D: ReactionCoefficient
    """Linear-temperature coefficient d, 1/K."""
    EnthalpyShift: ReactionCoefficient
    """Source reaction enthalpy increment in J/mol, with a required pivot temperature; preserves ln
    K at the pivot.
    """

class ReactionCoordinate:
    """Named source coefficient or enthalpy-shift coordinate.

    Notes
    -----
    Plain assignments precede an enthalpy-shift addition on a fresh observation copy. Duplicate assignments and two shifts on one reaction are refused. Sequence starts and minimax do not admit EnthalpyShift; see ALG-FIT.
    """
    name: str
    """Reaction name shared by every matching observation/source law."""
    coefficient: ReactionCoefficient
    """Selected a/b/c/d coefficient or pivoted enthalpy increment."""
    pivot_temperature: float
    """Positive pivot in K within the source interval, required for EnthalpyShift; otherwise
    unused.
    """

class CoordinateTarget:
    """One native EOS target of an affine optimizer coordinate."""
    family: ParameterFamily
    """Native parameter family, physical meaning and unit."""
    indices: list[int]
    """Native storage slots tied to the target; same flattened-slot convention as ActiveParameter."""
    coefficient: float
    """Multiplier mapping base coordinate q to this target; units are target units divided by q
    units, so cross-family targets need not share a unit.
    """

class Coordinate:
    """Affine physical coordinate q = origin + scale*v.

    Notes
    -----
    q uses the reaction-coefficient unit or the units declared by its EOS target map. Optimizer v is dimensionless. Native admission owns tied-slot consistency, target conflicts and valid bounds; see affine_parameter_coordinate_map in engine/docs/equations.md.
    """
    targets: list[CoordinateTarget]
    """Ordered EOS targets, each equal to coefficient*q; empty for a reaction coordinate."""
    reaction: Optional[ReactionCoordinate]
    """Optional named source coordinate; mutually exclusive with EOS targets."""
    origin: float
    """Physical base value of q at dimensionless optimizer v=0, in q units."""
    scale: float
    """Positive finite affine scale in q units per dimensionless optimizer increment."""
    lower: float
    """Inclusive lower physical bound on q, in q units; may be negative infinity when undeclared."""
    upper: float
    """Inclusive upper physical bound on q, in q units; may be positive infinity when undeclared."""

class CoordinateSet:
    """Ordered affine parameter coordinates for native estimation.

    Notes
    -----
    Physical rows are EOS targets in coordinate/target order, then held Born diameters where evaluated, then reaction rows. ActiveParameter lists contain EOS rows only; reaction coordinates are separate.
    """
    coordinates: list[Coordinate]
    """Optimizer coordinates in declaration order; this order owns all optimizer columns."""
    def active(self) -> list[ActiveParameter]:
        """Read EOS physical derivative directions from the coordinate declarations.

        Returns
        -------
        list[ActiveParameter]
            One ActiveParameter per declared EOS target, in coordinate/target order; excludes reaction coordinates and unevaluated holds.
        """
        ...
    def tie(self) -> list[int]:
        """Read the physical-row to optimizer-coordinate membership.

        Returns
        -------
        list[int]
            Optimizer-coordinate index per physical EOS/reaction row in native physical-row order; held rows arise only during evaluation.
        """
        ...
    def start(self, parameters: Parameters, observations: list[Observation] = ..., thermochemistry: Optional[ThermochemistryRecord] = ...) -> list[float]:
        """Read dimensionless optimizer starts from current physical parameters and source laws.

        Parameters
        ----------
        parameters : Parameters
            Admitted native EOS parameters.
        observations : list[Observation]
            Optional observation/source declarations needed to identify named reaction coordinates.
        thermochemistry : Optional[ThermochemistryRecord]
            Optional species thermal records in component order; absence leaves ideal/total consumers unavailable.

        Returns
        -------
        list[float]
            Length u optimizer vector in coordinate declaration order.

        Notes
        -----
        Tied slots must agree; invalid or conflicting reaction sources are refused. This reads a start and performs no fit.
        """
        ...

class FitControls:
    """Native estimator tolerances, budgets and state-group workers.

    Notes
    -----
    Least squares uses Ceres with exact native columns. Minimax uses its explicit request and native Ipopt controls; no failed check silently selects it. See ALG-FIT.
    """
    solver_options: list[dict[str, Any]]
    """Effective named numerical settings, including criterion provenance."""
    function_tolerance: float
    """Dimensionless Ceres relative cost-change tolerance; zero selects the named source default."""
    gradient_tolerance: float
    """Ceres projected-gradient criterion in normalized optimizer/objective coordinates; zero
    selects its named default.
    """
    parameter_tolerance: float
    """Dimensionless Ceres optimizer-step criterion; zero selects its named default."""
    maximum_iterations: int
    """Required positive iteration limit, a count."""
    maximum_elapsed_time_seconds: float
    """Required positive elapsed-time budget in seconds."""
    workers: int
    """Requested native state-group process count including the caller; zero selects hardware
    concurrency, bounded by eligible groups. Minimax/point readers and active trace capture use
    one process.
    """
    def __init__(self, *, maximum_iterations: int, maximum_elapsed_time_seconds: float, function_tolerance: float = 0.0,
                 gradient_tolerance: float = 0.0, parameter_tolerance: float = 0.0, workers: int = 0, solver_options: Mapping[str, float] | None = None) -> None:
        """Initialize FitControls from the declared inputs.

        Parameters
        ----------
        maximum_iterations : int
            Initial value of FitControls.maximum_iterations; see its attribute description for units and ordering.
        maximum_elapsed_time_seconds : float
            Initial value of FitControls.maximum_elapsed_time_seconds; see its attribute description for units and ordering.
        function_tolerance : float
            Initial value of FitControls.function_tolerance; see its attribute description for units and ordering.
        gradient_tolerance : float
            Initial value of FitControls.gradient_tolerance; see its attribute description for units and ordering.
        parameter_tolerance : float
            Initial value of FitControls.parameter_tolerance; see its attribute description for units and ordering.
        workers : int
            Initial value of FitControls.workers; see its attribute description for units and ordering.
        solver_options : Mapping[str, float] | None
            Initial value of FitControls.solver_options; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class TrialFailure:
    parameter_actions: list[ResidualEnergyAction]
    """Retained refusal at one optimizer trial."""
    trial: int
    """Native trial counter."""
    observation: int
    """Zero-based failing observation row; the observation count itself denotes a request/model
    exception without a specific row.
    """
    optimizer: list[float]
    """Dimensionless optimizer point, including active measurement corrections where declared."""
    status: ObservationStatus
    """Typed native refusal cause."""
    message: str
    """Original native explanation; not replaced by an apparent zero residual."""
    jacobian: bool
    """Whether the refusal occurred while derivative columns were requested."""

    first_numerical_failure: Optional[Convergence]
    """First cause in this trial's serial row prefix."""

class FitCovariance:
    """Local least-squares covariance and rank evidence.

    Notes
    -----
    The unscaled matrix is the inverse local weighted normal matrix; apply variance_factor only under the stated unknown-common-variance assumption. No covariance is claimed for unfinished, rank-deficient, active-bound or minimax estimates. Model discrepancy and physical validation are outside this local estimate.
    """
    left_vectors: list[float]
    """Optional flattened row-major thin left singular vectors (M, min(M,U)) of the complete
    objective Jacobian; objective_rows labels rows. Empty unless requested and evaluated.
    """
    right_vectors: list[float]
    """Optional flattened row-major thin right singular vectors (U, min(M,U)); U includes active
    input corrections, and columns pair with singular_values. Empty unless requested and
    evaluated.
    """
    projected_gradient_norm: float
    """Maximum absolute projected objective-gradient component in native normalized optimizer
    coordinates; may be nonfinite when unavailable.
    """
    status: str
    """Covariance availability: available, not_converged, rank_deficient, active_bound, unavailable
    or minimax_unavailable.
    """
    assumption: str
    """Explicit local linearization/error assumptions; these limit interpretation of any returned
    covariance.
    """
    rank: int
    """Numerical rank of the complete weighted/whitened objective Jacobian, a count."""
    singular_values: list[float]
    """Descending singular values of the complete objective Jacobian, length min(M,U), in
    normalized residual/optimizer units.
    """
    objective_rows: list[int]
    """Positive-weight observation row indices followed by observation_count+input_column labels
    for measurement penalties.
    """
    active_bounds: list[int]
    """Indices of complete optimizer coordinates lying within the native bound-activity distance."""
    optimizer: list[float]
    """Flattened symmetric (u, u) model-coordinate covariance block after eliminating input
    corrections; native column-major storage (symmetry gives the same flat order as row-major).
    Empty unless available.
    """
    physical: list[float]
    """Flattened symmetric (r, r) covariance D*optimizer*D^T in physical-row order; units are
    products of physical-row units. Empty unless available.
    """
    variance_factor: Optional[float]
    """Optional sum(weight*r^2)/(objective_row_count-rank) for an unknown common residual variance;
    absent with known measurement sigmas/covariance or insufficient degrees of freedom.
    """

class MinimaxRequest:
    """Explicit constrained estimate of the largest absolute decadic target error.

    Notes
    -----
    This selects the native constrained estimator explicitly. It supplies exact first columns and limited-memory curvature, not a parameter Hessian or global optimum certificate. See ALG-FIT, Guarded minimax.
    """
    targets: list[int]
    """Unique observation-row indices with positive log-ratio arguments; their maximum absolute
    log10 error is minimized.
    """
    group_scales: list[float]
    """One positive constraint scale per check group, in that group's error unit; report-only
    groups remain omitted from trial constraints.
    """
    q: float
    """Positive epigraph scale in log10 decades; the minimized dimensionless objective is t/q."""
    tolerance: float
    """Positive dimensionless tolerance for scaled feasibility and first-order stationarity
    diagnostics.
    """
    activity_distance: float
    """Positive dimensionless normalized distance used to identify active variable bounds."""
    margin: float
    """Positive dimensionless initial interior margin; epigraph starts at target maximum + margin*q
    and slacks at absolute error + margin*group scale.
    """
    def __init__(self) -> None:
        """Initialize MinimaxRequest from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...
class MinimaxRow:
    """One scaled inequality of the native minimax problem."""
    kind: str
    """Native inequality role: target, maximum, link or mean."""
    row: int
    """Observation index for a scalar error row; unused for a mean-group row."""
    group: int
    """Check-group index; target rows use the number of check groups as a sentinel."""
    metric: CheckMetric
    """Absolute-difference or absolute-log10-ratio error convention."""
    sign: float
    """Dimensionless +1 or -1 for paired signed error inequalities; zero where unused."""
    threshold: float
    """Frozen physical error threshold in the metric's units."""
    scale: float
    """Positive physical row scale in the same units as the inequality."""
    auxiliary: int | None
    """Epigraph/slack point index subtracted by this row, or None for no such auxiliary."""
    members: list[int]
    """Point indices of mean-group slacks; empty for a scalar row."""
class MinimaxResult:
    """Candidate inequalities, multipliers and local first-order diagnostics.

    Notes
    -----
    Availability, optimizer termination, feasibility, stationarity and final check acceptance are distinct. Interpret post-solve scalar diagnostics, physical_values, multiplier conversions, bound distances/activity and hold/root diagnostics only when availability is 'available'. Point-only evaluation sets 'values_only': point, bounds, scaled values, Jacobian, row identities and frozen checks are retained, while post-solve diagnostics keep uncomputed defaults. Its actual epigraph is point[len(FitResult.optimizer)], not the default t field. Local infeasibility is not proof of global infeasibility. Point evaluation supplies no optimization certificate; minimax supplies no covariance. See ALG-FIT.
    """
    availability: str
    """Whether complete candidate values and derivative diagnostics were evaluated; not an
    optimality claim.
    """
    application_status: str
    """Native application-level termination or point-reader outcome."""
    solver_status: str
    """Ipopt optimizer termination name; separate from feasibility and stationarity."""
    linear_solver: str
    """Actual native Ipopt linear-solver name."""
    feasible: bool
    """Whether recomputed scaled primal/bound violations meet the declared tolerance."""
    stationary: bool
    """Whether stationarity, dual-sign and complementarity residuals meet the declared tolerance;
    separate from feasibility.
    """
    t: float
    """Nonnegative physical epigraph in log10 decades."""
    target_maximum: float
    """Maximum actual absolute target log10 error, in decades."""
    gap: float
    """t minus target_maximum, in decades."""
    primal: float
    """Maximum scaled inequality or bound violation, dimensionless."""
    stationarity: float
    """Maximum absolute stationarity residual of the t/q objective in the native point coordinates."""
    dual_sign: float
    """Maximum wrong-sign inequality/bound multiplier violation for the scaled objective."""
    complementarity: float
    """Maximum scaled multiplier-times-constraint or bound-distance residual."""
    point: list[float]
    """Full point: dimensionless model optimizer coordinates, physical t in decades, then physical
    group slacks in declaration order.
    """
    lower: list[float]
    """Per-point lower bounds in each point entry's units."""
    upper: list[float]
    """Per-point upper bounds in each point entry's units; may include infinity."""
    values: list[float]
    """Length m of dimensionless scaled inequalities c=C/scale, rows order, feasible when <=0."""
    physical_values: list[float]
    """Length m of unscaled inequalities C in each row's error unit."""
    jacobian: list[float]
    """Flattened row-major (m, v) scaled-inequality Jacobian, v=len(point); column units are
    inverse point-entry units.
    """
    multipliers: list[float]
    """Scaled-objective inequality multipliers lambda, rows order."""
    physical_multipliers: list[float]
    """Physical-objective inequality multipliers q*lambda/row.scale, decades per row-error unit."""
    z_lower: list[float]
    """Lower-bound multipliers for t/q, point order, inverse point-entry units."""
    z_upper: list[float]
    """Upper-bound multipliers for t/q, point order, inverse point-entry units."""
    physical_z_lower: list[float]
    """Physical-objective lower-bound multipliers q*z_lower, decades per point-entry unit."""
    physical_z_upper: list[float]
    """Physical-objective upper-bound multipliers q*z_upper, decades per point-entry unit."""
    lower_distances: list[float]
    """Dimensionless point-minus-lower distances divided by the declared variable scales."""
    upper_distances: list[float]
    """Dimensionless upper-minus-point distances divided by the declared variable scales."""
    hold_multipliers: list[float]
    """Reconstructed physical-objective multiplier per hydration hold, decades/(J/mol); nonfinite
    for an unavailable full-space reconstruction.
    """
    full_stationarity: list[float]
    """Physical-objective stationarity in independent optimizer order followed by held-diameter
    rows; units are decades per corresponding optimizer/diameter unit.
    """
    root_errors: list[float]
    """Local dimensionless inequality error estimate from held-root diameter uncertainty, one per
    inequality; not an enclosure.
    """
    slack_groups: list[int]
    """Check-group index per mean-error slack, in slack declaration order."""
    slack_rows: list[int]
    """Observation-row index per mean-error slack, same order as slack_groups."""
    active_lower: list[int]
    """Point indices near their lower bounds under activity_distance."""
    active_upper: list[int]
    """Point indices near their upper bounds under activity_distance."""
    hold_causes: list[str]
    """Full-space multiplier availability explanation per hold."""
    rows: list[MinimaxRow]
    """Inequality identities/scales in values and Jacobian row order."""
    frozen_checks: list[CheckResult]
    """Check reports fixing baseline thresholds before optimization."""
    ipopt_diagnostics: Optional[IpoptDiagnostics]
    """Actual invocation summary; absent for a point reader."""

class MeasuredInputKind:
    """Physical input whose reported uncertainty is declared.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    Temperature: MeasuredInputKind
    """Absolute temperature, K."""
    Pressure: MeasuredInputKind
    """Pressure, Pa."""
    MoleFraction: MeasuredInputKind
    """Independent component mole fraction, dimensionless."""
    Molality: MeasuredInputKind
    """Formula-unit molality, mol/kg solvent."""
class MeasuredInput:
    """One bounded physical input and its known standard uncertainty."""
    def __init__(self) -> None:
        """Initialize MeasuredInput from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...
    kind: MeasuredInputKind
    """Input quantity and its required physical unit."""
    species: int
    """Component index for MoleFraction; unused for other input kinds."""
    sigma: float
    """Nonnegative known standard deviation in unit; zero eliminates this correction coordinate."""
    lower: float
    """Inclusive lower physical input bound in unit, not a standardized correction."""
    upper: float
    """Inclusive upper physical input bound in unit."""
    unit: str
    """Explicit input unit admitted for kind; no basis is inferred from a bare magnitude."""
class MeasurementGroup:
    """Rows sharing the same measured physical inputs.

    Notes
    -----
    Unsupported observation/input combinations, overlapping groups and nonphysical input boxes are refused. See declared-measurement-uncertainty in the native parameter-estimation formulation.
    """
    def __init__(self) -> None:
        """Initialize MeasurementGroup from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...
    rows: list[int]
    """Observation-row indices sharing this explicitly declared measurement; equal state values
    alone do not establish sharing.
    """
    inputs: list[MeasuredInput]
    """Measured inputs in declaration order. A positive-sigma input becomes a standardized correction column only when the group contains at least one positive-weight output.
    """
    dependent_species: list[int]
    """Component indices reconstructed by the admitted normalization/charge chart; not independent
    measurement coordinates.
    """
    feed_total_mol: float
    """Required positive total feed mol for an equilibrium feed chart; fixes the extensive basis
    while fractions move.
    """
class MeasurementUncertainty:
    """Known output/input Gaussian measurement-error declaration.

    Notes
    -----
    The matrix must be finite, symmetric under the native criterion and positive definite; correlations have unit diagonal. Native admission never repairs, jitters or pseudo-inverts it. Whitened residuals use one retained Cholesky factor. See errors_in_variables_objective in engine/docs/equations.md.
    """
    def __init__(self) -> None:
        """Initialize MeasurementUncertainty from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...
    groups: list[MeasurementGroup]
    """Explicit shared-input measurement groups in declaration order."""
    output_sigma: list[float]
    """Length-N output standard deviations in observation order and each observation's physical unit. Positive-weight outputs require finite positive values; excluded zero-weight entries are unused. Difference residuals and unit positive relative weights are required.
    """
    covariance: list[float] | None
    """Optional flattened row-major (M, M) known physical covariance: positive-weight outputs in observation order, then active input correction columns in group/input order. Entry units are products of the corresponding physical units.
    """
    correlation: list[float] | None
    """Optional flattened row-major (M, M) dimensionless known correlation in the same order;
    mutually exclusive with covariance.
    """
class InputReconciliation:
    """Candidate physical input corrections and their objective contribution.

    Notes
    -----
    Reconciled inputs and parameters are fitted candidates. Final checks at the original measured inputs retain a distinct meaning from reconciled-input predictions.
    """
    values: list[float]
    """Reconciled positive-sigma input values in group/input-column order, with each
    MeasuredInput.unit.
    """
    corrections: list[float]
    """Dimensionless standardized physical corrections (value-measured)/sigma, in input-column
    order; these remain unwhitened under a declared matrix.
    """
    penalties: list[float]
    """Squared standardized input residuals, or squared whitened input components with a matrix, in
    input-column order.
    """
    jacobian: list[float]
    """Flattened row-major (observation_count+input_count, complete_optimizer_count) standardized
    objective Jacobian; whitened under a declared matrix.
    """
    row_groups: list[int]
    """Measurement-group identity per output row followed by each input penalty row; native
    sentinel denotes an ungrouped output.
    """
    column_groups: list[int]
    """Measurement-group index per active positive-sigma input column."""
    column_inputs: list[int]
    """Within-group MeasuredInput index per active input column."""
    output_cost: float
    """Half the squared standardized/whitened output residual norm, dimensionless."""
    input_cost: float
    """Half the squared input residual norm, dimensionless; with joint covariance it is conditional
    on the outputs.
    """
class FitResult:
    """Native parameter-estimation candidate, objective and local evidence.

    Notes
    -----
    Only positive-weight objective rows enter least-squares trial values/rank/variance. Zero-weight rows receive final reporting and checks. Fitted parameters are candidates; neither convergence nor covariance proves independent physical validation. See ALG-FIT.
    """
    saturation_predictions: list[Optional[Prediction]]

    first_numerical_failure: Optional[Convergence]
    """First chronological producer refusal, retained after successful recovery."""
    rank_records: list[NumericalRank]
    """Numerical determination of the complete weighted optimizer Jacobian."""
    solver_options: list[dict[str, Any]]
    """Effective named numerical settings and provenance."""
    reconciliation: InputReconciliation | None
    """Candidate input corrections and standardized objective evidence, or None without a declared
    measurement model.
    """
    minimax: MinimaxResult | None
    """Constrained-estimator evidence, or None for least squares."""
    chart: EvaluatedCoordinates
    """Evaluated physical parameter chart, holds and completed reaction laws at the candidate."""
    checks: list[CheckResult]
    """Final native check reports at stored candidate parameters, without resolving a hydration
    hold.
    """
    acceptance: str
    """Native check acceptance classification, separate from optimizer termination and usable."""
    optimizer: list[float]
    """Dimensionless model coordinates in declaration order; fitted input corrections are published
    separately in reconciliation.corrections.
    """
    status: str
    """Native estimator termination classification; a finite candidate can be usable after reaching
    a limit.
    """
    message: str
    """Estimator conclusion and failure explanation."""
    usable: bool
    """Whether the native estimator marks the retained candidate usable. Least squares requires a published candidate with finite physical values and a termination other than failed; a limit may still be usable. This does not establish convergence, check acceptance or physical validity.
    """
    physical: list[float]
    """Physical rows in chart order: EOS targets, held Born diameters, then reaction coordinates.
    The first len(chart.active) entries match chart.active; units vary by row.
    """
    parameters: Parameters
    """Stored candidate native EOS parameters, including resolved held diameters."""
    predictions: list[float]
    """Length N of physical predictions in observation order, each ObservationKind unit; inspect
    statuses for unavailable entries.
    """
    residuals: list[float]
    """Length N of unweighted dimensionless residuals in the declared form; statuses retain
    missing-row causes.
    """
    weighted_residuals: list[float]
    """Length N of least-squares output objective residuals sqrt(weight)*residual, whitened when
    declared; input penalties are published separately in reconciliation. Excluded rows are zero
    and unavailable rows are not successful zeros.
    """
    physical_jacobian: list[float]
    """Flattened row-major (N, r) physical prediction partials in observation/physical-row order;
    row units divided by physical-coordinate units.
    """
    optimizer_jacobian: list[float]
    """Least squares: row-major weighted/whitened output residual columns (N, u) for model
    coordinates only; complete input/penalty columns are in reconciliation.jacobian. Minimax:
    raw prediction columns in model optimizer order, not Ceres weighted columns.
    """
    statuses: list[ObservationStatus]
    """One observation availability cause per row, including NotEvaluated after skipped/failed
    evaluation.
    """
    failures: list[TrialFailure]
    """Rejected native trials with original cause and optimizer point."""
    iterations: int
    """Estimator iteration count."""
    residual_evaluations: int
    """Number of residual evaluations."""
    jacobian_evaluations: int
    """Number of exact Jacobian evaluations."""
    trial_failures: int
    """Number of rejected native trial evaluations."""
    iteration_seconds: list[float]
    """Per-Ceres-iteration elapsed timing in seconds; entry zero is the starting iteration."""
    iteration_costs: list[float]
    """Per-Ceres-iteration dimensionless half squared weighted residual norm; entry zero is the
    start.
    """
    initial_cost: float
    """Least-squares starting cost 0.5*sum(weight*r^2), extended for measurement errors; NaN for
    minimax. Failed evaluations may leave zero, which is meaningful only with status/usable.
    """
    final_cost: float
    """Least-squares final cost in the same convention as initial_cost; NaN for minimax. A zero
    after refusal is not evidence of a perfect fit.
    """
    covariance: FitCovariance
    """Local least-squares uncertainty/rank evidence with explicit availability and assumptions;
    unavailable for minimax.
    """

def predict(model: HelmholtzModel, observation: Observation) -> Prediction:
    """Evaluate one declared physical observation natively.

    Parameters
    ----------
    model : HelmholtzModel
        Native Helmholtz model with declared parameters, active coordinates and optional thermochemistry.
    observation : Observation
        Reported quantity, physical inputs and declared residual convention.

    Returns
    -------
    Prediction
        Prediction with physical unit, unweighted residual and independent refusal evidence.

    Notes
    -----
    Native validation precedes state/reference/equilibrium evaluation. A failed physical state is unavailable, not a zero prediction or residual.
    """
    ...
def differentiate(model: HelmholtzModel, observation: Observation) -> PredictionGradient:
    """Evaluate exact first columns of one declared physical prediction.

    Parameters
    ----------
    model : HelmholtzModel
        Native Helmholtz model with declared parameters, active coordinates and optional thermochemistry.
    observation : Observation
        Physical row to evaluate and differentiate.

    Returns
    -------
    PredictionGradient
        PredictionGradient in active native coordinate order, or an empty gradient with cause.

    Notes
    -----
    Columns are prediction partials, not weighted residual columns. The native property/reference/solved-state owner supplies them.
    """
    ...
def fit(parameters: Parameters, thermochemistry: Optional[ThermochemistryRecord], coordinates: Sequence[Coordinate],
        observations: Sequence[Observation], weights: Sequence[float], controls: FitControls, holds: Sequence[HydrationHold] = ...,
        checks: Sequence[CheckGroup] = ..., baseline: Parameters | None = ..., singular_vectors: bool = ...,
        starts: Optional[Sequence[Sequence[float]]] = None, iterations: Optional[int] = None,
        seconds: Optional[float] = None, minimax: MinimaxRequest | None = ...,
        uncertainty: MeasurementUncertainty | None = ...) -> Union[FitResult, list[FitResult]]:
    """Estimate bounded parameters with native least squares or explicitly requested minimax.

    Parameters
    ----------
    parameters : Parameters
        Admitted native EOS parameters.
    thermochemistry : Optional[ThermochemistryRecord]
        Optional species thermal records in component order; absence leaves ideal/total consumers unavailable.
    coordinates : Sequence[Coordinate]
        Ordered affine physical parameter/source coordinates.
    observations : Sequence[Observation]
        Ordered reported quantities with physical inputs, source conventions and residual forms.
    weights : Sequence[float]
        Nonnegative relative row weights; empty means all one. Least-squares residuals are sqrt(weight)*r; zero-weight rows are final reports. Known measurement models require unit positive weights.
    controls : FitControls
        Required iteration/elapsed budgets and named native tolerances.
    holds : Sequence[HydrationHold]
        Ordered dependent Born-diameter constraints, resolved by exact native held-row solves.
    checks : Sequence[CheckGroup]
        Declared final error groups; these do not silently change the estimator.
    baseline : Parameters | None
        Optional separately stored source parameters used to freeze baseline check thresholds.
    singular_vectors : bool
        Request local objective-Jacobian singular vectors in the returned covariance evidence.
    starts : Optional[Sequence[Sequence[float]]]
        Optional ordered dimensionless model optimizer starts; None requests a single start. Sequence mode has shared budgets and refuses unsupported coordinate/measurement forms.
    iterations : Optional[int]
        Optional shared sequence iteration budget; None uses controls.maximum_iterations.
    seconds : Optional[float]
        Optional shared sequence wall-time budget in seconds; None uses controls.maximum_elapsed_time_seconds.
    minimax : MinimaxRequest | None
        Explicit constrained-estimation request; None selects Ceres least squares.
    uncertainty : MeasurementUncertainty | None
        Optional known Gaussian physical measurement-error declaration, including shared input corrections and optional covariance.

    Returns
    -------
    Union[FitResult, list[FitResult]]
        One FitResult, or ordered list for explicit sequence starts; each retains its own termination, checks and availability.

    Notes
    -----
    Ceres uses exact native columns. Minimax has distinct local feasibility/stationarity evidence and no covariance. Failed trials retain their original cause; returned candidates remain calibration results, not physical validation. See ALG-FIT.
    """
    ...

class CheckMetric:
    """Error convention for native parameter-estimation checks.

    Notes
    -----
    Native enum name/value properties follow the module convention in epcsaft._core.
    """
    AbsoluteDifference: CheckMetric
    """Absolute prediction-minus-observation difference in the observation unit; the retained
    signed error is prediction-observation.
    """
    AbsoluteLog10Ratio: CheckMetric
    """Absolute log10(prediction/observation), in decades; both arguments must be positive, and
    retained row errors keep their sign.
    """

class HydrationHold:
    """Eliminate one Born diameter by holding a single-ion solvation target.

    Notes
    -----
    Both inner criteria apply. The dependent derivative is formed from the exact held-row slope; a failed hold is not substituted by a fixed guessed diameter. See dependent_born_diameter_derivative in engine/docs/equations.md.
    """
    def __init__(self) -> None:
        """Initialize HydrationHold from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...
    reference: Observation
    """SolvationGibbsEnergy observation defining the target in observed, J/mol."""
    lower: float
    """Lower Born-diameter bracket in angstrom."""
    upper: float
    """Upper Born-diameter bracket in angstrom."""
    energy_tolerance: float
    """Positive target energy-residual criterion in J/mol."""
    diameter_tolerance: float
    """Positive achieved diameter-uncertainty criterion in angstrom."""
    maximum_iterations: int
    """Positive inner safeguarded-Newton iteration limit."""

class HoldResult:
    """Held-diameter residual and exact implicit derivatives."""
    energy_error: float
    """Signed held solvation-energy residual in J/mol."""
    correction: float
    """Magnitude of the last Newton correction abs(H/H_d), angstrom."""
    uncertainty: float
    """Achieved diameter uncertainty min(correction, bracket width), angstrom; a local numerical
    estimate.
    """
    derivative: float
    """Held-row partial H_d in J/(mol angstrom)."""
    parameter_derivative: list[float]
    """Held energy-row partials in EvaluatedCoordinates.active order: independent EOS targets followed by held Born diameters, excluding reaction-coordinate rows. Units are J/mol divided by each active parameter's native unit.
    """
    iterations: int
    """Inner safeguarded-Newton iteration count."""

    rank: NumericalRank
    """Scalar hold derivative rank and its actual cutoff."""

class EvaluatedCoordinates:
    """Mapped candidate parameters, dependent diameters and reaction laws."""
    reaction_laws: dict[str, ReactionLogPolynomial]
    """Completed native ln K(T) laws keyed by reaction name; replay these laws without applying the
    fitted coordinates again.
    """
    parameters: Parameters
    """Native EOS parameters after affine mapping and held-diameter resolution."""
    active: list[ActiveParameter]
    """EOS physical-row ActiveParameters, independent targets followed by held Born diameters;
    reaction rows are not ActiveParameters.
    """
    physical: list[float]
    """Physical target values in EOS-target, held-diameter, then reaction order, with each row's
    native unit.
    """
    map: list[float]
    """Flattened row-major (r, u) exact D=d(physical)/d(model optimizer), including holds and
    target multipliers; row units per dimensionless optimizer increment.
    """
    holds: list[HoldResult]
    """HoldResult per declared hold in input order."""

class CheckGroup:
    """Declared final error group and optional fixed or baseline thresholds.

    Notes
    -----
    Checks re-evaluate stored parameters at the original observations without resolving hydration holds. A failed group is not proof that no feasible fit exists.
    """
    def __init__(self) -> None:
        """Initialize CheckGroup from the declared inputs.

        Returns
        -------
        None
            None.
        """
        ...
    name: str
    """Caller check-group identity."""
    rows: list[int]
    """Observation indices in the group, in declared reporting order."""
    metric: CheckMetric
    """Absolute-difference or absolute-decades error convention."""
    maximum: float | None
    """Optional maximum-error threshold in the metric unit."""
    mean: float | None
    """Optional arithmetic mean-absolute-error threshold in the metric unit."""
    baseline: bool
    """Freeze thresholds from the separately supplied baseline parameter set under the native check
    rules.
    """
    report_only: bool
    """Report failures without changing candidate acceptance; omitted from minimax trial
    constraints.
    """

class CheckResult:
    """Signed row errors and absolute-error group statistics."""
    name: str
    """Declared check-group identity."""
    unit: str
    """Metric unit: physical observation unit for Difference or log10 decades for ratio errors."""
    status: str
    """Native availability and threshold outcome; unavailable rows remain explicit."""
    rows: list[int]
    """Observation indices in reporting order."""
    failed_rows: list[int]
    """Observation indices whose final check evaluation was unavailable or whose absolute error
    exceeded the maximum threshold; a mean-only failure need not identify an individual row.
    """
    errors: list[float]
    """Signed physical differences or signed log10 ratios, corresponding to rows; group statistics
    use absolute values.
    """
    slopes: list[float]
    """Signed-error derivative with respect to prediction: 1 for Difference or
    1/(prediction*ln(10)) for log10 ratio; inverse prediction units for the latter.
    """
    maximum: float
    """Maximum retained absolute error in unit; interpret only with status, since an unavailable
    row need not replace a partial maximum.
    """
    mean: float
    """Arithmetic mean absolute error in unit; interpret only with status and failed_rows."""
    argmax: int
    """Observation-row index attaining the maximum absolute error; meaningful only for an available
    group.
    """
    maximum_threshold: float | None
    """Effective maximum threshold in unit, or None without one."""
    mean_threshold: float | None
    """Effective mean threshold in unit, or None without one."""
    baseline_maximum: float | None
    """Frozen baseline maximum absolute error in unit; native binding returns None when no baseline
    value exists.
    """
    baseline_mean: float | None
    """Frozen baseline mean absolute error in unit; native binding returns None when no baseline
    value exists.
    """
    report_only: bool
    """Whether this group's failure is excluded from acceptance decisions."""

def evaluate_coordinates(parameters: Parameters, set: CoordinateSet, optimizer: Sequence[float], holds: Sequence[HydrationHold], thermochemistry: ThermochemistryRecord | None = ..., solver_options: Mapping[str, float] | None = None) -> EvaluatedCoordinates:
    """Map one optimizer point and solve its declared dependent diameters.

    Parameters
    ----------
    parameters : Parameters
        Admitted native EOS parameters.
    set : CoordinateSet
        Ordered native CoordinateSet.
    optimizer : Sequence[float]
        Length u of dimensionless model coordinates.
    holds : Sequence[HydrationHold]
        Ordered hydration holds, each with a declared bracket and two numerical criteria.
    thermochemistry : ThermochemistryRecord | None
        Optional species thermal records in component order; absence leaves ideal/total consumers unavailable.
    solver_options : Mapping[str, float] | None
        Optional mapping of admitted named numerical settings; None uses source defaults.

    Returns
    -------
    EvaluatedCoordinates
        EvaluatedCoordinates including physical target order, exact map and completed reaction laws.

    Notes
    -----
    The exact map includes target coefficients and implicit held-diameter derivatives; failed inner admission is not replaced by a guessed value.
    """
    ...
def evaluate_checks(stored: Parameters, observations: Sequence[Observation], checks: Sequence[CheckGroup], baseline: Parameters | None = ..., solver_options: Mapping[str, float] | None = None) -> list[CheckResult]:
    """Re-evaluate final check groups at stored physical parameters.

    Parameters
    ----------
    stored : Parameters
        Stored native parameter set; held diameters are consumed as stored, not resolved again.
    observations : Sequence[Observation]
        Original declared physical observations.
    checks : Sequence[CheckGroup]
        Ordered final check groups.
    baseline : Parameters | None
        Optional stored source parameters for native baseline threshold rules.
    solver_options : Mapping[str, float] | None
        Optional mapping of admitted named numerical settings; None uses source defaults.

    Returns
    -------
    list[CheckResult]
        CheckResult list in check-group order, preserving unavailable rows.

    Notes
    -----
    This is a point check and does not fit parameters. Report-only failures are retained independently of acceptance.
    """
    ...

def observations_at(parameters: Parameters, coordinates: Sequence[Coordinate], optimizer: Sequence[float], holds: Sequence[HydrationHold], observations: Sequence[Observation], errors: Sequence[float],
    thermochemistry: Optional[ThermochemistryRecord] = None, solver_options: Mapping[str, float] | None = None) -> list[Observation]:
    """Construct observation values at a declared native parameter point with signed errors.

    Parameters
    ----------
    parameters : Parameters
        Admitted native EOS parameters.
    coordinates : Sequence[Coordinate]
        Ordered affine physical coordinate declarations.
    optimizer : Sequence[float]
        Dimensionless model point in coordinate order.
    holds : Sequence[HydrationHold]
        Ordered native hydration holds.
    observations : Sequence[Observation]
        Template observations supplying physical states and source conventions.
    errors : Sequence[float]
        One signed physical difference or signed log10 error per observation in native check-error convention, not an uncertainty.
    thermochemistry : Optional[ThermochemistryRecord]
        Optional species thermal records in component order; absence leaves ideal/total consumers unavailable.
    solver_options : Mapping[str, float] | None
        Optional mapping of admitted named numerical settings; None uses source defaults.

    Returns
    -------
    list[Observation]
        Fresh Observation list in input order with manufactured reported values.

    Notes
    -----
    This makes numerical recovery cases, not new measurements or physical validation. Native code evaluates the point before changing observation values.
    """
    ...

class ObservationEvaluation:
    """Native point predictions, exact columns and optional measurement evidence.

    Notes
    -----
    For Difference, dr/dv=(dy/dv)/scale; for LogRatio, dr/dv=(dy/dv)/prediction. These point-reader prediction columns must not be confused with weighted Ceres residual columns. See ALG-FIT.
    """
    saturation_predictions: list[Optional[Prediction]]

    first_numerical_failure: Optional[Convergence]
    """First producer refusal in the eligible row prefix."""
    first_numerical_failure_observation: int
    """Source row, meaningful only when a first cause is present."""
    solver_options: list[dict[str, Any]]
    """Effective named numerical settings and provenance."""
    reconciliation: InputReconciliation | None
    """Standardized candidate input corrections and complete objective columns, or None without a
    measurement declaration.
    """
    state_densities: list[list[float]]
    """Ragged per-observation density evidence in mol/m^3 when retained for uncertainty evaluation:
    central state/declared phases first, then reference ladder and pure-solvent states where
    present; otherwise empty.
    """
    reference_values: list[list[float]]
    """Ragged per-observation dimensionless reference ln fugacity coefficients in cation/anion
    order, followed by pure-solvent ln fugacity coefficient when present; otherwise empty.
    """
    predictions: list[float]
    """Length N of physical predictions in observation order; units depend on ObservationKind and
    statuses retain refusals.
    """
    gradient: list[float]
    """Flattened row-major (N, r) exact raw physical prediction partials, observation/physical-row
    order; prediction units per native physical-row unit.
    """
    optimizer_gradient: list[float]
    """Flattened row-major prediction columns in model optimizer order, including physical
    input-motion columns when declared; prediction units per dimensionless optimizer increment.
    """
    residual_optimizer_gradient: list[float]
    """Unit-weight dimensionless residual columns in optimizer order; with uncertainty, the
    complete standardized and optionally whitened objective columns including input penalties.
    """
    chart: EvaluatedCoordinates
    """Evaluated physical coordinate map, held diameters and completed reaction laws."""
    statuses: list[ObservationStatus]
    """One availability or refusal cause per observation in input order."""
def evaluate_observations(parameters: Parameters, coordinates: Sequence[Coordinate], optimizer: Sequence[float],
    holds: Sequence[HydrationHold], observations: Sequence[Observation],
    thermochemistry: Optional[ThermochemistryRecord] = None,
    uncertainty: MeasurementUncertainty | None = ..., solver_options: Mapping[str, float] | None = None) -> ObservationEvaluation:
    """Evaluate declared predictions and exact point columns without optimizing.

    Parameters
    ----------
    parameters : Parameters
        Admitted native EOS parameters.
    coordinates : Sequence[Coordinate]
        Ordered affine physical parameter/source coordinates.
    optimizer : Sequence[float]
        Dimensionless model coordinates, extended by standardized input corrections when declared.
    holds : Sequence[HydrationHold]
        Ordered native hydration holds.
    observations : Sequence[Observation]
        Ordered physical observation declarations.
    thermochemistry : Optional[ThermochemistryRecord]
        Optional species thermal records in component order; absence leaves ideal/total consumers unavailable.
    uncertainty : MeasurementUncertainty | None
        Optional known measurement-error model for shared input corrections and whitening.
    solver_options : Mapping[str, float] | None
        Optional mapping of admitted named numerical settings; None uses source defaults.

    Returns
    -------
    ObservationEvaluation
        ObservationEvaluation with physical and optimizer prediction columns and unit-weight residual columns.

    Notes
    -----
    No optimizer success is inferred from point evaluation. Raw prediction columns and standardized/whitened residual columns have distinct meanings.
    """
    ...
def evaluate_minimax(parameters: Parameters, coordinates: Sequence[Coordinate], optimizer: Sequence[float],
    holds: Sequence[HydrationHold], observations: Sequence[Observation], checks: Sequence[CheckGroup],
    baseline: Parameters | None, minimax: MinimaxRequest, thermochemistry: Optional[ThermochemistryRecord] = None, solver_options: Mapping[str, float] | None = None) -> FitResult:
    """Evaluate the native minimax rows at the declared starting point without optimizing.

    Parameters
    ----------
    parameters : Parameters
        Admitted native EOS parameters.
    coordinates : Sequence[Coordinate]
        Ordered affine physical coordinates supported by minimax.
    optimizer : Sequence[float]
        Dimensionless model optimizer point.
    holds : Sequence[HydrationHold]
        Ordered native hydration holds.
    observations : Sequence[Observation]
        Ordered physical observations.
    checks : Sequence[CheckGroup]
        Declared final/error guard groups.
    baseline : Parameters | None
        Optional source parameter set for frozen check thresholds.
    minimax : MinimaxRequest
        Explicit positive epigraph, row-scale and feasibility controls.
    thermochemistry : Optional[ThermochemistryRecord]
        Optional species thermal records in component order; absence leaves ideal/total consumers unavailable.
    solver_options : Mapping[str, float] | None
        Optional mapping of admitted named numerical settings; None uses source defaults.

    Returns
    -------
    FitResult
        FitResult with minimax row values/columns and frozen checks; no optimization or stationarity certificate.

    Notes
    -----
    The native point reader shares the constrained estimator row algebra and retains its local availability limits.
    """
    ...
def formula_composition(solvent_fractions: Sequence[float], parameters: Parameters, ion_indices: tuple[int, int],
    stoichiometry: tuple[int, int], molality: float) -> list[float]:
    """Convert one charge-neutral salt formula molality into explicit mole fractions.

    Parameters
    ----------
    solvent_fractions : Sequence[float]
        Length n of normalized, charge-neutral solvent mole fractions in component_ids order; no formula ions are inserted here.
    parameters : Parameters
        Admitted native EOS parameters.
    ion_indices : tuple[int, int]
        Two component indices in (cation, anion) order.
    stoichiometry : tuple[int, int]
        Two positive formula-unit ion counts in (cation, anion) order; the resulting formula must be charge neutral.
    molality : float
        Finite positive formula-unit concentration in mol/kg declared solvent; zero is refused.

    Returns
    -------
    list[float]
        Length n normalized dimensionless mole fractions in component_ids order.

    Notes
    -----
    Uses declared solvent molar masses and positive ion counts; formula molality is mol salt/kg solvent. The native owner enforces solvent and formula electroneutrality. See formula_unit_composition in engine/docs/equations.md.
    """
    ...

class EquilibriumStart:
    """One seed for the same thermodynamic declaration.

    Notes
    -----
    Seeds carry no thermodynamic constraints. Native admission rejects incompatible or conflicting seed representations.
    """
    coordinates: list[float]
    """Optional full native compiled coordinate vector, including multipliers; empty selects
    physical seed construction.
    """
    amounts: list[float]
    """Optional fluid-phase amount guesses in mol, phase order; empty leaves them unspecified."""
    compositions: list[list[float]]
    """Optional ragged mole-fraction guesses in phase order. Each nonempty row follows that phase's declared support, with one entry per declared support component; an empty row leaves that phase unspecified.
    """
    packing_fractions: list[Optional[float]]
    """Optional per-phase dimensionless packing guesses; each None leaves that phase unspecified."""
    temperature: Optional[float]
    """Optional positive temperature guess in K for a free temperature."""
    pressure: Optional[float]
    """Optional positive pressure guess in Pa for a free pressure."""
    previous: Optional[EquilibriumResult]
    """Optional previous EquilibriumResult supplying a compatible starting state; not a replacement
    physical declaration.
    """
    def __init__(self, *, coordinates: Sequence[float] = ..., amounts: Sequence[float] = ...,
                 compositions: Sequence[Sequence[float]] = ..., packing_fractions: Sequence[Optional[float]] = ...,
                 temperature: Optional[float] = None, pressure: Optional[float] = None,
                 previous: Optional[EquilibriumResult] = None) -> None:
        """Initialize EquilibriumStart from the declared inputs.

        Parameters
        ----------
        coordinates : Sequence[float]
            Initial value of EquilibriumStart.coordinates; see its attribute description for units and ordering.
        amounts : Sequence[float]
            Initial value of EquilibriumStart.amounts; see its attribute description for units and ordering.
        compositions : Sequence[Sequence[float]]
            Initial value of EquilibriumStart.compositions; see its attribute description for units and ordering.
        packing_fractions : Sequence[Optional[float]]
            Initial value of EquilibriumStart.packing_fractions; see its attribute description for units and ordering.
        temperature : Optional[float]
            Initial value of EquilibriumStart.temperature; see its attribute description for units and ordering.
        pressure : Optional[float]
            Initial value of EquilibriumStart.pressure; see its attribute description for units and ordering.
        previous : Optional[EquilibriumResult]
            Initial value of EquilibriumStart.previous; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class EquilibriumWorkLimits:
    """Shared and per-attempt budgets for declared equilibrium starts.

    Notes
    -----
    None leaves that limit undeclared. Limits are checked by the native driver; a callback/evaluation in progress can overrun the elapsed deadline.
    """
    attempts: Optional[int]
    """Optional maximum attempted-start count."""
    iterations: Optional[int]
    """Optional shared iteration budget across attempts and seed sub-solves."""
    attempt_iterations: Optional[int]
    """Optional maximum iterations within one attempt."""
    seconds: Optional[float]
    """Optional total elapsed-time budget in seconds."""
    attempt_seconds: Optional[float]
    """Optional per-attempt elapsed-time budget in seconds."""
    def __init__(self, *, attempts: Optional[int] = None, iterations: Optional[int] = None,
                 attempt_iterations: Optional[int] = None, seconds: Optional[float] = None,
                 attempt_seconds: Optional[float] = None) -> None:
        """Initialize EquilibriumWorkLimits from the declared inputs.

        Parameters
        ----------
        attempts : Optional[int]
            Initial value of EquilibriumWorkLimits.attempts; see its attribute description for units and ordering.
        iterations : Optional[int]
            Initial value of EquilibriumWorkLimits.iterations; see its attribute description for units and ordering.
        attempt_iterations : Optional[int]
            Initial value of EquilibriumWorkLimits.attempt_iterations; see its attribute description for units and ordering.
        seconds : Optional[float]
            Initial value of EquilibriumWorkLimits.seconds; see its attribute description for units and ordering.
        attempt_seconds : Optional[float]
            Initial value of EquilibriumWorkLimits.attempt_seconds; see its attribute description for units and ordering.

        Returns
        -------
        None
            None.
        """
        ...

class EquilibriumStartsResult:
    """Ordered start outcomes and actual budget use."""
    selected: Optional[int]
    """Index of the first result with success=True and requested_tolerance_met=True, or None without such a result; a retained floor-admitted state need not be selected."""
    results: list[Optional[EquilibriumResult]]
    """One optional result per declared start in input order; None preserves an unrun/refused
    attempt with causes.
    """
    causes: list[str]
    """Explanation per start, including why an attempt was not run."""
    attempt_seconds: list[float]
    """Elapsed seconds per declared start."""
    iterations: int
    """Total consumed iterations across attempts."""
    seconds: float
    """Total elapsed wall time in seconds."""
    overrun_seconds: float
    """Elapsed deadline overrun in seconds, reported rather than hidden."""
    budget_exhausted: bool
    """Whether a declared shared work limit was exhausted."""

def solve_equilibrium_starts(model: HelmholtzModel, problem: EquilibriumProblem, starts: list[EquilibriumStart], limits: EquilibriumWorkLimits) -> EquilibriumStartsResult:
    """Evaluate declared seeds under shared native work limits.

    Parameters
    ----------
    model : HelmholtzModel
        Native Helmholtz model with declared parameters, active coordinates and optional thermochemistry.
    problem : EquilibriumProblem
        One thermodynamic declaration shared by every seed.
    starts : list[EquilibriumStart]
        Ordered native seeds; none selects new physical inputs.
    limits : EquilibriumWorkLimits
        Declared shared and per-attempt work limits.

    Returns
    -------
    EquilibriumStartsResult
        EquilibriumStartsResult preserving input order, unrun causes and selected-result identity.

    Notes
    -----
    All attempts use the same native solve and stopping admission. A finite seed search is not proof of global equilibrium.
    """
    ...

def reaction_log10_reciprocal(law: ReactionLogPolynomial) -> list[float]:
    """Convert a natural-log reaction law to a decadic reciprocal-temperature form.

    Parameters
    ----------
    law : ReactionLogPolynomial
        Native ln K=a+b/T+c*ln(T/Tref)+d*T source correlation.

    Returns
    -------
    list[float]
        Three coefficients [A, B, C] of -log10 K = A/T + B + C*T, in K, dimensionless and 1/K.

    Notes
    -----
    Returns [-b/ln(10), -a/ln(10), -d/ln(10)]. A nonzero logarithmic coefficient c is refused; no refit, standard-state change or validity-range extension is performed.
    """
    ...

class SpeciesActivities:
    standard_molality_mol_per_kg: float | None
    solvent_index: int | None
    log_activities: list[float | None]
    activities: list[float | None]
    causes: list[str]

class ReactionActivities:
    phase_index: int
    reaction_index: int
    basis: str
    cause: str
    reference_pressure_pa: float | None
    log_quotient: float | None
    species: SpeciesActivities

@overload
def reaction_activities(model: HelmholtzModel, result: EquilibriumResult) -> list[ReactionActivities]: ...
@overload
def reaction_activities(model: IdealMoleFraction, result: EquilibriumResult) -> list[ReactionActivities]: ...

class AffineLoading:
    """Declared n(alpha)=initial+alpha*direction, mol; positive loading bounds in mol/mol and target pressure in Pa."""
    initial: list[float]
    direction: list[float]
    lower: float
    upper: float
    seed: float
    pressure: float

class SaturationOperand:
    quantity: Quantity
    """Physical candidate comparison value and its native local numerical estimate."""
    value: float
    error: NumericalEstimate

class SaturationCandidate:
    """Composed state and comparison operands: loading, T, P, then density, packing, amount and all mole fractions per declared phase."""
    state: EquilibriumResult
    operands: list[SaturationOperand]

class SaturationComparison:
    """All pairwise candidate comparison judgments; no transitive chaining of numerical error estimates."""
    left: int
    right: int
    identity: str
    operands: list[EquilibriumComparison]

class SaturationNode:
    """One required grid node: loading in mol/mol, pressure mismatch in Pa, slope in Pa per mol/mol."""
    loading: float
    mismatch: float
    slope: float
    state: EquilibriumResult
    estimates: EquilibriumNumericalEstimates
    action_cause: str
    action_seconds: float
    estimate_seconds: float
    local_estimates: EquilibriumNumericalEstimates
    action_diagnostics: Optional[SolvedStateActionDiagnostics]

class SaturationSearch:
    """Finite grid/search evidence; unique within the declared search, without interval or global certification."""
    grid: list[float]
    nodes: list[SaturationNode]
    events: list[tuple[int, int]]
    candidates: list[SaturationCandidate]
    comparisons: list[SaturationComparison]
    causes: list[str]
    groups: list[list[int]]
    accepted: Optional[int]
    first_numerical_failure: Optional[Convergence]
    brackets: int
    iterations: int
    evaluations: int
    hessian_calls: int
    attempts: int
    complete: bool
    declaration: AffineLoading
    work_at_entry: tuple[int, int, int, int]
    """Iterations, evaluations, Hessian callbacks and driver attempts spent before this search."""
    uniqueness_scope: str
    detection_limits: str
