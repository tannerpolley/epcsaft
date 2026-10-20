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
EOS_STANDARD_STATE_ID: str
"""Native reaction activity identity: (rho*x/(1 mol/m^3))*exp(mu_res/(RT)), dimensionless."""

class ParameterError(ValueError):
    """A parameter value, unit, shape or declaration is inadmissible."""
    ...
class EosError(RuntimeError):
    """A native thermodynamic evaluation cannot supply the requested result."""
    ...
class DensityStructureError(EosError):
    """The density search cannot establish the requested root structure.

    Notes
    -----
    A failed or incomplete scan is not evidence of a unique stable root. See ALG-DENSITY in engine/docs/science/algorithms.md.
    """
    ...
class ConvergenceError(EosError):
    """A native iterative calculation stops without an admitted result.

    Notes
    -----
    The exception message retains the residual, stopping criterion and available convergence diagnostics.
    """
    ...
class DeclarationError(ParameterError):
    """A required physical declaration is absent or inadmissible."""
    ...

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
    @property
    def amplification_over_declared(self) -> Optional[float]:
        """Read estimated amplification relative to its declared factor.

        Returns
        -------
        Optional[float]
            Dimensionless amplification/criterion.factor, or None without either input.
        """
        ...

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
    """True if collecting trace evidence failed; the retained trace may omit observations."""

class Capture:
    """Capture native numerical traces in a context manager.

    Notes
    -----
    Capture observes calculations without changing their stopping criteria. Detailed probes are enabled by diagnostics; exceptions from the enclosed calculation are not suppressed.
    """
    def __init__(self, diagnostics: bool = False) -> None:
        """Initialize Capture from the declared inputs.

        Parameters
        ----------
        diagnostics : bool
            Enable detailed numerical probes in addition to operation spans.

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
    """Convert source ideal-gas coefficients to the Engine gas-constant convention.

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
    NASA coefficients and entropy pressure transfers follow the native thermochemistry owner; this operation does not fit or invent species calorics.
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
    """7 values in source order a0 through a6; temperature is in K. Heat capacity divided by R uses
    powers T**0 through T**4; a5/T and a6 are the enthalpy and entropy integration terms.
    Coefficient units follow these powers; retain the source gas constant until
    prepare_thermochemistry.
    """
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
    """Source equilibrium law ln K = a + b/T + c ln(T/Tref) + d T.

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
    def __init__(self, a: float, b: float, c: float, d: float, reference_temperature: float,
                 temperature_min: float, temperature_max: float, standard_state_id: str) -> None:
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
    """Ion suppression coefficient in kg/mol, component slots."""
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
    """Mole-fraction increments, length n or empty for no composition motion; dimensionless,
    preserving normalization and charge balance.
    """
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
    """Mole-fraction increments, length n or empty for no composition motion; dimensionless,
    preserving normalization and charge balance.
    """
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
    """Dimensionless association factorization pivot ratio."""
    association_residual: float
    """Dimensionless association mass-action residual magnitude."""
    pressure_density_derivative: float
    """Central dP/drho at fixed T and composition, Pa m^3/mol."""
    pressure_residual: float
    """Density closure pressure residual, Pa."""
    stability_determinant: float
    """Dimensionless pressure slope (dP/drho)/(RT); mechanical regularity evidence, not global
    stability.
    """
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

class PropertyActionBatch:
    """Property actions sharing one central state and implicit work."""
    solver_options: list[dict[str, Any]]
    """Effective numerical settings as named dictionaries with values, provenance and units where
    declared.
    """
    results: list[PropertyActionResult]
    """Results in observable order; log-fugacity requests expand into component_ids order."""
    maximum_energy_order: int
    """Maximum Helmholtz derivative order used across the batch."""
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
    """Density-root convergence record, or None without such a solve."""

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
    """Length n in component_ids order, ln(phi_i)=mu_res_i/(RT)-ln(Z), dimensionless. The list is
    None unless Z is positive and finite; individual entries preserve unavailable chemical
    potentials.
    """
    component_property_status: list[PropertyActionStatus]
    """Length n, component_ids order; explains component-property availability, including
    IonicInsertion.
    """
    log_fugacity_parameter_derivative: list[Optional[float]]
    """Flattened row-major (n, a) first partials of ln(phi_i) at fixed T/P/x, inverse native
    parameter units; unavailable entries are None.
    """
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
    """Length n in component_ids order, m^3/mol; None unless requested with positive pressure
    slope; individual unavailable insertions remain None.
    """
    partial_molar_enthalpy: Optional[list[Optional[float]]]
    """Length n in component_ids order, J/mol; None unless requested with positive pressure slope
    and admitted thermochemistry; unavailable component insertions remain None.
    """
    partial_molar_entropy: Optional[list[Optional[float]]]
    """Length n in component_ids order, J/(mol K); None unless requested with positive pressure
    slope and thermochemistry; zero-fraction species and unavailable insertions have None
    entries.
    """

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
    """Sign of dP/drho: positive for mechanical stability, negative for instability, zero at a zero
    slope or when unresolved.
    """
    linked: bool
    """Whether the diagnostic linked this point through its scan/refinement logic; not a global
    phase-stability statement.
    """
    declared: bool
    """Whether this point belongs to the declared search evidence rather than an unresolved
    auxiliary result.
    """
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
    """Refined finite-domain boundary point, or None if no such search was requested or retained."""
    geometric_cap: float
    """Packing-derived upper density bound in mol/m^3; not proof of finite evaluability below it."""
    declared_crossings: int
    """Number of pressure crossings declared by the search."""
    evaluations: int
    """Number of diagnostic EOS evaluations."""
    stop: str
    """Diagnostic termination explanation; an incomplete search does not establish uniqueness."""

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
            Ordered physical increments; mixed actions contract all directions, whereas each first density action is reported separately.

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
    Liquid: PhaseKind
    Vapor: PhaseKind

class AmountKind:
    Free: AmountKind
    Pinned: AmountKind

class FeedKind:
    MoleFractions: FeedKind
    Amounts: FeedKind
    MassFractions: FeedKind
    Molalities: FeedKind

class AmountSpec:
    kind: AmountKind
    fraction: float
    def __init__(self) -> None: ...
    @staticmethod
    def free() -> AmountSpec: ...
    @staticmethod
    def pinned(value: float) -> AmountSpec: ...

class FreeValue:
    free: bool
    value: float
    def __init__(self) -> None: ...
    @staticmethod
    def fixed(value: float) -> FreeValue: ...
    @staticmethod
    def variable(guess: Optional[float] = None) -> FreeValue: ...

class PhaseSpec:
    name: str
    kind: PhaseKind
    support: list[int]
    support_declared: bool
    amount: AmountSpec
    composition_guess: list[float]
    packing_guess: Optional[float]
    def __init__(self) -> None: ...

class FeedSpec:
    kind: FeedKind
    values: list[float]
    solvent: list[int]
    solvent_mass_kg: float
    def __init__(self) -> None: ...

class SpeciesReference:
    InfiniteDilution: SpeciesReference
    PureLiquid: SpeciesReference
    PureGas: SpeciesReference

class NeutralReferenceDeclaration:
    solvent_index: Optional[int]
    solute_molality_ratios: list[float]
    standard_molality_mol_per_kg: float
    species_reference: list[SpeciesReference]
    def __init__(self) -> None: ...

class ReferenceSourceBasis:
    CommonMolalityInfiniteDilution: ReferenceSourceBasis
    RawMoleFractionInfiniteDilution: ReferenceSourceBasis

class NeutralReferenceStatus:
    Available: NeutralReferenceStatus
    InvalidDeclaration: NeutralReferenceStatus
    NonNeutralDeclaration: NeutralReferenceStatus
    UnsupportedSourceBasis: NeutralReferenceStatus
    FreeTemperature: NeutralReferenceStatus
    FreePressure: NeutralReferenceStatus
    DensityUnavailable: NeutralReferenceStatus
    BranchChange: NeutralReferenceStatus
    Nonfinite: NeutralReferenceStatus
    TerminalChangeExceeded: NeutralReferenceStatus

class NeutralReferenceReactionDiagnostic:
    solver_options: list[dict[str, Any]]
    reaction_index: int
    terminal_contraction: Optional[float]
    observed_terminal_change: Optional[float]
    terminal_contraction_pressure_derivative: Optional[float]
    terminal_contraction_pressure_second_derivative: Optional[float]
    reference_pressure_pa: float
    status: NeutralReferenceStatus
    message: str

class ReactionReference:
    source_basis: ReferenceSourceBasis
    reference_pressure_pa: Optional[float]
    def __init__(self, source_basis: ReferenceSourceBasis,
                 reference_pressure_pa: Optional[float] = None) -> None: ...

class ReactionSpec:
    name: str
    stoichiometry: list[float]
    equilibrium_constant: Optional[float]
    affinity_j_per_mol: Optional[float]
    correlation: Optional[ReactionLogPolynomial]
    reference: Optional[ReactionReference]
    def __init__(self) -> None: ...

class SolidSpec:
    name: str
    stoichiometry: list[float]
    equilibrium_constant: Optional[float]
    correlation: Optional[ReactionLogPolynomial]
    reference: Optional[ReactionReference]
    molar_volume_m3_per_mol: float
    amount: AmountSpec
    def __init__(self) -> None: ...

class StabilityOptions:
    evaluate_supplied_trials: bool
    trials: Optional[list[PhaseSpec]]
    declaration_budget: Optional[int]
    trial_budget: Optional[int]
    standard_uncertainty_over_rt: list[float]
    def __init__(self) -> None: ...

class StabilityTrial:
    phase: PhaseSpec
    seed: str
    evidence: str
    message: str
    optimizer: EquilibriumResult
    raw_coordinates: list[float]
    raw_residuals: list[float]
    composition: list[float]
    split_coordinates: list[float]
    basis: list[list[float]]
    row_transformation: list[list[float]]
    quantities: dict[str, float]
    split_checks: list[dict[str, float]]
    donor: Optional[int]
    resolved_descent: bool

class StabilityReport:
    completed: bool
    message: str
    trials: list[StabilityTrial]
    candidate_coordinates: list[float]
    standard_over_rt: list[float]
    plane_over_rt: list[float]
    present_phases: list[dict[str, float]]
    work: dict[str, float]

class DiscoveryRecord:
    declaration: EquilibriumProblem
    candidate: EquilibriumResult
    rounds: list[StabilityReport]
    declarations: list[EquilibriumResult]
    declaration_roles: list[str]
    declaration_attempt_coordinates: list[list[float]]
    work: dict[str, float]


def cooling_onset(compiled: CompiledSystem) -> tuple[EquilibriumResult, list[StabilityReport]]: ...

def stability_test(compiled: CompiledSystem, central: EquilibriumResult,
                   options: StabilityOptions = ...) -> StabilityReport: ...

class CriticalDeclaration:
    density_guess_mol_per_m3: float
    mode_guess: list[float]
    def __init__(self) -> None: ...

class EnvelopeCoordinate:
    LogTemperature: EnvelopeCoordinate
    LogPressure: EnvelopeCoordinate
    LogK: EnvelopeCoordinate

class EnvelopeCondition:
    kind: EnvelopeCoordinate
    component: int
    target: float
    def __init__(self) -> None: ...

class CriticalRecord:
    eigenvalues: list[float]
    mode: list[float]
    cubic: float
    zero_band: float
    condition_estimate: float
    backward_error: float
    rank: int

class EquilibriumProblem:
    @property
    def solver_options(self) -> list[dict[str, Any]]: ...
    @solver_options.setter
    def solver_options(self, value: Mapping[str, float] | None) -> None: ...
    discovery: Optional[StabilityOptions]
    phases: list[PhaseSpec]
    temperature: FreeValue
    pressure: FreeValue
    feed: FeedSpec
    reactions: list[ReactionSpec]
    standard_gibbs_j_per_mol: list[float]
    topology_declared: bool
    row_tolerance: float
    enthalpy_j: float | None
    entropy_j_per_k: float | None
    critical: Optional[CriticalDeclaration]
    envelope: Optional[EnvelopeCondition]
    azeotrope: bool
    neutral_reference: Optional[NeutralReferenceDeclaration]
    solids: list[SolidSpec]
    def __init__(self) -> None: ...

class InventoryBound:
    species: int
    status: str
    termination: str
    maximum: Optional[float]
    feed_normalized_amount: Optional[float]
    log_fraction: Optional[float]
    balance_residual: float
    normalization_residual: float
    lower_bound_residual: float
    witness: list[float]
    amounts_mol: list[float]
    finite_witness: list[float]
    limiting: bool
    supported: bool
    floor: float
    lp_error: float
    finite_residual: float
    feasibility_termination: str

def inventory_bounds(compiled: CompiledSystem, species: Sequence[int], fraction: bool,
                     floor: float) -> list[InventoryBound]: ...

class RowKind:
    Chemical: RowKind
    Pressure: RowKind
    Moiety: RowKind
    Electro: RowKind
    Normalization: RowKind
    Enthalpy: RowKind
    Entropy: RowKind
    CriticalMode: RowKind
    CriticalNormalization: RowKind
    CriticalCubic: RowKind
    EnvelopeCondition: RowKind
    AzeotropeCondition: RowKind
    SolidAffinity: RowKind

class Row:
    kind: RowKind
    phase: int
    index: int
    eqid: str

class PointEvaluation:
    values: list[float]
    objective: float
    floor: float

class CompiledSystem:
    rows: list[Row]
    pinned_count: int
    free_intensive_count: int
    rank: int
    charge_row_count: int
    primal_size: int
    generic_nullity: int
    feed_total: float
    has_ions: bool
    has_reaction: bool
    conserves_total_moles: bool
    jacobian_rank: int
    removed_species: list[int]
    def generic_point(self) -> list[float]: ...
    def evaluate(self, point: Sequence[float]) -> PointEvaluation: ...
    def jacobian(self, point: Sequence[float]) -> list[float]: ...
    def lever_jacobian(self, point: Sequence[float]) -> list[float]: ...
    def objective_gradient(self, point: Sequence[float]) -> list[float]: ...
    def hessian(self, point: Sequence[float], constraint_multipliers: Sequence[float] = ...) -> list[float]: ...

class EquilibriumResult:
    solver_settings: dict[str, float]
    diagnostic_options_file: str
    linear_solver: str
    solver_options: list[dict[str, Any]]
    """Solved phases in declared order (with Discover(): discovery.declaration.phases). Units: K, Pa, mol, mol/m^3.
    Per-phase-per-component lists are phase-major: entry p*n + i is phase p, component i of component_ids
    (n components). State lists are empty and temperature and pressure are None when success is False."""
    discovery: Optional[DiscoveryRecord]
    neutral_reference_diagnostics: list[NeutralReferenceReactionDiagnostic]
    success: bool
    message: str
    solver_status: str
    temperature: Optional[float]  # K; None without a state
    pressure: Optional[float]  # Pa; None without a state
    amounts: list[float]  # phase amounts in mol on the feed's basis (1 mol total for MoleFractions), per phase
    mole_fractions: list[float]  # phase-major, n per phase
    molar_densities: list[float]  # mol/m^3, per phase
    packing_fractions: list[float]  # per phase
    multipliers: list[float]
    constraint_multipliers: list[float]
    residuals: list[float]
    classified_residuals: list[float]
    coordinates: list[float]
    phase_count: int
    iterations: int
    evaluations: int
    hessian_calls: int
    objective: float
    reduced_hessian_minimum: float
    convergence: Convergence
    rows: list[Row]
    unknown_count: int
    jacobian_rank: int
    generic_nullity: int
    requested_tolerance_met: Optional[bool]
    topology_event: TopologyEvent
    topology_phase: Optional[int]
    topology_pair: Optional[tuple[int, int]]
    critical: Optional[CriticalRecord]
    solid_amounts: list[float]  # mol, per declared solid
    solid_affinity_decades: list[float]  # SI_s in decades, per declared solid
    def phase_mole_fractions(self, phase: int) -> list[float]: ...  # one declared phase, component_ids order

class TopologyEvent:
    None_: TopologyEvent
    PhaseDisappeared: TopologyEvent
    PhasesCoalesced: TopologyEvent

class ContinuationStop:
    reached_target: ContinuationStop
    initial_failure: ContinuationStop
    invalid_path: ContinuationStop
    singular_jacobian: ContinuationStop
    branch_lost: ContinuationStop
    phase_disappeared: ContinuationStop
    phases_coalesced: ContinuationStop
    left_domain: ContinuationStop
    infeasible: ContinuationStop
    nonfinite: ContinuationStop
    stalled: ContinuationStop
    step_too_small: ContinuationStop
    budget_exhausted: ContinuationStop
    critical_endpoint: ContinuationStop

class ContinuationOptions:
    """Ordinary path: steps are fractions of the start-to-target path; a rejected step halves and a step never
    grows, so initial_step is the largest step. Envelope mode (Problem.envelope set) ignores initial_step and
    the target problem: an ordinary step is at most envelope_max_log_step (log-coordinate units), the step after
    an accepted one is sized from the measured predictor density error, and a rejected ordinary step (branch-window
    miss, or a corrector that ends on coincident phases) halves; a rejected critical-point crossing, closing step or
    terminal failure ends the trace at once with the states accepted so far. It needs envelope_straddle > 0,
    2 * envelope_straddle <= envelope_max_log_step and T and P guesses through Free(guess). Both modes stop with
    budget_exhausted after more than max_retries_per_step consecutive rejections of retryable steps (reason in
    cause_stop and cause_message) and with step_too_small when the halved step falls below minimum_step."""
    initial_step: float
    minimum_step: float
    max_steps: int
    max_retries_per_step: int
    envelope_max_log_step: float
    envelope_straddle: float
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, *, initial_step: float, minimum_step: float, max_steps: int, max_retries_per_step: int,
                 envelope_max_log_step: float = 0.0, envelope_straddle: float = 0.0) -> None: ...

class ContinuationRequest:
    start_problem: EquilibriumProblem
    target_problem: EquilibriumProblem
    start_result: Optional[EquilibriumResult]
    options: ContinuationOptions
    parameter_targets: Sequence[float]
    def __init__(self) -> None: ...

class ContinuationResult:
    success: bool
    message: str
    stop: ContinuationStop
    path_values: list[float]  # per state: path fraction, or in envelope mode the conditions[k] target (ln T, ln P or ln K)
    parameter_values: list[list[float]]
    states: list[EquilibriumResult]
    reached_path: float
    attempted_path: Optional[float]
    event_phase: Optional[int]
    event_pair: Optional[tuple[int, int]]
    attempts: int
    accepted_steps: int
    retries: int
    last_result: Optional[EquilibriumResult]
    cause_stop: Optional[ContinuationStop]
    cause_convergence: Optional[Convergence]
    cause_message: Optional[str]
    conditions: list[EnvelopeCondition]
    critical_points: list[EquilibriumResult]
    critical_seeds: list[list[float]]

def continue_equilibrium(model: HelmholtzModel, request: ContinuationRequest) -> ContinuationResult: ...

class SolvedStateObservableKind:
    Temperature: SolvedStateObservableKind
    Pressure: SolvedStateObservableKind
    PhaseAmount: SolvedStateObservableKind
    PhaseComponentAmount: SolvedStateObservableKind
    PhaseMoleFraction: SolvedStateObservableKind
    PhaseVolume: SolvedStateObservableKind
    PhaseMolarDensity: SolvedStateObservableKind
    PhasePackingFraction: SolvedStateObservableKind
    PhasePressure: SolvedStateObservableKind
    PhaseResidualChemicalPotentialOverRt: SolvedStateObservableKind
    PhaseProperty: SolvedStateObservableKind
    PhaseLogActivity: SolvedStateObservableKind
    PhaseGibbsMolarVolume: SolvedStateObservableKind
    SystemEnthalpy: SolvedStateObservableKind
    EquilibriumIsobaricHeatCapacity: SolvedStateObservableKind
    SolidAmount: SolvedStateObservableKind
    SolidSaturationIndex: SolvedStateObservableKind

class SolvedStateObservable:
    kind: SolvedStateObservableKind
    phase: int
    component: int
    property: PropertyObservable
    def __init__(self, kind: SolvedStateObservableKind, phase: int = 0, component: int = 0,
                 property: PropertyObservable = ...) -> None: ...

class SolvedStateActionDirection:
    temperature_k: float
    pressure_pa: float
    feed_amounts_mol: list[float]
    active_parameters: list[float]
    reaction_coefficients: list[list[float]]
    enthalpy_j: float
    entropy_j_per_k: float
    def __init__(self, temperature_k: float, pressure_pa: float, feed_amounts_mol: Sequence[float],
                 active_parameters: Sequence[float], reaction_coefficients: Sequence[Sequence[float]] = ...,
                 enthalpy_j: float = ..., entropy_j_per_k: float = ...) -> None: ...

class SolvedStateActionRequest:
    observables: list[SolvedStateObservable]
    directions: list[SolvedStateActionDirection]
    def __init__(self, observables: Sequence[SolvedStateObservable],
                 directions: Sequence[SolvedStateActionDirection]) -> None: ...

class SolvedStateActionStatus:
    Available: SolvedStateActionStatus
    CentralSolveUnavailable: SolvedStateActionStatus
    DefiningEquationsNotSatisfied: SolvedStateActionStatus
    PhaseBoundaryNotDifferentiable: SolvedStateActionStatus
    NonuniqueCoalescedPhases: SolvedStateActionStatus
    RankDeficient: SolvedStateActionStatus
    NonfiniteJacobian: SolvedStateActionStatus
    LinearSolveResidual: SolvedStateActionStatus
    SupportChange: SolvedStateActionStatus
    IonicInsertion: SolvedStateActionStatus
    EnergyOrderUnavailable: SolvedStateActionStatus
    DensityUnavailable: SolvedStateActionStatus
    PressureSingular: SolvedStateActionStatus
    AssociationSingular: SolvedStateActionStatus
    AssociationIllConditioned: SolvedStateActionStatus
    AssociationResidual: SolvedStateActionStatus
    MissingThermochemistry: SolvedStateActionStatus
    ReferenceUnavailable: SolvedStateActionStatus
    OutsideInterval: SolvedStateActionStatus
    ThermalBoundary: SolvedStateActionStatus
    CentralPropertyUnavailable: SolvedStateActionStatus
    Nonfinite: SolvedStateActionStatus
    IllConditioned: SolvedStateActionStatus

class SolvedStateActionResult:
    observable: SolvedStateObservable
    unit: str
    order: int
    value: Optional[float]
    action: Optional[float]
    status: SolvedStateActionStatus
    property_diagnostics: Optional[PropertyActionDiagnostics]
    diagnostic_message: Optional[str]

class SolvedStateActionDiagnostics:
    equation_count: int
    state_dimension: int
    maximum_energy_order: int
    state_factorizations: int
    linear_solves: int
    rank: Optional[int]
    condition_estimate: Optional[float]
    maximum_linear_residual: Optional[float]
    central_residual: Optional[float]
    central_floor: Optional[float]

class SolvedStateActionBatch:
    results: list[SolvedStateActionResult]
    diagnostics: SolvedStateActionDiagnostics

def solved_state_actions(compiled: CompiledSystem, central: EquilibriumResult,
                         request: SolvedStateActionRequest) -> SolvedStateActionBatch: ...
def solved_state_first_actions(compiled: CompiledSystem, central: EquilibriumResult,
                               observables: Sequence[SolvedStateObservable],
                               directions: Sequence[SolvedStateActionDirection]) -> list[SolvedStateActionBatch]: ...
def linear_solver() -> str: ...
def compile_problem(model: HelmholtzModel, problem: EquilibriumProblem) -> CompiledSystem: ...
def equilibrium_asl_record(model: Any, problem: Any, continuation: ContinuationOptions, feed_offset: Sequence[float],
                           feed_matrix: Sequence[Sequence[float]], anchor_controls: Sequence[float],
                           identity_component: int, minimum_gap: float) -> str: ...
def property_asl_record(model: Any) -> str: ...
def sensitivity(compiled: CompiledSystem, point: Sequence[float], lever: Sequence[float]) -> list[float]: ...
def solve_equilibrium(
    model: HelmholtzModel,
    problem: EquilibriumProblem,
    start: Optional[EquilibriumResult] = ...,
) -> EquilibriumResult: ...

class InfiniteDilutionReference:
    """Two-ion infinite-dilution reference at declared temperature and pressure.

    Notes
    -----
    See ALG-REFERENCE and reference_molality_sequence in the canonical Engine documents. Each derivative column retains its own availability; a missing column does not erase the reference value.
    """
    solver_options: list[dict[str, Any]]
    """Effective numerical settings for the reference calculation."""
    parameter_derivative: list[tuple[Optional[float], Optional[float]]]
    """Length a of (cation, anion) reference log-fugacity partials in inverse parameter units;
    individual columns are None when their limiting admission fails.
    """
    parameter_convergence_error: list[float]
    """Length a, maximum two-ion terminal-pair derivative change in inverse parameter units."""
    log_fugacity: tuple[float, float]
    """Dimensionless reference log fugacities in (cation, anion) order, using the native density
    standard.
    """
    molality: float
    """Retained reference molality in mol/kg solvent; this finite rung is not an exact
    zero-concentration state.
    """
    convergence_error: float
    """Dimensionless terminal change of the reference log fugacities."""
    pure_solvent_log_fugacity: Optional[float]
    """Same-T/P pure-solvent log fugacity, dimensionless, or None when unavailable."""
    pure_solvent_state: Optional[State]
    """Same-T/P pure-solvent liquid state, or None when unavailable."""
    solvent_status: PropertyActionStatus
    """Availability or refusal cause of the pure-solvent reference."""
    solvent_cause: str
    """Native pure-solvent refusal explanation, possibly empty."""
    solvent_convergence: Optional[Convergence]
    """Pure-solvent density convergence evidence, or None without a retained record."""

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
    """Finite-composition homogeneous state used for the activity calculation."""

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
        Length n of normalized, charge-neutral solvent mole fractions in component_ids order; no formula ions are inserted here.
    solute_indices : Sequence[int]
        Ordered solute component indices; each result retains its own reference availability.

    Returns
    -------
    list[HenryConstant]
        One HenryConstant per solute_indices entry, preserving requested order.

    Notes
    -----
    See henry_constant_px and ALG-REFERENCE. Each constant is in Pa, with solvent composition and liquid reference held at the supplied T/P.
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
        Infinite-dilution reference for the same model, T/P, solvent, ions and stoichiometry.

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
        Nonnegative formula-unit concentration in mol/kg declared solvent.
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
    MeanIonicActivityCoefficient: ObservationKind
    OsmoticCoefficient: ObservationKind
    MixtureOsmoticCoefficient: ObservationKind
    SolvationGibbsEnergy: ObservationKind
    TransferGibbsEnergy: ObservationKind
    BulkRelativePermittivity: ObservationKind
    LogEquilibriumConstant: ObservationKind
    IsobaricHeatCapacity: ObservationKind
    EquilibriumIsobaricHeatCapacity: ObservationKind
    Density: ObservationKind
    PhaseDensity: ObservationKind
    EquilibriumPressure: ObservationKind
    EquilibriumTemperature: ObservationKind
    PhaseMoleFraction: ObservationKind
    PhaseFraction: ObservationKind
    PartialPressure: ObservationKind
    PhaseAmountRatio: ObservationKind
    AdditionEnthalpy: ObservationKind

class ResidualForm:
    Difference: ResidualForm
    LogRatio: ResidualForm

class ObservationStatus:
    NotEvaluated: ObservationStatus
    ReferenceUnavailable: ObservationStatus
    ThermochemistryUnavailable: ObservationStatus
    Available: ObservationStatus
    DensityUnavailable: ObservationStatus
    CentralSolveUnavailable: ObservationStatus
    TopologyEvent: ObservationStatus
    MissingMolarMass: ObservationStatus
    NonpositiveLogArgument: ObservationStatus
    ActionUnavailable: ObservationStatus
    Nonfinite: ObservationStatus
    InadmissibleTrial: ObservationStatus

class Observation:
    ion_indices: tuple[int, int]
    stoichiometry: tuple[int, int]
    molality: float
    reference_composition: list[float]
    reaction: str
    kind: ObservationKind
    temperature: float
    pressure: float
    composition: list[float]
    branch: DensityBranch
    problem: EquilibriumProblem
    final_problem: EquilibriumProblem
    feed_amount: float
    phase: int
    species: list[int]
    denominator_phase: int
    denominator_species: list[int]
    mass_basis: bool
    observed: float
    form: ResidualForm
    scale: float

class Prediction:
    solver_options: list[dict[str, Any]]
    value: Optional[float]
    unit: str
    status: ObservationStatus
    message: str
    residual: Optional[float]
    convergence: Optional[Convergence]
    branch: Optional[DensityBranch]
    stable_root_count: Optional[int]
    requested_tolerance_met: Optional[bool]

class PredictionGradient:
    value: list[float]
    status: ObservationStatus
    message: str

class ReactionCoefficient:
    A: ReactionCoefficient
    B: ReactionCoefficient
    C: ReactionCoefficient
    D: ReactionCoefficient
    EnthalpyShift: ReactionCoefficient

class ReactionCoordinate:
    name: str
    coefficient: ReactionCoefficient
    pivot_temperature: float

class CoordinateTarget:
    family: ParameterFamily
    indices: list[int]
    coefficient: float

class Coordinate:
    targets: list[CoordinateTarget]
    reaction: Optional[ReactionCoordinate]
    origin: float
    scale: float
    lower: float
    upper: float

class CoordinateSet:
    coordinates: list[Coordinate]
    def active(self) -> list[ActiveParameter]: ...
    def tie(self) -> list[int]: ...
    def start(self, parameters: Parameters, observations: list[Observation] = ..., thermochemistry: Optional[ThermochemistryRecord] = ...) -> list[float]: ...

class FitControls:
    solver_options: list[dict[str, Any]]
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float
    maximum_iterations: int
    maximum_elapsed_time_seconds: float
    workers: int
    def __init__(self, *, maximum_iterations: int, maximum_elapsed_time_seconds: float, function_tolerance: float = 0.0,
                 gradient_tolerance: float = 0.0, parameter_tolerance: float = 0.0, workers: int = 0, solver_options: Mapping[str, float] | None = None) -> None: ...

class TrialFailure:
    trial: int
    observation: int
    optimizer: list[float]
    status: ObservationStatus
    message: str
    jacobian: bool

class FitCovariance:
    left_vectors: list[float]
    right_vectors: list[float]
    projected_gradient_norm: float
    status: str
    assumption: str
    rank: int
    singular_values: list[float]
    objective_rows: list[int]
    active_bounds: list[int]
    optimizer: list[float]
    physical: list[float]
    variance_factor: Optional[float]

class MinimaxRequest:
    targets: list[int]
    group_scales: list[float]
    q: float
    tolerance: float
    activity_distance: float
    margin: float
    def __init__(self) -> None: ...
class MinimaxRow:
    kind: str
    row: int
    group: int
    metric: CheckMetric
    sign: float
    threshold: float
    scale: float
    auxiliary: int | None
    members: list[int]
class MinimaxResult:
    availability: str
    application_status: str
    solver_status: str
    linear_solver: str
    feasible: bool
    stationary: bool
    t: float
    target_maximum: float
    gap: float
    primal: float
    stationarity: float
    dual_sign: float
    complementarity: float
    point: list[float]
    lower: list[float]
    upper: list[float]
    values: list[float]
    physical_values: list[float]
    jacobian: list[float]
    multipliers: list[float]
    physical_multipliers: list[float]
    z_lower: list[float]
    z_upper: list[float]
    physical_z_lower: list[float]
    physical_z_upper: list[float]
    lower_distances: list[float]
    upper_distances: list[float]
    hold_multipliers: list[float]
    full_stationarity: list[float]
    root_errors: list[float]
    slack_groups: list[int]
    slack_rows: list[int]
    active_lower: list[int]
    active_upper: list[int]
    hold_causes: list[str]
    rows: list[MinimaxRow]
    frozen_checks: list[CheckResult]
class MeasuredInputKind:
    Temperature: MeasuredInputKind
    Pressure: MeasuredInputKind
    MoleFraction: MeasuredInputKind
    Molality: MeasuredInputKind
class MeasuredInput:
    def __init__(self) -> None: ...
    kind: MeasuredInputKind
    species: int
    sigma: float
    lower: float
    upper: float
    unit: str
class MeasurementGroup:
    def __init__(self) -> None: ...
    rows: list[int]
    inputs: list[MeasuredInput]
    dependent_species: list[int]
    feed_total_mol: float
class MeasurementUncertainty:
    def __init__(self) -> None: ...
    groups: list[MeasurementGroup]
    output_sigma: list[float]
    covariance: list[float] | None
    correlation: list[float] | None
class InputReconciliation:
    values: list[float]
    corrections: list[float]
    penalties: list[float]
    jacobian: list[float]
    row_groups: list[int]
    column_groups: list[int]
    column_inputs: list[int]
    output_cost: float
    input_cost: float
class FitResult:
    solver_options: list[dict[str, Any]]
    reconciliation: InputReconciliation | None
    minimax: MinimaxResult | None
    chart: EvaluatedCoordinates
    checks: list[CheckResult]
    acceptance: str
    optimizer: list[float]
    status: str
    message: str
    usable: bool
    # Native chart order: one value per target of each non-reaction coordinate (declaration order), then each held
    # Born diameter, then each reaction coordinate. Entry k < len(chart.active) is the value of chart.active[k].
    physical: list[float]
    parameters: Parameters
    predictions: list[float]  # per observation, in the observation's unit
    residuals: list[float]
    weighted_residuals: list[float]
    physical_jacobian: list[float]
    optimizer_jacobian: list[float]
    statuses: list[ObservationStatus]
    failures: list[TrialFailure]
    iterations: int
    residual_evaluations: int
    jacobian_evaluations: int
    trial_failures: int
    iteration_seconds: list[float]
    iteration_costs: list[float]
    initial_cost: float
    final_cost: float
    covariance: FitCovariance

def predict(model: HelmholtzModel, observation: Observation) -> Prediction: ...
def differentiate(model: HelmholtzModel, observation: Observation) -> PredictionGradient: ...
def fit(parameters: Parameters, thermochemistry: Optional[ThermochemistryRecord], coordinates: Sequence[Coordinate],
        observations: Sequence[Observation], weights: Sequence[float], controls: FitControls, holds: Sequence[HydrationHold] = ...,
        checks: Sequence[CheckGroup] = ..., baseline: Parameters | None = ..., singular_vectors: bool = ...,
        starts: Optional[Sequence[Sequence[float]]] = None, iterations: Optional[int] = None,
        seconds: Optional[float] = None, minimax: MinimaxRequest | None = ...,
        uncertainty: MeasurementUncertainty | None = ...) -> Union[FitResult, list[FitResult]]: ...

class CheckMetric:
    AbsoluteDifference: CheckMetric
    AbsoluteLog10Ratio: CheckMetric

class HydrationHold:
    def __init__(self) -> None: ...
    reference: Observation
    lower: float
    upper: float
    energy_tolerance: float
    diameter_tolerance: float
    maximum_iterations: int

class HoldResult:
    energy_error: float
    correction: float
    uncertainty: float
    derivative: float
    parameter_derivative: list[float]
    iterations: int

class EvaluatedCoordinates:
    reaction_laws: dict[str, ReactionLogPolynomial]
    parameters: Parameters
    active: list[ActiveParameter]
    physical: list[float]
    map: list[float]
    holds: list[HoldResult]

class CheckGroup:
    def __init__(self) -> None: ...
    name: str
    rows: list[int]
    metric: CheckMetric
    maximum: float | None
    mean: float | None
    baseline: bool
    report_only: bool

class CheckResult:
    name: str
    unit: str
    status: str
    rows: list[int]
    failed_rows: list[int]
    errors: list[float]
    slopes: list[float]
    maximum: float
    mean: float
    argmax: int
    maximum_threshold: float | None
    mean_threshold: float | None
    baseline_maximum: float
    baseline_mean: float
    report_only: bool

def evaluate_coordinates(parameters: Parameters, set: CoordinateSet, optimizer: Sequence[float], holds: Sequence[HydrationHold], thermochemistry: ThermochemistryRecord | None = ..., solver_options: Mapping[str, float] | None = None) -> EvaluatedCoordinates: ...
def evaluate_checks(stored: Parameters, observations: Sequence[Observation], checks: Sequence[CheckGroup], baseline: Parameters | None = ..., solver_options: Mapping[str, float] | None = None) -> list[CheckResult]: ...

def observations_at(parameters: Parameters, coordinates: Sequence[Coordinate], optimizer: Sequence[float], holds: Sequence[HydrationHold], observations: Sequence[Observation], errors: Sequence[float],
    thermochemistry: Optional[ThermochemistryRecord] = None, solver_options: Mapping[str, float] | None = None) -> list[Observation]: ...

class ObservationEvaluation:
    solver_options: list[dict[str, Any]]
    reconciliation: InputReconciliation | None
    state_densities: list[list[float]]
    reference_values: list[list[float]]
    predictions: list[float]
    gradient: list[float]
    optimizer_gradient: list[float]
    residual_optimizer_gradient: list[float]
    chart: EvaluatedCoordinates
    statuses: list[ObservationStatus]
def evaluate_observations(parameters: Parameters, coordinates: Sequence[Coordinate], optimizer: Sequence[float],
    holds: Sequence[HydrationHold], observations: Sequence[Observation],
    thermochemistry: Optional[ThermochemistryRecord] = None,
    uncertainty: MeasurementUncertainty | None = ..., solver_options: Mapping[str, float] | None = None) -> ObservationEvaluation: ...
def evaluate_minimax(parameters: Parameters, coordinates: Sequence[Coordinate], optimizer: Sequence[float],
    holds: Sequence[HydrationHold], observations: Sequence[Observation], checks: Sequence[CheckGroup],
    baseline: Parameters | None, minimax: MinimaxRequest, thermochemistry: Optional[ThermochemistryRecord] = None, solver_options: Mapping[str, float] | None = None) -> FitResult: ...
def formula_composition(solvent_fractions: Sequence[float], parameters: Parameters, ion_indices: tuple[int, int],
    stoichiometry: tuple[int, int], molality: float) -> list[float]: ...

class EquilibriumStart:
    coordinates: list[float]
    amounts: list[float]
    compositions: list[list[float]]
    packing_fractions: list[Optional[float]]
    temperature: Optional[float]
    pressure: Optional[float]
    previous: Optional[EquilibriumResult]
    def __init__(self, *, coordinates: Sequence[float] = ..., amounts: Sequence[float] = ...,
                 compositions: Sequence[Sequence[float]] = ..., packing_fractions: Sequence[Optional[float]] = ...,
                 temperature: Optional[float] = None, pressure: Optional[float] = None,
                 previous: Optional[EquilibriumResult] = None) -> None: ...

class EquilibriumWorkLimits:
    attempts: Optional[int]
    iterations: Optional[int]
    attempt_iterations: Optional[int]
    seconds: Optional[float]
    attempt_seconds: Optional[float]
    def __init__(self, *, attempts: Optional[int] = None, iterations: Optional[int] = None,
                 attempt_iterations: Optional[int] = None, seconds: Optional[float] = None,
                 attempt_seconds: Optional[float] = None) -> None: ...

class EquilibriumStartsResult:
    selected: Optional[int]
    results: list[Optional[EquilibriumResult]]
    causes: list[str]
    attempt_seconds: list[float]
    iterations: int
    seconds: float
    overrun_seconds: float
    budget_exhausted: bool

def solve_equilibrium_starts(model: HelmholtzModel, problem: EquilibriumProblem, starts: list[EquilibriumStart], limits: EquilibriumWorkLimits) -> EquilibriumStartsResult: ...

def reaction_log10_reciprocal(law: ReactionLogPolynomial) -> list[float]: ...
