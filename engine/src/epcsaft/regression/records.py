from __future__ import annotations

import csv
import hashlib
import io
import math
from collections import Counter
from dataclasses import dataclass
from importlib.resources import files
from itertools import pairwise

EXPECTED_HEADER = ("species", "T_K", "p_sat_Pa", "rho_sat_liq_kg_m3", "source")
PROPANE_EXPECTED_HEADER = (
    "row_id",
    "component_id",
    "role",
    "T_K",
    "p_sat_Pa",
    "p_sat_expanded_uncertainty_Pa",
    "rho_sat_liq_kg_m3",
    "rho_sat_liq_expanded_uncertainty_kg_m3",
    "rho_sat_vap_kg_m3",
    "rho_sat_vap_expanded_uncertainty_kg_m3",
)
SOURCE_RETRIEVED_ON = "2026-07-17"
SOURCE_USE_BASIS = (
    "NIST Standard Reference Data retained as compact source-backed candidate evidence; "
    "redistribution and use remain subject to the NIST SRD terms"
)
SOURCE_UNITS = (
    ("temperature", "K"),
    ("pressure", "Pa"),
    ("saturated_liquid_mass_density", "kg/m3"),
)

METHANE_DATA_SHA256 = "a5e16df3bf8ec78483fc340782cddc89ab8b284a9f6dfaecd6cda3ffde579227"
METHANE_PACKAGED_DATA_SHA256 = "dec64d5a6cac414a4a92393a0d728fa27c02135c6a159d0d1881d7b6dde6d26c"
METHANE_SOURCE_ID = "nist-webbook-srd69-methane-saturation"
METHANE_SOURCE_URL = (
    "https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C74828&"
    "Type=SatP&Digits=8&THigh=180&TLow=100&TInc=10&RefState=DEF&TUnit=K&"
    "PUnit=Pa&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm"
)
METHANE_SOURCE_CITATION = "NIST Chemistry WebBook, SRD 69, methane (CAS 74-82-8) fluid properties"
METHANE_SOURCE_LOCATOR = "Saturation properties query, 100 K through 180 K in 10 K increments"
METHANE_SOURCE_TRANSFORMATION = (
    "Exact retained CSV fields and decimal strings; CRLF line endings normalized to LF. "
    f"Retained source SHA-256: {METHANE_DATA_SHA256}; packaged SHA-256: "
    f"{METHANE_PACKAGED_DATA_SHA256}."
)
METHANE_TEMPERATURES_K = tuple(float(value) for value in range(100, 181, 10))
METHANE_TRAINING_TEMPERATURES_K = (110.0, 130.0, 150.0, 170.0)
METHANE_HELD_OUT_TEMPERATURES_K = (100.0, 120.0, 140.0, 160.0, 180.0)

ETHANE_DATA_SHA256 = "ed09b8781acfb7025ca505878b884f6353ddd9f3f4bd7aae2e6df88bbe847a67"
ETHANE_PACKAGED_DATA_SHA256 = "b01333e827933c0a7148672c8ae3eef78393320c0d18f2c4d5a0fc40d9bef6b2"
ETHANE_SOURCE_ID = "nist-webbook-srd69-ethane-saturation"
ETHANE_SOURCE_URL = (
    "https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C74840&"
    "Type=SatP&Digits=8&THigh=280&TLow=100&TInc=20&RefState=DEF&TUnit=K&"
    "PUnit=Pa&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm"
)
ETHANE_SOURCE_CITATION = "NIST Chemistry WebBook, SRD 69, ethane (CAS 74-84-0) fluid properties"
ETHANE_SOURCE_LOCATOR = "Saturation properties query, 100 K through 280 K in 20 K increments"
ETHANE_SOURCE_TRANSFORMATION = (
    "Exact retained CSV fields and decimal strings; CRLF line endings normalized to LF. "
    f"Retained source SHA-256: {ETHANE_DATA_SHA256}; packaged SHA-256: "
    f"{ETHANE_PACKAGED_DATA_SHA256}."
)
ETHANE_TEMPERATURES_K = tuple(float(value) for value in range(100, 281, 20))
ETHANE_TRAINING_TEMPERATURES_K = (140.0, 180.0, 220.0, 260.0)
ETHANE_HELD_OUT_TEMPERATURES_K = (120.0, 160.0, 200.0, 240.0)
ETHANE_STRESS_TEMPERATURES_K = (100.0, 280.0)

PROPANE_PACKAGED_DATA_SHA256 = "ccd1cfa15ec44432b06cbf22316d168c61b282631c9b1e1591e497b8d48b5676"
PROPANE_PACKET_YAML_SHA256 = "ba31448989f565d05d63908076e836977780aa87199f208310e9b80b03f64697"
PROPANE_SOURCE_RECEIPT_SHA256 = "ed5eb703ccd3e6bb4c4cfa82ecd58c58f9da0c93ab07a204dee94d8b0ae8d081"
PROPANE_FIT_TARGET_CONTRACT_SHA256 = (
    "7f25259265dfa42f1de36bc04740baf6c78e09c8bc35a42392f06a4b8a32cb90"
)
PROPANE_SOURCE_VERIFICATION_CONTRACT_SHA256 = (
    "b0cb440613ec5fc764d1ccce7c40af371af208a129bb211fb1d749d34046020c"
)
PROPANE_COMPARISON_CONTRACT_SHA256 = (
    "522b55f8c9641bab7b572f1741fc24cf48b7a2df10706ade17064cd4c79ba2f2"
)
PROPANE_THERMOML_JSON_SHA256 = "322495c5a01c003e83376e5bad544c3abced330d5054ff0411a7a00b70a963c9"
PROPANE_THERMOML_XML_SHA256 = "1b2e47d4cafff0f21cf7779d8d01b522bc2fa8d885ce4d6ebc04c151e0504829"
PROPANE_SOURCE_ID = "glos-2004-propane-coexistence-experiment"
PROPANE_SOURCE_URL = "https://trc.nist.gov/ThermoML/10.1016/j.jct.2004.07.017.json"
PROPANE_SOURCE_CITATION = (
    "Glos, Kleinrahm, and Wagner, Journal of Chemical Thermodynamics 36 "
    "(2004) 1037-1059, doi:10.1016/j.jct.2004.07.017"
)
PROPANE_SOURCE_LOCATOR = (
    "Table 2 propane coexistence measurements; NIST ThermoML datasets 1, 2, and 3"
)
PROPANE_SOURCE_USE_BASIS = (
    "Direct primary experimental Glos 2004 measurements retained from the hash-bound "
    "Validation packet as source evidence, not model-acceptance cutoffs"
)
PROPANE_SOURCE_TRANSFORMATION = (
    "Exact target CSV bytes from Validation commit "
    "7e51590757f1cb85f51df98e9fe1f88cd4255a88, tree "
    "05af9e948c786ddfcf43dba701970f1cbb6435a2. Target SHA-256: "
    f"{PROPANE_PACKAGED_DATA_SHA256}; packet YAML SHA-256: "
    f"{PROPANE_PACKET_YAML_SHA256}; 63-row source receipt SHA-256: "
    f"{PROPANE_SOURCE_RECEIPT_SHA256}; fit-target contract SHA-256: "
    f"{PROPANE_FIT_TARGET_CONTRACT_SHA256}; source-verification contract SHA-256: "
    f"{PROPANE_SOURCE_VERIFICATION_CONTRACT_SHA256}; comparison contract SHA-256: "
    f"{PROPANE_COMPARISON_CONTRACT_SHA256}; ThermoML JSON SHA-256: "
    f"{PROPANE_THERMOML_JSON_SHA256}; ThermoML XML SHA-256: "
    f"{PROPANE_THERMOML_XML_SHA256}. Pressure converted exactly from kPa to Pa by "
    "Validation; density units unchanged."
)
PROPANE_TEMPERATURES_K = tuple(float(value) for value in range(110, 341, 10))
PROPANE_TRAINING_TEMPERATURES_K = (150.0, 210.0, 270.0, 330.0)
PROPANE_HELD_OUT_TEMPERATURES_K = tuple(
    temperature
    for temperature in PROPANE_TEMPERATURES_K
    if temperature not in (*PROPANE_TRAINING_TEMPERATURES_K, 110.0, 340.0)
)
PROPANE_STRESS_TEMPERATURES_K = (110.0, 340.0)
PROPANE_ROW_IDS = tuple(
    f"glos2004-propane-sat-{int(temperature)}-k" for temperature in PROPANE_TEMPERATURES_K
)

# Retained import names used only by the accepted methane receipt tooling.
EXPECTED_DATA_SHA256 = METHANE_DATA_SHA256
EXPECTED_PACKAGED_DATA_SHA256 = METHANE_PACKAGED_DATA_SHA256

FIGIEL_TARGETS_PACKAGED_SHA256 = "10a6f7a777458fe29285461cc24028298d488a053a7d9565dc51a6a5949a4c42"
FIGIEL_VALIDATION_LEDGER_SHA256 = "f405a3e48d21cd979a8dd480d5f8cb3be40754f5d6babf368b505b5f305607f0"
FIGIEL_VALIDATION_PARAMETER_PACKET_SHA256 = (
    "932e8baa90fcefbaa8c3a8730cdeadd83a4c01f0a3b109f4e4cd0319aee9312b"
)
FIGIEL_VALIDATION_METADATA_SHA256 = (
    "8ea06c6ca5452d01448a03f9a76cf7d0c35bb99c9abe23ccb1729d56c71d468f"
)
FIGIEL_TARGET_HEADER = (
    "target_id",
    "ion_label",
    "active_component_id",
    "counterion_component_id",
    "target_kj_per_mol",
    "target_j_per_mol",
    "published_born_diameter_angstrom",
    "expected_provider_fingerprint",
)


@dataclass(frozen=True, slots=True)
class BornDiameterTarget:
    target_id: str
    ion_label: str
    active_component_id: str
    counterion_component_id: str
    target_kj_per_mol: float
    target_j_per_mol: float
    published_diameter_angstrom: float
    expected_provider_fingerprint: str

    @property
    def component_order(self) -> tuple[str, str, str]:
        return ("water", self.active_component_id, self.counterion_component_id)


@dataclass(frozen=True, slots=True)
class BornDiameterTracerSpecification:
    specification_id: str
    targets: tuple[BornDiameterTarget, ...]
    source_validation_commit: str
    source_validation_tree: str
    source_ledger_sha256: str
    source_parameter_packet_sha256: str
    source_metadata_sha256: str
    packaged_targets_sha256: str
    source_doi: str
    source_si_doi: str
    source_locator: str
    source_basis: str
    temperature_k: float
    pressure_pa: float
    reference_molality_mol_per_kg: float
    reference_convergence_error_max: float
    diameter_origin_angstrom: float
    diameter_scale_angstrom: float
    diameter_bounds_angstrom: tuple[float, float]
    scaled_bounds: tuple[float, float]
    start_diameters_angstrom: tuple[tuple[float, ...], ...]
    max_num_iterations: int
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float
    scaled_residual_max: float
    confirmation_parameter_scaled_max_delta: float
    observable_round_trip_j_per_mol: float
    published_diameter_reporting_half_increment_angstrom: float
    rank_threshold_multiplier: float
    ceres_linear_solver: str
    ceres_num_threads: int
    ceres_logging: str

    def __post_init__(self) -> None:
        expected_targets = (
            (
                "figiel2025-s5-Lip-reported-average",
                "Li+",
                "lithium-cation",
                "chloride-anion",
                -486.2,
                -486_200.0,
                2.784,
                "sha256:9cef6373572a729c7c3519ee04c261b4f9d821d37ae2b74dc916a4a02e863149",
            ),
            (
                "figiel2025-s5-Nap-reported-average",
                "Na+",
                "sodium-cation",
                "chloride-anion",
                -381.1,
                -381_100.0,
                3.445,
                "sha256:e97c462061bf92d2d75622a8d025a4f7e94106edccf2c924c88142b1b6a03439",
            ),
            (
                "figiel2025-s5-Kp-reported-average",
                "K+",
                "potassium-cation",
                "chloride-anion",
                -309.1,
                -309_100.0,
                4.150,
                "sha256:1e885c57ebf7226d28a153ad4c3afb5e9ed8786b04aa0ff1b65e776b157d67d3",
            ),
            (
                "figiel2025-s5-Clm-reported-average",
                "Cl-",
                "chloride-anion",
                "sodium-cation",
                -314.9,
                -314_900.0,
                4.100,
                "sha256:ac5a3ee0bf044d9fc63086a41f04f7b9f34a506daf67ac9627fd082a519c62a3",
            ),
            (
                "figiel2025-s5-Brm-reported-average",
                "Br-",
                "bromide-anion",
                "sodium-cation",
                -290.9,
                -290_900.0,
                4.480,
                "sha256:7db8316e0c2352d992c7482d4c1218463020d749d5365ec8a8f5e491643eb9b2",
            ),
        )
        observed_targets = tuple(
            (
                target.target_id,
                target.ion_label,
                target.active_component_id,
                target.counterion_component_id,
                target.target_kj_per_mol,
                target.target_j_per_mol,
                target.published_diameter_angstrom,
                target.expected_provider_fingerprint,
            )
            for target in self.targets
        )
        if observed_targets != expected_targets:
            raise ValueError("Born tracer targets must match the exact five-target contract")
        exact_identity = (
            self.specification_id,
            self.source_validation_commit,
            self.source_validation_tree,
            self.source_ledger_sha256,
            self.source_parameter_packet_sha256,
            self.source_metadata_sha256,
            self.packaged_targets_sha256,
            self.source_doi,
            self.source_si_doi,
            self.source_locator,
            self.source_basis,
        )
        expected_identity = (
            "figiel-2025-five-ion-born-diameter-tracer-v1",
            "8944d34f7002cda1bb8760e606cc1f11696f58cd",
            "6c8fd350dcd6bfdd7be1918f73fd33a23e2070dd",
            FIGIEL_VALIDATION_LEDGER_SHA256,
            FIGIEL_VALIDATION_PARAMETER_PACKET_SHA256,
            FIGIEL_VALIDATION_METADATA_SHA256,
            FIGIEL_TARGETS_PACKAGED_SHA256,
            "10.1021/acs.iecr.5c00475",
            "10.1021/acs.iecr.5c00475.s001",
            "SI Table S5, PDF page 9 of 10, reported-average lit column",
            (
                "x-treatment at infinite dilution; gas at 1 bar to hypothetical "
                "dilute-ideal aqueous solution at 1 bar; negative is favorable"
            ),
        )
        if exact_identity != expected_identity:
            raise ValueError("Born tracer source identity must match the approved packet")
        numerical_contract = (
            self.temperature_k,
            self.pressure_pa,
            self.reference_molality_mol_per_kg,
            self.reference_convergence_error_max,
            self.diameter_origin_angstrom,
            self.diameter_scale_angstrom,
            self.diameter_bounds_angstrom,
            self.scaled_bounds,
            self.start_diameters_angstrom,
            self.max_num_iterations,
            self.function_tolerance,
            self.gradient_tolerance,
            self.parameter_tolerance,
            self.scaled_residual_max,
            self.confirmation_parameter_scaled_max_delta,
            self.observable_round_trip_j_per_mol,
            self.published_diameter_reporting_half_increment_angstrom,
            self.rank_threshold_multiplier,
            self.ceres_linear_solver,
            self.ceres_num_threads,
            self.ceres_logging,
        )
        expected_numerical_contract = (
            298.15,
            100_000.0,
            1.0e-12,
            5.0e-5,
            3.0,
            1.0,
            (1.0, 6.0),
            (-2.0, 3.0),
            ((3.0,) * 5, (2.0,) * 5, (5.0,) * 5),
            500,
            1.0e-10,
            1.0e-10,
            1.0e-10,
            1.0e-8,
            1.0e-5,
            50.0,
            0.0005,
            100.0,
            "DENSE_QR",
            1,
            "SILENT",
        )
        if numerical_contract != expected_numerical_contract:
            raise ValueError("Born tracer numerical controls must match the frozen design")

    @property
    def target_ids(self) -> tuple[str, ...]:
        return tuple(target.target_id for target in self.targets)

    @property
    def ion_labels(self) -> tuple[str, ...]:
        return tuple(target.ion_label for target in self.targets)

    @property
    def active_component_ids(self) -> tuple[str, ...]:
        return tuple(target.active_component_id for target in self.targets)

    @property
    def targets_j_per_mol(self) -> tuple[float, ...]:
        return tuple(target.target_j_per_mol for target in self.targets)

    @property
    def published_diameters_angstrom(self) -> tuple[float, ...]:
        return tuple(target.published_diameter_angstrom for target in self.targets)


def _load_figiel_born_targets() -> tuple[BornDiameterTarget, ...]:
    data = (
        files("epcsaft.regression.data").joinpath("figiel-born-diameter-targets.csv").read_bytes()
    )
    if hashlib.sha256(data).hexdigest() != FIGIEL_TARGETS_PACKAGED_SHA256:
        raise ValueError("packaged Figiel Born target hash does not match the approved contract")
    reader = csv.DictReader(io.StringIO(data.decode("utf-8"), newline=""))
    if tuple(reader.fieldnames or ()) != FIGIEL_TARGET_HEADER:
        raise ValueError("packaged Figiel Born target header does not match the approved contract")
    targets = tuple(
        BornDiameterTarget(
            target_id=row["target_id"],
            ion_label=row["ion_label"],
            active_component_id=row["active_component_id"],
            counterion_component_id=row["counterion_component_id"],
            target_kj_per_mol=float(row["target_kj_per_mol"]),
            target_j_per_mol=float(row["target_j_per_mol"]),
            published_diameter_angstrom=float(row["published_born_diameter_angstrom"]),
            expected_provider_fingerprint=row["expected_provider_fingerprint"],
        )
        for row in reader
    )
    if any(target.target_j_per_mol != 1000.0 * target.target_kj_per_mol for target in targets):
        raise ValueError("Figiel target kJ/mol to J/mol conversion must be exact")
    return targets


FIGIEL_BORN_DIAMETER_TRACER_V1 = BornDiameterTracerSpecification(
    specification_id="figiel-2025-five-ion-born-diameter-tracer-v1",
    targets=_load_figiel_born_targets(),
    source_validation_commit="8944d34f7002cda1bb8760e606cc1f11696f58cd",
    source_validation_tree="6c8fd350dcd6bfdd7be1918f73fd33a23e2070dd",
    source_ledger_sha256=FIGIEL_VALIDATION_LEDGER_SHA256,
    source_parameter_packet_sha256=FIGIEL_VALIDATION_PARAMETER_PACKET_SHA256,
    source_metadata_sha256=FIGIEL_VALIDATION_METADATA_SHA256,
    packaged_targets_sha256=FIGIEL_TARGETS_PACKAGED_SHA256,
    source_doi="10.1021/acs.iecr.5c00475",
    source_si_doi="10.1021/acs.iecr.5c00475.s001",
    source_locator="SI Table S5, PDF page 9 of 10, reported-average lit column",
    source_basis=(
        "x-treatment at infinite dilution; gas at 1 bar to hypothetical dilute-ideal "
        "aqueous solution at 1 bar; negative is favorable"
    ),
    temperature_k=298.15,
    pressure_pa=100_000.0,
    reference_molality_mol_per_kg=1.0e-12,
    reference_convergence_error_max=5.0e-5,
    diameter_origin_angstrom=3.0,
    diameter_scale_angstrom=1.0,
    diameter_bounds_angstrom=(1.0, 6.0),
    scaled_bounds=(-2.0, 3.0),
    start_diameters_angstrom=((3.0,) * 5, (2.0,) * 5, (5.0,) * 5),
    max_num_iterations=500,
    function_tolerance=1.0e-10,
    gradient_tolerance=1.0e-10,
    parameter_tolerance=1.0e-10,
    scaled_residual_max=1.0e-8,
    confirmation_parameter_scaled_max_delta=1.0e-5,
    observable_round_trip_j_per_mol=50.0,
    published_diameter_reporting_half_increment_angstrom=0.0005,
    rank_threshold_multiplier=100.0,
    ceres_linear_solver="DENSE_QR",
    ceres_num_threads=1,
    ceres_logging="SILENT",
)

FIGIEL_NABR_MIAC_PACKAGED_SHA256 = (
    "0591d90f4775d672eb830229cf99a925523b101dcf3f4881a8425318da25d97d"
)
FIGIEL_NABR_MIAC_HEADER = ("salt", "molality_mol_kg", "gamma_pm_m")
FIGIEL_FIXED_BORN_DIAMETERS_ANGSTROM = (
    2.7888130173797934,
    3.4524616464076425,
    4.147266741279482,
    4.101505615791675,
    4.476998527506598,
)
FIGIEL_FIXED_WATER_SOLVATION_FACTOR = 1.5590515389548207
FIGIEL_AQUEOUS_KIJ_COORDINATES = (
    ("water", "lithium-cation"),
    ("water", "sodium-cation"),
    ("water", "potassium-cation"),
    ("water", "chloride-anion"),
    ("water", "bromide-anion"),
    ("lithium-cation", "chloride-anion"),
    ("sodium-cation", "chloride-anion"),
    ("potassium-cation", "chloride-anion"),
    ("lithium-cation", "bromide-anion"),
    ("sodium-cation", "bromide-anion"),
    ("potassium-cation", "bromide-anion"),
)
FIGIEL_AQUEOUS_PUBLISHED_KIJ = (
    -0.4,
    -0.3,
    -0.1,
    -0.3,
    -0.3,
    0.8,
    0.8,
    0.0,
    0.5,
    0.65,
    -0.35,
)


@dataclass(frozen=True, slots=True)
class AqueousMiacObservation:
    row_id: str
    salt: str
    cation_component_id: str
    anion_component_id: str
    molality_mol_per_kg: float
    gamma_pm_m: float

    def __post_init__(self) -> None:
        if (
            not self.row_id.strip()
            or self.salt != "NaBr"
            or self.cation_component_id != "sodium-cation"
            or self.anion_component_id != "bromide-anion"
        ):
            raise ValueError("water-factor observation must identify one NaBr row")
        if not math.isfinite(self.molality_mol_per_kg) or self.molality_mol_per_kg <= 0.0:
            raise ValueError("water-factor molality must be positive and finite")
        if not math.isfinite(self.gamma_pm_m) or self.gamma_pm_m <= 0.0:
            raise ValueError("water-factor MIAC must be positive and finite")


@dataclass(frozen=True, slots=True)
class FigielWaterSolvationFactorSpecification:
    specification_id: str
    observations: tuple[AqueousMiacObservation, ...]
    source_validation_commit: str
    source_validation_tree: str
    source_ledger_sha256: str
    source_parameter_packet_sha256: str
    source_metadata_sha256: str
    source_si_extraction_sha256: str
    source_hamer_wu_csv_sha256: str
    packaged_nabr_sha256: str
    fixed_born_evidence_file_sha256: str
    fixed_born_evidence_subject_sha256: str
    fixed_born_diameters_angstrom: tuple[float, ...]
    fixed_aqueous_kij_coordinate_order: tuple[tuple[str, str], ...]
    fixed_aqueous_kij: tuple[float, ...]
    expected_provider_fingerprint: str
    temperature_k: float
    pressure_pa: float
    parameter_name: str
    parameter_bounds: tuple[float, float]
    starts: tuple[float, float]
    max_num_iterations: int
    start_wall_time_max_seconds: float
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float
    rank_threshold_multiplier: float
    start_agreement_max_abs: float
    residual_definition: str
    ceres_linear_solver: str
    ceres_num_threads: int
    ceres_logging: str

    def __post_init__(self) -> None:
        if len(self.observations) != 21 or len({row.row_id for row in self.observations}) != 21:
            raise ValueError("water-factor contract requires 21 unique NaBr rows")
        exact_contract = (
            self.specification_id,
            self.source_validation_commit,
            self.source_validation_tree,
            self.source_ledger_sha256,
            self.source_parameter_packet_sha256,
            self.source_metadata_sha256,
            self.source_si_extraction_sha256,
            self.source_hamer_wu_csv_sha256,
            self.packaged_nabr_sha256,
            self.fixed_born_evidence_file_sha256,
            self.fixed_born_evidence_subject_sha256,
            self.fixed_born_diameters_angstrom,
            self.fixed_aqueous_kij_coordinate_order,
            self.fixed_aqueous_kij,
            self.expected_provider_fingerprint,
            tuple(row.molality_mol_per_kg for row in self.observations),
            self.temperature_k,
            self.pressure_pa,
            self.parameter_name,
            self.parameter_bounds,
            self.starts,
            self.max_num_iterations,
            self.start_wall_time_max_seconds,
            self.function_tolerance,
            self.gradient_tolerance,
            self.parameter_tolerance,
            self.rank_threshold_multiplier,
            self.start_agreement_max_abs,
            self.residual_definition,
            self.ceres_linear_solver,
            self.ceres_num_threads,
            self.ceres_logging,
        )
        expected_contract = (
            "figiel-2025-water-solvation-factor-nabr-v1",
            "8944d34f7002cda1bb8760e606cc1f11696f58cd",
            "6c8fd350dcd6bfdd7be1918f73fd33a23e2070dd",
            FIGIEL_VALIDATION_LEDGER_SHA256,
            FIGIEL_VALIDATION_PARAMETER_PACKET_SHA256,
            FIGIEL_VALIDATION_METADATA_SHA256,
            "85bd39f727158d5a9d6eea6828c1673f73850e783a655b09660cc9b66d84321a",
            "2f63e13f06a5b0f4e8bca2980b6a8d9d7fb0f839153c43e3a71952daf9796595",
            FIGIEL_NABR_MIAC_PACKAGED_SHA256,
            "99d46eafbdae3428f690543364096fb414b818db201d1c35b0c0da8b03ae91d5",
            "55ea2cd69af62c45b26179cfab6939760de23058b5a7e8c880a79f67faa417ed",
            FIGIEL_FIXED_BORN_DIAMETERS_ANGSTROM,
            FIGIEL_AQUEOUS_KIJ_COORDINATES,
            FIGIEL_AQUEOUS_PUBLISHED_KIJ,
            "sha256:80ec4bde74ef7af307e44098e2e495e20a5be57d7cb23452ff321cb644816783",
            (
                0.001,
                0.002,
                0.005,
                0.01,
                0.02,
                0.05,
                0.1,
                0.6,
                0.7,
                0.8,
                0.9,
                1.0,
                1.2,
                1.4,
                3.0,
                3.5,
                4.0,
                4.5,
                5.0,
                5.5,
                6.0,
            ),
            298.15,
            100_000.0,
            "water_solvation_factor",
            (1.0, 2.0),
            (1.2, 1.8),
            500,
            180.0,
            1.0e-10,
            1.0e-10,
            1.0e-10,
            100.0,
            1.0e-5,
            "1 - gamma_model/gamma_observed",
            "DENSE_QR",
            1,
            "SILENT",
        )
        if exact_contract != expected_contract:
            raise ValueError("water-factor controls must match the frozen one-stage contract")


def _load_figiel_nabr_miac_targets() -> tuple[AqueousMiacObservation, ...]:
    data = files("epcsaft.regression.data").joinpath("figiel-nabr-miac-targets.csv").read_bytes()
    if hashlib.sha256(data).hexdigest() != FIGIEL_NABR_MIAC_PACKAGED_SHA256:
        raise ValueError("packaged NaBr MIAC hash does not match the frozen packet")
    reader = csv.DictReader(io.StringIO(data.decode("utf-8"), newline=""))
    if tuple(reader.fieldnames or ()) != FIGIEL_NABR_MIAC_HEADER:
        raise ValueError("packaged NaBr MIAC header does not match the frozen packet")
    return tuple(
        AqueousMiacObservation(
            row_id=f"hamer-wu-1972-nabr-{index:03d}",
            salt=row["salt"],
            cation_component_id="sodium-cation",
            anion_component_id="bromide-anion",
            molality_mol_per_kg=float(row["molality_mol_kg"]),
            gamma_pm_m=float(row["gamma_pm_m"]),
        )
        for index, row in enumerate(reader, start=1)
    )


FIGIEL_WATER_SOLVATION_FACTOR_V1 = FigielWaterSolvationFactorSpecification(
    specification_id="figiel-2025-water-solvation-factor-nabr-v1",
    observations=_load_figiel_nabr_miac_targets(),
    source_validation_commit="8944d34f7002cda1bb8760e606cc1f11696f58cd",
    source_validation_tree="6c8fd350dcd6bfdd7be1918f73fd33a23e2070dd",
    source_ledger_sha256=FIGIEL_VALIDATION_LEDGER_SHA256,
    source_parameter_packet_sha256=FIGIEL_VALIDATION_PARAMETER_PACKET_SHA256,
    source_metadata_sha256=FIGIEL_VALIDATION_METADATA_SHA256,
    source_si_extraction_sha256=(
        "85bd39f727158d5a9d6eea6828c1673f73850e783a655b09660cc9b66d84321a"
    ),
    source_hamer_wu_csv_sha256=("2f63e13f06a5b0f4e8bca2980b6a8d9d7fb0f839153c43e3a71952daf9796595"),
    packaged_nabr_sha256=FIGIEL_NABR_MIAC_PACKAGED_SHA256,
    fixed_born_evidence_file_sha256=(
        "99d46eafbdae3428f690543364096fb414b818db201d1c35b0c0da8b03ae91d5"
    ),
    fixed_born_evidence_subject_sha256=(
        "55ea2cd69af62c45b26179cfab6939760de23058b5a7e8c880a79f67faa417ed"
    ),
    fixed_born_diameters_angstrom=FIGIEL_FIXED_BORN_DIAMETERS_ANGSTROM,
    fixed_aqueous_kij_coordinate_order=FIGIEL_AQUEOUS_KIJ_COORDINATES,
    fixed_aqueous_kij=FIGIEL_AQUEOUS_PUBLISHED_KIJ,
    expected_provider_fingerprint=(
        "sha256:80ec4bde74ef7af307e44098e2e495e20a5be57d7cb23452ff321cb644816783"
    ),
    temperature_k=298.15,
    pressure_pa=100_000.0,
    parameter_name="water_solvation_factor",
    parameter_bounds=(1.0, 2.0),
    starts=(1.2, 1.8),
    max_num_iterations=500,
    start_wall_time_max_seconds=180.0,
    function_tolerance=1.0e-10,
    gradient_tolerance=1.0e-10,
    parameter_tolerance=1.0e-10,
    rank_threshold_multiplier=100.0,
    start_agreement_max_abs=1.0e-5,
    residual_definition="1 - gamma_model/gamma_observed",
    ceres_linear_solver="DENSE_QR",
    ceres_num_threads=1,
    ceres_logging="SILENT",
)

FIGIEL_AQUEOUS_MIAC_PACKAGED_SHA256 = (
    "2f63e13f06a5b0f4e8bca2980b6a8d9d7fb0f839153c43e3a71952daf9796595"
)
FIGIEL_AQUEOUS_SALT_CONTRACTS = (
    ("LiCl", "lithium-cation", "chloride-anion", (0, 3, 5)),
    ("NaCl", "sodium-cation", "chloride-anion", (1, 3, 6)),
    ("KCl", "potassium-cation", "chloride-anion", (2, 3, 7)),
    ("LiBr", "lithium-cation", "bromide-anion", (0, 4, 8)),
    ("NaBr", "sodium-cation", "bromide-anion", (1, 4, 9)),
    ("KBr", "potassium-cation", "bromide-anion", (2, 4, 10)),
)
FIGIEL_AQUEOUS_KIJ_PROVIDER_FINGERPRINTS = (
    "sha256:00d5bc709841bffd8cbde9df6afafaa9520d23a500be59aae54ecd63b4369ad9",
    "sha256:569f4ace0d776e8b631df0ff563a2686356b28ad5e24556c38909b371fd811ac",
    "sha256:17955d601749bc94802325dbc88cdc62b3ee954f1e5e1ad02c0547d8187b27d4",
    "sha256:31f40e1441358b8729c0963c808cb5290cde6d455f6f910a5c7abe95315a8aea",
    "sha256:8fedd3d395a20591b05ee7fff2c1e698872c97b07e82a0668923d656f818790e",
    "sha256:c0af15d984b8324ebc4b622a9ab68153c43eb13026bd74a6bcd08bba7d49dea5",
)


@dataclass(frozen=True, slots=True)
class AqueousKijObservation:
    row_id: str
    salt: str
    cation_component_id: str
    anion_component_id: str
    molality_mol_per_kg: float
    gamma_pm_m: float

    def __post_init__(self) -> None:
        expected = {
            salt: (cation, anion) for salt, cation, anion, _ in FIGIEL_AQUEOUS_SALT_CONTRACTS
        }
        if not self.row_id.strip() or expected.get(self.salt) != (
            self.cation_component_id,
            self.anion_component_id,
        ):
            raise ValueError("aqueous-kij observation has the wrong salt identity")
        if not math.isfinite(self.molality_mol_per_kg) or self.molality_mol_per_kg <= 0.0:
            raise ValueError("aqueous-kij molality must be positive and finite")
        if not math.isfinite(self.gamma_pm_m) or self.gamma_pm_m <= 0.0:
            raise ValueError("aqueous-kij MIAC must be positive and finite")


@dataclass(frozen=True, slots=True)
class FigielAqueousKijSpecification:
    specification_id: str
    observations: tuple[AqueousKijObservation, ...]
    source_validation_commit: str
    source_validation_tree: str
    source_hamer_wu_csv_sha256: str
    packaged_data_sha256: str
    fixed_born_evidence_file_sha256: str
    fixed_born_evidence_subject_sha256: str
    fixed_water_factor_regression_commit: str
    fixed_water_factor_regression_tree: str
    fixed_born_diameters_angstrom: tuple[float, ...]
    fixed_water_solvation_factor: float
    salt_contracts: tuple[tuple[str, str, str, tuple[int, int, int]], ...]
    coordinate_order: tuple[tuple[str, str], ...]
    published_parameters: tuple[float, ...]
    expected_provider_fingerprints: tuple[str, ...]
    temperature_k: float
    pressure_pa: float
    parameter_bounds: tuple[float, float]
    start_schedules: tuple[tuple[str, float, str], ...]
    max_num_iterations: int
    start_wall_time_max_seconds: float
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float
    rank_threshold_multiplier: float
    start_agreement_max_abs: float
    published_parameter_max_abs_delta: float
    residual_definition: str
    ceres_linear_solver: str
    ceres_num_threads: int
    ceres_logging: str

    def __post_init__(self) -> None:
        counts = Counter(row.salt for row in self.observations)
        if len(self.observations) != 164 or len({row.row_id for row in self.observations}) != 164:
            raise ValueError("aqueous-kij contract requires 164 unique rows")
        if counts != Counter(
            {"LiCl": 29, "NaCl": 29, "KCl": 28, "LiBr": 29, "NaBr": 21, "KBr": 28}
        ):
            raise ValueError("aqueous-kij salt counts do not match the frozen packet")
        if (
            self.specification_id != "figiel-2025-aqueous-kij-v1"
            or self.source_validation_commit != "8944d34f7002cda1bb8760e606cc1f11696f58cd"
            or self.source_validation_tree != "6c8fd350dcd6bfdd7be1918f73fd33a23e2070dd"
            or self.source_hamer_wu_csv_sha256 != FIGIEL_AQUEOUS_MIAC_PACKAGED_SHA256
            or self.packaged_data_sha256 != FIGIEL_AQUEOUS_MIAC_PACKAGED_SHA256
            or self.fixed_born_evidence_file_sha256
            != "99d46eafbdae3428f690543364096fb414b818db201d1c35b0c0da8b03ae91d5"
            or self.fixed_born_evidence_subject_sha256
            != "55ea2cd69af62c45b26179cfab6939760de23058b5a7e8c880a79f67faa417ed"
            or self.fixed_water_factor_regression_commit
            != "882e0735ed1b5586a591682da1fd3d78f46636d4"
            or self.fixed_water_factor_regression_tree != "9860f7bf3f0548a6a7a527da627e7f4d97a47797"
            or self.fixed_born_diameters_angstrom != FIGIEL_FIXED_BORN_DIAMETERS_ANGSTROM
            or self.fixed_water_solvation_factor != FIGIEL_FIXED_WATER_SOLVATION_FACTOR
            or self.salt_contracts != FIGIEL_AQUEOUS_SALT_CONTRACTS
            or self.coordinate_order != FIGIEL_AQUEOUS_KIJ_COORDINATES
            or self.published_parameters != FIGIEL_AQUEOUS_PUBLISHED_KIJ
            or self.expected_provider_fingerprints != FIGIEL_AQUEOUS_KIJ_PROVIDER_FINGERPRINTS
            or self.temperature_k != 298.15
            or self.pressure_pa != 100_000.0
            or self.parameter_bounds != (-1.0, 1.0)
            or self.start_schedules
            != (
                ("primary", 0.0, "forward"),
                ("confirmation", 0.25, "reverse"),
            )
            or self.max_num_iterations != 50
            or self.start_wall_time_max_seconds != 180.0
            or self.function_tolerance != 1.0e-10
            or self.gradient_tolerance != 1.0e-10
            or self.parameter_tolerance != 1.0e-10
            or self.rank_threshold_multiplier != 100.0
            or self.start_agreement_max_abs != 1.0e-5
            or self.published_parameter_max_abs_delta != 0.05
            or self.residual_definition != "1 - gamma_model/gamma_observed"
            or self.ceres_linear_solver != "DENSE_QR"
            or self.ceres_num_threads != 1
            or self.ceres_logging != "SILENT"
        ):
            raise ValueError("aqueous-kij controls do not match the frozen contract")


def _load_figiel_aqueous_miac_targets() -> tuple[AqueousKijObservation, ...]:
    data = files("epcsaft.regression.data").joinpath("figiel-aqueous-miac-targets.csv").read_bytes()
    if hashlib.sha256(data).hexdigest() != FIGIEL_AQUEOUS_MIAC_PACKAGED_SHA256:
        raise ValueError("packaged aqueous MIAC hash does not match the frozen packet")
    reader = csv.DictReader(io.StringIO(data.decode("utf-8"), newline=""))
    if tuple(reader.fieldnames or ()) != FIGIEL_NABR_MIAC_HEADER:
        raise ValueError("packaged aqueous MIAC header does not match the frozen packet")
    salt_contracts = {
        salt: (cation, anion) for salt, cation, anion, _ in FIGIEL_AQUEOUS_SALT_CONTRACTS
    }
    counters: Counter[str] = Counter()
    observations: list[AqueousKijObservation] = []
    for row in reader:
        counters[row["salt"]] += 1
        cation, anion = salt_contracts[row["salt"]]
        observations.append(
            AqueousKijObservation(
                row_id=(f"hamer-wu-1972-{row['salt'].lower()}-{counters[row['salt']]:03d}"),
                salt=row["salt"],
                cation_component_id=cation,
                anion_component_id=anion,
                molality_mol_per_kg=float(row["molality_mol_kg"]),
                gamma_pm_m=float(row["gamma_pm_m"]),
            )
        )
    return tuple(observations)


FIGIEL_AQUEOUS_KIJ_V1 = FigielAqueousKijSpecification(
    specification_id="figiel-2025-aqueous-kij-v1",
    observations=_load_figiel_aqueous_miac_targets(),
    source_validation_commit="8944d34f7002cda1bb8760e606cc1f11696f58cd",
    source_validation_tree="6c8fd350dcd6bfdd7be1918f73fd33a23e2070dd",
    source_hamer_wu_csv_sha256=FIGIEL_AQUEOUS_MIAC_PACKAGED_SHA256,
    packaged_data_sha256=FIGIEL_AQUEOUS_MIAC_PACKAGED_SHA256,
    fixed_born_evidence_file_sha256=(
        "99d46eafbdae3428f690543364096fb414b818db201d1c35b0c0da8b03ae91d5"
    ),
    fixed_born_evidence_subject_sha256=(
        "55ea2cd69af62c45b26179cfab6939760de23058b5a7e8c880a79f67faa417ed"
    ),
    fixed_water_factor_regression_commit=("882e0735ed1b5586a591682da1fd3d78f46636d4"),
    fixed_water_factor_regression_tree=("9860f7bf3f0548a6a7a527da627e7f4d97a47797"),
    fixed_born_diameters_angstrom=FIGIEL_FIXED_BORN_DIAMETERS_ANGSTROM,
    fixed_water_solvation_factor=FIGIEL_FIXED_WATER_SOLVATION_FACTOR,
    salt_contracts=FIGIEL_AQUEOUS_SALT_CONTRACTS,
    coordinate_order=FIGIEL_AQUEOUS_KIJ_COORDINATES,
    published_parameters=FIGIEL_AQUEOUS_PUBLISHED_KIJ,
    expected_provider_fingerprints=FIGIEL_AQUEOUS_KIJ_PROVIDER_FINGERPRINTS,
    temperature_k=298.15,
    pressure_pa=100_000.0,
    parameter_bounds=(-1.0, 1.0),
    start_schedules=(
        ("primary", 0.0, "forward"),
        ("confirmation", 0.25, "reverse"),
    ),
    max_num_iterations=50,
    start_wall_time_max_seconds=180.0,
    function_tolerance=1.0e-10,
    gradient_tolerance=1.0e-10,
    parameter_tolerance=1.0e-10,
    rank_threshold_multiplier=100.0,
    start_agreement_max_abs=1.0e-5,
    published_parameter_max_abs_delta=0.05,
    residual_definition="1 - gamma_model/gamma_observed",
    ceres_linear_solver="DENSE_QR",
    ceres_num_threads=1,
    ceres_logging="SILENT",
)


def _positive_finite(value: float, field: str) -> None:
    if not math.isfinite(value):
        raise ValueError(f"{field} must be finite")
    if value <= 0.0:
        raise ValueError(f"{field} must be positive")


def _source_fields(source_id: str) -> tuple[str, str, str, str, str, str]:
    if source_id == METHANE_SOURCE_ID:
        return (
            METHANE_SOURCE_CITATION,
            METHANE_SOURCE_LOCATOR,
            METHANE_SOURCE_URL,
            METHANE_SOURCE_TRANSFORMATION,
            METHANE_DATA_SHA256,
            METHANE_PACKAGED_DATA_SHA256,
        )
    if source_id == ETHANE_SOURCE_ID:
        return (
            ETHANE_SOURCE_CITATION,
            ETHANE_SOURCE_LOCATOR,
            ETHANE_SOURCE_URL,
            ETHANE_SOURCE_TRANSFORMATION,
            ETHANE_DATA_SHA256,
            ETHANE_PACKAGED_DATA_SHA256,
        )
    if source_id == PROPANE_SOURCE_ID:
        return (
            PROPANE_SOURCE_CITATION,
            PROPANE_SOURCE_LOCATOR,
            PROPANE_SOURCE_URL,
            PROPANE_SOURCE_TRANSFORMATION,
            PROPANE_PACKAGED_DATA_SHA256,
            PROPANE_PACKAGED_DATA_SHA256,
        )
    raise ValueError("source_id must identify an admitted pure-saturation table")


def _source_provenance(source_id: str) -> tuple[str, str, tuple[tuple[str, str], ...]]:
    if source_id in (METHANE_SOURCE_ID, ETHANE_SOURCE_ID):
        return SOURCE_RETRIEVED_ON, SOURCE_USE_BASIS, SOURCE_UNITS
    if source_id == PROPANE_SOURCE_ID:
        return SOURCE_RETRIEVED_ON, PROPANE_SOURCE_USE_BASIS, SOURCE_UNITS
    raise ValueError("source_id must identify an admitted pure-saturation table")


@dataclass(frozen=True, slots=True)
class SourceIdentity:
    source_id: str
    citation: str
    locator: str
    url: str
    retrieved_on: str
    use_basis: str
    transformation: str
    units: tuple[tuple[str, str], ...]
    data_sha256: str
    packaged_data_sha256: str

    def __post_init__(self) -> None:
        for field in (
            "source_id",
            "citation",
            "locator",
            "url",
            "retrieved_on",
            "use_basis",
            "transformation",
        ):
            if not getattr(self, field).strip():
                raise ValueError(f"source {field} must be nonblank")
        expected = _source_fields(self.source_id)
        if (
            self.citation,
            self.locator,
            self.url,
            self.transformation,
            self.data_sha256,
            self.packaged_data_sha256,
        ) != expected:
            raise ValueError("source fields must match the exact retained source identity")
        if (self.retrieved_on, self.use_basis, self.units) != _source_provenance(self.source_id):
            raise ValueError("source provenance must match the exact retained source identity")


METHANE_SOURCE_V1 = SourceIdentity(
    source_id=METHANE_SOURCE_ID,
    citation=METHANE_SOURCE_CITATION,
    locator=METHANE_SOURCE_LOCATOR,
    url=METHANE_SOURCE_URL,
    retrieved_on=SOURCE_RETRIEVED_ON,
    use_basis=SOURCE_USE_BASIS,
    transformation=METHANE_SOURCE_TRANSFORMATION,
    units=SOURCE_UNITS,
    data_sha256=METHANE_DATA_SHA256,
    packaged_data_sha256=METHANE_PACKAGED_DATA_SHA256,
)
ETHANE_SOURCE_V1 = SourceIdentity(
    source_id=ETHANE_SOURCE_ID,
    citation=ETHANE_SOURCE_CITATION,
    locator=ETHANE_SOURCE_LOCATOR,
    url=ETHANE_SOURCE_URL,
    retrieved_on=SOURCE_RETRIEVED_ON,
    use_basis=SOURCE_USE_BASIS,
    transformation=ETHANE_SOURCE_TRANSFORMATION,
    units=SOURCE_UNITS,
    data_sha256=ETHANE_DATA_SHA256,
    packaged_data_sha256=ETHANE_PACKAGED_DATA_SHA256,
)
PROPANE_SOURCE_V1 = SourceIdentity(
    source_id=PROPANE_SOURCE_ID,
    citation=PROPANE_SOURCE_CITATION,
    locator=PROPANE_SOURCE_LOCATOR,
    url=PROPANE_SOURCE_URL,
    retrieved_on=SOURCE_RETRIEVED_ON,
    use_basis=PROPANE_SOURCE_USE_BASIS,
    transformation=PROPANE_SOURCE_TRANSFORMATION,
    units=SOURCE_UNITS,
    data_sha256=PROPANE_PACKAGED_DATA_SHA256,
    packaged_data_sha256=PROPANE_PACKAGED_DATA_SHA256,
)


@dataclass(frozen=True, slots=True)
class SaturationObservation:
    row_id: str
    component_id: str
    temperature_k: float
    pressure_pa: float
    liquid_density_kg_m3: float
    source_id: str
    pressure_expanded_uncertainty_pa: float | None = None
    liquid_density_expanded_uncertainty_kg_m3: float | None = None
    vapor_density_kg_m3: float | None = None
    vapor_density_expanded_uncertainty_kg_m3: float | None = None

    def __post_init__(self) -> None:
        if not self.row_id.strip():
            raise ValueError("row_id must be nonblank")
        if self.component_id not in ("methane", "ethane", "propane"):
            raise ValueError("component_id must be 'methane', 'ethane', or 'propane'")
        _positive_finite(self.temperature_k, "temperature_k")
        _positive_finite(self.pressure_pa, "pressure_pa")
        _positive_finite(self.liquid_density_kg_m3, "liquid_density_kg_m3")
        expected_source = {
            "methane": METHANE_SOURCE_ID,
            "ethane": ETHANE_SOURCE_ID,
            "propane": PROPANE_SOURCE_ID,
        }[self.component_id]
        if self.source_id != expected_source:
            raise ValueError("source_id does not match the pure component")
        uncertainty_fields = (
            (self.pressure_expanded_uncertainty_pa, "pressure uncertainty"),
            (
                self.liquid_density_expanded_uncertainty_kg_m3,
                "liquid density uncertainty",
            ),
            (self.vapor_density_kg_m3, "vapor density"),
            (
                self.vapor_density_expanded_uncertainty_kg_m3,
                "vapor density uncertainty",
            ),
        )
        for value, field in uncertainty_fields:
            if value is not None:
                _positive_finite(value, field)
        if self.component_id != "propane" and any(
            value is not None for value, _ in uncertainty_fields
        ):
            raise ValueError("source uncertainty fields belong only to the propane packet")
        if (self.vapor_density_kg_m3 is None) != (
            self.vapor_density_expanded_uncertainty_kg_m3 is None
        ):
            raise ValueError("vapor density and uncertainty must be present together")


@dataclass(frozen=True, slots=True)
class PureSaturationDataset:
    dataset_id: str
    component_id: str
    temperature_unit: str
    pressure_unit: str
    liquid_density_unit: str
    source: SourceIdentity
    rows: tuple[SaturationObservation, ...]
    training_temperatures_k: tuple[float, ...]
    held_out_temperatures_k: tuple[float, ...]
    stress_temperatures_k: tuple[float, ...]

    def __post_init__(self) -> None:
        if self.component_id == "methane":
            expected = (
                "nist-webbook-methane-saturation-100-180-k-v1",
                METHANE_SOURCE_ID,
                METHANE_TEMPERATURES_K,
                METHANE_TRAINING_TEMPERATURES_K,
                METHANE_HELD_OUT_TEMPERATURES_K,
                (),
            )
        elif self.component_id == "ethane":
            expected = (
                "nist-webbook-ethane-saturation-100-280-k-v1",
                ETHANE_SOURCE_ID,
                ETHANE_TEMPERATURES_K,
                ETHANE_TRAINING_TEMPERATURES_K,
                ETHANE_HELD_OUT_TEMPERATURES_K,
                ETHANE_STRESS_TEMPERATURES_K,
            )
        elif self.component_id == "propane":
            expected = (
                "glos-2004-experimental-propane-saturation-110-340-k-v1",
                PROPANE_SOURCE_ID,
                PROPANE_TEMPERATURES_K,
                PROPANE_TRAINING_TEMPERATURES_K,
                PROPANE_HELD_OUT_TEMPERATURES_K,
                PROPANE_STRESS_TEMPERATURES_K,
            )
        else:
            raise ValueError("component_id must be 'methane', 'ethane', or 'propane'")
        dataset_id, source_id, temperatures, training, held_out, stress = expected
        if self.dataset_id != dataset_id:
            raise ValueError("dataset_id does not match the admitted component table")
        if self.source.source_id != source_id:
            raise ValueError("dataset source does not match its component")
        if (self.temperature_unit, self.pressure_unit, self.liquid_density_unit) != (
            "K",
            "Pa",
            "kg/m3",
        ):
            raise ValueError("dataset units must be K, Pa, and kg/m3")
        row_ids = tuple(row.row_id for row in self.rows)
        if len(set(row_ids)) != len(row_ids):
            raise ValueError("dataset contains a duplicate row_id")
        observed_temperatures = tuple(row.temperature_k for row in self.rows)
        if any(right <= left for left, right in pairwise(observed_temperatures)):
            raise ValueError("dataset temperatures must be strictly increasing")
        if observed_temperatures != temperatures:
            raise ValueError("dataset temperatures do not match the retained reporting grid")
        expected_row_ids = (
            PROPANE_ROW_IDS
            if self.component_id == "propane"
            else tuple(
                f"nist-{self.component_id}-sat-{int(temperature)}-k" for temperature in temperatures
            )
        )
        if row_ids != expected_row_ids:
            raise ValueError("dataset row IDs must match their retained temperatures")
        if any(
            row.component_id != self.component_id or row.source_id != source_id for row in self.rows
        ):
            raise ValueError("every row must share the dataset component and source identity")
        if (
            self.training_temperatures_k,
            self.held_out_temperatures_k,
            self.stress_temperatures_k,
        ) != (training, held_out, stress):
            raise ValueError("dataset row partition does not match the admitted specification")
        if (
            set(training) & set(held_out)
            or set(training) & set(stress)
            or set(held_out) & set(stress)
        ):
            raise ValueError("dataset row partitions must be disjoint")
        if set(training) | set(held_out) | set(stress) != set(temperatures):
            raise ValueError("dataset row partitions must cover every retained row")

    def _rows_at(self, temperatures: tuple[float, ...]) -> tuple[SaturationObservation, ...]:
        selected = frozenset(temperatures)
        return tuple(row for row in self.rows if row.temperature_k in selected)

    @property
    def training_rows(self) -> tuple[SaturationObservation, ...]:
        return self._rows_at(self.training_temperatures_k)

    @property
    def held_out_rows(self) -> tuple[SaturationObservation, ...]:
        return self._rows_at(self.held_out_temperatures_k)

    @property
    def stress_rows(self) -> tuple[SaturationObservation, ...]:
        return self._rows_at(self.stress_temperatures_k)

    @property
    def training_row_ids(self) -> tuple[str, ...]:
        return tuple(row.row_id for row in self.training_rows)


@dataclass(frozen=True, slots=True)
class PureSaturationFitSpecification:
    specification_id: str
    component_id: str
    dataset_id: str
    source_id: str
    expected_provider_fingerprint: str
    parameter_names: tuple[str, str, str]
    parameter_units: tuple[str, str, str]
    start: tuple[float, float, float]
    lower_bounds: tuple[float, float, float]
    upper_bounds: tuple[float, float, float]
    parameter_scales: tuple[float, float, float]
    fixed_amount_mol: float
    molar_mass_kg_per_mol: float
    residual_names: tuple[str, str, str, str]
    residual_weights: tuple[float, float, float, float]
    liquid_volume_bounds_m3: tuple[float, float]
    vapor_volume_bounds_m3: tuple[float, float]
    training_temperatures_k: tuple[float, float, float, float]
    max_num_iterations: int
    function_tolerance: float
    gradient_tolerance: float
    parameter_tolerance: float
    topology_relative_separation_min: float
    reporting_pressure_bounds_pa: tuple[float, float]
    confirmation_liquid_volume_start_multiplier: float
    confirmation_vapor_volume_start_multiplier: float
    confirmation_parameter_scaled_max_delta: float
    confirmation_cost_relative_delta: float
    reporting_pressure_scaled_residual_max: float
    reporting_chemical_potential_residual_max: float
    ceres_linear_solver: str
    ceres_num_threads: int
    ceres_logging: str

    def __post_init__(self) -> None:
        if self.component_id == "methane":
            expected_identity = (
                "pure-methane-saturation-lifted-volumes-v1",
                "nist-webbook-methane-saturation-100-180-k-v1",
                METHANE_SOURCE_ID,
                "sha256:4a5cde0fad05a150c7fa54bb5ac7db508424f0126cf512596dcc248c284dd9e0",
                (1.08, 3.555744, 157.5315),
                0.016043,
                (2.0e-5, 1.0e-4),
                (1.5e-4, 0.1),
                METHANE_TRAINING_TEMPERATURES_K,
                (1.0e3, 1.0e7),
            )
        elif self.component_id == "ethane":
            expected_identity = (
                "pure-ethane-saturation-lifted-volumes-v1",
                "nist-webbook-ethane-saturation-100-280-k-v1",
                ETHANE_SOURCE_ID,
                "sha256:6381ef30d4f25c63fee1fd098c7024dd59254f3021d0588ee46e4eecfb31619b",
                (1.6069, 3.5206, 191.42),
                0.030070,
                (2.0e-5, 1.0e-4),
                (1.5e-4, 100.0),
                ETHANE_TRAINING_TEMPERATURES_K,
                (1.0, 1.0e7),
            )
        elif self.component_id == "propane":
            expected_identity = (
                "pure-propane-saturation-lifted-volumes-v1",
                "glos-2004-experimental-propane-saturation-110-340-k-v1",
                PROPANE_SOURCE_ID,
                "sha256:7d893c35288fceec76b7dde4a16fbf2ef95830b7d18827d3205537ce54878140",
                (2.002, 3.6184, 208.11),
                0.044096,
                (2.0e-5, 1.2e-4),
                (1.5e-4, 2.0e3),
                PROPANE_TRAINING_TEMPERATURES_K,
                (0.1, 1.0e7),
            )
        else:
            raise ValueError("component_id must be 'methane', 'ethane', or 'propane'")
        if (
            self.specification_id,
            self.dataset_id,
            self.source_id,
            self.expected_provider_fingerprint,
            self.start,
            self.molar_mass_kg_per_mol,
            self.liquid_volume_bounds_m3,
            self.vapor_volume_bounds_m3,
            self.training_temperatures_k,
            self.reporting_pressure_bounds_pa,
        ) != expected_identity:
            raise ValueError("fit specification does not match the admitted component contract")
        if self.parameter_names != (
            "segment_count",
            "segment_diameter_angstrom",
            "dispersion_energy_over_k_kelvin",
        ) or self.parameter_units != ("1", "angstrom", "K"):
            raise ValueError(
                "parameter names and units must match the provider coordinate contract"
            )
        if any(
            not math.isfinite(value)
            for group in (
                self.start,
                self.lower_bounds,
                self.upper_bounds,
                self.parameter_scales,
            )
            for value in group
        ):
            raise ValueError("parameter contract values must be finite")
        if any(
            not lower < start < upper
            for lower, start, upper in zip(
                self.lower_bounds, self.start, self.upper_bounds, strict=True
            )
        ):
            raise ValueError("every parameter start must lie strictly inside its bounds")
        if self.lower_bounds != (0.5, 2.0, 50.0) or self.upper_bounds != (
            3.5,
            5.0,
            400.0,
        ):
            raise ValueError("parameter bounds do not match the pure-saturation contract")
        if self.parameter_scales != (0.1, 0.1, 10.0):
            raise ValueError("parameter scales do not match the pure-saturation contract")
        _positive_finite(self.fixed_amount_mol, "fixed_amount_mol")
        if self.fixed_amount_mol != 1.0:
            raise ValueError("the pure-saturation slice fixes amount at exactly one mole")
        _positive_finite(self.molar_mass_kg_per_mol, "molar mass")
        if self.residual_names != (
            "liquid_pressure",
            "vapor_pressure",
            "chemical_potential_equality",
            "liquid_density",
        ):
            raise ValueError("residual ordering must match the lifted-volume formulation")
        if self.residual_weights != (0.25, 0.25, 0.25, 0.25):
            raise ValueError("residual weights must be the declared equal row normalization")
        if self.liquid_volume_bounds_m3[1] >= self.vapor_volume_bounds_m3[0]:
            raise ValueError("phase volume bounds must enforce liquid below vapor")
        expected_max_num_iterations = 5000 if self.component_id == "propane" else 500
        if self.max_num_iterations != expected_max_num_iterations:
            raise ValueError("max_num_iterations does not match the pure-saturation contract")
        for value in (
            self.function_tolerance,
            self.gradient_tolerance,
            self.parameter_tolerance,
            self.topology_relative_separation_min,
            self.confirmation_liquid_volume_start_multiplier,
            self.confirmation_vapor_volume_start_multiplier,
            self.confirmation_parameter_scaled_max_delta,
            self.confirmation_cost_relative_delta,
            self.reporting_pressure_scaled_residual_max,
            self.reporting_chemical_potential_residual_max,
        ):
            _positive_finite(value, "solver tolerance")
        if (
            self.function_tolerance,
            self.gradient_tolerance,
            self.parameter_tolerance,
            self.topology_relative_separation_min,
        ) != (1.0e-10, 1.0e-10, 1.0e-10, 1.0e-3):
            raise ValueError("solver tolerances do not match the pure-saturation contract")
        if (
            self.confirmation_liquid_volume_start_multiplier,
            self.confirmation_vapor_volume_start_multiplier,
        ) != (1.01, 0.98):
            raise ValueError("confirmation start multipliers do not match the contract")
        if (
            self.confirmation_parameter_scaled_max_delta,
            self.confirmation_cost_relative_delta,
        ) != (1.0e-5, 1.0e-8):
            raise ValueError("confirmation acceptance thresholds do not match the contract")
        if (
            self.reporting_pressure_scaled_residual_max,
            self.reporting_chemical_potential_residual_max,
        ) != (1.0e-8, 1.0e-8):
            raise ValueError("reporting closure thresholds do not match the contract")
        if (self.ceres_linear_solver, self.ceres_num_threads, self.ceres_logging) != (
            "DENSE_QR",
            1,
            "SILENT",
        ):
            raise ValueError("Ceres execution controls do not match the contract")


def _fit_specification(
    *,
    component_id: str,
    specification_id: str,
    dataset_id: str,
    source_id: str,
    expected_provider_fingerprint: str,
    start: tuple[float, float, float],
    molar_mass_kg_per_mol: float,
    vapor_volume_bounds_m3: tuple[float, float],
    training_temperatures_k: tuple[float, float, float, float],
    reporting_pressure_bounds_pa: tuple[float, float],
    liquid_volume_bounds_m3: tuple[float, float] = (2.0e-5, 1.0e-4),
    max_num_iterations: int = 500,
) -> PureSaturationFitSpecification:
    return PureSaturationFitSpecification(
        specification_id=specification_id,
        component_id=component_id,
        dataset_id=dataset_id,
        source_id=source_id,
        expected_provider_fingerprint=expected_provider_fingerprint,
        parameter_names=(
            "segment_count",
            "segment_diameter_angstrom",
            "dispersion_energy_over_k_kelvin",
        ),
        parameter_units=("1", "angstrom", "K"),
        start=start,
        lower_bounds=(0.5, 2.0, 50.0),
        upper_bounds=(3.5, 5.0, 400.0),
        parameter_scales=(0.1, 0.1, 10.0),
        fixed_amount_mol=1.0,
        molar_mass_kg_per_mol=molar_mass_kg_per_mol,
        residual_names=(
            "liquid_pressure",
            "vapor_pressure",
            "chemical_potential_equality",
            "liquid_density",
        ),
        residual_weights=(0.25, 0.25, 0.25, 0.25),
        liquid_volume_bounds_m3=liquid_volume_bounds_m3,
        vapor_volume_bounds_m3=vapor_volume_bounds_m3,
        training_temperatures_k=training_temperatures_k,
        max_num_iterations=max_num_iterations,
        function_tolerance=1.0e-10,
        gradient_tolerance=1.0e-10,
        parameter_tolerance=1.0e-10,
        topology_relative_separation_min=1.0e-3,
        reporting_pressure_bounds_pa=reporting_pressure_bounds_pa,
        confirmation_liquid_volume_start_multiplier=1.01,
        confirmation_vapor_volume_start_multiplier=0.98,
        confirmation_parameter_scaled_max_delta=1.0e-5,
        confirmation_cost_relative_delta=1.0e-8,
        reporting_pressure_scaled_residual_max=1.0e-8,
        reporting_chemical_potential_residual_max=1.0e-8,
        ceres_linear_solver="DENSE_QR",
        ceres_num_threads=1,
        ceres_logging="SILENT",
    )


METHANE_SATURATION_FIT_V1 = _fit_specification(
    component_id="methane",
    specification_id="pure-methane-saturation-lifted-volumes-v1",
    dataset_id="nist-webbook-methane-saturation-100-180-k-v1",
    source_id=METHANE_SOURCE_ID,
    expected_provider_fingerprint=(
        "sha256:4a5cde0fad05a150c7fa54bb5ac7db508424f0126cf512596dcc248c284dd9e0"
    ),
    start=(1.08, 3.555744, 157.5315),
    molar_mass_kg_per_mol=0.016043,
    vapor_volume_bounds_m3=(1.5e-4, 0.1),
    training_temperatures_k=METHANE_TRAINING_TEMPERATURES_K,
    reporting_pressure_bounds_pa=(1.0e3, 1.0e7),
)
ETHANE_SATURATION_FIT_V1 = _fit_specification(
    component_id="ethane",
    specification_id="pure-ethane-saturation-lifted-volumes-v1",
    dataset_id="nist-webbook-ethane-saturation-100-280-k-v1",
    source_id=ETHANE_SOURCE_ID,
    expected_provider_fingerprint=(
        "sha256:6381ef30d4f25c63fee1fd098c7024dd59254f3021d0588ee46e4eecfb31619b"
    ),
    start=(1.6069, 3.5206, 191.42),
    molar_mass_kg_per_mol=0.030070,
    vapor_volume_bounds_m3=(1.5e-4, 100.0),
    training_temperatures_k=ETHANE_TRAINING_TEMPERATURES_K,
    reporting_pressure_bounds_pa=(1.0, 1.0e7),
)
PROPANE_SATURATION_FIT_V1 = _fit_specification(
    component_id="propane",
    specification_id="pure-propane-saturation-lifted-volumes-v1",
    dataset_id="glos-2004-experimental-propane-saturation-110-340-k-v1",
    source_id=PROPANE_SOURCE_ID,
    expected_provider_fingerprint=(
        "sha256:7d893c35288fceec76b7dde4a16fbf2ef95830b7d18827d3205537ce54878140"
    ),
    start=(2.002, 3.6184, 208.11),
    molar_mass_kg_per_mol=0.044096,
    liquid_volume_bounds_m3=(2.0e-5, 1.2e-4),
    vapor_volume_bounds_m3=(1.5e-4, 2.0e3),
    training_temperatures_k=PROPANE_TRAINING_TEMPERATURES_K,
    reporting_pressure_bounds_pa=(0.1, 1.0e7),
    max_num_iterations=5000,
)


def _load_dataset(
    filename: str,
    component_id: str,
    display_name: str,
    source: SourceIdentity,
    dataset_id: str,
    temperatures: tuple[float, ...],
    training: tuple[float, ...],
    held_out: tuple[float, ...],
    stress: tuple[float, ...],
) -> PureSaturationDataset:
    data = files("epcsaft.regression").joinpath(f"data/{filename}").read_bytes()
    if hashlib.sha256(data).hexdigest() != source.packaged_data_sha256:
        raise ValueError(f"packaged {component_id} data SHA-256 does not match the source record")
    parsed = csv.reader(io.StringIO(data.decode("utf-8"), newline=""))
    if tuple(next(parsed)) != EXPECTED_HEADER:
        raise ValueError(f"{component_id} source data header or units changed")
    rows: list[SaturationObservation] = []
    for values in parsed:
        if len(values) != len(EXPECTED_HEADER):
            raise ValueError(f"{component_id} source row has a missing field")
        species, temperature, pressure, density, source_url = values
        if species != display_name or source_url != source.url:
            raise ValueError(f"{component_id} source identity changed")
        temperature_k = float(temperature)
        rows.append(
            SaturationObservation(
                row_id=f"nist-{component_id}-sat-{int(temperature_k)}-k",
                component_id=component_id,
                temperature_k=temperature_k,
                pressure_pa=float(pressure),
                liquid_density_kg_m3=float(density),
                source_id=source.source_id,
            )
        )
    if tuple(row.temperature_k for row in rows) != temperatures:
        raise ValueError(f"{component_id} source rows changed")
    return PureSaturationDataset(
        dataset_id=dataset_id,
        component_id=component_id,
        temperature_unit="K",
        pressure_unit="Pa",
        liquid_density_unit="kg/m3",
        source=source,
        rows=tuple(rows),
        training_temperatures_k=training,
        held_out_temperatures_k=held_out,
        stress_temperatures_k=stress,
    )


def _load_propane_dataset() -> PureSaturationDataset:
    data = files("epcsaft.regression").joinpath("data/propane_saturation.csv").read_bytes()
    if hashlib.sha256(data).hexdigest() != PROPANE_PACKAGED_DATA_SHA256:
        raise ValueError("packaged propane data SHA-256 does not match the source record")
    parsed = csv.reader(io.StringIO(data.decode("utf-8"), newline=""))
    if tuple(next(parsed)) != PROPANE_EXPECTED_HEADER:
        raise ValueError("propane source data header or units changed")
    rows: list[SaturationObservation] = []
    for values in parsed:
        if len(values) != len(PROPANE_EXPECTED_HEADER):
            raise ValueError("propane source row has a missing field")
        (
            row_id,
            component_id,
            role,
            temperature,
            pressure,
            pressure_uncertainty,
            liquid_density,
            liquid_density_uncertainty,
            vapor_density,
            vapor_density_uncertainty,
        ) = values
        temperature_k = float(temperature)
        expected_role = (
            "training"
            if temperature_k in PROPANE_TRAINING_TEMPERATURES_K
            else "held_out"
            if temperature_k in PROPANE_HELD_OUT_TEMPERATURES_K
            else "stress"
        )
        if component_id != "propane" or role != expected_role:
            raise ValueError("propane source component or frozen role changed")
        rows.append(
            SaturationObservation(
                row_id=row_id,
                component_id=component_id,
                temperature_k=temperature_k,
                pressure_pa=float(pressure),
                liquid_density_kg_m3=float(liquid_density),
                source_id=PROPANE_SOURCE_ID,
                pressure_expanded_uncertainty_pa=float(pressure_uncertainty),
                liquid_density_expanded_uncertainty_kg_m3=float(liquid_density_uncertainty),
                vapor_density_kg_m3=float(vapor_density) if vapor_density else None,
                vapor_density_expanded_uncertainty_kg_m3=(
                    float(vapor_density_uncertainty) if vapor_density_uncertainty else None
                ),
            )
        )
    return PureSaturationDataset(
        dataset_id="glos-2004-experimental-propane-saturation-110-340-k-v1",
        component_id="propane",
        temperature_unit="K",
        pressure_unit="Pa",
        liquid_density_unit="kg/m3",
        source=PROPANE_SOURCE_V1,
        rows=tuple(rows),
        training_temperatures_k=PROPANE_TRAINING_TEMPERATURES_K,
        held_out_temperatures_k=PROPANE_HELD_OUT_TEMPERATURES_K,
        stress_temperatures_k=PROPANE_STRESS_TEMPERATURES_K,
    )


def load_pure_saturation_dataset(component_id: str) -> PureSaturationDataset:
    if type(component_id) is not str:
        raise TypeError("component_id must be an exact string")
    if component_id == "methane":
        return _load_dataset(
            "methane_saturation.csv",
            "methane",
            "Methane",
            METHANE_SOURCE_V1,
            "nist-webbook-methane-saturation-100-180-k-v1",
            METHANE_TEMPERATURES_K,
            METHANE_TRAINING_TEMPERATURES_K,
            METHANE_HELD_OUT_TEMPERATURES_K,
            (),
        )
    if component_id == "ethane":
        return _load_dataset(
            "ethane_saturation.csv",
            "ethane",
            "Ethane",
            ETHANE_SOURCE_V1,
            "nist-webbook-ethane-saturation-100-280-k-v1",
            ETHANE_TEMPERATURES_K,
            ETHANE_TRAINING_TEMPERATURES_K,
            ETHANE_HELD_OUT_TEMPERATURES_K,
            ETHANE_STRESS_TEMPERATURES_K,
        )
    if component_id == "propane":
        return _load_propane_dataset()
    raise ValueError("component_id must be 'methane', 'ethane', or 'propane'")
