#include "pure_saturation_fit_internal.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace epcsaft_engine::regression::internal {
namespace {

std::string text(PyObject* object, const char* label) {
    if (!PyUnicode_Check(object)) {
        throw std::invalid_argument(std::string(label) + " must be text");
    }
    Py_ssize_t size = 0;
    const char* value = PyUnicode_AsUTF8AndSize(object, &size);
    if (value == nullptr) {
        throw std::invalid_argument(std::string(label) + " must be UTF-8 text");
    }
    return std::string(value, static_cast<std::size_t>(size));
}

std::vector<std::string> texts(PyObject* object, std::size_t expected, const char* label) {
    OwnedPyObject sequence{PySequence_Fast(object, label)};
    if (sequence == nullptr) throw std::invalid_argument(label);
    if (PySequence_Fast_GET_SIZE(sequence.get()) != static_cast<Py_ssize_t>(expected)) {
        throw std::invalid_argument(std::string(label) + " has the wrong length");
    }
    std::vector<std::string> values;
    values.reserve(expected);
    for (std::size_t index = 0; index < expected; ++index) {
        values.push_back(text(
            PySequence_Fast_GET_ITEM(sequence.get(), static_cast<Py_ssize_t>(index)), label
        ));
    }
    return values;
}

template <std::size_t Size>
std::array<double, Size> fixed_doubles(PyObject* object, const char* label) {
    const std::vector<double> parsed = doubles(object, Size, label);
    std::array<double, Size> values{};
    std::copy(parsed.begin(), parsed.end(), values.begin());
    return values;
}

constexpr std::array<const char*, 9> methane_row_ids{
    "nist-methane-sat-100-k", "nist-methane-sat-110-k", "nist-methane-sat-120-k",
    "nist-methane-sat-130-k", "nist-methane-sat-140-k", "nist-methane-sat-150-k",
    "nist-methane-sat-160-k", "nist-methane-sat-170-k", "nist-methane-sat-180-k",
};
constexpr std::array<double, 9> methane_pressures{
    34375.892, 88130.038, 191430.08, 367319.94, 641181.43,
    1039961.3, 1592078.0, 2328348.8, 3285180.7,
};
constexpr std::array<double, 9> methane_liquid_densities{
    438.88524, 424.77725, 409.90234, 394.03734, 376.86505,
    357.89846, 336.31495, 310.50203, 276.22850,
};

constexpr std::array<const char*, 10> ethane_row_ids{
    "nist-ethane-sat-100-k", "nist-ethane-sat-120-k", "nist-ethane-sat-140-k",
    "nist-ethane-sat-160-k", "nist-ethane-sat-180-k", "nist-ethane-sat-200-k",
    "nist-ethane-sat-220-k", "nist-ethane-sat-240-k", "nist-ethane-sat-260-k",
    "nist-ethane-sat-280-k",
};
constexpr std::array<double, 10> ethane_pressures{
    11.080787, 352.30167, 3813.5564, 21405.224, 78638.137,
    217232.94, 492046.39, 966787.79, 1711835.4, 2806735.8,
};
constexpr std::array<double, 10> ethane_liquid_densities{
    640.94852, 618.94997, 596.58156, 573.55144, 549.50874,
    523.97698, 496.27145, 465.30887, 429.07617, 382.72712,
};

constexpr std::array<const char*, 24> propane_row_ids{
    "glos2004-propane-sat-110-k", "glos2004-propane-sat-120-k",
    "glos2004-propane-sat-130-k", "glos2004-propane-sat-140-k",
    "glos2004-propane-sat-150-k", "glos2004-propane-sat-160-k",
    "glos2004-propane-sat-170-k", "glos2004-propane-sat-180-k",
    "glos2004-propane-sat-190-k", "glos2004-propane-sat-200-k",
    "glos2004-propane-sat-210-k", "glos2004-propane-sat-220-k",
    "glos2004-propane-sat-230-k", "glos2004-propane-sat-240-k",
    "glos2004-propane-sat-250-k", "glos2004-propane-sat-260-k",
    "glos2004-propane-sat-270-k", "glos2004-propane-sat-280-k",
    "glos2004-propane-sat-290-k", "glos2004-propane-sat-300-k",
    "glos2004-propane-sat-310-k", "glos2004-propane-sat-320-k",
    "glos2004-propane-sat-330-k", "glos2004-propane-sat-340-k",
};
constexpr std::array<double, 24> propane_pressures{
    0.6, 3.2, 18.0, 78.0, 283.0, 851.0, 2205.0, 5068.0,
    10547.0, 20193.0, 36032.0, 60574.0, 96775.0, 148000.0,
    217964.0, 310685.0, 430425.0, 581684.0, 769143.0, 997682.0,
    1272430.0, 1598870.0, 1983000.0, 2431450.0,
};
constexpr std::array<double, 24> propane_liquid_densities{
    707.968, 697.825, 687.713, 677.601, 667.462, 657.272,
    647.004, 636.628, 626.123, 615.456, 604.594, 593.499,
    582.132, 570.444, 558.383, 545.88, 532.853, 519.198,
    504.796, 489.465, 472.968, 454.951, 434.869, 411.772,
};

std::vector<std::string> expected_identity(const std::string& component_id) {
    if (component_id == "methane") {
        return {
            "nist-webbook-methane-saturation-100-180-k-v1",
            "methane", "K", "Pa", "kg/m3",
            "nist-webbook-srd69-methane-saturation",
            "NIST Chemistry WebBook, SRD 69, methane (CAS 74-82-8) fluid properties",
            "Saturation properties query, 100 K through 180 K in 10 K increments",
            "https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C74828&Type=SatP&Digits=8&THigh=180&TLow=100&TInc=10&RefState=DEF&TUnit=K&PUnit=Pa&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm",
            "2026-07-17",
            "NIST Standard Reference Data retained as compact source-backed candidate evidence; redistribution and use remain subject to the NIST SRD terms",
            "Exact retained CSV fields and decimal strings; CRLF line endings normalized to LF. Retained source SHA-256: a5e16df3bf8ec78483fc340782cddc89ab8b284a9f6dfaecd6cda3ffde579227; packaged SHA-256: dec64d5a6cac414a4a92393a0d728fa27c02135c6a159d0d1881d7b6dde6d26c.",
            "temperature", "K", "pressure", "Pa", "saturated_liquid_mass_density", "kg/m3",
            "a5e16df3bf8ec78483fc340782cddc89ab8b284a9f6dfaecd6cda3ffde579227",
            "dec64d5a6cac414a4a92393a0d728fa27c02135c6a159d0d1881d7b6dde6d26c",
            "pure-methane-saturation-lifted-volumes-v1",
            "segment_count", "segment_diameter_angstrom", "dispersion_energy_over_k_kelvin",
            "1", "angstrom", "K",
            "liquid_pressure", "vapor_pressure", "chemical_potential_equality", "liquid_density",
            "mol", "m3/mol", "epcsaft.native_sdk.v1",
            "sha256:4a5cde0fad05a150c7fa54bb5ac7db508424f0126cf512596dcc248c284dd9e0",
            "p_j = start_j + parameter_scale_j * z_j",
            "V_liquid = (molar_mass / observed_liquid_density) * exp(u_liquid)",
            "V_vapor = (R * T / observed_pressure) * exp(u_vapor)",
            "P_report = observed_pressure * exp(u_pressure)",
            "DENSE_QR", "SILENT",
        };
    }
    if (component_id == "ethane") {
        return {
            "nist-webbook-ethane-saturation-100-280-k-v1",
            "ethane", "K", "Pa", "kg/m3",
            "nist-webbook-srd69-ethane-saturation",
            "NIST Chemistry WebBook, SRD 69, ethane (CAS 74-84-0) fluid properties",
            "Saturation properties query, 100 K through 280 K in 20 K increments",
            "https://webbook.nist.gov/cgi/fluid.cgi?Action=Data&Wide=on&ID=C74840&Type=SatP&Digits=8&THigh=280&TLow=100&TInc=20&RefState=DEF&TUnit=K&PUnit=Pa&DUnit=kg%2Fm3&HUnit=kJ%2Fmol&WUnit=m%2Fs&VisUnit=uPa*s&STUnit=N%2Fm",
            "2026-07-17",
            "NIST Standard Reference Data retained as compact source-backed candidate evidence; redistribution and use remain subject to the NIST SRD terms",
            "Exact retained CSV fields and decimal strings; CRLF line endings normalized to LF. Retained source SHA-256: ed09b8781acfb7025ca505878b884f6353ddd9f3f4bd7aae2e6df88bbe847a67; packaged SHA-256: b01333e827933c0a7148672c8ae3eef78393320c0d18f2c4d5a0fc40d9bef6b2.",
            "temperature", "K", "pressure", "Pa", "saturated_liquid_mass_density", "kg/m3",
            "ed09b8781acfb7025ca505878b884f6353ddd9f3f4bd7aae2e6df88bbe847a67",
            "b01333e827933c0a7148672c8ae3eef78393320c0d18f2c4d5a0fc40d9bef6b2",
            "pure-ethane-saturation-lifted-volumes-v1",
            "segment_count", "segment_diameter_angstrom", "dispersion_energy_over_k_kelvin",
            "1", "angstrom", "K",
            "liquid_pressure", "vapor_pressure", "chemical_potential_equality", "liquid_density",
            "mol", "m3/mol", "epcsaft.native_sdk.v1",
            "sha256:6381ef30d4f25c63fee1fd098c7024dd59254f3021d0588ee46e4eecfb31619b",
            "p_j = start_j + parameter_scale_j * z_j",
            "V_liquid = (molar_mass / observed_liquid_density) * exp(u_liquid)",
            "V_vapor = (R * T / observed_pressure) * exp(u_vapor)",
            "P_report = observed_pressure * exp(u_pressure)",
            "DENSE_QR", "SILENT",
        };
    }
    if (component_id == "propane") {
        return {
            "glos-2004-experimental-propane-saturation-110-340-k-v1",
            "propane", "K", "Pa", "kg/m3",
            "glos-2004-propane-coexistence-experiment",
            "Glos, Kleinrahm, and Wagner, Journal of Chemical Thermodynamics 36 (2004) 1037-1059, doi:10.1016/j.jct.2004.07.017",
            "Table 2 propane coexistence measurements; NIST ThermoML datasets 1, 2, and 3",
            "https://trc.nist.gov/ThermoML/10.1016/j.jct.2004.07.017.json",
            "2026-07-17",
            "Direct primary experimental Glos 2004 measurements retained from the hash-bound Validation packet as source evidence, not model-acceptance cutoffs",
            "Exact target CSV bytes from Validation commit 7e51590757f1cb85f51df98e9fe1f88cd4255a88, tree 05af9e948c786ddfcf43dba701970f1cbb6435a2. Target SHA-256: ccd1cfa15ec44432b06cbf22316d168c61b282631c9b1e1591e497b8d48b5676; packet YAML SHA-256: ba31448989f565d05d63908076e836977780aa87199f208310e9b80b03f64697; 63-row source receipt SHA-256: ed5eb703ccd3e6bb4c4cfa82ecd58c58f9da0c93ab07a204dee94d8b0ae8d081; fit-target contract SHA-256: 7f25259265dfa42f1de36bc04740baf6c78e09c8bc35a42392f06a4b8a32cb90; source-verification contract SHA-256: b0cb440613ec5fc764d1ccce7c40af371af208a129bb211fb1d749d34046020c; comparison contract SHA-256: 522b55f8c9641bab7b572f1741fc24cf48b7a2df10706ade17064cd4c79ba2f2; ThermoML JSON SHA-256: 322495c5a01c003e83376e5bad544c3abced330d5054ff0411a7a00b70a963c9; ThermoML XML SHA-256: 1b2e47d4cafff0f21cf7779d8d01b522bc2fa8d885ce4d6ebc04c151e0504829. Pressure converted exactly from kPa to Pa by Validation; density units unchanged.",
            "temperature", "K", "pressure", "Pa", "saturated_liquid_mass_density", "kg/m3",
            "ccd1cfa15ec44432b06cbf22316d168c61b282631c9b1e1591e497b8d48b5676",
            "ccd1cfa15ec44432b06cbf22316d168c61b282631c9b1e1591e497b8d48b5676",
            "pure-propane-saturation-lifted-volumes-v1",
            "segment_count", "segment_diameter_angstrom", "dispersion_energy_over_k_kelvin",
            "1", "angstrom", "K",
            "liquid_pressure", "vapor_pressure", "chemical_potential_equality", "liquid_density",
            "mol", "m3/mol", "epcsaft.native_sdk.v1",
            "sha256:7d893c35288fceec76b7dde4a16fbf2ef95830b7d18827d3205537ce54878140",
            "p_j = start_j + parameter_scale_j * z_j",
            "V_liquid = (molar_mass / observed_liquid_density) * exp(u_liquid)",
            "V_vapor = (R * T / observed_pressure) * exp(u_vapor)",
            "P_report = observed_pressure * exp(u_pressure)",
            "DENSE_QR", "SILENT",
        };
    }
    throw std::invalid_argument("component identity must be methane, ethane, or propane");
}

bool component_specific_contract_matches(const Payload& payload) {
    if (payload.identity[1] == "methane") {
        return payload.start == std::array<double, 3>{1.08, 3.555744, 157.5315}
            && payload.molar_mass == 0.016043
            && payload.liquid_volume_bounds == std::array<double, 2>{2.0e-5, 1.0e-4}
            && payload.vapor_volume_bounds == std::array<double, 2>{1.5e-4, 0.1}
            && payload.max_iterations == 500
            && payload.reporting_pressure_bounds == std::array<double, 2>{1.0e3, 1.0e7};
    }
    if (payload.identity[1] == "ethane") {
        return payload.start == std::array<double, 3>{1.6069, 3.5206, 191.42}
            && payload.molar_mass == 0.030070
            && payload.liquid_volume_bounds == std::array<double, 2>{2.0e-5, 1.0e-4}
            && payload.vapor_volume_bounds == std::array<double, 2>{1.5e-4, 100.0}
            && payload.max_iterations == 500
            && payload.reporting_pressure_bounds == std::array<double, 2>{1.0, 1.0e7};
    }
    if (payload.identity[1] == "propane") {
        return payload.start == std::array<double, 3>{2.002, 3.6184, 208.11}
            && payload.molar_mass == 0.044096
            && payload.liquid_volume_bounds == std::array<double, 2>{2.0e-5, 1.2e-4}
            && payload.vapor_volume_bounds == std::array<double, 2>{1.5e-4, 2.0e3}
            && payload.max_iterations == 5000
            && payload.reporting_pressure_bounds == std::array<double, 2>{0.1, 1.0e7};
    }
    return false;
}

}  // namespace

bool positive_finite(double value) {
    return std::isfinite(value) && value > 0.0;
}

std::vector<double> doubles(PyObject* object, std::size_t expected, const char* label) {
    OwnedPyObject sequence{PySequence_Fast(object, label)};
    if (sequence == nullptr) throw std::invalid_argument(label);
    const Py_ssize_t size = PySequence_Fast_GET_SIZE(sequence.get());
    if (size != static_cast<Py_ssize_t>(expected)) {
        throw std::invalid_argument(std::string(label) + " has the wrong length");
    }
    std::vector<double> values;
    values.reserve(expected);
    for (Py_ssize_t index = 0; index < size; ++index) {
        const double value = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(sequence.get(), index));
        if (PyErr_Occurred() != nullptr) {
            throw std::invalid_argument(std::string(label) + " must contain numbers");
        }
        values.push_back(value);
    }
    return values;
}

std::size_t reporting_row_count(const Payload& payload) {
    if (payload.identity[1] == "methane") return methane_row_ids.size();
    if (payload.identity[1] == "ethane") return ethane_row_ids.size();
    if (payload.identity[1] == "propane") return propane_row_ids.size();
    throw std::invalid_argument("component identity must be methane, ethane, or propane");
}

Row parse_row(PyObject* object, const std::string& component_id, std::size_t source_index) {
    OwnedPyObject sequence{PySequence_Fast(object, "source row must be a sequence")};
    if (sequence == nullptr || PySequence_Fast_GET_SIZE(sequence.get()) != 6) {
        throw std::invalid_argument("source row must contain exactly six fields");
    }
    PyObject** items = PySequence_Fast_ITEMS(sequence.get());
    Row row{
        text(items[0], "row_id"),
        text(items[1], "component_id"),
        PyFloat_AsDouble(items[2]),
        PyFloat_AsDouble(items[3]),
        PyFloat_AsDouble(items[4]),
        text(items[5], "source_id"),
    };
    bool matches = false;
    if (component_id == "methane" && source_index < methane_row_ids.size()) {
        matches = row.row_id == methane_row_ids[source_index]
            && row.component_id == "methane"
            && row.temperature == 100.0 + 10.0 * static_cast<double>(source_index)
            && row.pressure == methane_pressures[source_index]
            && row.liquid_density == methane_liquid_densities[source_index]
            && row.source_id == "nist-webbook-srd69-methane-saturation";
    } else if (component_id == "ethane" && source_index < ethane_row_ids.size()) {
        matches = row.row_id == ethane_row_ids[source_index]
            && row.component_id == "ethane"
            && row.temperature == 100.0 + 20.0 * static_cast<double>(source_index)
            && row.pressure == ethane_pressures[source_index]
            && row.liquid_density == ethane_liquid_densities[source_index]
            && row.source_id == "nist-webbook-srd69-ethane-saturation";
    } else if (component_id == "propane" && source_index < propane_row_ids.size()) {
        matches = row.row_id == propane_row_ids[source_index]
            && row.component_id == "propane"
            && row.temperature == 110.0 + 10.0 * static_cast<double>(source_index)
            && row.pressure == propane_pressures[source_index]
            && row.liquid_density == propane_liquid_densities[source_index]
            && row.source_id == "glos-2004-propane-coexistence-experiment";
    }
    if (PyErr_Occurred() != nullptr || !matches) {
        throw std::invalid_argument("source row violates the exact retained identity and values");
    }
    return row;
}

Payload parse_payload(PyObject* object) {
    OwnedPyObject sequence{
        PySequence_Fast(object, "pure-saturation fit payload must be a sequence")
    };
    if (sequence == nullptr) {
        throw std::invalid_argument("pure-saturation fit payload must be a sequence");
    }
    if (PySequence_Fast_GET_SIZE(sequence.get()) != 24) {
        throw std::invalid_argument("pure-saturation fit payload must contain exactly 24 fields");
    }
    PyObject** items = PySequence_Fast_ITEMS(sequence.get());
    Payload payload{};
    payload.identity = texts(items[0], 41, "compiled problem identity");
    if (payload.identity[1] != "methane" && payload.identity[1] != "ethane"
        && payload.identity[1] != "propane") {
        throw std::invalid_argument("component identity must be methane, ethane, or propane");
    }
    const std::vector<std::string> expected = expected_identity(payload.identity[1]);
    if (payload.identity != expected) {
        throw std::invalid_argument("compiled problem identity does not match an admitted component");
    }
    OwnedPyObject rows{PySequence_Fast(items[1], "training rows must be a sequence")};
    if (rows == nullptr
        || PySequence_Fast_GET_SIZE(rows.get()) != static_cast<Py_ssize_t>(row_count)) {
        throw std::invalid_argument("training rows must contain exactly four rows");
    }
    std::array<std::size_t, row_count> source_indices{};
    if (payload.identity[1] == "methane") {
        source_indices = {1, 3, 5, 7};
    } else if (payload.identity[1] == "ethane") {
        source_indices = {2, 4, 6, 8};
    } else {
        source_indices = {4, 10, 16, 22};
    }
    for (std::size_t index = 0; index < row_count; ++index) {
        payload.rows[index] = parse_row(
            PySequence_Fast_GET_ITEM(rows.get(), static_cast<Py_ssize_t>(index)),
            payload.identity[1],
            source_indices[index]
        );
    }
    payload.start = fixed_doubles<3>(items[2], "parameter start");
    payload.lower = fixed_doubles<3>(items[3], "parameter lower bounds");
    payload.upper = fixed_doubles<3>(items[4], "parameter upper bounds");
    payload.parameter_scale = fixed_doubles<3>(items[5], "parameter scales");
    payload.amount = PyFloat_AsDouble(items[6]);
    payload.molar_mass = PyFloat_AsDouble(items[7]);
    payload.weights = fixed_doubles<4>(items[8], "residual weights");
    payload.liquid_volume_bounds = fixed_doubles<2>(items[9], "liquid volume bounds");
    payload.vapor_volume_bounds = fixed_doubles<2>(items[10], "vapor volume bounds");
    payload.topology_separation = PyFloat_AsDouble(items[11]);
    const long max_iterations = PyLong_AsLong(items[12]);
    payload.function_tolerance = PyFloat_AsDouble(items[13]);
    payload.gradient_tolerance = PyFloat_AsDouble(items[14]);
    payload.parameter_tolerance = PyFloat_AsDouble(items[15]);
    payload.reporting_pressure_bounds = fixed_doubles<2>(items[16], "reporting pressure bounds");
    payload.confirmation_liquid_start_multiplier = PyFloat_AsDouble(items[17]);
    payload.confirmation_vapor_start_multiplier = PyFloat_AsDouble(items[18]);
    payload.confirmation_parameter_delta = PyFloat_AsDouble(items[19]);
    payload.confirmation_cost_delta = PyFloat_AsDouble(items[20]);
    payload.reporting_pressure_closure = PyFloat_AsDouble(items[21]);
    payload.reporting_mu_closure = PyFloat_AsDouble(items[22]);
    const long num_threads = PyLong_AsLong(items[23]);
    if (PyErr_Occurred() != nullptr) {
        throw std::invalid_argument("pure-saturation fit payload contains a nonnumeric scalar");
    }
    payload.max_iterations = static_cast<int>(max_iterations);
    payload.num_threads = static_cast<int>(num_threads);
    if (!component_specific_contract_matches(payload)
        || payload.lower != std::array<double, 3>{0.5, 2.0, 50.0}
        || payload.upper != std::array<double, 3>{3.5, 5.0, 400.0}
        || payload.parameter_scale != std::array<double, 3>{0.1, 0.1, 10.0}
        || payload.amount != 1.0
        || payload.weights != std::array<double, 4>{0.25, 0.25, 0.25, 0.25}
        || payload.topology_separation != 1.0e-3
        || payload.function_tolerance != 1.0e-10
        || payload.gradient_tolerance != 1.0e-10
        || payload.parameter_tolerance != 1.0e-10
        || payload.confirmation_liquid_start_multiplier != 1.01
        || payload.confirmation_vapor_start_multiplier != 0.98
        || payload.confirmation_parameter_delta != 1.0e-5
        || payload.confirmation_cost_delta != 1.0e-8
        || payload.reporting_pressure_closure != 1.0e-8
        || payload.reporting_mu_closure != 1.0e-8
        || payload.num_threads != 1) {
        throw std::invalid_argument("pure-saturation fit payload does not match its component contract");
    }
    return payload;
}

}  // namespace epcsaft_engine::regression::internal
