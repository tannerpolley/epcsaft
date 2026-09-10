#include "core/convergence.hpp"
#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"
#include "properties/helmholtz_model.hpp"
#include "properties_support.hpp"
#include "contributions/association.hpp"
#include "recording_support.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>
#include <string>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#define CHECK(condition) recording_support::check( \
    static_cast<bool>(condition), __FILE__, __LINE__, #condition)

namespace {
void core_constant_identities() {
    const double r = epcsaft::AVOGADRO_PER_MOL * epcsaft::BOLTZMANN_J_PER_K;
    const double r_source = 8.31446261815324;
    // R has four literal/product roundings and C_rho has six at most.  The
    // binary64 unit roundoff is epsilon/2, so eight epsilon is a conservative
    // source-scale envelope with room for higher-order terms.
    constexpr double source_roundoff_factor = 8.0;
    const double r_bound =
        source_roundoff_factor * std::numeric_limits<double>::epsilon()
        * std::abs(r_source);
    const double debye_product = epcsaft::DEBYE_C_M * epcsaft::SPEED_OF_LIGHT_M_PER_S;
    CHECK(epcsaft::GAS_CONSTANT_J_PER_MOL_K == r);
    CHECK(std::abs(epcsaft::GAS_CONSTANT_J_PER_MOL_K - r_source) <= r_bound);
    const double epsilon_source = 8.8541878128e-12;
    const double epsilon_bound =
        source_roundoff_factor * std::numeric_limits<double>::epsilon()
        * std::abs(epsilon_source);
    CHECK(std::abs(epcsaft::VACUUM_PERMITTIVITY_F_PER_M - epsilon_source)
        <= epsilon_bound);
    CHECK(std::abs(debye_product - 1.0e-21) <= 1.0e-36);
    CHECK(epcsaft::ANGSTROM_M == 1.0e-10);
    const double density_reference = epcsaft::AVOGADRO_PER_MOL * 1.0e-30;
    const double density_source = 6.02214076e-7;
    const double density_bound =
        source_roundoff_factor * std::numeric_limits<double>::epsilon()
        * std::abs(density_source);
    CHECK(std::abs(epcsaft::NUMBER_DENSITY_PER_ANGSTROM3 - density_reference)
        <= density_bound);
    CHECK(std::abs(epcsaft::NUMBER_DENSITY_PER_ANGSTROM3 - density_source)
        <= density_bound);
    recording_support::compare("core-constants/R", "R", "J/mol/K",
        epcsaft::GAS_CONSTANT_J_PER_MOL_K, r_source, {}, "Relationship", "source-decimal",
        "R=N_A*k_B", "strict_absolute", r_bound, 0.0,
        std::abs(epcsaft::GAS_CONSTANT_J_PER_MOL_K - r_source) <= r_bound,
        "c1-constants", "R-source", "R", "N_A,k_B", "R=8.31446261815324", "constant");
    recording_support::compare("core-constants/R-identity", "R_identity_error", "J/mol/K",
        epcsaft::GAS_CONSTANT_J_PER_MOL_K, r, {}, "Relationship", "SI-definitions",
        "R=N_A*k_B identity", "strict_absolute", r_bound, 0.0,
        std::abs(epcsaft::GAS_CONSTANT_J_PER_MOL_K - r) <= r_bound,
        "c1-constants", "R-identity", "R", "N_A,k_B", "R=N_A*k_B", "constant");
    recording_support::compare("core-constants/debye", "D*c", "C m^2/s",
        debye_product, 1.0e-21, {}, "Relationship", "SI-definitions",
        "D*c=1e-21", "strict_absolute", 1.0e-36, 0.0,
        std::abs(debye_product - 1.0e-21) <= 1.0e-36,
        "c1-constants", "Debye", "D*c", "D,c", "D*c=1e-21", "constant");
    recording_support::compare("core-constants/angstrom", "angstrom", "m",
        epcsaft::ANGSTROM_M, 1.0e-10, {}, "Relationship", "SI-definitions",
        "angstrom=1e-10 m", "absolute", 0.0, 0.0,
        epcsaft::ANGSTROM_M == 1.0e-10,
        "c1-constants", "angstrom", "angstrom", "angstrom", "angstrom=1e-10 m", "constant");
    recording_support::compare("core-constants/epsilon", "epsilon_0", "F/m",
        epcsaft::VACUUM_PERMITTIVITY_F_PER_M, epsilon_source, {},
        "Relationship", "CODATA-2018", "retained CODATA 2018 central value",
        "strict_absolute", epsilon_bound, 0.0,
        std::abs(epcsaft::VACUUM_PERMITTIVITY_F_PER_M - epsilon_source)
            <= epsilon_bound,
        "c1-constants", "epsilon_0", "epsilon_0", "CODATA-2018", "epsilon_0=CODATA-2018", "constant");
    recording_support::compare("core-constants/number-density", "C_rho",
        "(number/angstrom^3)/(mol/m^3)",
        epcsaft::NUMBER_DENSITY_PER_ANGSTROM3, density_source, {},
        "Relationship", "source-decimal", "C_rho=N_A*10^{-30}",
        "strict_absolute", density_bound, 0.0,
        std::abs(epcsaft::NUMBER_DENSITY_PER_ANGSTROM3 - density_source) <= density_bound,
        "c1-constants", "C_rho", "C_rho", "N_A,10^{-30}", "C_rho=6.02214076e-7", "constant");
}

using epcsaft::core::Stop;

int stop_rank(Stop stop) {
    switch (stop) {
    case Stop::converged: return 0;
    case Stop::converged_to_floor: return 1;
    case Stop::budget_exhausted: return 2;
    case Stop::stalled: return 3;
    default: return 4;
    }
}

void classifier_partition() {
    const double v = 2.0;
    const epcsaft::policy::Tolerance criterion{
        "criterion", v, epcsaft::policy::Category::AlgorithmTolerance, "test"
    };
    const std::vector<double> residuals = {0.0, v / 2.0, v, 2.0 * v, 5.0 * v, 20.0 * v, 200.0 * v};
    const std::vector<std::optional<double>> floors = {std::nullopt, v / 10.0, v, 10.0 * v};
    const std::vector<std::optional<double>> contractions = {std::nullopt, 0.5, 1.0, 1.5};
    std::size_t classifications = 0;
    bool all_classifications = true;
    for (std::size_t floor_index = 0; floor_index < floors.size(); ++floor_index) {
        const auto floor = floors[floor_index];
        for (std::size_t contraction_index = 0;
             contraction_index < contractions.size(); ++contraction_index) {
            const auto contraction = contractions[contraction_index];
            for (const bool budget : {false, true}) {
                int previous = -1;
                for (std::size_t residual_index = 0;
                     residual_index < residuals.size(); ++residual_index) {
                    const double residual = residuals[residual_index];
                    const Stop stop = epcsaft::core::classify(
                        residual, criterion, floor, contraction, budget
                    );
                    ++classifications;
                    const bool in_floor = floor.has_value() && residual > v
                        && residual <= epcsaft::policy::FLOOR_MARGIN.value * *floor;
                    const bool exhausted = !((residual <= v) || in_floor) && budget
                        && (!contraction.has_value() || *contraction < 1.0);
                    const int expected = residual <= v ? 0
                        : in_floor ? 1 : exhausted ? 2 : 3;
                    const int actual = stop_rank(stop);
                    const bool verdict = actual == expected && actual >= previous;
                    const auto point = std::to_string(floor_index) + "/"
                        + std::to_string(contraction_index) + "/"
                        + (budget ? "1/" : "0/") + std::to_string(residual_index);
                    const auto floor_value = floor.has_value()
                        ? recording_support::number(*floor) : std::string("none");
                    const auto contraction_value = contraction.has_value()
                        ? recording_support::number(*contraction) : std::string("none");
                    recording_support::compare_custom(
                        (std::string("classification/") + point).c_str(),
                        "stop_rank", "rank", static_cast<double>(actual),
                        static_cast<double>(expected),
                        recording_support::fields({
                            recording_support::field("criterion", v),
                            recording_support::field("floor_margin", epcsaft::policy::FLOOR_MARGIN.value),
                            recording_support::field("residual", residual),
                            recording_support::field("floor", floor_value),
                            recording_support::field("contraction", contraction_value),
                            recording_support::field("budget", budget ? 1U : 0U),
                            recording_support::field("previous_rank", previous),
                            recording_support::field("actual_rank", actual),
                            recording_support::field("expected_rank", expected),
                        }),
                        "Relationship", "partition", "classifier returns the expected stop state",
                        "actual_rank == expected_rank && actual_rank >= previous_rank", verdict,
                        "d1-classification", point.c_str(), "s", "criterion,residual,floor,contraction,budget",
                        "s=expected(s) and s_k\\geq s_{k-1}", "floor,contraction,budget,residual");
                    all_classifications = all_classifications && verdict;
                    previous = actual;
                }
            }
        }
    }
    CHECK(all_classifications);
    CHECK(classifications == 224);
}

void contraction_window() {
    epcsaft::core::ResidualHistory history;
    const double q = 0.7;
    bool all_contractions = true;
    for (std::size_t k = 0; k < 12; ++k) {
        history.push(std::pow(3.0, 1.0) * std::pow(q, static_cast<double>(k)));
        const auto contraction = history.contraction();
        if (k < static_cast<std::size_t>(epcsaft::policy::STALL_WINDOW.value)) {
            all_contractions = all_contractions && !contraction.has_value();
        } else {
            all_contractions = all_contractions && contraction.has_value();
            if (contraction.has_value())
                all_contractions = all_contractions
                    && std::abs(*contraction - q) <= 1e-12 * q;
        }
    }
    const auto contraction = history.contraction();
    const bool final_contraction = contraction.has_value()
        && std::abs(*contraction - q) <= 1e-12 * q;
    recording_support::compare(
        "contraction-limit", "contraction", "1", contraction.value_or(0.0), q,
        recording_support::fields({recording_support::field("q", q),
                                   recording_support::field("history_length", 12U),
                                   recording_support::field("window", epcsaft::policy::STALL_WINDOW.value)}),
        "Relationship", "convergence-history", "contraction estimate reaches q",
        "reference_or_one", 1.0e-12 * q, 0.0,
        final_contraction,
        "d2-contraction", "final", "q", "q,N,window", "q_{estimated}=q", "history");
    CHECK(all_contractions);
    CHECK(final_contraction);
}

void declaration_guards() {
    bool threw = false;
    try {
        const double factor = 0.5;
        (void)epcsaft::policy::derived(
            epcsaft::policy::ASSOCIATION_RESIDUAL, factor, "bad", "test"
        );
    } catch (const char*) { threw = true; }
    CHECK(threw);
    const auto derived = epcsaft::policy::derived(
        epcsaft::policy::ASSOCIATION_RESIDUAL, 3.0, "derived", "test"
    );
    CHECK(derived.value == 3.0 * epcsaft::policy::ASSOCIATION_RESIDUAL.value);
    CHECK(derived.upstream == &epcsaft::policy::ASSOCIATION_RESIDUAL);

    epcsaft::core::Declared<epcsaft::ElectrolyteForm> missing;
    CHECK(!missing.declared());
    threw = false;
    try { missing.require("x"); }
    catch (const epcsaft::core::MissingDeclaration& error) {
        threw = std::string(error.what()).find("x") != std::string::npos;
    }
    CHECK(threw);

    auto aqueous = properties_support::aqueous_electrolyte();
    aqueous.electrolyte_form = {};
    threw = false;
    try { epcsaft::HelmholtzModel(std::move(aqueous)); }
    catch (const std::invalid_argument& error) {
        threw = std::string(error.what()).find("electrolyte_form") != std::string::npos;
    }
    CHECK(threw);
}

void check_density_record(
    epcsaft::HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, epcsaft::DensityBranch branch,
    const char* point_id
) {
    const auto result = model.state_at_pressure(
        temperature, pressure, composition, branch
    );
    CHECK(result.density_convergence.has_value());
    const auto& record = *result.density_convergence;
    const auto returned = model.state(temperature, result.molar_density, composition);
    const double residual = std::abs(
        returned.compressibility_factor
        - pressure / (result.molar_density * epcsaft::GAS_CONSTANT_J_PER_MOL_K * temperature)
    );
    const bool residual_matches = std::abs(record.residual - residual) <= 1e-15;
    recording_support::compare(
        (std::string("density-residual/") + point_id).c_str(),
        "density_residual", "1", record.residual, residual,
        recording_support::fields({recording_support::field("T", temperature),
                                   recording_support::field("P", pressure),
                                   recording_support::field("x", composition),
                                   recording_support::field("rho", result.molar_density),
                                   recording_support::field("branch", branch == epcsaft::DensityBranch::Liquid ? "liquid" : "vapor")}),
        "Relationship", "convergence-record", "recorded residual equals recomputed residual",
        "absolute", 1e-15, 0.0,
        residual_matches, "d4-density-residual", point_id,
        "residual", "T,P,\\mathbf{x},\\rho,branch", "r_{record}=r_{state}", "state");
    CHECK(residual_matches);
    CHECK(record.criterion != nullptr);
    CHECK(std::string(record.criterion->name) == "DENSITY_ROOT_RESIDUAL");
    if (record.stop == Stop::converged) CHECK(record.residual <= record.criterion->value);
}

void density_records() {
    epcsaft::HelmholtzModel ethane(properties_support::ethane());
    for (const auto& row : std::array<std::tuple<double, double>, 2>{
             std::tuple{200.0, 217800.69494692178},
             std::tuple{240.0, 969152.4807171321}}) {
        check_density_record(ethane, std::get<0>(row), std::get<1>(row), {1.0},
                             epcsaft::DensityBranch::Vapor,
                             (std::string("ethane-") + std::to_string(std::get<0>(row))
                              + "-vapor").c_str());
        check_density_record(ethane, std::get<0>(row), std::get<1>(row), {1.0},
                             epcsaft::DensityBranch::Liquid,
                             (std::string("ethane-") + std::to_string(std::get<0>(row))
                              + "-liquid").c_str());
    }
    check_density_record(ethane, 320.0, 1.0e6, {1.0}, epcsaft::DensityBranch::Vapor,
                         "ethane-320-vapor");
    check_density_record(ethane, 320.0, 1.0e6, {1.0}, epcsaft::DensityBranch::Liquid,
                         "ethane-320-liquid");
    epcsaft::HelmholtzModel nine(properties_support::nine_species());
    check_density_record(nine, 313.15, 1.0e5,
                         {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0},
                         epcsaft::DensityBranch::Liquid, "nine-313-liquid");
    epcsaft::HelmholtzModel aqueous(properties_support::aqueous_electrolyte());
    check_density_record(aqueous, 298.15, 1.0e5, {1.0, 0.0, 0.0},
                         epcsaft::DensityBranch::Liquid, "aqueous-298-liquid");
}

void density_publication_boundaries() {
    epcsaft::HelmholtzModel model(properties_support::neutral_binary());
    const auto valid = model.state_at_pressure(
        300.0, 1.0e5, {0.4, 0.6}, epcsaft::DensityBranch::Vapor
    );
    const bool valid_root = valid.density_convergence.has_value()
        && valid.density_convergence->stop == Stop::converged
        && std::isfinite(valid.density_convergence->residual)
        && std::isfinite(valid.molar_density);

    recording_support::expect_group(
        recording_support::active_test_id(), "density-publication",
        "density-publication", 1, "C3 density publication boundaries", "check"
    );
    recording_support::compare(
        "density-publication", "valid-root", "1", valid_root ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("T", 300.0),
            recording_support::field("P", 1.0e5),
            recording_support::field("x", std::vector<double>{0.4, 0.6}),
            recording_support::field("stop", valid.density_convergence.has_value()
                ? epcsaft::core::stop_name(valid.density_convergence->stop) : "absent"),
            recording_support::field("residual", valid.density_convergence.has_value()
                ? valid.density_convergence->residual : std::numeric_limits<double>::infinity()),
        }),
        "Relationship", "density-root", "finite density root remains available",
        "absolute", 0.0, 0.0, valid_root, "native.diagnostics.c3-correction", "C3",
        "density-publication", "availability", "T,P,x", "valid root published"
    );
    CHECK(valid_root);
}

void association_record(
    const epcsaft::Parameters& parameters, double density,
    const std::vector<double>& composition, const char* point_id
) {
    const auto sigma = epcsaft::segment_sigma(300.0, parameters);
    const auto diameter = epcsaft::segment_diameter(300.0, parameters, sigma);
    epcsaft::core::Convergence record;
    const auto fractions = epcsaft::solve_site_fractions(
        300.0, density, composition, parameters, sigma, diameter, &record
    );
    const auto delta = epcsaft::association_delta_values(
        300.0, density, composition, parameters, sigma, diameter,
        parameters.association_energy_over_k, parameters.association_volume
    );
    const auto residuals = epcsaft::association_residuals(
        density, composition, parameters, fractions, delta
    );
    double norm = 0.0;
    for (double value : residuals) norm = std::max(norm, std::abs(value));
    recording_support::compare(
        (std::string("association-residual/") + point_id).c_str(),
        "association_residual", "1", record.residual, norm,
        recording_support::fields({recording_support::field("T", 300.0),
                                   recording_support::field("rho", density),
                                   recording_support::field("x", composition),
                                   recording_support::field("site_fractions", fractions),
                                   recording_support::field("residuals", residuals)}),
        "Relationship", "convergence-record", "association record stores residual norm",
        "absolute", 1e-15, 0.0, std::abs(record.residual - norm) <= 1e-15,
        "d5-association-residual", point_id, "r", "T,\\rho,\\mathbf{x},\\chi",
        "r_{record}=||f(\\chi)||_\\infty", "state");
    CHECK(std::abs(record.residual - norm) <= 1e-15);
}

void association_records() {
    association_record(properties_support::associating_binary(), 1.0e4, {0.1385, 0.8615}, "associating");
    association_record(properties_support::nine_species(), 1.2e4,
                      {0.795, 0.10, 0.05, 0.02, 0.02, 0.005, 0.001, 0.008, 0.001},
                      "nine-species");
}

void same_state(const epcsaft::State& left, const epcsaft::State& right) {
#define SAME(field) CHECK(left.field == right.field)
    SAME(hard_chain); SAME(dispersion); SAME(association); SAME(debye_huckel);
    SAME(born); SAME(bulk_relative_permittivity); SAME(residual_helmholtz);
    SAME(compressibility_factor); SAME(pressure); SAME(residual_entropy);
    SAME(residual_internal_energy); SAME(residual_enthalpy);
    SAME(residual_isochoric_heat_capacity); SAME(molar_density);
    SAME(pressure_density_derivative); SAME(residual_chemical_potential_over_rt);
    SAME(log_fugacity_coefficient); SAME(residual_isobaric_heat_capacity);
    SAME(site_fractions);
    CHECK(left.density_convergence.has_value() == right.density_convergence.has_value());
#undef SAME
}

void trace_case(epcsaft::Parameters parameters, double temperature, double pressure,
                std::vector<double> composition, epcsaft::DensityBranch branch) {
    epcsaft::HelmholtzModel model(std::move(parameters));
    const auto uncaptured = model.state_at_pressure(
        temperature, pressure, composition, branch
    );
    epcsaft::core::Capture capture;
    const auto captured = model.state_at_pressure(
        temperature, pressure, composition, branch
    );
    const auto trace = capture.take();
    same_state(uncaptured, captured);
    bool root_seen = false;
    for (std::size_t index = 0; index < trace.spans.size(); ++index) {
        const auto& span = trace.spans[index];
        if (!span.parent.has_value()) {
            CHECK(!root_seen);
            root_seen = true;
        } else {
            CHECK(*span.parent < index);
        }
    }
    CHECK(root_seen);
}

void trace_state_case(
    epcsaft::Parameters parameters, double temperature, double density,
    std::vector<double> composition
) {
    epcsaft::HelmholtzModel model(std::move(parameters));
    const auto uncaptured = model.state(temperature, density, composition);
    epcsaft::core::Capture capture;
    const auto captured = model.state(temperature, density, composition);
    const auto trace = capture.take();
    same_state(uncaptured, captured);
    CHECK(!trace.spans.empty());
    CHECK(std::string(trace.spans.front().name) == "state");
}

void trace_activity_case() {
    epcsaft::HelmholtzModel model(properties_support::aqueous_electrolyte());
    const std::vector<double> solvent = {1.0, 0.0, 0.0};
    const std::array<std::size_t, 2> ions = {1, 2};
    const std::array<int, 2> stoichiometry = {1, 1};
    const auto reference = epcsaft::infinite_dilution_reference(
        model, 298.15, 1.0e5, solvent, ions, stoichiometry
    );
    const auto uncaptured = epcsaft::electrolyte_activity(
        model, 298.15, 1.0e5, solvent, ions, stoichiometry, 1.0e-8, reference
    );
    epcsaft::core::Capture capture;
    const auto captured = epcsaft::electrolyte_activity(
        model, 298.15, 1.0e5, solvent, ions, stoichiometry, 1.0e-8, reference
    );
    const auto trace = capture.take();
    CHECK(uncaptured.single_ion == captured.single_ion);
    CHECK(uncaptured.mean_ionic_log_molality_scale
          == captured.mean_ionic_log_molality_scale);
    CHECK(!trace.spans.empty());
    CHECK(std::string(trace.spans.front().name) == "electrolyte_activity");
}

void trace_invariance() {
    trace_case(properties_support::neutral_binary(), 300.0, 1.0e5, {0.4, 0.6},
               epcsaft::DensityBranch::Vapor);
    trace_case(properties_support::ethane(), 200.0, 217800.69494692178, {1.0},
               epcsaft::DensityBranch::Liquid);
    trace_case(properties_support::ethane(), 320.0, 1.0e6, {1.0},
               epcsaft::DensityBranch::Vapor);
    trace_state_case(properties_support::neutral_binary(), 300.0, 1.0, {0.4, 0.6});
    trace_activity_case();
}

void termination_boundaries() {
    const epcsaft::policy::Tolerance criterion{"c3", 1.0, epcsaft::policy::Category::AlgorithmTolerance, "test"};
    const std::vector<double> nonfinite = {
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity(),
        -std::numeric_limits<double>::infinity()
    };
    bool valid = true;
    for (const double residual : nonfinite) {
        valid = valid && epcsaft::core::classify(residual, criterion, std::nullopt,
            std::nullopt, false) == Stop::nonfinite;
    }
    valid = valid && epcsaft::core::classify(2.0, criterion,
        std::numeric_limits<double>::infinity(), std::nullopt, false) == Stop::nonfinite;
    valid = valid && epcsaft::core::classify(2.0, criterion, std::nullopt,
        std::numeric_limits<double>::infinity(), false) == Stop::nonfinite;
    valid = valid && epcsaft::core::classify(2.0, criterion, std::nullopt,
        std::nullopt, false, std::numeric_limits<double>::infinity()) == Stop::nonfinite;
    epcsaft::core::ResidualHistory history;
    bool history_rejected = false;
    try { history.push(std::numeric_limits<double>::quiet_NaN()); }
    catch (const std::invalid_argument&) { history_rejected = true; }
    valid = valid && history_rejected;
    recording_support::expect_group(recording_support::active_test_id(),
        "termination-boundaries", "termination-boundaries", 1,
        "C3 nonfinite termination", "check");
    recording_support::compare("termination-boundaries", "valid", "1",
        valid ? 1.0 : 0.0, 1.0, {}, "Relationship", "classifier", "nonfinite boundaries",
        "absolute", 0.0, 0.0, valid, "native.diagnostics.c3", "C3",
        "termination", "stop", "residual", "nonfinite inputs rejected");
}

void producer_boundaries() {
    const auto parameters = properties_support::associating_binary();
    const auto sigma = epcsaft::segment_sigma(300.0, parameters);
    const auto diameter = epcsaft::segment_diameter(300.0, parameters, sigma);

    auto nonfinite_parameters = parameters;
    nonfinite_parameters.association_energy_over_k[1] = 1.0e308;
    epcsaft::core::Convergence nonfinite_record;
    bool nonfinite = false;
    try {
        (void)epcsaft::solve_site_fractions(
            300.0, 1.0e4, {0.4, 0.6}, nonfinite_parameters, sigma, diameter,
            &nonfinite_record);
    } catch (const epcsaft::core::ConvergenceFailure& error) {
        nonfinite_record = error.record;
        nonfinite = error.record.stop == Stop::nonfinite;
    }

    epcsaft::core::Convergence domain_record;
    bool domain = false;
    try {
        (void)epcsaft::solve_site_fractions(
            300.0, -1.0e8, {0.4, 0.6}, parameters, sigma, diameter,
            &domain_record);
    } catch (const epcsaft::core::ConvergenceFailure& error) {
        domain_record = error.record;
        domain = error.record.stop == Stop::left_domain;
    }

    const bool verdict = nonfinite && domain
        && !std::isfinite(nonfinite_record.residual)
        && std::isinf(domain_record.residual);
    recording_support::expect_group(
        recording_support::active_test_id(), "producer-boundaries", "producer-boundaries", 1,
        "C3 typed producer stops", "check");
    recording_support::compare(
        "producer-boundaries", "typed-stops", "1", verdict ? 1.0 : 0.0, 1.0,
        recording_support::fields({
            recording_support::field("nonfinite_stop",
                                     epcsaft::core::stop_name(nonfinite_record.stop)),
            recording_support::field("domain_stop",
                                     epcsaft::core::stop_name(domain_record.stop)),
            recording_support::field("nonfinite_residual", nonfinite_record.residual),
            recording_support::field("domain_residual", domain_record.residual),
            recording_support::field("T", 300.0),
            recording_support::field("rho_nonfinite", 1.0e4),
            recording_support::field("rho_domain", -1.0e8),
        }),
        "Relationship", "producer-failure", "producers retain typed stop reasons",
        "absolute", 0.0, 0.0, verdict, "native.diagnostics.c3", "C3",
        "producer-boundaries", "stop", "T,rho,x,parameters", "typed stop retained");
}

}  // namespace

int main(int argc, char** argv) {
    return recording_support::run(argc, argv, {
        {{"native.diagnostics.c1", "Core_constant_identities", "Numerical", "Core", "Explicit", "Numerical", nullptr, "CORE1", "core-constants/*", "engine/tests/native/diagnostics_checks.cpp", "diagnostics_checks", "C1"}, core_constant_identities},
        {{"native.diagnostics.d1", "Classifier_partition", "Numerical", "Core", "Explicit", "Numerical", nullptr, "CORE3", "classification/*", "engine/tests/native/diagnostics_checks.cpp", "diagnostics_checks", "D1", "classification/*:224"}, classifier_partition},
        {{"native.diagnostics.d2", "Contraction_window", "Numerical", "Core", "Explicit", "Numerical", nullptr, "CORE3", "contraction-limit", "engine/tests/native/diagnostics_checks.cpp", "diagnostics_checks", "D2"}, contraction_window},
        {{"native.diagnostics.d3", "Declaration_guards", "Package", "Surface", "Explicit", "Outcome_only", nullptr, "S1", nullptr, "engine/tests/native/diagnostics_checks.cpp", "diagnostics_checks", "D3"}, declaration_guards},
        {{"native.diagnostics.d4", "Density_convergence_records", "Numerical", "Properties", "Explicit", "Numerical", "Residual", "P2", "density-residual/*", "engine/tests/native/diagnostics_checks.cpp", "diagnostics_checks", "D4", "density-residual/*:8"}, density_records},
        {{"native.diagnostics.c3-density", "Density_publication_boundaries", "Numerical", "Properties", "Explicit", "Numerical", nullptr, "CORE3", "density-publication", "engine/tests/native/diagnostics_checks.cpp", "diagnostics_checks", "C3", "density-publication:1"}, density_publication_boundaries},
        {{"native.diagnostics.d5", "Association_convergence_records", "Numerical", "Contributions", "Explicit", "Numerical", "Residual", "C2", "association-residual/*", "engine/tests/native/diagnostics_checks.cpp", "diagnostics_checks", "D5", "association-residual/*:2"}, association_records},
        {{"native.diagnostics.c3", "Termination_boundaries", "Numerical", "Core", "Explicit", "Numerical", nullptr, "CORE3", "termination-boundaries", "engine/tests/native/diagnostics_checks.cpp", "diagnostics_checks", "C3", "termination-boundaries:1"}, termination_boundaries},
        {{"native.diagnostics.c3-producers", "Producer_boundaries", "Numerical", "Core", "Explicit", "Numerical", nullptr, "CORE3", "producer-boundaries", "engine/tests/native/diagnostics_checks.cpp", "diagnostics_checks", "C3", "producer-boundaries:1"}, producer_boundaries},
        {{"native.diagnostics.d6", "Trace_invariance", "Package", "Core", "Explicit", "Outcome_only", nullptr, "B1", nullptr, "engine/tests/native/diagnostics_checks.cpp", "diagnostics_checks", "D6"}, trace_invariance},
    });
}
