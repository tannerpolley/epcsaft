#include "core/convergence.hpp"
#include "core/numerical_policy.hpp"
#include "model/helmholtz_model.hpp"
#include "model_support.hpp"
#include "physics/association.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <optional>
#include <string>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#define CHECK(condition) do { if (!(condition)) { std::fprintf(stderr, "%s:%d: failed: %s\n", __FILE__, __LINE__, #condition); std::abort(); } } while (false)

namespace {

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
    for (const auto floor : floors) for (const auto contraction : contractions) {
        for (bool budget : {false, true}) {
            int previous = -1;
            for (const double residual : residuals) {
                const Stop stop = epcsaft::core::classify(
                    residual, criterion, floor, contraction, budget
                );
                CHECK(stop == Stop::converged || stop == Stop::converged_to_floor
                      || stop == Stop::budget_exhausted || stop == Stop::stalled);
                CHECK((stop == Stop::converged) == (residual <= v));
                const bool in_floor = floor.has_value() && residual > v
                    && residual <= epcsaft::policy::FLOOR_MARGIN.value * *floor;
                CHECK((stop == Stop::converged_to_floor) == in_floor);
                const bool exhausted = !((residual <= v) || in_floor) && budget
                    && (!contraction.has_value() || *contraction < 1.0);
                CHECK((stop == Stop::budget_exhausted) == exhausted);
                CHECK(stop_rank(stop) >= previous);
                previous = stop_rank(stop);
            }
        }
    }
}

void contraction_window() {
    epcsaft::core::ResidualHistory history;
    const double q = 0.7;
    for (std::size_t k = 0; k < 12; ++k) {
        history.push(std::pow(3.0, 1.0) * std::pow(q, static_cast<double>(k)));
        const auto contraction = history.contraction();
        if (k < static_cast<std::size_t>(epcsaft::policy::STALL_WINDOW.value)) {
            CHECK(!contraction.has_value());
        } else {
            CHECK(contraction.has_value());
            CHECK(std::abs(*contraction - q) <= 1e-12 * q);
        }
    }
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
    CHECK(!missing.is(epcsaft::ElectrolyteForm::None));
    threw = false;
    try { missing.require("x"); }
    catch (const epcsaft::core::MissingDeclaration& error) {
        threw = std::string(error.what()).find("x") != std::string::npos;
    }
    CHECK(threw);

    auto aqueous = model_support::aqueous_electrolyte();
    aqueous.electrolyte_form = {};
    threw = false;
    try { epcsaft::HelmholtzModel(std::move(aqueous)); }
    catch (const std::invalid_argument& error) {
        threw = std::string(error.what()).find("electrolyte_form") != std::string::npos;
    }
    CHECK(threw);
    auto neutral = model_support::neutral_binary();
    neutral.electrolyte_form = epcsaft::ElectrolyteForm::None;
    threw = false;
    try { epcsaft::HelmholtzModel(std::move(neutral)); }
    catch (const std::invalid_argument& error) {
        threw = std::string(error.what()).find("electrolyte_form") != std::string::npos;
    }
    CHECK(threw);
}

void check_density_record(
    epcsaft::HelmholtzModel& model, double temperature, double pressure,
    const std::vector<double>& composition, epcsaft::DensityBranch branch
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
    CHECK(std::abs(record.residual - residual) <= 1e-15);
    CHECK(record.criterion != nullptr);
    CHECK(std::string(record.criterion->name) == "DENSITY_ROOT_RESIDUAL");
    if (record.stop == Stop::converged) CHECK(record.residual <= record.criterion->value);
}

void density_records() {
    epcsaft::HelmholtzModel ethane(model_support::ethane());
    for (const auto& row : std::array<std::tuple<double, double>, 2>{
             std::tuple{200.0, 217800.69494692178},
             std::tuple{240.0, 969152.4807171321}}) {
        check_density_record(ethane, std::get<0>(row), std::get<1>(row), {1.0},
                             epcsaft::DensityBranch::Vapor);
        check_density_record(ethane, std::get<0>(row), std::get<1>(row), {1.0},
                             epcsaft::DensityBranch::Liquid);
    }
    check_density_record(ethane, 320.0, 1.0e6, {1.0}, epcsaft::DensityBranch::Vapor);
    check_density_record(ethane, 320.0, 1.0e6, {1.0}, epcsaft::DensityBranch::Liquid);
    epcsaft::HelmholtzModel nine(model_support::nine_species());
    check_density_record(nine, 313.15, 1.0e5,
                         {0.80, 0.10, 0.05, 0.02, 0.02, 0.005, 0.0, 0.005, 0.0},
                         epcsaft::DensityBranch::Liquid);
    epcsaft::HelmholtzModel aqueous(model_support::aqueous_electrolyte());
    check_density_record(aqueous, 298.15, 1.0e5, {1.0, 0.0, 0.0},
                         epcsaft::DensityBranch::Liquid);
}

void association_record(
    const epcsaft::Parameters& parameters, double density,
    const std::vector<double>& composition
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
    CHECK(std::abs(record.residual - norm) <= 1e-15);
}

void association_records() {
    association_record(model_support::associating_binary(), 1.0e4, {0.1385, 0.8615});
    association_record(model_support::nine_species(), 1.2e4,
                      {0.795, 0.10, 0.05, 0.02, 0.02, 0.005, 0.001, 0.008, 0.001});
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

void trace_invariance() {
    trace_case(model_support::neutral_binary(), 300.0, 1.0e5, {0.4, 0.6},
               epcsaft::DensityBranch::Vapor);
    trace_case(model_support::ethane(), 200.0, 217800.69494692178, {1.0},
               epcsaft::DensityBranch::Liquid);
    trace_case(model_support::ethane(), 320.0, 1.0e6, {1.0},
               epcsaft::DensityBranch::Vapor);
}

}  // namespace

int main() {
    classifier_partition();
    contraction_window();
    declaration_guards();
    density_records();
    association_records();
    trace_invariance();
    std::printf("diagnostics checks passed\n");
}
