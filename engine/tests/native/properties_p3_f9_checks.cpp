// P3, Issue26, energy and Born actions, F8, F9, reference limits and residual cv.
#include "equilibrium_support.hpp"
#include "properties_checks_shared.hpp"

namespace properties_checks {

std::string p3_inputs(const char* form, double temperature = 320.0) {
    return recording_support::fields({
        recording_support::field("form", form),
        recording_support::field("T", temperature),
        recording_support::field("rho", 1000.0),
        recording_support::field("P0", 1.0e5),
    });
}

void p3_compare(const std::string& id, const char* quantity, const char* unit,
                double model, double reference, const std::string& inputs,
                const char* group, const char* point, double tolerance = 2.0e-12) {
    const bool ok = close(model, reference, tolerance);
    recording_support::compare(
        id.c_str(), quantity, unit, model, reference, inputs, "Relationship",
        "independent-reference", "P3 property/reference relation", "reference_or_one",
        0.0, tolerance, ok, group, point, quantity, "T,rho,P0", "model=reference", "state"
    );
}

void properties_thermochemistry_p3() {
    const double R = epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const double temperature = 320.0;
    const double density = 1000.0;
    const double standard_entropy = 50.0 + 30.0 * std::log(temperature / 300.0);
    const double expected_h = 1000.0 + 30.0 * (temperature - 300.0);
    const double expected_u = expected_h - R * temperature;
    const double expected_s = standard_entropy
        - R * std::log(density * R * temperature / 1.0e5);
    const auto forms = p3_forms();


    for (const auto& [name, correlation] : forms) {
        const auto record = p3_record(correlation);
        const auto ideal = epcsaft::evaluate_ideal_properties(record, temperature, density, {1.0});
        CHECK(ideal.value.has_value() && ideal.status == epcsaft::ThermochemistryStatus::Available);
        const auto inputs = p3_inputs(name);
        p3_compare(std::string("p3-form-value/") + name + "/h", "ideal_enthalpy", "J/mol",
                   ideal.value->enthalpy, expected_h, inputs, "p3-form-values", name);
        p3_compare(std::string("p3-form-value/") + name + "/u", "ideal_internal_energy", "J/mol",
                   ideal.value->internal_energy, expected_u, inputs, "p3-form-values", name);
        p3_compare(std::string("p3-form-value/") + name + "/s", "ideal_entropy", "J/(mol K)",
                   ideal.value->entropy, expected_s, inputs, "p3-form-values", name);
        p3_compare(std::string("p3-form-value/") + name + "/cp", "ideal_heat_capacity_cp", "J/(mol K)",
                   ideal.value->isobaric_heat_capacity, 30.0, inputs, "p3-form-values", name);
        p3_compare(std::string("p3-form-value/") + name + "/cv", "ideal_heat_capacity_cv", "J/(mol K)",
                   ideal.value->isochoric_heat_capacity, 30.0 - R, inputs, "p3-form-values", name);


        epcsaft::HelmholtzModel model(properties_support::ethane(), {}, record);
        const auto state = model.state(temperature, density, {1.0});
        CHECK(state.thermochemistry_status == epcsaft::ThermochemistryStatus::Available);
        CHECK(state.ideal_enthalpy.has_value() && state.ideal_internal_energy.has_value()
              && state.ideal_entropy.has_value() && state.ideal_isobaric_heat_capacity.has_value()
              && state.ideal_isochoric_heat_capacity.has_value()
              && state.total_enthalpy.has_value() && state.total_internal_energy.has_value()
              && state.total_entropy.has_value() && state.total_isochoric_heat_capacity.has_value());
        // Recorded for the Python transport rows; total = ideal + residual is the Engine's own definition.
        for (const auto& [key, unit, value] : std::array<std::tuple<const char*, const char*, double>, 5>{{
                 {"h", "J/mol", *state.total_enthalpy}, {"u", "J/mol", *state.total_internal_energy},
                 {"s", "J/(mol K)", *state.total_entropy}, {"cp", "J/(mol K)", *state.total_isobaric_heat_capacity},
                 {"cv", "J/(mol K)", *state.total_isochoric_heat_capacity}}})
            recording_support::observe((std::string("p3-total/") + name + "/" + key).c_str(), key, unit, value,
                                       std::nullopt, inputs, "Relationship", "state-total", "P3 total property",
                                       "p3-total", name);
    }

    for (const auto& [name, correlation] : forms) {
        for (double endpoint_temperature : {250.0, 350.0}) {
            const auto endpoint_state = epcsaft::evaluate_ideal_properties(
                p3_record(correlation), endpoint_temperature, density, {1.0}
            );
            CHECK(endpoint_state.value.has_value());
            const double h = 1000.0 + 30.0 * (endpoint_temperature - 300.0);
            const double u = h - R * endpoint_temperature;
            const double s = 50.0 + 30.0 * std::log(endpoint_temperature / 300.0)
                - R * std::log(density * R * endpoint_temperature / 1.0e5);
            const auto point = std::string(name) + "/" + recording_support::number(endpoint_temperature);
            const auto endpoint_inputs = p3_inputs(name, endpoint_temperature);
            p3_compare(std::string("p3-endpoint/") + point + "/h", "ideal_enthalpy", "J/mol",
                       endpoint_state.value->enthalpy, h, endpoint_inputs, "p3-endpoints", point.c_str());
            p3_compare(std::string("p3-endpoint/") + point + "/u", "ideal_internal_energy", "J/mol",
                       endpoint_state.value->internal_energy, u, endpoint_inputs, "p3-endpoints", point.c_str());
            p3_compare(std::string("p3-endpoint/") + point + "/s", "ideal_entropy", "J/(mol K)",
                       endpoint_state.value->entropy, s, endpoint_inputs, "p3-endpoints", point.c_str());
            p3_compare(std::string("p3-endpoint/") + point + "/cp", "ideal_heat_capacity_cp", "J/(mol K)",
                       endpoint_state.value->isobaric_heat_capacity, 30.0, endpoint_inputs, "p3-endpoints", point.c_str());
            p3_compare(std::string("p3-endpoint/") + point + "/cv", "ideal_heat_capacity_cv", "J/(mol K)",
                       endpoint_state.value->isochoric_heat_capacity, 30.0 - R, endpoint_inputs, "p3-endpoints", point.c_str());
        }
    }

    const auto binary = epcsaft::evaluate_ideal_properties(
        p3_binary_record(), temperature, density, {0.25, 0.75}
    );
    CHECK(binary.value.has_value());
    const double binary_h = 0.25 * expected_h + 0.75 * (2000.0 + 40.0 * (temperature - 300.0));
    const double binary_u = binary_h - R * temperature;
    const double binary_cp = 0.25 * 30.0 + 0.75 * 40.0;
    const double binary_s = 0.25 * (standard_entropy - R * std::log(
        0.25 * density * R * temperature / 1.0e5))
        + 0.75 * (70.0 + 40.0 * std::log(temperature / 300.0) - R * std::log(
            0.75 * density * R * temperature / 1.0e5));
    p3_compare("p3-mixing/h", "ideal_enthalpy", "J/mol", binary.value->enthalpy, binary_h,
               "x=0.25,0.75,T=320", "p3-mixing", "h");
    p3_compare("p3-mixing/u", "ideal_internal_energy", "J/mol", binary.value->internal_energy, binary_u,
               "x=0.25,0.75,T=320", "p3-mixing", "u");
    p3_compare("p3-mixing/s", "ideal_entropy", "J/(mol K)", binary.value->entropy, binary_s,
               "x=0.25,0.75,T=320", "p3-mixing", "s");
    p3_compare("p3-mixing/cp", "ideal_heat_capacity_cp", "J/(mol K)", binary.value->isobaric_heat_capacity,
               binary_cp, "x=0.25,0.75,T=320", "p3-mixing", "cp");
    p3_compare("p3-mixing/cv", "ideal_heat_capacity_cv", "J/(mol K)", binary.value->isochoric_heat_capacity,
               binary_cp - R, "x=0.25,0.75,T=320", "p3-mixing", "cv");

    const auto endpoint_record = p3_record(forms.front().second);
    const auto outside = epcsaft::HelmholtzModel(
        properties_support::ethane(), {}, endpoint_record
    ).state(400.0, density, {1.0});
    const bool outside_partial = !outside.ideal_enthalpy.has_value()
        && std::isfinite(outside.residual_enthalpy);
    recording_support::compare_custom("p3-boundary/outside", "residual_only", "bool",
        outside_partial ? 1.0 : 0.0, 1.0, "T=400", "Availability", "interval-owner",
        "outside interval preserves residual-only state", "absolute", outside_partial,
        "p3-boundaries", "outside");
    bool gap_failed = false, overlap_failed = false, pressure_failed = false;
    try {
        auto gap = endpoint_record;
        gap.component_intervals.front().push_back(
            epcsaft::IdealInterval{351.0, 360.0, true, true, forms.front().second}
        );
        epcsaft::validate_thermochemistry(gap, 1);
    } catch (const std::invalid_argument&) { gap_failed = true; }
    try {
        auto overlap = endpoint_record;
        overlap.component_intervals.front().push_back(
            epcsaft::IdealInterval{250.0, 350.0, true, true, forms.front().second}
        );
        epcsaft::validate_thermochemistry(overlap, 1);
    } catch (const std::invalid_argument&) { overlap_failed = true; }
    try {
        auto invalid = endpoint_record;
        invalid.component_intervals.front().front().correlation.source_pressure = 0.0;
        epcsaft::validate_thermochemistry(invalid, 1);
    } catch (const std::invalid_argument&) { pressure_failed = true; }
    recording_support::compare_custom("p3-boundary/gap", "gap_rejected", "bool",
        gap_failed ? 1.0 : 0.0, 1.0, "351-360 after 250-350", "Availability", "interval-owner",
        "gapped intervals are rejected", "absolute", gap_failed, "p3-boundaries", "gap");
    recording_support::compare_custom("p3-boundary/overlap", "overlap_rejected", "bool",
        overlap_failed ? 1.0 : 0.0, 1.0, "250-350 plus 250-350", "Availability", "interval-owner",
        "overlapping intervals are rejected", "absolute", overlap_failed, "p3-boundaries", "overlap");
    recording_support::compare_custom("p3-boundary/pressure", "invalid_pressure_rejected", "bool",
        pressure_failed ? 1.0 : 0.0, 1.0, "Psrc=0", "Availability", "input-validation",
        "invalid source pressure is rejected", "absolute", pressure_failed, "p3-boundaries", "pressure");
    bool missing_component_failed = false;
    try {
        auto missing = endpoint_record;
        missing.component_intervals.clear();
        epcsaft::validate_thermochemistry(missing, 1);
    } catch (const std::invalid_argument&) { missing_component_failed = true; }
    recording_support::compare_custom("p3-boundary/missing-component", "missing_component_rejected", "bool",
        missing_component_failed ? 1.0 : 0.0, 1.0, "component thermal record absent", "Availability",
        "input-validation", "missing component thermochemistry is unavailable", "absolute",
        missing_component_failed, "p3-boundaries", "missing-component");

    const auto missing_record_state = epcsaft::HelmholtzModel(properties_support::ethane()).state(
        temperature, density, {1.0}
    );
    recording_support::compare_custom("p3-availability/missing-record", "thermochemistry_status", "enum",
        static_cast<double>(missing_record_state.thermochemistry_status), 0.0, "no record",
        "Availability", "status-owner", "missing record has a typed cause", "absolute",
        missing_record_state.thermochemistry_status == epcsaft::ThermochemistryStatus::MissingRecord,
        "p3-availability", "missing-record");
    recording_support::compare_custom("p3-availability/outside", "thermochemistry_status", "enum",
        static_cast<double>(outside.thermochemistry_status), 2.0, "T=400 outside [250,350]",
        "Availability", "status-owner", "outside interval has a typed cause", "absolute",
        outside.thermochemistry_status == epcsaft::ThermochemistryStatus::OutsideInterval,
        "p3-availability", "outside");
    auto nonfinite_record = endpoint_record;
    std::get<epcsaft::IdealPolynomial>(
        nonfinite_record.component_intervals.front().front().correlation.form
    ).coefficients[0] = 1.0e308;
    const auto nonfinite = epcsaft::HelmholtzModel(
        properties_support::ethane(), {}, nonfinite_record
    ).state(temperature, density, {1.0});
    recording_support::compare_custom("p3-availability/nonfinite", "thermochemistry_status", "enum",
        static_cast<double>(nonfinite.thermochemistry_status), 3.0, "finite record with nonfinite evaluation",
        "Availability", "status-owner", "nonfinite evaluation has a typed cause", "absolute",
        nonfinite.thermochemistry_status == epcsaft::ThermochemistryStatus::Nonfinite,
        "p3-availability", "nonfinite");
    const auto negative_d = epcsaft::HelmholtzModel(
        properties_support::ethane(), {}, endpoint_record
    ).state(280.0, 4000.0, {1.0});
    const bool negative_d_partial = negative_d.thermochemistry_status == epcsaft::ThermochemistryStatus::Available
        && !negative_d.total_isobaric_heat_capacity.has_value()
        && std::isfinite(negative_d.residual_enthalpy);
    recording_support::compare_custom("p3-availability/negative-D", "partial_state", "bool",
        negative_d_partial ? 1.0 : 0.0, 1.0, "T=280,rho=4000,D<0",
        "Availability", "partial-state", "negative-D preserves other fields and total cp unavailability",
        "absolute", negative_d_partial, "p3-availability", "negative-D");
    const auto available_state = epcsaft::HelmholtzModel(
        properties_support::ethane(), {}, endpoint_record
    ).state(temperature, density, {1.0});
    recording_support::compare_custom("p3-availability/available", "thermochemistry_status", "enum",
        static_cast<double>(available_state.thermochemistry_status), 1.0, "T=320 inside [250,350]",
        "Availability", "status-owner", "available record has a typed cause", "absolute",
        available_state.thermochemistry_status == epcsaft::ThermochemistryStatus::Available,
        "p3-availability", "available");

    auto high_pressure = endpoint_record;
    high_pressure.standard_pressure = 2.0e5;
    const auto low_state = epcsaft::evaluate_ideal_properties(endpoint_record, temperature, density, {1.0});
    const auto high_ideal = epcsaft::evaluate_ideal_properties(high_pressure, temperature, density, {1.0});
    CHECK(low_state.value.has_value() && high_ideal.value.has_value());
    p3_compare("p3-reference/h", "h_pressure_invariance", "J/mol",
               high_ideal.value->enthalpy, low_state.value->enthalpy, "P0=1e5 to 2e5",
               "p3-reference", "h");
    p3_compare("p3-reference/u", "u_pressure_invariance", "J/mol",
               high_ideal.value->internal_energy, low_state.value->internal_energy, "P0=1e5 to 2e5",
               "p3-reference", "u");
    p3_compare("p3-reference/s", "s_pressure_invariance", "J/(mol K)",
               high_ideal.value->entropy, low_state.value->entropy, "P0=1e5 to 2e5",
               "p3-reference", "s");
    p3_compare("p3-reference/cp", "cp_pressure_invariance", "J/(mol K)",
               high_ideal.value->isobaric_heat_capacity, low_state.value->isobaric_heat_capacity,
               "P0=1e5 to 2e5", "p3-reference", "cp");
    p3_compare("p3-reference/cv", "cv_pressure_invariance", "J/(mol K)",
               high_ideal.value->isochoric_heat_capacity, low_state.value->isochoric_heat_capacity,
               "P0=1e5 to 2e5", "p3-reference", "cv");
    const auto low_total = epcsaft::HelmholtzModel(
        properties_support::ethane(), {}, endpoint_record
    ).state(temperature, density, {1.0});
    const auto high_total = epcsaft::HelmholtzModel(
        properties_support::ethane(), {}, high_pressure
    ).state(temperature, density, {1.0});
    for (const auto& item : std::array<std::tuple<const char*, double, double, const char*>, 5>{
        std::tuple{"h", *low_total.total_enthalpy, *high_total.total_enthalpy, "J/mol"},
        {"u", *low_total.total_internal_energy, *high_total.total_internal_energy, "J/mol"},
        {"s", *low_total.total_entropy, *high_total.total_entropy, "J/(mol K)"},
        {"cp", *low_total.total_isobaric_heat_capacity, *high_total.total_isobaric_heat_capacity, "J/(mol K)"},
        {"cv", *low_total.total_isochoric_heat_capacity, *high_total.total_isochoric_heat_capacity, "J/(mol K)"}
    }) {
        const auto [quantity, low, high, unit] = item;
        const auto point = std::string("total-") + quantity;
        p3_compare(std::string("p3-reference/") + point, "total_pressure_invariance", unit,
                   high, low, "State P0=1e5 to 2e5", "p3-reference", point.c_str());
    }

    const epcsaft::ReactionLogPolynomial correlation{
        2.0, -1200.0, 0.4, -0.003, 300.0, 293.15, 393.15, "P0"
    };

    // PHREEQC v3.9.0 WATEQ4F, H+ + CO3-2 = HCO3-; pinned LK_SPECIES replay (#790).
    const double ln10 = std::log(10.0), tref_phreeqc = 298.15;
    const epcsaft::ReactionLogPolynomial bicarbonate{ln10 * 107.8871 - 38.92561 * std::log(tref_phreeqc),
        ln10 * -5151.79, -38.92561, ln10 * 0.03252849, tref_phreeqc, 298.15, 423.15,
        "PHREEQC-WATEQ4F", ln10 * 563713.9, 0.0};
    const std::array<double, 6> phreeqc{10.32885437851, 10.17440127629, 10.12525183557,
        10.15523535789, 10.24662341286, 10.38703684093};
    for (std::size_t i = 0; i < phreeqc.size(); ++i) {
        const double t = 298.15 + 25.0 * i;
        const double value = epcsaft::evaluate_reaction_properties(bicarbonate, t).ln_k / ln10;
        recording_support::compare(("p3-reaction/bicarbonate/" + std::to_string(i)).c_str(), "log10_K", "1",
            value, phreeqc[i], "T=" + std::to_string(t), "Reference", "PHREEQC-LK_SPECIES",
            "pinned WATEQ4F bicarbonate source law", "absolute", 1e-10, 0.0,
            std::abs(value - phreeqc[i]) <= 1e-10, "p3-reaction", "bicarbonate");
    }
    for (const auto terms : {std::pair{90000.0, 0.0}, std::pair{0.0, 1e-6}, std::pair{90000.0, 1e-6}}) {
        const epcsaft::ReactionLogPolynomial law{0, 0, 0, 0, 350, 300, 400, "basis-verification", terms.first, terms.second};
        for (std::size_t order = 1; order <= 3; ++order) {
            int successive = 0;
            for (const double h : {0.01, 0.005, 0.0025}) {
                const double exact = epcsaft::reaction_log_value(law, 350.0, order);
                const double fd = (epcsaft::reaction_log_value(law, 350.0 + h, order - 1)
                    - epcsaft::reaction_log_value(law, 350.0 - h, order - 1)) / (2.0 * h);
                const bool ok = std::abs(exact - fd) <= 1e-7 * std::abs(fd) + 1e-12;
                successive = ok ? successive + 1 : 0;
                recording_support::compare(("p3-reaction/terms/" + std::to_string(terms.first) + "/"
                    + std::to_string(terms.second) + "/" + std::to_string(order) + "/" + std::to_string(h)).c_str(),
                    "dlnK_dT_order", ("K^-" + std::to_string(order)).c_str(), exact, fd,
                    "T=350;h=" + std::to_string(h), "Derivative", "central-difference",
                    "difference of preceding native derivative", "reference_scale", 1e-12, 1e-7, ok,
                    "p3-reaction", "basis-terms");
            }
            CHECK(successive >= 2);
        }
    }
    for (const auto member : {&epcsaft::ReactionLogPolynomial::e, &epcsaft::ReactionLogPolynomial::f}) {
        auto invalid = correlation;
        invalid.*member = std::numeric_limits<double>::infinity();
        bool refused = false;
        try { (void)epcsaft::evaluate_reaction_properties(invalid, 350.0); }
        catch (const std::invalid_argument&) { refused = true; }
        CHECK(refused);
    }

    const auto reaction = epcsaft::evaluate_reaction_properties(correlation, 350.0);
    const double expected_ln_k = 2.0 - 1200.0 / 350.0
        + 0.4 * std::log(350.0 / 300.0) - 0.003 * 350.0;
    const double expected_derivative = 1200.0 / (350.0 * 350.0) + 0.4 / 350.0 - 0.003;
    p3_compare("p3-reaction/ln-k", "ln_k", "1", reaction.ln_k, expected_ln_k,
               "T=350", "p3-reaction", "ln-k");
    CHECK(reaction.k.has_value());
    p3_compare("p3-reaction/dk", "dln_k_dT", "1/K", reaction.derivative, expected_derivative,
               "T=350", "p3-reaction", "dk");
    p3_compare("p3-reaction/h", "reaction_enthalpy", "J/mol", reaction.enthalpy,
               R * 350.0 * 350.0 * expected_derivative, "T=350", "p3-reaction", "h");

    bool endpoints = false, finite_k = true, outside_k = false;
    try {
        const auto lower = epcsaft::evaluate_reaction_properties(correlation, correlation.temperature_min);
        const auto upper = epcsaft::evaluate_reaction_properties(correlation, correlation.temperature_max);
        endpoints = true;
        finite_k = lower.k.has_value() && upper.k.has_value();
    } catch (const std::exception&) {}
    try {
        epcsaft::evaluate_reaction_properties(correlation, 280.0);
    } catch (const std::invalid_argument&) { outside_k = true; }
    recording_support::compare_custom("p3-reaction/endpoints", "K_endpoint_available", "bool",
        endpoints ? 1.0 : 0.0, 1.0, "Tmin/Tmax inclusive", "Availability", "correlation-interval",
        "correlation endpoints are owned and available", "absolute", endpoints, "p3-reaction", "endpoints");
    recording_support::compare_custom("p3-reaction/finite-k", "K_finite", "bool",
        finite_k ? 1.0 : 0.0, 1.0, "Tmin/Tmax", "Availability", "finite-K",
        "finite ln K produces finite K", "absolute", finite_k, "p3-reaction", "finite-k");
    const auto overflow = epcsaft::evaluate_reaction_properties(
        {1000.0, 0.0, 0.0, 0.0, 300.0, 293.15, 393.15, "P0"}, 300.0
    );
    const auto underflow = epcsaft::evaluate_reaction_properties(
        {-1000.0, 0.0, 0.0, 0.0, 300.0, 293.15, 393.15, "P0"}, 300.0
    );
    const bool k_separate = !overflow.k.has_value() && !underflow.k.has_value()
        && std::isfinite(overflow.ln_k) && std::isfinite(underflow.ln_k)
        && std::isfinite(overflow.derivative) && std::isfinite(underflow.enthalpy);
    recording_support::compare_custom("p3-reaction/k-unavailable", "K_availability_separate", "bool",
        k_separate ? 1.0 : 0.0, 1.0, "ln K=+/-1000", "Availability", "finite-K",
        "finite ln K, derivative and enthalpy survive K overflow/underflow",
        "absolute", k_separate, "p3-reaction", "k-unavailable");
    recording_support::compare_custom("p3-reaction/outside", "K_outside_rejected", "bool",
        outside_k ? 1.0 : 0.0, 1.0, "T=280 outside [293.15,393.15]", "Availability",
        "correlation-interval", "out-of-domain K is unavailable", "absolute", outside_k,
        "p3-reaction", "outside");


    const double tref = correlation.reference_temperature;
    const double species_h_ref = R * (-correlation.b + correlation.c * tref + correlation.d * tref * tref);
    const auto species_h = [=](double t) {
        return species_h_ref + R * (correlation.c * (t - tref) + correlation.d * (t * t - tref * tref));
    };
    p3_compare("p3-reaction/species-vanthoff", "species_delta_h", "J/mol",
               species_h(350.0), reaction.enthalpy, "supplied h_i(Tref),cp_i(T)",
               "p3-reaction", "species-vanthoff");
    constexpr int integration_steps = 1000;
    const double integration_step = (350.0 - tref) / integration_steps;
    double species_integral = species_h(tref) / (R * tref * tref) + species_h(350.0) / (R * 350.0 * 350.0);
    for (int i = 1; i < integration_steps; ++i) {
        const double t = tref + i * integration_step;
        species_integral += (i % 2 == 0 ? 2.0 : 4.0) * species_h(t) / (R * t * t);
    }
    species_integral *= integration_step / 3.0;
    p3_compare("p3-reaction/species-integral", "ln_k_species_integral", "1",
               species_integral, reaction.ln_k
                   - epcsaft::evaluate_reaction_properties(correlation, tref).ln_k,
               "supplied h_i(Tref),cp_i(T),Simpson-1000", "p3-reaction", "species-integral", 1.0e-8);
    bool missing_standard_state = false;
    try {
        auto missing = correlation;
        missing.standard_state_id.clear();
        epcsaft::evaluate_reaction_properties(missing, 350.0);
    } catch (const std::invalid_argument&) { missing_standard_state = true; }
    recording_support::compare_custom("p3-reaction/missing-standard-state", "standard_state_required", "bool",
        missing_standard_state ? 1.0 : 0.0, 1.0, "correlation identity absent and no gauge substitute",
        "Availability", "input-validation", "missing source standard state remains unavailable",
        "absolute", missing_standard_state, "p3-reaction", "missing-standard-state");
}


// Issue 26: every native permittivity form on the aqueous 2B-water Na+/Cl- liquid.
// A density-dependent eps_r enters pressure, the closure and every pressure-chart chain.
void modular_permittivity_actions() {
    auto ion_suppression = properties_support::aqueous_electrolyte();
    auto component_mixing = ion_suppression;
    component_mixing.permittivity_form = epcsaft::PermittivityForm::ComponentMixing;
    component_mixing.relative_permittivity = {78.09, 8.0, 8.0};
    auto density_temperature = ion_suppression;  // Yang2024 water d_V, d_T (full EOS 6.2)
    density_temperature.permittivity_form = epcsaft::PermittivityForm::DensityTemperature;
    density_temperature.dielectric_volume_m3_per_mol = {3.777e-4, 0.0, 0.0};
    density_temperature.dielectric_temperature_k = {1403.0, 0.0, 0.0};
    auto dipolar = ion_suppression;  // Rueben2024 PT2 records (full EOS 6.4)
    dipolar.permittivity_form = epcsaft::PermittivityForm::DipolarPerturbation;
    dipolar.dipole_scaling_debye2 = {5.199, 5.407, 2.563};
    dipolar.correlation_integral_parameter = {0.1276, 0.215, 0.138};
    const std::array<std::pair<const char*, epcsaft::Parameters>, 4> forms = {{
        {"ion-suppression", ion_suppression}, {"component-mixing", component_mixing},
        {"density-temperature", density_temperature}, {"dipolar", dipolar}}};
    const std::vector<epcsaft::ActiveParameter> active = {
        {epcsaft::ParameterFamily::Kij, {1, 3}},
        {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 2}}};
    const std::vector<epcsaft::PropertyObservable> observables = {
        epcsaft::PropertyObservable::ResidualEnthalpy,
        epcsaft::PropertyObservable::ResidualEntropy,
        epcsaft::PropertyObservable::ResidualIsobaricHeatCapacity};
    const auto parameter = [](double kij, double energy) {
        return std::vector<epcsaft::ParameterActionDirection>{
            {epcsaft::ParameterFamily::Kij, {1, 3}, kij},
            {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 2}, energy}};
    };
    const std::array<std::vector<epcsaft::ParameterActionDirection>, 3> theta = {
        parameter(0.01, 20.0), parameter(-0.005, -10.0), parameter(0.002, 5.0)};
    const std::array<double, 3> temperature_direction = {1.0, -0.5, 0.2};
    const std::array<std::vector<double>, 3> composition_direction = {{
        {-0.02, 0.01, 0.01}, {0.01, -0.005, -0.005}, {0.004, -0.002, -0.002}}};
    const std::array<double, 3> density_direction = {200.0, -100.0, 50.0};
    const std::array<double, 3> pressure_direction = {3.0e5, -2.0e5, 1.0e5};  // every stencil P > 0
    const std::array<double, 3> tolerances = {5.0e-6, 5.0e-5, 5.0e-4};
    const std::vector<double> x = {0.96, 0.02, 0.02};
    // 2 MPa: at least 1 MPa from the ion-bearing vapor-root fold of every form (scan below).
    constexpr double temperature = 298.15, pressure = 2.0e6;
    const std::array<double, 12> scan = {
        2e3, 1e4, 2e4, 5e4, 9e4, 9.5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6};
    const auto actions_match = [&](const char* name, const epcsaft::Parameters& parameters,
            const std::vector<epcsaft::ActiveParameter>& active,
            const std::array<std::vector<epcsaft::ParameterActionDirection>, 3>& theta, const std::vector<double>& x,
            const std::array<std::vector<double>, 3>& composition_direction, double pressure) {
        const epcsaft::HelmholtzModel model(parameters, active);
        const auto root = model.state_at_pressure(
            temperature, pressure, x, epcsaft::DensityBranch::Liquid);
        const double rho = root.molar_density;
        // Z - 1 = rho d(a_res)/d(rho) at fixed T, x: a dropped d(eps_r)/d(rho) fails it.
        const double h = 1.0e-5 * rho;
        const double slope = (model.state(temperature, rho + h, x).residual_helmholtz
            - model.state(temperature, rho - h, x).residual_helmholtz) / (2.0 * h);
        const double identity_reference = 1.0 + rho * slope;
        const bool identity_ok = std::abs(root.compressibility_factor - identity_reference)
            <= 1.0e-8 * std::abs(identity_reference - 1.0);
        recording_support::compare_custom(
            (std::string("issue26/pressure/") + name).c_str(), "compressibility_factor", "1",
            root.compressibility_factor, identity_reference,
            recording_support::fields({recording_support::field("form", name),
                recording_support::field("rho", rho),
                recording_support::field("epsilon_r", root.bulk_relative_permittivity),
                recording_support::field("relative_step", 1.0e-5)}),
            "Relationship", "central-difference", "pressure carries every density chain",
            "relative", identity_ok, "issue26-pressure-identity", name, "Z",
            "T,\\rho,\\mathbf{x}", "Z=1+\\rho\\,\\partial a^{res}/\\partial\\rho", "form");
        for (const bool pressure_chart : {false, true}) {
            for (std::size_t order = 1; order <= 3; ++order) {
                std::vector<epcsaft::DensityActionDirection> density_request;
                std::vector<epcsaft::PressureActionDirection> pressure_request;
                for (std::size_t k = 0; k < order; ++k) {
                    density_request.push_back({temperature_direction[k], density_direction[k],
                        composition_direction[k], theta[k]});
                    pressure_request.push_back({temperature_direction[k], pressure_direction[k],
                        composition_direction[k], theta[k]});
                }
                const auto batch = pressure_chart
                    ? model.pressure_actions(temperature, pressure, x,
                        epcsaft::DensityBranch::Liquid, observables, pressure_request)
                    : model.density_actions(temperature, rho, x, observables, density_request);
                const auto action_reference = pressure_chart
                    ? property_action_reference(parameters, active, std::nullopt, temperature, pressure, x,
                        observables, pressure_request, epcsaft::DensityBranch::Liquid)
                    : property_action_reference(parameters, active, std::nullopt, temperature, rho, x,
                        observables, density_request);
                for (std::size_t index = 0; index < observables.size(); ++index) {
                    const double fine = action_reference.value[index];
                    const double middle = action_reference.middle[index];
                    const double allowance = tolerances[order - 1] * action_reference.scale[index];
                    const auto& result = batch.results[index];
                    const bool ok = result.status == epcsaft::PropertyActionStatus::Available
                        && result.action.has_value()
                        && std::abs(*result.action - fine) <= allowance
                        && action_reference.uncertainty[index] <= allowance;
                    const std::string point = std::string(name) + (pressure_chart ? "/P/" : "/rho/")
                        + std::to_string(order) + "/" + std::to_string(index);
                    recording_support::compare_custom(("issue26/action/" + point).c_str(),
                        "action", "native",
                        result.action.value_or(std::numeric_limits<double>::quiet_NaN()), fine,
                        recording_support::fields({recording_support::field("form", name),
                            recording_support::field("chart", pressure_chart ? "T,P,x,theta" : "T,rho,x,theta"),
                            recording_support::field("order", order),
                            recording_support::field("observable", index),
                            recording_support::field("status", static_cast<int>(result.status)),
                            recording_support::field("reference_fine", fine),
                            recording_support::field("reference_middle", middle),
                            recording_support::field("tolerance", tolerances[order - 1]),
                            recording_support::field("term_scale", action_reference.scale[index]),
                            recording_support::field("reference_uncertainty", action_reference.uncertainty[index]),
                            recording_support::field("allowance", allowance)}),
                        "Relationship", "central-difference-refinement",
                        "liquid electrolyte action retains every permittivity chain",
                        "tolerance_times_chart_term_scale", ok, "issue26-actions", point.c_str(), "action",
                        "T,\\rho|P,\\mathbf{x},\\theta", "available", "form,chart,order,observable");
                }
            }
        }
        return root.bulk_relative_permittivity;
    };
    for (const auto& [name, parameters] : forms) {
        const epcsaft::HelmholtzModel model(parameters, active);
        double previous_liquid = 0.0;
        std::size_t scan_index = 0;
        // Independent native 4096/16384-cell diagnostics retained for #513, at the same declared inputs.
        const std::array<double, 12> liquid_reference = {56247.49942614573, 56247.597994176765, 56247.72120232484,
            56248.09081416554, 56248.583600548125, 56248.64519648327, 56248.70679189342, 56249.93858986418,
            56253.63272468506, 56259.7854243374, 56272.075143997776, 56308.81964538034};
        const std::array<double, 3> vapor_reference = {0.8503966976383335, 4.685632278617003, 10.471760503385712};
        for (const double scan_pressure : scan) {
            const auto liquid = model.state_at_pressure(
                temperature, scan_pressure, x, epcsaft::DensityBranch::Liquid);
            const auto vapor = model.state_at_pressure(
                temperature, scan_pressure, x, epcsaft::DensityBranch::Vapor);
            const int roots = liquid.stable_root_count.value_or(-1);
            // One stable root: both requests return it unlabelled. Two: liquid is the denser.
            // The liquid root rises monotonically in P, so no request switched branch.
            const bool ok = liquid.molar_density > previous_liquid
                && vapor.stable_root_count == liquid.stable_root_count
                && (roots == 1 ? close(vapor.molar_density, liquid.molar_density, 1e-10)
                        && liquid.density_branch == epcsaft::DensityBranch::Unspecified
                    : roots == 2 && vapor.molar_density < 1.0e-2 * liquid.molar_density
                        && liquid.density_branch == epcsaft::DensityBranch::Liquid)
                && std::abs(liquid.pressure - scan_pressure) <= epcsaft::policy::DENSITY_ROOT_RESIDUAL.value
                    * liquid.molar_density * epcsaft::GAS_CONSTANT_J_PER_MOL_K * temperature;
            if (parameters.permittivity_form.is(epcsaft::PermittivityForm::DensityTemperature)) {
                CHECK(close(liquid.molar_density, liquid_reference[scan_index], 1e-10));
                CHECK(close(vapor.molar_density, scan_index < 3 ? vapor_reference[scan_index] : liquid_reference[scan_index], 1e-10));
                CHECK(roots == (scan_index < 3 ? 2 : 1));
                CHECK(vapor.density_branch == (scan_index < 3 ? epcsaft::DensityBranch::Vapor : epcsaft::DensityBranch::Unspecified));
                CHECK(vapor.pressure_density_derivative > 0.0 && liquid.pressure_density_derivative > 0.0);
            }
            ++scan_index;
            previous_liquid = liquid.molar_density;
            const std::string point = std::string(name) + "/" + recording_support::number(scan_pressure);
            recording_support::compare_custom(("issue26/branch/" + point).c_str(),
                "liquid_molar_density", "mol/m3", liquid.molar_density, vapor.molar_density,
                recording_support::fields({recording_support::field("form", name),
                    recording_support::field("P", scan_pressure),
                    recording_support::field("stable_roots", roots),
                    recording_support::field("epsilon_r_liquid", liquid.bulk_relative_permittivity),
                    recording_support::field("epsilon_r_vapor", vapor.bulk_relative_permittivity),
                    recording_support::field("liquid_pressure_residual", liquid.pressure - scan_pressure)}),
                "Relationship", "declared-branch-scan", "liquid root persists; sole roots are unlabelled",
                "custom", ok, "issue26-branch", point.c_str(), "\\rho_L", "T,P,\\mathbf{x}",
                "branch", "form,P");
        }
        actions_match(name, parameters, active, theta, x, composition_direction, pressure);
    }
    // Issue 409: Bulow 2021 Fig. 9 [C4mim][TfO] + NaCl, every species charged; 1 bar liquid.
    auto il = properties_support::parameters(4);
    il.segment_count = {2.4805, 3.7432, 1.0, 1.0};
    il.sigma_constant = {3.6371, 3.8771, 2.8232, 2.756};
    il.dispersion_energy_over_k = {218.1441, 509.3113, 230.0, 170.0};
    il.molar_mass_kg_per_mol = {0.13922, 0.14907, 0.02298977, 0.035453};
    il.charge_number = {1, -1, 1, -1};
    il.packing_diameter_angstrom = {3.588496, 3.874333, 2.8232, 2.756};  // IL ions: PC-SAFT d(298.15 K)
    il.debye_huckel_diameter_angstrom = il.born_diameter_angstrom = il.sigma_constant;
    il.relative_permittivity = {13.2, 13.2, 8.0, 8.0};  // Bulow 2021 Table 5
    il.k_ij[11] = il.k_ij[14] = 0.317;  // Na+-Cl-, Bulow 2021 Table 4
    il.electrolyte_form = epcsaft::ElectrolyteForm::DebyeHuckelBorn;
    il.ionic_dispersion = epcsaft::IonicDispersion::AllIonPairs;
    il.permittivity_form = epcsaft::PermittivityForm::ComponentMixing;
    il.ionic_region_relative_permittivity = 1.0;
    const double salt = 0.001509;  // n_salt / (n_salt + n_IL)
    const auto kij = [](double value) {  // Na+-TfO-
        return std::vector<epcsaft::ParameterActionDirection>{{epcsaft::ParameterFamily::Kij, {6, 9}, value}};
    };
    // Normalized, electroneutral, and every salt fraction stays positive over the reference steps.
    const double epsilon = actions_match("all-ion-component-mixing", il, {{epcsaft::ParameterFamily::Kij, {6, 9}}},
        {kij(0.01), kij(-0.005), kij(0.002)}, {(1 - salt) / 2, (1 - salt) / 2, salt / 2, salt / 2},
        {{{-1e-3, -1e-3, 1e-3, 1e-3}, {5e-4, 5e-4, -5e-4, -5e-4}, {-2e-4, -2e-4, 2e-4, 2e-4}}}, 1.0e5);
    const double expected = 13.2 * (1 - salt) + 8.0 * salt;
    recording_support::compare_custom("issue409/epsilon", "bulk_relative_permittivity", "1", epsilon, expected,
        recording_support::fields({recording_support::field("x_salt", salt)}), "Reference", "source-equation",
        "Bulow 2019 Eq. 7 on ions alone", "relative", std::abs(epsilon - expected) <= 1e-12 * expected,
        "issue409-epsilon", "C4mimTfO-NaCl", "\\epsilon_r", "\\mathbf{x}", "13.2(1-x_s)+8x_s", "x_salt");
    // The forms that average over neutral solvents, or whose ion records are aqueous, still refuse.
    il.ion_suppression_coefficient = {7.01, 7.01, 7.01, 7.01};
    for (const auto form : {epcsaft::PermittivityForm::IonSuppression, epcsaft::PermittivityForm::DensityTemperature,
                            epcsaft::PermittivityForm::DipolarPerturbation}) {
        il.permittivity_form = form;
        bool refused = false;
        try { const epcsaft::HelmholtzModel model(il); }
        catch (const std::invalid_argument& error) { refused = std::string(error.what()) == "invalid dielectric solvent"; }
        CHECK(refused);
    }
}

// Issue 26: implicit-association E3/E4 energy actions against a centered difference of the
// next-lower exact order (u^T H v for E3, E3 for E4) along the last direction, Richardson h, h/2, h/4.
struct EnergyCase final {
    const char* name;
    epcsaft::Parameters parameters;
    std::vector<epcsaft::ActiveParameter> active;
    double temperature, molar_density;
    std::vector<double> composition;
    std::vector<std::vector<double>> directions;  // width 2 + C + A: T, V, n, theta
};

double shifted_energy(const EnergyCase& item, const std::vector<std::vector<double>>& lower,
                      const std::vector<double>& along, double t) {
    const std::size_t components = item.composition.size();
    std::vector<double> shift(item.active.size()), amounts(components);
    for (std::size_t a = 0; a < shift.size(); ++a) shift[a] = t * along[2 + components + a];
    for (std::size_t i = 0; i < components; ++i) amounts[i] = item.composition[i] + t * along[2 + i];
    const epcsaft::HelmholtzModel model(
        p4_shifted_parameters(item.parameters, item.active, shift), item.active);
    const double temperature = item.temperature + t * along[0];
    const double volume = 1.0 / item.molar_density + t * along[1];
    if (lower.size() == 3) return model.residual_energy_action(temperature, volume, amounts, lower).action.value();
    const auto derivatives = model.evaluate(temperature, volume, amounts, epcsaft::DerivativeOrder::Hessian);
    const std::size_t width = along.size();
    CHECK(derivatives.gradient.size() == width);
    double sum = 0.0;
    for (std::size_t i = 0; i < width; ++i)
        for (std::size_t j = 0; j < width; ++j)
            sum += lower[0][i] * derivatives.hessian[i * width + j] * lower[1][j];
    return sum;
}

void energy_action_ladders() {
    auto density_temperature = properties_support::aqueous_electrolyte();
    density_temperature.permittivity_form = epcsaft::PermittivityForm::DensityTemperature;
    density_temperature.dielectric_volume_m3_per_mol = {3.777e-4, 0.0, 0.0};
    density_temperature.dielectric_temperature_k = {1403.0, 0.0, 0.0};
    auto dipolar = properties_support::aqueous_electrolyte();
    dipolar.permittivity_form = epcsaft::PermittivityForm::DipolarPerturbation;
    dipolar.dipole_scaling_debye2 = {5.199, 5.407, 2.563};
    dipolar.correlation_integral_parameter = {0.1276, 0.215, 0.138};
    auto component_mixing = properties_support::aqueous_electrolyte();
    component_mixing.permittivity_form = epcsaft::PermittivityForm::ComponentMixing;
    component_mixing.relative_permittivity = {78.09, 8.0, 8.0};
    const std::vector<epcsaft::ActiveParameter> aqueous_active = {
        {epcsaft::ParameterFamily::Kij, {1, 3}},
        {epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 2}}};
    const double aqueous_volume = 1.0 / 56000.0;
    // Charge-neutral amount directions; V, T, k_ij and association-energy components nonzero.
    const std::vector<std::vector<double>> aqueous_directions = {
        {0.5, 1.0e-2 * aqueous_volume, -0.02, 0.01, 0.01, 0.01, 10.0},
        {-0.3, -0.5e-2 * aqueous_volume, 0.01, -0.005, -0.005, -0.004, -6.0},
        {0.2, 0.3e-2 * aqueous_volume, 0.004, -0.002, -0.002, 0.002, 4.0},
        {-0.1, 0.2e-2 * aqueous_volume, -0.006, 0.003, 0.003, -0.003, -3.0}};
    const std::vector<double> aqueous_x = {0.96, 0.02, 0.02};
    const double binary_volume = 1.0 / 42441.09427799721;
    const std::vector<EnergyCase> cases = {
        {"ion-suppression", properties_support::aqueous_electrolyte(), aqueous_active, 298.15, 56000.0, aqueous_x, aqueous_directions},
        {"component-mixing", component_mixing, aqueous_active, 298.15, 56000.0, aqueous_x, aqueous_directions},
        {"density-temperature", density_temperature, aqueous_active, 298.15, 56000.0, aqueous_x, aqueous_directions},
        {"dipolar", dipolar, aqueous_active, 298.15, 56000.0, aqueous_x, aqueous_directions},
        {"cross-associating-binary", properties_support::associating_binary(),
         {{epcsaft::ParameterFamily::AssociationEnergyOverK, {1, 4}},
          {epcsaft::ParameterFamily::AssociationVolume, {1, 4}}},
         313.15, 42441.09427799721, {0.25, 0.75}, {
            {0.5, 1.0e-2 * binary_volume, -0.01, 0.01, 10.0, 1.0e-3},
            {-0.3, -0.5e-2 * binary_volume, 0.005, -0.004, -6.0, -5.0e-4},
            {0.2, 0.3e-2 * binary_volume, 0.003, -0.002, 4.0, 3.0e-4},
            {-0.1, 0.2e-2 * binary_volume, -0.004, 0.003, -3.0, 2.0e-4}}},
    };
    for (const auto& item : cases) {
        const auto& d = item.directions;
        const std::array<std::pair<const char*, std::vector<std::vector<double>>>, 3> requests = {{
            {"E3", {d[0], d[1], d[2]}}, {"E3-repeated", {d[0], d[0], d[2]}}, {"E4", {d[0], d[1], d[2], d[3]}}}};
        const epcsaft::HelmholtzModel model(item.parameters, item.active);
        for (const auto& [label, request] : requests) {
            const auto action = model.residual_energy_action(
                item.temperature, 1.0 / item.molar_density, item.composition, request);
            const std::vector<std::vector<double>> lower(request.begin(), request.end() - 1);
            const auto difference = [&](double h) {
                // NOLINTNEXTLINE(clang-analyzer-core.CallAndMessage): clang-tidy 18 misreads a lambda's reference capture of a structured binding as uninitialized.
                return (shifted_energy(item, lower, request.back(), h)
                    - shifted_energy(item, lower, request.back(), -h)) / (2.0 * h);
            };
            // Two Richardson extrapolations (h, h/2) and (h/2, h/4): their difference is an
            // asymptotic, conservative estimate of the finer one's error (about 15x it for an h^4
            // remainder), not a bound. It must resolve 1e-7 relative, a decade below the 1e-6
            // agreement criterion that this falsifier uses to expose a missing chain.
            const double d1 = difference(1.0e-2), d2 = difference(5.0e-3), d3 = difference(2.5e-3);
            const double r1 = (4.0 * d2 - d1) / 3.0, reference = (4.0 * d3 - d2) / 3.0;
            const double scale = std::abs(reference);
            const double uncertainty = std::abs(reference - r1);
            const bool resolved = uncertainty <= 1.0e-7 * scale;
            const bool ok = action.status == epcsaft::PropertyActionStatus::Available
                && action.action.has_value() && resolved
                && std::abs(*action.action - reference) <= 1.0e-6 * scale;
            const std::string point = std::string(item.name) + "/" + label;
            recording_support::compare_custom(("issue26/energy/" + point).c_str(), "energy_action",
                "extensive residual scalar", action.action.value_or(std::numeric_limits<double>::quiet_NaN()),
                reference, recording_support::fields({recording_support::field("case", item.name),
                    recording_support::field("action", label),
                    recording_support::field("status", static_cast<int>(action.status)),
                    recording_support::field("difference_h", d1),
                    recording_support::field("difference_h_over_2", d2),
                    recording_support::field("difference_h_over_4", d3),
                    recording_support::field("reference_uncertainty", uncertainty),
                    recording_support::field("association_residual", action.diagnostics.association_residual)}),
                "Relationship", "centered-difference-richardson",
                "implicit association and permittivity chains complete through fourth order",
                "relative", ok, "issue26-energy", point.c_str(), label, "T,V,\\mathbf{n},\\theta",
                "D^kA=\\partial_t D^{k-1}A", "case,action");
        }
    }
}

// C4 (issue #224 B): exact E5 F5[w,w,w,a,b] against the centered difference of exact E4 F4[w,w,w,a] along b,
// w = (0, 0, sqrt(n0), -sqrt(n1)); directions scaled by T, V and each amount, plus one mixed pair.
void energy_e5_ladder() {
    struct Case { const char* name; epcsaft::Parameters parameters; double temperature; std::vector<double> n; };
    for (const auto& item : std::vector<Case>{{"nonassociating", properties_support::neutral_binary(), 260.0, {0.27, 0.73}},
                                              {"associating", properties_support::associating_binary(), 376.88, {0.1385, 0.8615}}}) {
        const epcsaft::HelmholtzModel model(item.parameters);
        const double T = item.temperature, V = 1.0e-4;
        const std::vector<double> w = {0, 0, std::sqrt(item.n[0]), -std::sqrt(item.n[1])};
        std::vector<std::vector<double>> axes = {{T, 0, 0, 0}, {0, V, 0, 0}, {0, 0, item.n[0], 0}, {0, 0, 0, item.n[1]}};
        std::vector<std::pair<std::vector<double>, std::vector<double>>> pairs;
        for (std::size_t a = 0; a < 4; ++a) for (std::size_t b = a; b < 4; ++b) pairs.push_back({axes[a], axes[b]});
        pairs.push_back({{0.3 * T, 0.2 * V, 0.5 * item.n[0], -0.4 * item.n[1]}, {-0.2 * T, 0.4 * V, -0.3 * item.n[0], 0.6 * item.n[1]}});
        const auto E = [&](const std::vector<std::vector<double>>& d, double t, const std::vector<double>& b) {
            std::vector<double> n = {item.n[0] + t * b[2], item.n[1] + t * b[3]};
            return model.residual_energy_action(T + t * b[0], V + t * b[1], n, d);
        };
        for (std::size_t k = 0; k < pairs.size(); ++k) {
            const auto& [a, b] = pairs[k];
            const auto exact = E({w, w, w, a, b}, 0.0, b);
            CHECK(exact.action.has_value());
            double magnitude = 0.0;  // fixed chart sub-direction contractions, the scale near a zero action
            for (std::size_t i = 0; i < 4; ++i) {
                auto unit = std::vector<double>(4, 0.0);
                unit[i] = b[i];
                if (b[i] != 0.0) magnitude += std::abs(E({w, w, w, a, unit}, 0.0, b).action.value());
            }
            const double scale = std::max(std::abs(*exact.action), magnitude);
            std::vector<double> errors;
            for (const double h : {1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5})
                errors.push_back(std::abs(*exact.action - (E({w, w, w, a}, h, b).action.value() - E({w, w, w, a}, -h, b).action.value()) / (2 * h)) / scale);
            const double best = *std::min_element(errors.begin(), errors.end());
            const std::string id = std::string("C4/") + item.name + "/" + std::to_string(k);
            recording_support::compare_custom(id.c_str(), "e5_energy_action", "extensive residual scalar", *exact.action, scale,
                recording_support::fields({recording_support::field("errors_over_scale", errors), recording_support::field("error", best),
                    recording_support::field("criterion", 1e-5), recording_support::field("site_residual", exact.diagnostics.association_residual)}),
                "Relationship", "issue-224", "E5 equals the centered difference of exact E4", "action or sub-direction magnitude",
                best <= 1e-5, "critical-envelope", id.c_str());
        }
    }
}

// Copied anchor compositions from 2026-born-unavailable-states/inputs, reordered
// CO2/MEA/water -> properties_support's water/MEA/CO2. EqIDs own the radial work.
void born_parameter_columns() {
    const std::array<std::vector<double>, 2> anchors = {{{0.8863734508600574, 0.07827110880066358, 3.8480681348098447e-08, 0.01818111244399814, 0.015766643324345936, 0.0003959814177106597, 0.0010068230563202352, 1.3460805650294977e-11, 4.841602761876763e-06},
        {0.8448866085505489, 0.004900981383765356, 0.0002933071241251476, 0.07818889597821514, 0.029096076931579815, 0.03617530531138919, 0.006458693289371563, 2.1372510655078036e-09, 1.2929375400277957e-07}}};
    for (std::size_t case_index = 0; case_index < anchors.size(); ++case_index)
    for (const double c : {0.0, 1.0}) for (const double epsilon : {2.0, 8.0}) {
        auto p = properties_support::nine_species();
        p.born_c_shell = 1.0; p.born_c_dielectric = c;
        p.ionic_region_relative_permittivity = epsilon;
        const auto& x = anchors[case_index];
        const epcsaft::HelmholtzModel model(p, {{epcsaft::ParameterFamily::BornCDielectric, {0}},
            {epcsaft::ParameterFamily::IonicRegionRelativePermittivity, {0}}});
        const double t = 313.15, rho = 45000.0;
        const auto state = model.state(t, rho, x);
        const auto derivatives = model.evaluate(t, 1.0 / rho, x, epcsaft::DerivativeOrder::Gradient);
        double solvation = 0.0, total = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) { solvation += x[i] * p.solvation_factor[i]; total += x[i]; }
        solvation /= total;
        double shell_sum = 0.0, correction = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) if (p.charge_number[i]) {
            const double z = p.charge_number[i], d = p.born_diameter_angstrom[i] * epcsaft::ANGSTROM_M;
            const double outer = d * (1.0 + (solvation - 1.0) / std::abs(z));
            shell_sum += x[i] * z * z / outer;
            correction += x[i] * z * z * (1.0 / d - 1.0 / outer);
        }
        const double prefactor = -epcsaft::ELEMENTARY_CHARGE_C * epcsaft::ELEMENTARY_CHARGE_C
            / (4.0 * epcsaft::PI * epcsaft::VACUUM_PERMITTIVITY_F_PER_M * epcsaft::BOLTZMANN_J_PER_K * t);
        const std::array<double, 3> reference = {
            prefactor * ((1.0 - 1.0 / state.bulk_relative_permittivity) * shell_sum + c * (1.0 - 1.0 / epsilon) * correction),
            prefactor * (1.0 - 1.0 / epsilon) * correction,
            prefactor * c * correction / (epsilon * epsilon)};
        const std::array<double, 3> actual = {state.born, derivatives.gradient[2 + x.size()], derivatives.gradient[3 + x.size()]};
        for (std::size_t column = 0; column < actual.size(); ++column) {
            const auto id = std::to_string(case_index + 1) + "/" + std::to_string(c) + "/" + std::to_string(epsilon) + "/" + std::to_string(column);
            const double atol = reference[column] == 0.0 ? 1e-12 : 0.0;
            const bool ok = std::abs(actual[column] - reference[column]) <= atol + 1e-12 * std::abs(reference[column]);
            recording_support::compare(("born-parameter/" + id).c_str(), "Born radial work and first columns", "dimensionless",
                actual[column], reference[column], recording_support::fields({recording_support::field("epsilon", epsilon), recording_support::field("c", c)}),
                "Relationship", "born_radial_charging_work", "frozen-state Born value, c coefficient and dq/depsilon chain", "reference_scale",
                atol, 1e-12, ok, "born-parameter", id.c_str(), "a_Born", "c,epsilon", "radial work", "case,c,epsilon,column");
        }
    }
}

using PO = epcsaft::PropertyObservable;
using PD = epcsaft::PressureActionDirection;
using PF = epcsaft::ParameterFamily;
constexpr double reference_nan = std::numeric_limits<double>::quiet_NaN();
// Source-pressure liquid root at a fixed terminal composition, then the runtime reference partial.
epcsaft::ResidualEnergyAction anchored_reference_partial(const epcsaft::HelmholtzModel& model, double T, double P,
    const std::vector<double>& x, const std::vector<double>& nu, std::size_t a, std::size_t b,
    std::optional<std::size_t> parameter, std::optional<double> anchor = std::nullopt) {
    return epcsaft::reference_partial_from_state(model, T, x, nu, a, b, parameter,
        model.state_at_pressure(T, P, x, epcsaft::DensityBranch::Liquid, anchor),
        [&] { for (std::size_t i=0;i<x.size();++i) if (x[i]>0 && model.parameters().charge_number[i]) return false; return true; }());
}
double factorial(unsigned n) { double v = 1; for (unsigned k=2; k<=n; ++k) v *= k; return v; }
double binomial(unsigned n, unsigned k) { return factorial(n)/(factorial(k)*factorial(n-k)); }

// Independent P4 generator: C_T=(h_res+RT)/(RT^2).
double pure_reference_thermal(const epcsaft::HelmholtzModel& model, double T, double P,
    unsigned a, unsigned b, bool theta, double rho) {
    const double ps=model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid,rho).pressure_density_derivative
        *model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid,rho).molar_density;
    const auto h = [&](unsigned thermal) {
        std::vector<PD> d;
        for (unsigned j=0; j<thermal; ++j) { PD q; q.composition={0.0}; q.temperature=T; d.push_back(q); }
        for (unsigned j=0; j<b; ++j) { PD q; q.composition={0.0}; q.pressure=ps; d.push_back(q); }
        if (theta) { PD q; q.composition={0.0}; q.parameter={{PF::DispersionEnergyOverK,{0},191.42}}; d.push_back(q); }
        if (d.empty()) return model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid,rho).residual_enthalpy;
        // Heat capacity first/second actions independently supply h thermal derivatives.
        const auto observable = thermal ? PO::ResidualIsobaricHeatCapacity : PO::ResidualEnthalpy;
        if (thermal) d.erase(d.begin());
        if (d.empty()) return *model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid,rho).residual_isobaric_heat_capacity;
        return model.pressure_actions(T,P,{1},epcsaft::DensityBranch::Liquid,{observable},d,rho)
            .results.front().action.value_or(reference_nan) / (std::pow(T,int(thermal?thermal-1:0))*std::pow(ps,int(b))*(theta?191.42:1));
    };
    double value=0;
    for (unsigned j=0; j<a; ++j)
        value += binomial(a-1,j) * (j%2 ? -1:1) * factorial(j+1)
            * h(a-1-j)/std::pow(T,int(j+2))/epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    if (!b && !theta) value += ((a-1)%2 ? -1:1)*factorial(a-1)/std::pow(T,int(a));
    return value;
}

// Independent pressure-volume inverse through E4; no reference Taylor path.
double pure_reference_pressure(const epcsaft::HelmholtzModel& model, double T, double P,
    unsigned b, bool theta, double rho) {
    rho=model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid,rho).molar_density;
    const double V=1/rho, rt=epcsaft::GAS_CONSTANT_J_PER_MOL_K*T;
    const auto e=model.evaluate(T,V,{1},epcsaft::DerivativeOrder::Hessian);
    const auto w=e.gradient.size();
    std::vector<double> v(w,0), k(w,0); v[1]=1; k[w-1]=1;
    const auto energy=[&](std::vector<std::vector<double>> d) {
        return model.residual_energy_action(T,V,{1},d).action.value_or(reference_nan);
    };
    const double pv=rt*(-1/(V*V)-e.hessian[w+1]);
    const double pvv=rt*(2/(V*V*V)-energy({v,v,v}));
    const double pvvv=rt*(-6/std::pow(V,4)-energy({v,v,v,v}));
    double value=V;
    if (theta) {
        const double vk=rt*e.hessian[w+w-1]/pv;
        const double pvk=-rt*energy({v,v,k})+pvv*vk;
        const double pvvk=-rt*energy({v,v,v,k})+pvvv*vk;
        if (b==1) value=vk;
        if (b==2) value=-pvk/(pv*pv);
        if (b==3) value=-pvvk/std::pow(pv,3)+3*pvv*pvk/std::pow(pv,4);
    } else {
        if (b==2) value=1/pv;
        if (b==3) value=-pvv/std::pow(pv,3);
        if (b==4) value=3*pvv*pvv/std::pow(pv,5)-pvvv/std::pow(pv,4);
    }
    return -value/rt;
}

void reference_f8() {
    epcsaft::HelmholtzModel model(properties_support::ethane(),{{PF::DispersionEnergyOverK,{0}}});
    const double T=260, R=epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    for (double P : {3e6,5e6}) {
        const auto state=model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid);
        for (unsigned c=0; c<=1; ++c) for (unsigned degree=1; degree<=4; ++degree)
            for (unsigned b=0; b<=degree-c; ++b) {
                const unsigned a=degree-b-c;
                const double scale=R*std::pow(T,int(a))*std::pow(P,int(b))*(c?191.42:1);
                const auto exact=anchored_reference_partial(model,T,P,{1},{-1},a,b,c?std::optional<std::size_t>(0):std::nullopt,state.molar_density);
                const double oracle=(c && !a && !b) ? -model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid,state.molar_density).log_fugacity_parameter_derivative[0].value() : a ? pure_reference_thermal(model,T,P,a,b,c,state.molar_density)
                    : pure_reference_pressure(model,T,P,b,c,state.molar_density);
                const double x=scale*exact.action.value_or(reference_nan), y=scale*oracle;
                const double target=1e-7+1e-9*std::abs(y);
                const bool ok=exact.status==epcsaft::PropertyActionStatus::Available && std::abs(x-y)<=target;
                const std::string id="F8/"+std::to_string(int(P))+"/"+std::to_string(a)+std::to_string(b)+std::to_string(c);
                recording_support::compare_custom(id.c_str(),"scaled reference partial","J/mol/K",x,y,
                    recording_support::fields({recording_support::field("T",T),recording_support::field("P_source",P),
                    recording_support::field("a",double(a)),recording_support::field("b",double(b)),recording_support::field("theta",double(c)),
                    recording_support::field("target",target),recording_support::field("error",std::abs(x-y)),recording_support::field("error_over_target",std::abs(x-y)/target),
                    recording_support::field("site_residual",exact.diagnostics.association_residual),
                    recording_support::field("pressure_series_residual",exact.diagnostics.pressure_residual)}),
                    "Relationship","pure-P4-identities","finite source-pressure reference derivatives","absolute+relative",ok,"F8");
            }
    }
    auto parameters=properties_support::parameters(2);
    parameters.segment_count={1.6069,1.6069};parameters.sigma_constant={3.5206,3.5206};
    parameters.dispersion_energy_over_k={191.42,191.42};parameters.molar_mass_kg_per_mol={0.03007,0.03007};
    epcsaft::HelmholtzModel pair(parameters,{{PF::DispersionEnergyOverK,{0}},{PF::DispersionEnergyOverK,{1}}});
    for (double P:{3e6,5e6}) for(double molality:{1e-10,1e-12,1e-14}) {
        const double solute=molality*0.03007;
        const std::vector<double> composition{1/(1+solute),solute/(1+solute)};
        const auto state=pair.state_at_pressure(T,P,composition,epcsaft::DensityBranch::Liquid);
        for (const auto& column:std::vector<std::array<unsigned,3>>{{1,0,0},{2,0,0},{0,0,1},{1,0,1},{1,1,0}}) {
            const auto a=column[0],b=column[1],c=column[2];
            auto exact=anchored_reference_partial(pair,T,P,composition,{-2,1},a,b,c?std::optional<std::size_t>(0):std::nullopt,state.molar_density);
            if (c) {
                const auto other=anchored_reference_partial(pair,T,P,composition,{-2,1},a,b,1,state.molar_density);
                if(exact.action && other.action) *exact.action+=*other.action;else exact.action.reset();
            }
            const auto pure=model.state_at_pressure(T,P,{1},epcsaft::DensityBranch::Liquid);
            const double oracle=a?pure_reference_thermal(model,T,P,a,b,c,pure.molar_density)
                :-pure.log_fugacity_parameter_derivative[0].value();
            const double scale=R*std::pow(T,int(a))*std::pow(P,int(b))*(c?191.42:1);
            const double value=scale*exact.action.value_or(reference_nan), expected=scale*oracle;
            const double target=1e-7+1e-9*std::abs(expected);
            const bool ok=exact.action && std::abs(value-expected)<=target;
            const std::string id="F8/lower-terminal/"+std::to_string(int(P))+"/"+std::to_string(int(-std::log10(molality)))
                +"/"+std::to_string(a)+std::to_string(b)+std::to_string(c);
            recording_support::compare_custom(id.c_str(),"scaled reference partial","J/mol/K",value,expected,
                recording_support::fields({recording_support::field("T",T),recording_support::field("P_source",P),
                    recording_support::field("molality",molality),recording_support::field("a",a),recording_support::field("b",b),
                    recording_support::field("theta",c),recording_support::field("target",target),
                    recording_support::field("error_over_target",std::abs(value-expected)/target),
                    recording_support::field("density",state.molar_density),recording_support::field("pressure_root_residual",state.pressure-P),
                    recording_support::field("site_residual",exact.diagnostics.association_residual)}),
                "Relationship","pure-P4-identities","lower slopes at fixed E terminal compositions; tied EOS direction",
                "absolute+relative",ok,"F8-lower");
        }
    }

}

// E2/E3 Hessian-chain oracle at a literal fixed composition. No limiting selection.
double reference_lower(const epcsaft::HelmholtzModel& m, double T,
    const std::vector<double>& x, unsigned a, unsigned b,
    const epcsaft::State& s, const epcsaft::HelmholtzDerivatives& e,
    const std::vector<double>& reaction = {-1,1,1}) {
    const double V=1/s.molar_density, R=epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const auto w=e.gradient.size();
    std::vector<double> nu(w,0), ev(w,0); ev[1]=1;
    double delta=0, value=0;
    for (std::size_t i=0; i<reaction.size(); ++i) {
        nu[2+i]=reaction[i];delta+=reaction[i];
        value+=reaction[i]*(std::log(s.molar_density)+s.residual_chemical_potential_over_rt[i].value());
    }
    if (a+b==0) return value;
    const double pv=R*T*(-1/(V*V)-e.hessian[w+1]);
    const double pt=R*(1/V-e.gradient[1])-R*T*e.hessian[w];
    double cv=-delta/V;
    for (std::size_t i=2; i<2+reaction.size(); ++i) cv+=nu[i]*e.hessian[i*w+1];
    std::vector<std::vector<double>> path;
    for (unsigned j=0; j<a+b; ++j) {
        std::vector<double> d(w,0);
        if (j<a) { d[0]=1;d[1]=-pt/pv; } else d[1]=1/pv;
        path.push_back(d);
    }
    if (a+b==1) {
        double val=cv*path[0][1];
        for (std::size_t i=2; i<2+reaction.size(); ++i) val+=nu[i]*e.hessian[i*w]*path[0][0];
        return val;
    }
    const auto& u=path[0];const auto& v=path[1];
    const auto mixed=[&](const std::vector<double>& q) {
        // Size both noncomposition motions before half-sign polarization.
        const auto norm=[&](const auto& d) { return std::max(std::abs(d[0])/T,std::abs(d[1])/V); };
        const double alpha=norm(v)/norm(u);
        auto plus=v,minus=v;
        for(std::size_t i=0;i<w;++i){plus[i]+=alpha*u[i];minus[i]-=alpha*u[i];}
        const auto p=m.residual_energy_action(T,V,x,{q,plus,plus});
        const auto n=m.residual_energy_action(T,V,x,{q,minus,minus});
        return (p.action.value_or(reference_nan)-n.action.value_or(reference_nan))/(4*alpha);
    };
    const auto hv=[&](const auto& d){double v=0;for(std::size_t i=0;i<w;++i)v+=e.hessian[w+i]*d[i];return v;};
    const double p2=R*(u[0]*(-v[1]/(V*V)-hv(v))+v[0]*(-u[1]/(V*V)-hv(u)))
        +R*T*(2*u[1]*v[1]/std::pow(V,3)-mixed(ev));
    return mixed(nu)+delta*u[1]*v[1]/(V*V)-cv*p2/pv;
}

void reference_f9() {
    const double T=298.15,R=epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    auto params=properties_support::aqueous_electrolyte();
    for (const auto family : {PF::DispersionEnergyOverK, PF::BornDiameter}) {
    const std::size_t index = family == PF::BornDiameter ? 1 : 0;
    const double theta = index ? 3.445 : 353.95;
    epcsaft::HelmholtzModel model(params,{{family,{index}}});
    // Each parameter offset owns one model; T/P shifts use the unchanged model.
    std::map<double, std::unique_ptr<epcsaft::HelmholtzModel>> shifted_models;
    const auto shifted_model = [&](double step) -> const epcsaft::HelmholtzModel& {
        auto& cached = shifted_models[step];
        if (!cached) {
            auto p = params;
            if (index) p.born_diameter_angstrom[index] += step;
            else p.dispersion_energy_over_k[index] += step;
            cached = std::make_unique<epcsaft::HelmholtzModel>(p,
                std::vector<epcsaft::ActiveParameter>{{family,{index}}});
        }
        return *cached;
    };
    for(double P:{1e5,2e5}) for(double ratio:{1.0,2.0}) {
        std::optional<double> anchor;
        for(double molality:{0.0,1e-10,1e-12,1e-14}) {
            const double salt=molality*params.molar_mass_kg_per_mol[0]*ratio;
            const std::vector<double> x{1/(1+2*salt),salt/(1+2*salt),salt/(1+2*salt)};
            const auto state=model.state_at_pressure(T,P,x,epcsaft::DensityBranch::Liquid,anchor);
            anchor=state.molar_density;
            // The lower-column oracle uses the same state and Hessian for every
            // column at a given shifted model, temperature and pressure.
            using LowerPoint = std::tuple<double, double, double>;
            std::map<LowerPoint, std::pair<epcsaft::State, epcsaft::HelmholtzDerivatives>> lower_states;
            for(unsigned c=0;c<=1;++c) for(unsigned d=1;d<=4;++d) for(unsigned b=0;b<=d-c;++b) {
                const unsigned a=d-b-c;
                if(d==1 && b)continue;
                if(d==2 && !(a==2 || (a==1 && b==1) || (a==1 && c==1) || (c && b)))continue;
                const auto exact=anchored_reference_partial(model,T,P,x,{-1,1,1},a,b,
                    c?std::optional<std::size_t>(0):std::nullopt,*anchor);
                const double scale=R*std::pow(T,int(a))*std::pow(P,int(b))*(c?theta:1);
                const double value=scale*exact.action.value_or(reference_nan);
                const auto lower=[&](double step) {
                    double t=T,pressure=P;
                    unsigned la=a,lb=b;
                    if(!c && a){t+=step;--la;}else if(!c){pressure+=step;--lb;}
                    const auto& perturbed = c ? shifted_model(step) : model;
                    if(molality != 0.0 && la+lb<=2) {
                        const LowerPoint key{c ? step : 0.0, t, pressure};
                        auto found = lower_states.find(key);
                        if (found == lower_states.end()) {
                            const auto state = perturbed.state_at_pressure(t, pressure, x,
                                epcsaft::DensityBranch::Liquid, *anchor);
                            auto energy = perturbed.evaluate(t, 1/state.molar_density, x,
                                epcsaft::DerivativeOrder::Hessian);
                            found = lower_states.emplace(key, std::make_pair(state, std::move(energy))).first;
                        }
                        return reference_lower(perturbed,t,x,la,lb,
                            found->second.first, found->second.second);
                    }
                    return anchored_reference_partial(perturbed,t,pressure,x,{-1,1,1},la,lb,std::nullopt,*anchor).action.value_or(reference_nan);
                };
                double h=c?(index?0.002:0.02):a?0.2:2e4;
                bool accepted=false;
                double coarse=scale*(lower(h)-lower(-h))/(2*h);
                for(unsigned halving=0;halving<=6;++halving) {
                    const double plus=lower(h/2),minus=lower(-h/2);
                    const double fine=scale*(plus-minus)/h;
                    const double rich=(4*fine-coarse)/3,refine=std::abs(fine-coarse)/3;
                    const double sigma=std::abs(value); // One reference-chart basis monomial.
                    const double target=1e-7+1e-6*std::max(std::abs(rich),sigma);
                    const bool resolved=refine<=target/4;
                    accepted=resolved&&std::abs(value-rich)<=target;
                    const std::string id="F9/partial/"+std::to_string(int(P))+"/"+std::to_string(int(ratio))+"/"+std::to_string(molality==0?0:int(-std::log10(molality)))+"/"+std::to_string(index)
                        +"/"+std::to_string(a)+std::to_string(b)+std::to_string(c);
                    const auto inputs=recording_support::fields({recording_support::field("T",T),recording_support::field("P_source",P),
                        recording_support::field("molality",molality),recording_support::field("ratio",ratio),
                        recording_support::field("a",double(a)),recording_support::field("b",double(b)),recording_support::field("theta",double(c)),
                        recording_support::field("h",h),recording_support::field("coarse",coarse),recording_support::field("fine",fine),
                        recording_support::field("refinement",refine),recording_support::field("target",target),recording_support::field("sigma",sigma),
                        recording_support::field("error",std::abs(value-rich)),recording_support::field("error_over_target",std::abs(value-rich)/target),recording_support::field("resolved",resolved?1:0),
                        recording_support::field("site_residual",exact.diagnostics.association_residual),
                        recording_support::field("density",state.molar_density),recording_support::field("pressure_root_residual",state.pressure-P),
                        recording_support::field("pressure_density_slope",state.pressure_density_derivative),
                        recording_support::field("site_pivot_ratio",exact.diagnostics.site_pivot_ratio),
                        recording_support::field("rounding_over_step",64*std::numeric_limits<double>::epsilon()*scale*(std::abs(plus)+std::abs(minus))/h),
                        recording_support::field("pressure_series_residual",exact.diagnostics.pressure_residual)});
                    if (!resolved && halving < 6)
                        recording_support::observe(("F9/refinement/"+id.substr(11)+"/"+std::to_string(halving)).c_str(),"scaled reference partial","J/mol/K",
                            value,rich,inputs,"Numerical","fixed-composition-Richardson","coarse step requires refinement");
                    else recording_support::compare_custom(id.c_str(),"scaled reference partial","J/mol/K",value,rich,inputs,
                        "Relationship","fixed-composition-Richardson","finite terminal ionic derivatives","absolute+relative",accepted,"F9");
                    if(resolved)break;
                    h/=2;coarse=fine;
                }
                CHECK(accepted);
            }
        }
    }
}
}
void reference_limits() {
    const double T=298.15,P=1e5,R=epcsaft::GAS_CONSTANT_J_PER_MOL_K;
    const auto params=properties_support::aqueous_electrolyte();
    epcsaft::HelmholtzModel model(params,{{PF::DispersionEnergyOverK,{0}}});
    for(bool fixed:{false,true}) {
        const double src=fixed?2e5:P;
        epcsaft::NeutralReferenceDeclaration declaration{0,{0,1,1},1};
        epcsaft::ReactionReference reference{fixed?std::optional<double>(src):std::nullopt,
            epcsaft::ReferenceSourceBasis::CommonMolalityInfiniteDilution};
        const auto jet=epcsaft::reaction_references(model,T,P,declaration,{{-1,1,1}},{reference},true,std::nullopt,2,true,4)[0];
        CHECK(jet.status==epcsaft::NeutralReferenceStatus::Available);
        const std::vector<std::array<unsigned,3>> lower_columns{{1,0,0},{2,0,0},{0,0,1},{1,0,1},{1,1,0}};
        for(const auto& column:lower_columns) {
            const auto a=column[0],b=column[1],c=column[2];
            std::array<double,3> values{};
            unsigned rung=0;
            for(double molality:{1e-10,1e-12,1e-14}) {
                const double salt=molality*params.molar_mass_kg_per_mol[0];
                const std::vector<double> terminal{1/(1+2*salt),salt/(1+2*salt),salt/(1+2*salt)};
                values[rung++]=anchored_reference_partial(model,T,src,terminal,{-1,1,1},a,b,
                    c?std::optional<std::size_t>(0):std::nullopt).action.value_or(reference_nan);
            }
            const double scale=std::pow(T,int(a))*std::pow(src,int(b))*(c?353.95:1)*(a==2?R:1);
            const double tolerance=a==2?0.01:5e-5;
            const auto admitted=a==2?jet.terminal_contraction_temperature_second_derivative
                :c?(a?jet.terminal_contraction_temperature_parameter_derivative[0]:jet.terminal_contraction_parameter_derivative[0])
                :b?jet.terminal_contraction_temperature_pressure_derivative:jet.terminal_contraction_temperature_derivative;
            CHECK(admitted);
            for(unsigned k=0;k<2;++k) {
                const double change=scale*std::abs(values[k+1]-values[k]);
                const std::string id="F9/lower-limit/"+std::to_string(fixed)+"/"+std::to_string(a)+std::to_string(b)+std::to_string(c)+"/"+std::to_string(k);
                const bool ok=change<=tolerance;
                recording_support::compare_custom(id.c_str(),"scaled lower reference change",a==2?"J/mol/K":"dimensionless",change,0,
                    recording_support::fields({recording_support::field("target",tolerance),recording_support::field("source_pressure",src),
                        recording_support::field("lower_molality",k?1e-12:1e-10),recording_support::field("upper_molality",k?1e-14:1e-12)}),
                    "Relationship","terminal-ladder","unchanged lower-column admission scales","absolute",ok,"F9-lower-limit");
            }
            const double limit=(10*values[2]-values[1])/9;  // sqrt(m) Richardson limit of the terminal pair (#353)
            CHECK(scale*std::abs(*admitted-(fixed && b?0:limit))<=tolerance);
        }
        auto doubled=declaration;doubled.solute_molality_ratios={0,2,2};
        const auto jet2=epcsaft::reaction_references(model,T,P,doubled,{{-1,1,1}},{reference},true,std::nullopt,2,true,4)[0];
        CHECK(jet2.status==epcsaft::NeutralReferenceStatus::Available);
        for(std::size_t i=0;i<jet.derivatives.size();++i) {
            const auto& column=jet.derivatives[i];const auto& second=jet2.derivatives[i];
            CHECK(column.value && second.value);
            const double scale=R*std::pow(T,int(column.temperature_order))*std::pow(src,int(column.pressure_order))
                *(column.parameter?353.95:1);
            const double change=scale*std::abs(*column.value-*second.value);
            const std::string id="F9/doubled/"+std::to_string(fixed)+"/"+std::to_string(i);
            recording_support::compare_custom(id.c_str(),"doubled terminal concentration change","J/mol/K",change,0,
                recording_support::fields({recording_support::field("target",0.01)}),"Relationship","same-neutral-curve",
                "doubled terminal concentration, no different ionic-ratio claim","absolute",change<=0.01,"F9-limit");
            if(fixed && column.pressure_order) {
                const std::string zero="F8/system-P-zero/"+std::to_string(i);
                recording_support::compare_custom(zero.c_str(),"fixed-source system-pressure derivative","raw SI",*column.value,0,
                    "{}","Relationship","analytic-zero","fixed-source system-P terms vanish exactly","absolute",*column.value==0,"F8");
            }
        }

        for(const auto& column:jet.derivatives) {
            CHECK(column.value);
            if(fixed && column.pressure_order) continue;
            std::array<double,3> values{};
            unsigned rung=0;
            for(double molality:{1e-10,1e-12,1e-14}) {
                const double salt=molality*params.molar_mass_kg_per_mol[0];
                const std::vector<double> x{1/(1+2*salt),salt/(1+2*salt),salt/(1+2*salt)};
                values[rung++]=anchored_reference_partial(model,T,src,x,{-1,1,1},column.temperature_order,
                    column.pressure_order,column.parameter).action.value_or(reference_nan);
            }
            const double limit=(10*values[2]-values[1])/9;
            const double scale=R*std::pow(T,int(column.temperature_order))*std::pow(src,int(column.pressure_order))
                *(column.parameter?353.95:1);
            const double gap=scale*(*column.value-limit);
            const std::string id="F9/limit/"+std::to_string(fixed)+"/"+std::to_string(column.temperature_order)
                +std::to_string(column.pressure_order)+std::to_string(column.parameter.has_value());
            recording_support::compare_custom(id.c_str(),"scaled analytic-minus-ladder reference","J/mol/K",gap,0,
                recording_support::fields({recording_support::field("target",0.01)}),"Relationship","terminal-ladder",
                "analytic limit compared with retained finite calculation","absolute",std::abs(gap)<=0.01,"F9-limit");
        }
        // Existing equilibrium driver supplies only the case-center composition.
        // This independent verification completion never supplies runtime calorics.
        auto problem=equilibrium_support::flash(3);
        problem.temperature=epcsaft::equilibrium::FreeValue::fixed(T);
        problem.pressure=epcsaft::equilibrium::FreeValue::fixed(P);
        problem.feed.values={1,1e-6,1e-6};problem.source_reference=declaration;
        problem.reactions={{"H2O=Na++Cl-",{-1,1,1},1e-8,std::nullopt,std::nullopt,reference}};
        const auto solved=epcsaft::equilibrium::solve_equilibrium(model,problem);
        CHECK(solved.success);
        const auto& x=solved.mole_fractions;
        const auto center=model.state_at_pressure(T,P,x,epcsaft::DensityBranch::Liquid);
        PD dt,dp;dt.temperature=1;dp.pressure=1;dt.composition=dp.composition={0,0,0};
        const auto thermal=model.pressure_actions(T,P,x,epcsaft::DensityBranch::Liquid,{}, {dt},center.molar_density);
        const auto pressure=model.pressure_actions(T,P,x,epcsaft::DensityBranch::Liquid,{}, {dp},center.molar_density);
        const double v=1/center.molar_density, vt=-thermal.density_action[0]*v*v,vp=-pressure.density_action[0]*v*v;
        const double y=(x[1]+x[2])/2,anchor_weight=x[0]+y;
        std::array<double,3> cp{},cv{},u{};
        for(unsigned k=0;k<3;++k) {
            const double molality=k==0?1e-10:k==1?1e-12:1e-14;
            const double salt=molality*params.molar_mass_kg_per_mol[0];
            const std::vector<double> base{1/(1+2*salt),salt/(1+2*salt),salt/(1+2*salt)};
            const auto partial=[&](unsigned a,unsigned b) {
                return *anchored_reference_partial(model,T,src,base,{-1,1,1},a,b,std::nullopt).action;
            };
            const double ct=partial(1,0),ctt=partial(2,0);
            const double cpref=fixed?0:partial(0,1),ctpref=fixed?0:partial(1,1),cppref=fixed?0:partial(0,2);
            const double vg=v-y*R*T*cpref, vgt=vt-y*R*(cpref+T*ctpref),vgp=vp-y*R*T*cppref;
            const double h=anchor_weight*(-241826-R*T)+y*R*T*T*ct+R*T+center.residual_enthalpy;
            cp[k]=anchor_weight*(33.59-R)+y*R*(2*T*ct+T*T*ctt)+R+*center.residual_isobaric_heat_capacity;
            cv[k]=cp[k]+T*vgt*vgt/vgp;u[k]=h-P*vg;
            const std::string id="F9/completion/"+std::to_string(fixed)+"/"+std::to_string(k);
            recording_support::observe(id.c_str(),"completed frozen cv","J/mol/K",cv[k],std::nullopt,
                recording_support::fields({recording_support::field("cp",cp[k]),recording_support::field("u_G",u[k]),
                recording_support::field("v_G",vg),recording_support::field("v_G_T",vgt),recording_support::field("v_G_P",vgp),
                recording_support::field("quotient_magnitude",std::abs(T*vgt*vgt/vgp)),
                recording_support::field("x_water",x[0]),recording_support::field("x_Na",x[1]),recording_support::field("x_Cl",x[2]),
                recording_support::field("equilibrium_residual",solved.residual_norm()),
                recording_support::field("molality",molality),recording_support::field("source_pressure",src)}),
                "Numerical","independent-anchored-completion","case-center quotient conditioning");
            CHECK(vg>0 && vgp<0 && std::isfinite(cv[k]));
        }
        for(unsigned k=0;k<2;++k)for(bool is_cv:{false,true}) {
            const double change=std::abs((is_cv?cv:cp)[k+1]-(is_cv?cv:cp)[k]);
            const std::string id="F9/completed-change/"+std::to_string(fixed)+"/"+std::to_string(is_cv)+"/"+std::to_string(k);
            const bool ok=change<=0.01;
            recording_support::compare_custom(id.c_str(),is_cv?"completed cv change":"completed cp change","J/mol/K",change,0,
                recording_support::fields({recording_support::field("u_G_change",std::abs(u[k+1]-u[k])),recording_support::field("target",0.01)}),
                "Relationship","completed-ladder","completed heat-capacity convergence","absolute",ok,"F9-limit");
        }
    }
}
// Frozen N-P-fix vapor center and normalized tangent from the independently
// differentiated binary128 cv action in f6-a1-extraction/states.tsv.
void residual_cv_temperature_extraction() {
    using namespace epcsaft;
    constexpr double T = 260.0, V = 0.000390185267693186768959;
    const std::vector<double> x = {0.517913560605011924487,
        0.4700137196266039874, 0.012072719768384160971};
    const HelmholtzModel model(properties_support::neutral_ternary(),
        {{ParameterFamily::DispersionEnergyOverK, {1}}});
    const DensityActionDirection direction{T, 0.00371339890091887002074 / (V * V),
        {-0.194749819139095670106, 0.0800169602891089848384,
         0.114732858849986685268}, {}};
    const auto batch = model.density_actions(T, 1.0 / V, x,
        {PropertyObservable::ResidualIsochoricHeatCapacity}, {direction});
    CHECK(batch.results.front().action);
    const double action = *batch.results.front().action / T; // restore the 1 K action
    recording_support::compare("residual-cv/temperature-extraction/A1", "residual cv A1", "J/(mol K)",
        action, 0.04758556315383267,
        recording_support::fields({recording_support::field("T_K", T),
            recording_support::field("V_m3_per_mol", V), recording_support::field("x", x),
            recording_support::field("normalized_delta_T_K", direction.temperature),
            recording_support::field("normalized_delta_rho_mol_per_m3", direction.molar_density),
            recording_support::field("normalized_delta_x", direction.composition),
            recording_support::field("physical_delta_T_K", 1.0)}),
        "Numerical", "binary128 nested exact T derivatives", "Frozen residual cv first action",
        "absolute", 1e-12, 0.0, std::abs(action - 0.04758556315383267) <= 1e-12);
}

}  // namespace properties_checks
