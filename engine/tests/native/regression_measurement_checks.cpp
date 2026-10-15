#include "regression_checks_shared.hpp"
#include "core/trace.hpp"
#include <Eigen/Dense>

namespace regression_checks {
Parameters coupled_lengths(double);
CoordinateSet coupled_coordinate();
Observation log_constant(const EquilibriumProblem&, const char*);
namespace {
MeasuredInput input(MeasuredInputKind kind, double sigma, double lower, double upper, std::size_t species = 0) {
    return {kind, species, sigma, lower, upper, kind == MeasuredInputKind::Temperature ? "K"
        : kind == MeasuredInputKind::Pressure ? "Pa" : kind == MeasuredInputKind::Molality ? "mol/kg" : "mol/mol"};
}
FitRequest ionic_measurements() {
    FitRequest request;
    request.parameters = coupled_lengths(3.2); request.coordinates = coupled_coordinate(); request.controls = budgeted_controls();
    for (const auto kind : {ObservationKind::MeanIonicActivityCoefficient, ObservationKind::OsmoticCoefficient}) {
        Observation row;
        row.kind = kind; row.temperature = 298.15; row.pressure = 1e5; row.molality = 0.5;
        row.composition = {1, 0, 0}; row.ion_indices = {1, 2}; row.stoichiometry = {1, 1}; row.branch = DensityBranch::Liquid;
        row.observed = 1.0; request.observations.push_back(row);
    }
    request.uncertainty = MeasurementUncertainty{{{{0, 1}, {
        input(MeasuredInputKind::Temperature, 0.5, 288.15, 308.15), input(MeasuredInputKind::Pressure, 2000, 5e4, 1.5e5),
        input(MeasuredInputKind::Molality, 0.05, 0.001, 8.0)}, {}, 0.0}}, {0.005, 0.005}, {}, {}};
    return request;
}
ObservationEvaluation point(const FitRequest& request, const std::vector<double>& z, bool columns = true) {
    std::vector<std::size_t> rows(request.observations.size()); std::iota(rows.begin(), rows.end(), 0);
    FitResult failures; ObservationEvaluator evaluator(request, failures, rows); ObservationEvaluation out;
    if (!evaluator.evaluate(z, columns, out)) {
        std::cerr << "point failure: " << failures.failures.back().message << '\n'; CHECK(false);
    }
    return out;
}
void derivative_probe(const std::string& id, const FitRequest& request, const std::vector<double>& z) {
    constexpr double h = 0.02;
    const auto exact = point(request, z);
    std::vector<std::vector<double>> directions(z.size(), std::vector<double>(z.size(), 0.0));
    for (std::size_t c = 0; c < z.size(); ++c) directions[c][c] = 1.0;
    directions.push_back(std::vector<double>(z.size(), 0.37));
    for (std::size_t c = 1; c < z.size(); c += 2) directions.back()[c] = -0.23;
    for (std::size_t c = 0; c < directions.size(); ++c) {
        std::vector<std::vector<double>> differences;
        for (const double step : {h, h/2, h/4}) {
            auto plus = z, minus = z;
            for (std::size_t k = 0; k < z.size(); ++k) { plus[k] += step*directions[c][k]; minus[k] -= step*directions[c][k]; }
            const auto a = point(request, plus, false), b = point(request, minus, false);
            std::vector<double> difference(a.weighted.size());
            for (std::size_t i = 0; i < difference.size(); ++i) difference[i] = (a.weighted[i] - b.weighted[i])/(2*step);
            differences.push_back(std::move(difference));
        }
        for (std::size_t row = 0; row < exact.weighted.size(); ++row) {
            const double reference = (4*differences[2][row] - differences[1][row])/3;
            const double coarse = (4*differences[1][row] - differences[0][row])/3;
            double action = 0.0;
            for (std::size_t k = 0; k < z.size(); ++k) action += exact.residual_optimizer_gradient[row*z.size()+k]*directions[c][k];
            const double allowance = 1e-7 + 5e-6*std::abs(reference);
            const auto name = id + "/" + std::to_string(c) + "/" + std::to_string(row);
            reference_scale_agrees(name, "standardized_residual_action", "1", action, reference, 1e-7, 5e-6,
                fields({field("h", h), field("refinement_change", std::abs(reference-coarse))}),
                "exact physical input and mixed columns agree with rebuilt states", "measurement-columns");
            CHECK(std::abs(reference-coarse) <= allowance/4);
        }
    }
}
FitRequest linear_measurements(bool shared = false) {
    FitRequest request;
    request.parameters = properties_support::neutral_ternary(); request.controls = budgeted_controls();
    request.controls.solver_options.values[epcsaft::policy::function] = request.controls.solver_options.values[epcsaft::policy::gradient] = 1e-12;
    request.controls.solver_options.values[epcsaft::policy::parameter] = 1e-10;
    Coordinate slope; slope.reaction = ReactionCoordinate{"reaction_a", ReactionCoefficient::D};
    slope.origin = 0.0; slope.scale = 0.001; slope.lower = 0.001; slope.upper = 0.01;
    request.coordinates.coordinates = {slope}; request.uncertainty.emplace();
    auto problem = equilibrium_support::lattice_problem(equilibrium_support::PhaseSet::OneVapor,
        equilibrium_support::AmountPattern::AllFree, equilibrium_support::IntensivePattern::BothFixed,
        equilibrium_support::Chemistry::TwoReactions);
    problem.reactions[0].equilibrium_constant.reset();
    problem.reactions[0].correlation = ReactionLogPolynomial{0, 0, 0, 0.0045, 300, 250, 450, equilibrium::EOS_STANDARD_STATE_ID};
    const std::vector<double> temperatures{300, 320, 340}, values{1.52, 1.59, 1.72};
    for (std::size_t i = 0; i < temperatures.size(); ++i) {
        problem.temperature = FreeValue::fixed(shared ? 320 : temperatures[i]);
        auto row = log_constant(problem, "reaction_a"); row.observed = values[i]; request.observations.push_back(row);
        request.uncertainty->output_sigma.push_back(0.02);
        if (!shared) request.uncertainty->groups.push_back({{i}, {input(MeasuredInputKind::Temperature, 2, 250, 450)}, {}, 0});
    }
    if (shared) request.uncertainty->groups = {{{0, 1, 2}, {input(MeasuredInputKind::Temperature, 2, 250, 450)}, {}, 0}};
    return request;
}
}  // namespace

void check_measurement_columns() {
    auto ionic = ionic_measurements(); derivative_probe("ionic", ionic, {1, 0, 0, 0});
    FitRequest density_request = ionic;
    auto& row = density_request.observations[0]; row.kind = ObservationKind::Density; row.mass_basis = true;
    row.composition = {0.98, 0.01, 0.01}; density_request.observations.resize(1);
    density_request.uncertainty->output_sigma = {1.0};
    density_request.uncertainty->groups = {{{0}, {input(MeasuredInputKind::Pressure, 2000, 5e4, 1.5e5),
        input(MeasuredInputKind::MoleFraction, 0.001, 0.005, 0.02, 1)}, {0, 2}, 0}};
    derivative_probe("density-mass", density_request, {1, 0, 0});
    auto cp_request = density_request;
    cp_request.observations[0].kind = ObservationKind::IsobaricHeatCapacity; cp_request.observations[0].mass_basis = false;
    auto thermo = properties_support::p3_record(properties_support::p3_forms()[0].second);
    thermo.component_intervals.resize(3, thermo.component_intervals.front()); cp_request.thermochemistry = thermo;
    cp_request.uncertainty->groups[0].inputs.push_back(input(MeasuredInputKind::Temperature, 0.5, 288.15, 308.15));
    derivative_probe("pressure-cp", cp_request, {1,0,0,0});
    auto apparent = cp_request;
    auto& cp = apparent.observations[0]; cp.kind = ObservationKind::EquilibriumIsobaricHeatCapacity;
    cp.problem.temperature = FreeValue::fixed(cp.temperature); cp.problem.pressure = FreeValue::fixed(cp.pressure);
    cp.problem.feed = {equilibrium::FeedKind::Amounts, cp.composition, {}, 0};
    equilibrium::PhaseSpec liquid; liquid.name = "liquid"; liquid.kind = equilibrium::PhaseKind::Liquid;
    cp.problem.phases = {liquid}; apparent.uncertainty->groups[0].feed_total_mol = 1;
    derivative_probe("equilibrium-cp", apparent, {1,0,0,0});
    FitRequest equilibrium_request;
    equilibrium_request.parameters = properties_support::neutral_binary(); equilibrium_request.controls = budgeted_controls();
    equilibrium_request.coordinates.coordinates = {{{{ParameterFamily::Kij, {1, 2}, 1}}, 0, 0.01, -0.5, 0.5}};
    Observation equilibrium_row; equilibrium_row.kind = ObservationKind::PartialPressure;
    equilibrium_row.problem = equilibrium_support::bubble_at(180.0); equilibrium_row.phase = 1; equilibrium_row.species = {0};
    equilibrium_request.observations = {equilibrium_row};
    const auto feed = equilibrium::resolve_feed(HelmholtzModel(equilibrium_request.parameters), equilibrium_row.problem.feed);
    equilibrium_request.uncertainty = MeasurementUncertainty{{{{0}, {input(MeasuredInputKind::Temperature, 0.2, 165, 195),
        input(MeasuredInputKind::MoleFraction, 0.002, 0.05, 0.95, 0)}, {1}, std::accumulate(feed.begin(), feed.end(), 0.0)}}, {10000}, {}, {}};
    derivative_probe("equilibrium", equilibrium_request, {0, 0, 0});
    derivative_probe("source-law", linear_measurements(), {4.5, 0, 0, 0});
}

void check_measurement_fit() {
    for (const bool shared : {false, true}) for (const auto& controls :
        std::vector<std::array<double,3>>{{1e-12,1e-12,1e-10}, {1e-13,1e-13,1e-11}, {1e-14,1e-14,1e-12}, {1e-16,1e-14,1e-12}}) {
        auto request = linear_measurements(shared);
        request.controls.solver_options.values[epcsaft::policy::function] = controls[0]; request.controls.solver_options.values[epcsaft::policy::gradient] = controls[1];
        request.controls.solver_options.values[epcsaft::policy::parameter] = controls[2];
        const bool verify = controls[0] == (shared ? 1e-16 : 1e-14);
        core::Capture capture;
        const auto result = fit(request);
        CHECK(result.usable && result.status == "converged" && result.reconciliation && result.covariance.status == "available");
        const auto& rec = *result.reconciliation;
        const std::vector<double> temperatures = shared ? std::vector<double>{320, 320, 320} : std::vector<double>{300, 320, 340};
        const std::vector<double> y{1.52, 1.59, 1.72};
        double xx = 0, xy = 0, yy = 0;
        for (std::size_t i = 0; i < y.size(); ++i) { xx += temperatures[i]*temperatures[i]; xy += temperatures[i]*y[i]; yy += y[i]*y[i]; }
        const double lambda = 0.0001, d = yy - lambda*xx;
        const double slope = shared ? std::accumulate(y.begin(), y.end(), 0.0)/960.0 : (d + std::sqrt(d*d + 4*lambda*xy*xy))/(2*xy);
        const std::string id = shared ? "linear-shared" : "linear";
        double max_q_error = 0, max_j_error = 0;
        std::vector<double> expected_q;
        const std::size_t latent = shared ? 1 : 3, width = latent + 1;
        Eigen::MatrixXd j = Eigen::MatrixXd::Zero(3 + latent, width);
        double cost = 0.0;
        for (std::size_t c = 0; c < latent; ++c) {
            double numerator = temperatures[c]/4, denominator = 0.25;
            for (std::size_t i = 0; i < 3; ++i) if (shared || i == c) { numerator += slope*y[i]/0.0004; denominator += slope*slope/0.0004; }
            const double q = numerator/denominator, u = (q-temperatures[c])/2;
            max_q_error = std::max(max_q_error, std::abs(rec.values[c]-q));
            expected_q.push_back(q);
            cost += 0.5*u*u; j(3+c, 1+c) = 1;
            for (std::size_t i = 0; i < 3; ++i) if (shared || i == c) {
                const double residual = (slope*q-y[i])/0.02; cost += 0.5*residual*residual;
                j(i, 0) = 0.001*q/0.02; j(i, 1+c) = slope*2/0.02;
            }
        }
        for (Eigen::Index i = 0; i < j.rows(); ++i) for (Eigen::Index c = 0; c < j.cols(); ++c) {
            const double error = std::abs(rec.jacobian[static_cast<std::size_t>(i)*width+static_cast<std::size_t>(c)] - j(i,c));
            max_j_error = std::max(max_j_error, error);
        }
        const Eigen::MatrixXd h = j.transpose()*j;
        const double covariance = h.ldlt().solve(Eigen::VectorXd::Unit(width, 0))[0];
        Eigen::VectorXd r(3+latent);
        for (std::size_t i = 0; i < 3; ++i) r[static_cast<Eigen::Index>(i)] = result.weighted_residuals[i];
        for (std::size_t c = 0; c < latent; ++c) r[static_cast<Eigen::Index>(3+c)] = rec.corrections[c];
        Eigen::MatrixXd actual(3+latent, width);
        for (std::size_t i = 0; i < 3+latent; ++i) for (std::size_t c = 0; c < width; ++c) actual(i,c) = rec.jacobian[i*width+c];
        const Eigen::VectorXd gradient = actual.transpose()*r;
        const Eigen::VectorXd step = (actual.transpose()*actual).ldlt().solve(gradient);
        std::cout << std::setprecision(17) << "LINEAR_REFINEMENT shared=" << shared << " function_tolerance=" << controls[0]
            << " slope_error=" << std::abs(result.physical[0]-slope) << " max_q_error=" << max_q_error
            << " cost_error=" << std::abs(result.final_cost-cost) << " max_J_error=" << max_j_error
            << " covariance_error=" << std::abs(result.covariance.optimizer[0]-covariance)
            << " gradient_norm=" << gradient.cwiseAbs().maxCoeff() << " newton_step_norm=" << step.cwiseAbs().maxCoeff() << '\n';
        std::cout << "CERES_TERMINATION shared=" << shared << " function_tolerance=" << controls[0] << " type="
            << result.termination << " message=" << result.message << '\n';
        for (const auto& event : capture.trace().events) if (std::string(capture.trace().spans[event.span].name) == "fit_iteration") {
            std::cout << "CERES_ITERATION shared=" << shared << " function_tolerance=" << controls[0] << " " << event.name;
            for (const auto value : event.values) std::cout << ' ' << value;
            std::cout << '\n';
        }
        if (shared) {
            Eigen::MatrixXd analytic = Eigen::MatrixXd::Zero(4,2); Eigen::VectorXd residual(4);
            for (std::size_t i = 0; i < 3; ++i) {
                analytic(i,0) = 0.001*rec.values[0]/0.02; analytic(i,1) = 2*result.physical[0]/0.02;
                residual[i] = (result.physical[0]*rec.values[0]-y[i])/0.02;
            }
            analytic(3,1) = 1; residual[3] = (rec.values[0]-320)/2;
            const Eigen::VectorXd analytic_step = (analytic.transpose()*analytic).ldlt().solve(analytic.transpose()*residual);
            const std::vector<double> z{result.optimizer[0],rec.corrections[0]};
            const auto center = point(request,z);
            for (std::size_t i = 0; i < 4; ++i) {
                CHECK(center.weighted[i] == r[i]);
                CHECK(std::abs(r[i]-residual[i]) <= 1e-10);
                for (std::size_t c = 0; c < 2; ++c) {
                    CHECK(std::abs(actual(i,c)-analytic(i,c)) <= 1e-10*std::max(1.0,std::abs(analytic(i,c))));
                    auto p = z, m = z; p[c] += 0.02; m[c] -= 0.02;
                    const double difference = (point(request,p,false).weighted[i]-point(request,m,false).weighted[i])/0.04;
                    CHECK(std::abs(actual(i,c)-difference) <= 1e-7+5e-6*std::abs(difference));
                }
            }
            CHECK((step-analytic_step).cwiseAbs().maxCoeff() <= 1e-10);
            const auto corrected = point(request,{z[0]-step[0],z[1]-step[1]});
            std::cout << "SHARED_GN code_dq_K=" << 2*step[1] << " analytic_dq_K=" << 2*analytic_step[1]
                << " corrected_q_K=" << corrected.reconciliation->values[0] << " corrected_slope=" << corrected.chart.physical[0] << '\n';
            if (verify) {
                CHECK(std::abs(corrected.reconciliation->values[0]-320) <= 320e-10);
                CHECK(std::abs(corrected.chart.physical[0]-slope) <= 1e-10);
                CHECK(std::abs(corrected.cost-cost) <= 1e-10);
                for (std::size_t i = 0; i < 4; ++i) for (std::size_t c = 0; c < 2; ++c)
                    CHECK(std::abs(corrected.reconciliation->jacobian[i*2+c]-j(i,c)) <= 1e-10*std::max(1.0,std::abs(j(i,c))));
                const Eigen::Map<const Eigen::Matrix<double,4,2,Eigen::RowMajor>> jc(corrected.reconciliation->jacobian.data());
                const double corrected_cov = (jc.transpose()*jc).ldlt().solve(Eigen::Vector2d::Unit(0))[0];
                CHECK(std::abs(corrected_cov-covariance) <= 1e-10*std::max(1.0,std::abs(covariance)));
            }
        }
        if (verify) {
            for (std::size_t c = 0; c < expected_q.size(); ++c)
                reference_scale_agrees(id+"/q/"+std::to_string(c), "temperature", "K", rec.values[c], expected_q[c], 0, 1e-10, "", "analytic eliminated temperature", "measurement-fit");
            reference_scale_agrees(id, "slope", "1/K", result.physical[0], slope, 1e-10, 0, "", "analytic errors-in-variables slope", "measurement-fit");
            reference_scale_agrees(id+"/cost", "objective", "1", result.final_cost, cost, 1e-10, 0, "", "complete analytic objective", "measurement-fit");
            CHECK(max_j_error <= 1e-10*std::max(1.0,j.cwiseAbs().maxCoeff()));
            CHECK(std::abs(result.covariance.optimizer[0]-covariance) <= 1e-10*std::max(1.0,std::abs(covariance)));
            CHECK(!result.covariance.variance_factor);
        }
    }
}

void check_measurement_reuse() {
    const auto request = ionic_measurements();
    const std::vector<double> base{1, 0, 0, 0};
    const auto compare = [&](const std::string& id, const ObservationEvaluation& a, const ObservationEvaluation& b, double molality) {
        CHECK(a.statuses == b.statuses);
        for (std::size_t row = 0; row < a.predictions.size(); ++row) {
            const double floor = row%2 ? 1e-12/(2*0.0180153*molality) : 1e-12;
            CHECK(agree(id+"/value/"+std::to_string(row), "prediction", "1", a.predictions[row], b.predictions[row], floor, 1e-10,
                "", "changed-input reuse agrees with cold evaluation", "measurement-reuse"));
            for (std::size_t k = 0; k < a.state_densities[row].size(); ++k)
                CHECK(agree(id+"/rho/"+std::to_string(row)+"/"+std::to_string(k), "density", "mol/m3",
                    a.state_densities[row][k], b.state_densities[row][k], 1e-12, 1e-10, "", "finite and reference densities", "measurement-reuse"));
            for (std::size_t k = 0; k < a.reference_values[row].size(); ++k)
                CHECK(agree(id+"/reference/"+std::to_string(row)+"/"+std::to_string(k), "ln_phi_reference", "1",
                    a.reference_values[row][k], b.reference_values[row][k], 1e-12, 1e-10, "", "terminal and pure reference values", "measurement-reuse"));
        }
        for (std::size_t i = 0; i < a.residual_optimizer_gradient.size(); ++i)
            CHECK(reference_scale_agrees(id+"/column/"+std::to_string(i), "standardized_column", "1",
                a.residual_optimizer_gradient[i], b.residual_optimizer_gradient[i], 1e-7, 5e-6, "",
                "same-point promotion retains exact parameter and input columns", "measurement-reuse"));
    };
    for (const auto& shift : std::vector<std::vector<double>>{{0,1,0,0}, {0,0,1,0}, {0,0,0,1}, {0,1,1,1}}) {
        FitResult failures; ObservationEvaluator reused(request, failures, {0,1}); ObservationEvaluation initial, value, promoted;
        CHECK(reused.evaluate(base, true, initial)); auto z = base;
        for (std::size_t c = 0; c < z.size(); ++c) z[c] += shift[c];
        CHECK(reused.evaluate(z, false, value)); CHECK(reused.evaluate(z, true, promoted));
        const auto cold = point(request, z);
        compare("changed-"+std::to_string(shift[1])+std::to_string(shift[2])+std::to_string(shift[3]), promoted, cold, 0.5+0.05*shift[3]);
        CHECK(value.predictions == promoted.predictions);
    }
    auto two = request;
    two.observations.insert(two.observations.end(), request.observations.begin(), request.observations.end());
    two.uncertainty->output_sigma = {0.005,0.005,0.005,0.005};
    auto second = two.uncertainty->groups[0]; second.rows = {2,3};
    second.inputs[0].sigma = 5; second.inputs[1].sigma = 20000;
    two.uncertainty->groups.push_back(second);
    auto separate = request; separate.uncertainty->groups[0].inputs = second.inputs;
    const std::vector<ObservationEvaluation> cold{point(request, base), point(separate, base)};
    for (int order = 0; order < 2; ++order) {
        if (order) std::reverse(two.uncertainty->groups.begin(), two.uncertainty->groups.end());
        const auto shared = point(two, {1,0,0,0,0,0,0});
        for (std::size_t row = 0; row < 4; ++row) {
            CHECK(shared.predictions[row] == cold[row/2].predictions[row%2]);
            CHECK(shared.reference_values[row] == cold[row/2].reference_values[row%2]);
            for (std::size_t c = 0; c < 6; ++c) {
                const auto g = shared.reconciliation->column_groups[c], j = shared.reconciliation->column_inputs[c];
                if (two.uncertainty->groups[g].rows.front()/2 != row/2) continue;
                CHECK(reference_scale_agrees("sigma-order/"+std::to_string(order)+"/"+std::to_string(row)+"/"+std::to_string(j),
                    "standardized_input_column", "1", shared.residual_optimizer_gradient[row*7+1+c],
                    cold[row/2].residual_optimizer_gradient[(row%2)*4+1+j], 1e-7, 5e-6, "",
                    "same physical reference admits each consuming sigma independently", "measurement-reuse"));
                CHECK(shared.input_gradient[row*6+c] == cold[row/2].input_gradient[(row%2)*3+j]);
            }
        }
    }
    const HelmholtzModel model(request.parameters, request.coordinates.active_parameters());
    const auto prediction = predict(model, request.observations[0]); CHECK(prediction.reference && prediction.activity);
    const auto action = [&](double sigma) { return electrolyte_input_action(model, 298.15, 1e5, {1,0,0}, {1,2}, {1,1}, 0.5,
        *prediction.reference, *prediction.activity, {1,0,0}, sigma); };
    const auto raw = action(0);
    const double delta = std::max(std::abs(raw.reference_difference[0]), std::abs(raw.reference_difference[1])); CHECK(delta > 0);
    const auto low = action(0.5*5e-5/delta), high = action(2*5e-5/delta);
    CHECK(low.status == PropertyActionStatus::Available && !low.witness_evaluated && high.witness_evaluated);
    CHECK(low.reference_column == raw.reference_column && high.reference_column == raw.reference_column);
    std::cout << std::setprecision(17) << "REFERENCE_INPUT_ADMISSION raw_difference=" << delta
        << " sigma_low=" << 0.5*5e-5/delta << " sigma_high=" << 2*5e-5/delta
        << " high_available=" << (high.status == PropertyActionStatus::Available) << " cause=" << high.message << '\n';
}

void check_measurement_admission() {
    FitRequest roundoff;
    roundoff.parameters = properties_support::neutral_ternary();
    roundoff.coordinates.coordinates = {{{{ParameterFamily::SegmentCount,{0},1}},1,0.1,0.5,5}};
    Observation face; face.kind = ObservationKind::Density; face.temperature = 400; face.pressure = 1e5;
    face.branch = DensityBranch::Vapor; face.composition = {0.30000000000000004,0.2,0.49999999999999994};
    roundoff.observations = {face};
    roundoff.uncertainty = MeasurementUncertainty{{{{0},{input(MeasuredInputKind::MoleFraction,0.002,0,0.8,0)},{2},0}},{1},{},{}};
    CHECK(face.composition[2]-(0.8-face.composition[0]) < 0);
    const auto admitted = point(roundoff,{roundoff.coordinates.start(roundoff.parameters)[0],0},false);
    CHECK(admitted.statuses[0] == ObservationStatus::Available);
    const auto request = linear_measurements();
    const auto reject = [&](FitRequest candidate, const std::string& cause, const std::string& order = "") {
        const auto result = fit(candidate);
        truth("measurement-refusal/"+cause+order, "invalid measurement declarations retain their cause",
            result.status == "invalid_setup" && result.message.find(cause) != std::string::npos,
            fields({field("message", result.message)}), "measurement-admission");
    };
    auto item = request; item.uncertainty->output_sigma.clear(); reject(item, "output sigma");
    item = request; item.uncertainty->groups[0].inputs[0].sigma = -1; reject(item, "nonnegative");
    item = request; item.uncertainty->groups[0].inputs[0].unit = "C"; reject(item, "unit");
    item = request; item.uncertainty->correlation = std::vector<double>{}; reject(item, "correlated");
    item = request; item.uncertainty->covariance = std::vector<double>{}; reject(item, "correlated");
    item = request; item.weights = {1,2,1}; reject(item, "relative weights");
    item = ionic_measurements(); item.observations[0].form = ResidualForm::LogRatio; reject(item, "LogRatio");
    item = request; item.uncertainty->groups[1].rows = {0}; reject(item, "overlapping");
    item = request; item.uncertainty->groups[0].inputs.push_back(item.uncertainty->groups[0].inputs[0]); reject(item, "duplicate");
    item = ionic_measurements(); item.observations[1].pressure += 1; reject(item, "conflicting");
    item = ionic_measurements();
    item.uncertainty->groups[0].inputs = {input(MeasuredInputKind::Molality, 0.05, 0, 8)};
    reject(item, "molality box leaves the admitted domain", "/gamma-first");
    item.uncertainty->groups[0].rows = {1,0};
    reject(item, "molality box leaves the admitted domain", "/osmotic-first");
    item = ionic_measurements(); item.observations.resize(1); item.observations[0].kind = ObservationKind::Density;
    item.observations[0].composition = {0.98,0.01,0.01}; item.uncertainty->output_sigma = {1};
    item.uncertainty->groups = {{{0}, {input(MeasuredInputKind::MoleFraction, 0.001, 0, 0.6, 1)}, {0,2}, 0}};
    reject(item, "nonphysical feed box");
    item = FitRequest{}; item.parameters = properties_support::neutral_binary(); item.controls = budgeted_controls();
    item.coordinates.coordinates = {{{{ParameterFamily::Kij,{1,2},1}},0,0.01,-0.5,0.5}};
    Observation free_pressure; free_pressure.kind = ObservationKind::EquilibriumPressure; free_pressure.problem = equilibrium_support::bubble_at(180);
    item.observations = {free_pressure};
    item.uncertainty = MeasurementUncertainty{{{{0}, {input(MeasuredInputKind::Pressure, 1000, 1e4, 1e6)}, {}, 0}}, {1}, {}, {}};
    reject(item, "free or absent intensive");
    item = ionic_measurements(); item.uncertainty->groups[0].inputs[0].sigma = 1e20;
    FitResult failures; ObservationEvaluator failed_action(item, failures, {0,1}); ObservationEvaluation unavailable;
    CHECK(!failed_action.evaluate({1,0,0,0}, true, unavailable));
    CHECK(failures.failures.back().status == ObservationStatus::ActionUnavailable);
    CHECK(failures.failures.back().message.find("reference input action") != std::string::npos);
    std::cout << "REFERENCE_ACTION_REFUSAL " << failures.failures.back().message << '\n';
    auto weighted = request; weighted.uncertainty.reset(); weighted.weights = {2500,2500,2500};
    auto zero = request;
    for (auto& group : zero.uncertainty->groups) for (auto& input : group.inputs) input.sigma = 0;
    const auto a = fit(weighted), b = fit(zero);
    CHECK(a.usable && b.usable && a.status == b.status && b.reconciliation && b.reconciliation->corrections.empty());
    CHECK(a.physical == b.physical && a.predictions == b.predictions && a.weighted_residuals == b.weighted_residuals);
    CHECK(a.optimizer_jacobian == b.optimizer_jacobian && a.covariance.optimizer == b.covariance.optimizer);
    CHECK(a.covariance.variance_factor && !b.covariance.variance_factor);
    auto ordinary = weighted; ordinary.weights.clear();
    auto empty = ordinary; empty.uncertainty.emplace();
    const auto plain = fit(ordinary), canonical = fit(empty);
    CHECK(plain.physical == canonical.physical && plain.predictions == canonical.predictions && plain.optimizer_jacobian == canonical.optimizer_jacobian);
    CHECK(plain.residual_evaluations == canonical.residual_evaluations && plain.jacobian_evaluations == canonical.jacobian_evaluations);
    CHECK(plain.covariance.optimizer == canonical.covariance.optimizer && plain.covariance.variance_factor == canonical.covariance.variance_factor);
    auto report = request; report.weights = {1,1,0};
    const auto reports = fit(report); CHECK(reports.usable && reports.reconciliation && reports.reconciliation->corrections.size() == 2);
    CHECK(reports.statuses.back() == ObservationStatus::Available && std::isfinite(reports.predictions.back()));
    auto bound = request; bound.uncertainty->groups[0].inputs[0].upper = 300.1;
    const auto active = fit(bound); CHECK(active.usable && active.covariance.status == "active_bound");
    auto deficient = request; deficient.observations.resize(1); deficient.uncertainty->groups.resize(1); deficient.uncertainty->output_sigma.resize(1);
    Coordinate intercept; intercept.reaction = ReactionCoordinate{"reaction_a", ReactionCoefficient::A};
    intercept.scale = 0.1; intercept.lower = -1; intercept.upper = 1;
    deficient.coordinates.coordinates.push_back(intercept);
    const auto rank = fit(deficient); CHECK(rank.usable && rank.covariance.status == "rank_deficient");
    std::cout << "MEASUREMENT_LIMITS zero_sigma_equal=1 empty_equal=1 report_latent=" << reports.reconciliation->corrections.size()
        << " active=" << active.covariance.active_bounds.size() << " deficient_rank=" << rank.covariance.rank << '\n';
}
}  // namespace regression_checks
