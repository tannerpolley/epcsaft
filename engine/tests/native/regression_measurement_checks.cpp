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
constexpr double D_LAW = 0.0045, RHO = 0.5;  // #572 joint case: R_ij = rho^|i-j| over y, then T
Eigen::MatrixXd kac(Eigen::Index m, double rho, bool inverse) {  // rho^|i-j|, or its tridiagonal inverse
    Eigen::MatrixXd r(m, m);
    for (Eigen::Index i = 0; i < m; ++i) for (Eigen::Index j = 0; j < m; ++j)
        r(i, j) = !inverse ? std::pow(rho, double(std::abs(i - j)))
            : (i == j ? (i == 0 || i == m - 1 ? 1.0 : 1 + rho*rho) : std::abs(i - j) == 1 ? -rho : 0.0) / (1 - rho*rho);
    return r;
}
std::vector<double> row_major(const Eigen::MatrixXd& m) {
    std::vector<double> out; for (Eigen::Index i = 0; i < m.rows(); ++i) for (Eigen::Index j = 0; j < m.cols(); ++j) out.push_back(m(i, j));
    return out;
}
bool scaled(const std::string& id, const char* quantity, const char* unit, double model, double reference, double s) {
    return reference_scale_agrees(id, quantity, unit, model, reference, 1e-10*std::max(s, std::abs(reference)), 0, "",
        "independent joint covariance reference", "measurement-correlation");
}
}  // namespace

void check_measurement_correlation() {
    const Eigen::Vector3d t(300, 320, 340), y = Eigen::Vector3d::Constant(0.1) + D_LAW*t + Eigen::Vector3d(0.02, -0.03, 0.01);
    Eigen::VectorXd sd(6); sd << 0.02, 0.02, 0.02, 2, 2, 2;
    const auto joint = [&](double rho, bool covariance) {
        auto request = intercept_request({300, 320, 340}, {y[0], y[1], y[2]}, {0.02, 0.02, 0.02});
        for (std::size_t i = 0; i < 3; ++i) request.uncertainty->groups.push_back({{i}, {input(MeasuredInputKind::Temperature, 2, 250, 450)}, {}, 0});
        const Eigen::MatrixXd r = kac(6, rho, false);
        if (covariance) request.uncertainty->covariance = row_major(sd.asDiagonal() * r * sd.asDiagonal());
        else request.uncertainty->correlation = row_major(r);
        return request;
    };
    for (const double rho : {0.0, RHO}) for (const bool covariance : {false, true}) {
        const auto id = std::string(covariance ? "covariance" : "correlation") + "/rho=" + std::to_string(rho);
        // Independent normal equations: beta = (A, dT), H = [[1, d I], [0, I]], b = (y - d T_M, 0), C^-1 = D^-1 R^-1 D^-1.
        const Eigen::MatrixXd w = sd.cwiseInverse().asDiagonal() * kac(6, rho, true) * sd.cwiseInverse().asDiagonal();
        Eigen::MatrixXd h = Eigen::MatrixXd::Zero(6, 4); Eigen::VectorXd b = Eigen::VectorXd::Zero(6);
        for (Eigen::Index i = 0; i < 3; ++i) { h(i, 0) = 1; h(i, 1 + i) = D_LAW; h(3 + i, 1 + i) = 1; b[i] = y[i] - D_LAW*t[i]; }
        const Eigen::MatrixXd information = h.transpose()*w*h;
        const Eigen::VectorXd beta = information.ldlt().solve(h.transpose()*w*b), e = h*beta - b;
        const Eigen::Vector3d ey = e.head(3) / 0.02;
        const double cost = 0.5*e.dot(w*e), output_cost = 0.5*ey.dot(kac(3, rho, true)*ey);
        const auto result = fit(joint(rho, covariance));
        CHECK(result.usable && result.status == "converged" && result.reconciliation && result.covariance.status == "available");
        CHECK(result.covariance.assumption.find("joint") != std::string::npos && !result.covariance.variance_factor);
        const auto& rec = *result.reconciliation;
        const double variance = information.inverse()(0, 0);
        std::cout << std::setprecision(17) << "JOINT_REFERENCE " << id << " A=" << result.physical[0] << " reference_A=" << beta[0]
            << " variance=" << result.covariance.optimizer[0] << " reference_variance=" << variance << " cost=" << result.final_cost
            << " reference_cost=" << cost << " iterations=" << result.iterations << " projected_gradient=" << result.covariance.projected_gradient_norm
            << " message=" << result.message << '\n';
        // Lead-approved allowance (2026-10-16, cost resolution): each fitted coordinate adds its component of one
        // Gauss-Newton step from the returned endpoint, z = L^-1 (r/sigma_y, u) and the returned complete whitened columns.
        Eigen::VectorXd ee(6); Eigen::MatrixXd jr(6, 4);
        for (Eigen::Index i = 0; i < 6; ++i) {
            ee[i] = i < 3 ? result.residuals[std::size_t(i)] / 0.02 : rec.corrections[std::size_t(i - 3)];
            for (Eigen::Index c = 0; c < 4; ++c) jr(i, c) = rec.jacobian[std::size_t(i*4 + c)];
        }
        Eigen::VectorXd ez = ee;
        for (Eigen::Index i = 5; i > 0; --i) ez[i] = (ez[i] - rho*ez[i - 1]) / std::sqrt(1 - rho*rho);
        const Eigen::VectorXd step = (jr.transpose()*jr).ldlt().solve(jr.transpose()*ez);
        const double decrease = 0.5*step.dot(jr.transpose()*ez);
        const auto within = [&](const std::string& name, const char* quantity, const char* unit, double model, double reference, double s, double component) {
            std::cout << std::setprecision(17) << "JOINT_STEP_BOUND " << id << '/' << name << " raw_error=" << std::abs(model - reference)
                << " step_component=" << component << " predicted_decrease=" << decrease << '\n';
            reference_scale_agrees(id + "/" + name, quantity, unit, model, reference, 1e-10*std::max(s, std::abs(reference)) + component, 0,
                fields({field("step_component", component)}), "independent joint reference within the one-step stationarity allowance",
                "measurement-correlation");
        };
        within("A", "A", "1", result.physical[0], beta[0], 1, std::abs(step[0]));
        for (std::size_t c = 0; c < 3; ++c) {
            const auto k = Eigen::Index(1 + c);
            within("dT/" + std::to_string(c), "temperature_correction", "K", rec.values[c] - t[k - 1], beta[k], 2, 2*std::abs(step[k]));
            within("u/" + std::to_string(c), "standardized_correction", "1", rec.corrections[c], beta[k] / 2, 1, std::abs(step[k]));
        }
        // Lead decision 2026-10-16: the split is not stationary, so it is checked against the independent decomposition at
        // the returned endpoint; the reference-point differences are first order in the endpoint offset and only printed.
        const Eigen::Vector3d ye = ee.head(3);
        const double endpoint_cost = 0.5*ee.dot(kac(6, rho, true)*ee), endpoint_output = 0.5*ye.dot(kac(3, rho, true)*ye);
        std::cout << std::setprecision(17) << "COST_SPLIT " << id << " output=" << rec.output_cost << " endpoint_reference=" << endpoint_output
            << " reference_point=" << output_cost << " input=" << rec.input_cost << " endpoint_reference=" << endpoint_cost - endpoint_output
            << " reference_point=" << cost - output_cost << '\n';
        scaled(id + "/cost", "objective", "1", result.final_cost, cost, 1);
        scaled(id + "/output_cost", "output_cost", "1", rec.output_cost, endpoint_output, 1);
        scaled(id + "/input_cost", "input_cost", "1", rec.input_cost, endpoint_cost - endpoint_output, 1);
        scaled(id + "/variance/optimizer", "A_variance", "1", result.covariance.optimizer[0], variance, 1);
        scaled(id + "/variance/physical", "A_variance", "1", result.covariance.physical[0], variance, 1);
    }
    // Nonstationary point: raw prediction columns [1, 2d I] and the independently whitened complete objective.
    const std::vector<double> z0{0, 0.25, -0.5, 0.75};
    const auto request = joint(RHO, false);
    const auto at = point(request, z0);
    Eigen::VectorXd e(6); Eigen::MatrixXd j = Eigen::MatrixXd::Zero(6, 4);
    for (Eigen::Index i = 0; i < 3; ++i) {
        e[i] = (D_LAW*(t[i] + 2*z0[std::size_t(1 + i)]) - y[i]) / 0.02; e[3 + i] = z0[std::size_t(1 + i)];
        j(i, 0) = 50; j(i, 1 + i) = 2*D_LAW/0.02; j(3 + i, 1 + i) = 1;
    }
    const auto whiten = [](Eigen::MatrixXd v) {  // z_0 = e_0, z_i = (e_i - rho e_(i-1)) / sqrt(1 - rho^2)
        for (Eigen::Index i = v.rows() - 1; i > 0; --i) v.row(i) = (v.row(i) - RHO*v.row(i - 1)) / std::sqrt(1 - RHO*RHO);
        return v;
    };
    const Eigen::VectorXd z = whiten(e); const Eigen::MatrixXd jz = whiten(j);
    for (std::size_t i = 0; i < 3; ++i) for (std::size_t c = 0; c < 4; ++c)
        scaled("point/raw/" + std::to_string(i) + "/" + std::to_string(c), "raw_prediction_column", "1", at.optimizer_gradient[i*4 + c],
            c == 0 ? 1.0 : c == 1 + i ? 2*D_LAW : 0.0, c == 0 ? 1 : 2*D_LAW);
    for (std::size_t i = 0; i < 6; ++i) {
        scaled("point/whitened_residual/" + std::to_string(i), "whitened_residual", "1", at.weighted[i], z[Eigen::Index(i)], 1);
        for (std::size_t c = 0; c < 4; ++c) {
            const double reference = jz(Eigen::Index(i), Eigen::Index(c));
            scaled("point/whitened_column/" + std::to_string(i) + "/" + std::to_string(c), "whitened_column", "1", at.residual_optimizer_gradient[i*4 + c], reference, 1);
            CHECK(at.reconciliation->jacobian[i*4 + c] == at.residual_optimizer_gradient[i*4 + c]);
        }
    }
    scaled("point/whitened_A/0", "whitened_column", "1", at.residual_optimizer_gradient[0], 50, 1);
    scaled("point/whitened_A/1", "whitened_column", "1", at.residual_optimizer_gradient[4], 50/std::sqrt(3.0), 1);
    const Eigen::Vector3d ey = e.head(3);
    const double cost = 0.5*e.dot(kac(6, RHO, true)*e), output_cost = 0.5*ey.dot(kac(3, RHO, true)*ey);
    scaled("point/cost", "objective", "1", at.cost, cost, 1);
    scaled("point/output_cost", "output_cost", "1", at.reconciliation->output_cost, output_cost, 1);
    scaled("point/input_cost", "input_cost", "1", at.reconciliation->input_cost, cost - output_cost, 1);
    for (std::size_t c = 0; c < 3; ++c) {
        scaled("point/penalty/" + std::to_string(c), "penalty", "1", at.reconciliation->penalties[c], z[Eigen::Index(3 + c)]*z[Eigen::Index(3 + c)], 1);
        CHECK(at.reconciliation->corrections[c] == z0[1 + c]);
    }
    std::cout << std::setprecision(17) << "CORRELATED_POINT raw_A=" << at.optimizer_gradient[0] << ',' << at.optimizer_gradient[4] << ','
        << at.optimizer_gradient[8] << " whitened_A=" << at.residual_optimizer_gradient[0] << ',' << at.residual_optimizer_gradient[4] << '\n';
    derivative_probe("correlated", request, z0);
    // Admission: refused before fitting, never repaired; roundoff-level asymmetry within 64 m eps is admitted.
    const auto refused = [&](const std::string& cause, const std::function<void(FitRequest&)>& edit) {
        auto item = joint(RHO, false); edit(item); const auto result = fit(item);
        truth("correlation-refusal/" + cause, "invalid measurement matrices retain their cause",
            result.status == "invalid_setup" && result.message.find(cause) != std::string::npos, fields({field("message", result.message)}), "measurement-correlation");
    };
    refused("both", [&](auto& item) { item.uncertainty->covariance = item.uncertainty->correlation; });
    refused("dimension", [&](auto& item) { item.uncertainty->correlation->resize(25); });
    refused("finite", [&](auto& item) { (*item.uncertainty->correlation)[1] = NAN; });
    refused("symmetric", [&](auto& item) { (*item.uncertainty->correlation)[1] = 0.6; });
    refused("positive definite", [&](auto& item) { (*item.uncertainty->correlation)[1] = (*item.uncertainty->correlation)[6] = 1; });
    refused("diagonal", [&](auto& item) { item.uncertainty->covariance = row_major(sd.asDiagonal() * kac(6, RHO, false) * sd.asDiagonal());
        item.uncertainty->correlation.reset(); (*item.uncertainty->covariance)[0] *= 1.01; });
    auto roundoff = joint(RHO, false); (*roundoff.uncertainty->correlation)[1] += 4e-14;
    CHECK(fit(roundoff).status == "converged");
    // Finite triangular operands whose squared norm overflows, after an admitted LLT.
    FitRequest overflow; overflow.parameters = properties_support::ethane(); overflow.controls = budgeted_controls();
    overflow.controls.maximum_iterations = 1;
    overflow.coordinates.coordinates = {{{{ParameterFamily::DispersionEnergyOverK, {0}}}, 191.42, 10, 180, 200}};
    for (double t : {300., 320.}) {
        auto row = density(t, 1e5, {1}, DensityBranch::Vapor); row.form = ResidualForm::Difference; row.scale = 1; row.observed = t == 300 ? 1e153 : -1e153;
        overflow.observations.push_back(row);
    }
    overflow.uncertainty.emplace(); overflow.uncertainty->output_sigma = {1, 1};
    const double r = 1 - 1e-12; overflow.uncertainty->correlation = {1, r, r, 1};
    for (bool detailed : {false, true}) {
        core::Capture capture(detailed); const auto result = fit(overflow);
        CHECK(result.status == "initial_evaluation_failed" && result.trial_failures == 1 && result.jacobian_evaluations == 1);
        CHECK(result.first_numerical_failure && result.first_numerical_failure->stop == core::Stop::nonfinite);
        CHECK(result.first_numerical_failure->cause == "whitened cost not finite");
        CHECK(!result.first_numerical_failure->rank);
        if (detailed) {
            bool solved = false, cost = false;
            for (const auto& event : capture.trace().events) {
                if (event.name == std::string("solved_operands")) { solved = true; for (double v : event.values) CHECK(std::isfinite(v)); }
                if (event.name == std::string("cost")) { cost = true; CHECK(!std::isfinite(event.values[0])); }
            }
            CHECK(solved && cost);
        }
    }
    overflow.uncertainty->correlation = {1, 1, 1, 1};
    const auto llt = fit(overflow);
    CHECK(llt.status == "invalid_setup" && llt.first_numerical_failure && llt.first_numerical_failure->stop == core::Stop::infeasible);
    CHECK(!llt.first_numerical_failure->rank && std::isnan(llt.first_numerical_failure->residual));
    // Permuted groups with an intervening zero-weight output and a fixed zero-sigma input equal the canonical order.
    auto canonical = intercept_request({300, 320, 340}, {y[0], y[1], y[2]}, {0.02, 0.02, 0.02});
    canonical.uncertainty->groups = {{{0}, {input(MeasuredInputKind::Temperature, 2, 250, 450)}, {}, 0},
        {{1}, {input(MeasuredInputKind::Temperature, 0, 320, 320)}, {}, 0}, {{2}, {input(MeasuredInputKind::Temperature, 2, 250, 450)}, {}, 0}};
    const Eigen::MatrixXd rc = kac(5, RHO, false);  // y300, y320, y340, T300, T340
    canonical.uncertainty->correlation = row_major(rc);
    auto permuted = intercept_request({340, 330, 300, 320}, {y[2], 1.6, y[0], y[1]}, {0.02, 0.02, 0.02, 0.02});
    permuted.weights = {1, 0, 1, 1};
    permuted.uncertainty->groups = {{{2}, {input(MeasuredInputKind::Temperature, 2, 250, 450)}, {}, 0},
        {{3}, {input(MeasuredInputKind::Temperature, 0, 320, 320)}, {}, 0}, {{0}, {input(MeasuredInputKind::Temperature, 2, 250, 450)}, {}, 0}};
    const std::array<Eigen::Index, 5> order{2, 0, 1, 3, 4};  // compact rows y340, y300, y320, T300, T340
    Eigen::MatrixXd rp(5, 5);
    for (Eigen::Index i = 0; i < 5; ++i) for (Eigen::Index k = 0; k < 5; ++k) rp(i, k) = rc(order[std::size_t(i)], order[std::size_t(k)]);
    permuted.uncertainty->correlation = row_major(rp);
    const auto a = fit(canonical), p = fit(permuted);
    CHECK(a.usable && p.usable && a.status == "converged" && p.status == "converged" && p.reconciliation && a.reconciliation);
    const auto& pr = *p.reconciliation;
    CHECK(p.weighted_residuals.size() == 4 && p.weighted_residuals[1] == 0.0 && pr.jacobian.size() == 6*3);
    CHECK(pr.jacobian[3] == 0.0 && pr.jacobian[4] == 0.0 && pr.jacobian[5] == 0.0);
    CHECK((pr.column_groups == std::vector<std::size_t>{0, 2}) && (pr.column_inputs == std::vector<std::size_t>{0, 0}));
    CHECK((pr.row_groups == std::vector<std::size_t>{2, 4, 0, 1, 0, 2}) && p.covariance.objective_rows.size() == 5);
    scaled("order/A", "A", "1", p.physical[0], a.physical[0], 1);
    scaled("order/cost", "objective", "1", p.final_cost, a.final_cost, 1);
    scaled("order/variance", "A_variance", "1", p.covariance.optimizer[0], a.covariance.optimizer[0], 1);
    for (std::size_t c = 0; c < 2; ++c) {
        scaled("order/u/" + std::to_string(c), "standardized_correction", "1", pr.corrections[c], a.reconciliation->corrections[c], 1);
        scaled("order/penalty/" + std::to_string(c), "penalty", "1", pr.penalties[c], a.reconciliation->penalties[c], 1);
    }
    scaled("order/output_cost", "output_cost", "1", pr.output_cost, a.reconciliation->output_cost, 1);
    const std::array<std::size_t, 3> rows{2, 3, 0};  // canonical output i at permuted row
    for (std::size_t i = 0; i < 3; ++i)
        scaled("order/residual/" + std::to_string(i), "ln_K_residual", "1", p.residuals[rows[i]], a.residuals[i], 0.02);
}

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
    item = request; item.uncertainty->correlation = std::vector<double>{}; reject(item, "dimension");
    item = request; item.uncertainty->covariance = std::vector<double>{}; reject(item, "dimension");
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
