// #572 frozen repeated-sample comparison: G (declared correlation R) against D (independent sigmas) on the same
// 256 correlated ln K draws; retained draws and bootstrap indices, then 512 native fits and the paired statistics.
//   tools/native-driver regression_checks analyses/2026-correlated-measurement-errors/scripts/campaign.cpp
//   driver <leaf> [prepare]   (prepare: retain draws and bootstrap indices only)
#include "regression_checks_shared.hpp"
#include <ceres/version.h>
#include <random>
#include <sstream>

using namespace regression_checks;
namespace {
constexpr std::size_t N = 256, B = 4096;
constexpr double A_TRUE = 0.1, D_LAW = 0.0045;
const Eigen::Vector3d T(280, 320, 380), SIGMA(0.01, 0.04, 0.02);
#ifdef __clang__
constexpr const char* COMPILER = "clang ";
#else
constexpr const char* COMPILER = "gcc ";
#endif
Eigen::Matrix3d correlation() { Eigen::Matrix3d r; r << 1, .8, .3, .8, 1, .2, .3, .2, 1; return r; }
std::vector<std::vector<double>> read(const std::string& path) {
    std::vector<std::vector<double>> rows;
    std::ifstream in(path);
    for (std::string line; std::getline(in, line);) {
        if (line.empty() || !std::isdigit(static_cast<unsigned char>(line[0]))) continue;
        std::vector<double> row; std::stringstream fields(line);
        for (std::string item; std::getline(fields, item, ',');) row.push_back(std::stod(item));
        rows.push_back(std::move(row));
    }
    return rows;
}
// Draws and indices are generated once, before any fit; an existing file is the experiment and is only read.
std::vector<std::vector<double>> retained(const std::string& path, const std::function<void(std::ofstream&)>& generate) {
    if (!std::ifstream(path)) { std::ofstream out(path); out << std::setprecision(17); generate(out); }
    return read(path);
}
double quantile(const std::vector<double>& sorted, double p) {  // Hyndman-Fan type 7, zero-based
    const double h = double(sorted.size() - 1) * p, j = std::floor(h);
    const auto k = static_cast<std::size_t>(j);
    return sorted[k] + (h - j) * (sorted[std::min(k + 1, sorted.size() - 1)] - sorted[k]);
}
}  // namespace
int main(int argc, char** argv) {
    const std::string leaf = argv[1];
    const Eigen::Matrix3d r = correlation(), c = SIGMA.asDiagonal() * r * SIGMA.asDiagonal();
    const Eigen::Matrix3d lr = r.llt().matrixL();
    const auto draws = retained(leaf + "/inputs/draws.csv", [&](std::ofstream& out) {
        std::mt19937_64 engine(0x572); std::normal_distribution<double> normal(0.0, 1.0);
        out << "replicate,n0,n1,n2,error0,error1,error2,y0,y1,y2\n";
        for (std::size_t k = 0; k < N; ++k) {
            Eigen::Vector3d n; for (int i = 0; i < 3; ++i) n[i] = normal(engine);
            const Eigen::Vector3d e = SIGMA.asDiagonal() * (lr * n), y = Eigen::Vector3d::Constant(A_TRUE) + D_LAW * T + e;
            out << k << ',' << n[0] << ',' << n[1] << ',' << n[2] << ',' << e[0] << ',' << e[1] << ',' << e[2] << ',' << y[0] << ',' << y[1] << ',' << y[2] << '\n';
        }
    });
    const auto indices = retained(leaf + "/inputs/bootstrap-indices.csv", [&](std::ofstream& out) {
        std::mt19937_64 engine(0x572b); std::uniform_int_distribution<std::size_t> pick(0, N - 1);
        for (std::size_t b = 0; b < B; ++b) for (std::size_t k = 0; k < N; ++k) out << pick(engine) << (k + 1 < N ? ',' : '\n');
    });
    if (draws.size() != N || indices.size() != B) throw std::runtime_error("retained draws or bootstrap indices incomplete");
    if (argc > 2 && std::string(argv[2]) == "prepare") return 0;
    // Independent linear estimators from the declared values: a_G = C^-1 1 / 1^T C^-1 1, a_D = D^-2 1 / 1^T D^-2 1.
    const Eigen::Vector3d ones = Eigen::Vector3d::Ones(), cg = c.ldlt().solve(ones), cd = SIGMA.cwiseAbs2().cwiseInverse();
    const Eigen::Vector3d ag = cg / cg.sum(), ad = cd / cd.sum();
    const double covariance_ref[2] = {1.0 / cg.sum(), 1.0 / cd.sum()}, population[2] = {ag.dot(c * ag), ad.dot(c * ad)};
    std::ofstream endpoints(leaf + "/results/endpoints.csv"), failures(leaf + "/results/failures.csv");
    endpoints << std::setprecision(17) << "replicate,method,converged,usable,interior_full_rank,within_tolerance,A,A_reference,raw_error,"
                 "step_bound,allowance,predicted_decrease,rounding_noise,covariance,covariance_reference,final_cost,projected_gradient,"
                 "iterations,residual_evaluations,jacobian_evaluations,trial_failures\n";
    failures << "replicate,method,kind,detail,message\n";
    std::vector<double> error[2], reference_error[2], tau[2];
    bool complete = true, frozen_complete = true;
    for (std::size_t k = 0; k < N; ++k) {
        const std::vector<double> y(draws[k].begin() + 7, draws[k].end());
        const Eigen::Vector3d z = Eigen::Map<const Eigen::Vector3d>(y.data()) - D_LAW * T;
        for (int m = 0; m < 2; ++m) {
            auto request = intercept_request({T[0], T[1], T[2]}, y, {SIGMA[0], SIGMA[1], SIGMA[2]});
            if (m == 0) request.uncertainty->correlation = std::vector<double>(r.data(), r.data() + 9);  // symmetric: order-free
            const auto fit_result = fit(request);
            if (k == 0 && m == 0) {  // the resolved policy every fit uses, and the toolchain that compiled this driver
                std::ofstream provenance(leaf + "/results/campaign-provenance.txt");
                provenance << std::setprecision(17) << "compiler=" << COMPILER << __VERSION__ << "\nstandard_library=libstdc++ "
                           << _GLIBCXX_RELEASE << ' ' << __GLIBCXX__ << "\ncplusplus=" << __cplusplus << "\neigen="
                           << EIGEN_WORLD_VERSION << '.' << EIGEN_MAJOR_VERSION << '.' << EIGEN_MINOR_VERSION << "\nceres="
                           << CERES_VERSION_STRING << "\nmaximum_iterations=" << request.controls.maximum_iterations
                           << "\nmaximum_elapsed_time_seconds=" << request.controls.maximum_elapsed_time_seconds
                           << "\nworkers=" << request.controls.workers << '\n';
                for (std::size_t i = 0; i < policy::SOLVER_POLICIES.size(); ++i)
                    provenance << policy::SOLVER_POLICIES[i]->name << '=' << fit_result.solver_options.value(i) << '\n';
            }
            // Lead-approved allowance (2026-10-16, cost resolution): add |one Gauss-Newton step| from the returned endpoint,
            // and record its predicted cost decrease beside a one-roundoff-per-standardized-output cost change.
            double jj = 0.0, g = 0.0, noise = 0.0;
            for (std::size_t i = 0; i < 3; ++i) {
                const double j = fit_result.optimizer_jacobian.at(i), w = fit_result.weighted_residuals.at(i);
                jj += j * j; g += j * w; noise += std::abs(w) * std::numeric_limits<double>::epsilon() * std::abs(y[i]) / SIGMA[Eigen::Index(i)];
            }
            const double step = std::abs(g / jj), frozen = 1e-10 * std::max(1.0, std::abs((m == 0 ? ag : ad).dot(z)));
            const double reference = (m == 0 ? ag : ad).dot(z), allowance = frozen + step;
            const bool usable = fit_result.usable && !fit_result.physical.empty() && std::isfinite(fit_result.physical[0]);
            const bool regular = fit_result.covariance.status == "available";
            const double a = usable ? fit_result.physical[0] : NAN, variance = regular ? fit_result.covariance.physical[0] : NAN;
            const bool covariance_ok = std::abs(variance - covariance_ref[m]) <= 1e-10 * std::max(1.0, covariance_ref[m]);
            const bool within = usable && regular && covariance_ok && std::abs(a - reference) <= allowance;
            frozen_complete &= usable && fit_result.status == "converged" && regular && covariance_ok && std::abs(a - reference) <= frozen;
            complete &= usable && fit_result.status == "converged" && regular && within;
            endpoints << k << ',' << (m ? "D" : "G") << ',' << (fit_result.status == "converged") << ',' << usable << ',' << regular << ','
                      << within << ',' << a << ',' << reference << ',' << std::abs(a - reference) << ',' << step << ',' << allowance << ','
                      << 0.5 * g * g / jj << ',' << noise << ',' << variance << ',' << covariance_ref[m] << ',' << fit_result.final_cost
                      << ',' << fit_result.covariance.projected_gradient_norm << ',' << fit_result.iterations << ','
                      << fit_result.residual_evaluations << ',' << fit_result.jacobian_evaluations << ',' << fit_result.trial_failures << '\n';
            if (fit_result.status != "converged" || !regular)
                failures << k << ',' << (m ? "D" : "G") << ",fit," << fit_result.status << ' ' << fit_result.termination << ' '
                         << fit_result.covariance.status << ",\"" << fit_result.message << "\"\n";
            for (const auto& trial : fit_result.failures)
                failures << k << ',' << (m ? "D" : "G") << ",trial," << trial.trial << ",\"" << trial.message << "\"\n";
            error[m].push_back(a - A_TRUE); reference_error[m].push_back(reference - A_TRUE); tau[m].push_back(allowance);
        }
    }

    // Paired statistic over all 256 pairs: d_r = e_G^2 - e_D^2, its 4096 percentile resamples and the numerical margin eta.
    std::ofstream summary(leaf + "/results/summary.csv"), paired(leaf + "/results/paired-errors.csv"), means(leaf + "/results/bootstrap-means.csv");
    summary << std::setprecision(17) << "quantity,value\n";
    paired << std::setprecision(17) << "replicate,error_G,error_D,difference,reference_error_G,reference_error_D\n";
    std::vector<double> d(N);
    double eta = 0.0, frozen_eta = 0.0, delta = 0.0;
    for (std::size_t k = 0; k < N; ++k) {
        d[k] = error[0][k] * error[0][k] - error[1][k] * error[1][k]; delta += d[k] / N;
        double bound = 0.0;
        double frozen_bound = 0.0;
        for (int m = 0; m < 2; ++m) {
            bound += tau[m][k] * (2 * std::abs(reference_error[m][k]) + tau[m][k]);
            const double t0 = 1e-10 * std::max(1.0, std::abs(reference_error[m][k] + A_TRUE));
            frozen_bound += t0 * (2 * std::abs(reference_error[m][k]) + t0);
        }
        eta = std::max(eta, bound); frozen_eta = std::max(frozen_eta, frozen_bound);
        paired << k << ',' << error[0][k] << ',' << error[1][k] << ',' << d[k] << ',' << reference_error[0][k] << ',' << reference_error[1][k] << '\n';
    }
    for (int m = 0; m < 2; ++m) {
        double bias = 0.0, variance = 0.0, mse = 0.0;
        for (const double e : error[m]) { bias += e / N; mse += e * e / N; }
        for (const double e : error[m]) variance += (e - bias) * (e - bias) / N;
        const std::string name = m ? "D" : "G";
        summary << "bias_" << name << ',' << bias << "\nvariance_" << name << ',' << variance << "\nmse_" << name << ',' << mse
                << "\nmse_identity_residual_" << name << ',' << mse - (bias * bias + variance) << "\npopulation_mse_" << name << ','
                << population[m] << "\ncovariance_reference_" << name << ',' << covariance_ref[m] << '\n';
    }
    std::vector<double> resampled(B);
    means << std::setprecision(17) << "resample,delta\n";
    for (std::size_t b = 0; b < B; ++b) {
        for (const double k : indices[b]) resampled[b] += d[static_cast<std::size_t>(k)] / N;
        means << b << ',' << resampled[b] << '\n';
    }
    std::sort(resampled.begin(), resampled.end());
    const double lower = quantile(resampled, 0.025), upper = quantile(resampled, 0.975);
    const char* decision = !complete ? "incomplete" : upper < -eta ? "reduction" : lower > eta ? "higher_mse_G" : "inconclusive";
    summary << "delta," << delta << "\ninterval_lower," << lower << "\ninterval_upper," << upper << "\neta," << eta
            << "\nall_numerical_conditions," << complete << "\ndecision," << decision << "\nfrozen_rule_numerical_conditions," << frozen_complete
            << "\nfrozen_rule_eta," << frozen_eta << '\n';
    std::cout << "decision=" << decision << " delta=" << delta << " interval=[" << lower << ',' << upper << "] eta=" << eta << '\n';
    return 0;
}
