// Frozen synthetic measurement experiment, issue #510, design sections 8–9.
#include <iostream>
#ifdef EPCSAFT_EIV
#include "binary128_support.hpp"
#endif
#include "regression/fit.hpp"
#include "surface/parameter_builder.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "core/trace.hpp"
#include "properties_support.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <random>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
using namespace epcsaft;
using namespace epcsaft::regression;
using namespace epcsaft::equilibrium;
namespace fs = std::filesystem;
#include "mea-declarations.hpp"
std::ofstream output(const fs::path& path) {
    fs::create_directories(path.parent_path()); std::ofstream stream(path); stream << std::setprecision(17); return stream;
}
template<class T> void cell(std::ostream& stream, const T& value) { stream << value; }
void cell(std::ostream& stream, const std::string& value) {
    stream << '"'; for (char c : value) { if (c == '"') stream << '"'; stream << c; } stream << '"';
}
template<class... T> void csv(std::ostream& stream, const T&... values) {
    std::size_t n = 0; ((stream << (n++ ? "," : ""), cell(stream, values)), ...); stream << '\n';
}
FitRequest ionic_request() {
    FitRequest request; request.parameters = properties_support::aqueous_electrolyte();
    auto& p = request.parameters; p.sigma_constant[1] = 3.2;
    p.packing_diameter_angstrom[1] = p.debye_huckel_diameter_angstrom[1] = 0.88*3.2; p.born_diameter_angstrom[1] = 3.5;
    request.coordinates = {{{{{ParameterFamily::PackingDiameter,{1},0.88}, {ParameterFamily::SigmaConstant,{1},1},
        {ParameterFamily::DebyeHuckelDiameter,{1},0.88}},3.0,0.2,2.4,4.4}}};
    for (int j = 0; j < 16; ++j) for (auto kind : {ObservationKind::MeanIonicActivityCoefficient, ObservationKind::OsmoticCoefficient}) {
        Observation row; row.kind = kind; row.temperature = 298.15; row.pressure = 1e5;
        row.molality = 0.1*std::pow(45.0,j/15.0); row.composition = {1,0,0};
        row.ion_indices = {1,2}; row.stoichiometry = {1,1}; row.branch = DensityBranch::Liquid;
        row.observed = 1; request.observations.push_back(row);
    }
    return request;
}
void controls(FitRequest& request, bool replay = false) {
    request.controls = {}; auto& c = request.controls; c.solver_options.values[policy::function] = c.solver_options.values[policy::gradient] = replay ? 1e-13 : 1e-12;
    c.solver_options.values[policy::parameter] = replay ? 1e-11 : 1e-10; c.maximum_iterations = 100; c.maximum_elapsed_time_seconds = 600; c.workers = 1;
    request.singular_vectors = true;
}
ObservationEvaluation evaluate(const FitRequest& request, const std::vector<double>& z) {
    std::vector<std::size_t> rows(request.observations.size()); std::iota(rows.begin(),rows.end(),0);
    FitResult failure; ObservationEvaluator evaluator(request,failure,rows); ObservationEvaluation result;
    if (!evaluator.evaluate(z,true,result)) throw std::runtime_error(failure.failures.empty() ? "point unavailable" : failure.failures.back().message);
    return result;
}
void vector(std::ostream& stream, const std::string& quantity, const std::vector<double>& values) {
    for (std::size_t i = 0; i < values.size(); ++i) csv(stream,quantity,i,values[i]);
}
void retain(const fs::path& directory, const FitResult& fit, const core::Trace& trace, double seconds) {
    auto values = output(directory/"values.csv"); csv(values,"quantity","index","value");
    vector(values,"optimizer",fit.optimizer); vector(values,"physical",fit.physical);
    vector(values,"prediction",fit.predictions); vector(values,"residual",fit.residuals);
    vector(values,"weighted_residual",fit.weighted_residuals); vector(values,"physical_jacobian",fit.physical_jacobian);
    vector(values,"optimizer_jacobian",fit.optimizer_jacobian); vector(values,"covariance_optimizer",fit.covariance.optimizer);
    vector(values,"covariance_physical",fit.covariance.physical); vector(values,"singular_values",fit.covariance.singular_values);
    vector(values,"right_vectors",fit.covariance.right_vectors); vector(values,"iteration_cost",fit.iteration_costs);
    vector(values,"left_vectors",fit.covariance.left_vectors);
    vector(values,"iteration_seconds",fit.iteration_seconds);
    csv(values,"cost",0,fit.final_cost); csv(values,"projected_gradient_norm",0,fit.covariance.projected_gradient_norm);
    auto details = output(directory/"details.txt"); details << fit.status << '\n' << fit.termination << '\n' << fit.message << '\n'
        << fit.covariance.status << '\n' << fit.covariance.assumption << '\n';
    auto failures = output(directory/"failures.csv"); csv(failures,"trial","row","status","jacobian","message","optimizer");
    for (const auto& failure : fit.failures) {
        std::ostringstream point; for (double value : failure.optimizer) point << std::setprecision(17) << value << ' ';
        csv(failures,failure.trial,failure.observation,observation_status_name(failure.status),failure.jacobian,failure.message,point.str());
    }
    auto work = output(directory/"work.csv"); csv(work,"quantity","count");
    std::map<std::string,std::size_t> counts;
    for (const auto& span : trace.spans) {
        ++counts[span.name];
        if (span.convergence) { counts[std::string(span.name)+"/iterations"] += span.convergence->iterations;
            counts[std::string(span.name)+"/evaluations"] += span.convergence->evaluations; }
    }
    for (const auto& entry : counts) csv(work,entry.first,entry.second);
    csv(work,"fit_iterations",fit.iterations); csv(work,"fit_residual_evaluations",fit.residual_evaluations);
    csv(work,"fit_jacobian_evaluations",fit.jacobian_evaluations); csv(work,"trial_failures",fit.trial_failures);
    csv(work,"trace_incomplete",trace.incomplete);
    auto diagnostics = output(directory/"diagnostics.csv"); csv(diagnostics,"seconds","rank","active_bounds","usable");
    csv(diagnostics,seconds,fit.covariance.rank,fit.covariance.active_bounds.size(),fit.usable);
}
#ifdef EPCSAFT_EIV
struct Experiment {
    std::string name; FitRequest truth; ObservationEvaluation values;
    std::vector<double> sigma; std::vector<FitRequest> datasets; std::vector<bool> valid;
    double truth_coordinate;
};
MeasuredInput measured(MeasuredInputKind kind, double sigma, double lower, double upper, std::size_t species = 0) {
    return {kind,species,sigma,lower,upper,kind == MeasuredInputKind::Temperature ? "K" : kind == MeasuredInputKind::Pressure
        ? "Pa" : kind == MeasuredInputKind::Molality ? "mol/kg" : "mol/mol"};
}
double quantity(const Observation& row, const MeasuredInput& input) {
    const bool ionic = electrolyte_observation(row.kind);
    if (input.kind == MeasuredInputKind::Temperature) return ionic ? row.temperature : row.problem.temperature.value;
    if (input.kind == MeasuredInputKind::Pressure) return ionic ? row.pressure : row.problem.pressure.value;
    if (input.kind == MeasuredInputKind::Molality) return row.molality;
    return row.problem.feed.values[input.species]/std::accumulate(row.problem.feed.values.begin(),row.problem.feed.values.end(),0.0);
}
// Physical declarations for the frozen neutral CO2/water exchange and pure-solvent salt.
void declare(Observation& row, const MeasuredInput& input, double value, double total) {
    const bool ionic = electrolyte_observation(row.kind);
    if (input.kind == MeasuredInputKind::Temperature) (ionic ? row.temperature : row.problem.temperature.value) = value;
    if (input.kind == MeasuredInputKind::Pressure) (ionic ? row.pressure : row.problem.pressure.value) = value;
    if (input.kind == MeasuredInputKind::Molality) row.molality = value;
    if (input.kind == MeasuredInputKind::MoleFraction) {
        const double change = value*total-row.problem.feed.values[input.species];
        row.problem.feed.values[input.species] += change; row.problem.feed.values[2] -= change;
    }
}
Experiment experiment(bool ionic, const fs::path& directory, bool generate) {
    Experiment experiment; auto& request = experiment.truth; experiment.name = ionic ? "ionic" : "mea";
    request = ionic ? ionic_request() : mea_request(); controls(request); request.uncertainty.emplace();
    experiment.truth_coordinate = ionic ? 1.0 : 0.0;
    for (std::size_t g = 0; g < 16; ++g) {
        MeasurementGroup group;
        if (ionic) {
            group.rows = {2*g,2*g+1}; group.inputs = {measured(MeasuredInputKind::Temperature,0.5,288.15,308.15),
                measured(MeasuredInputKind::Pressure,2000,5e4,1.5e5),
                measured(MeasuredInputKind::Molality,0.1*request.observations[2*g].molality,0.001,8)};
        } else {
            auto& row = request.observations[g]; const HelmholtzModel model(request.parameters);
            row.problem.feed.values = resolve_feed(model,row.problem.feed); row.problem.feed.kind = FeedKind::Amounts;
            group.feed_total_mol = std::accumulate(row.problem.feed.values.begin(),row.problem.feed.values.end(),0.0);
            group.rows = {g}; group.dependent_species = {2}; group.inputs = {measured(MeasuredInputKind::Temperature,1,303.15,323.15),
                measured(MeasuredInputKind::MoleFraction,0.002,0,row.problem.feed.values[0]/group.feed_total_mol+row.problem.feed.values[2]/group.feed_total_mol)};
        }
        request.uncertainty->groups.push_back(group);
    }
    request.uncertainty->output_sigma.assign(request.observations.size(),1);
    auto z = std::vector<double>(1+(ionic ? 48 : 32),0); z[0] = experiment.truth_coordinate;
    experiment.values = evaluate(request,z);
    std::ofstream truths;
    if (generate) { truths = output(directory/experiment.name/"truth.csv"); csv(truths,"row","value","sigma_y"); }
    for (std::size_t i = 0; i < request.observations.size(); ++i) {
        const double value = experiment.values.predictions[i], sigma = (ionic ? 0.005 : 0.01)*std::abs(value);
        if (!std::isfinite(sigma) || sigma == 0) throw std::runtime_error("unavailable/zero true signal");
        request.observations[i].observed = value; request.uncertainty->output_sigma[i] = sigma;
        experiment.sigma.push_back(sigma); if (generate) csv(truths,i,value,sigma);
    }
    experiment.values = evaluate(request,z);
    if (!generate) {
        experiment.datasets.assign(128,request); experiment.valid.assign(128,true);
        std::ifstream retained(directory/experiment.name/"draws.csv"); std::string line; std::getline(retained,line);
        std::size_t count = 0;
        while (std::getline(retained,line)) {
            std::replace(line.begin(),line.end(),',',' '); std::istringstream row(line);
            std::size_t r,g,k; std::string kind; double draw,truth,sigma,value; bool admitted;
            row >> r >> std::quoted(kind) >> g >> k >> draw >> truth >> sigma >> value >> admitted;
            if (!row || r >= 128) throw std::runtime_error("invalid retained draw");
            auto& data = experiment.datasets[r]; experiment.valid[r] = experiment.valid[r] && admitted;
            if (kind == "input") {
                const auto& group = data.uncertainty->groups.at(g); const auto& input = group.inputs.at(k);
                if (truth != quantity(request.observations[group.rows[0]],input) || sigma != input.sigma) throw std::runtime_error("retained input identity changed");
                for (auto i : group.rows) declare(data.observations[i],input,value,group.feed_total_mol);
            } else {
                if (truth != request.observations.at(g).observed || sigma != experiment.sigma[g]) throw std::runtime_error("retained true output changed");
                data.observations[g].observed = value;
            }
            ++count;
        }
        if (count != 128*(16*(ionic ? 3 : 2)+request.observations.size())) throw std::runtime_error("retained draws incomplete");
        for (auto& data : experiment.datasets) data.parameters = data.coordinates.apply(data.parameters,{ionic ? 0.0 : 1.0});
        return experiment;
    }
    auto draws = output(directory/experiment.name/"draws.csv"); csv(draws,"replicate","kind","group_or_row","input","normal","truth","sigma","measured","in_domain");
    std::mt19937_64 generator(ionic ? 51002 : 51001); std::normal_distribution<double> normal;
    for (std::size_t r = 0; r < 128; ++r) {
        auto data = request; bool valid = true;
        for (std::size_t g = 0; g < data.uncertainty->groups.size(); ++g) {
            const auto& group = data.uncertainty->groups[g];
            for (std::size_t k = 0; k < group.inputs.size(); ++k) {
                const auto& input = group.inputs[k]; const double truth = quantity(request.observations[group.rows[0]],input);
                const double draw = normal(generator), value = truth+input.sigma*draw;
                const bool admitted = value >= input.lower && value <= input.upper; valid &= admitted;
                csv(draws,r,"input",g,k,draw,truth,input.sigma,value,admitted);
                for (auto row : group.rows) declare(data.observations[row],input,value,group.feed_total_mol);
            }
        }
        for (std::size_t i = 0; i < data.observations.size(); ++i) {
            const double draw = normal(generator), truth = request.observations[i].observed;
            data.observations[i].observed = truth+experiment.sigma[i]*draw;
            csv(draws,r,"output",i,0,draw,truth,experiment.sigma[i],data.observations[i].observed,true);
        }
        data.parameters = data.coordinates.apply(data.parameters,{ionic ? 0.0 : 1.0});
        experiment.datasets.push_back(std::move(data)); experiment.valid.push_back(valid);
    }
    return experiment;
}
std::vector<double> complete_point(const FitResult& fit) {
    auto z = fit.optimizer;
    if (fit.reconciliation) z.insert(z.end(),fit.reconciliation->corrections.begin(),fit.reconciliation->corrections.end());
    return z;
}
std::vector<Observation> endpoint(const FitRequest& request, const FitResult& fit) {
    auto observations = request.observations;
    if (fit.reconciliation) {
        const auto& rec = *fit.reconciliation;
        for (std::size_t c = 0; c < rec.values.size(); ++c) {
            const auto& group = request.uncertainty->groups[rec.column_groups[c]];
            for (auto row : group.rows) declare(observations[row],group.inputs[rec.column_inputs[c]],rec.values[c],group.feed_total_mol);
        }
    }
    return observations;
}
void reconciliation(const fs::path& directory, const FitResult& fit) {
    if (!fit.reconciliation) return;
    const auto& rec = *fit.reconciliation; auto values = output(directory/"reconciliation.csv");
    csv(values,"column","group","input","physical_value","standardized_correction","squared_penalty");
    for (std::size_t c = 0; c < rec.values.size(); ++c)
        csv(values,c,rec.column_groups[c],rec.column_inputs[c],rec.values[c],rec.corrections[c],rec.penalties[c]);
    auto matrix = output(directory/"augmented-jacobian.csv"); csv(matrix,"row","column","value");
    const auto width = complete_point(fit).size();
    for (std::size_t i = 0; i < rec.jacobian.size(); ++i) csv(matrix,i/width,i%width,rec.jacobian[i]);
}
std::vector<double> wide_ionic_values(const FitRequest& request, const std::vector<double>& z) {
    auto row = request.observations[0]; const auto& group = request.uncertainty->groups[0];
    for (std::size_t k = 0; k < group.inputs.size(); ++k) declare(row,group.inputs[k],quantity(row,group.inputs[k])+group.inputs[k].sigma*z[k+1],0);
    auto parameters = request.coordinates.apply(request.parameters,{z[0]}); const HelmholtzModel model(parameters);
    const auto prediction = predict(model,row); if (!prediction.activity) throw std::runtime_error(prediction.message);
    EquilibriumProblem problem; problem.temperature = FreeValue::fixed(row.temperature); problem.pressure = FreeValue::fixed(row.pressure);
    ThermochemistryRecord record;
    binary128_support::WideIonic oracle(parameters,record,problem,
        formula_composition(row.composition,parameters,row.ion_indices,row.stoichiometry,row.molality)[1]);
    const auto values = oracle.activity_values(row.molality,prediction.activity->state.molar_density);
    std::vector<double> residuals;
    for (std::size_t i = 0; i < 2; ++i) residuals.push_back(double((values[i]-binary128_support::Wide(request.observations[i].observed))
        /binary128_support::Wide(request.uncertainty->output_sigma[i])));
    residuals.insert(residuals.end(),z.begin()+1,z.end()); return residuals;
}
std::string family(const Experiment& experiment, std::size_t c) {
    if (c == 0) return "model";
    const std::vector<std::string> names = experiment.name == "mea" ? std::vector<std::string>{"T","feed"} : std::vector<std::string>{"T","P","molality"};
    return names[(c-1)%names.size()];
}
bool needs_observed_estimate(const Experiment& e, std::size_t c) { return e.name == "mea" || c == 0; }
void verify(const Experiment& experiment, const fs::path& directory) {
    auto stream = output(directory/experiment.name/"verification.csv");
    csv(stream,"group","column","family","row","h","exact","difference","difference_half","error","allowance","refinement_change");
    std::map<std::string,double> rho;
    for (std::size_t g = 0; g < 16; ++g) {
        auto request = experiment.truth; const auto original = request.uncertainty->groups[g];
        request.observations.clear(); request.uncertainty->output_sigma.clear();
        for (auto row : original.rows) { request.observations.push_back(experiment.truth.observations[row]); request.uncertainty->output_sigma.push_back(experiment.sigma[row]); }
        auto group = original; std::iota(group.rows.begin(),group.rows.end(),0); request.uncertainty->groups = {group};
        std::vector<double> z(1+group.inputs.size(),0); z[0] = experiment.truth_coordinate;
        const auto exact = evaluate(request,z); const auto width = z.size();
        for (std::size_t c = 0; c < width; ++c) {
            std::vector<std::vector<double>> differences; double norm = 0;
            for (const double h : {0.02,0.01,0.005,0.0025}) {
                auto plus = z, minus = z; plus[c] += h; minus[c] -= h;
                const auto a = experiment.name == "ionic" ? wide_ionic_values(request,plus) : evaluate(request,plus).weighted;
                const auto b = experiment.name == "ionic" ? wide_ionic_values(request,minus) : evaluate(request,minus).weighted;
                std::vector<double> difference(a.size());
                for (std::size_t i = 0; i < difference.size(); ++i) difference[i] = (a[i]-b[i])/(2*h);
                differences.push_back(difference);
            }
            double best = INFINITY;
            for (std::size_t step = 0; step < 3; ++step) {
                double discrepancy = 0;
                for (std::size_t i = 0; i < exact.weighted.size(); ++i) {
                    const double j = exact.residual_optimizer_gradient[i*width+c]; norm = std::max(norm,std::abs(j));
                    discrepancy = std::max(discrepancy,std::abs(j-differences[step][i])+std::abs(differences[step][i]-differences[step+1][i]));
                    const double richardson = (4*differences[2][i]-differences[1][i])/3;
                    const double previous = (4*differences[1][i]-differences[0][i])/3;
                    const double allowance = 1e-7+5e-6*std::abs(richardson);
                    csv(stream,g,c,family(experiment,c),i,0.02/std::pow(2,step),j,differences[step][i],differences[step+1][i],
                        std::abs(j-richardson),allowance,std::abs(richardson-previous));
                    if (std::abs(j-richardson) > allowance || std::abs(richardson-previous) > allowance/4) throw std::runtime_error("campaign derivative falsifier failed");
                }
                best = std::min(best,discrepancy);
            }
            const double relative = norm ? best/norm : best ? INFINITY : 0;
            rho[family(experiment,c)] = std::max(rho[family(experiment,c)],relative);
        }
    }
    auto estimates = output(directory/experiment.name/"rho.csv"); csv(estimates,"family","relative_estimate");
    for (const auto& entry : rho) csv(estimates,entry.first,entry.second);
}
struct Errors { std::vector<std::array<double,4>> value; std::vector<std::array<double,5>> action; };
// #510/#520's finite-ladder allowance belongs to its retained Engine source, not the analytic limit.
Errors errors(const Experiment&, const FitRequest&, const FitResult&,
              const fs::path&, const std::map<std::string,double>&) {
    throw std::runtime_error("FrozenReferenceAllowanceUnavailable: analytic references have no terminal-change "
        "error estimate; reproduce #510/#520 with their recorded Engine and campaign source");
}

std::array<double,9> propagate(const Experiment& experiment, const FitResult& result, const Errors& point, const Errors& truth) {
    const auto width = complete_point(result).size(), n = experiment.sigma.size();
    const auto& j = result.reconciliation ? result.reconciliation->jacobian : result.optimizer_jacobian;
    const auto& svd = result.covariance; std::vector<double> inverse_row(width,0); std::array<double,9> h{};
    if (svd.rank != width || svd.right_vectors.size() != width*width) { h.fill(NAN); return h; }
    for (std::size_t c = 0; c < width; ++c) {
        for (std::size_t k = 0; k < width; ++k) inverse_row[c] += svd.right_vectors[k]*svd.right_vectors[c*width+k]/std::pow(svd.singular_values[k],2);
        inverse_row[c] = std::abs(inverse_row[c]);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t term = 0; term < 4; ++term)
                h[term] += inverse_row[c]*std::abs(j[i*width+c])*(point.value[i][term]+truth.value[i][term])/experiment.sigma[i];
            const double residual = std::abs(result.weighted_residuals[i]);
            for (std::size_t term = 0; term < 5; ++term) h[4+term] += inverse_row[c]*point.action[i*width+c][term]*residual;
        }
    }
    return h;
}
std::map<std::string,double> read_rho(const fs::path& path) {
    std::ifstream file(path); std::string line; std::getline(file,line); std::map<std::string,double> result;
    while (std::getline(file,line)) { std::replace(line.begin(),line.end(),',',' '); std::istringstream row(line); std::string name; double value;
        row >> std::quoted(name) >> value; if (!row || !std::isfinite(value)) throw std::runtime_error("invalid frozen rho"); result[name] = value; }
    return result;
}
void replicate(const Experiment& experiment, std::size_t r, const fs::path& directory, const Errors& truth, const std::map<std::string,double>& rho) {
    for (bool augmented : {false,true}) for (bool replay : {false,true}) {
        if (replay && r%32) continue;
        const auto path = directory/experiment.name/std::to_string(r)/(augmented ? "eiv" : "fixed")/(replay ? "replay" : "main");
        auto request = experiment.datasets[r]; controls(request,replay);
        if (!augmented) { request.uncertainty.reset(); for (double sigma : experiment.sigma) request.weights.push_back(1/(sigma*sigma)); }
        std::array<double,9> h; h.fill(NAN); double error = NAN; bool regular = false;
        try {
            if (!experiment.valid[r]) throw std::runtime_error("invalid experiment: measured input outside qualification box");
            core::Capture capture; const auto begin = std::chrono::steady_clock::now(); const auto result = fit(request);
            retain(path,result,capture.trace(),std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()); reconciliation(path,result);
            if (result.usable && !result.optimizer.empty()) error = result.optimizer[0]-experiment.truth_coordinate;
            regular = result.status == "converged" && result.covariance.status == "available" && std::isfinite(error);
            if (result.usable) h = propagate(experiment,result,errors(experiment,request,result,path,rho),truth);
        } catch (const std::exception& failure) { auto message = output(path/"exception.txt"); message << failure.what() << '\n'; }
        auto summary = output(path/"endpoint.csv"); csv(summary,"regular","signed_error","pressure","reference_value","observed_value","rounding_value","linear_action","reference_action","observed_action","rounding_action","verification_action");
        csv(summary,regular,error,h[0],h[1],h[2],h[3],h[4],h[5],h[6],h[7],h[8]);
    }
}
void correct_point(const Experiment&, std::size_t, const fs::path&, const Errors&, const std::map<std::string,double>&);
void run(Experiment& experiment, const fs::path& inputs, const fs::path& directory, const std::string& mode) {
    const auto rho = read_rho(inputs/experiment.name/"rho.csv");
    FitResult generated; generated.optimizer = {experiment.truth_coordinate}; generated.chart = experiment.values.chart;
    generated.predictions = experiment.values.predictions; generated.optimizer_jacobian = experiment.values.residual_optimizer_gradient;
    generated.reconciliation = experiment.values.reconciliation;
    const auto truth = errors(experiment,experiment.truth,generated,directory/experiment.name/"truth",rho);
    std::map<pid_t,std::size_t> children; std::size_t next = 0; const auto count = std::max(1u,std::thread::hardware_concurrency());
    while (next < 128 || !children.empty()) {
        while (next < 128 && children.size() < count) {
            const auto index = next++; const auto pid = fork();
            if (pid < 0) throw std::runtime_error("campaign fork failed");
            if (!pid) {
                if (mode == "run") replicate(experiment,index,directory,truth,rho);
                else for (auto method : {"fixed","eiv"}) for (auto control : {"main","replay"}) if (std::string(control) == "main" || index%32 == 0)
                    correct_point(experiment,index,directory/experiment.name/std::to_string(index)/method/control,truth,rho);
                _exit(0); }
            children[pid] = index;
        }
        int status; const auto pid = waitpid(-1,&status,0); if (pid < 0) throw std::runtime_error("campaign wait failed");
        std::cout << experiment.name << " replicate " << children.at(pid) << " process status " << status << '\n' << std::flush; children.erase(pid);
    }
}
std::vector<double> endpoint_values(const fs::path& path) {
    std::ifstream file(path); std::vector<double> values; std::string value,line; std::getline(file,line);
    if (line.find("regular") == 0) std::getline(file,line);
    std::replace(line.begin(),line.end(),',',' '); std::istringstream row(line);
    while (row >> value) values.push_back(std::stod(value));
    if (values.size() != 11) values.assign(11,NAN); return values;
}
std::map<std::string,std::vector<double>> retained_values(const fs::path& path) {
    std::ifstream file(path); std::string line,key,index; std::getline(file,line); std::map<std::string,std::vector<double>> values;
    while (std::getline(file,line)) { std::istringstream row(line); std::getline(row,key,','); std::getline(row,index,',');
        key.erase(std::remove(key.begin(),key.end(),'"'),key.end()); std::getline(row,line); values[key].push_back(std::stod(line)); }
    return values;
}
void correct_point(const Experiment& experiment, std::size_t r, const fs::path& path, const Errors& truth, const std::map<std::string,double>& rho) {
    auto point = endpoint_values(path/"endpoint.csv"); if (!std::isfinite(point[10])) return;
    const auto values = retained_values(path/"values.csv"); FitResult result; result.optimizer = values.at("optimizer");
    result.covariance.singular_values = values.at("singular_values"); result.covariance.right_vectors = values.at("right_vectors");
    const auto width = result.covariance.singular_values.size(), n = experiment.sigma.size(); result.covariance.rank = width;
    result.weighted_residuals = values.at("weighted_residual"); result.predictions = values.at("prediction"); auto request = experiment.datasets[r];
    if (width == 1) { request.uncertainty.reset(); for (double sigma : experiment.sigma) request.weights.push_back(1/(sigma*sigma)); result.optimizer_jacobian = values.at("optimizer_jacobian"); }
    else { const auto rows = retained_values(path/"augmented-jacobian.csv"); auto& rec = result.reconciliation.emplace();
        for (std::size_t i = 0; i < n+width-1; ++i) { const auto& row = rows.at(std::to_string(i)); rec.jacobian.insert(rec.jacobian.end(),row.begin(),row.end()); }
        std::ifstream file(path/"reconciliation.csv"); std::string line; std::getline(file,line);
        while (std::getline(file,line)) { std::replace(line.begin(),line.end(),',',' '); std::istringstream row(line); std::size_t c,g,k; double q,u,penalty; row >> c >> g >> k >> q >> u >> penalty;
            rec.column_groups.push_back(g); rec.column_inputs.push_back(k); rec.values.push_back(q); rec.corrections.push_back(u); rec.penalties.push_back(penalty); }
    }
    result.chart = evaluate_coordinates(request.parameters,request.coordinates,result.optimizer,request.holds,request.thermochemistry);
    for (const auto* name : {"endpoint.csv","action-errors.csv","prediction-errors.csv","raw-errors.csv"}) if (!fs::exists(path/(std::string("before-cold-scaling-")+name)))
        fs::copy_file(path/name,path/(std::string("before-cold-scaling-")+name));
    const auto h = propagate(experiment,result,errors(experiment,request,result,path,rho),truth);
    auto record = output(path/"endpoint.csv"); csv(record,"regular","signed_error","pressure","reference_value","observed_value","rounding_value","linear_action","reference_action","observed_action","rounding_action","verification_action");
    csv(record,point[0],point[1],h[0],h[1],h[2],h[3],h[4],h[5],h[6],h[7],h[8]);
}
void statistics(const std::string& name, const fs::path& inputs, const fs::path& directory) {
    std::array<std::vector<double>,2> errors, means; std::array<std::array<double,9>,2> components{};
    std::array<double,2> bias{},mcse{},rmse{},lower{},upper{},numerical{},replay{}; std::array<int,2> failures{};
    auto paired = output(directory/name/"paired-errors.csv"); csv(paired,"replicate","fixed","eiv","fixed_regular","eiv_regular");
    for (std::size_t r = 0; r < 128; ++r) {
        std::array<bool,2> regular;
        for (std::size_t a = 0; a < 2; ++a) {
            const auto path = directory/name/std::to_string(r)/(a ? "eiv" : "fixed"); const auto value = endpoint_values(path/"main/endpoint.csv");
            regular[a] = value[0] == 1; failures[a] += !regular[a]; errors[a].push_back(value[1]);
            for (std::size_t t = 0; t < 9; ++t) components[a][t] += value[2+t]/128;
            if (r%32 == 0) { const auto tight = endpoint_values(path/"replay/endpoint.csv");
                if (tight[0] != 1 || !std::isfinite(tight[1])) replay[a] = NAN;
                else if (std::isfinite(replay[a])) replay[a] = std::max(replay[a],std::abs(tight[1]-value[1])); }
        }
        csv(paired,r,errors[0].back(),errors[1].back(),regular[0],regular[1]);
    }
    auto resamples = output(directory/name/"bootstrap-means.csv"); csv(resamples,"resample","fixed_signed_mean","eiv_signed_mean");
    std::ifstream indices(inputs/name/"bootstrap-indices.csv"); std::string line; std::getline(indices,line); std::array<double,2> sum{}; std::size_t count = 0;
    while (std::getline(indices,line)) {
        std::replace(line.begin(),line.end(),',',' '); std::istringstream row(line); std::size_t b,k,r; row >> b >> k >> r;
        if (!row || b != count/128 || k != count%128 || r >= 128) throw std::runtime_error("invalid bootstrap indices");
        for (std::size_t a = 0; a < 2; ++a) sum[a] += errors[a][r]/128;
        if (++count%128 == 0) { csv(resamples,b,sum[0],sum[1]); for (std::size_t a = 0; a < 2; ++a) means[a].push_back(sum[a]); sum = {}; }
    }
    if (count != 2000*128) throw std::runtime_error("incomplete bootstrap indices");
    auto summary = output(directory/name/"summary.csv"); csv(summary,"method","signed_mean","mcse","absolute_bias","rmse","signed_lower","signed_upper","n_A","t_A","failed_main_fits");
    auto contributions = output(directory/name/"numerical-components.csv"); csv(contributions,"method","pressure","reference_value","observed_value","rounding_value","linear_action","reference_action","observed_action","rounding_action","verification_action");
    for (std::size_t a = 0; a < 2; ++a) {
        bias[a] = std::accumulate(errors[a].begin(),errors[a].end(),0.0)/128; double variance = 0,squares = 0;
        for (double value : errors[a]) { variance += std::pow(value-bias[a],2); squares += value*value; }
        mcse[a] = std::sqrt(variance/(127*128)); rmse[a] = std::sqrt(squares/128);
        if (std::isfinite(bias[a])) { std::sort(means[a].begin(),means[a].end()); lower[a] = means[a][12]; upper[a] = means[a][1987]; }
        else lower[a] = upper[a] = NAN;
        numerical[a] = std::accumulate(components[a].begin(),components[a].end(),0.0);
        csv(summary,a ? "eiv" : "fixed",bias[a],mcse[a],std::abs(bias[a]),rmse[a],lower[a],upper[a],numerical[a],replay[a],failures[a]);
        const auto& v = components[a]; csv(contributions,a ? "eiv" : "fixed",v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8]);
    }
    const double delta = std::abs(bias[1])-std::abs(bias[0]), df = lower[0] <= 0 && upper[0] >= 0 ? 0 : std::min(std::abs(lower[0]),std::abs(upper[0]));
    const double bound = std::max(std::abs(lower[1]),std::abs(upper[1]))-df, epsilon = 2*(numerical[0]+numerical[1]+replay[0]+replay[1]);
    auto decision = output(directory/name/"decision.csv"); csv(decision,"delta","U_Delta","epsilon_num","failed_main_fits","accepted");
    csv(decision,delta,bound,epsilon,failures[0]+failures[1],!failures[0] && !failures[1] && std::isfinite(epsilon) && bound < -epsilon);
}
#endif
int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::runtime_error("usage: campaign prepare|bootstrap|run|allowance|statistics DIRECTORY");
#ifdef EPCSAFT_EIV
        if (std::string(argv[1]) == "prepare" || std::string(argv[1]) == "bootstrap") {
            std::mt19937_64 generator(51003); std::uniform_int_distribution<std::size_t> index(0,127);
            for (bool ionic : {false,true}) {
                const std::string name = ionic ? "ionic" : "mea";
                if (std::string(argv[1]) == "prepare") {
                    auto data = experiment(ionic,argv[2],!fs::exists(fs::path(argv[2])/name/"draws.csv"));
                    verify(data,argv[2]); std::cout << name << " prepared and verified\n" << std::flush;
                }
                const auto path = fs::path(argv[2])/name/"bootstrap-indices.csv"; std::ofstream indices;
                if (!fs::exists(path)) { indices = output(path); csv(indices,"resample","draw","replicate"); }
                for (std::size_t b = 0; b < 2000; ++b) for (std::size_t k = 0; k < 128; ++k) {
                    const auto r = index(generator); if (indices) csv(indices,b,k,r);
                }
            }
            return 0;
        }
        if (std::string(argv[1]) == "run" || std::string(argv[1]) == "statistics" || std::string(argv[1]) == "allowance") {
            const fs::path leaf(argv[2]);
            for (bool ionic : {false,true}) {
                if (std::string(argv[1]) == "run") { auto data = experiment(ionic,leaf/"inputs",false); run(data,leaf/"inputs",leaf/"results","run"); }
                if (std::string(argv[1]) == "allowance") {
                    auto data = experiment(ionic,leaf/"inputs",false);
                    const auto path = leaf/"results"/data.name; if (!fs::exists(path/"truth-before-complete-norm")) fs::rename(path/"truth",path/"truth-before-complete-norm");
                    run(data,leaf/"inputs",leaf/"results","allowance");
                }
                statistics(ionic ? "ionic" : "mea",leaf/"inputs",leaf/"results");
            }
            return 0;
        }
#endif
        throw std::runtime_error("unknown mode");
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
