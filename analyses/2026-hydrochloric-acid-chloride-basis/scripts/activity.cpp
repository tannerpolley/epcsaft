#include "regression/fit.hpp"
#include <algorithm>
#include <cstring>
#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
using namespace epcsaft;
namespace r = epcsaft::regression;
using Rows = std::vector<std::map<std::string, std::string>>;
Rows read(const std::filesystem::path& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("missing input " + path.string());
    std::string line, field; std::getline(f, line); if (!line.empty() && line.back() == '\r') line.pop_back(); std::stringstream header(line);
    std::vector<std::string> names; while (std::getline(header, field, ',')) names.push_back(field);
    Rows rows;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::stringstream stream(line); std::map<std::string, std::string> row;
        for (const auto& name : names) { if (!std::getline(stream, field, ',') && !(name == names.back() && line.back() == ',')) throw std::runtime_error("short input row"); row[name] = field; }
        rows.push_back(row);
    }
    return rows;
}
std::string quote(std::string value) {
    std::string out = "\"";
    for (char c : value) out += c == '"' ? "\"\"" : std::string(1, c);
    return out + "\"";
}
int main(int argc, char** argv) try {
    if (argc != 3 && argc != 4 && argc != 6) return 2;
    const bool supplied_endpoint = argc == 6;
    if (supplied_endpoint) argc = 4;
    const std::filesystem::path leaf(argv[1]), inputs = leaf / "inputs", results(argv[2]);
    std::ifstream maps("/proc/self/maps"); std::string mapping; bool loaded = false;
    std::ofstream identity(results / "native-mappings.txt");
    while (std::getline(maps, mapping)) if (mapping.find("libepcsaft_native.so") != std::string::npos) {
        const auto path = mapping.substr(mapping.find('/')); identity << path << '\n';
        if (!std::filesystem::equivalent(path, leaf / "../../build/engine-meson/shared/libepcsaft_native.so")) throw std::runtime_error("unexpected native library mapping");
        loaded = true;
    }
    if (!loaded) throw std::runtime_error("missing native library mapping");
    auto p = declared_parameters();
    const bool replacement_reference = mixture_reference == "Hamer-Wu1972";
    const bool trace_required = verify_lengths || replacement_reference;
    const std::size_t count = p.component_count(), H = ion, Cl = chloride;
    if (supplied_endpoint) {
        const double u = std::stod(argv[4]), c = std::stod(argv[5]);
        if (!std::isfinite(u) || !std::isfinite(c)) throw std::runtime_error("nonfinite endpoint");
        p.k_ij[H*count+water] = p.k_ij[water*count+H] = u;
        p.k_ij[H*count+Cl] = p.k_ij[Cl*count+H] = c;
    }
    const double k = pair_correction(temperature_K, p, temperature_law_values(temperature_K, p), TemperatureLaw::Kij, H, Cl);
    const double energy = std::sqrt(p.dispersion_energy_over_k[H] * p.dispersion_energy_over_k[Cl]) * (1 - k);
    if (k == 1 && energy != 0) throw std::runtime_error("nonzero H-Cl cross dispersion energy at k=1");
    if (!sulfate_training) std::ofstream(results / (ion_label + "-Cl-dispersion.csv")) << "k_ion_Cl,cross_energy_over_k_K\n" << std::setprecision(17) << k << ',' << energy << '\n';
    std::vector<double> solvent(count, 0.0); solvent[water] = 1.0;
    if (argc == 4 && std::string(argv[3]) == "binary-water") {
        if (!supplied_endpoint && !packet_readback) p.k_ij[H*count+water] = p.k_ij[water*count+H] = -1.0706252324767123;
        if (!supplied_endpoint && !packet_readback && k != 0) throw std::runtime_error("selected H-Cl convention");
        HelmholtzModel model(p);
        const auto reference = solvent_reference(model, temperature_K, pressure_Pa, solvent);
        std::ofstream out(results / "water.csv");
        out << std::setprecision(17) << "molality,K_measured,K_adjusted,aw_unadjusted,phi_unadjusted,aw_adjusted,phi_adjusted,aw_prediction,error,status,cause,density_mol_m3,pressure_slope\n";
        const auto rows = read(inputs / "messnaoui2020-hcl.csv");
        if (rows.size() != 14) throw std::runtime_error("Messnaoui inventory");
        for (const auto& row : rows) {
            const double m = std::stod(row.at("molality"));
            auto x = solvent; x[water] = 1/p.molar_mass_kg_per_mol[water]; x[H] = x[Cl] = m;
            const double total = std::accumulate(x.begin(), x.end(), 0.0); for (auto& v : x) v /= total;
            const auto value = solvent_activity(model, temperature_K, pressure_Pa, x, reference);
            const double aw = value.solvent_activity.value_or(NAN);
            for (const auto& field : {"molality","K_measured","K_adjusted","aw_unadjusted","phi_unadjusted","aw_adjusted","phi_adjusted"}) out << row.at(field) << ',';
            out << aw << ',' << aw-std::stod(row.at("aw_adjusted")) << ',' << int(value.status) << ',' << quote(value.cause) << ',' << value.state.molar_density << ',' << value.state.pressure_density_derivative << '\n';
        }
        return 0;
    }
    if (argc == 4 && std::string(argv[3]) == "numerical") {
        if (!supplied_endpoint && !packet_readback) p.k_ij[H*count+water] = p.k_ij[water*count+H] = -1.0706252324767123;
        HelmholtzModel model(p);
        const auto reference = infinite_dilution_reference(model, temperature_K, pressure_Pa, solvent, {H,Cl}, {1,1});
        const auto wr = solvent_reference(model, temperature_K, pressure_Pa, solvent);
        std::ofstream out(results / "numerical.csv");
        out << std::setprecision(17) << "background,m_HCl,m_background,log_gamma,expected,error,phi,aw,density,pressure_slope,density_residual,reference_density,reference_residual,reference_status,reference_cause,passes\n";
        bool complete = true;
        for (std::size_t background : {std::size_t(water),std::size_t(lithium),std::size_t(magnesium),std::size_t(4)})
            for (const auto& a : background == std::size_t(water) ? std::vector<std::array<double,2>>{{.001,0},{.2,0},{1,0},{3,0},{7,0}} : std::vector<std::array<double,2>>{{.2,.4},{.7,1.1},{1e-15,.4}}) {
                std::vector<double> m(count,0); m[H]=a[0]; if (background != std::size_t(water)) m[background]=a[1];
                m[Cl]=a[0]+p.charge_number[background]*a[1];
                const auto full = electrolyte_activity(model,temperature_K,pressure_Pa,solvent,{H,Cl},{1,1},a[0],reference,std::nullopt,m);
                const auto scalar = electrolyte_activity(model,temperature_K,pressure_Pa,solvent,{H,Cl},{1,1},a[0],reference);
                const auto& logs = *full.state.log_fugacity_coefficient;
                const double total = std::accumulate(m.begin(),m.end(),0.0), mass=p.molar_mass_kg_per_mol[water];
                const double expected = background == std::size_t(water) ? scalar.mean_ionic_log_molality_scale : .5*(require_component_property(logs[H])-reference.log_fugacity[0]+require_component_property(logs[Cl])-reference.log_fugacity[1])-std::log(1+mass*total);
                auto x=m; x[water]=1/mass; for(auto& v:x) v /= 1/mass+total;
                const auto aw=solvent_activity(model,temperature_K,pressure_Pa,x,wr,full.state);
                const double phi=full.osmotic_coefficient && std::holds_alternative<double>(*full.osmotic_coefficient) ? std::get<double>(*full.osmotic_coefficient) : NAN;
                const double error=std::abs(full.mean_ionic_log_molality_scale-expected);
                const bool ok=std::isfinite(phi) && aw.solvent_activity && error<=1e-10 && std::abs(std::inner_product(m.begin(),m.end(),p.charge_number.begin(),0.0))<=1e-10 && reference.solvent_state && reference.solvent_state->density_convergence && full.state.density_convergence && full.state.pressure_density_derivative>0;
                out << background << ',' << a[0] << ',' << a[1] << ',' << full.mean_ionic_log_molality_scale << ',' << expected << ',' << error << ',' << phi << ',' << aw.solvent_activity.value_or(NAN) << ',' << full.state.molar_density << ',' << full.state.pressure_density_derivative << ',' << (full.state.density_convergence ? full.state.density_convergence->residual : NAN) << ',' << (reference.solvent_state ? reference.solvent_state->molar_density : NAN) << ',' << (reference.solvent_state->density_convergence ? reference.solvent_state->density_convergence->residual : NAN) << ',' << int(reference.solvent_status) << ',' << quote(reference.solvent_cause) << ',' << ok << '\n';
                complete=complete && ok;
            }
        std::ofstream training_out(results / "training.csv"); training_out << std::setprecision(17) << "molality,phi,gamma\n";
        for (const auto& row:read(inputs / training)) if(std::stod(row.at("molality"))<=7) {
            const double m=std::stod(row.at("molality")); const auto a=electrolyte_activity(model,temperature_K,pressure_Pa,solvent,{H,Cl},{1,1},m,reference);
            training_out << m << ',' << std::get<double>(*a.osmotic_coefficient) << ',' << a.mean_ionic_molality_scale << '\n';
        }
        return complete ? 0 : 1;
    }
    if (argc == 4 && std::string(argv[3]) == "report-binary") {
        if (!supplied_endpoint && !packet_readback) p.k_ij[H*count+water] = p.k_ij[water*count+H] = -1.0706252324767123;
        HelmholtzModel model(p); const auto ref=solvent_reference(model,temperature_K,pressure_Pa,solvent);
        std::ofstream out(results / "report.csv"); out << std::setprecision(17) << "source,molality,aw,phi,cause\n";
        for (const auto& file:{"el-guendouzi2001-hcl.csv","hamer-wu-hcl.csv"}) for(const auto& row:read(inputs/file)) {
            const double m=std::stod(row.at("molality")); if(std::string(file)==training && m<=7) continue;
            auto x=solvent; x[water]=1/p.molar_mass_kg_per_mol[water]; x[H]=x[Cl]=m; const double sum=std::accumulate(x.begin(),x.end(),0.0); for(auto& v:x) v/=sum;
            const auto a=solvent_activity(model,temperature_K,pressure_Pa,x,ref);
            out << file << ',' << m << ',' << a.solvent_activity.value_or(NAN) << ',' << (a.osmotic_coefficient && std::holds_alternative<double>(*a.osmotic_coefficient) ? std::get<double>(*a.osmotic_coefficient) : NAN) << ',' << quote(a.cause) << '\n';
        }
        return 0;
    }
    if (argc == 4 && std::string(argv[3]).rfind("source-",0)==0) {
        const std::string mode(argv[3]), name=mode.substr(7,mode.rfind(':')-7);
        const double pressure=std::stod(mode.substr(mode.rfind(':')+1));
        const std::size_t background=name=="jiang1996-hcl-li" ? lithium : name=="khoo1977-hcl-ca" ? 4 : magnesium;
        if (!supplied_endpoint && !packet_readback) p.k_ij[H*count+water] = p.k_ij[water*count+H] = -1.0706252324767123;
        std::optional<r::EvaluatedCoordinates> source_chart;
        if (supplied_endpoint) {
            r::CoordinateSet coordinates; coordinates.coordinates=coordinate_declarations;
            source_chart=r::evaluate_coordinates(p,coordinates,{std::stod(argv[4]),std::stod(argv[5])},{});
        }
        HelmholtzModel model(p,source_chart ? source_chart->active : std::vector<ActiveParameter>{}); const auto ref=infinite_dilution_reference(model,temperature_K,pressure,solvent,{H,Cl},{1,1});
        std::ofstream out(results / "source.csv");
        out << std::setprecision(17) << "ionic_strength,y_background,m_HCl,m_background,gamma_prediction,log_gamma,phi,density,pressure_slope,reference_density,density_residual,reference_residual,density_stop,reference_stop,reference_status,reference_cause,identity_error,passes,locator,gamma_observed,status,cause,density_convergence,reference_convergence,reference_log_H,reference_log_Cl\n";
        std::ofstream derivatives(results/"source-columns.csv");
        derivatives << std::setprecision(17) << "row,coordinate,prediction_partial,log_partial,status,cause\n";
        bool complete=true; std::size_t source_row=0;
        for(const auto& row:read(inputs / (name+".csv"))) {
            const double I=std::stod(row.at("ionic_strength")), acid=background==lithium ? std::stod(row.at("m_HCl")) : I*(1-std::stod(row.at("y_background"))), salt=background==lithium ? std::stod(row.at("m_LiCl")) : I*std::stod(row.at("y_background"))/3;
            std::vector<double> m(count,0); m[H]=acid; m[background]=salt; m[Cl]=acid+p.charge_number[background]*salt;
            double charge=0, strength=0;
            for (std::size_t i=0; i<count; ++i) { charge += m[i]*p.charge_number[i]; strength += .5*m[i]*p.charge_number[i]*p.charge_number[i]; }
            if (std::abs(charge)>1e-10 || std::abs(strength-I)>1e-10) throw std::runtime_error("source composition balance");
            r::Observation observation; observation.kind=r::ObservationKind::MeanIonicActivityCoefficient;
            observation.temperature=temperature_K; observation.pressure=pressure; observation.composition=solvent;
            observation.branch=DensityBranch::Liquid; observation.ion_indices={H,Cl}; observation.stoichiometry={1,1};
            observation.molality=acid; observation.finite_solute_molalities=m;
            observation.observed=std::pow(10.0, std::stod(row.at(background==lithium ? "log10_gamma" : "negative_log10_gamma"))*(background==lithium ? 1 : -1));
            const auto prediction=r::predict(model,observation);
            const auto state=prediction.activity ? std::optional<State>(prediction.activity->state) : prediction.state; const auto& rs=ref.solvent_state;
            const auto dc=state ? state->density_convergence : std::nullopt;
            const auto rc=rs ? rs->density_convergence : std::nullopt;
            const double gamma=prediction.value.value_or(NAN), log_gamma=gamma>0 ? std::log(gamma) : NAN;
            double expected=NAN;
            if (state && state->log_fugacity_coefficient && state->log_fugacity_coefficient->size()>std::max(H,Cl)) {
                const auto& logs=*state->log_fugacity_coefficient;
                if (logs[H] && logs[Cl]) expected=.5*(*logs[H]-ref.log_fugacity[0]+*logs[Cl]-ref.log_fugacity[1])-std::log(1+p.molar_mass_kg_per_mol[water]*std::accumulate(m.begin(),m.end(),0.0));
            }
            const double phi=prediction.activity && prediction.activity->osmotic_coefficient && std::holds_alternative<double>(*prediction.activity->osmotic_coefficient) ? std::get<double>(*prediction.activity->osmotic_coefficient) : NAN;
            const auto converged=[](const auto& c) { return c && c->criterion && std::isfinite(c->residual) && std::abs(c->residual)<=c->criterion->value && c->stop==core::Stop::converged; };
            const double error=std::abs(log_gamma-expected);
            const bool ok=prediction.status==r::ObservationStatus::Available && std::isfinite(gamma) && gamma>0 && state && state->molar_density>0 && state->pressure_density_derivative>0 && ref.solvent_cause.empty() && rs && rs->molar_density>0 && rs->pressure_density_derivative>0 && converged(dc) && converged(rc) && error<=1e-10;
            out << I << ',' << (I-acid)/I << ',' << acid << ',' << salt << ',' << gamma << ',' << log_gamma << ',' << phi << ',' << (state ? state->molar_density : NAN) << ',' << (state ? state->pressure_density_derivative : NAN) << ',' << (rs ? rs->molar_density : NAN) << ',' << (dc ? dc->residual : NAN) << ',' << (rc ? rc->residual : NAN) << ',' << quote(dc ? core::stop_name(dc->stop) : "absent") << ',' << quote(rc ? core::stop_name(rc->stop) : "absent") << ',' << int(ref.solvent_status) << ',' << quote(ref.solvent_cause) << ',' << error << ',' << ok;
            out << ',' << quote(row.at("locator")) << ',' << observation.observed << ',' << r::observation_status_name(prediction.status) << ',' << quote(prediction.message) << ',' << quote(dc ? core::describe(*dc) : "absent") << ',' << quote(rc ? core::describe(*rc) : "absent") << ',' << ref.log_fugacity[0] << ',' << ref.log_fugacity[1] << '\n';
            if (source_chart) {
                const auto gradient=r::differentiate(model,observation,prediction);
                for(std::size_t a=0;a<2;++a) {
                    double partial=0;
                    if(gradient.value.size()==source_chart->active.size()) for(std::size_t b=0;b<gradient.value.size();++b) partial+=gradient.value[b]*source_chart->map.at(2*b+a);
                    else partial=NAN;
                    derivatives << source_row << ',' << a << ',' << partial << ',' << partial/gamma << ',' << r::observation_status_name(gradient.status) << ',' << quote(gradient.message) << '\n';
                }
            }
            ++source_row; complete=complete && ok;
        }
        return complete ? 0 : 1;
    }
    std::ofstream status(results / "fit-status.csv"), values(results / "fit-values.csv"), failures(results / "trial-failures.csv");
    status << "stage,start,status,usable,cost,iterations,residual_evaluations,jacobian_evaluations,covariance,rank,message\n";
    values << std::setprecision(17) << "stage,start,field,index,value\n";
    failures << "stage,start,trial,observation,status,jacobian,optimizer,message\n";
    auto dump = [&](int q, std::size_t start, const std::string& field, const auto& vector) {
        for (std::size_t i = 0; i < vector.size(); ++i) values << q << ',' << start << ',' << field << ',' << i << ',' << vector[i] << '\n';
    };
    r::FitRequest request; request.parameters = p; request.singular_vectors = true;
    request.controls.maximum_iterations = maximum_iterations; request.controls.maximum_elapsed_time_seconds = maximum_seconds;
    struct Row { std::string source, cell, observable; double mNa, mSr; std::string salt{}, bucket{}; };
    std::vector<Row> metadata;
    if (sulfate_training) {
        std::ofstream admission(results / "mixture-admission.csv");
        admission << std::setprecision(17) << "observation,source,cell,salt,bucket,molality,m_LiCl,I,D,aw,observed,x_water,x_Na,x_Mg,x_Li,x_Cl,x_SO4,weight,passes\n";
        std::map<std::string, int> groups; std::map<std::string, bool> seen;
        for (const auto& row : read(inputs / training)) {
            const auto& salt = row.at("salt"); const auto& bucket = row.at("bucket");
            const double m = std::stod(row.at("molality")), a = std::stod(row.at("m_LiCl")), aw = std::stod(row.at("aw"));
            const bool mixed = salt == "LiCl-MgSO4", na = salt == "Na2SO4";
            const double D = mixed ? 2*a+2*m : (na ? 3*m : 2*m), I = (na ? 3*m : 4*m)+a;
            const double observed = mixed ? -std::log(aw)/(p.molar_mass_kg_per_mol[water]*D) : std::stod(row.at("phi"));
            if (!(m > 0 && aw > 0 && aw < 1 && std::isfinite(observed) && observed > 0) || (na ? m > 1 : m > 2) || (mixed ? !(a > 0 && I <= 8.3) : a != 0) || !(na || mixed || salt == "MgSO4")) throw std::runtime_error("sulfate conditions");
            const std::string cell = row.at("source")+"/"+row.at("cell");
            if (seen[cell]) throw std::runtime_error("duplicate sulfate cell");
            seen[cell] = true; ++groups[bucket];
            r::Observation o; o.temperature = temperature_K; o.pressure = pressure_Pa; o.branch = DensityBranch::Liquid;
            o.composition.assign(count, 0); o.composition[water] = 1/p.molar_mass_kg_per_mol[water];
            o.composition[sodium] = na ? 2*m : 0; o.composition[magnesium] = na ? 0 : m; o.composition[lithium] = a; o.composition[Cl] = a; o.composition[H] = m;
            double charge = 0; for (std::size_t i = 0; i < count; ++i) charge += o.composition[i]*p.charge_number[i];
            const double total = std::accumulate(o.composition.begin(), o.composition.end(), 0.0); for (auto& x : o.composition) x /= total;
            if (std::abs(charge) > 1e-10 || std::abs(std::accumulate(o.composition.begin(),o.composition.end(),0.0)-1) > 1e-10) throw std::runtime_error("sulfate balance");
            const auto composition = o.composition;
            o.kind = mixed ? r::ObservationKind::MixtureOsmoticCoefficient : r::ObservationKind::OsmoticCoefficient;
            o.ion_indices = {na ? std::size_t(sodium) : std::size_t(magnesium), H}; o.stoichiometry = {na ? 2 : 1, 1};
            if (!mixed) o.composition = solvent;
            o.form = r::ResidualForm::Difference; o.scale = .01; o.observed = observed; o.molality = m;
            const double weight = 1.0/((mixed ? 2 : 8)*std::stoi(row.at("bucket_states")));
            request.observations.push_back(o); request.weights.push_back(weight); metadata.push_back({mixed ? "Zhang2016" : "binary", cell, mixed ? "mixture_phi" : "phi", a, m, salt, bucket});
            admission << metadata.size()-1 << ',' << row.at("source") << ',' << row.at("cell") << ',' << salt << ',' << bucket << ',' << m << ',' << a << ',' << I << ',' << D << ',' << aw << ',' << observed;
            for (auto index : {water,sodium,magnesium,lithium,chloride,ion}) admission << ',' << composition[index];
            admission << ',' << weight << ",1\n";
        }
        if (groups != std::map<std::string,int>{{"Na-low",5},{"Na-high",4},{"Mg-low",9},{"Mg-high",5},{"mixture",21}} || std::abs(std::accumulate(request.weights.begin(),request.weights.end(),0.0)-1) > 1e-14) throw std::runtime_error("sulfate source groups");
    } else {
    const auto binary = read(inputs / training);
    for (const auto& row : binary) {
        const double m = std::stod(row.at("molality")); if (m > cutoff) continue;
        for (bool gamma : {false, true}) {
            r::Observation o; o.temperature = temperature_K; o.pressure = pressure_Pa; o.composition = solvent;
            o.branch = DensityBranch::Liquid;
            o.ion_indices = {H, Cl}; o.stoichiometry = {nu_cation, nu_anion}; o.molality = m;
            o.kind = gamma ? r::ObservationKind::MeanIonicActivityCoefficient : r::ObservationKind::OsmoticCoefficient;
            o.form = gamma ? r::ResidualForm::LogRatio : r::ResidualForm::Difference;
            o.observed = std::stod(row.at(gamma ? "gamma" : "phi")); o.scale = gamma ? 1.0 : 0.01;
            request.observations.push_back(o);
            metadata.push_back({"binary", std::to_string(metadata.size()/2+1), gamma ? "gamma" : "phi", 0, m});
            request.weights.push_back(1.0 / ((mixture_states ? 8.0 : 4.0) * (m <= split ? low_states : training_states - low_states) * (gamma ? std::pow(std::log(1.02), 2) : 1.0)));
        }
    }
    if (request.observations.size() != 2 * training_states) throw std::runtime_error("training count");
    if (mixture_states) {
        const auto mixtures = read(inputs / mixture_training);
        std::ofstream admission;
        if (!verify_lengths || (argc == 4 && std::string(argv[3]) == "prepare")) admission.open(results / "mixture-admission.csv");
        admission << std::setprecision(17) << "observation,cell,y,m_total,m_NaCl,m_" << ion_label << "Cl2,I,D,observed,x_water,x_Na,x_" << ion_label << ",x_Cl,weight,passes\n";
        std::map<std::string, bool> seen; std::vector<int> counts(3, 0);
        const std::vector<double> fractions{0.17066, 0.47366, 0.82682};
        Rows knots, references, original, ratios; std::ofstream audit, deviation;
        if (replacement_reference) {
            knots = read(inputs / "hamer-wu1972-nacl-knots.csv"); references = read(inputs / "robinson-bower1965-references.csv"); original = read(inputs / "robinson-bower1965-table1.csv"); ratios = read(inputs / "robinson-bower1965-ratios.csv");
            if (knots.size() != 13 || references.size() != 15 || original.size() != 60 || ratios.size() != 60 || p.molar_mass_kg_per_mol[water] != .0180153) throw std::runtime_error("replacement input inventory");
            audit.open(results / "reference-reproduction.csv"); audit << std::setprecision(17) << "observation,Table1_row,run,m_ref,phi_ref_historical,phi_ref_native,phi_mix_historical,phi_mix_native,target_difference,ratio_interval_passes,conversion\n";
            deviation.open(results / "reference-deviation.csv"); deviation << std::setprecision(17) << "observation,aw,aw_printed,delta_aw,delta_phi_ref_historical,delta_phi_ref_native,delta_phi_mix_native,table_rounding_bound\n";
        }
        for (const auto& row : mixtures) {
            const std::string cell = replacement_reference ? "T1-row" + row.at("Table1_row") + "-run" + row.at("run") : "T" + row.at("table") + "-S" + row.at("series") + "-E" + row.at("experiment") + "-y" + row.at("y");
            const auto require = [&](bool ok, const char* why) { if (!ok) throw std::runtime_error(cell + ": " + why); };
            const auto number = [&](const char* key) { const double v = std::stod(row.at(key)); require(std::isfinite(v), key); return v; };
            double y, m, ref, phi, I, mNa, mSr, D, observed;
            if (replacement_reference) {
                const auto equal = [](double a, double b) { return std::isfinite(a) && std::isfinite(b) && std::abs(a-b) <= 1e-10*std::max(1.0, std::abs(b)); };
                const auto index = std::stoul(row.at("Table1_row")); const auto& source = original.at(index-1);
                ref = number("m_ref"); mNa = number("m_NaCl"); mSr = number("m_BaCl2"); m = mNa+mSr; I = mNa+3*mSr; D = 2*mNa+3*mSr; y = mNa/I;
                require(!seen[cell] && mNa > 0 && mSr > 0 && number("observation_index") == metadata.size() && equal(number("I"), I) && equal(number("D"), D), "mixture inventory"); seen[cell] = true;
                require(ref == std::stod(source.at("NaCl_reference_mol_per_kg_water")) && mNa == std::stod(source.at("NaCl_mol_per_kg_water")) && mSr == std::stod(source.at("BaCl2_mol_per_kg_water")), "source composition");
                const auto run = std::find_if(references.begin(), references.end(), [&](const auto& r) { return r.at("run") == row.at("run"); }); require(run != references.end() && ref == std::stod(run->at("m_ref_mol_kg_water")) && source.at("set") == run->at("printed_set"), "reference link");
                phi = NAN;
                for (std::size_t a = 1; a < knots.size(); ++a) { const double l = std::stod(knots[a-1].at("m_NaCl_mol_kg_water")), h = std::stod(knots[a].at("m_NaCl_mol_kg_water")); if (ref >= l && ref <= h) { const double low = std::stod(knots[a-1].at("phi_HW_Table16")); phi = low+(ref-l)/(h-l)*(std::stod(knots[a].at("phi_HW_Table16"))-low); break; } }
                const double conversion = 1/(55.51*p.molar_mass_kg_per_mol[water]), ratio = 2*ref/D, historical = ratio*phi; observed = conversion*historical;
                require(equal(phi, std::stod(run->at("Phi_ref_historical_basis"))) && equal(conversion*phi, std::stod(run->at("Phi_ref_native"))) && equal(historical, number("Phi_mix_historical_basis")) && equal(observed, number("Phi_mix_native")), "817-REF or 817-MAP");
                const double h = .00005, printed = std::stod(ratios.at(index-1).at("R_obs_printed"));
                require(printed+h >= (ref-h)/(mNa+h+1.5*(mSr+h)) && printed-h <= (ref+h)/(mNa-h+1.5*(mSr-h)), "printed ratio interval");
                audit << metadata.size() << ',' << index << ',' << row.at("run") << ',' << ref << ',' << phi << ',' << conversion*phi << ',' << historical << ',' << observed << ',' << observed-number("Phi_mix_native") << ",1," << conversion << '\n';
                const double aw = std::exp(-2*ref*phi/55.51), shift = phi+55.51*std::log(.9453)/(2*ref);
                if (row.at("run") == "5") deviation << metadata.size() << ',' << aw << ",0.9453," << aw-.9453 << ',' << shift << ',' << conversion*shift << ',' << conversion*ratio*shift << ',' << conversion*ratio*.0005 << '\n';
            } else {
            y = number("y"); m = number("m_total"); ref = number("m_reference"); phi = number("phi_reference");
            require(!seen[cell] && row.at("table") == "2" && m > 0 && ref > 0 && phi > 0 && y > 0 && y < 1, "cell inventory or positive inputs"); seen[cell] = true;
            require(number("m_total_plus_minus") >= 0 && number("m_reference_plus_minus") >= 0 && number("time_days") > 0, "printed entries");
            const auto group = std::find(fractions.begin(), fractions.end(), y); require(group != fractions.end(), "composition group");
            const auto g = group - fractions.begin(); require(std::size_t(g) == (metadata.size()-2*training_states < 17 ? 0 : metadata.size()-2*training_states < 35 ? 1 : 2), "source order"); ++counts[g];
            if (row.at("series") == "2") {
                const std::map<std::string, double> corrected{{"1", 1.62565}, {"2", 1.30315}, {"3", 1.16184}};
                require(corrected.count(row.at("experiment")) && ref == corrected.at(row.at("experiment")) && row.at("footnotes").find('c') != std::string::npos, "corrected Series 2 reference");
            }
            I = 3*m/(1+2*y); mNa = y*I; mSr = (1-y)*I/3; D = 2*mNa+3*mSr; observed = 2*ref*phi/D;
            }
            const auto equal = [](double a, double b) { return std::isfinite(a) && std::abs(a-b) <= 1e-10*std::max(1.0, std::abs(b)); };
            require(equal(mNa+mSr, m) && equal(mNa+3*mSr, I) && equal(D*observed, 2*ref*phi/(replacement_reference ? 55.51*p.molar_mass_kg_per_mol[water] : 1)), "molality or osmotic balance");
            r::Observation o; o.temperature = temperature_K; o.pressure = pressure_Pa; o.branch = DensityBranch::Liquid;
            o.composition.assign(count, 0); o.composition[water] = 1/p.molar_mass_kg_per_mol[water];
            o.composition[sodium] = mNa; o.composition[H] = mSr; o.composition[Cl] = mNa+2*mSr;
            double charge = 0; for (std::size_t i = 0; i < count; ++i) charge += o.composition[i]*p.charge_number[i];
            require(equal(charge, 0), "electroneutral inventory");
            const double amount = std::accumulate(o.composition.begin(), o.composition.end(), 0.0);
            for (auto& x : o.composition) x /= amount;
            require(equal(std::accumulate(o.composition.begin(), o.composition.end(), 0.0), 1), "normalized composition");
            o.kind = r::ObservationKind::MixtureOsmoticCoefficient; o.form = r::ResidualForm::Difference; o.scale = 0.01; o.observed = observed; o.molality = m;
            request.observations.push_back(o); request.weights.push_back(1.0/(2.0*mixture_states)); metadata.push_back({replacement_reference ? "Robinson-Bower" : "Clegg", cell, "mixture_phi", mNa, mSr});
            admission << metadata.size()-1 << ',' << cell << ',' << y << ',' << m << ',' << mNa << ',' << mSr << ',' << I << ',' << D << ',' << observed << ',' << o.composition[water] << ',' << o.composition[sodium] << ',' << o.composition[H] << ',' << o.composition[Cl] << ',' << request.weights.back() << ",1\n";
        }
        if (mixtures.size() != std::size_t(mixture_states) || (!replacement_reference && counts != std::vector<int>{17, 18, 14})) throw std::runtime_error("mixture occupied-cell counts");
    }
    }
    if (hcl_mixture_training) {
        const auto roy = read(inputs / "roy1980-hcl-mg.csv");
        std::map<double, int> counts;
        for (const auto& row : roy) if (std::stod(row.at("y_background")) > 0) ++counts[std::stod(row.at("ionic_strength"))];
        if (roy.size() != 58 || counts.size() != 11) throw std::runtime_error("Roy inventory");
        for (const auto& row : roy) {
            const double I = std::stod(row.at("ionic_strength")), y = std::stod(row.at("y_background"));
            const double acid = I*(1-y), salt = I*y/3;
            r::Observation o; o.temperature = temperature_K; o.pressure = pressure_Pa; o.composition = solvent;
            o.branch = DensityBranch::Liquid; o.ion_indices = {H,Cl}; o.stoichiometry = {1,1}; o.molality = acid;
            o.kind = r::ObservationKind::MeanIonicActivityCoefficient; o.form = r::ResidualForm::LogRatio;
            o.observed = std::pow(10., -std::stod(row.at("negative_log10_gamma")));
            o.finite_solute_molalities = std::vector<double>(count, 0);
            auto& m = *o.finite_solute_molalities; m[H] = acid; m[magnesium] = salt; m[Cl] = acid+2*salt;
            request.observations.push_back(o);
            request.weights.push_back(y > 0 ? 1/(11.*counts.at(I)*std::pow(std::log(1.05),2)) : 0);
            metadata.push_back({"Roy", row.at("locator"), "gamma", salt, acid});
        }
    }
    if (profile_enabled && argc == 4 && (std::string(argv[3]) == "profile" || std::string(argv[3]) == "profile-replay")) {
        for (auto& weight : request.weights) weight *= .8;
        for (const auto& row : read(inputs / "jiang1996-hcl-li.csv")) if (std::stod(row.at("ionic_strength")) == 5) {
            const double acid=std::stod(row.at("m_HCl")), salt=std::stod(row.at("m_LiCl"));
            if (std::abs(acid+salt-5)>1e-10) throw std::runtime_error("Jiang I=5 inventory");
            auto o=request.observations[1]; o.molality=acid;
            o.observed=std::pow(10.0,std::stod(row.at("log10_gamma")));
            o.finite_solute_molalities=std::vector<double>(count,0);
            (*o.finite_solute_molalities)[H]=acid; (*o.finite_solute_molalities)[lithium]=salt; (*o.finite_solute_molalities)[Cl]=acid+salt;
            request.observations.push_back(o); request.weights.push_back(1/(50*std::pow(std::log(1.05),2)));
            metadata.push_back({"Jiang-I5",row.at("locator"),"gamma",salt,acid});
        }
        if (request.observations.size()!=70) throw std::runtime_error("profile observation count");
        std::ofstream rows(results/"training.csv"), states(results/"states.csv"), attempts(results/"attempts.csv"), columns(results/"columns.csv");
        rows << std::setprecision(17) << "point,start,c,u,row,source,locator,molality,background,observable,observed,prediction,residual,weight,weighted_residual,status,cause,verified\n";
        states << std::setprecision(17) << "point,start,row,role,density,pressure,slope,branch,convergence,cause,log_H,log_Cl,identity_error\n";
        attempts << std::setprecision(17) << "point,start,c,u,converged,verified,cost,B,V,mean_I5,worst_I5,sequence_wall_seconds\n";
        std::ofstream column_summary(results/"column-summary.csv"); column_summary << std::setprecision(17) << "point,start,norm_u,norm_c,inverse_normal_correlation\n";
        columns << std::setprecision(17) << "point,start,row,coordinate,prediction_partial,log_partial,residual_partial,objective_partial,status,cause\n";
        struct Endpoint { double c,u,J,B,V,mean,worst; bool verified,converged; std::size_t point,start; };
        std::vector<Endpoint> endpoints;
        const auto trace=[&](std::size_t point,std::size_t start,std::size_t row,const std::string& role,const auto& state,const std::string& cause) {
            states << point << ',' << start << ',' << row << ',' << role << ',' << (state ? state->molar_density : NAN) << ',' << (state ? state->pressure : NAN) << ',' << (state ? state->pressure_density_derivative : NAN) << ',' << (state ? int(state->density_branch) : -1) << ',' << quote(state && state->density_convergence ? core::describe(*state->density_convergence) : "absent") << ',' << quote(cause) << ",nan,nan,nan\n";
            return state && state->molar_density>0 && std::isfinite(state->molar_density) && state->pressure_density_derivative>0 && state->density_convergence && state->density_convergence->stop==core::Stop::converged && state->density_convergence->criterion && std::isfinite(state->density_convergence->residual) && std::abs(state->density_convergence->residual)<=state->density_convergence->criterion->value;
        };
        const auto evaluate=[&](std::size_t point,std::size_t start,double c,double u,bool converged,double J,double wall) {
            r::CoordinateSet coordinates; coordinates.coordinates=coordinate_declarations;
            auto chart=r::evaluate_coordinates(p,coordinates,{u,c},{});
            if (packet_readback) chart.parameters=p;
            HelmholtzModel model(chart.parameters,chart.active); r::ReferenceCache references;
            Endpoint e{c,u,J,0,0,0,0,true,converged,point,start};
            double gram[3]{};
            double errors[2][3]{}, maxima[2][3]{}; std::size_t counts[2][3]{};
            for (std::size_t i=0; i<request.observations.size(); ++i) {
                const auto& o=request.observations[i]; const auto prediction=r::predict(model,o,std::nullopt,std::nullopt,&references);
                const auto state=prediction.activity ? std::optional<State>(prediction.activity->state) : prediction.state;
                const double y=prediction.value.value_or(NAN); const auto gradient=r::differentiate(model,o,prediction);
                bool verified=prediction.status==r::ObservationStatus::Available && std::isfinite(y) && trace(point,start,i,"state",state,prediction.message);
                double identity_error=NAN;
                if (prediction.reference) {
                    const auto& ref=*prediction.reference;
                    verified=trace(point,start,i,"reference",ref.solvent_state,ref.solvent_cause) && ref.solvent_cause.empty() && verified;
                    if (prediction.activity && state && state->log_fugacity_coefficient) {
                        const auto& logs=*state->log_fugacity_coefficient;
                        const double solutes=o.finite_solute_molalities ? std::accumulate(o.finite_solute_molalities->begin(),o.finite_solute_molalities->end(),0.0) : 2*o.molality;
                        if (logs.size()>std::max(H,Cl) && logs[H] && logs[Cl]) identity_error=std::abs(prediction.activity->mean_ionic_log_molality_scale-(.5*(*logs[H]-ref.log_fugacity[0]+*logs[Cl]-ref.log_fugacity[1])-std::log(1+p.molar_mass_kg_per_mol[water]*solutes)));
                        states << point << ',' << start << ',' << i << ",ionic-reference,nan,nan,nan,-1," << quote(ref.solvent_convergence ? core::describe(*ref.solvent_convergence) : "absent") << ',' << quote(ref.solvent_cause) << ',' << ref.log_fugacity[0] << ',' << ref.log_fugacity[1] << ',' << identity_error << '\n';
                    }
                    verified=verified && std::isfinite(identity_error) && identity_error<=1e-10;
                } else verified=false;
                const bool gamma=o.kind==r::ObservationKind::MeanIonicActivityCoefficient;
                verified=verified && (!gamma || y>0) && gradient.status==r::ObservationStatus::Available;
                double A[2]{NAN,NAN};
                for (std::size_t a=0; a<2; ++a) {
                    double partial=0;
                    if (gradient.value.size()==chart.active.size()) for (std::size_t b=0; b<gradient.value.size(); ++b) partial+=gradient.value[b]*chart.map.at(2*b+a);
                    else partial=NAN;
                    const double residual_partial=gamma ? partial/y : partial/o.scale;
                    verified=verified && std::isfinite(partial); A[a]=std::sqrt(request.weights[i])*residual_partial;
                    columns << point << ',' << start << ',' << i << ',' << a << ',' << partial << ',' << (gamma ? partial/y : NAN) << ',' << residual_partial << ',' << std::sqrt(request.weights[i])*residual_partial << ',' << r::observation_status_name(gradient.status) << ',' << quote(gradient.message) << '\n';
                }
                gram[0]+=A[0]*A[0]; gram[1]+=A[0]*A[1]; gram[2]+=A[1]*A[1];
                const auto& m=metadata[i]; const double residual=prediction.residual.value_or(NAN);
                rows << point << ',' << start << ',' << c << ',' << u << ',' << i << ',' << m.source << ',' << quote(m.cell) << ',' << o.molality << ',' << m.mNa << ',' << m.observable << ',' << o.observed << ',' << y << ',' << residual << ',' << request.weights[i] << ',' << std::sqrt(request.weights[i])*residual << ',' << r::observation_status_name(prediction.status) << ',' << quote(prediction.message) << ',' << verified << '\n';
                e.verified=e.verified && verified;
                const double error=gamma ? 100*(y/o.observed-1) : y-o.observed;
                if (i<60) {
                    const auto bucket=o.molality>3 ? 1 : 0, quantity=gamma ? 1 : 0;
                    errors[bucket][quantity]+=gamma ? std::abs(error) : error*error;
                    maxima[bucket][quantity]=std::max(maxima[bucket][quantity],std::abs(error)); ++counts[bucket][quantity];
                } else if (m.mNa>0) { e.mean+=std::abs(error)/9; e.worst=std::max(e.worst,std::abs(error)); }
            }
            const auto ref=solvent_reference(model,temperature_K,pressure_Pa,solvent);
            std::size_t water_row=70;
            for (const auto& row:read(inputs/"messnaoui2020-hcl.csv")) {
                const double m=std::stod(row.at("molality")); auto x=solvent;
                x[water]=1/p.molar_mass_kg_per_mol[water]; x[H]=x[Cl]=m;
                const double total=std::accumulate(x.begin(),x.end(),0.0); for(auto& v:x) v/=total;
                const auto value=solvent_activity(model,temperature_K,pressure_Pa,x,ref);
                const double y=value.solvent_activity.value_or(NAN), observed=std::stod(row.at("aw_adjusted")), error=y-observed;
                const bool verified=std::isfinite(y) && value.cause.empty() && trace(point,start,water_row,"water",std::optional<State>(value.state),value.cause) && trace(point,start,water_row,"water-reference",ref.pure_solvent_state,ref.solvent_cause);
                rows << point << ',' << start << ',' << c << ',' << u << ',' << water_row << ",Messnaoui," << quote("Messnaoui Table 2 m="+row.at("molality")) << ',' << m << ",0,aw," << observed << ',' << y << ',' << error << ",0,0," << int(value.status) << ',' << quote(value.cause) << ',' << verified << '\n';
                for (std::size_t a=0; a<2; ++a) {
                    double partial=0; bool available=value.solvent_activity_parameter_derivative.size()==chart.active.size();
                    for (std::size_t b=0; b<value.solvent_activity_parameter_derivative.size(); ++b) {
                        available=available && value.solvent_activity_parameter_derivative[b].has_value();
                        partial+=value.solvent_activity_parameter_derivative[b].value_or(NAN)*chart.map.at(2*b+a);
                    }
                    columns << point << ',' << start << ',' << water_row << ',' << a << ',' << (available ? partial : NAN) << ",nan,nan,0," << int(value.status) << ',' << quote(value.cause) << '\n';
                }
                const auto bucket=m>3 ? 1 : 0; errors[bucket][2]+=error*error;
                maxima[bucket][2]=std::max(maxima[bucket][2],std::abs(error)); ++counts[bucket][2]; e.verified=e.verified && verified; ++water_row;
            }
            const auto& mean_limits=profile_mean_limits; const auto& max_limits=profile_max_limits;
            for (int b=0; b<2; ++b) for (int q=0; q<3; ++q) {
                const double mean=errors[b][q]/counts[b][q];
                e.B=std::max({e.B,(q==1 ? mean : std::sqrt(mean))/mean_limits[q],maxima[b][q]/max_limits[q]});
            }
            column_summary << point << ',' << start << ',' << std::sqrt(gram[0]) << ',' << std::sqrt(gram[2]) << ',' << -gram[1]/std::sqrt(gram[0]*gram[2]) << '\n'; column_summary.flush();
            e.V=std::max(e.mean/profile_target_mean,e.worst/profile_target_max);
            attempts << point << ',' << start << ',' << c << ',' << u << ',' << converged << ',' << e.verified << ',' << J << ',' << e.B << ',' << e.V << ',' << e.mean << ',' << e.worst << ',' << wall << '\n'; attempts.flush(); rows.flush(); states.flush(); columns.flush();
            return e;
        };
        if (supplied_endpoint || packet_readback) {
            const auto replay=evaluate(0,0,packet_readback ? p.k_ij[H*count+Cl] : std::stod(argv[5]),packet_readback ? p.k_ij[H*count+water] : std::stod(argv[4]),true,0,0);
            return replay.verified ? 0 : 1;
        }
        const auto baseline=evaluate(0,0,0,profile_anchor,true,0,0);
        if (!baseline.verified) throw std::runtime_error("baseline numerical admission failed");
        std::vector<double> previous;
        for(const auto& row:read(leaf/"results/issue-798-r2/fit-values.csv")) if(row.at("stage")=="1" && row.at("start")=="3" && row.at("field")=="predictions") previous.push_back(std::stod(row.at("value")));
        const auto baseline_rows=read(results/"training.csv");
        if(previous.size()!=60 || baseline_rows.size()!=84) throw std::runtime_error("baseline replay inventory");
        for(std::size_t i=0;i<60;++i) {
            const double now=std::stod(baseline_rows[i].at("prediction")), old=previous[i], m=std::stod(baseline_rows[i].at("molality"));
            const double tolerance=std::max(1e-10*std::abs(old),i%2==0 ? 1e-12/(2*p.molar_mass_kg_per_mol[water]*m) : 1e-12);
            if(std::abs(now-old)>tolerance) throw std::runtime_error("binary baseline replay mismatch at "+std::to_string(i));
        }
        const auto water_baseline=read(leaf/"results/binary-water-qualification/native/water.csv");
        if(water_baseline.size()!=14) throw std::runtime_error("water baseline inventory");
        for(std::size_t i=0;i<14;++i) if(std::abs(std::stod(baseline_rows[i+70].at("prediction"))-std::stod(water_baseline[i].at("aw_prediction")))>1e-10) throw std::runtime_error("water baseline replay mismatch");
        const auto jiang_baseline=read(inputs/"issue-901-jiang-baseline.csv");
        if(jiang_baseline.size()!=10) throw std::runtime_error("Jiang baseline inventory");
        for(const auto& row:jiang_baseline) {
            const double y=std::stod(row.at("y")), old=std::stod(row.at("baseline"));
            const auto found=std::find_if(baseline_rows.begin()+60,baseline_rows.begin()+70,[&](const auto& r){return std::abs(std::stod(r.at("background"))/5-y)<1e-12;});
            if(found==baseline_rows.begin()+70 || std::abs(std::stod(found->at("prediction"))-old)>std::max(1e-12,1e-10*std::abs(old))) throw std::runtime_error("Jiang baseline replay mismatch");
        }
        std::ofstream policies(results/"solver-policy.csv"); policies << std::setprecision(17) << "name,value,unit\n";
        for(std::size_t i=0;i<policy::SOLVER_POLICIES.size();++i) policies << policy::SOLVER_POLICIES[i]->name << ',' << request.controls.solver_options.value(i) << ',' << policy::SOLVER_UNITS[i] << '\n';
        std::vector<double> grid=profile_coarse;
        std::vector<Endpoint> ordinary;
        const auto best_point=[&]() {
            std::vector<Endpoint> complete;
            for(const auto& e:ordinary) {
                bool ok=e.verified && e.converged;
                std::size_t count_starts=0;
                for(const auto& a:endpoints) if(a.point==e.point) {
                    ++count_starts; ok=ok && a.verified && a.converged && a.u>coordinate_declarations[0].lower+1e-6 && a.u<coordinate_declarations[0].upper-1e-6;
                    if(a.J<=e.J+std::max(1e-10,1e-4*std::abs(e.J))) ok=ok && std::abs(a.u-e.u)<=1e-6;
                }
                if(ok && count_starts==3) complete.push_back(e);
            }
            if(complete.empty()) throw std::runtime_error("no numerically complete coarse point");
            const bool eligible=std::any_of(complete.begin(),complete.end(),[](const auto& e){return e.B<=1;});
            std::stable_sort(complete.begin(),complete.end(),[&](const auto& a,const auto& b) {
                return std::make_tuple(eligible && a.B>1,eligible ? a.V : a.B,a.V,a.J,std::abs(a.c),a.point)<std::make_tuple(eligible && b.B>1,eligible ? b.V : b.B,b.V,b.J,std::abs(b.c),b.point);
            });
            return complete.front();
        };
        for(std::size_t point=0;point<grid.size();++point) {
            const double c=grid[point];
            std::cout << "profile point " << point+1 << " c=" << c << std::endl;
            request.parameters=p; request.parameters.k_ij[H*count+Cl]=request.parameters.k_ij[Cl*count+H]=c;
            request.coordinates.coordinates={coordinate_declarations[0]};
            const std::vector<std::vector<double>> allocated{{profile_anchor+profile_compensation*c},{profile_starts[0]},{profile_starts[1]}};
            const auto began=std::chrono::steady_clock::now();
            const auto fits=r::fit_sequence(request,allocated,total_maximum_iterations,total_maximum_seconds);
            const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
            std::optional<Endpoint> best;
            for(std::size_t j=0;j<3;++j) {
                if(j>=fits.size()) {
                    attempts << point+1 << ',' << j+1 << ',' << c << ",nan,0,0,nan,nan,nan,nan,nan," << seconds << '\n';
                    continue;
                }
                dump(point+1,j+1,"start",allocated[j]);
                const auto& f=fits[j]; const auto& covariance=f.covariance;
                status << std::setprecision(17) << point+1 << ',' << j+1 << ',' << f.status << ',' << f.usable << ',' << f.final_cost << ',' << f.iterations << ',' << f.residual_evaluations << ',' << f.jacobian_evaluations << ',' << covariance.status << ',' << covariance.rank << ',' << quote(f.message) << '\n';
                for(const auto& item:std::vector<std::pair<std::string,std::vector<double>>>{{"optimizer",f.optimizer},{"physical",f.physical},{"predictions",f.predictions},{"residuals",f.residuals},{"weighted_residuals",f.weighted_residuals},{"optimizer_jacobian",f.optimizer_jacobian},{"physical_jacobian",f.physical_jacobian},{"singular_values",covariance.singular_values},{"right_vectors",covariance.right_vectors},{"left_vectors",covariance.left_vectors},{"C",covariance.optimizer},{"iteration_seconds",f.iteration_seconds},{"iteration_costs",f.iteration_costs},{"projected_gradient",{covariance.projected_gradient_norm}},{"initial_cost",{f.initial_cost}},{"chart_physical",f.chart.physical},{"chart_map",f.chart.map}}) dump(point+1,j+1,item.first,item.second);
                dump(point+1,j+1,"active_bounds",covariance.active_bounds);
                for(std::size_t i=0;i<f.statuses.size();++i) values << point+1 << ',' << j+1 << ",observation_status," << i << ',' << int(f.statuses[i]) << '\n';
                for(const auto& failure:f.failures) {
                    std::ostringstream optimizer; optimizer << std::setprecision(17); for(double v:failure.optimizer) optimizer << v << ' ';
                    failures << point+1 << ',' << j+1 << ',' << failure.trial << ',' << failure.observation << ',' << r::observation_status_name(failure.status) << ',' << failure.jacobian << ',' << quote(optimizer.str()) << ',' << quote(failure.message+(failure.first_numerical_failure ? "; "+core::describe(*failure.first_numerical_failure) : "")) << '\n';
                }
                std::ofstream detail(results/("fit-"+std::to_string(point+1)+"-"+std::to_string(j+1)+".txt"));
                detail << f.termination << '\n' << (f.first_numerical_failure ? core::describe(*f.first_numerical_failure) : "no first numerical failure") << '\n';
                for(const auto& rank:f.rank_records) detail << rank.basis << ',' << rank.dimension << ',' << rank.rank << ',' << rank.threshold << '\n';
                const bool converged=f.status=="converged" && f.usable && std::isfinite(f.final_cost) && f.optimizer.size()==1 && std::isfinite(f.optimizer[0]);
                if(!converged) {
                    attempts << point+1 << ',' << j+1 << ',' << c << ",nan,0,0,nan,nan,nan,nan,nan," << seconds << '\n';
                    continue;
                }
                auto e=evaluate(point+1,j+1,c,f.optimizer[0],converged,f.final_cost,seconds); endpoints.push_back(e);
                dump(point+1,j+1,"bound_distances",std::vector<double>{e.u-coordinate_declarations[0].lower,coordinate_declarations[0].upper-e.u});
                if(!best || e.J<best->J) best=e;
            }
            if(best) ordinary.push_back(*best);
            values.flush(); failures.flush(); status.flush();
            if(point+1==profile_coarse.size()) {
                const auto center=best_point();
                for(double offset:profile_refinement_offsets) {
                    const double c=std::clamp(center.c+offset,profile_bounds[0],profile_bounds[1]);
                    if(std::none_of(grid.begin(),grid.end(),[&](double old){return std::abs(c-old)<1e-12;})) grid.push_back(c);
                }
            }
        }
        const auto diagnostic=best_point();
        std::vector<Endpoint> successful;
        for(const auto& e:endpoints) if(e.converged && e.verified && e.B<=1 && e.V<=1) successful.push_back(e);
        auto frozen=diagnostic;
        if(!successful.empty()) {
            const double minimum=std::min_element(successful.begin(),successful.end(),[](const auto& a,const auto& b){return a.J<b.J;})->J;
            successful.erase(std::remove_if(successful.begin(),successful.end(),[&](const auto& e){return e.J>minimum+std::max(1e-10,1e-4*std::abs(minimum));}),successful.end());
            std::stable_sort(successful.begin(),successful.end(),[](const auto& a,const auto& b){return std::make_tuple(std::abs(a.c),a.worst,a.point,a.start)<std::make_tuple(std::abs(b.c),b.worst,b.point,b.start);});
            frozen=successful.front();
        }
        std::ofstream frozen_out(results/"frozen.csv");
        frozen_out << std::setprecision(17) << "point,start,c,u,cost,B,V,success\n" << frozen.point << ',' << frozen.start << ',' << frozen.c << ',' << frozen.u << ',' << frozen.J << ',' << frozen.B << ',' << frozen.V << ',' << !successful.empty() << '\n';
        std::ofstream(results/"outcome.txt") << (successful.empty() ? "Profile evaluated; completeness and interval edge require scoring." : "Successful diagnostic candidate; transfer remains unevaluated.") << '\n';
        return 0;
    }
    if (request.observations.size() != std::size_t(training_observations)) throw std::runtime_error("total observation count");
    if (argc == 4 && std::string(argv[3]) == "prepare") return 0;
    std::ofstream decision(results / "identification.csv");
    decision << std::setprecision(17) << "stage,start,selected,comparable,J_min,delta_J,rank_ratio,max_correlation,max_resolution,max_endpoint_distance,passes,reason\n";
    std::ofstream resolution(results / "local-resolution.csv");
    resolution << std::setprecision(17) << "stage,start,coordinate,C_aa,span,J_min,delta_J,h,R\n";
    std::optional<r::FitResult> selected;
    std::size_t offset = 0;
    std::string outcome = "Not identified under the frozen free-coordinate rule; qualification remains unopened.";
    for (int q : {static_cast<int>(coordinate_declarations.size())}) {
        request.coordinates.coordinates = coordinate_declarations;
        const std::size_t n = starts.size();
        const std::vector<std::vector<double>> allocated(starts.begin() + offset, starts.begin() + offset + n);
        r::validate_fit_request(request);
        std::ofstream policies(results / "solver-policy.csv"); policies << std::setprecision(17) << "name,value,unit\n";
        for (std::size_t i = 0; i < policy::SOLVER_POLICIES.size(); ++i) policies << policy::SOLVER_POLICIES[i]->name << ',' << request.controls.solver_options.value(i) << ',' << policy::SOLVER_UNITS[i] << '\n';
        if (argc == 4) {
            std::ofstream diagnostics; if (!trace_required) { diagnostics.open(results / "prefit-diagnostics.csv"); diagnostics << std::setprecision(17) << "start,field,state,index,value\n"; }
            std::ofstream sensitivities(results / "sensitivity.csv"), summary(results / "sensitivity-summary.csv");
            sensitivities << std::setprecision(17) << (verify_lengths && !sulfate_training ? "ion,start,sigma,k_water,coordinate,coordinate_name,source_row,molality,observable,prediction,residual,prediction_derivative,b,weight,A,status,cause" : "ion,start,u,k_Cl,coordinate,coordinate_name,source_row,molality,observable,prediction,residual,prediction_derivative,b,weight,A,status,cause") << (sulfate_training ? ",source,salt,bucket,observed\n" : mixture_states ? (",source,observed,m_NaCl,m_" + ion_label + "Cl2,x_water,x_Na,x_" + ion_label + ",x_Cl\n") : "\n");
            summary << std::setprecision(17) << (verify_lengths && !sulfate_training ? "start,sigma,k_water,entries,sigma_weighted_norm,water_weighted_norm,objective_cost,correlation,passes,reason,evaluation_seconds" : "start,u,k_Cl,entries,weighted_norm,chloride_weighted_norm,objective_cost,correlation,passes,reason,evaluation_seconds" ) << (sulfate_training ? ",C_11,delta_J,h,lower_distance,upper_distance\n" : "\n");
            std::vector<std::size_t> rows(request.observations.size()); std::iota(rows.begin(), rows.end(), 0);
            bool complete = true;
            std::ofstream trace, verification, lengths, previous;
            Rows inherited; std::map<std::pair<std::size_t, std::string>, std::pair<DensityBranch, std::optional<std::size_t>>> branches;
            if (verify_lengths) {
                if (!sulfate_training) { inherited = read(leaf / "results/sr-mixture/sensitivity.csv"); previous.open(results / "length-previous.csv"); previous << "point,observation,field,previous,current,passes\n" << std::setprecision(17); }
                verification.open(results / "length-verification.csv"); lengths.open(results / "length-coordinates.csv");
                verification << std::setprecision(17) << "evaluation,observation,coordinate,exact,physical_partial_sum,step,centered_difference,passes\n";
                lengths << std::setprecision(17) << (sulfate_training ? "evaluation,k_water,reciprocal,active,map,passes\n" : "evaluation,sigma,packing,DH,Born,k_water,k_Cl,lower_sigma,map0,map1,map2,map3,map4,passes\n");
            }
            if (trace_required) { trace.open(results / "prefit-diagnostics.csv"); trace << std::setprecision(17) << "evaluation,observation,state,density_mol_m3,pressure_Pa,pressure_slope,branch,stable_roots,convergence,status,cause,reference0,reference1,reference_error,pure_solvent_lnphi\n"; }
            const auto equal = [&](double a, double b, double tolerance) { return std::isfinite(a) && std::isfinite(b) && std::abs(a-b) <= tolerance*std::max(1.0, std::abs(b)); };
            const auto verify = [&](const std::vector<double>& point, std::size_t evaluation, const r::ObservationEvaluation* screen) {
                const auto e = r::evaluate_coordinates(p, request.coordinates, point, {}, {}, request.controls.solver_options);
                HelmholtzModel model(e.parameters, e.active, {}, request.controls.solver_options); r::ReferenceCache references;
                std::vector<double> predictions;
                if (sulfate_training) {
                    const bool mapped = e.active.size() == 1 && e.map == std::vector<double>{1} && e.parameters.k_ij[H*count+water] == point[0] && e.parameters.k_ij[water*count+H] == point[0];
                    complete = complete && mapped;
                    lengths << evaluation << ',' << e.parameters.k_ij[H*count+water] << ',' << e.parameters.k_ij[water*count+H] << ',' << e.active.size() << ',' << e.map.at(0) << ',' << mapped << '\n';
                } else if (verify_lengths) {
                const auto& t = coordinate_declarations[0].targets;
                const std::vector<double> actual{e.parameters.sigma_constant[H], e.parameters.packing_diameter_angstrom[H], e.parameters.debye_huckel_diameter_angstrom[H], e.parameters.born_diameter_angstrom[H]};
                bool mapped = e.active.size() == 5 && e.map.size() == 10 && e.parameters.k_ij[H*count+Cl] == 0;
                for (std::size_t a = 0; a < actual.size(); ++a) mapped = mapped && equal(actual[a], t.at(a).coefficient*point[0], verification_equivalence) && e.map.at(a*q) == t[a].coefficient && e.map.at(a*q+1) == 0;
                mapped = mapped && e.map.at(8) == 0 && e.map.at(9) == 1 && e.parameters.k_ij[H*count+water] == point[1] && e.parameters.k_ij[water*count+H] == point[1]; complete = complete && mapped;
                lengths << evaluation; for (double v : actual) lengths << ',' << v;
                lengths << ',' << e.parameters.k_ij[H*count+water] << ',' << e.parameters.k_ij[H*count+Cl] << ',' << coordinate_declarations[0].lower;
                for (int a = 0; a < 5; ++a) lengths << ',' << e.map.at(a*q+(a == 4));
                lengths << ',' << mapped << '\n';
                }
                for (auto i : rows) {
                    const auto prediction = r::predict(model, request.observations[i], {}, {}, &references);
                    const auto* ref = prediction.reference ? &*prediction.reference : nullptr;
                    std::vector<std::pair<std::string, const State*>> states{{"activity", prediction.activity ? &prediction.activity->state : nullptr}, {"pure", ref && ref->pure_solvent_state ? &*ref->pure_solvent_state : nullptr}};
                    if (ref && ref->solvent_state && metadata[i].source == "binary") states.push_back({"solvent_reference", &*ref->solvent_state});
                    for (const auto& [name, s] : states) {
                        const std::string cause = prediction.message + (ref ? ref->solvent_cause : "reference unavailable") + (s ? "" : "; state unavailable");
                        trace << evaluation << ',' << i << ',' << name << ',' << (s ? s->molar_density : NAN) << ',' << (s ? s->pressure : NAN) << ',' << (s ? s->pressure_density_derivative : NAN) << ',' << (s ? int(s->density_branch) : -1) << ',' << (s && s->stable_root_count ? std::to_string(*s->stable_root_count) : "unavailable") << ',' << quote(s && s->density_convergence ? core::describe(*s->density_convergence) : "unavailable") << ',' << r::observation_status_name(prediction.status) << ',' << quote(cause) << ',' << (ref ? ref->log_fugacity[0] : NAN) << ',' << (ref ? ref->log_fugacity[1] : NAN) << ',' << (ref ? ref->convergence_error.value_or(NAN) : NAN) << ',' << (ref ? ref->pure_solvent_log_fugacity.value_or(NAN) : NAN) << '\n';
                        complete = complete && s && std::isfinite(s->molar_density) && s->molar_density > 0 && s->density_convergence.has_value() && std::isfinite(s->pressure_density_derivative) && s->pressure_density_derivative > 0;
                        if (verify_lengths && s && evaluation == std::size_t(verification_point+1)) branches[{i,name}] = {s->density_branch, s->stable_root_count};
                        if (verify_lengths && s && evaluation >= screen_points.size()+1) complete = complete && branches.at({i,name}) == std::make_pair(s->density_branch, s->stable_root_count);
                    }
                    const double y = prediction.value.value_or(NAN); predictions.push_back(y);
                    complete = complete && prediction.status == r::ObservationStatus::Available && std::isfinite(y) && ref;
                    if (!screen || !verify_lengths) continue;
                    const auto gradient = r::differentiate(model, request.observations[i], prediction);
                    for (int a = 0; a < q; ++a) {
                        double partial = 0; for (std::size_t b = 0; b < gradient.value.size(); ++b) partial += gradient.value[b]*e.map.at(b*q+a);
                        const double exact = screen->optimizer_gradient.size() > i*q+a ? screen->optimizer_gradient[i*q+a] : NAN;
                        const bool ok = gradient.status == r::ObservationStatus::Available && gradient.value.size() == (sulfate_training ? 1 : 5) && equal(partial, exact, verification_equivalence) && screen->predictions.size() > i && equal(y, screen->predictions[i], verification_equivalence);
                        verification << evaluation << ',' << i << ',' << a << ',' << exact << ',' << partial << ",0,nan," << ok << '\n'; complete = complete && ok;
                    }
                }
                return predictions;
            };
            for (std::size_t j = 0; j < screen_points.size(); ++j) {
                r::FitResult failed; r::ObservationEvaluation out;
                const auto begin = std::chrono::steady_clock::now();
                const bool evaluated = r::ObservationEvaluator(request, failed, rows).evaluate(screen_points[j], true, out);
                bool traced = true;
                if (replacement_reference) { const bool prior = complete; complete = true; verify(screen_points[j], j+1, nullptr); traced = complete; complete = prior && traced; }
                if (verify_lengths) {
                    verify(screen_points[j], j+1, &out);
                    if (j == std::size_t(verification_point)) for (int a = 0; a < q; ++a) for (std::size_t h = 0; h < verification_steps.size(); ++h) {
                        auto plus = screen_points[j], minus = plus; plus[a] += verification_steps[h]; minus[a] -= verification_steps[h];
                        const auto yp = verify(plus, screen_points.size()+1+2*verification_steps.size()*a+2*h, nullptr), ym = verify(minus, screen_points.size()+2+2*verification_steps.size()*a+2*h, nullptr);
                        for (auto i : rows) {
                            const double difference = (yp[i]-ym[i])/(2*verification_steps[h]), exact = out.optimizer_gradient.size() > i*q+a ? out.optimizer_gradient[i*q+a] : NAN;
                            const bool ok = equal(difference, exact, verification_difference_tolerance); complete = complete && ok;
                            verification << j+1 << ',' << i << ',' << a << ',' << exact << ",nan," << verification_steps[h] << ',' << difference << ',' << ok << '\n';
                        }
                    }
                }
                const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
                for (const auto& field : std::vector<std::pair<std::string, std::vector<std::vector<double>>>>{{"density", out.state_densities}, {"reference", out.reference_values}}) for (std::size_t i = 0; i < field.second.size(); ++i) for (std::size_t a = 0; a < field.second[i].size(); ++a) diagnostics << j+1 << ',' << field.first << ',' << i << ',' << a << ',' << field.second[i][a] << '\n';
                double cost = 0, correlation = 0; std::vector<double> norms(q, 0), weighted(rows.size()*q);
                bool available = evaluated && traced && out.predictions.size() == rows.size() && out.residuals.size() == rows.size()
                    && out.statuses.size() == rows.size() && out.optimizer_gradient.size() == rows.size()*q && out.residual_optimizer_gradient.size() == rows.size()*q;
                const auto value = [](const auto& v, std::size_t i) { return i < v.size() ? v[i] : NAN; };
                for (auto i : rows) {
                    const auto& o = request.observations[i]; const double w = request.weights[i];
                    const double y = value(out.predictions, i), residual = value(out.residuals, i);
                    if (verify_lengths && !sulfate_training && j%3 == 0) for (const auto& [field, current] : std::vector<std::pair<std::string, double>>{{"prediction",y}, {"prediction_derivative",value(out.optimizer_gradient,i*q+1)}, {"b",value(out.residual_optimizer_gradient,i*q+1)}, {"A",std::sqrt(w)*value(out.residual_optimizer_gradient,i*q+1)}}) {
                        const double old = std::stod(inherited.at(j*rows.size()*q+i*q).at(field)); const bool same = equal(current,old,verification_equivalence); complete = complete && same;
                        previous << j+1 << ',' << i << ',' << field << ',' << old << ',' << current << ',' << same << '\n';
                    }
                    const auto state = i < out.statuses.size() ? out.statuses[i] : r::ObservationStatus::NotEvaluated;
                    std::string cause; for (const auto& f : failed.failures) if (f.observation == i) cause += f.message + "; ";
                    cost += .5*w*residual*residual;
                    for (int a = 0; a < q; ++a) {
                        const double derivative = value(out.optimizer_gradient, i*q+a), b = value(out.residual_optimizer_gradient, i*q+a), A = std::sqrt(w)*b;
                        available = available && state == r::ObservationStatus::Available && std::isfinite(y) && std::isfinite(residual) && std::isfinite(derivative) && std::isfinite(b) && std::isfinite(A);
                        weighted[i*q+a] = A; norms[a] = std::hypot(norms[a], A);
                        const auto& row = metadata[i];
                        auto composition = o.composition;
                        if (mixture_states && !sulfate_training && row.source == "binary") { composition[water] = 1/p.molar_mass_kg_per_mol[water]; composition[H] = row.mSr; composition[Cl] = 2*row.mSr; const double total = std::accumulate(composition.begin(), composition.end(), 0.0); for (auto& x : composition) x /= total; }
                        sensitivities << ion_label << ',' << j+1 << ',' << screen_points[j][0] << ',' << (q == 2 ? screen_points[j][1] : k) << ',' << a << ',' << coordinate_names[a] << ',' << row.cell << ',' << o.molality << ',' << row.observable << ',' << y << ',' << residual << ',' << derivative << ',' << b << ',' << w << ',' << A << ',' << r::observation_status_name(state) << ',' << quote(cause);
                        if (sulfate_training) sensitivities << ',' << row.source << ',' << row.salt << ',' << row.bucket << ',' << o.observed;
                        if (mixture_states) sensitivities << ',' << row.source << ',' << o.observed << ',' << row.mNa << ',' << row.mSr << ',' << composition[water] << ',' << composition[sodium] << ',' << composition[H] << ',' << composition[Cl];
                        sensitivities << '\n';
                    }
                }
                available = available && std::isfinite(cost) && std::all_of(norms.begin(), norms.end(), [](double v) { return std::isfinite(v) && v > 0; });
                std::string reason = available ? "" : "incomplete_or_nonfinite_columns_or_norms";
                if (q == 2) for (auto i : rows) correlation -= (weighted[i*q]/norms[0])*(weighted[i*q+1]/norms[1]);
                if (reason.empty() && (!std::isfinite(correlation) || std::abs(correlation) > screening_max_correlation)) reason = "correlation";
                const double C11 = 1/(norms[0]*norms[0]), delta = std::max(identification_cost_floor, identification_cost_fraction*std::abs(cost)), h = std::sqrt(2*delta*C11);
                const double lower = screen_points[j][0]-coordinate_declarations[0].lower, upper = coordinate_declarations[0].upper-screen_points[j][0];
                if (sulfate_training && reason.empty() && (!(lower > identification_bound_distance_strict && upper > identification_bound_distance_strict) || !std::isfinite(C11) || !(C11 > 0) || !std::isfinite(h) || h > identification_max_h)) reason = "local_resolution";
                available = reason.empty(); complete = complete && available;
                summary << j+1 << ',' << screen_points[j][0] << ',' << (q == 2 ? screen_points[j][1] : k) << ',' << rows.size() << ',' << norms[0] << ',' << (q == 2 ? norms[1] : NAN) << ',' << cost << ',' << (sulfate_training ? NAN : correlation) << ',' << available << ',' << quote(reason) << ',' << seconds;
                if (sulfate_training) summary << ',' << C11 << ',' << delta << ',' << h << ',' << lower << ',' << upper;
                summary << '\n';
            }
            std::ofstream(results / "outcome.txt") << (complete ? "Prefit sensitivities available; fit not yet run." : "Prefit sensitivity failure; fitting and 2a not reached.") << '\n';
            return complete ? 0 : 1;
        }
        if (q == 2 || sulfate_training) {
            if (verify_lengths) { const auto checked = read(results / "length-verification.csv"); if (checked.size() != request.observations.size()*(screen_points.size()*q+q*verification_steps.size())) throw std::runtime_error("missing derivative verification"); for (const auto& row : checked) if (row.at("passes") != "1") throw std::runtime_error("derivative verification does not permit fitting"); }
            if (verify_lengths) for (const auto& file : (sulfate_training ? std::vector<std::string>{"length-coordinates.csv"} : std::vector<std::string>{"length-coordinates.csv", "length-previous.csv"})) { const auto checked = read(results / file); if (checked.empty()) throw std::runtime_error("empty length verification"); for (const auto& row : checked) if (row.at("passes") != "1") throw std::runtime_error("length verification does not permit fitting"); }
            if (verify_lengths) { std::map<std::pair<int,int>, bool> coverage; for (const auto& row : read(results / "prefit-diagnostics.csv")) coverage[{std::stoi(row.at("evaluation")),std::stoi(row.at("observation"))}] = true; if (coverage.size() != (screen_points.size()+2*q*verification_steps.size())*request.observations.size()) throw std::runtime_error("incomplete length diagnostic coverage"); }
            const auto screened = read(results / "sensitivity-summary.csv");
            if (screened.size() != screen_points.size()) throw std::runtime_error("missing complete screen");
            for (std::size_t j = 0; j < screened.size(); ++j) if (screened[j].at("passes") != "1" || std::stod(screened[j].at(verify_lengths && !sulfate_training ? "sigma" : "u")) != screen_points[j][0] || (!sulfate_training && std::stod(screened[j].at(verify_lengths ? "k_water" : "k_Cl")) != screen_points[j][1]) || (sulfate_training && (!std::isfinite(std::stod(screened[j].at("h"))) || std::stod(screened[j].at("h")) > identification_max_h || std::stoul(screened[j].at("entries")) != request.observations.size()))) throw std::runtime_error("screen does not permit fitting");
        }
        std::cout << "native Ceres stage " << q << ": " << n << " starts" << std::endl;
        const auto finite_columns = [&](const auto& column) {
            if (column.size() != request.observations.size()*q || !std::all_of(column.begin(), column.end(), [](double v) { return std::isfinite(v); })) return false;
            for (int a = 0; a < q; ++a) {
                double norm = 0; for (std::size_t i = 0; i < request.observations.size(); ++i) norm = std::hypot(norm, column[i*q+a]);
                if (!std::isfinite(norm) || norm == 0) return false;
            }
            return true;
        };
        const auto fits = r::fit_sequence(request, allocated, total_maximum_iterations, total_maximum_seconds);
        bool incomplete = fits.size() != n;
        double minimum = INFINITY; std::size_t best = 0;
        for (std::size_t j = 0; j < fits.size(); ++j) {
            const auto& f = fits[j]; const auto& c = f.covariance;
            status << std::setprecision(17) << q << ',' << offset+j+1 << ',' << f.status << ',' << f.usable << ',' << f.final_cost << ',' << f.iterations << ',' << f.residual_evaluations << ',' << f.jacobian_evaluations << ',' << c.status << ',' << c.rank << ',' << quote(f.message) << '\n'; status.flush();
            for (const auto& item : std::vector<std::pair<std::string, std::vector<double>>>{{"optimizer", f.optimizer}, {"physical", f.physical}, {"predictions", f.predictions}, {"residuals", f.residuals}, {"weighted_residuals", f.weighted_residuals}, {"optimizer_jacobian", f.optimizer_jacobian}, {"physical_jacobian", f.physical_jacobian}, {"singular_values", c.singular_values}, {"right_vectors", c.right_vectors}, {"left_vectors", c.left_vectors}, {"C", c.optimizer}, {"iteration_seconds", f.iteration_seconds}, {"iteration_costs", f.iteration_costs}, {"projected_gradient", {c.projected_gradient_norm}}}) dump(q, offset+j+1, item.first, item.second);
            dump(q, offset+j+1, "active_bounds", c.active_bounds);
            dump(q, offset+j+1, "active_bound_count", std::vector<double>{double(c.active_bounds.size())});
            for (std::size_t i = 0; i < f.statuses.size(); ++i) values << q << ',' << offset+j+1 << ",observation_status," << i << ',' << int(f.statuses[i]) << '\n';
            for (const auto& failure : f.failures) {
                std::ostringstream point; point << std::setprecision(17); for (double v : failure.optimizer) point << v << ' ';
                failures << q << ',' << offset+j+1 << ',' << failure.trial << ',' << failure.observation << ',' << r::observation_status_name(failure.status) << ',' << failure.jacobian << ',' << quote(point.str()) << ',' << quote(failure.message) << '\n';
            }
            const auto finite = [](const auto& x) { return std::all_of(x.begin(), x.end(), [](double v) { return std::isfinite(v); }); };
            incomplete = incomplete || !finite_columns(f.optimizer_jacobian) || f.statuses.size() != request.observations.size()
                || !std::all_of(f.statuses.begin(), f.statuses.end(), [](auto v) { return v == r::ObservationStatus::Available; }) || f.predictions.size() != request.observations.size() || !finite(f.predictions) || f.optimizer.size() != std::size_t(q) || !finite(f.optimizer) || f.status != "converged" || !f.usable || !std::isfinite(f.final_cost);
            if (f.status == "converged" && f.usable && std::isfinite(f.final_cost) && f.final_cost < minimum) { minimum = f.final_cost; best = j; }
        }
        const double delta = std::max(identification_cost_floor, identification_cost_fraction * std::abs(minimum));
        bool passes = !incomplete && std::isfinite(minimum);
        for (std::size_t j = 0; j < fits.size(); ++j) {
            const auto& f = fits[j]; const auto& c = f.covariance;
            const bool comparable = std::isfinite(minimum) && f.status == "converged" && f.final_cost <= minimum + delta;
            double ratio = c.singular_values.size() == std::size_t(q) && c.singular_values.front() > 0 ? c.singular_values.back()/c.singular_values.front() : 0;
            double correlation = 0, maximum = 0, distance = 0;
            std::string reason;
            const bool sensitivity = finite_columns(f.optimizer_jacobian);
            if (f.status != "converged" || !f.usable) reason = "incomplete_start";
            else if (!sensitivity) reason = "zero_or_unavailable_sensitivity";
            else if (c.rank != std::size_t(q) || !std::isfinite(ratio) || !(ratio > identification_rank_ratio_strict)) reason = "rank";
            else if (c.status != "available" || !c.active_bounds.empty() || c.optimizer.size() != std::size_t(q*q) || !std::all_of(c.optimizer.begin(), c.optimizer.end(), [](double v) { return std::isfinite(v); })) reason = "covariance_or_bound";
            const auto& coordinates = request.coordinates.coordinates;
            for (int a = 0; a < q && f.optimizer.size() == std::size_t(q); ++a) {
                auto coord = coordinates[a]; coord.lower = (coord.lower - coord.origin)/coord.scale; coord.upper = (coord.upper - coord.origin)/coord.scale; const double span = coord.upper - coord.lower;
                dump(q, offset+j+1, "bound_distances", std::vector<double>{f.optimizer[a]-coord.lower, coord.upper-f.optimizer[a]});
                if (reason.empty() && (f.optimizer[a] - coord.lower <= identification_bound_distance_strict/coord.scale || coord.upper - f.optimizer[a] <= identification_bound_distance_strict/coord.scale)) reason = "near_bound";
                if (std::isfinite(minimum) && fits[best].optimizer.size() == std::size_t(q)) distance = std::max(distance, std::abs(f.optimizer[a]-fits[best].optimizer[a])/span);
                if (c.optimizer.size() != std::size_t(q*q)) continue;
                const double C = c.optimizer[a*q+a], h = std::sqrt(2*delta*C), R = h/span;
                resolution << q << ',' << offset+j+1 << ',' << a << ',' << C << ',' << span << ',' << minimum << ',' << delta << ',' << h << ',' << R << '\n';
                if (reason.empty() && (!(C > 0) || !std::isfinite(R))) reason = "unavailable_resolution";
                if (reason.empty() && (verify_lengths ? R > identification_relative_resolution : h > identification_max_h/coord.scale)) reason = "local_resolution";
                maximum = std::max(maximum, R);
                for (int b = 0; b < a; ++b) {
                    const double rho = c.optimizer[a*q+b]/std::sqrt(C*c.optimizer[b*q+b]);
                    if (reason.empty() && !std::isfinite(rho)) reason = "unavailable_correlation";
                    correlation = std::max(correlation, std::abs(rho));
                }
            }
            if (reason.empty() && correlation > screening_max_correlation) reason = "correlation";
            if (reason.empty() && distance > (verify_lengths ? identification_relative_resolution : identification_max_endpoint_distance/(coordinates.front().upper - coordinates.front().lower)) && comparable) reason = "endpoint_distance";
            if (comparable && !reason.empty()) passes = false;
            decision << q << ',' << offset+j+1 << ',' << (j == best && std::isfinite(minimum)) << ',' << comparable << ',' << minimum << ',' << delta << ',' << ratio << ',' << (sulfate_training ? NAN : correlation) << ',' << maximum << ',' << distance << ',' << reason.empty() << ',' << reason << '\n';
        }
        values.flush(); decision.flush(); resolution.flush(); failures.flush();
        if (incomplete) { outcome = "Not identified: incomplete parameter estimation within the frozen work limits; qualification remains unopened."; break; }
        if (passes) { selected = fits[best]; outcome = "Identified binary endpoint under the frozen free-coordinate rule; accuracy and qualification remain required and unopened."; break; }
        offset += n;
    }
    std::ofstream(results / "outcome.txt") << outcome << '\n';
    std::cout << outcome << std::endl;
    return 0;
} catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
