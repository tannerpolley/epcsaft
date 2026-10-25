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
    if (argc != 3 && argc != 4) return 2;
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
    const double k = pair_correction(temperature_K, p, temperature_law_values(temperature_K, p), TemperatureLaw::Kij, H, Cl);
    const double energy = std::sqrt(p.dispersion_energy_over_k[H] * p.dispersion_energy_over_k[Cl]) * (1 - k);
    if (k == 1 && energy != 0) throw std::runtime_error("nonzero H-Cl cross dispersion energy at k=1");
    if (!sulfate_training) std::ofstream(results / (ion_label + "-Cl-dispersion.csv")) << "k_ion_Cl,cross_energy_over_k_K\n" << std::setprecision(17) << k << ',' << energy << '\n';
    std::vector<double> solvent(count, 0.0); solvent[water] = 1.0;
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
