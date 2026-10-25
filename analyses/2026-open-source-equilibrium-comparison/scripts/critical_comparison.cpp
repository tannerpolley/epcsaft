// Issue 224 B: C1/C3 critical-equation candidates against FeOs v0.10.1 and teqp v0.23.2, and E2 envelope
// points against FeOs bubble/dew points at common temperatures. The Engine solves natively; Python declares calls.
#include "properties/helmholtz_model.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "properties_support.hpp"
#include "equilibrium_support.hpp"
#include "caloric_specification_inputs.hpp"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
using namespace epcsaft;
using namespace epcsaft::equilibrium;
namespace {
std::vector<std::map<std::string, std::string>> read(const std::filesystem::path& path) {
    std::ifstream in(path);std::string line;std::getline(in,line);
    std::vector<std::string> names;std::istringstream header(line);
    for(std::string item;std::getline(header,item,',');)names.push_back(item);
    std::vector<std::map<std::string, std::string>> rows;
    while(std::getline(in,line)) {std::istringstream s(line);auto& row=rows.emplace_back();for(const auto& name:names)std::getline(s,row[name],',');}
    return rows;
}
EquilibriumProblem cell(const std::vector<double>& z) {
    auto p=equilibrium_support::flash(z.size());p.feed.values=z;p.solver_options.values[epcsaft::policy::row]=1e-11;
    p.temperature=FreeValue::variable(z.size()==1?375:260);p.pressure=FreeValue::variable(1e6);
    p.critical=CriticalDeclaration{z.size()==1?5000.:6000.,{}};
    return p;
}
EquilibriumResult critical(const HelmholtzModel& model,const std::vector<double>& z) {
    auto r=solve_equilibrium(model,cell(z));
    if(!r.has_state()||!r.critical)throw std::runtime_error("critical candidate: "+r.message);
    return r;
}
}  // namespace

int main(int argc,char** argv) {
    if(argc!=3)return 2;
    const std::filesystem::path root=argv[1];std::filesystem::create_directories(root);
    const HelmholtzModel propane(caloric_specification_checks::propane_parameters()),binary(properties_support::neutral_binary());
    // E1's level-h trace supplies the branch; each common T is corrected natively from its bracketing state.
    auto start=equilibrium_support::flash(2);start.feed.values={0.5,0.5};start.solver_options.values[epcsaft::policy::row]=1e-11;
    start.phases.push_back(equilibrium_support::phase("V",PhaseKind::Vapor,AmountSpec::pinned(0)));
    start.temperature=FreeValue::variable(120);start.pressure=FreeValue::variable(1e5);
    start.envelope=EnvelopeCondition{EnvelopeCoordinate::LogPressure,0,std::log(1e5)};
    ContinuationRequest request{start,start,std::nullopt,{1,1e-7,4000,6},{}};
    request.options.envelope_max_log_step=0.02;request.options.envelope_straddle=0.01;
    const auto trace=continue_equilibrium(binary,request);
    if(!trace.success)throw std::runtime_error("E1 trace: "+trace.message);
    std::vector<EquilibriumResult> points;
    for(const double T:{130.,180.,230.,260.,268.})
        for(std::size_t k=0;k+1<trace.states.size();++k) {
            const auto &a=trace.states[k],&b=trace.states[k+1];
            if(a.declaration->phases[0].kind!=b.declaration->phases[0].kind||(*a.temperature-T)*(*b.temperature-T)>0)continue;
            auto problem=*a.declaration;problem.envelope=EnvelopeCondition{EnvelopeCoordinate::LogTemperature,0,std::log(T)};
            points.push_back(solve_equilibrium_from_coordinates(binary,problem,a.coordinates));
            if(!points.back().has_state())throw std::runtime_error("E2 corrected cell at "+std::to_string(T)+" K: "+points.back().message);
        }
    {
        std::ofstream f(root/"reference-inputs.json");f<<std::setprecision(17)<<"{\"components\":[";
        const std::vector<std::pair<std::string,Parameters>> sources={{"propane",propane.parameters()},{"methane",binary.parameters()},{"ethane",binary.parameters()}};
        for(std::size_t k=0;k<sources.size();++k) {
            const auto& p=sources[k].second;const std::size_t i=k==2?1:0;
            f<<(k?",":"")<<"{\"name\":\""<<sources[k].first<<"\",\"mass\":"<<p.molar_mass_kg_per_mol[i]*1000<<",\"residual\":{\"m\":"<<p.segment_count[i]
             <<",\"sigma\":"<<p.sigma_constant[i]<<",\"epsilon_k\":"<<p.dispersion_energy_over_k[i]<<"}}";
        }
        f<<"],\"feos_critical_z0\":[0.5,0.27],\"teqp_locus_near_z0\":[0.5,0.27],\"teqp_trace_options\":{\"polish\":true,"
            "\"polish_exception_on_fail\":true,\"verbosity\":0,\"init_dt\":10.0,\"max_dt\":1000.0},\"envelope_points\":[";
        for(std::size_t k=0;k<points.size();++k)
            f<<(k?",":"")<<"{\"bubble\":"<<(points[k].declaration->phases[0].kind==PhaseKind::Liquid?"true":"false")
             <<",\"T\":"<<*points[k].temperature<<",\"P\":"<<*points[k].pressure<<",\"x0\":"<<points[k].mole_fractions[2]<<"}";
        f<<"]}";
    }
    const std::string command=std::string(argv[2])+" -I analyses/2026-open-source-equilibrium-comparison/scripts/critical_reference.py "+root.string();
    if(std::system(command.c_str())!=0)return 3;
    std::ofstream out(root/"comparisons.csv");out<<std::setprecision(17)<<"falsifier,quantity,engine,reference,scale,error,criterion,pass\n";
    bool all=true;
    const auto compare=[&](const std::string& id,const std::string& quantity,double engine,double reference,double scale,double criterion) {
        const double error=std::abs(engine-reference)/scale;const bool pass=std::isfinite(error)&&error<=criterion;all&=pass;
        out<<id<<","<<quantity<<","<<engine<<","<<reference<<","<<scale<<","<<error<<","<<criterion<<","<<pass<<"\n";
    };
    const auto candidate=[&](const std::string& id,const HelmholtzModel& model,const std::map<std::string,std::string>& row) {
        const double z0=std::stod(row.at("z0"));
        const auto r=critical(model,model.parameters().component_count()==1?std::vector<double>{1}:std::vector<double>{z0,1-z0});
        const std::string source=row.at("source")+"/z0="+row.at("z0");
        for(const auto& [quantity,engine]:std::vector<std::pair<std::string,double>>{{"T_K",*r.temperature},{"rho_mol_m3",r.molar_densities[0]},{"P_Pa",*r.pressure}})
            compare(id,source+"/"+quantity,engine,std::stod(row.at(quantity)),std::stod(row.at(quantity)),1e-8);
    };
    for(const auto& row:read(root/"critical-pure.csv"))candidate("C1",propane,row);
    for(const auto& row:read(root/"critical-binary.csv"))candidate("C3",binary,row);
    const auto references=read(root/"envelope-points.csv");
    for(std::size_t k=0;k<points.size();++k) {
        const auto& r=points[k];const auto& row=references.at(k);
        const std::string source=std::string(row.at("bubble")=="1"?"bubble":"dew")+"/T="+row.at("T_K")+"/";
        compare("E2",source+"P_Pa",*r.pressure,std::stod(row.at("P_Pa")),std::stod(row.at("P_Pa")),1e-7);
        compare("E2",source+"rho_feed",r.molar_densities[0],std::stod(row.at("rho_feed")),std::stod(row.at("rho_feed")),1e-7);
        compare("E2",source+"rho_incipient",r.molar_densities[1],std::stod(row.at("rho_incipient")),std::stod(row.at("rho_incipient")),1e-7);
        compare("E2",source+"x0_incipient",r.mole_fractions[2],std::stod(row.at("x0_incipient")),1,1e-7);
    }
    return all?0:1;
}
