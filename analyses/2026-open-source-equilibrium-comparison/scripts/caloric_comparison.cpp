// Issue 224 PH6/PS6: one common caloric increment; FeOs v0.10.1 native flashes.
#include "properties/helmholtz_model.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "properties_support.hpp"
#include "equilibrium_support.hpp"
#include "caloric_specification_inputs.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <cstdlib>
#include <numeric>
#include <map>
using namespace epcsaft;
using namespace epcsaft::equilibrium;
std::vector<double> calorics(const HelmholtzModel& model,const EquilibriumResult& result) {
    std::vector<double> value(2);
    for(std::size_t p=0;p<result.phase_count;++p) {
        const auto first=result.mole_fractions.begin()+p*2;
        const auto state=model.state(*result.temperature,result.molar_densities[p],{first,first+2});
        value[0]+=result.amounts[p]* *state.total_enthalpy;
        value[1]+=result.amounts[p]* *state.total_entropy;
    }
    return value;
}
std::vector<std::vector<double>> read(const std::filesystem::path& path) {
    std::ifstream f(path);std::string line;std::getline(f,line);std::vector<std::vector<double>> rows;
    while(std::getline(f,line)) {std::istringstream s(line);std::string item;std::vector<double> row;
        while(std::getline(s,item,','))row.push_back(std::stod(item));
        rows.push_back(row);}
    if(rows.size()!=2)throw std::runtime_error("reference must retain both rows");
    return rows;
}
int retain(const std::filesystem::path& root,const char* log) {
    struct Group {std::size_t count=0;double maximum=0;bool pass=true;};
    std::map<std::string,Group> groups;
    std::ifstream input(log);std::string line;
    std::ofstream raw(root/"native-comparisons.csv");raw<<std::setprecision(17)<<"comparison_id,quantity,unit,model,reference,scale,error,criterion,pass,inputs\n";
    // Read scalar fields from the Engine's fixed recording format, retaining the full input object.
    const auto field=[](const std::string& record,const char* key) {
        const auto marker=std::string("\"")+key+"\":";const auto found=record.find(marker);
        if(found==std::string::npos)return std::string();
        const auto start=found+marker.size();
        if(record[start]=='"')return record.substr(start+1,record.find('"',start+1)-start-1);
        const auto value=record.substr(start,record.find_first_of(",}",start)-start);
        return value=="null"?std::string():value;
    };
    while(std::getline(input,line)) {
        const std::string prefix="EPCSAFT_COMPARE_JSON ";if(line.rfind(prefix,0)!=0)continue;
        const auto id=field(line,"comparison_id");
        const auto start=line.find("\"inputs\":")+9,end=line.find(",\"basis\":",start);
        auto context=line.substr(start,end-start);std::string escaped;
        for(char c:context){if(c=='"')escaped+='"';escaped+=c;}
        raw<<id<<","<<field(line,"quantity")<<",\""<<field(line,"unit")<<"\","<<field(line,"model")<<","<<field(line,"reference")
           <<","<<field(context,"scale")<<","<<field(context,"error")<<","<<field(context,"criterion")<<","<<field(line,"verdict")<<",\""<<escaped<<"\"\n";
        if(field(line,"verdict").empty())continue;
        auto group=id.substr(0,id.find('/'));
        if(id.find("/F5/J/")!=std::string::npos)group="Assembled Jacobian";
        else if(id.find("/F5/H/")!=std::string::npos)group="Weighted Hessian";
        else if(group=="admission")group="Caloric admission and join actions";
        auto& result=groups[group];++result.count;result.pass&=field(line,"verdict")=="true";
        result.maximum=std::max(result.maximum,std::stod(field(context,"error"))/std::stod(field(context,"criterion")));
    }
    std::ifstream reference(root/"reference/comparisons.csv");std::getline(reference,line);
    while(std::getline(reference,line)) {
        std::istringstream row(line);std::vector<std::string> fields;std::string field;
        while(std::getline(row,field,','))fields.push_back(field);
        auto& result=groups[fields[0]];++result.count;result.pass&=fields[7]=="1";
        result.maximum=std::max(result.maximum,std::stod(fields[5])/std::stod(fields[6]));
    }
    std::ofstream csv(root/"qualification.csv"),md(root/"qualification-table.md");
    csv<<std::setprecision(17)<<"falsifier,count,maximum_error_over_criterion,pass\n";
    md<<std::setprecision(9)<<"| Falsifier | Numerical comparisons | Maximum error / criterion (limit: one) | Outcome |\n|---|---:|---:|---|\n";
    bool all=true;
    for(const auto& [name,r]:groups) {
        csv<<name<<","<<r.count<<","<<r.maximum<<","<<(r.pass?"True":"False")<<"\n";
        md<<"| "<<name<<" | "<<r.count<<" | "<<r.maximum<<" | "<<(r.pass?"Pass":"Fail")<<" |\n";all&=r.pass;
    }
    return input.eof()&&!groups.empty()&&all?0:1;
}

int main(int argc,char** argv) {
    if(argc==4 && std::string(argv[3])=="retain")return retain(argv[1],argv[2]);
    if(argc!=3)return 2;
    const std::filesystem::path root=argv[1];std::filesystem::create_directories(root);
    const HelmholtzModel model(properties_support::neutral_binary(),{},caloric_specification_checks::records());
    auto problem=equilibrium_support::flash(2);problem.phases.push_back(equilibrium_support::phase("V",PhaseKind::Vapor));
    problem.pressure=FreeValue::fixed(1e6);problem.solver_options.values[epcsaft::policy::row]=1e-11;
    std::vector<EquilibriumResult> tp;std::vector<std::vector<double>> hs;
    for(double T:{180.,185.}) {problem.temperature=FreeValue::fixed(T);tp.push_back(solve_equilibrium(model,problem));
        if(!tp.back().has_state())throw std::runtime_error(tp.back().message);
        hs.push_back(calorics(model,tp.back()));}
    auto inputs=[&](const std::vector<std::vector<double>>& reference,bool flashes) {
        std::ofstream f(root/"reference-inputs.json");f<<std::setprecision(17)<<"{\"components\":[";
        const auto& params=model.parameters();const auto& records=*model.thermochemistry();
        for(std::size_t i=0;i<2;++i) {
            auto coefficients=std::get<IdealDIPPR107>(records.component_intervals[i][0].correlation.form).coefficients;
            coefficients[0]*=1000;coefficients[1]*=1000;coefficients[3]*=1000;
            if(i)f<<",";
            f<<"{\"name\":\""<<(i?"ethane":"methane")<<"\",\"mass\":"<<params.molar_mass_kg_per_mol[i]*1000
             <<",\"residual\":{\"m\":"<<params.segment_count[i]<<",\"sigma\":"<<params.sigma_constant[i]<<",\"epsilon_k\":"<<params.dispersion_energy_over_k[i]
             <<"},\"ideal\":{\"DIPPR107\":[";
            for(std::size_t j=0;j<5;++j){if(j)f<<",";f<<coefficients[j];}f<<"]}}";
        }
        f<<"],\""<<(flashes?"flashes":"baselines")<<"\":[";
        for(std::size_t i=0;i<2;++i) {
            if(i)f<<",";
            f<<"{\"kind\":\""<<(flashes?(i?"PS":"PH"):"TP")<<"\",\"T\":"<<(flashes?180:180+5*i)<<",\"P\":1000000,\"feed\":[0.5,0.5]";
            if(flashes)f<<",\"target\":"<<reference[0][i+1]+hs[1][i]-hs[0][i];
            f<<"}";
        }
        f<<"]}";
    };
    const std::string command=std::string(argv[2])+" analyses/2026-open-source-equilibrium-comparison/scripts/caloric_reference.py "+root.string();
    inputs({},false);if(std::system((command+" baselines").c_str())!=0)return 3;
    const auto baseline=read(root/"baselines.csv");inputs(baseline,true);if(std::system((command+" flashes").c_str())!=0)return 4;
    const auto flashed=read(root/"flashes.csv");
    std::ofstream out(root/"comparisons.csv");out<<std::setprecision(17)<<"falsifier,quantity,engine,reference,scale,error,criterion,pass\n";
    bool all=true;auto compare=[&](const char* id,const std::string& name,double actual,double ref,double scale) {
        const double error=std::abs(actual-ref)/scale;const bool pass=std::isfinite(error)&&error<=1e-8;all&=pass;
        out<<id<<","<<name<<","<<actual<<","<<ref<<","<<scale<<","<<error<<",1e-8,"<<pass<<"\n";};
    for(std::size_t i=0;i<2;++i) {
        const char* id=i?"PS6":"PH6";const double increment=hs[1][i]-hs[0][i];
        compare(id,"increment",increment,baseline[1][i+1]-baseline[0][i+1],GAS_CONSTANT_J_PER_MOL_K*(i?1:185));
        auto caloric=problem;caloric.temperature=FreeValue::variable(180);if(i)caloric.entropy_j_per_k=hs[0][i]+increment;else caloric.enthalpy_j=hs[0][i]+increment;
        const auto result=solve_equilibrium(model,caloric);if(!result.has_state())throw std::runtime_error(result.message);
        compare(id,"temperature",*result.temperature,flashed[i][0],185);
        compare(id,"vapor-fraction",result.amounts[1],flashed[i][3],1);
        for(std::size_t j=0;j<4;++j)compare(id,"composition-"+std::to_string(j),result.mole_fractions[j],flashed[i][4+j],1);
    }
    std::ofstream engine(root/"engine-baselines.csv");engine<<std::setprecision(17)<<"T_K,H_J,S_J_per_K,beta,xL0,xL1,xV0,xV1\n";
    for(std::size_t i=0;i<2;++i){engine<<*tp[i].temperature<<","<<hs[i][0]<<","<<hs[i][1]<<","<<tp[i].amounts[1];for(auto x:tp[i].mole_fractions)engine<<","<<x;engine<<"\n";}
    return all?0:1;
}
