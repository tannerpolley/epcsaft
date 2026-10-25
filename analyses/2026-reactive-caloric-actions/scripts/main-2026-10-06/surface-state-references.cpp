#include "properties_checks_shared.hpp"
using namespace properties_checks;
#define main retained_surface_main
#include "surface_example.cpp"
#undef main
#include <Eigen/Dense>
#include <fstream>
#include <map>
std::map<std::string,double> surface_rows() {
    std::ifstream f("build/current-calorics/main-2026-10-06/native-rows/native-rows.tsv");
    std::map<std::string,double> rows;std::string line;
    while(std::getline(f,line)) if(line.rfind("surface_example\t",0)==0) {
        const auto a=line.find('\t'),b=line.find('\t',a+1);
        rows[line.substr(a+1,b-a-1)]=std::stod(line.substr(line.rfind('\t')+1));
    }return rows;
}
int main(){try {
    const auto rows=surface_rows();const auto value=[&](const std::string& id){return rows.at(id);};
    const double T=353.15,R=GAS_CONSTANT_J_PER_MOL_K,rho=value("state/root/0"),V=1/rho;
    const HelmholtzModel model(build_parameters(associating_binary()),{},ideal_record());
    const std::vector<double> x{.25,.75};
    const auto d=model.evaluate(T,V,x,DerivativeOrder::PhaseHessian);
    const auto width=d.gradient.size();
    const double P=R*T*(rho-d.gradient[1]),Prho=R*T*(1+V*V*d.hessian[width+1]);
    const double PT=R*(rho-d.gradient[1])-R*T*d.hessian[1];
    const double ih=.25*(1000+80*(T-300))+.75*(2000+35*(T-300));
    const double is=.25*(50+80*std::log(T/300)-R*std::log(.25*rho*R*T/1e5))
        +.75*(70+35*std::log(T/300)-R*std::log(.75*rho*R*T/1e5));
    const double cvres=-R*(2*T*d.gradient[0]+T*T*d.hessian[0]);
    const double cp=46.25+cvres+T*PT*PT/(rho*rho*Prho)-R;
    const double references[10]={ih,ih-R*T,is,46.25,46.25-R,
        ih-R*T*T*d.gradient[0]+P/rho-R*T,ih-R*T-R*T*T*d.gradient[0],
        is-R*(d.value+T*d.gradient[0]),cp,46.25-R+cvres};
    for(std::size_t k=0;k<10;++k) {
        const auto id="state/calorics/"+std::to_string(k);const double a=value(id),b=references[k];
        CHECK(recording_support::compare(id.c_str(),"caloric value","J/mol or J/(mol K)",a,b,
            "literal constant-cp anchors; T=353.15 K, x=(0.25,0.75), recorded rho",
            "Identity","primitive Helmholtz gradient/Hessian and ideal-gas identities",
            "original surface symmetric value criterion","symmetric_close",1e-12,1e-12,same(a,b)));
    }
    const double reduced=std::abs(P-1e5)/(rho*R*T),root_target=policy::DENSITY_ROOT_RESIDUAL.value;
    for(const auto& id:{"state/root/0","state/root/1","state/root/3"}) {
        const double a=value(id),b=std::string(id)=="state/root/0"?rho:std::string(id)=="state/root/1"?1e5:0;
        CHECK(recording_support::compare_custom(id,"density-root quantity","mixed",a,b,
            recording_support::fields({recording_support::field("metric_error",reduced),recording_support::field("target",root_target)}),
            "Numerical","native primitive EOS pressure","original density root admission",
            "positive liquid density and original reduced-pressure criterion",rho>0&&Prho>0&&reduced<=root_target&&value("state/root/3")<=root_target));
    }
    CHECK(recording_support::compare("state/root/2","pressure-density slope","Pa m3/mol",value("state/root/2"),Prho,"primitive EOS Hessian",
        "Identity","pressure gradient","original surface symmetric value criterion","symmetric_close",1e-12,1e-12,same(value("state/root/2"),Prho)));
    EquilibriumProblem bubble;
    PhaseSpec liquid,vapor;liquid.name="L";liquid.kind=PhaseKind::Liquid;vapor.name="V";vapor.kind=PhaseKind::Vapor;vapor.amount=AmountSpec::pinned(0);
    bubble.phases={liquid,vapor};bubble.temperature=FreeValue::fixed(T);bubble.pressure=FreeValue::variable(3e4);bubble.feed.kind=FeedKind::Amounts;bubble.feed.values=x;
    const auto boiling=solve_equilibrium(model,bubble);CHECK(boiling.success);
    auto flash=bubble;flash.phases[1].amount=AmountSpec::free();flash.pressure=FreeValue::fixed(1.8e4);
    const auto split=solve_equilibrium(model,flash);CHECK(split.success&&*split.requested_tolerance_met);
    std::vector<double> n{value("flash/amounts/0"),value("flash/amounts/1")},fractions;
    for(std::size_t k=0;k<4;++k)fractions.push_back(value("flash/mole_fractions/"+std::to_string(k)));
    double balance=std::abs(n[0]+n[1]-1),chemical=0;
    std::vector<double> mu;
    for(std::size_t p=0;p<2;++p) {
        const std::vector<double> y{fractions[2*p],fractions[2*p+1]};
        const auto state=model.state_at_pressure(T,1.8e4,y,p?DensityBranch::Vapor:DensityBranch::Liquid);
        for(std::size_t i=0;i<2;++i)mu.push_back(std::log(y[i]*state.molar_density)+*state.residual_chemical_potential_over_rt[i]);
    }
    for(std::size_t i=0;i<2;++i) {balance=std::max(balance,std::abs(n[0]*fractions[i]+n[1]*fractions[2+i]-x[i]));chemical=std::max(chemical,std::abs(mu[i]-mu[2+i]));}
    for(std::size_t k=0;k<2;++k) {
        const auto id="flash/amounts/"+std::to_string(k);const double lever=(x[0]-fractions[2])/(fractions[0]-fractions[2]);
        CHECK(recording_support::compare_custom(id.c_str(),"phase amount","mol",n[k],k?1-lever:lever,
            recording_support::fields({recording_support::field("metric_error",balance),recording_support::field("target",1e-12),recording_support::field("chemical_residual",chemical)}),
            "Relationship","species balances and independent EOS chemical potentials","original surface balances and driver chemical rows",
            "balance<=1e-12 mol and chemical rows<=row_tolerance",balance<=1e-12&&chemical<=flash.solver_options.value(policy::row)));
    }
    for(std::size_t k=0;k<4;++k) {
        const auto id="flash/mole_fractions/"+std::to_string(k);const auto p=k/2,i=k%2;
        const double b=(x[i]-n[1-p]*fractions[2*(1-p)+i])/n[p];
        CHECK(recording_support::compare_custom(id.c_str(),"mole fraction","1",fractions[k],b,
            recording_support::fields({recording_support::field("metric_error",balance),recording_support::field("target",1e-12),recording_support::field("chemical_residual",chemical)}),
            "Relationship","species balances and independent EOS chemical potentials","original surface balances and driver chemical rows",
            "balance<=1e-12 mol and chemical rows<=row_tolerance",balance<=1e-12&&chemical<=flash.solver_options.value(policy::row)));
    }
    for(std::size_t i=0;i<2;++i) {
        const auto id="flash/multipliers/"+std::to_string(i);
        CHECK(recording_support::compare(id.c_str(),"chemical multiplier","1",value(id),mu[i],"native EOS ln(x rho)+mu_res/RT",
            "Relationship","independent EOS chemical row","original equilibrium row criterion","absolute",flash.solver_options.value(policy::row),0,std::abs(value(id)-mu[i])<=flash.solver_options.value(policy::row)));
        const auto dual="flash/constraint_multipliers/"+std::to_string(i);
        CHECK(recording_support::compare(dual.c_str(),"balance objective dual","1",value(dual),-value(id),"unit absolute feed and reduced objective",
            "Identity","constraint dual sign and objective basis","original surface value criterion","symmetric_close",1e-12,1e-12,same(value(dual),-value(id))));
    }
    double bubble_chemical=0;
    for(std::size_t i=0;i<2;++i) {
        double potentials[2];for(std::size_t p=0;p<2;++p) {
            const std::vector<double> y{boiling.mole_fractions[2*p],boiling.mole_fractions[2*p+1]};
            const auto state=model.state_at_pressure(T,value("bubble/pressure"),y,p?DensityBranch::Vapor:DensityBranch::Liquid);
            potentials[p]=std::log(y[i]*state.molar_density)+*state.residual_chemical_potential_over_rt[i];
        }bubble_chemical=std::max(bubble_chemical,std::abs(potentials[0]-potentials[1]));
    }
    CHECK(recording_support::compare_custom("bubble/pressure","bubble pressure","Pa",value("bubble/pressure"),*boiling.pressure,
        recording_support::fields({recording_support::field("metric_error",bubble_chemical),recording_support::field("target",bubble.solver_options.value(policy::row))}),
        "Numerical","native EOS chemical-potential coexistence","original bubble chemical-row criterion","chemical residual<=row_tolerance",bubble_chemical<=bubble.solver_options.value(policy::row)));
    const auto compiled=compile_problem(model,flash);const auto prepared=compiled.prepare(split.coordinates,true);
    const auto data=compiled.jacobian_at(prepared);const auto count=split.coordinates.size();
    Eigen::MatrixXd J(count,count);for(std::size_t i=0;i<count;++i)for(std::size_t j=0;j<count;++j)J(i,j)=data[i*count+j];
    const auto power=[](double v){return v>0?std::exp2(-std::round(std::log2(v))):1.;};
    Eigen::VectorXd column=J.cwiseAbs().colwise().maxCoeff().transpose().unaryExpr(power);
    Eigen::MatrixXd A=J*column.asDiagonal();Eigen::VectorXd row_scale=A.cwiseAbs().rowwise().maxCoeff().unaryExpr(power);A=row_scale.asDiagonal()*A;
    const Eigen::MatrixXd inverse=A.colPivHouseholderQr().solve(Eigen::MatrixXd::Identity(count,count));
    const double condition=A.cwiseAbs().rowwise().sum().maxCoeff()*inverse.cwiseAbs().rowwise().sum().maxCoeff();
    CHECK(recording_support::compare("solved/condition","equilibrated infinity condition","1",value("solved/condition"),condition,
        "ALG-ACTIONS power-of-two R J C; independent QR inverse","Numerical","condition norm from native defining Jacobian",
        "original surface symmetric criterion","symmetric_close",1e-12,1e-12,same(value("solved/condition"),condition)));
    const double radial=n[0]*fractions[1];
    CHECK(recording_support::compare("solved/radial/0","radial component amount","mol",value("solved/radial/0"),radial,
        "Euler homogeneous amount relation","Identity","absolute-feed extensivity","original surface Euler criterion","symmetric_close",5e-6,5e-6,same(value("solved/radial/0"),radial,5e-6)));
    const double floor=value("floor/residual");
    CHECK(recording_support::compare_custom("floor/residual","achieved row residual","1",floor,1e-18,"original explicit unattainable row request",
        "Numerical","original unmet-tolerance discrimination","returned finite residual exceeds declared tolerance","finite residual>1e-18",std::isfinite(floor)&&floor>1e-18));
    return 0;
}catch(const std::exception& e){std::fprintf(stderr,"%s\n",e.what());return 1;}}
