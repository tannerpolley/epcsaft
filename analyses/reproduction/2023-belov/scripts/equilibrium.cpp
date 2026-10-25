#include "equilibrium/problem.hpp"
#include "surface/parameter_builder.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace epcsaft;
using namespace epcsaft::equilibrium;
int main(int argc, char** argv) {
    if (argc!=2) return 2;
    double T,P,P0,co,co2,tolerance;
    std::cin >> T >> P >> P0 >> co >> co2 >> tolerance;
    std::vector<double> g(8), expected(8),feed(8);
    std::vector<std::vector<double>> atoms(2,std::vector<double>(8));
    ParameterDocument doc;
    doc.base_choice="ideal-gas";
    doc.association_choice=doc.electrolyte_choice=doc.permittivity_choice=doc.polar_choice="none";
    ThermochemistryRecord record; record.standard_pressure=P0;
    for (std::size_t i=0;i<8;++i) {
        std::string name;
        std::cin >> name >> g[i] >> atoms[0][i] >> atoms[1][i] >> feed[i] >> expected[i];
        doc.component_ids.push_back(name);
        const double mass=(atoms[0][i]*(2*co-co2)+atoms[1][i]*(co2-co))/1000;
        doc.entries.push_back({name+"/mass","molar_mass",{i},mass,"kilogram / mole"});
        IdealPolynomial anchor{{0},T,GAS_CONSTANT_J_PER_MOL_K*T*g[i],0};
        record.component_intervals.push_back({{T,T+1,true,true,{anchor,P0}}});
    }
    EquilibriumProblem problem;
    problem.temperature=FreeValue::fixed(T); problem.pressure=FreeValue::fixed(P);
    problem.feed={FeedKind::Amounts,feed};
    PhaseSpec vapor; vapor.name="V"; vapor.kind=PhaseKind::Vapor;
    problem.phases={vapor}; problem.solver_options.set("EQUILIBRIUM_ROW_RESIDUAL",1e-10);
    for (int k=0;k<6;++k) {
        ReactionSpec reaction; reaction.stoichiometry.resize(8);
        for (auto& value:reaction.stoichiometry) std::cin >> value;
        problem.reactions.push_back(reaction);
    }
    if (!std::cin) throw std::runtime_error("incomplete Belov source declaration");
    HelmholtzModel model(build_parameters(doc),{},record);
    std::ofstream comparison(std::string(argv[1])+"/comparison.csv");
    comparison << std::setprecision(17) << "quantity,value,reference,error,tolerance,outcome,evidence_label\n";
    auto compare=[&](std::string quantity,double value,double reference,double error,double limit) {
        comparison << quantity << ',' << value << ',' << reference << ',' << error << ',' << limit << ','
                   << (std::isfinite(value) && std::isfinite(error) && std::abs(error)<=limit ? "pass":"fail")
                   << ",numerical verification of eight-gas restriction\n";
    };
    for (std::size_t i=0;i<8;++i) {
        std::vector<double> unit(8); unit[i]=1;
        auto anchor=reaction_thermochemistry_density(record,unit,T);
        const double value=anchor.pressure_gibbs/(GAS_CONSTANT_J_PER_MOL_K*T);
        compare(doc.component_ids[i]+"_G_over_RT",value,g[i],value-g[i],1e-10);
        if (anchor.status!=ThermochemistryStatus::Available || !std::isfinite(value) || std::abs(value-g[i])>1e-10)
            throw std::runtime_error("one-state Gibbs encoding failed");
    }
    const auto result=solve_equilibrium(model,problem);
    std::ofstream diagnostic(std::string(argv[1])+"/diagnostics.txt");
    diagnostic << std::setprecision(17) << "success " << result.success << " message " << result.message
               << "\nsolver " << result.solver_status << " iterations " << result.iterations << " evaluations " << result.evaluations << '\n';
    diagnostic << "kind phase index equation raw classified\n";
    for (std::size_t i=0;i<result.rows.size();++i) {
        const auto& row=result.rows[i];
        diagnostic << int(row.kind) << ' ' << row.phase << ' ' << row.index << ' ' << row_kind_name(row.kind) << ' ' << result.residuals[i] << ' ' << result.classified_residuals[i] << '\n';
    }
    if (!result.has_state()) return 1;
    const auto state=model.state(T,result.molar_densities[0],result.mole_fractions);
    diagnostic << "T_K P_Pa rho_mol_m3 phase_amount_mol dP_drho\n" << T << ' ' << state.pressure << ' ' << result.molar_densities[0] << ' ' << result.amounts[0] << ' ' << state.pressure_density_derivative << '\n';
    diagnostic << "species mole_fraction\n";
    for (std::size_t i=0;i<8;++i) diagnostic << doc.component_ids[i] << ' ' << result.mole_fractions[i] << '\n';
    std::vector<double> amounts(8);
    for (std::size_t i=0;i<8;++i) {
        amounts[i]=result.amounts[0]*result.mole_fractions[i];
        const double error=amounts[i]>0?std::abs(std::log(amounts[i])-std::log(expected[i])):INFINITY;
        compare(doc.component_ids[i]+"_amount_mol",amounts[i],expected[i],error,tolerance);
        const double mu=require_component_property(state.residual_chemical_potential_over_rt[i]);
        compare(doc.component_ids[i]+"_mu_res_RT",mu,0,mu,0);
    }
    for (std::size_t j=0;j<2;++j) {
        double defect=0;
        for (std::size_t i=0;i<8;++i) defect+=atoms[j][i]*(amounts[i]-feed[i]);
        compare(j?"oxygen_balance_mol":"carbon_balance_mol",defect,0,defect,1e-10);
    }
    compare("raw_residual",maximum_abs(result.residuals),0,maximum_abs(result.residuals),1e-10);
    compare("classified_residual",maximum_abs(result.classified_residuals),0,maximum_abs(result.classified_residuals),1e-10);
    compare("pressure_relative",state.pressure/P-1,0,state.pressure/P-1,1e-10);
    const double slope_defect=std::isfinite(state.pressure_density_derivative) && state.pressure_density_derivative>0?0:1;
    compare("mechanical_slope_defect",slope_defect,0,slope_defect,0);
}
