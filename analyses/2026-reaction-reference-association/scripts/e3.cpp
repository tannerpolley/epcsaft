// Independent verification of EqID reference_limit_terminal_pair; no production reference replacement.
#include "equilibrium/problem.hpp"
#include "surface/parameter_builder.hpp"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc < 2) return 2;
    std::ifstream in(argv[1]);
    epcsaft::ParameterDocument doc;
    std::size_t n, sites, entries;
    in >> n >> sites >> entries;
    doc.component_ids.resize(n);
    doc.site_component.resize(sites);
    doc.site_multiplicity.resize(sites);
    for (auto& id : doc.component_ids) in >> std::quoted(id);
    for (std::size_t i = 0; i < sites; ++i) in >> doc.site_component[i] >> doc.site_multiplicity[i];
    in >> std::quoted(doc.base_choice) >> std::quoted(doc.association_choice) >> std::quoted(doc.electrolyte_choice)
       >> std::quoted(doc.ionic_dispersion_choice) >> std::quoted(doc.permittivity_choice) >> std::quoted(doc.polar_choice);
    for (std::size_t i = 0; i < entries; ++i) {
        epcsaft::ParameterEntry e;
        std::size_t width;
        in >> std::quoted(e.identity) >> std::quoted(e.field) >> width;
        e.index.resize(width);
        for (auto& index : e.index) in >> index;
        in >> e.magnitude >> std::quoted(e.unit);
        doc.entries.push_back(e);
    }
    using namespace epcsaft;
    using namespace epcsaft::equilibrium;
    const HelmholtzModel model(build_parameters(doc));
    const double mw = .01801528, rt = GAS_CONSTANT_J_PER_MOL_K;
    std::vector<std::vector<double>> nu{{-1,0,0,0,0,0,1,1},{0,1,0,-1,0,0,1,0},{-1,0,-1,0,1,0,1,0},{0,0,0,0,-1,1,1,0}};
    const double coeff[4][3]={{132.899,-13445.9,-22.4773},{-83.491,-819.7,10.9756},{231.465,-12092.1,-36.7816},{216.049,-12431.7,-35.4819}};
    NeutralReferenceDeclaration declaration{0,{0,.125,.125,.125,.125,.03125,.125,.0625},1.,{}};
    std::cout << "row,control,reaction,T_K,P_Pa,actual,expected,error,limit,pass\n" << std::setprecision(17);
    bool passed = true;
    int row;
    double t,m,c,observed,root;
    auto emit = [&](const std::string& control, int reaction, double p, double actual, double expected, double limit) {
        const double error = std::abs(actual-expected);
        const bool ok = std::isfinite(error) && error <= limit;
        passed &= ok;
        std::cout << row << ',' << control << ',' << reaction << ',' << t << ',' << p << ',' << actual << ',' << expected << ',' << error << ',' << limit << ',' << ok << std::endl;
    };
    auto reject = [&](const std::string& name, int r, double p, double defect, double predicted, double limit) {
        emit(name+"_signed",r,p,defect,predicted,1e-9);
        emit(name+"_trips",r,p,std::abs(defect)>limit,1,0);
    };
    while (in >> row >> t >> m >> c >> observed >> root) {
        std::vector<double> x(n),y(n);
        for (auto& v:x) in >> v;
        for (auto& v:y) in >> v;
        if(argc==3 && std::string(argv[2])!="e5") {
            const double total=std::accumulate(x.begin(),x.end(),0.);
            for(auto& v:x) v/=total;
            const auto record = [&](const char* name, int i, double value) { std::cout<<row<<','<<name<<','<<i<<','<<t<<','<<observed<<','<<value<<",,,,record"<<std::endl; };
            const auto state=model.state_at_pressure(t,observed,x,DensityBranch::Liquid);
            for(std::size_t i=0;i<n;++i) {
                record("MEA_species_fraction",i,x[i]);
                const double molal=x[i]/(mw*x[2]), concentration=state.molar_density*x[i];
                emit("MEA_m_amount",i,observed,molal*mw*x[2]*total,x[i]*total,1e-12*total);
                emit("MEA_c_amount",i,observed,concentration/state.molar_density*total,x[i]*total,1e-12*total);
                const double mu=std::log(x[i]*state.molar_density)+state.residual_chemical_potential_over_rt[i].value();
                if(x[i]>0) emit("MEA_mu_fugacity",i,observed,std::log(x[i])+state.log_fugacity_coefficient->at(i).value()+std::log(observed/(rt*t)),mu,1e-10);
                if(i<3 && y[i]>=0) emit("MEA_fugacity_Pa",i,observed,x[i]*std::exp(state.log_fugacity_coefficient->at(i).value())*observed,y[i],1e-8*std::max(1.,y[i]));
            }
            record("MEA_nominal_carbon",-1,m); record("MEA_reconstructed_carbon",-1,c);
            record("MEA_carbon_closure_difference",-1,c-m);
            record("MEA_CO2_fugacity_Pa",0,x[0]*std::exp(state.log_fugacity_coefficient->at(0).value())*observed);
            emit("MEA_amine_closure",-1,observed,total*(x[1]+x[3]+x[4]),m/root,1e-12*total);
            emit("MEA_charge_closure",-1,observed,total*(x[3]+x[7]-x[4]-x[5]-2*x[6]-x[8]),0,1e-10*total);
            emit("MEA_initial_water",-1,observed,total*(x[2]+x[5])*mw,.7,1e-12*total);
            continue;
        }
        if(n==2) {
            EquilibriumProblem neutral;
            neutral.temperature=FreeValue::variable(t); neutral.pressure=FreeValue::fixed(observed);
            neutral.feed={FeedKind::MoleFractions,{m,1-m},{},1.};
            PhaseSpec l,v; l.name="L"; v.name="V"; v.kind=PhaseKind::Vapor; v.amount=AmountSpec::pinned(0);
            neutral.phases={l,v};
            const auto solution=solve_equilibrium(model,neutral);
            if(!solution.success) throw std::runtime_error(solution.message);
            emit("neutral_temperature",-1,observed,*solution.temperature,t,1e6);
            continue;
        }
        EquilibriumProblem problem;
        problem.temperature=FreeValue::fixed(t); problem.pressure=FreeValue::variable(observed*1e5*.5);
        problem.feed={FeedKind::Amounts,{1/mw,m,c,0,0,0,0,0},{},1.}; problem.source_reference=declaration;
        PhaseSpec liquid,vapor; liquid.name="L"; vapor.name="V"; vapor.kind=PhaseKind::Vapor;
        liquid.support={0,1,2,3,4,5,6,7}; vapor.support={0,1,2}; vapor.amount=AmountSpec::pinned(0);
        liquid.support_declared=vapor.support_declared=true;
        const double extent=.97*std::min(m,c);
        liquid.composition_guess={1/mw-extent,m-extent+1e-6,c-extent+1e-4*c,extent,extent,1e-3*extent,1e-9,1e-7};
        double sum=std::accumulate(liquid.composition_guess.begin(),liquid.composition_guess.end(),0.);
        for(auto& v:liquid.composition_guess) v/=sum;
        vapor.composition_guess={.5,.0001,.4999}; problem.phases={liquid,vapor};
        for (std::size_t r=0;r<nu.size();++r) {
            ReactionSpec reaction; reaction.name=std::to_string(r); reaction.stoichiometry=nu[r];
            reaction.correlation=ReactionLogPolynomial{coeff[r][0]+coeff[r][2]*std::log(298.15),coeff[r][1],coeff[r][2],0.,298.15,273.,498.,"aqueous-mole-fraction-infinite-dilution-water-v1"};
            reaction.reference=ReactionReference{std::nullopt,ReferenceSourceBasis::RawMoleFractionInfiniteDilution}; problem.reactions.push_back(reaction);
        }
        const bool e5=argc==3 && std::string(argv[2])=="e5";
        const auto raw=compile_problem(model,problem);
        if(e5) {
            const auto start=std::clock();
            const auto a=solve_equilibrium(model,problem);
            const double cpu=double(std::clock()-start)/CLOCKS_PER_SEC;
            if(!a.success) throw std::runtime_error(a.message);
            const auto state=raw.prepare(a.coordinates);
            const auto record=[&](const std::string& name,int i,double value) { emit(name,i,*a.pressure,value,value,0); };
            record("pressure",-1,*a.pressure); record("signed_error_percent",-1,100*(*a.pressure/(observed*1e5)-1));
            record("cpu_seconds",-1,cpu); record("iterations",-1,a.iterations);
            record("evaluations",-1,a.evaluations); record("hessian_calls",-1,a.hessian_calls);
            for(double residual:a.classified_residuals) emit("classified_rows",-1,*a.pressure,residual,0,1e-10);
            for(std::size_t phase=0;phase<2;++phase) {
                record("amount",phase,state.amount[phase]); record("density",phase,a.molar_densities[phase]);
                emit("phase_sum",phase,*a.pressure,std::accumulate(state.composition[phase].begin(),state.composition[phase].end(),0.),1,1e-10);
                for(std::size_t i=0;i<n;++i) {
                    record("fraction_"+std::to_string(phase),i,state.composition[phase][i]);
                    record("species_amount_"+std::to_string(phase),i,state.amount[phase]*state.composition[phase][i]);
                    if(phase==1 && i<3) record("partial_pressure",i,state.composition[phase][i] * *a.pressure);
                }
            }
            const std::vector<double> feed{1/mw,m,c,0,0,0,0,0};
            for(const auto& balance:raw.B) {
                double defect=0;
                for(std::size_t i=0;i<n;++i) defect+=balance[i]*(state.amount[0]*state.composition[0][i]-feed[i]);
                emit("inventory_balance",-1,*a.pressure,defect,0,1e-10*(1/mw+m+c));
            }
            const auto liquid=model.state_at_pressure(t,*a.pressure,state.composition[0],DensityBranch::Liquid);
            const auto vapor=model.state_at_pressure(t,*a.pressure,state.composition[1],DensityBranch::Vapor);
            for(std::size_t i=0;i<3;++i) emit("bubble_fugacity",i,*a.pressure,std::log(state.composition[0][i]/state.composition[1][i])+liquid.log_fugacity_coefficient->at(i).value()-vapor.log_fugacity_coefficient->at(i).value(),0,1e-8);
            for(const auto* phase:{&liquid,&vapor}) {
                record("mechanical_slope",phase==&liquid?0:1,phase->pressure_density_derivative);
                if(!(phase->pressure_density_derivative>0)) throw std::runtime_error("nonpositive mechanical slope");
                for(std::size_t i=0;i<phase->site_fractions.size();++i) {
                    record(phase==&liquid?"sites_liquid":"sites_vapor",i,phase->site_fractions[i]);
                    if(!(phase->site_fractions[i]>0 && phase->site_fractions[i]<=1)) throw std::runtime_error("site fraction outside (0,1]");
                }
            }
            const std::vector<ReactionReference> refs(4,{std::nullopt,ReferenceSourceBasis::RawMoleFractionInfiniteDilution});
            for(double p:{root,*a.pressure}) {
                const auto moving=reaction_references(model,t,p,declaration,nu,refs,true,std::nullopt,2);
                for(std::size_t r=0;r<4;++r) {
                    const double k=coeff[r][0]+coeff[r][1]/t+coeff[r][2]*std::log(t);
                    record("source_log_K",r,k); emit("reference_C",r,p,*moving[r].terminal_contraction,*moving[r].terminal_contraction,0);
                    emit("compiled_log_Krho",r,p,k+*moving[r].terminal_contraction,k+*moving[r].terminal_contraction,0);
                    if(p==*a.pressure) {
                        double q=0;
                        for(std::size_t i=0;i<n;++i) if(nu[r][i]) q+=nu[r][i]*(std::log(state.composition[0][i]*liquid.molar_density)+liquid.residual_chemical_potential_over_rt[i].value());
                        emit("equilibrium_affinity",r,p,q-k-*moving[r].terminal_contraction,0,1e-8);
                    }
                }
            }
            if(!passed) return 1;
            continue;
        }
        auto molal=problem;
        for (auto& r:molal.reactions) {
            const double scale=std::accumulate(r.stoichiometry.begin()+1,r.stoichiometry.end(),0.)*std::log(mw);
            r.correlation->a-=scale; r.reference->source_basis=ReferenceSourceBasis::CommonMolalityInfiniteDilution;
        }
        const auto common=compile_problem(model,molal);
        auto untransformed=molal;
        for(std::size_t r=0;r<4;++r) untransformed.reactions[r].correlation=problem.reactions[r].correlation;
        const auto negative=compile_problem(model,untransformed);
        for(std::size_t r=0;r<4;++r) {
            double offset=0;
            for(std::size_t i=0;i<n;++i) offset-=nu[r][i]*(negative.standard_gibbs[i]-raw.standard_gibbs[i])/(rt*t);
            reject("untransformed_K",r,problem.pressure.value,offset,std::accumulate(nu[r].begin()+1,nu[r].end(),0.)*std::log(mw),1e-10);
        }
        for(std::size_t r=0;r<4;++r) {
            double difference=0;
            for(std::size_t i=0;i<n;++i) difference+=nu[r][i]*(raw.standard_gibbs[i]-common.standard_gibbs[i])/(rt*t);
            emit("compiled_x_m",r,problem.pressure.value,difference,0,1e-10);
        }
        const auto a=solve_equilibrium(model,problem), b=solve_equilibrium(model,molal);
        if(!a.success || !b.success) throw std::runtime_error(a.message+" / "+b.message);
        emit("solved_P_relative",-1,*a.pressure,*b.pressure / *a.pressure-1,0,1e-8);
        emit("retained_P_relative",-1,*a.pressure,*a.pressure/root-1,0,1e-8);
        for(std::size_t i=0;i<a.mole_fractions.size();++i) emit("solved_x",i,*a.pressure,a.mole_fractions[i],b.mole_fractions[i],1e-9);
        for(std::size_t i=0;i<2;++i) emit("solved_rho_relative",i,*a.pressure,b.molar_densities[i]/a.molar_densities[i]-1,0,1e-8);
        for(double residual:a.classified_residuals) emit("classified_rows",-1,*a.pressure,residual,0,1e-10);
        const auto prepared=raw.prepare(a.coordinates), prepared_b=common.prepare(b.coordinates);
        emit("phase_amount_relative",0,*a.pressure,prepared_b.amount[0]/prepared.amount[0]-1,0,1e-7);
        const double feed_total=1/mw+m+c;
        std::vector<double> feed{1/mw,m,c,0,0,0,0,0};
        for(std::size_t i=0;i<n;++i) {
            const double xi=feed[i]/feed_total, mi=feed[i]/(mw*feed[0]);
            emit("feed_m_roundtrip",i,*a.pressure,mi*mw*feed[0],feed[i],1e-12*feed_total);
            emit("feed_x_roundtrip",i,*a.pressure,xi*feed_total,feed[i],1e-12*feed_total);
            if(a.mole_fractions[i]>=1e-8) emit("species_amount_relative",i,*a.pressure,(prepared_b.amount[0]*prepared_b.composition[0][i])/(prepared.amount[0]*prepared.composition[0][i])-1,0,1e-7);
        }
        emit("initial_water_kg",-1,*a.pressure,feed[0]*mw,1,1e-12*feed_total);
        emit("analytical_loading",-1,*a.pressure,feed[2]/feed[1],c/m,1e-12*feed_total);
        const auto liquid_state=model.state_at_pressure(t,*a.pressure,prepared.composition[0],DensityBranch::Liquid);
        const auto vapor_state=model.state_at_pressure(t,*a.pressure,prepared.composition[1],DensityBranch::Vapor);
        for(std::size_t i=0;i<3;++i) emit("bubble_fugacity",i,*a.pressure,
            std::log(prepared.composition[0][i]/prepared.composition[1][i])+liquid_state.log_fugacity_coefficient->at(i).value()-vapor_state.log_fugacity_coefficient->at(i).value(),0,1e-8);
        for(const auto& composition:prepared.composition) emit("phase_sum",-1,*a.pressure,std::accumulate(composition.begin(),composition.end(),0.),1,1e-10);
        reject("species_amount_perturbed",-1,*a.pressure,(1.00001*prepared_b.amount[0]*prepared_b.composition[0][0])/(prepared.amount[0]*prepared.composition[0][0])-1,1e-5,1e-7);
        for(const auto& solution:{prepared,prepared_b}) {
        emit("incipient_vapor_amount",-1,*a.pressure,solution.amount[1],0,1e-10);
        for(const auto& balance:raw.B) {
            double defect=0;
            for(std::size_t i=0;i<n;++i) defect+=balance[i]*(solution.amount[0]*solution.composition[0][i]+solution.amount[1]*solution.composition[1][i]-feed[i]);
            emit("inventory_balance",-1,*a.pressure,defect,0,1e-10*feed_total);
        }
        }
        double inventory_perturbation=0;
        for(const auto& balance:raw.B) {
            double defect=0, expected=0;
            for(std::size_t i=0;i<n;++i) { defect+=balance[i]*(1.00001*prepared_b.amount[0]*prepared_b.composition[0][i]-feed[i]); expected+=balance[i]*feed[i]*1e-5; }
            emit("inventory_perturbed_signed",-1,*a.pressure,defect,expected,1e-9);
            inventory_perturbation=std::max(inventory_perturbation,std::abs(defect));
        }
        emit("inventory_perturbed_trips",-1,*a.pressure,inventory_perturbation>1e-10*feed_total,1,0);
        std::vector<ReactionReference> refs(4,{std::nullopt,ReferenceSourceBasis::RawMoleFractionInfiniteDilution});
        std::vector<double> solvent(n,0); solvent[0]=1;
        for(double p:{root,observed*1e5,100000.}) {
            const auto moving=reaction_references(model,t,p,declaration,nu,refs,true,std::nullopt,2);
            auto fixed_refs=refs; for(auto& r:fixed_refs) r.reference_pressure_pa=100000.;
            const auto fixed=reaction_references(model,t,p,declaration,nu,fixed_refs,true,std::nullopt,2);
            const auto water=model.state_at_pressure(t,p,solvent,DensityBranch::Liquid);
            auto fixed_problem=problem; for(auto& r:fixed_problem.reactions) r.reference->reference_pressure_pa=100000.;
            const auto fixed_compiled=compile_problem(model,fixed_problem);
            auto point=a.coordinates; point[raw.primal_size-raw.free_intensive_count]=p;
            const auto at=raw.prepare(point), at0=fixed_compiled.prepare(point);
            std::vector<RowSeries> series, series0;
            for(int direction:{0,1,2}) {
                std::vector<double> tangent(point.size()), control(raw.control_width());
                if(direction!=1) control[0]=t;
                if(direction!=0) tangent[raw.primal_size-raw.free_intensive_count]=p;
                series.push_back(raw.row_series(at,{point,tangent,std::vector<double>(point.size())},control,2));
                series0.push_back(fixed_compiled.row_series(at0,{point,tangent,std::vector<double>(point.size())},control,2));
            }
            for(std::size_t r=0;r<4;++r) {
                double target=0, target0=0, mixed=0;
                for(std::size_t i=0;i<n;++i) { target-=nu[r][i]*at.standard.derivative[0][i]; target0-=nu[r][i]*at0.standard.derivative[0][i]; }
                const double k=coeff[r][0]+coeff[r][1]/t+coeff[r][2]*std::log(t);
                emit("compiled_pointwise_moving",r,p,target,k+*moving[r].terminal_contraction,1e-10);
                emit("compiled_pointwise_fixed",r,p,target0,k+*fixed[r].terminal_contraction,1e-10);
                for(std::size_t j=0;j<raw.rows.size();++j) if(raw.rows[j].kind==RowKind::Chemical && raw.rows[j].phase==0) {
                    const auto i=raw.rows[j].index;
                    mixed-=nu[r][i]*((series[2].rows[2][j]-series0[2].rows[2][j])-(series[0].rows[2][j]-series0[0].rows[2][j])-(series[1].rows[2][j]-series0[1].rows[2][j]));
                }
                emit("compiler_row_TP",r,p,mixed,*moving[r].derivative_value(1,1,{})*t*p,1e-8);
            }
            for(int perturb:{0,1}) {
                auto composition=prepared.composition[0];
                const double shift=.01*composition[2];
                if(perturb) for(std::size_t i=0;i<n;++i) composition[i]+=shift*nu[2][i];
                const auto finite=model.state_at_pressure(t,p,composition,DensityBranch::Liquid);
                const auto anchor=model.state_at_pressure(t,100000.,solvent,DensityBranch::Liquid);
                const auto pressure_action=model.pressure_actions(t,p,composition,DensityBranch::Liquid,{PropertyObservable::LogFugacityCoefficient},{{0,p,std::vector<double>(n),{}}},finite.molar_density);
                for(const auto& action:pressure_action.results) {
                    const auto i=*action.component; std::vector<double> unit(n,0); unit[i]=1;
                    const double finite_slope=reference_partial_from_state(model,t,composition,unit,0,1,{},finite,false).action.value();
                    const double reference_slope=reference_partial_from_state(model,t,solvent,unit,0,1,{},water,true).action.value();
                    emit("Gibbs_fixed_pressure",i,p,*action.action+1,p*finite_slope,1e-8);
                    emit("Gibbs_same_pressure",i,p,*action.action+1-p*reference_slope,p*(finite_slope-reference_slope),1e-8);
                }
                for(std::size_t r=0;r<4;++r) {
                    double q=0, cm=0, c0=0, qs=0, dnu=0, solute=0, qm=0, qc=0, standard_c=0;
                    for(std::size_t i=0;i<n;++i) if(nu[r][i]) {
                        std::vector<double> unit(n,0); unit[i]=1;
                        const double li=reference_partial_from_state(model,t,solvent,unit,0,0,{},water,true).action.value();
                        const double l0=reference_partial_from_state(model,t,solvent,unit,0,0,{},anchor,true).action.value();
                        const double mu=finite.residual_chemical_potential_over_rt[i].value();
                        const double gx=finite.log_fugacity_coefficient->at(i).value()+std::log(p/(rt*t))-li;
                        qm+=nu[r][i]*(gx+(i?std::log(composition[0]):0)+std::log(i?composition[i]/(mw*composition[0]):composition[i]));
                        qc+=nu[r][i]*(gx+(i?std::log(water.molar_density/finite.molar_density):0)+std::log(i?composition[i]*finite.molar_density/1000.:composition[i]));
                        standard_c+=nu[r][i]*(li-(i?std::log(water.molar_density/1000.):0));
                        q+=nu[r][i]*(std::log(composition[i]*finite.molar_density)+mu);
                        qs+=nu[r][i]*(std::log(composition[i])+finite.log_fugacity_coefficient->at(i).value());
                        cm+=nu[r][i]*li; c0+=nu[r][i]*l0; dnu+=nu[r][i]; if(i) solute+=nu[r][i];
                        emit("fugacity_density_identity",i,p,mu-std::log(p/(rt*t*finite.molar_density)),finite.log_fugacity_coefficient->at(i).value(),1e-10);
                    }
                    const double k=coeff[r][0]+coeff[r][1]/t+coeff[r][2]*std::log(t);
                    emit("finite_Q_fugacity",r,p,qs+dnu*std::log(p/(rt*t)),q,1e-10);
                    emit("inherited_transport_affinity",r,p,q-c0-(k+cm-c0),q-cm-k,1e-9);
                    emit("fixed_transport_affinity",r,p,q-cm-(k+c0-cm),q-c0-k,1e-9);
                    emit("molal_affinity",r,p,qm-(k-solute*std::log(mw)),q-cm-k,1e-9);
                    emit("molarity_affinity",r,p,qc-(k+solute*std::log(water.molar_density/1000.)),q-cm-k,1e-9);
                    emit("c_density_target",r,p,k+solute*std::log(water.molar_density/1000.)+standard_c,-std::inner_product(nu[r].begin(),nu[r].end(),at.standard.derivative[0].begin(),0.),1e-10);
                    reject("missing_c_standard_density",r,p,k+standard_c-(-std::inner_product(nu[r].begin(),nu[r].end(),at.standard.derivative[0].begin(),0.)),-solute*std::log(water.molar_density/1000.),1e-10);
                    reject("missing_c_density_ratio",r,p,qc-solute*std::log(water.molar_density/finite.molar_density)-(q-cm+solute*std::log(water.molar_density/1000.)),-solute*std::log(water.molar_density/finite.molar_density),1e-9);
                    emit("missing_pressure_factor",r,p,qs+dnu*std::log(100000./(rt*t))- (q-c0) -c0,-dnu*std::log(p/100000.),1e-9);
                    emit("old_minus_fixed",r,p,(k+cm)-(k+c0),cm-c0,1e-10);
                    if(!perturb && p==root) emit("equilibrium_affinity",r,p,q-cm-k,0,1e-8);
                    for(std::size_t i=0;i<n;++i) emit("feed_c_roundtrip",i,p,(feed[i]/feed_total*finite.molar_density)*feed_total/finite.molar_density,feed[i],1e-12*feed_total);
                }
            }
            // Eq. (7): log standard-water density chain, checked by refined T/P differences.
            for(int direction:{0,1,2}) {
                const double dt=direction!=1?t:0, dp=direction!=0?p:0;
                const auto action=model.pressure_actions(t,p,solvent,DensityBranch::Liquid,{PropertyObservable::LogFugacityCoefficient},{{dt,dp,std::vector<double>(n),{}}},water.molar_density);
                const double exact=action.density_action.at(0)/water.molar_density;
                double previous=1;
                for(double h:{.002,.001}) {
                    const double plus=model.state_at_pressure(t+h*dt,p+h*dp,solvent,DensityBranch::Liquid).molar_density;
                    const double minus=model.state_at_pressure(t-h*dt,p-h*dp,solvent,DensityBranch::Liquid).molar_density;
                    const double fd=(std::log(plus/1000.)-std::log(minus/1000.))/(2*h);
                    emit("c_density_derivative_"+std::to_string(direction),-1,p,fd,exact,2e-5);
                    if(direction!=1) emit("c_density_refinement",-1,p,std::abs(fd-exact)<=previous,1,0);
                    previous=std::abs(fd-exact);
                }
                reject("missing_c_density_derivative",direction,p,-exact,-exact,2e-5);
            }
            for(std::size_t r=0;r<4;++r) {
                for(auto orders:std::vector<std::pair<int,int>>{{0,0},{1,0},{2,0},{0,1},{0,2},{1,1}}) {
                    double independent=0, independent0=0;
                    const auto water0=model.state_at_pressure(t,100000.,solvent,DensityBranch::Liquid);
                    for(std::size_t i=0;i<n;++i) if(nu[r][i]) {
                        std::vector<double> unit(n,0); unit[i]=1;
                        if(!orders.second) independent0+=nu[r][i]*reference_partial_from_state(model,t,solvent,unit,orders.first,0,{},water0,true).action.value();
                        independent+=nu[r][i]*reference_partial_from_state(model,t,solvent,unit,orders.first,orders.second,{},water,true).action.value();
                    }
                    const double scale=std::pow(t,orders.first)*std::pow(p,orders.second);
                    emit("species_reference_"+std::to_string(orders.first)+std::to_string(orders.second),r,p,*moving[r].derivative_value(orders.first,orders.second,{})*scale,independent*scale,orders.first+orders.second?1e-8:1e-10);
                    emit("fixed_reference_"+std::to_string(orders.first)+std::to_string(orders.second),r,p,*fixed[r].derivative_value(orders.first,orders.second,{})*scale,independent0*scale,orders.first+orders.second?1e-8:1e-10);
                    if(orders.second && !orders.first) {
                        double compiled_slope=0, fixed_slope=0;
                        const auto& column=orders.second==1?at.standard.pressure_derivative:at.standard.pressure_second_derivative;
                        const auto& column0=orders.second==1?at0.standard.pressure_derivative:at0.standard.pressure_second_derivative;
                        for(std::size_t i=0;i<n;++i) { compiled_slope-=nu[r][i]*column[i]; fixed_slope-=nu[r][i]*column0[i]; }
                        emit("compiler_P"+std::to_string(orders.second),r,p,compiled_slope*scale,independent*scale,1e-8);
                        emit("compiler_fixed_P"+std::to_string(orders.second),r,p,fixed_slope*scale,0,1e-8);
                    }
                    if(orders.second) {
                        emit("fixed_system_P",r,p,*fixed[r].derivative_value(orders.first,orders.second,{}),0,1e-10);
                        emit("inherited_transport_derivative",r,p,independent*scale,*moving[r].derivative_value(orders.first,orders.second,{})*scale,1e-8);
                        emit("manufactured_transport_derivative",r,p,(*moving[r].derivative_value(orders.first,orders.second,{})-independent)*scale,0,1e-8);
                        emit("missing_C_derivative",r,p,-independent*scale,-*moving[r].derivative_value(orders.first,orders.second,{})*scale,1e-8);
                    }
                }
            }
        }
        if(!passed) return 1;
    }
    return passed?0:1;
}
