// Issue 224: manufactured calorics on the Gross (2001) methane/ethane and propane packets.
#include "equilibrium_checks_shared.hpp"
#include "caloric_specification_inputs.hpp"
#include "caloric_specification_shared.hpp"
#include <cstdlib>
#include <filesystem>
namespace caloric_specification_checks {
using namespace epcsaft;
using namespace epcsaft::equilibrium;
inline double total(const HelmholtzModel& m, const EquilibriumResult& r, bool h) {
    double value=0; const std::size_t c=m.parameters().component_count();
    for(std::size_t p=0;p<r.phase_count;++p) {
        std::vector<double> x(r.mole_fractions.begin()+p*c,r.mole_fractions.begin()+(p+1)*c);
        const auto s=m.state(*r.temperature,r.molar_densities[p],x);
        value+=r.amounts[p]*(h?*s.total_enthalpy:*s.total_entropy);
    }
    return value;
}
inline EquilibriumProblem condition(EquilibriumProblem p, bool h, double target, double seed) {
    p.temperature=FreeValue::variable(seed); p.solver_options.values[epcsaft::policy::row] = 1e-11;
    if(h) p.enthalpy_j=target; else p.entropy_j_per_k=target;
    return p;
}
inline EquilibriumResult recovered(const HelmholtzModel& m,const EquilibriumProblem& p,const EquilibriumResult& ref,
                      bool h,double seed,const std::string& id) {
    const double target=total(m,ref,h), feed=std::accumulate(p.feed.values.begin(),p.feed.values.end(),0.0);
    const auto request=condition(p,h,target,seed); const auto result=solve_equilibrium(m,request);
    if(!result.has_state()) throw std::runtime_error(id+": "+result.message);
    compare(id+"/T",*result.temperature,*ref.temperature,*ref.temperature,1e-8);
    compare(id+"/closure",total(m,result,h),target,feed*GAS_CONSTANT_J_PER_MOL_K*(h?*result.temperature:1),1e-10);
    for(std::size_t i=0;i<result.amounts.size();++i) compare(id+"/split-"+std::to_string(i),result.amounts[i],ref.amounts[i],feed,1e-9);
    for(std::size_t i=0;i<result.mole_fractions.size();++i) compare(id+"/x-"+std::to_string(i),result.mole_fractions[i],ref.mole_fractions[i],1,1e-8);
    return result;
}
inline EquilibriumProblem binary() {
    auto p=equilibrium_support::flash(2); p.phases.push_back(equilibrium_support::phase("V",PhaseKind::Vapor));
    p.temperature=FreeValue::fixed(180);p.pressure=FreeValue::fixed(1e6);p.solver_options.values[epcsaft::policy::row] = 1e-11;
    return p;
}
void roundtrips() {
    const HelmholtzModel m(properties_support::neutral_binary(),{},records());
    auto p=binary();const auto base=solve_equilibrium(m,p);CHECK(base.has_state());
    const auto missing = std::filesystem::current_path() / "build/issue-1164/missing-ipopt-options.opt";
    CHECK(!std::filesystem::exists(missing));
    const char* setting = std::getenv("EPCSAFT_IPOPT_OPTIONS_FILE");
    const auto previous = setting ? std::optional<std::string>(setting) : std::nullopt;
    CHECK(setenv("EPCSAFT_IPOPT_OPTIONS_FILE", missing.c_str(), 1) == 0);
    EquilibriumStart complete; complete.coordinates = base.coordinates;
    const auto failed = solve_equilibrium_starts(m, p, {complete}, {});
    if (previous) CHECK(setenv("EPCSAFT_IPOPT_OPTIONS_FILE", previous->c_str(), 1) == 0);
    else CHECK(unsetenv("EPCSAFT_IPOPT_OPTIONS_FILE") == 0);
    CHECK(!failed.selected && failed.results[0]);
    const auto& diagnostic = *failed.results[0];
    const bool retained = diagnostic.diagnostic_options_file == missing.string()
        && diagnostic.solver_status == "UNASSIGNED"
        && diagnostic.ipopt_diagnostics && diagnostic.solver_settings.empty() && diagnostic.linear_solver.empty()
        && diagnostic.message.find("EPCSAFT_IPOPT_OPTIONS_FILE must name a readable regular file: " + missing.string()) != std::string::npos;
    recording_support::compare("ipopt-options/missing", "retained-diagnostics", "1", retained ? 1.0 : 0.0, 1.0,
        recording_support::fields({recording_support::field("cause", diagnostic.message),
            recording_support::field("options_file", diagnostic.diagnostic_options_file),
            recording_support::field("solver_status", diagnostic.solver_status)}),
        "Relationship", "driver-diagnostics", "a configuration refusal retains the driver observations", "absolute", 0.0, 0.0, retained);
    for(bool h:{true,false}) for(double n:{1.0,2.5}) {
        auto q=p;q.feed.values={0.5*n,0.5*n};const auto ref=solve_equilibrium(m,q);CHECK(ref.has_state());
        recovered(m,q,ref,h,165,(h?"PH1/":"PS1/")+std::to_string(n));
        for(std::size_t i=0;i<ref.amounts.size();++i) compare((h?"PH1/":"PS1/")+std::to_string(n)+"/amount-"+std::to_string(i),ref.amounts[i],n*base.amounts[i],n,1e-9);
    }
    // #610: the fixed-T seed at 240 K fails and the feed has only a vapor root there (phases_coincide before).
    recovered(m,p,base,true,240,"PH1/240");
    const HelmholtzModel pure(propane_parameters(),{},records(true));
    for(bool vapor:{false,true}) {
        auto q=equilibrium_support::flash(1);q.phases[0].kind=vapor?PhaseKind::Vapor:PhaseKind::Liquid;
        q.temperature=FreeValue::fixed(vapor?350:250);q.pressure=FreeValue::fixed(vapor?1e5:1e6);q.solver_options.values[epcsaft::policy::row] = 1e-11;
        const auto ref=solve_equilibrium(pure,q);CHECK(ref.has_state());
        for(bool h:{true,false}) recovered(pure,q,ref,h,vapor?330:240,(h?"PH2/":"PS2/")+std::string(vapor?"vapor":"liquid"));
    }
    auto sat=equilibrium_support::flash(1);sat.feed.values={2.5};sat.temperature=FreeValue::variable(300);
    sat.pressure=FreeValue::fixed(1e6);sat.solver_options.values[epcsaft::policy::row] = 1e-11;
    sat.phases.push_back(equilibrium_support::phase("V",PhaseKind::Vapor,AmountSpec::pinned(0)));
    auto ref=solve_equilibrium(pure,sat);CHECK(ref.has_state());ref.amounts={2.5*0.65,2.5*0.35};
    auto q=sat;q.phases[1].amount=AmountSpec{};
    const auto liquid=pure.state(*ref.temperature,ref.molar_densities[0],{1});
    const auto vapor=pure.state(*ref.temperature,ref.molar_densities[1],{1});
    const double hl=*liquid.total_enthalpy,hv=*vapor.total_enthalpy,sl=*liquid.total_entropy,sv=*vapor.total_entropy;
    for(bool h:{true,false}) {
        const auto result=recovered(pure,q,ref,h,300,h?"PH3":"PS3");
        for(std::size_t i=0;i<2;++i) CHECK(pure.state(*result.temperature,result.molar_densities[i],{1}).pressure_density_derivative>0);
        compare(std::string(h?"PH3":"PS3")+"/H-lever",(total(pure,result,true)/2.5-hl)/(hv-hl),0.35,1,1e-8);
        compare(std::string(h?"PH3":"PS3")+"/S-lever",(total(pure,result,false)/2.5-sl)/(sv-sl),0.35,1,1e-8);
    }
    // #579: Gross (2001) pure benzene and toluene, k_ij -0.003 from hda-pc-saft/1 (Yu 2026, Table 2), NASA Glenn (2002)
    // NASA9 calorics at 1 bar. The start 1.05 x 370 K lies above the dew point, so the fixed-T seed has no liquid; it
    // stops at its own budget and is discarded (66 -> 24).
    auto bt=properties_support::parameters(2);
    bt.segment_count={2.4653,2.8149};bt.sigma_constant={3.6478,3.7169};bt.dispersion_energy_over_k={287.35,285.69};
    bt.molar_mass_kg_per_mol={0.078114,0.092141};bt.k_ij[1]=bt.k_ij[2]=-0.003;
    ThermochemistryRecord nasa;nasa.standard_pressure=100000.0;
    for(auto a:std::vector<std::array<double,9>>{
            {-167734.0902,4404.50004,-37.1737791,0.1640509559,-0.0002020812374,1.307915264e-07,-3.4442841e-11,-10354.55401,216.9853345},
            {-287796.222,6133.94152,-45.7470676,0.1936895724,-0.0002304305304,1.459301178e-07,-3.7907961e-11,-23084.02499,269.3915042}}) {
        for(auto& v:a) v*=8.31451/GAS_CONSTANT_J_PER_MOL_K;  // the record's gas constant to the Engine's
        nasa.component_intervals.push_back({{200,1000,true,true,{IdealNASA9{a},100000}}});
    }
    const HelmholtzModel aromatics(bt,{},nasa);
    auto flash=binary();flash.temperature=FreeValue::fixed(370);flash.pressure=FreeValue::fixed(1e5);
    const auto at370=solve_equilibrium(aromatics,flash);CHECK(at370.has_state());
    const auto ph=recovered(aromatics,flash,at370,true,1.05*370,"PH5/benzene-toluene");
    CHECK(ph.iterations<=26);
}
// The identical ethane pair A=B at 300 K, 1 MPa: h_B-h_A = 5000 J/mol and equal c_p; the reaction law by kind
// (0 affinity, 1 constant K, 2 ln K polynomial on the EOS standard state).
std::pair<ThermochemistryRecord,EquilibriumProblem> ethane_reaction(int kind) {
    auto r=records();const double R=GAS_CONSTANT_J_PER_MOL_K;
    for(auto& species:r.component_intervals) species[0].correlation={IdealPolynomial{{40},298.15,0,0},100000};
    r.component_intervals[1][0].correlation={IdealPolynomial{{40},298.15,kind==1?0.:5000.,kind==0?0.:R*std::log(2.)},100000};
    auto p=equilibrium_support::flash(2);p.temperature=FreeValue::fixed(300);p.pressure=FreeValue::fixed(1e6);
    p.phases[0].kind=PhaseKind::Vapor;  // vapor-side sole root at 300 K, 1 MPa (#423)
    p.feed.values={2.5,0};p.solver_options.values[epcsaft::policy::row] = 1e-11;
    ReactionSpec rx;rx.name="A=B";rx.stoichiometry={-1,1};
    if(kind==0) rx.affinity_j_per_mol=5000;else if(kind==1) rx.equilibrium_constant=2;
    else rx.correlation=ReactionLogPolynomial{std::log(2.),-5000/R,0,0,298.15,100,550,EOS_STANDARD_STATE_ID};
    p.reactions={rx};return {r,p};
}
void reactions() {
    for(int kind=0;kind<3;++kind) {
        auto [r,p]=ethane_reaction(kind);const double R=GAS_CONSTANT_J_PER_MOL_K;const auto rx=p.reactions[0];
        const HelmholtzModel m(e2::identical_ethane_pair(),{},r);
        const auto ref=solve_equilibrium(m,p);CHECK(ref.has_state());
        for(bool h:{true,false}) recovered(m,p,ref,h,280,(h?"PH5/":"PS5/")+std::to_string(kind));
        const auto target=reaction_target_series(300,p.reactions,2);
        const auto thermal=reaction_thermochemistry_density(r,rx.stoichiometry,300);CHECK(thermal.status==ThermochemistryStatus::Available);
        compare("PH5/"+std::to_string(kind)+"/gibbs",thermal.gibbs,target[0][0],5000,1e-8);
        compare("PH5/"+std::to_string(kind)+"/enthalpy",thermal.enthalpy,target[0][0]-300*target[1][0],5000,1e-8);
        compare("PH5/"+std::to_string(kind)+"/cp",thermal.heat_capacity,-300*target[2][0],40,1e-8);
        auto& polynomial=std::get<IdealPolynomial>(r.component_intervals[1][0].correlation.form);polynomial.reference_enthalpy+=100;
        const HelmholtzModel wrong(e2::identical_ethane_pair(),{},r);
        bool rejected=false;try{(void)compile_problem(wrong,condition(p,true,total(m,ref,true),280));}
        catch(const std::invalid_argument& e){rejected=std::string(e.what()).find("reference inconsistency")!=std::string::npos;}
        compare("PH5/"+std::to_string(kind)+"/wrong-anchor",rejected?1:0,1,1,1e-10);
    }
}
inline SolvedStateActionDirection direction(const CompiledSystem& c) {
    return {0,0,std::vector<double>(c.n_feed.size()),std::vector<double>(c.model->active_parameters().size())};
}
inline double action(const CompiledSystem& c,const EquilibriumResult& r,const std::vector<SolvedStateActionDirection>& d,
                     SolvedStateObservableKind kind=SolvedStateObservableKind::Temperature) {
    const auto batch=solved_state_actions(c,r,{{{kind}},d});
    const auto& value=batch.results[0];
    if(!value.action) throw std::runtime_error("PH/PS action unavailable: "+value.diagnostic_message.value_or("no diagnostic"));
    // Numbered within the running case, so a row's id does not depend on which cases ran before it.
    static const char* numbered_case=nullptr;static std::size_t call=0;
    if(numbered_case!=recording_support::active_test_id()) {numbered_case=recording_support::active_test_id();call=0;}
    recording_support::observe(("PH4/action-call-"+std::to_string(++call)).c_str(),kind==SolvedStateObservableKind::SystemEnthalpy?"system-enthalpy-action":"solved-temperature-action",value.unit.c_str(),*value.action,std::nullopt,
        recording_support::fields({recording_support::field("order",d.size()),recording_support::field("T_K",*r.temperature),
            recording_support::field("P_Pa",*r.pressure),recording_support::field("condition",batch.diagnostics.condition_estimate.value_or(0)),
            recording_support::field("backward_error",batch.diagnostics.maximum_linear_residual.value_or(0)),
            recording_support::field("energy_order",batch.diagnostics.maximum_energy_order)}),
        "Observation","implicit-exact-action","conditioning and recurrence diagnostics","caloric-specification");
    return *value.action;
}
void actions() {
    const HelmholtzModel m(properties_support::neutral_binary(),{},records());
    auto p=binary();const auto tp=solve_equilibrium(m,p);const auto tc=compile_problem(m,p);
    auto endpoint=p;endpoint.temperature.value=185;const auto end_tp=solve_equilibrium(m,endpoint);CHECK(end_tp.has_state());
    for(bool h:{true,false}) {
        auto q=condition(p,h,total(m,tp,h),165);const auto c=compile_problem(m,q);const auto r=solve_equilibrium(m,q);
        auto d=direction(c);if(h)d.enthalpy_j=c.feed_total*GAS_CONSTANT_J_PER_MOL_K*180;else d.entropy_j_per_k=c.feed_total*GAS_CONSTANT_J_PER_MOL_K;
        auto pressure=direction(c);pressure.pressure_pa=1e5;auto feed=direction(c);feed.feed_amounts_mol={0.13,0.09};
        if(h) {
            auto thermal=direction(tc);thermal.temperature_k=1;
            compare("PH4/inverse",action(c,r,{d})/d.enthalpy_j,1/action(tc,tp,{thermal},SolvedStateObservableKind::SystemEnthalpy),1/std::abs(action(tc,tp,{thermal},SolvedStateObservableKind::SystemEnthalpy)),1e-8);
        }
        for(int order=2;order<=3;++order) {
            std::vector<SolvedStateActionDirection> ds={d,pressure};if(order==3)ds.push_back(feed);
            const double exact=action(c,r,ds);double chart_scale=std::abs(exact);
            if(order==3) {
                chart_scale=0;
                for(std::size_t i=0;i<feed.feed_amounts_mol.size();++i) {
                    auto basis=direction(c);basis.feed_amounts_mol[i]=feed.feed_amounts_mol[i];
                    chart_scale+=std::abs(action(c,r,{d,pressure,basis}));
                }
            }
            const auto moving=ds.back();ds.pop_back();
            double previous=0,previous_richardson=0,best=std::numeric_limits<double>::infinity();bool resolved=false;
            for(int rung=0;rung<6;++rung) {
                const double step=std::ldexp(0.02,-rung);double shifted[2];
                for(int sign:{-1,1}) {
                    auto v=q;v.pressure.value+=sign*step*moving.pressure_pa;
                    for(std::size_t i=0;i<v.feed.values.size();++i)v.feed.values[i]+=sign*step*moving.feed_amounts_mol[i];
                    if(h)v.enthalpy_j=*q.enthalpy_j+sign*step*moving.enthalpy_j;else v.entropy_j_per_k=*q.entropy_j_per_k+sign*step*moving.entropy_j_per_k;
                    const auto cc=compile_problem(m,v);const auto rr=solve_equilibrium_from_coordinates(m,v,r.coordinates);
                    shifted[sign==1]=action(cc,rr,ds);
                }
                const double reference=(shifted[1]-shifted[0])/(2*step),scale=std::max(chart_scale,std::abs(reference));
                const double error=std::abs(exact-reference)/scale;best=std::min(best,error);
                recording_support::observe(((h?"PH4/":"PS4/")+std::to_string(order)+"/rung-"+std::to_string(rung)).c_str(),"action-ladder","K",exact,reference,
                    recording_support::fields({recording_support::field("step",step),recording_support::field("error",error),recording_support::field("scale",scale)}),"Observation","lower-exact-action-difference","mixed H/S,P,feed refinement","caloric-specification");
                if(rung>0) {
                    const double richardson=(4*reference-previous)/3;
                    best=std::min(best,std::abs(exact-richardson)/scale);
                    if(rung>1&&std::abs(richardson-previous_richardson)/scale<2.5e-7)resolved=true;
                    recording_support::observe(((h?"PH4/":"PS4/")+std::to_string(order)+"/Richardson-"+std::to_string(rung)).c_str(),"refined-action","K",exact,richardson,recording_support::fields({recording_support::field("step",step),recording_support::field("scale",scale),recording_support::field("refinement",std::abs(richardson-previous_richardson)/scale)}),"Observation","centered-second-order-extrapolation","lower exact action refinement","caloric-specification");
                    previous_richardson=richardson;
                }
                previous=reference;
            }
            compare((h?"PH4/":"PS4/")+std::to_string(order)+"/refined",best,0,1,1e-6);CHECK(resolved);
        }
        auto end=q;
        if(h)end.enthalpy_j=total(m,end_tp,true);else end.entropy_j_per_k=total(m,end_tp,false);
        ContinuationRequest path{q,end,r,{0.5,1e-5,10,3},{}};
        const auto traced=continue_equilibrium(m,path);
        if(!traced.success || traced.states.empty())throw std::runtime_error("PH/PS target continuation: "+traced.message);
        const auto& last=traced.states.back();
        compare(h?"PH4/continuation-T":"PS4/continuation-T",*last.temperature,185,185,1e-8);
        compare(h?"PH4/continuation-split":"PS4/continuation-split",last.amounts[1],end_tp.amounts[1],1,1e-8);
        compare(h?"PH4/continuation-closure":"PS4/continuation-closure",total(m,last,h),total(m,end_tp,h),GAS_CONSTANT_J_PER_MOL_K*(h?185:1),1e-10);
    }
    const HelmholtzModel pure(propane_parameters(),{},records(true));
    auto vapor=equilibrium_support::flash(1);vapor.phases[0].kind=PhaseKind::Vapor;vapor.temperature.value=350;vapor.pressure.value=1e5;
    const auto ref=solve_equilibrium(pure,vapor);const auto q=condition(vapor,false,total(pure,ref,false),330);const auto c=compile_problem(pure,q);const auto r=solve_equilibrium(pure,q);
    auto pd=direction(c);pd.pressure_pa=1;
    const auto e=pure.evaluate(350,1/ref.molar_densities[0],{1},DerivativeOrder::Hessian);
    // alpha/cp identity: P_T = R/v - R F_V - RT F_TV, alpha = P_T/(rho P_rho).
    const auto state=pure.state(350,ref.molar_densities[0],{1});
    const double pt=GAS_CONSTANT_J_PER_MOL_K*(ref.molar_densities[0]-e.gradient[1]-350*e.hessian[3]);
    const double expected=350*pt/(ref.molar_densities[0]*ref.molar_densities[0]*state.pressure_density_derivative* *state.total_isobaric_heat_capacity);
    compare("PS4/isentropic",action(c,r,{pd}),expected,std::abs(expected),1e-8);
    auto sat=vapor;sat.phases={equilibrium_support::phase("L",PhaseKind::Liquid),equilibrium_support::phase("V",PhaseKind::Vapor,AmountSpec::pinned(0))};
    sat.pressure.value=1e6;sat.temperature=FreeValue::variable(300);sat.feed.values={2.5};sat.solver_options.values[epcsaft::policy::row] = 1e-11;
    auto sr=solve_equilibrium(pure,sat);const auto sc=compile_problem(pure,sat);auto sp=direction(sc);sp.pressure_pa=1;
    const double tangent=action(sc,sr,{sp});sr.amounts={1.625,0.875};sat.phases[1].amount={};
    const auto ps=condition(sat,false,total(pure,sr,false),300);const auto pc=compile_problem(pure,ps);const auto pr=solve_equilibrium(pure,ps);
    compare("PS4/saturation-tangent",action(pc,pr,{sp}),tangent,std::abs(tangent),1e-8);
}
// Assembled off-solution Jacobian and weighted-Hessian ladder (F5): centered row and weighted-gradient
// differences, the alternating weights and each condition row alone, objective factors 0, 0.5 and 1
// (only 0 for a feasibility cell, whose objective gradient is zero).
void curvature() {
    const HelmholtzModel m(properties_support::neutral_binary(),{},records());
    auto p=binary();const auto reference=solve_equilibrium(m,p);CHECK(reference.has_state());p.feed.values={0.27,0.73};
    for(bool h:{true,false}) {
        auto q=condition(p,h,total(m,reference,h),180);q.pressure.value*=1.03;const auto c=compile_problem(m,q);
        auto y=c.generic_point();const auto prepared=c.prepare(y);
        y[c.primal_size-1]*=1.02;
        for(const auto& l:prepared.shape)y[l.volume]+=0.01;
        y[prepared.shape[0].x[0]]+=0.02;
        for(std::size_t j=c.primal_size;j<y.size();++j)y[j]=0.11+0.13*j;
        std::vector<double> scale(y.size(),1);scale[c.primal_size-1]=y[c.primal_size-1];
        ladder(c,y,scale,h?"PH":"PS");
    }
}

void admission() {
    const auto params=propane_parameters();auto p=equilibrium_support::flash(1);p.temperature.value=250;p.pressure.value=1e6;
    const HelmholtzModel m(params,{},records(true));const auto ref=solve_equilibrium(m,p);CHECK(ref.has_state());
    for(int cause=0;cause<2;++cause) {
        const HelmholtzModel broken(params,{},cause?std::optional<ThermochemistryRecord>(records(true)):std::nullopt);
        auto q=condition(p,true,total(m,ref,true),cause?600:240);bool rejected=false;
        try{(void)compile_problem(broken,q);}catch(const std::invalid_argument& e){
            rejected=std::string(e.what()).find(cause?"outside interval":"missing thermochemistry")!=std::string::npos;}
        compare("admission/"+std::to_string(cause),rejected?1:0,1,1,1e-10);
        CHECK(broken.state(250,ref.molar_densities[0],{1}).pressure>0);
        CHECK(solve_equilibrium(broken,p).has_state());
    }
    // Saturation at 1 MPa (300.053 K, PH3) lies outside [100,300] K while the 300 K seed lies inside.
    auto narrow=records(true);narrow.component_intervals[0][0].upper=300;const HelmholtzModel outside(params,{},narrow);
    auto two=p;two.feed.values={2.5};two.phases.push_back(equilibrium_support::phase("V",PhaseKind::Vapor));
    bool refused=false;std::string cause;
    try { solve_equilibrium(outside,condition(two,true,0,300)); }
    catch(const std::invalid_argument& error) { refused=true;cause=error.what(); }
    CHECK(refused);
    compare("admission/saturation-outside",cause.find("outside interval: species 0 at 300.05")!=std::string::npos?1:0,1,1,1e-10,recording_support::field("cause",cause));
    ThermochemistryRecord joined;joined.standard_pressure=100000;
    joined.component_intervals={{{100,250,true,true,{IdealPolynomial{{40},200,0,0},100000}},
        {250,550,false,true,{IdealPolynomial{{50},300,100,1},100000}}}};
    const HelmholtzModel discontinuous(params,{},joined);const auto tp=solve_equilibrium(discontinuous,p);CHECK(tp.has_state());
    const auto q=condition(p,true,total(discontinuous,tp,true),240);const auto r=solve_equilibrium(discontinuous,q);CHECK(r.has_state());
    const auto c=compile_problem(discontinuous,q);auto d=direction(c);d.enthalpy_j=1;
    // Use the exact owned join coordinates so proximity cannot change its availability meaning.
    auto at=r.coordinates;const auto layout=c.prepare(at).shape[0];
    at[layout.x[0]]=0;at[*layout.amount]=1;at[layout.volume]=std::log(discontinuous.packing_fraction_coefficient(250,{1})*tp.molar_densities[0]);at[c.primal_size-1]=250;
    for(std::size_t j=0;j<c.rank;++j)at[c.primal_size+j]=tp.multipliers[j];
    auto join=r;join.coordinates=at;const auto evaluated=c.evaluate(at);join.residuals=evaluated.values;
    const auto actions=solved_state_actions(c,join,{{{SolvedStateObservableKind::Temperature}},{d}});
    compare("admission/join",actions.results[0].status==SolvedStateActionStatus::ThermalBoundary?1:0,1,1,1e-10);
    CHECK(actions.results[0].value.has_value());
}

inline double volume(const EquilibriumResult& r) {
    double v=0;for(std::size_t p=0;p<r.phase_count;++p)v+=r.amounts[p]/r.molar_densities[p];return v;
}
// Paterson (2023), Table 2: conserved total volume, pressure as its solved multiplier.
void volume_specification() {
    const HelmholtzModel model(properties_support::neutral_binary());
    const std::array<std::array<double,6>,3> references = {{{180,0.00047042881945383514,1e6,0.2945375955329262,0.9244914722077252,0.32615467905619605},
        {190,0.0006523810823353041,1e6,0.21537559277082224,0.8647803491384229,0.43828506711469967},
        {220,0.00044242974279205897,2e6,0.2160681158366867,0.7337137610815471,0.5485062740728838}}};
    for(const auto& ref:references) {
        auto p=binary();p.temperature=FreeValue::fixed(ref[0]);p.pressure=FreeValue::variable(ref[2]);p.volume_m3=ref[1];
        auto r=solve_equilibrium(model,p);
        if(!r.has_state())throw std::runtime_error("TV: "+r.message);
        const std::string id="TV/"+std::to_string(ref[0]);
        compare(id+"/P",*r.pressure,ref[2],ref[2],1e-8);
        compare(id+"/x",r.mole_fractions[0],ref[3],ref[3],1e-8);
        compare(id+"/y",r.mole_fractions[2],ref[4],ref[4],1e-8);
        compare(id+"/N-vapor",r.amounts[1],ref[5],1,1e-8);
        compare(id+"/V",volume(r),ref[1],ref[1],1e-10);
        for(std::size_t i=0;i<2;++i)compare(id+"/n"+std::to_string(i),r.amounts[0]*r.mole_fractions[i]+r.amounts[1]*r.mole_fractions[2+i],0.5,0.5,1e-10);
        CHECK(r.reduced_hessian_minimum>0);
        recording_support::observe((id+"/curvature").c_str(),"reduced-Helmholtz-curvature","chart",r.reduced_hessian_minimum,std::nullopt,"","Observation","exact-Lagrangian-Hessian","positive on the total-volume tangent","volume-specification");
        const auto c=compile_problem(model,p);auto d=direction(c);d.volume_m3=ref[1];
        std::vector<SolvedStateObservable> o={{SolvedStateObservableKind::Pressure},{SolvedStateObservableKind::PhaseAmount,1},
            {SolvedStateObservableKind::PhaseMoleFraction,0,0},{SolvedStateObservableKind::PhaseMoleFraction,1,0},
            {SolvedStateObservableKind::PhaseMolarDensity,0},{SolvedStateObservableKind::PhaseMolarDensity,1}};
        const auto a=solved_state_actions(c,r,{o,{d}});
        for(double h:{1e-4,5e-5}) {
            auto plus=p,minus=p;plus.volume_m3=ref[1]*(1+h);minus.volume_m3=ref[1]*(1-h);
            const auto rp=solve_equilibrium(model,plus),rm=solve_equilibrium(model,minus);
            CHECK(rp.has_state()&&rm.has_state());
            const std::array<double,6> vp={*rp.pressure,rp.amounts[1],rp.mole_fractions[0],rp.mole_fractions[2],rp.molar_densities[0],rp.molar_densities[1]};
            const std::array<double,6> vm={*rm.pressure,rm.amounts[1],rm.mole_fractions[0],rm.mole_fractions[2],rm.molar_densities[0],rm.molar_densities[1]};
            for(std::size_t j=0;j<a.results.size();++j){CHECK(a.results[j].action);const double fd=(vp[j]-vm[j])/(2*h);
                compare(id+"/action/"+std::to_string(j)+"/"+std::to_string(h),*a.results[j].action,fd,std::abs(fd),1e-6);}
        }
        if(ref[0]==180) {
            const auto tp=solve_equilibrium(model,binary());CHECK(tp.has_state());
            auto inverse=p;inverse.volume_m3=volume(tp);const auto back=solve_equilibrium(model,inverse);CHECK(back.has_state());
            compare("TV/TP-roundtrip/P",*back.pressure,*tp.pressure,*tp.pressure,1e-8);
            for(std::size_t i=0;i<tp.mole_fractions.size();++i)compare("TV/TP-roundtrip/x"+std::to_string(i),back.mole_fractions[i],tp.mole_fractions[i],1,1e-8);
            // The total-volume row and its optimizer-weighted Hessian away from the root.
            auto y=r.coordinates;y[0]+=0.01;auto scale=std::vector<double>(y.size(),1.0);scale[c.primal_size-1]=ref[2];
            ladder(c,y,scale,"TV/curvature",{0.});
            auto twice=p;twice.feed.values={1,1};twice.volume_m3=2*ref[1];const auto rr=solve_equilibrium(model,twice);CHECK(rr.has_state());
            compare("TV/extensive/P",*rr.pressure,*r.pressure,*r.pressure,1e-10);
            for(std::size_t i=0;i<2;++i)compare("TV/extensive/N"+std::to_string(i),rr.amounts[i],2*r.amounts[i],2,1e-10);
        }
    }
    auto p=binary();p.temperature=FreeValue::fixed(220);p.volume_m3=0.002;p.pressure=FreeValue::variable();p.phases.erase(p.phases.begin());
    // Clapeyron 0.6.29, matched Gross (2001) parameters, one mole and rho=500 mol/m3.
    const auto r=solve_equilibrium(model,p);CHECK(r.has_state());compare("TV/homogeneous/P",*r.pressure,830502.7179639882,830502.7179639882,1e-8);
    const auto c=compile_problem(model,p);auto d=direction(c);d.volume_m3=*p.volume_m3;
    const auto state=model.state(220,500,{0.5,0.5});
    compare("TV/homogeneous/dP",action(c,r,{d},SolvedStateObservableKind::Pressure),-500*state.pressure_density_derivative,500*state.pressure_density_derivative,1e-10);
    for(bool dew:{false,true}) {
        const double vb=5.2524821398533176e-5,vd=0.00557698412667614;
        auto q=binary();q.temperature=FreeValue::fixed(190);q.pressure=FreeValue::variable(dew?270991.0143801206:2095823.2840843329);
        q.volume_m3=dew?vd:vb;q.phases={q.phases[dew?1:0]};
        const auto endpoint=solve_equilibrium(model,q);CHECK(endpoint.has_state());
        compare(dew?"TV/dew/P":"TV/bubble/P",*endpoint.pressure,q.pressure.value,q.pressure.value,1e-8);
        double previous=1.0;
        for(double f:{1e-3,1e-4}) {
            auto near=binary();near.temperature=q.temperature;near.pressure=FreeValue::variable(q.pressure.value*(dew?1.01:0.99));near.volume_m3=dew?vd-f*(vd-vb):vb+f*(vd-vb);
            const auto rr=solve_equilibrium(model,near);if(!rr.has_state())throw std::runtime_error("TV boundary: "+rr.message);
            compare("TV/boundary/V",volume(rr),*near.volume_m3,*near.volume_m3,1e-10);
            const double absent=rr.amounts[dew?0:1];CHECK(absent>0&&absent<previous);previous=absent;
            recording_support::observe(dew?"TV/dew/limit":"TV/bubble/limit","absent-phase-amount","mol",absent,std::nullopt,
                recording_support::fields({recording_support::field("fractional_volume_distance",f),recording_support::field("P_Pa",*rr.pressure),recording_support::field("V_m3",*near.volume_m3),recording_support::field("mole_fractions",rr.mole_fractions)}),
                "Observation","volume-limit","positive amount decreases toward the phase boundary","volume-specification");
        }
    }

    const HelmholtzModel pure(propane_parameters());
    auto sat=equilibrium_support::flash(1);sat.temperature=FreeValue::fixed(300);sat.pressure=FreeValue::variable(1e6);
    sat.phases.push_back(equilibrium_support::phase("V",PhaseKind::Vapor,AmountSpec::pinned(0)));
    const auto s=solve_equilibrium(pure,sat);CHECK(s.has_state());
    sat.phases[1].amount=AmountSpec::free();sat.volume_m3=0.65/s.molar_densities[0]+0.35/s.molar_densities[1];
    const auto tv=solve_equilibrium(pure,sat);if(!tv.has_state())throw std::runtime_error("pure TV: "+tv.message);
    compare("TV/pure/P",*tv.pressure,*s.pressure,*s.pressure,1e-8);
    compare("TV/pure/split",tv.amounts[1],0.35,1,1e-8);CHECK(tv.reduced_hessian_minimum>0);
    recording_support::observe("TV/pure/curvature","reduced-Helmholtz-curvature","chart",tv.reduced_hessian_minimum,std::nullopt,"","Observation","exact-Lagrangian-Hessian","pure coexistence retains volume and both free amounts","volume-specification");
}
void internal_energy_specification() {
    const double R=GAS_CONSTANT_J_PER_MOL_K;
    ThermochemistryRecord thermal; thermal.standard_pressure=1e5;
    const double s0=R*(2.5+2.5*std::log(298.15)-std::log(1e5/R));
    for(int i=0;i<2;++i)thermal.component_intervals.push_back({{100,600,false,false,{IdealPolynomial{{2.5*R},298.15,2.5*R*298.15,s0},1e5}}});
    const HelmholtzModel model(properties_support::neutral_binary(),{},thermal);
    const auto energy=[](const HelmholtzModel& m,const EquilibriumResult& r) {
        double value=0;const auto c=m.parameters().component_count();
        for(std::size_t p=0;p<r.phase_count;++p) {
            const std::vector<double> x(r.mole_fractions.begin()+p*c,r.mole_fractions.begin()+(p+1)*c);
            value+=r.amounts[p]* *m.state(*r.temperature,r.molar_densities[p],x).total_internal_energy;
        }
        return value;
    };
    const auto recover=[&](const HelmholtzModel& m,EquilibriumProblem p,EquilibriumResult ref,double start,const std::string& id) {
        const double u=energy(m,ref),v=volume(ref),n=std::accumulate(p.feed.values.begin(),p.feed.values.end(),0.0);
        p.temperature=FreeValue::variable(start);p.pressure=FreeValue::variable(1.05* *ref.pressure);
        p.internal_energy_j=u;p.volume_m3=v;
        const auto r=solve_equilibrium(m,p);
        if(!r.has_state())throw std::runtime_error(id+": "+r.message);
        compare(id+"/T",*r.temperature,*ref.temperature,*ref.temperature,1e-8);
        compare(id+"/P",*r.pressure,*ref.pressure,*ref.pressure,1e-8);
        compare(id+"/U",energy(m,r),u,n*R* *r.temperature,1e-10);
        compare(id+"/V",volume(r),v,v,1e-10);
        for(std::size_t j=0;j<r.amounts.size();++j)compare(id+"/N"+std::to_string(j),r.amounts[j],ref.amounts[j],ref.amounts[j],1e-8);
        for(std::size_t j=0;j<r.mole_fractions.size();++j)compare(id+"/x"+std::to_string(j),r.mole_fractions[j],ref.mole_fractions[j],ref.mole_fractions[j],1e-8);
        const auto c=m.parameters().component_count();
        for(std::size_t i=0;i<c;++i){double total=0;for(std::size_t q=0;q<r.phase_count;++q)total+=r.amounts[q]*r.mole_fractions[q*c+i];compare(id+"/balance"+std::to_string(i),total,p.feed.values[i],n,1e-10);}
        const double entropy=-total(m,r,false)/R;
        compare(id+"/negative-entropy",r.objective,entropy,std::max(1.0,std::abs(entropy)),1e-8);
        for(std::size_t j=0;j<r.molar_densities.size();++j)compare(id+"/rho"+std::to_string(j),r.molar_densities[j],ref.molar_densities[j],ref.molar_densities[j],1e-8);
        CHECK(r.requested_tolerance_met.value_or(false));CHECK(r.reduced_hessian_minimum>0);
        return std::make_pair(p,r);
    };
    const std::array<std::array<double,4>,3> refs={{{180,1e6,-5269.669334777402,0.00047042881945383514},{190,1e6,-4108.403820583579,0.0006523810823353041},{220,2e6,-2131.692379591598,0.00044242974279205897}}};
    for(const auto& row:refs) {
        auto p=binary();p.temperature.value=row[0];p.pressure.value=row[1];
        const auto tp=solve_equilibrium(model,p);CHECK(tp.has_state());
        compare("UV/reference/U",energy(model,tp),row[2],std::abs(row[2]),1e-8);
        compare("UV/reference/V",volume(tp),row[3],row[3],1e-8);
        const auto solved=recover(model,p,tp,1.02*row[0],"UV/"+std::to_string(row[0]));
        if(row[0]==180) {
            auto doubled=p;doubled.feed.values={1,1};auto ref=tp;for(auto& n:ref.amounts)n*=2;
            recover(model,doubled,ref,1.02*row[0],"UV/extensive");
            const auto c=compile_problem(model,solved.first);auto y=solved.second.coordinates;y[0]+=0.01;
            std::vector<double> scales(y.size(),1);scales[c.primal_size-2]=180;scales[c.primal_size-1]=1e6;
            ladder(c,y,scales,"UV/curvature",{0.},1e-6);
            // In physical amounts and volumes, L_Tz=0, L_TT=Cv/(R T²), L_zz=A_zz/(RT).
            const auto at=c.prepare(solved.second.coordinates);const auto h=c.entropy_hessian_at(at);
            const std::size_t width=at.point.size();std::vector<double> tangent(width,0.0);const double dt=1.8;
            tangent[c.primal_size-2]=dt;double expected=0.0;
            for(std::size_t phase=0;phase<at.shape.size();++phase) {
                const auto& shape=at.shape[phase];const double amount=at.amount[phase],vol=amount*at.volume[phase];
                std::vector<double> amounts(2);for(std::size_t i=0;i<2;++i)amounts[i]=amount*at.composition[phase][i];
                const double dn[2]={0.03,-0.01},dtotal=dn[0]+dn[1],dv=0.02*vol;
                const auto derivative=model.evaluate(180,vol,amounts,DerivativeOrder::Hessian);
                const std::vector<double> move={0,dv,dn[0],dn[1]};
                for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j)expected+=move[i]*derivative.hessian[i*4+j]*move[j];
                expected+=dn[0]*dn[0]/amounts[0]+dn[1]*dn[1]/amounts[1]-2*dv*dtotal/vol+amount*dv*dv/(vol*vol);
                const auto state=model.state(180,1/at.volume[phase],at.composition[phase]);
                expected+=amount* *state.total_isochoric_heat_capacity/(R*180*180)*dt*dt;
                for(std::size_t i=0;i<2;++i)tangent[shape.x[i]]=dn[i]/amounts[i]-dtotal/amount;
                tangent[*shape.amount]=dtotal;
                tangent[shape.volume]=(at.packing[phase].component[0]*dn[0]/amounts[0]+at.packing[phase].component[1]*dn[1]/amounts[1])/at.packing[phase].coefficient
                    +at.packing[phase].temperature_derivative/at.packing[phase].coefficient*dt-dv/vol;
            }
            double actual=0.0;for(std::size_t i=0;i<width;++i)for(std::size_t j=0;j<width;++j)actual+=tangent[i]*h[i*width+j]*tangent[j];
            compare("UV/entropy-physical-curvature",actual,expected,std::abs(expected),1e-8);
            auto du=direction(c),dv_control=direction(c);du.internal_energy_j=R*180;dv_control.volume_m3=*solved.first.volume_m3;
            std::vector<SolvedStateObservable> observables;
            for(std::size_t phase=0;phase<2;++phase) {
                observables.push_back({SolvedStateObservableKind::PhaseAmount,phase});
                observables.push_back({SolvedStateObservableKind::PhaseVolume,phase});
                observables.push_back({SolvedStateObservableKind::PhaseProperty,phase,0,PropertyObservable::TotalInternalEnergy});
                for(std::size_t i=0;i<2;++i)observables.push_back({SolvedStateObservableKind::PhaseMoleFraction,phase,i});
            }
            const auto actions=solved_state_first_actions(c,solved.second,observables,{du,dv_control});
            for(std::size_t d=0;d<2;++d) {
                double energy_action=0,volume_action=0;std::array<double,2> inventory{};
                for(std::size_t phase=0;phase<2;++phase) {
                    const auto& a=actions[d].results;const std::size_t offset=5*phase;
                    for(std::size_t i=0;i<5;++i)CHECK(a[offset+i].action.has_value());
                    energy_action+=*a[offset].action* *a[offset+2].value+*a[offset].value* *a[offset+2].action;
                    volume_action+=*a[offset+1].action;
                    for(std::size_t i=0;i<2;++i)inventory[i]+=*a[offset].action* *a[offset+3+i].value+*a[offset].value* *a[offset+3+i].action;
                }
                compare("UV/pair-action/U"+std::to_string(d),energy_action,d==0?R*180:0,R*180,1e-8);
                compare("UV/pair-action/V"+std::to_string(d),volume_action,d==1?*solved.first.volume_m3:0,*solved.first.volume_m3,1e-8);
                for(std::size_t i=0;i<2;++i)compare("UV/pair-action/inventory"+std::to_string(d)+std::to_string(i),inventory[i],0,1,1e-8);
            }


        }
    }
    auto p=binary();p.phases.erase(p.phases.begin());p.temperature.value=300;p.pressure.value=1e5;
    const auto tp=solve_equilibrium(model,p);CHECK(tp.has_state());
    const auto solved=recover(model,p,tp,294,"UV/homogeneous");
    compare("UV/reference/homogeneous-U",energy(model,tp),3720.0048586268294,3720.0048586268294,1e-8);
    const auto c=compile_problem(model,solved.first);auto u=direction(c),v=direction(c);u.internal_energy_j=1;v.volume_m3=*solved.first.volume_m3;
    const auto state=model.state(300,tp.molar_densities[0],{0.5,0.5});const auto eos=model.evaluate(300,1/tp.molar_densities[0],{0.5,0.5},DerivativeOrder::Hessian);
    const double cv=*state.total_isochoric_heat_capacity,pt=R*(tp.molar_densities[0]-eos.gradient[1]-300*eos.hessian[1]);
    compare("UV/action/TU",action(c,solved.second,{u}),1/cv,1/cv,1e-8);
    compare("UV/action/PU",action(c,solved.second,{u},SolvedStateObservableKind::Pressure),pt/cv,std::abs(pt/cv),1e-8);
    const double tv=(1e5-300*pt)/cv* *solved.first.volume_m3;
    compare("UV/action/TV",action(c,solved.second,{v}),tv,std::abs(tv),1e-8);
    compare("UV/action/PV",action(c,solved.second,{v},SolvedStateObservableKind::Pressure),pt*tv-tp.molar_densities[0]*state.pressure_density_derivative,1e5,1e-8);
    thermal.component_intervals.resize(1);const HelmholtzModel pure(propane_parameters(),{},thermal);
    auto sat=equilibrium_support::flash(1);sat.temperature=FreeValue::fixed(300);sat.pressure=FreeValue::variable(1e6);
    sat.phases.push_back(equilibrium_support::phase("V",PhaseKind::Vapor,AmountSpec::pinned(0)));
    auto sr=solve_equilibrium(pure,sat);CHECK(sr.has_state());sr.amounts={0.6,0.4};sat.phases[1].amount=AmountSpec::free();
    recover(pure,sat,sr,306,"UV/propane");
    recover(pure,sat,sr,350,"UV/propane-exterior-seed");
    auto edge=thermal;edge.component_intervals[0][0].upper=300;edge.component_intervals[0][0].upper_inclusive=true;
    const HelmholtzModel bounded(propane_parameters(),{},edge);
    auto q=sat;q.temperature=FreeValue::variable(300);q.internal_energy_j=energy(pure,sr);q.volume_m3=volume(sr);
    bool refused=false;try{(void)solve_equilibrium(bounded,q);}catch(const std::invalid_argument& e){refused=std::string(e.what()).find("interval endpoint")!=std::string::npos;}CHECK(refused);
}
// Apparent Cp under PH (free T, fixed P, enthalpy target) on the PH5/2 state, read from solved_state_actions of the
// PH-compiled problem, against central differences of the driven equilibrium enthalpy along T (F1) and of the driven
// apparent Cp at the re-solved T over epsilon/k at fixed H, P and feed (F2); each extrapolated ladder within its criterion.
void apparent_heat_capacity() {
    const double R=GAS_CONSTANT_J_PER_MOL_K;const auto reaction=ethane_reaction(2);const auto& r=reaction.first;const auto& tp=reaction.second;
    const std::vector<ActiveParameter> active={{ParameterFamily::DispersionEnergyOverK,{0}}};
    const auto model=[&](double shift){auto q=e2::identical_ethane_pair();q.dispersion_energy_over_k[0]+=shift;return HelmholtzModel(q,active,r);};
    const auto m=model(0);
    const auto read=[](const CompiledSystem& c,const EquilibriumResult& x,SolvedStateObservableKind kind,std::vector<SolvedStateActionDirection> d) {
        const auto item=solved_state_actions(c,x,{{{kind}},d}).results[0];
        if(!item.value) throw std::runtime_error("apparent Cp read: "+item.diagnostic_message.value_or("unavailable"));
        return item;
    };
    const auto driven=[&](const HelmholtzModel& model,double T,SolvedStateObservableKind kind) {
        auto q=tp;q.temperature=FreeValue::fixed(T);const auto c=compile_problem(model,q);const auto x=solve_equilibrium(model,q);CHECK(x.has_state());
        return *read(c,x,kind,{direction(c)}).value;
    };
    constexpr auto Cp=SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity,H=SolvedStateObservableKind::SystemEnthalpy;
    const auto ph=condition(tp,true,driven(m,300,H),280);
    const auto c=compile_problem(m,ph);const auto x=solve_equilibrium(m,ph);CHECK(x.has_state());
    const double T=*x.temperature,feed=c.feed_total;auto wq=direction(c);wq.active_parameters[0]=1;
    const double cp=*read(c,x,Cp,{direction(c)}).value;
    const auto parameter=read(c,x,Cp,{wq});CHECK(parameter.status==SolvedStateActionStatus::Available&&parameter.action);
    const double dcp=*parameter.action;
    compare("PH6/T",T,300,300,1e-8);
    // Cross-checks: the driven-T/P apparent Cp at this T, and c_p,frozen + (h_B-h_A)^2 K/((1+K)^2 R T^2) for identical species.
    compare("PH6/driven-cp",cp,driven(m,T,Cp),cp,1e-12);
    const double K=std::exp(std::log(2.)-5000/(R*T));SolvedStateObservable frozen{SolvedStateObservableKind::PhaseProperty};
    frozen.property=PropertyObservable::TotalIsobaricHeatCapacity;
    const double analytic=*solved_state_actions(c,x,{{frozen},{direction(c)}}).results[0].value+25e6*K/((1+K)*(1+K))/(R*T*T);
    compare("PH6/analytic-cp",cp,analytic,cp,1e-12);
    const auto ladder=[&](const std::string& id,double exact,double h0,double tolerance,const auto& difference) {
        double fd[3];
        for(int i=0;i<3;++i) {
            fd[i]=difference(h0*std::ldexp(1.,-i));
            recording_support::observe((id+"/step-"+std::to_string(i)).c_str(),"apparent-cp-difference","J/(mol-feed K)",exact,fd[i],
                recording_support::fields({recording_support::field("step",h0*std::ldexp(1.,-i)),recording_support::field("error",std::abs(exact-fd[i])/std::abs(exact))}),
                "Observation","centered-difference","fixed-H apparent Cp","caloric-specification");
        }
        for(int i=0;i<2;++i) compare(id+"/ratio-"+std::to_string(i),std::abs(exact-fd[i])/std::abs(exact-fd[i+1]),4,1,0.5);
        compare(id+"/Richardson",exact,(4*fd[2]-fd[1])/3,exact,tolerance);
    };
    ladder("PH6/F1",cp,0.4,1e-10,[&](double h){return (driven(m,T+h,H)-driven(m,T-h,H))/(2*h*feed);});
    ladder("PH6/F2",dcp,0.2,1e-9,[&](double d){
        double v[2];
        for(int side=0;side<2;++side) {
            const auto shifted=model(side?d:-d);const auto xs=solve_equilibrium_from_coordinates(shifted,ph,x.coordinates);CHECK(xs.has_state());
            v[side]=driven(shifted,*xs.temperature,Cp);
        }
        return (v[1]-v[0])/(2*d);
    });
    bool refused=false;try{(void)solved_state_actions(c,x,{{{Cp}},{wq,wq}});}
    catch(const std::invalid_argument& e){refused=std::string(e.what()).find("fixed-P enthalpy target")!=std::string::npos;}
    CHECK(refused);
}
} // namespace caloric_specification_checks
