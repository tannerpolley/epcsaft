// Issue 224: manufactured calorics on the Gross (2001) methane/ethane and propane packets.
#include "equilibrium_checks_shared.hpp"
#include "caloric_specification_inputs.hpp"
#include "caloric_specification_shared.hpp"
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
    for(bool h:{true,false}) {
        recovered(pure,q,ref,h,300,h?"PH3":"PS3");
        const auto result=solve_equilibrium(pure,condition(q,h,total(pure,ref,h),300));
        for(std::size_t i=0;i<2;++i) {
            const auto state=pure.state(*result.temperature,result.molar_densities[i],{1});
            CHECK(state.pressure_density_derivative>0);
            const auto liquid=pure.state(*ref.temperature,ref.molar_densities[0],{1});
            const auto vapor=pure.state(*ref.temperature,ref.molar_densities[1],{1});
            const double hl=*liquid.total_enthalpy,hv=*vapor.total_enthalpy,sl=*liquid.total_entropy,sv=*vapor.total_entropy;
            compare(std::string(h?"PH3":"PS3")+"/H-lever",(total(pure,result,true)/2.5-hl)/(hv-hl),0.35,1,1e-8);
            compare(std::string(h?"PH3":"PS3")+"/S-lever",(total(pure,result,false)/2.5-sl)/(sv-sl),0.35,1,1e-8);
        }
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
void reactions() {
    for(int kind=0;kind<3;++kind) {
        auto r=records();const double R=GAS_CONSTANT_J_PER_MOL_K;
        for(auto& species:r.component_intervals) species[0].correlation={IdealPolynomial{{40},298.15,0,0},100000};
        r.component_intervals[1][0].correlation={IdealPolynomial{{40},298.15,kind==1?0.:5000.,kind==0?0.:R*std::log(2.)},100000};
        auto p=equilibrium_support::flash(2);p.temperature=FreeValue::fixed(300);p.pressure=FreeValue::fixed(1e6);
        p.phases[0].kind=PhaseKind::Vapor;  // vapor-side sole root at 300 K, 1 MPa (#423)
        p.feed.values={2.5,0};p.solver_options.values[epcsaft::policy::row] = 1e-11;
        ReactionSpec rx;rx.name="A=B";rx.stoichiometry={-1,1};
        if(kind==0) rx.affinity_j_per_mol=5000;else if(kind==1) rx.equilibrium_constant=2;
        else rx.correlation=ReactionLogPolynomial{std::log(2.),-5000/R,0,0,298.15,100,550,EOS_STANDARD_STATE_ID};
        p.reactions={rx};const HelmholtzModel m(e2::identical_ethane_pair(),{},r);
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
    std::string cause;try{cause=solve_equilibrium(outside,condition(two,true,0,300)).message;}catch(const std::exception& e){cause=e.what();}
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

} // namespace caloric_specification_checks
