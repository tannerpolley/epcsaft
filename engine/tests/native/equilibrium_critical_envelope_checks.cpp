// Issue 224 B: critical-equation candidates and phase envelopes on the Gross (2001) propane and
// methane/ethane packets; numerical verification only.
#include "equilibrium_checks_shared.hpp"
#include "caloric_specification_inputs.hpp"
#include "caloric_specification_shared.hpp"
namespace critical_envelope_checks {
using namespace epcsaft;
using namespace epcsaft::equilibrium;
using caloric_specification_checks::compare;
inline EquilibriumProblem critical_cell(const std::vector<double>& z, double t, double rho) {
    auto p=equilibrium_support::flash(z.size());p.feed.values=z;p.solver_options.values[epcsaft::policy::row] = 1e-11;
    p.temperature=FreeValue::variable(t);p.pressure=FreeValue::variable(1e6);p.critical=CriticalDeclaration{rho,{}};
    return p;
}
inline EquilibriumResult solved(const HelmholtzModel& m,const EquilibriumProblem& p,const std::string& id) {
    const auto r=solve_equilibrium(m,p);
    if(!r.has_state()||!r.critical) throw std::runtime_error(id+": "+r.message);
    const auto& k=*r.critical;
    recording_support::observe((id+"/candidate").c_str(),"critical-equation-candidate","K",*r.temperature,std::nullopt,
        recording_support::fields({recording_support::field("P_Pa",*r.pressure),recording_support::field("rho_mol_m3",r.molar_densities[0]),
            recording_support::field("eigenvalues",k.eigenvalues),recording_support::field("zero_band",k.zero_band),recording_support::field("cubic",k.cubic),
            recording_support::field("condition",k.condition_estimate),recording_support::field("backward_error",k.backward_error)}),
        "Observation","root identity and spectrum of M","candidate only, no stability certificate","critical-envelope");
    return r;
}
inline EquilibriumProblem envelope_cell(const std::vector<double>& z,double t,double pressure,EnvelopeCondition condition) {
    auto p=equilibrium_support::flash(z.size());p.feed.values=z;p.solver_options.values[epcsaft::policy::row] = 1e-11;
    p.phases.push_back(equilibrium_support::phase("V",PhaseKind::Vapor,AmountSpec::pinned(0)));
    p.temperature=FreeValue::variable(t);p.pressure=FreeValue::variable(pressure);p.envelope=condition;
    return p;
}
inline ContinuationResult trace(const HelmholtzModel& m,const EquilibriumProblem& p,double step,double straddle) {
    ContinuationRequest r{p,p,std::nullopt,{1,1e-7,4000,6},{}};r.options.envelope_max_log_step=step;r.options.envelope_straddle=straddle;
    return continue_equilibrium(m,r);
}
inline void same_state(const std::string& id,const EquilibriumResult& a,const EquilibriumResult& b) {
    compare(id+"/T",*a.temperature,*b.temperature,*b.temperature,1e-8);
    compare(id+"/P",*a.pressure,*b.pressure,*b.pressure,1e-8);
    for(std::size_t i=0;i<a.molar_densities.size();++i)compare(id+"/rho-"+std::to_string(i),a.molar_densities[i],b.molar_densities[i],b.molar_densities[i],1e-8);
    for(std::size_t i=0;i<a.mole_fractions.size();++i)compare(id+"/x-"+std::to_string(i),a.mole_fractions[i],b.mole_fractions[i],1,1e-8);
}
// C1/C2: the pure propane candidate is admitted directly; pressure-closed cp is unavailable while T, P and
// density survive; P_rho/(RT) and rho P_rhorho/(RT) from pressure_state_hessian vanish within a derived floor.
void pure_critical() {
    const HelmholtzModel m(caloric_specification_checks::propane_parameters());
    const auto p=critical_cell({1},375,5000);const auto r=solved(m,p,"C1");const auto c=compile_problem(m,p);
    const auto actions=solved_state_actions(c,r,{{{SolvedStateObservableKind::Temperature},{SolvedStateObservableKind::PhaseMolarDensity},
        {SolvedStateObservableKind::PhaseProperty,0,0,PropertyObservable::TotalIsobaricHeatCapacity}},{{0,0,{0.1},{}}}});
    compare("C1/direct-values",actions.results[0].action&&actions.results[1].action&&actions.results[1].value?1:0,1,1,1e-10);
    compare("C1/pressure-closed-cp",actions.results[2].status==SolvedStateActionStatus::PressureSingular&&!actions.results[2].value?1:0,1,1,1e-10);
    // Floor: the exact Newton correction of the critical system moves each identity, plus rounding of its terms.
    auto y=r.coordinates;const auto rows=c.evaluate(y).values;const auto J=c.jacobian(y);const auto n=static_cast<Eigen::Index>(y.size());
    const Eigen::VectorXd delta=Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(J.data(),n,n)
        .fullPivLu().solve(-Eigen::Map<const Eigen::VectorXd>(rows.data(),n));
    const auto identities=[&](const std::vector<double>& at) {
        const auto prepared=c.prepare(at);const double T=prepared.temperature,rho=prepared.composition[0][0]/prepared.volume[0];
        const auto P=m.neutral_property(T,rho,NeutralProperty::Pressure,DerivativeOrder::Hessian);
        const auto a=m.neutral_property(T,rho,NeutralProperty::ResidualHelmholtz,DerivativeOrder::Hessian);
        const double rt=GAS_CONSTANT_J_PER_MOL_K*T,i1=P.gradient[1]/rt,i2=rho*P.hessian[3]/rt,t1=2*rho*a.gradient[1],t2=rho*rho*a.hessian[3];
        return std::array<double,4>{i1,i2,8*std::numeric_limits<double>::epsilon()*(1+std::abs(t1)+std::abs(t2)),
            8*std::numeric_limits<double>::epsilon()*(std::abs(t1)+4*std::abs(t2)+std::abs(i2-t1-4*t2))};
    };
    const auto at=identities(y);for(Eigen::Index i=0;i<n;++i)y[i]+=delta(i);const auto corrected=identities(y);
    for(std::size_t k=0;k<2;++k) {
        const double floor=std::abs(corrected[k]-at[k])+at[2+k];const std::string id=k?"C2/rho_P_rhorho":"C2/P_rho";
        recording_support::observe((id+"/floor").c_str(),"identity-floor","1",floor,std::nullopt,
            recording_support::fields({recording_support::field("newton_change",std::abs(corrected[k]-at[k])),recording_support::field("rounding",at[2+k])}),
            "Observation","propagated floor","Newton correction of the critical system plus term rounding","critical-envelope");
        CHECK(floor<1e-9);
        compare(id,at[k],0,1,1e-9);
    }
}
// C3: binary candidates and the first composition action along (0.1, -0.1) against Richardson-refined
// differences of central critical solves; the second action against differences of first actions; residual and total
// cp (P_rho > 0 off a pure fluid) equal the HelmholtzModel::state route.
void binary_critical() {
    const HelmholtzModel m(properties_support::neutral_binary(),{},caloric_specification_checks::records());
    for(const double z0:{0.5,0.27}) {
        const std::string id="C3/"+std::to_string(z0);
        const auto p=critical_cell({z0,1-z0},260,6000);const auto r=solved(m,p,id);const auto c=compile_problem(m,p);
        const SolvedStateActionDirection d{0,0,{0.1,-0.1},{}},d0{0,0,{0.1,0},{}},d1{0,0,{0,-0.1},{}};
        for(const auto kind:{SolvedStateObservableKind::Temperature,SolvedStateObservableKind::Pressure,SolvedStateObservableKind::PhaseMolarDensity}) {
            const auto action=[&](const SolvedStateActionDirection& v) {
                const auto a=solved_state_actions(c,r,{{{kind}},{v}}).results[0];
                if(!a.action)throw std::runtime_error(id+": critical first action: "+a.diagnostic_message.value_or(""));
                return *a.action;
            };
            const double exact=action(d),scale=std::abs(action(d0))+std::abs(action(d1));
            double best=INFINITY,previous=0;
            for(int rung=0;rung<5;++rung) {
                const double h=std::ldexp(0.02,-rung);double side[2];
                for(int s:{0,1}) {
                    auto q=p;q.feed.values={z0+(2*s-1)*h*0.1,1-z0-(2*s-1)*h*0.1};
                    const auto rr=solve_equilibrium_from_coordinates(m,q,r.coordinates);if(!rr.has_state())throw std::runtime_error(id+": "+rr.message);
                    side[s]=kind==SolvedStateObservableKind::Temperature?*rr.temperature:kind==SolvedStateObservableKind::Pressure?*rr.pressure:rr.molar_densities[0];
                }
                const double reference=(side[1]-side[0])/(2*h);
                if(rung)best=std::min(best,std::abs(exact-(4*reference-previous)/3)/scale);
                previous=reference;
            }
            compare(id+"/first-action/"+std::to_string(static_cast<int>(kind)),best,0,1,1e-6);
        }
        // A2 along (d, d) against Richardson-refined central differences of A1 at the perturbed critical solves.
        for(const auto kind:{SolvedStateObservableKind::Temperature,SolvedStateObservableKind::Pressure,SolvedStateObservableKind::PhaseMolarDensity}) {
            const auto second=solved_state_actions(c,r,{{{kind}},{d,d}}).results[0];
            if(!second.action)throw std::runtime_error(id+": critical second action: "+std::string(solved_state_status_name(second.status)));
            double best=INFINITY,previous=0;
            for(int rung=0;rung<5;++rung) {
                const double h=std::ldexp(0.02,-rung);double side[2];
                for(int s:{0,1}) {
                    auto q=p;q.feed.values={z0+(2*s-1)*h*0.1,1-z0-(2*s-1)*h*0.1};
                    const auto rr=solve_equilibrium_from_coordinates(m,q,r.coordinates);if(!rr.has_state())throw std::runtime_error(id+": "+rr.message);
                    const auto a=solved_state_actions(compile_problem(m,q),rr,{{{kind}},{d}}).results[0];
                    if(!a.action)throw std::runtime_error(id+": shifted first action: "+a.diagnostic_message.value_or(""));
                    side[s]=*a.action;
                }
                const double reference=(side[1]-side[0])/(2*h);
                if(rung)best=std::min(best,std::abs(*second.action-(4*reference-previous)/3)/std::abs(*second.action));
                previous=reference;
            }
            compare(id+"/second-action/"+std::to_string(static_cast<int>(kind)),best,0,1,1e-6);
        }
        const auto cp=solved_state_actions(c,r,{{{SolvedStateObservableKind::PhaseProperty,0,0,PropertyObservable::ResidualIsobaricHeatCapacity},
            {SolvedStateObservableKind::PhaseProperty,0,0,PropertyObservable::TotalIsobaricHeatCapacity}},{d}}).results;
        const auto state=m.state(*r.temperature,r.molar_densities[0],{r.mole_fractions[0],r.mole_fractions[1]});
        compare(id+"/cp-residual",cp[0].value.value_or(NAN),*state.residual_isobaric_heat_capacity,std::abs(*state.residual_isobaric_heat_capacity),1e-10,
            recording_support::field("P_rho_over_RT",state.pressure_density_derivative/(GAS_CONSTANT_J_PER_MOL_K * *r.temperature)));
        compare(id+"/cp-total",cp[1].value.value_or(NAN),*state.total_isobaric_heat_capacity,std::abs(*state.total_isobaric_heat_capacity),1e-10);
    }
}
// B-curvature: the shared F5 ladder on critical (nonassociating and associating) and envelope callbacks.
inline void off_solution(const HelmholtzModel& m,const EquilibriumProblem& p,const std::string& id) {
    const auto c=compile_problem(m,p);auto y=c.generic_point();const auto prepared=c.prepare(y);
    const std::size_t t=c.primal_size-2;y[t]*=1.02;y[t+1]*=1.03;
    for(const auto& l:prepared.shape)y[l.volume]+=0.01;
    y[prepared.shape[0].x[0]]+=0.02;
    for(std::size_t i=0;i<c.mode_count;++i)y[c.mode_offset()+i]*=1.1;
    for(std::size_t j=c.primal_size;j<y.size();++j)y[j]=0.11+0.13*j;
    std::vector<double> scale(y.size(),1);scale[t]=y[t];scale[t+1]=y[t+1];
    caloric_specification_checks::ladder(c,y,scale,id,{0.});
}
void curvature() {
    const HelmholtzModel m(properties_support::neutral_binary()),a(properties_support::associating_binary());
    off_solution(m,critical_cell({0.27,0.73},260,6000),"B/critical");
    off_solution(a,critical_cell({0.1385,0.8615},376.88,5000),"B/critical-associating");
    for(const auto kind:{EnvelopeCoordinate::LogTemperature,EnvelopeCoordinate::LogPressure,EnvelopeCoordinate::LogK})
        off_solution(m,envelope_cell({0.27,0.73},180,1e5,{kind,0,0.3}),"B/envelope-"+std::to_string(static_cast<int>(kind)));
}
// E1: bubble at 0.1 MPa through the passage to dew at 0.1 MPa: closure on the other branch and the candidate
// against C3; E4: the predicted-anchor step budget.
void mixture() {
    const HelmholtzModel m(properties_support::neutral_binary());
    const auto start=envelope_cell({0.5,0.5},120,1e5,{EnvelopeCoordinate::LogPressure,0,std::log(1e5)});
    auto dew=envelope_cell({0.5,0.5},170,1e5,{});dew.envelope.reset();dew.pressure=FreeValue::fixed(1e5);
    dew.phases[0].kind=PhaseKind::Vapor;dew.phases[1].kind=PhaseKind::Liquid;
    const auto cold=solve_equilibrium(m,dew);if(!cold.has_state())throw std::runtime_error("E1 cold dew: "+cold.message);
    const auto c3=solved(m,critical_cell({0.5,0.5},260,6000),"E1/C3");
    {
        const std::string id="E1/"+std::to_string(1.0);
        const auto t=trace(m,start,0.02,0.01);
        if(!t.success||t.critical_points.size()!=1)throw std::runtime_error(id+" trace: "+t.message);
        same_state(id+"/dew",t.states.back(),cold);
        same_state(id+"/candidate",t.critical_points[0],c3);
        double tmax=0,pmax=0;std::size_t switches=0,straddle=0;
        for(std::size_t k=0;k<t.states.size();++k) {
            tmax=std::max(tmax,*t.states[k].temperature);pmax=std::max(pmax,*t.states[k].pressure);
            switches+=k&&(t.conditions[k].kind!=t.conditions[k-1].kind||t.conditions[k].component!=t.conditions[k-1].component);
            if(k&&t.states[k].declaration->phases[0].kind!=t.states[k-1].declaration->phases[0].kind)straddle=t.conditions[k].component;
        }
        recording_support::observe((id+"/trace").c_str(),"envelope-trace","K",tmax,std::nullopt,
            recording_support::fields({recording_support::field("max_P_Pa",pmax),recording_support::field("states",t.states.size()),
                recording_support::field("attempts",t.attempts),recording_support::field("switches",switches),recording_support::field("straddle_component",straddle),
                recording_support::field("seed",t.critical_seeds[0])}),
            "Observation","accepted-state extrema and condition switches","declared envelope branch only","critical-envelope");
    }
    // E4 (#571, Opus P05): bubble at 180 K through the passage to the dew point at the initial pressure in the
    // predicted-anchor step budget, with one identity swap and the candidate and closing dew state of the dense trace.
    const auto p05=envelope_cell({0.5,0.5},180,1.6e6,{EnvelopeCoordinate::LogTemperature,0,std::log(180.)});
    for(const auto [cap,budget]:std::vector<std::pair<double,std::size_t>>{{0.05,120},{0.2,80}}) {
        const std::string id="E4/"+std::to_string(cap);
        ContinuationRequest r{p05,p05,std::nullopt,{0.1,1e-6,3000,8},{}};r.options.envelope_max_log_step=cap;r.options.envelope_straddle=0.01;
        r.start_problem.solver_options.values[epcsaft::policy::row]=1e-10;
        const auto t=continue_equilibrium(m,r);
        if(t.stop!=ContinuationStop::reached_target||t.critical_points.size()!=1)throw std::runtime_error(id+" trace: "+t.message);
        std::size_t swaps=0;
        for(std::size_t k=1;k<t.states.size();++k)swaps+=t.states[k].declaration->phases[0].kind!=t.states[k-1].declaration->phases[0].kind;
        recording_support::observe((id+"/budget").c_str(),"envelope-budget","1",static_cast<double>(t.accepted_steps),std::nullopt,
            recording_support::fields({recording_support::field("attempts",t.attempts),recording_support::field("retries",t.retries),
                recording_support::field("swaps",swaps)}),"Observation","accepted steps to closure","declared envelope branch only","critical-envelope");
        CHECK(swaps==1&&t.accepted_steps<=budget&&(cap>0.1||t.attempts<=150));
        compare(id+"/critical-T",*t.critical_points[0].temperature,265.9667440343549,1,1e-8);
        compare(id+"/critical-P",*t.critical_points[0].pressure,6842019.677409215,1,0.1);
        compare(id+"/dew-T",*t.states.back().temperature,232.31429972874673,1,1e-8);
    }
}
// E3: propane saturation; Clapeyron dP/dT = dh_res/(T dv) at 250, 300, 350 K, and the R1 pure endpoint:
// last state 2W < g <= 4W and the candidate equal to C1.
void pure() {
    const HelmholtzModel m(caloric_specification_checks::propane_parameters());
    const auto start=envelope_cell({1},250,1e5,{EnvelopeCoordinate::LogTemperature,0,std::log(250.)});
    for(const auto [t,seed]:std::vector<std::pair<double,double>>{{250.,2e5},{300.,1e6},{350.,2.5e6}}) {  // P seeds in the two-root region
        auto sat=start;sat.envelope.reset();sat.temperature=FreeValue::fixed(t);sat.pressure=FreeValue::variable(seed);
        const auto r=solve_equilibrium(m,sat);if(!r.has_state())throw std::runtime_error("E3 saturation: "+r.message);
        const auto slope=solved_state_actions(compile_problem(m,sat),r,{{{SolvedStateObservableKind::Pressure}},{{1,0,{0},{}}}}).results[0];CHECK(slope.action);
        const auto liquid=m.state(t,r.molar_densities[0],{1}),vapor=m.state(t,r.molar_densities[1],{1});
        const double clapeyron=(vapor.residual_enthalpy-liquid.residual_enthalpy)/(t*(1/r.molar_densities[1]-1/r.molar_densities[0]));
        compare("E3/clapeyron-"+std::to_string(static_cast<int>(t)),*slope.action,clapeyron,std::abs(clapeyron),1e-8);
    }
    const auto c1=solved(m,critical_cell({1},375,5000),"E3/C1");const double W=policy::DENSITY_ANCHOR_LOG_WINDOW.value;
    {
        const std::string id="E3/"+std::to_string(1.0);
        const auto t=trace(m,start,0.02,0.01);
        if(t.stop!=ContinuationStop::critical_endpoint||t.critical_points.size()!=1)throw std::runtime_error(id+" trace: "+t.message);
        const auto& last=t.states.back();const double g=std::abs(std::log(last.molar_densities[0]/last.molar_densities[1]));
        compare(id+"/last-gap",g>2*W&&g<=4*W?1:0,1,1,1e-10,recording_support::field("gap",g));
        same_state(id+"/candidate",t.critical_points[0],c1);
    }
}
// Gross & Sadowski (2002) Table 1/2 binaries with 2B self sites and Wolbach-Sandler cross edges, as FeOs combines them.
inline Parameters azeotropic_binary(std::array<double,2> m_,std::array<double,2> s,std::array<double,2> e,
                                    std::array<double,2> eab,std::array<double,2> kab,double kij) {
    auto p=properties_support::parameters(2);p.segment_count={m_[0],m_[1]};p.sigma_constant={s[0],s[1]};
    p.dispersion_energy_over_k={e[0],e[1]};p.molar_mass_kg_per_mol={0.05,0.05};p.k_ij[1]=p.k_ij[2]=kij;
    std::vector<int> owner;for(int i=0;i<2;++i)if(kab[i]>0)owner.insert(owner.end(),{i,i});
    const std::size_t n=owner.size();p.site_component=owner;p.site_multiplicity.assign(n,1);
    p.association_active.assign(n*n,0);p.association_energy_over_k.assign(n*n,0.);p.association_volume.assign(n*n,0.);
    for(std::size_t a=0;a<n;a+=2)for(std::size_t b=1;b<n;b+=2) {  // donor a, acceptor b
        const auto i=static_cast<std::size_t>(owner[a]),j=static_cast<std::size_t>(owner[b]);const double f=std::pow(std::sqrt(s[i]*s[j])/(0.5*(s[i]+s[j])),3);
        for(const auto [u,v]:{std::pair{a,b},std::pair{b,a}}) {
            p.association_active[u*n+v]=1;p.association_energy_over_k[u*n+v]=0.5*(eab[i]+eab[j]);
            p.association_volume[u*n+v]=std::sqrt(kab[i]*kab[j])*f;
        }
    }
    return p;
}
inline EquilibriumProblem azeotrope_cell(FreeValue t,FreeValue pressure) {
    auto p=equilibrium_support::flash(2);p.feed.values={0.5,0.5};p.solver_options.values[epcsaft::policy::row]=1e-12;
    p.phases.push_back(equilibrium_support::phase("V",PhaseKind::Vapor,AmountSpec::pinned(0)));
    p.temperature=t;p.pressure=pressure;p.azeotrope=true;
    return p;
}
// A1: azeotrope cells from the z = 0.5 seed against FeOs 0.10.1 binary_azeotrope at identical parameters (#571);
// methane/ethane, which has none, and its mixture critical point (x = y trivially) are not admitted.
void azeotrope() {
    const HelmholtzModel pb(azeotropic_binary({2.9997,2.4653},{3.2522,3.6478},{233.40,287.35},{2276.8,0},{0.015268,0},0.020)),
        ew(azeotropic_binary({2.3827,1.0656},{3.1771,3.0007},{198.24,366.51},{2653.4,2500.7},{0.032384,0.034868},0.0));
    struct Case { std::string id; const HelmholtzModel* m; EquilibriumProblem p; double x, value; bool fixed_t; };
    const std::vector<Case> cases{
        {"A1/propanol-benzene-313",&pb,azeotrope_cell(FreeValue::fixed(313.15),FreeValue::variable()),0.12076600330935693,26070.101425286382,true},
        {"A1/ethanol-water-343",&ew,azeotrope_cell(FreeValue::fixed(343.15),FreeValue::variable()),0.7912786135777201,73377.95467480272,true},
        {"A1/ethanol-water-101325",&ew,azeotrope_cell(FreeValue::variable(),FreeValue::fixed(101325)),0.802916090415992,351.2092202913194,false}};
    for(const auto& k:cases) {
        const auto r=solve_equilibrium(*k.m,k.p);if(!r.success)throw std::runtime_error(k.id+": "+r.message);
        compare(k.id+"/x",r.mole_fractions[0],k.x,1,1e-9);compare(k.id+"/y",r.mole_fractions[2],k.x,1,1e-9);
        if(k.fixed_t)compare(k.id+"/P",*r.pressure,k.value,k.value,1e-8);else compare(k.id+"/T",*r.temperature,k.value,1,1e-6);
    }
    const HelmholtzModel me(properties_support::neutral_binary());
    for(const auto& [id,p]:std::vector<std::pair<std::string,EquilibriumProblem>>{
            {"A1/methane-ethane-200",azeotrope_cell(FreeValue::fixed(200),FreeValue::variable())},
            {"A1/methane-ethane-critical",azeotrope_cell(FreeValue::fixed(265.96674403),FreeValue::variable(6.84e6))}}) {
        const auto r=solve_equilibrium(me,p);
        recording_support::observe((id+"/refusal").c_str(),"azeotrope-refusal","1",r.has_state()?1.0:0.0,std::nullopt,
            recording_support::field("message",r.message),"Observation","typed driver refusal","no certificate that none exists","critical-envelope");
        CHECK(!r.success&&!r.has_state());
    }
}
}  // namespace critical_envelope_checks
