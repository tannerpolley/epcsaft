#pragma once
// Issue 224 B: critical-equation candidates and phase envelopes on the Gross (2001) propane and
// methane/ethane packets; numerical verification only.
#include "caloric_specification_checks.hpp"
namespace critical_envelope_checks {
using namespace epcsaft;
using namespace epcsaft::equilibrium;
using caloric_specification_checks::compare;
inline EquilibriumProblem critical_cell(const std::vector<double>& z, double t, double rho) {
    auto p=equilibrium_support::flash(z.size());p.feed.values=z;p.row_tolerance=1e-11;
    p.temperature=FreeValue::variable(t);p.pressure=FreeValue::variable(1e6);p.critical=CriticalDeclaration{rho,{}};
    return p;
}
inline EquilibriumResult solved(const HelmholtzModel& m,const EquilibriumProblem& p,const std::string& id) {
    const auto r=solve_equilibrium(m,p);
    if(!r.has_state()||!r.critical) throw std::runtime_error(id+": "+r.message);
    const auto& k=*r.critical;
    recording_support::observe((id+"/candidate").c_str(),"critical-equation-candidate","K",r.temperature[0],std::nullopt,
        recording_support::fields({recording_support::field("P_Pa",r.pressure[0]),recording_support::field("rho_mol_m3",r.molar_densities[0]),
            recording_support::field("eigenvalues",k.eigenvalues),recording_support::field("zero_band",k.zero_band),recording_support::field("cubic",k.cubic),
            recording_support::field("condition",k.condition_estimate),recording_support::field("backward_error",k.backward_error)}),
        "Observation","root identity and spectrum of M","candidate only, no stability certificate","critical-envelope");
    return r;
}
inline EquilibriumProblem envelope_cell(const std::vector<double>& z,double t,double pressure,EnvelopeCondition condition) {
    auto p=equilibrium_support::flash(z.size());p.feed.values=z;p.row_tolerance=1e-11;
    p.phases.push_back(equilibrium_support::phase("V",PhaseKind::Vapor,AmountSpec::pinned(0)));
    p.temperature=FreeValue::variable(t);p.pressure=FreeValue::variable(pressure);p.envelope=condition;
    return p;
}
inline ContinuationResult trace(const HelmholtzModel& m,const EquilibriumProblem& p,double step,double straddle) {
    ContinuationRequest r{p,p,std::nullopt,{1,1e-7,4000,6},{}};r.options.envelope_max_log_step=step;r.options.envelope_straddle=straddle;
    return continue_equilibrium(m,r);
}
inline void same_state(const std::string& id,const EquilibriumResult& a,const EquilibriumResult& b) {
    compare(id+"/T",a.temperature[0],b.temperature[0],b.temperature[0],1e-8);
    compare(id+"/P",a.pressure[0],b.pressure[0],b.pressure[0],1e-8);
    for(std::size_t i=0;i<a.molar_densities.size();++i)compare(id+"/rho-"+std::to_string(i),a.molar_densities[i],b.molar_densities[i],b.molar_densities[i],1e-8);
    for(std::size_t i=0;i<a.mole_fractions.size();++i)compare(id+"/x-"+std::to_string(i),a.mole_fractions[i],b.mole_fractions[i],1,1e-8);
}
// C1/C2: the pure propane candidate is admitted directly; pressure-closed cp is unavailable while T, P and
// density survive; P_rho/(RT) and rho P_rhorho/(RT) from pressure_state_hessian vanish within a derived floor.
inline void pure_critical() {
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
inline void binary_critical() {
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
                    side[s]=kind==SolvedStateObservableKind::Temperature?rr.temperature[0]:kind==SolvedStateObservableKind::Pressure?rr.pressure[0]:rr.molar_densities[0];
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
        const auto state=m.state(r.temperature[0],r.molar_densities[0],{r.mole_fractions[0],r.mole_fractions[1]});
        compare(id+"/cp-residual",cp[0].value.value_or(NAN),*state.residual_isobaric_heat_capacity,std::abs(*state.residual_isobaric_heat_capacity),1e-10,
            recording_support::field("P_rho_over_RT",state.pressure_density_derivative/(GAS_CONSTANT_J_PER_MOL_K*r.temperature[0])));
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
    caloric_specification_checks::ladder(c,y,scale,id);
}
inline void curvature() {
    const HelmholtzModel m(properties_support::neutral_binary()),a(properties_support::associating_binary());
    off_solution(m,critical_cell({0.27,0.73},260,6000),"B/critical");
    off_solution(a,critical_cell({0.1385,0.8615},376.88,5000),"B/critical-associating");
    for(const auto kind:{EnvelopeCoordinate::LogTemperature,EnvelopeCoordinate::LogPressure,EnvelopeCoordinate::LogK})
        off_solution(m,envelope_cell({0.27,0.73},180,1e5,{kind,0,0.3}),"B/envelope-"+std::to_string(static_cast<int>(kind)));
}
// E1: bubble at 0.1 MPa through the passage to dew at 0.1 MPa at three bracketing levels: closure on the other
// branch, the candidate against C3, and Hermite seed refinement.
inline void mixture() {
    const HelmholtzModel m(properties_support::neutral_binary());
    const auto start=envelope_cell({0.5,0.5},120,1e5,{EnvelopeCoordinate::LogPressure,0,std::log(1e5)});
    auto dew=envelope_cell({0.5,0.5},170,1e5,{});dew.envelope.reset();dew.pressure=FreeValue::fixed(1e5);
    dew.phases[0].kind=PhaseKind::Vapor;dew.phases[1].kind=PhaseKind::Liquid;
    const auto cold=solve_equilibrium(m,dew);if(!cold.has_state())throw std::runtime_error("E1 cold dew: "+cold.message);
    const auto c3=solved(m,critical_cell({0.5,0.5},260,6000),"E1/C3");
    std::vector<std::vector<double>> seeds;
    for(const double f:{1.0,0.5,0.25}) {
        const std::string id="E1/"+std::to_string(f);
        const auto t=trace(m,start,0.02*f,0.01*f);
        if(!t.success||t.critical_points.size()!=1)throw std::runtime_error(id+" trace: "+t.message);
        seeds.push_back(t.critical_seeds[0]);
        same_state(id+"/dew",t.states.back(),cold);
        same_state(id+"/candidate",t.critical_points[0],c3);
        double tmax=0,pmax=0;std::size_t switches=0,straddle=0;
        for(std::size_t k=0;k<t.states.size();++k) {
            tmax=std::max(tmax,t.states[k].temperature[0]);pmax=std::max(pmax,t.states[k].pressure[0]);
            switches+=k&&(t.conditions[k].kind!=t.conditions[k-1].kind||t.conditions[k].component!=t.conditions[k-1].component);
            if(k&&t.states[k].declaration->phases[0].kind!=t.states[k-1].declaration->phases[0].kind)straddle=t.conditions[k].component;
        }
        recording_support::observe((id+"/trace").c_str(),"envelope-trace","K",tmax,std::nullopt,
            recording_support::fields({recording_support::field("max_P_Pa",pmax),recording_support::field("states",t.states.size()),
                recording_support::field("attempts",t.attempts),recording_support::field("switches",switches),recording_support::field("straddle_component",straddle),
                recording_support::field("seed",t.critical_seeds[0])}),
            "Observation","accepted-state extrema and condition switches","declared envelope branch only","critical-envelope");
    }
    const std::array<double,3> solved{std::log(c3.temperature[0]),std::log(c3.pressure[0]),std::log(c3.molar_densities[0])};
    for(std::size_t k=0;k<3;++k) {  // successive differences only; no convergence order is claimed
        const double d1=std::abs(seeds[1][k]-seeds[0][k]),d2=std::abs(seeds[2][k]-seeds[1][k]),ratio=d1/d2;
        recording_support::observe(("E1/hermite-"+std::to_string(k)).c_str(),"seed-refinement","1",seeds[2][k],solved[k],
            recording_support::fields({recording_support::field("h_to_h2",d1),recording_support::field("h2_to_h4",d2),recording_support::field("ratio",ratio)}),"Observation","bracketing-step refinement","seed need not meet the solved tolerance","critical-envelope");
    }
}
// E3: propane saturation; Clapeyron dP/dT = dh_res/(T dv) at 250, 300, 350 K, and the R1 pure endpoint at
// three step levels: last state 2W < g <= 4W and the candidate equal to C1.
inline void pure() {
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
    for(const double f:{1.0,0.5,0.25}) {
        const std::string id="E3/"+std::to_string(f);
        const auto t=trace(m,start,0.02*f,0.01*f);
        if(t.stop!=ContinuationStop::critical_endpoint||t.critical_points.size()!=1)throw std::runtime_error(id+" trace: "+t.message);
        const auto& last=t.states.back();const double g=std::abs(std::log(last.molar_densities[0]/last.molar_densities[1]));
        compare(id+"/last-gap",g>2*W&&g<=4*W?1:0,1,1,1e-10,recording_support::field("gap",g));
        same_state(id+"/candidate",t.critical_points[0],c1);
    }
}
}  // namespace critical_envelope_checks
