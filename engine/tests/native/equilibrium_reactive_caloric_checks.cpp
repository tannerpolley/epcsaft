#include "reactive_caloric_shared.hpp"
#include "binary128_support.hpp"
#include "equilibrium_support.hpp"

#include <chrono>
#include <cstring>
#include <map>
#include <memory>
#include <numeric>

// Issue #225 F1/F3/F4/F5/F6. The oracles are native value calculations, except the F6
// free-intensive pair ladders, which difference centered exact first actions.
namespace reactive_caloric_checks {
using namespace epcsaft;
using namespace epcsaft::equilibrium;
using Direction = SolvedStateActionDirection;
using Observable = SolvedStateObservable;
constexpr double R = GAS_CONSTANT_J_PER_MOL_K;
using binary128_support::NamedDirection;
using binary128_support::Quantity;
#ifndef _WIN32
using binary128_support::Wide;
using binary128_support::WideIonic;
#endif

struct Case {
    std::string name;
    Parameters parameters;
    std::vector<ActiveParameter> active;
    ThermochemistryRecord records;
    EquilibriumProblem problem;
    bool ionic = false;
};

inline Case input(bool ionic, bool fixed, bool unreactive = false) {
    Case c;
    c.ionic = ionic;
    c.name = unreactive ? "U" : std::string(ionic ? "I" : "E") + (fixed ? "-fix" : "-sys");
    c.parameters = ionic ? properties_support::aqueous_electrolyte() : properties_support::parameters(2);
    if (!ionic) {
        c.parameters.segment_count = {1.6069, 1.6069};
        c.parameters.sigma_constant = {3.5206, 3.5206};
        c.parameters.dispersion_energy_over_k = {191.42, 191.42};
        c.parameters.molar_mass_kg_per_mol = {0.03007, 0.03007};
    }
    c.active = {{ParameterFamily::DispersionEnergyOverK, {0}}};
    if (!ionic) c.active.push_back({ParameterFamily::DispersionEnergyOverK, {1}});
    c.records.standard_pressure = 1e5;
    c.records.component_intervals.resize(ionic ? 3 : 2);
    c.records.component_intervals[0] = {{200.0, ionic ? 500.0 : 400.0, true, true,
        {ionic ? IdealPolynomial{{33.59}, 298.15, -241826.0, 188.83}
               : IdealPolynomial{{52.0, 0.12}, 300.0, 0.0, 0.0}, 1e5}}};
    if (unreactive) c.records.component_intervals[1] = c.records.component_intervals[0];
    c.problem = equilibrium_support::flash(ionic ? 3 : 2);
    c.problem.temperature = FreeValue::fixed(ionic ? 298.15 : 260.0);
    c.problem.pressure = FreeValue::fixed(ionic ? 1e5 : 3e6);
    c.problem.feed.values = ionic ? std::vector<double>{1.0, 1e-6, 1e-6} : std::vector<double>{2.0, 0.2};
    c.problem.source_reference = NeutralReferenceDeclaration{0,
        ionic ? std::vector<double>{0.0, 1.0, 1.0} : std::vector<double>{0.0, 1.0}, 1.0};
    if (!unreactive) c.problem.reactions = {{ionic ? "H2O=Na++Cl-" : "2A=B",
        ionic ? std::vector<double>{-1.0, 1.0, 1.0} : std::vector<double>{-2.0, 1.0},
        std::nullopt, std::nullopt,
        ionic ? ReactionLogPolynomial{std::log(1e-8), 0, 0, 0, 298.15, 250, 350, "analytic-source"}
              : ReactionLogPolynomial{std::log(0.3), 200, 3, -0.004, 300, 250, 350, "analytic-source"},
        ReactionReference{fixed ? std::optional<double>(ionic ? 2e5 : 5e6) : std::nullopt,
            ionic ? ReferenceSourceBasis::CommonMolalityInfiniteDilution
                  : ReferenceSourceBasis::RawMoleFractionInfiniteDilution}}};
    return c;
}

inline Direction zero(const Case& c) {
    Direction d;
    d.feed_amounts_mol.assign(c.problem.feed.values.size(), 0.0);
    d.active_parameters.assign(c.active.size(), 0.0);
    d.reaction_coefficients.assign(c.problem.reactions.size(), {});
    return d;
}
inline Direction add(Direction a, const Direction& b, double scale = 1.0) {
    a.temperature_k += scale*b.temperature_k;
    a.pressure_pa += scale*b.pressure_pa;
    for (std::size_t i=0;i<a.feed_amounts_mol.size();++i) a.feed_amounts_mol[i] += scale*b.feed_amounts_mol[i];
    for (std::size_t i=0;i<a.active_parameters.size();++i) a.active_parameters[i] += scale*b.active_parameters[i];
    for (std::size_t r=0;r<a.reaction_coefficients.size();++r)
        for (std::size_t k=0;k<4;++k) a.reaction_coefficients[r][k] += scale*b.reaction_coefficients[r][k];
    return a;
}
inline std::vector<NamedDirection> directions(const Case& c) {
    std::vector<NamedDirection> result;
    auto put = [&](std::string name, double step, Direction d, bool parameter=false) {
        result.push_back({name,step,d,{d},parameter});
    };
    auto d=zero(c); d.temperature_k=1; put("T",0.2,d);
    d=zero(c); d.pressure_pa=1e5; put("P",0.2,d);
    d=zero(c); d.feed_amounts_mol[0]=1; put("water-feed",c.ionic?0.001:0.002,d);
    if (c.ionic) {d=zero(c);d.feed_amounts_mol[1]=d.feed_amounts_mol[2]=1e-6;put("salt-feed",0.01,d);}
    d=zero(c); for (auto& value:d.active_parameters) value=1; put("epsilon",0.02,d,true);
    result.back().chart.clear();
    for (std::size_t a=0;a<c.active.size();++a) {auto part=zero(c);part.active_parameters[a]=1;result.back().chart.push_back(part);}
    if (!c.problem.reactions.empty()) {
        d=zero(c);d.reaction_coefficients[0][c.ionic?0:1]=1;put(c.ionic?"reaction-a":"reaction-b",c.ionic?0.001:2.0,d,true);
    }
    return result;
}

inline std::vector<Quantity> quantities(std::size_t n) {
    std::vector<Quantity> q;
    auto put=[&](const char* name,const char* unit,double atol,SolvedStateObservableKind kind,std::size_t component=0) {
        Observable o;o.kind=kind;o.component=component;q.push_back({name,unit,atol,o});
    };
    put("T","K",1e-8,SolvedStateObservableKind::Temperature);
    put("P","Pa",1e-3,SolvedStateObservableKind::Pressure);
    put("n","mol",1e-10,SolvedStateObservableKind::PhaseAmount);
    put("V_EOS","m^3",1e-12,SolvedStateObservableKind::PhaseVolume);
    put("rho","mol/m^3",1e-6,SolvedStateObservableKind::PhaseMolarDensity);
    put("eta","1",1e-11,SolvedStateObservableKind::PhasePackingFraction);
    put("P_EOS","Pa",1e-3,SolvedStateObservableKind::PhasePressure);
    for(std::size_t i=0;i<n;++i) {
        for (auto item: {std::tuple{"x", "1",1e-11,SolvedStateObservableKind::PhaseMoleFraction},
                        std::tuple{"n_i","mol",1e-10,SolvedStateObservableKind::PhaseComponentAmount},
                        std::tuple{"mu_res_RT","1",1e-11,SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt},
                        std::tuple{"ln_a","1",1e-11,SolvedStateObservableKind::PhaseLogActivity}}) {
            put(std::get<0>(item),std::get<1>(item),std::get<2>(item),std::get<3>(item),i);
            q.back().name += "-"+std::to_string(i);
        }
    }
    const char* names[]={"h","u_G","s","cp","cv_G"};
    for(std::size_t family:{0U,10U}) for(std::size_t k=0;k<5;++k) {
        Observable o;o.kind=SolvedStateObservableKind::PhaseProperty;o.property=static_cast<PropertyObservable>(family+k);
        q.push_back({std::string(family==0?"ideal-":"total-")+names[k],k<2?"J/mol":"J/(mol K)",k<2?1e-5:1e-7,o});
    }
    put("v_G","m^3/mol",1e-12,SolvedStateObservableKind::PhaseGibbsMolarVolume);
    // Last four quantities are native extensive products; no second runtime owner.
    for(auto item:{std::tuple{"H","J",1e-5},std::tuple{"U_G","J",1e-5},
                   std::tuple{"S","J/K",1e-7},std::tuple{"V_G","m^3",1e-12}})
        q.push_back({std::get<0>(item),std::get<1>(item),std::get<2>(item),{}});
    return q;
}
inline std::size_t index(const std::vector<Quantity>& q,const std::string& name) {
    auto it=std::find_if(q.begin(),q.end(),[&](const auto& item){return item.name==name;});
    CHECK(it!=q.end());return it-q.begin();
}
inline double column(const ReferenceCalorics& ref, std::size_t a,std::size_t b,std::size_t i) {
    for(const auto& d:ref.gibbs) if(d.temperature_order==a&&d.pressure_order==b&&!d.parameter) {
        if(d.column.status!=SolvedStateActionStatus::Available) throw std::runtime_error(d.column.message);
        return d.column.values.at(i);
    }
    throw std::runtime_error("missing raw anchored Gibbs column "+std::to_string(a)+","+std::to_string(b));
}
struct Local {
    std::array<double,5> ideal{}, total{};
    double g=0,v=0,vt=0,vp=0,quotient=0,rho=0;
    PropertyActionDiagnostics diagnostics;
};
inline Local local(const HelmholtzModel& model,const CompiledSystem& compiled,double T,double P,
                   const std::vector<double>& x,double anchor,DensityBranch branch=DensityBranch::Liquid) {
    Local result;
    const auto state=model.state_at_pressure(T,P,x,branch,anchor);
    result.rho=state.molar_density;
    const auto raw=model.pressure_actions(T,P,x,branch,
        {PropertyObservable::MolarVolumeTemperatureDerivative,PropertyObservable::MolarVolumePressureDerivative},
        {{0,0,std::vector<double>(x.size(),0),{}}},anchor);
    for(const auto& item:raw.results) CHECK(item.status==PropertyActionStatus::Available&&item.value);
    result.diagnostics=raw.diagnostics;
    result.v=1/result.rho;result.vt=*raw.results[0].value;result.vp=*raw.results[1].value;
    if (!compiled.problem.reactions.empty()) {
        const auto ref=reference_calorics(compiled,T,P,true,true,false,2);
        if(ref.status!=SolvedStateActionStatus::Available) throw std::runtime_error(ref.message);
        for(std::size_t i=0;i<x.size();++i) {
            const double g=column(ref,0,0,i),gt=column(ref,1,0,i);
            result.ideal[0]+=x[i]*(g-T*gt);
            result.ideal[2]-=x[i]*(gt+R*std::log(x[i]*result.rho));
            result.ideal[3]-=x[i]*T*column(ref,2,0,i);
            result.v+=x[i]*column(ref,0,1,i);
            result.vt+=x[i]*column(ref,1,1,i);
            result.vp+=x[i]*column(ref,0,2,i);
            result.g+=x[i]*g;
        }
        result.ideal[0]+=R*T;result.ideal[2]+=R;result.ideal[3]+=R;
        result.ideal[1]=result.ideal[0]-R*T-P*(result.v-1/result.rho);
        const double eos_quotient=T*(*raw.results[0].value)*(*raw.results[0].value)/(*raw.results[1].value);
        result.ideal[4]=result.ideal[3]-R+T*result.vt*result.vt/result.vp-eos_quotient;
    } else {
        CHECK(state.ideal_enthalpy&&state.ideal_internal_energy&&state.ideal_entropy
              &&state.ideal_isobaric_heat_capacity&&state.ideal_isochoric_heat_capacity);
        result.ideal={*state.ideal_enthalpy,*state.ideal_internal_energy,*state.ideal_entropy,
                      *state.ideal_isobaric_heat_capacity,*state.ideal_isochoric_heat_capacity};
        for(std::size_t i=0;i<x.size();++i) {
            const auto& correlation=model.thermochemistry()->component_intervals[i][0].correlation;
            const auto thermal=evaluate_ideal_correlation(correlation,T,1e5);
            // g_rho = g_P0 + RT ln(RT rho0/P0); rho0 = 1 mol/m^3.
            result.g+=x[i]*(thermal.h-T*thermal.s+R*T*std::log(R*T/1e5));
        }
    }
    CHECK(state.residual_isobaric_heat_capacity);
    const std::array<double,5> residual={state.residual_enthalpy,state.residual_internal_energy,state.residual_entropy,
        *state.residual_isobaric_heat_capacity,state.residual_isochoric_heat_capacity};
    for(std::size_t k=0;k<5;++k) result.total[k]=result.ideal[k]+residual[k];
    const auto energy=model.evaluate(T,1/result.rho,x,DerivativeOrder::Value);
    double mixing=-1.0;for(std::size_t i=0;i<x.size();++i) mixing+=x[i]*std::log(x[i]*result.rho);
    result.g+=R*T*energy.value+P/result.rho+R*T*mixing;
    result.quotient=T*result.vt*result.vt/result.vp;
    CHECK(std::isfinite(result.v)&&std::isfinite(result.vp)&&result.vp!=0);
    return result;
}

struct Sample {std::vector<double> values;Local phase;EquilibriumResult solved;std::vector<Local> phases;};
inline Sample values(const HelmholtzModel& model,const CompiledSystem& compiled,const EquilibriumResult& solved,
                     const std::vector<Quantity>& q) {
    CHECK(solved.has_state());
    const std::size_t components=model.parameters().component_count();
    const double T=*solved.temperature,P=*solved.pressure;
    Sample sample;sample.solved=solved;
    std::vector<Local> phases;
    std::vector<State> states;
    for(std::size_t phase=0;phase<compiled.problem.phases.size();++phase) {
        const std::vector<double> x(solved.mole_fractions.begin()+phase*components,
                                   solved.mole_fractions.begin()+(phase+1)*components);
        const double rho=solved.molar_densities[phase];
        const auto branch=compiled.problem.phases[phase].kind==PhaseKind::Vapor?DensityBranch::Vapor:DensityBranch::Liquid;
        phases.push_back(local(model,compiled,T,P,x,rho,branch));
        states.push_back(model.state(T,rho,x));
        CHECK(states.back().pressure_density_derivative>0);
        CHECK(std::abs(states.back().pressure-P)<=1e-3+1e-9*P);
        CHECK(std::abs(phases.back().rho-rho)<=1e-6+1e-10*rho);
    }
    sample.phase=phases.front();
    CHECK(solved.residual_norm()<=compiled.problem.solver_options.value(epcsaft::policy::row));
    for(const auto& item:q) {
        const auto& o=item.observable;const auto phase=o.phase;
        const double n=solved.amounts[phase],rho=solved.molar_densities[phase];
        const auto& state=states[phase];const auto& caloric=phases[phase];
        const auto x=solved.mole_fractions.begin()+phase*components;double value=0;
        switch(o.kind) {
        case SolvedStateObservableKind::Temperature:value=T;break;
        case SolvedStateObservableKind::Pressure:value=P;break;
        case SolvedStateObservableKind::PhaseAmount:value=n;break;
        case SolvedStateObservableKind::PhaseVolume:value=n/rho;break;
        case SolvedStateObservableKind::PhaseMolarDensity:value=rho;break;
        case SolvedStateObservableKind::PhasePackingFraction:value=solved.packing_fractions[phase];break;
        case SolvedStateObservableKind::PhasePressure:value=state.pressure;break;
        case SolvedStateObservableKind::PhaseMoleFraction:value=x[o.component];break;
        case SolvedStateObservableKind::PhaseComponentAmount:value=n*x[o.component];break;
        case SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt:value=require_component_property(state.residual_chemical_potential_over_rt[o.component]);break;
        case SolvedStateObservableKind::PhaseLogActivity:value=std::log(x[o.component]*rho)+require_component_property(state.residual_chemical_potential_over_rt[o.component]);break;
        case SolvedStateObservableKind::PhaseGibbsMolarVolume:value=caloric.v;break;
        case SolvedStateObservableKind::PhaseProperty:{auto k=std::size_t(o.property);value=k<5?caloric.ideal[k]:caloric.total[k-10];break;}
        }
        if(item.name=="H") value=n*sample.phase.total[0];
        if(item.name=="U_G") value=n*sample.phase.total[1];
        if(item.name=="S") value=n*sample.phase.total[2];
        if(item.name=="V_G") value=n*sample.phase.v;
        sample.values.push_back(value);
    }
    sample.phases=std::move(phases);
    return sample;
}

inline std::string context(const Case& c,const EquilibriumResult& solved,const Local& l) {
    std::string law;
    if(!c.problem.reactions.empty()) {
        const auto& polynomial=*c.problem.reactions[0].correlation;
        law=recording_support::fields({recording_support::field("source_a",polynomial.a),
            recording_support::field("source_b_K",polynomial.b),recording_support::field("source_c",polynomial.c),
            recording_support::field("source_d_per_K",polynomial.d),recording_support::field("source_Tref_K",polynomial.reference_temperature),
            recording_support::field("source_Tmin_K",polynomial.temperature_min),recording_support::field("source_Tmax_K",polynomial.temperature_max),
            recording_support::field("source_P_Pa",c.problem.reactions[0].reference->reference_pressure_pa.value_or(*solved.pressure))});
    }
    return recording_support::fields({recording_support::field("case",c.name),law,
        recording_support::field("record_cp_J_per_mol_K",c.ionic?std::vector<double>{33.59}:std::vector<double>{52.0,0.12}),
        recording_support::field("record_Tref_K",c.ionic?298.15:300.0),
        recording_support::field("record_h_J_per_mol",c.ionic?-241826.0:0.0),
        recording_support::field("record_s_J_per_mol_K",c.ionic?188.83:0.0),
        recording_support::field("T_K",*solved.temperature),recording_support::field("P_Pa",*solved.pressure),
        recording_support::field("feed_mol",c.problem.feed.values),recording_support::field("x",solved.mole_fractions),
        recording_support::field("n_mol",solved.amounts[0]),recording_support::field("rho",l.rho),
        recording_support::field("solver_status",solved.solver_status),recording_support::field("solver_residual",solved.residual_norm()),
        recording_support::field("iterations",solved.iterations),recording_support::field("evaluations",solved.evaluations),
        recording_support::field("v_G",l.v),recording_support::field("v_G_T_raw",l.vt),recording_support::field("v_G_P_raw",l.vp),
        recording_support::field("v_G_positive",l.v>0),recording_support::field("v_G_P_negative",l.vp<0),
        recording_support::field("quotient_J_per_mol_K",l.quotient),
        recording_support::field("quotient_sensitivity_vT_Pa",std::abs(2 * *solved.temperature*l.vt/l.vp)),
        recording_support::field("quotient_sensitivity_vP_Pa2_per_K",std::abs(l.quotient/l.vp)),
        recording_support::field("pressure_density_derivative",l.diagnostics.pressure_density_derivative),
        recording_support::field("pressure_residual",l.diagnostics.pressure_residual),
        recording_support::field("association_residual",l.diagnostics.association_residual),
        recording_support::field("site_rank",l.diagnostics.site_rank),
        recording_support::field("site_pivot_ratio",l.diagnostics.site_pivot_ratio)});
}
inline bool compare(const std::string& id,const Quantity& q,double exact,double reference,
                    const std::string& inputs,double rtol,const char* method="native-identity") {
    const double target=q.atol+rtol*std::abs(reference);
    const bool good=std::isfinite(exact)&&std::isfinite(reference)&&std::abs(exact-reference)<=target;
    recording_support::compare_custom(id.c_str(),q.name.c_str(),q.unit.c_str(),exact,reference,
        recording_support::fields({inputs,recording_support::field("target",target),recording_support::field("error",std::abs(exact-reference)),
            recording_support::field("error_over_target",std::abs(exact-reference)/target)}),
        "Relationship",method,"Issue 225 manufactured branch numerical verification",
        "abs(exact-reference)<=atol+rtol*abs(reference)",good,"reactive-caloric",id.c_str());
    return good;
}

class Actions {
    const CompiledSystem& compiled;
    const EquilibriumResult& solved;
    const std::vector<Quantity>& q;
    std::map<std::string,std::vector<double>> cache;
public:
    SolvedStateActionDiagnostics diagnostics;
    Actions(const CompiledSystem& c,const EquilibriumResult& s,const std::vector<Quantity>& quantities):compiled(c),solved(s),q(quantities){}
    const std::vector<double>& get(const std::vector<Direction>& ds) {
        std::string key;
        for(const auto& d:ds) {
            key+=recording_support::number(d.temperature_k)+","+recording_support::number(d.pressure_pa)+"/";
            for(double v:d.feed_amounts_mol) key+=recording_support::number(v)+",";
            key+="/";for(double v:d.active_parameters) key+=recording_support::number(v)+",";
            key+="/";for(const auto& coefficients:d.reaction_coefficients) for(double v:coefficients) key+=recording_support::number(v)+",";
            key+=";";
        }
        auto found=cache.find(key);if(found!=cache.end()) return found->second;
        std::vector<Observable> os;for(std::size_t i=0;i+4<q.size();++i) os.push_back(q[i].observable);
        Direction z;z.feed_amounts_mol.assign(compiled.n_feed.size(),0);z.active_parameters.assign(compiled.model->active_parameters().size(),0);
        const auto batch=solved_state_actions(compiled,solved,{os,ds.empty()?std::vector<Direction>{z}:ds});
        diagnostics=batch.diagnostics;
        std::vector<double> result;
        for(std::size_t i=0;i<batch.results.size();++i) {
            const auto& entry=batch.results[i];const auto& v=ds.empty()?entry.value:entry.action;
            if(entry.status!=SolvedStateActionStatus::Available||!v) throw std::runtime_error(
                "Issue225 "+q[i].name+" order "+std::to_string(ds.size())+" "+solved_state_status_name(entry.status)+" "+entry.diagnostic_message.value_or(""));
            result.push_back(*v);
        }
        const std::size_t n=index(q,"n");
        for(const auto& name:{"total-h","total-u_G","total-s","v_G"}) {
            const auto k=index(q,name);double product=0;
            if(ds.empty()) {result.push_back(result[n]*result[k]);continue;}
            for(std::size_t mask=0;mask<(1U<<ds.size());++mask) {
                std::vector<Direction> left,right;
                for(std::size_t j=0;j<ds.size();++j) ((mask>>j)&1?left:right).push_back(ds[j]);
                if(left.size()==ds.size()) product+=result[n]*get({})[k];
                else if(right.size()==ds.size()) product+=get({})[n]*result[k];
                else product+=get(left)[n]*get(right)[k];
            }
            result.push_back(product);
        }
        return cache.emplace(key,std::move(result)).first->second;
    }
};

// A binary128 oracle rebuilds its own values; it reads only the solved state and density (rebuilt=false).
inline Sample shifted(const Case& c,const Direction& delta,const std::vector<Quantity>& q,const EquilibriumResult& anchor,const Direction* first_action=nullptr,bool rebuilt=true) {
    auto p=c.parameters;auto problem=c.problem;
    problem.temperature.value+=delta.temperature_k;problem.pressure.value+=delta.pressure_pa;
    for(std::size_t i=0;i<problem.feed.values.size();++i) {problem.feed.values[i]+=delta.feed_amounts_mol[i];CHECK(problem.feed.values[i]>=0);}
    double charge=0;for(std::size_t i=0;i<problem.feed.values.size();++i) charge+=problem.feed.values[i]*p.charge_number[i];CHECK(std::abs(charge)<=1e-20);
    for(std::size_t a=0;a<c.active.size();++a) p.dispersion_energy_over_k[c.active[a].indices[0]]+=delta.active_parameters[a];
    for(std::size_t r=0;r<problem.reactions.size();++r) {
        auto& law=*problem.reactions[r].correlation;const auto& d=delta.reaction_coefficients[r];
        law.a+=d[0];law.b+=d[1];law.c+=d[2];law.d+=d[3];
        if(!problem.temperature.free) CHECK(problem.temperature.value>=law.temperature_min&&problem.temperature.value<=law.temperature_max);
    }
    HelmholtzModel model(p,c.active,c.records);
    const auto compiled=compile_problem(model,problem);
    const auto solved=solve_compiled(compiled,anchor);
    if(!solved.has_state()) throw std::runtime_error(c.name+" oracle solve: "+solved.message+" residual="+recording_support::number(solved.residual_norm()));
    if(!rebuilt) {Sample sample;sample.solved=solved;sample.phase.rho=solved.molar_densities[0];return sample;}
    auto sample=values(model,compiled,solved,q);
    if(first_action) {Actions actions(compiled,solved,q);sample.values=actions.get({*first_action});}
    return sample;
}
#ifndef _WIN32
#include "reactive_caloric_wide.hpp"
struct Difference {
    std::vector<Wide> value;
    std::vector<double> rounding;
    double precise_pressure_residual=0,precise_pressure_floor=0;
    double residual=0,pressure_residual=0,site_residual=0,minimum_v=std::numeric_limits<double>::infinity();
    double maximum_vp=-std::numeric_limits<double>::infinity(),denominator=1;
};
inline Difference difference(const Case& c,const std::vector<NamedDirection>& ds,double factor,
                             const std::vector<Quantity>& q,const EquilibriumResult& anchor,WideNeutral* wide,WideIonic* ionic) {
    const bool first_difference=(c.problem.temperature.free||c.problem.pressure.free)&&ds.size()==2;
    const auto perturbed=first_difference?std::vector<NamedDirection>{ds.back()}:ds;
    Difference out;out.value.assign(q.size(),0);out.rounding.assign(q.size(),0);
    std::vector<Wide> sum(q.size(),0);Wide denominator=1;
    for(const auto& d:perturbed) denominator*=Wide(2)*Wide(d.step)*Wide(factor);out.denominator=double(denominator);
    for(std::size_t mask=0;mask<(1U<<perturbed.size());++mask) {
        auto delta=zero(c);int sign=1;
        for(std::size_t k=0;k<perturbed.size();++k) {const double s=(mask>>k)&1?1:-1;sign*=int(s);delta=add(delta,perturbed[k].value,s*perturbed[k].step*factor);}
        const auto sample=shifted(c,delta,q,anchor,first_difference?&ds.front().value:nullptr,!wide&&!ionic);out.residual=std::max(out.residual,sample.solved.residual_norm());
        if(!wide&&!ionic) {
            out.pressure_residual=std::max(out.pressure_residual,std::abs(sample.phase.diagnostics.pressure_residual));
            out.site_residual=std::max(out.site_residual,std::abs(sample.phase.diagnostics.association_residual));
            out.minimum_v=std::min(out.minimum_v,sample.phase.v);out.maximum_vp=std::max(out.maximum_vp,sample.phase.vp);
        }
        const auto precise=wide?wide->values(q,ds,factor,mask,sample):ionic?ionic->values(q,ds,factor,mask,sample.solved,sample.phase.rho)
            :std::vector<Wide>(sample.values.begin(),sample.values.end());
        out.precise_pressure_residual=std::max(out.precise_pressure_residual,ionic?ionic->pressure_residual():wide?wide->pressure_residual():0.0);
        out.precise_pressure_floor=std::max(out.precise_pressure_floor,ionic?ionic->pressure_limit():wide?wide->pressure_limit():0.0);
        for(std::size_t k=0;k<q.size();++k) {sum[k]+=sign*precise[k];out.rounding[k]+=double(CppAD::abs(precise[k]));}
    }
    const double epsilon=wide||ionic?double(CppAD::numeric_limits<Wide>::epsilon()):std::numeric_limits<double>::epsilon();
    for(std::size_t k=0;k<q.size();++k) {out.value[k]=sum[k]/denominator;out.rounding[k]*=epsilon/double(denominator);}
    return out;
}
inline std::vector<double> scale(Actions& actions,const std::vector<NamedDirection>& ds,std::size_t size) {
    std::vector<double> sigma(size,0);std::vector<Direction> parts;
    std::function<void(std::size_t)> visit=[&](std::size_t k) {
        if(k==ds.size()) {const auto& v=actions.get(parts);for(std::size_t i=0;i<size;++i) sigma[i]+=std::abs(v[i]);return;}
        for(const auto& part:ds[k].chart) {parts.push_back(part);visit(k+1);parts.pop_back();}
    };visit(0);return sigma;
}
inline void ladder(const Case& c,Actions& actions,const std::vector<NamedDirection>& ds,
                   const std::vector<Quantity>& q,const Sample& center,const std::string& inputs,
                   WideNeutral* wide=nullptr,WideIonic* ionic=nullptr) {
    std::vector<Direction> request;std::string label;
    std::string literal_directions;
    for(const auto& d:ds) {
        request.push_back(d.value);if(!label.empty()) label+="+";label+=d.name;
        literal_directions=recording_support::fields({literal_directions,
            recording_support::field((d.name+"_step").c_str(),d.step),recording_support::field((d.name+"_delta_T_K").c_str(),d.value.temperature_k),
            recording_support::field((d.name+"_delta_P_Pa").c_str(),d.value.pressure_pa),
            recording_support::field((d.name+"_delta_feed_mol").c_str(),d.value.feed_amounts_mol),
            recording_support::field((d.name+"_delta_epsilon_K").c_str(),d.value.active_parameters)});
        for(std::size_t r=0;r<d.value.reaction_coefficients.size();++r) {
            const auto& v=d.value.reaction_coefficients[r];
            literal_directions=recording_support::fields({literal_directions,
                recording_support::field((d.name+"_delta_abcd").c_str(),std::vector<double>(v.begin(),v.end()))});
        }
    }
    const auto exact=actions.get(request);const auto diagnostics=actions.diagnostics;
    const auto sigma=scale(actions,ds,q.size());
    const bool free_chart=c.problem.temperature.free||c.problem.pressure.free;
    CHECK(free_chart==(!wide&&!ionic));
    std::vector<std::vector<NamedDirection>> orders{ds};
    if(free_chart&&ds.size()==2&&ds[0].name!=ds[1].name) orders.push_back({ds[1],ds[0]});
    struct Reference {double rich=0,refinement=0,target=0;std::size_t rung=0;bool resolved=false;std::string row;};
    std::vector<std::vector<Reference>> references(orders.size(),std::vector<Reference>(q.size()));
    const auto id_for=[&](std::size_t i,std::size_t rung) {return std::string(free_chart?"F6/":"F1/")+c.name+"/A"
        +std::to_string(ds.size())+"/"+label+"/"+q[i].name+"/rung-"+std::to_string(rung);};
    const char* method=free_chart&&ds.size()==2?"native-centered-first-actions":"native-centered-rebuilt-values";
    for(std::size_t order=0;order<orders.size();++order) {
        auto coarse=difference(c,orders[order],1.0,q,center.solved,wide,ionic);
        for(std::size_t rung=0;rung<6;++rung) {
            const double factor=std::ldexp(1.0,-int(rung+1));
            const auto fine=difference(c,orders[order],factor,q,center.solved,wide,ionic);
            for(std::size_t i=0;i<q.size();++i) if(!references[order][i].resolved) {
                const double rich=double((4*fine.value[i]-coarse.value[i])/3),refinement=double(CppAD::abs(fine.value[i]-coarse.value[i])/3);
                const double target=q[i].atol+1e-6*std::max(std::abs(rich),sigma[i]);
                const bool resolved=refinement<=target/4;
                const auto row=recording_support::fields({inputs,literal_directions,recording_support::field("directions",label),
                    recording_support::field("fine_step_factor",factor),recording_support::field("oracle_precision_bits",free_chart?53:113),recording_support::field("D_h",double(coarse.value[i])),
                    recording_support::field("D_h2",double(fine.value[i])),recording_support::field("Richardson",rich),
                    recording_support::field("refinement",refinement),recording_support::field("sigma",sigma[i]),
                    recording_support::field("atol",q[i].atol),recording_support::field("target",target),
                    recording_support::field("resolved",resolved),
                    recording_support::field("rounding_over_step",fine.rounding[i]),recording_support::field("oracle_residual",fine.residual),
                    free_chart?recording_support::fields({recording_support::field("oracle_pressure_residual",fine.pressure_residual),
                        recording_support::field("oracle_site_residual",fine.site_residual),
                        recording_support::field("oracle_minimum_v_G",fine.minimum_v),recording_support::field("oracle_maximum_v_G_P",fine.maximum_vp)})
                    :recording_support::fields({recording_support::field("wide_pressure_residual_Pa",fine.precise_pressure_residual),
                        recording_support::field("wide_pressure_closure_limit_Pa",fine.precise_pressure_floor)}),
                    recording_support::field("oracle_residual_over_step_product",fine.residual/fine.denominator),
                    recording_support::field("condition_times_residual_over_step_product",diagnostics.condition_estimate.value_or(0)*fine.residual/fine.denominator),
                    recording_support::field("action_condition",diagnostics.condition_estimate.value_or(0)),
                    recording_support::field("linear_residual",diagnostics.maximum_linear_residual.value_or(0)),
                    recording_support::field("central_residual",diagnostics.central_residual.value_or(0))});
                references[order][i]={rich,refinement,target,rung,resolved,row};
                const std::string id=id_for(i,rung)+(orders.size()==2?(order?"/reversed":"/original"):"");
                if(orders.size()==2) recording_support::observe(id.c_str(),q[i].name.c_str(),q[i].unit.c_str(),rich,std::nullopt,row,
                    "Refinement",method,"Reference uncertainty before production agreement","reactive-caloric",id.c_str());
                else if(!resolved&&rung<5) recording_support::observe(id.c_str(),q[i].name.c_str(),q[i].unit.c_str(),exact[i],rich,row,
                    "Refinement",method,"Unresolved coarse step; retained and refined","reactive-caloric",id.c_str());
            }
            if(std::all_of(references[order].begin(),references[order].end(),[](const auto& r){return r.resolved;})) break;
            coarse=fine;
        }
    }
    for(std::size_t i=0;i<q.size();++i) {
        // Select only by each order's first-resolved refinement, before comparing with production.
        std::size_t selected=0;
        if(orders.size()==2&&references[1][i].resolved&&(!references[0][i].resolved||references[1][i].refinement<references[0][i].refinement)) selected=1;
        const auto& r=references[selected][i];
        std::vector<double> estimates,richardson;std::vector<int> rungs,resolved;
        for(const auto& order:references) {estimates.push_back(order[i].refinement);richardson.push_back(order[i].rich);
            rungs.push_back(int(order[i].rung));resolved.push_back(order[i].resolved);}
        const bool good=r.resolved&&std::isfinite(exact[i])&&std::isfinite(r.rich)&&std::abs(exact[i]-r.rich)<=r.target;
        const auto row=recording_support::fields({r.row,recording_support::field("reference_orders",orders.size()==2?"original,reversed":"original"),
            recording_support::field("reference_refinements",estimates),recording_support::field("reference_Richardson",richardson),
            recording_support::field("reference_rungs",rungs),recording_support::field("reference_resolved",resolved),
            recording_support::field("selected_reference_order",r.resolved?(selected?"reversed":"original"):"none"),
            recording_support::field("reference_selection",r.resolved?"minimum first-resolved refinement with original-order tie break":"no resolved reference"),
            recording_support::field("error",std::abs(exact[i]-r.rich)),recording_support::field("error_over_target",std::abs(exact[i]-r.rich)/r.target)});
        const auto id=id_for(i,r.rung);
        recording_support::compare_custom(id.c_str(),q[i].name.c_str(),q[i].unit.c_str(),exact[i],r.rich,row,
            "Relationship",method,free_chart&&ds.size()==2?"Issue 225 A2 compared with centered exact A1 differences":"Issue 225 rebuilt values independent of action recorder",
            "abs(exact-R)<=atol+1e-6*max(abs(R),sigma) && refinement<=target/4",good,"reactive-caloric",id.c_str());
    }
}

#endif

// Frozen-composition oracle: P4 with complete declared A records supplies only
// the EOS/ideal mixing chain. Subtract those density standards and insert the
// already qualified anchored raw Gibbs columns; no equilibrium path is used.
class Frozen {
    const Case& c;
    double T,P;
    std::vector<double> x;
    HelmholtzModel baseline;
    ReferenceCalorics ref;
    std::array<double,5> baseline_g{};
    std::map<std::string,std::array<double,8>> cache;
    static ThermochemistryRecord complete(const Case& c) {
        auto records=c.records;for(auto& rows:records.component_intervals) if(rows.empty()) rows=records.component_intervals[0];return records;
    }
    double standard(const std::vector<Direction>& ds,std::size_t a,std::size_t b,bool correction) const {
        double result=0;
        std::function<void(std::size_t,std::size_t,std::size_t,std::optional<std::size_t>,double)> contract;
        contract=[&](std::size_t k,std::size_t ta,std::size_t pb,std::optional<std::size_t> param,double weight) {
            if(k==ds.size()) {
                for(std::size_t i=0;i<x.size();++i) {
                    double value=gibbs_contraction(ref,i,ta,pb,param);
                    if(correction&&!param&&pb==0) value-=baseline_g.at(ta);
                    result+=x[i]*weight*value;
                }return;
            }
            if(ds[k].temperature_k!=0) contract(k+1,ta+1,pb,param,weight*ds[k].temperature_k);
            if(ds[k].pressure_pa!=0) contract(k+1,ta,pb+1,param,weight*ds[k].pressure_pa);
            if(!param) for(std::size_t j=0;j<ds[k].active_parameters.size();++j)
                if(ds[k].active_parameters[j]!=0) contract(k+1,ta,pb,j,weight*ds[k].active_parameters[j]);
        };contract(0,a,b,std::nullopt,1.0);return result;
    }
public:
    // h,s,cp,v,u,cv,vT,vP, all at fixed x; increment actions retain the unit.
    Frozen(const Case& source,const CompiledSystem& compiled,const Sample& center):
        c(source),T(*center.solved.temperature),P(*center.solved.pressure),x(center.solved.mole_fractions),
        baseline(c.parameters,c.active,complete(c)),ref(reference_calorics(compiled,T,P,true,true,true,4)) {
        CHECK(ref.status==SolvedStateActionStatus::Available);
        using AD=CppAD::AD<double>;std::vector<AD> independent{T};CppAD::Independent(independent);
        const AD t=independent[0];const auto thermal=evaluate_ideal_correlation(c.records.component_intervals[0][0].correlation,t,1e5);
        std::vector<AD> dependent{thermal.h-t*thermal.s+R*t*CppAD::log(R*t/1e5)};
        CppAD::ADFun<double> tape(independent,dependent);baseline_g[0]=tape.Forward(0,std::vector<double>{T})[0];
        double factorial=1;
        for(std::size_t k=1;k<5;++k) {factorial*=k;baseline_g[k]=factorial*tape.Forward(k,std::vector<double>{k==1?1.0:0.0})[0];}
    }
    const std::array<double,8>& get(const std::vector<Direction>& ds) {
        std::string key;
        for(const auto& d:ds) {key+=recording_support::number(d.temperature_k)+","+recording_support::number(d.pressure_pa);
            for(double v:d.active_parameters) key+=","+recording_support::number(v);key+=';';}
        auto found=cache.find(key);if(found!=cache.end()) return found->second;
        std::vector<PressureActionDirection> pd;
        for(const auto& d:ds) {PressureActionDirection p;p.temperature=d.temperature_k;p.pressure=d.pressure_pa;p.composition.assign(x.size(),0);
            for(std::size_t a=0;a<c.active.size();++a) p.parameter.push_back({c.active[a].family,c.active[a].indices,d.active_parameters[a]});pd.push_back(p);}
        if(pd.empty()) pd.push_back({0,0,std::vector<double>(x.size(),0),{}});
        std::vector<PropertyObservable> observables={PropertyObservable::TotalEnthalpy,PropertyObservable::TotalEntropy,
            PropertyObservable::TotalIsobaricHeatCapacity,PropertyObservable::ResidualEnthalpy,PropertyObservable::ResidualInternalEnergy};
        if(ds.size()<3) {observables.push_back(PropertyObservable::MolarVolumeTemperatureDerivative);observables.push_back(PropertyObservable::MolarVolumePressureDerivative);}
        const auto batch=baseline.pressure_actions(T,P,x,DensityBranch::Liquid,observables,pd);
        std::vector<double> v;
        for(const auto& item:batch.results) {const auto& value=ds.empty()?item.value:item.action;CHECK(item.status==PropertyActionStatus::Available&&value);v.push_back(*value);}
        if(ds.size()>=3) v.resize(7,0);
        double mech=(v[3]-v[4]+(ds.empty()?R*T:ds.size()==1?R*ds[0].temperature_k:0))/P;
        double h=v[0]+standard(ds,0,0,true)-T*standard(ds,1,0,true);
        double s=v[1]-standard(ds,1,0,true),cp=ds.size()<3?v[2]-T*standard(ds,2,0,true):0;
        for(std::size_t j=0;j<ds.size();++j) {
            auto lower=ds;lower.erase(lower.begin()+j);
            mech-=ds[j].pressure_pa*(get(lower)[3]-standard(lower,0,1,false))/P;
            h-=ds[j].temperature_k*standard(lower,1,0,true);
            if(ds.size()<3) cp-=ds[j].temperature_k*standard(lower,2,0,true);
        }
        const double volume=mech+standard(ds,0,1,false);
        double u=h-P*volume;
        for(std::size_t j=0;j<ds.size();++j) {auto lower=ds;lower.erase(lower.begin()+j);u-=ds[j].pressure_pa*get(lower)[3];}
        std::array<double,8> out={h,s,cp,volume,u,0,ds.size()<3?v[5]+standard(ds,1,1,false):0,ds.size()<3?v[6]+standard(ds,0,2,false):0};
        return cache.emplace(key,out).first->second;
    }
};

inline void identities(const Case& c,const HelmholtzModel& model,const CompiledSystem& compiled,
                       Actions& actions,const std::vector<Quantity>& q,const Sample& center,const std::string& inputs) {
    const auto& value=actions.get({});
    const double T=*center.solved.temperature,P=*center.solved.pressure,n=center.solved.amounts[0];
    auto check=[&](std::string id,std::string quantity,double a,double b,double relative) {
        compare("reference-calorics/"+id+"/"+c.name,q[index(q,quantity)],a,b,inputs,relative);
    };
    check("F4/G","H",n*center.phase.g,value[index(q,"H")]-T*value[index(q,"S")],1e-9);
    check("F5/u","total-u_G",value[index(q,"total-u_G")],value[index(q,"total-h")]-P*value[index(q,"v_G")],1e-10);
    check("F5/cv-quotient","total-cv_G",value[index(q,"total-cv_G")],value[index(q,"total-cp")]+center.phase.quotient,1e-10);
    auto basis=directions(c);const auto& dt=basis[0].value;const auto& dp=basis[1].value;
    const auto t=actions.get({dt}),p=actions.get({dp}),tt=actions.get({dt,dt});
    check("F4/T-dS","H",T*t[index(q,"S")],t[index(q,"H")],1e-9);
    check("F4/TT-dS","H",T*tt[index(q,"S")]+dt.temperature_k*t[index(q,"S")],tt[index(q,"H")],1e-9);
    check("F4/P-dH","H",p[index(q,"H")],dp.pressure_pa*value[index(q,"V_G")]+T*p[index(q,"S")],1e-9);
    check("F4/Maxwell","S",p[index(q,"S")],-dp.pressure_pa*t[index(q,"V_G")]/dt.temperature_k,1e-9);
    if(!c.ionic&&!c.problem.reactions[0].reference->reference_pressure_pa) {
        const auto references=reaction_references(model,T,P,*c.problem.source_reference,
            {c.problem.reactions[0].stoichiometry},{*c.problem.reactions[0].reference},true,std::nullopt,1);
        CHECK(references[0].terminal_contraction_pressure_derivative);
        const double rtcp=R*T**references[0].terminal_contraction_pressure_derivative;
        const auto ref=reference_calorics(compiled,T,P,true,true,false,2);
        const double anchored_a=column(ref,0,1,0),anchored_b=column(ref,0,1,1);
        check("F4/gauge-A","v_G",anchored_a,0,1e-10);
        check("F4/gauge-B","v_G",anchored_b,-rtcp,1e-10);
        const double wrong=n/center.phase.rho+n*rtcp*(0.4*center.solved.mole_fractions[0]-0.2*center.solved.mole_fractions[1]);
        const double shift=wrong-n*center.phase.v;
        const double target=1e-12+1e-9*std::abs(n*center.phase.v);
        const bool discriminates=std::abs(shift)>10*target;
        recording_support::compare_custom(("reference-calorics/F4/gauge-discrimination/"+c.name).c_str(),"V_row_minus_V_anchored","m^3",shift,0,
            recording_support::fields({inputs,recording_support::field("V_row",wrong),recording_support::field("V_anchored",n*center.phase.v),
                recording_support::field("anchored_coefficients","0,-1"),recording_support::field("row_coefficients","0.4,-0.2"),recording_support::field("target",target)}),
            "Discrimination","conserved-moiety-gauge","Row representative cannot replace the anchored volume",
            "abs(V_row-V_anchored)>10*target",discriminates,"reactive-caloric");
    }
    Frozen frozen(c,compiled,center);
    std::vector<std::vector<Direction>> requests={{},{dt},{dp}};
    const auto& theta=basis[c.ionic?4:3].value;requests.push_back({theta});
    for(const auto& a:{dt,dp}) for(const auto& b:{dt,dp,theta}) requests.push_back({a,b});
    std::size_t row=0;
    for(const auto& ds:requests) {
        const auto v=frozen.get(ds);auto tds=ds,pds=ds;tds.push_back(dt);pds.push_back(dp);
        const auto tjet=frozen.get(tds),pjet=frozen.get(pds);
        double thermal=T*tjet[1],pressure=v[3]+T*pjet[1]/dp.pressure_pa;
        double duT=T*tjet[1]-P*v[6]*dt.temperature_k,duP=T*pjet[1]-P*v[7]*dp.pressure_pa;
        for(std::size_t j=0;j<ds.size();++j) {
            auto lower=ds;lower.erase(lower.begin()+j);auto lowerT=lower,lowerP=lower;lowerT.push_back(dt);lowerP.push_back(dp);
            const auto lt=frozen.get(lowerT),lp=frozen.get(lowerP),lv=frozen.get(lower);
            thermal+=ds[j].temperature_k*lt[1];pressure+=ds[j].temperature_k*lp[1]/dp.pressure_pa;
            duT+=ds[j].temperature_k*lt[1]-ds[j].pressure_pa*lv[6]*dt.temperature_k;
            duP+=ds[j].temperature_k*lp[1]-ds[j].pressure_pa*lv[7]*dp.pressure_pa;
        }
        const auto suffix=std::to_string(row++);
        check("F4/frozen-hT/"+suffix,"total-cp",tjet[0]/dt.temperature_k,v[2],1e-9);
        check("F4/frozen-TsT/"+suffix,"total-cp",thermal/dt.temperature_k,v[2],1e-9);
        check("F4/frozen-hP/"+suffix,"v_G",pjet[0]/dp.pressure_pa,pressure,1e-9);
        check("F4/frozen-Maxwell/"+suffix,"v_G",pjet[1]/dp.pressure_pa,-v[6],1e-9);
        check("F5/frozen-duT/"+suffix,"total-u_G",tjet[4],duT,1e-10);
        check("F5/frozen-duP/"+suffix,"total-u_G",pjet[4],duP,1e-10);
    }
    if(c.problem.reactions[0].reference->reference_pressure_pa) {
        check("F5/fixed-mechanical-u","total-u_G",value[index(q,"total-u_G")],value[index(q,"total-h")]-P/center.phase.rho,1e-10);
        check("F5/fixed-mechanical-v","v_G",value[index(q,"v_G")],1/center.phase.rho,1e-10);
        const auto state=model.state(T,center.phase.rho,center.solved.mole_fractions);
        check("F5/fixed-cv","total-cv_G",value[index(q,"total-cv_G")],center.phase.ideal[3]-R+state.residual_isochoric_heat_capacity,1e-10);
    }
}

inline bool inversion(const Case& c,const HelmholtzModel& model,const CompiledSystem& compiled,const Sample& center,
                      const std::vector<Quantity>& q,const std::string& inputs,double exact_cv) {
    const double T=*center.solved.temperature,P=*center.solved.pressure,targetVolume=center.phase.v;
    bool good=true;
    const auto u_at=[&](double delta) {
        double pressure=P;Local value;double residual=0;std::size_t iterations=0;
        const double limit=32*std::numeric_limits<double>::epsilon()*std::abs(targetVolume);
        for(;iterations<12;++iterations) {
            value=local(model,compiled,T+delta,pressure,center.solved.mole_fractions,center.phase.rho);
            residual=value.v-targetVolume;if(std::abs(residual)<=limit) break;
            const double trial=pressure-residual/value.vp;
            if(!(trial>0)&&std::isfinite(trial))
                throw OracleOutsideChart("delta_T_K="+recording_support::number(delta)+", trial_P_Pa="+recording_support::number(trial)
                    +", volume_residual="+recording_support::number(residual));
            pressure=trial;CHECK(pressure>0&&std::isfinite(pressure));
        }
        const bool closed=std::abs(residual)<=limit;
        const auto row=recording_support::fields({inputs,recording_support::field("delta_T_K",delta),recording_support::field("inverted_P_Pa",pressure),
            recording_support::field("inversion_residual",residual),recording_support::field("inversion_limit",limit),
            recording_support::field("inversion_iterations",iterations),recording_support::field("v_G_P_raw_at_inversion",value.vp),
            recording_support::field("v_G_at_inversion",value.v),recording_support::field("pressure_residual",value.diagnostics.pressure_residual),
            recording_support::field("site_residual",value.diagnostics.association_residual)});
        const std::string id="F5/"+c.name+"/inversion/"+recording_support::number(delta);
        recording_support::compare_custom(id.c_str(),"fixed_v_G_residual","m^3/mol",residual,0,row,"Closure","native-Newton-vG",
            "Frozen x local volume inversion uses exact raw v_G,P","abs(volume_residual)<=32*epsilon*abs(v_G_center)",closed,"reactive-caloric");return value.total[1];
    };
    good=scalar_ladder("F5/"+c.name+"/fixed-vG-cv",q[index(q,"total-cv_G")],exact_cv,0.2,inputs,1e-6,
        [&](double h){return (u_at(h)-u_at(-h))/(2*h);})&&good;
    for(std::size_t axis=0;axis<2;++axis) {
        const auto basis=directions(c);const auto d=basis[axis].value;
        const double exact=axis==0?-center.values[index(q,"S")]:d.pressure_pa*center.values[index(q,"V_G")];
        good=scalar_ladder("F4/"+c.name+(axis==0?"/G_T":"/G_P"),q[index(q,"H")],exact,0.2,inputs,1e-9,
            [&](double h){const auto plus=shifted(c,add(zero(c),d,h),q,center.solved),minus=shifted(c,add(zero(c),d,-h),q,center.solved);
                return (plus.solved.amounts[0]*plus.phase.g-minus.solved.amounts[0]*minus.phase.g)/(2*h);})&&good;
    }
    return good;
}

inline void product_identities(const Case& c,Actions& actions,const std::vector<NamedDirection>& ds,
                               const std::vector<Quantity>& q,const std::string& inputs) {
    std::vector<Direction> request;for(const auto& d:ds) request.push_back(d.value);
    const auto value=actions.get(request);const double P=c.problem.pressure.value;
    double expected=value[index(q,"total-h")]-P*value[index(q,"v_G")];
    for(std::size_t j=0;j<ds.size();++j) {
        auto lower=request;lower.erase(lower.begin()+j);
        expected-=request[j].pressure_pa*actions.get(lower)[index(q,"v_G")];
    }
    std::string label;for(const auto& d:ds) label+=d.name+"+";
    compare("F5/"+c.name+"/u-product/A"+std::to_string(ds.size())+"/"+label,
        q[index(q,"total-u_G")],value[index(q,"total-u_G")],expected,inputs,1e-10);
}

inline std::vector<std::vector<NamedDirection>> requests(const Case& c,const std::vector<NamedDirection>& basis) {
    std::vector<std::vector<NamedDirection>> result;
    for(const auto& a:basis) result.push_back({a});
    for(std::size_t i=0;i<basis.size();++i) for(std::size_t j=i;j<basis.size();++j)
        if(!(basis[i].parameter&&basis[j].parameter)) result.push_back({basis[i],basis[j]});
    const std::size_t theta=c.ionic?4:3,kappa=theta+1;
    auto combined=[&](std::size_t a,std::size_t b,const char* label) {
        NamedDirection d{label,std::min(basis[a].step,basis[b].step),add(basis[a].value,basis[b].value),basis[a].chart,
            basis[a].parameter||basis[b].parameter};
        d.chart.insert(d.chart.end(),basis[b].chart.begin(),basis[b].chart.end());return d;
    };
    const auto tf=combined(0,2,"T+feed"),pf=combined(1,2,"P+feed"),tp=combined(0,1,"T+P");
    result.push_back({tf});result.push_back({pf});result.push_back({tp});result.push_back({tf,tp});result.push_back({pf,tp});
    if(kappa<basis.size()) {
        const auto parameters=combined(theta,kappa,"epsilon+coefficient");
        result.push_back({parameters});result.push_back({tf,parameters});result.push_back({pf,parameters});
        result.push_back({tp,combined(0,theta,"T+epsilon")});
    }
    std::vector<NamedDirection> feeds={basis[2]};if(c.ionic) feeds.push_back(basis[3]);
    for(std::size_t a=0;a<feeds.size();++a) for(std::size_t b=a;b<feeds.size();++b) for(std::size_t d=b;d<feeds.size();++d)
        result.push_back({feeds[a],feeds[b],feeds[d]});
    return result;
}

#ifndef _WIN32
void run(bool ionic,bool fixed) {
    const auto c=input(ionic,fixed);HelmholtzModel model(c.parameters,c.active,c.records);
    const auto compiled=compile_problem(model,c.problem);const auto solved=solve_compiled(compiled);
    if(!solved.has_state()) throw std::runtime_error(c.name+" center: "+solved.message);
    const auto q=quantities(c.parameters.component_count());const auto basis=directions(c);
    const auto center=values(model,compiled,solved,q);const auto inputs=recording_support::fields({equilibrium_support::weighted_hessian_context(compiled),context(c,solved,center.phase)});
    Actions actions(compiled,solved,q);
    {   // #449: standalone first H actions request standard Gibbs degree two; a batch with Cp keeps degree three.
        // Neither demand may move H or Cp: each must equal its standalone request bit for bit.
        Observable h;h.kind=SolvedStateObservableKind::SystemEnthalpy;const auto cp=q[index(q,"total-cp")].observable;
        std::vector<Direction> ds;for(const auto& d:basis) ds.push_back(d.value);
        const auto alone=solved_state_first_actions(compiled,solved,{h},ds),cp_alone=solved_state_first_actions(compiled,solved,{cp},ds),
            both=solved_state_first_actions(compiled,solved,{h,cp},ds);
        for(std::size_t k=0;k<ds.size();++k) for(std::size_t i=0;i<2;++i) for(bool action:{false,true}) {
            const auto& a=(i?cp_alone:alone)[k].results[0];const auto& b=both[k].results[i];
            const auto& x=action?a.action:a.value;const auto& y=action?b.action:b.value;
            const bool same=a.status==SolvedStateActionStatus::Available&&b.status==a.status&&x&&y&&std::memcmp(&*x,&*y,sizeof(double))==0;
            const std::string id="449/"+c.name+"/H+cp-batch/"+basis[k].name+"/"+(i?"total-cp":"H")+(action?"/action":"/value");
            recording_support::compare_custom(id.c_str(),i?"total-cp":"H",i?"J/(mol K)":"J",y.value_or(NAN),x.value_or(NAN),inputs,
                "Relationship","standalone-request","Issue 449: the reference-column demand leaves H and Cp unchanged",
                "batch bits == standalone bits",same,"reactive-caloric");
        }
    }
    const auto exact=actions.get({});
    // Order 0 against the binary128 oracle, which shares only the EOS kernels with production; one tape per case.
    std::unique_ptr<WideNeutral> wide;std::unique_ptr<WideIonic> ionic_oracle;
    if(ionic) ionic_oracle=std::make_unique<WideIonic>(c.parameters,c.records,c.problem,center.solved.mole_fractions[1]);else wide=std::make_unique<WideNeutral>(c);
    const auto oracle=ionic?ionic_oracle->values(q,{},1,0,center.solved,center.phase.rho):wide->values(q,{},1,0,center);
    for(std::size_t i=0;i<q.size();++i) compare("F1/"+c.name+"/value/"+q[i].name,q[i],exact[i],double(oracle[i]),inputs,1e-10,"binary128-values");
    if(ionic||fixed) identities(c,model,compiled,actions,q,center,inputs);  // E-sys identities run in cheap()
    // A ladder whose every step leaves the chart records no comparison; it must still fail.
    CHECK(inversion(c,model,compiled,center,q,inputs,exact[index(q,"total-cv_G")]));
    for(const auto& request:requests(c,basis)) {
        ladder(c,actions,request,q,center,inputs,wide.get(),ionic_oracle.get());
        product_identities(c,actions,request,q,inputs);
        if(!ionic&&!fixed) {
            const bool zero_pressure=request.size()==1&&request[0].name=="P";
            const bool zero_theta=request.size()==1&&request[0].name=="epsilon";
            if(zero_pressure||zero_theta) {
                const auto sigma=scale(actions,request,q.size());const auto action=actions.get({request[0].value});
                for(std::size_t i=0;i<q.size();++i) if(q[i].observable.kind==SolvedStateObservableKind::PhaseAmount
                    ||q[i].observable.kind==SolvedStateObservableKind::PhaseMoleFraction) {
                    const double target=q[i].atol+64*std::numeric_limits<double>::epsilon()*sigma[i];const bool zero_ok=std::abs(action[i])<=target;
                    recording_support::compare_custom(("F1/"+c.name+"/analytic-zero/"+request[0].name+"/"+q[i].name).c_str(),q[i].name.c_str(),q[i].unit.c_str(),action[i],0,
                        recording_support::fields({inputs,recording_support::field("sigma",sigma[i]),recording_support::field("target",target)}),
                        "Relationship","identical-ethane-analytic-zero","Identical system-reference chemistry has no P or tied EOS composition response",
                        "abs(action)<=atol+64*epsilon*sigma",zero_ok,"reactive-caloric");
                }
            }
        }
    }
}

// P4 composition chain for U: x=n/N has an analytic derivative for every
// requested subset. Set partitions compose the curved feed path exactly.
inline std::vector<double> composition_action(const Case& c,const std::vector<Direction>& ds) {
    const auto& feed=c.problem.feed.values;const double N=std::accumulate(feed.begin(),feed.end(),0.0);
    std::vector<double> out(feed.size(),0);const std::size_t k=ds.size();
    if(k==0) {for(std::size_t i=0;i<feed.size();++i) out[i]=feed[i]/N;return out;}
    std::vector<double> Ns;for(const auto& d:ds) Ns.push_back(std::accumulate(d.feed_amounts_mol.begin(),d.feed_amounts_mol.end(),0.0));
    double factorial=1;for(std::size_t j=1;j<=k;++j) factorial*=j;
    for(std::size_t i=0;i<feed.size();++i) {
        double product=1;for(double v:Ns) product*=v;
        out[i]=std::pow(-1.0,int(k))*factorial*feed[i]*product/std::pow(N,int(k+1));
        for(std::size_t j=0;j<k;++j) {product=1;for(std::size_t l=0;l<k;++l) if(l!=j) product*=Ns[l];
            out[i]+=std::pow(-1.0,int(k-1))*(factorial/k)*ds[j].feed_amounts_mol[i]*product/std::pow(N,int(k));}
    }return out;
}
void unreactive() {
    const auto c=input(false,false,true);HelmholtzModel model(c.parameters,c.active,c.records);
    const auto compiled=compile_problem(model,c.problem);const auto solved=solve_compiled(compiled);
    CHECK(solved.has_state());const auto q=quantities(2);const auto center=values(model,compiled,solved,q);
    const auto inputs=recording_support::fields({equilibrium_support::weighted_hessian_context(compiled),context(c,solved,center.phase)});
    Actions actions(compiled,solved,q);auto basis=directions(c);WideNeutral wide(c);
    auto checks=requests(c,basis);checks.insert(checks.begin(),{});
    const std::vector<PropertyObservable> os={PropertyObservable::IdealEnthalpy,PropertyObservable::IdealInternalEnergy,PropertyObservable::IdealEntropy,
        PropertyObservable::IdealIsobaricHeatCapacity,PropertyObservable::IdealIsochoricHeatCapacity,
        PropertyObservable::TotalEnthalpy,PropertyObservable::TotalInternalEnergy,PropertyObservable::TotalEntropy,
        PropertyObservable::TotalIsobaricHeatCapacity,PropertyObservable::TotalIsochoricHeatCapacity};
    const char* names[]={"ideal-h","ideal-u_G","ideal-s","ideal-cp","ideal-cv_G","total-h","total-u_G","total-s","total-cp","total-cv_G"};
    std::size_t row=0;
    for(const auto& requested:checks) {
        std::vector<Direction> ds;for(const auto& d:requested) ds.push_back(d.value);
        const auto exact=actions.get(ds);std::vector<double> oracle(os.size(),0);
        std::vector<std::vector<std::size_t>> partition;
        std::function<void(std::size_t)> visit=[&](std::size_t k) {
            if(k==ds.size()) {
                std::vector<PressureActionDirection> pd;
                for(const auto& block:partition) {
                    std::vector<Direction> subset;for(auto i:block) subset.push_back(ds[i]);
                    PressureActionDirection d;d.composition=composition_action(c,subset);
                    if(block.size()==1) {const auto& w=ds[block[0]];d.temperature=w.temperature_k;d.pressure=w.pressure_pa;
                        for(std::size_t a=0;a<c.active.size();++a) d.parameter.push_back({c.active[a].family,c.active[a].indices,w.active_parameters[a]});}
                    pd.push_back(d);
                }
                if(pd.empty()) pd.push_back({0,0,std::vector<double>(2,0),{}});
                const auto p4=model.pressure_actions(c.problem.temperature.value,c.problem.pressure.value,center.solved.mole_fractions,DensityBranch::Liquid,os,pd,center.phase.rho);
                for(std::size_t i=0;i<os.size();++i) {const auto& v=ds.empty()?p4.results[i].value:p4.results[i].action;
                    CHECK(p4.results[i].status==PropertyActionStatus::Available&&v);oracle[i]+=*v;}return;
            }
            for(std::size_t i=0;i<partition.size();++i) {partition[i].push_back(k);visit(k+1);partition[i].pop_back();}
            partition.push_back({k});visit(k+1);partition.pop_back();
        };visit(0);
        const auto id=std::to_string(row++);
        for(std::size_t i=0;i<os.size();++i) compare("F3/U/P4/"+id+"/"+names[i],q[index(q,names[i])],exact[index(q,names[i])],oracle[i],inputs,1e-10,"native-P4-curved-composition-chain");
        product_identities(c,actions,requested,q,inputs);
        if(!requested.empty()) ladder(c,actions,requested,q,center,inputs,&wide);
    }
    const auto& base=actions.get({});
    compare("F3/U/v_G",q[index(q,"v_G")],base[index(q,"v_G")],1/center.phase.rho,inputs,1e-10);
}
#endif

inline void dependency_checks() {
    auto c = input(false, true);
    HelmholtzModel model(c.parameters, c.active, c.records);
    const auto compiled = compile_problem(model, c.problem);
    const auto solved = solve_compiled(compiled);
    CHECK(solved.has_state());
    const auto q = quantities(2);
    const auto ref = reference_calorics(compiled, 260.0, 3e6, true, true, false, 2);
    const double hB = column(ref, 0, 0, 1) - 260.0 * column(ref, 1, 0, 1) + R * 260.0;
    const double sB = -column(ref, 1, 0, 1) + R * (1 + std::log(R * 260.0 / 1e5));
    const double cpB = -260.0 * column(ref, 2, 0, 1) + R;
    auto check = [&](const std::string& id, double actual, double expected, double atol) {
        Quantity quantity{id, "J/mol", atol, {}};
        compare("reference-calorics/dependencies/" + id, quantity, actual, expected, "{}", 1e-10);
    };
    const auto zero_direction = zero(c);
    const std::vector<Observable> outputs = {q[index(q, "total-h")].observable,
        q[index(q, "total-cp")].observable,
        {SolvedStateObservableKind::PhaseProperty, 0, 0, PropertyObservable::ResidualEnthalpy},
        q[index(q, "x-1")].observable};
    const auto baseline = solved_state_actions(compiled, solved, {outputs, {zero_direction}});
    for (bool inconsistent_cp : {false, true}) {
        auto records = c.records;
        records.component_intervals[1] = {{200, 400, true, true,
            {IdealPolynomial{{cpB + double(inconsistent_cp)}, 260.0, hB, sB}, 1e5}}};
        HelmholtzModel anchored(c.parameters, c.active, records);
        const auto complete = compile_problem(anchored, c.problem);
        const auto state = solve_compiled(complete);
        CHECK(state.has_state());
        auto source = zero_direction;
        source.reaction_coefficients[0][0] = 1.0;
        const auto values = solved_state_actions(complete, state, {outputs, {zero_direction}});
        const std::string record = inconsistent_cp ? "inconsistent-cp-record/" : "consistent-cp-record/";
        check(record + "anchored-h", *values.results[0].value, *baseline.results[0].value, 1e-5);
        CHECK(values.results[2].status == SolvedStateActionStatus::Available);
        check(record + "independent-residual", *values.results[2].value, *baseline.results[2].value, 1e-5);
        CHECK(values.results[3].status == SolvedStateActionStatus::Available);
        if (inconsistent_cp) {
            CHECK(values.results[1].status == SolvedStateActionStatus::ReferenceUnavailable && !values.results[1].value);
        } else {
            check("consistent-cp", *values.results[1].value, *baseline.results[1].value, 1e-7);
            const auto action = solved_state_actions(complete, state, {outputs, {source}});
            CHECK(action.results[0].status == SolvedStateActionStatus::ReferenceUnavailable
                  && action.results[0].value && !action.results[0].action);
            check("source-retains-center-h", *action.results[0].value, *baseline.results[0].value, 1e-5);
            CHECK(action.results[2].status == SolvedStateActionStatus::Available);
            CHECK(action.results[3].status == SolvedStateActionStatus::Available && action.results[3].action);
            CHECK(std::abs(*action.results[3].action) > 1e-11);
        }
    }
    auto endpoint_records = c.records;
    auto& endpoint_interval = endpoint_records.component_intervals[0][0];
    const auto endpoint_thermal = evaluate_ideal_correlation(endpoint_interval.correlation, 260.0, 1e5);
    // Rebase the same linear Cp law: its thermal anchor must belong to the interval.
    endpoint_interval.upper = 260.0;
    endpoint_interval.correlation.form = IdealPolynomial{{endpoint_thermal.cp, 0.12},
        260.0, endpoint_thermal.h, endpoint_thermal.s};
    HelmholtzModel endpoint(c.parameters, c.active, endpoint_records);
    const auto endpoint_compiled = compile_problem(endpoint, c.problem);
    const auto endpoint_state = solve_compiled(endpoint_compiled);
    CHECK(endpoint_state.has_state());
    auto basis = directions(c);
    for (std::size_t k = 0; k < 3; ++k) {
        const auto& d = basis[k].value;
        const auto action = solved_state_actions(endpoint_compiled, endpoint_state, {outputs, {d}});
        check("endpoint-center-h/" + basis[k].name, *action.results[0].value, *baseline.results[0].value, 1e-5);
        const auto interior = solved_state_actions(compiled, solved, {outputs, {d}});
        if (d.temperature_k != 0.0)
            CHECK(action.results[0].status == SolvedStateActionStatus::ThermalBoundary && !action.results[0].action);
        else {
            CHECK(action.results[0].status == SolvedStateActionStatus::Available && action.results[0].action);
            check("endpoint-isothermal-action/" + basis[k].name, *action.results[0].action, *interior.results[0].action, 1e-5);
        }
        CHECK(action.results[2].status == SolvedStateActionStatus::Available);
    }
}

inline void volume_partial_values() {
    const auto c=input(false,false,true);HelmholtzModel model(c.parameters,c.active,c.records);
    const std::vector<PropertyObservable> os={PropertyObservable::MolarVolumeTemperatureDerivative,
        PropertyObservable::MolarVolumePressureDerivative};
    PressureActionDirection zero_direction{0,0,{0,0},{}};
    const auto center=model.pressure_actions(260,3e6,{1,0},DensityBranch::Liquid,os,{zero_direction});
    auto insertion=zero_direction;insertion.composition={-1,1};
    // A neutral insertion is its exact infinite-dilution limit (#73): f(0) = 2 f(eps) - f(2 eps) + O(eps^2).
    const auto at=[&](double eps){return model.pressure_actions(260,3e6,{1-eps,eps},DensityBranch::Liquid,os,{insertion});};
    const auto refused=at(0),one=at(1e-7),two=at(2e-7);
    for(std::size_t k=0;k<os.size();++k) {
        CHECK(center.results[k].value && refused.results[k].value);
        CHECK(refused.results[k].status==PropertyActionStatus::Available && refused.results[k].action);
        const double limit=2*one.results[k].action.value()-two.results[k].action.value();
        // Identical components: the limit is zero; the floor is 1e-10 of the quantity itself.
        CHECK(std::abs(*refused.results[k].action-limit)<=std::max(1e-9*std::abs(limit),1e-10*std::abs(*center.results[k].value)));
        Quantity q{k==0?"vT":"vP",k==0?"m^3/(mol K)":"m^3/(mol Pa)",k==0?1e-12:1e-18,{}};
        compare("reference-calorics/dependencies/retained-"+q.name,q,*refused.results[k].value,*center.results[k].value,
            "T_K:260;P_Pa:3000000;x:1,0;insertion:-1,1",1e-10);
    }
}
inline void isothermal_row_columns() {
    auto c=input(true,true);c.active={{ParameterFamily::BornDiameter,{1}}};
    // A steep declared dielectric slope no longer has a finite-dilution truncation gate.
    c.parameters.permittivity_term_component={0};c.parameters.permittivity_term_power={1};
    c.parameters.permittivity_term_coefficient={1e10};
    HelmholtzModel model(c.parameters,c.active);const auto compiled=compile_problem(model,c.problem);
    const auto solved=solve_compiled(compiled);CHECK(solved.has_state());
    const auto diagnostics=reaction_references(model,298.15,1e5,*c.problem.source_reference,
        {c.problem.reactions[0].stoichiometry},{*c.problem.reactions[0].reference},false,std::nullopt,1,true);
    CHECK(diagnostics[0].terminal_contraction_parameter_derivative[0]);
    CHECK(diagnostics[0].terminal_contraction_temperature_parameter_derivative[0]);
    const auto lower_entry=std::find_if(diagnostics[0].lower_derivatives.begin(),diagnostics[0].lower_derivatives.end(),
        [](const auto& d){return d.temperature_order==1 && d.pressure_order==0 && d.parameter==0;});
    CHECK(lower_entry!=diagnostics[0].lower_derivatives.end());
    const auto& lower=*lower_entry;
    const auto inputs=recording_support::fields({equilibrium_support::weighted_hessian_context(compiled),
        recording_support::field("dielectric_T_slope_per_K",1e10),
        recording_support::field("C_theta",*diagnostics[0].terminal_contraction_parameter_derivative[0]),
        recording_support::field("C_Ttheta",*lower.value),
        recording_support::field("site_residual",lower.diagnostics.association_residual)});
    Quantity theta{"C_theta","1/angstrom",1e-10,{}};
    compare("reference-calorics/dependencies/lower-column/available-theta",theta,
        *diagnostics[0].terminal_contraction_parameter_derivative[0],46.619570441218094,inputs,1e-10);
    const auto shifted_rows=[&](double delta) {
        auto p=c.parameters;auto problem=c.problem;p.born_diameter_angstrom[1]+=delta;problem.feed.values[0]+=delta;
        HelmholtzModel shifted_model(p,c.active);const auto shifted_compiled=compile_problem(shifted_model,problem);
        return shifted_compiled.evaluate(solved.coordinates).values;
    };
    for(bool warm:{false,true}) {
        const auto prepared=compiled.prepare(solved.coordinates);
        const std::vector<double> z(prepared.point.size(),0);
        const std::vector<double> control{0,0,1,0,0,1};
        if(warm) (void)compiled.row_series(prepared,{prepared.point,z},control,1);
        const auto second=compiled.row_series(prepared,{prepared.point,z,z},control,2);
        for(std::size_t row=0;row<3;++row) {
            Quantity q{"chemical-row-"+std::to_string(row),"1",1e-11,{}};
            CHECK(scalar_ladder("reference-calorics/dependencies/isothermal-row/"+std::to_string(warm)+"/"+q.name,
                q,second.rows[1][row],0.001,inputs,1e-6,
                [&](double h){return (shifted_rows(h)[row]-shifted_rows(-h)[row])/(2*h);}));
        }
    }
}
#ifndef _WIN32
inline Case boundary_input(bool nonidentical,bool free_temperature,bool fixed) {
    auto c=input(false,fixed);
    c.name=std::string(nonidentical?"N":"E")+"-bubble-"+(free_temperature?"T":"P")+(fixed?"-fix":"-sys");
    if(nonidentical) {
        c.parameters=properties_support::neutral_ternary();
        c.active={{ParameterFamily::DispersionEnergyOverK,{1}}};
        c.records.component_intervals.resize(3);
        c.records.component_intervals[2]=c.records.component_intervals[0];
        c.problem.feed.values={2,0.2,0.2};
        c.problem.source_reference->solute_molality_ratios={0,1,1};
        c.problem.reactions[0].name="A=2B";
        c.problem.reactions[0].stoichiometry={-1,2,0};
        c.problem.reactions[0].correlation=ReactionLogPolynomial{std::log(0.3),0,0,0,300,250,350,"analytic-source"};
    }
    c.problem.phases.push_back(equilibrium_support::phase("V",PhaseKind::Vapor,AmountSpec::pinned(0)));
    if(nonidentical) {
        c.problem.phases[0].composition_guess={0.21145628710236675,0.7334775180293439,0.05506619486828941};
        c.problem.phases[1].composition_guess={0.5275400084331457,0.46022444490546427,0.012235546661389953};
    }
    c.problem.temperature=free_temperature?FreeValue::variable(nonidentical?260:270):FreeValue::fixed(nonidentical?260:265);
    c.problem.pressure=free_temperature?FreeValue::fixed(nonidentical?3.953695785877087e6:1.952823e6)
        :FreeValue::variable(nonidentical?3.953695785877087e6:2e6);
    return c;
}

inline void pure_boundary(const Case& c,Actions& actions,const Sample& center,const std::vector<Quantity>& q,
                          const std::vector<NamedDirection>& basis,const std::string& inputs) {
    const HelmholtzModel pure(properties_support::ethane());
    const double P=*center.solved.pressure;
    const auto states=[&](double T) {
        return std::array<State,2>{pure.state_at_pressure(T,P,{1},DensityBranch::Liquid),
                                   pure.state_at_pressure(T,P,{1},DensityBranch::Vapor)};
    };
    const auto fugacity=[&](double T) {
        const auto s=states(T);
        CHECK(s[0].molar_density>s[1].molar_density);
        return require_component_property(s[0].log_fugacity_coefficient->at(0))
            -require_component_property(s[1].log_fugacity_coefficient->at(0));
    };
    double lo=250,hi=280,flo=fugacity(lo),fhi=fugacity(hi);
    CHECK(flo*fhi<0);
    for(std::size_t k=0;k<80&&hi-lo>1e-9;++k) {
        const double mid=(lo+hi)/2,f=fugacity(mid);
        if(f*flo>0) {lo=mid;flo=f;} else hi=mid;
    }
    const double T=(lo+hi)/2;const auto s=states(T);
    compare("F6/"+c.name+"/pure-saturation",q[index(q,"T")],*center.solved.temperature,T,
  recording_support::fields({inputs,recording_support::field("pure_fugacity_residual",fugacity(T)),
      recording_support::field("pure_T_bracket",std::vector<double>{lo,hi})}),1e-8,"native-pure-fugacity-bisection");
    if(c.problem.temperature.free) {
        const double clapeyron=T*(1/s[1].molar_density-1/s[0].molar_density)
            /(s[1].residual_enthalpy-s[0].residual_enthalpy);
        auto quantity=q[index(q,"T")];quantity.name="dT_dP";quantity.unit="K/Pa";quantity.atol/=1e5;
        compare("F6/"+c.name+"/Clapeyron",quantity,actions.get({basis[0].value})[index(q,"T")]/1e5,
      clapeyron,inputs,1e-8,"pure-Properties-Clapeyron");
        const auto sigma=scale(actions,{basis[1]},q.size());
        const auto k=index(q,"T");const double tangent=actions.get({basis[1].value})[k];
        const bool good=std::abs(tangent)<=q[k].atol+64*std::numeric_limits<double>::epsilon()*sigma[k];
        recording_support::compare_custom(("F6/"+c.name+"/pure-feed-zero").c_str(),"DT_feed","K",tangent,0,
            inputs,"Relationship","identical-species","identical ethane saturation is independent of feed",
            "abs(DT)<=atol+64*epsilon*sigma",good,"reactive-caloric");
    } else {
        const double clapeyron=(s[1].residual_enthalpy-s[0].residual_enthalpy)
            /(T*(1/s[1].molar_density-1/s[0].molar_density));
        auto quantity=q[index(q,"P")];quantity.name="dP_dT";quantity.unit="Pa/K";
        compare("F6/"+c.name+"/Clapeyron",quantity,actions.get({basis[0].value})[index(q,"P")],
      clapeyron,inputs,1e-8,"pure-Properties-Clapeyron");
    }
}

void boundary(bool nonidentical,bool free_temperature,bool fixed) {
    const auto c=boundary_input(nonidentical,free_temperature,fixed);
    HelmholtzModel model(c.parameters,c.active,c.records);
    const auto solve=[&](const EquilibriumProblem& problem,const char* label) {
        core::Capture capture;const auto start=std::chrono::steady_clock::now();
        const auto compiled=compile_problem(model,problem);const auto solved=solve_compiled(compiled);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        const auto references=std::count_if(capture.trace().spans.begin(),capture.trace().spans.end(),
            [](const auto& span){return std::string(span.name)=="reaction_references";});
        recording_support::observe(("F6/"+c.name+"/"+label+"-work").c_str(),"wall_time","s",seconds,std::nullopt,
            recording_support::fields({recording_support::field("iterations",solved.iterations),
                recording_support::field("evaluations",solved.evaluations),recording_support::field("reference_evaluations",references),
                recording_support::field("solver_status",solved.solver_status),recording_support::field("solver_message",solved.message),
                recording_support::field("residual",solved.residual_norm()),recording_support::field("T",solved.temperature.value_or(NAN)),
                recording_support::field("P",solved.pressure.value_or(NAN)),recording_support::field("roots",solved.molar_densities),
                recording_support::field("seed_T",problem.temperature.value),
                recording_support::field("seed_P",problem.pressure.value)}),
            "Observation","native-one-driver","Equivalent work; no performance ceiling","reactive-caloric");
        if(!solved.success||!solved.has_state()) throw std::runtime_error(c.name+" "+label+": "+solved.message);
        return solved;
    };
    const auto solved=solve(c.problem,"free-intensive");const auto compiled=compile_problem(model,c.problem);
    CHECK(solved.topology_event==TopologyEvent::None);
    CHECK(solved.molar_densities[0]>solved.molar_densities[1]);
    auto q=quantities(c.parameters.component_count());
    auto vapor=q;vapor.resize(vapor.size()-4);
    vapor.erase(vapor.begin(),vapor.begin()+2); // T/P are common controls/outputs.
    for(auto& quantity:vapor) {quantity.name+="-V";quantity.observable.phase=1;}
    q.insert(q.end()-4,vapor.begin(),vapor.end());
    const auto center=values(model,compiled,solved,q);const auto inputs=recording_support::fields({equilibrium_support::weighted_hessian_context(compiled),context(c,solved,center.phase)});
    for(std::size_t phase=0;phase<center.phases.size();++phase) {
        const auto& local=center.phases[phase];const auto id="F6/"+c.name+"/volume-chart/"+std::to_string(phase);
        recording_support::observe(id.c_str(),"v_G","m^3/mol",local.v,std::nullopt,
            recording_support::fields({context(c,solved,local),recording_support::field("phase",phase)}),
            "Observation","native-local-chart","Regular volume chart; stable signs and quotient sensitivities recorded","reactive-caloric");
    }
    Actions actions(compiled,solved,q);
    auto basis=directions(c);
    basis.erase(basis.begin()+(free_temperature?0:1));
    std::vector<std::vector<NamedDirection>> combinations;
    for(const auto& a:basis) combinations.push_back({a});
    for(std::size_t i=0;i<basis.size();++i) for(std::size_t j=i;j<basis.size();++j)
        if(!(basis[i].parameter&&basis[j].parameter)) combinations.push_back({basis[i],basis[j]});
    const auto combine=[&](std::size_t a,std::size_t b,const char* name) {
        auto d=basis[a];d.name=name;d.value=add(d.value,basis[b].value);d.step=std::min(d.step,basis[b].step);
        d.chart.insert(d.chart.end(),basis[b].chart.begin(),basis[b].chart.end());d.parameter|=basis[b].parameter;return d;
    };
    const auto controls=combine(0,1,free_temperature?"P+feed":"T+feed");
    const auto parameters=combine(2,3,"epsilon+coefficient");
    combinations.push_back({controls});combinations.push_back({parameters});
    combinations.push_back({controls,controls});combinations.push_back({controls,parameters});
    auto mixed=combine(0,2,free_temperature?"P+epsilon":"T+epsilon");
    combinations.push_back({controls,mixed});
    for(const auto& ds:combinations) ladder(c,actions,ds,q,center,inputs);
    const auto intensive=index(q,free_temperature?"T":"P");
    if(nonidentical) {
        for(const auto& ds:std::vector<std::vector<Direction>>{{basis[1].value},{basis[1].value,basis[2].value},
                                                              {basis[0].value,basis[2].value}}) {
            const auto response=actions.get(ds)[intensive];
            recording_support::compare_custom(("F6/"+c.name+"/nonzero-intensive-"+std::to_string(ds.size())
                +(ds.front().feed_amounts_mol[0]!=0?"-feed":"-control")).c_str(),"solved-intensive-response",
                free_temperature?"K":"Pa",response,0,inputs,"Relationship","nonidentical-witness",
                "N exposes actual feed and mixed parameter motion","finite(response) && response!=0",
                std::isfinite(response)&&response!=0,"reactive-caloric");
        }
    } else pure_boundary(c,actions,center,q,basis,inputs);
    auto driven=c.problem;driven.temperature=FreeValue::fixed(*solved.temperature);
    driven.pressure=FreeValue::fixed(*solved.pressure);driven.phases.resize(1);
    const auto fixed_center=solve(driven,"fixed-intensive");
    const auto fixed_compiled=compile_problem(model,driven);
    const auto liquid_quantities=quantities(c.parameters.component_count());
    const auto fixed_values=values(model,fixed_compiled,fixed_center,liquid_quantities);
    for(const auto& name:std::vector<std::string>{"rho","x-0","total-h"}) {
        const auto k=index(liquid_quantities,name);
        compare("F6/"+c.name+"/fixed-center/"+name,liquid_quantities[k],center.values[index(q,name)],
            fixed_values.values[k],inputs,1e-10,"same-liquid-state-on-fixed-chart");
    }
}

#endif

void cheap() {
    dependency_checks();
    volume_partial_values();
    isothermal_row_columns();
    const auto c=input(false,false);HelmholtzModel model(c.parameters,c.active,c.records);
    const auto compiled=compile_problem(model,c.problem);const auto solved=solve_compiled(compiled);CHECK(solved.has_state());
    const auto q=quantities(2);const auto center=values(model,compiled,solved,q);
    Actions actions(compiled,solved,q);identities(c,model,compiled,actions,q,center,context(c,solved,center.phase));
}
#include "reactive_cp_checks.hpp"
} // namespace reactive_caloric_checks
