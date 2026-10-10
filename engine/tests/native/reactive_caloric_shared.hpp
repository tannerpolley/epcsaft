#pragma once
#include "equilibrium_checks_shared.hpp"

// S2 and the reactive-caloric checks share the frozen Gibbs contraction and the scalar ladder.
namespace reactive_caloric_checks {
using namespace epcsaft;
using namespace epcsaft::equilibrium;
using binary128_support::Quantity;

inline double gibbs_contraction(const ReferenceCalorics& ref,std::size_t component,
                               std::size_t a,std::size_t b,std::optional<std::size_t> parameter) {
    for(const auto& term:ref.gibbs) if(term.temperature_order==a&&term.pressure_order==b&&term.parameter==parameter) {
        if(term.column.status!=SolvedStateActionStatus::Available) throw std::runtime_error(term.column.message);
        return term.column.values.at(component);
    }
    throw std::runtime_error("missing frozen Gibbs contraction");
}

struct OracleOutsideChart : std::runtime_error { using std::runtime_error::runtime_error; };
inline bool scalar_ladder(const std::string& id,const Quantity& q,double exact,double step,const std::string& inputs,
                          double relative,const std::function<double(double)>& centered) {
    auto evaluate=[&](double h)->std::optional<double> {
        try { return centered(h); }
        catch(const OracleOutsideChart& error) {
            const auto key=id+"/outside-chart/"+recording_support::number(h);
            recording_support::observe(key.c_str(),q.name.c_str(),q.unit.c_str(),exact,std::nullopt,
                recording_support::fields({inputs,recording_support::field("requested_step",h),
                    recording_support::field("resolved",0),recording_support::field("cause",error.what())}),
                "Refinement","native-independent-value-difference","No derivative reference outside the positive-pressure chart; halve the prescribed step","reactive-caloric");
            return std::nullopt;
        }
    };
    auto coarse_value=evaluate(step);bool good=false;
    for(std::size_t rung=0;rung<6;++rung) {
        const double h=std::ldexp(step,-int(rung+1));const auto fine_value=evaluate(h);
        if(!coarse_value||!fine_value) {coarse_value=fine_value;continue;}
        const double coarse=*coarse_value,fine=*fine_value,rich=(4*fine-coarse)/3,refinement=std::abs(fine-coarse)/3;
        const double sigma=std::abs(exact),target=q.atol+relative*std::max(std::abs(rich),sigma);
        const bool resolved=refinement<=target/4;good=resolved&&std::abs(exact-rich)<=target;
        const auto row=recording_support::fields({inputs,recording_support::field("coarse_step",2*h),recording_support::field("fine_step",h),
            recording_support::field("D_h",coarse),recording_support::field("D_h2",fine),recording_support::field("Richardson",rich),
            recording_support::field("refinement",refinement),recording_support::field("sigma",sigma),recording_support::field("target",target),
            recording_support::field("error",std::abs(exact-rich)),recording_support::field("error_over_target",std::abs(exact-rich)/target),
            recording_support::field("resolved",resolved)});
        const auto key=id+"/rung-"+std::to_string(rung);
        if(!resolved&&rung<5) recording_support::observe(key.c_str(),q.name.c_str(),q.unit.c_str(),exact,rich,row,
            "Refinement","native-independent-value-difference","Coarse unresolved step retained and refined","reactive-caloric");
        else recording_support::compare_custom(key.c_str(),q.name.c_str(),q.unit.c_str(),exact,rich,row,
            "Relationship","native-independent-value-difference","Independent Gibbs or fixed-volume verification",
            "abs(exact-R)<=atol+relative*max(abs(R),sigma) && refinement<=target/4",good,"reactive-caloric");
        if(resolved) break;coarse_value=fine_value;
    }
    return good;
}

}  // namespace reactive_caloric_checks
