#pragma once
#include "equilibrium_checks_shared.hpp"

// The issue-224 caloric and critical-envelope checks share the recorded comparison and the
// derivative ladder.
namespace caloric_specification_checks {
using namespace epcsaft;
using namespace epcsaft::equilibrium;
using namespace equilibrium_checks;
inline void compare(const std::string& id, double actual, double reference, double scale, double tolerance,
                    const std::string& inputs = "") {
    const double error = std::abs(actual-reference)/scale;
    recording_support::compare_custom(id.c_str(),"caloric-specification","SI",actual,reference,
        recording_support::fields({inputs,recording_support::field("scale",scale),
            recording_support::field("error",error),recording_support::field("criterion",tolerance)}),
        "Relationship","issue-224","PH/PS value or exact derivative falsifier","positive quantity scale",
        std::isfinite(error)&&error<=tolerance,"caloric-specification",id.c_str());
}
inline void ladder(const CompiledSystem& c,const std::vector<double>& y,const std::vector<double>& scale,const std::string& id,
                   const std::vector<double>& objectives={0.,0.5,1.}) {
    const std::size_t n=y.size(),m=c.rows.size();
    std::vector<std::vector<double>> weights(1,std::vector<double>(m));
    for(std::size_t j=0;j<m;++j)weights[0][j]=(j%2?-1:1)*(0.3+0.07*j);
    for(std::size_t r=0;r<m;++r)if(c.rows[r].kind>=RowKind::Enthalpy){weights.emplace_back(m,0.);weights.back()[r]=1;}
    recording_support::observe((id+"/F5/point").c_str(),"callback-coordinates","chart",y[c.primal_size-1],std::nullopt,
        recording_support::fields({recording_support::field("point",y),recording_support::field("feed_mol",c.problem.feed.values)}),
        "Observation","off-solution-chart","declared callback coordinates and independent chemical duals","caloric-specification");
    const auto J=c.jacobian(y);std::vector<std::vector<double>> H;
    for(const auto& w:weights)H.push_back(c.hessian(y,w));  // feasibility objective is constant for every factor
    const std::size_t sets=weights.size()*objectives.size();
    for(std::size_t column=0;column<n;++column) {
        std::vector<double> bestJ(m,1e100),previous_change(sets,0),largest(sets,0),noise(sets,0);
        std::vector<std::vector<double>> bestH(sets,std::vector<double>(n,1e100)),previousH(sets);std::vector<bool> refinement(sets,false);
        std::vector<std::array<double,2>> atJ(m);std::vector<std::vector<std::array<double,2>>> atH(sets,std::vector<std::array<double,2>>(n));  // exact, reference at the best step
        for(int rung=0;rung<5;++rung) {
            const double step=std::ldexp(0.002,-rung)*scale[column];auto plus=y,minus=y;plus[column]+=step;minus[column]-=step;
            const auto fp=c.evaluate(plus).values,fm=c.evaluate(minus).values;
            for(std::size_t row=0;row<m;++row) {
                const double exact=J[row*n+column]*scale[column],ref=(fp[row]-fm[row])/(2*step)*scale[column];
                const double error=std::abs(exact-ref)/std::max(1.,std::abs(ref));
                if(error<bestJ[row]){bestJ[row]=error;atJ[row]={exact,ref};}
            }
            for(std::size_t k=0;k<sets;++k) {
                const auto& w=weights[k/objectives.size()];const double objective=objectives[k%objectives.size()];
                const auto gp=weighted_gradient(c,plus,w,objective),gm=weighted_gradient(c,minus,w,objective);
                double change=0,floor=0;std::vector<double> refs(n);
                for(std::size_t row=0;row<n;++row) {
                    refs[row]=(gp[row]-gm[row])/(2*step);const double factor=scale[row]*scale[column];
                    const double exact=H[k/objectives.size()][row*n+column]*factor,error=std::abs(exact-refs[row]*factor)/std::max(1.,std::abs(refs[row]*factor));
                    if(error<bestH[k][row]){bestH[k][row]=error;atH[k][row]={exact,refs[row]*factor};}
                    if(rung)change=std::max(change,std::abs(refs[row]-previousH[k][row])*factor);
                    floor=std::max(floor,std::numeric_limits<double>::epsilon()*(std::abs(gp[row])+std::abs(gm[row]))*factor/step);
                }
                if(rung) {
                    largest[k]=std::max(largest[k],change);noise[k]=std::max(noise[k],floor);
                    const double ratio=change?previous_change[k]/change:0;
                    if(rung>1&&change>32*floor&&ratio>3.5&&ratio<4.5)refinement[k]=true;
                    previous_change[k]=change;
                }
                previousH[k]=refs;
            }
        }
        const auto entry=[&](const std::array<double,2>& at) {
            return recording_support::fields({recording_support::field("exact",at[0]),recording_support::field("reference",at[1]),recording_support::field("coordinate_scale",scale[column])});
        };
        for(std::size_t row=0;row<m;++row)compare(id+"/F5/J/"+std::to_string(column)+"/"+std::to_string(row),bestJ[row],0,1,1e-6,entry(atJ[row]));
        for(std::size_t k=0;k<sets;++k) {
            const std::string label=id+"/F5/H/"+std::to_string(k/objectives.size())+"/"+std::to_string(objectives[k%objectives.size()])+"/"+std::to_string(column);
            for(std::size_t row=0;row<n;++row)compare(label+"/"+std::to_string(row),bestH[k][row],0,1,1e-5,entry(atH[k][row]));
            recording_support::observe((label+"/refinement").c_str(),"second-order-refinement","chart",largest[k],noise[k],
                recording_support::fields({recording_support::field("refined",refinement[k]?1:0)}),"Observation","centered-refinement","truncation before the rounding floor","caloric-specification");
            CHECK(refinement[k]||largest[k]<=32*noise[k]);
        }
    }
}
}  // namespace caloric_specification_checks
