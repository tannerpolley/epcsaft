#pragma once
// Includes and helpers shared by the equilibrium_checks translation units.
#include "quad_base.hpp"
#include "properties/helmholtz_tape.hpp"
#include "equilibrium/gibbs_minimizer.hpp"
#include "equilibrium/stoichiometry.hpp"
#include "properties/thermochemistry.hpp"
#include "equilibrium_support.hpp"
#include "core/numerical_policy.hpp"
#include "core/physical_constants.hpp"
#include "recording_support.hpp"
#include "binary128_support.hpp"
#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iterator>
#include <memory>
#include <limits>
#include <string>
#include <stdexcept>
#include <tuple>
#include <map>
#include <numeric>

namespace equilibrium_checks {

inline double convergence_limit(const epcsaft::core::Convergence& convergence) {
    CHECK(convergence.criterion != nullptr);
    return std::max(convergence.criterion->value, convergence.floor.value_or(0.0));
}

inline bool close(double left, double right) {
    return std::abs(left - right) <= 1.0e-6 * std::max(1.0, std::abs(right));
}

inline bool close_with_tolerance(double left, double right, double tolerance) {
    return std::abs(left - right) <= tolerance * std::max(1.0, std::abs(right));
}

inline double centered_step(double value) {
    return std::sqrt(std::numeric_limits<double>::epsilon())
        * std::max(std::abs(value), 1.0);
}

// Old normalized correlation ln K = a + 1200/T + 0.4 ln T + 0.002 T; P3 carries the
// implicit 1 K reference through its legacy normalization at 300 K on [250, 450] K.
// The retained legacy offset a = 0.3 gives K ~ 870-1010 (dilute A, x_A ~ 1e-3); the
// owner-selected qualification offset makes K(350 K) = 2 with unchanged b, c, d.
constexpr double DILUTE_ISOMERISATION_A = 0.3;
inline double moderate_isomerisation_a() {
    return std::log(2.0) - (1200.0 / 350.0 + 0.4 * std::log(350.0) + 0.002 * 350.0);
}
inline epcsaft::ReactionLogPolynomial legacy_isomerisation(double sign, double a) {
    return {sign * a + sign * 0.4 * std::log(300.0), sign * 1200.0, sign * 0.4, sign * 0.002, 300.0, 250.0, 450.0,
        epcsaft::equilibrium::EOS_STANDARD_STATE_ID};
}

inline void check_finite_difference(const epcsaft::equilibrium::CompiledSystem& compiled,
                             const std::vector<double>& point) {
    const auto analytic = compiled.jacobian(point);
    const auto base = compiled.evaluate(point).values;
    const auto context = equilibrium_support::weighted_hessian_context(compiled);
    CHECK(analytic.size() == base.size() * point.size());
    for (std::size_t column = 0; column < point.size(); ++column) {
        auto plus = point;
        auto minus = point;
        // Centred differences balance truncation and rounding at eps^(1/3), not eps^(1/2):
        // the rows cancel O(10) potentials, so the sqrt(eps) step carries ~5e-7 rounding noise.
        const double delta = std::cbrt(std::numeric_limits<double>::epsilon()) * std::max(std::abs(point[column]), 1.0);
        plus[column] += delta;
        minus[column] -= delta;
        const auto up = compiled.evaluate(plus).values;
        const auto down = compiled.evaluate(minus).values;
        for (std::size_t row = 0; row < base.size(); ++row) {
            const double centered = (up[row] - down[row]) / (2.0 * delta);
            const double model_value = analytic[row * point.size() + column];
            const bool verdict = close(model_value, centered);
            const auto point_id = std::to_string(row) + "-" + std::to_string(column);
            recording_support::compare(
                (std::string("jacobian/") + point_id).c_str(),
                "equilibrium_jacobian", "1", model_value, centered,
                recording_support::fields({
                    context,
                    recording_support::field("point", point),
                    recording_support::field("row", row),
                    recording_support::field("column", column),
                    recording_support::field("delta", delta),
                }),
                "Relationship", "finite-difference", "equilibrium Jacobian matches centered residual derivative",
                "reference_or_one", 0.0, 1.0e-6, verdict,
                "f5-jacobian", point_id.c_str(), "J", "T,P,\\mathbf{x},\\Delta z",
                "J_{ij}=\\Delta r_i/\\Delta z_j", "row,column");
        }
    }
}

inline std::vector<std::size_t> callback_rows(
    const epcsaft::equilibrium::CompiledSystem& compiled
) {
    std::vector<std::size_t> result;
    for (std::size_t row = 0; row < compiled.rows.size(); ++row) {
        if (compiled.feasibility_cell()
            || (compiled.rows[row].kind != epcsaft::equilibrium::RowKind::Chemical
                && compiled.rows[row].kind != epcsaft::equilibrium::RowKind::Pressure)) {
            result.push_back(row);
        }
    }
    return result;
}

inline std::vector<double> weighted_gradient(
    const epcsaft::equilibrium::CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& multipliers,
    double objective_factor
) {
    auto result = compiled.objective_gradient(point);
    const auto jacobian = compiled.jacobian(point);
    const auto rows = callback_rows(compiled);
    if (result.size() != point.size()
        || jacobian.size() != compiled.rows.size() * point.size()
        || multipliers.size() != rows.size()) {
        throw std::runtime_error("weighted-gradient shape mismatch");
    }
    for (double& value : result) value *= objective_factor;
    for (std::size_t multiplier = 0; multiplier < rows.size(); ++multiplier) {
        for (std::size_t column = 0; column < point.size(); ++column) {
            result[column] += multipliers[multiplier]
                * jacobian[rows[multiplier] * point.size() + column];
        }
    }
    return result;
}

inline double maximum_abs(const std::vector<double>& values) {
    double result = 0.0;
    for (const double value : values) result = std::max(result, std::abs(value));
    return result;
}

inline void check_weighted_hessian_gradient_ladder(
    const epcsaft::equilibrium::CompiledSystem& compiled,
    const std::vector<double>& point,
    const std::vector<double>& multipliers,
    const char* label,
    double tolerance
) {
    const auto rows = callback_rows(compiled);
    epcsaft::equilibrium::GibbsMinimizer nlp(compiled, point);
    Ipopt::Index variables = 0, constraints = 0, jacobian_entries = 0, hessian_entries = 0;
    Ipopt::TNLP::IndexStyleEnum index_style;
    CHECK(nlp.get_nlp_info(variables, constraints, jacobian_entries, hessian_entries, index_style));
    CHECK(variables >= 0 && constraints >= 0);
    CHECK(point.size() == static_cast<std::size_t>(variables));
    CHECK(multipliers.size() == rows.size());
    CHECK(multipliers.size() == static_cast<std::size_t>(constraints));

    const auto objective_factors = {0.0, 0.5, 1.0};
    const auto step_factors = {10.0, 3.0, 1.0};
    const auto context = equilibrium_support::weighted_hessian_context(compiled);

    std::vector<Ipopt::Index> hessian_rows(static_cast<std::size_t>(hessian_entries));
    std::vector<Ipopt::Index> hessian_columns(static_cast<std::size_t>(hessian_entries));
    CHECK(nlp.eval_h(variables, point.data(), true, 1.0, constraints,
                     multipliers.data(), true, hessian_entries,
                     hessian_rows.data(), hessian_columns.data(), nullptr));
    for (const double objective_factor : objective_factors) {
        std::vector<Ipopt::Number> hessian(
            static_cast<std::size_t>(hessian_entries), 0.0
        );
        CHECK(nlp.eval_h(variables, point.data(), false, objective_factor,
                         constraints, multipliers.data(), false,
                         hessian_entries, hessian_rows.data(),
                         hessian_columns.data(), hessian.data()));
        for (const double step_factor : step_factors) {
            for (std::size_t column = 0; column < point.size(); ++column) {
                auto plus = point;
                auto minus = point;
                const double step = step_factor * centered_step(point[column]);
                plus[column] += step;
                minus[column] -= step;
                const auto grad_plus = weighted_gradient(
                    compiled, plus, multipliers, objective_factor
                );
                const auto grad_minus = weighted_gradient(
                    compiled, minus, multipliers, objective_factor
                );
                for (std::size_t row = column; row < point.size(); ++row) {
                    const double centered = (grad_plus[row] - grad_minus[row]) / (2.0 * step);
                    const std::size_t index = row * (row + 1) / 2 + column;
                    const auto comparison_id = std::string("hessian/") + label
                        + "/objective-" + recording_support::number(objective_factor)
                        + "/step-" + recording_support::number(step_factor)
                        + "/" + std::to_string(row) + "-" + std::to_string(column);
                    const auto point_id = std::string("case-") + label
                        + "/objective-" + recording_support::number(objective_factor)
                        + "/step-" + recording_support::number(step_factor)
                        + "/row-" + std::to_string(row)
                        + "/column-" + std::to_string(column);
                    const auto coordinates = recording_support::fields({
                        recording_support::field("case", label),
                        recording_support::field("objective_factor", objective_factor),
                        recording_support::field("step_factor", step_factor),
                        recording_support::field("row", row),
                        recording_support::field("column", column),
                    });
                    const auto inputs = recording_support::fields({
                        context,
                        recording_support::field("case", label),
                        recording_support::field("objective_factor", objective_factor),
                        recording_support::field("step_factor", step_factor),
                        recording_support::field("point", point),
                        recording_support::field("multipliers", multipliers),
                        recording_support::field("row", row),
                        recording_support::field("column", column),
                        recording_support::field("step", step),
                        recording_support::field("variables", point.size()),
                        recording_support::field("constraints", multipliers.size()),
                    });
                    const bool pressure_entry = compiled.problem.pressure.free
                        && (row == compiled.primal_size - 1 || column == compiled.primal_size - 1);
                    const double coordinate_scale = pressure_entry ? compiled.problem.pressure.value : 1.0;
                    const double absolute = pressure_entry ? 1.0e-12 / coordinate_scale : 0.0;
                    const double allowance = pressure_entry ? absolute + tolerance * std::abs(centered)
                        : tolerance * std::max(1.0, std::abs(centered));
                    const bool verdict = std::abs(hessian[index] - centered) <= allowance;
                    recording_support::compare(
                        comparison_id.c_str(), "weighted_lagrangian_hessian", "coordinate",
                        hessian[index], centered, inputs, "Relationship",
                        "weighted_hessian_gradient_finite_difference",
                        "nabla2_L=d(weighted_gradient)/d(coordinate)", pressure_entry ? "reference_scale" : "reference_or_one",
                        absolute, tolerance, verdict, "f5-hessian-ladder", point_id.c_str(),
                        "\\nabla^2\\mathcal L", "\\mathbf{z},\\lambda,\\Delta z",
                        "\\nabla^2\\mathcal L_{ij}=\\Delta(\\nabla\\mathcal L)_i/\\Delta z_j",
                        coordinates.c_str());
                }
            }
        }
    }
}
// E2 and E5 share the identical-ethane pair and its ideal-gas records.
namespace e2 {

inline epcsaft::Parameters identical_ethane_pair() {
    auto result = properties_support::parameters(2);
    const auto ethane = properties_support::ethane();
    for (std::size_t i = 0; i < 2; ++i) {
        result.segment_count[i] = ethane.segment_count[0];
        result.sigma_constant[i] = ethane.sigma_constant[0];
        result.dispersion_energy_over_k[i] = ethane.dispersion_energy_over_k[0];
        result.molar_mass_kg_per_mol[i] = ethane.molar_mass_kg_per_mol[0];
    }
    return result;
}

inline epcsaft::ThermochemistryRecord species_record(double cp_a, double cp_b) {
    epcsaft::ThermochemistryRecord record;
    record.standard_pressure = 1.0e5;
    const double anchor_entropy = std::log(350.0 / 300.0);
    record.component_intervals = {
        {epcsaft::IdealInterval{250.0, 450.0, true, true, epcsaft::IdealCorrelation{
            epcsaft::IdealPolynomial{{cp_a}, 350.0, 1500.0, 30.0 * anchor_entropy}, 1.0e5}}},
        {epcsaft::IdealInterval{250.0, 450.0, true, true, epcsaft::IdealCorrelation{
            epcsaft::IdealPolynomial{{cp_b}, 350.0, 1000.0, 20.0 * anchor_entropy}, 1.0e5}}}};
    return record;
}

}  // namespace e2

}  // namespace equilibrium_checks

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
