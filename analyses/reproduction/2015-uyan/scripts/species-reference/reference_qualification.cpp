#include "properties_support.hpp"
#include "equilibrium_support.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace epcsaft;
constexpr double R = GAS_CONSTANT_J_PER_MOL_K;
std::ofstream out;
int failures = 0;
void compare(std::string name, double t, double p, double step, double actual, double expected,
             double scale, double tolerance) {
    const bool accepted = step != .1 && step != 100. && (name.find("gas-limit-") != 0 || p == .01);
    const double error = std::abs(actual-expected)*scale/std::max(1.0,std::abs(expected)*scale);
    out << name << '\t' << t << '\t' << p << '\t' << step << '\t' << actual << '\t'
        << expected << '\t' << scale << '\t' << error << '\t' << tolerance << '\t' << accepted << std::endl;
    if (accepted && !(error <= tolerance)) { ++failures; std::cerr << name << " error " << error << '\n'; }
}
Parameters ethane_pair() {
    auto p = properties_support::parameters(2);
    p.segment_count={1.6069,1.6069};p.sigma_constant={3.5206,3.5206};
    p.dispersion_energy_over_k={191.42,191.42};p.molar_mass_kg_per_mol={.03007,.03007};
    return p;
}
NeutralReferenceDeclaration declaration() {
    NeutralReferenceDeclaration d{0,{0,1},1};
    d.species_reference={SpeciesReference::PureLiquid,SpeciesReference::PureGas};
    return d;
}
NeutralReferenceReactionDiagnostic reference(const Parameters& parameters,double t,double p,bool fixed=false,
    const std::vector<double>& nu={-2,1}, std::size_t order=4) {
    HelmholtzModel model(parameters,{{ParameterFamily::DispersionEnergyOverK,{1}}});
    ReactionReference ref{fixed ? std::optional<double>(2e5) : std::nullopt,ReferenceSourceBasis::RawMoleFractionInfiniteDilution};
    auto d = reaction_references(model,t,p,declaration(),{nu},{ref},true,std::nullopt,2,true,order)[0];
    if (!d.terminal_contraction) throw std::runtime_error(d.message);
    return d;
}
double partial(const NeutralReferenceReactionDiagnostic& d,int a,int b,bool theta=false) {
    for (const auto* columns : {&d.lower_derivatives,&d.derivatives}) for (const auto& c : *columns)
        if (c.temperature_order==std::size_t(a) && c.pressure_order==std::size_t(b)
            && c.parameter.has_value()==theta && c.value) return *c.value;
    if (!theta && a==0 && b==1 && d.terminal_contraction_pressure_derivative) return *d.terminal_contraction_pressure_derivative;
    if (!theta && a==0 && b==2 && d.terminal_contraction_pressure_second_derivative) return *d.terminal_contraction_pressure_second_derivative;
    if (!theta && a==1 && b==1 && d.terminal_contraction_temperature_pressure_derivative) return *d.terminal_contraction_temperature_pressure_derivative;
    throw std::runtime_error(d.derivative_cause(a,b,theta ? std::optional<std::size_t>(0) : std::nullopt));
}
double lambda(const State& s,double t,double p,int i) {
    return require_component_property((*s.log_fugacity_coefficient)[i])+std::log(p/(R*t));
}
int main(int argc,char** argv) {
    if (argc!=3) return 2;
    out.open(argv[1]);out << std::setprecision(17);
    out << "quantity\tT_K\tP_Pa\tstep\tactual\texpected\tscale\tscaled_error\ttolerance\taccepted\n";
    auto co2=ethane_pair();co2.segment_count[1]=2.0729;co2.sigma_constant[1]=2.7852;
    co2.dispersion_energy_over_k[1]=169.21;co2.molar_mass_kg_per_mol[1]=.0440095;
    for(double t:{313.15,353.15,393.15}) {
        HelmholtzModel m(co2);
        auto d=reference(co2,t,1e5,false,{0,1},0);
        const auto vapor=m.state_at_pressure(t,1e5,{0,1},DensityBranch::Vapor);
        compare("co2-vapor-identity",t,1e5,0,*d.terminal_contraction,lambda(vapor,t,1e5,1),1,1e-10);
        compare("co2-vapor-T",t,1e5,0,partial(d,1,0),-(vapor.residual_enthalpy+R*t)/(R*t*t),t,1e-9);
        compare("co2-vapor-P",t,1e5,0,partial(d,0,1),1/(vapor.molar_density*R*t),1e5,1e-9);
        for(double p:{100.,10.,1.,.1,.01}) {
            d=reference(co2,t,p,false,{0,1},0);
            compare("gas-limit-L",t,p,0,*d.terminal_contraction-std::log(p/(R*t)),0,1,1e-7);
            compare("gas-limit-T",t,p,0,t*partial(d,1,0)+1,0,1,1e-7);
            compare("gas-limit-P",t,p,0,p*partial(d,0,1)-1,0,1,1e-7);
        }
    }
    auto p=ethane_pair();const double t=260,pressure=1e5;
    for(bool fixed:{false,true}) {
        const auto center=reference(p,t,pressure,fixed);
        for(double step:{.1,.05}) {
            const auto lo=reference(p,t-step,pressure,fixed);
            const auto hi=reference(p,t+step,pressure,fixed);
            compare(fixed?"fixed-T":"system-T",t,pressure,step,partial(center,1,0),(*hi.terminal_contraction-*lo.terminal_contraction)/(2*step),t,1e-6);
            for(int order:{2,3,4})
                compare((fixed?"fixed-T":"system-T")+std::to_string(order),t,pressure,step,
                    partial(center,order,0),(partial(hi,order-1,0)-partial(lo,order-1,0))/(2*step),std::pow(t,order),order>2?1e-5:1e-6);
            compare(fixed?"fixed-TP":"system-TP",t,pressure,step,partial(center,1,1),
                (partial(hi,0,1)-partial(lo,0,1))/(2*step),t*pressure,1e-6);
            auto minus=p,plus=p;minus.dispersion_energy_over_k[1]-=step;plus.dispersion_energy_over_k[1]+=step;
            const auto l=reference(minus,t,pressure,fixed),h=reference(plus,t,pressure,fixed);
            compare(fixed?"fixed-theta":"system-theta",t,pressure,step,partial(center,0,0,true),
                (*h.terminal_contraction-*l.terminal_contraction)/(2*step),191.42,1e-6);
            compare(fixed?"fixed-Ttheta":"system-Ttheta",t,pressure,step,partial(center,1,0,true),
                (partial(h,1,0)-partial(l,1,0))/(2*step),t*191.42,1e-6);
        }
        for(double step:{100.,50.}) {
            const auto lo=reference(p,t,pressure-step,fixed),hi=reference(p,t,pressure+step,fixed);
            compare(fixed?"fixed-P":"system-P",t,pressure,step,partial(center,0,1),
                (*hi.terminal_contraction-*lo.terminal_contraction)/(2*step),pressure,1e-6);
            compare(fixed?"fixed-PP":"system-PP",t,pressure,step,partial(center,0,2),
                (partial(hi,0,1)-partial(lo,0,1))/(2*step),pressure*pressure,1e-6);
        }
        if(fixed) for(const auto* columns:{&center.lower_derivatives,&center.derivatives}) for(const auto& c:*columns)
            if(c.pressure_order && c.value && *c.value!=0) throw std::runtime_error("nonzero fixed-source P column");
    }
    HelmholtzModel model(p);
    auto decl=declaration();ReactionReference source{std::nullopt,ReferenceSourceBasis::RawMoleFractionInfiniteDilution};
    // Direct native pure branch identity at the existing 180 K saturation state.
    const double ps=std::stod(argv[2]);
    const auto liquid=model.state_at_pressure(180,ps,{0,1},DensityBranch::Liquid);
    const auto vapor=model.state_at_pressure(180,ps,{0,1},DensityBranch::Vapor);
    auto d=reaction_references(model,180,ps,decl,{{0,1}},{source})[0];
    compare("ethane-vapor-branch-P",180,ps,0,*d.terminal_contraction_pressure_derivative,1/(vapor.molar_density*R*180),ps,1e-9);
    if (!(vapor.molar_density < liquid.molar_density/2)) throw std::runtime_error("branches are not distinct");
    out << "ethane-liquid-P" << '\t' << 180 << '\t' << ps << '\t' << 0 << '\t'
        << 1/(liquid.molar_density*R*180) << "\t0\t1\t0\t0\t0\n";
    decl.species_reference={SpeciesReference::PureGas,SpeciesReference::PureLiquid};
    d=reaction_references(model,180,ps,decl,{{-2,1}},{source},true,std::nullopt,2)[0];
    const auto pure_solvent=model.state_at_pressure(180,ps,{1,0},DensityBranch::Vapor);
    compare("all-pure-value",180,ps,0,*d.terminal_contraction,
        -2*lambda(pure_solvent,180,ps,0)+lambda(liquid,180,ps,1),1,1e-10);
    compare("all-pure-P",180,ps,0,*d.terminal_contraction_pressure_derivative,
        -2/(pure_solvent.molar_density*R*180)+1/(liquid.molar_density*R*180),ps,1e-9);
    decl=declaration();
    source.source_basis=ReferenceSourceBasis::CommonMolalityInfiniteDilution;
    d=reaction_references(model,t,pressure,decl,{{-2,1}},{source})[0];
    if(d.status!=NeutralReferenceStatus::UnsupportedSourceBasis || d.terminal_contraction || d.message.find("molality")==std::string::npos)
        throw std::runtime_error("mixed molality not refused with cause");
    std::cout << "mixed molality: " << d.message << '\n';
    auto ions=properties_support::aqueous_electrolyte();HelmholtzModel ionic(ions);
    NeutralReferenceDeclaration charged{0,{0,1,1},1};charged.species_reference={SpeciesReference::PureLiquid,SpeciesReference::PureGas,SpeciesReference::InfiniteDilution};
    source.source_basis=ReferenceSourceBasis::RawMoleFractionInfiniteDilution;
    d=reaction_references(ionic,298.15,1e5,charged,{{-1,1,1}},{source})[0];
    if(d.status!=NeutralReferenceStatus::NonNeutralDeclaration || d.terminal_contraction || d.message.find("charged")==std::string::npos)
        throw std::runtime_error("charged pure reference not refused");
    std::cout << "charged pure: " << d.message << '\n';
    auto oversized=p;oversized.sigma_constant[1]=1e6;
    HelmholtzModel bounded(oversized);auto unavailable=declaration();unavailable.solute_molality_ratios[1]=1e-30;
    d=reaction_references(bounded,260,1e5,unavailable,{{-2,1}},{source})[0];
    if(d.terminal_contraction || d.status!=NeutralReferenceStatus::DensityUnavailable
        || d.message.find("pure reference species 1")==std::string::npos)
        throw std::runtime_error("unavailable pure root lost its cause: "+d.message);
    std::cout << "unavailable pure: " << d.message << '\n';
    // Existing record-anchored identical-ethane caloric case, with only the reference changed.
    // At 260 K and 3 or 5 MPa pure ethane has only its liquid root, so s1 declares the pure liquid.
    auto liquid_reference=declaration();liquid_reference.species_reference[1]=SpeciesReference::PureLiquid;
    using namespace epcsaft::equilibrium;
    ThermochemistryRecord record;record.standard_pressure=1e5;record.component_intervals.resize(2);
    record.component_intervals[0]={{200,400,true,true,{IdealPolynomial{{52,.12},300,0,0},1e5}}};
    HelmholtzModel caloric(p,{},record);
    for(bool fixed:{false,true}) {
        auto problem=equilibrium_support::flash(2);problem.temperature=FreeValue::fixed(260);
        problem.pressure=FreeValue::fixed(3e6);problem.feed.values={2,.2};problem.source_reference=liquid_reference;
        problem.reactions={{"2A=B",{-2,1},std::nullopt,std::nullopt,
            ReactionLogPolynomial{std::log(.3),200,3,-.004,300,250,350,"analytic-source"},
            ReactionReference{fixed?std::optional<double>(5e6):std::nullopt,ReferenceSourceBasis::RawMoleFractionInfiniteDilution}}};
        const auto compiled=compile_problem(caloric,problem);
        const auto column=[](const ReferenceCalorics& ref,int a) {
            for(const auto& c:ref.gibbs) if(c.temperature_order==std::size_t(a)&&c.pressure_order==0&&!c.parameter
                &&c.column.status==SolvedStateActionStatus::Available) return c.column.values[1];
            throw std::runtime_error(ref.message);
        };
        const auto center=reference_calorics(compiled,260,3e6,true,true,false,4);
        for(double step:{.1,.05}) {
            const auto lo=reference_calorics(compiled,260-step,3e6,true,true,false,4);
            const auto hi=reference_calorics(compiled,260+step,3e6,true,true,false,4);
            for(int order:{3,4}) compare((fixed?"caloric-fixed-g":"caloric-system-g")+std::to_string(order),260,3e6,step,
                column(center,order),(column(hi,order-1)-column(lo,order-1))/(2*step),std::pow(260,order)/R,1e-5);
        }
    }
    return failures?1:0;
}
