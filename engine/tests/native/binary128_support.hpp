#pragma once
#include "quad_base.hpp"
#include "equilibrium_support.hpp"
#include "properties/helmholtz_tape.hpp"
#include "properties/thermochemistry.hpp"
#include "core/physical_constants.hpp"

#include <array>
#include <string>
#include <vector>

// Binary128 value oracles for verification: a named observable, a named Richardson
// direction, and the water/Na/Cl oracle built from the EOS kernels alone.
namespace binary128_support {
using namespace epcsaft;
using namespace epcsaft::equilibrium;


struct Quantity { std::string name, unit; double atol; SolvedStateObservable observable; };
struct NamedDirection {
    std::string name;
    double step;
    SolvedStateActionDirection value;
    std::vector<SolvedStateActionDirection> chart; // Feasible chart sub-directions, including both tied EOS slots.
    bool parameter = false;
};

#ifndef _WIN32
using Wide = CppAD::Quad;

// Manufactured electroneutral water/Na/Cl value oracle. The two symmetric
// water sites have the exact mass-action root X=2/(1+sqrt(1+4 rho_N xw Delta)).
// Chemical equilibrium is independently re-solved on x=(1-2y,y,y); no
// production action or standard-state derivative supplies a perturbed value.
// `values` takes the double solve at the same inputs only to seed and check it.
class WideIonic {
    static constexpr double R = GAS_CONSTANT_J_PER_MOL_K;
    struct Inputs {const Parameters& parameters;const ThermochemistryRecord& records;const EquilibriumProblem& problem;} c;
    CppAD::ADFun<Wide> fluid,reference,reference_coarse;
    Wide max_pressure_residual=0,max_pressure_floor=0;
    struct State {
        Wide v,a,at,av,att,vt,vp,C,ct,ctt,cp,ctp,cpp,cy;
        std::array<Wide,3> mu;
        Wide pressure;
    };
    CppAD::ADFun<Wide> record(double y) {
        using AD=CppAD::AD<Wide>;
        const std::vector<double> rx{1-2*y,y,y};
        std::vector<AD> independent{Wide(c.problem.temperature.value),Wide(1.8e-5),Wide(rx[0]),Wide(y),Wide(y),Wide(c.parameters.dispersion_energy_over_k[0])};
        CppAD::Independent(independent);
        ParameterRecord<AD> p{{}, AD(Wide(0)), AD(Wide(0))};static_cast<FixedParameters&>(p)=c.parameters;
        for(std::size_t family=0;family<detail::PARAMETER_FAMILY_COUNT;++family) {
            const auto f=static_cast<detail::ParameterFamily>(family);const auto& src=detail::family_values(c.parameters,f);
            const auto dst=detail::family_values(p,f,src.second);
            for(std::size_t i=0;i<src.second;++i) dst.first[i]=AD(Wide(src.first[i]));
        }
        p.dispersion_energy_over_k[0]=independent[5];
        const AD N=independent[2]+independent[3]+independent[4],rho=N/independent[1];
        const std::vector<AD> x{independent[2]/N,independent[3]/N,independent[4]/N};
        const auto sigma=segment_sigma(independent[0],p),diameter=segment_diameter(independent[0],p,sigma);
        const auto delta=association_delta_values(independent[0],rho,x,p,sigma,diameter,p.association_energy_over_k,p.association_volume);
        const AD site=2/(1+CppAD::sqrt(1+4*rho*NUMBER_DENSITY_PER_ANGSTROM3*x[0]*delta[1]));
        const std::vector<AD> sites{site,site};
        const auto electro=electrolyte_contribution(independent[0],rho,x,rx,p);
        const auto born=born_contribution(independent[0],x,p,electro.bulk_relative_permittivity);
        const auto terms=residual_helmholtz(independent[0],rho,x,p,sigma,diameter,sites,electro.debye_huckel,born);
        std::vector<AD> dependent{N*(terms.hard_chain+terms.dispersion+terms.association+terms.debye_huckel+terms.born+terms.dipolar)};
        CppAD::ADFun<Wide> energy(independent,dependent);auto nested=energy.base2ad();
        std::vector<AD> point{Wide(c.problem.temperature.value),Wide(1.8e-5),Wide(rx[0]),Wide(y),Wide(y),Wide(c.parameters.dispersion_energy_over_k[0])};
        CppAD::Independent(point);
        const auto a=nested.Forward(0,point),j=nested.Jacobian(point);
        const AD total=point[2]+point[3]+point[4];
        dependent={a[0],AD(Wide(R))*point[0]*(total/point[1]-j[1]),j[3]+j[4]-j[2]+CppAD::log(total/point[1]),j[2],j[3],j[4]};
        return CppAD::ADFun<Wide>(point,dependent);
    }
    State state(CppAD::ADFun<Wide>& tape,Wide T,Wide P,Wide epsilon,Wide y,Wide anchor,bool jets=true) {
        std::vector<Wide> point{T,1/anchor,1-2*y,y,y,epsilon};
        // At liquid density, ideal and residual pressures cancel. The measured
        // binary128 floor is governed by RT/v, not the small net pressure P.
        auto limit=[&] {const Wide floor=128*CppAD::numeric_limits<Wide>::epsilon()*Wide(R)*T/point[1];
            const Wide requested=Wide(ldexpq(1,-100))*P;return floor>requested?floor:requested;};
        for(std::size_t iteration=0;iteration<12;++iteration) {
            const auto f=tape.Forward(0,point),j=tape.Jacobian(point);const Wide residual=f[1]-P;
            if(CppAD::abs(residual)<=limit()) break;
            CHECK(j[7]<0);point[1]-=residual/j[7];
        }
        const auto f=tape.Forward(0,point),j=tape.Jacobian(point);
        // Volume forward error |dv/v|<=2^-100 (#272): it moves the values ~1e-30 relative, far below the
        // 1e-6 Richardson targets; the pressure limit above is the evaluation noise floor, not an accuracy.
        const Wide closure=Wide(ldexpq(1,-100))*CppAD::abs(j[7])*point[1];
        CHECK(CppAD::abs(f[1]-P)<=closure && point[1]>0 && j[7]<0);
        if(CppAD::abs(f[1]-P)>max_pressure_residual) max_pressure_residual=CppAD::abs(f[1]-P);
        if(closure>max_pressure_floor) max_pressure_floor=closure;
        const auto ah=jets?tape.Hessian(point,0):std::vector<Wide>(36),
            ph=jets?tape.Hessian(point,1):std::vector<Wide>(36),ch=jets?tape.Hessian(point,2):std::vector<Wide>(36);
        const Wide vt=-j[6]/j[7],vp=1/j[7];
        const Wide vtt=-(ph[0]+2*ph[1]*vt+ph[7]*vt*vt)/j[7];
        const Wide vtp=-(ph[1]+ph[7]*vt)*vp/j[7],vpp=-ph[7]*vp*vp/j[7];
        const Wide py=-2*j[8]+j[9]+j[10],vy=-py/j[7];
        return {point[1],f[0],j[0],j[1],ah[0],vt,vp,f[2],j[12]+j[13]*vt,
            ch[0]+2*ch[1]*vt+ch[7]*vt*vt+j[13]*vtt,j[13]*vp,
            (ch[1]+ch[7]*vt)*vp+j[13]*vtp,ch[7]*vp*vp+j[13]*vpp,
            -2*j[14]+j[15]+j[16]+j[13]*vy,{f[3],f[4],f[5]},f[1]};
    }
public:
    // EqIDs: mean_ionic_activity_coefficient, molality_activity_conversion,
    // practical_osmotic_coefficient. Verification values for the same symmetric
    // water/Na/Cl fixture, at fixed salt molality, without a reaction solve.
    std::array<Wide,2> activity_values(double molality, double anchor) {
        const Wide T=c.problem.temperature.value,P=c.problem.pressure.value;
        const Wide mass=c.parameters.molar_mass_kg_per_mol[0],s=2*mass*Wide(molality);
        const Wide epsilon=c.parameters.dispersion_energy_over_k[0];
        const auto finite=state(fluid,T,P,epsilon,s/(2*(1+s)),Wide(anchor),false);
        auto pure_tape=record(0);
        const auto pure=state(pure_tape,T,P,epsilon,Wide(0),Wide(anchor),false);
        const auto log_z=[&](const State& state){return CppAD::log(P*state.v/(Wide(R)*T));};
        // EqID reference_limit_terminal_pair: the ion reference is the limit of the 1e-12 and 1e-14 mol/kg rungs, as in values().
        const auto ion_sum=[&](CppAD::ADFun<Wide>& tape,Wide m){const Wide y=mass*m/(1+2*mass*m);
            const auto trace=state(tape,T,P,epsilon,y,Wide(anchor),false);return trace.mu[1]+trace.mu[2]-2*log_z(trace);};
        const Wide limit=(10*ion_sum(reference,Wide(1e-14))-ion_sum(reference_coarse,Wide(1e-12)))/9;
        const Wide mean=(finite.mu[1]+finite.mu[2]-limit)/2-log_z(finite)-CppAD::log(1+s);
        const Wide solvent=-CppAD::log(1+s)+finite.mu[0]-pure.mu[0]-log_z(finite)+log_z(pure);
        return {CppAD::exp(mean),-solvent/s};
    }
    // Solvent molar enthalpy at the declared (T, P): the reference tape at its 1e-14 mol/kg trace state
    // (ion fraction 1.8e-16) with the record's ideal-gas enthalpy.
    Wide solvent_enthalpy(double anchor) {
        const Wide T=c.problem.temperature.value,mass=c.parameters.molar_mass_kg_per_mol[0],m=Wide(1e-14);
        const auto s=state(reference,T,Wide(c.problem.pressure.value),Wide(c.parameters.dispersion_energy_over_k[0]),mass*m/(1+2*mass*m),Wide(anchor));
        const auto thermal=evaluate_ideal_correlation(c.records.component_intervals[0][0].correlation,CppAD::AD<Wide>(T),1e5);
        return CppAD::Value(thermal.h)-Wide(R)*T*T*s.at-Wide(R)*T*s.v*s.av;
    }
    double pressure_residual() const {return double(max_pressure_residual);}
    double pressure_limit() const {return double(max_pressure_floor);}
    // `salt_fraction` is the ion mole fraction y of the double solve that seeds the fluid tape.
    WideIonic(const Parameters& parameters,const ThermochemistryRecord& records,const EquilibriumProblem& problem,double salt_fraction)
        :c{parameters,records,problem} {
        CHECK(c.parameters.site_component==std::vector<int>({0,0}));
        CHECK(c.parameters.site_multiplicity==std::vector<int>({1,1}));
        CHECK(c.parameters.association_active==std::vector<int>({0,1,1,0}));
        CHECK(c.parameters.association_energy_over_k[1]==c.parameters.association_energy_over_k[2]);
        CHECK(c.parameters.association_volume[1]==c.parameters.association_volume[2]);
        fluid=record(salt_fraction);reference=record(c.parameters.molar_mass_kg_per_mol[0]*1e-14);
        reference_coarse=record(c.parameters.molar_mass_kg_per_mol[0]*1e-12);
    }
    std::vector<Wide> values(const std::vector<Quantity>& q,const std::vector<NamedDirection>& ds,double factor,std::size_t mask,
                             const EquilibriumResult& solved,double rho) {
        max_pressure_residual=0;max_pressure_floor=0;
        Wide T=c.problem.temperature.value,P=c.problem.pressure.value,epsilon=c.parameters.dispersion_energy_over_k[0];
        std::array<Wide,3> feed{Wide(c.problem.feed.values[0]),Wide(c.problem.feed.values[1]),Wide(c.problem.feed.values[2])};
        const auto& law=*c.problem.reactions[0].correlation;
        std::array<Wide,4> coeff{Wide(law.a),Wide(law.b),Wide(law.c),Wide(law.d)};
        for(std::size_t k=0;k<ds.size();++k) {
            const Wide scale=Wide((mask>>k)&1?1:-1)*Wide(ds[k].step)*Wide(factor);const auto& d=ds[k].value;
            T+=scale*d.temperature_k;P+=scale*d.pressure_pa;epsilon+=scale*d.active_parameters[0];
            for(std::size_t i=0;i<3;++i) feed[i]+=scale*d.feed_amounts_mol[i];
            for(std::size_t i=0;i<4;++i) coeff[i]+=scale*d.reaction_coefficients[0][i];
        }
        CHECK(feed[1]==feed[2]);
        const Wide mass=c.parameters.molar_mass_kg_per_mol[0],m=Wide(1e-14);
        const Wide ys=mass*m/(1+2*mass*m);const auto& source=*c.problem.reactions[0].reference;const bool fixed=source.reference_pressure_pa.has_value();
        const Wide sourceP=fixed?Wide(*source.reference_pressure_pa):P;
        // The production reference is the sqrt(m) Richardson limit of its 1e-12 and 1e-14 mol/kg points (#353).
        auto src=state(reference,T,sourceP,epsilon,ys,rho);
        const Wide yc=mass*Wide(1e-12)/(1+2*mass*Wide(1e-12));const auto coarse=state(reference_coarse,T,sourceP,epsilon,yc,rho);
        for(auto field:{&State::C,&State::ct,&State::ctt,&State::cp,&State::ctp,&State::cpp}) src.*field=(10*(src.*field)-coarse.*field)/9;
        const Wide lnK=coeff[0]+coeff[1]/T+coeff[2]*CppAD::log(T/Wide(law.reference_temperature))+coeff[3]*T;
        const Wide lnKT=-coeff[1]/(T*T)+coeff[2]/T+coeff[3],lnKTT=2*coeff[1]/(T*T*T)-coeff[2]/(T*T);
        // The raw mole-fraction basis has no common-molality scale 2 ln M_w.
        const Wide scale=source.source_basis==ReferenceSourceBasis::RawMoleFractionInfiniteDilution?Wide(0):2*CppAD::log(mass);
        const Wide target=lnK+scale+src.C;
        Wide y=solved.mole_fractions[1];State sys{};
        Wide chemical=0;
        for(std::size_t iteration=0;iteration<8;++iteration) {
            sys=state(fluid,T,P,epsilon,y,rho,false);
            chemical=2*CppAD::log(y)-CppAD::log(1-2*y)+sys.C-target;
            if(CppAD::abs(chemical)<ldexpq(1,-96)) break;
            y-=chemical/(2/y+2/(1-2*y)+sys.cy);CHECK(y>0&&y<Wide(0.5));
        }
        CHECK(CppAD::abs(chemical)<ldexpq(1,-96));
        sys=state(fluid,T,P,epsilon,y,rho);
        const std::array<Wide,3> x{1-2*y,y,y};const Wide n=(feed[0]+(feed[1]+feed[2])/2)/(1-y),rt=Wide(R)*T;
        CHECK(CppAD::abs(n-Wide(solved.amounts[0]))<Wide(1e-10));
        for(std::size_t i=0;i<3;++i) CHECK(CppAD::abs(x[i]-Wide(solved.mole_fractions[i]))<Wide(1e-10));
        CHECK(CppAD::abs(1/sys.v-Wide(rho))<Wide(1e-6));
        const auto thermal=evaluate_ideal_correlation(c.records.component_intervals[0][0].correlation,CppAD::AD<Wide>(T),1e5);
        const Wide h=CppAD::Value(thermal.h),s=CppAD::Value(thermal.s),cp=CppAD::Value(thermal.cp);
        std::array<Wide,3> g{},gt{},gtt{},gp{},gtp{},gpp{};
        g[0]=h-T*s+rt*CppAD::log(rt/Wide(1e5));gt[0]=-s+Wide(R)*(CppAD::log(rt/Wide(1e5))+1);gtt[0]=(Wide(R)-cp)/T;
        for(std::size_t i=1;i<3;++i) {
            g[i]=(g[0]-rt*target)/2;gt[i]=(gt[0]-Wide(R)*(target+T*(lnKT+src.ct)))/2;
            gtt[i]=(gtt[0]-Wide(R)*(2*(lnKT+src.ct)+T*(lnKTT+src.ctt)))/2;
            if(!fixed) {gp[i]=-rt*src.cp/2;gtp[i]=-Wide(R)*(src.cp+T*src.ctp)/2;gpp[i]=-rt*src.cpp/2;}
        }
        std::array<Wide,5> ideal{},total{};Wide volume=sys.v,vt=sys.vt,vp=sys.vp;
        for(std::size_t i=0;i<3;++i) {
            ideal[0]+=x[i]*(g[i]-T*gt[i]);ideal[2]-=x[i]*(gt[i]+Wide(R)*CppAD::log(x[i]/sys.v));ideal[3]-=x[i]*T*gtt[i];
            volume+=x[i]*gp[i];vt+=x[i]*gtp[i];vp+=x[i]*gpp[i];
        }
        ideal[0]+=rt;ideal[2]+=R;ideal[3]+=R;ideal[1]=ideal[0]-rt-P*(volume-sys.v);
        ideal[4]=ideal[3]-R+T*vt*vt/vp-T*sys.vt*sys.vt/sys.vp;
        const Wide cvres=-Wide(R)*(2*T*sys.at+T*T*sys.att);
        const std::array<Wide,5> residual{-Wide(R)*T*T*sys.at-rt*sys.v*sys.av,-Wide(R)*T*T*sys.at,-Wide(R)*(sys.a+T*sys.at),cvres-T*sys.vt*sys.vt/sys.vp-R,cvres};
        for(std::size_t i=0;i<5;++i) total[i]=ideal[i]+residual[i];
        const std::vector<Wide> wx{x.begin(),x.end()},eps{epsilon,Wide(c.parameters.dispersion_energy_over_k[1]),Wide(c.parameters.dispersion_energy_over_k[2])};
        ParameterRecord<Wide> geometry{{},Wide(0),Wide(0)};
        static_cast<FixedParameters&>(geometry)=c.parameters;
        geometry.dispersion_energy_over_k=eps;
        geometry.packing_diameter_angstrom.assign(c.parameters.packing_diameter_angstrom.begin(),c.parameters.packing_diameter_angstrom.end());
        const auto sigma=segment_sigma(T,c.parameters),diameter=segment_diameter(T,geometry,sigma);
        const Wide eta=packing_fraction_from_geometry(1/sys.v,wx,c.parameters,diameter);
        std::vector<Wide> out;
        for(const auto& item:q) {
            const auto& o=item.observable;Wide v=0;
            switch(o.kind) {
            case SolvedStateObservableKind::Temperature:v=T;break;
            case SolvedStateObservableKind::Pressure:v=P;break;
            case SolvedStateObservableKind::PhaseAmount:v=n;break;
            case SolvedStateObservableKind::PhaseVolume:v=n*sys.v;break;
            case SolvedStateObservableKind::PhaseMolarDensity:v=1/sys.v;break;
            case SolvedStateObservableKind::PhasePackingFraction:v=eta;break;
            case SolvedStateObservableKind::PhasePressure:v=sys.pressure;break;
            case SolvedStateObservableKind::PhaseMoleFraction:v=x[o.component];break;
            case SolvedStateObservableKind::PhaseComponentAmount:v=n*x[o.component];break;
            case SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt:v=sys.mu[o.component];break;
            case SolvedStateObservableKind::PhaseLogActivity:v=CppAD::log(x[o.component]/sys.v)+sys.mu[o.component];break;
            case SolvedStateObservableKind::PhaseGibbsMolarVolume:v=volume;break;
            case SolvedStateObservableKind::PhaseProperty:{const auto k=std::size_t(o.property);v=k<5?ideal[k]:total[k-10];break;}
            }
            if(item.name=="H") v=n*total[0];if(item.name=="U_G") v=n*total[1];
            if(item.name=="S") v=n*total[2];if(item.name=="V_G") v=n*volume;out.push_back(v);
        }
        return out;
    }
};
#endif
} // namespace binary128_support
