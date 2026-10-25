#pragma once
// Included inside reactive_caloric_checks. A manufactured neutral value oracle,
// not a second general equilibrium driver. Only the EOS kernels are shared;
// neither production reference jets nor solved-state actions supply its values.
class WideNeutral {
    const Case& c;
    CppAD::ADFun<Wide> thermodynamic;
    Wide max_pressure_residual=0,max_pressure_limit=0;
    struct State {
        Wide v, pressure_residual, a, at, av, att, mu, vt, vp, lambda, lt, ltt, lp, ltp, lpp;
    };
    State state(Wide T, Wide P, Wide epsilon, Wide anchor) {
        std::vector<Wide> point{T, 1/anchor, epsilon};
        Wide residual=0;
        for (std::size_t iteration=0;iteration<12;++iteration) {
            const auto values=thermodynamic.Forward(0,point);
            const auto jacobian=thermodynamic.Jacobian(point);
            residual=values[1]-P;
            if (CppAD::abs(residual)<=ldexpq(1,-100)*P) break;
            CHECK(jacobian[4]<0);
            point[1]-=residual/jacobian[4];
        }
        const auto value=thermodynamic.Forward(0,point),j=thermodynamic.Jacobian(point);
        residual=value[1]-P;CHECK(CppAD::abs(residual)<=ldexpq(1,-100)*P && point[1]>0 && j[4]<0);
        if(CppAD::abs(residual)>max_pressure_residual) max_pressure_residual=CppAD::abs(residual);
        const Wide limit=Wide(ldexpq(1,-100))*P;if(limit>max_pressure_limit) max_pressure_limit=limit;
        const auto ah=thermodynamic.Hessian(point,0),ph=thermodynamic.Hessian(point,1),lh=thermodynamic.Hessian(point,2);
        const Wide vt=-j[3]/j[4],vp=1/j[4];
        const Wide vtt=-(ph[0]+2*ph[1]*vt+ph[4]*vt*vt)/j[4];
        const Wide vtp=-(ph[1]+ph[4]*vt)*vp/j[4],vpp=-ph[4]*vp*vp/j[4];
        return {point[1],residual,value[0],j[0],j[1],ah[0],value[2]+CppAD::log(point[1]),vt,vp,value[2],
            j[6]+j[7]*vt,lh[0]+2*lh[1]*vt+lh[4]*vt*vt+j[7]*vtt,
            j[7]*vp,(lh[1]+lh[4]*vt)*vp+j[7]*vtp,lh[4]*vp*vp+j[7]*vpp};
    }
public:
    double pressure_residual() const {return double(max_pressure_residual);}
    double pressure_limit() const {return double(max_pressure_limit);}
    explicit WideNeutral(const Case& source):c(source) {
        CHECK(!c.ionic && c.parameters.site_component.empty());
        using AD=CppAD::AD<Wide>;
        std::vector<AD> independent{Wide(c.problem.temperature.value),Wide(1e-4),Wide(c.parameters.dispersion_energy_over_k[0])};
        CppAD::Independent(independent);
        ParameterRecord<AD> p{{}, AD(Wide(0)), AD(Wide(0))};
        static_cast<FixedParameters&>(p)=c.parameters;
        for(std::size_t family=0;family<detail::PARAMETER_FAMILY_COUNT;++family) {
            const auto f=static_cast<detail::ParameterFamily>(family);
            const auto& src=detail::family_values(c.parameters,f);
            const auto dst=detail::family_values(p,f,src.second);
            for(std::size_t i=0;i<src.second;++i) dst.first[i]=AD(Wide(src.first[i]));
        }
        for(auto& e:p.dispersion_energy_over_k) e=independent[2];
        const std::vector<AD> x{AD(1),AD(0)};
        const auto sigma=segment_sigma(independent[0],p),diameter=segment_diameter(independent[0],p,sigma);
        const auto terms=residual_helmholtz(independent[0],AD(1)/independent[1],x,p,sigma,diameter,{},AD(0),AD(0));
        std::vector<AD> dependent{terms.hard_chain+terms.dispersion};
        CppAD::ADFun<Wide> energy(independent,dependent);
        auto nested=energy.base2ad();
        std::vector<AD> point{Wide(c.problem.temperature.value),Wide(1e-4),Wide(c.parameters.dispersion_energy_over_k[0])};
        CppAD::Independent(point);
        const auto a=nested.Forward(0,point),j=nested.Jacobian(point);
        dependent={a[0],AD(Wide(R))*point[0]*(1/point[1]-j[1]),a[0]-point[1]*j[1]-CppAD::log(point[1])};
        thermodynamic=CppAD::ADFun<Wide>(point,dependent);
    }
    std::vector<Wide> values(const std::vector<Quantity>& q,const std::vector<NamedDirection>& ds,
                             double factor,std::size_t mask,const Sample& production) {
        // Construct affine coordinates without a double-rounded perturbed input.
        max_pressure_residual=0;max_pressure_limit=0;
        Wide T=c.problem.temperature.value,P=c.problem.pressure.value,epsilon=c.parameters.dispersion_energy_over_k[0];
        std::array<Wide,2> feed{Wide(c.problem.feed.values[0]),Wide(c.problem.feed.values[1])};
        std::array<Wide,4> coefficients{};
        if(!c.problem.reactions.empty()) {const auto& law=*c.problem.reactions[0].correlation;coefficients={Wide(law.a),Wide(law.b),Wide(law.c),Wide(law.d)};}
        for(std::size_t k=0;k<ds.size();++k) {
            const Wide scale=Wide((mask>>k)&1?1:-1)*Wide(ds[k].step)*Wide(factor);const auto& d=ds[k].value;
            T+=scale*d.temperature_k;P+=scale*d.pressure_pa;epsilon+=scale*d.active_parameters[0];
            for(std::size_t i=0;i<2;++i) feed[i]+=scale*d.feed_amounts_mol[i];
            if(!c.problem.reactions.empty()) for(std::size_t i=0;i<4;++i) coefficients[i]+=scale*d.reaction_coefficients[0][i];
            if(d.active_parameters.size()==2) CHECK(d.active_parameters[0]==d.active_parameters[1]);
        }
        const auto sys=state(T,P,epsilon,production.phase.rho);
        const Wide rt=Wide(R)*T;
        std::array<Wide,2> x{},g{},gt{},gtt{},gp{},gtp{},gpp{};
        Wide n=0;
        const auto thermal=evaluate_ideal_correlation(c.records.component_intervals[0][0].correlation,CppAD::AD<Wide>(T),1e5);
        const Wide h=CppAD::Value(thermal.h),s=CppAD::Value(thermal.s),cp=CppAD::Value(thermal.cp);
        g[0]=h-T*s+rt*CppAD::log(rt/Wide(1e5));gt[0]=-s+Wide(R)*(CppAD::log(rt/Wide(1e5))+1);gtt[0]=(Wide(R)-cp)/T;
        if(c.problem.reactions.empty()) {
            n=feed[0]+feed[1];x={feed[0]/n,feed[1]/n};g[1]=g[0];gt[1]=gt[0];gtt[1]=gtt[0];
        } else {
            const auto& law=*c.problem.reactions[0].correlation;
            const bool fixed=c.problem.reactions[0].reference->reference_pressure_pa.has_value();
            const auto src=fixed?state(T,Wide(*c.problem.reactions[0].reference->reference_pressure_pa),epsilon,production.phase.rho):sys;
            const Wide lnK=coefficients[0]+coefficients[1]/T+coefficients[2]*CppAD::log(T/Wide(law.reference_temperature))+coefficients[3]*T;
            const Wide lnKT=-coefficients[1]/(T*T)+coefficients[2]/T+coefficients[3];
            const Wide lnKTT=2*coefficients[1]/(T*T*T)-coefficients[2]/(T*T);
            const Wide K=CppAD::exp(lnK+sys.lambda-src.lambda);
            x[0]=2/(1+CppAD::sqrt(1+4*K));x[1]=K*x[0]*x[0];n=(feed[0]+2*feed[1])/(x[0]+2*x[1]);
            g[1]=2*g[0]-rt*lnK+rt*src.lambda;
            gt[1]=2*gt[0]-Wide(R)*(lnK+T*lnKT)+Wide(R)*(src.lambda+T*src.lt);
            gtt[1]=2*gtt[0]-Wide(R)*(2*lnKT+T*lnKTT)+Wide(R)*(2*src.lt+T*src.ltt);
            if(!fixed) {gp[1]=rt*src.lp;gtp[1]=Wide(R)*(src.lp+T*src.ltp);gpp[1]=rt*src.lpp;}
        }
        CHECK(CppAD::abs(x[0]+x[1]-1)<ldexpq(1,-100));
        CHECK(CppAD::abs(n-Wide(production.solved.amounts[0]))<Wide(1e-10));
        for(std::size_t i=0;i<2;++i) CHECK(CppAD::abs(x[i]-Wide(production.solved.mole_fractions[i]))<Wide(1e-10));
        CHECK(CppAD::abs(1/sys.v-Wide(production.phase.rho))<Wide(1e-6));
        std::array<Wide,5> ideal{},total{};
        Wide volume=sys.v,vt=sys.vt,vp=sys.vp;
        for(std::size_t i=0;i<2;++i) {
            ideal[0]+=x[i]*(g[i]-T*gt[i]);ideal[2]-=x[i]*(gt[i]+Wide(R)*CppAD::log(x[i]/sys.v));ideal[3]-=x[i]*T*gtt[i];
            volume+=x[i]*gp[i];vt+=x[i]*gtp[i];vp+=x[i]*gpp[i];
        }
        ideal[0]+=rt;ideal[2]+=R;ideal[3]+=R;
        ideal[1]=ideal[0]-rt-P*(volume-sys.v);
        ideal[4]=ideal[3]-R+T*vt*vt/vp-T*sys.vt*sys.vt/sys.vp;
        const Wide cvres=-Wide(R)*(2*T*sys.at+T*T*sys.att);
        const std::array<Wide,5> residual{-Wide(R)*T*T*sys.at-rt*sys.v*sys.av,-Wide(R)*T*T*sys.at,
            -Wide(R)*(sys.a+T*sys.at),cvres-T*sys.vt*sys.vt/sys.vp-R,cvres};
        for(std::size_t i=0;i<5;++i) total[i]=ideal[i]+residual[i];
        const std::vector<Wide> pure_x{Wide(1),Wide(0)},eps{epsilon,epsilon};
        ParameterRecord<Wide> geometry{{},Wide(0),Wide(0)};
        static_cast<FixedParameters&>(geometry)=c.parameters;
        geometry.dispersion_energy_over_k=eps;
        geometry.packing_diameter_angstrom.assign(c.parameters.packing_diameter_angstrom.begin(),c.parameters.packing_diameter_angstrom.end());
        const auto sigma=segment_sigma(T,c.parameters),diameter=segment_diameter(T,geometry,sigma);
        const Wide eta=packing_fraction_from_geometry(1/sys.v,pure_x,c.parameters,diameter);
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
            case SolvedStateObservableKind::PhasePressure:v=P+sys.pressure_residual;break;
            case SolvedStateObservableKind::PhaseMoleFraction:v=x[o.component];break;
            case SolvedStateObservableKind::PhaseComponentAmount:v=n*x[o.component];break;
            case SolvedStateObservableKind::PhaseResidualChemicalPotentialOverRt:v=sys.mu;break;
            case SolvedStateObservableKind::PhaseLogActivity:v=CppAD::log(x[o.component]/sys.v)+sys.mu;break;
            case SolvedStateObservableKind::PhaseGibbsMolarVolume:v=volume;break;
            case SolvedStateObservableKind::PhaseProperty:{const auto k=std::size_t(o.property);v=k<5?ideal[k]:total[k-10];break;}
            }
            if(item.name=="H") v=n*total[0];if(item.name=="U_G") v=n*total[1];
            if(item.name=="S") v=n*total[2];if(item.name=="V_G") v=n*volume;
            out.push_back(v);
        }
        return out;
    }
};
