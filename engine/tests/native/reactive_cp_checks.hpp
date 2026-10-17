#pragma once
// #225 F1: c_p against the qualified H actions; F2: quadrature against endpoint H and c_p; F3: against pressure actions.
inline void cp_record(const std::string& id, double exact, double reference, double target,
                      const std::string& inputs, bool resolved=true, const char* unit="J/(mol-feed K)") {
    const bool good=resolved&&std::isfinite(exact)&&std::isfinite(reference)&&std::abs(exact-reference)<=target;
    recording_support::compare_custom(id.c_str(),"apparent-Cp",unit,exact,reference,
        recording_support::fields({inputs,recording_support::field("target",target),
            recording_support::field("error",std::abs(exact-reference)),
            recording_support::field("error_over_target",std::abs(exact-reference)/target)}),
        "Relationship","native-identity","Issue 225 apparent Cp numerical verification",
        "abs(exact-reference)<=target",good,"reactive-cp",id.c_str());
}
inline double cp_feed(const Case& c) {return std::accumulate(c.problem.feed.values.begin(),c.problem.feed.values.end(),0.0);}
inline SolvedStateActionBatch cp_batch(const CompiledSystem& compiled,const EquilibriumResult& solved,const Direction& d) {
    const auto batch=solved_state_actions(compiled,solved,{{{SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity}}, {d}});
    CHECK(batch.results[0].value&&batch.results[0].action&&batch.results[0].status==SolvedStateActionStatus::Available);
    return batch;
}
// F1: c_p = H_T/N_f by definition, so its action along d is (H_Td - H_T N_f,d/N_f)/N_f. The H actions
// are the ones the E/I ladders qualify against binary128 values at the same centres.
void cp_f1(bool ionic_case,bool fixed) {
    const auto c=input(ionic_case,fixed);HelmholtzModel model(c.parameters,c.active,c.records);
    const auto compiled=compile_problem(model,c.problem);const auto solved=solve_compiled(compiled);CHECK(solved.success);
    const auto q=quantities(c.problem.feed.values.size());const std::size_t H=index(q,"H");
    auto ds=directions(c);const auto thermal=ds.front().value;
    auto mixed=ds[0];mixed.name="2T+P+feed+epsilon+source";mixed.value=add(ds[0].value,ds[0].value);
    for(std::size_t k=1;k<ds.size();++k) mixed.value=add(mixed.value,ds[k].value);
    ds.push_back(mixed);
    Actions old(compiled,solved,q);const double feed=cp_feed(c),h_t=old.get({thermal})[H];
    for(std::size_t test=0;test<=ds.size();++test) {
        const bool value=test==0;const auto d=value?zero(c):ds[test-1].value;
        const auto batch=cp_batch(compiled,solved,d);const auto& item=batch.results[0];
        const double exact=value?*item.value:*item.action;
        const double feed_d=std::accumulate(d.feed_amounts_mol.begin(),d.feed_amounts_mol.end(),0.0);
        const double h_td=value?h_t:old.get({thermal,d})[H],coupling=value?0.0:h_t*feed_d/feed;
        // 1e-10 of the terms, with a floor of 1e-16 c_p (the largest error seen is 4e-18 c_p). The floor binds only
        // on reaction-a, whose action is 1e-9 of c_p, so that row is held to about 1e-7 of its action.
        const double reference=(h_td-coupling)/feed,target=(1e-10*(std::abs(h_td)+std::abs(coupling))+1e-16*std::abs(h_t))/feed;
        cp_record("F1/"+c.name+"/apparent-Cp/"+(value?"value":ds[test-1].name),exact,reference,target,
            recording_support::fields({recording_support::field("H_T",h_t),recording_support::field("H_Td",h_td),
                recording_support::field("N_f_mol",feed),recording_support::field("dN_f_mol",feed_d),
                recording_support::field("action_condition",batch.diagnostics.condition_estimate.value_or(0))}));
        const auto h=solved_state_actions(compiled,solved,{{{SolvedStateObservableKind::SystemEnthalpy}}, {d}}).results[0];
        CHECK(h.value&&h.action);const auto reference_h=old.get(value?std::vector<Direction>{}:std::vector<Direction>{d})[H];
        CHECK(compare("preservation/"+c.name+"/H/"+std::to_string(test),q[H],value?*h.value:*h.action,reference_h,"",1e-10));
    }
}
void cp_f3() {
    const auto c=input(false,false,true);HelmholtzModel model(c.parameters,c.active,c.records);
    const auto compiled=compile_problem(model,c.problem);const auto solved=solve_compiled(compiled);CHECK(solved.success);
    std::vector<double> x=solved.mole_fractions;const auto ds=directions(c);
    for(const auto k:{0U,1U,3U}) {
        const auto batch=cp_batch(compiled,solved,ds[k].value);
        PressureActionDirection physical{ds[k].value.temperature_k,ds[k].value.pressure_pa,{0,0},{}};
        for(std::size_t a=0;a<c.active.size();++a) physical.parameter.push_back({c.active[a].family,c.active[a].indices,ds[k].value.active_parameters[a]});
        const auto p=model.pressure_actions(c.problem.temperature.value,c.problem.pressure.value,x,DensityBranch::Liquid,
            {PropertyObservable::TotalIsobaricHeatCapacity},{physical},solved.molar_densities[0]);
        CHECK(p.results[0].value&&p.results[0].action);
        cp_record("F3/U/apparent-Cp/value/"+ds[k].name,*batch.results[0].value,*p.results[0].value,1e-7+1e-10*std::abs(*p.results[0].value),"");
        cp_record("F3/U/apparent-Cp/"+ds[k].name,*batch.results[0].action,*p.results[0].action,1e-7+1e-10*std::abs(*p.results[0].action),"");
    }
}
inline std::vector<std::pair<double,double>> cp_gauss(std::size_t n) {
    std::vector<std::pair<double,double>> result;
    for(std::size_t i=0;i<n;++i) {
        double z=std::cos(3.14159265358979323846*(i+0.75)/(n+0.5)),derivative=0;
        for(int iteration=0;iteration<30;++iteration) {
            double p=1,previous=0;for(std::size_t k=1;k<=n;++k) {const double next=((2*k-1)*z*p-(k-1)*previous)/k;previous=p;p=next;}
            derivative=n*(z*p-previous)/(z*z-1);const double step=p/derivative;z-=step;
            if(std::abs(step)<=4*std::numeric_limits<double>::epsilon()) break;
            CHECK(iteration<29);
        }
        result.push_back({z,2/((1-z*z)*derivative*derivative)});
    }
    return result;
}
void cp_f2() {
    for(bool ionic:{false,true}) for(bool fixed:{false,true}) {
        auto c=input(ionic,fixed);HelmholtzModel model(c.parameters,c.active,c.records);
        const auto q=quantities(c.problem.feed.values.size());const double left=ionic?288.15:250,right=ionic?308.15:275;
        auto thermal=zero(c);thermal.temperature_k=1;std::array<double,2> endpoints_h{},endpoints_cp{};
        for(std::size_t i=0;i<2;++i) {
            c.problem.temperature.value=i?right:left;const auto compiled=compile_problem(model,c.problem);const auto solved=solve_compiled(compiled);CHECK(solved.success);
            const auto sample=values(model,compiled,solved,q);endpoints_h[i]=sample.values[index(q,"H")];endpoints_cp[i]=*cp_batch(compiled,solved,thermal).results[0].value;
            recording_support::observe(("F2/"+c.name+"/endpoint-"+std::to_string(i)).c_str(),"H","J",endpoints_h[i],endpoints_cp[i],context(c,solved,sample.phase),
                "Observation","independent-endpoint","H and Cp endpoint","reactive-cp",c.name.c_str());
        }
        // One Gauss-Legendre panel at 8 and 10 nodes (6 fails the gate on E-sys); the quarter-target refinement gate rejects a coarse rule.
        std::array<std::array<double,2>,2> integral{};
        for(std::size_t level=0;level<2;++level) {
            std::size_t node=0;for(const auto& [z,weight]:cp_gauss(level?10:8)) {
                c.problem.temperature.value=(left+right)/2+(right-left)*z/2;const auto compiled=compile_problem(model,c.problem);const auto solved=solve_compiled(compiled);CHECK(solved.success);
                const auto batch=cp_batch(compiled,solved,thermal);
                integral[level][0]+=(right-left)*weight/2 * *batch.results[0].value*cp_feed(c);
                integral[level][1]+=(right-left)*weight/2 * *batch.results[0].action;
                const auto id="F2/"+c.name+"/nodes-"+std::to_string(level?10:8)+"/node-"+std::to_string(node++);
                recording_support::observe(id.c_str(),"Cp","J/(mol-feed K)",*batch.results[0].value,*batch.results[0].action,
                    recording_support::fields({recording_support::field("T_K",*solved.temperature),recording_support::field("z",z),recording_support::field("weight",weight),
                        recording_support::field("N_f_mol",cp_feed(c)),recording_support::field("condition",batch.diagnostics.condition_estimate.value_or(0)),
                        recording_support::field("linear_residual",batch.diagnostics.maximum_linear_residual.value_or(0))}),
                    "Observation","native-Gauss-Legendre","Quadrature node Cp and thermal action","reactive-cp",id.c_str());
            }
        }
        for(std::size_t k=0;k<2;++k) {
            const double increment=k?endpoints_cp[1]-endpoints_cp[0]:endpoints_h[1]-endpoints_h[0];
            const double target=(k?1e-7:1e-5)+1e-8*std::max(std::abs(increment),std::abs(integral[1][k]));
            const double refinement=std::abs(integral[1][k]-integral[0][k]);
            cp_record("F2/"+c.name+(k?"/Cp-integral":"/H-integral"),integral[1][k],increment,target,
                recording_support::fields({recording_support::field("integral_8",integral[0][k]),recording_support::field("integral_10",integral[1][k]),
                    recording_support::field("refinement",refinement),recording_support::field("refinement_over_quarter_target",4*refinement/target),
                    recording_support::field("N_f_mol",cp_feed(c)),recording_support::field("identity_unit",k?"J/(mol-feed K)":"J")}),refinement<=target/4,k?"J/(mol-feed K)":"J");
        }
    }
}
