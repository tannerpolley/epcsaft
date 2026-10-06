#pragma once
// #225 F1/F2/F3: references reconstruct H values; no production H A1 is used.
inline void cp_record(const std::string& id, double exact, double reference, double target,
                      const std::string& inputs, bool resolved=true, const char* unit="J/(mol-feed K)") {
    const bool good=resolved&&std::isfinite(exact)&&std::isfinite(reference)&&std::abs(exact-reference)<=target;
    recording_support::compare_custom(id.c_str(),"apparent-Cp",unit,exact,reference,
        recording_support::fields({inputs,recording_support::field("target",target),
            recording_support::field("error",std::abs(exact-reference)),
            recording_support::field("error_over_target",std::abs(exact-reference)/target)}),
        "Relationship","native-rebuilt-values","Issue 225 apparent Cp numerical verification",
        "absolute plus relative; quarter-target refinement",good,"reactive-cp",id.c_str());
    CHECK(good);
}
inline double cp_feed(const Case& c) {return std::accumulate(c.problem.feed.values.begin(),c.problem.feed.values.end(),0.0);}
inline SolvedStateActionBatch cp_batch(const CompiledSystem& compiled,const EquilibriumResult& solved,const Direction& d) {
    const auto batch=solved_state_actions(compiled,solved,{{{SolvedStateObservableKind::EquilibriumIsobaricHeatCapacity}}, {d}});
    CHECK(batch.results[0].value&&batch.results[0].action&&batch.results[0].status==SolvedStateActionStatus::Available);
    return batch;
}
inline Wide cp_difference(const Case& c,const Sample& center,const std::vector<Quantity>& q,
                          const std::vector<NamedDirection>& ds,double factor,WideNeutral* wide,WideIonic* ionic,
                          const std::string& id) {
    Wide sum=0,denominator=1;
    for(const auto& d:ds) denominator*=Wide(2)*Wide(d.step)*Wide(factor);
    for(std::size_t mask=0;mask<(1U<<ds.size());++mask) {
        auto delta=zero(c);Wide feed=cp_feed(c);int sign=1;
        for(std::size_t k=0;k<ds.size();++k) {
            const int s=(mask>>k)&1?1:-1;sign*=s;
            delta=add(delta,ds[k].value,s*ds[k].step*factor);
            for(double v:ds[k].value.feed_amounts_mol) feed+=Wide(s)*Wide(ds[k].step)*Wide(factor)*Wide(v);
        }
        const auto sample=shifted(c,delta,q,center.solved);
        const auto precise=wide?wide->values(q,ds,factor,mask,sample):ionic->values(q,ds,factor,mask,sample);
        const Wide h=precise[index(q,"H")];sum+=sign*h/feed;
        const auto label=id+"/factor-"+recording_support::number(factor)+"/sign-"+std::to_string(mask);
        recording_support::observe(label.c_str(),"independent-H","J",double(h),std::nullopt,
            recording_support::fields({context(c,sample.solved,sample.phase),recording_support::field("N_f_mol",double(feed)),
                recording_support::field("delta_T_K",delta.temperature_k),recording_support::field("delta_P_Pa",delta.pressure_pa),
                recording_support::field("delta_abcd",std::vector<double>(delta.reaction_coefficients[0].begin(),delta.reaction_coefficients[0].end())),
                recording_support::field("delta_feed_mol",delta.feed_amounts_mol),recording_support::field("delta_epsilon_K",delta.active_parameters),
                recording_support::field("wide_pressure_residual_Pa",wide?wide->pressure_residual():ionic->pressure_residual()),
                recording_support::field("wide_pressure_limit_Pa",wide?wide->pressure_limit():ionic->pressure_limit())}),
            "Observation","native-independent-H","Rebuilt normalized H reference","reactive-cp",label.c_str());
    }
    return sum/denominator;
}
inline void cp_f1(bool ionic_case,bool fixed) {
    const auto c=input(ionic_case,fixed);HelmholtzModel model(c.parameters,c.active,c.records);
    const auto compiled=compile_problem(model,c.problem);const auto solved=solve_compiled(compiled);CHECK(solved.success);
    const auto q=quantities(c.problem.feed.values.size());const auto center=values(model,compiled,solved,q);
    std::unique_ptr<WideNeutral> wide;std::unique_ptr<WideIonic> ionic;
    if(c.ionic) ionic=std::make_unique<WideIonic>(c,center);else wide=std::make_unique<WideNeutral>(c);
    auto ds=directions(c);const auto thermal=ds.front();
    auto mixed=ds[0];mixed.name="2T+P+feed+epsilon+source";mixed.value=add(ds[0].value,ds[0].value);
    mixed.chart={ds[0].value,ds[0].value};
    for(std::size_t k=1;k<ds.size();++k) {mixed.value=add(mixed.value,ds[k].value);mixed.chart.insert(mixed.chart.end(),ds[k].chart.begin(),ds[k].chart.end());}
    mixed.step=0.0001;ds.push_back(mixed);
    Actions old(compiled,solved,q);
    for(std::size_t test=0;test<=ds.size();++test) {
        const bool value=test==0;const auto d=value?zero(c):ds[test-1].value;
        const auto batch=cp_batch(compiled,solved,d);const auto& item=batch.results[0];
        const double exact=value?*item.value:*item.action;double sigma=std::abs(exact);
        if(!value) {sigma=0;for(const auto& part:ds[test-1].chart) sigma+=std::abs(*cp_batch(compiled,solved,part).results[0].action);}
        const std::vector<NamedDirection> stencil=value?std::vector<NamedDirection>{thermal}:std::vector<NamedDirection>{thermal,ds[test-1]};
        const auto id="F1/"+c.name+"/apparent-Cp/"+(value?"value":ds[test-1].name);
        auto coarse=cp_difference(c,center,q,stencil,1,wide.get(),ionic.get(),id);bool resolved=false;
        for(int rung=0;rung<6;++rung) {
            const double factor=std::ldexp(1.0,-rung-1);const auto fine=cp_difference(c,center,q,stencil,factor,wide.get(),ionic.get(),id);
            const double rich=double((4*fine-coarse)/3),refinement=double(CppAD::abs(fine-coarse)/3);
            const double target=1e-7+1e-6*std::max(std::abs(rich),sigma);resolved=refinement<=target/4;
            const auto row=recording_support::fields({context(c,solved,center.phase),recording_support::field("N_f_mol",cp_feed(c)),
                recording_support::field("D_h",double(coarse)),recording_support::field("D_h2",double(fine)),
                recording_support::field("Richardson",rich),recording_support::field("refinement",refinement),
                recording_support::field("sigma",sigma),recording_support::field("rung",rung),recording_support::field("resolved",resolved),
                recording_support::field("outer_step",value?0:ds[test-1].step),recording_support::field("thermal_step",thermal.step),
                recording_support::field("action_condition",batch.diagnostics.condition_estimate.value_or(0)),
                recording_support::field("linear_residual",batch.diagnostics.maximum_linear_residual.value_or(0))});
            if(resolved||rung==5) {cp_record(id,exact,rich,target,row,resolved);break;}
            recording_support::observe((id+"/rung-"+std::to_string(rung)).c_str(),"refinement","J/(mol-feed K)",refinement,target/4,row,
                "Refinement","native-independent-H","Unresolved rung retained","reactive-cp",id.c_str());coarse=fine;
        }
        const auto h=solved_state_actions(compiled,solved,{{{SolvedStateObservableKind::SystemEnthalpy}}, {d}}).results[0];
        CHECK(h.value&&h.action);const auto reference=old.get(value?std::vector<Direction>{}:std::vector<Direction>{d})[index(q,"H")];
        CHECK(compare("preservation/"+c.name+"/H/"+std::to_string(test),q[index(q,"H")],value?*h.value:*h.action,reference,"",1e-10));
    }
}
inline void cp_f3() {
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
inline void cp_f2() {
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
        std::array<std::array<double,2>,2> integral{};
        for(std::size_t level=0;level<2;++level) for(std::size_t panel=0;panel<4;++panel) {
            const double a=left+(right-left)*panel/4,b=left+(right-left)*(panel+1)/4;
            std::size_t node=0;for(const auto& [z,weight]:cp_gauss(level?16:8)) {
                c.problem.temperature.value=(a+b)/2+(b-a)*z/2;const auto compiled=compile_problem(model,c.problem);const auto solved=solve_compiled(compiled);CHECK(solved.success);
                const auto batch=cp_batch(compiled,solved,thermal);const auto sample=values(model,compiled,solved,q);
                integral[level][0]+=(b-a)*weight/2 * *batch.results[0].value*cp_feed(c);
                integral[level][1]+=(b-a)*weight/2 * *batch.results[0].action;
                const auto id="F2/"+c.name+"/nodes-"+std::to_string(level?16:8)+"/panel-"+std::to_string(panel)+"/node-"+std::to_string(node++);
                recording_support::observe(id.c_str(),"Cp","J/(mol-feed K)",*batch.results[0].value,*batch.results[0].action,
                    recording_support::fields({context(c,solved,sample.phase),recording_support::field("z",z),recording_support::field("weight",weight),
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
                recording_support::fields({recording_support::field("integral_8",integral[0][k]),recording_support::field("integral_16",integral[1][k]),
                    recording_support::field("refinement",refinement),recording_support::field("refinement_over_quarter_target",4*refinement/target),
                    recording_support::field("N_f_mol",cp_feed(c)),recording_support::field("identity_unit",k?"J/(mol-feed K)":"J")}),refinement<=target/4,k?"J/(mol-feed K)":"J");
        }
    }
}
