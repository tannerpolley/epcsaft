#pragma once
// #225 F7 uses the existing literal E declarations and native fitter.
void cp_column(const std::string& id,double exact,double sigma,double step,
               const std::function<double(double)>& value,const std::string& inputs) {
    const auto difference=[&](double h){return (value(h)-value(-h))/(2*h);};
    double coarse=difference(step);bool passed=false;
    for(std::size_t rung=0;rung<6;++rung) {
        step/=2;const double fine=difference(step),rich=(4*fine-coarse)/3,refinement=std::abs(fine-coarse)/3;
        const double target=1e-7+1e-6*std::max(std::abs(rich),sigma);const bool resolved=refinement<=target/4;
        passed=resolved&&std::isfinite(exact)&&std::abs(exact-rich)<=target;
        const auto row=fields({inputs,field("step",step),field("D_h",coarse),field("D_h2",fine),field("Richardson",rich),
            field("refinement",refinement),field("sigma",sigma),field("target",target),field("resolved",resolved),
            field("error",std::abs(exact-rich)),field("error_over_target",std::abs(exact-rich)/target)});
        const auto name=id+"/rung-"+std::to_string(rung);
        if(!resolved&&rung<5) recording_support::observe(name, "Cp-column", "native", exact, rich,
            "rebuilt-native-prediction", row);
        else recording_support::compare_custom(name, "Cp-column", "native", exact, rich,
            "abs(exact-R)<=1e-7+1e-6 max(abs(R),sigma); refinement<=target/4", passed, "rebuilt-native-prediction", row);
        if(resolved) break;coarse=fine;
    }
}
void check_reactive_cp() {
    FitRequest request;request.parameters=addition_parameters();request.thermochemistry=addition_record();
    request.controls=budgeted_controls();
    auto epsilon=addition_coordinates().coordinates[1];epsilon.origin=191.42;epsilon.scale=2;epsilon.lower=180;epsilon.upper=200;
    auto b=reaction_coordinate("2A=B",ReactionCoefficient::B);b.origin=200;b.scale=10;b.lower=100;b.upper=300;
    request.coordinates.coordinates={epsilon,b};
    const HelmholtzModel truth_model(request.parameters,request.coordinates.active_parameters(),request.thermochemistry);
    for(double T:{255.,260.,265.,270.,275.}) {
        auto ob=solved(ObservationKind::EquilibriumIsobaricHeatCapacity,addition_problem(T,3e6,{2,0.2},5e6));
        ob.form=ResidualForm::Difference;ob.scale=1;const auto p=predict(truth_model,ob);
        CHECK(p.status==ObservationStatus::Available&&p.value);ob.observed=*p.value;request.observations.push_back(ob);
    }
    {
        using namespace epcsaft::equilibrium;
        using Limit = core::WorkLimit;
        const auto& row = request.observations[2];
        const std::vector<std::string> ids{"A", "B"};
        const auto base = predict_validated(truth_model, row, {}, nullptr, nullptr, ids);
        CHECK(base.solved && base.value);
        const auto same = predict_validated(truth_model, row, {}, &base, nullptr, ids);
        CHECK(same.value && std::abs(*same.value / *base.value - 1) <= 1e-10);
        EquilibriumWorkLimits limits;
        limits.attempt_iterations = 1;
        const auto polished = predict_validated(truth_model, row, {}, nullptr, &limits);
        CHECK(polished.value && polished.solved && polished.solved->iterations == 1
            && polished.solved->attempt == Attempt::Cold && !polished.refusal_convergence
            && std::abs(*polished.value / *base.value - 1) <= 1e-10);
        limits.attempt_iterations = 800;
        const auto larger = predict_validated(truth_model, row, {}, nullptr, &limits);
        CHECK(larger.solved && larger.solved->solver_settings.at("max_iter") == 800 && larger.value == base.value);
        limits = {}; limits.seconds = 0;
        check_prediction_limit(predict_validated(truth_model, row, {}, nullptr, &limits),
            RefusedOperation::central_solve, Limit::deadline, false);
        const auto no_calorics = predict_validated(HelmholtzModel(request.parameters), row);
        CHECK(no_calorics.status == ObservationStatus::ThermochemistryUnavailable
            && no_calorics.refused_operation == RefusedOperation::thermochemistry && no_calorics.refusal_convergence
            && std::isnan(no_calorics.refusal_convergence->residual) && !no_calorics.refusal_convergence->criterion);
        auto changed = row;
        changed.problem.reactions[0].correlation->b = 220;
        CHECK(!same_reaction(row.problem.reactions[0], changed.problem.reactions[0]));
        auto parameters = request.parameters; parameters.dispersion_energy_over_k = {195.25, 195.25};
        const HelmholtzModel other(parameters, request.coordinates.active_parameters(), request.thermochemistry);
        const auto cold = predict_validated(other, changed);
        const auto started = predict_validated(other, changed, {}, &base, nullptr, ids);
        CHECK(cold.value && started.value && std::abs(*cold.value / *started.value - 1) <= 1e-10
            && started.solved->attempt != Attempt::PhysicalAmounts);
        auto reject = [&](const Observation& item, const Prediction& seed, const auto& components) {
            IpoptLog log;
            CHECK(recording_support::refusal([&] { (void)predict_validated(other, item, {}, &seed, nullptr, components); })
                .find("start:") == 0 && log.rows.empty());
        };
        reject(changed, base, std::vector<std::string>{"B", "A"});
        auto unavailable = base; unavailable.status = ObservationStatus::CentralSolveUnavailable;
        reject(changed, unavailable, ids);
        auto zero = base; zero.solved->mole_fractions.front() = 0;
        reject(changed, zero, ids);
        auto altered = base;
        auto declaration = *base.solved->declaration;
        for (auto change : {0, 1, 2, 3, 4}) {
            auto incompatible = declaration;
            if (change == 0) incompatible.feed.values.front() *= 2;
            if (change == 1) incompatible.temperature.value += 1;
            if (change == 2) incompatible.reactions[0].stoichiometry.front() -= 1;
            if (change == 3) incompatible.reactions[0].reference->reference_pressure_pa = 3e6;
            if (change == 4) incompatible.phases[0].support.front() = incompatible.phases[0].support.back();
            altered.solved->declaration = std::make_shared<EquilibriumProblem>(std::move(incompatible));
            reject(changed, altered, ids);
        }
    }
    request.weights.assign(5,1);request.parameters.dispersion_energy_over_k={195.25,195.25};
    for(auto& ob:request.observations) ob.problem.reactions[0].correlation->b=220;
    const auto start=request.coordinates.start(request.parameters,request.observations,request.thermochemistry);
    const auto result=fit(request);
    const auto evidence=fields({summary(result),field("true_epsilon_K",191.42),field("true_b_K",200),field("start_epsilon_K",195.25),field("start_b_K",220),
        field("bounds_epsilon_K",std::vector<double>{180,200}),field("bounds_b_K",std::vector<double>{100,300}),field("weights",request.weights),field("Difference_scale",1),
        field("physical",result.physical),field("predictions",result.predictions),field("residuals",result.residuals),
        field("epsilon_relative_error",std::abs(result.physical[0]-191.42)/191.42),
        field("b_relative_error",std::abs(result.physical[2]-200)/200),
        field("singular_values",result.covariance.singular_values),field("rank",result.covariance.rank),field("active_bounds",std::vector<int>(result.covariance.active_bounds.begin(),result.covariance.active_bounds.end())),
        field("covariance_status",covariance_status(result.covariance)),field("covariance_assumption",result.covariance.assumption),
        field("conditional_covariance",result.covariance.physical),field("variance_factor",result.covariance.variance_factor.value_or(-1)),
        field("optimizer_jacobian",result.optimizer_jacobian),field("physical_jacobian",result.physical_jacobian),
        field("iteration_seconds",result.iteration_seconds),field("iteration_costs",result.iteration_costs)});
    truth("F7/recovery/rank",result.covariance.rank==2
        &&result.covariance.singular_values.size()==2&&result.covariance.active_bounds.empty()&&result.covariance.outcome == FitCovarianceStatus::Available,evidence);
    agree("F7/recovery/epsilon","epsilon","K",result.physical[0],191.42,0,1e-6,evidence);
    agree("F7/recovery/b","b","K",result.physical[2],200,0,1e-6,evidence);
    agree("F7/recovery/cost","cost","1",result.final_cost,0,1e-20,0,evidence);
    auto columns=request.coordinates;
    for(const auto k:{ReactionCoefficient::A,ReactionCoefficient::C,ReactionCoefficient::D}) columns.coordinates.push_back(reaction_coordinate("2A=B",k));
    for(std::size_t point=0;point<2;++point) {
        const auto optimizer=point?result.optimizer:start;auto rows=request.observations;request.coordinates.apply(rows,optimizer);
        const auto parameters=request.coordinates.apply(request.parameters,optimizer);
        const HelmholtzModel model(parameters,columns.active_parameters(),request.thermochemistry);
        const auto mapped=evaluated_columns(request,optimizer);
        for(std::size_t row=0;row<rows.size();++row) {
            const auto& ob=rows[row];const auto prediction=predict(model,ob);CHECK(prediction.value&&prediction.solved&&prediction.compiled);
            regular_reaction_state(model,ob,"F7/chart");
            Observation h=ob;h.kind=ObservationKind::SystemEnthalpy;auto hp=prediction;project_prediction(model,h,hp);
            const auto grouped=differentiate_group(model,{&h,&ob},prediction,columns);
            const auto& gradient=grouped[1];CHECK(gradient.status==ObservationStatus::Available&&gradient.value.size()==6);
            const auto diagnostics=fields({field("T_K",ob.problem.temperature.value),field("P_Pa",ob.problem.pressure.value),field("source_P_Pa",5e6),
                field("feed_mol",ob.problem.feed.values),field("x",prediction.solved->mole_fractions),field("rho",prediction.solved->molar_densities),
                field("solver_residual",prediction.solved->residual_norm()),field("iterations",prediction.solved->iterations),field("evaluations",prediction.solved->evaluations),
                field("prediction",*prediction.value),field("observed",ob.observed),field("physical_columns",gradient.value),field("optimizer_point",optimizer)});
            for(std::size_t k=0;k<7;++k) {
                const bool tied=k==6;const double exact=tied?gradient.value[0]+gradient.value[1]:gradient.value[k];
                const double sigma=tied?std::abs(gradient.value[0])+std::abs(gradient.value[1]):std::abs(exact);
                const double step=k<2||tied?0.02:k==2?2:k==5?1e-5:1e-3;
                const auto value=[&](double shift) {
                    auto p=parameters;auto moved=ob;
                    if(k<2) p.dispersion_energy_over_k[k]+=shift;
                    else if(tied) for(auto& v:p.dispersion_energy_over_k) v+=shift;
                    else {auto& law=*moved.problem.reactions[0].correlation;
                        if(k==2) law.b+=shift;else if(k==3) law.a+=shift;else if(k==4) law.c+=shift;else law.d+=shift;}
                    const auto prediction=predict(HelmholtzModel(p,columns.active_parameters(),request.thermochemistry),moved);
                    CHECK(prediction.status==ObservationStatus::Available&&prediction.value&&prediction.solved);
                    return *prediction.value;
                };
                cp_column("F7/"+std::string(point?"solution":"start")+"/T-"+number(ob.problem.temperature.value)+"/physical-"+std::to_string(k),exact,sigma,step,value,diagnostics);
            }
            for(std::size_t k=0;k<2;++k) {
                const auto value=[&](double shift) {auto v=optimizer;v[k]+=shift;return evaluated_columns(request,v).predictions[row];};
                const double sigma=k==0?epsilon.scale*(std::abs(gradient.value[0])+std::abs(gradient.value[1])):b.scale*std::abs(gradient.value[2]);
                cp_column("F7/"+std::string(point?"solution":"start")+"/T-"+number(ob.problem.temperature.value)+"/optimizer-"+std::to_string(k),
                    mapped.residual_optimizer_gradient[2*row+k],sigma,k?2/b.scale:0.02/epsilon.scale,value,diagnostics);
            }
            agree("F7/grouped/H/"+std::to_string(point)+"/"+std::to_string(row),"H","J",*hp.value,
                endpoint_enthalpy(model,h.problem,*hp.compiled,*hp.solved),1e-5,1e-10,diagnostics);
        }
    }
    auto signed_row=request.observations.front();signed_row.form=ResidualForm::LogRatio;
    auto negative_record=*request.thermochemistry;
    negative_record.component_intervals[0][0].correlation.form=IdealPolynomial{{-1000.0},300,0,0};
    const auto negative=predict(HelmholtzModel(request.parameters,{},negative_record),signed_row);
    truth("F7/nonpositive-log",
        negative.value&&*negative.value<0&&negative.status==ObservationStatus::NonpositiveLogArgument,
        fields({field("value",negative.value.value_or(0)),field("message",negative.message)}));
}
