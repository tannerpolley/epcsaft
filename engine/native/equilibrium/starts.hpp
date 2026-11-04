// Included once, by gibbs_minimizer.cpp in namespace epcsaft::equilibrium: the multi-start equilibrium solve.
EquilibriumStartsResult solve_equilibrium_starts(const HelmholtzModel& model,
    const EquilibriumProblem& problem, const std::vector<EquilibriumStart>& starts,
    const EquilibriumWorkLimits& limits) {
    core::Scope scope("solve_equilibrium_starts");
    using Clock = std::chrono::steady_clock;
    const auto begun = Clock::now();
    const auto deadline_at = [](Clock::time_point now, const std::optional<double>& seconds) {
        if (seconds && (!std::isfinite(*seconds) || *seconds < 0.0))
            throw std::invalid_argument("work seconds must be finite and nonnegative");
        if (!seconds || *seconds >= std::chrono::duration<double>(Clock::time_point::max() - now).count())
            return Clock::time_point::max();
        return now + std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(*seconds));
    };
    const auto deadline = deadline_at(begun, limits.seconds);
    deadline_at(begun, limits.attempt_seconds);  // validate even an unattempted set
    EquilibriumStartsResult out;
    out.results.resize(starts.size()); out.causes.assign(starts.size(), "unattempted");
    out.attempt_seconds.resize(starts.size(), 0.0);
    const auto exhausted = [&] { return Clock::now() >= deadline
        || (limits.iterations && out.iterations >= *limits.iterations); };
    std::optional<CompiledSystem> retained;
    if (!exhausted()) {
        auto& compiled = retained.emplace(compile_problem(model, problem));
        for (const auto& seed : starts) if (seed.previous) {
            const auto& old = *seed.previous;
            if (!old.declaration || old.model_identity != &model)
                throw std::invalid_argument("previous start model/parameters or declaration differs");
            if (const auto difference = declared_state_difference(*old.declaration, compiled.problem); !difference.empty())
                throw std::invalid_argument("previous start: " + difference);
        }
        for (std::size_t i = 0; i < starts.size(); ++i) {
            if (exhausted() || (limits.attempts && i >= *limits.attempts)) { out.budget_exhausted = true; break; }
            core::Scope scope("equilibrium_start_attempt");
            scope.probe("declared_start_index", static_cast<double>(i));
            const auto& seed = starts[i];
            const auto attempt_begun = Clock::now();
            DriverLimits work{std::min(limits.attempt_iterations.value_or(static_cast<std::size_t>(policy::EQUILIBRIUM_MAX_ITERATIONS.value)),
                limits.iterations ? *limits.iterations - out.iterations : std::numeric_limits<std::size_t>::max()),
                std::min(deadline, deadline_at(Clock::now(), limits.attempt_seconds))};
            if (work.iterations == 0 || work.expired()) { out.causes[i] = "attempt budget exhausted"; continue; }
            std::optional<EquilibriumResult> audit;
            StageWork stages;
            try {
                auto initial = seed.coordinates;
                const auto offsets = compiled.phase_offsets();
                const auto check = [](bool ok, const std::string& field = "invalid declared seed") {
                    if (!ok) throw std::invalid_argument(field);
                };
                check((initial.empty() && !seed.previous) || (seed.compositions.empty() && seed.packing_fractions.empty()
                    && seed.amounts.empty() && !seed.temperature && !seed.pressure && !(seed.previous && !initial.empty())));
                if (seed.previous) { check(seed.previous->has_state()); initial = seed.previous->coordinates; }
                else if (initial.empty()) {
                    auto seeded = compiled;
                    check(seed.compositions.empty() || seed.compositions.size() == compiled.problem.phases.size());
                    check(seed.packing_fractions.empty() || seed.packing_fractions.size() == compiled.problem.phases.size());
                    check(seed.amounts.empty() || seed.amounts.size() == compiled.problem.phases.size());
                    // Complete: the seed supplies every coordinate the plain solve's cold start would generate.
                    bool uncomposed = true, complete = (!seed.amounts.empty() || compiled.problem.phases.size() == compiled.pinned_count)
                        && (seed.temperature || !problem.temperature.free) && (seed.pressure || !problem.pressure.free);
                    const auto composed = [&](std::size_t p) { return !seed.compositions.empty() && !seed.compositions[p].empty(); };
                    const auto packed = [&](std::size_t p) { return !seed.packing_fractions.empty() && seed.packing_fractions[p]; };
                    for (std::size_t p = 0; p < seeded.problem.phases.size(); ++p) {
                        auto& phase = seeded.problem.phases[p];
                        if (composed(p)) {
                            check(seed.compositions[p].size() == compiled.declared_supports[p].size()
                                  && std::all_of(seed.compositions[p].begin(), seed.compositions[p].end(), [](double x) { return std::isfinite(x) && x >= 0.0 && x <= 1.0; }));
                            phase.composition_guess = seed.compositions[p];
                            phase.declared_composition_guess = false;
                            restrict_composition_guess(phase, compiled.declared_supports[p]);
                        }
                        check(!model.parameters().ideal_gas || !packed(p), "ideal-gas packing seeds are undefined");
                        if (packed(p)) phase.packing_guess = seed.packing_fractions[p];
                        check(phase.composition_guess.empty() || phase.composition_guess.size() == phase.support.size());
                        for (const double x : phase.composition_guess)
                            check(std::isfinite(x) && x > 0.0 && x <= 1.0,
                                  "composition for phase " + phase.name + " is outside the admitted range");
                        check(!phase.packing_guess || admitted_packing_fraction(*phase.packing_guess),
                              "packing for phase " + phase.name + " is outside the admitted range");
                        if (!seed.amounts.empty() && (!std::isfinite(seed.amounts[p]) || seed.amounts[p] < 0.0
                            || (phase.is_pinned() && seed.amounts[p] != compiled.phase_amount(phase))))
                            throw std::invalid_argument("amount for phase " + phase.name + " is nonfinite, negative or conflicts with its pin");
                        complete = complete && composed(p) && packed(p);
                        uncomposed = uncomposed && !composed(p);
                    }
                    if (seed.temperature) { check(problem.temperature.free); seeded.problem.temperature = FreeValue::variable(*seed.temperature); }
                    if (seed.pressure) { check(problem.pressure.free); seeded.problem.pressure = FreeValue::variable(*seed.pressure); }
                    if (work.expired()) throw std::runtime_error("elapsed work budget exhausted in initialization");
                    auto speciated = complete || !uncomposed ? std::nullopt : speciated_start(seeded, &work, stages);
                    std::vector<bool> packings(seeded.problem.phases.size());
                    for (std::size_t p = 0; p < packings.size(); ++p) packings[p] = packed(p);
                    initial = complete ? seeded.generic_point() : speciated ? std::move(*speciated)
                        : boundary_start(seeded, uncomposed, &work, &stages, packings);
                    for (std::size_t p = 0; p < seeded.problem.phases.size(); ++p) {
                        const auto& phase = seeded.problem.phases[p];
                        if (composed(p)) for (std::size_t j = 0; j < phase.composition_guess.size(); ++j)
                            initial[offsets[p] + j] = std::log(phase.composition_guess[j]);
                        auto cursor = offsets[p] + phase.support.size();
                        if (phase.is_free()) { if (!seed.amounts.empty()) initial[cursor] = seed.amounts[p]; ++cursor; }
                        if (packed(p)) initial[cursor] = std::log(*phase.packing_guess);
                    }
                    auto cursor = compiled.primal_size - compiled.free_intensive_count;
                    if (problem.temperature.free) { if (seed.temperature) initial[cursor] = *seed.temperature; ++cursor; }
                    if (problem.pressure.free && seed.pressure) initial[cursor] = *seed.pressure;
                }
                if (work.consumed >= work.iterations) throw std::runtime_error("iteration budget exhausted by the seed solves");
                check(initial.size() == compiled.unknown_count() && all_numeric(initial));
                for (std::size_t p = 0; p < offsets.size(); ++p) {
                    const auto& phase = compiled.problem.phases[p]; auto cursor = offsets[p];
                    for (std::size_t j = 0; j < phase.support.size(); ++j) {
                        const double x = std::exp(initial[cursor++]);
                        check(x > 0.0);
                    }
                    if (phase.is_free()) check(initial[cursor++] >= 0.0);
                    check(model.parameters().ideal_gas ? std::isfinite(std::exp(initial[cursor])) && std::exp(initial[cursor]) > 0.0
                          : admitted_packing_fraction(std::exp(initial[cursor])));
                }
                for (auto j = compiled.primal_size - compiled.free_intensive_count; j < compiled.primal_size; ++j) check(initial[j] > 0.0);
                audit = structural_result(compiled, scope);
                audit->convergence.residual = std::numeric_limits<double>::quiet_NaN();
                auto result = compiled.feasible ? run_driver(model, compiled, std::move(initial), *audit, scope, nullptr, &*audit, &work)
                                                : *audit;
                stages.report(result);
                if (work.expired()) { result.message += "; elapsed work budget exhausted";
                    reject_elapsed(result, scope); }
                if (result.success) {
                    BranchCheck branch;
                    try { branch = compiled.problem.critical ? admit_critical(compiled, result)
                              : check_density_branches(model, compiled, result, {}, true, scope, &work); }
                    catch (const std::exception& error) { branch.message = error.what(); }
                    if (!branch.accepted || work.expired()) {
                        result.success = false;
                        result.message = (work.expired() ? "elapsed work budget exhausted in admission" : branch.message) + "; driver: " + result.message;
                        if (work.expired()) reject_elapsed(result, scope); else clear_state(result);
                    }
                }
                scope.finish(result.convergence);
                out.causes[i] = result.message; out.results[i] = std::move(result);
                if (out.results[i]->success && out.results[i]->requested_tolerance_met.value_or(false)) out.selected = i;
            } catch (const std::exception& error) {
                out.causes[i] = std::string(audit ? "driver/post-driver: " : "initialization: ") + error.what();
                if (work.expired()) out.causes[i] += "; elapsed deadline exhausted";
                if (audit) { audit->message = out.causes[i]; stages.report(*audit);
                    if (work.expired()) reject_elapsed(*audit, scope); else clear_state(*audit);
                    out.results[i] = std::move(audit); }
            }
            if (out.results[i]) for (auto& rank : out.results[i]->rank_records)
                if (rank.first_conditioning_warning && !rank.first_conditioning_warning->attempt) (rank.first_conditioning_warning = std::make_shared<core::ConditioningObservation>(*rank.first_conditioning_warning))->attempt = i;
            out.iterations += work.consumed;
            out.attempt_seconds[i] = std::chrono::duration<double>(Clock::now() - attempt_begun).count();
            if (out.selected) break;
        }
    } else out.budget_exhausted = true;
    out.budget_exhausted = out.budget_exhausted || (!out.selected && exhausted());
    out.seconds = std::chrono::duration<double>(Clock::now() - begun).count();
    out.overrun_seconds = limits.seconds ? std::max(0.0, out.seconds - *limits.seconds) : 0.0;
    if (limits.attempt_seconds) for (const double t : out.attempt_seconds)
        out.overrun_seconds = std::max(out.overrun_seconds, t - *limits.attempt_seconds);
    if (scope.diagnostics_on() && retained) for (std::size_t i = 0; i < out.results.size(); ++i)
        if (out.results[i]) observe_endpoint(model, &*retained, *out.results[i], scope, i);
    return out;
}
