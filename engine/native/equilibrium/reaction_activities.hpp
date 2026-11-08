// Included once, by gibbs_minimizer.cpp in namespace epcsaft::equilibrium: species activities and ln Q at a returned state.
namespace {
// EqID reaction_activity_quotient: products-positive sum of nu_i ln a_i; a species cause or opposing infinities leave `cause`.
std::optional<double> log_quotient(const std::vector<double>& nu, const SpeciesActivities& species, std::string& cause) {
    double sum = 0.0;
    for (std::size_t i = 0; cause.empty() && i < nu.size(); ++i) {
        if (nu[i] == 0) continue;
        if (!species.log_activities[i]) cause = species.causes[i];
        else sum += nu[i] * *species.log_activities[i];
    }
    if (std::isnan(sum)) cause = "conflicting infinite log activity terms";
    return cause.empty() ? std::optional<double>(sum) : std::nullopt;
}

// EqIDs reaction_activity_quotient (pressure basis b = P), standard_pressure_density_conversion:
// ln Q_P = ln Q_rho + sum(nu) ln(RT/P0). Actions are partial derivatives at the phase state along T, P and e_j - x
// (amount j at fixed T, P), from the exact log-fugacity actions plus the ln x_i and ln P terms.
void pressure_basis(const HelmholtzModel& model, const EquilibriumResult& result, std::size_t phase,
    const std::vector<double>& x, const std::vector<double>& nu, const SpeciesActivities& density, ReactionActivities& row) {
    const auto& record = model.thermochemistry();
    row.standard_pressure_pa = row.reference_pressure_pa ? row.reference_pressure_pa
        : record ? std::optional<double>(record->standard_pressure) : std::nullopt;
    if (!row.standard_pressure_pa) { row.pressure_cause = "standard pressure undeclared"; return; }
    const auto ln_q = log_quotient(nu, density, row.pressure_cause);
    if (!ln_q) return;
    const double T = *result.temperature, P = *result.pressure;
    double dnu = 0.0;
    for (const double v : nu) dnu += v;
    row.log_pressure_quotient = *ln_q + dnu * std::log(GAS_CONSTANT_J_PER_MOL_K * T / *row.standard_pressure_pa);
    if (!std::isfinite(*ln_q)) { row.pressure_cause = "pressure-basis actions undefined at a structurally zero activity"; return; }
    const auto branch = result.declaration->phases[phase].kind == PhaseKind::Vapor ? DensityBranch::Vapor : DensityBranch::Liquid;
    row.log_pressure_quotient_actions.assign(x.size() + 2, std::nullopt);
    for (std::size_t k = 0; k < x.size() + 2; ++k) try {
        PressureActionDirection direction{double(k == 0), double(k == 1), std::vector<double>(x.size(), 0.0), {}};
        if (k > 1) for (std::size_t i = 0; i < x.size(); ++i) direction.composition[i] = (i == k - 2) - x[i];
        const auto batch = model.pressure_actions(T, P, x, branch, {PropertyObservable::LogFugacityCoefficient}, {direction});
        double action = k == 1 ? dnu / P : 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) if (nu[i] != 0)
            action += nu[i] * (batch.results.at(i).action.value() + (k > 1 ? double(i == k - 2) / x[i] - 1.0 : 0.0));
        row.log_pressure_quotient_actions[k] = action;
    } catch (const std::exception& error) { row.pressure_cause = error.what(); }
}

std::vector<ReactionActivities> read_reaction_activities(const HelmholtzModel* model,
    const IdealMoleFraction* ideal, const EquilibriumResult& result) {
    core::Scope scope("reaction_activities");
    if ((model && result.model_identity != model) || (ideal && result.ideal_identity != ideal))
        throw std::invalid_argument("reaction activities: model/result mismatch");
    const auto n = model ? model->parameters().component_count() : ideal->component_ids.size();
    if (!result.declaration || !result.temperature || !result.pressure || !std::isfinite(*result.temperature)
        || *result.temperature <= 0 || !std::isfinite(*result.pressure) || *result.pressure <= 0
        || result.phase_count == 0 || result.phase_count != result.declaration->phases.size()
        || result.amounts.size() != result.phase_count || result.mole_fractions.size() != n * result.phase_count
        || (model && result.molar_densities.size() != result.phase_count))
        throw std::invalid_argument("reaction activities: physical state or declaration unavailable");
    require_solution_shape(result, result.phase_count, n);
    const auto& problem = *result.declaration;
    std::vector<std::vector<double>> rows;
    std::vector<ReactionReference> references;
    std::vector<std::size_t> referenced;
    for (std::size_t r = 0; r < problem.reactions.size(); ++r) {
        const auto& reaction = problem.reactions[r];
        if (reaction.stoichiometry.size() != n || !all_numeric(reaction.stoichiometry))
            throw std::invalid_argument("reaction activities: stoichiometry");
        if (reaction.reference) {
            if (!model || !problem.source_reference) throw std::invalid_argument("reaction activities: source declaration unavailable");
            rows.push_back(reaction.stoichiometry); references.push_back(*reaction.reference); referenced.push_back(r);
        }
    }
    std::vector<std::vector<NeutralReferenceReactionDiagnostic::Derivative>> species_references;
    const auto diagnostics = rows.empty() ? std::vector<NeutralReferenceReactionDiagnostic>{}
        : reaction_references(*model, *result.temperature, *result.pressure, *problem.source_reference,
            rows, references, false, std::nullopt, 0, false, 0, &species_references);
    std::vector<ReactionActivities> output;
    for (std::size_t p = 0; p < result.phase_count; ++p) {
        std::vector<double> x(result.mole_fractions.begin() + p * n, result.mole_fractions.begin() + (p + 1) * n);
        std::optional<State> state;
        if (ideal && ideal->solvent_index && x[*ideal->solvent_index] * result.amounts[p] <= 0.0)
            throw std::invalid_argument("unit-activity solvent absent from physical state");
        if (model) state = model->state(*result.temperature, result.molar_densities[p], x);
        const auto density = model ? species_activities(x, &*state) : SpeciesActivities{};
        for (std::size_t r = 0; r < problem.reactions.size(); ++r) {
            ReactionActivities row; row.phase_index = p; row.reaction_index = r;
            row.basis = model ? EOS_STANDARD_STATE_ID : ideal->standard_state_id();
            const auto found = std::find(referenced.begin(), referenced.end(), r);
            if (found != referenced.end()) {
                const auto k = static_cast<std::size_t>(found - referenced.begin());
                const auto& declaration = *problem.source_reference;
                row.reference_pressure_pa = references[k].reference_pressure_pa.value_or(*result.pressure);
                row.basis = references[k].source_basis == ReferenceSourceBasis::CommonMolalityInfiniteDilution
                    ? "common-molality-infinite-dilution; EOS single-ion convention" : "raw-mole-fraction-species-reference; EOS single-ion convention";
                if (diagnostics[k].status != NeutralReferenceStatus::Available) row.cause = diagnostics[k].message;
                row.species = species_activities(x, &*state, species_references[k], &model->parameters(), &declaration, references[k].source_basis);
            } else row.species = model ? density : species_activities(x, nullptr);
            if (ideal && ideal->activity_basis == "molarity") {
                row.species.solvent_index = ideal->solvent_index;
                for (std::size_t i = 0; i < n; ++i) {
                    const double value = ideal->molarity_log_activity(std::log(result.amounts[p]) + std::log(x[i]), i);
                    row.species.log_activities[i] = value;
                    const double activity = std::exp(value);
                    row.species.activities[i].reset();
                    if (!std::isfinite(activity) || (activity == 0.0 && std::isfinite(value)))
                        row.species.causes[i] = "activity exponentiation overflow or underflow";
                    else row.species.activities[i] = activity;
                }
            }
            row.log_quotient = log_quotient(problem.reactions[r].stoichiometry, row.species, row.cause);
            if (model) pressure_basis(*model, result, p, x, problem.reactions[r].stoichiometry, density, row);
            output.push_back(std::move(row));
        }
    }
    return output;
}
} // namespace
std::vector<ReactionActivities> reaction_activities(const HelmholtzModel& model, const EquilibriumResult& result) {
    return read_reaction_activities(&model, nullptr, result);
}
std::vector<ReactionActivities> reaction_activities(const IdealMoleFraction& model, const EquilibriumResult& result) {
    return read_reaction_activities(nullptr, &model, result);
}
// EqID solid_standard_gibbs: a present pure solid (positive amount) has unit activity on its pressure-corrected standard state.
std::vector<std::optional<double>> solid_activities(const EquilibriumResult& result) {
    core::Scope scope("solid_activities");
    std::vector<std::optional<double>> activities;
    for (const double amount : result.solid_amounts) activities.push_back(amount > 0.0 ? std::optional<double>(1.0) : std::nullopt);
    return activities;
}
