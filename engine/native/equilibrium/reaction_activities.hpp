// Included once, by gibbs_minimizer.cpp in namespace epcsaft::equilibrium: species activities and ln Q at a returned state.
namespace {
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
            } else row.species = species_activities(x, state ? &*state : nullptr);
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
            row.log_quotient = 0.0;
            for (std::size_t i = 0; row.cause.empty() && i < n; ++i) {
                const double nu = problem.reactions[r].stoichiometry[i];
                if (nu == 0) continue;
                if (!row.species.log_activities[i]) row.cause = row.species.causes[i];
                else *row.log_quotient += nu * *row.species.log_activities[i];
            }
            if (std::isnan(*row.log_quotient)) row.cause = "conflicting infinite log activity terms";
            if (!row.cause.empty()) row.log_quotient.reset();
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
