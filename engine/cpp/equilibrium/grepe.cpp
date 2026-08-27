#include "grepe.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace epcsaft_engine::equilibrium {
namespace {

constexpr double kGasConstant = 8.31446261815324;

CompiledReactionSystem expand_species_face(
    CompiledReactionSystem compiled,
    const ReactionSystemInput& input
) {
    if (compiled.species_count == input.species_ids.size()) {
        return compiled;
    }
    const std::size_t full_count = input.species_ids.size();
    const std::size_t reduced_count = compiled.species_count;
    const auto expand_matrix = [&](const DenseMatrix& reduced) {
        DenseMatrix full{
            reduced.rows,
            full_count,
            std::vector<double>(reduced.rows * full_count, 0.0),
        };
        for (std::size_t row = 0; row < reduced.rows; ++row) {
            for (std::size_t retained = 0; retained < reduced_count; ++retained) {
                full(row, compiled.retained_species_indices[retained]) =
                    reduced(row, retained);
            }
        }
        return full;
    };
    const auto expand_species = [&](const std::vector<double>& reduced, std::size_t width) {
        std::vector<double> full(full_count * width, 0.0);
        for (std::size_t retained = 0; retained < reduced_count; ++retained) {
            const std::size_t species = compiled.retained_species_indices[retained];
            std::copy_n(
                reduced.begin() + static_cast<std::ptrdiff_t>(retained * width),
                width,
                full.begin() + static_cast<std::ptrdiff_t>(species * width)
            );
        }
        return full;
    };
    compiled.balance_matrix = expand_matrix(compiled.balance_matrix);
    compiled.reaction_matrix = expand_matrix(compiled.reaction_matrix);
    compiled.g_ref = expand_species(compiled.g_ref, 1);
    if (!compiled.g_ref_pressure_derivatives_per_pa.empty()) {
        compiled.g_ref_pressure_derivatives_per_pa = expand_species(
            compiled.g_ref_pressure_derivatives_per_pa,
            1
        );
    }
    if (!compiled.g_ref_parameter_derivatives.empty()) {
        compiled.g_ref_parameter_derivatives = expand_species(
            compiled.g_ref_parameter_derivatives,
            compiled.g_ref_parameter_derivatives.size() / reduced_count
        );
    }
    compiled.species_count = full_count;
    compiled.charges = input.charges;
    compiled.molar_masses_kg_per_mol = input.molar_masses_kg_per_mol;
    compiled.feed_amounts = input.feed_amounts;
    return compiled;
}

SourceReferenceTransferEvidence source_evidence(
    const SourceReferenceResult& source,
    const GrepeSourceStandardState& standard,
    const ProviderActiveParameterSet& parameters
) {
    if (source.status != SOURCE_REFERENCE_OK) {
        throw std::invalid_argument(
            source.error.empty() ? "source-reference evaluation failed" : source.error
        );
    }
    SourceReferenceTransferEvidence result;
    result.status = "ok";
    result.domain_status = "admitted-domain";
    result.convergence_status = "converged-limit";
    result.reference_state_id = standard.reference.reference_state_id;
    result.activity_convention_id = standard.reference.activity_convention_id;
    result.component_ids = standard.reference.component_ids;
    result.neutral_basis_row_count = source.neutral_basis_row_count;
    result.neutral_basis = source.neutral_basis;
    result.log_fugacity_contractions = source.log_fugacity_contractions;
    result.activity_scale_log_contractions = source.activity_scale_log_contractions;
    result.transfer_log_contractions = source.transfer_log_contractions;
    result.source_composition = source.source_composition;
    if ((source.derivative_availability & 2u) != 0u) {
        result.derivative_availability.push_back("pressure");
    }
    if ((source.derivative_availability & 8u) != 0u) {
        result.derivative_availability.push_back("parameters");
    }
    result.pressure_derivatives_per_pa = source.pressure_derivatives_per_pa;
    result.active_parameter_identities = parameters.identities;
    result.parameter_derivatives = source.parameter_derivatives;
    result.temperature_k = source.temperature_k;
    result.pressure_pa = source.pressure_pa;
    result.source_reference_pressure_pa = source.source_reference_pressure_pa;
    result.standard_molality_mol_per_kg = source.standard_molality_mol_per_kg;
    result.reference_convergence_error = source.reference_convergence_error;
    result.source_temperature_min_k = source.source_temperature_min_k;
    result.source_temperature_max_k = source.source_temperature_max_k;
    result.source_pressure_min_pa = source.source_pressure_min_pa;
    result.source_pressure_max_pa = source.source_pressure_max_pa;
    result.parameter_fingerprint = source.parameter_fingerprint;
    result.topology_fingerprint = source.topology_fingerprint;
    result.component_order_fingerprint = source.component_order_fingerprint;
    result.reference_state_fingerprint = source.reference_state_fingerprint;
    result.domain_fingerprint = source.domain_fingerprint;
    result.artifact_fingerprint = source.artifact_fingerprint;
    result.helmholtz_basis_id = source.helmholtz_basis_id;
    return result;
}

ReactionSystemInput effective_reaction(
    const GrepeRequest& request,
    const std::optional<double>& coordinate = std::nullopt
) {
    ReactionSystemInput input = request.reaction;
    if (!request.feed.has_value()) {
        return input;
    }
    const GrepeFeed& feed = *request.feed;
    if (feed.affine) {
        if (!coordinate.has_value() || feed.coordinates.size() != 1
            || feed.component_coordinate_matrix.rows != input.species_ids.size()
            || feed.component_coordinate_matrix.columns != 1) {
            throw std::invalid_argument("affine GREPE feed requires one resolved coordinate");
        }
        input.feed_amounts.resize(feed.amounts_or_base_mol.size());
        for (std::size_t component = 0; component < input.feed_amounts.size(); ++component) {
            input.feed_amounts[component] = feed.amounts_or_base_mol[component]
                + feed.component_coordinate_matrix(component, 0) * *coordinate;
        }
    } else {
        input.feed_amounts = feed.amounts_or_base_mol;
    }
    input.conserved_totals.assign(input.balance_matrix.rows, 0.0);
    for (std::size_t row = 0; row < input.balance_matrix.rows; ++row) {
        for (std::size_t component = 0; component < input.balance_matrix.columns; ++component) {
            input.conserved_totals[row] +=
                input.balance_matrix(row, component) * input.feed_amounts[component];
        }
    }
    return input;
}

struct CompiledGrepeFinite {
    CompiledReactionSystem system;
    double source_residual = 0.0;
};

CompiledGrepeFinite compile_finite(
    const ProviderContext& provider,
    ReactionSystemInput input,
    const GrepeRequest& request,
    bool ideal_gas
) {
    input.temperature_k = request.temperature.value;
    input.pressure_pa = request.pressure.value;
    input.provider_fingerprint = provider.fingerprint();
    const ProviderActiveParameterSet parameters = provider.active_parameters();
    const bool configured_directions = request.directions.has_value();
    CompiledGrepeFinite result;
    if (!request.source.has_value()) {
        for (EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
            record.temperature_k = input.temperature_k;
            record.pressure_pa = input.pressure_pa;
        }
        result.system = compile_reaction_system(input, !configured_directions);
    } else {
        const GrepeSourceStandardState& standard = *request.source;
        for (const EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
            if (record.reference_id != standard.id
                || record.conversion_id
                    != "source-standard-state-to-provider-neutral-reference"
                || record.pressure_pa != standard.reference_pressure_pa) {
                throw std::invalid_argument(
                    "source equilibrium-constant declaration is incompatible"
                );
            }
        }
        const ProviderContext source_provider(provider.model());
        const SourceReferenceResult source = source_provider.evaluate_source_reference_transfer(
            input.temperature_k,
            input.pressure_pa,
            standard.reference
        );
        const SourceReferenceTransferEvidence evidence = source_evidence(
            source,
            standard,
            parameters
        );
        const SourceStandardStateResult transformed = transform_source_standard_state(
            input.reaction_matrix,
            input.ln_k,
            standard.log_activity_scale_factors,
            input.species_ids,
            input.charges,
            input.provider_fingerprint,
            input.temperature_k,
            input.pressure_pa,
            evidence
        );
        result.source_residual = transformed.representation_residual_inf_norm;
        input.ln_k = transformed.ln_k_provider_basis;
        for (EquilibriumConstantRecord& record : input.equilibrium_constant_records) {
            record.reference_id = "provider-helmholtz-coordinate-basis";
            record.conversion_id = "already-provider-basis";
            record.pressure_pa = input.pressure_pa;
        }
        result.system = compile_reaction_system(input, !configured_directions);
        if (!transformed.ln_k_pressure_derivatives_per_pa.empty()) {
            ReactionSystemInput derivative = input;
            derivative.ln_k = transformed.ln_k_pressure_derivatives_per_pa;
            const CompiledReactionSystem compiled = compile_reaction_system(
                derivative,
                !configured_directions
            );
            if (compiled.retained_species_indices != result.system.retained_species_indices) {
                throw std::invalid_argument(
                    "source-reference pressure derivative species face changed"
                );
            }
            result.system.g_ref_pressure_derivatives_per_pa = compiled.g_ref;
        }
        const std::size_t parameter_count = parameters.identities.size();
        if (parameter_count != 0) {
            if (transformed.ln_k_parameter_derivatives.size()
                != input.reaction_matrix.rows * parameter_count) {
                throw std::invalid_argument(
                    "source-reference parameter derivatives are unavailable"
                );
            }
            result.system.g_ref_parameter_derivatives.assign(
                result.system.species_count * parameter_count,
                0.0
            );
            for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
                ReactionSystemInput derivative = input;
                for (std::size_t reaction = 0; reaction < input.reaction_matrix.rows; ++reaction) {
                    derivative.ln_k[reaction] = transformed.ln_k_parameter_derivatives[
                        reaction * parameter_count + parameter
                    ];
                }
                const CompiledReactionSystem compiled = compile_reaction_system(
                    derivative,
                    !configured_directions
                );
                if (compiled.retained_species_indices
                    != result.system.retained_species_indices) {
                    throw std::invalid_argument(
                        "source-reference derivative species face changed"
                    );
                }
                for (std::size_t species = 0; species < result.system.species_count; ++species) {
                    result.system.g_ref_parameter_derivatives[
                        species * parameter_count + parameter
                    ] = compiled.g_ref[species];
                }
            }
        }
    }
    result.system = expand_species_face(std::move(result.system), input);
    result.system.ideal_gas = ideal_gas;
    return result;
}

std::optional<PhaseDirectionCompilation> compile_directions(
    const GrepeRequest& request,
    const std::vector<ReactiveFinitePhaseStart>& starts
) {
    if (!request.directions.has_value()) {
        return std::nullopt;
    }
    const GrepeDirectionSystem& declaration = *request.directions;
    const std::size_t component_count = request.reaction.species_ids.size();
    DenseMatrix references{
        starts.size(),
        component_count,
        std::vector<double>(starts.size() * component_count),
    };
    for (std::size_t phase = 0; phase < starts.size(); ++phase) {
        for (std::size_t component = 0; component < component_count; ++component) {
            references(phase, component) = starts[phase].amounts_mol[component];
        }
    }
    PhaseDirectionCompilation compiled = compile_phase_direction_system(
        declaration.phase_ids,
        declaration.species_ids,
        declaration.charges,
        declaration.balance_matrix,
        references,
        declaration.direction_matrix,
        declaration.driving_constants
    );
    for (std::size_t direction = 0;
         direction < compiled.direction_matrix.rows;
         ++direction) {
        double change = 0.0;
        double scale = 0.0;
        for (std::size_t phase = 0; phase < starts.size(); ++phase) {
            for (std::size_t species = 0; species < component_count; ++species) {
                const double contribution = compiled.direction_matrix(
                    direction, phase * component_count + species
                ) * request.reaction.molar_masses_kg_per_mol[species];
                change += contribution;
                scale += std::abs(contribution);
            }
        }
        if (std::abs(change)
            > kPhaseDirectionConservationTolerance * std::max(1.0, scale)) {
            throw std::invalid_argument(
                "phase direction violates global molar-mass conservation"
            );
        }
    }
    return compiled;
}

std::vector<std::vector<std::vector<double>>> generated_allocations(
    const GrepeRequest& request,
    const ReactionSystemInput& reaction,
    const std::vector<GrepePhase>& phases
) {
    const std::size_t phase_count = phases.size();
    const std::size_t component_count = reaction.feed_amounts.size();
    if (phase_count == 1) {
        return {{{std::vector<double>(component_count, 1.0)}}};
    }
    if (std::any_of(
            reaction.feed_amounts.begin(),
            reaction.feed_amounts.end(),
            [&](double amount) {
                return amount <= phase_count * 0.1 * request.trace_floor;
            })) {
        throw std::invalid_argument(
            "automatic finite-phase starts require a strict-interior species feed"
        );
    }
    std::vector<std::vector<std::vector<double>>> candidates;
    if (phase_count == 2) {
        for (std::size_t reaction_index = 0;
             reaction_index < reaction.reaction_matrix.rows;
             ++reaction_index) {
            std::vector<double> first(component_count, 0.5);
            bool high = false;
            bool low = false;
            for (std::size_t component = 0; component < component_count; ++component) {
                const double coefficient = reaction.reaction_matrix(reaction_index, component);
                if (reaction.charges[component] == 0 && coefficient != 0.0) {
                    first[component] = coefficient > 0.0 ? 0.8 : 0.2;
                    high = high || coefficient > 0.0;
                    low = low || coefficient < 0.0;
                }
            }
            if (high && low) {
                std::vector<double> second(component_count);
                std::transform(
                    first.begin(), first.end(), second.begin(), [](double value) {
                        return 1.0 - value;
                    }
                );
                candidates.push_back({std::move(first), std::move(second)});
            }
        }
        for (std::size_t anchor = 0; anchor < component_count; ++anchor) {
            if (reaction.charges[anchor] != 0) {
                continue;
            }
            std::vector<double> first(component_count, 0.2);
            first[anchor] = 0.8;
            for (std::size_t component = 0; component < component_count; ++component) {
                if (reaction.charges[component] != 0) {
                    first[component] = 0.5;
                }
            }
            std::vector<double> second(component_count);
            std::transform(first.begin(), first.end(), second.begin(), [](double value) {
                return 1.0 - value;
            });
            candidates.push_back({std::move(first), std::move(second)});
        }
    }
    candidates.push_back(std::vector<std::vector<double>>(
        phase_count,
        std::vector<double>(component_count, 1.0 / static_cast<double>(phase_count))
    ));
    std::vector<std::vector<std::vector<double>>> projected;
    for (const auto& candidate : candidates) {
        std::vector<std::vector<double>> rows(
            phase_count,
            std::vector<double>(component_count, 0.0)
        );
        for (std::size_t component = 0; component < component_count; ++component) {
            double scale = 0.0;
            for (std::size_t phase = 0; phase < phase_count; ++phase) {
                if (std::find(
                        phases[phase].support_indices.begin(),
                        phases[phase].support_indices.end(),
                        component
                    ) != phases[phase].support_indices.end()) {
                    scale += candidate[phase][component];
                }
            }
            if (!(scale > 0.0)) {
                throw std::invalid_argument(
                    "every fed species must belong to a declared finite-phase support"
                );
            }
            for (std::size_t phase = 0; phase < phase_count; ++phase) {
                if (std::find(
                        phases[phase].support_indices.begin(),
                        phases[phase].support_indices.end(),
                        component
                    ) != phases[phase].support_indices.end()) {
                    rows[phase][component] = candidate[phase][component] / scale;
                }
            }
        }
        if (std::find(projected.begin(), projected.end(), rows) == projected.end()) {
            projected.push_back(std::move(rows));
        }
    }
    return projected;
}

ReactiveFinitePhaseStart start_from_allocation(
    const epcsaft::NativeModel& model,
    const GrepePhase& phase,
    const std::vector<double>& feed,
    const std::vector<double>& allocation,
    double temperature_k,
    double pressure_pa
) {
    ReactiveFinitePhaseStart start;
    start.amounts_mol.resize(feed.size());
    for (std::size_t component = 0; component < feed.size(); ++component) {
        start.amounts_mol[component] = feed[component] * allocation[component];
    }
    const double total = std::accumulate(
        start.amounts_mol.begin(), start.amounts_mol.end(), 0.0
    );
    if (!(total > 0.0)) {
        throw std::invalid_argument("finite phase start has zero total amount");
    }
    std::vector<double> composition(start.amounts_mol.size());
    std::transform(
        start.amounts_mol.begin(),
        start.amounts_mol.end(),
        composition.begin(),
        [total](double amount) { return amount / total; }
    );
    const auto state = model.evaluate_pressure(
        temperature_k,
        pressure_pa,
        composition,
        phase.role == PhaseRole::Liquid ? "liquid" : "vapor"
    );
    start.molar_volume_m3_per_mol = 1.0 / state.molar_density;
    return start;
}

std::optional<ReactiveFinitePhaseStart> continuation_start(
    const GrepeRequest& request,
    const ReactionSystemInput& reaction,
    const std::vector<GrepePhase>& phases
) {
    if (phases.size() != 1 || request.continuation_phases.size() != 1) {
        return std::nullopt;
    }
    const GrepeContinuationPhase& continuation = request.continuation_phases[0];
    if (continuation.mole_fractions.size() != reaction.balance_matrix.columns) {
        return std::nullopt;
    }
    std::vector<double> projected(reaction.balance_matrix.rows, 0.0);
    for (std::size_t row = 0; row < reaction.balance_matrix.rows; ++row) {
        for (std::size_t component = 0; component < reaction.balance_matrix.columns;
             ++component) {
            projected[row] += reaction.balance_matrix(row, component)
                * continuation.mole_fractions[component];
        }
    }
    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t row = 0; row < projected.size(); ++row) {
        numerator += projected[row] * reaction.conserved_totals[row];
        denominator += projected[row] * projected[row];
    }
    const double total = numerator / denominator;
    if (!std::isfinite(total) || total <= 0.0) {
        return std::nullopt;
    }
    ReactiveFinitePhaseStart start;
    start.amounts_mol.resize(continuation.mole_fractions.size());
    std::transform(
        continuation.mole_fractions.begin(),
        continuation.mole_fractions.end(),
        start.amounts_mol.begin(),
        [total](double fraction) { return total * fraction; }
    );
    start.molar_volume_m3_per_mol = continuation.molar_volume_m3_per_mol;
    return start;
}

std::vector<std::vector<ReactiveFinitePhaseStart>> start_sets(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    const ReactionSystemInput& reaction,
    const std::vector<GrepePhase>& phases
) {
    std::vector<std::vector<ReactiveFinitePhaseStart>> result;
    if (std::all_of(phases.begin(), phases.end(), [](const GrepePhase& phase) {
            return phase.start.has_value();
        })) {
        std::vector<ReactiveFinitePhaseStart> declared;
        for (const GrepePhase& phase : phases) {
            declared.push_back(*phase.start);
        }
        result.push_back(std::move(declared));
    }
    if (request.declared_starts_only) {
        return result;
    }
    if (const auto continuation = continuation_start(request, reaction, phases)) {
        result.push_back({*continuation});
    }
    const std::optional<std::vector<ReactiveFinitePhaseStart>> anchor =
        result.empty()
        ? std::nullopt
        : std::optional<std::vector<ReactiveFinitePhaseStart>>(result.front());
    for (const auto& allocation : generated_allocations(request, reaction, phases)) {
        std::vector<std::size_t> order(phases.size());
        std::iota(order.begin(), order.end(), 0);
        std::vector<std::vector<std::size_t>> orders;
        do {
            bool supported = true;
            for (std::size_t phase = 0; phase < phases.size() && supported; ++phase) {
                for (std::size_t component = 0; component < reaction.feed_amounts.size();
                     ++component) {
                    if (allocation[order[phase]][component] != 0.0
                        && std::find(
                            phases[phase].support_indices.begin(),
                            phases[phase].support_indices.end(),
                            component
                        ) == phases[phase].support_indices.end()) {
                        supported = false;
                        break;
                    }
                }
            }
            if (!supported) {
                continue;
            }
            orders.push_back(order);
        } while (std::next_permutation(order.begin(), order.end()));
        if (anchor.has_value()) {
            std::stable_sort(
                orders.begin(),
                orders.end(),
                [&](const auto& left, const auto& right) {
                    const auto score = [&](const auto& candidate) {
                        double value = 0.0;
                        for (std::size_t phase = 0; phase < phases.size(); ++phase) {
                            for (std::size_t component = 0;
                                 component < reaction.feed_amounts.size();
                                 ++component) {
                                if (reaction.feed_amounts[component] > 0.0) {
                                    const double difference = allocation[candidate[phase]][component]
                                        - (*anchor)[phase].amounts_mol[component]
                                            / reaction.feed_amounts[component];
                                    value += difference * difference;
                                }
                            }
                        }
                        return value;
                    };
                    return score(left) < score(right);
                }
            );
        } else if (!orders.empty()) {
            std::vector<std::size_t> canonical(phases.size());
            std::vector<std::size_t> sorted(phases.size());
            std::iota(sorted.begin(), sorted.end(), 0);
            std::stable_sort(sorted.begin(), sorted.end(), [&](std::size_t left, std::size_t right) {
                return phases[left].identity < phases[right].identity;
            });
            for (std::size_t phase = 0; phase < phases.size(); ++phase) {
                canonical[phase] = static_cast<std::size_t>(std::find(
                    sorted.begin(), sorted.end(), phase
                ) - sorted.begin());
            }
            orders = {std::move(canonical)};
        }
        for (const auto& selected_order : orders) {
            std::vector<ReactiveFinitePhaseStart> starts;
            try {
                for (std::size_t phase = 0; phase < phases.size(); ++phase) {
                    starts.push_back(start_from_allocation(
                        model,
                        phases[phase],
                        reaction.feed_amounts,
                        allocation[selected_order[phase]],
                        request.temperature.value,
                        request.pressure.value
                    ));
                }
            } catch (const std::exception&) {
                starts.clear();
            }
            if (!starts.empty()) {
                result.push_back(std::move(starts));
            }
        }
    }
    return result;
}

double failure_score(const ReactiveFinitePhaseResult& result) {
    const double values[] = {
        result.balance_inf_norm,
        result.pressure_relative_inf_norm,
        result.reaction_affinity_inf_norm,
        result.transfer_chemical_potential_inf_norm,
        result.kkt_stationarity_inf_norm,
    };
    double score = 0.0;
    for (double value : values) {
        score += std::isfinite(value) ? std::abs(value) : 1.0e100;
    }
    return score;
}

GrepeResult unsupported(const std::string& reason) {
    GrepeResult result;
    result.formulation = "unsupported";
    result.certificate_scope = "unavailable_formulation";
    result.failure_code = "unsupported_topology";
    result.equilibrium.solver_status = "not_started";
    result.equilibrium.numerical_status = "unavailable";
    result.equilibrium.physical_status = "not_adjudicated";
    result.equilibrium.provider_domain_status = "not_adjudicated";
    result.equilibrium.failure_reason = reason;
    return result;
}

struct FiniteGrepeContext final {
    std::vector<GrepePhase> phases;
    bool ideal;
    ProviderContext provider;
    ReactionSystemInput reaction;
    CompiledGrepeFinite compiled;
    ProviderActiveParameterSet parameters;
    std::vector<std::array<double, 2>> packing;
    std::vector<std::vector<std::size_t>> supports;
    GrepeResult result;
};

std::optional<FiniteGrepeContext> prepare_finite_grepe(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request
) {
    std::vector<GrepePhase> phases = [&] {
        std::vector<GrepePhase> selected;
        for (const GrepePhase& phase : request.phases) {
            if (phase.finite) {
                selected.push_back(phase);
            }
        }
        return selected;
    }();
    const bool ideal = phases.size() == 1 && phases[0].ideal_gas;
    if (phases.empty() || (!ideal && std::any_of(
            phases.begin(), phases.end(), [](const GrepePhase& phase) {
                return phase.ideal_gas;
            }))) {
        return std::nullopt;
    }
    ProviderContext provider(model);
    ReactionSystemInput reaction = effective_reaction(request);
    CompiledGrepeFinite compiled = compile_finite(provider, reaction, request, ideal);
    ProviderActiveParameterSet parameters = provider.active_parameters();
    std::vector<std::array<double, 2>> packing;
    std::vector<std::vector<std::size_t>> supports;
    for (const GrepePhase& phase : phases) {
        packing.push_back(ideal ? std::array<double, 2>{0.0, 0.0}
                                : phase.packing_fraction_interval);
        supports.push_back(phase.support_indices);
        supports.back().erase(
            std::remove_if(
                supports.back().begin(),
                supports.back().end(),
                [&compiled](std::size_t species) {
                    return std::binary_search(
                        compiled.system.removed_species_indices.begin(),
                        compiled.system.removed_species_indices.end(),
                        species
                    );
                }
            ),
            supports.back().end()
        );
    }
    GrepeResult result;
    result.formulation = "finite";
    result.certificate_scope = phases.size() == 1
        ? "fixed_TP_single_finite_reacting_phase"
        : "fixed_TP_configured_finite_phases_local_minimum";
    result.temperature_k = request.temperature.value;
    result.pressure_pa = request.pressure.value;
    result.source_reference_representation_residual_inf_norm = compiled.source_residual;
    return FiniteGrepeContext{
        std::move(phases),
        ideal,
        std::move(provider),
        std::move(reaction),
        std::move(compiled),
        std::move(parameters),
        std::move(packing),
        std::move(supports),
        std::move(result),
    };
}

GrepeResult solve_finite_result(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request
) {
    std::optional<FiniteGrepeContext> prepared = prepare_finite_grepe(model, request);
    if (!prepared.has_value()) {
        return unsupported(
            "declared reactive phase/state topology has no implemented GREPE formulation"
        );
    }
    FiniteGrepeContext& context = *prepared;
    const std::vector<GrepePhase>& phases = context.phases;
    const bool ideal = context.ideal;
    const ReactionSystemInput& reaction = context.reaction;
    const CompiledGrepeFinite& compiled = context.compiled;
    const ProviderActiveParameterSet& parameters = context.parameters;
    const std::vector<std::array<double, 2>>& packing = context.packing;
    const std::vector<std::vector<std::size_t>>& supports = context.supports;
    GrepeResult& result = context.result;
    std::vector<std::vector<ReactiveFinitePhaseStart>> candidates;
    if (ideal) {
        ReactiveFinitePhaseStart start;
        if (phases[0].start.has_value()) {
            start = *phases[0].start;
        } else {
            start.amounts_mol = reaction.feed_amounts;
            start.molar_volume_m3_per_mol =
                kGasConstant * request.temperature.value / request.pressure.value;
        }
        candidates = {{std::move(start)}};
    } else {
        candidates = start_sets(model, request, reaction, phases);
    }
    if (candidates.empty()) {
        throw std::invalid_argument("Provider could not construct a finite-phase candidate start");
    }
    bool selected = false;
    std::size_t total_iterations = 0;
    std::size_t total_evaluations = 0;
    double selected_score = 0.0;
    for (const auto& starts : candidates) {
        ProviderContext attempt_provider(model);
        const std::optional<PhaseDirectionCompilation> directions =
            compile_directions(request, starts);
        ReactiveFinitePhaseResult attempt = solve_provider_reactive_finite_phases(
            compiled.system,
            attempt_provider,
            request.temperature.value,
            request.pressure.value,
            packing,
            starts,
            request.trace_floor,
            parameters.identities.empty() ? nullptr : &parameters,
            500,
            request.perform_derivative_audits,
            directions.has_value() ? &*directions : nullptr,
            supports,
            request.outputs
        );
        total_iterations += attempt.optimizer_iterations;
        total_evaluations += attempt.provider_evaluation_count;
        const double score = failure_score(attempt);
        if (!selected || attempt.accepted || score < selected_score) {
            selected = true;
            selected_score = score;
            result.equilibrium = std::move(attempt);
        }
        if (result.equilibrium.accepted) {
            break;
        }
    }
    result.equilibrium.optimizer_iterations = total_iterations;
    result.equilibrium.provider_evaluation_count = total_evaluations;
    result.branch_kind = "finite_phases_local_minimum";
    result.phase_mechanical_classes.assign(
        phases.size(), ideal ? "analytical_ideal_gas" : "strict_stable_local_minimum"
    );
    return result;
}

GrepeResult evaluate_finite_result(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    const std::vector<double>& local_variables,
    const std::vector<double>& local_lower_bounds,
    const std::vector<double>& local_upper_bounds
) {
    std::optional<FiniteGrepeContext> prepared = prepare_finite_grepe(model, request);
    if (!prepared.has_value()) {
        return unsupported(
            "declared reactive phase/state topology has no implemented GREPE formulation"
        );
    }
    FiniteGrepeContext& context = *prepared;
    context.result.equilibrium = evaluate_provider_reactive_finite_phases_local(
        context.compiled.system,
        context.provider,
        request.temperature.value,
        request.pressure.value,
        context.packing,
        local_variables,
        local_lower_bounds,
        local_upper_bounds,
        context.parameters.identities.empty() ? nullptr : &context.parameters,
        nullptr,
        context.supports,
        request.outputs
    );
    return std::move(context.result);
}

GrepeResult finalize_finite_result(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    const std::vector<double>& local_variables,
    const std::vector<double>& local_lower_bounds,
    const std::vector<double>& local_upper_bounds
) {
    std::optional<FiniteGrepeContext> prepared = prepare_finite_grepe(model, request);
    if (!prepared.has_value()) {
        return unsupported(
            "declared reactive phase/state topology has no implemented GREPE formulation"
        );
    }
    FiniteGrepeContext& context = *prepared;
    context.result.equilibrium = finalize_provider_reactive_finite_phases_local(
        context.compiled.system,
        context.provider,
        request.temperature.value,
        request.pressure.value,
        context.packing,
        local_variables,
        local_lower_bounds,
        local_upper_bounds,
        context.parameters.identities.empty() ? nullptr : &context.parameters,
        nullptr,
        context.supports,
        request.outputs,
        request.trace_floor,
        request.perform_derivative_audits
    );
    if (context.result.equilibrium.accepted) {
        context.result.branch_kind = "finite_phases_local_minimum";
        context.result.phase_mechanical_classes.assign(
            context.phases.size(),
            context.ideal ? "analytical_ideal_gas" : "strict_stable_local_minimum"
        );
    }
    return std::move(context.result);
}

GrepeResult initialize_finite_result(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request
) {
    std::optional<FiniteGrepeContext> prepared = prepare_finite_grepe(model, request);
    if (!prepared.has_value()) {
        return unsupported(
            "declared reactive phase/state topology has no implemented GREPE formulation"
        );
    }
    FiniteGrepeContext& context = *prepared;
    if (request.continuation_phases.size() != 1
        || request.continuation_phases[0].identity != context.phases[0].identity) {
        throw std::invalid_argument(
            "full-space GREPE requires one matching certified continuation anchor"
        );
    }
    const auto start = continuation_start(request, context.reaction, context.phases);
    if (!start.has_value()) {
        throw std::invalid_argument(
            "full-space GREPE requires one matching certified continuation anchor"
        );
    }
    const std::vector<ReactiveFinitePhaseStart> starts{*start};
    const std::optional<PhaseDirectionCompilation> directions =
        compile_directions(request, starts);
    context.result.equilibrium = solve_provider_reactive_finite_phases(
        context.compiled.system,
        context.provider,
        request.temperature.value,
        request.pressure.value,
        context.packing,
        starts,
        request.trace_floor,
        context.parameters.identities.empty() ? nullptr : &context.parameters,
        0,
        false,
        directions.has_value() ? &*directions : nullptr,
        context.supports,
        request.outputs
    );
    return std::move(context.result);
}

class GrepeFailure final : public std::runtime_error {
public:
    GrepeFailure(std::string code, const std::string& diagnostic)
        : std::runtime_error(diagnostic), code_(std::move(code)) {}

    [[nodiscard]] const std::string& code() const noexcept { return code_; }

private:
    std::string code_;
};

std::size_t input_index(
    const ReactiveFinitePhaseResult& result,
    const std::string& identity
) {
    const auto found = std::find(
        result.input_order.begin(), result.input_order.end(), identity
    );
    if (found == result.input_order.end()) {
        throw GrepeFailure(
            "parameter_topology_mismatch",
            "reactive finite-phase sensitivity order changed"
        );
    }
    return static_cast<std::size_t>(found - result.input_order.begin());
}

std::vector<double> solve_dense(
    std::vector<double> matrix,
    std::vector<double> right_hand_side
) {
    const std::size_t dimension = right_hand_side.size();
    if (matrix.size() != dimension * dimension) {
        throw GrepeFailure(
            "owner_contract_failure", "incipient tangent matrix dimensions changed"
        );
    }
    if (dimension == 0) {
        return {};
    }
    double scale = 1.0;
    for (double value : matrix) {
        scale = std::max(scale, std::abs(value));
    }
    double minimum_pivot = std::numeric_limits<double>::infinity();
    double maximum_pivot = 0.0;
    for (std::size_t column = 0; column < dimension; ++column) {
        std::size_t pivot = column;
        for (std::size_t row = column + 1; row < dimension; ++row) {
            if (std::abs(matrix[row * dimension + column])
                > std::abs(matrix[pivot * dimension + column])) {
                pivot = row;
            }
        }
        const double value = std::abs(matrix[pivot * dimension + column]);
        if (!std::isfinite(value) || value <= 1.0e-13 * scale) {
            throw GrepeFailure(
                "singular_sensitivity", "incipient-vapor composition Jacobian is singular"
            );
        }
        minimum_pivot = std::min(minimum_pivot, value);
        maximum_pivot = std::max(maximum_pivot, value);
        if (pivot != column) {
            for (std::size_t entry = column; entry < dimension; ++entry) {
                std::swap(
                    matrix[column * dimension + entry],
                    matrix[pivot * dimension + entry]
                );
            }
            std::swap(right_hand_side[column], right_hand_side[pivot]);
        }
        for (std::size_t row = column + 1; row < dimension; ++row) {
            const double multiplier = matrix[row * dimension + column]
                / matrix[column * dimension + column];
            for (std::size_t entry = column; entry < dimension; ++entry) {
                matrix[row * dimension + entry] -=
                    multiplier * matrix[column * dimension + entry];
            }
            right_hand_side[row] -= multiplier * right_hand_side[column];
        }
    }
    if (maximum_pivot / minimum_pivot > 1.0e10) {
        throw GrepeFailure(
            "ill_conditioned_sensitivity",
            "incipient-vapor composition Jacobian is ill-conditioned"
        );
    }
    std::vector<double> solution(dimension, 0.0);
    for (std::size_t row = dimension; row-- > 0;) {
        double value = right_hand_side[row];
        for (std::size_t column = row + 1; column < dimension; ++column) {
            value -= matrix[row * dimension + column] * solution[column];
        }
        solution[row] = value / matrix[row * dimension + row];
    }
    if (std::any_of(solution.begin(), solution.end(), [](double value) {
            return !std::isfinite(value);
        })) {
        throw GrepeFailure(
            "singular_sensitivity", "incipient-vapor linear solve returned nonfinite values"
        );
    }
    return solution;
}

std::vector<double> softmax_reference(const std::vector<double>& log_ratios) {
    const double maximum = std::max(
        0.0,
        log_ratios.empty() ? 0.0 : *std::max_element(log_ratios.begin(), log_ratios.end())
    );
    std::vector<double> fractions(log_ratios.size() + 1);
    double total = std::exp(-maximum);
    fractions.back() = total;
    for (std::size_t index = 0; index < log_ratios.size(); ++index) {
        fractions[index] = std::exp(log_ratios[index] - maximum);
        total += fractions[index];
    }
    for (double& fraction : fractions) {
        fraction /= total;
        if (!std::isfinite(fraction) || fraction <= 0.0) {
            throw GrepeFailure(
                "active_bound_ambiguity",
                "declared vapor support left its strict positive composition chart"
            );
        }
    }
    return fractions;
}

struct LiquidEvaluation {
    GrepeResult solved;
    std::vector<int> output_positions;
    std::size_t log_fugacity_offset = 0;
};

GrepeRequest resolved_request(const GrepeRequest& request, double coordinate) {
    GrepeRequest resolved = request;
    resolved.reaction = effective_reaction(request, coordinate);
    resolved.feed.reset();
    return resolved;
}

LiquidEvaluation evaluate_liquid(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    double pressure_pa
) {
    if (std::any_of(
            request.reaction.feed_amounts.begin(),
            request.reaction.feed_amounts.end(),
            [&](double amount) { return amount <= request.trace_floor; }
        )) {
        throw GrepeFailure(
            "active_bound_ambiguity",
            "the declared strict-interior floor reaches the reactive liquid feed"
        );
    }
    const auto finite = std::find_if(
        request.phases.begin(), request.phases.end(), [](const GrepePhase& phase) {
            return phase.finite;
        }
    );
    if (finite == request.phases.end()) {
        throw GrepeFailure("unsupported_topology", "incipient solve requires one finite phase");
    }
    const std::size_t finite_index = static_cast<std::size_t>(
        finite - request.phases.begin()
    );
    GrepeRequest liquid = request;
    liquid.pressure = {false, pressure_pa, {}, {}};
    liquid.phases = {*finite};
    liquid.boundary.reset();
    liquid.directions.reset();
    liquid.outputs.clear();
    LiquidEvaluation result;
    result.output_positions.assign(request.outputs.size(), -1);
    for (std::size_t index = 0; index < request.outputs.size(); ++index) {
        const ReactiveFinitePhaseOutput& output = request.outputs[index];
        if (output.selector != "system.pressure"
            && output.selector != "phase.partial_pressure"
            && output.phase_index == finite_index) {
            ReactiveFinitePhaseOutput copied = output;
            copied.phase_index = 0;
            result.output_positions[index] = static_cast<int>(liquid.outputs.size());
            liquid.outputs.push_back(std::move(copied));
        }
    }
    result.log_fugacity_offset = liquid.outputs.size();
    for (std::size_t component = 0;
         component < request.reaction.species_ids.size();
         ++component) {
        std::vector<double> coefficient(request.reaction.species_ids.size(), 0.0);
        coefficient[component] = 1.0;
        liquid.outputs.push_back({"phase.log_fugacity", 0, std::move(coefficient), {}});
    }
    liquid.continuation_phases.erase(
        std::remove_if(
            liquid.continuation_phases.begin(),
            liquid.continuation_phases.end(),
            [&](const GrepeContinuationPhase& phase) {
                return phase.identity != finite->identity;
            }
        ),
        liquid.continuation_phases.end()
    );
    result.solved = solve_finite_result(model, liquid);
    if (!result.solved.equilibrium.accepted) {
        const ReactiveFinitePhaseResult& failed = result.solved.equilibrium;
        const std::string code = failed.provider_domain_status == "failed"
            ? "provider_domain_rejection"
            : "numerical_convergence_failure";
        throw GrepeFailure(
            code,
            failed.failure_reason.empty() ? failed.solver_status : failed.failure_reason
        );
    }
    if (std::any_of(
            result.solved.equilibrium.phase_amounts_mol[0].begin(),
            result.solved.equilibrium.phase_amounts_mol[0].end(),
            [&](double amount) { return amount <= request.trace_floor; }
        )) {
        throw GrepeFailure(
            "active_bound_ambiguity",
            "the reactive liquid contains a trace-bound species at the coupled state"
        );
    }
    return result;
}

struct VaporPoint {
    std::vector<double> fractions;
    std::vector<double> chemical_potentials;
    std::vector<double> mu_log_ratio_at_fixed_volume;
    std::vector<double> mu_log_volume;
    std::vector<double> mu_log_ratio;
    std::vector<double> pressure_log_ratio;
    std::vector<double> mu_pressure;
    std::vector<double> mu_parameters;
    double volume_m3 = 0.0;
    double pressure_pa = 0.0;
    double packing_fraction = 0.0;
    double pressure_log_volume = 0.0;
};

VaporPoint evaluate_vapor(
    const ProviderContext& provider,
    const GrepeRequest& request,
    const GrepePhase& vapor,
    const std::vector<double>& log_ratios,
    double log_volume,
    const std::vector<std::string>& parameter_identities
) {
    VaporPoint result;
    result.fractions = softmax_reference(log_ratios);
    result.volume_m3 = std::exp(log_volume);
    const auto bounds = provider.evaluate_molar_volume_bounds(
        request.temperature.value,
        result.fractions,
        vapor.packing_fraction_interval[0],
        vapor.packing_fraction_interval[1]
    );
    if (!(bounds[0] < result.volume_m3 && result.volume_m3 < bounds[1])) {
        throw GrepeFailure(
            "active_bound_ambiguity",
            "declared vapor volume left its packing-fraction chart"
        );
    }
    const ProviderActiveParameterSet parameters = provider.active_parameters();
    const ParameterizedPhaseEvaluation block = provider.evaluate_parameterized_phase(
        request.temperature.value,
        result.fractions,
        result.volume_m3,
        vapor.packing_fraction_interval[0],
        vapor.packing_fraction_interval[1],
        parameters
    );
    const std::size_t count = result.fractions.size();
    const std::size_t reference = count - 1;
    const std::size_t state_count = count + 1;
    const std::size_t provider_parameter_count = parameters.identities.size();
    const std::size_t parameter_count = parameter_identities.size();
    bool ordered_support = block.supported_component_indices.size() == count;
    for (std::size_t index = 0; ordered_support && index < count; ++index) {
        ordered_support = block.supported_component_indices[index] == index;
    }
    if (!ordered_support
        || block.phase.gradient.size() != state_count
        || block.phase.hessian.size() != state_count * state_count
        || block.chemical_potential_parameter_derivatives_over_rt.size()
            != count * provider_parameter_count
        || block.pressure_parameter_derivatives_pa.size() != provider_parameter_count) {
        throw GrepeFailure(
            "owner_contract_failure", "vapor Provider derivative layout changed"
        );
    }
    const double h_vv = block.phase.hessian[count * state_count + count];
    double h_scale = 1.0;
    for (std::size_t index = 0; index < state_count; ++index) {
        h_scale = std::max(h_scale, std::abs(block.phase.hessian[count * state_count + index]));
    }
    if (!std::isfinite(h_vv) || h_vv <= 1.0e-14 * h_scale) {
        throw GrepeFailure(
            "second_order_failure", "declared vapor root lacks strict mechanical stability"
        );
    }
    std::vector<double> chart(count * reference, 0.0);
    for (std::size_t component = 0; component < count; ++component) {
        for (std::size_t coordinate = 0; coordinate < reference; ++coordinate) {
            chart[component * reference + coordinate] = result.fractions[component]
                * ((component == coordinate ? 1.0 : 0.0) - result.fractions[coordinate]);
        }
    }
    result.chemical_potentials.assign(block.phase.gradient.begin(), block.phase.gradient.begin() + static_cast<std::ptrdiff_t>(count));
    result.mu_log_ratio_at_fixed_volume.assign(count * reference, 0.0);
    result.mu_log_volume.assign(count, 0.0);
    result.mu_log_ratio.assign(count * reference, 0.0);
    result.pressure_log_ratio.assign(reference, 0.0);
    result.mu_pressure.assign(count, 0.0);
    const double rt = kGasConstant * request.temperature.value;
    for (std::size_t component = 0; component < count; ++component) {
        result.mu_log_volume[component] =
            block.phase.hessian[component * state_count + count] * result.volume_m3;
        for (std::size_t coordinate = 0; coordinate < reference; ++coordinate) {
            for (std::size_t species = 0; species < count; ++species) {
                result.mu_log_ratio_at_fixed_volume[component * reference + coordinate] +=
                    block.phase.hessian[component * state_count + species]
                    * chart[species * reference + coordinate];
            }
        }
    }
    std::vector<double> volume_log_ratio(reference, 0.0);
    for (std::size_t coordinate = 0; coordinate < reference; ++coordinate) {
        double row = 0.0;
        for (std::size_t component = 0; component < count; ++component) {
            row += block.phase.hessian[count * state_count + component]
                * chart[component * reference + coordinate];
        }
        result.pressure_log_ratio[coordinate] = -rt * row;
        volume_log_ratio[coordinate] = -row / h_vv;
        for (std::size_t component = 0; component < count; ++component) {
            result.mu_log_ratio[component * reference + coordinate] =
                result.mu_log_ratio_at_fixed_volume[component * reference + coordinate]
                + block.phase.hessian[component * state_count + count]
                    * volume_log_ratio[coordinate];
        }
    }
    result.pressure_log_volume = -rt * h_vv * result.volume_m3;
    const double volume_pressure = -1.0 / (rt * h_vv);
    for (std::size_t component = 0; component < count; ++component) {
        result.mu_pressure[component] =
            block.phase.hessian[component * state_count + count] * volume_pressure;
    }
    result.mu_parameters.assign(count * parameter_count, 0.0);
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        const auto found = std::find(
            parameters.identities.begin(),
            parameters.identities.end(),
            parameter_identities[parameter]
        );
        if (found == parameters.identities.end()) {
            continue;
        }
        const std::size_t provider_parameter = static_cast<std::size_t>(
            found - parameters.identities.begin()
        );
        const double volume_parameter =
            block.pressure_parameter_derivatives_pa[provider_parameter] / (rt * h_vv);
        for (std::size_t component = 0; component < count; ++component) {
            result.mu_parameters[component * parameter_count + parameter] =
                block.chemical_potential_parameter_derivatives_over_rt[
                    component * provider_parameter_count + provider_parameter
                ]
                + block.phase.hessian[component * state_count + count]
                    * volume_parameter;
        }
    }
    result.pressure_pa = block.phase.pressure_pa;
    result.packing_fraction = block.packing.value;
    return result;
}

VaporPoint solve_vapor(
    const ProviderContext& provider,
    const GrepeRequest& request,
    const GrepePhase& vapor,
    double pressure_pa,
    const std::vector<double>& liquid_potentials,
    const std::vector<std::string>& parameter_identities
) {
    const std::size_t count = vapor.support_indices.size();
    const std::size_t reference = count - 1;
    std::vector<std::vector<double>> starts;
    const auto continuation = std::find_if(
        request.continuation_phases.begin(),
        request.continuation_phases.end(),
        [&](const GrepeContinuationPhase& phase) { return phase.identity == vapor.identity; }
    );
    const bool continued = continuation != request.continuation_phases.end();
    double continuation_volume = 0.0;
    if (continued) {
        if (continuation->mole_fractions.size() != count
            || std::any_of(
                continuation->mole_fractions.begin(),
                continuation->mole_fractions.end(),
                [](double value) { return value <= 1.0e-200; }
            )) {
            throw GrepeFailure("branch_loss", "continued vapor state left its positive chart");
        }
        starts.push_back(continuation->mole_fractions);
        continuation_volume = continuation->molar_volume_m3_per_mol;
    } else {
        const double maximum = *std::max_element(
            liquid_potentials.begin(), liquid_potentials.end()
        );
        std::vector<double> weights(count);
        double total = 0.0;
        for (std::size_t component = 0; component < count; ++component) {
            weights[component] = std::exp(std::max(-700.0, liquid_potentials[component] - maximum));
            total += weights[component];
        }
        for (double& value : weights) {
            value /= total;
        }
        starts.push_back(std::move(weights));
        starts.push_back(std::vector<double>(count, 1.0 / static_cast<double>(count)));
        const double minor = 1.0e-4 / static_cast<double>(std::max<std::size_t>(1, count - 1));
        for (std::size_t dominant = 0; dominant < count; ++dominant) {
            std::vector<double> start(count, minor);
            start[dominant] = 1.0 - minor * static_cast<double>(count - 1);
            starts.push_back(std::move(start));
        }
    }
    std::vector<std::pair<VaporPoint, double>> accepted;
    std::optional<GrepeFailure> last_failure;
    for (const std::vector<double>& fractions : starts) {
        std::vector<double> coordinates(reference + 1, 0.0);
        for (std::size_t index = 0; index < reference; ++index) {
            coordinates[index] = std::log(fractions[index] / fractions[reference]);
        }
        const auto bounds = provider.evaluate_molar_volume_bounds(
            request.temperature.value,
            fractions,
            vapor.packing_fraction_interval[0],
            vapor.packing_fraction_interval[1]
        );
        const double proposed = continued
            ? continuation_volume
            : kGasConstant * request.temperature.value / pressure_pa;
        const double volume = std::min(
            std::nextafter(bounds[1], bounds[0]),
            std::max(std::nextafter(bounds[0], bounds[1]), proposed)
        );
        coordinates.back() = std::log(volume);
        std::optional<VaporPoint> point;
        double final_norm = std::numeric_limits<double>::infinity();
        try {
            for (int iteration = 0; iteration < 36; ++iteration) {
                point = evaluate_vapor(
                    provider,
                    request,
                    vapor,
                    std::vector<double>(coordinates.begin(), coordinates.end() - 1),
                    coordinates.back(),
                    parameter_identities
                );
                std::vector<double> residual(count, 0.0);
                for (std::size_t row = 0; row < reference; ++row) {
                    residual[row] = point->chemical_potentials[row]
                        - point->chemical_potentials[reference]
                        - liquid_potentials[row] + liquid_potentials[reference];
                }
                if (!(point->pressure_pa > 0.0)) {
                    throw GrepeFailure(
                        "no_incipient_phase", "declared vapor pressure is not positive"
                    );
                }
                residual.back() = std::log(point->pressure_pa / pressure_pa);
                final_norm = 0.0;
                for (double value : residual) {
                    final_norm = std::max(final_norm, std::abs(value));
                }
                if (final_norm <= 2.0e-10) {
                    break;
                }
                std::vector<double> jacobian(count * count, 0.0);
                for (std::size_t row = 0; row < reference; ++row) {
                    for (std::size_t column = 0; column < reference; ++column) {
                        jacobian[row * count + column] =
                            point->mu_log_ratio_at_fixed_volume[row * reference + column]
                            - point->mu_log_ratio_at_fixed_volume[
                                reference * reference + column
                            ];
                    }
                    jacobian[row * count + reference] = point->mu_log_volume[row]
                        - point->mu_log_volume[reference];
                }
                for (std::size_t column = 0; column < reference; ++column) {
                    jacobian[reference * count + column] =
                        point->pressure_log_ratio[column] / point->pressure_pa;
                }
                jacobian[reference * count + reference] =
                    point->pressure_log_volume / point->pressure_pa;
                for (double& value : residual) {
                    value = -value;
                }
                const std::vector<double> step = solve_dense(
                    std::move(jacobian), std::move(residual)
                );
                double composition_step = 0.0;
                for (std::size_t index = 0; index < reference; ++index) {
                    composition_step = std::max(composition_step, std::abs(step[index]));
                }
                double step_scale = std::min(
                    1.0,
                    composition_step > 2.0 ? 2.0 / composition_step : 1.0
                );
                if (std::abs(step.back()) > 0.1) {
                    step_scale = std::min(step_scale, 0.1 / std::abs(step.back()));
                }
                bool improved = false;
                for (int reduction = 0; reduction < 10; ++reduction) {
                    const double factor = step_scale * std::pow(0.5, reduction);
                    std::vector<double> candidate = coordinates;
                    for (std::size_t index = 0; index < count; ++index) {
                        candidate[index] += factor * step[index];
                    }
                    try {
                        std::vector<double> midpoint(count);
                        for (std::size_t index = 0; index < count; ++index) {
                            midpoint[index] = 0.5 * (coordinates[index] + candidate[index]);
                        }
                        static_cast<void>(evaluate_vapor(
                            provider,
                            request,
                            vapor,
                            std::vector<double>(midpoint.begin(), midpoint.end() - 1),
                            midpoint.back(),
                            parameter_identities
                        ));
                        const VaporPoint trial = evaluate_vapor(
                            provider,
                            request,
                            vapor,
                            std::vector<double>(candidate.begin(), candidate.end() - 1),
                            candidate.back(),
                            parameter_identities
                        );
                        double trial_norm = std::abs(std::log(trial.pressure_pa / pressure_pa));
                        for (std::size_t row = 0; row < reference; ++row) {
                            trial_norm = std::max(
                                trial_norm,
                                std::abs(
                                    trial.chemical_potentials[row]
                                    - trial.chemical_potentials[reference]
                                    - liquid_potentials[row]
                                    + liquid_potentials[reference]
                                )
                            );
                        }
                        if (trial_norm < final_norm) {
                            coordinates = std::move(candidate);
                            improved = true;
                            break;
                        }
                    } catch (const GrepeFailure&) {
                    } catch (const std::exception&) {
                    }
                }
                if (!improved) {
                    throw GrepeFailure(
                        "numerical_convergence_failure",
                        "incipient composition-volume Newton step did not reduce the residual"
                    );
                }
            }
            if (!point.has_value() || final_norm > 2.0e-10) {
                throw GrepeFailure(
                    "numerical_convergence_failure",
                    "incipient composition-volume Newton iteration limit reached"
                );
            }
            if (reference != 0) {
                std::vector<double> chart(count * reference, 0.0);
                for (std::size_t component = 0; component < count; ++component) {
                    for (std::size_t coordinate = 0; coordinate < reference; ++coordinate) {
                        chart[component * reference + coordinate] = point->fractions[component]
                            * ((component == coordinate ? 1.0 : 0.0)
                               - point->fractions[coordinate]);
                    }
                }
                std::vector<double> hessian(reference * reference, 0.0);
                for (std::size_t row = 0; row < reference; ++row) {
                    for (std::size_t column = 0; column < reference; ++column) {
                        for (std::size_t component = 0; component < count; ++component) {
                            hessian[row * reference + column] +=
                                chart[component * reference + row]
                                * point->mu_log_ratio[component * reference + column];
                        }
                    }
                }
                double scale = 1.0;
                for (double value : hessian) {
                    scale = std::max(scale, std::abs(value));
                }
                std::vector<double> lower(reference * reference, 0.0);
                for (std::size_t row = 0; row < reference; ++row) {
                    for (std::size_t column = 0; column <= row; ++column) {
                        if (std::abs(hessian[row * reference + column]
                                     - hessian[column * reference + row])
                            > 2.0e-9 * scale) {
                            throw GrepeFailure(
                                "first_order_failure",
                                "incipient tangent-plane Hessian is not symmetric"
                            );
                        }
                        double value = 0.5 * (
                            hessian[row * reference + column]
                            + hessian[column * reference + row]
                        );
                        for (std::size_t inner = 0; inner < column; ++inner) {
                            value -= lower[row * reference + inner]
                                * lower[column * reference + inner];
                        }
                        if (row == column) {
                            if (!std::isfinite(value) || value <= 1.0e-12 * scale) {
                                throw GrepeFailure(
                                    "second_order_failure",
                                    "incipient composition is not a strict local tangent-plane minimum"
                                );
                            }
                            lower[row * reference + column] = std::sqrt(value);
                        } else {
                            lower[row * reference + column] =
                                value / lower[column * reference + column];
                        }
                    }
                }
            }
            bool duplicate = false;
            for (const auto& existing : accepted) {
                double distance = 0.0;
                for (std::size_t component = 0; component < count; ++component) {
                    distance = std::max(
                        distance,
                        std::abs(std::log(point->fractions[component]
                                          / existing.first.fractions[component]))
                    );
                }
                duplicate = duplicate || distance <= 2.0e-8;
            }
            if (!duplicate) {
                accepted.emplace_back(std::move(*point), final_norm);
            }
        } catch (const GrepeFailure& error) {
            last_failure = error;
        } catch (const std::exception& error) {
            last_failure = GrepeFailure("provider_domain_rejection", error.what());
        }
    }
    if (accepted.empty()) {
        if (continued) {
            throw GrepeFailure(
                "branch_loss",
                last_failure.has_value()
                    ? last_failure->what()
                    : "continued incipient composition-volume corrector failed"
            );
        }
        if (last_failure.has_value()
            && (last_failure->code() == "second_order_failure"
                || last_failure->code() == "singular_sensitivity"
                || last_failure->code() == "ill_conditioned_sensitivity")) {
            throw *last_failure;
        }
        throw GrepeFailure(
            "no_incipient_phase",
            last_failure.has_value()
                ? last_failure->what()
                : "no declared strict-stable vapor composition satisfies fugacity equality"
        );
    }
    if (continued) {
        return std::move(accepted.front().first);
    }
    std::sort(accepted.begin(), accepted.end(), [&](const auto& left, const auto& right) {
        const double left_distance = left.first.chemical_potentials[reference]
            - liquid_potentials[reference];
        const double right_distance = right.first.chemical_potentials[reference]
            - liquid_potentials[reference];
        return left_distance == right_distance
            ? left.first.fractions < right.first.fractions
            : left_distance < right_distance;
    });
    if (accepted.size() > 1) {
        const double first = accepted[0].first.chemical_potentials[reference]
            - liquid_potentials[reference];
        const double second = accepted[1].first.chemical_potentials[reference]
            - liquid_potentials[reference];
        if (std::abs(second - first)
            <= 2.0e-10 + std::max(accepted[0].second, accepted[1].second)) {
            throw GrepeFailure(
                "multiple_unresolved_branches",
                "multiple strict local incipient minima have unresolved tangent-plane distance"
            );
        }
    }
    return std::move(accepted.front().first);
}

struct IncipientClosure {
    LiquidEvaluation liquid;
    VaporPoint vapor;
    std::vector<double> vapor_fractions;
    std::vector<double> boundary_directional;
    std::vector<double> boundary_direct;
    std::vector<double> vapor_directional;
    std::vector<double> vapor_direct;
    double pressure_pa = 0.0;
    double vapor_pressure_pa = 0.0;
    double vapor_volume_m3 = 0.0;
    double vapor_packing_fraction = 0.0;
    double residual = 0.0;
    double residual_log_pressure_derivative = 0.0;
    std::size_t vapor_evaluations = 0;
    bool ideal = false;
};

IncipientClosure incipient_closure(
    const epcsaft::NativeModel& model,
    const epcsaft::NativeModel& incipient_model,
    const GrepeRequest& request,
    double pressure_pa
) {
    IncipientClosure closure;
    closure.pressure_pa = pressure_pa;
    closure.liquid = evaluate_liquid(model, request, pressure_pa);
    const ReactiveFinitePhaseResult& liquid = closure.liquid.solved.equilibrium;
    const auto vapor_position = std::find_if(
        request.phases.begin(), request.phases.end(), [](const GrepePhase& phase) {
            return !phase.finite;
        }
    );
    if (vapor_position == request.phases.end()) {
        throw GrepeFailure("unsupported_topology", "incipient vapor declaration is missing");
    }
    const GrepePhase& vapor = *vapor_position;
    const std::size_t count = vapor.support_indices.size();
    const std::size_t reference = count - 1;
    const std::size_t width = liquid.input_order.size();
    const std::size_t parameter_count = liquid.parameter_order.size();
    std::vector<double> liquid_potentials(count);
    const double log_rt = std::log(kGasConstant * request.temperature.value);
    for (std::size_t supported = 0; supported < count; ++supported) {
        liquid_potentials[supported] = liquid.phase_log_fugacities[0][
            vapor.support_indices[supported]
        ] - log_rt;
    }
    closure.boundary_directional.assign(width, 0.0);
    closure.boundary_direct.assign(parameter_count, 0.0);
    closure.vapor_directional.assign(width * count, 0.0);
    closure.vapor_direct.assign(parameter_count * count, 0.0);
    closure.ideal = vapor.ideal_gas;
    if (closure.ideal) {
        std::vector<double> fugacity(count);
        double boundary_pressure = 0.0;
        for (std::size_t supported = 0; supported < count; ++supported) {
            const std::size_t component = vapor.support_indices[supported];
            fugacity[supported] = std::exp(liquid.phase_log_fugacities[0][component]);
            boundary_pressure += fugacity[supported];
        }
        if (!std::isfinite(boundary_pressure) || boundary_pressure <= 0.0) {
            throw GrepeFailure(
                "no_incipient_phase", "ideal comparison has no positive incipient pressure"
            );
        }
        closure.vapor_fractions.resize(count);
        for (std::size_t supported = 0; supported < count; ++supported) {
            closure.vapor_fractions[supported] = fugacity[supported] / boundary_pressure;
        }
        for (std::size_t input = 0; input < width; ++input) {
            for (std::size_t supported = 0; supported < count; ++supported) {
                const std::size_t component = vapor.support_indices[supported];
                closure.boundary_directional[input] += fugacity[supported]
                    * liquid.input_log_fugacity_derivatives[component * width + input];
            }
            for (std::size_t supported = 0; supported < count; ++supported) {
                const std::size_t component = vapor.support_indices[supported];
                const double derivative = liquid.input_log_fugacity_derivatives[
                    component * width + input
                ];
                closure.vapor_directional[input * count + supported] = fugacity[supported]
                    * (derivative * boundary_pressure - closure.boundary_directional[input])
                    / (boundary_pressure * boundary_pressure);
            }
        }
        for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
            for (std::size_t supported = 0; supported < count; ++supported) {
                const std::size_t component = vapor.support_indices[supported];
                const std::size_t output = closure.liquid.log_fugacity_offset + component;
                const double partial = liquid.output_parameter_partials[
                    output * parameter_count + parameter
                ];
                closure.boundary_direct[parameter] += fugacity[supported] * partial;
            }
            for (std::size_t supported = 0; supported < count; ++supported) {
                const std::size_t component = vapor.support_indices[supported];
                const std::size_t output = closure.liquid.log_fugacity_offset + component;
                const double partial = liquid.output_parameter_partials[
                    output * parameter_count + parameter
                ];
                closure.vapor_direct[parameter * count + supported] = fugacity[supported]
                    * (partial * boundary_pressure - closure.boundary_direct[parameter])
                    / (boundary_pressure * boundary_pressure);
            }
        }
        closure.vapor_pressure_pa = boundary_pressure;
        closure.vapor_volume_m3 = kGasConstant * request.temperature.value / boundary_pressure;
        closure.vapor_packing_fraction = 0.0;
        closure.residual = std::log(boundary_pressure / pressure_pa);
        const std::size_t pressure = input_index(liquid, "pressure_pa");
        closure.residual_log_pressure_derivative = pressure_pa / boundary_pressure
                * closure.boundary_directional[pressure]
            - 1.0;
        closure.vapor.chemical_potentials = liquid_potentials;
        return closure;
    }

    ProviderContext provider(incipient_model);
    closure.vapor = solve_vapor(
        provider,
        request,
        vapor,
        pressure_pa,
        liquid_potentials,
        liquid.parameter_order
    );
    closure.vapor_evaluations = static_cast<std::size_t>(provider.evaluation_count());
    closure.vapor_fractions = closure.vapor.fractions;
    closure.vapor_pressure_pa = closure.vapor.pressure_pa;
    closure.vapor_volume_m3 = closure.vapor.volume_m3;
    closure.vapor_packing_fraction = closure.vapor.packing_fraction;
    std::vector<double> g_jacobian(reference * reference, 0.0);
    for (std::size_t row = 0; row < reference; ++row) {
        for (std::size_t column = 0; column < reference; ++column) {
            g_jacobian[row * reference + column] =
                closure.vapor.mu_log_ratio[row * reference + column]
                - closure.vapor.mu_log_ratio[reference * reference + column];
        }
    }
    std::vector<double> chart(count * reference, 0.0);
    for (std::size_t component = 0; component < count; ++component) {
        for (std::size_t coordinate = 0; coordinate < reference; ++coordinate) {
            chart[component * reference + coordinate] = closure.vapor.fractions[component]
                * ((component == coordinate ? 1.0 : 0.0)
                   - closure.vapor.fractions[coordinate]);
        }
    }
    const std::size_t pressure_input = input_index(liquid, "pressure_pa");
    for (std::size_t input = 0; input < width; ++input) {
        std::vector<double> potential(count, 0.0);
        const auto parameter = std::find(
            liquid.parameter_order.begin(),
            liquid.parameter_order.end(),
            liquid.input_order[input]
        );
        const bool active = parameter != liquid.parameter_order.end();
        const std::size_t parameter_index = static_cast<std::size_t>(
            parameter - liquid.parameter_order.begin()
        );
        for (std::size_t supported = 0; supported < count; ++supported) {
            const std::size_t component = vapor.support_indices[supported];
            const double vapor_base = input == pressure_input
                ? closure.vapor.mu_pressure[supported]
                : active
                ? closure.vapor.mu_parameters[supported * parameter_count + parameter_index]
                : 0.0;
            potential[supported] = vapor_base
                - liquid.input_log_fugacity_derivatives[component * width + input];
        }
        std::vector<double> right(reference, 0.0);
        for (std::size_t row = 0; row < reference; ++row) {
            right[row] = -(potential[row] - potential[reference]);
        }
        const std::vector<double> log_ratio = solve_dense(g_jacobian, std::move(right));
        double closure_direction = potential[reference];
        for (std::size_t coordinate = 0; coordinate < reference; ++coordinate) {
            closure_direction += closure.vapor.mu_log_ratio[
                reference * reference + coordinate
            ] * log_ratio[coordinate];
        }
        closure.boundary_directional[input] = pressure_pa * closure_direction;
        if (input == pressure_input) {
            closure.boundary_directional[input] += 1.0;
        }
        for (std::size_t component = 0; component < count; ++component) {
            for (std::size_t coordinate = 0; coordinate < reference; ++coordinate) {
                closure.vapor_directional[input * count + component] +=
                    chart[component * reference + coordinate] * log_ratio[coordinate];
            }
        }
    }
    closure.residual = closure.vapor.chemical_potentials[reference]
        - liquid_potentials[reference];
    closure.residual_log_pressure_derivative =
        closure.boundary_directional[pressure_input] - 1.0;
    return closure;
}

double incipient_distance(
    const IncipientClosure& left,
    const IncipientClosure& right
) {
    double distance = std::abs(std::log(left.pressure_pa / right.pressure_pa));
    for (std::size_t component = 0; component < left.vapor_fractions.size(); ++component) {
        distance = std::max(
            distance,
            std::abs(std::log(left.vapor_fractions[component]
                              / right.vapor_fractions[component]))
        );
    }
    return distance;
}

struct IncipientSolve {
    IncipientClosure closure;
    std::size_t closure_count = 0;
    std::size_t provider_evaluations = 0;
    std::size_t accepted_count = 0;
};

IncipientSolve solve_incipient_pressure(
    const epcsaft::NativeModel& model,
    const epcsaft::NativeModel& incipient_model,
    const GrepeRequest& request
) {
    const std::vector<double> starts = request.continuation_phases.empty()
        ? request.pressure.starts
        : std::vector<double>{request.continuation_pressure_pa};
    std::vector<IncipientClosure> accepted;
    std::vector<GrepeFailure> failures;
    std::size_t closure_count = 0;
    std::size_t provider_evaluations = 0;
    std::size_t accepted_count = 0;
    for (double pressure : starts) {
        std::optional<IncipientClosure> final;
        std::optional<GrepeFailure> failure;
        for (int iteration = 0; iteration < 18; ++iteration) {
            try {
                IncipientClosure current = incipient_closure(
                    model, incipient_model, request, pressure
                );
                ++closure_count;
                provider_evaluations += current.liquid.solved.equilibrium.provider_evaluation_count
                    + current.vapor_evaluations;
                if (std::abs(current.residual) <= 2.0e-10) {
                    final = std::move(current);
                    break;
                }
                const double slope = current.residual_log_pressure_derivative;
                if (std::abs(slope) < 1.0e-10) {
                    throw GrepeFailure(
                        "singular_sensitivity", "bubble closure slope is singular"
                    );
                }
                const double step = std::clamp(-current.residual / slope, -1.0, 1.0);
                const double candidate = std::exp(std::clamp(
                    std::log(pressure) + step,
                    std::log(request.pressure.bounds[0]),
                    std::log(request.pressure.bounds[1])
                ));
                if (candidate == pressure) {
                    throw GrepeFailure(
                        "no_incipient_phase",
                        "declared pressure interval excludes the Newton step"
                    );
                }
                pressure = candidate;
            } catch (const GrepeFailure& error) {
                failure = error;
                break;
            } catch (const std::exception& error) {
                failure = GrepeFailure("provider_domain_rejection", error.what());
                break;
            }
        }
        if (final.has_value()) {
            ++accepted_count;
            const auto duplicate = std::find_if(
                accepted.begin(), accepted.end(), [&](const IncipientClosure& existing) {
                    return incipient_distance(*final, existing) <= 2.0e-8;
                }
            );
            if (duplicate == accepted.end()) {
                accepted.push_back(std::move(*final));
            } else {
                const auto key = [](const IncipientClosure& value) {
                    return std::pair(value.pressure_pa, value.vapor_fractions);
                };
                if (key(*final) < key(*duplicate)) {
                    *duplicate = std::move(*final);
                }
            }
        } else {
            failures.push_back(
                failure.value_or(GrepeFailure(
                    "numerical_convergence_failure", "bubble Newton iteration limit reached"
                ))
            );
        }
    }
    if (accepted.empty()) {
        const std::vector<std::string_view> priority = {
            "multiple_unresolved_branches", "branch_loss", "active_bound_ambiguity",
            "parameter_topology_mismatch", "provider_domain_rejection",
            "ill_conditioned_sensitivity", "singular_sensitivity",
            "first_order_failure", "second_order_failure", "physical_rejection",
            "no_incipient_phase",
            "numerical_convergence_failure",
        };
        for (std::string_view code : priority) {
            const auto found = std::find_if(
                failures.begin(), failures.end(), [&](const GrepeFailure& failure) {
                    return failure.code() == code;
                }
            );
            if (found != failures.end()) {
                throw *found;
            }
        }
        throw GrepeFailure("no_incipient_phase", "no declared pressure start found a branch");
    }
    if (accepted.size() != 1) {
        throw GrepeFailure(
            "multiple_unresolved_branches", "declared starts found multiple bubble branches"
        );
    }
    return {
        std::move(accepted.front()), closure_count, provider_evaluations, accepted_count
    };
}

GrepeResult incipient_result(
    const GrepeRequest& request,
    IncipientSolve solved,
    const std::vector<double>* coordinate_direction = nullptr,
    double coordinate_value = 0.0
) {
    IncipientClosure& closure = solved.closure;
    ReactiveFinitePhaseResult& liquid = closure.liquid.solved.equilibrium;
    const auto vapor_position = std::find_if(
        request.phases.begin(), request.phases.end(), [](const GrepePhase& phase) {
            return !phase.finite;
        }
    );
    const std::size_t count = vapor_position->support_indices.size();
    const std::size_t component_count = request.reaction.species_ids.size();
    const std::size_t width = liquid.input_order.size();
    const std::size_t parameter_count = liquid.parameter_order.size();
    const std::size_t pressure_input = input_index(liquid, "pressure_pa");
    const bool fixed_pressure = coordinate_direction != nullptr;
    double scalar_slope = closure.boundary_directional[pressure_input] - 1.0;
    double residual_scale = 1.0;
    if (fixed_pressure) {
        if (coordinate_direction->size() != width) {
            throw GrepeFailure(
                "parameter_topology_mismatch", "fixed-pressure balance direction changed"
            );
        }
        residual_scale = closure.ideal ? closure.vapor_pressure_pa : closure.pressure_pa;
        scalar_slope = 0.0;
        for (std::size_t input = 0; input < width; ++input) {
            scalar_slope += closure.boundary_directional[input]
                * (*coordinate_direction)[input] / residual_scale;
        }
    }
    const double slope_scale = std::max(
        1.0,
        fixed_pressure
            ? std::abs(scalar_slope)
            : std::abs(closure.boundary_directional[pressure_input])
    );
    if (!std::isfinite(scalar_slope) || std::abs(scalar_slope) < 1.0e-12 * slope_scale) {
        throw GrepeFailure(
            "singular_sensitivity",
            fixed_pressure
                ? "fixed-pressure affine-feed closure is singular"
                : "coupled bubble Schur complement is singular"
        );
    }
    if (std::abs(scalar_slope) < 1.0e-8 * slope_scale) {
        throw GrepeFailure(
            "ill_conditioned_sensitivity",
            fixed_pressure
                ? "fixed-pressure affine-feed closure is ill-conditioned"
                : "coupled bubble Schur complement is ill-conditioned"
        );
    }
    std::vector<double> pressure_derivatives(parameter_count, 0.0);
    std::vector<double> vapor_derivatives(parameter_count * count, 0.0);
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        const std::size_t input = input_index(liquid, liquid.parameter_order[parameter]);
        const double coordinate_parameter = fixed_pressure
            ? -(closure.boundary_direct[parameter] + closure.boundary_directional[input])
                / residual_scale / scalar_slope
            : 0.0;
        pressure_derivatives[parameter] = fixed_pressure
            ? 0.0
            : -(closure.boundary_direct[parameter] + closure.boundary_directional[input])
                / scalar_slope;
        for (std::size_t component = 0; component < count; ++component) {
            double implicit = closure.vapor_directional[
                pressure_input * count + component
            ] * pressure_derivatives[parameter];
            if (fixed_pressure) {
                implicit = 0.0;
                for (std::size_t direction = 0; direction < width; ++direction) {
                    implicit += closure.vapor_directional[direction * count + component]
                        * (*coordinate_direction)[direction] * coordinate_parameter;
                }
            }
            vapor_derivatives[parameter * count + component] =
                closure.vapor_direct[parameter * count + component]
                + closure.vapor_directional[input * count + component] + implicit;
        }
    }

    GrepeResult result;
    result.formulation = fixed_pressure ? "intensive_boundary" : "finite_incipient";
    result.certificate_scope = fixed_pressure
        ? "fixed_T_one_finite_liquid_one_incipient_vapor_affine_feed_imposed_boundaries"
        : "fixed_T_one_finite_liquid_one_incipient_vapor";
    result.temperature_k = request.temperature.value;
    result.pressure_pa = closure.pressure_pa;
    result.source_reference_representation_residual_inf_norm =
        closure.liquid.solved.source_reference_representation_residual_inf_norm;
    result.branch_kind = fixed_pressure ? "intensive_boundary" : "finite_incipient_local_branch";
    result.branch_coordinates = fixed_pressure
        ? std::vector<double>{coordinate_value, closure.pressure_pa}
        : std::vector<double>{closure.pressure_pa};
    result.phase_mechanical_classes = {
        "strict_stable_local_minimum",
        closure.ideal ? "analytical_ideal_gas" : "strict_stable_vapor",
    };
    result.equilibrium = liquid;
    ReactiveFinitePhaseResult& equilibrium = result.equilibrium;
    for (std::size_t component = 0; component < component_count; ++component) {
        equilibrium.phase_chemical_potentials_over_rt[0][component] =
            equilibrium.phase_log_fugacities[0][component]
            - std::log(kGasConstant * request.temperature.value);
    }
    equilibrium.phase_amounts_mol.resize(2);
    equilibrium.phase_volumes_m3.resize(2);
    equilibrium.phase_pressures_pa.resize(2);
    equilibrium.phase_packing_fractions.resize(2);
    equilibrium.phase_chemical_potentials_over_rt.resize(2);
    equilibrium.phase_log_fugacities.resize(2);
    equilibrium.phase_amounts_mol[1].assign(component_count, 0.0);
    equilibrium.phase_chemical_potentials_over_rt[1].assign(component_count, 0.0);
    equilibrium.phase_log_fugacities[1].assign(component_count, 0.0);
    for (std::size_t supported = 0; supported < count; ++supported) {
        const std::size_t component = vapor_position->support_indices[supported];
        equilibrium.phase_amounts_mol[1][component] = closure.vapor_fractions[supported];
        equilibrium.phase_chemical_potentials_over_rt[1][component] =
            closure.vapor.chemical_potentials[supported];
        equilibrium.phase_log_fugacities[1][component] =
            closure.vapor.chemical_potentials[supported]
            + std::log(kGasConstant * request.temperature.value);
    }
    equilibrium.phase_volumes_m3[1] = closure.vapor_volume_m3;
    equilibrium.phase_pressures_pa[1] = closure.vapor_pressure_pa;
    equilibrium.phase_packing_fractions[1] = closure.vapor_packing_fraction;
    equilibrium.output_values.assign(request.outputs.size(), 0.0);
    equilibrium.output_input_derivatives.assign(
        request.outputs.size() * parameter_count, 0.0
    );
    for (std::size_t output = 0; output < request.outputs.size(); ++output) {
        const ReactiveFinitePhaseOutput& declared = request.outputs[output];
        const bool vapor_output = declared.phase_index
            == static_cast<std::size_t>(vapor_position - request.phases.begin());
        if (declared.selector == "system.pressure") {
            equilibrium.output_values[output] = closure.pressure_pa;
            for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
                equilibrium.output_input_derivatives[output * parameter_count + parameter] =
                    pressure_derivatives[parameter];
            }
        } else if (vapor_output && (
                declared.selector == "phase.mole_fraction"
                || declared.selector == "phase.partial_pressure")) {
            double fraction = 0.0;
            for (std::size_t component = 0; component < count; ++component) {
                fraction += declared.coefficients[component]
                    * closure.vapor_fractions[component];
            }
            equilibrium.output_values[output] = declared.selector == "phase.partial_pressure"
                ? closure.pressure_pa * fraction
                : fraction;
            for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
                double derivative = 0.0;
                for (std::size_t component = 0; component < count; ++component) {
                    derivative += declared.coefficients[component]
                        * vapor_derivatives[parameter * count + component];
                }
                if (declared.selector == "phase.partial_pressure") {
                    derivative = pressure_derivatives[parameter] * fraction
                        + closure.pressure_pa * derivative;
                }
                equilibrium.output_input_derivatives[output * parameter_count + parameter] =
                    derivative;
            }
        } else {
            const int position = closure.liquid.output_positions[output];
            if (position < 0) {
                throw GrepeFailure(
                    "observation_out_of_support",
                    "finite-incipient output is unavailable on its declared support"
                );
            }
            equilibrium.output_values[output] = liquid.output_values[
                static_cast<std::size_t>(position)
            ];
            for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
                const std::size_t input = input_index(liquid, liquid.parameter_order[parameter]);
                double implicit = liquid.output_input_derivatives[
                    static_cast<std::size_t>(position) * width + pressure_input
                ] * pressure_derivatives[parameter];
                if (fixed_pressure) {
                    double coordinate_output = 0.0;
                    for (std::size_t direction = 0; direction < width; ++direction) {
                        coordinate_output += liquid.output_input_derivatives[
                            static_cast<std::size_t>(position) * width + direction
                        ] * (*coordinate_direction)[direction];
                    }
                    const double coordinate_parameter = -(
                        closure.boundary_direct[parameter]
                        + closure.boundary_directional[input]
                    ) / residual_scale / scalar_slope;
                    implicit = coordinate_output * coordinate_parameter;
                }
                equilibrium.output_input_derivatives[output * parameter_count + parameter] =
                    liquid.output_input_derivatives[
                        static_cast<std::size_t>(position) * width + input
                    ] + implicit;
            }
        }
    }
    equilibrium.input_order = equilibrium.parameter_order;
    equilibrium.provider_evaluation_count = solved.provider_evaluations;
    result.evidence = {
        {"bubble_log_closure_abs", std::abs(closure.residual)},
        {"bubble_log_pressure_derivative", closure.residual_log_pressure_derivative},
        {"balance_inf_norm", equilibrium.balance_inf_norm},
        {"reaction_affinity_inf_norm", equilibrium.reaction_affinity_inf_norm},
        {"pressure_relative_inf_norm", equilibrium.pressure_relative_inf_norm},
        {"phase_pressure_residual_inf_pa", std::abs(
            equilibrium.phase_pressures_pa[0] - equilibrium.phase_pressures_pa[1]
        )},
        {"closure_evaluation_count", static_cast<double>(solved.closure_count)},
        {"provider_evaluation_count", static_cast<double>(solved.provider_evaluations)},
        {"accepted_candidate_count", static_cast<double>(solved.accepted_count)},
        {"distinct_branch_count", 1.0},
    };
    return result;
}

GrepeResult failed_result(
    const GrepeRequest& request,
    const std::string& formulation,
    const std::string& certificate_scope,
    const GrepeFailure& failure
) {
    GrepeResult result;
    result.formulation = formulation;
    result.certificate_scope = certificate_scope;
    result.failure_code = failure.code();
    result.temperature_k = request.temperature.value;
    result.pressure_pa = request.pressure.value;
    result.equilibrium.solver_status = "failed";
    result.equilibrium.numerical_status = "failed";
    result.equilibrium.physical_status = failure.code() == "physical_rejection"
        || failure.code() == "second_order_failure"
        ? "failed" : "not_adjudicated";
    result.equilibrium.provider_domain_status = failure.code() == "provider_domain_rejection"
        ? "failed" : "not_adjudicated";
    result.equilibrium.failure_reason = failure.what();
    return result;
}

std::vector<double> feed_input_direction(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    const ReactiveFinitePhaseResult& equilibrium
) {
    if (!request.feed.has_value() || !request.feed->affine
        || request.feed->component_coordinate_matrix.columns != 1) {
        throw GrepeFailure(
            "parameter_topology_mismatch", "intensive boundary requires one affine feed coordinate"
        );
    }
    std::vector<double> result(equilibrium.input_order.size(), 0.0);
    ReactionSystemInput input = request.reaction;
    input.provider_fingerprint = ProviderContext(model).fingerprint();
    input.temperature_k = request.temperature.value;
    input.pressure_pa = request.pressure.value;
    const CompiledReactionSystem compiled = expand_species_face(
        compile_reaction_system(input, true), input
    );
    for (std::size_t balance = 0; balance < compiled.balance_matrix.rows; ++balance) {
        double direction = 0.0;
        for (std::size_t component = 0;
             component < compiled.balance_matrix.columns;
             ++component) {
            direction += compiled.balance_matrix(balance, component)
                * request.feed->component_coordinate_matrix(component, 0);
        }
        result[input_index(
            equilibrium, "balance_total[" + std::to_string(balance) + "]"
        )] = direction;
    }
    return result;
}

GrepeResult solve_fixed_pressure_incipient_boundary(
    const epcsaft::NativeModel& model,
    const epcsaft::NativeModel& incipient_model,
    const GrepeRequest& request
) {
    const GrepeIntensiveBoundary& boundary = *request.boundary;
    const GrepeFeed& feed = *request.feed;
    const double lower = feed.coordinates[0].bounds[0];
    const double upper = feed.coordinates[0].bounds[1];
    std::size_t evaluations = 0;
    std::size_t provider_evaluations = 0;
    auto evaluate = [&](double coordinate) {
        GrepeRequest resolved = resolved_request(request, coordinate);
        IncipientClosure closure = incipient_closure(
            model, incipient_model, resolved, boundary.target
        );
        ++evaluations;
        provider_evaluations += closure.liquid.solved.equilibrium.provider_evaluation_count
            + closure.vapor_evaluations;
        const std::vector<double> direction = feed_input_direction(
            model, request, closure.liquid.solved.equilibrium
        );
        const double scale = closure.ideal ? closure.vapor_pressure_pa : closure.pressure_pa;
        double slope = 0.0;
        for (std::size_t input = 0; input < direction.size(); ++input) {
            slope += closure.boundary_directional[input] * direction[input] / scale;
        }
        return std::tuple(
            std::move(resolved), std::move(closure), std::move(direction), slope
        );
    };
    double coordinate = boundary.initial;
    std::optional<std::tuple<GrepeRequest, IncipientClosure, std::vector<double>, double>> state;
    for (int iteration = 0; iteration < 64; ++iteration) {
        state = evaluate(coordinate);
        IncipientClosure& closure = std::get<1>(*state);
        const double residual = closure.residual;
        const double slope = std::get<3>(*state);
        if (std::abs(residual) <= 2.0e-10) {
            break;
        }
        if (!std::isfinite(slope) || std::abs(slope) <= 1.0e-13) {
            throw GrepeFailure(
                "singular_sensitivity", "fixed-pressure affine-feed closure is singular"
            );
        }
        const double candidate = std::clamp(coordinate - residual / slope, lower, upper);
        if (candidate == coordinate) {
            throw GrepeFailure(
                "no_incipient_phase", "fixed-pressure boundary is not bracketed"
            );
        }
        bool improved = false;
        double length = 1.0;
        for (int reduction = 0; reduction < 20; ++reduction) {
            const double trial_coordinate = coordinate + length * (candidate - coordinate);
            auto trial = evaluate(trial_coordinate);
            if (std::abs(std::get<1>(trial).residual) < std::abs(residual)) {
                coordinate = trial_coordinate;
                state = std::move(trial);
                improved = true;
                break;
            }
            length *= 0.5;
        }
        if (!improved) {
            throw GrepeFailure(
                "numerical_convergence_failure",
                "fixed-pressure boundary step did not reduce the residual"
            );
        }
    }
    if (!state.has_value() || std::abs(std::get<1>(*state).residual) > 2.0e-10) {
        throw GrepeFailure(
            "numerical_convergence_failure", "fixed-pressure boundary iteration limit reached"
        );
    }
    IncipientSolve solved{
        std::move(std::get<1>(*state)), evaluations, provider_evaluations, 1
    };
    GrepeResult result = incipient_result(
        std::get<0>(*state),
        std::move(solved),
        &std::get<2>(*state),
        coordinate
    );
    return result;
}

GrepeResult solve_homogeneous_boundary(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request
) {
    const GrepeIntensiveBoundary& boundary = *request.boundary;
    const GrepeFeed& feed = *request.feed;
    const auto phase = std::find_if(
        request.phases.begin(), request.phases.end(), [&](const GrepePhase& candidate) {
            return candidate.identity == boundary.phase_identity;
        }
    );
    if (phase == request.phases.end()) {
        throw GrepeFailure(
            "observation_out_of_support", "intensive boundary names an undeclared phase"
        );
    }
    const std::size_t phase_index = static_cast<std::size_t>(phase - request.phases.begin());
    const double lower = feed.coordinates[0].bounds[0];
    const double upper = feed.coordinates[0].bounds[1];
    auto evaluate = [&](double coordinate) {
        GrepeRequest resolved = resolved_request(request, coordinate);
        ReactiveFinitePhaseOutput boundary_output{
            boundary.selector,
            phase_index,
            boundary.coefficients,
            {},
        };
        resolved.outputs.push_back(std::move(boundary_output));
        GrepeResult owner = solve_finite_result(model, resolved);
        if (!owner.equilibrium.accepted) {
            throw GrepeFailure(
                owner.equilibrium.provider_domain_status == "failed"
                    ? "provider_domain_rejection" : "numerical_convergence_failure",
                owner.equilibrium.failure_reason.empty()
                    ? owner.equilibrium.solver_status
                    : owner.equilibrium.failure_reason
            );
        }
        const ReactiveFinitePhaseResult& native = owner.equilibrium;
        const std::vector<double> direction = feed_input_direction(model, request, native);
        const std::size_t boundary_row = request.outputs.size();
        const std::size_t width = native.input_order.size();
        double derivative = 0.0;
        for (std::size_t input = 0; input < width; ++input) {
            derivative += native.output_input_derivatives[boundary_row * width + input]
                * direction[input];
        }
        const double value = native.output_values[boundary_row];
        const double residual = boundary.selector == "system.pressure"
            ? std::log(value / boundary.target) : value - boundary.target;
        if (boundary.selector == "system.pressure") {
            derivative /= value;
        }
        return std::tuple(
            std::move(owner), std::move(direction), residual, derivative
        );
    };
    double coordinate = boundary.initial;
    std::optional<std::tuple<GrepeResult, std::vector<double>, double, double>> selected;
    for (int iteration = 0; iteration < 32; ++iteration) {
        selected = evaluate(coordinate);
        const double residual = std::get<2>(*selected);
        const double derivative = std::get<3>(*selected);
        if (std::abs(residual) <= 2.0e-11) {
            break;
        }
        if (std::abs(derivative) <= 1.0e-13) {
            throw GrepeFailure(
                "singular_sensitivity", "intensive-boundary Jacobian is singular"
            );
        }
        const double step = -residual / derivative;
        double maximum_step = 1.0;
        if (step > 0.0) {
            maximum_step = std::min(maximum_step, 0.95 * (upper - coordinate) / step);
        } else if (step < 0.0) {
            maximum_step = std::min(maximum_step, 0.95 * (lower - coordinate) / step);
        }
        bool improved = false;
        double length = std::max(0.0, maximum_step);
        for (int reduction = 0; reduction < 20; ++reduction) {
            auto trial = evaluate(coordinate + length * step);
            if (std::abs(std::get<2>(trial)) < std::abs(residual)) {
                coordinate += length * step;
                selected = std::move(trial);
                improved = true;
                break;
            }
            length *= 0.5;
        }
        if (!improved) {
            throw GrepeFailure(
                "numerical_convergence_failure",
                "intensive-boundary Newton step did not reduce the residual"
            );
        }
    }
    if (!selected.has_value() || std::abs(std::get<2>(*selected)) > 2.0e-11) {
        throw GrepeFailure(
            "numerical_convergence_failure", "intensive-boundary solve exceeded its iteration budget"
        );
    }
    GrepeResult result = std::move(std::get<0>(*selected));
    ReactiveFinitePhaseResult& equilibrium = result.equilibrium;
    const std::vector<double>& direction = std::get<1>(*selected);
    const std::size_t width = equilibrium.input_order.size();
    const std::size_t boundary_row = request.outputs.size();
    const double boundary_derivative = std::get<3>(*selected);
    const std::size_t parameter_count = equilibrium.parameter_order.size();
    std::vector<double> total(request.outputs.size() * parameter_count, 0.0);
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        const std::size_t input = input_index(equilibrium, equilibrium.parameter_order[parameter]);
        double boundary_parameter = equilibrium.output_input_derivatives[
            boundary_row * width + input
        ];
        if (boundary.selector == "system.pressure") {
            boundary_parameter /= equilibrium.output_values[boundary_row];
        }
        const double coordinate_parameter = -boundary_parameter / boundary_derivative;
        for (std::size_t output = 0; output < request.outputs.size(); ++output) {
            double coordinate_output = 0.0;
            for (std::size_t direction_index = 0; direction_index < width; ++direction_index) {
                coordinate_output += equilibrium.output_input_derivatives[
                    output * width + direction_index
                ] * direction[direction_index];
            }
            total[output * parameter_count + parameter] =
                equilibrium.output_input_derivatives[output * width + input]
                + coordinate_output * coordinate_parameter;
        }
    }
    equilibrium.output_values.resize(request.outputs.size());
    equilibrium.output_input_derivatives = std::move(total);
    equilibrium.input_order = equilibrium.parameter_order;
    result.formulation = "intensive_boundary";
    result.certificate_scope = "fixed_TP_affine_feed_imposed_intensive_boundaries";
    result.branch_kind = "intensive_boundary";
    result.branch_coordinates = {coordinate};
    return result;
}

PhaseSupport grepe_phase_support(const GrepePhase& phase) {
    if (phase.support_name == "all_components") {
        return PhaseSupport::AllComponents;
    }
    if (phase.support_name == "neutral_components") {
        return PhaseSupport::NeutralComponents;
    }
    return PhaseSupport::DeclaredNeutralComponents;
}

double sensitivity_value(
    const SolvedStateSensitivity& sensitivity,
    const std::string& identity,
    std::size_t parameter
) {
    const auto found = std::find(
        sensitivity.output_identities.begin(),
        sensitivity.output_identities.end(),
        identity
    );
    if (found == sensitivity.output_identities.end()) {
        throw GrepeFailure(
            "owner_contract_failure", "coexistence sensitivity output layout changed"
        );
    }
    const std::size_t output = static_cast<std::size_t>(
        found - sensitivity.output_identities.begin()
    );
    return sensitivity.output_parameter_jacobian[
        output * sensitivity.parameter_identities.size() + parameter
    ];
}

GrepeResult solve_temperature_boundary(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request
) {
    const auto finite = std::find_if(
        request.phases.begin(), request.phases.end(), [](const GrepePhase& phase) {
            return phase.finite;
        }
    );
    if (request.phases.size() != 2 || finite == request.phases.end()) {
        throw GrepeFailure(
            "unsupported_topology", "temperature boundary requires two declared phases"
        );
    }
    const std::size_t reference_index = static_cast<std::size_t>(
        finite - request.phases.begin()
    );
    ReactionSystemInput reaction = effective_reaction(request);
    std::vector<double> reference = reaction.feed_amounts;
    const double total = std::accumulate(reference.begin(), reference.end(), 0.0);
    for (double& fraction : reference) {
        fraction /= total;
    }
    FixedTCoexistenceProblem problem;
    problem.phases.reserve(2);
    for (const GrepePhase& phase : request.phases) {
        problem.phases.push_back({phase.role, grepe_phase_support(phase)});
    }
    problem.reference_phase_index = reference_index;
    problem.reference_mole_fractions = reference;
    for (std::size_t index = 0; index < request.continuation_phases.size(); ++index) {
        const GrepeContinuationPhase& continued = request.continuation_phases[index];
        const auto declared = std::find_if(
            request.phases.begin(), request.phases.end(), [&](const GrepePhase& phase) {
                return phase.identity == continued.identity;
            }
        );
        if (declared != request.phases.end()) {
            problem.hints.push_back({
                static_cast<std::size_t>(declared - request.phases.begin()),
                continued.mole_fractions,
                1.0 / continued.molar_volume_m3_per_mol,
            });
        }
    }
    const double target_pressure = request.pressure.value;
    std::size_t evaluation_count = 0;
    std::size_t provider_evaluations = 0;
    auto evaluate = [&](double temperature_k) {
        ProviderContext reference_provider(model);
        ProviderContext partner_provider(model);
        problem.temperature_k = temperature_k;
        FixedTCoexistenceResult solved = solve_fixed_t_coexistence(
            reference_provider, partner_provider, problem
        );
        ++evaluation_count;
        provider_evaluations += static_cast<std::size_t>(
            std::max(0, solved.diagnostics.provider_evaluation_count)
        );
        if (solved.diagnostics.physical_status != "accepted"
            || solved.phases.size() != 2) {
            const std::string reason = solved.diagnostics.failure_reason.empty()
                ? "fixed-temperature coexistence was not certified"
                : solved.diagnostics.failure_reason;
            throw GrepeFailure(
                reason.find("provider") != std::string::npos
                    || reason.find("domain") != std::string::npos
                    ? "provider_domain_rejection" : "physical_rejection",
                reason
            );
        }
        if (!solved.temperature_sensitivity.has_value()
            || solved.temperature_sensitivity->status != "available") {
            throw GrepeFailure(
                "owner_contract_failure",
                "fixed-temperature coexistence temperature derivative is unavailable"
            );
        }
        const double slope = sensitivity_value(
            *solved.temperature_sensitivity, "boundary.pressure_pa", 0
        );
        const double residual = solved.pressure_pa - target_pressure;
        return std::tuple(std::move(solved), residual, slope);
    };
    const double lower = request.temperature.bounds[0];
    const double upper = request.temperature.bounds[1];
    auto lower_state = evaluate(lower);
    auto upper_state = evaluate(upper);
    if (std::get<1>(lower_state) * std::get<1>(upper_state) > 0.0) {
        throw GrepeFailure(
            "physical_rejection", "fixed_pressure_temperature_not_bracketed"
        );
    }
    double left = lower;
    double right = upper;
    double coordinate = std::clamp(request.temperature.value, left, right);
    std::optional<std::tuple<FixedTCoexistenceResult, double, double>> selected;
    const double tolerance = 1.0e-10 * std::max(1.0, target_pressure);
    for (int iteration = 0; iteration < 64; ++iteration) {
        selected = evaluate(coordinate);
        const double residual = std::get<1>(*selected);
        if (std::abs(residual) <= tolerance) {
            break;
        }
        if (residual * std::get<1>(lower_state) <= 0.0) {
            right = coordinate;
            upper_state = *selected;
        } else {
            left = coordinate;
            lower_state = *selected;
        }
        const double slope = std::get<2>(*selected);
        double candidate = coordinate - residual / slope;
        if (!std::isfinite(candidate) || !(left < candidate && candidate < right)) {
            candidate = 0.5 * (left + right);
        }
        coordinate = candidate;
    }
    if (!selected.has_value() || std::abs(std::get<1>(*selected)) > tolerance) {
        throw GrepeFailure(
            "numerical_convergence_failure", "fixed_pressure_temperature_iteration_limit"
        );
    }
    FixedTCoexistenceResult& solved = std::get<0>(*selected);
    for (std::size_t phase = 0; phase < solved.phases.size(); ++phase) {
        const auto interval = request.phases[phase].packing_fraction_interval;
        if (!(interval[0] < solved.phases[phase].packing_fraction
              && solved.phases[phase].packing_fraction < interval[1])) {
            throw GrepeFailure(
                "provider_domain_rejection",
                "fixed-pressure coexistence left a declared packing-fraction interval"
            );
        }
    }
    const std::size_t parameter_count = ProviderContext(model).active_parameters().identities.size();
    if (parameter_count != 0
        && (!solved.sensitivity.has_value() || solved.sensitivity->status != "available")) {
        throw GrepeFailure(
            "owner_contract_failure",
            "fixed-pressure coexistence did not provide exact total derivatives"
        );
    }
    const double pressure_temperature = std::get<2>(*selected);
    std::vector<double> temperature_parameters(parameter_count, 0.0);
    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
        temperature_parameters[parameter] = -sensitivity_value(
            *solved.sensitivity, "boundary.pressure_pa", parameter
        ) / pressure_temperature;
    }
    GrepeResult result;
    result.formulation = "temperature_boundary";
    result.certificate_scope = "fixed_p_declared_coexistence";
    result.temperature_k = coordinate;
    result.pressure_pa = target_pressure;
    result.branch_kind = "fixed_p_declared_coexistence";
    result.branch_coordinates = {coordinate, static_cast<double>(reference_index)};
    result.phase_mechanical_classes.reserve(2);
    ReactiveFinitePhaseResult& equilibrium = result.equilibrium;
    equilibrium.accepted = true;
    equilibrium.solver_status = "converged";
    equilibrium.numerical_status = "passed";
    equilibrium.physical_status = "passed";
    equilibrium.provider_domain_status = "passed";
    equilibrium.parameter_order = ProviderContext(model).active_parameters().identities;
    equilibrium.input_order = equilibrium.parameter_order;
    for (const PhaseEquilibriumState& phase : solved.phases) {
        equilibrium.phase_amounts_mol.push_back(phase.mole_fractions);
        equilibrium.phase_volumes_m3.push_back(phase.molar_volume_m3_per_mol);
        equilibrium.phase_pressures_pa.push_back(phase.pressure_pa);
        equilibrium.phase_packing_fractions.push_back(phase.packing_fraction);
        equilibrium.phase_chemical_potentials_over_rt.push_back(
            phase.chemical_potential_over_rt
        );
        std::vector<double> log_fugacity = phase.chemical_potential_over_rt;
        for (double& value : log_fugacity) {
            value += std::log(kGasConstant * coordinate);
        }
        equilibrium.phase_log_fugacities.push_back(std::move(log_fugacity));
        result.phase_mechanical_classes.push_back(phase.mechanical_class);
    }
    equilibrium.output_values.assign(request.outputs.size(), 0.0);
    equilibrium.output_input_derivatives.assign(
        request.outputs.size() * parameter_count, 0.0
    );
    for (std::size_t output = 0; output < request.outputs.size(); ++output) {
        const ReactiveFinitePhaseOutput& declared = request.outputs[output];
        std::string identity;
        std::vector<double> fixed(parameter_count, 0.0);
        double temperature_partial = 0.0;
        if (declared.selector == "system.temperature") {
            equilibrium.output_values[output] = coordinate;
            for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
                equilibrium.output_input_derivatives[output * parameter_count + parameter] =
                    temperature_parameters[parameter];
            }
            continue;
        }
        if (declared.selector == "system.pressure" || declared.selector == "phase.pressure") {
            equilibrium.output_values[output] = target_pressure;
            continue;
        }
        const PhaseEquilibriumState& phase = solved.phases[declared.phase_index];
        if (declared.selector == "phase.mole_fraction") {
            for (std::size_t component = 0; component < phase.mole_fractions.size(); ++component) {
                equilibrium.output_values[output] += declared.coefficients[component]
                    * phase.mole_fractions[component];
                identity = "phase[" + std::to_string(declared.phase_index)
                    + "].mole_fraction[" + phase.supported_component_ids[component] + "]";
                if (parameter_count != 0) {
                    for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
                        fixed[parameter] += declared.coefficients[component]
                            * sensitivity_value(*solved.sensitivity, identity, parameter);
                    }
                }
                temperature_partial += declared.coefficients[component]
                    * sensitivity_value(*solved.temperature_sensitivity, identity, 0);
            }
        } else if (declared.selector == "phase.molar_density") {
            equilibrium.output_values[output] = phase.molar_density_mol_m3;
            identity = "phase[" + std::to_string(declared.phase_index)
                + "].molar_volume_m3_per_mol";
            temperature_partial = -sensitivity_value(
                *solved.temperature_sensitivity, identity, 0
            ) * phase.molar_density_mol_m3 * phase.molar_density_mol_m3;
            if (parameter_count != 0) {
                for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
                    fixed[parameter] = -sensitivity_value(
                        *solved.sensitivity, identity, parameter
                    ) * phase.molar_density_mol_m3 * phase.molar_density_mol_m3;
                }
            }
        } else {
            throw GrepeFailure(
                "observation_out_of_support", "temperature-boundary output is unavailable"
            );
        }
        for (std::size_t parameter = 0; parameter < parameter_count; ++parameter) {
            equilibrium.output_input_derivatives[output * parameter_count + parameter] =
                fixed[parameter] + temperature_partial * temperature_parameters[parameter];
        }
    }
    equilibrium.provider_evaluation_count = provider_evaluations;
    result.evidence = {
        {"temperature_k", coordinate},
        {"pressure_residual_inf", std::abs(std::get<1>(*selected))},
        {"optimizer_solve_count", static_cast<double>(evaluation_count)},
        {"provider_evaluation_count", static_cast<double>(provider_evaluations)},
    };
    return result;
}

}  // namespace

GrepeResult solve_grepe(
    const epcsaft::NativeModel& model,
    const epcsaft::NativeModel& incipient_model,
    const GrepeRequest& request
) {
    const std::size_t incipient_count = static_cast<std::size_t>(std::count_if(
        request.phases.begin(),
        request.phases.end(),
        [](const GrepePhase& phase) { return !phase.finite; }
    ));
    try {
        if (request.temperature.solved) {
            return solve_temperature_boundary(model, request);
        }
        if (request.boundary.has_value()) {
            if (incipient_count == 1
                && request.boundary->selector == "system.pressure") {
                return solve_fixed_pressure_incipient_boundary(
                    model, incipient_model, request
                );
            }
            if (incipient_count == 0) {
                return solve_homogeneous_boundary(model, request);
            }
            return unsupported("declared intensive boundary topology is unavailable");
        }
        if (incipient_count == 1 && request.pressure.solved) {
            return incipient_result(
                request, solve_incipient_pressure(model, incipient_model, request)
            );
        }
        if (incipient_count != 0 || request.pressure.solved) {
            return unsupported("declared finite-incipient topology is unavailable");
        }
        return solve_finite_result(model, request);
    } catch (const GrepeFailure& failure) {
        return failed_result(
            request,
            request.temperature.solved
                ? "temperature_boundary"
                : request.boundary.has_value()
                ? "intensive_boundary"
                : incipient_count == 1 ? "finite_incipient" : "finite",
            request.temperature.solved
                ? "fixed_p_declared_coexistence"
                : request.boundary.has_value()
                ? incipient_count == 1
                    ? "fixed_T_one_finite_liquid_one_incipient_vapor_affine_feed_imposed_boundaries"
                    : "fixed_TP_affine_feed_imposed_intensive_boundaries"
                : incipient_count == 1
                ? "fixed_T_one_finite_liquid_one_incipient_vapor"
                : "fixed_TP_configured_finite_phases_local_minimum",
            failure
        );
    } catch (const std::invalid_argument& error) {
        const std::string diagnostic = error.what();
        if (diagnostic.find("Provider") == std::string::npos
            && diagnostic.find("source domain") == std::string::npos
            && diagnostic.find("packing interval") == std::string::npos) {
            throw;
        }
        GrepeResult result;
        result.formulation = "finite";
        result.certificate_scope = "fixed_TP_configured_finite_phases_local_minimum";
        result.failure_code = "provider_domain_rejection";
        result.temperature_k = request.temperature.value;
        result.pressure_pa = request.pressure.value;
        result.equilibrium.solver_status = "not_started";
        result.equilibrium.numerical_status = "not_adjudicated";
        result.equilibrium.physical_status = "not_adjudicated";
        result.equilibrium.provider_domain_status = "failed";
        result.equilibrium.failure_reason = diagnostic;
        return result;
    }
}

GrepeResult evaluate_grepe_local(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    const std::vector<double>& local_variables,
    const std::vector<double>& local_lower_bounds,
    const std::vector<double>& local_upper_bounds
) {
    return evaluate_finite_result(
        model, request, local_variables, local_lower_bounds, local_upper_bounds
    );
}

GrepeResult finalize_grepe_local(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request,
    const std::vector<double>& local_variables,
    const std::vector<double>& local_lower_bounds,
    const std::vector<double>& local_upper_bounds
) {
    return finalize_finite_result(
        model, request, local_variables, local_lower_bounds, local_upper_bounds
    );
}

GrepeResult initialize_grepe_local(
    const epcsaft::NativeModel& model,
    const GrepeRequest& request
) {
    return initialize_finite_result(model, request);
}

}  // namespace epcsaft_engine::equilibrium
