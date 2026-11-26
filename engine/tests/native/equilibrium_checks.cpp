#include "recording_support.hpp"

// The case registry; meson reads it for the per-case tests. Each case function is defined in
// one equilibrium_*_checks.cpp unit.
namespace equilibrium_checks {
void check_discovery_work();
void check_concurrent_first_discovery();
void check_s1();
void check_ideal_gas();
void check_species_thermochemistry();
void check_pressure_basis_quotient();
void check_ideal_mole_fraction();
void check_ideal_molarity();
void check_s2();
void check_reaction_quotient_pressure();
void check_reaction_quotient_actions();
void check_neutral_reference_permutations();
void check_mixed_reference_work();
void check_trace_reaction_balance();
void check_s6();
void check_s7();
void check_c1_multiplier_scale();
void check_c1_ideal_caller_chart();
void check_f5();
void check_s11();
void check_s13();
void check_c3_availability();
void check_s16();
void check_e2_thermal();
void check_e2_solve();
void check_e2_reuse();
void check_e5_neutral();
void check_e5_pins();
void check_e5_reactive();
void check_e5_associating();
void check_e5_electrolyte();
void check_e5_roots();
void check_e5_thermal();
void check_e5_refusals();
void check_e5_identities();
void check_phase_removal_inventory();
void check_tangent_plane_curvature();
void check_density_correlated_ion_suppression_admission();
void check_rank_probe_permittivity_denominator();
void check_ion_suppression_reacting_liquid_actions();
void check_tangent_plane_admission();
void check_prepared_supplied_trials();
void check_feed_charge_permutation();
void check_exact_zero_trace_dot();
void check_liquid_trial_charge_permutation();
void check_trace_ion_curvature();
void check_candidate_publication();
void check_action_scaled_jacobian_overflow();
void check_budget_retry();
void check_physical_amount_work();
void check_reaction_relaxed_split();
void check_capped_gibbs_descent();
void check_mdea_endpoint_gibbs();
void check_associating_isomer_gibbs();
void check_mdea_full_stability();
void check_piperazine_a12s_bubble_stability();
void check_mdea_bubble_work_limit();
void check_mdea_bubble_recovery();
void check_mdea_bubble_refusal_work();
void check_bubble_branch_end();
void check_reaction_relaxed_search();
void check_reaction_relaxed_pz_row();
void check_mea_census_state_3();
void check_mea_trace_ion_discovery();
void check_reference_parameter_coverage();
void check_mea_discovery_omitted_species();
void check_electric_potential_difference();
void check_prescribed_liquid_coincidence();
void check_electric_potential_difference_phase_order();
void check_vanished_water_liquid_seed();
void check_mea_below_bubble_hold();
void check_mea_packet3_trace_root();
void check_stability_grid();
void check_condensed_trial_recovery();
void check_solve_error_curvature();
void check_stability_inventory_face();
void check_fixed_state_plane_uncertainty();
void check_stability_phase_comparison();
void check_stability_amount_correction();
void check_trial_refinement();
void check_normalized_reaction_affinity();
void check_reacting_liquid_split();
void check_declared_reacting_liquid_split();
void check_gross_2002_bubble();
void check_gross_2002_cooling_onset();
void check_density_branch_admission();
void check_tangent_plane_jacobian();
}  // namespace equilibrium_checks
namespace caloric_specification_checks {
void admission();
void actions();
void curvature();
void roundtrips();
void volume_specification();
void internal_energy_specification();
void reactions();
void apparent_heat_capacity();
}  // namespace caloric_specification_checks
namespace critical_envelope_checks {
void pure_critical();
void binary_critical();
void spinodal();
void curvature();
void mixture();
void pure();
void azeotrope();
}  // namespace critical_envelope_checks
namespace reactive_caloric_checks {
void cheap();
void cp_f1(bool ionic_case, bool fixed);
void cp_f2();
void cp_f3();
void run(bool ionic, bool fixed);
void unreactive();
void boundary(bool nonidentical, bool free_temperature, bool fixed);
namespace caloric_ph_ps { void run(); }
}  // namespace reactive_caloric_checks

int main(int argc, char** argv) {
    using namespace equilibrium_checks;
    return recording_support::run(argc, argv, {
        {{"native.equilibrium.capped-gibbs-descent", "capped-gibbs-descent", "Automatic_fast"}, check_capped_gibbs_descent},
        {{"native.equilibrium.mdea-endpoint-gibbs", "mdea-endpoint-gibbs", "Automatic_fast"}, check_mdea_endpoint_gibbs},
        {{"native.equilibrium.discovery-work", "discovery-work", "Automatic_fast"}, check_discovery_work},
        {{"native.equilibrium.concurrent-first-discovery", "concurrent-first-discovery", "Automatic_fast"}, check_concurrent_first_discovery},
        {{"native.equilibrium.species-thermochemistry", "species-thermochemistry", "Automatic_fast"}, check_species_thermochemistry},
        {{"native.equilibrium.pressure-basis-quotient", "pressure-basis-quotient", "Automatic_fast"}, check_pressure_basis_quotient},
        {{"native.equilibrium.ideal-gas", "ideal-gas", "Automatic_fast"}, check_ideal_gas},
        {{"native.equilibrium.caloric-admission", "caloric-admission", "Automatic_fast"}, caloric_specification_checks::admission, "caloric-specification"},
        {{"native.equilibrium.caloric-actions", "caloric-actions", "Automatic_fast"}, caloric_specification_checks::actions, "caloric-specification"},
        {{"native.equilibrium.caloric-curvature", "caloric-curvature", "Automatic_fast"}, caloric_specification_checks::curvature, "caloric-specification"},
        {{"native.equilibrium.internal-energy-specification", "internal-energy-specification", "Automatic_fast"}, caloric_specification_checks::internal_energy_specification, "internal-energy-specification"},
        {{"native.equilibrium.volume-specification", "volume-specification", "Automatic_fast"}, caloric_specification_checks::volume_specification, "volume-specification"},
        {{"native.equilibrium.caloric-roundtrips", "caloric-roundtrips", "Automatic_fast"}, caloric_specification_checks::roundtrips, "caloric-specification"},
        {{"native.equilibrium.caloric-reactions", "caloric-reactions", "Automatic_fast"}, caloric_specification_checks::reactions, "caloric-specification"},
        {{"native.equilibrium.caloric-apparent-heat-capacity", "caloric-apparent-heat-capacity", "Automatic_fast"}, caloric_specification_checks::apparent_heat_capacity, "caloric-specification"},
        {{"native.equilibrium.diffusive-spinodal", "diffusive-spinodal", "Automatic_fast"}, critical_envelope_checks::spinodal, "critical-envelope"},
        {{"native.equilibrium.critical-pure", "critical-pure", "Automatic_fast"}, critical_envelope_checks::pure_critical, "critical-envelope"},
        {{"native.equilibrium.critical-binary", "critical-binary", "Automatic_fast"}, critical_envelope_checks::binary_critical, "critical-envelope"},
        {{"native.equilibrium.critical-curvature", "critical-curvature", "Automatic_fast"}, critical_envelope_checks::curvature, "critical-envelope"},
        {{"native.equilibrium.envelope-mixture", "envelope-mixture", "Explicit"}, critical_envelope_checks::mixture, "critical-envelope"},
        {{"native.equilibrium.envelope-pure", "envelope-pure", "Explicit"}, critical_envelope_checks::pure, "critical-envelope"},
        {{"native.equilibrium.azeotrope", "azeotrope", "Automatic_fast"}, critical_envelope_checks::azeotrope, "critical-envelope"},
        {{"native.equilibrium.removal", "TPD", "Automatic_fast"}, check_phase_removal_inventory},
        {{"native.equilibrium.density-correlated-ion-suppression", "DensityCorrelatedIonSuppression", "Automatic_fast"}, check_density_correlated_ion_suppression_admission},
        {{"native.equilibrium.rank-probe-permittivity-denominator", "RankProbePermittivityDenominator", "Automatic_fast"}, check_rank_probe_permittivity_denominator},
        {{"native.equilibrium.ion-suppression-reacting-actions", "IonSuppressionReactingActions", "Automatic_fast"}, check_ion_suppression_reacting_liquid_actions},
        {{"native.equilibrium.plane-curvature", "TPD", "Automatic_fast"}, check_tangent_plane_curvature},
        {{"native.equilibrium.exact-zero-trace-dot", "exact-zero-trace-dot", "Automatic_fast"}, check_exact_zero_trace_dot},
        {{"native.equilibrium.feed-charge-permutation", "feed-charge-permutation", "Automatic_fast"}, check_feed_charge_permutation},
        {{"native.equilibrium.liquid-trial-charge-permutation", "liquid-trial-charge-permutation", "Automatic_fast"}, check_liquid_trial_charge_permutation},
        {{"native.equilibrium.prepared-supplied-trials", "prepared-supplied-trials", "Automatic_fast"}, check_prepared_supplied_trials},
        {{"native.equilibrium.plane-admission", "TPD", "Automatic_fast"}, check_tangent_plane_admission},
        {{"native.equilibrium.trace-ion-curvature", "trace-ion-curvature", "Automatic_fast"}, check_trace_ion_curvature},
        {{"native.equilibrium.candidate-publication", "candidate-publication", "Automatic_fast"}, check_candidate_publication},
        {{"native.equilibrium.action-scaled-jacobian-overflow", "action-scaled-jacobian-overflow", "Automatic_fast"}, check_action_scaled_jacobian_overflow},
        {{"native.equilibrium.budget-retry", "budget-retry", "Automatic_fast"}, check_budget_retry},
        {{"native.equilibrium.associating-isomer-gibbs", "associating-isomer-gibbs", "Automatic_fast"}, check_associating_isomer_gibbs},
        {{"native.equilibrium.reaction-relaxed-split", "reaction-relaxed-split", "Automatic_fast"}, check_reaction_relaxed_split},
        {{"native.equilibrium.mdea-full-stability", "mdea-full-stability", "Automatic_fast"}, check_mdea_full_stability},
        {{"native.equilibrium.piperazine-a12s-bubble-stability", "piperazine-a12s-bubble-stability", "Automatic_fast"}, check_piperazine_a12s_bubble_stability},
        {{"native.equilibrium.mdea-bubble-work-limit", "mdea-bubble-work-limit", "Automatic_fast"}, check_mdea_bubble_work_limit},
        {{"native.equilibrium.mdea-bubble-recovery", "mdea-bubble-recovery", "Automatic_fast"}, check_mdea_bubble_recovery},
        {{"native.equilibrium.mdea-bubble-refusal-work", "mdea-bubble-refusal-work", "Automatic_fast"}, check_mdea_bubble_refusal_work},
        {{"native.equilibrium.bubble-branch-end", "bubble-branch-end", "Automatic_fast"}, check_bubble_branch_end},
        {{"native.equilibrium.reaction-relaxed-search", "reaction-relaxed-search", "Automatic_fast"}, check_reaction_relaxed_search},
        {{"native.equilibrium.reaction-relaxed-pz-row", "reaction-relaxed-pz-row", "Automatic_fast"}, check_reaction_relaxed_pz_row},
        {{"native.equilibrium.mea-census-state-3", "mea-census-state-3", "Automatic_fast"}, check_mea_census_state_3},
        {{"native.equilibrium.mea-trace-ion-discovery", "mea-trace-ion-discovery", "Automatic_fast"}, check_mea_trace_ion_discovery},
        {{"native.equilibrium.reference-parameter-coverage", "reference-parameter-coverage", "Automatic_fast"}, check_reference_parameter_coverage},
        {{"native.equilibrium.mea-discovery-omitted-species", "mea-discovery-omitted-species", "Automatic_fast"}, check_mea_discovery_omitted_species},
        {{"native.equilibrium.prescribed-liquid-coincidence", "prescribed-liquid-coincidence", "Automatic_fast"}, check_prescribed_liquid_coincidence},
        {{"native.equilibrium.electric-potential-difference", "electric-potential-difference", "Automatic_fast"}, check_electric_potential_difference},
        {{"native.equilibrium.electric-potential-difference-phase-order", "electric-potential-difference-phase-order", "Automatic_fast"}, check_electric_potential_difference_phase_order},
        {{"native.equilibrium.vanished-water-liquid-seed", "vanished-water-liquid-seed", "Automatic_fast"}, check_vanished_water_liquid_seed},
        {{"native.equilibrium.mea-below-bubble-hold", "mea-below-bubble-hold", "Automatic_fast"}, check_mea_below_bubble_hold},
        {{"native.equilibrium.mea-packet3-trace-root", "mea-packet3-trace-root", "Automatic_fast"}, check_mea_packet3_trace_root},
        {{"native.equilibrium.stability-grid", "stability-grid", "Automatic_fast"}, check_stability_grid},
        {{"native.equilibrium.condensed-trial-recovery", "condensed-trial-recovery", "Automatic_fast"}, check_condensed_trial_recovery},
        {{"native.equilibrium.physical-amount-work", "physical-amount-work", "Automatic_fast"}, check_physical_amount_work},
        {{"native.equilibrium.solve-error-curvature", "solve-error-curvature", "Automatic_fast"}, check_solve_error_curvature},
        {{"native.equilibrium.stability-inventory-face", "stability-inventory-face", "Automatic_fast"}, check_stability_inventory_face},
        {{"native.equilibrium.fixed-state-plane-uncertainty", "fixed-state-plane-uncertainty", "Automatic_fast"}, check_fixed_state_plane_uncertainty},
        {{"native.equilibrium.stability-phase-comparison", "stability-phase-comparison", "Automatic_fast"}, check_stability_phase_comparison},
        {{"native.equilibrium.stability-amount-correction", "stability-amount-correction", "Automatic_fast"}, check_stability_amount_correction},
        {{"native.equilibrium.trial-refinement", "trial-refinement", "Automatic_fast"}, check_trial_refinement},
        {{"native.equilibrium.normalized-reaction", "TPD", "Automatic_fast"}, check_normalized_reaction_affinity},
        {{"native.equilibrium.reacting-liquid-split", "TPD", "Automatic_fast"}, check_reacting_liquid_split},
        {{"native.equilibrium.declared-reacting-liquid-split", "TPD", "Automatic_fast"}, check_declared_reacting_liquid_split},
        {{"native.equilibrium.gross-2002-bubble", "TPD", "Automatic_fast"}, check_gross_2002_bubble},
        {{"native.equilibrium.gross-2002-cooling-onset", "TPD", "Automatic_fast"}, check_gross_2002_cooling_onset},
        {{"native.equilibrium.density-branch-admission", "TPD", "Automatic_fast"}, check_density_branch_admission},
        {{"native.equilibrium.tangent-plane-jacobian", "TPD", "Automatic_fast"}, check_tangent_plane_jacobian},
        {{"native.equilibrium.f5", "F5", "Automatic_fast"}, check_f5},
        {{"native.equilibrium.ideal-molarity", "ideal-molarity", "Automatic_fast"}, check_ideal_molarity},
        {{"native.equilibrium.ideal-mole-fraction", "ideal-mole-fraction", "Automatic_fast"}, check_ideal_mole_fraction},
        {{"native.equilibrium.s1", "S1", "Automatic_fast"}, check_s1},
        {{"native.equilibrium.reference-caloric-identities", "reference-caloric-identities", "Automatic_fast"}, reactive_caloric_checks::cheap, "S2"},
        {{"native.equilibrium.reactive-caloric-F2", "reactive-caloric-F2", "Explicit"}, reactive_caloric_checks::cp_f2, "reactive-cp"},
        {{"native.equilibrium.reactive-cp-F1", "reactive-cp-F1", "Explicit"}, [] {for(bool ionic:{false,true}) for(bool fixed:{false,true}) reactive_caloric_checks::cp_f1(ionic,fixed);}, "reactive-cp"},
        {{"native.equilibrium.reactive-cp-F3", "reactive-cp-F3", "Explicit"}, reactive_caloric_checks::cp_f3, "reactive-cp"},
        {{"native.equilibrium.reactive-caloric-ph-ps", "reactive-caloric-ph-ps", "Explicit"}, reactive_caloric_checks::caloric_ph_ps::run, "caloric-specification"},
#if defined(__linux__)
        {{"native.equilibrium.reactive-caloric-E-sys", "reactive-caloric-E-sys", "Explicit"}, [] {reactive_caloric_checks::run(false,false);}, "reactive-caloric"},
        {{"native.equilibrium.reactive-caloric-E-fix", "reactive-caloric-E-fix", "Explicit"}, [] {reactive_caloric_checks::run(false,true);}, "reactive-caloric"},
        {{"native.equilibrium.reactive-caloric-I-sys", "reactive-caloric-I-sys", "Explicit"}, [] {reactive_caloric_checks::run(true,false);}, "reactive-caloric"},
        {{"native.equilibrium.reactive-caloric-I-fix", "reactive-caloric-I-fix", "Explicit"}, [] {reactive_caloric_checks::run(true,true);}, "reactive-caloric"},
        {{"native.equilibrium.reactive-caloric-U", "reactive-caloric-U", "Explicit"}, reactive_caloric_checks::unreactive, "reactive-caloric"},
        {{"native.equilibrium.reactive-caloric-F6-E-T-sys", "reactive-caloric-F6-E-T-sys", "Explicit"}, [] {reactive_caloric_checks::boundary(false,true,false);}, "reactive-caloric"},
        {{"native.equilibrium.reactive-caloric-F6-E-T-fix", "reactive-caloric-F6-E-T-fix", "Explicit"}, [] {reactive_caloric_checks::boundary(false,true,true);}, "reactive-caloric"},
        {{"native.equilibrium.reactive-caloric-F6-E-P-sys", "reactive-caloric-F6-E-P-sys", "Explicit"}, [] {reactive_caloric_checks::boundary(false,false,false);}, "reactive-caloric"},
        {{"native.equilibrium.reactive-caloric-F6-E-P-fix", "reactive-caloric-F6-E-P-fix", "Explicit"}, [] {reactive_caloric_checks::boundary(false,false,true);}, "reactive-caloric"},
        {{"native.equilibrium.reactive-caloric-F6-N-T-sys", "reactive-caloric-F6-N-T-sys", "Explicit"}, [] {reactive_caloric_checks::boundary(true,true,false);}, "reactive-caloric"},
        {{"native.equilibrium.reactive-caloric-F6-N-T-fix", "reactive-caloric-F6-N-T-fix", "Explicit"}, [] {reactive_caloric_checks::boundary(true,true,true);}, "reactive-caloric"},
        {{"native.equilibrium.reactive-caloric-F6-N-P-sys", "reactive-caloric-F6-N-P-sys", "Explicit"}, [] {reactive_caloric_checks::boundary(true,false,false);}, "reactive-caloric"},
        {{"native.equilibrium.reactive-caloric-F6-N-P-fix", "reactive-caloric-F6-N-P-fix", "Explicit"}, [] {reactive_caloric_checks::boundary(true,false,true);}, "reactive-caloric"},
#endif
        {{"native.equilibrium.reaction-quotient-pressure", "reaction-quotient-pressure", "Automatic_fast"}, check_reaction_quotient_pressure},
        {{"native.equilibrium.reaction-quotient-actions", "reaction-quotient-actions", "Automatic_fast"}, check_reaction_quotient_actions},
        {{"native.equilibrium.neutral-reference-permutations", "neutral-reference-permutations", "Automatic_fast"}, check_neutral_reference_permutations},
        {{"native.equilibrium.mixed-reference-work", "mixed-reference-work", "Automatic_fast"}, check_mixed_reference_work},
        {{"native.equilibrium.s2", "S2", "Automatic_fast"}, check_s2},
        {{"native.equilibrium.trace-reaction-balance", "trace-reaction-balance", "Automatic_fast"}, check_trace_reaction_balance},
        {{"native.equilibrium.s6", "S6", "Automatic_fast"}, check_s6},
        {{"native.equilibrium.s7", "S7", "Automatic_fast"}, check_s7},
        {{"native.equilibrium.c1-multipliers", "C1", "Automatic_fast"}, check_c1_multiplier_scale},
        {{"native.equilibrium.c1-ideal-caller", "C1", "Automatic_fast"}, check_c1_ideal_caller_chart},
        {{"native.equilibrium.s11", "S11", "Automatic_fast"}, check_s11},
        {{"native.equilibrium.s13", "S13", "Automatic_fast"}, check_s13},
        {{"native.equilibrium.c3", "C3", "Automatic_fast"}, check_c3_availability},
        {{"native.equilibrium.s16", "S16", "Automatic_fast"}, check_s16},
        {{"native.equilibrium.e2-thermal", "E2-2", "Automatic_fast"}, check_e2_thermal},
        {{"native.equilibrium.e2-solve", "E2-6", "Automatic_fast"}, check_e2_solve},
        {{"native.equilibrium.e2-reuse", "E2-R", "Automatic_fast"}, check_e2_reuse},
        {{"native.equilibrium.e5-neutral", "E5-1", "Automatic_fast"}, check_e5_neutral},
        {{"native.equilibrium.e5-pins", "E5-2", "Automatic_fast"}, check_e5_pins},
        {{"native.equilibrium.e5-reactive", "E5-3", "Automatic_fast"}, check_e5_reactive},
        {{"native.equilibrium.e5-associating", "E5-4", "Automatic_fast"}, check_e5_associating},
        {{"native.equilibrium.e5-electrolyte", "E5-5", "Automatic_fast"}, check_e5_electrolyte},
        {{"native.equilibrium.e5-roots", "E5-6", "Automatic_fast"}, check_e5_roots},
        {{"native.equilibrium.e5-thermal", "E5-7", "Automatic_fast"}, check_e5_thermal},
        {{"native.equilibrium.e5-refusals", "E5-8", "Automatic_fast"}, check_e5_refusals},
        {{"native.equilibrium.e5-identities", "E5-10", "Automatic_fast"}, check_e5_identities},
    });
}
