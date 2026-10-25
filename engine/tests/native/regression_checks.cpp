#include "recording_support.hpp"

// The case registry. Each case function is defined in one regression_*_checks.cpp unit.
namespace regression_checks {
void check_propane_inventory_boundary();
void check_inventory_face();
void check_neutral_electrostatics();
void check_coupled_columns();
void check_coupled_recovery();
void check_hydration_failures();
void check_ionic_addition_enthalpy();
void check_free_temperature_addition();
void check_addition_enthalpy();
void check_log_constant();
void check_f1();
void check_f2();
void check_f3_f4();
void check_f5_f6();
void check_f7_f8();
void check_f9_f10();
void check_d86();
void check_j1();
void check_j2();
void check_minimax_consistency();
void check_d7_minimax();
void check_reaction_derivatives();
void check_reaction_recovery();
void check_reaction_rejections();
void check_reaction_enthalpy_shift();
void check_reactive_cp();
void check_cp_derivatives();
void check_cp_recovery();
void check_r3();
void check_measurement_columns();
void check_measurement_fit();
void check_measurement_reuse();
void check_measurement_admission();
void check_measurement_correlation();
}  // namespace regression_checks

int main(int argc, char** argv) {
    using namespace regression_checks;
    return recording_support::run(argc, argv, {
        {{"native.regression.measurement-columns", "measurement-columns"}, check_measurement_columns, "R2"},
        {{"native.regression.measurement-fit", "measurement-fit"}, check_measurement_fit, "R2"},
        {{"native.regression.measurement-reuse", "measurement-reuse"}, check_measurement_reuse, "R2"},
        {{"native.regression.measurement-admission", "measurement-admission"}, check_measurement_admission, "R2"},
        {{"native.regression.measurement-correlation", "measurement-correlation"}, check_measurement_correlation, "R2"},
        {{"native.regression.inventory-boundary", "inventory-boundary"}, check_propane_inventory_boundary, "R1"},
        {{"native.regression.inventory-face", "inventory-face"}, check_inventory_face, "R1"},
        {{"native.regression.neutral-electrostatics", "neutral-electrostatics"}, check_neutral_electrostatics, "R1"},
        {{"native.regression.affine-columns", "affine-columns"}, check_coupled_columns, "R2"},
        {{"native.regression.affine-recovery", "affine-recovery"}, check_coupled_recovery, "R2"},
        {{"native.regression.hydration-failures", "hydration-failures"}, check_hydration_failures, "R1"},
        {{"native.regression.ionic-addition", "ionic-addition"}, check_ionic_addition_enthalpy, "R2"},
        {{"native.regression.reactive-caloric-F6-addition", "reactive-caloric-F6-addition"}, check_free_temperature_addition, "R2"},
        {{"native.regression.addition-enthalpy", "addition-enthalpy"}, check_addition_enthalpy, "R2"},
        {{"native.regression.ln-k", "ln-k"}, check_log_constant, "R2"},
        {{"native.regression.r1-f1", "R1-F1"}, check_f1, "R1"},
        {{"native.regression.r1-f2", "R1-F2"}, check_f2, "R1"},
        {{"native.regression.r1-f3", "R1-F3"}, check_f3_f4, "R1"},
        {{"native.regression.r1-f5", "R1-F5"}, check_f5_f6, "R1"},
        {{"native.regression.r1-f7", "R1-F7"}, check_f7_f8, "R1"},
        {{"native.regression.r1-f9", "R1-F9"}, check_f9_f10, "R1"},
        {{"native.regression.d86", "D86"}, check_d86, "R1"},
        {{"native.regression.r2-j1", "R2-J1"}, check_j1, "R2"},
        {{"native.regression.r2-j2", "R2-J2"}, check_j2, "R2"},
        {{"native.regression.r2-minimax-consistency", "minimax-consistency"}, check_minimax_consistency, "R2"},
        {{"native.regression.r2-d7-minimax-analytic", "R2-D7"}, check_d7_minimax, "R2"},
        {{"native.regression.reaction-derivatives", "reaction-derivatives"}, check_reaction_derivatives, "R2"},
        {{"native.regression.reaction-recovery", "reaction-recovery"}, check_reaction_recovery, "R2"},
        {{"native.regression.reaction-enthalpy-shift", "reaction-enthalpy-shift"}, check_reaction_enthalpy_shift, "R2"},
        {{"native.regression.reaction-rejections", "reaction-rejections"}, check_reaction_rejections, "R2"},
        {{"native.regression.reactive-cp", "reactive-cp"}, check_reactive_cp},
        {{"native.regression.cp-derivatives", "cp-derivatives"}, check_cp_derivatives, "R2"},
        {{"native.regression.cp-recovery", "cp-recovery"}, check_cp_recovery, "R2"},
        {{"native.regression.r3", "R3"}, check_r3, "R2"},
    });
}
