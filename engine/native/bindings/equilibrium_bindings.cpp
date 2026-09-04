#include "equilibrium/problem.hpp"
#include "equilibrium/gibbs_minimizer.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <optional>

namespace py = pybind11;
namespace epcsaft::equilibrium {

void register_equilibrium(py::module_& module) {
    py::enum_<PhaseKind>(module, "PhaseKind")
        .value("Liquid", PhaseKind::Liquid)
        .value("Vapor", PhaseKind::Vapor);
    py::enum_<AmountKind>(module, "AmountKind")
        .value("Free", AmountKind::Free)
        .value("Pinned", AmountKind::Pinned);
    py::enum_<FeedKind>(module, "FeedKind")
        .value("MoleFractions", FeedKind::MoleFractions)
        .value("Amounts", FeedKind::Amounts)
        .value("MassFractions", FeedKind::MassFractions)
        .value("Molalities", FeedKind::Molalities);
    py::class_<AmountSpec>(module, "AmountSpec")
        .def(py::init<>())
        .def_readwrite("kind", &AmountSpec::kind)
        .def_readwrite("fraction", &AmountSpec::fraction)
        .def_static("free", &AmountSpec::free)
        .def_static("pinned", &AmountSpec::pinned);
    py::class_<FreeValue>(module, "FreeValue")
        .def(py::init<>())
        .def_readwrite("free", &FreeValue::free)
        .def_readwrite("value", &FreeValue::value)
        .def_static("fixed", &FreeValue::fixed)
        .def_static("variable", [](std::optional<double> guess) {
            return guess.has_value() ? FreeValue::variable(*guess) : FreeValue::variable();
        }, py::arg("guess") = py::none());
    py::class_<PhaseSpec>(module, "PhaseSpec")
        .def(py::init<>())
        .def_readwrite("name", &PhaseSpec::name)
        .def_readwrite("kind", &PhaseSpec::kind)
        .def_readwrite("support", &PhaseSpec::support)
        .def_readwrite("amount", &PhaseSpec::amount)
        .def_readwrite("composition_guess", &PhaseSpec::composition_guess)
        .def_readwrite("packing_guess", &PhaseSpec::packing_guess);
    py::class_<FeedSpec>(module, "FeedSpec")
        .def(py::init<>())
        .def_readwrite("kind", &FeedSpec::kind)
        .def_readwrite("values", &FeedSpec::values)
        .def_readwrite("solvent", &FeedSpec::solvent)
        .def_readwrite("solvent_mass_kg", &FeedSpec::solvent_mass_kg);
    py::class_<ReactionSpec>(module, "ReactionSpec")
        .def(py::init<>())
        .def_readwrite("name", &ReactionSpec::name)
        .def_readwrite("stoichiometry", &ReactionSpec::stoichiometry)
        .def_readwrite("equilibrium_constant", &ReactionSpec::equilibrium_constant)
        .def_readwrite("affinity_j_per_mol", &ReactionSpec::affinity_j_per_mol);
    py::class_<EquilibriumProblem>(module, "EquilibriumProblem")
        .def(py::init<>())
        .def_readwrite("phases", &EquilibriumProblem::phases)
        .def_readwrite("temperature", &EquilibriumProblem::temperature)
        .def_readwrite("pressure", &EquilibriumProblem::pressure)
        .def_readwrite("feed", &EquilibriumProblem::feed)
        .def_readwrite("reactions", &EquilibriumProblem::reactions)
        .def_readwrite("standard_gibbs_j_per_mol", &EquilibriumProblem::standard_gibbs_j_per_mol)
        .def_readwrite("topology_declared", &EquilibriumProblem::topology_declared);
    py::enum_<RowKind>(module, "RowKind")
        .value("Chemical", RowKind::Chemical)
        .value("Pressure", RowKind::Pressure)
        .value("Moiety", RowKind::Moiety)
        .value("Electro", RowKind::Electro)
        .value("Normalization", RowKind::Normalization);
    py::class_<Row>(module, "Row")
        .def_readonly("kind", &Row::kind)
        .def_readonly("phase", &Row::phase)
        .def_readonly("index", &Row::index)
        .def_readonly("eqid", &Row::eqid);
    py::class_<PointEvaluation>(module, "PointEvaluation")
        .def_readonly("values", &PointEvaluation::values)
        .def_readonly("jacobian", &PointEvaluation::jacobian)
        .def_readonly("objective", &PointEvaluation::objective)
        .def_readonly("floor", &PointEvaluation::floor);
    py::class_<CompiledSystem>(module, "CompiledSystem")
        .def_property_readonly("rows", [](const CompiledSystem& system) { return system.rows; })
        .def_readonly("pinned_count", &CompiledSystem::pinned_count)
        .def_readonly("free_intensive_count", &CompiledSystem::free_intensive_count)
        .def_readonly("rank", &CompiledSystem::rank)
        .def_readonly("charge_row_count", &CompiledSystem::charge_row_count)
        .def_readonly("primal_size", &CompiledSystem::primal_size)
        .def_readonly("rows_minus_unknowns", &CompiledSystem::rows_minus_unknowns)
        .def_readonly("generic_nullity", &CompiledSystem::generic_nullity)
        .def_readonly("feed_total", &CompiledSystem::feed_total)
        .def_readonly("has_ions", &CompiledSystem::has_ions)
        .def_readonly("has_reaction", &CompiledSystem::has_reaction)
        .def_readonly("jacobian_rank", &CompiledSystem::jacobian_rank)
        .def_readonly("removed_species", &CompiledSystem::removed_species)
        .def("generic_point", &CompiledSystem::generic_point)
        .def("evaluate", &CompiledSystem::evaluate)
        .def("jacobian", &CompiledSystem::jacobian)
        .def("lever_jacobian", &CompiledSystem::lever_jacobian)
        .def("objective_gradient", &CompiledSystem::objective_gradient)
        .def("hessian", &CompiledSystem::hessian);
    py::class_<EquilibriumResult>(module, "EquilibriumResult")
        .def_readonly("success", &EquilibriumResult::success)
        .def_readonly("message", &EquilibriumResult::message)
        .def_readonly("temperature", &EquilibriumResult::temperature)
        .def_readonly("pressure", &EquilibriumResult::pressure)
        .def_readonly("amounts", &EquilibriumResult::amounts)
        .def_readonly("mole_fractions", &EquilibriumResult::mole_fractions)
        .def_readonly("molar_densities", &EquilibriumResult::molar_densities)
        .def_readonly("packing_fractions", &EquilibriumResult::packing_fractions)
        .def_readonly("multipliers", &EquilibriumResult::multipliers)
        .def_readonly("residuals", &EquilibriumResult::residuals)
        .def_readonly("sensitivity", &EquilibriumResult::sensitivity)
        .def_readonly("coordinates", &EquilibriumResult::coordinates)
        .def_readonly("phase_count", &EquilibriumResult::phase_count)
        .def_readonly("iterations", &EquilibriumResult::iterations)
        .def_readonly("evaluations", &EquilibriumResult::evaluations)
        .def_readonly("hessian_calls", &EquilibriumResult::hessian_calls)
        .def_readonly("objective", &EquilibriumResult::objective)
        .def_readonly("reduced_hessian_minimum", &EquilibriumResult::reduced_hessian_minimum)
        .def_readonly("convergence", &EquilibriumResult::convergence)
        .def_readonly("rows", &EquilibriumResult::rows);
    py::class_<SensitivityInfo>(module, "SensitivityInfo")
        .def_readonly("equations", &SensitivityInfo::equations)
        .def_readonly("parameters", &SensitivityInfo::parameters)
        .def_readonly("rank", &SensitivityInfo::rank)
        .def_readonly("condition_estimate", &SensitivityInfo::condition_estimate)
        .def_readonly("nonsingular", &SensitivityInfo::nonsingular);
    module.def("compile_problem", &compile_problem, py::keep_alive<0, 1>());
    module.def("sensitivity", &sensitivity);
    module.def("solve_equilibrium", &solve_equilibrium, py::arg("model"), py::arg("problem"),
               py::arg("start") = std::nullopt);
}

}  // namespace epcsaft::equilibrium
