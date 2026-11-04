"""IDAES pressure properties and classified-root, neutral smooth VLE.

FTPx, LogBubbleDew and SmoothVLE own the process equations (Burgard et al.,
PSE 2018). EOS values and exact derivatives come from the native ASL adapter.
"""

from __future__ import annotations

from ctypes import byref
from pathlib import Path
from typing import Any
import warnings

from idaes.core import Component, LiquidPhase, VaporPhase  # type: ignore[import-untyped]
from idaes.models.properties.modular_properties.eos.eos_base import EoSBase  # type: ignore[import-untyped]
from idaes.models.properties.modular_properties.phase_equil import SmoothVLE  # type: ignore[import-untyped]
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import LogBubbleDew  # type: ignore[import-untyped]
from idaes.models.properties.modular_properties.phase_equil.forms import log_fugacity  # type: ignore[import-untyped]
from idaes.models.properties.modular_properties.state_definitions import FTPx  # type: ignore[import-untyped]
from pyomo.core.base.external import AMPLExternalFunction, _ARGLIST  # type: ignore[import-untyped]
from pyomo.environ import check_optimal_termination, log, value, units as u  # type: ignore[import-untyped]

from . import _core
from .eos import DensityBranch, PressureActionDirection, property_asl_record


class PropertyExternalFunction(AMPLExternalFunction):  # type: ignore[misc]
    """Read the native ASL error that Pyomo 6.10.1 leaves unchecked."""

    def _evaluate(self, args: Any, fixed: Any, fgh: int) -> Any:
        if self._so is None:
            self.load_library()
        arguments = _ARGLIST(args, fgh, fixed)
        result = self._known_functions[self._function][0](byref(arguments))
        if arguments.Errmsg:
            message = arguments.Errmsg.decode()
            if message[:1] in ("'", '"'):
                message = message[1:]
            raise RuntimeError(message)
        n = len(args)
        gradient = [float("nan")] * n if fgh >= 1 else None
        hessian = [float("nan")] * (n * (n + 1) // 2) if fgh >= 2 else None
        for j in range(n):
            rj = arguments.at[j]
            if rj < 0:
                continue
            if gradient is not None:
                gradient[j] = arguments.derivs[rj]
            if hessian is not None:
                for i in range(j + 1):
                    ri = arguments.at[i]
                    if ri >= 0:
                        hessian[i + j * (j + 1) // 2] = arguments.hes[ri + rj * (rj + 1) // 2]
        return result, gradient, hessian


def asl_library() -> Path:
    """Locate the adapter belonging to the loaded Engine extension."""
    folder = Path(_core.__file__).resolve().parent
    for name in ("libepcsaft_asl.so", "libepcsaft_asl.dylib", "epcsaft_asl.dll"):
        if (folder / name).is_file():
            return folder / name
    raise FileNotFoundError(f"ASL property library missing beside {_core.__file__}")


class PressureProperties(EoSBase):  # type: ignore[misc]
    """Generic-property EOS methods supplied by native pressure actions."""

    @staticmethod
    def build_parameters(phase: Any) -> None:
        """The immutable native record already contains the physical parameters."""

    @staticmethod
    def common(b: Any, phase: Any) -> None:
        options = phase.config.equation_of_state_options
        if b.component("_pressure_property") is None:
            for name, unit in (
                ("_pressure_property", u.dimensionless),
                ("_density_property", u.mol / u.m**3),
                ("_enthalpy_property", u.J / u.mol),
            ):
                b.add_component(
                    name,
                    PropertyExternalFunction(
                        library=options["library"],
                        function="epcsaft_pressure_property",
                        units=unit,
                        arg_units=[u.K, u.Pa, *[u.dimensionless for _ in b.component_list], None, None, None],
                    ),
                )

    @staticmethod
    def calculate_scaling_factors(b: Any, phase: Any) -> None:
        """Use generic-state scaling; log-fugacity equations have unit scale."""

    @staticmethod
    def _at(b: Any, p: str, row: int, T: Any, P: Any, x: Any, strict: bool) -> Any:
        phase = b.params.get_phase(p)
        options = phase.config.equation_of_state_options
        function = b._density_property if row == 0 else b._enthalpy_property if row == 1 else b._pressure_property
        return function(
            T,
            P,
            *[x[j] for j in b.component_list],
            options["equilibrium_record" if strict else "record"],
            str(row),
            "0" if phase.is_liquid_phase() else "1",
        )

    @staticmethod
    def dens_mol_phase(b: Any, p: str) -> Any:
        return PressureProperties._at(
            b, p, 0, b.temperature, b.pressure, {j: b.mole_frac_phase_comp[p, j] for j in b.component_list}, False
        )

    @staticmethod
    def enth_mol_phase(b: Any, p: str) -> Any:
        return PressureProperties._at(
            b, p, 1, b.temperature, b.pressure, {j: b.mole_frac_phase_comp[p, j] for j in b.component_list}, False
        )

    @staticmethod
    def _log_fugacity(b: Any, p: str, j: str, T: Any, P: Any, x: Any) -> Any:
        row = 2 + list(b.component_list).index(j)
        return log(x[j]) + log(P / b.params.pressure_ref) + PressureProperties._at(b, p, row, T, P, x, True)

    @staticmethod
    def log_fug_phase_comp_eq(b: Any, p: str, j: str, pp: Any) -> Any:
        return PressureProperties._log_fugacity(
            b, p, j, b._teq[pp], b.pressure, {k: b.mole_frac_phase_comp[p, k] for k in b.component_list}
        )

    @staticmethod
    def _boundary(b: Any, p: str, j: str, pp: Any, kind: str) -> Any:
        phase = b.params.get_phase(p)
        saturated = phase.is_liquid_phase() if kind.endswith("bub") else phase.is_vapor_phase()
        x = {
            k: b.mole_frac_comp[k] if saturated else b.component("_mole_frac_" + kind)[pp, k] for k in b.component_list
        }
        T = (
            b.component("temperature_bubble" if kind == "tbub" else "temperature_dew")[pp]
            if kind[0] == "t"
            else b.temperature
        )
        P = b.component("pressure_bubble" if kind == "pbub" else "pressure_dew")[pp] if kind[0] == "p" else b.pressure
        return PressureProperties._log_fugacity(b, p, j, T, P, x)

    @staticmethod
    def log_fug_phase_comp_Tbub(b: Any, p: str, j: str, pp: Any) -> Any:
        return PressureProperties._boundary(b, p, j, pp, "tbub")

    @staticmethod
    def log_fug_phase_comp_Tdew(b: Any, p: str, j: str, pp: Any) -> Any:
        return PressureProperties._boundary(b, p, j, pp, "tdew")

    @staticmethod
    def log_fug_phase_comp_Pbub(b: Any, p: str, j: str, pp: Any) -> Any:
        return PressureProperties._boundary(b, p, j, pp, "pbub")

    @staticmethod
    def log_fug_phase_comp_Pdew(b: Any, p: str, j: str, pp: Any) -> Any:
        return PressureProperties._boundary(b, p, j, pp, "pdew")


def property_configuration(
    model: Any, *, phases: tuple[str, ...], state_bounds: dict[str, Any], temperature_ref: float, pressure_ref: float
) -> dict[str, Any]:
    """Declare a neutral GenericParameterBlock; starts and state bounds are caller inputs.

    Phases are liquid, vapor, or both. Each VLE phase mixture density root must
    classify as its requested liquid or vapor side on its own (T, P, x) isotherm,
    or the call is refused. A component above its own critical temperature is
    admitted when both mixture roots classify; no pure-component critical test
    is applied. Generic initialization requiring pure saturation correlations
    is unavailable; supply state starts.
    """
    if not phases or len(set(phases)) != len(phases) or set(phases) - {"liquid", "vapor"}:
        raise ValueError("phases must select liquid, vapor, or both once")
    if any(e.magnitude != 0 for e in model.parameters._document.entries if e.field == "charge_number"):
        raise ValueError("IDAES pressure properties qualify neutral components only")
    record = property_asl_record(model)
    equilibrium_record = property_asl_record(model, require_phase_identity=True) if len(phases) == 2 else record
    options = {
        "library": str(asl_library()),
        "record": record,
        "equilibrium_record": equilibrium_record,
        "model": model,
    }
    pair = ("vapor", "liquid")
    config = dict(
        components={
            j: {"type": Component, **({"phase_equilibrium_form": {pair: log_fugacity}} if len(phases) == 2 else {})}
            for j in model.component_ids
        },
        phases={
            p: {
                "type": LiquidPhase if p == "liquid" else VaporPhase,
                "equation_of_state": PressureProperties,
                "equation_of_state_options": options,
            }
            for p in phases
        },
        base_units={"time": u.s, "length": u.m, "mass": u.kg, "amount": u.mol, "temperature": u.K},
        state_definition=FTPx,
        state_bounds=state_bounds,
        temperature_ref=temperature_ref,
        pressure_ref=pressure_ref,
    )
    if len(phases) == 2:
        config.update(
            phases_in_equilibrium=[pair], phase_equilibrium_state={pair: SmoothVLE}, bubble_dew_method=LogBubbleDew
        )
    return config


def check_state(state: Any, phase: str) -> Any:
    """Read native endpoint diagnostics; warn when the single root has no classified side."""
    pobj = state.params.get_phase(phase)
    model = pobj.config.equation_of_state_options["model"]
    batch = model.pressure_actions(
        value(state.temperature),
        value(state.pressure),
        x=[value(state.mole_frac_phase_comp[phase, j]) for j in state.component_list],
        phase=phase,
        observables=[],
        directions=[PressureActionDirection(composition=[0.0] * len(model.component_ids))],
    )
    if batch.molar_density is None:
        raise RuntimeError(batch.density_message)
    diagnostics = batch.diagnostics
    if diagnostics.returned_branch == DensityBranch.Unspecified:
        warnings.warn(
            f"Requested {phase}; returned branch Unspecified; stable_root_count={diagnostics.stable_root_count}; "
            "phase identity is not established by the finite search",
            RuntimeWarning,
            stacklevel=2,
        )
    elif diagnostics.returned_branch != (DensityBranch.Liquid if pobj.is_liquid_phase() else DensityBranch.Vapor):
        raise RuntimeError(f"branch_mismatch: requested {phase}, returned {diagnostics.returned_branch}")
    return batch


def solve_with_diagnostics(
    solver: Any, model: Any, *, logfile: Path, endpoints: tuple[tuple[Any, str], ...] | None = None
) -> Any:
    """Keep subprocess diagnostics in a caller-owned log and surface failure in Python."""
    logfile = Path(logfile).resolve()
    logfile.parent.mkdir(parents=True, exist_ok=True)
    solver._log = ""
    try:
        result = solver.solve(model, logfile=str(logfile), symbolic_solver_labels=True)
        if not check_optimal_termination(result):
            raise RuntimeError(str(result.solver))
    except Exception as error:
        if solver._log:
            logfile.write_text(solver._log)
        messages = (
            [line for line in logfile.read_text().splitlines() if "epcsaft_property:" in line]
            if logfile.exists()
            else []
        )
        raise RuntimeError(f"{error}; ASL diagnostics: {'; '.join(messages)}; solver log: {logfile}") from error
    for state, phase in endpoints or ():
        check_state(state, phase)
    return result
