"""Run from the checkout root: .venv/bin/python docs/site/examples/idaes-vle-flash.py."""
import os
from pathlib import Path
import idaes
from idaes.core import EnergyBalanceType, FlowsheetBlock
from idaes.models.properties.modular_properties.base.generic_property import GenericParameterBlock
from idaes.models.unit_models import Flash
from pyomo.environ import ConcreteModel, SolverFactory, value
from epcsaft import Mixture, Parameters, equilibrium as eq
from epcsaft.idaes import asl_library, check_state, property_configuration, solve_with_diagnostics

root = Path.cwd()
out = root / "build/docs/idaes-vle"
out.mkdir(parents=True, exist_ok=True)
model = Mixture(Parameters.from_json(
    root / "data/packets/gross-2001-methane-ethane/3/parameters/parameters.json",
    components=["methane", "ethane"]))
feed = eq.Amounts({"methane": 0.5, "ethane": 0.5})
boundaries = []
for zero, guess in (("vapor", 170), ("liquid", 210)):
    boundaries.append(eq.solve_equilibrium(model, eq.Problem(
        T=eq.Free(guess), P=1e6, feed=feed, row_tolerance=1e-11,
        phases=[eq.Phase(p, kind=p, amount=eq.Pinned(0) if p == zero else eq.Free())
                for p in ("liquid", "vapor")])))
m = ConcreteModel()
m.fs = FlowsheetBlock(dynamic=False)
m.fs.properties = GenericParameterBlock(**property_configuration(
    model, phases=("liquid", "vapor"), temperature_ref=200, pressure_ref=1e5,
    state_bounds={"flow_mol": (0, 1, 2), "temperature": (140, 180, 260), "pressure": (1e5, 1e6, 5e6)}))
m.fs.flash = Flash(property_package=m.fs.properties,
                   energy_balance_type=EnergyBalanceType.none, has_pressure_change=False)
f = m.fs.flash
f.inlet.flow_mol.fix(1)
f.inlet.temperature.fix(180)
f.inlet.pressure.fix(1e6)
for j in model.component_ids:
    f.inlet.mole_frac_comp[0, j].fix(0.5)
state = f.control_volume.properties_out[0]
state.temperature.fix(180)
for s in (f.control_volume.properties_in[0], state):
    s.temperature_bubble["vapor", "liquid"].set_value(boundaries[0].temperature)
    s.temperature_dew["vapor", "liquid"].set_value(boundaries[1].temperature)
    s._t1_vapor_liquid.set_value(180)
    s._teq["vapor", "liquid"].set_value(180)
    s.eps_1_vapor_liquid.set_value(1e-6)
    s.eps_2_vapor_liquid.set_value(5e-8)
    for k, j in enumerate(model.component_ids):
        s._mole_frac_tbub["vapor", "liquid", j].set_value(boundaries[0].mole_fractions[2+k])
        s._mole_frac_tdew["vapor", "liquid", j].set_value(boundaries[1].mole_fractions[k])
        s.mole_frac_phase_comp["liquid", j].set_value((0.3, 0.7)[k])
        s.mole_frac_phase_comp["vapor", j].set_value((0.92, 0.08)[k])
    for p, amount in (("liquid", 0.67), ("vapor", 0.33)):
        s.flow_mol_phase[p].set_value(amount)
        s.phase_frac[p].set_value(amount)
options = dict(tol=1e-11, constr_viol_tol=1e-11, acceptable_iter=0, linear_solver="mumps",
               hessian_approximation="exact", print_user_options="yes", print_level=5)
option_file = out / "ipopt.opt"
option_file.write_text("".join(f"{key} {v}\n" for key, v in options.items()))
os.environ["AMPLFUNC"] = "\n".join(filter(None, (os.environ.get("AMPLFUNC"), str(asl_library()))))
solver = SolverFactory("ipopt", executable=str(Path(idaes.bin_directory) / "ipopt"))
solver.options["option_file_name"] = str(option_file.resolve())
solve_with_diagnostics(solver, m, logfile=out / "ipopt.log")
reference = eq.solve_equilibrium(model, eq.Problem(T=180, P=1e6, feed=feed,
    phases=[eq.Phase("liquid", kind="liquid"), eq.Phase("vapor", kind="vapor")], row_tolerance=1e-11))
for index, p in enumerate(state.phase_list):
    print(p, "mol/s:", value(state.flow_mol_phase[p]), "composition:",
          [value(state.mole_frac_phase_comp[p, j]) for j in model.component_ids],
          "density mol/m3:", check_state(state, p).molar_density,
          "native fraction:", reference.amounts[index])
print("component residuals mol/s:", [value(c.body) for c in state.component_flow_balances.values()])
print("equilibrium residuals:", [value(c.body) for c in state.equilibrium_constraint.values()])
