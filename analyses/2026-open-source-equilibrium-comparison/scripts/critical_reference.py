"""Issue 224 B: FeOs v0.10.1 / teqp v0.23.2 calls on the declared Gross (2001) inputs; C++ owns every comparison."""
import csv, json, sys
from pathlib import Path
import numpy as np, teqp
from feos import EquationOfState, Identifier, Parameters, PhaseEquilibrium, PureRecord, State, set_num_threads
from si_units import KELVIN, METER, MOL, PASCAL
set_num_threads(1)
root = Path(sys.argv[1])
spec = json.loads((root / 'reference-inputs.json').read_text())
rho_unit = MOL / METER**3
def components(names): return [c for c in spec['components'] if c['name'] in names]
def feos_eos(names):
    return EquationOfState.pcsaft(Parameters.from_records([PureRecord(Identifier(name=c['name']), c['mass'], **c['residual']) for c in components(names)]))
def teqp_model(names):
    return teqp.make_model({'kind': 'PCSAFT', 'model': {'coeffs': [{'name': c['name'], 'm': c['residual']['m'], 'sigma_Angstrom': c['residual']['sigma'],
        'epsilon_over_k': c['residual']['epsilon_k'], 'BibTeXKey': 'Gross2001'} for c in components(names)]}})
def table(name, header, rows):
    with (root / name).open('w') as f:
        writer = csv.writer(f, lineterminator='\n'); writer.writerow(header); writer.writerows(rows)
def teqp_pressure(model, T, rhovec): return float(rhovec.sum() * model.get_R(rhovec / rhovec.sum()) * T + model.get_pr(T, rhovec))
c = State.critical_point(feos_eos(['propane']), initial_temperature=375 * KELVIN, max_iter=200, tol=1e-12)
m = teqp_model(['propane']); T, r = m.solve_pure_critical(375.0, 5000.0)
table('critical-pure.csv', ['source', 'z0', 'T_K', 'rho_mol_m3', 'P_Pa'],
      [['feos', 1, c.temperature / KELVIN, c.density / rho_unit, c.pressure() / PASCAL], ['teqp', 1, T, r, teqp_pressure(m, T, np.array([r]))]])
binary, rows = feos_eos(['methane', 'ethane']), []
for z0 in spec['feos_critical_z0']:
    c = State.critical_point(binary, np.array([z0, 1 - z0]), initial_temperature=260 * KELVIN, max_iter=200, tol=1e-12)
    rows.append(['feos', z0, c.temperature / KELVIN, c.density / rho_unit, c.pressure() / PASCAL])
mb, options = teqp_model(['methane', 'ethane']), teqp.TCABOptions()
for key, value in spec['teqp_trace_options'].items(): setattr(options, key, value)
T0, r0 = teqp_model(['methane']).solve_pure_critical(190.0, 10000.0)
locus = [(p['rho0 / mol/m^3'] / (p['rho0 / mol/m^3'] + p['rho1 / mol/m^3']), p) for p in mb.trace_critical_arclength_binary(T0, np.array([r0, 0.0]), None, options)]
for target in spec['teqp_locus_near_z0']:  # teqp has no fixed-composition critical solve: its own nearest locus point
    z0, p = min(locus, key=lambda item: abs(item[0] - target))
    rhovec = np.array([p['rho0 / mol/m^3'], p['rho1 / mol/m^3']])
    rows.append(['teqp', z0, p['T / K'], rhovec.sum(), teqp_pressure(mb, p['T / K'], rhovec)])
table('critical-binary.csv', ['source', 'z0', 'T_K', 'rho_mol_m3', 'P_Pa'], rows)
points = []
for q in spec['envelope_points']:  # Engine P and incipient x seed the branch, e.g. one of two retrograde dew points
    seed, solve = dict(tp_init=q['P'] * PASCAL, max_iter_inner=1000, max_iter_outer=1000, tol_inner=1e-11, tol_outer=1e-11), np.array([q['x0'], 1 - q['x0']])
    if q['bubble']:
        e = PhaseEquilibrium.bubble_point(binary, q['T'] * KELVIN, np.array([0.5, 0.5]), vapor_molefracs=solve, **seed)
        feed, incipient = e.liquid, e.vapor
    else:
        e = PhaseEquilibrium.dew_point(binary, q['T'] * KELVIN, np.array([0.5, 0.5]), liquid_molefracs=solve, **seed)
        feed, incipient = e.vapor, e.liquid
    points.append([int(q['bubble']), q['T'], feed.pressure() / PASCAL, feed.density / rho_unit, incipient.density / rho_unit, incipient.molefracs[0]])
table('envelope-points.csv', ['bubble', 'T_K', 'P_Pa', 'rho_feed', 'rho_incipient', 'x0_incipient'], points)
