"""#223 case 5: FeOs 0.10.1 reference refinement and native Engine comparison."""
import argparse
import csv
import json
import math
import struct
from pathlib import Path

import numpy as np
import baseline as b
import qualify as q
import costs

SOURCE = b.HERE/'inputs/feos-test-parameters.json'
SOURCE_SHA256 = 'bd09af1cb49e4b525f4c589d66b10579ac13600634b889f551de9d533f843580'


def engine_parameters():
    """Transport the FeOs two-site association and PC-SAFT quantities unchanged."""
    data = json.loads(SOURCE.read_text())
    selected = [next(row for row in data if row['identifier']['name'] == name) for name in ('water_np', 'hexane')]
    doc = b._core.ParameterDocument()
    doc.component_ids = ['water_np', 'hexane']
    doc.base_choice, doc.association_choice = 'pc-saft', 'general-site'
    doc.electrolyte_choice, doc.ionic_dispersion_choice, doc.permittivity_choice = 'none', '', 'none'
    doc.polar_choice = 'none'
    site = selected[0]['association_sites'][0]
    doc.site_component, doc.site_multiplicity = [0, 0], [int(site['na']), int(site['nb'])]
    rows = []
    for i, row in enumerate(selected):
        rows += [('molar_mass', [i], row['molarweight'], 'gram / mole'),
                 ('segment_count', [i], row['m'], 'dimensionless'),
                 ('sigma_constant', [i], row['sigma'], 'angstrom'),
                 ('dispersion_energy_over_k', [i], row['epsilon_k'], 'kelvin'),
                 ('charge_number', [i], 0., 'elementary-charge')]
    rows += [('k_ij', [0, 1], 0., 'dimensionless'),
             ('association_energy_over_k', [0, 1], site['epsilon_k_ab'], 'kelvin'),
             ('association_volume', [0, 1], site['kappa_ab'], 'dimensionless')]
    entries = []
    for field, index, value, unit in rows:
        entry = b._core.ParameterEntry()
        entry.identity, entry.field, entry.index, entry.magnitude, entry.unit = f'{field}/{index}', field, index, value, unit
        entries.append(entry)
    doc.entries = entries
    return b.Parameters(b._core.build_parameters(doc), tuple(doc.component_ids), (('water_np', 'A'), ('water_np', 'B')))


def feos_record(case, state, reference_plane, si, *, reduced_trial=False):
    x = np.asarray(state.molefracs, dtype=float)
    rho = float(state.density/(si.MOL/si.METER**3))
    mu = np.log(x)+np.asarray(state.ln_phi(), dtype=float)
    gradient = mu-reference_plane
    tpd = float(x@gradient)
    stationary = float(np.max(np.abs(gradient-tpd)))
    pressure = float(state.pressure()/si.PASCAL)
    values = dict(x=x.tolist(), density_mol_m3=rho,
                  pressure_Pa=pressure, pressure_relative_residual=abs(pressure/1e5-1),
                  tpd_over_rt=tpd)
    values['normalized_stationarity_max_abs' if reduced_trial else 'feed_plane_stationarity_max_abs'] = stationary
    if reduced_trial:
        # FeOs stability_analysis.rs uses Moles::from_reduced(y), not mol-valued y.
        y = np.asarray(state.moles*si.NAV, dtype=float)
        values['trial_coordinates_reduced'] = y.tolist()
        values['unnormalized_stationarity_max_abs'] = float(np.max(np.abs(np.log(y)+np.asarray(state.ln_phi(), dtype=float)-reference_plane)))
    q.quantities(case, values)
    return values



def retained_reference(path):
    metadata = json.loads(path.read_text())
    q.require(metadata['failure'] is None and metadata['feos_version'] == '0.10.1', 'FeOs reference is incomplete or has the wrong version')
    q.require(metadata['sources'] == json.loads((b.HERE/'inputs/feos-sources.json').read_text()), 'FeOs reference inputs/options changed')
    source = path.with_suffix('.csv')
    groups = {}
    for row in csv.DictReader(source.open()):
        if not row['case'].startswith('feos/refinement-') or '/trial-' not in row['case'] and '/flash-' not in row['case']:
            continue
        value = float(row['value'])
        q.require(struct.pack('>d', value).hex() == row['binary64'], 'FeOs retained binary64 value differs')
        group = groups.setdefault(row['case'], {})
        quantity = row['quantity']
        if '/' in quantity:
            field, index = quantity.rsplit('/', 1)
            group.setdefault(field, {})[int(index)] = value
        else:
            group[quantity] = value
    refs = []
    for run in (0, 1):
        reference = {}
        for kind, label in (('trials', 'trial'), ('phases', 'flash')):
            rows = [row for case, row in groups.items() if case.startswith(f'feos/refinement-{run}/{label}-')]
            for row in rows:
                for field, value in tuple(row.items()):
                    if isinstance(value, dict):
                        q.require(sorted(value) == list(range(len(value))), 'FeOs reference coordinate indices are incomplete')
                        row[field] = [value[i] for i in range(len(value))]
            q.require(len(rows) == 2, f'FeOs retained {kind} are incomplete')
            reference[kind] = sorted(rows, key=lambda row: row['x'][0])
        refs.append(reference)
    return refs, dict(metadata=str(path.resolve()), metadata_sha256=b.digest(path), operands_csv=str(source.resolve()),
                      operands_sha256=b.digest(source), method='Unchanged independently refined FeOs 0.10.1 reference reused; no FeOs calls or copied reference rows')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference', type=Path, help='Reuse a completed pinned FeOs JSON/CSV reference without recalculating it')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    target = args.output/'feos.csv'
    if target.exists() or (args.output/'feos.json').exists():
        raise RuntimeError('output already contains retained values')
    q.require(b.digest(SOURCE) == SOURCE_SHA256, 'FeOs parameter input changed')
    refs = []
    reference_reuse = None
    feos_import = None
    failure = None
    direct = json.loads(b.metadata.distribution('epcsaft').read_text('direct_url.json'))
    wheel = Path(direct['url'].removeprefix('file://'))
    identity = dict(wheel=str(wheel), wheel_sha256=b.digest(wheel), executed_script_sha256=b.digest(__file__),
                    baseline_script_sha256=b.digest(b.__file__), qualifier_script_sha256=b.digest(q.__file__),
                    source_mtimes={str(p.relative_to(b.ROOT)): p.stat().st_mtime_ns for p in (b.ROOT/'engine/native/equilibrium/gibbs_minimizer.cpp', b.ROOT/'engine/native/equilibrium/request_compile.cpp')},
                    thread_environment={name: b.os.environ.get(name) for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS')})
    try:
        if args.reference:
            refs, reference_reuse = retained_reference(args.reference)
        else:
            import feos
            import si_units as si
            from feos import EquationOfState, Parameters, PhaseEquilibrium, State, Verbosity
            feos_import = feos.__file__
            q.require(b.metadata.version('feos') == '0.10.1', 'FeOs version must be 0.10.1')
            identity['si_units_version'] = b.metadata.version('si-units')
            eos = EquationOfState.pcsaft(Parameters.from_json(['water_np', 'hexane'], str(SOURCE)))
            feed = State(eos, temperature=300*si.KELVIN, pressure=1e5*si.PASCAL,
                         composition=np.array([.5, .5])*si.MOL, density_initialization='liquid')
            plane = np.log(np.asarray(feed.molefracs))+np.asarray(feed.ln_phi())
            options = json.loads((b.HERE/'inputs/feos-sources.json').read_text())['comparison']
            feos_record('feos/feed-liquid', feed, plane, si)
            for run, option in enumerate(options['refinements']):
                trials = feed.stability_analysis(max_iter=option['max_iter'], tol=option['tol'], verbosity=Verbosity.Iter)
                records = [feos_record(f'feos/refinement-{run}/trial-{i}', state, plane, si, reduced_trial=True) for i, state in enumerate(trials)]
                flash = PhaseEquilibrium.tp_flash(eos, 300*si.KELVIN, 1e5*si.PASCAL, np.array([.5, .5])*si.MOL,
                                                 max_iter=option['max_iter'], tol=option['tol'], verbosity=Verbosity.Iter)
                phases = [feos_record(f'feos/refinement-{run}/flash-{i}', state, plane, si) for i, state in enumerate((flash.liquid, flash.vapor))]
                gap = np.log(np.asarray(flash.liquid.molefracs))+np.asarray(flash.liquid.ln_phi())-np.log(np.asarray(flash.vapor.molefracs))-np.asarray(flash.vapor.ln_phi())
                q.quantities(f'feos/refinement-{run}/flash', dict(chemical_potential_gap_max_abs=float(np.max(np.abs(gap)))))
                refs.append(dict(trials=sorted(records, key=lambda row: row['x'][0]), phases=sorted(phases, key=lambda row: row['x'][0])))
        for kind in ('trials', 'phases'):
            q.require(len(refs[0][kind]) == len(refs[1][kind]), f'FeOs {kind}: reference phase count changed under refinement')
            for i, (coarse, fine) in enumerate(zip(refs[0][kind], refs[1][kind])):
                dx = max(abs(a-bb) for a, bb in zip(coarse['x'], fine['x']))
                drho = abs(coarse['density_mol_m3']/fine['density_mol_m3']-1)
                if reference_reuse is None:
                    q.quantities(f'feos/refinement-comparison/{kind}-{i}', dict(composition_max_abs=dx, density_max_relative=drho))
                q.require(dx < 1e-8 and drho < 1e-7, f'FeOs {kind}: reference refinement changes exceed issue 5 limits')
        model = b.Mixture(engine_parameters())
        problem = b.eq.Problem(T=300., P=1e5, row_tolerance=1e-10, phases=[b.eq.Phase('L', kind='liquid')],
                              feed=b.eq.Amounts({'water_np': .5, 'hexane': .5}))
        central = q.timed('case-5/FeOs', 'declared', lambda: b.eq.solve_equilibrium(model, problem))
        q.require(central.success, f'case 5 feed failed: {central.message}')
        compiled = b.eq.compile_problem(model, problem)
        report = q.timed('case-5/FeOs', 'explicit-test', lambda: b.eq.stability_test(compiled, central, b.eq.StabilityOptions()))
        q.check_report('case-5/FeOs', report, True)
        problem.phases = b.eq.Discover()
        problem.topology_declared = False
        found = q.timed('case-5/FeOs', 'discover', lambda: b.eq.solve_equilibrium(model, problem))
        q.require(found.success and found.phase_count == 2, f'case 5 discovery failed: {found.message}')
        n = len(model.component_ids)
        phases = sorted([dict(x=found.mole_fractions[i*n:(i+1)*n], density_mol_m3=found.molar_densities[i]) for i in range(found.phase_count)], key=lambda row: row['x'][0])
        for i, (engine, external) in enumerate(zip(phases, refs[1]['phases'])):
            dx = max(abs(a-bb) for a, bb in zip(engine['x'], external['x']))
            drho = abs(engine['density_mol_m3']/external['density_mol_m3']-1)
            q.quantities(f'case-5/flash-comparison-{i}', dict(composition_max_abs=dx, density_max_relative=drho))
            q.require(dx <= 1e-7 and drho <= 1e-6, 'Engine/FeOs flash comparison exceeds issue 5 tolerance')
        engine_trials = sorted([t for t in report.trials if t.seed in ('near-pure-0', 'near-pure-1') and 'Liquid' in str(t.phase.kind)], key=lambda t: t.composition[0])
        q.require(len(engine_trials) == len(refs[1]['trials']), 'Engine/FeOs near-pure trial enumeration differs')
        for i, (engine, external) in enumerate(zip(engine_trials, refs[1]['trials'])):
            dx = max(abs(a-bb) for a, bb in zip(engine.composition, external['x']))
            drho = abs((1/engine.quantities['volume_m3_per_mol'])/external['density_mol_m3']-1)
            q.quantities(f'case-5/stationary-comparison-{i}', dict(composition_max_abs=dx, density_max_relative=drho))
            q.require(dx <= 1e-7 and drho <= 1e-6, 'Engine/FeOs stationary-trial comparison exceeds issue 5 tolerance')
    except Exception as error:
        failure = f'{type(error).__name__}: {error}'
        raise
    finally:
        b.csv_write(target, b.ROWS)
        costs.write(args.output/'feos.json', dict(identity=identity, feos_version='0.10.1', feos_import=feos_import, reference_reuse=reference_reuse,
            sources=json.loads((b.HERE/'inputs/feos-sources.json').read_text()), failure=failure, costs=q.COSTS,
            limits=['independent implementation comparison, not physical validation', 'FeOs negative-TPD stopping relaxation is retained in achieved residuals', 'reused FeOs iteration output remains with the identified original reference; fresh reference output belongs in its run log']))


if __name__ == '__main__':
    main()
