"""Bounded inverse-isobar and label-invariance checks of a converged disagreement."""
import copy
import gzip
import json
from pathlib import Path
import time

from epcsaft import Mixture, Parameters, provenance
from epcsaft.equilibrium import Free, MoleFractions, Phase, Pinned, Problem, solve_equilibrium

from reproduce import LEAF, RESULTS, declarations, dump, sha

run = json.loads(gzip.decompress((RESULTS / 'native-results.json.gz').read_bytes()))
original = next(p for p in run['points'] if p['row_id'] == '6a/bubble/1' and p['offset'] == 0 and p['requested_tolerance'] == 1e-10)
model = declarations()
document = json.loads((LEAF / 'inputs/parameters/parameters.json').read_text())
swapped = copy.deepcopy(document)
for site in swapped['topology']['sites']:
    if site['component_id'] == 'MDEA':
        site['multiplicity'] = 3 - site['multiplicity']
identity = dict(parent_run_id=run['execution_identity']['run_id'], engine_provenance=provenance(),
                extension_sha256=sha(Path(__import__('epcsaft')._core.__file__)),
                producer_sha256=sha(Path(__file__)),
                parent_result_sha256=sha(RESULTS / 'native-results.json.gz'))
records = []
started = time.perf_counter()
for name, candidate, temperature, pressure in [
    ('inverse-bubble', model, original['temperature_K'], Free(original['pressure_Pa'])),
    ('MDEA-site-label-swap', Mixture(Parameters.from_mapping(swapped, components=['MDEA', 'water'])), Free(original['temperature_K']), original['pressure_Pa'])]:
    x = original['feed_mdea']
    result = solve_equilibrium(candidate, Problem(T=temperature, P=pressure,
        feed=MoleFractions({'MDEA': x, 'water': 1-x}),
        phases=[Phase('L', kind='liquid'), Phase('V', kind='vapor', amount=Pinned(0))], row_tolerance=1e-11))
    records.append(dict(name=name, source_row_id=original['row_id'], execution_identity=identity,
                        temperature_K=result.temperature, pressure_Pa=result.pressure,
                        success=result.success, residuals=result.classified_residuals,
                        message=result.message, requested_tolerance_met=result.requested_tolerance_met))
dump(RESULTS / 'diagnosis.json', dict(records=records, wall_seconds=time.perf_counter()-started,
     conclusion='Tighter row tolerance, inverse pressure solve and site-label invariance do not explain the bubble-curve disagreement; cause undetermined. No parameter change adopted.'))
