"""Read native results and completed captures; no numerical reconstruction."""
import csv
import json
import sys
from contextlib import nullcontext
from pathlib import Path
sys.path.insert(0, str(Path.cwd() / 'analyses/2026-tpd-phase-discovery/scripts'))
import baseline as b
from epcsaft import trace

def call(name, model, operation, mode, output):
    b.ROWS = []
    model._model.reset_counters()
    capture = None if mode == 'off' else trace() if mode == 'ordinary' else trace(diagnostics=True)
    with capture if capture else nullcontext():
        result = operation()
    if hasattr(result, 'convergence'):
        b.state(name, result)
        b.record(name, 'message', result.message)
    else:
        b.record(name, 'molar_density', result.molar_density)
    b.fields(name + '/sweeps', model._model.counters(), ('forward', 'reverse', 'hessian_calls', 'recordings'))
    if getattr(result, 'discovery', None) is not None:
        for key, value in sorted(result.discovery.work.items()): b.record(name + '/search', key, value)
    output.mkdir(parents=True, exist_ok=True)
    b.csv_write(output / f'{name}-state.csv', b.ROWS)
    if capture is not None:
        b.csv_write(output / f'{name}-ordinary.csv', b.trace_rows(capture))
        events, spans = [], []
        width = max(12, max((len(e.values) for e in capture.events), default=0))
        for e in capture.events:
            span = capture.spans[e.span]
            events.append([e.name, e.span, span.parent, span.name, e.iteration, e.reason if hasattr(e, 'reason') else '', *e.values])
        for i, s in enumerate(capture.spans):
            c = s.convergence
            spans.append([i, s.parent, s.name, s.wall_us, *([c.residual, c.criterion.name if c.criterion else '',
                c.criterion.value if c.criterion else '', c.floor, c.amplification, c.contraction,
                c.iterations, c.evaluations, c.stop.name, repr(c)] if c else [''] * 10)])
        for suffix, header, rows in [('events', ['event', 'span', 'parent', 'scope', 'iteration', 'reason', *[f'v{i}' for i in range(width)]], events),
            ('spans', ['span', 'parent', 'scope', 'wall_us', 'residual', 'criterion', 'criterion_value', 'floor', 'amplification', 'contraction', 'iterations', 'evaluations', 'stop', 'reason'], spans)]:
            with (output / f'{name}-{suffix}.csv').open('w', newline='') as stream:
                writer = csv.writer(stream)
                writer.writerow(header)
                writer.writerows(row + [''] * (len(header) - len(row)) for row in rows)
        (output / f'{name}-capture.json').write_text(json.dumps({'incomplete': getattr(capture.trace, 'incomplete', False)}) + '\n')
    return result, capture
