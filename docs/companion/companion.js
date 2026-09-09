const area = document.querySelector('#companion');
const connection = document.querySelector('#connection');
const page = area?.dataset.page || 'tests';
const REFRESH_MS = 2000;
const REQUEST_TIMEOUT_MS = 5000;
const BADGE_EXPIRY_MS = 10000;

const state = {
  snapshot: null,
  revision: null,
  online: null,
  lastSuccess: 0,
  request: null,
  requestNumber: 0,
  timer: null,
  expiryTimer: null,
  query: '',
  kind: 'All',
  invalidated: false,
  rendered: false,
  searchFocused: false,
  searchSelection: [0, 0],
};

export function make(tag, className, value) {
  const element = document.createElement(tag);
  if (className) element.className = className;
  if (value !== undefined && value !== null) element.textContent = String(value);
  return element;
}

function text(value, fallback = '—') {
  return value === undefined || value === null || value === '' ? fallback : String(value);
}

function scopeLabel(declaration) {
  const raw = String(declaration?.scope || '').trim().replaceAll('_', ' ');
  if (!raw || raw.toLowerCase() === 'real') return null;
  return raw.charAt(0).toUpperCase() + raw.slice(1);
}

function declarationLabels(declaration) {
  const labels = [make('span', 'kind-label', text(declaration.kind, 'Unclassified'))];
  const scope = scopeLabel(declaration);
  if (scope) labels.push(make('span', 'kind-label scope-label', scope));
  return labels;
}

function array(value) {
  return Array.isArray(value) ? value : [];
}

function math(value, className = 'math', explicit = false) {
  const raw = text(value);
  const latex = explicit || raw.includes('\\') || raw.includes('$') || raw.includes('^') || raw.includes('_')
    ? (raw.startsWith('$') || raw.startsWith('\\(') ? raw : `\\(${raw}\\)`)
    : raw;
  return make('span', className, latex);
}

function unit(value) {
  if (!value) return null;
  const raw = String(value);
  const wrapped = raw.startsWith('$') || raw.startsWith('\\(') ? raw : `\\(${raw}\\)`;
  return make('span', 'units', ` (${wrapped})`);
}

function scientific(value) {
  if (typeof value !== 'number' || !Number.isFinite(value)) return '—';
  if (value === 0) return '0';
  const [coefficient, exponent] = value.toExponential(5).split('e');
  const power = Number(exponent);
  return power === 0 ? coefficient : `${coefficient} × 10^${power}`;
}

function displayValue(value) {
  if (typeof value === 'number') return scientific(value);
  if (typeof value === 'boolean') return value ? 'true' : 'false';
  if (value === null || value === undefined) return '—';
  if (Array.isArray(value)) return `(${value.map(displayValue).join(', ')})`;
  if (typeof value === 'object') {
    return Object.entries(value)
      .filter(([key]) => !key.startsWith('_') && !internalKey(key))
      .map(([key, nested]) => `${key} = ${displayValue(nested)}`)
      .join(', ');
  }
  return String(value);
}

function displayInputs(inputs) {
  const entries = inputEntries(inputs);
  return entries.length ? entries.map(([key, value]) => `${key} = ${displayValue(value)}`).join('; ') : '—';
}

function internalKey(key) {
  return /^(?:source|provenance|path|filename|function|command|identity|digest|sha256|artifact|input_generation|schema_version)$/i.test(key);
}

function formatCriterion(criterion) {
  if (!criterion || typeof criterion !== 'object') return null;
  if (criterion.type && criterion.type !== 'scalar') {
    const parts = [text(criterion.rule || criterion.type, 'Criterion')];
    for (const key of ['threshold', 'scale', 'norm']) if (criterion[key] !== undefined) parts.push(`${key} ${displayValue(criterion[key])}`);
    return parts.join(' · ');
  }
  const parts = [];
  if (typeof criterion.atol === 'number') parts.push(`atol ${scientific(criterion.atol)}`);
  if (typeof criterion.rtol === 'number') parts.push(`rtol ${scientific(criterion.rtol)}`);
  if (criterion.rule) parts.push(String(criterion.rule).replaceAll('_', ' '));
  return parts.join(' · ') || 'Criterion recorded';
}

function badge(label, className = '') {
  const normalized = String(label || 'Unavailable');
  const lower = normalized.toLowerCase();
  let tone = '';
  if (/(passed|current|observed|synced|on)/.test(lower) && !/(outdated|unverified|off)/.test(lower)) tone = 'passed';
  if (/(failed|error|offline|unavailable|outdated|unverified|interrupted|skipped)/.test(lower)) tone = 'failed';
  if (/(running|queued|updating)/.test(lower)) tone = 'active';
  return make('span', `badge ${tone} ${className}`.trim(), normalized);
}

function rowHeader() {
  const row = make('div', 'test-columns');
  for (const label of ['Quantity / Check', 'Inputs', 'Relationship / Comparison', 'Status']) row.append(make('span', '', label));
  return row;
}

function statusFromJudgment(judgment) {
  if (judgment === 'Criterion satisfied') return 'Passed';
  if (judgment === 'Criterion failed') return 'Failed';
  if (judgment === 'Observation only') return 'Observed';
  return 'Unavailable';
}

function freshness(run) {
  return run?.freshness === 'Current' ? 'Current' : run?.freshness === 'Outdated' ? 'Outdated' : 'Unverified';
}

function freshnessOf(runs) {
  const values = runs.map(freshness);
  return values.length && values.every(value => value === 'Current') ? 'Current'
    : values.includes('Unverified') ? 'Unverified' : 'Outdated';
}

function judgmentOf(comparisons) {
  return comparisons.some(row => row.judgment === 'Criterion failed') ? 'Failed'
    : comparisons.some(row => row.judgment === 'Criterion satisfied') ? 'Passed'
    : comparisons.some(row => row.judgment === 'Observation only') ? 'Observed' : 'Unavailable';
}

function executionOf(runs) {
  const executions = runs.map(run => run?.execution);
  if (executions.includes('Running')) return 'Running';
  if (executions.includes('Queued')) return 'Queued';
  if (executions.includes('Error')) return 'Error';
  if (executions.includes('Interrupted')) return 'Interrupted';
  if (executions.includes('Failed')) return 'Failed';
  if (executions.includes('Skipped')) return 'Skipped';
  if (executions.length && executions.every(value => value === 'Passed')) return 'Passed';
  return 'Unavailable';
}

function executionSuffix(runs) {
  const execution = executionOf(runs);
  return ['Passed', 'Running', 'Queued'].includes(execution) ? '' : ` · test ${execution.toLowerCase()}`;
}

function aggregate(declaration, runs) {
  if (declaration.lifecycle === 'Planned') return {execution: 'Planned', freshness: null, judgment: null};
  if (!runs.length) return {execution: 'Not run', freshness: 'Unverified', judgment: 'Unavailable'};
  const execution = executionOf(runs);
  if (execution === 'Running') return {execution, freshness: null, judgment: null};
  if (execution === 'Queued') return {execution, freshness: null, judgment: null};
  const aggregateFreshness = freshnessOf(runs);
  const comparisons = runs.flatMap(run => array(run.comparisons));
  const judgment = comparisons.length ? judgmentOf(comparisons) : null;
  return {execution, freshness: aggregateFreshness, judgment};
}

function statusLabel(declaration, runs, snapshot) {
  const summary = aggregate(declaration, runs);
  if (summary.execution === 'Planned') return 'Planned';
  if (summary.execution === 'Running' || summary.execution === 'Queued') return summary.execution;
  if (!runs.length) return 'Not run';
  const invalid = state.invalidated || snapshot.state !== 'Synced' || Date.now() - state.lastSuccess > BADGE_EXPIRY_MS;
  const inventoryComplete = snapshot.inventory?.state === 'Complete';
  if (invalid || !inventoryComplete) return `Unverified · last ${summary.execution.toLowerCase()}`;
  if (summary.freshness !== 'Current') return `${summary.freshness} · last ${summary.execution.toLowerCase()}`;
  if (summary.judgment === 'Observed') return 'Observed · Current';
  if (summary.judgment === 'Unavailable') return `Unavailable · ${summary.execution}`;
  if (summary.judgment === 'Failed') return 'Failed · Current';
  return `${summary.execution} · Current`;
}

function comparisonStatusLabel(comparison, run, snapshot) {
  if (run.execution === 'Running' || run.execution === 'Queued') return run.execution;
  const judgment = statusFromJudgment(comparison.judgment);
  const invalid = state.invalidated || snapshot.state !== 'Synced' || !recordCurrent(run);
  const stale = snapshot.state === 'Synced' && !state.invalidated ? freshness(run) : 'Unverified';
  const suffix = executionSuffix([run]);
  if (invalid) return `${stale} · last ${judgment.toLowerCase()}${suffix}`;
  if (suffix) return `${judgment}${suffix}`;
  return `${judgment} · Current`;
}

function comparisonTitle(comparison) {
  const title = make('span', 'comparison-title');
  title.append(math(comparison.quantity_latex || comparison.quantity, 'math quantity-name', Boolean(comparison.quantity_latex)));
  const quantityUnit = unit(comparison.unit);
  if (quantityUnit) title.append(quantityUnit);
  if (comparison.basis || comparison.method) title.append(make('small', 'comparison-meta', [comparison.basis, comparison.method].filter(Boolean).join(' · ')));
  return title;
}

function observations(comparison) {
  return array(comparison.observations);
}

function recordCurrent(run) {
  return freshness(run) === 'Current'
    && !state.invalidated
    && state.snapshot?.state === 'Synced'
    && Date.now() - state.lastSuccess <= BADGE_EXPIRY_MS
    && state.snapshot?.inventory?.state === 'Complete';
}

function inputObject(value) {
  return value && typeof value === 'object' && !Array.isArray(value) ? value : {};
}

function sharedInputsFor(run, comparison) {
  const groupRef = comparison?.group_ref;
  if (groupRef === undefined || groupRef === null) return {};
  const group = array(run?.input_groups).find(item => item && item.group_ref === groupRef);
  return inputObject(group?.shared_inputs);
}

function effectiveInputs(run, comparison, override) {
  return {...sharedInputsFor(run, comparison), ...inputObject(comparison?.inputs), ...inputObject(override)};
}

function comparisonRows(records) {
  return records.flatMap(record => {
    const comparison = record.comparison;
    const points = observations(comparison);
    const inputs = record.inputs || effectiveInputs(record.run, comparison);
    if (!points.length) return [{...comparison, inputs, _run: record.run, _comparison: comparison}];
    return points.map(point => ({...comparison, ...point, inputs: effectiveInputs(record.run, comparison, point.inputs), coordinates: point.coordinates || comparison.coordinates,
      model: point.model ?? comparison.model, reference: point.reference ?? comparison.reference,
      criterion: point.criterion || comparison.criterion, judgment: point.judgment || comparison.judgment,
      _run: record.run, _comparison: comparison}));
  });
}

function commonInputText(rows) {
  const maps = rows.map(row => Object.fromEntries(inputEntries(row.inputs)));
  const keys = [...new Set(maps.flatMap(map => Object.keys(map)))];
  return keys.filter(key => {
    const values = maps.map(map => map[key]).filter(value => value !== undefined);
    return values.length === maps.length && new Set(values.map(displayValue)).size === 1;
  }).map(key => `${key} = ${displayValue(maps[0][key])}`).join('; ');
}

function groupedNumericalTable(records) {
  const rows = comparisonRows(records);
  const hasNumericOperands = rows.some(row => 'model' in row || 'reference' in row);
  if (!rows.length || !hasNumericOperands) return null;
  const coordinateKeys = [...new Set(rows.flatMap(row => Object.keys(row.coordinates || {})))];
  const inputKeys = [...new Set(rows.flatMap(row => inputEntries(row.inputs).map(([key]) => key)))];
  const varyingInputKeys = inputKeys.filter(key => {
    const values = rows.map(row => Object.fromEntries(inputEntries(row.inputs))[key]).filter(value => value !== undefined);
    return values.length !== rows.length || new Set(values.map(displayValue)).size > 1;
  });
  const current = records.every(record => recordCurrent(record.run));
  const wrapper = make('div', 'table-scroll');
  const table = make('table', 'numerics');
  const head = make('thead');
  const header = make('tr');
  if (rows.some(row => row.point_id)) header.append(make('th', '', 'Point'));
  for (const key of coordinateKeys) header.append(make('th', '', key));
  for (const key of varyingInputKeys) header.append(make('th', '', key));
  for (const label of ['Reference', 'Model', 'Error', 'Criterion', 'Verdict']) header.append(make('th', '', label));
  head.append(header);
  table.append(head);
  const body = make('tbody');
  for (const row of rows) {
    const tr = make('tr');
    if (rows.some(item => item.point_id)) tr.append(make('td', '', text(row.point_id)));
    for (const key of coordinateKeys) tr.append(make('td', '', displayValue(row.coordinates?.[key])));
    const inputs = Object.fromEntries(inputEntries(row.inputs));
    for (const key of varyingInputKeys) tr.append(make('td', '', displayValue(inputs[key])));
    const reference = row.reference;
    const model = row.model;
    const error = row.absolute_error ?? (typeof model === 'number' && typeof reference === 'number' ? Math.abs(model - reference) : null);
    const values = current
      ? [reference, model, error, formatCriterion(row.criterion) || '—', statusFromJudgment(row.judgment)]
      : ['', '', '', '', ''];
    for (const value of values) tr.append(make('td', '', typeof value === 'number' ? displayValue(value) : value));
    body.append(tr);
  }
  table.append(body);
  wrapper.append(table);
  if (!current) wrapper.append(make('p', 'table-note', 'Current result cells are empty because this record is no longer applicable.'));
  return wrapper;
}

function comparisonDetailsGroup(group, snapshot) {
  const records = group.records;
  const representative = records[0].comparison;
  const rows = comparisonRows(records);
  const block = make('div', 'comparison');
  const heading = make('div', 'comparison-heading');
  heading.append(comparisonTitle(representative), badge(comparisonGroupStatus(group, snapshot)));
  block.append(heading);
  const label = inputLabel(records.map(record => record.comparison));
  const common = commonInputText(rows);
  if (label || common) {
    const line = make('p', 'detail-line');
    if (label) line.append(math(label, 'math input-label', true));
    if (common) line.append(make('span', 'actual-inputs', `${label ? ' · ' : ''}${common}`));
    block.append(line);
  }
  const relation = representative.relationship_latex || representative.claim || representative.method;
  if (relation) {
    const line = make('p', 'detail-line');
    line.append(make('span', 'detail-label', 'Comparison: '), math(relation, 'math', Boolean(representative.relationship_latex)));
    block.append(line);
  }
  const table = groupedNumericalTable(records);
  if (table) block.append(table);
  else if (rows.length > 1) block.append(make('p', 'detail-line', `${rows.length} state points recorded without a numerical table.`));
  else if (representative.judgment === 'Observation only') block.append(make('p', 'detail-line', 'Outcome recorded as an observation; no numerical table was emitted.'));
  if (records.some(record => record.comparison.evidence_state === 'Unavailable')) block.append(make('p', 'detail-line', 'Evidence unavailable for this comparison.'));
  const executions = [...new Set(records.map(record => record.run.execution).filter(execution => execution !== 'Passed'))];
  if (executions.length) block.append(make('p', 'detail-line', `Execution: ${executions.join(', ')}`));
  return block;
}

function runDetails(run) {
  const details = make('div', 'run-details');
  details.append(make('p', 'run-heading', run.parameter_case_id && run.parameter_case_id !== 'single' ? `Parameter case ${run.parameter_case_id}` : 'Execution'));
  const records = comparisonRecords([run]);
  const groups = comparisonGroups(records);
  if (!groups.length) details.append(make('p', 'detail-line', run.execution === 'Passed' ? 'Outcome recorded without numerical comparisons.' : 'No complete comparison was recorded.'));
  else for (const group of groups) details.append(comparisonDetailsGroup(group, state.snapshot || {}));
  return details;
}

function declarationsFor(snapshot) {
  return array(snapshot.declarations).filter(declaration => declaration && typeof declaration === 'object');
}

function runsFor(snapshot, testId) {
  return array(snapshot.runs).filter(run => run && run.test_id === testId);
}

function comparisonRecords(runs) {
  return runs.flatMap(run => array(run.comparisons).map(comparison => ({comparison, run, inputs: effectiveInputs(run, comparison)})));
}

function comparisonGroups(records) {
  const groups = new Map();
  for (const record of records) {
    const key = record.comparison.comparison_group_id || record.comparison.comparison_id || record.comparison.point_id;
    if (!key) continue;
    if (!groups.has(key)) groups.set(key, {key, records: []});
    groups.get(key).records.push(record);
  }
  return [...groups.values()];
}

function relationshipCell(comparison) {
  const value = comparison.relationship_latex || comparison.claim || comparison.method;
  if (!value) return make('span', 'relationship', '—');
  const cell = make('span', 'relationship');
  cell.append(math(value, 'math', Boolean(comparison.relationship_latex)));
  return cell;
}

function inputEntries(inputs) {
  if (!inputs || typeof inputs !== 'object' || Array.isArray(inputs)) return [];
  return Object.entries(inputs).filter(([key]) => !key.startsWith('_') && !internalKey(key));
}

function recordInputs(record) {
  return record?.comparison ? record.inputs || effectiveInputs(record.run, record.comparison) : record?.inputs;
}

function compactInputs(records) {
  const maps = records.map(record => Object.fromEntries(inputEntries(recordInputs(record))));
  const keys = [...new Set(maps.flatMap(map => Object.keys(map)))];
  const common = [];
  const varying = [];
  for (const key of keys) {
    const values = maps.map(map => map[key]).filter(value => value !== undefined);
    const displays = [...new Set(values.map(displayValue))];
    if (values.length === maps.length && displays.length === 1) common.push(`${key} = ${displays[0]}`);
    else if (displays.length) varying.push(`${key}: ${displays.join(' → ')}`);
  }
  return [...common, varying.length ? `varying ${varying.join('; ')}` : ''].filter(Boolean).join('; ');
}

function inputLabel(comparisons) {
  const labels = [...new Set(comparisons.map(comparison => comparison.inputs_latex || comparison.input_labels).filter(Boolean))];
  return labels.length === 1 ? labels[0] : null;
}

function inputsCell(records) {
  const comparisons = records.map(record => record.comparison || record);
  const cell = make('span', 'inputs');
  const label = inputLabel(comparisons);
  if (label) cell.append(math(label, 'math input-label', true));
  const actual = compactInputs(records);
  if (actual) cell.append(make('span', 'actual-inputs', `${label ? ' · ' : ''}${comparisons.length > 1 ? `${comparisons.length} state points · ` : ''}${actual}`));
  else if (comparisons.length > 1) cell.append(make('span', 'actual-inputs', `${comparisons.length} state points`));
  else if (!label) cell.textContent = '—';
  return cell;
}

function comparisonGroupStatus(group, snapshot) {
  const records = group.records;
  const runs = records.map(record => record.run);
  const execution = executionOf(runs);
  if (execution === 'Running' || execution === 'Queued') return execution;
  const judgment = judgmentOf(records.map(record => record.comparison));
  const current = freshnessOf(runs);
  const invalid = state.invalidated || snapshot.state !== 'Synced' || !runs.every(record => recordCurrent(record));
  const stale = snapshot.state === 'Synced' && !state.invalidated ? current : 'Unverified';
  const suffix = executionSuffix(runs);
  if (invalid || current !== 'Current') return `${stale} · last ${judgment.toLowerCase()}${suffix}`;
  if (suffix) return `${judgment}${suffix}`;
  return `${judgment} · Current`;
}

function comparisonRow(declaration, group, snapshot) {
  const row = make('details', 'test-row');
  row.id = `comparison-${declaration.test_id}-${group.key}`;
  const representative = group.records[0].comparison;
  const summary = make('summary');
  const quantity = make('span', 'quantity');
  quantity.append(...declarationLabels(declaration), comparisonTitle(representative));
  summary.append(quantity, inputsCell(group.records), relationshipCell(representative), badge(comparisonGroupStatus(group, snapshot)));
  row.append(summary);
  const inside = make('div', 'row-details');
  const scope = scopeLabel(declaration);
  if (scope) inside.append(make('p', 'detail-line', `Scope: ${scope}`));
  inside.append(comparisonDetailsGroup(group, snapshot));
  row.append(inside);
  return row;
}

function checkRow(declaration, runs, snapshot) {
  const row = make('details', 'test-row');
  row.id = `test-${declaration.test_id}`;
  const summary = make('summary');
  const quantity = make('span', 'quantity');
  quantity.append(...declarationLabels(declaration), make('span', 'check-name', text(declaration.title, 'Unclassified check')));
  summary.append(quantity);
  const inputs = new Set(comparisonRecords(runs).map(record => displayInputs(record.inputs)).filter(value => value !== '—'));
  summary.append(make('span', 'inputs', inputs.size === 1 ? [...inputs][0] : inputs.size ? `${inputs.size} recorded input cases` : '—'));
  const comparisons = runs.flatMap(run => array(run.comparisons));
  const relationships = comparisons.map(comparison => comparison.relationship_latex || comparison.claim || comparison.method).filter(Boolean);
  const relationship = make('span', 'relationship');
  if (relationships.length === 1) relationship.append(math(relationships[0], 'math', Boolean(comparisons[0].relationship_latex)));
  else relationship.textContent = relationships.length ? `${relationships.length} recorded comparisons` : declaration.evidence_expectation === 'Outcome only' ? 'Outcome only' : text(declaration.unmet_prerequisite);
  summary.append(relationship);
  summary.append(badge(statusLabel(declaration, runs, snapshot)));
  row.append(summary);
  const inside = make('div', 'row-details');
  const scope = scopeLabel(declaration);
  if (scope) inside.append(make('p', 'detail-line', `Scope: ${scope}`));
  if (declaration.lifecycle === 'Planned') inside.append(make('p', 'detail-line', `Planned: ${text(declaration.unmet_prerequisite, 'Prerequisite not recorded')}`));
  else if (!runs.length) inside.append(make('p', 'detail-line', 'No execution record is available for this declaration.'));
  else for (const run of runs) inside.append(runDetails(run));
  row.append(inside);
  return row;
}

function categoryOrder(snapshot, declarations) {
  const plan = array(snapshot.plan).map(group => group.name).filter(Boolean);
  const categories = [...new Set(declarations.map(declaration => declaration.category).filter(Boolean))];
  if (array(snapshot.studies).length) categories.push('Studies');
  return [...plan, ...categories.filter(category => !plan.includes(category))];
}

function searchable(declaration, runs) {
  return [declaration.test_id, declaration.title, declaration.kind, declaration.scope, declaration.category, declaration.requirement_id, declaration.unmet_prerequisite,
    ...comparisonRecords(runs).flatMap(record => {
      const comparison = record.comparison;
      return [comparison.quantity_latex, comparison.quantity, comparison.claim, comparison.method, comparison.basis,
        comparison.relationship_latex, displayInputs(record.inputs)];
    })]
    .filter(Boolean).join(' ').toLowerCase();
}

function controls(snapshot, render) {
  const wrapper = make('div', 'test-controls');
  const search = make('input', 'search');
  search.id = 'test-search';
  search.type = 'search';
  search.placeholder = 'Find a quantity or scientific check';
  search.setAttribute('aria-label', 'Find a quantity or scientific check');
  search.value = state.query;
  search.addEventListener('input', () => { state.query = search.value; render(); });
  const select = make('select', 'kind-filter');
  select.setAttribute('aria-label', 'Filter by kind');
  for (const kind of ['All', ...new Set(declarationsFor(snapshot).map(declaration => declaration.kind).filter(Boolean))]) {
    const option = make('option', '', kind);
    option.value = kind;
    option.selected = state.kind === kind;
    select.append(option);
  }
  select.addEventListener('change', () => { state.kind = select.value; render(); });
  wrapper.append(search, select);
  return wrapper;
}

function planView(snapshot) {
  const open = new Set([...area.querySelectorAll('details[open]')].map(element => element.id));
  area.replaceChildren();
  const plan = array(snapshot.plan);
  const rows = plan.flatMap(group => array(group.rows));
  const qualified = rows.filter(row => ['Qualified', '[x]'].includes(row.status)).length;
  const metrics = make('div', 'metrics');
  for (const [value, label] of [[`${qualified} / ${rows.length}`, 'Qualified sublayers'], [String(plan.length), 'Capability groups']]) {
    const metric = make('div', 'metric');
    metric.append(make('b', '', value), make('span', '', label));
    metrics.append(metric);
  }
  area.append(metrics);
  for (const [index, group] of plan.entries()) {
    const layer = make('details', 'layer');
    layer.id = `layer-${index}`;
    layer.open = open.has(layer.id);
    const summary = make('summary');
    summary.append(make('small', '', String(index + 1).padStart(2, '0')), make('strong', '', text(group.name, 'Capability group')), badge(`${array(group.rows).filter(row => ['Qualified', '[x]'].includes(row.status)).length} / ${array(group.rows).length}`));
    layer.append(summary);
    const inside = make('div', 'inside');
    for (const row of array(group.rows)) {
      const work = make('div', 'work');
      const head = make('div', 'work-head');
      head.append(make('small', '', text(row.id)), make('strong', '', text(row.name)), badge(text(row.status)));
      work.append(head, make('p', '', text(row.required)));
      if (row.remaining) work.append(make('p', '', `${row.remaining} h remaining · ${row.actual === '—' ? 'Time not measured' : text(row.actual) + ' h measured'}`));
      inside.append(work);
    }
    layer.append(inside);
    area.append(layer);
  }
  area.append(make('p', 'page-note', 'Qualification reflects accepted evidence recorded by the scientific owners.'));
  if (array(snapshot.preparation).length) {
    const preparation = make('details', 'layer');
    preparation.id = 'test-preparation';
    preparation.open = open.has(preparation.id);
    const summary = make('summary');
    summary.append(make('strong', '', 'Test-system preparation'));
    preparation.append(summary);
    const inside = make('div', 'inside');
    for (const item of array(snapshot.preparation)) {
      const work = make('div', 'work');
      const head = make('div', 'work-head');
      head.append(make('strong', '', text(item.name, 'Preparation')), badge(text(item.status)));
      work.append(head, make('p', '', text(item.required || item.completion)));
      inside.append(work);
    }
    preparation.append(inside);
    area.append(preparation);
  }
}

function testsView(snapshot) {
  const open = new Set([...area.querySelectorAll('details[open]')].map(element => element.id));
  const searchInput = area.querySelector('#test-search');
  state.searchFocused = document.activeElement === searchInput;
  state.searchSelection = searchInput ? [searchInput.selectionStart ?? state.query.length, searchInput.selectionEnd ?? state.query.length] : state.searchSelection;
  const declarations = declarationsFor(snapshot);
  const render = () => testsView(snapshot);
  area.replaceChildren();
  const banner = make('div', 'snapshot-banner');
  const inventory = snapshot.inventory || {};
  banner.append(badge(inventory.state === 'Complete' ? 'Inventory complete' : 'Inventory incomplete'));
  banner.append(badge(`Automatic execution ${snapshot.automatic_execution === 'on' ? 'on' : 'off'}`, snapshot.automatic_execution === 'on' ? 'automatic-on' : 'automatic-off'));
  if (snapshot.worker?.state && snapshot.worker.state !== 'Idle') banner.append(badge(`Worker ${snapshot.worker.state}`, 'active'));
  if (array(snapshot.runs).some(run => ['Queued', 'Running'].includes(run.execution))) banner.append(badge('Checks queued or running', 'active'));
  if (inventory.state !== 'Complete') banner.append(make('span', 'banner-note', 'Current results are unverified until discovery completes.'));
  area.append(banner, controls(snapshot, render));
  const search = state.query.trim().toLowerCase();
  for (const [index, category] of categoryOrder(snapshot, declarations).entries()) {
    const selected = declarations.filter(declaration => declaration.category === category && (state.kind === 'All' || declaration.kind === state.kind));
    const visible = selected.filter(declaration => !search || searchable(declaration, runsFor(snapshot, declaration.test_id)).includes(search));
    const categoryStudies = category === 'Studies' ? array(snapshot.studies).filter(study => !search || JSON.stringify(study).toLowerCase().includes(search)) : [];
    if (!visible.length && !categoryStudies.length) continue;
    const layer = make('details', `layer ${category === 'Tooling' ? 'tooling-layer' : ''}`.trim());
    layer.id = `tests-category-${category}-${index}`;
    layer.open = open.has(layer.id);
    const summary = make('summary');
    summary.append(make('strong', '', category));
    layer.append(summary);
    const inside = make('div', 'inside');
    if (visible.length) {
      const group = make('section', 'test-group');
      group.append(rowHeader());
      for (const declaration of visible) {
        const runs = runsFor(snapshot, declaration.test_id);
        const groups = comparisonGroups(comparisonRecords(runs));
        if (declaration.lifecycle === 'Planned' || declaration.evidence_expectation === 'Outcome only' || !groups.length) group.append(checkRow(declaration, runs, snapshot));
        else for (const comparison of groups) group.append(comparisonRow(declaration, comparison, snapshot));
      }
      inside.append(group);
    }
    if (categoryStudies.length) {
      const group = make('section', 'test-group');
      group.append(make('h3', '', 'Studies'), rowHeader());
      for (const study of categoryStudies) {
        const declaration = {test_id: study.study_id || study.id || `study-${categoryStudies.indexOf(study)}`, title: study.title || 'Explicit study', kind: study.kind || 'Numerical', category: 'Studies', lifecycle: study.lifecycle || 'Implemented', evidence_expectation: study.evidence_expectation || 'Outcome only'};
        const runs = study.run ? [study.run] : [];
        const groups = comparisonGroups(comparisonRecords(runs));
        if (declaration.evidence_expectation === 'Outcome only' || !groups.length) group.append(checkRow(declaration, runs, snapshot));
        else for (const comparison of groups) group.append(comparisonRow(declaration, comparison, snapshot));
      }
      inside.append(group);
    }
    layer.append(inside);
    area.append(layer);
  }
  if (!area.querySelector('.layer')) area.append(make('p', 'empty', 'No matching checks.'));
  for (const detail of area.querySelectorAll('details')) if (open.has(detail.id)) detail.open = true;
  if (state.searchFocused) {
    const input = area.querySelector('#test-search');
    if (input) {
      input.focus();
      input.setSelectionRange(...state.searchSelection);
    }
  }
}

function renderSnapshot(snapshot) {
  if (page === 'plan') planView(snapshot);
  else testsView(snapshot);
  if (window.MathJax?.typesetPromise) window.MathJax.typesetPromise([area]);
}

function paintConnection() {
  if (!connection) return;
  if (state.online === false || state.snapshot?.state === 'Offline') connection.textContent = '○ Offline';
  else if (state.snapshot?.state === 'Updating') connection.textContent = '◌ Updating';
  else if (state.invalidated || Date.now() - state.lastSuccess > BADGE_EXPIRY_MS) connection.textContent = '○ Updating';
  else if (state.snapshot?.inventory?.state !== 'Complete') connection.textContent = '◌ Updating';
  else connection.textContent = '● Synced';
}

function schedule() {
  clearTimeout(state.timer);
  state.timer = setTimeout(refresh, REFRESH_MS);
}

function scheduleExpiry() {
  clearTimeout(state.expiryTimer);
  if (!state.lastSuccess) return;
  const delay = Math.max(1, BADGE_EXPIRY_MS - (Date.now() - state.lastSuccess) + 1);
  state.expiryTimer = setTimeout(() => {
    if (Date.now() - state.lastSuccess < BADGE_EXPIRY_MS) return scheduleExpiry();
    state.invalidated = true;
    if (state.snapshot && state.rendered) renderSnapshot(state.snapshot);
    paintConnection();
  }, delay);
}

async function refresh() {
  const requestNumber = ++state.requestNumber;
  if (state.request) state.request.abort();
  const controller = new AbortController();
  state.request = controller;
  paintConnection();
  const timeout = setTimeout(() => controller.abort(), REQUEST_TIMEOUT_MS);
  try {
    const response = await fetch('/companion', {cache: 'no-store', signal: controller.signal});
    if (!response.ok) throw new Error('snapshot unavailable');
    const snapshot = await response.json();
    const updating = snapshot?.state === 'Updating';
    const requiresRevision = snapshot?.state === 'Synced';
    if (requestNumber !== state.requestNumber || !snapshot || typeof snapshot !== 'object' || (requiresRevision && !snapshot.revision) || snapshot.presentation_available === false) throw new Error('incoherent snapshot');
    const previous = state.snapshot;
    const wasInvalidated = state.invalidated;
    const changed = !previous
      || (!updating && state.revision !== snapshot.revision)
      || (wasInvalidated && !updating)
      || previous.state !== snapshot.state
      || previous.automatic_execution !== snapshot.automatic_execution
      || previous.inventory?.state !== snapshot.inventory?.state
      || previous.worker?.state !== snapshot.worker?.state;
    state.snapshot = snapshot;
    state.online = true;
    if (updating) {
      state.invalidated = true;
    } else {
      state.revision = snapshot.revision;
      state.lastSuccess = Date.now();
      state.invalidated = false;
    }
    if (changed) renderSnapshot(snapshot);
    state.rendered = true;
    if (!updating) scheduleExpiry();
    paintConnection();
  } catch (_) {
    if (requestNumber === state.requestNumber) {
      state.online = false;
      paintConnection();
    }
  } finally {
    clearTimeout(timeout);
    if (requestNumber === state.requestNumber) {
      state.request = null;
      schedule();
    }
  }
}

document.addEventListener('visibilitychange', () => {
  if (document.visibilityState === 'hidden') {
    clearTimeout(state.timer);
    clearTimeout(state.expiryTimer);
    if (state.request) state.request.abort();
    return;
  }
  state.invalidated = true;
  paintConnection();
  if (state.snapshot && state.rendered) renderSnapshot(state.snapshot);
  refresh();
});

refresh();
