const { test } = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

const html = fs.readFileSync(path.join(__dirname, '..', 'index.html'), 'utf8');
const script = [...html.matchAll(/<script>([\s\S]*?)<\/script>/g)].at(-1)[1];
const examples = [
  { values: [1, 5, 4.5, 1.3, 80, 300, 6, 5, 4.8, 42], disco: -2.118126394406057, raw: 2.152268597777491 },
  { values: [8, 7, 6, 0.8, 130, 480, 12, 9, 3.4, 30], disco: 4.921159259470046, raw: 10.030433484658058 },
  { values: [0, 3.9, 3.6, 1.6, 60, 200, 4, 2.5, 4, 50], disco: 1.4551748543547456, raw: 3.9506639677639126 }
];
// Independent R rank-one covariance calculations, using the published matrices.
const subsetScores = {
  CHARLS: [-4.7177596820757026, 3.2975198908337378, 0.85220305782945438],
  RuLAS: [-2.3242376471311776, 4.054018365019223, 0.30870613426190235]
};
const close = (actual, expected, tolerance = 1e-11) => assert.ok(Math.abs(actual - expected) < tolerance, `${actual} != ${expected}`);

function calculator() {
  const nodes = new Map();
  function getElementById(id) {
    if (!nodes.has(id)) {
      const classes = new Set();
      nodes.set(id, {
        value: '', textContent: '', innerHTML: '', hidden: true, attributes: {}, handlers: {},
        get valueAsNumber() { return this.value.trim() === '' ? NaN : Number(this.value); },
        classList: { add: c => classes.add(c), remove: c => classes.delete(c), contains: c => classes.has(c) },
        addEventListener(event, callback) { this.handlers[event] = callback; },
        setAttribute(name, value) { this.attributes[name] = value; },
        removeAttribute(name) { delete this.attributes[name]; },
        focus() { this.focused = true; }
      });
    }
    return nodes.get(id);
  }
  const context = vm.createContext({ document: { getElementById } });
  vm.runInContext(script, context);
  const api = vm.runInContext('({ MARKERS, COHORT_STATS, COHORT_PANELS, referenceMean: STANDARD_REF.mu, zVector, computeDisco, computeDM, estimatePercentile, percentileLabel, ordinal })', context);
  function fill(values = examples[0].values, cohort = 'UKB') {
    api.MARKERS.forEach((m, i) => { getElementById('in_' + m).value = String(values[i]); });
    getElementById('cohortSelect').value = cohort;
  }
  function submit() { getElementById('discoForm').handlers.submit({ preventDefault() {} }); }
  fill();
  return { ...api, node: getElementById, fill, submit };
}

test('fixed examples preserve full-panel scores and match cohort distributions to their panels', () => {
  const c = calculator();
  for (const [index, example] of examples.entries()) {
    const z = c.zVector(example.values);
    const dm = c.computeDM(z);
    close(c.computeDisco(z), example.disco);
    close(dm, example.raw);
    for (const [cohort, stats] of Object.entries(c.COHORT_STATS)) {
      c.fill(example.values, cohort);
      c.submit();
      assert.equal(c.node('discoResult').textContent, example.disco.toFixed(3));
      assert.equal(c.node('dmResult').textContent, example.raw.toFixed(3));
      const comparison = subsetScores[cohort]?.[index] ?? example.disco;
      close(c.computeDisco(z, c.COHORT_PANELS[cohort]), comparison);
      assert.equal(c.node('zResult').textContent, ((comparison - stats.mean) / stats.sd).toFixed(2));
      assert.ok(c.node('zLabel').textContent.includes(`${c.COHORT_PANELS[cohort].length} markers`));
      assert.equal(c.node('resultSection').classList.contains('show'), true);
      assert.equal(c.node('formError').hidden, true);
    }
  }
});

for (const cohort of ['UKB', 'NHANES', 'CHARLS', 'CLHLS', 'RuLAS']) {
  test(`${cohort}: percentile knots, round-trip endpoints, midpoints and tails`, () => {
    const c = calculator();
    const stats = c.COHORT_STATS[cohort];
    const q = stats.pct.slice(1);
    const p = [5, 10, 25, 50, 75, 90, 95];
    q.forEach((z, i) => {
      const score = stats.mean + stats.sd * z;
      const estimate = c.estimatePercentile((score - stats.mean) / stats.sd, stats);
      close(estimate.value, p[i]);
      assert.equal(estimate.bound, null);
    });
    for (let i = 0; i < q.length - 1; i++) {
      close(c.estimatePercentile((q[i] + q[i + 1]) / 2, stats).value, (p[i] + p[i + 1]) / 2);
    }
    for (const delta of [1e-6, 1, 10]) {
      assert.equal(c.percentileLabel(c.estimatePercentile(q[0] - delta, stats)), '<5th');
      assert.equal(c.percentileLabel(c.estimatePercentile(q[6] + delta, stats)), '>95th');
    }
    assert.ok(c.estimatePercentile(q[0] + 1e-6, stats).value > 5);
    assert.ok(c.estimatePercentile(q[6] - 1e-6, stats).value < 95);
  });
}

test('percentile estimates use correct ordinals and agree between card and note', () => {
  const c = calculator();
  for (const [n, expected] of [[1, '1st'], [2, '2nd'], [3, '3rd'], [11, '11th'], [12, '12th'], [13, '13th'], [21, '21st'], [22, '22nd'], [23, '23rd'], [71, '71st']]) {
    assert.equal(c.ordinal(n), expected);
  }
  c.submit();
  assert.equal(c.node('pctResult').textContent, '≈19th');
  assert.ok(c.node('zNote').textContent.includes('estimated percentile: ≈19th'));
  c.fill(examples[1].values); c.submit();
  assert.equal(c.node('pctResult').textContent, '>95th');
  assert.ok(c.node('zNote').textContent.includes('estimated percentile: >95th'));
  for (const z of [NaN, Infinity, -Infinity]) assert.throws(() => c.estimatePercentile(z, c.COHORT_STATS.UKB));
});

test('every biomarker rejects blank, negative and non-finite inputs without showing stale results', () => {
  const c = calculator();
  for (const marker of c.MARKERS) {
    for (const value of ['', '-1', '-0.5', 'Infinity', 'NaN', '1e309']) {
      c.fill(); c.submit();
      const input = c.node('in_' + marker);
      input.value = value;
      c.submit();
      assert.equal(c.node('resultSection').classList.contains('show'), false);
      assert.equal(c.node('formError').hidden, false);
      assert.match(c.node('formError').textContent, /finite, non-negative/);
      assert.equal(input.attributes['aria-invalid'], 'true');
      assert.equal(input.focused, true);
    }
  }
});

test('finite inputs that overflow the calculation produce an error', () => {
  const c = calculator();
  c.fill(examples[0].values.map(() => 1e308));
  c.submit();
  assert.equal(c.node('resultSection').classList.contains('show'), false);
  assert.match(c.node('formError').textContent, /finite score could not be calculated/);
});

test('input and cohort edits invalidate results; a corrected submission recovers', () => {
  const c = calculator();
  for (const [id, event] of [['discoForm', 'input'], ['cohortSelect', 'change']]) {
    c.fill(); c.submit();
    c.node(id).handlers[event]();
    assert.equal(c.node('resultSection').classList.contains('show'), false);
  }
  c.node('in_CRP').value = '-1'; c.submit();
  c.fill(); c.node('discoForm').handlers.input(); c.submit();
  assert.equal(c.node('formError').hidden, true);
  assert.equal(c.node('in_CRP').attributes['aria-invalid'], undefined);
  assert.equal(c.node('resultSection').classList.contains('show'), true);
});

test('comparison label follows the corresponding panel score and cohort mean', () => {
  const c = calculator();
  const values = [1, 5, 4.5, 1.3, 80, 300, 6, 5, 4.8, 35];
  const seen = new Set();
  for (const [cohort, stats] of Object.entries(c.COHORT_STATS)) {
    const disco = c.computeDisco(c.zVector(values), c.COHORT_PANELS[cohort]);
    c.fill(values, cohort); c.submit();
    const expected = disco < stats.mean ? 'Below cohort mean' : disco > stats.mean ? 'Above cohort mean' : 'At cohort mean';
    assert.equal(c.node('cohortComparison').textContent, expected);
    seen.add(expected);
  }
  assert.equal(seen.size, 2);
});

test('DM at the reference mean is zero and finite without a logarithm', () => {
  const c = calculator();
  assert.equal(c.computeDM(c.referenceMean), 0);
});

test('excluded markers cannot change CHARLS or RuLAS relative estimates', () => {
  const c = calculator();
  for (const [cohort, excluded] of [['CHARLS', ['RBC', 'ALB']], ['RuLAS', ['UREA']]]) {
    c.fill(examples[0].values, cohort); c.submit();
    const initial = ['zResult', 'pctResult', 'cohortComparison', 'zNote'].map(id => c.node(id).textContent);
    const full = c.node('discoResult').textContent;
    for (const marker of excluded) c.node('in_' + marker).value = '20';
    c.submit();
    assert.deepEqual(['zResult', 'pctResult', 'cohortComparison', 'zNote'].map(id => c.node(id).textContent), initial);
    assert.notEqual(c.node('discoResult').textContent, full);
  }
});
