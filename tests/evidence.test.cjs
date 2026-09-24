const { test } = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const root = path.join(__dirname, '..');
const data = JSON.parse(fs.readFileSync(path.join(root, 'data/cohort-evidence.json')));
const html = fs.readFileSync(path.join(root, 'index.html'), 'utf8');

test('coverage and mortality estimates are complete and internally consistent', () => {
  assert.deepEqual(data.cohorts.map(c => [c.cohort, c.marker_count]),
    [['UKB',10], ['NHANES',10], ['CHARLS',8], ['CLHLS',10], ['RuLAS',9]]);
  assert.equal((html.match(/class="heat-cell"/g) || []).length, 47);
  assert.equal((html.match(/class="heat-cell unavailable"/g) || []).length, 3);
  assert.equal((html.match(/<tr data-cohort=/g) || []).length, 5);
  for (const c of data.cohorts) {
    assert.equal(c.markers.length, c.marker_count);
    assert.equal(new Set(c.markers).size, c.marker_count);
    assert.ok(c.markers.every(m => data.markers.includes(m)));
    assert.ok(c.events > 0 && c.events < c.n && c.n <= c.selection.source_n);
    assert.equal(c.n, c.selection.complete_case_n);
    for (const e of Object.values(c.estimates)) {
      assert.ok(Object.values(e).every(Number.isFinite));
      assert.ok(0 <= e.c_low && e.c_low <= e.c_index && e.c_index <= e.c_high && e.c_high <= 1);
      assert.ok(0 < e.hr_low && e.hr_low <= e.hr && e.hr <= e.hr_high);
      assert.ok(e.score_sd > 0 && e.ph_test_p >= 0 && e.ph_test_p <= 1);
    }
    for (const cov of ['age', 'sex', 'bmi', 'education', 'smoking', 'alcohol']) assert.ok(c.covariates.includes(cov));
  }
});

test('published metadata distinguishes apparent adjusted performance and raw DM', () => {
  assert.match(data.analysis.c_index, /full adjusted model; apparent/);
  assert.equal(data.analysis.dm_scale, 'Untransformed Mahalanobis distance');
  assert.match(data.analysis.sampling_weights, /Unweighted/);
  assert.match(html, /proportional hazards/);
  assert.match(html, /Download aggregate results/);
  assert.doesNotMatch(html, /pending confirmation|will be documented/);
});
