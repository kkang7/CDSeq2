"""Independent CSV verification of Figure S4; no R or inference dependencies.

Usage: python3 verify_expression_prioritization.py OUTPUT_FOLDER [PREVIOUS_VERSION]
The optional previous-version comparison also requires Pillow for pixel checks.
"""
import csv
import hashlib
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

p = Path(sys.argv[1])
td = p / 'supplementary_tables'
def read(name):
    with (td / name).open(newline='') as f:
        return list(csv.DictReader(f))
def key(d):
    return d['scenario'], d['cell_type'], d['gene']
def number(s):
    return float(s) if s not in ('NA', '', 'NaN') else math.nan
def close(a, b):
    assert math.isclose(float(a), float(b), rel_tol=1e-11, abs_tol=1e-11), (a, b)
def digest(f):
    return hashlib.sha256(f.read_bytes()).hexdigest()

m = read('S23_S6_GAM6_membership.csv')
assert len({key(d) for d in m}) == len(m)
groups = defaultdict(list)
for d in m:
    groups[key(d)[:2]].append(d)
selected = set()
for group, rows in groups.items():
    eligible = [d for d in rows if math.isfinite(number(d['GAM6']))]
    eligible.sort(key=lambda d: (-number(d['GAM6']), d['gene']))
    n = math.floor(len(eligible) * 0.2)
    selected.update(key(d) for d in eligible[:n])
    for rank, d in enumerate(eligible, 1):
        assert int(d['GAM6_rank']) == rank
    for d in rows:
        assert int(d['n_common_panel']) == len(rows)
        assert int(d['n_score_eligible']) == len(eligible)
        assert (d['selected'] == 'TRUE') == (key(d) in selected)
        assert (d['rank_eligible'] == 'TRUE') == math.isfinite(number(d['GAM6']))
assert len(groups) == 63
scores = {key(d): number(d['GAM6']) for d in m}
matched = 0
for d in read('S8_eqtl_tests.csv'):
    if d['engine'] == 'CDSeq2' and key(d) in scores:
        close(d['scenario_GAM6_score'], scores[key(d)])
        matched += 1

profiles = read('S14_expression_accuracy.csv')
base = {e: {key(d): number(d['Pearson']) for d in profiles
            if d['scenario'] != 'direct' and d['oracle'] == 'observed'
            and d['engine'] == e} for e in ('BayesPrism', 'CDSeq2')}
assert set(base['BayesPrism']) == set(base['CDSeq2']) == set(scores)
# Compare the panel itself against native BayesPrism preprocessing records.
qc = {(d['scenario'], d['gene']) for d in read('S19_BayesPrism_preprocessing.csv')
      if d['retained'] == 'TRUE'}
for (sc, ct), rows in groups.items():
    assert {d['gene'] for d in rows} == {g for s, g in qc if s == sc}

summary = read('S24_S6_recovery_counts.csv')
assert len(summary) == 210
for d in summary:
    keys = [k for k in scores if k[1] == d['cell_type'] and
            (d['scenario'] == 'ALL' or k[0] == d['scenario'])]
    keep = [k for k in keys if k in selected] if d['workflow'] == 'CDSeq2+GAM6' else keys
    engine = 'BayesPrism' if d['workflow'] == 'BayesPrism' else 'CDSeq2'
    values = [base[engine][k] for k in keep]
    finite = sorted(v for v in values if math.isfinite(v))
    assert int(d['n_common_panel']) == len(keys)
    assert int(d['n_GAM6_eligible']) == sum(math.isfinite(scores[k]) for k in keys)
    assert int(d['n_retained']) == len(keep)
    assert int(d['n_defined']) == len(finite)
    assert int(d['n_undefined']) == len(keep) - len(finite)
    close(d['retained_proportion'], len(keep) / len(keys))
    close(d['median_Pearson'], statistics.median(finite))
    # R type-7 quantile interpolation, computed independently.
    for field, fraction in [('q25_Pearson', .25), ('q75_Pearson', .75)]:
        index = (len(finite) - 1) * fraction
        low = math.floor(index)
        value = finite[low] + (index-low) * (finite[math.ceil(index)]-finite[low])
        close(d[field], value)

print(f'PASS: all {len(m):,} profile memberships, 63 ranking groups and {len(summary)} count/summary rows.')
print(f'PASS: common-panel genes match BayesPrism QC; {matched:,} scores match archived nominated-test scores.')
print('PASS: ranks use only archived score/gene ties; undefined truth correlations do not affect membership.')
for d in summary:
    if d['scenario']=='ALL' and d['workflow']=='CDSeq2+GAM6':
        print(f"PASS: {d['cell_type']}: {d['n_retained']} / {d['n_common_panel']} profiles retained ({100*float(d['retained_proportion']):.6f}%).")
if len(sys.argv) > 2:
    old = Path(sys.argv[2])
    # Compare the seven source/result tables used by this expression benchmark;
    # other applications can change independently in later manuscript revisions.
    numbers = {8, 12, 13, 14, 19, 23, 24}
    original = [f for f in (old / 'supplementary_tables').glob('S[0-9]*_*.csv')
                if int(f.name.split('_')[0][1:]) in numbers]
    assert len(original) == len(numbers)
    assert all(digest(f) == digest(td / f.name) for f in original)
    print('PASS: the seven expression-benchmark source/result tables are byte-identical to the previous version.')
    from PIL import Image, ImageChops
    a = Image.open(next((old / 'figures').glob('Figure_S*_ground_truth_comparison.png'))).convert('RGB')
    b = Image.open(p / 'figures/Figure_S4_ground_truth_comparison.png').convert('RGB')
    assert a.size == b.size == (3000, 2520)
    for box in [(0,0,3000,600),(0,840,3000,1680)]:
        assert ImageChops.difference(a.crop(box), b.crop(box)).getbbox() is None
    print('PASS: S4B and the upper plot region of S4A are pixel-identical; the panel-A workflow legend was renamed.')
