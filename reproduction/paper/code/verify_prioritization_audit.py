"""Independent checks for the additional diagnostics (no R implementation reuse)."""
from pathlib import Path
import csv, math, sys
from collections import defaultdict
p=Path(sys.argv[1]); td=p/'supplementary_tables'
def rows(name):return list(csv.DictReader((td/name).open()))
groups=defaultdict(list)
for r in rows('S37_matched_retention_controls.csv'):
    assert int(r['selected'])==80
    tp,fp,n=map(float,(r['TP'],r['FP'],r['targets']))
    assert tp+fp==float(r['calls'])
    if tp+fp: assert math.isclose(float(r['precision']),tp/(tp+fp),abs_tol=1e-12)
    if n: assert math.isclose(float(r['recall']),tp/n,abs_tol=1e-12)
    groups[(r['id'],r['cell_type'],r['method'])].append(r)
assert len(groups)==30*7*4
assert all(len(g)==(100 if k[2]=='Random' else 1) for k,g in groups.items())
primary=rows('S29_LD_accuracy_coverage.csv')
matched=0
for r in primary:
    if r['workflow']!='CDSeq2_with' or float(r['alpha'])!=.05 or r['scope']=='all' or float(r['h2'])==0:continue
    z=groups[(r['id'],r['scope'],'GAM6')][0]
    for field in ['TP','FP','calls','selected']:assert float(r[field])==float(z[field]),(r['id'],field)
    matched+=1
assert matched==210
decs=defaultdict(list)
for r in rows('S36_score_deciles.csv'):decs[(r['id'],r['cell_type'])].append(r)
assert len(decs)==210
for k,g in decs.items():
    assert len(g)==10 and sum(int(x['selected']) for x in g)==200
    top=[x for x in g if int(x['decile'])<=4]
    z=groups[(*k,'GAM6')][0]
    assert sum(int(x['TP']) for x in top)==int(z['TP'])
    assert sum(int(x['FP']) for x in top)==int(z['FP'])
n=rows('S39_null_stability.csv')
primary=[r for r in n if r['record']=='null_trait' and r['engine']=='CDSeq2' and r['normalization']=='logCPM' and r['n']=='39']
assert len(primary)==200
assert sum(int(r['calls']) for r in primary)==13183
assert sum(int(r['calls'])>0 for r in primary)==5
assert {int(r['trait']):int(r['calls']) for r in primary if int(r['calls'])>1000}=={10:4011,55:4702,112:4467}
omitted=[r for r in n if r['record']=='leave_one_more_donor' and r['omitted_donor']=='38']
assert len(omitted)==3 and all(int(r['calls'])==0 for r in omitted)
print('PASS: 210 matched primary selections; 21,000 random subsets; score groups; null counts and donor sensitivity.')
