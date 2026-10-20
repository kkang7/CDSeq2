"""Synchronize book assets from a reviewed manuscript bundle: python3 tools/sync_manuscript.py BUNDLE."""
from pathlib import Path
import csv, hashlib, json, re, shutil, subprocess, sys, tarfile, zipfile
repo=Path(__file__).resolve().parents[1];paper=Path(sys.argv[1]).resolve();dest=repo/'book/downloads';dest.mkdir(exist_ok=True)
tables=sorted((paper/'supplementary_tables').glob('S[0-9]*_*.csv'),key=lambda p:int(re.match(r'S(\d+)',p.name)[1]))
figs=list((paper/'figures').glob('*.png'));assert len(tables)==42 and len(figs)==25
for p in figs:shutil.copy2(p,repo/'book/assets/figures'/p.name)
with zipfile.ZipFile(dest/'supplementary_tables.zip','w',zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for p in tables:z.write(p,'supplementary_tables/'+p.name)
with (dest/'table_manifest.csv').open('w') as f:
    w=csv.writer(f);w.writerow(['table','file','bytes','sha256'])
    for p in tables:w.writerow([re.match(r'S\d+',p.name)[0],p.name,p.stat().st_size,hashlib.sha256(p.read_bytes()).hexdigest()])
for p in list(paper.glob('*.R'))+list(paper.glob('*.py')):shutil.copy2(p,repo/'reproduction/paper/code'/p.name)
for p in ['manuscript.pdf','supplement.pdf']:shutil.copy2(paper/p,dest/p)
with tarfile.open(dest/'independent_LD_benchmark.tar.gz','w:gz',compresslevel=9) as z:z.add(paper/'CDSeq2_independent_LD_benchmark_v2',arcname='independent_ld')
with zipfile.ZipFile(dest/'reconstructed_brain_inputs.zip','w',zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for p in sorted((paper/'ReconstructedData').rglob('*')):
        if p.is_file():z.write(p,'ReconstructedData/'+str(p.relative_to(paper/'ReconstructedData')))
def balanced(s,start):
    depth=1;i=start
    while depth:
        if s[i]=='{' and s[i-1]!='\\':depth+=1
        elif s[i]=='}' and s[i-1]!='\\':depth-=1
        i+=1
    return s[start:i-1]
pages={'figures':['# Main figures\n\nAll seven main figures and captions match manuscript version 16. [Download the manuscript](../downloads/manuscript.pdf) or [underlying tables](data.qmd).'], 'supplementary':['# Supplementary figures\n\nAll eighteen supplementary figures match manuscript version 16. [Download the supplement](../downloads/supplement.pdf).']};mapping=[]
for doc in ['manuscript.tex','supplement.tex']:
    for b in re.findall(r'\\begin\{figure\}.*?\\end\{figure\}',(paper/doc).read_text(),re.S):
        m=re.search(r'\{figures/(Figure_([^_]+)_[^}]+)\.pdf\}',b)
        if not m:continue
        fn,num=m.groups();cap=balanced(b,b.index('\\caption{')+9)
        pandoc=[shutil.which('pandoc')] if shutil.which('pandoc') else ['quarto','pandoc']
        md=subprocess.run(pandoc+['--from=latex','--to=markdown','--wrap=none'],input=cap,text=True,capture_output=True,check=True).stdout.strip()
        tail=cap.split('Underlying data are provided')[-1]; refs=[]
        for x,y in re.findall(r'S(\d+)\s*--\s*S?(\d+)',tail):refs.extend('S'+str(i) for i in range(int(x),int(y)+1))
        refs.extend(re.findall(r'\bS\d+\b',tail));refs=sorted(set(refs),key=lambda x:int(x[1:]));assert refs and all(int(x[1:])<=42 for x in refs)
        mapping.append(dict(figure=num,supplementary_tables=';'.join(refs)))
        page='supplementary' if num.startswith('S')else'figures'
        links=', '.join(f'[{x}](data.qmd#{x.lower()})' for x in refs)
        pages[page].append(f'## Figure {num} {{#figure-{num.lower()}}}\n\n![](../assets/figures/{fn}.png){{fig-alt="Figure {num}; detailed caption follows." loading="lazy"}}\n\n{md}\n\n**Table links:** {links}.')
for page,parts in pages.items():(repo/f'book/chapters/{page}.qmd').write_text('\n\n'.join(parts)+'\n')
for path in [dest/'figure_data_map.csv',paper/'figure_data_map.csv']:
    with path.open('w') as f:w=csv.DictWriter(f,fieldnames=['figure','supplementary_tables']);w.writeheader();w.writerows(mapping)
data='''# Data and downloads

All **42 supplementary tables**, **7 main figures** and **18 supplementary figures** match manuscript version **16**. The primary results use the reconstructed ROSMAP inputs; the new prioritization and stability diagnostics retain those fits and primary analyses.

- [Supplementary tables](../downloads/supplementary_tables.zip), [SHA-256 table manifest](../downloads/table_manifest.csv) and [figure-to-table map](../downloads/figure_data_map.csv).
- [Manuscript](../downloads/manuscript.pdf) and [supplement](../downloads/supplement.pdf).
- [Verified reconstructed brain inputs and workflow](../downloads/reconstructed_brain_inputs.zip).
- [Independent LD benchmark](../downloads/independent_LD_benchmark.tar.gz): compact inputs, archived inference source, code and completed results.
- [Additional diagnostic inputs](../downloads/revision16_analysis_inputs.rds) and [source manifest](../downloads/revision16_input_sources.csv).
- [Package-compatible brain scorer](../downloads/brain_GAM6_release_0.1.0.rds) and [model description](../downloads/brain_scorer_model_card.txt).
- [Tutorial scripts](../downloads/tutorial_scripts.zip), [reproduction scripts](../downloads/reproduction_scripts.zip) and [validation record](../downloads/book_validation.txt).

Rendering figures repeats presentation from saved data; it does not rerun native inference. The brain input archive contains the verified reconstructed numerical inputs and nucleus selections. Source verification separately reads the nine downloaded ROSMAP files. Native fit caches are excluded. See [reproduction](reproduce.qmd).

Tables preserve original machine identifiers: `Excit` means excitatory neurons, `Inhib` inhibitory neurons, `Astro` astrocytes, `Oligo` oligodendrocytes, `OPC` oligodendrocyte precursor cells, `Immune` immune cells and `Vascular` vascular cells. Undefined precision is not zero. The additional diagnostic table definitions are in Supplementary Section 5.

## Table inventory

'''
for p in tables:
    num=re.match(r'S\d+',p.name)[0];title=re.sub(r'^S\d+_|\.csv$','',p.name).replace('_',' ')
    refs=[x['figure'] for x in mapping if num in x['supplementary_tables'].split(';')]
    data+=f'### {num} — {title} {{#{num.lower()}}}\n\n`{p.name}` · {p.stat().st_size/1024/1024:.2f} MB. '+('Supports figures '+', '.join(refs)+'.' if refs else 'Additional analysis and provenance data.')+'\n\n'
(repo/'book/chapters/data.qmd').write_text(data)
release={'manuscript_version':16,'package_version':'0.1.0','table_count':len(tables),'figure_count':len(figs),'primary_input_version':'reconstructed_20261018'}
(repo/'reproducibility.json').write_text(json.dumps(release,indent=2)+'\n')
print('Synchronized version 16: 42 tables and 25 figures.')
