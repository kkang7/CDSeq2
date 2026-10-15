"""Portable entry point. From repository root: python3 reproduction/reproduce.py --help"""
from pathlib import Path
import argparse,subprocess,zipfile,tarfile,shutil,csv,hashlib,sys
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('part',choices=['tables','figures','independent-prepare','independent-run','independent-figures'])
p.add_argument('--out',type=Path,default=root/'work')
p.add_argument('--workers',type=int,default=4)
p.add_argument('--python',default=sys.executable,help='Python with reportlab and pypdfium2 for schematic')
a=p.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
if a.workers<1:p.error('--workers must be positive')
def run(*cmd):subprocess.run([str(x) for x in cmd],check=True)
def tables():
    dest=out/'paper';dest.mkdir(exist_ok=True)
    with zipfile.ZipFile(root/'book/downloads/supplementary_tables.zip') as z:
        for item in z.infolist():
            path=dest/item.filename;payload=z.read(item)
            if path.exists() and path.read_bytes()!=payload:
                raise SystemExit(f'Refusing to overwrite changed table: {path}; choose another --out.')
            path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(payload)
    for row in csv.DictReader((root/'book/downloads/table_manifest.csv').open()):
        f=dest/'supplementary_tables'/row['file']
        assert hashlib.sha256(f.read_bytes()).hexdigest()==row['sha256'],f
    shutil.copy2(root/'reproduction/paper/code/cell_type_labels.R',dest/'cell_type_labels.R')
    (dest/'figures').mkdir(exist_ok=True)
    print('Verified all 35 supplementary tables:',dest)
    return dest
def independent():
    dest=out/'independent_ld'
    if not dest.exists():
        with tarfile.open(root/'book/downloads/independent_LD_benchmark.tar.gz') as t:
            # All archive members must be ordinary files/directories within destination.
            for member in t.getmembers():
                target=(out/member.name).resolve()
                if not target.is_relative_to(out) or not(member.isfile() or member.isdir()):
                    raise SystemExit('Invalid benchmark archive member')
            t.extractall(out)
    print('Independent benchmark:',dest)
    return dest
if a.part=='tables':tables()
elif a.part=='figures':
    dest=tables();code=root/'reproduction/paper/code'
    run('Rscript','--vanilla',code/'prepare_expression_prioritization.R','-',dest)
    run('Rscript','--vanilla',code/'eqtl_accuracy_metrics.R',dest)
    for name in ['revised_figures.R','additional_figures.R','application_figures.R',
                 'eqtl_accuracy_figures.R','figure_eqtl_rankings.R']:
        run('Rscript','--vanilla',code/name,dest)
    run(a.python,code/'draw_workflow.py',dest)
    ld=independent();run('Rscript','--vanilla',ld/'code/figures.R',ld,dest/'figures')
    run(a.python,code/'verify_eqtl_accuracy.py',dest)
    run(a.python,code/'verify_expression_prioritization.py',dest)
    assert len(list((dest/'figures').glob('*.png')))==21
    print('Rebuilt all 6 main and 15 supplementary figures:',dest/'figures')
elif a.part=='independent-prepare':independent()
elif a.part=='independent-figures':
    ld=independent();run('Rscript','--vanilla',ld/'code/figures.R',ld)
elif a.part=='independent-run':
    ld=independent();run(sys.executable,ld/'code/run_all.py',ld,a.workers)
    for name in ['analyze.R','export_design.R','summarize_null.R','verify.R']:
        run('Rscript','--vanilla',ld/'code'/name,ld)
    run(sys.executable,ld/'code/verify_tables.py',ld)
    run('Rscript','--vanilla',ld/'code/verify_pooled_curves.R',ld)
    run('Rscript','--vanilla',ld/'code/figures.R',ld)
