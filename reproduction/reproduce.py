"""Portable entry point. From repository root: python3 reproduction/reproduce.py --help"""
from pathlib import Path
import argparse,subprocess,zipfile,tarfile,shutil,csv,hashlib,sys,os,json
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('part',choices=['tables','figures','independent-prepare','independent-run','independent-figures','diagnostics','brain-prepare','brain-run','benchmark-library'])
p.add_argument('--out',type=Path,default=root/'work')
p.add_argument('--workers',type=int,default=4)
p.add_argument('--python',default=sys.executable,help='Python with reportlab and pypdfium2 for schematic')
a=p.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
if a.workers<1:p.error('--workers must be positive')
env=os.environ.copy()
def run(*cmd):subprocess.run([str(x) for x in cmd],check=True,env=env)
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
    print('Verified supplementary tables:',dest)
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
def benchmark_library():
    ld=independent();lib=out/'benchmark-library';lib.mkdir(exist_ok=True)
    run('R','CMD','INSTALL','--library='+str(lib),ld/'source_snapshot/CDSeq2_0.0.0.9000.tar.gz')
    env['R_LIBS_USER']=str(lib)+os.pathsep+env.get('R_LIBS_USER','')
    return lib

def brain():
    dest=out/'ReconstructedData'
    if not dest.exists():
        with zipfile.ZipFile(root/'book/downloads/reconstructed_brain_inputs.zip') as z:
            for name in z.namelist():
                if not (out/name).resolve().is_relative_to(out):raise SystemExit('Invalid brain archive path')
            z.extractall(out)
    # The plotting code sits next to the copied brain directory.
    code=dest/'reproduction/paper/code';shutil.copytree(root/'reproduction/paper/code',code,dirs_exist_ok=True)
    print('Verified reconstructed-input workflow:',dest)
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
    run('Rscript','--vanilla',code/'prioritization_audit_figures.R',dest,'--figures-only')
    assert len(list((dest/'figures').glob('*.png')))==25
    print('Rebuilt all 7 main and 18 supplementary figures:',dest/'figures')
elif a.part=='independent-prepare':independent()
elif a.part=='independent-figures':
    ld=independent();run('Rscript','--vanilla',ld/'code/figures.R',ld)
elif a.part=='independent-run':
    benchmark_library();ld=independent();run(sys.executable,ld/'code/run_all.py',ld,a.workers)
    for name in ['analyze.R','export_design.R','summarize_null.R','verify.R']:
        run('Rscript','--vanilla',ld/'code'/name,ld)
    run(sys.executable,ld/'code/verify_tables.py',ld)
    run('Rscript','--vanilla',ld/'code/verify_pooled_curves.R',ld)
    run('Rscript','--vanilla',ld/'code/figures.R',ld)

elif a.part=='benchmark-library':benchmark_library()
elif a.part=='brain-prepare':brain()
elif a.part=='brain-run':
    benchmark_library();b=brain()
    for name in ['00_preflight.R','01_prepare.R','02_deconvolution_and_gam7.R','03_gam6.R','04_continuous_phenotype.R','05_eqtl.R','06_standard_bayesprism.R','07_paper_outputs.R']:
        run('Rscript','--vanilla',b/'reproduction/brain'/name,b/'config.R')
elif a.part=='diagnostics':
    dest=tables();ld=independent();code=root/'reproduction/paper/code';inp=root/'book/downloads/revision16_analysis_inputs.rds'
    run('Rscript','--vanilla',code/'prioritization_audit.R',inp,ld/'results',dest/'supplementary_tables')
    run('Rscript','--vanilla',code/'null_stability_audit.R',inp,dest/'supplementary_tables')
    run('Rscript','--vanilla',code/'prioritization_audit_figures.R',dest)
    run(sys.executable,code/'verify_prioritization_audit.py',dest)
