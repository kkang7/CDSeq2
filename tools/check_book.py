"""Check the source asset inventory and rendered local navigation."""
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit,unquote
import re,csv,zipfile,hashlib,argparse,json
root=Path(__file__).resolve().parents[1];book=root/'book'
p=argparse.ArgumentParser();p.add_argument('--source-only',action='store_true');a=p.parse_args()
errors=[]
release=json.loads((root/'reproducibility.json').read_text())
def need(ok,msg):
    if not ok:errors.append(msg)
manifest=list(csv.DictReader((book/'downloads/table_manifest.csv').open()))
need(len(manifest)==release['table_count'],'Supplementary table count differs from release metadata')
with zipfile.ZipFile(book/'downloads/supplementary_tables.zip') as z:
    for row in manifest:
        payload=z.read('supplementary_tables/'+row['file'])
        need(hashlib.sha256(payload).hexdigest()==row['sha256'],'Changed table '+row['file'])
need(len(list((book/'assets/figures').glob('*.png')))==release['figure_count'],'Figure count differs from release metadata')
for f in root.rglob('*'):
    if f.is_file() and not any(x in f.relative_to(root).parts for x in ['work','_book','.quarto','.git']):
        need(f.stat().st_size<100*1024*1024,'File exceeds 100 MiB: '+str(f.relative_to(root)))
class Links(HTMLParser):
    def __init__(self):super().__init__();self.links=[];self.ids=set()
    def handle_starttag(self,tag,attrs):
        at=dict(attrs)
        if at.get('id'):self.ids.add(at['id'])
        for attr in ['href','src']:
            if at.get(attr):self.links.append(at[attr])
if not a.source_only:
    rendered=book/'_book';need((rendered/'index.html').exists(),'Render book first')
    parsed={}
    for f in rendered.rglob('*.html'):
        x=Links();x.feed(f.read_text());parsed[f.resolve()]=x
    need(len(parsed)>=17,'Expected all book chapters')
    for f,x in parsed.items():
        for link in x.links:
            u=urlsplit(link)
            if u.scheme or u.netloc or link.startswith('javascript:'):continue
            path=(f.parent/unquote(u.path)).resolve() if u.path else f
            if path.is_dir():path=path/'index.html'
            need(path.exists(),f'Missing local link in {f.name}: {link}')
            if u.fragment and path in parsed:
                need(unquote(u.fragment) in parsed[path].ids,f'Missing anchor in {f.name}: {link}')
    for row in manifest:
        need(row['table'].lower() in parsed.get((rendered/'chapters/data.html').resolve(),Links()).ids,'Missing table anchor '+row['table'])
if errors:raise SystemExit('\n'.join(errors))
print(f"PASS: {len(manifest)} verified tables, {release['figure_count']} figure assets, GitHub file-size limits"+
      ('' if a.source_only else ', rendered pages, local links, downloads and table anchors'))
