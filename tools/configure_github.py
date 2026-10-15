"""Set public links. Does not create a repository or publish anything."""
import argparse,re,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('owner');p.add_argument('--repo',default='CDSeq2');a=p.parse_args()
if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9-]*',a.owner):p.error('Use a GitHub user or organization name.')
if not re.fullmatch(r'[A-Za-z0-9_.-]+',a.repo):p.error('Invalid repository name.')
root=Path(__file__).resolve().parents[1]
url=f'https://github.com/{a.owner}/{a.repo}'
f=root/'book/_quarto.yml';s=f.read_text();s=re.sub(r'^  repo-url:.*\n','',s,flags=re.M)
s=s.replace('  repo-actions: [issue]\n','')
s=s.replace('  title: "CDSeq2"',f'  title: "CDSeq2"\n  repo-url: {json.dumps(url)}\n  repo-actions: [issue]')
f.write_text(s)
v=root/'book/_variables.yml';v.write_text(f'github_owner: {json.dumps(a.owner)}\nrepository: {json.dumps(a.repo)}\nrepository_url: {json.dumps(url)}\n')
f=root/'DESCRIPTION';s=f.read_text();s=re.sub(r'^(URL|BugReports):.*\n','',s,flags=re.M)
s+=f'URL: {url}\nBugReports: {url}/issues\n';f.write_text(s)
print(f'Configured {url}; no files uploaded. Enable GitHub Pages with source GitHub Actions after pushing.')
