"""Vector biological-cartoon workflow. Schematic; no data encoded."""
from pathlib import Path
import sys,math,random
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor,Color,white
out=Path(sys.argv[1] if len(sys.argv)>1 else '.')/'figures';out.mkdir(exist_ok=True)
W,H=1000,650;c=canvas.Canvas(str(out/'Figure_1_workflow.pdf'),pagesize=(W,H))
ink='#173343';muted='#546B78';blue='#237EAB';orange='#DC8E4B';green='#469E89';purple='#9675B2'
def color(x):return HexColor(x)
def text(x,y,s,size=13,bold=False,col=ink,center=False):
 c.setFillColor(color(col));c.setFont('Helvetica-Bold' if bold else 'Helvetica',size);(c.drawCentredString if center else c.drawString)(x,y,s)
def lines(x,y,ss,size=12,gap=17,**kw):
 for i,s in enumerate(ss):text(x,y-i*gap,s,size,**kw)
def line(x1,y1,x2,y2,col=muted,w=1.3):c.setStrokeColor(color(col));c.setLineWidth(w);c.line(x1,y1,x2,y2)
def arrow(x1,y1,x2,y2):
 line(x1,y1,x2,y2,blue,2);a=math.atan2(y2-y1,x2-x1);p=c.beginPath();p.moveTo(x2,y2);p.lineTo(x2-8*math.cos(a-.5),y2-8*math.sin(a-.5));p.lineTo(x2-8*math.cos(a+.5),y2-8*math.sin(a+.5));p.close();c.setFillColor(color(blue));c.drawPath(p,fill=1,stroke=0)
def cell(x,y,r,col,kind=0):
 c.setFillColor(color(col));c.setStrokeColor(color(col));c.setLineWidth(1.5)
 if kind==0: # neuron with branching neurites
  for ang in [0.2,1.1,2.3,3.4,4.4,5.3]:
   dx,dy=math.cos(ang),math.sin(ang);end=(x+dx*r*2.1,y+dy*r*2.1);c.line(x,y,*end)
   for side in [-.5,.5]:c.line(*end,end[0]+r*.65*math.cos(ang+side),end[1]+r*.65*math.sin(ang+side))
  p=c.beginPath();p.moveTo(x,y+r);p.lineTo(x-r,y-r*.7);p.lineTo(x+r,y-r*.7);p.close();c.drawPath(p,fill=1,stroke=0)
 elif kind==1: # stellate glial cartoon
  p=c.beginPath()
  for j in range(20):
   a=j*math.pi/10;rr=r*(1.75 if j%2==0 else .6);(p.moveTo if j==0 else p.lineTo)(x+rr*math.cos(a),y+rr*math.sin(a))
  p.close();c.drawPath(p,fill=1,stroke=0)
 else:c.circle(x,y,r,fill=1,stroke=0)
 c.setFillColor(Color(1,1,1,.75));c.circle(x,y,r*.32,fill=1,stroke=0)
def donor(x,y,col=muted):
 c.setFillColor(color(col));c.circle(x,y+15,5.2,fill=1,stroke=0);c.roundRect(x-8,y-4,16,13,4,fill=1,stroke=0)
def heatmap(x,y,nr,nc,size=9):
 rng=random.Random(32+int(x))
 for r in range(nr):
  for k in range(nc):
   z=rng.random();c.setFillColor(Color(.90-z*.62,.96-z*.37,.98-z*.22));c.rect(x+k*size,y+r*size,size-1,size-1,stroke=0,fill=1)
def axes(x,y,w,h):line(x,y,x+w,y,'#8BA0AD');line(x,y,x,y+h,'#8BA0AD')

# Three-stage workflow, with downstream applications retained below.
text(500,620,'CDSeq2: from human brain tissue to prioritized expression profiles',22,True,center=True)
for x,title,label,col in [(20,'Input','A',purple),(350,'CDSeq2','B',blue),(680,'Output','C',green)]:
 c.setFillColor(color('#F5F8FA'));c.setStrokeColor(color('#D8E3EA'));c.roundRect(x,235,300,350,12,fill=1,stroke=1)
 c.setFillColor(color(col));c.circle(x+23,560,12,fill=1,stroke=0);text(x+23,555,label,14,True,col='#FFFFFF',center=True);text(x+44,555,title,17,True)
arrow(322,410,347,410);arrow(652,410,677,410)
# Brain cartoon with cortical folds and a magnified mixed-cell tissue sample.
c.saveState();c.translate(65,475);c.scale(.78,.78)
brain=c.beginPath();brain.moveTo(-39,-13)
brain.curveTo(-52,-4,-49,13,-42,19);brain.curveTo(-45,35,-29,46,-17,44)
brain.curveTo(-7,57,12,53,20,45);brain.curveTo(37,47,48,34,44,23)
brain.curveTo(57,10,49,-6,40,-10);brain.curveTo(41,-24,23,-32,10,-25)
brain.curveTo(-6,-35,-18,-24,-20,-24);brain.curveTo(-33,-28,-45,-24,-39,-13);brain.close()
c.setFillColor(color('#E5C8DB'));c.setStrokeColor(color(purple));c.setLineWidth(2);c.drawPath(brain,fill=1,stroke=1)
# Cerebellum and short brainstem.
c.setFillColor(color('#D2B2CC'));c.ellipse(14,-38,39,-19,fill=1,stroke=1)
q=c.beginPath();q.moveTo(8,-24);q.curveTo(6,-34,9,-42,15,-45);q.lineTo(21,-43);q.curveTo(15,-35,15,-29,19,-25);c.drawPath(q,fill=1,stroke=1)
for pts in [(-36,14,-23,24,-18,12,-27,3),(-28,35,-13,41,-8,28,-15,21),(-11,43,3,32,-4,21,5,15),(14,37,27,33,17,20,29,16),(33,27,43,16,29,6,37,-3),(-35,-10,-24,-19,-12,-7,-19,1),(-8,11,4,6,-7,-4,3,-14),(9,20,18,9,8,-1,18,-9),(20,-15,25,-24,34,-16,34,-11)]:
 q=c.beginPath();q.moveTo(*pts[:2]);q.curveTo(*pts[2:]);c.drawPath(q,stroke=1,fill=0)
for yy in [-26,-30,-34]:line(21,yy,33,yy,purple,1)
c.restoreState()
line(98,497,120,509,purple,1);line(99,460,120,438,purple,1)
# Organic tissue outline, filled with colored neurons and glia.
pth=c.beginPath();pth.moveTo(124,441);pth.curveTo(103,464,117,512,145,522);pth.curveTo(171,537,194,523,219,529);pth.curveTo(252,534,280,512,287,487);pth.curveTo(301,464,275,438,251,437);pth.curveTo(223,423,204,436,181,429);pth.curveTo(158,420,130,424,124,441);pth.close()
c.setFillColor(color('#F3EAF1'));c.setStrokeColor(color('#C9ABBE'));c.drawPath(pth,fill=1,stroke=1)
cell_colors=[blue,orange,purple,green]
for xx,yy,co,kind in [(146,487,blue,0),(197,501,purple,1),(253,496,orange,0),(164,451,green,2),(224,458,blue,0),(271,460,purple,1)]:cell(xx,yy,7,co,kind)
text(196,405,'Human brain tissue',16,True,center=True);text(196,385,'Mixture of different cell types',13,center=True)
arrow(170,370,170,351)
# Bulk assay and external annotated reference as two complementary inputs.
# Individual short reads interleave colors from different cell-type origins.
for row in range(5):
 for col in range(3):
  xx=48+col*23+(row%2)*3;yy=333-row*7;co=cell_colors[(row+col)%4]
  line(xx,yy,xx+15,yy,co,2.8)
  for offset in [4,8,12]:line(xx+offset,yy-1.6,xx+offset,yy+1.6,co,.8)
lines(82,280,['Bulk RNA','counts'],14,center=True)
for i,(co,kind) in enumerate(zip(cell_colors,[0,0,1,2])):cell(196+i*27,319,4.4,co,kind)
lines(235,280,['Annotated single-cell /','single-nucleus reference'],12.5,center=True)
# Combined package: inference followed by optional reliability scoring.
c.setFillColor(color('#E5F0F7'));c.setStrokeColor(color(blue));c.roundRect(375,424,250,99,9,fill=1,stroke=1)
text(500,498,'CDSeq2 inference',18,True,center=True)
lines(500,475,['Reference-informed Bayesian model','Individual expression + RNA fractions'],13,center=True)
arrow(500,418,500,392)
c.setFillColor(color('#E3F1EC'));c.setStrokeColor(color(green));c.roundRect(375,283,250,103,9,fill=1,stroke=1)
text(500,360,'Gene-type prioritization (optional)',14,True,center=True)
lines(500,337,['GAM6: six observable support features','Rank gene-cell-type recoverability','Choose the proportion to retain'],13,center=True)
text(500,252,'No phenotype or genotype input to scoring',12,center=True)
# Colors match each cell type to the corresponding matrix column.
for i,(co,nm,kind) in enumerate(zip(cell_colors,['Excitatory\nneurons','Inhibitory\nneurons','Astrocytes','Oligodendrocytes'],[0,0,1,2])):
 cell(715+i*73,510,6,co,kind);lines(715+i*73,487,nm.split('\n'),9,gap=10,bold=True,center=True)
text(830,458,'Individual-level cell-type expression',14,True,center=True)
rng=random.Random(42)
for k in range(3):
 x=710+k*83
 text(x+28,432,'Donor '+str(k+1),12,center=True)
 for row in range(6):
  for col,co in enumerate(cell_colors):
   base=color(co);z=.22+.72*rng.random();c.setFillColor(Color(1-z*(1-base.red),1-z*(1-base.green),1-z*(1-base.blue)))
   c.rect(x+col*14,336+row*14,13,13,stroke=0,fill=1)
 for col,co in enumerate(cell_colors):c.setFillColor(color(co));c.rect(x+col*14,423,13,4,stroke=0,fill=1)
text(830,317,'Rows: genes | Columns: cell types',12,center=True)
text(830,296,'Full estimates + prioritized gene-type pairs',12.5,center=True)
for k in range(3):
 xx=710+k*83
 for width,co in zip([23,10,13,10],cell_colors):c.setFillColor(color(co));c.rect(xx,270,width,9,fill=1,stroke=0);xx+=width
text(830,250,'RNA fractions for each individual',12,center=True)
arrow(830,231,830,207)
# Downstream cartoons: continuous phenotype scatter and genotype-expression plot.
c.setFillColor(color('#F5F8FA'));c.setStrokeColor(color('#D8E3EA'));c.roundRect(20,35,960,168,12,fill=1,stroke=1)
c.setFillColor(color(orange));c.circle(43,178,12,fill=1,stroke=0);text(43,173,'D',14,True,col='#FFFFFF',center=True);text(65,173,'Downstream analysis',16,True)
text(168,143,'Continuous phenotypes',15,True,center=True);lines(168,120,['Relate expression to','quantitative traits'],13,center=True)
axes(320,70,140,79);rng=random.Random(3)
for i in range(23):
 xx=324+i*5.6;yy=82+i*2.2+rng.uniform(-12,12);c.setFillColor(color(blue));c.circle(xx,yy,2.3,fill=1,stroke=0)
line(324,82,451,133,blue,1.7);text(390,50,'Phenotype',12,center=True)
text(630,143,'Cell-type-specific eQTLs',15,True,center=True)
for i in range(14):
 yy=122-i*4;x1=562+10*math.sin(i*.65);x2=562-10*math.sin(i*.65);line(x1,yy,x2,yy,purple,1)
 if i:line(prev[0],prev[1],x1,yy,blue,1.5);line(prev[2],prev[1],x2,yy,orange,1.5)
 prev=(x1,yy,x2)
lines(653,104,['Genetic regulation','within cell types'],13,center=True)
axes(780,70,151,79)
for g in range(3):
 for i in range(7):
  xx=805+g*49+rng.uniform(-5,5);yy=84+g*22+rng.uniform(-9,9);c.setFillColor(color(blue));c.circle(xx,yy,2.3,fill=1,stroke=0)
 text(805+g*49,56,str(g),10,center=True)
text(855,42,'Genotype dosage',11,center=True)
text(500,14,'Schematic: complete expression estimates and optional gene-type priorities support downstream analysis.',11,col=muted,center=True)
c.showPage();c.save()
try:
 import pypdfium2 as pdfium
 doc=pdfium.PdfDocument(str(out/'Figure_1_workflow.pdf'));doc[0].render(scale=3).to_pil().save(out/'Figure_1_workflow.png')
except ImportError:pass
