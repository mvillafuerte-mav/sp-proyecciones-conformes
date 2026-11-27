"""Geometría de la proyección estereográfica y panel de retículas existentes."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pymupdf

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family': 'DejaVu Serif', 'font.size': 11})
phi = np.deg2rad(28)
lam = np.deg2rad(25)
S = np.array([0., 0., -1.])
O = np.zeros(3)
P = np.array([np.cos(phi)*np.cos(lam), np.cos(phi)*np.sin(lam), np.sin(phi)])
A = np.array([0., 0., P[2]])
Q = np.array([P[0]/(1+P[2]), P[1]/(1+P[2]), 0])
assert np.isclose(np.linalg.norm(Q), np.cos(phi)/(1+np.sin(phi)))
assert np.linalg.norm(np.cross(Q-S, P-S)) < 1e-12
fig = plt.figure(figsize=(10, 4.8), facecolor='white')
ax = fig.add_subplot(121, projection='3d')
u, v = np.meshgrid(np.linspace(0, 2*np.pi, 40), np.linspace(0, np.pi, 25))
ax.plot_wireframe(np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v),
                  color='#c6cdd4', linewidth=.35, alpha=.35)
for triangle,color in [([S,A,P], '#e66d24'), ([S,O,Q], '#1769ab')]:
    ax.add_collection3d(Poly3DCollection([triangle],facecolors=color,alpha=.27,edgecolors=color,linewidths=1.5))
for pts,color in [([S,P], '#c84917'),([O,P], '#4a254f'),([O,A], '#222222')]:
    points=np.array(pts);ax.plot(*points.T,color=color,lw=2)
for name,point in [('S',S),('O',O),('P',P),('A',A),('Q',Q)]:
    ax.scatter(*point,color='#111111',s=18)
    ax.text(*(point+np.array([.03,.03,.02])),name,fontsize=12,fontweight='bold')
ax.text(*((O+Q)/2+np.array([0,0,-.1])),'r',color='#1769ab')
ax.text(*((S+O)/2+np.array([-.13,0,0])),'R',color='#1769ab')
ax.text(*((A+P)/2+np.array([0,0,.08])),r'$R\cos\varphi$',color='#a94315')
ax.text(-.37,0,.15,r'$R\sin\varphi$',color='#4a254f',fontsize=10)
ax.set_box_aspect((1,1,1),zoom=1.4);ax.set_proj_type('ortho');ax.view_init(elev=18,azim=-63)
ax.set_axis_off();ax.set_title('(a) Esfera y plano ecuatorial',pad=0)

ax2=fig.add_subplot(122)
p=np.array([np.cos(phi),np.sin(phi)]);s=np.array([0,-1]);o=np.array([0,0]);a=np.array([0,p[1]]);q=np.array([np.cos(phi)/(1+np.sin(phi)),0])
t=np.linspace(-np.pi/2,np.pi/2,300)
ax2.plot(np.cos(t),np.sin(t),color='#b3bdc8',lw=1.2)
ax2.axhline(0,color='#a3adb7',lw=1);ax2.axvline(0,color='#a3adb7',lw=1)
ax2.fill(*np.array([s,a,p]).T,color='#e66d24',alpha=.17)
ax2.fill(*np.array([s,o,q]).T,color='#1769ab',alpha=.27)
for pts,color in [([s,p],'#c84917'),([o,p],'#4a254f'),([a,p],'#c84917'),([s,a],'#c84917'),([o,q],'#1769ab')]:
    points=np.array(pts);ax2.plot(*points.T,color=color,lw=1.8)
for name,point,off in [('S',s,(-.08,-.09)),('O',o,(-.08,.04)),('A',a,(-.08,.04)),('Q',q,(.02,-.11)),('P',p,(.03,.02))]:
    ax2.scatter(*point,color='#111111',s=18);ax2.text(*(point+off),name,fontweight='bold')
arc=np.linspace(0,phi,50);ax2.plot(.25*np.cos(arc),.25*np.sin(arc),color='#4a254f');ax2.text(.27,.06,r'$\varphi$',color='#4a254f')
ax2.text(q[0]/2,-.12,'r',color='#1769ab');ax2.text(-.14,-.53,'R',color='#1769ab')
ax2.text(p[0]/2,p[1]+.08,r'$R\cos\varphi$',color='#a94315',ha='center')
ax2.text(-.52,.21,r'$R\sin\varphi$',color='#4a254f')
ax2.text(.48,.31,'R',color='#4a254f')
ax2.text(.45,-.79,r'$SO=R,\quad SA=R+R\sin\varphi$',ha='center',fontsize=10)
ax2.set_aspect('equal');ax2.set_xlim(-.58,1.15);ax2.set_ylim(-1.16,1.12);ax2.axis('off');ax2.set_title('(b) Corte por O, S y P')
fig.tight_layout()
geometry=ROOT/'wiki/recursos/construccion_estereografica.pdf'
fig.savefig(geometry,bbox_inches='tight');plt.close(fig)

g=pymupdf.open(geometry)
ret=pymupdf.open(ROOT/'wiki/figuras/fig1_reticula.pdf')
out=pymupdf.open();page=out.new_page(width=720,height=620)
page.show_pdf_page(pymupdf.Rect(0,0,720,350),g,0)
page.show_pdf_page(pymupdf.Rect(35,355,685,620),ret,0)
out.save(ROOT/'wiki/figuras/fig1_geometria_reticula.pdf')

page.get_pixmap(matrix=pymupdf.Matrix(3, 3)).save(ROOT/"wiki/figuras/fig1_geometria_reticula.png")
