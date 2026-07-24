import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import glob, os, re
for p in glob.glob(os.path.expanduser('~/.fonts/Arial*.ttf')):
    try: fm.fontManager.addfont(p)
    except Exception: pass
plt.rcParams['font.family']='sans-serif'
plt.rcParams['font.sans-serif']=['Arial','Liberation Sans','DejaVu Sans']
plt.rcParams['mathtext.default']='regular'   # superscripts render in Arial
plt.rcParams['pdf.fonttype']=42; plt.rcParams['ps.fonttype']=42
plt.rcParams['svg.fonttype']='none'
CORAL='#F4633A'

def _supplus(txt):
    # turn '+' into a superscript plus at upper-right of the preceding char
    return txt.replace('+', r'$^{+}$')

def umap_right_legend(adata, key, outpath, title='Cell types', size=6, dpi=140, square=False):
    import scanpy as sc, numpy as np, matplotlib
    from matplotlib.lines import Line2D
    a=adata
    vc=a.obs[key].value_counts()
    cats=list(a.obs[key].cat.categories) if hasattr(a.obs[key],'cat') else sorted(a.obs[key].unique())
    uns_key=f'{key}_colors'
    fig,ax=plt.subplots(figsize=(6.6,5.4) if square else (7.6,5.2))
    xy=a.obsm['X_umap']
    if uns_key in a.uns: palette=list(a.uns[uns_key])
    else:
        cmap=matplotlib.cm.get_cmap('tab20'); palette=[matplotlib.colors.to_hex(cmap(i%20)) for i in range(len(cats))]
    col={c:palette[i%len(palette)] for i,c in enumerate(cats)}
    for c in cats:
        m=(a.obs[key]==c).values
        ax.scatter(xy[m,0],xy[m,1],s=size,c=col[c],linewidths=0,rasterized=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(''); ax.set_ylabel('')
    for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(1.0); s.set_color('black')
    if square: ax.set_box_aspect(1)   # force square plotting box
    handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor=col[c],markersize=9,
             label=f'{_supplus(str(c))} ({int(vc.get(c,0))})') for c in cats]
    leg=ax.legend(handles=handles,loc='center left',bbox_to_anchor=(1.02,0.5),
              frameon=False,fontsize=11,handletextpad=0.4,labelspacing=0.7,
              title=title,title_fontsize=14)
    leg.get_title().set_color(CORAL); leg.get_title().set_ha('left')
    plt.tight_layout(); plt.savefig(outpath,dpi=dpi,bbox_inches='tight'); plt.close()
