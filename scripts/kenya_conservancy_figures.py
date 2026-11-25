"""Render the three 200 m conservancy comparison figures from the results CSV."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import patheffects as pe
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--stats", type=Path,
                    default=root / "outputs/kenya_conservancies/kenya_conservancy_stats.csv")
parser.add_argument("--outdir", type=Path, default=root / "outputs/kenya_conservancies/figs")
args = parser.parse_args()
output = args.outdir
data = pd.read_csv(args.stats)
sample = data[data.category.str.contains("conservancy", case=False, na=False)].copy()
stats=["200m_mean_p_km2","200m_p99_p_km2","200m_max_p_km2"]
names=["Mean","99th percentile","Maximum"]
ring="200m_mean_5km_outside_p_km2"
if sample.empty or sample[[*stats, ring]].isna().any().any():
    parser.error("The conservancy sample must have complete 200 m statistics.")
output.mkdir(parents=True, exist_ok=True)
counts=[int(sample[c].lt(50).sum()) for c in stats]
rho=spearmanr(sample[stats[0]],sample[ring])[0]
opc=sample.loc[sample.name.eq("Ol Pejeta Conservancy")].iloc[0]

ink="#243332"
slate="#515e5a"
teal="#3c716d"
dark_teal="#245d59"
light_teal="#afc9c2"
mist="#ebf0ec"
grid="#e1e7e2"
accent="#ae541f"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8.5,"pdf.fonttype":42,
                     "text.color":ink,"axes.labelcolor":ink,"xtick.color":slate,"ytick.color":slate,
                     "axes.spines.top":False,"axes.spines.right":False,"axes.edgecolor":"#86918b",
                     "axes.linewidth":.55})

def figure(h):
    return plt.figure(figsize=(125/25.4,h/25.4),facecolor="white")

def place(fig,x,y,text,**kwargs):
    w,h=fig.get_size_inches()*25.4
    return fig.text(x/w,y/h,text,**kwargs)

def save(fig,stem):
    fig.savefig(output/(stem+".pdf"))
    fig.savefig(output/(stem+".png"),dpi=300)
    plt.close(fig)

def labels(a,axis,ticks):
    values=np.log1p(ticks)
    names=[f"{x:,}" for x in ticks]
    (a.set_xticks if axis=="x" else a.set_yticks)(values,names)

# Share of conservancies below the 50 people/km² reference threshold.
fig=figure(34)
ax=fig.add_axes([.295,.05,.530,.90])
ax.set(xlim=(0,100),ylim=(-.6,2.6))
ax.set_axis_off()
for i,(name,n) in enumerate(zip(names,counts)):
    y=2-i
    pct=100*n/len(sample)
    ax.barh(y,100,height=.48,color=mist,linewidth=0,zorder=1)
    ax.barh(y,pct,height=.48,color=teal,linewidth=0,zorder=2)
    ax.text(-5.2,y,name,ha="right",va="center",fontsize=9.2,weight="medium")
    ax.text(106,y+.040,f"{pct:.0f}%",ha="left",va="center",
            fontsize=12,weight="bold",color=dark_teal)
    ax.text(106,y-.29,f"({n}/{len(sample)})",ha="left",va="center",color=slate,fontsize=8)
save(fig,"01_threshold_comparison")

# Median and interquartile range for each density statistic.
fig=figure(48)
ax=fig.add_axes([.245,.195,.552,.755])
lim=np.log1p(max(430, 1.1 * sample[stats].quantile(.75).max(), 1.1 * opc[stats].max()))
ax.set(xlim=(0,lim),ylim=(-.57,2.68))
labels(ax,"x",[0,5,10,25,50,100,250])
ax.tick_params(axis="x",length=2,pad=2)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_axisbelow(True)
ax.grid(axis="x",color=grid,lw=.48)
t=np.log1p(50)
ax.axvline(t,color="#6a756d",ls=(0,(3,2)),lw=.72,zorder=2)
for i,(name,col) in enumerate(zip(names,stats)):
    y=2-i
    q1,med,q3=sample[col].quantile([.25,.5,.75])
    a,m,b=np.log1p([q1,med,q3])
    ax.plot([a,b],[y,y],lw=7.8,color=light_teal,solid_capstyle="round",zorder=3)
    ax.plot(m,y,marker="o",markersize=5,color=dark_teal,mec="white",mew=.6,zorder=4)
    ax.scatter(np.log1p(opc[col]),y-.235,marker="D",s=21,facecolor=accent,edgecolor="white",lw=.55,zorder=5)
    ax.text(-.145,y,name,ha="right",va="center",fontsize=9.2,weight="medium",clip_on=False)
    ax.annotate(f"{med:.0f}",(m,y),xytext=(0,6),textcoords="offset points",
                ha="center",va="bottom",weight="bold",fontsize=8.5,color=dark_teal,zorder=6)
place(fig,65.1,2.4,"Inside density (people/km²)",ha="center",fontsize=8.4)
legend = [
    Line2D([], [], color=light_teal, lw=5.0, solid_capstyle="round", label="Middle 50%"),
    Line2D([], [], linestyle="", marker="o", color=dark_teal, markersize=4.7, label="Median"),
    Line2D([], [], linestyle="", marker="D", color=accent, markersize=4.5, label="Ol Pejeta"),
]
fig.legend(handles=legend, loc="center left", bbox_to_anchor=(.815,.59), frameon=False,
           ncol=1, fontsize=8, handlelength=1.05, handletextpad=.39, labelspacing=1.05,
           borderaxespad=0)
save(fig,"02_typical_density_ranges")

# Inside mean versus the 5 km outside-ring mean.
fig=figure(59)
ax=fig.add_axes([.168,.177,.795,.76])
ax.set_ylabel("Mean inside density (people/km²)",color="#000000",rotation=90,labelpad=5)
ax.text(.027,.948,f"Spearman ρ = {rho:.2f}",transform=ax.transAxes,va="top",
        fontsize=8.4,color=slate)
xmax, ymax = max(280, 1.1 * sample[ring].max()), max(1300, 1.1 * sample[stats[0]].max())
ax.set(xlim=(0,np.log1p(xmax)),ylim=(0,np.log1p(ymax)))
labels(ax,"x",[0,1,5,10,25,50,100,250])
labels(ax,"y",[0,5,25,100,500,1000])
ax.tick_params(axis="both",length=2,pad=2)
ax.set_axisbelow(True)
ax.grid(axis="y",color=grid,lw=.55)
p=np.log1p(xmax)
ax.plot([0, min(p, np.log1p(ymax))], [0, min(p, np.log1p(ymax))],color="#7c8982",lw=.85,ls=(0,(3,3)),zorder=2)
rest=sample[~sample.name.eq(opc["name"])]
ax.scatter(np.log1p(rest[ring]),np.log1p(rest[stats[0]]),s=13.5,
           facecolors=teal,edgecolors="white",lw=.26,alpha=.62,zorder=3)
halo=[pe.withStroke(linewidth=2,foreground="white")]
x,y=np.log1p([opc[ring],opc[stats[0]]])
ax.scatter(x,y,marker="D",s=36,facecolor=accent,edgecolor="white",lw=.72,zorder=6)
ax.annotate("Ol Pejeta",(x,y),xytext=(6,-4),textcoords="offset points",
            ha="left",va="top",weight="bold",fontsize=8.5,color=accent,
            path_effects=halo,zorder=7)
place(fig,70.7,2.4,"Mean density in the 5 km ring outside (people/km²)",
      ha="center",fontsize=8.7)
save(fig,"03_inside_and_outside")
