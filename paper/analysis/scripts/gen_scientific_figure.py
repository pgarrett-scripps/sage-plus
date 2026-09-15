"""Generate public timing and expanded entrapment release comparisons."""
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from _assets import record
from _scientific import PAPER, INPUTS, ENGINE, load

source=PAPER.parent/'benchmarks/scientific-results/20260914/figures/public-timing-pilot.png'
target=PAPER/'figures/scientific-public-timing-pilot.png'
shutil.copyfile(source,target)
record('fig.scientific-public-timing-pilot',str(target.relative_to(PAPER)),kind='figure',inputs=['../benchmarks/scientific-results/20260914/figures/public-timing-pilot.png','../benchmarks/scientific-results/20260914/figures/figures.json'],desc='Frozen public release timing comparison')

pilot,_=load()
colors={'upstream':'#32658a','plus':'#d87532'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False})
fig,axes=plt.subplots(2,2,figsize=(8,6.5),layout='constrained')
for col,study in enumerate(('human','hye')):
    ax=axes[0,col]
    for e in ('upstream','plus'):
        cells=[r for r in pilot['calibration'] if r['suite'].startswith(f'entrapment-{study}-') and r['job'].endswith('-'+e)]
        x=np.array([p['nominal_q'] for p in cells[0]['thresholds']])*100
        values=np.array([[p['paired_fdp_tie_max']*100 for p in r['thresholds']] for r in cells])
        for row in values:
            ax.plot(x,row,color=colors[e],alpha=.18,lw=.65)
        ax.plot(x,values.mean(axis=0),'o-',color=colors[e],label=ENGINE[e])
    ax.plot([.1,5],[.1,5],ls='--',color='#777',lw=.8,label='FDP = nominal q')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xticks([.1,.5,1,2,5],['0.1','0.5','1','2','5'])
    ax.set_yticks([.1,.5,1,2,5],['0.1','0.5','1','2','5'])
    ax.set_xlabel('Nominal peptide q (%)')
    ax.set_ylabel('Conservative paired FDP (%)')
    ax.set_title(f"{'AB'[col]}  {'HEK' if study=='human' else 'Mixture'} calibration",loc='left',fontweight='bold')
    ax=axes[1,col]
    deltas=[]
    for file in (0,1):
        for seed in (20260914,20260915,20260916):
            values={e:next(p['paired_fdp_tie_max'] for r in pilot['calibration'] if r['suite']==f'entrapment-{study}-{seed}' and r['job']==f'file-{file}-{e}' for p in r['thresholds'] if p['nominal_q']==.01) for e in ('upstream','plus')}
            deltas.append(100*(values['plus']-values['upstream']))
    ax.scatter(deltas,range(6),color='#287b79',s=35)
    summary=next(r for r in pilot['calibration_uncertainty'] if r['study']==study)
    mean=100*(summary['means']['plus']-summary['means']['upstream'])
    lo,hi=np.array(summary['percentile_intervals']['paired_difference'])*100
    ax.errorbar(mean,6,xerr=[[mean-lo],[hi-mean]],fmt='D',color='#804f87',capsize=3)
    ax.set_yticks(range(7),[f'File {f+1}, seed {i+1}' for f in (0,1) for i in range(3)]+['Mean and interval'])
    ax.invert_yaxis()
    ax.axvline(0,color='#777',ls='--',lw=.8)
    ax.set_xlim(-.14,.08)
    ax.set_xticks([-.10,-.05,0,.05],['−0.10','−0.05','0.00','0.05'])
    ax.set_xlabel('Sage Plus minus Sage FDP (pp)')
    ax.set_title(f"{'CD'[col]}  Paired differences at 1%",loc='left',fontweight='bold')
axes[0,0].legend(fontsize=8)
target=PAPER/'figures/scientific-entrapment-pilot.png'
fig.savefig(target,dpi=320,bbox_inches='tight',metadata={'Software':None})
plt.close(fig)
record('fig.scientific-entrapment-pilot',str(target.relative_to(PAPER)),kind='figure',inputs=INPUTS,desc='Threshold calibration and file-by-seed paired differences for the frozen release pair')
