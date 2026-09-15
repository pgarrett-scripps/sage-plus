"""Generate comparative release figures from the audited report snapshots."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from _assets import record
from _scientific import PAPER, INPUTS, ENGINE, load
from _report import report, REPORT_INPUTS, workloads

COLORS = {'upstream':'#32658a','plus':'#d87532'}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,'figure.facecolor':'white','savefig.facecolor':'white'})


def finish(fig,name,desc):
    target=PAPER/'figures'/f'report-{name}.png'
    fig.savefig(target,dpi=320,bbox_inches='tight',metadata={'Software':None})
    plt.close(fig)
    sections={'identification':'public','disagreement':'public','scaling':'scaling','ptm':'ptm','lfq':'lfq','control':'lfq'}
    inputs=INPUTS if name=='workloads' else [f'analysis/data/report-extension/{sections[name]}.json']
    record(f'fig.report-{name}',str(target.relative_to(PAPER)),kind='figure',inputs=inputs,desc=desc)


def label(ax,letter,title):
    ax.set_title(f'{letter}  {title}',loc='left',fontweight='bold',pad=12)
    ax.grid(axis='y',color='#e6eaee',linewidth=.7,zorder=0)
    ax.set_axisbelow(True)


def identification():
    rows=report('public')
    names=['HEK 1','HEK 2','A Alpha','B Alpha','A Beta','B Beta']
    fig,axes=plt.subplots(2,2,figsize=(8,6.6),layout='constrained')
    for ax,metric,title,letter in zip(axes[0],('target_psms','target_peptidoforms'),('Accepted spectrum matches','Accepted peptidoforms'),('A','B')):
        for e,offset in (('upstream',-.18),('plus',.18)):
            ax.bar(np.arange(6)+offset,[r['engines'][e]['0.01'][metric]/1000 for r in rows],.34,color=COLORS[e],label=ENGINE[e])
        ax.set_xticks(range(6),names,rotation=35,ha='right')
        ax.set_ylabel('Accepted targets (thousands)')
        label(ax,letter,title)
    axes[0,0].legend(fontsize=9)
    for ax,study,letter in zip(axes[1],('PXD001468','PXD028735'),('C','D')):
        selected=[r for r in rows if study in r['pair']]
        qs=(.001,.005,.01,.02,.05)
        for metric,marker,text in (('target_psms','o','PSMs'),('target_peptidoforms','s','Peptidoforms')):
            curves=np.array([[100*(r['engines']['plus'][str(q)][metric]/r['engines']['upstream'][str(q)][metric]-1) for q in qs] for r in selected])
            color='#804f87' if metric=='target_peptidoforms' else '#287b79'
            ax.plot(np.array(qs)*100,curves.mean(axis=0),marker=marker,color=color,label=text)
            ax.fill_between(np.array(qs)*100,curves.min(axis=0),curves.max(axis=0),color=color,alpha=.12)
        ax.axhline(0,color='#555',lw=.8)
        ax.set_xscale('log')
        ax.set_xticks([.1,.5,1,2,5],['0.1','0.5','1','2','5'])
        ax.set_xlabel('Reported q threshold (%)')
        ax.set_ylabel('Sage Plus change from Sage (%)')
        label(ax,letter,'HEK threshold response' if study=='PXD001468' else 'Mixture threshold response')
    axes[1,0].legend(fontsize=9)
    finish(fig,'identification','Public PSM and peptidoform yields and nominal-threshold response')


def disagreement():
    rows=report('public')
    names=['HEK 1','HEK 2','A Alpha','B Alpha','A Beta','B Beta']
    fig,axes=plt.subplots(1,2,figsize=(8,3.8),layout='constrained')
    cats=[('same_assignment_above_threshold','Same assignment, above q','#77a9b5'),('different_assignment','Different rank-one assignment','#d3a34b'),('no_matching_rank_one_spectrum_charge','No matched spectrum and charge','#9b86ad')]
    for ax,e,letter in zip(axes,('upstream','plus'),('A','B')):
        left=np.zeros(6)
        for key,title,color in cats:
            vals=np.array([r['disagreement'][e].get(key,0) for r in rows])
            ax.barh(range(6),vals,left=left,color=color,label=title)
            left+=vals
        ax.set_yticks(range(6),names)
        ax.invert_yaxis()
        ax.set_xlabel('Engine-only accepted PSMs')
        ax.set_xlim(0,1150)
        label(ax,letter,f'{ENGINE[e]} only')
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=1,fontsize=9)
    finish(fig,'disagreement','Decomposition of accepted PSM disagreement using raw rank-one output')


def tradeoff():
    rows=workloads()
    fig,axes=plt.subplots(1,2,figsize=(8,3.8),layout='constrained')
    for ax,key,title,letter in zip(axes,('seconds','rss'),('Wall time ratio','Peak memory ratio'),('A','B')):
        values=[r['engines']['plus'][key]/r['engines']['upstream'][key] for r in rows]
        for i,(row,value) in enumerate(zip(rows,values)):
            color=['#287b79','#804f87','#737b83'][i//2]
            ax.plot([1,value],[i,i],color=color,lw=2)
            ax.scatter(value,i,color=color,s=45,zorder=3)
            ax.annotate(f'{value:.2f}',(value,i),xytext=(5,5),textcoords='offset points',fontsize=9)
        ax.set_yticks(range(len(rows)),[r['label'] for r in rows])
        ax.invert_yaxis()
        ax.axvline(1,color='#555',lw=.8,ls='--')
        ax.set_xlim(.6,1.35)
        ax.set_xlabel('Sage Plus / Sage')
        label(ax,letter,title)
    fig.supxlabel('Green: repeated public input    Purple: file and seed summaries    Gray: local context',fontsize=8)
    finish(fig,'workloads','Within-workload performance ratios with distinct replication designs')


def scaling():
    rows=report('scaling')
    fig,axes=plt.subplots(1,3,figsize=(8.4,3.2),layout='constrained')
    for e in ('upstream','plus'):
        baseline=np.median([r['seconds'] for r in rows if r['engine']==e and r['threads']==1 and not r['warmup']])
        for ax,key in zip(axes,('seconds','rss','speedup')):
            groups=[[r for r in rows if r['engine']==e and r['threads']==t and not r['warmup']] for t in (1,2,4,8)]
            values=[np.median([r['seconds' if key=='speedup' else key] for r in g]) for g in groups]
            if key=='speedup':
                values=baseline/np.array(values)
            elif key=='rss':
                values=np.array(values)/1024
            ax.plot((1,2,4,8),values,'o-',color=COLORS[e],label=ENGINE[e])
            if key!='speedup':
                for t,g in zip((1,2,4,8),groups):
                    ax.scatter([t]*len(g),[r[key]/(1024 if key=='rss' else 1) for r in g],s=10,color=COLORS[e],alpha=.5)
    for ax,letter,title,y in zip(axes,'ABC',('Wall time','Peak resident memory','Parallel speedup'),('Seconds','GiB','One-worker time / time')):
        label(ax,letter,title)
        ax.set_xticks([1,2,4,8])
        ax.set_xlabel('Workers')
        ax.set_ylabel(y)
    axes[0].legend(fontsize=9)
    finish(fig,'scaling','Matched release thread scaling with warmups excluded')


def ptm():
    rows=report('ptm')
    fig,axes=plt.subplots(1,2,figsize=(8,3.5),layout='constrained')
    for ax,lib,letter in zip(axes,(1,2),'AB'):
        for e,offset in (('upstream',-.18),('plus',.18)):
            r=next(r for r in rows if r['engine']==e and r['library']==lib)
            vals=[r['spectrum_accepted'],r['joint_accepted']]
            ax.bar(np.arange(2)+offset,vals,.34,color=COLORS[e],label=ENGINE[e])
            for x,y in zip(np.arange(2)+offset,vals):
                ax.text(x,y+60,f'{y:,}',ha='center',fontsize=10,color=COLORS[e])
        ax.set_xticks([0,1],['Spectrum q ≤ 1%','Spectrum and peptide\nq ≤ 1%'])
        ax.set_ylabel('Accepted target PSMs')
        ax.set_ylim(0,max(r['spectrum_accepted'] for r in rows)*1.23)
        label(ax,letter,f'Synthetic HCD {lib}')
    axes[0].legend(fontsize=9,loc='upper right')
    finish(fig,'ptm','Synthetic PTM acceptance under spectrum-only and joint peptide filters')


def quantification():
    q=report('lfq')
    fig,axes=plt.subplots(2,3,figsize=(8.5,6),layout='constrained')
    species=('human','yeast','ecoli')
    expected=(0,-1,2)
    for ax,sp,exp,letter in zip(axes[0],species,expected,'ABC'):
        for row in q['engines']:
            values=list(row['ratios'][sp].values())
            ax.hist(values,bins=np.linspace(-5,6,89),density=False,weights=np.ones(len(values))/len(values),histtype='step',linewidth=1.5,color=COLORS[row['engine']],label=ENGINE[row['engine']])
        ax.axvline(exp,color='#555',ls='--',lw=.8)
        ax.set_xlim(exp-2,exp+2)
        ax.set_xlabel('Observed log₂(B/A)')
        ax.set_ylabel('Fraction per bin')
        label(ax,letter,{'human':'Human','yeast':'Yeast','ecoli':'E. coli'}[sp])
    axes[0,0].legend(fontsize=9)
    for ax,key,title,letter in zip(axes[1],('median_absolute_log2_error','median_preparation_cv','missing_fraction_observed_union'),('Absolute ratio error','Preparation CV','Feature missingness'),'DEF'):
        for row,offset in zip(q['engines'],(-.18,.18)):
            vals=[row['species'][sp][key]*(100 if key!='median_absolute_log2_error' else 1) for sp in species]
            ax.bar(np.arange(3)+offset,vals,.34,color=COLORS[row['engine']])
        ax.set_xticks(range(3),['Human','Yeast','E. coli'])
        ax.set_ylabel('Median |log₂ ratio error|' if key=='median_absolute_log2_error' else ('Median CV (%)' if key=='median_preparation_cv' else 'Missing file slots (%)'))
        label(ax,letter,title)
    finish(fig,'lfq','Matched LFQ distributions, ratio error, preparation variability and missingness')


def control():
    rows=report('lfq')['engines']
    fig,axes=plt.subplots(1,3,figsize=(8.5,3.6),layout='constrained')
    for e,offset in (('upstream',-.18),('plus',.18)):
        row=next(r for r in rows if r['engine']==e)
        axes[0].bar(np.arange(2)+offset,[row['thresholds'][q]['target_precursors']/1000 for q in ('0.01','0.05')],.34,color=COLORS[e],label=ENGINE[e])
    axes[0].set_xticks([0,1],['1%','5%'])
    axes[0].set_xlabel('LFQ q threshold')
    axes[0].set_ylabel('Target precursors (thousands)')
    label(axes[0],'A','Quantification yield')
    for i,row in enumerate(rows):
        c=row['control']
        for ax,total,foreign in ((axes[1],c['quantified'],c['foreign']),(axes[2],c['without_strict_ms2'],c['foreign_without_strict_ms2'])):
            ax.bar(i,total-foreign,color='#94b8bf',label='Human' if i==0 else None)
            ax.bar(i,foreign,bottom=total-foreign,color='#a86674',label='Foreign' if i==0 else None)
            ax.text(i,total+max(200,total*.03),f'{foreign/total:.1%}',ha='center',fontsize=9)
    for ax,letter,title in zip(axes[1:],'BC',('Human-only control','No direct MS2 support')):
        ax.set_xticks([0,1],['Sage','Sage Plus'])
        ax.set_ylabel('Positive quantified file rows')
        label(ax,letter,title)
        ax.margins(y=.2)
    axes[0].legend(fontsize=8)
    axes[1].set_ylim(0,max(r['control']['quantified'] for r in rows)*1.38)
    axes[1].legend(fontsize=8,loc='upper right')
    finish(fig,'control','LFQ threshold yield and species-absent control diagnostics for both releases')


if __name__ == '__main__':
    for function in (identification,disagreement,tradeoff,scaling,ptm,quantification,control):
        function()
