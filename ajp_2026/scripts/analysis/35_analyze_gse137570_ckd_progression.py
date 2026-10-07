from pathlib import Path
import os
import pandas as pd, numpy as np
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt

base=Path(os.environ.get("AJP_WORKDIR", Path(__file__).resolve().parents[2] / "work"))
base.mkdir(parents=True, exist_ok=True)
out=base/"35_gse137570_ckd_progression_screen"
xlsx=Path(os.environ.get("GSE137570_XLSX", out/"GSE137570_Doyle_et_al.Normalized_read_counts.xlsx"))
sig=["HAVCR1","LCN2","VCAM1","CLU","COL1A1","FN1","ACTA2","TGFB1","TLR4","NLRP3"]

def load_sheet(sheet, symbol_col):
    df=pd.read_excel(xlsx,sheet_name=sheet)
    df=sanitize(df)
    return df

def sanitize(df):
    # strip column strings
    df.columns=[str(c).strip() for c in df.columns]
    return df

def compute_scores(df, symbol_col, sample_cols, meta_rows):
    expr=df.copy()
    expr[symbol_col]=expr[symbol_col].astype(str).str.upper()
    expr=expr.groupby(symbol_col,as_index=True)[sample_cols].mean()
    present=[g for g in sig if g in expr.index]
    missing=[g for g in sig if g not in expr.index]
    mat=expr.loc[present].T
    z=mat.apply(lambda x:(x-x.mean())/x.std(ddof=0),axis=0)
    score=z.mean(axis=1)
    sdf=pd.DataFrame({'sample':score.index,'candidate_score':score.values})
    for k,row in meta_rows.items():
        vals=df[df[symbol_col].astype(str).str.upper()==k.upper()]
        if not vals.empty:
            r=vals.iloc[0]
            sdf[k]=sdf['sample'].map({c:r[c] for c in sample_cols})
    return sdf, present, missing, expr

# Cohort 2 progression.
df2=pd.read_excel(xlsx,sheet_name='COHORT 2_NORMALIZED READ COUNTS')
df2.columns=[str(c).strip() for c in df2.columns]
samples2=[c for c in df2.columns if c!='Symbol']
scores2,present2,missing2,expr2=compute_scores(df2,'Symbol',samples2,{'Age':'Age','Gender':'Gender','CKD_Progression':'CKD_Progression'})
scores2['CKD_Progression']=pd.to_numeric(scores2['CKD_Progression'],errors='coerce')
scores2['progression_status']=scores2['CKD_Progression'].map({0:'non-progressive',1:'progressive'})
scores2.to_csv(out/'GSE137570_cohort2_candidate_scores.csv',index=False,encoding='utf-8-sig')

# Components cohort 2.
sets={
 'candidate_10_gene':sig,
 'injury_HAVCR1_LCN2_VCAM1':['HAVCR1','LCN2','VCAM1'],
 'CLU':['CLU'],
 'fibrosis_ECM_COL1A1_FN1_ACTA2_TGFB1':['COL1A1','FN1','ACTA2','TGFB1'],
 'inflammation_TLR4_NLRP3':['TLR4','NLRP3'],
}
for g in sig: sets['gene_'+g]=[g]
rows=[]; tests=[]
for name,genes in sets.items():
    present=[g for g in genes if g in expr2.index]
    if not present: continue
    mat=expr2.loc[present].T
    z=mat.apply(lambda x:(x-x.mean())/x.std(ddof=0),axis=0)
    score=z.mean(axis=1)
    tmp=pd.DataFrame({'sample':score.index,'score':score.values,'score_name':name,'n_genes':len(present),'genes_present':';'.join(present)})
    tmp=tmp.merge(scores2[['sample','progression_status','CKD_Progression']],on='sample',how='left')
    rows.append(tmp)
    da=tmp.loc[tmp['CKD_Progression']==1,'score'].dropna(); db=tmp.loc[tmp['CKD_Progression']==0,'score'].dropna()
    if len(da)>0 and len(db)>0:
        u,p=stats.mannwhitneyu(da,db,alternative='two-sided')
        tests.append({'score_name':name,'n_genes':len(present),'genes_present':';'.join(present),'n_progressive':len(da),'n_nonprogressive':len(db),'mean_progressive':da.mean(),'mean_nonprogressive':db.mean(),'median_progressive':da.median(),'median_nonprogressive':db.median(),'delta_mean_progressive_minus_nonprogressive':da.mean()-db.mean(),'mannwhitney_p':p})
long2=pd.concat(rows,ignore_index=True)
tests2=pd.DataFrame(tests).sort_values('mannwhitney_p')
if not tests2.empty:
    pvals=tests2['mannwhitney_p'].values; m=len(pvals); order=np.argsort(pvals); adj=np.empty(m); prev=1.0
    for i in range(m-1,-1,-1):
        idx=order[i]; val=pvals[idx]*m/(i+1); prev=min(prev,val); adj[idx]=prev
    tests2['bh_fdr']=adj
long2.to_csv(out/'GSE137570_cohort2_component_gene_scores_long.csv',index=False,encoding='utf-8-sig')
tests2.to_csv(out/'GSE137570_cohort2_component_gene_progression_tests.csv',index=False,encoding='utf-8-sig')

# Cohort 1 GFR/TIF.
df1=pd.read_excel(xlsx,sheet_name='COHORT 1_NORMALIZED READ COUNTS')
df1.columns=[str(c).strip() for c in df1.columns]
samples1=[c for c in df1.columns if c!='SAMPLE ID']
scores1,present1,missing1,expr1=compute_scores(df1,'SAMPLE ID',samples1,{'Age':'Age','Gender':'Gender','GFR':'GFR','TIF':'TIF'})
for col in ['GFR','TIF']:
    scores1[col]=pd.to_numeric(scores1[col],errors='coerce')
scores1.to_csv(out/'GSE137570_cohort1_candidate_scores.csv',index=False,encoding='utf-8-sig')
cor_rows=[]
for var in ['GFR','TIF']:
    sub=scores1[['candidate_score',var]].dropna()
    if len(sub)>2:
        rho,p=stats.spearmanr(sub['candidate_score'],sub[var])
        cor_rows.append({'cohort':'cohort1','variable':var,'n':len(sub),'spearman_rho':rho,'p_value':p})
cor=pd.DataFrame(cor_rows)
cor.to_csv(out/'GSE137570_cohort1_candidate_correlations.csv',index=False,encoding='utf-8-sig')

# Gene presence.
pd.DataFrame({'gene':sig,'present_cohort1':[g in present1 for g in sig],'present_cohort2':[g in present2 for g in sig]}).to_csv(out/'GSE137570_signature_gene_presence.csv',index=False,encoding='utf-8-sig')

print('cohort2 present',present2,'missing',missing2)
print(tests2.head(15).to_string(index=False))
print('cohort1 cor')
print(cor.to_string(index=False))

# Plot.
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],'svg.fonttype':'none','pdf.fonttype':42,'font.size':8,'axes.spines.right':False,'axes.spines.top':False})
fig,axes=plt.subplots(1,3,figsize=(8.6,3.2),gridspec_kw={'width_ratios':[0.9,1,1]})
# progression box/scatter
ax=axes[0]
order=['non-progressive','progressive']; colors={'non-progressive':'#6C757D','progressive':'#D55E00'}
for i,g in enumerate(order):
    y=scores2.loc[scores2['progression_status']==g,'candidate_score'].dropna()
    rng=np.random.default_rng(200+i); x=np.full(len(y),i)+rng.uniform(-0.08,0.08,len(y))
    ax.scatter(x,y,color=colors[g],s=32,edgecolor='white',linewidth=0.5,zorder=3)
    ax.plot([i-0.18,i+0.18],[y.mean(),y.mean()],color='black',lw=1.4)
ax.set_xticks([0,1]); ax.set_xticklabels(['Non-prog.','Prog.'])
ax.set_ylabel('10-gene candidate score')
main=tests2[tests2['score_name']=='candidate_10_gene']
if not main.empty:
    r=main.iloc[0]
    ax.text(0.03,0.96,f"delta={r['delta_mean_progressive_minus_nonprogressive']:.2f}\nP={r['mannwhitney_p']:.3g}",transform=ax.transAxes,ha='left',va='top')
ax.set_title('Cohort 2 progression',loc='left',fontweight='bold')
ax.grid(axis='y',color='#E9ECEF')
# GFR
ax=axes[1]
ax.scatter(scores1['GFR'],scores1['candidate_score'],color='#0072B2',s=32,edgecolor='white',linewidth=0.5)
sub=scores1[['GFR','candidate_score']].dropna()
if len(sub)>2:
    coef=np.polyfit(sub['GFR'],sub['candidate_score'],1); xs=np.linspace(sub['GFR'].min(),sub['GFR'].max(),100); ax.plot(xs,coef[0]*xs+coef[1],color='#333333',lw=1)
    hit=cor[cor['variable']=='GFR'].iloc[0]
    ax.text(0.04,0.96,f"rho={hit['spearman_rho']:.2f}\nP={hit['p_value']:.3g}",transform=ax.transAxes,ha='left',va='top')
ax.set_xlabel('GFR'); ax.set_ylabel('10-gene candidate score')
ax.set_title('Cohort 1 kidney function',loc='left',fontweight='bold')
ax.grid(color='#E9ECEF')
# TIF
ax=axes[2]
ax.scatter(scores1['TIF'],scores1['candidate_score'],color='#009E73',s=32,edgecolor='white',linewidth=0.5)
sub=scores1[['TIF','candidate_score']].dropna()
if len(sub)>2:
    coef=np.polyfit(sub['TIF'],sub['candidate_score'],1); xs=np.linspace(sub['TIF'].min(),sub['TIF'].max(),100); ax.plot(xs,coef[0]*xs+coef[1],color='#333333',lw=1)
    hit=cor[cor['variable']=='TIF'].iloc[0]
    ax.text(0.04,0.96,f"rho={hit['spearman_rho']:.2f}\nP={hit['p_value']:.3g}",transform=ax.transAxes,ha='left',va='top')
ax.set_xlabel('Tubulointerstitial fibrosis (%)'); ax.set_ylabel('10-gene candidate score')
ax.set_title('Cohort 1 histology',loc='left',fontweight='bold')
ax.grid(color='#E9ECEF')
fig.suptitle('GSE137570 patient-level CKD progression and clinicopathologic screen',x=0.01,ha='left',fontweight='bold',fontsize=10)
fig.tight_layout(rect=[0,0,1,0.92])
for ext in ['png','pdf','svg','tiff']:
    if ext in ['png','tiff']:
        fig.savefig(base/f'35_gse137570_candidate_screen.{ext}',dpi=600,bbox_inches='tight')
    else:
        fig.savefig(base/f'35_gse137570_candidate_screen.{ext}',bbox_inches='tight')
plt.close(fig)

# Component plot.
plot=tests2.copy().head(10)
if 'candidate_10_gene' not in set(plot['score_name']):
    plot=pd.concat([tests2[tests2['score_name']=='candidate_10_gene'],plot.head(9)],ignore_index=True)
plot=plot.iloc[::-1]
fig,ax=plt.subplots(figsize=(6.8,3.6))
cols=['#D55E00' if d>0 else '#6C757D' for d in plot['delta_mean_progressive_minus_nonprogressive']]
ax.barh(plot['score_name'],plot['delta_mean_progressive_minus_nonprogressive'],color=cols)
for y,(_,r) in enumerate(plot.iterrows()):
    d=r['delta_mean_progressive_minus_nonprogressive']
    ax.text(d+(0.03 if d>=0 else -0.03),y,f"P={r['mannwhitney_p']:.2g}",va='center',ha='left' if d>=0 else 'right',fontsize=7)
ax.axvline(0,color='#333333',lw=0.8)
ax.set_xlabel('Mean score delta: progressive minus non-progressive')
ax.set_title('GSE137570 Cohort 2 component/gene progression screen',loc='left',fontweight='bold')
fig.tight_layout()
for ext in ['png','pdf','svg','tiff']:
    if ext in ['png','tiff']:
        fig.savefig(base/f'35_gse137570_component_gene_screen.{ext}',dpi=600,bbox_inches='tight')
    else:
        fig.savefig(base/f'35_gse137570_component_gene_screen.{ext}',bbox_inches='tight')
plt.close(fig)

# Report.
lines=[]
lines.append('# Step 35. GSE137570 CKD progression and clinicopathologic screen')
lines.append('')
lines.append('## Dataset')
lines.append('')
lines.append('GSE137570 contains human kidney biopsy transcriptomic data from two cohorts. Cohort 2 includes binary progressive CKD status; Cohort 1 includes GFR and tubulointerstitial fibrosis (TIF).')
lines.append('')
lines.append('## Gene recovery')
lines.append('')
lines.append('Cohort 2 detected candidate genes: '+', '.join(present2))
lines.append('Cohort 2 missing candidate genes: '+(', '.join(missing2) if missing2 else 'none'))
lines.append('Cohort 1 detected candidate genes: '+', '.join(present1))
lines.append('Cohort 1 missing candidate genes: '+(', '.join(missing1) if missing1 else 'none'))
lines.append('')
lines.append('## Primary 10-gene candidate')
lines.append('')
if not main.empty:
    r=main.iloc[0]
    lines.append(f"- Cohort 2 progressive vs non-progressive CKD: n={int(r['n_progressive'])} vs {int(r['n_nonprogressive'])}; delta={r['delta_mean_progressive_minus_nonprogressive']:.3g}; Mann-Whitney P={r['mannwhitney_p']:.3g}; FDR={r['bh_fdr']:.3g}.")
for _,r in cor.iterrows():
    lines.append(f"- Cohort 1 candidate score vs {r['variable']}: n={int(r['n'])}; Spearman rho={r['spearman_rho']:.3g}; P={r['p_value']:.3g}.")
lines.append('')
lines.append('## Component/gene screen, Cohort 2')
lines.append('')
for _,r in tests2.head(10).iterrows():
    lines.append(f"- {r['score_name']}: delta={r['delta_mean_progressive_minus_nonprogressive']:.3g}; P={r['mannwhitney_p']:.3g}; FDR={r['bh_fdr']:.3g}; genes={r['genes_present']}")
lines.append('')
lines.append('## Interpretation')
lines.append('')
if not main.empty and float(main.iloc[0]['delta_mean_progressive_minus_nonprogressive'])>0 and float(main.iloc[0]['mannwhitney_p'])<0.05:
    lines.append('The 10-gene candidate is higher in progressive CKD samples in Cohort 2, providing patient-level outcome support in an external CKD biopsy dataset. Cohort 1 also provides cross-sectional clinicopathologic context through GFR and TIF correlations. This can be added as a supplementary validation layer, with care not to claim individualized prognostic validation.')
else:
    lines.append('The primary 10-gene candidate does not provide a strong progression-group validation in Cohort 2. Use only as a screened dataset unless components show a biologically useful pattern.')
lines.append('')
lines.append('## Output files')
lines.append('')
for f in ['35_gse137570_ckd_progression_screen/GSE137570_cohort2_candidate_scores.csv','35_gse137570_ckd_progression_screen/GSE137570_cohort2_component_gene_progression_tests.csv','35_gse137570_ckd_progression_screen/GSE137570_cohort1_candidate_correlations.csv','35_gse137570_candidate_screen.png','35_gse137570_component_gene_screen.png']:
    lines.append(f'- `{f}`')
(base/'35_gse137570_ckd_progression_screen.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
print(base/'35_gse137570_ckd_progression_screen.md')
# AJP 2026 deposition copy. Configure local data/output paths as described in
# ajp_2026/docs/run_order.md; raw third-party matrices are not redistributed.
