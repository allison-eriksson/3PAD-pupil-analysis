"""Read-only audit of saved raw pupil samples, interpolation and retained trials.
No exclusions are applied and no model inputs/results are changed.
"""
import argparse
import csv,json,time,hashlib
from pathlib import Path
import h5py,numpy as np
from scipy.io import loadmat
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--est-root", type=Path, required=True)
parser.add_argument("--retention-csv", type=Path, required=True)
parser.add_argument("--output-dir", type=Path, required=True)
args = parser.parse_args()
BASE, OUT = args.est_root, args.output_dir
OUT.mkdir(parents=True, exist_ok=False)
with args.retention_csv.open() as stream:
    COHORT = list(csv.DictReader(stream))
def write(name,rows):
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in rows for k in row)));w.writeheader();w.writerows(rows)
def string(f,ref): return ''.join(chr(int(x)) for x in f[ref][...].ravel())
def getrun(f,pid,visit):
 sub=f[f[f'AE_EYE_STRUCT_EST_R{visit}/sub'][pid-1,0]]
 run=sub['run']
 if isinstance(run,h5py.Dataset):
  refs=run[...].ravel();run=f[next(x for x in refs if x)]
 return run
runs=[];trials=[]
for visit in [1,2]:
 bev=loadmat(BASE/f'Results_EST/AE_BEV_MAT_EST_R{visit}.mat')[f'AE_BEV_MAT_EST_R{visit}']
 saved=loadmat(BASE/f'Results_EST/EYE_ANAL/RUN{visit}/allSubsAllBev_LowPass_R{visit}.mat')['allSubsAllBev']
 with h5py.File(BASE/f'Results_EST/AE_EYE_STRUCT_EST_R{visit}.mat') as processed:
  for rec in [r for r in COHORT if int(r['visit'])==visit and r['record_available']=='True']:
   pid=int(rec['pid']);row={'pid':pid,'visit':visit}
   try:
    folder=BASE/f'{pid:03d}'/f'r{visit}'
    exps=sorted(p for p in folder.glob('EST*') if p.is_dir())
    files=[e.with_suffix('.mat') for exp in exps for e in sorted(exp.glob('*.edf')) if e.with_suffix('.mat').exists()]
    row['converted_file_candidates']=len(files)
    if len(files)!=1:raise ValueError(f'Expected one converted file; found {len(files)}')
    path=files[0];row['source']=str(path)
    with h5py.File(path) as f:
     mc=f['#subsystem#/MCOS'];s=f[mc[0,7]];ev=f[mc[0,8]]
     assert 'pupilSize' in s
     t=s['time'][...].ravel();v=s['pupilSize'][0,:]
     m=ev['Messages'];mt=m['time'][...].ravel();msgs=[string(f,x) for x in m['info'][...].ravel()]
     starts=[mt[i] for i,x in enumerate(msgs) if x.startswith('EXPERIMENT START')]
     assert len(starts)==1,starts
     start=starts[0]
    step=int(np.median(np.diff(t)));row['sample_step_ms']=step
    width=100 if visit==1 else 200
    row['gaussian_support_samples']=width
    raw=v<=100;mask=np.convolve(raw.astype(np.int32),np.ones(width,dtype=np.int32),'full')[width//2:width//2+len(raw)]>0
    ons=np.where(np.diff(mask.astype(int))==1)[0]+1;ends=np.where(np.diff(mask.astype(int))==-1)[0]+1
    at_start=bool(mask[0]);at_end=bool(mask[-1])
    if at_start:ons=np.r_[0,ons]
    if at_end:ends=np.r_[ends,len(v)-1]
    w=v.copy();applied=np.zeros(len(v),bool)
    for a,z in zip(ons,ends):w[a:z+1]=np.linspace(v[a],v[z],z-a+1);applied[a:z+1]=True
    if at_start:w[ons[0]:ends[0]+1]=v[ends[0]]
    if at_end:w[ons[-1]:ends[-1]+1]=v[ons[-1]]
    stored=getrun(processed,pid,visit)['currP'][...].ravel()
    rebuilt=np.repeat(w[t>=start],step)
    row['stored_samples']=len(stored);row['reconstructed_samples']=len(rebuilt)
    if len(stored)!=len(rebuilt):raise ValueError('Stored signal length differs from raw reconstruction')
    row['interpolation_max_abs_difference']=float(np.nanmax(np.abs(rebuilt-stored)))
    row['interpolation_matches_saved']=bool(np.allclose(rebuilt,stored,rtol=0,atol=1e-8,equal_nan=True))
    b=bev[bev[:,16]==pid];sv=saved[saved[:,16]==pid]
    assert len(b)==len(sv)==int(rec['recorded'])==60
    assert set(b[:,0])==set(sv[:,0])
    row['recorded_trials']=len(b);row['retained_trials']=int(rec['retained'])
    row['raw_nonfinite_samples']=int(np.sum(~np.isfinite(v)))
    # Diagnostic windows cover the fixed pupil windows; they are not exclusion rules.
    rel=t-start
    for trial in b:
     onset=trial[7]-trial[14];rt=trial[4];included=bool(trial[1]!=99)
     q={'pid':pid,'visit':visit,'trial':int(trial[0]),'response_status':int(trial[1]),'retained':included,'stimulus_onset_ms':float(onset),'response_time_ms':float(rt)}
     windows={'fixed_0_to_4138':(0,4138),'early_mean_500_to_850':(500,850),'early_peak_500_to_2000':(500,2000),'late_2128_to_4138':(2128,4138)}
     if included and np.isfinite(rt):windows.update({'response_minus250_to_plus250':(rt-250,rt+250),'onset_to_response':(0,rt)})
     for name,(a,z) in windows.items():
      ix=(rel>=onset+a)&(rel<=onset+z);n=int(ix.sum());q[name+'_samples']=n
      q[name+'_raw_low_samples']=int(raw[ix].sum());q[name+'_nonfinite_samples']=int((~np.isfinite(v[ix])).sum());q[name+'_interpolated_samples']=int(applied[ix].sum())
      q[name+'_raw_low_fraction']=float(raw[ix].mean()) if n else None
      q[name+'_interpolated_fraction']=float(applied[ix].mean()) if n else None
     trials.append(q)
    row['status']='ok'
   except Exception as exc:row['status']='error';row['error']=str(exc)
   runs.append(row)
   write('raw_signal_recovery_by_assessment.csv',runs);write('raw_signal_quality_by_trial.csv',trials)
   print(f"R{visit} pid {pid}: {row['status']} matches={row.get('interpolation_matches_saved')} completed={len(runs)}/171",flush=True)
write('raw_signal_recovery_by_assessment.csv',runs);write('raw_signal_quality_by_trial.csv',trials)
summary={}
for visit in [1,2]:
 r=[x for x in runs if x['visit']==visit];q=[x for x in trials if x['visit']==visit and x['retained']]
 s={'assessments':len(r),'errors':[x for x in r if x['status']!='ok'],'exact_interpolation_matches':sum(x.get('interpolation_matches_saved',False) for x in r),'retained_trials_audited':len(q)}
 for stem in ['fixed_0_to_4138','early_mean_500_to_850','early_peak_500_to_2000','late_2128_to_4138','response_minus250_to_plus250','onset_to_response']:
  a=np.array([x[stem+'_raw_low_fraction'] for x in q if x.get(stem+'_raw_low_fraction') is not None]);im=np.array([x[stem+'_interpolated_fraction'] for x in q if x.get(stem+'_interpolated_fraction') is not None])
  if len(a):s[stem]={'raw_low_percent_quantiles':(100*np.quantile(a,[0,.25,.5,.75,1])).tolist(),'interpolated_percent_quantiles':(100*np.quantile(im,[0,.25,.5,.75,1])).tolist(),'raw_low_at_least_20pct_trials':int(sum(a>=.2)),'raw_low_at_least_50pct_trials':int(sum(a>=.5)),'interpolated_at_least_50pct_trials':int(sum(im>=.5))}
 summary[str(visit)]=s
(OUT/'quality_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2),flush=True)
