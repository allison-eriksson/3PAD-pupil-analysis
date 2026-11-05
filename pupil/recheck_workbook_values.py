from pathlib import Path
import argparse
import csv,json,gc
import numpy as np
from scipy.io import loadmat
from openpyxl import load_workbook
parser = argparse.ArgumentParser(description="Calculate the eight corrected workbook columns.")
parser.add_argument("--original-workbook", type=Path, required=True)
parser.add_argument("--est-results", type=Path, required=True)
parser.add_argument("--timing-output", type=Path, required=True)
parser.add_argument("--audit-dir", type=Path, required=True)
args = parser.parse_args()
source, est, audit = args.original_workbook, args.est_results, args.audit_dir
audit.mkdir(parents=True, exist_ok=False)
w=load_workbook(source,data_only=False);s=w.worksheets[0];heads=[c.value for c in s[1]]
raw={int(row[0]):dict(zip(heads,row)) for row in s.iter_rows(min_row=2,values_only=True)}
model=list(csv.DictReader((args.timing_output/'analysis_unimputed_corrected.csv').open()))
model={int(x['pid']):x for x in model};updates=[];examples=[];stats=[]
for visit in (1,2):
 print('Rechecking saved pupil traces, visit',visit,flush=True)
 t=loadmat(est/f'EYE_ANAL/RUN{visit}/allSubsAllBev_LowPass_R{visit}.mat')['allSubsAllBev']
 p=loadmat(est/f'EYE_ANAL/RUN{visit}/allSubsAllPup_LowPass_R{visit}.mat')['allSubsAllPup']
 suffix='' if visit==1 else '2';fixed=t.copy();n=0;nonii=0;old_err=0.;rt_err=0.;example_count=0
 for i in np.flatnonzero(np.isin(t[:,16],list(raw)) & (t[:,1]!=99)):
  rt=int(t[i,4]);condition=int(t[i,3]);z=10000;trace=p[i]
  true_rod=float(np.sum(np.diff(trace[z:z+rt+1])))
  true_par=float(np.mean(trace[z+rt-250:z+rt+251]))
  if int(t[i,15])==4:
   bad_rod=float(np.sum(np.diff(trace[z:z+condition+1])));bad_par=float(np.mean(trace[z+condition-250:z+condition+251]))
   assert np.isclose(bad_rod,t[i,20],rtol=0,atol=1e-10)
   assert np.isclose(bad_par,t[i,19],rtol=0,atol=1e-10)
   old_err=max(old_err,abs(bad_rod-t[i,20]));rt_err=max(rt_err,abs(true_rod-t[i,20]));n+=1
   fixed[i,19]=true_par;fixed[i,20]=true_rod
   if example_count<4 and abs(true_rod-t[i,20])>.01:
    examples.append(dict(visit=visit,participant=int(t[i,16]),saved_trial_row=int(i+1),condition_code=condition,actual_response_ms=rt,stored_ROD=float(t[i,20]),ROD_using_condition_code=bad_rod,ROD_using_response_time=true_rod,stored_PAR=float(t[i,19]),PAR_using_condition_code=bad_par,PAR_using_response_time=true_par));example_count+=1
  else:
   assert np.isclose(true_rod,t[i,20],rtol=0,atol=1e-10)
   assert np.isclose(true_par,t[i,19],rtol=0,atol=1e-10)
   nonii+=1
 for pid,row in raw.items():
  old=t[(t[:,16]==pid)&(t[:,1]!=99)];new=fixed[(fixed[:,16]==pid)&(fixed[:,1]!=99)]
  for stem,col in [('PARDiff',19),('RodDiff',20),('PARMean',19),('RodMean',20)]:
   v=stem+suffix
   if row[v] is None:continue
   assert len(old), (pid,v,'no trial records')
   before=float(old[:,col].mean()) if stem.endswith('Mean') else float(old[old[:,2]==2,col].mean()-old[old[:,2]==1,col].mean())
   after=float(new[:,col].mean()) if stem.endswith('Mean') else float(new[new[:,2]==2,col].mean()-new[new[:,2]==1,col].mean())
   assert abs(before-float(row[v]))<1e-8,(pid,v,before,row[v])
   if pid in model and v in model[pid]:assert abs(after-float(model[pid][v]))<1e-8,(pid,v,'model mismatch')
   address=s.cell(row=list(raw).index(pid)+2,column=heads.index(v)+1).coordinate
   updates.append(dict(pid=pid,visit=visit,variable=v,cell=address,original=float(row[v]),corrected=after,analytic_cohort=pid in model,model_input=v in model[next(iter(model))]))
 stats.append(dict(visit=visit,II_trials_verified=n,other_trials_correct_RT_verified=nonii,max_error_matching_condition_code_ROD=old_err,max_difference_from_actual_RT_ROD=rt_err))
 del p,t,fixed;gc.collect()
(audit/'workbook_updates.json').write_text(json.dumps(updates,indent=2))
for name,rows in [('timing_trial_examples.csv',examples),('merged_workbook_correction_log.csv',updates)]:
 with (base/'audit'/name).open('w',newline='') as stream:
  writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
report=dict(source=str(source),rows=len(raw),columns=len(heads),updated_observed_values=len(updates),analytic_cohort_values=sum(x['analytic_cohort'] for x in updates),model_input_values=sum(x['analytic_cohort'] and x['model_input'] for x in updates),extra_nonanalytic_values=sum(not x['analytic_cohort'] for x in updates),all_original_observed_values_reconstructed=True,corrected_analytic_values_match_bootstrap_input=True,trial_checks=stats)
(audit/'timing_recheck.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2),flush=True)
