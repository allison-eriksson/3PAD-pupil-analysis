# 3PAD PLOS ONE revision analysis code

3 November 2026 · v01

The original eye-tracking preprocessing used the study's MATLAB scripts. The timing correction uses Python to recompute Dilation at Response and Net Pupil Change for consecutive incongruent trials from the saved MATLAB low-pass pupil traces, using the recorded reaction time. Other trial measures and the original interpolation and trial-retention rules carry forward. The merged workbook then enters the R analysis. The R workflow can be run from RStudio or the command line with the appropriate study inputs and software.

## Regenerate results in R from the corrected workbook and saved paired fits

Unpack this archive and open a terminal in the unpacked code directory. Set `PROJECT` to the authorized 3PAD revision project directory and `NEW_OUTPUT` to a new directory outside the existing results.

```sh
Rscript --vanilla reproduce.R "$PROJECT" "$NEW_OUTPUT"
```

This rebuilds the cohort and observed-data Tables 1a and 1b, recalculates the 90 observed-pair correlations, draws S2 and S3 Figs, regenerates Tables 2–4 and Fig 2 from the saved paired fits, and recalculates the sample-size sensitivity analysis. It checks the reproduced tables against the study project. It also exports cohort-flow and outcome-overlap counts. Fig 1 and S1 Fig are existing study illustrations supplied with the manuscript. The code exports the counts underlying the outcome overlap.

To rerun the separate 50-imputation correlation sensitivity analysis, use a different new output directory and add `mi`.

```sh
Rscript --vanilla reproduce.R "$PROJECT" "$NEW_MI_OUTPUT" mi
```

The archive includes the revision Quarto notebook and its R scripts. Use the command above for the isolated reproduction run. The Quarto notebook is the full working analysis record and uses the study project's usual paths and saved inputs.

### Inputs for the R run

Keep the study inputs in the authorized project. The paths below are relative to that project.

| Input | Path |
|---|---|
| Corrected merged workbook | `data/raw/3PADCombined_R_allSubs_25Oct26.xlsx` |
| Verified model-ready data | `data/derived/analysis_unimputed.rds` |
| Saved paired model results and apparent fits | `output/derived/<outcome>_final_bootstrap_1000x10.rds`, `<outcome>_apparent_baseline.rds`, `<outcome>_apparent_full.rds`, `<outcome>_apparent_idi.rds` for `epds`, `stais`, `stait` |
| Saved comparison outputs | `output/tables/master_*.csv`, `output/tables/pmsampsize_*.csv`, `output/table1_observed_audit/Table_1{a,b}_corrected.csv`, and `pd_correlations_all_90.csv` under both `output/correlations_observed/` and `output/correlations_mi_sensitivity/` |

## Refit the corrected models in R

The refit route uses the original saved baseline fits and refits the full models on the identical bootstrap draws. `ORIGINAL` is the archived pre-correction project. The runner reads its original typed RDS and paired fits. The original full EPDS fit is used as a numerical reproduction check before any corrected models are run. `CORRECTED_RDS` is the corrected model-ready RDS, usually `$PROJECT/data/derived/analysis_unimputed.rds`. Create a new private `REFIT` directory for the results. The per-resample outputs contain participant-level predictions.

```sh
mkdir -p "$REFIT/audit"
Rscript --vanilla refit/run_corrected.R "$REFIT" "$ORIGINAL" "$CORRECTED_RDS" check
Rscript --vanilla refit/run_corrected.R "$REFIT" "$ORIGINAL" "$CORRECTED_RDS" pilot epds
for outcome in epds stais stait; do
  Rscript --vanilla refit/run_corrected.R "$REFIT" "$ORIGINAL" "$CORRECTED_RDS" run "$outcome" || exit 1
done
Rscript --vanilla refit/repair_flagged_fits.R "$REFIT" "$ORIGINAL" "$CORRECTED_RDS" epds
Rscript --vanilla refit/repair_flagged_fits.R "$REFIT" "$ORIGINAL" "$CORRECTED_RDS" stais
Rscript --vanilla refit/check_calibration_convergence.R "$REFIT" "$ORIGINAL" "$CORRECTED_RDS"
Rscript --vanilla refit/prepare_calibration_repair.R "$REFIT" "$ORIGINAL" "$CORRECTED_RDS"
Rscript --vanilla refit/verify_numerical_repairs.R "$REFIT" "$ORIGINAL" "$CORRECTED_RDS"
THREEPAD_REFIT_DIR="$REFIT" Rscript --vanilla reproduce.R "$PROJECT" "$NEW_REFIT_OUTPUT" refit
```

The full fit uses 1,000 paired resamples per outcome, 10 imputations per resample and seed 26092026. Valid completed resamples are reused on restart. The three numerical safeguards reproduce the checked EPDS 220 and STAI-S 295 path refits and STAI-T 616 calibration-slope repair.

## Rebuild the corrected input from saved MATLAB outputs

`EST_RESULTS` contains `EYE_ANAL/RUN1/` and `EYE_ANAL/RUN2/`, each with the matching `allSubsAllBev_LowPass_R*.mat` and `allSubsAllPup_LowPass_R*.mat` files. The correction and workbook audit need Python with NumPy, SciPy and openpyxl, plus R with jsonlite. The raw-quality audit also needs h5py. Use new private output locations.

```sh
python3 pupil/rebuild_timing_corrected_data.py --est-results "$EST_RESULTS" \
  --original-data "$ORIGINAL/data/derived/analysis_unimputed.rds" --output-dir "$TIMING"
python3 pupil/recheck_workbook_values.py --est-results "$EST_RESULTS" \
  --original-workbook "$ORIGINAL/data/raw/3PADCombined_R_allSubs_14Mar25.xlsx" \
  --timing-output "$TIMING" --audit-dir "$WORKBOOK_AUDIT"
python3 pupil/write_corrected_workbook.py \
  --original-workbook "$ORIGINAL/data/raw/3PADCombined_R_allSubs_14Mar25.xlsx" \
  --audit-dir "$WORKBOOK_AUDIT" --output-workbook "$NEW_WORKBOOK"
```

For S3 Table, `pupil/audit_raw_quality.py --help` describes the raw audit inputs. It uses the original converted EyeLink files and saved processed structures, plus the retention CSV emitted by the timing script. `pupil/summarize_trial_quality.py RETENTION_CSV QUALITY_AUDIT_DIR NEW_PRIVATE_OUTPUT` exports the two sets of values behind S3 Table. These audits describe the recordings used in the analysis and impose no new trial exclusions.

## Software and access

See `environment/bootstrap_sessionInfo.txt` for the actual refit environment and `environment/packages.md` for the package versions used in verification. The included `renv.lock` is the working project's historical lock. Use the recorded package versions for numerical reproduction.

This code archive contains source code and documentation. Participant workbooks, recordings, typed data, per-resample predictions, and audit logs are held separately under the study's access restrictions. The original MATLAB snapshots, local path inventory, and upstream dependency notes are maintained with the author's private replication materials.
