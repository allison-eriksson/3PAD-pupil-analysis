# 3PAD MATLAB preprocessing

These are the 3PAD study's early- and late-pregnancy MATLAB preprocessing scripts, adapted from the workflow supplied by Marcus Grueschow. The repository copies use configurable local paths, require local study inputs, and omit original comments and hard-coded participant identifiers. The archival originals are held with the study records.

| Visit | Command script | Raw processing | Eye analysis |
| --- | --- | --- | --- |
| Early pregnancy | `AE_EST_COMM_R1.m` | `AE_EST_RAW_R1.m` | `AE_EST_EYE_ANAL_R1.m` |
| Late pregnancy | `AE_EST_COMM_R2.m` | `AE_EST_RAW_R2.m` | `AE_EST_EYE_ANAL_R2.m` |

## Local setup for authorized users

Set `THREEPAD_EST_ROOT` to an authorized working copy of the study's EST directory. Set `THREEPAD_MATLAB_TOOLS` to its MATLAB tools directory. The scripts expect the `analysisScriptEyetrackingHelenaMarcB`, `edf-converter-master` and `eeglab2024.0` subdirectories there. The tools include Edf2Mat with the appropriate EyeLink EDF API for the platform, study helpers such as `VisAngMarc` and `createGauss`, and EEGLAB's `eegfilt`. The code also calls MATLAB toolbox functions. Obtain the required tools separately in the authorized environment.

For an authorized raw rerun, use a private working copy of the input tree because the historical scripts write converted files and intermediate outputs. To use the low-pass extraction branch, set `doExtraction` and `doLowPassFilteringOnly` to 1 and `doBandPassFiltering` to 0 in your local working copies of both eye-analysis scripts. The repository copies retain the archived flag values. The MATLAB scripts were reviewed as source, and the archived low-pass inputs were verified by checksum.

## Revision workflow

The MATLAB scripts document the historical preprocessing. For the revised analyses, use [the Python timing-correction workflow](../README.md#rebuild-the-corrected-input-from-saved-matlab-outputs) on the saved low-pass pupil and behavioral arrays. It corrects response timing for consecutive incongruent trials before the corrected workbook enters the R analysis. Participant lists, recordings, workbooks, and intermediate outputs belong in the authorized study environment.
