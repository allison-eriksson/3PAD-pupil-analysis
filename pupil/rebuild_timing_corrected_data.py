"""Reconstruct the four response-timing-corrected PD columns without MATLAB.

Requires Python with numpy/scipy and R with jsonlite. Source files are read only.
Only the output directory receives new files. No model fitting is performed.
"""
import argparse
import csv
import gc
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
from scipy.io import loadmat


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--est-results", type=Path, required=True)
parser.add_argument("--original-data", type=Path, required=True)
parser.add_argument("--rscript", default=shutil.which("Rscript"))
parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
args = parser.parse_args()
out = args.output_dir.resolve()
out.mkdir(parents=True, exist_ok=True)
for name in ["analysis_unimputed_corrected.csv", "analysis_unimputed_corrected.rds"]:
    if (out / name).exists():
        raise FileExistsError(f"Refusing to overwrite {out / name}; use a new output directory.")

# Read the original model-ready R data, including types, missingness and outcomes.
r_read = '''
d <- readRDS(commandArgs(TRUE)[1])
cat(jsonlite::toJSON(d, dataframe="rows", na="null", digits=17))
'''
raw = subprocess.check_output(
    [args.rscript, "--vanilla", "-e", r_read, str(args.original_data)], text=True
)
original = json.loads(raw)
cohort = {int(row["pid"]): row for row in original}
assert len(original) == len(cohort) == 91
corrections = {pid: {"pid": pid} for pid in cohort}
log, trial_counts = [], []
manifest = [{"file": str(args.original_data.resolve()), "sha256": sha256(args.original_data)}]
zero = 10000  # Zero-based stimulus-onset position in the saved 20,001-point epochs.

for visit in (1, 2):
    folder = args.est_results / "EYE_ANAL" / f"RUN{visit}"
    table_path = folder / f"allSubsAllBev_LowPass_R{visit}.mat"
    pupil_path = folder / f"allSubsAllPup_LowPass_R{visit}.mat"
    print(f"Reading R{visit} saved trial data and pupil epochs.", flush=True)
    table = loadmat(table_path)["allSubsAllBev"]
    pupil = loadmat(pupil_path)["allSubsAllPup"]
    assert table.shape[0] == pupil.shape[0] and pupil.shape[1] == 20001
    corrected_table = table.copy()
    # MATLAB columns: 2=response status, 3=congruency, 4=condition code,
    # 5=RT in ms, 16=sequence, 17=participant, 20=PAR, 21=ROD.
    affected = np.where(
        (table[:, 15] == 4) & (table[:, 1] != 99)
        & np.isin(table[:, 16], list(cohort))
    )[0]
    for row in affected:
        rt = int(table[row, 4])
        condition = int(table[row, 3])
        assert rt >= 250 and zero + rt + 250 < pupil.shape[1]
        trace = pupil[row]
        # Confirm the error in the original stored values before correcting it.
        original_par = np.mean(trace[zero + condition - 250:zero + condition + 251])
        original_rod = trace[zero + condition] - trace[zero]
        assert np.isclose(original_par, table[row, 19], atol=1e-10, rtol=0)
        assert np.isclose(original_rod, table[row, 20], atol=1e-10, rtol=0)
        corrected_table[row, 19] = np.mean(trace[zero + rt - 250:zero + rt + 251])
        # Preserve the original net-change quantity; do not divide by RT.
        corrected_table[row, 20] = trace[zero + rt] - trace[zero]

    for pid, original_row in cohort.items():
        all_rows = table[table[:, 16] == pid]
        if not len(all_rows):
            trial_counts.append({"pid": pid, "visit": visit, "record_available": False})
            continue
        old_rows = all_rows[all_rows[:, 1] != 99]
        new_rows = corrected_table[(corrected_table[:, 16] == pid) & (corrected_table[:, 1] != 99)]
        corrected_ii = int(np.sum(old_rows[:, 15] == 4))
        trial_counts.append({
            "pid": pid, "visit": visit, "record_available": True,
            "recorded": len(all_rows), "missed_response_excluded": int(np.sum(all_rows[:, 1] == 99)),
            "retained": len(old_rows), "incorrect_responses_retained": int(np.sum(old_rows[:, 1] == 0)),
            "congruent_retained": int(np.sum(old_rows[:, 2] == 1)),
            "incongruent_retained": int(np.sum(old_rows[:, 2] == 2)),
            "ii_trials_timing_corrected": corrected_ii,
        })
        for stem, column, name in [("PARDiff", 19, "Dilation at Response"), ("RodDiff", 20, "Net Pupil Change")]:
            variable = stem + ("2" if visit == 2 else "")
            before = np.mean(old_rows[old_rows[:, 2] == 2, column]) - np.mean(old_rows[old_rows[:, 2] == 1, column])
            after = np.mean(new_rows[new_rows[:, 2] == 2, column]) - np.mean(new_rows[new_rows[:, 2] == 1, column])
            assert abs(before - original_row[variable]) < 1e-8
            assert np.isfinite(after)
            corrections[pid][variable] = float(after)
            log.append({"pid": pid, "visit": visit, "variable": variable, "measure_name": name,
                        "original": original_row[variable], "corrected": float(after),
                        "change": float(after - original_row[variable]),
                        "ii_trials_timing_corrected": corrected_ii})
    del table, corrected_table, pupil
    gc.collect()
    for path in (table_path, pupil_path):
        manifest.append({"file": str(path.resolve()), "sha256": sha256(path)})

write_csv(out / "pd_timing_correction_log.csv", log,
          ["pid", "visit", "variable", "measure_name", "original", "corrected", "change", "ii_trials_timing_corrected"])
write_csv(out / "trial_retention_by_participant_visit.csv", trial_counts,
          ["pid", "visit", "record_available", "recorded", "missed_response_excluded", "retained",
           "incorrect_responses_retained", "congruent_retained", "incongruent_retained", "ii_trials_timing_corrected"])

r_write = r'''
args <- commandArgs(TRUE)
original <- readRDS(args[1]); corrected <- original
updates <- jsonlite::fromJSON(paste(readLines(file("stdin"),warn=FALSE),collapse="\n"))
changed <- c("PARDiff","PARDiff2","RodDiff","RodDiff2")
for(v in changed){
  ii <- match(corrected$pid,updates$pid)
  stopifnot(identical(is.na(original[[v]]),is.na(updates[[v]][ii])))
  corrected[[v]] <- updates[[v]][ii]
  stopifnot(identical(class(original[[v]]),class(corrected[[v]])))
}
unchanged <- setdiff(names(original),changed)
stopifnot(identical(original[unchanged],corrected[unchanged]),identical(is.na(original),is.na(corrected)))
stopifnot(identical(original$pid,corrected$pid),nrow(corrected)==91,ncol(corrected)==46)
saveRDS(corrected,file.path(args[2],"analysis_unimputed_corrected.rds"))
write.csv(corrected,file.path(args[2],"analysis_unimputed_corrected.csv"),row.names=FALSE,na="")
csv <- read.csv(file.path(args[2],"analysis_unimputed_corrected.csv"),check.names=FALSE,na.strings=c("","NA"))
for(v in names(corrected)){
 stopifnot(identical(is.na(csv[[v]]),is.na(corrected[[v]])))
 if(is.factor(corrected[[v]]))stopifnot(identical(as.character(csv[[v]]),as.character(corrected[[v]])))
 else stopifnot(isTRUE(all.equal(as.numeric(csv[[v]]),as.numeric(corrected[[v]]),tolerance=1e-12)))
}
description <- c(
 pid="Original study participant identifier; retained to link correction records. Not a newly assigned public-release identifier.",
 gb1_age="Age at early-pregnancy assessment, in years.",
 country="Country-of-birth grouping from the original R analysis: Sweden versus Other.",
 bmi="Early-pregnancy weight in kg divided by squared height in metres.",
 edu="University education versus Less, using original REDCap coding.",
 first_preg="First pregnancy (Yes/No); this is not parity.",
 workFT="Any full-time employment, self-employment or study versus Other, using original R coding.",
 pcos="PCOS self-report; Yes/No using the original R recode.",
 exercise1="Early-pregnancy exercise: Low to none versus Moderate to high, using the original R recode.",
 exercise2="Late-pregnancy exercise: Low to none versus Moderate to high, using the original R recode.",
 SSRI="Early-pregnancy SSRI use: Yes/No.",
 pos_LE1="Positive life experience in early pregnancy: Yes/No using original R coding.",
 pos_LE2="Positive life experience in late pregnancy: Yes/No using original R coding.",
 neg_LE1="Negative life experience in early pregnancy: Yes/No using original R coding.",
 neg_LE2="Negative life experience in late pregnancy: Yes/No using original R coding.",
 outcome_epds="Six-week postpartum EPDS: 1 if score >=12, otherwise 0; outcome not imputed.",
 outcome_stai_s="Six-week postpartum STAI-S: 1 if score >=40, otherwise 0; outcome not imputed.",
 outcome_stai_t="Six-week postpartum STAI-T: 1 if score >=40, otherwise 0; outcome not imputed."
)
scales <- c(epds="Edinburgh Postnatal Depression Scale",stai_s="State-Trait Anxiety Inventory State",stai_t="State-Trait Anxiety Inventory Trait",ders16="Difficulties in Emotion Regulation Scale 16-item",saaq="Swedish Acceptance and Action Questionnaire",rs14="Resilience Scale 14-item",pss="Perceived Stress Scale",soc13="Sense of Coherence Scale 13-item",gad7="Generalized Anxiety Disorder 7-item scale")
for(v in grep("^gb[12]_.*_sum$",names(corrected),value=TRUE)){
 code <- sub("_sum$","",sub("^gb[12]_","",v))
 visit <- if(startsWith(v,"gb1_"))"Early pregnancy"else"Late pregnancy"
 description[v] <- paste0(visit," ",scales[[code]]," sum score, preserved from original analysis data.")
}
pd_stems <- c(PARDiff="Dilation at Response",RodDiff="Net Pupil Change",PCUDiff="Early Peak Pupil Change",PTRDiff="Early Mean Pupil Size",PAADiff="Late Sustained Dilation")
pd_defs <- c(PARDiff="Mean z-scored pupil size from 250 ms before to 250 ms after the recorded response. II-trial response timing corrected.",RodDiff="Net z-scored pupil change from stimulus onset to the recorded response, equivalent to the sum of successive changes. Not divided by response time. II-trial response timing corrected.",PCUDiff="Mean z-scored pupil size over 500-2000 ms after stimulus onset. Original name retained; this is not a maximum or a baseline-subtracted peak.",PTRDiff="Mean z-scored pupil size over 500-850 ms after stimulus onset.",PAADiff="Mean z-scored pupil size over 2128-4138 ms after stimulus onset.")
for(stem in names(pd_stems))for(suffix in c("","2")){
 v <- paste0(stem,suffix);visit <- if(suffix=="")"Early pregnancy"else"Late pregnancy"
 description[v] <- paste0(visit,": ",pd_stems[[stem]],". ",pd_defs[[stem]]," Participant-level incongruent mean minus congruent mean; missed responses (status 99) excluded, incorrect responses retained.")
}
stopifnot(all(names(corrected)%in%names(description)))
dictionary <- data.frame(variable=names(corrected),description=unname(description[names(corrected)]),R_class=vapply(corrected,function(x)paste(class(x),collapse="/"),character(1)),factor_levels=vapply(corrected,function(x)if(is.factor(x))paste(levels(x),collapse="; ")else"",character(1)),n_missing=vapply(corrected,function(x)sum(is.na(x)),integer(1)),timing_corrected=names(corrected)%in%changed,stringsAsFactors=FALSE)
write.csv(dictionary,file.path(args[2],"data_dictionary.csv"),row.names=FALSE,na="")
verification <- list(rows=nrow(corrected),columns=ncol(corrected),changed_columns=changed,unchanged_columns_verified=length(unchanged),participant_order_preserved=TRUE,missingness_preserved=TRUE,factor_levels_preserved=TRUE,CSV_round_trip_verified=TRUE,observed_values_updated=vapply(changed,function(v)sum(!is.na(corrected[[v]])),integer(1)),outcome_events=c(EPDS=sum(corrected$outcome_epds),STAI_S=sum(corrected$outcome_stai_s),STAI_T=sum(corrected$outcome_stai_t)),models_rerun=FALSE)
cat(jsonlite::toJSON(verification,auto_unbox=TRUE,pretty=TRUE))
'''
verified = subprocess.check_output(
    [args.rscript, "--vanilla", "-e", r_write, str(args.original_data), str(out)],
    input=json.dumps(list(corrections.values())), text=True,
)
(out / "verification.json").write_text(verified + "\n", encoding="utf-8")
manifest.append({"file": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__))})
(out / "source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
print(verified, flush=True)
print(f"New corrected files saved in {out}. Original sources were read only.", flush=True)
