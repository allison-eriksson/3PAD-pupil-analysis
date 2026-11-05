additional_library <- Sys.getenv("THREEPAD_ADDITIONAL_R_LIB")
if (nzchar(additional_library)) .libPaths(c(.libPaths(), additional_library))
# Run from the unpacked code directory with an isolated output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: Rscript --vanilla reproduce.R STUDY_PROJECT NEW_OUTPUT_DIR [refit|mi]")
}
project <- normalizePath(args[1], mustWork = TRUE)
output <- args[2]
mode <- if (length(args) == 3L) args[3] else "saved"
stopifnot(mode %in% c("saved", "refit", "mi"), !file.exists(output))
dir.create(output, recursive = TRUE)
output <- normalizePath(output, mustWork = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(stringr)
})
options(dplyr.summarise.inform = FALSE)
p <- file.path(output, "output", "tables")
dir.create(p, recursive = TRUE)
b <- file.path(output, "output", "figures")
dir.create(b, recursive = TRUE)

compare_csv <- function(actual, reference) {
  read <- function(path) read.csv(path, check.names = FALSE)
  delta <- all.equal(read(actual), read(reference), tolerance = 1e-9,
                     check.attributes = FALSE)
  if (!isTRUE(delta)) stop(basename(actual), ": ", paste(delta, collapse = "; "))
  message("Verified ", basename(actual))
}

source("refit/original_workflow_functions.R")
df_raw <- readxl::read_excel(file.path(project, "data/raw/3PADCombined_R_allSubs_25Oct26.xlsx"))
source("code/prepare_data.R")
stopifnot(isTRUE(all.equal(analysis_unimputed,
  readRDS(file.path(project, "data/derived/analysis_unimputed.rds")), tolerance = 1e-12)))
write.csv(cohort_flow, file.path(p, "cohort_flow.csv"), row.names = FALSE)
overlap <- df |> count(outcome_epds, outcome_stai_s, outcome_stai_t)
write.csv(overlap, file.path(p, "outcome_overlap_counts.csv"), row.names = FALSE)
source("code/rebuild_tables_1a_1b.R")
table1 <- file.path(output, "output", "table1_observed_audit")
invisible(rebuild_tables_1a_1b(df, table1))
for (name in c("Table_1a_corrected.csv", "Table_1b_corrected.csv")) {
  compare_csv(file.path(table1, name), file.path(project, "output/table1_observed_audit", name))
}

source("code/supplementary_pd_correlations.R")
observed <- file.path(output, "output", "correlations_observed")
mi <- file.path(output, "output", "correlations_mi_sensitivity")
run_pd_correlations(df_raw, observed, method = "observed_pairs")
correlation_csv <- "pd_correlations_all_90.csv"
compare_csv(file.path(observed, correlation_csv),
            file.path(project, "output/correlations_observed", correlation_csv))
if (mode == "mi") {
  run_pd_correlations(df_raw, mi, method = "multiple_imputation", m = 50L,
                      maxit = 20L, seed = 26102026L)
  compare_csv(file.path(mi, correlation_csv),
              file.path(project, "output/correlations_mi_sensitivity", correlation_csv))
} else {
  dir.create(mi, recursive = TRUE)
  stopifnot(file.copy(file.path(project, "output/correlations_mi_sensitivity", correlation_csv), mi))
}
source("code/rebuild_pd_presentation.R")
rebuild_pd_presentation(output, file.path(output, "output"), normalizePath("code"))

source("code/load_reporting_inputs.R")
load_reporting_inputs(project, if (mode == "refit") Sys.getenv("THREEPAD_REFIT_DIR") else NULL)
source("code/reporting_assignments.R")
for (name in c("master_model_performance", "master_incremental_value", "master_pd_stability_long",
               "master_pd_stability_table", "master_calibration_plot_data", "master_calibration_labels")) {
  filename <- paste0(name, ".csv")
  write.csv(get(name), file.path(p, filename), row.names = FALSE)
  compare_csv(file.path(p, filename), file.path(project, "output/tables", filename))
}
source("code/figure2.R")
source("code/sample_size.R")
for (name in c("pmsampsize_sensitivity_analysis.csv", "pmsampsize_primary_summary.csv")) {
  compare_csv(file.path(p, name), file.path(project, "output/tables", name))
}
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
message("Finished writing and verifying outputs in ", output)
