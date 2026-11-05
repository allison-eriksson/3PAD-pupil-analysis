load_reporting_inputs <- function(project, refit = NULL) {
  if (!is.null(refit)) {
    stopifnot(nzchar(refit))
    refit <- normalizePath(refit, mustWork = TRUE)
  }
  derived <- file.path(project, "output", "derived")
  targets <- c(epds = 220L, stais = 295L, stait = 616L)
  for (key in names(targets)) {
    get_saved <- function(suffix) readRDS(file.path(derived, paste0(key, "_", suffix, ".rds")))
    baseline <- get_saved("apparent_baseline")
    if (is.null(refit)) {
      full <- get_saved("apparent_full")
      idi <- get_saved("apparent_idi")
      combined <- get_saved("final_bootstrap_1000x10")
    } else {
      full <- readRDS(file.path(refit, paste0(key, "_apparent_full_corrected.rds")))
      idi <- calculate_idi(baseline$predictions$outcome, baseline$predictions$predicted_prob,
                           full$predictions$predicted_prob)
      resamples <- lapply(seq_len(1000L), function(b) {
        folder <- if (b == targets[[key]]) "numerical_repairs" else "bootstrap"
        result <- readRDS(file.path(refit, folder, key, sprintf("%s_bootstrap_%04d.rds", key, b)))
        stopifnot(result$comparison$bootstrap == b)
        result
      })
      parts <- c("performance", "comparison", "diagnostics", "tuning", "coefficients")
      combined <- setNames(lapply(parts, function(n) bind_rows(lapply(resamples, "[[", n))), parts)
      combined$predictions <- map_dfr(seq_along(resamples), function(b) {
        map_dfr(c("baseline", "full"), function(model) {
          map_dfr(c("bootstrap", "validation"), function(sample) {
            resamples[[b]]$predictions[[model]][[sample]] |>
              mutate(bootstrap = b, model = model, sample = sample)
          })
        })
      })
      combined$bootstrap_indices <- tibble(bootstrap = seq_along(resamples),
                                         bootstrap_index = lapply(resamples, "[[", "bootstrap_index"))
    }
    stopifnot(identical(sort(unique(combined$performance$bootstrap)), seq_len(1000L)))
    for (suffix in c("apparent_baseline", "apparent_full", "apparent_idi", "final_bootstrap")) {
      value <- switch(suffix, apparent_baseline = baseline, apparent_full = full,
                      apparent_idi = idi, final_bootstrap = combined)
      assign(paste0(key, "_", suffix), value, envir = parent.frame())
    }
  }
}
