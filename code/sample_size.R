# 14. Formal sample-size adequacy assessment

sample_size_inputs <- tibble::tribble(
  ~outcome, ~events, ~prevalence,
  "EPDS",   13L,     13/91,
  "STAI-S", 12L,     12/91,
  "STAI-T", 21L,     21/91
)

candidate_cstats <- c(
  0.70,
  0.75,
  0.80
)

sample_size_scenarios <-
  tidyr::crossing(
    sample_size_inputs,
    model = c("baseline", "full"),
    cstatistic = candidate_cstats
  ) %>%
  mutate(
    parameters = if_else(
      model == "baseline",
      32L,
      42L
    )
  )

sample_size_scenarios

pmsampsize_results <-
  purrr::pmap(
    sample_size_scenarios,
    function(
      outcome,
      events,
      prevalence,
      model,
      cstatistic,
      parameters
    ) {

      pmsampsize::pmsampsize(
        type = "b",
        cstatistic = cstatistic,
        parameters = parameters,
        prevalence = prevalence,
        shrinkage = 0.90
      )
    }
  )

pmsampsize_summary <-
  purrr::map2_dfr(
    seq_len(nrow(sample_size_scenarios)),
    pmsampsize_results,
    function(i, result) {

      scenario <-
        sample_size_scenarios[i, ]

      tibble(
        outcome =
          scenario$outcome,

        model =
          scenario$model,

        parameters =
          scenario$parameters,

        prevalence =
          scenario$prevalence,

        observed_events =
          scenario$events,

        assumed_cstatistic =
          scenario$cstatistic,

        estimated_cox_snell_r2 =
          result$rsquared,

        minimum_sample_size =
          result$sample_size,

        required_events =
          ceiling(result$events),

        required_EPP =
          result$EPP,

        final_shrinkage =
          result$final_shrinkage,

        observed_sample_size =
          91L,

        sample_size_shortfall =
          result$sample_size - 91L,

        observed_EPP =
          scenario$events /
          scenario$parameters
      )
    }
  )

pmsampsize_summary


readr::write_csv(
  pmsampsize_summary,
  file.path(p,
    "pmsampsize_sensitivity_analysis.csv"
  )
)

pmsampsize_primary_summary <-
  pmsampsize_summary %>%
  filter(assumed_cstatistic == 0.75) %>%
  select(
    outcome,
    model,
    parameters,
    prevalence,
    observed_events,
    observed_EPP,
    minimum_sample_size,
    required_events,
    required_EPP
  )

pmsampsize_primary_summary




readr::write_csv(
  pmsampsize_primary_summary,
  file.path(p,
    "pmsampsize_primary_summary.csv"
  )
)
