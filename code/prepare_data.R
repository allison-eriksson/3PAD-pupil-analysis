
#| label: define-analysis-cohort

# Define participants with all three observed postpartum outcomes
df_raw$has_postpartum_outcomes <-
  !is.na(df_raw$ppv6_epds_sum) &
  !is.na(df_raw$ppv6_stai_s_sum) &
  !is.na(df_raw$ppv6_stai_t_sum)

# Participant-flow summary
cohort_flow <- tibble(
  category = c(
    "Raw records",
    "Postpartum outcomes unavailable",
    "Final analytic cohort"
  ),
  n = c(
    nrow(df_raw),
    sum(!df_raw$has_postpartum_outcomes),
    sum(df_raw$has_postpartum_outcomes)
  )
)

cohort_flow

# Create analytic cohort
df <- df_raw %>%
  filter(has_postpartum_outcomes)

# Reproducibility checks
stopifnot(nrow(df_raw) == 99)
stopifnot(nrow(df) == 91)

cat("Analytic cohort:", nrow(df), "participants\n")


#| label: construct-background-variables

df <- df %>%
  mutate(

    # BMI
    gb1_weight = as.numeric(gb1_weight),
    bmi = gb1_weight / ((gb1_height / 100)^2),

    # Country of birth: Sweden vs other
    country = case_when(
      is.na(gb1_country) ~ NA_character_,
      str_detect(gb1_country, regex("sverige", ignore_case = TRUE)) ~ "Sweden",
      TRUE ~ "Other"
    ),
    country = factor(country, levels = c("Other", "Sweden")),

    # Education
    edu = case_when(
      is.na(gb1_edu___4) ~ NA_character_,
      gb1_edu___4 == 1 ~ "University",
      TRUE ~ "Less"
    ),
    edu = factor(edu, levels = c("Less", "University")),

    # First pregnancy - this is NOT parity
    first_preg = case_when(
      is.na(gb1_first_preg) ~ NA_character_,
      gb1_first_preg == 1 ~ "Yes",
      gb1_first_preg == 0 ~ "No"
    ),
    first_preg = factor(first_preg, levels = c("No", "Yes")),

    # Full-time employment, self-employment, or study
    workFT = if_else(
      if_any(
        c(
          gb1_job_work___1,
          gb1_job_business___1,
          gb1_job_student___1
        ),
        ~ . > 0
      ),
      "Full-time",
      "Other",
      missing = NA_character_
    ),
    workFT = factor(workFT, levels = c("Other", "Full-time")),

    # PCOS
    # REDCap code 1 = yes; other observed response categories collapsed to no,
    # matching the original analysis
    pcos = case_when(
      is.na(gb1_pcos) ~ NA_character_,
      gb1_pcos == 1 ~ "Yes",
      TRUE ~ "No"
    ),
    pcos = factor(pcos, levels = c("No", "Yes")),

    # Exercise
    # Original coding: 1-3 = moderate/high; 4-5 = low/none
    exercise1 = case_when(
      is.na(gb1_exercise) ~ NA_character_,
      gb1_exercise >= 4 ~ "Low to none",
      TRUE ~ "Moderate to high"
    ),
    exercise1 = factor(
      exercise1,
      levels = c("Low to none", "Moderate to high")
    ),

    exercise2 = case_when(
      is.na(gb2_exercise) ~ NA_character_,
      gb2_exercise >= 4 ~ "Low to none",
      TRUE ~ "Moderate to high"
    ),
    exercise2 = factor(
      exercise2,
      levels = c("Low to none", "Moderate to high")
    ),

    # SSRI use during early pregnancy
    SSRI = case_when(
      is.na(gb1_meds_ja___6) ~ NA_character_,
      gb1_meds_ja___6 == 1 ~ "Yes",
      TRUE ~ "No"
    ),
    SSRI = factor(SSRI, levels = c("No", "Yes"))
  )



df <- df %>%
  mutate(

    # Positive life experiences
    pos_LE1 = case_when(
      is.na(gb1_pos_lifeevents) ~ NA_character_,
      gb1_pos_lifeevents == 1 ~ "No",
      TRUE ~ "Yes"
    ),
    pos_LE1 = factor(pos_LE1, levels = c("No", "Yes")),

    pos_LE2 = case_when(
      is.na(gb2_pos_lifeevents) ~ NA_character_,
      gb2_pos_lifeevents == 1 ~ "No",
      TRUE ~ "Yes"
    ),
    pos_LE2 = factor(pos_LE2, levels = c("No", "Yes")),

    # Negative life experiences
    neg_LE1 = case_when(
      is.na(gb1_neg_lifeevents) ~ NA_character_,
      gb1_neg_lifeevents == 1 ~ "No",
      TRUE ~ "Yes"
    ),
    neg_LE1 = factor(neg_LE1, levels = c("No", "Yes")),

    neg_LE2 = case_when(
      is.na(gb2_neg_lifeevents) ~ NA_character_,
      gb2_neg_lifeevents == 1 ~ "No",
      TRUE ~ "Yes"
    ),
    neg_LE2 = factor(neg_LE2, levels = c("No", "Yes"))
  )



#| label: construct-outcomes

df <- df %>%
  mutate(
    outcome_epds = as.integer(ppv6_epds_sum >= 12),
    outcome_stai_s = as.integer(ppv6_stai_s_sum >= 40),
    outcome_stai_t = as.integer(ppv6_stai_t_sum >= 40)
  )

# Confirm outcome counts
outcome_check <- tibble(
  outcome = c(
    "EPDS >= 12",
    "STAI-S >= 40",
    "STAI-T >= 40"
  ),
  low = c(
    sum(df$outcome_epds == 0),
    sum(df$outcome_stai_s == 0),
    sum(df$outcome_stai_t == 0)
  ),
  high = c(
    sum(df$outcome_epds == 1),
    sum(df$outcome_stai_s == 1),
    sum(df$outcome_stai_t == 1)
  )
)

outcome_check


#| label: create-unimputed-analysis-data

analysis_unimputed <- df %>%
  select(
    pid,
    all_of(full_vars),
    outcome_epds,
    outcome_stai_s,
    outcome_stai_t
  )

stopifnot(nrow(analysis_unimputed) == 91)
stopifnot(length(full_vars) == 42)



