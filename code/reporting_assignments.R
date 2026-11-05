epds_final_results_new <- calculate_optimism_corrected(apparent_baseline = epds_apparent_baseline, apparent_full = epds_apparent_full, 
    apparent_idi = epds_apparent_idi, paired_bootstrap = epds_final_bootstrap)

epds_selection_by_bootstrap_new <- epds_final_bootstrap$coefficients %>% group_by(model, bootstrap, variable) %>% 
    summarise(n_selected = sum(selected), proportion_imputations_selected = mean(selected), n_positive = sum(selected & 
        coefficient > 0), n_negative = sum(selected & coefficient < 0), .groups = "drop") %>% mutate(majority_selected = n_selected >= 
    5, predominant_direction = case_when(n_positive > n_negative ~ "positive", n_negative > n_positive ~ "negative", 
    TRUE ~ "mixed"))

epds_bootstrap_stability_new <- epds_selection_by_bootstrap_new %>% group_by(model, variable) %>% summarise(mean_imputation_selection = mean(proportion_imputations_selected), 
    bootstrap_majority_frequency = mean(majority_selected), n_bootstraps_majority_selected = sum(majority_selected), 
    positive_when_majority_selected = sum(majority_selected & predominant_direction == "positive"), negative_when_majority_selected = sum(majority_selected & 
        predominant_direction == "negative"), .groups = "drop") %>% mutate(predominant_bootstrap_direction = case_when(positive_when_majority_selected > 
    negative_when_majority_selected ~ "positive", negative_when_majority_selected > positive_when_majority_selected ~ 
    "negative", TRUE ~ "mixed"), bootstrap_direction_consistency = if_else(n_bootstraps_majority_selected > 0, 
    pmax(positive_when_majority_selected, negative_when_majority_selected)/n_bootstraps_majority_selected, NA_real_)) %>% 
    arrange(model, desc(bootstrap_majority_frequency))

epds_pd_bootstrap_stability_new <- epds_bootstrap_stability_new %>% filter(model == "full", variable %in% pd_vars) %>% 
    arrange(desc(bootstrap_majority_frequency))

epds_incremental_variability <- epds_final_bootstrap$comparison %>% summarise(delta_auc_validation_mean = mean(delta_auc_validation), 
    delta_auc_validation_median = median(delta_auc_validation), delta_auc_validation_p025 = quantile(delta_auc_validation, 
        0.025), delta_auc_validation_p975 = quantile(delta_auc_validation, 0.975), proportion_delta_auc_positive = mean(delta_auc_validation > 
        0), idi_validation_mean = mean(idi_validation), idi_validation_median = median(idi_validation), idi_validation_p025 = quantile(idi_validation, 
        0.025), idi_validation_p975 = quantile(idi_validation, 0.975), proportion_idi_positive = mean(idi_validation > 
        0))

epds_optimism_variability <- epds_final_bootstrap$comparison %>% summarise(optimism_delta_auc_mean = mean(optimism_delta_auc), 
    optimism_delta_auc_p025 = quantile(optimism_delta_auc, 0.025), optimism_delta_auc_p975 = quantile(optimism_delta_auc, 
        0.975), optimism_idi_mean = mean(optimism_idi), optimism_idi_p025 = quantile(optimism_idi, 0.025), optimism_idi_p975 = quantile(optimism_idi, 
        0.975))

epds_calibration_grids <- list(baseline = seq(quantile(epds_apparent_baseline$predictions$predicted_prob, 0.025), 
    quantile(epds_apparent_baseline$predictions$predicted_prob, 0.975), length.out = 101), full = seq(quantile(epds_apparent_full$predictions$predicted_prob, 
    0.025), quantile(epds_apparent_full$predictions$predicted_prob, 0.975), length.out = 101))

epds_apparent_calibration_curves <- bind_rows(estimate_calibration_curve(outcome = epds_apparent_baseline$predictions$outcome, 
    predicted_prob = epds_apparent_baseline$predictions$predicted_prob, grid = epds_calibration_grids$baseline) %>% 
    mutate(model = "baseline", .before = 1), estimate_calibration_curve(outcome = epds_apparent_full$predictions$outcome, 
    predicted_prob = epds_apparent_full$predictions$predicted_prob, grid = epds_calibration_grids$full) %>% mutate(model = "full", 
    .before = 1))

epds_calibration_bootstrap <- purrr::map_dfr(c("baseline", "full"), function(model_name) {
    grid <- epds_calibration_grids[[model_name]]
    purrr::map_dfr(1:1000, function(b) {
        bootstrap_predictions <- epds_final_bootstrap$predictions %>% filter(bootstrap == b, model == model_name, 
            sample == "bootstrap")
        validation_predictions <- epds_final_bootstrap$predictions %>% filter(bootstrap == b, model == model_name, 
            sample == "validation")
        bootstrap_curve <- estimate_calibration_curve(outcome = bootstrap_predictions$outcome, predicted_prob = bootstrap_predictions$predicted_prob, 
            grid = grid)
        validation_curve <- estimate_calibration_curve(outcome = validation_predictions$outcome, predicted_prob = validation_predictions$predicted_prob, 
            grid = grid)
        tibble(bootstrap = b, model = model_name, predicted_prob = grid, calibration_bootstrap = bootstrap_curve$observed_prob, 
            calibration_validation = validation_curve$observed_prob, optimism = calibration_bootstrap - calibration_validation)
    })
})

epds_calibration_optimism <- epds_calibration_bootstrap %>% group_by(model, predicted_prob) %>% summarise(mean_optimism = mean(optimism, 
    na.rm = TRUE), n_contributing = sum(!is.na(optimism)), .groups = "drop")

epds_calibration_corrected <- epds_apparent_calibration_curves %>% left_join(epds_calibration_optimism, by = c("model", 
    "predicted_prob")) %>% mutate(corrected_observed_prob_raw = observed_prob - mean_optimism, corrected_observed_prob = pmin(pmax(corrected_observed_prob_raw, 
    0), 1), proportion_bootstraps_contributing = n_contributing/1000)

epds_calibration_plot_data <- epds_calibration_corrected %>% filter(proportion_bootstraps_contributing >= 0.8) %>% 
    mutate(model = factor(model, levels = c("baseline", "full"), labels = c("Baseline model", "Full model (+ PD)")))

epds_calibration_labels <- epds_final_results_new$model_performance %>% transmute(model = factor(model, levels = c("baseline", 
    "full"), labels = c("Baseline model", "Full model (+ PD)")), label = sprintf("Calibration intercept = %.2f\nCalibration slope = %.2f", 
    calibration_intercept_optimism_corrected, calibration_slope_optimism_corrected))

epds_final_results <- epds_final_results_new

epds_bootstrap_stability <- epds_bootstrap_stability_new

epds_pd_bootstrap_stability <- epds_pd_bootstrap_stability_new

stais_final_results_new <- calculate_optimism_corrected(apparent_baseline = stais_apparent_baseline, apparent_full = stais_apparent_full, 
    apparent_idi = stais_apparent_idi, paired_bootstrap = stais_final_bootstrap)

stais_final_results <- stais_final_results_new

stais_incremental_variability <- stais_final_bootstrap$comparison %>% summarise(delta_auc_validation_mean = mean(delta_auc_validation), 
    delta_auc_validation_median = median(delta_auc_validation), delta_auc_validation_p025 = quantile(delta_auc_validation, 
        0.025), delta_auc_validation_p975 = quantile(delta_auc_validation, 0.975), proportion_delta_auc_positive = mean(delta_auc_validation > 
        0), idi_validation_mean = mean(idi_validation), idi_validation_median = median(idi_validation), idi_validation_p025 = quantile(idi_validation, 
        0.025), idi_validation_p975 = quantile(idi_validation, 0.975), proportion_idi_positive = mean(idi_validation > 
        0))

stais_optimism_variability <- stais_final_bootstrap$comparison %>% summarise(optimism_delta_auc_mean = mean(optimism_delta_auc), 
    optimism_delta_auc_p025 = quantile(optimism_delta_auc, 0.025), optimism_delta_auc_p975 = quantile(optimism_delta_auc, 
        0.975), optimism_idi_mean = mean(optimism_idi), optimism_idi_p025 = quantile(optimism_idi, 0.025), optimism_idi_p975 = quantile(optimism_idi, 
        0.975))

stais_selection_by_bootstrap <- stais_final_bootstrap$coefficients %>% group_by(model, bootstrap, variable) %>% 
    summarise(n_selected = sum(selected), proportion_imputations_selected = mean(selected), n_positive = sum(selected & 
        coefficient > 0), n_negative = sum(selected & coefficient < 0), .groups = "drop") %>% mutate(majority_selected = n_selected >= 
    5, predominant_direction = case_when(n_positive > n_negative ~ "positive", n_negative > n_positive ~ "negative", 
    TRUE ~ "mixed"))

stais_bootstrap_stability <- stais_selection_by_bootstrap %>% group_by(model, variable) %>% summarise(mean_imputation_selection = mean(proportion_imputations_selected), 
    bootstrap_majority_frequency = mean(majority_selected), n_bootstraps_majority_selected = sum(majority_selected), 
    positive_when_majority_selected = sum(majority_selected & predominant_direction == "positive"), negative_when_majority_selected = sum(majority_selected & 
        predominant_direction == "negative"), .groups = "drop") %>% mutate(predominant_bootstrap_direction = case_when(positive_when_majority_selected > 
    negative_when_majority_selected ~ "positive", negative_when_majority_selected > positive_when_majority_selected ~ 
    "negative", TRUE ~ "mixed"), bootstrap_direction_consistency = if_else(n_bootstraps_majority_selected > 0, 
    pmax(positive_when_majority_selected, negative_when_majority_selected)/n_bootstraps_majority_selected, NA_real_)) %>% 
    arrange(model, desc(bootstrap_majority_frequency))

stais_pd_bootstrap_stability <- stais_bootstrap_stability %>% filter(model == "full", variable %in% pd_vars) %>% 
    arrange(desc(bootstrap_majority_frequency))

stais_calibration_grids <- list(baseline = seq(quantile(stais_apparent_baseline$predictions$predicted_prob, 0.025), 
    quantile(stais_apparent_baseline$predictions$predicted_prob, 0.975), length.out = 101), full = seq(quantile(stais_apparent_full$predictions$predicted_prob, 
    0.025), quantile(stais_apparent_full$predictions$predicted_prob, 0.975), length.out = 101))

stais_apparent_calibration_curves <- bind_rows(estimate_calibration_curve(outcome = stais_apparent_baseline$predictions$outcome, 
    predicted_prob = stais_apparent_baseline$predictions$predicted_prob, grid = stais_calibration_grids$baseline) %>% 
    mutate(model = "baseline", .before = 1), estimate_calibration_curve(outcome = stais_apparent_full$predictions$outcome, 
    predicted_prob = stais_apparent_full$predictions$predicted_prob, grid = stais_calibration_grids$full) %>% mutate(model = "full", 
    .before = 1))

stais_calibration_bootstrap <- purrr::map_dfr(c("baseline", "full"), function(model_name) {
    grid <- stais_calibration_grids[[model_name]]
    purrr::map_dfr(1:1000, function(b) {
        bootstrap_predictions <- stais_final_bootstrap$predictions %>% filter(bootstrap == b, model == model_name, 
            sample == "bootstrap")
        validation_predictions <- stais_final_bootstrap$predictions %>% filter(bootstrap == b, model == model_name, 
            sample == "validation")
        bootstrap_curve <- estimate_calibration_curve(outcome = bootstrap_predictions$outcome, predicted_prob = bootstrap_predictions$predicted_prob, 
            grid = grid)
        validation_curve <- estimate_calibration_curve(outcome = validation_predictions$outcome, predicted_prob = validation_predictions$predicted_prob, 
            grid = grid)
        tibble(bootstrap = b, model = model_name, predicted_prob = grid, calibration_bootstrap = bootstrap_curve$observed_prob, 
            calibration_validation = validation_curve$observed_prob, optimism = calibration_bootstrap - calibration_validation)
    })
})

stais_calibration_optimism <- stais_calibration_bootstrap %>% group_by(model, predicted_prob) %>% summarise(mean_optimism = mean(optimism, 
    na.rm = TRUE), n_contributing = sum(!is.na(optimism)), .groups = "drop")

stais_calibration_corrected <- stais_apparent_calibration_curves %>% left_join(stais_calibration_optimism, by = c("model", 
    "predicted_prob")) %>% mutate(corrected_observed_prob_raw = observed_prob - mean_optimism, corrected_observed_prob = pmin(pmax(corrected_observed_prob_raw, 
    0), 1), proportion_bootstraps_contributing = n_contributing/1000)

stais_calibration_plot_data <- stais_calibration_corrected %>% filter(proportion_bootstraps_contributing >= 0.8) %>% 
    mutate(model = factor(model, levels = c("baseline", "full"), labels = c("Baseline model", "Full model (+ PD)")))

stais_calibration_labels <- stais_final_results$model_performance %>% transmute(model = factor(model, levels = c("baseline", 
    "full"), labels = c("Baseline model", "Full model (+ PD)")), label = sprintf("Calibration intercept = %.2f\nCalibration slope = %.2f", 
    calibration_intercept_optimism_corrected, calibration_slope_optimism_corrected))

stait_final_results <- calculate_optimism_corrected(apparent_baseline = stait_apparent_baseline, apparent_full = stait_apparent_full, 
    apparent_idi = stait_apparent_idi, paired_bootstrap = stait_final_bootstrap)

stait_incremental_variability <- stait_final_bootstrap$comparison %>% summarise(delta_auc_validation_mean = mean(delta_auc_validation), 
    delta_auc_validation_median = median(delta_auc_validation), delta_auc_validation_p025 = quantile(delta_auc_validation, 
        0.025), delta_auc_validation_p975 = quantile(delta_auc_validation, 0.975), proportion_delta_auc_positive = mean(delta_auc_validation > 
        0), idi_validation_mean = mean(idi_validation), idi_validation_median = median(idi_validation), idi_validation_p025 = quantile(idi_validation, 
        0.025), idi_validation_p975 = quantile(idi_validation, 0.975), proportion_idi_positive = mean(idi_validation > 
        0))

stait_optimism_variability <- stait_final_bootstrap$comparison %>% summarise(optimism_delta_auc_mean = mean(optimism_delta_auc), 
    optimism_delta_auc_p025 = quantile(optimism_delta_auc, 0.025), optimism_delta_auc_p975 = quantile(optimism_delta_auc, 
        0.975), optimism_idi_mean = mean(optimism_idi), optimism_idi_p025 = quantile(optimism_idi, 0.025), optimism_idi_p975 = quantile(optimism_idi, 
        0.975))

stait_selection_by_bootstrap <- stait_final_bootstrap$coefficients %>% group_by(model, bootstrap, variable) %>% 
    summarise(n_selected = sum(selected), proportion_imputations_selected = mean(selected), n_positive = sum(selected & 
        coefficient > 0), n_negative = sum(selected & coefficient < 0), .groups = "drop") %>% mutate(majority_selected = n_selected >= 
    5, predominant_direction = case_when(n_positive > n_negative ~ "positive", n_negative > n_positive ~ "negative", 
    TRUE ~ "mixed"))

stait_bootstrap_stability <- stait_selection_by_bootstrap %>% group_by(model, variable) %>% summarise(mean_imputation_selection = mean(proportion_imputations_selected), 
    bootstrap_majority_frequency = mean(majority_selected), n_bootstraps_majority_selected = sum(majority_selected), 
    positive_when_majority_selected = sum(majority_selected & predominant_direction == "positive"), negative_when_majority_selected = sum(majority_selected & 
        predominant_direction == "negative"), .groups = "drop") %>% mutate(predominant_bootstrap_direction = case_when(positive_when_majority_selected > 
    negative_when_majority_selected ~ "positive", negative_when_majority_selected > positive_when_majority_selected ~ 
    "negative", TRUE ~ "mixed"), bootstrap_direction_consistency = if_else(n_bootstraps_majority_selected > 0, 
    pmax(positive_when_majority_selected, negative_when_majority_selected)/n_bootstraps_majority_selected, NA_real_)) %>% 
    arrange(model, desc(bootstrap_majority_frequency))

stait_pd_bootstrap_stability <- stait_bootstrap_stability %>% filter(model == "full", variable %in% pd_vars) %>% 
    arrange(desc(bootstrap_majority_frequency))

stait_calibration_grids <- list(baseline = seq(quantile(stait_apparent_baseline$predictions$predicted_prob, 0.025), 
    quantile(stait_apparent_baseline$predictions$predicted_prob, 0.975), length.out = 101), full = seq(quantile(stait_apparent_full$predictions$predicted_prob, 
    0.025), quantile(stait_apparent_full$predictions$predicted_prob, 0.975), length.out = 101))

stait_apparent_calibration_curves <- bind_rows(estimate_calibration_curve(outcome = stait_apparent_baseline$predictions$outcome, 
    predicted_prob = stait_apparent_baseline$predictions$predicted_prob, grid = stait_calibration_grids$baseline) %>% 
    mutate(model = "baseline", .before = 1), estimate_calibration_curve(outcome = stait_apparent_full$predictions$outcome, 
    predicted_prob = stait_apparent_full$predictions$predicted_prob, grid = stait_calibration_grids$full) %>% mutate(model = "full", 
    .before = 1))

stait_calibration_bootstrap <- purrr::map_dfr(c("baseline", "full"), function(model_name) {
    grid <- stait_calibration_grids[[model_name]]
    purrr::map_dfr(1:1000, function(b) {
        bootstrap_predictions <- stait_final_bootstrap$predictions %>% filter(bootstrap == b, model == model_name, 
            sample == "bootstrap")
        validation_predictions <- stait_final_bootstrap$predictions %>% filter(bootstrap == b, model == model_name, 
            sample == "validation")
        bootstrap_curve <- estimate_calibration_curve(outcome = bootstrap_predictions$outcome, predicted_prob = bootstrap_predictions$predicted_prob, 
            grid = grid)
        validation_curve <- estimate_calibration_curve(outcome = validation_predictions$outcome, predicted_prob = validation_predictions$predicted_prob, 
            grid = grid)
        tibble(bootstrap = b, model = model_name, predicted_prob = grid, calibration_bootstrap = bootstrap_curve$observed_prob, 
            calibration_validation = validation_curve$observed_prob, optimism = calibration_bootstrap - calibration_validation)
    })
})

stait_calibration_optimism <- stait_calibration_bootstrap %>% group_by(model, predicted_prob) %>% summarise(mean_optimism = mean(optimism, 
    na.rm = TRUE), n_contributing = sum(!is.na(optimism)), .groups = "drop")

stait_calibration_corrected <- stait_apparent_calibration_curves %>% left_join(stait_calibration_optimism, by = c("model", 
    "predicted_prob")) %>% mutate(corrected_observed_prob_raw = observed_prob - mean_optimism, corrected_observed_prob = pmin(pmax(corrected_observed_prob_raw, 
    0), 1), proportion_bootstraps_contributing = n_contributing/1000)

stait_calibration_plot_data <- stait_calibration_corrected %>% filter(proportion_bootstraps_contributing >= 0.8) %>% 
    mutate(model = factor(model, levels = c("baseline", "full"), labels = c("Baseline model", "Full model (+ PD)")))

stait_calibration_labels <- stait_final_results$model_performance %>% transmute(model = factor(model, levels = c("baseline", 
    "full"), labels = c("Baseline model", "Full model (+ PD)")), label = sprintf("Calibration intercept = %.2f\nCalibration slope = %.2f", 
    calibration_intercept_optimism_corrected, calibration_slope_optimism_corrected))

master_model_performance <- bind_rows(epds_final_results$model_performance %>% mutate(outcome = "EPDS"), stais_final_results$model_performance %>% 
    mutate(outcome = "STAI-S"), stait_final_results$model_performance %>% mutate(outcome = "STAI-T")) %>% select(outcome, 
    model, apparent_auc = auc, mean_optimism_auc, optimism_corrected_auc = auc_optimism_corrected, apparent_calibration_intercept = calibration_intercept, 
    mean_optimism_calibration_intercept, optimism_corrected_calibration_intercept = calibration_intercept_optimism_corrected, 
    apparent_calibration_slope = calibration_slope, mean_optimism_calibration_slope, optimism_corrected_calibration_slope = calibration_slope_optimism_corrected) %>% 
    mutate(outcome = factor(outcome, levels = c("EPDS", "STAI-S", "STAI-T")), model = factor(model, levels = c("baseline", 
        "full"))) %>% arrange(outcome, model)

master_pd_stability_long <- bind_rows(epds_pd_bootstrap_stability %>% mutate(outcome = "EPDS"), stais_pd_bootstrap_stability %>% 
    mutate(outcome = "STAI-S"), stait_pd_bootstrap_stability %>% mutate(outcome = "STAI-T")) %>% mutate(outcome = factor(outcome, 
    levels = c("EPDS", "STAI-S", "STAI-T")), timepoint = if_else(grepl("2$", variable), "Late pregnancy", "Early pregnancy"), 
    measure_code = sub("2$", "", variable), selection_frequency_percent = 100 * bootstrap_majority_frequency, direction_consistency_percent = 100 * 
        bootstrap_direction_consistency) %>% select(outcome, variable, measure_code, timepoint, mean_imputation_selection, 
    selection_frequency_percent, n_bootstraps_majority_selected, predominant_bootstrap_direction, direction_consistency_percent) %>% 
    arrange(measure_code, timepoint, outcome)

master_pd_stability_table <- master_pd_stability_long %>% mutate(pd_measure = recode(measure_code, PARDiff = "Dilation at Response", 
    RodDiff = "Net Pupil Change", PCUDiff = "Early Peak Pupil Change", PAADiff = "Late Sustained Dilation", PTRDiff = "Early Mean Pupil Size")) %>% 
    select(pd_measure, timepoint, outcome, selection_frequency_percent, predominant_bootstrap_direction, direction_consistency_percent) %>% 
    pivot_wider(names_from = outcome, values_from = c(selection_frequency_percent, predominant_bootstrap_direction, 
        direction_consistency_percent), names_glue = "{outcome}_{.value}") %>% arrange(factor(pd_measure, levels = c("Early Mean Pupil Size", 
    "Early Peak Pupil Change", "Late Sustained Dilation", "Dilation at Response", "Net Pupil Change")), factor(timepoint, 
    levels = c("Early pregnancy", "Late pregnancy")))

master_calibration_plot_data <- bind_rows(epds_calibration_plot_data %>% mutate(outcome = "EPDS"), stais_calibration_plot_data %>% 
    mutate(outcome = "STAI-S"), stait_calibration_plot_data %>% mutate(outcome = "STAI-T")) %>% mutate(outcome = factor(outcome, 
    levels = c("EPDS", "STAI-S", "STAI-T")), model = factor(model, levels = c("Baseline model", "Full model (+ PD)")))

master_calibration_labels <- bind_rows(epds_calibration_labels %>% mutate(outcome = "EPDS"), stais_calibration_labels %>% 
    mutate(outcome = "STAI-S"), stait_calibration_labels %>% mutate(outcome = "STAI-T")) %>% mutate(outcome = factor(outcome, 
    levels = c("EPDS", "STAI-S", "STAI-T")), model = factor(model, levels = c("Baseline model", "Full model (+ PD)")))

master_calibration_plot <- ggplot(master_calibration_plot_data, aes(x = predicted_prob)) + geom_abline(intercept = 0, 
    slope = 1, linetype = "dotted", linewidth = 0.7) + geom_line(aes(y = observed_prob, linetype = "Apparent", 
    color = "Apparent"), linewidth = 0.9) + geom_line(aes(y = corrected_observed_prob, linetype = "Optimism-corrected", 
    color = "Optimism-corrected"), linewidth = 1.1) + facet_grid(outcome ~ model) + geom_text(data = master_calibration_labels, 
    aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.05, vjust = 1.2, size = 3.2) + coord_equal() + 
    scale_linetype_manual(name = NULL, values = c(Apparent = "dashed", `Optimism-corrected` = "solid")) + scale_color_manual(name = NULL, 
    values = c(Apparent = "#D55E00", `Optimism-corrected` = "#0072B2")) + labs(x = "Predicted probability", y = "Observed probability", 
    caption = "Curves are displayed where ≥80% of bootstrap resamples contributed to the optimism estimate. Dotted line indicates perfect calibration.") + 
    theme_classic() + theme(legend.position = "bottom", strip.text.x = element_text(face = "bold"), strip.text.y = element_text(face = "bold"), 
    strip.background = element_rect(fill = "white", colour = "black"))

combine_incremental_results <- function(outcome_label, final_results, variability) {
    bind_cols(tibble(outcome = outcome_label), final_results$incremental_value, variability)
}

master_incremental_value <- bind_rows(combine_incremental_results("EPDS", epds_final_results, epds_incremental_variability), 
    combine_incremental_results("STAI-S", stais_final_results, stais_incremental_variability), combine_incremental_results("STAI-T", 
        stait_final_results, stait_incremental_variability))

