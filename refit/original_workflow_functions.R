background_vars <- c("gb1_age", "country", "bmi", "edu", "first_preg", "workFT", "pcos", "exercise1", "exercise2", 
    "SSRI", "pos_LE1", "pos_LE2", "neg_LE1", "neg_LE2")

psychological_vars <- c("gb1_epds_sum", "gb2_epds_sum", "gb1_stai_s_sum", "gb2_stai_s_sum", "gb1_stai_t_sum", "gb2_stai_t_sum", 
    "gb1_ders16_sum", "gb2_ders16_sum", "gb1_saaq_sum", "gb2_saaq_sum", "gb1_rs14_sum", "gb2_rs14_sum", "gb1_pss_sum", 
    "gb2_pss_sum", "gb1_soc13_sum", "gb2_soc13_sum", "gb1_gad7_sum", "gb2_gad7_sum")

pd_vars <- c("PARDiff", "PARDiff2", "RodDiff", "RodDiff2", "PCUDiff", "PCUDiff2", "PAADiff", "PAADiff2", "PTRDiff", 
    "PTRDiff2")

baseline_vars <- c(background_vars, psychological_vars)

full_vars <- c(baseline_vars, pd_vars)

binary_vars <- c("country", "edu", "first_preg", "workFT", "pcos", "exercise1", "exercise2", "SSRI", "pos_LE1", 
    "pos_LE2", "neg_LE1", "neg_LE2")

continuous_vars <- c("gb1_age", "bmi", psychological_vars, pd_vars)

fill_single_level_missing <- function(data, predictor_vars) {
    fill_log <- tibble(variable = character(), n_missing_filled = integer(), fill_value = character())
    for (v in predictor_vars) {
        x <- data[[v]]
        n_missing <- sum(is.na(x))
        if (n_missing == 0) {
            next
        }
        observed <- x[!is.na(x)]
        n_unique_observed <- dplyr::n_distinct(observed)
        if (n_unique_observed == 0) {
            stop(paste0("Predictor '", v, "' is completely missing in this bootstrap sample."))
        }
        if (n_unique_observed == 1) {
            fill_value <- observed[1]
            x[is.na(x)] <- fill_value
            data[[v]] <- x
            fill_log <- bind_rows(fill_log, tibble(variable = v, n_missing_filled = n_missing, fill_value = as.character(fill_value)))
        }
    }
    list(data = data, fill_log = fill_log)
}

make_mice_spec_development <- function(data, outcome_var, predictor_vars = full_vars) {
    imp_data <- data %>% select(all_of(predictor_vars), all_of(outcome_var))
    single_level_fix <- fill_single_level_missing(data = imp_data, predictor_vars = predictor_vars)
    imp_data <- single_level_fix$data
    single_level_fill_log <- single_level_fix$fill_log
    ini <- mice::mice(imp_data, maxit = 0, printFlag = FALSE)
    method <- ini$method
    predictor_matrix <- ini$predictorMatrix
    method[] <- ""
    predictor_matrix[, ] <- 0
    continuous_present <- intersect(continuous_vars, predictor_vars)
    continuous_missing <- continuous_present[map_lgl(imp_data[continuous_present], ~any(is.na(.x)))]
    method[continuous_missing] <- "pmm"
    binary_present <- intersect(binary_vars, predictor_vars)
    binary_missing <- binary_present[map_lgl(imp_data[binary_present], ~any(is.na(.x)))]
    method[binary_missing] <- "cart"
    method[outcome_var] <- ""
    core_vars <- c("gb1_age", "bmi", "gb1_epds_sum", "gb2_epds_sum", "gb1_stai_s_sum", "gb2_stai_s_sum", "gb1_stai_t_sum", 
        "gb2_stai_t_sum", "gb1_gad7_sum", "gb2_gad7_sum", "gb1_pss_sum", "gb2_pss_sum", outcome_var)
    core_vars <- intersect(core_vars, names(imp_data))
    vars_to_impute <- names(method)[method != ""]
    predictor_matrix[vars_to_impute, core_vars] <- 1
    paired_vars <- list(c("exercise1", "exercise2"), c("pos_LE1", "pos_LE2"), c("neg_LE1", "neg_LE2"), c("gb1_epds_sum", 
        "gb2_epds_sum"), c("gb1_stai_s_sum", "gb2_stai_s_sum"), c("gb1_stai_t_sum", "gb2_stai_t_sum"), c("gb1_ders16_sum", 
        "gb2_ders16_sum"), c("gb1_saaq_sum", "gb2_saaq_sum"), c("gb1_rs14_sum", "gb2_rs14_sum"), c("gb1_pss_sum", 
        "gb2_pss_sum"), c("gb1_soc13_sum", "gb2_soc13_sum"), c("gb1_gad7_sum", "gb2_gad7_sum"), c("PARDiff", "PARDiff2"), 
        c("RodDiff", "RodDiff2"), c("PCUDiff", "PCUDiff2"), c("PAADiff", "PAADiff2"), c("PTRDiff", "PTRDiff2"))
    for (pair in paired_vars) {
        a <- pair[1]
        b <- pair[2]
        if (a %in% vars_to_impute && b %in% names(imp_data)) {
            predictor_matrix[a, b] <- 1
        }
        if (b %in% vars_to_impute && a %in% names(imp_data)) {
            predictor_matrix[b, a] <- 1
        }
    }
    pd_present <- intersect(pd_vars, names(imp_data))
    for (v in intersect(vars_to_impute, pd_present)) {
        other_pd <- setdiff(pd_present, v)
        predictor_matrix[v, other_pd] <- 1
    }
    diag(predictor_matrix) <- 0
    predictor_matrix[outcome_var, ] <- 0
    list(data = imp_data, method = method, predictor_matrix = predictor_matrix, single_level_fill_log = single_level_fill_log)
}

make_mice_spec_validation <- function(data, outcome_var, predictor_vars = full_vars) {
    imp_data <- data %>% select(all_of(predictor_vars), all_of(outcome_var))
    single_level_fix <- fill_single_level_missing(data = imp_data, predictor_vars = predictor_vars)
    imp_data <- single_level_fix$data
    single_level_fill_log <- single_level_fix$fill_log
    ini <- mice::mice(imp_data, maxit = 0, printFlag = FALSE)
    method <- ini$method
    predictor_matrix <- ini$predictorMatrix
    method[] <- ""
    predictor_matrix[, ] <- 0
    continuous_present <- intersect(continuous_vars, predictor_vars)
    continuous_missing <- continuous_present[map_lgl(imp_data[continuous_present], ~any(is.na(.x)))]
    method[continuous_missing] <- "pmm"
    binary_present <- intersect(binary_vars, predictor_vars)
    binary_missing <- binary_present[map_lgl(imp_data[binary_present], ~any(is.na(.x)))]
    method[binary_missing] <- "cart"
    method[outcome_var] <- ""
    core_vars <- c("gb1_age", "bmi", "gb1_epds_sum", "gb2_epds_sum", "gb1_stai_s_sum", "gb2_stai_s_sum", "gb1_stai_t_sum", 
        "gb2_stai_t_sum", "gb1_gad7_sum", "gb2_gad7_sum", "gb1_pss_sum", "gb2_pss_sum", outcome_var)
    core_vars <- intersect(core_vars, names(imp_data))
    vars_to_impute <- names(method)[method != ""]
    predictor_matrix[vars_to_impute, core_vars] <- 1
    paired_vars <- list(c("exercise1", "exercise2"), c("pos_LE1", "pos_LE2"), c("neg_LE1", "neg_LE2"), c("gb1_epds_sum", 
        "gb2_epds_sum"), c("gb1_stai_s_sum", "gb2_stai_s_sum"), c("gb1_stai_t_sum", "gb2_stai_t_sum"), c("gb1_ders16_sum", 
        "gb2_ders16_sum"), c("gb1_saaq_sum", "gb2_saaq_sum"), c("gb1_rs14_sum", "gb2_rs14_sum"), c("gb1_pss_sum", 
        "gb2_pss_sum"), c("gb1_soc13_sum", "gb2_soc13_sum"), c("gb1_gad7_sum", "gb2_gad7_sum"), c("PARDiff", "PARDiff2"), 
        c("RodDiff", "RodDiff2"), c("PCUDiff", "PCUDiff2"), c("PAADiff", "PAADiff2"), c("PTRDiff", "PTRDiff2"))
    for (pair in paired_vars) {
        a <- pair[1]
        b <- pair[2]
        if (a %in% vars_to_impute && b %in% names(imp_data)) {
            predictor_matrix[a, b] <- 1
        }
        if (b %in% vars_to_impute && a %in% names(imp_data)) {
            predictor_matrix[b, a] <- 1
        }
    }
    pd_present <- intersect(pd_vars, names(imp_data))
    for (v in intersect(vars_to_impute, pd_present)) {
        other_pd <- setdiff(pd_present, v)
        predictor_matrix[v, other_pd] <- 1
    }
    diag(predictor_matrix) <- 0
    predictor_matrix[outcome_var, ] <- 0
    predictor_matrix[, outcome_var] <- 0
    list(data = imp_data, method = method, predictor_matrix = predictor_matrix, single_level_fill_log = single_level_fill_log)
}

make_stratified_group_folds <- function(data, outcome_var, k = 10, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    unique_people <- data %>% transmute(pid = pid, outcome = .data[[outcome_var]]) %>% distinct()
    fold_map <- unique_people %>% group_by(outcome) %>% mutate(fold = sample(rep(seq_len(k), length.out = n()))) %>% 
        ungroup()
    foldid <- fold_map$fold[match(data$pid, fold_map$pid)]
    pid_fold_check <- tibble(pid = data$pid, fold = foldid) %>% group_by(pid) %>% summarise(n_folds = n_distinct(fold), 
        .groups = "drop")
    stopifnot(all(pid_fold_check$n_folds == 1))
    foldid
}

make_model_matrices <- function(train_data, validation_data, predictor_vars) {
    x_train <- model.matrix(~., data = train_data %>% select(all_of(predictor_vars)))
    x_validation <- model.matrix(~., data = validation_data %>% select(all_of(predictor_vars)))
    x_train <- x_train[, colnames(x_train) != "(Intercept)", drop = FALSE]
    x_validation <- x_validation[, colnames(x_validation) != "(Intercept)", drop = FALSE]
    all_columns <- union(colnames(x_train), colnames(x_validation))
    align_matrix <- function(x, columns) {
        out <- matrix(0, nrow = nrow(x), ncol = length(columns), dimnames = list(NULL, columns))
        out[, colnames(x)] <- x
        out
    }
    list(train = align_matrix(x_train, all_columns), validation = align_matrix(x_validation, all_columns))
}

fit_elastic_net <- function(train_data, validation_data, outcome_var, predictor_vars, alpha_grid = c(0.1, 0.25, 
    0.5, 0.75, 1), k = 10, seed = 1) {
    matrices <- make_model_matrices(train_data = train_data, validation_data = validation_data, predictor_vars = predictor_vars)
    x_train <- matrices$train
    x_validation <- matrices$validation
    y_train <- train_data[[outcome_var]]
    stopifnot(!anyNA(x_train))
    stopifnot(!anyNA(x_validation))
    stopifnot(!anyNA(y_train))
    foldid <- make_stratified_group_folds(data = train_data, outcome_var = outcome_var, k = k, seed = seed)
    cv_fold_diagnostics <- map_dfr(sort(unique(foldid)), function(f) {
        training_rows <- foldid != f
        tibble(fold = f, n_train = sum(training_rows), n_events_train = sum(y_train[training_rows] == 1), n_non_events_train = sum(y_train[training_rows] == 
            0), n_unique_events_train = n_distinct(train_data$pid[training_rows & y_train == 1]), n_unique_non_events_train = n_distinct(train_data$pid[training_rows & 
            y_train == 0]))
    })
    glmnet_warnings <- character()
    capture_glmnet_warnings <- function(expr) {
        withCallingHandlers(expr, warning = function(w) {
            glmnet_warnings <<- c(glmnet_warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    }
    cv_models <- map(alpha_grid, function(a) {
        capture_glmnet_warnings(glmnet::cv.glmnet(x = x_train, y = y_train, family = "binomial", alpha = a, foldid = foldid, 
            type.measure = "deviance", grouped = FALSE, standardize = TRUE))
    })
    cv_deviance <- map_dbl(cv_models, ~min(.x$cvm, na.rm = TRUE))
    best_index <- which.min(cv_deviance)
    best_alpha <- alpha_grid[best_index]
    best_cv <- cv_models[[best_index]]
    best_lambda <- best_cv$lambda.min
    final_fit <- capture_glmnet_warnings(glmnet::glmnet(x = x_train, y = y_train, family = "binomial", alpha = best_alpha, 
        lambda = best_lambda, standardize = TRUE))
    train_prediction <- as.numeric(predict(final_fit, newx = x_train, type = "response"))
    validation_prediction <- as.numeric(predict(final_fit, newx = x_validation, type = "response"))
    coef_matrix <- as.matrix(coef(final_fit))
    coef_table <- tibble(variable = rownames(coef_matrix), coefficient = coef_matrix[, 1]) %>% filter(variable != 
        "(Intercept)") %>% mutate(selected = coefficient != 0, direction = case_when(coefficient > 0 ~ "positive", 
        coefficient < 0 ~ "negative", TRUE ~ "zero"))
    list(alpha = best_alpha, lambda = best_lambda, cv_deviance = cv_deviance[best_index], train_prediction = train_prediction, 
        validation_prediction = validation_prediction, coefficients = coef_table, cv_fold_diagnostics = cv_fold_diagnostics, 
        glmnet_warning_n = length(glmnet_warnings), glmnet_warning_text = paste(unique(glmnet_warnings), collapse = " | "))
}

calculate_calibration <- function(outcome, predicted_prob) {
    eps <- 1e-06
    p <- pmin(pmax(as.numeric(predicted_prob), eps), 1 - eps)
    lp <- qlogis(p)
    intercept_fit <- glm(outcome ~ offset(lp), family = binomial(), method = brglm2::brglmFit, type = "AS_mean")
    calibration_intercept <- unname(coef(intercept_fit)[1])
    lp_range <- diff(range(lp))
    slope_estimable <- is.finite(lp_range) && lp_range > 1e-12
    if (slope_estimable) {
        slope_fit <- glm(outcome ~ lp, family = binomial(), method = brglm2::brglmFit, type = "AS_mean")
        calibration_slope <- unname(coef(slope_fit)["lp"])
    }
    else {
        calibration_slope <- NA_real_
    }
    tibble(calibration_intercept = calibration_intercept, calibration_slope = calibration_slope, calibration_slope_estimable = slope_estimable)
}

run_one_bootstrap <- function(data, outcome_var, predictor_vars, m = 3, bootstrap_id = 1, seed = 26092026) {
    set.seed(seed + bootstrap_id)
    repeat {
        bootstrap_index <- sample.int(n = nrow(data), size = nrow(data), replace = TRUE)
        boot_data <- data[bootstrap_index, , drop = FALSE]
        if (length(unique(boot_data[[outcome_var]])) == 2) {
            break
        }
    }
    dev_spec <- make_mice_spec_development(data = boot_data, outcome_var = outcome_var, predictor_vars = predictor_vars)
    dev_imp <- mice::mice(data = dev_spec$data, m = m, maxit = 20, method = dev_spec$method, predictorMatrix = dev_spec$predictor_matrix, 
        printFlag = FALSE, seed = seed + bootstrap_id * 100)
    development_for_validation <- boot_data
    validation_for_imputation <- data
    validation_for_imputation[[outcome_var]] <- NA_integer_
    combined_for_validation <- bind_rows(development_for_validation, validation_for_imputation)
    validation_ignore <- c(rep(FALSE, nrow(development_for_validation)), rep(TRUE, nrow(validation_for_imputation)))
    val_spec <- make_mice_spec_validation(data = combined_for_validation, outcome_var = outcome_var, predictor_vars = predictor_vars)
    val_imp <- mice::mice(data = val_spec$data, m = m, maxit = 20, method = val_spec$method, predictorMatrix = val_spec$predictor_matrix, 
        ignore = validation_ignore, printFlag = FALSE, seed = seed + bootstrap_id * 100 + 50)
    train_pred_matrix <- matrix(NA_real_, nrow = nrow(boot_data), ncol = m)
    validation_pred_matrix <- matrix(NA_real_, nrow = nrow(data), ncol = m)
    tuning_list <- vector("list", m)
    coefficient_list <- vector("list", m)
    for (j in seq_len(m)) {
        train_j <- mice::complete(dev_imp, j)
        train_j$pid <- boot_data$pid
        validation_combined_j <- mice::complete(val_imp, j)
        first_validation_row <- nrow(boot_data) + 1
        validation_j <- validation_combined_j[first_validation_row:nrow(validation_combined_j), , drop = FALSE]
        validation_j$pid <- data$pid
        validation_j[[outcome_var]] <- data[[outcome_var]]
        model_j <- fit_elastic_net(train_data = train_j, validation_data = validation_j, outcome_var = outcome_var, 
            predictor_vars = predictor_vars, k = 5, seed = seed + bootstrap_id * 1000)
        train_pred_matrix[, j] <- model_j$train_prediction
        validation_pred_matrix[, j] <- model_j$validation_prediction
        tuning_list[[j]] <- tibble(bootstrap = bootstrap_id, imputation = j, alpha = model_j$alpha, lambda = model_j$lambda, 
            cv_deviance = model_j$cv_deviance, n_selected = sum(model_j$coefficients$selected), min_cv_events = min(model_j$cv_fold_diagnostics$n_events_train), 
            min_cv_unique_events = min(model_j$cv_fold_diagnostics$n_unique_events_train), glmnet_warning_n = model_j$glmnet_warning_n, 
            glmnet_warning_text = model_j$glmnet_warning_text)
        coefficient_list[[j]] <- model_j$coefficients %>% mutate(bootstrap = bootstrap_id, imputation = j)
    }
    mean_train_prediction <- rowMeans(train_pred_matrix)
    mean_validation_prediction <- rowMeans(validation_pred_matrix)
    stopifnot(!anyNA(mean_train_prediction), !anyNA(mean_validation_prediction))
    auc_bootstrap <- as.numeric(pROC::auc(response = boot_data[[outcome_var]], predictor = mean_train_prediction, 
        quiet = TRUE))
    auc_validation <- as.numeric(pROC::auc(response = data[[outcome_var]], predictor = mean_validation_prediction, 
        quiet = TRUE))
    calibration_bootstrap <- calculate_calibration(outcome = boot_data[[outcome_var]], predicted_prob = mean_train_prediction)
    calibration_validation <- calculate_calibration(outcome = data[[outcome_var]], predicted_prob = mean_validation_prediction)
    performance <- tibble(bootstrap = bootstrap_id, auc_bootstrap = auc_bootstrap, auc_validation = auc_validation, 
        optimism_auc = auc_bootstrap - auc_validation, calibration_intercept_bootstrap = calibration_bootstrap$calibration_intercept, 
        calibration_intercept_validation = calibration_validation$calibration_intercept, optimism_calibration_intercept = calibration_bootstrap$calibration_intercept - 
            calibration_validation$calibration_intercept, calibration_slope_bootstrap = calibration_bootstrap$calibration_slope, 
        calibration_slope_validation = calibration_validation$calibration_slope, optimism_calibration_slope = calibration_bootstrap$calibration_slope - 
            calibration_validation$calibration_slope)
    diagnostics <- tibble(bootstrap = bootstrap_id, n_bootstrap_rows = nrow(boot_data), n_unique_participants = n_distinct(boot_data$pid), 
        n_events = sum(boot_data[[outcome_var]] == 1), n_non_events = sum(boot_data[[outcome_var]] == 0), n_unique_events = boot_data %>% 
            distinct(pid, .data[[outcome_var]]) %>% summarise(n = sum(.data[[outcome_var]] == 1)) %>% pull(n), 
        dev_mice_logged_events = ifelse(is.null(dev_imp$loggedEvents), 0, nrow(dev_imp$loggedEvents)), val_mice_logged_events = ifelse(is.null(val_imp$loggedEvents), 
            0, nrow(val_imp$loggedEvents)))
    list(performance = performance, tuning = bind_rows(tuning_list), coefficients = bind_rows(coefficient_list), 
        diagnostics = diagnostics, predictions = list(bootstrap = tibble(pid = boot_data$pid, outcome = boot_data[[outcome_var]], 
            predicted_prob = mean_train_prediction), validation = tibble(pid = data$pid, outcome = data[[outcome_var]], 
            predicted_prob = mean_validation_prediction)), bootstrap_index = bootstrap_index)
}

calculate_idi <- function(outcome, predicted_baseline, predicted_full) {
    stopifnot(length(outcome) == length(predicted_baseline), length(outcome) == length(predicted_full))
    baseline_event_mean <- mean(predicted_baseline[outcome == 1])
    baseline_nonevent_mean <- mean(predicted_baseline[outcome == 0])
    full_event_mean <- mean(predicted_full[outcome == 1])
    full_nonevent_mean <- mean(predicted_full[outcome == 0])
    discrimination_slope_baseline <- baseline_event_mean - baseline_nonevent_mean
    discrimination_slope_full <- full_event_mean - full_nonevent_mean
    idi <- discrimination_slope_full - discrimination_slope_baseline
    tibble(discrimination_slope_baseline = discrimination_slope_baseline, discrimination_slope_full = discrimination_slope_full, 
        idi = idi)
}

run_one_paired_bootstrap <- function(data, outcome_var, baseline_vars, full_vars, m = 3, bootstrap_id = 1, seed = 26092026) {
    baseline_result <- run_one_bootstrap(data = data, outcome_var = outcome_var, predictor_vars = baseline_vars, 
        m = m, bootstrap_id = bootstrap_id, seed = seed)
    full_result <- run_one_bootstrap(data = data, outcome_var = outcome_var, predictor_vars = full_vars, m = m, 
        bootstrap_id = bootstrap_id, seed = seed)
    stopifnot(identical(baseline_result$bootstrap_index, full_result$bootstrap_index))
    stopifnot(identical(baseline_result$predictions$bootstrap$pid, full_result$predictions$bootstrap$pid))
    stopifnot(identical(baseline_result$predictions$bootstrap$outcome, full_result$predictions$bootstrap$outcome))
    baseline_auc_bootstrap <- baseline_result$performance$auc_bootstrap
    full_auc_bootstrap <- full_result$performance$auc_bootstrap
    baseline_auc_validation <- baseline_result$performance$auc_validation
    full_auc_validation <- full_result$performance$auc_validation
    delta_auc_bootstrap <- full_auc_bootstrap - baseline_auc_bootstrap
    delta_auc_validation <- full_auc_validation - baseline_auc_validation
    optimism_delta_auc <- delta_auc_bootstrap - delta_auc_validation
    idi_bootstrap <- calculate_idi(outcome = baseline_result$predictions$bootstrap$outcome, predicted_baseline = baseline_result$predictions$bootstrap$predicted_prob, 
        predicted_full = full_result$predictions$bootstrap$predicted_prob)
    stopifnot(identical(baseline_result$predictions$validation$pid, full_result$predictions$validation$pid))
    idi_validation <- calculate_idi(outcome = baseline_result$predictions$validation$outcome, predicted_baseline = baseline_result$predictions$validation$predicted_prob, 
        predicted_full = full_result$predictions$validation$predicted_prob)
    optimism_idi <- idi_bootstrap$idi - idi_validation$idi
    comparison <- tibble(bootstrap = bootstrap_id, baseline_auc_bootstrap = baseline_auc_bootstrap, full_auc_bootstrap = full_auc_bootstrap, 
        delta_auc_bootstrap = delta_auc_bootstrap, baseline_auc_validation = baseline_auc_validation, full_auc_validation = full_auc_validation, 
        delta_auc_validation = delta_auc_validation, optimism_delta_auc = optimism_delta_auc, baseline_discrimination_slope_bootstrap = idi_bootstrap$discrimination_slope_baseline, 
        full_discrimination_slope_bootstrap = idi_bootstrap$discrimination_slope_full, idi_bootstrap = idi_bootstrap$idi, 
        baseline_discrimination_slope_validation = idi_validation$discrimination_slope_baseline, full_discrimination_slope_validation = idi_validation$discrimination_slope_full, 
        idi_validation = idi_validation$idi, optimism_idi = optimism_idi)
    list(comparison = comparison, performance = bind_rows(baseline_result$performance %>% mutate(model = "baseline"), 
        full_result$performance %>% mutate(model = "full")), tuning = bind_rows(baseline_result$tuning %>% mutate(model = "baseline"), 
        full_result$tuning %>% mutate(model = "full")), coefficients = bind_rows(baseline_result$coefficients %>% 
        mutate(model = "baseline"), full_result$coefficients %>% mutate(model = "full")), diagnostics = bind_rows(baseline_result$diagnostics %>% 
        mutate(model = "baseline"), full_result$diagnostics %>% mutate(model = "full")), predictions = list(baseline = baseline_result$predictions, 
        full = full_result$predictions), bootstrap_index = baseline_result$bootstrap_index)
}

fit_apparent_model <- function(data, outcome_var, predictor_vars, m = 10, seed = 26092026) {
    apparent_spec <- make_mice_spec_development(data = data, outcome_var = outcome_var, predictor_vars = predictor_vars)
    apparent_imp <- mice::mice(data = apparent_spec$data, m = m, maxit = 20, method = apparent_spec$method, predictorMatrix = apparent_spec$predictor_matrix, 
        printFlag = FALSE, seed = seed)
    prediction_matrix <- matrix(NA_real_, nrow = nrow(data), ncol = m)
    tuning_list <- vector("list", m)
    coefficient_list <- vector("list", m)
    for (j in seq_len(m)) {
        data_j <- mice::complete(apparent_imp, j)
        data_j$pid <- data$pid
        model_j <- fit_elastic_net(train_data = data_j, validation_data = data_j, outcome_var = outcome_var, predictor_vars = predictor_vars, 
            k = 10, seed = seed)
        prediction_matrix[, j] <- model_j$train_prediction
        tuning_list[[j]] <- tibble(imputation = j, alpha = model_j$alpha, lambda = model_j$lambda, n_selected = sum(model_j$coefficients$selected), 
            glmnet_warning_n = model_j$glmnet_warning_n)
        coefficient_list[[j]] <- model_j$coefficients %>% mutate(imputation = j)
    }
    mean_prediction <- rowMeans(prediction_matrix)
    auc_apparent <- as.numeric(pROC::auc(response = data[[outcome_var]], predictor = mean_prediction, quiet = TRUE))
    calibration_apparent <- calculate_calibration(outcome = data[[outcome_var]], predicted_prob = mean_prediction)
    list(performance = tibble(auc = auc_apparent, calibration_intercept = calibration_apparent$calibration_intercept, 
        calibration_slope = calibration_apparent$calibration_slope), predictions = tibble(pid = data$pid, outcome = data[[outcome_var]], 
        predicted_prob = mean_prediction), tuning = bind_rows(tuning_list), coefficients = bind_rows(coefficient_list), 
        mice_logged_events = ifelse(is.null(apparent_imp$loggedEvents), 0, nrow(apparent_imp$loggedEvents)))
}

calculate_optimism_corrected <- function(apparent_baseline, apparent_full, apparent_idi, paired_bootstrap) {
    model_optimism <- paired_bootstrap$performance %>% group_by(model) %>% summarise(mean_optimism_auc = mean(optimism_auc, 
        na.rm = TRUE), mean_optimism_calibration_intercept = mean(optimism_calibration_intercept, na.rm = TRUE), 
        mean_optimism_calibration_slope = mean(optimism_calibration_slope, na.rm = TRUE), .groups = "drop")
    apparent_models <- bind_rows(apparent_baseline$performance %>% mutate(model = "baseline"), apparent_full$performance %>% 
        mutate(model = "full"))
    corrected_models <- apparent_models %>% left_join(model_optimism, by = "model") %>% mutate(auc_optimism_corrected = auc - 
        mean_optimism_auc, calibration_intercept_optimism_corrected = calibration_intercept - mean_optimism_calibration_intercept, 
        calibration_slope_optimism_corrected = calibration_slope - mean_optimism_calibration_slope)
    apparent_delta_auc <- apparent_full$performance$auc - apparent_baseline$performance$auc
    apparent_idi_value <- apparent_idi$idi
    mean_optimism_delta_auc <- mean(paired_bootstrap$comparison$optimism_delta_auc, na.rm = TRUE)
    mean_optimism_idi <- mean(paired_bootstrap$comparison$optimism_idi, na.rm = TRUE)
    corrected_comparison <- tibble(apparent_delta_auc = apparent_delta_auc, mean_optimism_delta_auc = mean_optimism_delta_auc, 
        optimism_corrected_delta_auc = apparent_delta_auc - mean_optimism_delta_auc, apparent_idi = apparent_idi_value, 
        mean_optimism_idi = mean_optimism_idi, optimism_corrected_idi = apparent_idi_value - mean_optimism_idi)
    list(model_performance = corrected_models, incremental_value = corrected_comparison)
}

estimate_calibration_curve <- function(outcome, predicted_prob, grid, f = 2/3) {
    keep <- !is.na(outcome) & is.finite(predicted_prob)
    x <- as.numeric(predicted_prob[keep])
    y <- as.numeric(outcome[keep])
    if (dplyr::n_distinct(x) < 2) {
        return(tibble(predicted_prob = grid, observed_prob = NA_real_))
    }
    smooth_fit <- stats::lowess(x = x, y = y, f = f, iter = 0)
    observed_grid <- stats::approx(x = smooth_fit$x, y = smooth_fit$y, xout = grid, rule = 1, ties = mean)$y
    tibble(predicted_prob = grid, observed_prob = pmin(pmax(observed_grid, 0), 1))
}

