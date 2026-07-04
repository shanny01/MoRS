mors_read_table <- function(path) {
  data.frame(data.table::fread(path, header = TRUE, data.table = FALSE))
}

mors_stop_if_missing_cols <- function(df, cols, df_name) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "Missing columns in %s: %s",
        df_name,
        paste(missing_cols, collapse = ", ")
      )
    )
  }
}

mors_residualize_outcome <- function(df, outcome_col, covariate_cols) {
  form <- stats::as.formula(
    paste(outcome_col, "~", paste(covariate_cols, collapse = " + "))
  )
  stats::resid(stats::lm(form, data = df))
}

mors_conditional_r2 <- function(y_resid, score) {
  keep <- stats::complete.cases(y_resid, score)
  if (sum(keep) < 3) {
    return(NA_real_)
  }
  
  corr <- stats::cor(y_resid[keep], score[keep])
  if (is.na(corr)) {
    return(NA_real_)
  }
  
  corr^2
}

mors_fit_lasso <- function(df, feature_cols, y_resid, family, nfolds, seed) {
  keep <- stats::complete.cases(df[, feature_cols, drop = FALSE], y_resid)
  x <- as.matrix(df[keep, feature_cols, drop = FALSE])
  y <- y_resid[keep]
  
  if (ncol(x) == 0) {
    stop("No feature columns supplied to lasso model.")
  }
  
  set.seed(seed)
  glmnet::cv.glmnet(
    x = x,
    y = y,
    family = family,
    alpha = 1,
    nfolds = nfolds,
    standardize = TRUE
  )
}

mors_predict_lasso_score <- function(model, df, feature_cols) {
  x <- as.matrix(df[, feature_cols, drop = FALSE])
  as.numeric(stats::predict(model, newx = x, s = "lambda.min"))
}

mors_tune_linear_model <- function(y_resid, prs_df, extra_scores = NULL, prs_cols) {
  best <- list(r2 = -Inf, prs = NA_character_, fit = NULL)
  
  for (prs_col in prs_cols) {
    design <- data.frame(prs = prs_df[[prs_col]])
    names(design)[1] <- prs_col
    
    if (!is.null(extra_scores)) {
      design <- cbind(design, extra_scores)
    }
    
    dat <- data.frame(y = y_resid, design)
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]
    if (nrow(dat) < 3) {
      next
    }
    
    fit <- stats::lm(y ~ ., data = dat)
    pred <- stats::predict(fit, newdata = dat)
    r2 <- mors_conditional_r2(dat$y, pred)
    
    if (!is.na(r2) && r2 > best$r2) {
      best <- list(r2 = r2, prs = prs_col, fit = fit)
    }
  }
  
  if (!is.finite(best$r2)) {
    stop("No valid PRS model could be tuned. Check missingness and column names.")
  }
  
  best
}

mors_tune_step_model <- function(y_resid, prs_df, omics_scores, prs_cols) {
  best <- list(r2 = -Inf, prs = NA_character_, fit = NULL)
  
  for (prs_col in prs_cols) {
    dat <- data.frame(y = y_resid, prs = prs_df[[prs_col]], omics_scores)
    names(dat)[2] <- prs_col
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]
    if (nrow(dat) < 3) {
      next
    }
    
    lower_form <- stats::as.formula(paste("y ~", prs_col))
    upper_form <- stats::as.formula(
      paste("y ~", paste(c(prs_col, names(omics_scores)), collapse = " + "))
    )
    
    lower_fit <- stats::lm(lower_form, data = dat)
    step_fit <- stats::step(
      lower_fit,
      scope = list(lower = lower_form, upper = upper_form),
      direction = "forward",
      trace = 0
    )
    
    pred <- stats::predict(step_fit, newdata = dat)
    r2 <- mors_conditional_r2(dat$y, pred)
    
    if (!is.na(r2) && r2 > best$r2) {
      best <- list(r2 = r2, prs = prs_col, fit = step_fit)
    }
  }
  
  if (!is.finite(best$r2)) {
    stop("No valid stepwise model could be tuned. Check missingness and input columns.")
  }
  
  best
}

mors_bootstrap_r2 <- function(y_resid, score_list, B, seed) {
  model_names <- names(score_list)
  out <- matrix(NA_real_, nrow = B, ncol = length(model_names))
  colnames(out) <- model_names
  
  set.seed(seed)
  n <- length(y_resid)
  
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    y_b <- y_resid[idx]
    
    for (model_name in model_names) {
      out[b, model_name] <- mors_conditional_r2(y_b, score_list[[model_name]][idx])
    }
  }
  
  out
}

mors_bootstrap_ci <- function(boot_mat) {
  ci <- t(apply(boot_mat, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  data.frame(
    model = rownames(ci),
    lower = ci[, 1],
    upper = ci[, 2],
    row.names = NULL
  )
}

mors_bootstrap_pvalue <- function(diff_vec) {
  diff_vec <- diff_vec[is.finite(diff_vec)]
  if (length(diff_vec) < 2) {
    return(NA_real_)
  }
  
  mu <- mean(diff_vec)
  sigma <- stats::sd(diff_vec)
  
  if (is.na(sigma) || sigma == 0) {
    return(ifelse(mu > 0, 0, 1))
  }
  
  stats::pnorm(0, mean = mu, sd = sigma)
}

mors_prepare_inputs <- function(
    phenotype_train_file,
    phenotype_tuning_file,
    phenotype_test_file,
    gpf_train_file,
    gpf_tuning_file,
    gpf_test_file,
    prs_tuning_file,
    prs_test_file,
    outcome_col,
    covariate_cols,
    id_col,
    trs_cols,
    prors_cols,
    pmrs_cols,
    smrs_cols,
    n_prs_cols
) {
  train_pheno <- mors_read_table(phenotype_train_file)
  tuning_pheno <- mors_read_table(phenotype_tuning_file)
  test_pheno <- mors_read_table(phenotype_test_file)
  
  train_gpf <- mors_read_table(gpf_train_file)
  tuning_gpf <- mors_read_table(gpf_tuning_file)
  test_gpf <- mors_read_table(gpf_test_file)
  
  tuning_prs <- mors_read_table(prs_tuning_file)
  test_prs <- mors_read_table(prs_test_file)
  
  required_pheno_cols <- c(id_col, outcome_col, covariate_cols)
  mors_stop_if_missing_cols(train_pheno, required_pheno_cols, "phenotype_train_file")
  mors_stop_if_missing_cols(tuning_pheno, required_pheno_cols, "phenotype_tuning_file")
  mors_stop_if_missing_cols(test_pheno, required_pheno_cols, "phenotype_test_file")
  
  omics_groups <- list(
    trs = trs_cols,
    prors = prors_cols,
    pmrs = pmrs_cols,
    smrs = smrs_cols
  )
  all_omics_cols <- unique(unlist(omics_groups, use.names = FALSE))
  
  mors_stop_if_missing_cols(train_gpf, c(id_col, all_omics_cols), "gpf_train_file")
  mors_stop_if_missing_cols(tuning_gpf, c(id_col, all_omics_cols), "gpf_tuning_file")
  mors_stop_if_missing_cols(test_gpf, c(id_col, all_omics_cols), "gpf_test_file")
  
  train_dat <- merge(
    train_pheno,
    train_gpf[, c(id_col, all_omics_cols), drop = FALSE],
    by = id_col
  )
  tuning_dat <- merge(
    tuning_pheno,
    tuning_gpf[, c(id_col, all_omics_cols), drop = FALSE],
    by = id_col
  )
  test_dat <- merge(
    test_pheno,
    test_gpf[, c(id_col, all_omics_cols), drop = FALSE],
    by = id_col
  )
  
  tuning_prs <- merge(tuning_dat[, id_col, drop = FALSE], tuning_prs, by = id_col)
  test_prs <- merge(test_dat[, id_col, drop = FALSE], test_prs, by = id_col)
  
  prs_cols <- setdiff(names(tuning_prs), id_col)
  if (!is.null(n_prs_cols)) {
    prs_cols <- prs_cols[seq_len(min(n_prs_cols, length(prs_cols)))]
  }
  if (length(prs_cols) == 0) {
    stop("No PRS columns found after applying n_prs_cols.")
  }
  
  mors_stop_if_missing_cols(test_prs, c(id_col, prs_cols), "prs_test_file")
  
  list(
    train_dat = train_dat,
    tuning_dat = tuning_dat,
    test_dat = test_dat,
    tuning_prs = tuning_prs,
    test_prs = test_prs,
    prs_cols = prs_cols,
    omics_groups = omics_groups,
    all_omics_cols = all_omics_cols
  )
}

mors_fit_all_lasso_models <- function(
    train_dat,
    tuning_dat,
    test_dat,
    omics_groups,
    all_omics_cols,
    y_train_resid,
    family,
    nfolds,
    seed
) {
  lasso_models <- list()
  tuning_scores <- list()
  test_scores <- list()
  
  for (group_name in names(omics_groups)) {
    feature_cols <- omics_groups[[group_name]]
    fit <- mors_fit_lasso(train_dat, feature_cols, y_train_resid, family, nfolds, seed)
    
    lasso_models[[group_name]] <- fit
    tuning_scores[[group_name]] <- mors_predict_lasso_score(fit, tuning_dat, feature_cols)
    test_scores[[group_name]] <- mors_predict_lasso_score(fit, test_dat, feature_cols)
  }
  
  combined_fit <- mors_fit_lasso(train_dat, all_omics_cols, y_train_resid, family, nfolds, seed)
  lasso_models$combined <- combined_fit
  tuning_scores$lasso <- mors_predict_lasso_score(combined_fit, tuning_dat, all_omics_cols)
  test_scores$lasso <- mors_predict_lasso_score(combined_fit, test_dat, all_omics_cols)
  
  list(
    lasso_models = lasso_models,
    tuning_scores = tuning_scores,
    test_scores = test_scores
  )
}

mors_tune_models <- function(y_tuning_resid, tuning_prs, prs_cols, tuning_scores) {
  tuning_omics_df <- data.frame(
    trs = tuning_scores$trs,
    prors = tuning_scores$prors,
    pmrs = tuning_scores$pmrs,
    smrs = tuning_scores$smrs
  )
  
  list(
    base = mors_tune_linear_model(y_tuning_resid, tuning_prs, NULL, prs_cols),
    trs = mors_tune_linear_model(y_tuning_resid, tuning_prs, data.frame(trs = tuning_scores$trs), prs_cols),
    prors = mors_tune_linear_model(y_tuning_resid, tuning_prs, data.frame(prors = tuning_scores$prors), prs_cols),
    pmrs = mors_tune_linear_model(y_tuning_resid, tuning_prs, data.frame(pmrs = tuning_scores$pmrs), prs_cols),
    smrs = mors_tune_linear_model(y_tuning_resid, tuning_prs, data.frame(smrs = tuning_scores$smrs), prs_cols),
    lasso = mors_tune_linear_model(y_tuning_resid, tuning_prs, data.frame(lasso = tuning_scores$lasso), prs_cols),
    step = mors_tune_step_model(y_tuning_resid, tuning_prs, tuning_omics_df, prs_cols)
  )
}

mors_build_test_predictions <- function(tuned_models, test_prs, test_scores) {
  build_prediction <- function(prs_col, fit, extra_scores = NULL) {
    design <- data.frame(prs = test_prs[[prs_col]])
    names(design)[1] <- prs_col
    
    if (!is.null(extra_scores)) {
      design <- cbind(design, extra_scores)
    }
    
    stats::predict(fit, newdata = design)
  }
  
  test_omics_df <- data.frame(
    trs = test_scores$trs,
    prors = test_scores$prors,
    pmrs = test_scores$pmrs,
    smrs = test_scores$smrs
  )
  
  step_model <- tuned_models$step$fit
  step_terms <- attr(stats::terms(step_model), "term.labels")
  step_newdata <- data.frame(test_prs[[tuned_models$step$prs]], test_omics_df)
  names(step_newdata)[1] <- tuned_models$step$prs
  step_newdata <- step_newdata[, step_terms, drop = FALSE]
  
  list(
    base = build_prediction(tuned_models$base$prs, tuned_models$base$fit),
    trs = build_prediction(tuned_models$trs$prs, tuned_models$trs$fit, data.frame(trs = test_scores$trs)),
    prors = build_prediction(tuned_models$prors$prs, tuned_models$prors$fit, data.frame(prors = test_scores$prors)),
    pmrs = build_prediction(tuned_models$pmrs$prs, tuned_models$pmrs$fit, data.frame(pmrs = test_scores$pmrs)),
    smrs = build_prediction(tuned_models$smrs$prs, tuned_models$smrs$fit, data.frame(smrs = test_scores$smrs)),
    step = stats::predict(step_model, newdata = step_newdata),
    lasso = build_prediction(tuned_models$lasso$prs, tuned_models$lasso$fit, data.frame(lasso = test_scores$lasso))
  )
}

mors_summarize_results <- function(y_test_resid, model_scores) {
  data.frame(
    model = c("base", "trs", "prors", "pmrs", "smrs", "step", "lasso"),
    conditional_r2 = c(
      mors_conditional_r2(y_test_resid, model_scores$base),
      mors_conditional_r2(y_test_resid, model_scores$trs),
      mors_conditional_r2(y_test_resid, model_scores$prors),
      mors_conditional_r2(y_test_resid, model_scores$pmrs),
      mors_conditional_r2(y_test_resid, model_scores$smrs),
      mors_conditional_r2(y_test_resid, model_scores$step),
      mors_conditional_r2(y_test_resid, model_scores$lasso)
    ),
    row.names = NULL
  )
}

mors_summarize_bootstrap <- function(results, y_test_resid, model_scores, B, seed) {
  boot_mat <- mors_bootstrap_r2(y_test_resid, model_scores, B, seed)
  results_boot_ci <- mors_bootstrap_ci(boot_mat)
  
  single_omics_models <- c("trs", "prors", "pmrs", "smrs")
  single_omics_r2 <- results$conditional_r2[match(single_omics_models, results$model)]
  single_omics_best <- single_omics_models[which.max(single_omics_r2)]
  
  results_boot_pv <- data.frame(
    comparison = c(
      "lasso_step",
      "lasso_best_single_omics",
      "step_best_single_omics",
      "lasso_base",
      "step_base",
      "trs_base",
      "prors_base",
      "pmrs_base",
      "smrs_base"
    ),
    p_value = c(
      mors_bootstrap_pvalue(boot_mat[, "lasso"] - boot_mat[, "step"]),
      mors_bootstrap_pvalue(boot_mat[, "lasso"] - boot_mat[, single_omics_best]),
      mors_bootstrap_pvalue(boot_mat[, "step"] - boot_mat[, single_omics_best]),
      mors_bootstrap_pvalue(boot_mat[, "lasso"] - boot_mat[, "base"]),
      mors_bootstrap_pvalue(boot_mat[, "step"] - boot_mat[, "base"]),
      mors_bootstrap_pvalue(boot_mat[, "trs"] - boot_mat[, "base"]),
      mors_bootstrap_pvalue(boot_mat[, "prors"] - boot_mat[, "base"]),
      mors_bootstrap_pvalue(boot_mat[, "pmrs"] - boot_mat[, "base"]),
      mors_bootstrap_pvalue(boot_mat[, "smrs"] - boot_mat[, "base"])
    ),
    row.names = NULL
  )
  
  list(
    results_boot_ci = results_boot_ci,
    results_boot_pv = results_boot_pv
  )
}

MoRS <- function(
    phenotype_train_file,
    phenotype_tuning_file,
    phenotype_test_file,
    gpf_train_file,
    gpf_tuning_file,
    gpf_test_file,
    prs_tuning_file,
    prs_test_file,
    outcome_col,
    covariate_cols,
    id_col = "IID",
    n_prs_cols = NULL,
    trs_cols,
    prors_cols,
    pmrs_cols,
    smrs_cols,
    family = "gaussian",
    nfolds = 5,
    B = 1000,
    seed = 123
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required.")
  }
  
  inputs <- mors_prepare_inputs(
    phenotype_train_file = phenotype_train_file,
    phenotype_tuning_file = phenotype_tuning_file,
    phenotype_test_file = phenotype_test_file,
    gpf_train_file = gpf_train_file,
    gpf_tuning_file = gpf_tuning_file,
    gpf_test_file = gpf_test_file,
    prs_tuning_file = prs_tuning_file,
    prs_test_file = prs_test_file,
    outcome_col = outcome_col,
    covariate_cols = covariate_cols,
    id_col = id_col,
    trs_cols = trs_cols,
    prors_cols = prors_cols,
    pmrs_cols = pmrs_cols,
    smrs_cols = smrs_cols,
    n_prs_cols = n_prs_cols
  )
  
  y_train_resid <- mors_residualize_outcome(inputs$train_dat, outcome_col, covariate_cols)
  y_tuning_resid <- mors_residualize_outcome(inputs$tuning_dat, outcome_col, covariate_cols)
  y_test_resid <- mors_residualize_outcome(inputs$test_dat, outcome_col, covariate_cols)
  
  fitted_models <- mors_fit_all_lasso_models(
    train_dat = inputs$train_dat,
    tuning_dat = inputs$tuning_dat,
    test_dat = inputs$test_dat,
    omics_groups = inputs$omics_groups,
    all_omics_cols = inputs$all_omics_cols,
    y_train_resid = y_train_resid,
    family = family,
    nfolds = nfolds,
    seed = seed
  )
  
  tuned_models <- mors_tune_models(
    y_tuning_resid = y_tuning_resid,
    tuning_prs = inputs$tuning_prs,
    prs_cols = inputs$prs_cols,
    tuning_scores = fitted_models$tuning_scores
  )
  
  model_scores <- mors_build_test_predictions(
    tuned_models = tuned_models,
    test_prs = inputs$test_prs,
    test_scores = fitted_models$test_scores
  )
  
  results <- mors_summarize_results(y_test_resid, model_scores)
  
  boot_summary <- mors_summarize_bootstrap(
    results = results,
    y_test_resid = y_test_resid,
    model_scores = model_scores,
    B = B,
    seed = seed
  )
  
  list(
    results = results,
    results_boot_ci = boot_summary$results_boot_ci,
    results_boot_pv = boot_summary$results_boot_pv,
    lasso_models = fitted_models$lasso_models,
    step_model = tuned_models$step$fit
  )
}