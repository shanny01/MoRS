## Example:
## 1. phenotype files contain: IID, phenotype, age, sex, PC1-PC10
## 2. gpf files contain: IID and genetically predicted multi-omics features
## 3. prs files contain: IID and candidate PRSs
## 4. trs_cols / prors_cols / pmrs_cols / smrs_cols must be provided as feature names

## Read the training GPF file to define omics feature groups
gpf_train <- data.frame(
  data.table::fread("data/gpf_train.txt", header = TRUE, data.table = FALSE)
)

## Replace these patterns with your real column naming rules
trs_cols <- grep("^trs_", names(gpf_train), value = TRUE)
prors_cols <- grep("^prors_", names(gpf_train), value = TRUE)
pmrs_cols <- grep("^pmrs_", names(gpf_train), value = TRUE)
smrs_cols <- grep("^smrs_", names(gpf_train), value = TRUE)

fit <- MoRS(
  phenotype_train_file = "data/phenotype_train.txt",
  phenotype_tuning_file = "data/phenotype_tuning.txt",
  phenotype_test_file = "data/phenotype_test.txt",
  gpf_train_file = "data/gpf_train.txt",
  gpf_tuning_file = "data/gpf_tuning.txt",
  gpf_test_file = "data/gpf_test.txt",
  prs_tuning_file = "data/prs_tuning.txt",
  prs_test_file = "data/prs_test.txt",
  outcome_col = "phenotype",
  covariate_cols = c("age", "sex", paste0("PC", 1:10)),
  id_col = "IID",
  n_prs_cols = 9,
  trs_cols = trs_cols,
  prors_cols = prors_cols,
  pmrs_cols = pmrs_cols,
  smrs_cols = smrs_cols,
  family = "gaussian",
  nfolds = 5,
  B = 1000,
  seed = 123
)

output <- list(
  results = fit$results,
  results_boot_ci = fit$results_boot_ci,
  results_boot_pv = fit$results_boot_pv,
  lasso_models = fit$lasso_models,
  step_model = fit$step_model
)

print(output$results)
print(output$results_boot_ci)
print(output$results_boot_pv)

saveRDS(output, file = "example/mors_output.rds")