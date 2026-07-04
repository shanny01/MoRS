# MoRS (Multi-omics Risk Score)

**MoRS** is a multi-omics integration method for improving polygenic risk prediction in the absence of individual-level omics data.
This package contains two main function: step-MoRS, which integrates single-omics scores using stepwise regression, and Lasso-MoRS, which directly models all predicted features across omics layers using LASSO regression. 

# Tutorial
## Getting Started
Clone this repository by 
```
git clone https://github.com/shanny01/MoRS.git

# necessary R packages
Rscript -e 'install.packages(c('data.table', 'dplyr', 'glmnet'))'
```

## Prepare the baseline PRSs

The clumping and thresholding (CT) approach was implemented using PLINK 1.9.0 with a window size of 250 kb and a linkage disequilibrium cutoff of 0.1 in the clumping step. The optimal P-value threshold was selected from a range of candidate values (1.0, 5×10-1, ⋯ , 5×10-8). PRSCS was implemented using Python 3.12.2 following the pipeline https://github.com/getian107/PRScs. The optimal global shrinkage parameter was determined from the default candidate values, 1.0, 10-2 ,10-4, and 10-6.


## Prepare genetically predicted multi-omics features

These levels were derived as a linear sum of genotype dosages weighted by effect sizes retrieved from the OmicsPred portal (www.omicspred.org; accession codes OPGS000001-OPGS017227). The calculations were implemented using the pgsc_calc pipline (https://github.com/PGScatalog/pgsc_calc), applying a default variant matching threshold of 0.75 between the weight files and target genomic data. 


## Inputs of MoRS

1. `phenotype_train_file`, `phenotype_tuning_file`, `phenotype_test_file` :
It contains the information of IID, outcome (Y) and covariates. 

2. 'gpf_train_file', 'gpf_tuning_file', 'gpf_test_file':
It contains the information of IID and genetically predicted multi-omics features.

3. 'prs_tuning_file', 'prs_test_file':
It contains the information of IID and candidate PRSs.

Note that we use the training data to select genetically predicted multi-omics features and train weights for calculating single-omics scores and Lasso-MoRS. The tuning data are used to choose the best tuning parameter in PRSs. The test data are used to assess the model performance. 

## Outputs of MoRS

1. `results`: 
the adjusted $R^2$ of all competing methods.

2. `results_boot.ci`: 
95% confidence intervals of all competing methods.

3. `results_boot.pv`: 
This component containing p values by one-sided paired bootstrap Z tests for paired differences in adjusted $R^2$ between two models for comparison. 

4. `lasso_model`: 
the list containing features selected in single-omics score and Lasso-MoRS.  

5. `step_model`: 
the list containing features selected in step-MoRS. 





## File Requirements

MoRS expects three types of input files.

### 1. Phenotype files

These files are used for train, tuning, and test datasets:
  
- `phenotype_train_file`
- `phenotype_tuning_file`
- `phenotype_test_file`

They must contain:
  
- an ID column, such as `IID`
- the outcome column
- all covariates used for residualization

Example columns:
  
```text
IID phenotype age sex PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10
```

### 2. GPF files

These files contain genetically predicted multi-omics features for the train, tuning, and test datasets:
  
- `gpf_train_file`
- `gpf_tuning_file`
- `gpf_test_file`

They must contain:
  
- the same ID column
- omics feature columns grouped into:
- `trs_cols`
- `prors_cols`
- `pmrs_cols`
- `smrs_cols`

These groups must be passed to `MoRS()` as vectors of column names.

### 3. PRS files

These files contain candidate PRSs for tuning and test datasets:
  
- `prs_tuning_file`
- `prs_test_file`

They must contain:
  
- the same ID column
- one or more PRS columns

Example columns:
  
```text
IID PRS1 PRS2 PRS3 PRS4 PRS5
```

## Main Function

```r
MoRS(
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
)
```

## Arguments

- `phenotype_train_file`  
Path to the phenotype training file.

- `phenotype_tuning_file`  
Path to the phenotype tuning file.

- `phenotype_test_file`  
Path to the phenotype test file.

- `gpf_train_file`  
Path to the genetically predicted feature training file.

- `gpf_tuning_file`  
Path to the genetically predicted feature tuning file.

- `gpf_test_file`  
Path to the genetically predicted feature test file.

- `prs_tuning_file`  
Path to the PRS tuning file.

- `prs_test_file`  
Path to the PRS test file.

- `outcome_col`  
Name of the outcome column in phenotype files.

- `covariate_cols`  
Character vector of covariate column names used to residualize the outcome.

- `id_col`  
Name of the subject ID column. Default is `"IID"`.

- `n_prs_cols`  
Optional integer limiting how many PRS columns from `prs_tuning_file` and `prs_test_file` are used.  
If `NULL`, all PRS columns are used.

- `trs_cols`  
Character vector of feature names for the `trs` omics layer.

- `prors_cols`  
Character vector of feature names for the `prors` omics layer.

- `pmrs_cols`  
Character vector of feature names for the `pmrs` omics layer.

- `smrs_cols`  
Character vector of feature names for the `smrs` omics layer.

- `family`  
Model family passed to `glmnet::cv.glmnet()`. Default is `"gaussian"`.

- `nfolds`  
Number of folds for lasso cross-validation. Default is `5`.

- `B`  
Number of bootstrap replicates for confidence intervals and p-values. Default is `1000`.

- `seed`  
Random seed for lasso cross-validation and bootstrap. Default is `123`.

## Output

`MoRS()` returns a list with the following components:
  
  - `results`  
A data frame of test-set conditional `R^2` values for:
  `base`, `trs`, `prors`, `pmrs`, `smrs`, `step`, and `lasso`.

- `results_boot_ci`  
A data frame containing bootstrap 95% confidence intervals for test-set conditional `R^2`.

- `results_boot_pv`  
A data frame containing bootstrap-based p-values for key comparisons:
- `lasso_step`
- `lasso_best_single_omics`
- `step_best_single_omics`
- `lasso_base`
- `step_base`
- `trs_base`
- `prors_base`
- `pmrs_base`
- `smrs_base`

- `lasso_models`  
A list of fitted `cv.glmnet` models for:
- `trs`
- `prors`
- `pmrs`
- `smrs`
- `combined`

- `step_model`  
The final forward stepwise `lm` model selected in the tuning dataset.

## Example

```r
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
```

## Support
Please send questions and issues related to MoRS software to Nayang Shan (shanny01@foxmail.com)