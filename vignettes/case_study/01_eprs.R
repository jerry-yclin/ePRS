# 01_eprs.R
#
# Description:
# This script builds a polygenic score model using ePRS
# An ePRS model that enhances the EN model with the IGE summary statistics.

# --- 1. Load Libraries & Data ---
library(glmnet)
library(pROC)
library(data.table)
library(dplyr)
library(bigsnpr) 

# Load data objects 
# Load "big_SNP" file containing genotype data (derived from library "bigsnpr")
#rds <- snp_attach("genotype.rds")
#G <- rds$genotype
# Load Phenotype data
#pheno <- fread("pheno.txt")

# --- Placeholder Simulation ---
n_samples <- 3624; n_snps <- 100000
fake_snps <- snp_fake(n_samples, n_snps)
G <- fake_snps$genotypes
G[] <- rbinom(n=n_samples * n_snps, size=2, prob=0.3)
pheno <- data.frame(STATUS = c(rep(1, 624), rep(0, 3000)), sex = rbinom(n_samples, 1, 0.5), pc = rnorm(n_samples)))
# --- End Placeholder ---

# Load train/test indices
train_indices <- c(sample(which(pheno$STATUS == 1), size = 312),  # 50% of cases
                   sample(which(pheno$STATUS == 0), size = 1500)) # 50% of controls
test_indices <- setdiff(1:nrow(pheno), train_indices)

# Load external IGE summary statistics
# IGE_summary <- fread("generalised_epilepsy_METAL_processed.tsv")
# --- Placeholder/Simulated IGE Summary Stats ---
IGE_summary <- data.frame(MarkerName = 1:n_snps, "P-value" = runif(n_snps)^2, beta = rnorm(n_snps, 0, 0.05))
# --- End Placeholder ---

# --- 2. Feature Selection & Data Preparation ---
# Optional In-sample GWAS
gwas_pval <- array(dim=100000)
for(i in ncol(G)){
  mod <- summary(glm(pheno_train$STATUS ~ G[train_indices,i] + pheno_train$sex + pheno_train$pc, family="binomial"))       # Only 1 PC included for illustrative purposes 
  gwas_pval[i] <- mod$coef[2,4]
}

# Optional pre-filtering step reduces the computational burden and noise.
pval_thresh_filter <- 1
snps_to_include <- which(gwas_pval < pval_thresh_filter)

# Create model matrix (genotypes + covariates)
covariates <- as.matrix(pheno[, c("sex", paste0("pc", 1:10))])
x_data <- cbind(G[, snps_to_include], covariates)
y_data <- pheno$STATUS

x_train <- x_data[train_indices, ]
y_train <- y_data[train_indices]
x_test <- x_data[test_indices, ]
y_test <- y_data[test_indices]

# The covariates should not be penalized
n_snps_selected <- length(snps_to_include)
penalty_factors <- c(rep(1, n_snps_selected), rep(0, ncol(covariates)))

# --- ePRS with Elastic Net ---
# Define the evidence term (E_j) from the external GWAS p-values.
ige_pvals <- IGE_summary$"P-value"[snps_to_include]
#E_j <- 10 * (1 - (1 - ige_pvals)^8)
E_j <- -1 / log10(ige_pvals) 

# Define the genetic correlation (r_g). Can be estimated using LD Score Regression on the training set of the target cohort.
r_g <- 0.5

# Construct the final penalty vector for ePRS
eprs_penalty_factors <- E_j * r_g + (1 - r_g)
eprs_penalty_factors <- c(eprs_penalty_factors, rep(0, ncol(covariates)))     # Do not penalize covariates

# Train the ePRS model
cv_eprs <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.8,      
                     penalty.factor = eprs_penalty_factors)
pred_eprs <- predict(cv_eprs, x_test, s = "lambda.min", type = "response")
auc_eprs <- roc(y_test, as.vector(pred_eprs), quiet = TRUE)$auc
print(paste("ePRS model AUC:", round(auc_eprs, 4)))


# --- Save Model Predictions for Analysis ---
# Store the predicted probabilities/scores
model_predictions <- data.frame(
  IID = pheno$IID[test_indices],
  y_true = y_test,
  pred_eprs = as.vector(pred_eprs),
)


