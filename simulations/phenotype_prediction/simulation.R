# run_simulation.R
#
# Description:
# This script simulates genotype-phenotype data to compare the predictive
# performance of different methods. The core comparison is between
# standard Elastic Net and ePRS, which leverages p-values from an external
# source GWAS to improve prediction in a smaller target cohort.
# The simulation varies the genetic correlation (r_g) between the source and
# target phenotypes to assess robustness.

# --- 1. Load Libraries---
library(glmnet)
library(MASS)
library(ggplot2)

# --- 2. Simulation Parameters ---
set.seed(100)
sim_count <- 40           # Number of full simulation runs
iter_count <- 40          # Number of different r_g values per run
n_row <- 300              # Sample size of target cohort
n_col <- 1000             # Number of simulated SNPs

# Define the SNP correlation structure (AR1 with rho=0)
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - (1:n - 1))
  rho^exponent
}
corr_matrix <- ar1_cor(n_col, 0)

# --- 3. Initialize Result Matrices ---
# We will compare three main methods:
# - p_t: P+T using target data for beta estimation
# - en: Standard Elastic Net
# - eprs_en: ePRS with Elastic Net

r_g_values <- p_t <- en <- eprs <- matrix(nrow = sim_count, ncol = iter_count)

# --- 4. Main Simulation Loop ---
for (sim in 1:sim_count) {

  # Source data (larger N)
  df_source <- mvrnorm(n = n_row * 5, mu = rep(0, n_col), Sigma = corr_matrix)
  # Target data (smaller N)
  df_target <- mvrnorm(n = n_row, mu = rep(0, n_col), Sigma = corr_matrix)

  # True effect sizes for source phenotype
  beta_source <- c(rep(1, iter_count), rep(0, n_col - iter_count))

  for (iter in 1:iter_count) {
    # --- Define Target Effect Sizes and Genetic Correlation ---
    # Slide the block of causal variants to control overlap with beta_source
    beta_target <- c(rep(0, iter), rep(1, iter_count), rep(0, n_col - iter_count - iter))
    r_g_values[sim, iter] <- sum(beta_source * beta_target) / sqrt((sum(beta_source^2) * sum(beta_target^2)))

    # --- Simulate Phenotypes ---
    y_source <- df_source %*% beta_source + rnorm(nrow(df_source), 0, 0)   
    y_target <- df_target %*% beta_target + rnorm(n_row, 0, 4)              

    # --- Run "External" Source GWAS ---
    # This generates the p-values that ePRS will use as external evidence.
    pval <- sapply(1:n_col, function(i) {
      summary(lm(y_source ~ df_source[, i]))$coef[2, 4]
    })

    # --- Split Target Data into Training/Validation/Test Sets ---
    # For tuning alpha and final evaluation
    train_idx <- 1:160
    valid_idx <- 161:200
    test_idx <- 201:300

    y_train <- y_target[train_idx]
    y_valid <- y_target[valid_idx]
    y_test <- y_target[test_idx]

    # --- Method 1: P+T (Clumping and Thresholding) using Target Betas ---
    effect_size = pval_target = array(dim=n_col)
    for(i in 1:n_col){
      model = summary(lm(y_train ~ df_target[train_idx,i]))$coef
      effect_size[i] = model[2,1]; pval_target[i] = model[2,4]
    }

    p_threshold = c(1, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01, 0.001); to_compare = NULL
    for(pthr in p_threshold){
      selected_snps = which(pval < pthr)
      y_pred = df_target[valid_idx, selected_snps] %*% as.matrix(effect_size[selected_snps])
      to_compare = append(to_compare, (cor(y_valid, y_pred)))
    }
    pthr = p_threshold[which.max(to_compare)]
    selected_snps = which(pval < pthr)
    y_pred = df_target[test_idx, selected_snps] %*% as.matrix(effect_size[selected_snps])
    p_t[sim, iter] = cor(y_test, y_pred)

    # --- Method 2: Standard Elastic Net (EN) ---
    # Tune alpha on a validation set
    alpha_thr <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
    validation_cors <- sapply(alpha_thr, function(alpha_tune) {
      m <- cv.glmnet(df_target[train_idx, ], y_train, nfolds = 3, alpha = alpha_tune)
      y_pred_valid <- predict(m, df_target[valid_idx, ], s = "lambda.min")
      cor(y_valid, y_pred_valid)
    })
    best_alpha_en <- alpha_thr[which.max(validation_cors)]
    
    # Train final model on combined train+valid set and evaluate on test set
    final_en_model <- cv.glmnet(df_target[c(train_idx, valid_idx), ], y_target[c(train_idx, valid_idx)], nfolds = 3, alpha = best_alpha_en)
    y_pred_en <- predict(final_en_model, df_target[test_idx, ], s = "lambda.min")
    en[sim, iter] <- cor(y_test, y_pred_en)

    # --- Method 3: ePRS with Elastic Net ---
    # Define the evidence term E_j from source p-values.
    # We use a transformation that is more sensitive to p-value changes than -log10(p).
    E_j <- 1/-log10(pval)
    # E_j <- 10 * (1 - (1 - pval)^8)

    # The penalty factor combines external evidence (E_j) with genetic correlation (r_g)
    pf <- E_j * r_g_values[sim, iter] + (1 - r_g_values[sim, iter])

    # Tune alpha for ePRS
    validation_cors_eprs <- sapply(alpha_thr, function(alpha_tune) {
      m <- cv.glmnet(df_target[train_idx, ], y_train, nfolds = 3, alpha = alpha_tune, penalty.factor = pf)
      y_pred_valid <- predict(m, df_target[valid_idx, ], s = "lambda.min")
      cor(y_valid, y_pred_valid)
    })
    best_alpha_eprs <- alpha_thr[which.max(validation_cors_eprs)]

    # Train final ePRS model and evaluate on test set
    final_eprs_model <- cv.glmnet(df_target[c(train_idx, valid_idx), ], y_target[c(train_idx, valid_idx)], nfolds = 3, alpha = best_alpha_eprs, penalty.factor = pf)
    y_pred_eprs <- predict(final_eprs_model, df_target[test_idx, ], s = "lambda.min")
    eprs[sim, iter] <- cor(y_test, y_pred_eprs)
    
    print(paste("sim:", sim, "| Iter:", iter, "| r_g:", round(r_g_values[sim, iter], 2)))
  }
}

# --- 5. Plotting Results ---
df = data.frame(r_sq=c(p_t,en,eprs),
                rg=rep(r_g, 3),
                method=rep(c("P+T","EN","ePRS"), each=1600))


ggplot(df, aes(x=rg, y=r_sq)) + geom_point(aes(col=method)) +
  geom_smooth(col="black", lty="dashed") +
  facet_grid(method~.)

ggplot(df, aes(x=rg, y=r_sq)) +
  geom_smooth(aes(col=method), lty="dashed", se=F)
