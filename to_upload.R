## Introduction to ePRS

library(glmnet) 
library(MASS)

set.seed(100)

# AR1 Correlation Matrix
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
corr_matrix_source <- ar1_cor(n_col, 0); corr_matrix_target <- ar1_cor(n_col, 0)

##########################################################################

# Simulating source/target datasets
n_row <- 300; n_col <- 1000
df_source <- mvrnorm(n=n_row*4, mu=rep(0, n_col), Sigma=corr_matrix_source)
df_target <- mvrnorm(n=n_row, mu=rep(0, n_col), Sigma=corr_matrix_target)

overlap <- 10                                                                             # Controls strength overlapping genetic signals: overlap=[0,40]
beta_source <- c(rep(1, 40), rep(0, n_col-40)) 
beta_target <- c(rep(0, overlap), rep(1, 40), rep(0, n_col-40-overlap))         
r_g <- sum(beta_source*beta_target) / sqrt((sum(beta_source^2) * sum(beta_target^2)))     # Genetic correlation between source/target phenotypes (r_g=0.75 when overlap=10)
y_source <- df_source %*% beta_source + rnorm(n_row, 0, 5)                                # Source phenotype
y_target <- df_target %*% beta_target + rnorm(n_row, 0, 5)                                # Target phenotype


# Derive p-values and the subsequent E_j values using the source phenotype
pval <- array(dim=n_col)
for(i in 1:n_col){
  model <- summary(lm(y_source ~ df_source[,i]))$coef
  pval[i] <- model[2,4]
}
E_j <- 1/-log10(pval)


####################################################################################
## ePRS
# Define alpha values to test
alpha_grid <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
cor_results <- numeric(length(alpha_grid))

# Define training and test indices
train_idx <- 1:160
test_idx <- 161:200 

# Loop over alpha values
for (i in seq_along(alpha_grid)) {
  alpha_val <- alpha_grid[i]
  
  # Fit penalized regression model with given alpha
  fit <- cv.glmnet(
    x = df_target[train_idx, ],
    y = y_target[train_idx],
    nfolds = 3,
    alpha = alpha_val,
    penalty.factor = E_j * r_g + (1 - r_g)
  )
  
  # Predict and compute correlation with true test labels
  preds <- as.vector(predict(fit, newx = df_target[test_idx, ], s = "lambda.min"))
  cor_results[i] <- cor(y_target[test_idx], preds)
}

# Select best alpha
best_alpha <- alpha_grid[which.max(cor_results)]
if (length(best_alpha) == 0) best_alpha <- 0  # fallback

# Refit on full training set
final_train_idx <- 1:200
final_test_idx <- 201:300

final_fit <- cv.glmnet(
  x = df_target[final_train_idx, ],
  y = y_target[final_train_idx],
  nfolds = 3,
  alpha = best_alpha,
  penalty.factor = E_j * r_g + (1 - r_g)
)

# Predict on final test set
y_pred <- as.vector(predict(final_fit, newx = df_target[final_test_idx, ], s = "lambda.min"))

# ePRS model performance; measured in R^2
eprs_en <- cor(y_target[201:300], y_pred)
print(eprs_en)

