# External PRS (ePRS): Enhancing Clinical Utility of Polygenic Scores with Small, Phenotypically Refined Cohorts

**ePRS** is a method to improve the accuracy and stability of polygenic risk scores (PRS) in clinical cohorts / phenotypes with limited sample sizes. It leverages summary statistics from large-scale, external Genome-Wide Association Studies (GWAS) of related phenotypes to inform a weighted, penalized regression model.

This repository provides the R code to implement ePRS and replicate the analyses from our research.

## Key Features

- **Improved Predictive Power**: Boosts prediction accuracy by borrowing information from a well-powered "source" GWAS.
- **Enhanced Stability**: Produces more stable coefficient estimates and risk predictions, especially when the target dataset is small.
- **Robustness**: Provides gains even when the genetic correlation (`r_g`) between the source and target phenotypes is moderate, and gracefully converges to a standard elastic net model when `r_g` is zero.
- **Flexibility**: Built on the powerful `glmnet` framework, allowing for both Ridge (`alpha=0`) and Elastic Net (`alpha` > 0) regularization.

## How it Works: The Core Idea

Standard penalized regression (like Ridge or Elastic Net) treats all genetic variants equally. ePRS introduces a simple but effective modification: it penalizes variants differently based on their statistical significance in an external "source" GWAS.

The penalty for each variant *j* is weighted using a `penalty.factor` in `glmnet`:

`penalty.factor_j = E_j * r_g + (1 - r_g)`

- **`E_j`**: A measure of evidence from the source GWAS. A common choice is `1 / -log10(p-value)`. Variants with smaller p-values in the source study receive a smaller `E_j`, and thus a smaller penalty, encouraging the model to include them.
- **`r_g`**: The estimated genetic correlation between the source and target phenotypes. This acts as a fail-safe mechanism, controlling how much we trust the external evidence. A high `r_g` puts more weight on the source GWAS p-values.

## Usage and Examples

This repository is organized into two main sections: simulations to demonstrate the method's properties and a real-world vignette.

### 1. Simulation: Phenotype Prediction

This simulation demonstrates how ePRS outperforms standard methods (P+T, Elastic Net) in predicting a quantitative trait across a range of genetic correlations.

**Goal**: Predict a target phenotype using a small training set, leveraging a GWAS of a related source phenotype.
**Location**: [`simulations/phenotype_prediction/`](simulations/phenotype_prediction/)


### 2. Simulation: Disease Subtype Differentiation

This simulation shows how ePRS can better differentiate between two disease subtypes, even when the training data only contains one subtype versus controls.

**Goal**: Differentiate Subtype A from Subtype B, using a general GWAS (of A+B combined) as the source and a training set of (Subtype A vs. Healthy) as the target.
**Location**: [`simulations/subtype_differentiation/`](simulations/subtype_differentiation/)

### 3. Case Study: ePRS through a toy study

This vignette provides a step-by-step workflow for applying ePRS using a toy example: building a target-specific polygenic score for a simulated phenotype.

**Location**: [`vignettes/case_study/`](vignettes/case_study/)


