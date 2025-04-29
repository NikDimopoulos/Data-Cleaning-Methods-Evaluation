# Data-Cleaning-Methods-Evaluation
Exploring the impact of outlier removal on multivariate testing 

This repository contains the R code for a simulation study that evaluates how robust statistical methods affect multivariate location testing, correlation structure testing, and linear regression inference in a dataset that has been contaminated by outliers.

---

##Overview

This paper explores how different outlier detection methods affect the outcomes of the following multivariate statistical tests:

- Hotelling's T² Test
- Likelihood Ratio (LR) Tests for location and correlation
- F-Test for multiple linear regression

 A multivariate dataset is simulated and contaminated with outliers. Three different outlier detection methods are then applied:

1. **Univariate Screening using Standard Deviation**
2. **Univariate Screening using Median Absolute Deviation (MAD)**
3. **Minimum Covariance Determinant (MCD)**

The test statistics mentioned above are computed for both the contaminated and cleaned datasets, and the results are then compared to evaluate the influence of these data cleaning techniques on multivariate testing.

  ## Contents

- **Functions Included**:
  - `generate_data()`: Simulates multivariate normal data with contamination.
  - `univariate_screening()`, `univariate_screening_mad()`, `mcd_cleaning()`: Data cleaning functions.
  - `hotelling_test()`, `lr_test_location()`, `lr_test_correlation()`: Hypothesis testing functions.
  - `response_variable()`, `ftest()`: Generate response variable and run F-test for linear regression.
  - `simulate()`: Runs the simulation loop for a defined number of replications (`R`).
  - Data visualization using `ggplot2`.
      
