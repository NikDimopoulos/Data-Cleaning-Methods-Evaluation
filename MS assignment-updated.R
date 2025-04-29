library(ggplot2)
library(mvtnorm)
library(robustbase)
library(dplyr)
library(Hotelling)
library(DescTools)
library(tidyverse)
library(car)

# Generate simulations
generate_data <- function(n, mu, sigma, mu_con, sigma_con, eps) {
  x <- rmvnorm(n, mu, sigma)
  x_con <- rmvnorm(n, mu_con, sigma_con)
  contamination <- rbinom(n, 1, eps)
  for (i in 1:ncol(x)) {
    x[, i] <- (1 - contamination) * x[, i] + contamination * x_con[, i]
  }
  as.data.frame(x)
}

# Univariate screening methods
univariate_screening <- function(data, c) {
  z <- scale(data)
  rows_to_keep <- apply(z, 1, function(row) all(abs(row) <= c))
  data[rows_to_keep, , drop = FALSE]
}

univariate_screening_mad <- function(data, c) {
  medians <- apply(data, 2, median)
  mads <- apply(data, 2, mad)
  z <- scale(data, center = medians, scale = mads)
  rows_to_keep <- apply(abs(z), 1, function(row) all(row <= c))
  data[rows_to_keep, , drop = FALSE]
}

# MCD-based cleaning
mcd_cleaning <- function(data, alpha) {
  mcd <- covMcd(data)
  mahal_distances <- mahalanobis(data, mcd$center, mcd$cov)
  threshold <- qchisq((1 - alpha / 100), ncol(data))
  rows_to_keep <- mahal_distances <= threshold
  data[rows_to_keep, , drop = FALSE]
}

# Hotelling's T2 test using Hotelling package
hotelling_test <- function(data, mu0) {
  #test <- hotelling.test(as.matrix(data),y=NULL, mu = mu0)
  test<-HotellingsT2Test(data, y = NULL, mu = mu0, test = "f")
  statistic=test$statistic
  p_value=test$p.value
  list(mean_statistic = statistic, mean_p_value =p_value,mean_significant = p_value < 0.05)
}



# Likelihood ratio tests
# Likelihood ratio tests
lr_test_location <- function(data, mu0) {
  p <- ncol(data)
  n <- nrow(data)
  Sigma_null_hat <- cov(data - mu0)
  xbar <- colMeans(data)
  Sigma_hat <- cov(data - xbar)
  lambda <- exp((n/2)* (log(det(Sigma_hat))-log(det(Sigma_null_hat))))
    statistic <- -2 * log(lambda)
    p_value <- pchisq(statistic, df = p, lower.tail = FALSE)
    list(loc_statistic = statistic, loc_p_value = p_value, loc_significant = p_value < 0.05)
}

lr_test_correlation <- function(data) {
  rho_hat <- cor(data)
  n <- nrow(data)
  p <- ncol(data)
  statistic <- -n * log(det(rho_hat))
  df <- p * (p - 1) / 2
  p_value <- pchisq(statistic, df, lower.tail = FALSE)
  list(cor_statistic = statistic, cor_p_value = p_value, cor_significant = p_value < 0.05)
}


#Generate a response variable y
response_variable<-function(X,beta){
  n=nrow(X)
  epsilon <-rnorm(n, mean=0, sd=3)
  y <- as.matrix(X) %*% beta + epsilon
  X$y <- y 
  return(X)
}

#Ftest for the regression
ftest<-function(data){
  colnames(data) <- c("X1", "X2", "X3", "X4","y")
  linear_model=lm(y~ ., data = data)
  result<-linearHypothesis(linear_model, c("X1 = 1", "X2 = 1", "X3 = 1", "X4 = 1"),
                           vcov = vcov(linear_model), 
                           test = "F")
  return(result)
}

# Descriptive statistics
print_descriptive_stats <- function(data, title) {
  cat("\nDescriptive statistics for:", title, "\n")
  print(summary(data))
  print(cor(data))
}

# Plot data
plot_data <- function(data, title) {
  pairs(data, main = title, pch = 19, col = "blue")
}

# Simulation study
simulate <- function(R, n, mu, sigma, mu_con, sigma_con, eps) {
  results <- replicate(R, {
    data <- generate_data(n, mu, sigma, mu_con, sigma_con, eps)
    
    # Data cleaning methods
    clean_sd <- univariate_screening(data, 3)
    #print_descriptive_stats(clean_sd, "Cleaned Data (Standard Deviation)")
    #plot_data(clean_sd, "Cleaned Data (Standard Deviation)")
    
    clean_mad <- univariate_screening_mad(data, 3)
    #print_descriptive_stats(clean_mad, "Cleaned Data (MAD)")
    #plot_data(clean_mad, "Cleaned Data (MAD)")
    
    clean_mcd <- mcd_cleaning(data, 2.5)
    #print_descriptive_stats(clean_mcd, "Cleaned Data (MCD)")
    #plot_data(clean_mcd, "Cleaned Data (MCD)")
    
    mu0 <- rep(0, ncol(data))
    
    
    # Test statistics
    hotelling_results <- list(
      Original = hotelling_test(data, mu0),
      SD = hotelling_test(clean_sd, mu0),
      MAD = hotelling_test(clean_mad, mu0),
      MCD = hotelling_test(clean_mcd, mu0)
    )
    
    lr_location_results <- list(
      Original = lr_test_location(data, mu0),
      SD = lr_test_location(clean_sd, mu0),
      MAD = lr_test_location(clean_mad, mu0),
      MCD = lr_test_location(clean_mcd, mu0)
    )
    
    lr_correlation_results <- list(
      Original = lr_test_correlation(data),
      SD = lr_test_correlation(clean_sd),
      MAD = lr_test_correlation(clean_mad),
      MCD = lr_test_correlation(clean_mcd)
    )
    
    #generate the response variable y and add it to the data
    beta=rep(1,ncol(data))
    regression_data<-response_variable(data,beta)
    cleaned_regression_data<-univariate_screening_mad(regression_data,2)
    
    regression_results<-list(
      Original=ftest(regression_data),
      cleaned=ftest(cleaned_regression_data)
      )
    
    
    listed_results<-c(
      hotelling_results,
      lr_location_results,
      lr_correlation_results,
      regression_results
    )
    
    unlisted=unlist(listed_results,use.names = TRUE)
    unlisted
  })
  
  return(t(results))
}

# Parameters
set.seed(1)
mu <- rep(0, 4)
sigma <- diag(4)
mu_con <- rep(10, 4)
sigma_con <- matrix(0.5^abs(row(matrix(1:4, 4, 4)) - col(matrix(1:4, 4, 4))), 4, 4)
print(sigma_con)
n <- 1000
eps <- 0.05
R <- 500

# Run simulation
results <- simulate(R, n, mu, sigma, mu_con, sigma_con, eps)
results <- as.data.frame(results)

# Display column names
colnames(results)

#Extract the results of Hotelling's T2 square for all 500 simulations
hotelling_statistics=results[,grep("mean_statistic.T.2", colnames(results))]
hotelling_pvalues <- results[, grep("mean_p_value", colnames(results))]
hotelling_significance <- results[, grep("mean_significant", colnames(results))]

#plot the hotelling statistics for 500 simulations
summary(hotelling_statistics)

hotelling_long=hotelling_statistics%>%
  pivot_longer(cols=everything(), names_to="Data",values_to = "Tsquared")

hotelling_long$Data <- factor(hotelling_long$Data, 
                              levels = c("Original", setdiff(unique(hotelling_long$Data),
                                                             "Original")))
  
ggplot(hotelling_long, aes(x=Data,y=Tsquared, fill=Data))+
  geom_boxplot() +
  labs(title = "Hotelling's T^2 Statistics by Method", x = "Method", y = "T^2 Statistic") +
  theme_minimal()

#significance of the statistics
#print(hotelling_pvalues)
#print(hotelling_significance)
summary(hotelling_significance)


#Do the same for the LR statistics
LR_location_statistics=results[,grep("loc_statistic",colnames(results))]
LR_location_significance=results[,grep("loc_significant",colnames(results))]
location_long=LR_location_statistics%>%
  pivot_longer(cols=everything(), names_to="Data",values_to = "LR")

location_long$Data <- factor(location_long$Data, 
                              levels = c("Original", setdiff(unique(location_long$Data),
                                                             "Original")))

ggplot(location_long, aes(x=Data,y=LR, fill=Data))+
  geom_boxplot() +
  labs(title = "Likelihood Ratio Statistics by Method (location)", 
       x = "Method", y = "LR") +
  theme_minimal()

summary(LR_location_statistics)
summary(LR_location_significance)


LR_correlation_statistics=results[,grep("cor_statistic",colnames(results))]
LR_correlation_significance=results[,grep("cor_significant",colnames(results))]
correlation_long=LR_correlation_statistics%>%
  pivot_longer(cols=everything(), names_to="Data",values_to = "LR")

correlation_long$Data <- factor(correlation_long$Data, 
                              levels = c("Original", setdiff(unique(correlation_long$Data),
                                                             "Original")))

ggplot(correlation_long, aes(x=Data,y=LR, fill=Data))+
  geom_boxplot() +
  labs(title = "Likelihood Ratio Statistics by Method (correlation)", 
       x = "Method", y = "LR") +
  theme_minimal()

summary(LR_correlation_statistics)
summary(LR_correlation_significance)

# Extract the f-test results
f_test_statistics <- results[, grep("F2", colnames(results))]
f_test_long <- f_test_statistics %>%
  pivot_longer(cols = everything(), names_to = "Data", values_to = "F_Statistic")
print(f_test_statistics)
summary(f_test_statistics)
# Boxplot for F-test statistics
ggplot(f_test_long, aes(x = Data, y = F_Statistic, fill = Data)) +
  geom_boxplot() +
  labs(title = "F-Test Statistics by Method", x = "Method", y = "F Statistic") +
  theme_minimal()



