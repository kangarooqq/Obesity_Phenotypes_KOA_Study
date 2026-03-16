# ==============================================================================
# Script Name: 12_NRI_IDI_Analysis.R
# Project: Obesity Phenotypes and Incident Knee Osteoarthritis (KOA)
# Purpose: Calculate continuous Net Reclassification Improvement (NRI) and 
#          Integrated Discrimination Improvement (IDI) for model comparison.
# ==============================================================================

library(survival)
library(rms)
library(dplyr)

cat("\n================ Pipeline: NRI & IDI Analysis ================\n")

# 1. Define base model and new model
# f_base: Model including BMI + baseline covariates
# f_new:  Model including BMI + WHtR + baseline covariates
f_base <- fit0 
f_new  <- fit1 

# 2. Calculate 5-year predicted risks
# Using the survival function derived from 'rms' models
S_base <- Survival(f_base)
S_new  <- Survival(f_new)

pred_base <- 1 - S_base(5, predict(f_base))
pred_new  <- 1 - S_new(5, predict(f_new))

# 3. Define event indices
ev_mask  <- (data_cox_nomo2$follow_up_time <= 5 & data_cox_nomo2$status == 1)
nev_mask <- !ev_mask

# 4. Bootstrap Estimation for Standard Errors and P-values
set.seed(2025)
B <- 1000
n <- nrow(data_cox_nomo2)
cNRI_boot <- numeric(B)
IDI_boot  <- numeric(B)

for (b in 1:B) {
  idx <- sample(1:n, n, replace = TRUE)
  pb <- pred_base[idx]
  pn <- pred_new[idx]
  
  # Define event status within the bootstrap sample
  ev_b  <- ev_mask[idx]
  nev_b <- nev_mask[idx]
  
  # Calculate cNRI (Continuous NRI)
  # Event group: proportion of individuals with increased risk
  cNRI_event_b <- mean(pn[ev_b] > pb[ev_b]) - mean(pn[ev_b] < pb[ev_b])
  # Non-event group: proportion of individuals with decreased risk
  cNRI_nonevent_b <- mean(pb[nev_b] > pn[nev_b]) - mean(pb[nev_b] < pn[nev_b])
  
  cNRI_boot[b] <- cNRI_event_b + cNRI_nonevent_b
  
  # Calculate IDI
  IDI_event_b <- mean(pn[ev_b]) - mean(pb[ev_b])
  IDI_nonevent_b <- mean(pn[nev_b]) - mean(pb[nev_b])
  IDI_boot[b] <- IDI_event_b - IDI_nonevent_b
}

# 5. Point Estimates and Statistical Inference
cNRI_total <- mean(cNRI_boot)
IDI_total  <- mean(IDI_boot)

cNRI_se <- sd(cNRI_boot)
IDI_se  <- sd(IDI_boot)

# Z-statistics and P-values
cNRI_p <- 2 * (1 - pnorm(abs(cNRI_total / cNRI_se)))
IDI_p  <- 2 * (1 - pnorm(abs(IDI_total / IDI_se)))

# 6. Report Results
cat("\n===== 5-Year Continuous NRI / IDI Results =====\n")
cat("Total cNRI: ", round(cNRI_total, 4), " (p =", signif(cNRI_p, 3), ")\n")
cat("Total IDI:  ", round(IDI_total, 4), " (p =", signif(IDI_p, 3), ")\n")

cat("\n================ NRI_IDI_Analysis Complete ================\n")