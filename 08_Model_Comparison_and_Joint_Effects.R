# ==============================================================================
# Script Name: 08_Model_Comparison_and_Joint_Effects.R
# Project: Obesity Phenotypes and Incident Knee Osteoarthritis (KOA)
# Purpose: Compare predictive performance (C-index/AIC) of four obesity indices 
#          and evaluate the "Dual-hit" hypothesis using a joint model (BMI + ABSI).
# ==============================================================================

library(survival)
library(dplyr)
library(gtsummary)
library(flextable)
library(officer)

cat("\n================ Pipeline: Joint Modeling and Performance Comparison ================\n")

# ------------------------------------------------------------------------------
# 1. Joint Model Analysis (BMI + ABSI)
# ------------------------------------------------------------------------------
# Rationale: BMI and ABSI are nearly orthogonal (r ≈ -0.06), allowing for a joint 
# model to test the "Dual-hit" hypothesis (mechanical vs. metabolic pathways).

cox_joint <- coxph(Surv(follow_up_time, status) ~ BMI_z + ABSI_z + 
                     age_base + gender + rural + smokev + drinkev + hibpe + diabe, 
                   data = data_cox)

cat("--- Joint Model Results (BMI and ABSI mutually adjusted) ---\n")
print(round(summary(cox_joint)$coefficients[c("BMI_z", "ABSI_z"), c("exp(coef)", "Pr(>|z|)")], 4))

# ------------------------------------------------------------------------------
# 2. Predictive Performance Comparison (C-index & AIC)
# ------------------------------------------------------------------------------

# Define models for performance comparison
models <- list(
  BMI  = coxph(Surv(follow_up_time, status) ~ BMI_z  + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = data_cox),
  WHtR = coxph(Surv(follow_up_time, status) ~ WHtR_z + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = data_cox),
  ABSI = coxph(Surv(follow_up_time, status) ~ ABSI_z + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = data_cox),
  BRI  = coxph(Surv(follow_up_time, status) ~ BRI_z  + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = data_cox)
)

# Extract metrics
performance_metrics <- data.frame(
  Model = c("Model 1: BMI", "Model 2: WHtR", "Model 3: ABSI", "Model 4: BRI"),
  C_index = sapply(models, function(m) summary(m)$concordance["C"]),
  AIC = sapply(models, function(m) extractAIC(m)[2])
)

cat("\n--- Predictive Performance Summary ---\n")
print(performance_metrics, digits = 4)

# ------------------------------------------------------------------------------
# 3. Generate Comparative Table (JAMA Style)
# ------------------------------------------------------------------------------

# Create regression tables for each model
tbls <- lapply(names(models), function(n) {
  tbl_regression(models[[n]], exponentiate = TRUE) %>% 
    modify_header(label = paste0("**Model: ", n, "**"))
})

# Merge and format as Flextable
ft_final <- tbl_merge(tbls = tbls, tab_spanner = paste0("**Model ", 1:4, ": ", names(models), "**")) %>%
  as_flex_table() %>%
  add_header_lines("Table: Comparative Predictive Performance of Obesity Indices for Incident KOA") %>%
  # Append performance metrics rows
  add_body_row(values = list(label = "C-index", 
                             stat_0_1 = as.character(round(performance_metrics$C_index[1], 3)),
                             stat_0_2 = as.character(round(performance_metrics$C_index[2], 3)),
                             stat_0_3 = as.character(round(performance_metrics$C_index[3], 3)),
                             stat_0_4 = as.character(round(performance_metrics$C_index[4], 3))), top = FALSE) %>%
  add_body_row(values = list(label = "AIC", 
                             stat_0_1 = as.character(round(performance_metrics$AIC[1], 0)),
                             stat_0_2 = as.character(round(performance_metrics$AIC[2], 0)),
                             stat_0_3 = as.character(round(performance_metrics$AIC[3], 0)),
                             stat_0_4 = as.character(round(performance_metrics$AIC[4], 0))), top = FALSE)

# Aesthetic styling per JAMA requirements
ft_final <- ft_final %>%
  fontsize(size = 9, part = "all") %>% 
  font(fontname = "Times New Roman", part = "all") %>% 
  hline_top(part = "header", border = fp_border(width = 2)) %>% 
  hline_bottom(part = "body", border = fp_border(width = 2))

save_as_docx(ft_final, path = "Table_Comparative_Performance.docx")

cat("\n================ Pipeline Complete ================\n")