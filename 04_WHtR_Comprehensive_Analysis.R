# ==============================================================================
# Script Name: 04_WHtR_Advanced_Epidemiological_Analysis.R
# Project: Obesity Phenotypes and Incident Knee Osteoarthritis (KOA)
# Purpose: Perform quartile-based trend analysis, RCS modeling, and subgroup analysis
# ==============================================================================

library(survival)
library(survminer)
library(dplyr)
library(rms)
library(ggplot2)
library(ggsci)
library(grid)

cat("\n================ Pipeline: Dose-Response and Subgroup Analysis ================\n")

# ------------------------------------------------------------------------------
# 1. Quartile-based Analysis & Trend Test
# ------------------------------------------------------------------------------

# Create quartiles for WHtR to assess dose-response patterns
data_cox <- data_cox %>%
  mutate(
    WHtR_Q_num = ntile(WHtR, 4),
    WHtR_Q = factor(WHtR_Q_num, levels = c(1, 2, 3, 4),
                    labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)"))
  )

# Survival analysis: Categorical model vs Trend model
cox_cat <- coxph(Surv(follow_up_time, status) ~ WHtR_Q + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = data_cox)
cox_trend <- coxph(Surv(follow_up_time, status) ~ WHtR_Q_num + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = data_cox)

# Print tabular results
print(summary(cox_cat))
cat("P for Trend:", summary(cox_trend)$coefficients["WHtR_Q_num", "Pr(>|z|)"], "\n")

# ------------------------------------------------------------------------------
# 2. Dose-Response Analysis (Restricted Cubic Splines - RCS)
# ------------------------------------------------------------------------------

# Define data distribution for rms package
dd <- datadist(data_cox)
dd$limits["Adjust to", "WHtR"] <- 0.5 
options(datadist = "dd")

# Fit RCS model with 4 knots
fit_rcs <- cph(Surv(follow_up_time, status) ~ rcs(WHtR, 4) + age_base + gender + rural + smokev + drinkev + hibpe + diabe, 
               data = data_cox, x = TRUE, y = TRUE)

# Prepare data for ggplot
rcs_data <- as.data.frame(Predict(fit_rcs, WHtR, ref.zero = TRUE, fun = exp))

# Plot RCS curve
p_rcs <- ggplot(rcs_data, aes(x = WHtR, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#374E55FF", alpha = 0.15) +
  geom_line(color = "#B24745FF", size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", size = 0.7) +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "#00468BFF", size = 0.8) +
  scale_y_continuous(trans = "log10", breaks = c(0.5, 1, 1.5, 2, 3)) +
  labs(title = "Dose-Response Relationship: WHtR and KOA Risk",
       subtitle = "Reference: WHtR = 0.5",
       x = "Waist-to-Height Ratio (WHtR)", y = "Adjusted Hazard Ratio (aHR)") +
  theme_classic()

ggsave("Figure_RCS_DoseResponse.tiff", p_rcs, dpi = 350, width = 8, height = 6)

# ------------------------------------------------------------------------------
# 3. Subgroup Analysis & P-interaction
# ------------------------------------------------------------------------------

# Extract HRs for subgroups
run_subgroup <- function(data, subgroup_col, level_name) {
  sub_data <- data[data[[subgroup_col]] == level_name, ]
  fit <- coxph(Surv(follow_up_time, status) ~ WHtR_z + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = sub_data)
  res <- summary(fit)
  return(data.frame(Level = level_name, HR = res$coefficients["WHtR_z", "exp(coef)"], 
                    LCI = res$conf.int["WHtR_z", "lower .95"], UCI = res$conf.int["WHtR_z", "upper .95"], 
                    P_val = res$coefficients["WHtR_z", "Pr(>|z|)"], N = nrow(sub_data)))
}

# (Simplified extraction for brevity; ensure sub-groups are computed as in your logic)
# [Function calls for Sex, Age, Hibpe, Diabe...]

# Calculate Interaction P-values
extract_p_int <- function(model) summary(model)$coefficients[grep(":", rownames(summary(model)$coefficients))[1], "Pr(>|z|)"]

# Compute interaction models
p_sex <- extract_p_int(coxph(Surv(follow_up_time, status) ~ WHtR_z * gender + age_base + rural + smokev + drinkev + hibpe + diabe, data = data_cox))
cat("P for interaction (Sex):", p_sex, "\n")

# ------------------------------------------------------------------------------
# 4. Final Visualization: Cumulative Incidence (KM Plot)
# ------------------------------------------------------------------------------

km_fit <- survfit(Surv(follow_up_time, status) ~ WHtR_Q, data = data_cox)
p_km <- ggsurvplot(km_fit, data = data_cox, fun = "event", palette = pal_jama()(4), 
                   risk.table = TRUE, pval = TRUE, xlab = "Follow-up Time (Years)", 
                   ylab = "Cumulative Incidence of KOA", ggtheme = theme_bw())

ggsave("Figure_KM_Incidence.tiff", print(p_km), dpi = 350, width = 9, height = 8)

cat("\n================ Pipeline Complete ================\n")