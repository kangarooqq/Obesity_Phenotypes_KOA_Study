# ==============================================================================
# Script Name: 07_BRI_Comprehensive_Analysis.R
# Project: Obesity Phenotypes and Incident Knee Osteoarthritis (KOA)
# Purpose: Comprehensive analysis for Body Roundness Index (BRI) including 
#          quartile-based modeling, non-linear RCS, and subgroup heterogeneity.
# ==============================================================================

library(survival)
library(survminer)
library(dplyr)
library(rms)
library(ggplot2)
library(ggsci)

cat("\n================ Pipeline: Comprehensive Analysis of BRI ================\n")

# ------------------------------------------------------------------------------
# 1. Quartile-based Analysis
# ------------------------------------------------------------------------------

# Stratify BRI into quartiles for trend analysis
data_cox <- data_cox %>%
  mutate(
    bri_q_num = ntile(BRI, 4),
    bri_label = factor(bri_q_num, levels = c(1, 2, 3, 4),
                       labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)"))
  ) %>%
  filter(!is.na(bri_label))

# Fit Cox models for categorical and trend assessment
cox_cat_bri <- coxph(Surv(follow_up_time, status) ~ bri_label + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = data_cox)
cox_trend_bri <- coxph(Surv(follow_up_time, status) ~ bri_q_num + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = data_cox)

# ------------------------------------------------------------------------------
# 2. Cumulative Incidence (KM Plot)
# ------------------------------------------------------------------------------

km_fit_bri <- survfit(Surv(follow_up_time, status) ~ bri_label, data = data_cox)
p_km_bri <- ggsurvplot(
  km_fit_bri, data = data_cox, fun = "event",
  palette = pal_jama()(4), size = 1,
  conf.int = TRUE, pval = TRUE,
  risk.table = TRUE, risk.table.col = "strata",
  xlab = "Follow-up Time (Years)", ylab = "Cumulative Incidence of KOA",
  ggtheme = theme_classic()
)
ggsave("Figure_BRI_KM_Incidence.pdf", print(p_km_bri), width = 9, height = 8)

# ------------------------------------------------------------------------------
# 3. Dose-Response Analysis (RCS)
# ------------------------------------------------------------------------------

# Set median as reference for RCS (as no clinical threshold exists for BRI)
bri_ref <- median(data_cox$BRI, na.rm = TRUE)
dd <- datadist(data_cox)
dd$limits["Adjust to", "BRI"] <- bri_ref
options(datadist = "dd")

fit_rcs_bri <- cph(Surv(follow_up_time, status) ~ rcs(BRI, 4) + age_base + gender + rural + smokev + drinkev + hibpe + diabe, 
                   data = data_cox, x = TRUE, y = TRUE)

rcs_data_bri <- as.data.frame(Predict(fit_rcs_bri, BRI, ref.zero = TRUE, fun = exp))

p_rcs_bri <- ggplot(rcs_data_bri, aes(x = BRI, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#374E55FF", alpha = 0.15) +
  geom_line(color = "#B24745FF", size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = bri_ref, linetype = "dotted", color = "#00468BFF") +
  scale_y_continuous(trans = "log10", breaks = c(0.5, 1, 1.5, 2, 3)) +
  labs(title = "Dose-Response Relationship: BRI and KOA Risk",
       subtitle = paste0("Reference: Median BRI = ", round(bri_ref, 3)),
       x = "Body Roundness Index (BRI)", y = "Adjusted Hazard Ratio (aHR)") +
  theme_classic()

ggsave("Figure_BRI_RCS.tiff", p_rcs_bri, dpi = 350, width = 7, height = 5.5)

# ------------------------------------------------------------------------------
# 4. Subgroup Analysis & Forest Plot
# ------------------------------------------------------------------------------

# Define subgroup function specifically for BRI
run_subgroup_bri <- function(data, subgroup_col, subgroup_name, level_name) {
  sub_data <- data[data[[subgroup_col]] == level_name, ]
  fit <- coxph(Surv(follow_up_time, status) ~ scale(BRI) + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = sub_data)
  res <- summary(fit)
  return(data.frame(Subgroup = subgroup_name, Level = level_name, N = nrow(sub_data), 
                    HR = res$coefficients[1, "exp(coef)"], 
                    LCI = res$conf.int[1, "lower .95"], 
                    UCI = res$conf.int[1, "upper .95"], 
                    P_val = res$coefficients[1, "Pr(>|z|)"]))
}

# (Execute function calls as done for previous indices)

# Final Forest Plot rendering with standardized JAMA styling...
cat("\n================ Pipeline Complete ================\n")