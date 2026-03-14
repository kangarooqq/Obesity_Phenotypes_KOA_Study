# ==============================================================================
# Script Name: 05_ABSI_Comprehensive_Analysis.R
# Project: Obesity Phenotypes and Incident Knee Osteoarthritis (KOA)
# Purpose: Comprehensive analysis for A Body Shape Index (ABSI) including 
#          quartile-based modeling, non-linear RCS curve, and subgroup heterogeneity.
# ==============================================================================

library(survival)
library(survminer)
library(dplyr)
library(rms)
library(ggplot2)
library(ggsci)
library(grid)

cat("\n================ Pipeline: Comprehensive Analysis of ABSI ================\n")

# ------------------------------------------------------------------------------
# 1. Quartile-based Analysis
# ------------------------------------------------------------------------------

# Define ABSI quartiles for trend analysis
data_cox <- data_cox %>%
  mutate(
    ABSI_Q_num = ntile(ABSI, 4),
    ABSI_Q = factor(ABSI_Q_num, levels = c(1, 2, 3, 4),
                    labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)"))
  )

# Fit categorical and trend models
cox_cat <- coxph(Surv(follow_up_time, status) ~ ABSI_Q + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = data_cox)
cox_trend <- coxph(Surv(follow_up_time, status) ~ ABSI_Q_num + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = data_cox)

cat("\n--- Cox Regression Results (ABSI Quartiles) ---\n")
print(summary(cox_cat))
cat("P for Trend:", summary(cox_trend)$coefficients["ABSI_Q_num", "Pr(>|z|)"], "\n")

# ------------------------------------------------------------------------------
# 2. Dose-Response Relationship (RCS)
# ------------------------------------------------------------------------------

# Use median as reference for RCS, as ABSI lacks a universal clinical threshold
absi_median <- median(data_cox$ABSI, na.rm = TRUE)
dd <- datadist(data_cox)
dd$limits["Adjust to", "ABSI"] <- absi_median
options(datadist = "dd")

fit_rcs <- cph(Surv(follow_up_time, status) ~ rcs(ABSI, 4) + age_base + gender + rural + smokev + drinkev + hibpe + diabe, 
               data = data_cox, x = TRUE, y = TRUE)

rcs_data <- as.data.frame(Predict(fit_rcs, ABSI, ref.zero = TRUE, fun = exp))

p_rcs_absi <- ggplot(rcs_data, aes(x = ABSI, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#374E55FF", alpha = 0.15) +
  geom_line(color = "#B24745FF", size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = absi_median, linetype = "dotted", color = "#00468BFF") +
  scale_y_continuous(trans = "log10", breaks = c(0.5, 1, 1.5, 2, 3)) +
  labs(title = "Dose-Response Relationship: ABSI and KOA Risk",
       subtitle = paste0("Reference: Median ABSI (", round(absi_median, 4), ")"),
       x = "A Body Shape Index (ABSI)", y = "Adjusted Hazard Ratio (aHR)") +
  theme_classic()

ggsave("Figure_ABSI_RCS.tiff", p_rcs_absi, dpi = 350, width = 7, height = 5)

# ------------------------------------------------------------------------------
# 3. Subgroup Analysis & Interaction
# ------------------------------------------------------------------------------

# Define subgroup variables
data_sub <- data_cox %>% mutate(
  Sub_Sex = ifelse(as.numeric(as.character(gender)) == 1, "Male", "Female"),
  Sub_Age = ifelse(age_base >= 60, ">= 60 years", "< 60 years"),
  Sub_Hibpe = ifelse(as.numeric(as.character(hibpe)) == 1, "Yes", "No"),
  Sub_Diabe = ifelse(as.numeric(as.character(diabe)) == 1, "Yes", "No")
)

# Automated subgroup function (using ABSI_z for standardization)
run_subgroup <- function(data, subgroup_col, subgroup_name, level_name) {
  sub_data <- data[data[[subgroup_col]] == level_name, ]
  fit <- coxph(Surv(follow_up_time, status) ~ ABSI_z + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = sub_data)
  res <- summary(fit)
  return(data.frame(Subgroup = subgroup_name, Level = level_name, N = nrow(sub_data), 
                    HR = res$coefficients["ABSI_z", "exp(coef)"], 
                    LCI = res$conf.int["ABSI_z", "lower .95"], 
                    UCI = res$conf.int["ABSI_z", "upper .95"], 
                    P_val = res$coefficients["ABSI_z", "Pr(>|z|)"]))
}

# (Compute list of results and plot...)
# [Visualisation code follows the forest plot pattern established in previous steps]

# ------------------------------------------------------------------------------
# 4. Cumulative Incidence Plot (KM)
# ------------------------------------------------------------------------------

km_fit_absi <- survfit(Surv(follow_up_time, status) ~ ABSI_Q, data = data_cox)
p_km_absi <- ggsurvplot(km_fit_absi, data = data_cox, fun = "event", palette = pal_jama()(4),
                        risk.table = TRUE, pval = TRUE, xlab = "Follow-up Time (Years)",
                        ylab = "Cumulative Incidence of KOA", ggtheme = theme_bw())

ggsave("Figure_ABSI_KM.tiff", print(p_km_absi), dpi = 350, width = 9, height = 8)

cat("\n================ Pipeline Complete ================\n")