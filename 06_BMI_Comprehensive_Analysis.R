# ==============================================================================
# Script Name: 06_BMI_Comprehensive_Analysis.R
# Project: Obesity Phenotypes and Incident Knee Osteoarthritis (KOA)
# Purpose: BMI-specific analysis including categorization, survival curves,
#          non-linear RCS modeling, and subgroup heterogeneity assessment.
# ==============================================================================

library(survival)
library(survminer)
library(dplyr)
library(rms)
library(ggplot2)
library(ggsci)

cat("\n================ Pipeline: Comprehensive Analysis of BMI ================\n")

# ------------------------------------------------------------------------------
# 1. Categorical Analysis (KM Plot)
# ------------------------------------------------------------------------------

# Define BMI categories based on Chinese guidelines
data_cox <- data_cox %>%
  mutate(
    bmi_group = cut(
      bmi_base, 
      breaks = c(-Inf, 18.5, 24.0, 28.0, Inf), 
      labels = c("Underweight (<18.5)", "Normal (18.5-23.9)", "Overweight (24.0-27.9)", "Obese (≥28.0)"),
      right = FALSE
    )
  ) %>%
  filter(!is.na(bmi_group))

# Fit KM model
km_fit_bmi <- survfit(Surv(follow_up_time, status) ~ bmi_group, data = data_cox)

# Plot Cumulative Incidence
p_km_bmi <- ggsurvplot(
  km_fit_bmi, data = data_cox, fun = "event",
  palette = pal_jama()(4), size = 1,
  conf.int = TRUE, pval = TRUE, pval.method = TRUE,
  risk.table = TRUE, risk.table.col = "strata",
  legend.title = "BMI Categories",
  xlab = "Follow-up Time (Years)", ylab = "Cumulative Incidence of KOA",
  title = "Cumulative Incidence of KOA by BMI Categories",
  ggtheme = theme_classic()
)

ggsave("Figure_BMI_KM_Incidence.pdf", print(p_km_bmi), width = 8, height = 7)

# ------------------------------------------------------------------------------
# 2. Non-linear Dose-Response Analysis (RCS)
# ------------------------------------------------------------------------------

dd <- datadist(data_cox)
bmi_ref <- 24 
dd$limits["Adjust to", "bmi_base"] <- bmi_ref 
options(datadist = "dd")

fit_rcs <- cph(Surv(follow_up_time, status) ~ rcs(bmi_base, 4) + age_base + gender + rural + smokev + drinkev + hibpe + diabe, 
               data = data_cox, x = TRUE, y = TRUE)

rcs_data <- as.data.frame(Predict(fit_rcs, bmi_base, ref.zero = TRUE, fun = exp))

p_rcs <- ggplot(rcs_data, aes(x = bmi_base, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#374E55FF", alpha = 0.15) +
  geom_line(color = "#B24745FF", size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = bmi_ref, linetype = "dotted", color = "#00468BFF") +
  scale_y_continuous(trans = "log10", breaks = c(0.5, 1, 1.5, 2, 3, 5)) + 
  labs(title = "Dose-Response Relationship: BMI and KOA Risk",
       subtitle = paste0("Reference: BMI = ", bmi_ref, " kg/m²"),
       x = "Body Mass Index (BMI, kg/m²)", y = "Adjusted Hazard Ratio (aHR)") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", size = 14), axis.title = element_text(face = "bold"))

ggsave("Figure_BMI_RCS.tiff", p_rcs, width = 7, height = 5.5, dpi = 350, compression = "lzw")

# ------------------------------------------------------------------------------
# 3. Subgroup Heterogeneity Analysis (Forest Plot)
# ------------------------------------------------------------------------------

# Define subgroup function
run_subgroup <- function(data, subgroup_col, subgroup_name, level_name) {
  sub_data <- data[data[[subgroup_col]] == level_name, ]
  fit <- coxph(Surv(follow_up_time, status) ~ BMI_z + age_base + gender + rural + smokev + drinkev + hibpe + diabe, data = sub_data)
  res <- summary(fit)
  return(data.frame(Subgroup = subgroup_name, Level = level_name, N = nrow(sub_data), 
                    HR = res$coefficients["BMI_z", "exp(coef)"], 
                    LCI = res$conf.int["BMI_z", "lower .95"], 
                    UCI = res$conf.int["BMI_z", "upper .95"], P_val = res$coefficients["BMI_z", "Pr(>|z|)"]))
}

# Run subgroup analysis for Sex, Age, Hypertension, Diabetes
# (Assuming data_sub is prepared as in previous steps)
# ... [Insert subgroup result extraction logic here] ...

# Final Forest Plot rendering
p_forest_bmi <- ggplot(df_forest, aes(x = HR, y = Label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_pointrange(aes(xmin = LCI, xmax = UCI, color = Subgroup), size = 0.7, shape = 15) +
  geom_text(aes(x = 1.35, label = sprintf("%.2f (%.2f - %.2f)", HR, LCI, UCI)), size = 3.5, hjust = 0) +
  facet_grid(Subgroup ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_color_jama() + 
  scale_x_continuous(limits = c(0.75, 1.55)) +
  labs(title = "Subgroup Analysis of BMI on Incident KOA Risk",
       subtitle = "Models adjusted for age, sex, rural status, smoking, alcohol, hypertension, and diabetes.",
       x = "Adjusted Hazard Ratio (95% CI) per 1-SD Increase", y = "") +
  theme_bw() +
  annotate("text", x = 1.35, y = Inf, label = "HR (95% CI)", vjust = -1, hjust = 0, size = 4, fontface = "bold")

ggsave("Figure_BMI_Subgroup_Forest.tiff", p_forest_bmi, width = 10, height = 7, dpi = 350, compression = "lzw")

cat("\n================ Pipeline Complete ================\n")