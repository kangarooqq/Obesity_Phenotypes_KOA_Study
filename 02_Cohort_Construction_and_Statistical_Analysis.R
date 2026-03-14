# ==============================================================================
# Script Name: 02_Cohort_Construction_and_Statistical_Analysis.R
# Project: Obesity Phenotypes and Incident Knee Osteoarthritis (KOA)
# Purpose: Construction of incident cohort, survival analysis, and visualizations
# ==============================================================================

# 1. Load required libraries
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(gtsummary)
library(flextable)
library(forestplot)
library(broom)
library(ggsci)

# ------------------------------------------------------------------------------
# 2. Cohort Construction and Data Preparation
# ------------------------------------------------------------------------------

cat("\n================ Pipeline: Constructing Incident Cohort ================\n")

# Step 1: Pre-processing and cleaning
charls_master <- charls_data %>%
  mutate(
    mheight = ifelse(mheight < 1.0 | mheight > 2.2, NA, mheight),
    mwaist  = ifelse(mwaist < 40 | mwaist > 150, NA, mwaist),
    bmi     = ifelse(bmi < 10 | bmi > 60, NA, bmi),
    bl_crp  = ifelse(bl_crp > 10, NA, bl_crp),
    year    = as.numeric(iwy),
    koa     = case_when(
      arthre == 1 & da042s12 == 1 ~ 1,
      arthre == 0 | da042s12 == 0 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(age >= 45) %>%
  select(ID, wave, year, koa, bmi, mheight, mwaist, bl_crp, tyg_bmi, 
         age, gender, rural, edu, smokev, drinkev, hibpe, diabe)

total_individuals_age_ok <- length(unique(charls_master$ID))

# Step 2: Define incident cohort (excluding baseline prevalent cases)
cohort_long <- charls_master %>%
  arrange(ID, year) %>%
  group_by(ID) %>%
  mutate(
    baseline_koa = first(koa),
    total_visits = n()
  ) %>%
  filter(baseline_koa == 0, total_visits > 1) %>%
  ungroup()

incident_cohort_ids <- length(unique(cohort_long$ID))
excluded_prevalent_or_single <- total_individuals_age_ok - incident_cohort_ids

# Step 3: Create analytic dataset (data_cox) and calculate obesity indices
data_cox <- cohort_long %>%
  group_by(ID) %>%
  summarise(
    status = ifelse(any(koa == 1, na.rm = TRUE), 1, 0),
    baseline_year  = first(year),
    end_year       = if_else(status == 1, year[which(koa == 1)[1]], last(year)),
    follow_up_time = end_year - baseline_year,
    bmi_base       = first(bmi),
    waist_base     = first(mwaist),
    mheight_base   = first(mheight),
    crp_base       = first(bl_crp),
    tyg_bmi_base   = first(tyg_bmi),
    age_base       = first(age),
    gender         = as.factor(first(gender)),
    rural          = as.factor(first(rural)),
    smokev         = as.factor(first(smokev)),
    drinkev        = as.factor(first(drinkev)),
    hibpe          = as.factor(first(hibpe)),
    diabe          = as.factor(first(diabe))
  ) %>%
  filter(follow_up_time > 0) %>%
  ungroup() %>%
  mutate(
    mheight_m = ifelse(mheight_base > 10, mheight_base / 100, mheight_base),
    mwaist_m  = waist_base / 100,
    WHtR      = mwaist_m / mheight_m,
    ABSI      = mwaist_m / ((bmi_base^(2/3)) * (mheight_m^(1/2))),
    BRI       = 364.2 - 365.5 * sqrt(1 - (mwaist_m / (pi * mheight_m))^2)
  ) %>%
  filter(is.finite(WHtR), is.finite(ABSI), is.finite(BRI))

final_cox_n <- nrow(data_cox)
excluded_missing_covariates <- incident_cohort_ids - final_cox_n

cat("Cohort construction complete. Final N:", final_cox_n, "\n")

# ------------------------------------------------------------------------------
# 3. Descriptive Statistics (Table 1)
# ------------------------------------------------------------------------------

table1_data <- data_cox %>%
  mutate(
    bmi_category = cut(
      bmi_base, 
      breaks = c(-Inf, 18.5, 24.0, 28.0, Inf), 
      labels = c("Underweight (<18.5)", "Normal (18.5-23.9)", "Overweight (24.0-27.9)", "Obese (≥28.0)"),
      right = FALSE
    ),
    gender = factor(gender, levels = c(0, 1), labels = c("Female", "Male"))
  ) %>%
  select(age_base, gender, bmi_base, waist_base, crp_base, bmi_category) %>%
  filter(!is.na(bmi_category))

table1 <- table1_data %>%
  tbl_summary(
    by = bmi_category,
    label = list(
      age_base ~ "Age, years",
      gender ~ "Sex",
      bmi_base ~ "Body Mass Index, kg/m²",
      waist_base ~ "Waist circumference, cm",
      crp_base ~ "C-Reactive Protein, mg/L"
    ),
    statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
    missing = "no"
  ) %>%
  add_p(test = list(all_continuous() ~ "aov", all_categorical() ~ "chisq.test")) %>% 
  add_overall() %>%
  bold_labels() %>%
  modify_header(label = "**Characteristic**") %>%
  modify_caption("**Table 1. Baseline Characteristics of the Incident KOA Cohort (n=14,637)**")

# Save table
table1 %>% as_flex_table() %>% flextable::save_as_docx(path = "Table1_Baseline_Characteristics.docx")

# ------------------------------------------------------------------------------
# 4. Survival Analysis and Visualization
# ------------------------------------------------------------------------------

# 4.1 Cox Proportional Hazards Model
cox_final <- coxph(Surv(follow_up_time, status) ~ bmi_base + age_base + gender + rural + hibpe + diabe, data = data_cox)

# 4.2 Forest Plot
tidy_data <- tidy(cox_final, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    res_label = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    term = case_when(
      term == "bmi_base" ~ "Baseline BMI",
      term == "age_base" ~ "Age (years)",
      term == "gender1"  ~ "Sex (Male vs Female)",
      term == "rural0"   ~ "Residence (Urban vs Rural)",
      term == "hibpe1"   ~ "Hypertension (Yes vs No)",
      term == "diabe0"   ~ "Diabetes (No vs Yes)",
      TRUE ~ term
    )
  )

table_text <- list(c("Variables", tidy_data$term), c("HR (95% CI)", tidy_data$res_label))

p_forest <- forestplot(
  labeltext = table_text,
  mean = c(NA, tidy_data$estimate),
  lower = c(NA, tidy_data$conf.low),
  upper = c(NA, tidy_data$conf.high),
  is.summary = c(TRUE, rep(FALSE, nrow(tidy_data))),
  zero = 1,
  col = fpColors(box = "#374E55FF", line = "#374E55FF", summary = "#374E55FF"),
  xlab = "Hazard Ratio (95% CI)",
  title = "Hazard Ratios for Incident KOA"
)

# 4.3 Cumulative Incidence Plots (Example: BMI Categories)
data_cox <- data_cox %>% mutate(
  bmi_group = cut(bmi_base, breaks = c(-Inf, 18.5, 24.0, 28.0, Inf), 
                  labels = c("Underweight (<18.5)", "Normal (18.5-23.9)", "Overweight (24.0-27.9)", "Obese (≥28.0)"), right = FALSE)
)
km_fit_bmi <- survfit(Surv(follow_up_time, status) ~ bmi_group, data = data_cox)

p_km_bmi <- ggsurvplot(
  km_fit_bmi, data = data_cox, fun = "event",
  palette = c("#374E55FF", "#8FADBCFF", "#DF8F44FF", "#B24745FF"),
  pval = TRUE, risk.table = TRUE,
  xlab = "Follow-up Time (Years)", ylab = "Cumulative Incidence of KOA",
  ggtheme = theme_classic()
)

# [Note: Repeat the ggsurvplot logic for WHtR, ABSI, and BRI as per your original code modules]

cat("\n================ Pipeline Complete ================\n")