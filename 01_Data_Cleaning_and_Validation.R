# ==============================================================================
# Script Name: 01_Data_Cleaning_and_Validation.R
# Project: Obesity Phenotypes and Incident Knee Osteoarthritis
# Purpose: Data import, cleaning, outlier management, and KOA case definition
# ==============================================================================

# 1. Load required packages
library(haven)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggsci)

cat("\n================ Pipeline Start: Data Cleaning & Definition ================\n")

# 2. Import Data

file_path <- "charls.dta" # The consolidated raw data from five waves of CHARLS were used for analysis
charls_data <- read_dta(file_path)

# 3. Step 1: Outlier Cleaning and Outcome Definition
cat("Cleaning extreme values and defining KOA cases...\n")

charls_master <- charls_data %>%
  # Recode physiologically implausible anthropometric and lab values to NA
  mutate(
    mheight = ifelse(mheight < 1.0 | mheight > 2.2, NA, mheight),
    mwaist  = ifelse(mwaist < 40 | mwaist > 150, NA, mwaist),
    bmi     = ifelse(bmi < 10 | bmi > 60, NA, bmi),
    bl_crp  = ifelse(bl_crp > 10, NA, bl_crp), # Exclude acute inflammation
    year    = as.numeric(iwy),
    # KOA defined as physician-diagnosed arthritis combined with chronic knee pain
    koa     = case_when(
      arthre == 1 & da042s12 == 1 ~ 1,
      arthre == 0 | da042s12 == 0 ~ 0, 
      TRUE ~ NA_real_
    )
  ) %>%
  # Restrict to middle-aged and older population
  filter(age >= 45) %>%
  # Select essential variables for longitudinal tracking
  select(ID, wave, year, koa, arthre, da042s12, 
         bmi, mheight, mwaist, bl_crp, tyg_bmi, 
         age, gender, rural, edu, smokev, drinkev, hibpe, diabe)

cat("Observations after cleaning:", nrow(charls_master), "\n")

# 4. Data Quality Visualization (JAMA-style)
theme_jama_pub <- function() {
  theme_bw() + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      text = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      legend.position = "bottom"
    )
}

jama_cols <- pal_jama()(7)

# Generate Quality Control Plots
p1 <- ggplot(charls_master, aes(x = row_number(), y = bmi)) +
  geom_point(color = jama_cols[3], alpha = 0.5, size = 1.5) +
  theme_jama_pub() +
  labs(title = "A. BMI Index Plot (Cleaned)", x = "Sample Index", y = "BMI (kg/m²)")

p2 <- ggplot(charls_master, aes(x = mheight, y = mweight)) +
  geom_point(color = jama_cols[1], alpha = 0.4, size = 1.5) +
  theme_jama_pub() +
  labs(title = "B. Height vs. Weight Scatterplot", x = "Height (m)", y = "Weight (kg)")

p3 <- ggplot(charls_master, aes(x = row_number(), y = bl_crp)) +
  geom_point(color = jama_cols[4], alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "black", size = 0.8) +
  theme_jama_pub() +
  labs(title = "C. C-Reactive Protein (CRP) Distribution", x = "Sample Index", y = "CRP (mg/L)")

# Combined plot for quality control
combined_plot <- p1 / p2 / p3 + plot_annotation(tag_levels = '1')
print(combined_plot)

# Save high-quality publication files
ggsave("Figure_QC_Data_Cleaning.pdf", combined_plot, width = 10, height = 14)
ggsave("Figure_QC_Data_Cleaning.tiff", combined_plot, width = 10, height = 14, dpi = 350, compression = "lzw")

cat("\n================ Pipeline Step 1 Complete ================\n")