# ==============================================================================
# Script Name: 03_Obesity_Indices_Comparison_and_Visualization.R
# Project: Obesity Phenotypes and Incident Knee Osteoarthritis (KOA)
# Purpose: Compare predictive performance of multiple obesity indices (BMI, WHtR, ABSI, BRI)
# ==============================================================================

library(survival)
library(dplyr)
library(forestplot)
library(broom)
library(grid)

cat("\n================ Pipeline: Comparison of Obesity Indices ================\n")

# 1. Standardize obesity indices (Z-score transformation)
# Z-score enables comparison across different scales (kg/m2 vs ratio/index)
data_cox <- data_cox %>% mutate(
  BMI_z  = as.vector(scale(bmi_base)),
  WHtR_z = as.vector(scale(WHtR)),
  ABSI_z = as.vector(scale(ABSI)),
  BRI_z  = as.vector(scale(BRI))
)

# Ensure categorical covariates are properly encoded
data_cox <- data_cox %>% mutate(
  smokev = as.factor(smokev),
  drinkev = as.factor(drinkev),
  hibpe = as.factor(hibpe),
  diabe = as.factor(diabe)
)

# 2. Fit Cox Proportional Hazards Models
# Formula includes standard adjustments for potential confounders
formula_base <- Surv(follow_up_time, status) ~ age_base + gender + rural + smokev + drinkev + hibpe + diabe

cox_bmi  <- coxph(update(formula_base, . ~ . + BMI_z),  data = data_cox)
cox_whtr <- coxph(update(formula_base, . ~ . + WHtR_z), data = data_cox)
cox_absi <- coxph(update(formula_base, . ~ . + ABSI_z), data = data_cox)
cox_bri  <- coxph(update(formula_base, . ~ . + BRI_z),  data = data_cox)

# Print results
cat("\n--- Hazard Ratios per 1-SD Increase for each Index ---\n")
print(summary(cox_bmi)$coefficients["BMI_z", c("exp(coef)", "Pr(>|z|)")])
print(summary(cox_whtr)$coefficients["WHtR_z", c("exp(coef)", "Pr(>|z|)")])
print(summary(cox_absi)$coefficients["ABSI_z", c("exp(coef)", "Pr(>|z|)")])
print(summary(cox_bri)$coefficients["BRI_z", c("exp(coef)", "Pr(>|z|)")])

# 3. Model Performance Comparison (C-index)
cat("\n--- Comparison of Model Predictive Performance (C-index) ---\n")
cat("BMI Model C-index:", summary(cox_bmi)$concordance["C"], "\n")
cat("WHtR Model C-index:", summary(cox_whtr)$concordance["C"], "\n")
cat("ABSI Model C-index:", summary(cox_absi)$concordance["C"], "\n")
cat("BRI Model C-index:", summary(cox_bri)$concordance["C"], "\n")

# 4. Generate Forest Plot for Model Comparison
# Combine model outputs into a single dataset
results_bmi  <- tidy(cox_bmi, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "BMI_z")
results_whtr <- tidy(cox_whtr, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "WHtR_z")
results_absi <- tidy(cox_absi, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "ABSI_z")
results_bri  <- tidy(cox_bri, exponentiate = TRUE, conf.int = TRUE) %>% filter(term == "BRI_z")

forest_data <- bind_rows(results_bmi, results_whtr, results_absi, results_bri) %>%
  mutate(
    label = c("Body Mass Index (BMI)", "Waist-to-Height Ratio (WHtR)", "A Body Shape Index (ABSI)", "Body Roundness Index (BRI)"),
    p_val_text = ifelse(p.value < 0.001, "<.001", sprintf("%.3f", p.value))
  )

# Define labels for the forest plot
table_text <- list(
  c("Obesity Indices (Z-score)", forest_data$label),
  c("HR (95% CI)", sprintf("%.2f (%.2f–%.2f)", forest_data$estimate, forest_data$conf.low, forest_data$conf.high)),
  c("P Value", forest_data$p_val_text)
)

# Plotting with JAMA style
jama_blue <- "#374E55FF"
model_note <- "Models were adjusted for age, gender, rural status, smoking, alcohol, hypertension, and diabetes."

p_forest_final <- forestplot(
  labeltext = table_text,
  mean  = c(NA, forest_data$estimate),
  lower = c(NA, forest_data$conf.low),
  upper = c(NA, forest_data$conf.high),
  is.summary = c(TRUE, rep(FALSE, nrow(forest_data))),
  clip = c(0.9, 1.2),
  xlog = FALSE, zero = 1,
  boxsize = 0.25,
  colgap = unit(5, "mm"),
  col = fpColors(box = jama_blue, lines = jama_blue, zero = "grey50"),
  txt_gp = fpTxtGp(label = gpar(fontfamily = "sans", cex = 0.95),
                   ticks = gpar(cex = 0.85),
                   xlab  = gpar(cex = 1.0, fontface = "bold")),
  xlab = "Hazard Ratio (per 1-SD Increase)"
)

# Post-processing to add model notes
final_plot <- forestplot(
  labeltext = table_text,
  mean  = c(NA, forest_data$estimate),
  lower = c(NA, forest_data$conf.low),
  upper = c(NA, forest_data$conf.high),
  is.summary = c(TRUE, rep(FALSE, nrow(forest_data))),
  fn.postform = function(x) {
    grid.text(model_note, x = unit(0.05, "npc"), y = unit(0.05, "npc"), 
              just = "left", gp = gpar(fontsize = 9, col = "grey30", fontface = "italic"))
    return(x)
  },
  xlab = "Hazard Ratio (per 1-SD Increase)"
)

# Save high-quality images
pdf("Figure_Model_Comparison_Forest.pdf", width = 10, height = 4)
print(final_plot)
dev.off()

tiff("Figure_Model_Comparison_Forest.tiff", width = 10, height = 4, units = "in", res = 350, compression = "lzw")
print(final_plot)
dev.off()

# ==============================================================================
# Correlation_Analysis.R
# Purpose: Generate a panoramic correlation heatmap to visualize the associations 
#          between obesity phenotypes and metabolic biomarkers.
# ==============================================================================
library(dplyr)
library(ggcorrplot)
library(ggplot2)
library(ggsci)
library(haven)

cat("\n================ Pipeline: Panoramic Correlation Analysis ================\n")

# 1. Data Preparation: Extracting variables from Wave 1 (2011) for cross-sectional correlation
full_cor_data <- charls_data %>%
  mutate(wave = zap_labels(wave)) %>% 
  filter(wave == 1) %>%
  mutate(
    mheight_m = ifelse(mheight > 10, mheight / 100, mheight),
    mheight_m = ifelse(mheight_m < 1.0 | mheight_m > 2.5, NA, mheight_m),
    mwaist_m  = mwaist / 100,
    mwaist_m  = ifelse(mwaist_m < 0.4 | mwaist_m > 1.6, NA, mwaist_m),
    bmi_clean = ifelse(bmi < 10 | bmi > 60, NA, bmi)
  ) %>%
  filter(!is.na(mheight_m), !is.na(mwaist_m), !is.na(bmi_clean)) %>%
  mutate(
    WHtR = mwaist_m / mheight_m,
    ABSI = mwaist_m / ((bmi_clean^(2/3)) * (mheight_m^(1/2))),
    BRI  = 364.2 - 365.5 * sqrt(1 - (mwaist_m / (pi * mheight_m))^2)
  ) %>%
  filter(is.finite(ABSI)) %>%
  dplyr::select(
    Weight = mweight, Waist = mwaist, BMI = bmi_clean, WHtR, ABSI, BRI,
    `TyG-BMI` = tyg_bmi, TyG = tyg, TG = bl_tg, HDL = bl_hdl, LDL = bl_ldl, CRP = bl_crp
  ) %>%
  na.omit()

cat("Clean sample size for correlation analysis:", nrow(full_cor_data), "\n")

# 2. Correlation Matrix Calculation (Spearman for non-parametric metabolic data)
corr_matrix_full <- cor(full_cor_data, method = "spearman")

# Helper function for P-value matrix
cor_pmat_full <- function(mat) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = "spearman", exact = FALSE)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}
p_mat_full <- cor_pmat_full(full_cor_data)

# 3. Visualization: Panoramic Correlation Heatmap
jama_colors <- c("#374E55FF", "white", "#B24745FF")

p_heatmap <- ggcorrplot(
  corr_matrix_full, 
  method = "circle", type = "lower", hc.order = TRUE,
  p.mat = p_mat_full, sig.level = 0.05, insig = "pch",
  pch.cex = 4, pch.col = "grey30",
  lab = TRUE, lab_size = 3.0,
  colors = jama_colors, outline.color = "white",
  title = "Spearman Correlation: Obesity Phenotypes & Metabolic Indicators",
  ggtheme = theme_minimal()
) +
  theme(
    text = element_text(family = "sans"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black")
  )

# 4. Save High-Quality Plots
ggsave("Figure_Correlation_Heatmap.pdf", plot = p_heatmap, width = 8, height = 7)
ggsave("Figure_Correlation_Heatmap.tiff", plot = p_heatmap, width = 8, height = 7, dpi = 350, compression = "lzw")

cat("\n================ Pipeline Complete ================\n")
