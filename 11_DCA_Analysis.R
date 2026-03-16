# ==============================================================================
# Script Name: 11_DCA_Analysis.R
# Project: Obesity Phenotypes and Incident Knee Osteoarthritis (KOA)
# Purpose: Decision Curve Analysis (DCA) for the WHtR-based nomogram
# ==============================================================================

# 1. Load necessary libraries
library(survival)
library(rms)
library(dcurves)
library(ggplot2)
library(ggsci)
library(dplyr)

cat("\n================ Pipeline: Decision Curve Analysis (DCA) ================\n")

# 2. Check for required model objects
if (!exists("f_nomo") || !exists("data_cox_nomo")) {
  stop("Error: Required model 'f_nomo' or data 'data_cox_nomo' not found in environment.")
}

# 3. Calculate 5-year KOA Risk
# S is the survival function object derived from 'f_nomo'
S <- Survival(f_nomo)

# Convert linear predictor to 5-year incident risk (1 - Survival probability)
data_cox_nomo$risk_5yr <- 1 - S(5, predict(f_nomo))

# 4. Define outcome variable for DCA
# DCA requires a binary status (1 if KOA occurred within 5 years, 0 otherwise)
data_cox_nomo <- data_cox_nomo %>%
  mutate(status_5yr = ifelse(follow_up_time <= 5 & status == 1, 1, 0))

# 5. Execute Decision Curve Analysis
dca_res <- dca(status_5yr ~ risk_5yr, data = data_cox_nomo)

# 6. Generate Publication-Quality Plot
p_dca <- plot(dca_res) + 
  theme_bw() + 
  scale_color_jama() + 
  labs(
    title = "Decision Curve Analysis for 5-Year KOA Risk",
    subtitle = "Comparing Clinical Net Benefit of the WHtR-based Nomogram",
    x = "Threshold Probability",
    y = "Net Benefit"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, color = "grey30", hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = c(0.7, 0.75), 
    legend.title = element_blank(),
    text = element_text(family = "sans")
  )

# 7. Print and Export Plots
print(p_dca)

# Save as PDF (Vector format for manuscript)
ggsave("Figure_DCA_5yr.pdf", p_dca, width = 8, height = 6)

# Save as TIFF (350 DPI for submission)
ggsave("Figure_DCA_5yr.tiff", p_dca, width = 8, height = 6, dpi = 350, compression = "lzw")


cat("\n================ DCA Analysis Complete ================\n")