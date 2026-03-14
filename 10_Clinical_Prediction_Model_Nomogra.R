# ==============================================================================
# Script Name: 10_Clinical_Prediction_Model_Nomogram.R
# Project: Obesity Phenotypes and Incident Knee Osteoarthritis (KOA)
# Purpose: Develop and validate a clinical nomogram using WHtR, assess 
#          calibration, and compute Brier scores for predictive accuracy.
# ==============================================================================

library(rms)
library(survival)
library(ggplot2)
library(riskRegression)
library(flextable)

cat("\n================ Pipeline: Clinical Nomogram Construction ================\n")

# ------------------------------------------------------------------------------
# 1. Model Preparation and Fitting
# ------------------------------------------------------------------------------

# Define factor variables for nomogram readability
data_cox_nomo <- data_cox %>%
  mutate(
    gender = factor(gender, levels=c(0,1), labels=c("Female", "Male")),
    rural  = factor(rural,  levels=c(0,1), labels=c("Urban", "Rural")),
    hibpe  = factor(hibpe,  levels=c(0,1), labels=c("No", "Yes"))
  )

# Label variables for nomogram display
label(data_cox_nomo$WHtR)     <- "Waist-to-Height Ratio"
label(data_cox_nomo$age_base) <- "Age (years)"
label(data_cox_nomo$gender)   <- "Sex"
label(data_cox_nomo$rural)    <- "Residence"
label(data_cox_nomo$hibpe)    <- "Hypertension"


# Setup datadist for rms package
dd <- datadist(data_cox_nomo)
options(datadist = "dd")

# Fit Cox model for nomogram
f_nomo <- cph(Surv(follow_up_time, status) ~ WHtR + age_base + gender + rural + hibpe, 
              data = data_cox_nomo, x = TRUE, y = TRUE, surv = TRUE)

# ------------------------------------------------------------------------------
# 2. Nomogram Plotting
# ------------------------------------------------------------------------------

nom <- nomogram(f_nomo, lp = FALSE, maxscale = 100, fun = function(x) surv(5, x), 
                funlabel = "5-Year Probability of KOA")

# Save high-resolution Nomogram
tiff("Figure_Nomogram_WHtR.tiff", width = 10, height = 7, units = "in", 
     res = 350, compression = "lzw", family = "sans")
plot(nom, xfrac = 0.25, cex.axis = 0.85, cex.var = 1.0, lmgp = 0.2)
dev.off()

# ------------------------------------------------------------------------------
# 3. Calibration Curve (1000 Bootstraps)
# ------------------------------------------------------------------------------

# Fit model with time.inc for calibration
data_clean <- na.omit(data_cox_nomo)
f_cal <- cph(Surv(follow_up_time, status) ~ WHtR + age_base + gender + rural + hibpe + diabe, 
             data = data_clean, x = TRUE, y = TRUE, surv = TRUE, time.inc = 5)

cal_plot <- calibrate(f_cal, method="boot", B=1000, u=5, m=round(nrow(data_clean)/3))

tiff("Figure_Calibration_5yr.tiff", width = 6, height = 6, units = "in", 
     res = 350, compression = "lzw", family = "sans")
plot(cal_plot, lwd=2, lty=1, errbar.col="gray50", xlab="Predicted 5-Year Probability", 
     ylab="Observed 5-Year Probability", subtitles=FALSE)
abline(0, 1, lty=2, lwd=1.5, col="gray40")
dev.off()

# ------------------------------------------------------------------------------
# 4. Comprehensive Accuracy: Brier Score
# ------------------------------------------------------------------------------

# Calculate Brier score for 5-year prediction
bs <- Score(list("WHtR_Model" = f_cal), 
            formula = Surv(follow_up_time, status) ~ 1, 
            data = data_clean, 
            times = 5, 
            metrics = "brier")

brier_score <- bs$Brier$score$Brier[2]
cat("\nBrier Score at 5 years:", round(brier_score, 4), "\n")

cat("\n================ Pipeline Complete ================\n")