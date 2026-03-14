# ==============================================================================
# Script Name: 09_Sensitivity_Analysis_E-value.R
# Project: Obesity Phenotypes and Incident Knee Osteoarthritis (KOA)
# Purpose: Perform E-value analysis to quantify the robustness of the association 
#          between WHtR and KOA against unmeasured confounding.
# ==============================================================================

# 1. Load necessary libraries
library(EValue)
library(survival)

cat("\n================ Pipeline: E-value Analysis for Sensitivity ================\n")

# ------------------------------------------------------------------------------
# 1. Extract HR and 95% CI Lower Bound
# ------------------------------------------------------------------------------
# Extract the point estimate (aHR) and the lower bound of the 95% CI 
# from the primary multivariable Cox model for WHtR_z.

# Get coefficients and confidence intervals
coef_summary <- summary(cox_whtr)$coefficients
conf_int <- confint(cox_whtr)

# Extract values for WHtR_z
hr_val  <- as.numeric(exp(coef_summary["WHtR_z", "coef"]))
lci_val <- as.numeric(exp(conf_int["WHtR_z", 1]))

cat("Primary association (WHtR_z per 1-SD):\n")
cat("Adjusted Hazard Ratio (aHR):", round(hr_val, 3), "\n")
cat("95% CI Lower Bound:", round(lci_val, 3), "\n")

# ------------------------------------------------------------------------------
# 2. Compute E-value
# ------------------------------------------------------------------------------
# Note: 'rare = FALSE' is used because the incidence of KOA in this cohort 
# is relatively high, which is more appropriate for standard Cox models.

e_val_res <- evalues.HR(est = hr_val, lo = lci_val, rare = FALSE)

cat("\n================ E-value Results ================\n")
print(e_val_res)

# ------------------------------------------------------------------------------
# 3. Interpretation Guidance
# ------------------------------------------------------------------------------
cat("\nInterpretation Note:\n")
cat("The E-value indicates the minimum strength of association that an unmeasured \n",
    "confounder would need to have with both the obesity index and the risk of KOA \n",
    "to explain away the observed effect. \n",
    "Thresholds: E-value > 1.25 is generally considered robust; > 2.0 is very robust.\n")

cat("\n================ Pipeline Complete ================\n")