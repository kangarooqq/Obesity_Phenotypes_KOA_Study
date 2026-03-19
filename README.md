**Obesity Phenotypes and Incident Knee Osteoarthritis: A Longitudinal Cohort Study**

This repository contains the R analysis scripts and documentation for the study titled "Obesity Phenotypes and Incident Knee Osteoarthritis: A Longitudinal Cohort Study".
This research investigates the predictive performance of four obesity indices—Body Mass Index (**BMI**), Waist-to-Height Ratio (**WHtR**), A Body Shape Index (**ABSI**), and Body Roundness Index (**BRI**)—on incident knee osteoarthritis (**KOA**) using data from the China Health and Retirement Longitudinal Study (CHARLS).

**Project Structure**
The analysis is organized into a modular pipeline to ensure transparency and reproducibility:
Script	Description

**01_Data_Cleaning**.R	Data import, outlier detection, and initial quality control.

**02_Cohort_Construction**.R	Definition of the incident cohort and exclusion criteria.

**03_Obesity_Comparison**.R	Comparative analysis of the four obesity indices.

**04_WHtR_Analysis**.R	Deep-dive into WHtR: non-linear RCS modeling & subgroup analysis.

**05_ABSI_Analysis**.R	Deep-dive into ABSI: non-linear RCS modeling & subgroup analysis.

**06_BMI_Analysis**.R	Deep-dive into BMI: non-linear RCS modeling & subgroup analysis.

**07_BRI_Analysis**.R	Deep-dive into BRI: non-linear RCS modeling & subgroup analysis.

**08_Joint_Modeling**.R	Joint model analysis to test the "Dual-hit" hypothesis.

**09_Sensitivity_Evalue**.R	Sensitivity analysis using E-values to assess robustness against confounding.

**10_Clinical_Prediction_Model_Nomogram**.R	Development of the clinical nomogram, calibration, and Brier score validation.

**11_DCA_Analysis**.R Decision Curve Analysis (DCA) for the WHtR-based nomogram.

**12_NRI_IDI_Analysis**.R Calculate continuous Net Reclassification Improvement (NRI) and Integrated Discrimination Improvement (IDI) for model comparison.

**Prerequisites**
To run these scripts, you will need the following R packages installed:

code
**R**
install.packages(c("dplyr", "tidyr", "survival", "survminer", "rms", "gtsummary", 
                   "flextable", "officer", "forestplot", "broom", "ggplot2", 
                   "ggsci", "riskRegression"))
                   
Running the Pipeline
Data Preparation: Place your raw charls.dta file in the project root directory.
Execution: Execute the scripts in numeric order from 01 to 12.
Reproducibility: Ensure your working directory is set to the project root. The scripts are designed to save high-quality publication-ready plots (PDF/TIFF) directly to the output folder.

Authors
**Cheng-Wei Kang, MD (First Author)**

Ethical Declarations
This study was approved by the Ethics Review Board of Peking University (IRB00001052–11015). All participants provided written informed consent.
Acknowledgments
We thank the CHARLS research team and study participants for their valuable time and effort.
