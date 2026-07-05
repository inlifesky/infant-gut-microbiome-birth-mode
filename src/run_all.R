# Run the full pipeline end-to-end and log the environment.
proj <- getwd()  # set working directory to repo root before running
src <- file.path(proj, "src")

steps <- c(
  # --- Cross-sectional dataset + noise-filter methods ---
  "01_load_clean.R", "02_diversity.R", "03_differential_abundance.R",
  "04_ml_classifier.R", "05_pathway.R",
  "06_DA_aldex2.R", "07_DA_ancombc2.R", "08_DA_maaslin2.R",
  "09_DA_concordance.R", "10_classifier_robustness.R",
  "12_pathway_aldex2.R", "13_pathway_ancombc2.R",
  "14_pathway_maaslin2.R", "15_pathway_concordance.R",
  "11_clinical_context.R",
  # --- Longitudinal + clinical-stratification analyses (Phase A) ---
  "01b_load_full_longitudinal.R",
  "17_longitudinal_recovery.R",
  "18_antibiotic_deconvolution.R",
  "19_cs_type_stratification.R",
  "20_rf_as_noise_filter.R",
  "21_intervention_report.R",
  "23_verify_blongum_reversal.R",
  "24_pathway_longitudinal_antibiotic.R",
  "25_pathway_intervention_report.R",
  "26_per_pathway_longitudinal.R",
  # --- Composite figures ---
  "16_publication_figure.R",
  "22_intervention_figure.R"
)

for (s in steps) {
  cat("\n==== running", s, "====\n")
  source(file.path(src, s))
}

writeLines(capture.output(sessionInfo()), file.path(proj, "results", "sessionInfo.txt"))
cat("\n[run_all done] sessionInfo written to results/sessionInfo.txt\n")
