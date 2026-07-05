# ─────────────────────────────────────────────────────────────────────────────
# 01b — Load BBS with FULL metadata, ALL timepoints
#
# Why this exists: the original step 01 collapsed BBS to one-sample-per-subject
# (day 4) and kept only 5 metadata columns. That design served a method demo
# ("can methods discriminate birth mode?") — but it cannot answer the clinical
# questions that motivate this project:
#
#   - When does the CS-vaginal gap close? (needs infant_age + all timepoints)
#   - Is the CS effect actually a peripartum antibiotic effect?
#     (needs antibiotics_current_use)
#   - Do elective and emergency CS disrupt the microbiome equally?
#     (needs c_section_type)
#   - Does feeding modulate the birth-mode effect?
#     (needs feeding_practice)
#
# We keep step 01's single-timepoint subset for the cross-sectional methods
# concordance (step 03-15) as the noise-filter foundation. THIS script adds the
# longitudinal+covariate-rich dataset for stratified clinical analysis (17-20).
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(dplyr)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running

# ── Pull BBS metadata + abundance ──
meta_all <- sampleMetadata |> filter(study_name == "ShaoY_2019")
cat("Total BBS samples:", nrow(meta_all), "/ subjects:",
    length(unique(meta_all$subject_id)), "\n")

# Keep infant samples only (exclude mothers; family_role tracks this)
meta_inf <- meta_all |> filter(family_role == "child" | is.na(family_role))
cat("Infant samples (family_role=child or NA):", nrow(meta_inf), "\n")

# Variables we need for the clinical questions
keep_cols <- c("sample_id", "subject_id", "born_method", "age_category",
               "country", "gender", "infant_age", "days_from_first_collection",
               "antibiotics_current_use", "c_section_type",
               "feeding_practice", "premature", "birth_weight",
               "gestational_age", "family_role", "visit_number",
               "number_reads", "number_bases")
meta_inf <- meta_inf[, intersect(keep_cols, colnames(meta_inf))]

# Non-NA coverage report (do not drop NA — let downstream stratifications decide)
cat("\nCoverage of clinical variables:\n")
for (v in setdiff(colnames(meta_inf), c("sample_id", "subject_id"))) {
  n_ok <- sum(!is.na(meta_inf[[v]]))
  cat(sprintf("  %-30s %d / %d  (%.0f%%)\n", v, n_ok, nrow(meta_inf),
              100 * n_ok / nrow(meta_inf)))
}
cat("\nTimepoints per subject (median, IQR):\n")
tps <- table(meta_inf$subject_id)
print(quantile(tps, c(0.25, 0.5, 0.75)))

# ── Pull abundance for the kept samples ──
cat("\nFetching MetaPhlAn3 species abundance for ShaoY_2019...\n")
tse_full <- curatedMetagenomicData("ShaoY_2019.relative_abundance",
                                   dryrun = FALSE, counts = FALSE,
                                   rownames = "short")[[1]]
cat("TSE assay dim:", dim(assay(tse_full)), "\n")

# Subset assay to infant samples we kept
samples_keep <- intersect(meta_inf$sample_id, colnames(tse_full))
tse_inf <- tse_full[, samples_keep]
meta_inf <- meta_inf[match(samples_keep, meta_inf$sample_id), ]
stopifnot(identical(meta_inf$sample_id, colnames(tse_inf)))

# Update colData with the rich metadata
colData(tse_inf) <- DataFrame(meta_inf, row.names = meta_inf$sample_id)

# Save
saveRDS(list(tse = tse_inf, meta = meta_inf,
             abundance = as.matrix(assay(tse_inf))),
        file.path(proj, "results", "01b_processed_longitudinal.rds"))

cat(sprintf("\n[01b done] saved %d samples × %d species, %d subjects\n",
            ncol(tse_inf), nrow(tse_inf), length(unique(meta_inf$subject_id))))
