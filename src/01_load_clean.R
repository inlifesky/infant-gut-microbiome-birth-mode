# ─────────────────────────────────────────────────────────────────────────────
# 01 — Load + clean: ShaoY_2019 (Baby Biome Study) infant gut metagenome
# Outcome: born_method (c_section vs vaginal) — a canonical early-life determinant
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(TreeSummarizedExperiment)
  library(SummarizedExperiment)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running
dir.create(file.path(proj, "results"), showWarnings = FALSE, recursive = TRUE)

# --- pull species-level relative abundance (cached in ExperimentHub after 1st run) ---
tse <- curatedMetagenomicData("2021-03-31.ShaoY_2019.relative_abundance",
                              dryrun = FALSE, rownames = "short")[[1]]
cat("raw TSE:", paste(dim(tse), collapse = " x "), "(features x samples)\n")

# --- restrict to INFANT samples with a clear binary birth mode ---
# (study contains 174 adult/mother samples — exclude; outcome is the infant's birth mode)
keep <- !is.na(tse$born_method) & tse$born_method %in% c("c_section", "vaginal") &
        tse$age_category == "newborn"
tse <- tse[, keep]
cat("infant samples with birth mode:", ncol(tse), "\n")

# --- collapse to ONE sample per subject: earliest infant_age (neonatal) ---
# birth mode is a between-subject variable + the cohort is longitudinal -> avoid
# pseudoreplication by keeping each infant's earliest timepoint.
ord <- order(tse$subject_id, tse$infant_age, na.last = TRUE)
tse <- tse[, ord]
first <- !duplicated(tse$subject_id)
tse <- tse[, first]
tse$born_method <- factor(tse$born_method, levels = c("vaginal", "c_section"))
cat("after one-per-subject (earliest timepoint):", ncol(tse), "independent infants\n")
cat("by birth mode:\n"); print(table(tse$born_method))
cat("earliest-timepoint infant_age (days) summary:\n"); print(summary(tse$infant_age))

# --- abundance matrix: features x samples, relative abundance (0-100) ---
ab <- assay(tse)
# convert to proportions if scaled to 100
if (max(colSums(ab), na.rm = TRUE) > 1.5) ab <- ab / 100

# --- prevalence + abundance filter (standard microbiome QC) ---
prev <- rowMeans(ab > 0)
maxab <- apply(ab, 1, max)
keepf <- prev >= 0.10 & maxab >= 1e-4        # present in >=10% samples, reaches 0.01%
cat("\nfeatures: ", nrow(ab), " -> kept after prevalence/abundance filter: ", sum(keepf), "\n", sep = "")
ab <- ab[keepf, , drop = FALSE]

# drop any sample left with zero total across kept features (avoids meaningless Bray-Curtis rows)
nz <- colSums(ab) > 0
if (any(!nz)) { cat("dropping", sum(!nz), "sample(s) with zero retained abundance\n"); ab <- ab[, nz, drop = FALSE]; tse <- tse[, nz] }

# --- assemble tidy objects ---
want <- c("subject_id", "born_method", "age_category", "country", "feeding_practice", "gender")
have <- intersect(want, colnames(colData(tse)))
cat("\nmetadata columns available:", paste(have, collapse = ", "), "\n")
meta <- as.data.frame(colData(tse))[, have, drop = FALSE]
meta$born_method <- factor(meta$born_method, levels = c("vaginal", "c_section"))

saveRDS(list(abundance = ab, meta = meta, tse = tse, sample_ids = colnames(ab)),
        file.path(proj, "results", "01_processed.rds"))

# provenance
sink(file.path(proj, "results", "dataset_provenance.txt"))
cat("DATASET PROVENANCE — infant microbiome birth-mode demo\n")
cat("source: curatedMetagenomicData, resource 2021-03-31.ShaoY_2019.relative_abundance\n")
cat("cohort: Shao et al. 2019, Baby Biome Study (infant gut shotgun metagenome)\n")
cat("outcome: born_method (vaginal vs c_section) — early-life determinant\n")
cat("cMD package version:", as.character(packageVersion("curatedMetagenomicData")), "\n\n")
cat("samples (after filter):", ncol(tse), "\n")
print(table(tse$born_method))
cat("unique subjects:", length(unique(tse$subject_id)), "\n")
cat("features kept (prevalence>=10% & maxRA>=1e-4):", nrow(ab), "\n")
cat("\nNote: secondary analysis of public uniformly-processed MetaPhlAn profiles;\n")
cat("no raw-read processing performed.\n")
sink()

cat("\n[01 done] saved results/01_processed.rds + dataset_provenance.txt\n")
