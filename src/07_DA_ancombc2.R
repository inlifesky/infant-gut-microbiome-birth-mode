# ─────────────────────────────────────────────────────────────────────────────
# 07 — Differential abundance: ANCOM-BC2
# Bias-corrected log-linear model; estimates sampling fraction per sample.
# Recommended by Nearing 2022 as conservative with well-controlled FDR.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ANCOMBC); library(TreeSummarizedExperiment)
  library(SummarizedExperiment); library(ggplot2)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running
d <- readRDS(file.path(proj, "results", "01_processed.rds"))
ab <- d$abundance; meta <- d$meta

# ANCOM-BC2 expects a TreeSummarizedExperiment with count-like integers
pseudo_counts <- round(ab * 1e6)
storage.mode(pseudo_counts) <- "integer"

# Build a minimal TSE
tse <- TreeSummarizedExperiment(
  assays = list(counts = pseudo_counts),
  colData = DataFrame(meta)
)

cat("ANCOM-BC2 input:", nrow(pseudo_counts), "features x", ncol(pseudo_counts), "samples\n")

# Run ANCOM-BC2
ancom_out <- ancombc2(
  data        = tse,
  assay_name  = "counts",
  fix_formula = "born_method",
  p_adj_method = "BH",
  prv_cut     = 0.10,
  verbose     = TRUE
)

res <- ancom_out$res
# The coefficient column for born_method c_section vs vaginal (reference)
# Column names: lfc_born_methodc_section, p_born_methodc_section, q_born_methodc_section, etc.
lfc_col <- grep("^lfc_", colnames(res), value = TRUE)[1]
p_col   <- grep("^p_born", colnames(res), value = TRUE)[1]
q_col   <- grep("^q_born", colnames(res), value = TRUE)[1]
diff_col <- grep("^diff_born", colnames(res), value = TRUE)[1]

out <- data.frame(
  taxon   = res$taxon,
  lfc     = res[[lfc_col]],
  p       = res[[p_col]],
  fdr     = res[[q_col]],
  diff    = if (!is.null(diff_col)) res[[diff_col]] else (res[[q_col]] < 0.05),
  stringsAsFactors = FALSE
)
out <- out[order(out$fdr), ]
write.csv(out, file.path(proj, "results", "DA_ancombc2.csv"), row.names = FALSE)

n_sig <- sum(out$fdr < 0.05, na.rm = TRUE)
cat("ANCOM-BC2 significant (FDR<0.05):", n_sig, "of", nrow(out), "\n")

# Volcano plot
out$sig <- ifelse(out$fdr < 0.05, "FDR<0.05", "ns")
p <- ggplot(out, aes(lfc, -log10(fdr), colour = sig)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70") +
  scale_colour_manual(values = c("FDR<0.05" = "#E7298A", "ns" = "grey70")) +
  labs(title = "ANCOM-BC2 differential abundance by birth mode",
       subtitle = sprintf("Bias-corrected — %d / %d FDR<0.05", n_sig, nrow(out)),
       x = "Log fold change (C-section vs vaginal)", y = "-log10(FDR)", colour = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(proj, "figures", "da_ancombc2.png"), p, width = 8, height = 5, dpi = 200)

cat("[07 done] results/DA_ancombc2.csv + figures/da_ancombc2.png\n")
