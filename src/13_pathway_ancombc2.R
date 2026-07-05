# ─────────────────────────────────────────────────────────────────────────────
# 13 — Pathway differential abundance: ANCOM-BC2
# Uses pathway matrix saved by Step 12.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ANCOMBC); library(TreeSummarizedExperiment)
  library(SummarizedExperiment); library(ggplot2)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running
d <- readRDS(file.path(proj, "results", "05_pathway_processed.rds"))
M <- d$abundance; meta <- d$meta

pseudo_counts <- round(M * 1e6)
storage.mode(pseudo_counts) <- "integer"

tse <- TreeSummarizedExperiment(
  assays  = list(counts = pseudo_counts),
  colData = DataFrame(meta)
)
cat("ANCOM-BC2 pathway input:", nrow(pseudo_counts), "features x", ncol(pseudo_counts), "samples\n")

ancom_out <- ancombc2(
  data         = tse,
  assay_name   = "counts",
  fix_formula  = "born_method",
  p_adj_method = "BH",
  prv_cut      = 0.10,
  verbose      = TRUE
)

res <- ancom_out$res
lfc_col  <- grep("^lfc_", colnames(res), value = TRUE)[1]
p_col    <- grep("^p_born", colnames(res), value = TRUE)[1]
q_col    <- grep("^q_born", colnames(res), value = TRUE)[1]

out <- data.frame(
  pathway = res$taxon,
  lfc     = res[[lfc_col]],
  p       = res[[p_col]],
  fdr     = res[[q_col]],
  stringsAsFactors = FALSE
)
out <- out[order(out$fdr), ]
write.csv(out, file.path(proj, "results", "pathway_ancombc2.csv"), row.names = FALSE)

n_sig <- sum(out$fdr < 0.05, na.rm = TRUE)
cat("Pathway ANCOM-BC2 significant (FDR<0.05):", n_sig, "of", nrow(out), "\n")

out$sig <- ifelse(out$fdr < 0.05, "FDR<0.05", "ns")
p <- ggplot(out, aes(lfc, -log10(fdr), colour = sig)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c("FDR<0.05" = "#E7298A", "ns" = "grey70")) +
  labs(title = "ANCOM-BC2: metabolic pathway DA by birth mode",
       subtitle = sprintf("Bias-corrected — %d / %d FDR<0.05", n_sig, nrow(out)),
       x = "Log fold change (C-section vs vaginal)", y = "-log10(FDR)", colour = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(proj, "figures", "pathway_ancombc2.png"), p, width = 8, height = 5, dpi = 200)

cat("[13 done] results/pathway_ancombc2.csv + figures/pathway_ancombc2.png\n")
