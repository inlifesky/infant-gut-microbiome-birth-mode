# ─────────────────────────────────────────────────────────────────────────────
# 12 — Pathway differential abundance: ALDEx2
# Mirrors Step 06 (species ALDEx2) but on HUMAnN MetaCyc pathway abundances.
# Also saves processed pathway matrix for Steps 13–14 (avoid re-download).
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(curatedMetagenomicData); library(SummarizedExperiment)
  library(ALDEx2); library(ggplot2)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running

# --- Load + filter pathway data (same logic as Step 05) ---
pw <- curatedMetagenomicData("2021-03-31.ShaoY_2019.pathway_abundance",
                             dryrun = FALSE, rownames = "short")[[1]]
ids <- readRDS(file.path(proj, "results", "01_processed.rds"))$sample_ids
pw <- pw[, intersect(colnames(pw), ids)]
grp <- factor(pw$born_method, levels = c("vaginal", "c_section"))
meta_pw <- data.frame(
  sample_id   = colnames(pw),
  born_method = as.character(grp),
  gender      = pw$gender,
  stringsAsFactors = FALSE
)

M <- assay(pw)
rn <- rownames(M)
M <- M[!grepl("\\|", rn) & !grepl("UNMAPPED|UNINTEGRATED", rn), , drop = FALSE]
if (max(colSums(M), na.rm = TRUE) > 1.5) M <- sweep(M, 2, colSums(M), "/")
prev <- rowMeans(M > 0); M <- M[prev >= 0.10, , drop = FALSE]
cat("Pathways after filter:", nrow(M), " Samples:", ncol(M), "\n")

# Save for Steps 13–14
saveRDS(list(abundance = M, meta = meta_pw),
        file.path(proj, "results", "05_pathway_processed.rds"))

# --- ALDEx2 ---
pseudo_counts <- round(M * 1e6)
storage.mode(pseudo_counts) <- "integer"

aldex_res <- aldex(pseudo_counts, conditions = as.character(grp),
                   mc.samples = 128, test = "t", effect = TRUE, verbose = TRUE)

res <- data.frame(
  pathway  = rownames(aldex_res),
  effect   = aldex_res$effect,
  wi_p     = aldex_res$wi.ep,
  wi_fdr   = aldex_res$wi.eBH,
  we_p     = aldex_res$we.ep,
  we_fdr   = aldex_res$we.eBH,
  stringsAsFactors = FALSE
)
res <- res[order(res$wi_fdr), ]
write.csv(res, file.path(proj, "results", "pathway_aldex2.csv"), row.names = FALSE)

n_sig <- sum(res$wi_fdr < 0.05, na.rm = TRUE)
cat("Pathway ALDEx2 significant (FDR<0.05):", n_sig, "of", nrow(res), "\n")

res$sig <- ifelse(res$wi_fdr < 0.05, "FDR<0.05", "ns")
p <- ggplot(res, aes(effect, -log10(wi_fdr), colour = sig)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c("FDR<0.05" = "#7570B3", "ns" = "grey70")) +
  labs(title = "ALDEx2: metabolic pathway DA by birth mode",
       subtitle = sprintf("Wilcoxon (128 MC) — %d / %d FDR<0.05", n_sig, nrow(res)),
       x = "ALDEx2 effect size", y = "-log10(FDR)", colour = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(proj, "figures", "pathway_aldex2.png"), p, width = 8, height = 5, dpi = 200)

cat("[12 done] results/pathway_aldex2.csv + figures/pathway_aldex2.png\n")
