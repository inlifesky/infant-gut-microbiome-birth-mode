# ─────────────────────────────────────────────────────────────────────────────
# 06 — Differential abundance: ALDEx2
# Dirichlet-multinomial Monte Carlo on proportions → CLR instances → Wilcoxon
# Rationale: ALDEx2 handles zero-replacement + compositional uncertainty
#   internally (128 MC instances), more rigorous than single-pseudocount CLR.
#   Top-recommended by Nearing et al. 2022 Nat Commun; Yang et al. 2025 Brief Bioinform.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(ALDEx2); library(ggplot2) })
set.seed(42)
proj <- getwd()  # set working directory to repo root before running
d <- readRDS(file.path(proj, "results", "01_processed.rds"))
ab <- d$abundance; meta <- d$meta

# ALDEx2 expects a features × samples matrix of COUNTS (integers).
# cMD gives proportions (0–1). Multiply by a large library size to create
# pseudo-counts that preserve relative structure.
# (This is the standard approach when only proportions are available;
#  see ALDEx2 vignette "working with relative abundance data".)
pseudo_counts <- round(ab * 1e6)
storage.mode(pseudo_counts) <- "integer"

grp <- as.character(meta$born_method)
cat("ALDEx2 input:", nrow(pseudo_counts), "features x", ncol(pseudo_counts), "samples\n")
cat("groups:", paste(names(table(grp)), table(grp), sep = "=", collapse = ", "), "\n")

# Run ALDEx2: 128 MC instances, Wilcoxon test (matches our Step 03 test choice)
aldex_res <- aldex(pseudo_counts, conditions = grp, mc.samples = 128,
                   test = "t", effect = TRUE, verbose = TRUE)

# Extract results
res <- data.frame(
  taxon     = rownames(aldex_res),
  effect    = aldex_res$effect,        # effect size (CLR-based)
  wi_p      = aldex_res$wi.ep,         # Wilcoxon p (expected)
  wi_fdr    = aldex_res$wi.eBH,        # Wilcoxon BH-FDR
  we_p      = aldex_res$we.ep,         # Welch t p (expected)
  we_fdr    = aldex_res$we.eBH,        # Welch t BH-FDR
  stringsAsFactors = FALSE
)
res <- res[order(res$wi_fdr), ]

write.csv(res, file.path(proj, "results", "DA_aldex2.csv"), row.names = FALSE)

n_sig_wi <- sum(res$wi_fdr < 0.05, na.rm = TRUE)
n_sig_we <- sum(res$we_fdr < 0.05, na.rm = TRUE)
cat("ALDEx2 significant (FDR<0.05): Wilcoxon =", n_sig_wi, ", Welch t =", n_sig_we,
    "of", nrow(res), "\n")

# Effect-size plot
res$sig <- ifelse(res$wi_fdr < 0.05, "FDR<0.05", "ns")
p <- ggplot(res, aes(effect, -log10(wi_fdr), colour = sig)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70") +
  scale_colour_manual(values = c("FDR<0.05" = "#7570B3", "ns" = "grey70")) +
  labs(title = "ALDEx2 differential abundance by birth mode",
       subtitle = sprintf("Wilcoxon (128 MC) — %d / %d FDR<0.05", n_sig_wi, nrow(res)),
       x = "ALDEx2 effect size", y = "-log10(FDR)", colour = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(proj, "figures", "da_aldex2.png"), p, width = 8, height = 5, dpi = 200)

cat("[06 done] results/DA_aldex2.csv + figures/da_aldex2.png\n")
