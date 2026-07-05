# ─────────────────────────────────────────────────────────────────────────────
# 03 — Differential abundance: c_section vs vaginal
# Primary: CLR transform + per-taxon Wilcoxon + BH FDR (compositionally aware,
#          appropriate for MetaPhlAn relative-abundance proportions).
# (ANCOM-BC is the count-native alternative; cMD provides proportions, so a
#  CLR-based test is the methodologically clean choice here.)
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(ggplot2) })
set.seed(42)
proj <- getwd()  # set working directory to repo root before running
d <- readRDS(file.path(proj, "results", "01_processed.rds"))
ab <- d$abundance; meta <- d$meta

# CLR transform (centred log-ratio) with a small pseudocount
clr <- function(mat) {
  m <- mat + min(mat[mat > 0]) / 2          # pseudocount = half the smallest nonzero
  logm <- log(m)
  sweep(logm, 2, colMeans(logm), "-")        # subtract per-sample geometric mean (in log)
}
cl <- clr(ab)                                 # features x samples (CLR)

grp <- meta$born_method
res <- data.frame(taxon = rownames(cl),
                  clr_vaginal = rowMeans(cl[, grp == "vaginal", drop = FALSE]),
                  clr_csection = rowMeans(cl[, grp == "c_section", drop = FALSE]))
res$clr_diff <- res$clr_csection - res$clr_vaginal   # +ve = higher in C-section
res$p <- apply(cl, 1, function(x) wilcox.test(x ~ grp)$p.value)
res$fdr <- p.adjust(res$p, method = "BH")
res <- res[order(res$fdr), ]
write.csv(res, file.path(proj, "results", "diff_abundance.csv"), row.names = FALSE)

n_sig <- sum(res$fdr < 0.05)
cat("significant taxa (FDR<0.05):", n_sig, "of", nrow(res), "\n")

# volcano-style: effect (CLR diff) vs -log10 FDR
res$sig <- ifelse(res$fdr < 0.05, "FDR<0.05", "ns")
res$short <- gsub("_", " ", gsub("^[a-z]+:", "", sub(".*s__", "", res$taxon)))
top <- head(res[res$fdr < 0.05, ], 12)
top$h <- ifelse(top$clr_diff > 0, 1, 0)   # right-side labels left-anchored, left-side right-anchored
p <- ggplot(res, aes(clr_diff, -log10(fdr), colour = sig)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70") +
  geom_text(data = top, aes(label = short, hjust = h), size = 2.5, vjust = -0.7,
            colour = "grey20", show.legend = FALSE) +
  scale_colour_manual(values = c("FDR<0.05" = "#D7301F", "ns" = "grey70")) +
  scale_x_continuous(expand = expansion(mult = c(0.30, 0.30))) +
  labs(title = "Differential abundance by birth mode (CLR + Wilcoxon, BH-FDR)",
       subtitle = "x>0 = enriched in C-section; x<0 = enriched in vaginal",
       x = "CLR mean difference (C-section - vaginal)", y = "-log10(FDR)", colour = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(proj, "figures", "da_volcano.png"), p, width = 8.5, height = 5.2, dpi = 200)
cat("[03 done] results/diff_abundance.csv + figures/da_volcano.png\n")
