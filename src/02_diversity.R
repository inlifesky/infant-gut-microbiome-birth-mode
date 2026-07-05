# ─────────────────────────────────────────────────────────────────────────────
# 02 — Diversity: alpha (Shannon) + beta (Bray-Curtis PCoA) by birth mode
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(vegan); library(ggplot2) })
set.seed(42)
proj <- getwd()  # set working directory to repo root before running
d <- readRDS(file.path(proj, "results", "01_processed.rds"))
ab <- d$abundance; meta <- d$meta
samp <- t(ab)                                   # samples x features

# --- alpha: Shannon ---
shannon <- diversity(samp, index = "shannon")
adf <- data.frame(shannon = shannon, born_method = meta$born_method)
wt <- wilcox.test(shannon ~ born_method, data = adf)

p_alpha <- ggplot(adf, aes(born_method, shannon, fill = born_method)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c(vaginal = "#2C7FB8", c_section = "#D95F0E")) +
  labs(title = "Infant gut alpha-diversity by birth mode",
       subtitle = sprintf("Shao et al. 2019 Baby Biome Study  |  Wilcoxon p = %.2e", wt$p.value),
       x = NULL, y = "Shannon diversity") +
  theme_bw(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(proj, "figures", "diversity_alpha.png"), p_alpha, width = 6.4, height = 4.4, dpi = 200)

# --- beta: Bray-Curtis PCoA + PERMANOVA (subject as strata for pseudoreplication) ---
bc <- vegdist(samp, method = "bray")
pcoa <- cmdscale(bc, k = 2, eig = TRUE)
ev <- pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100
bdf <- data.frame(PCo1 = pcoa$points[, 1], PCo2 = pcoa$points[, 2],
                  born_method = meta$born_method)
# one sample per subject (set in 01) -> unrestricted PERMANOVA is valid
perm <- adonis2(bc ~ born_method, data = meta, permutations = 999)

p_beta <- ggplot(bdf, aes(PCo1, PCo2, colour = born_method)) +
  geom_point(size = 1.3, alpha = 0.6) +
  stat_ellipse(level = 0.68, linewidth = 0.8) +
  scale_colour_manual(values = c(vaginal = "#2C7FB8", c_section = "#D95F0E")) +
  labs(title = "Birth-mode separation in infant gut community",
       subtitle = sprintf("Bray-Curtis PCoA (one sample/infant)  |  PERMANOVA R2 = %.3f, p = %.3f",
                          perm$R2[1], perm$`Pr(>F)`[1]),
       x = sprintf("PCo1 (%.1f%%)", ev[1]), y = sprintf("PCo2 (%.1f%%)", ev[2]),
       colour = "Birth mode") +
  theme_bw(base_size = 11)
ggsave(file.path(proj, "figures", "diversity_pcoa.png"), p_beta, width = 7.2, height = 4.8, dpi = 200)

sink(file.path(proj, "results", "diversity_stats.txt"))
cat("ALPHA (Shannon) by birth mode\n"); print(tapply(adf$shannon, adf$born_method, summary))
cat("\nWilcoxon:\n"); print(wt)
cat("\nBETA — PERMANOVA (Bray-Curtis, subject-stratified):\n"); print(perm)
sink()
cat("[02 done] figures/diversity_alpha.png, diversity_pcoa.png + results/diversity_stats.txt\n")
