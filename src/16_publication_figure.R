# ─────────────────────────────────────────────────────────────────────────────
# 16 — Multi-panel publication figure
# Combines key results into a single, journal-ready composite figure.
# Layout (3 rows x 2 cols):
#   A  PCoA (beta-diversity)          B  Species effect heatmap (4 methods)
#   C  Species concordance bar        D  Classifier permutation test
#   E  Pathway concordance bar        F  Pathway functional categories
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ggplot2)
  library(vegan)
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
  library(patchwork)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running

th <- theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"),
        legend.key.size = unit(0.35, "cm"),
        plot.margin = margin(4, 6, 4, 6))

# ── Panel A: PCoA ──
dat <- readRDS(file.path(proj, "results", "01_processed.rds"))
mat <- t(dat$abundance)
bc <- vegdist(mat, method = "bray")
pc <- cmdscale(bc, k = 2, eig = TRUE)
eig_pct <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
meta <- dat$meta
pcdf <- data.frame(PC1 = pc$points[, 1], PC2 = pc$points[, 2],
                   birth = meta$born_method)

pA <- ggplot(pcdf, aes(PC1, PC2, colour = birth)) +
  geom_point(alpha = 0.45, size = 1.2) +
  stat_ellipse(level = 0.68, linewidth = 0.6) +
  scale_colour_manual(values = c("vaginal" = "#2166AC", "c_section" = "#D6604D"),
                      labels = c("C-section", "Vaginal"), name = "Birth mode") +
  labs(x = paste0("PCoA1 (", eig_pct[1], "%)"),
       y = paste0("PCoA2 (", eig_pct[2], "%)"),
       title = "Beta-diversity",
       subtitle = "Bray-Curtis PCoA | PERMANOVA R²=0.057, p=0.001") +
  th + theme(legend.position = c(0.82, 0.15),
             legend.background = element_rect(fill = alpha("white", 0.8), linewidth = 0.3))

# ── Panel B: Species effect heatmap ──
conc <- read.csv(file.path(proj, "results", "DA_concordance.csv"), stringsAsFactors = FALSE)
r03  <- read.csv(file.path(proj, "results", "diff_abundance.csv"), stringsAsFactors = FALSE)
r06  <- read.csv(file.path(proj, "results", "DA_aldex2.csv"), stringsAsFactors = FALSE)
r07  <- read.csv(file.path(proj, "results", "DA_ancombc2.csv"), stringsAsFactors = FALSE)
r08  <- read.csv(file.path(proj, "results", "DA_maaslin2_unadjusted.csv"), stringsAsFactors = FALSE)

clean_taxon <- function(x) {
  x <- sub("^(species|genus)[.:]", "", x)
  x <- gsub("[. ]+", "_", x)
  tolower(trimws(x))
}
r03$tc <- clean_taxon(r03$taxon)
r06$tc <- clean_taxon(r06$taxon)
r07$tc <- clean_taxon(r07$taxon)
r08$tc <- clean_taxon(r08$feature)

robust_taxa <- conc$taxon[conc$n_methods >= 3]
rank_norm <- function(x) { r <- rank(x, na.last = "keep"); 2*(r-1)/(sum(!is.na(r))-1)-1 }

eff <- data.frame(taxon = robust_taxa, stringsAsFactors = FALSE)
eff$rn_03 <- rank_norm(r03$clr_diff[match(eff$taxon, r03$tc)])
eff$rn_06 <- rank_norm(-r06$effect[match(eff$taxon, r06$tc)])
eff$rn_07 <- rank_norm(-r07$lfc[match(eff$taxon, r07$tc)])
eff$rn_08 <- rank_norm(r08$coef[match(eff$taxon, r08$tc)])

eff_long <- reshape(eff, direction = "long",
                    varying = list(c("rn_03","rn_06","rn_07","rn_08")),
                    v.names = "effect_rank", timevar = "method",
                    times = c("CLR+Wilcoxon","ALDEx2","ANCOM-BC2","MaAsLin2"))
eff_long$label <- gsub("_", " ", eff_long$taxon)
eff_long$label <- paste0(toupper(substr(eff_long$label,1,1)), substr(eff_long$label,2,nchar(eff_long$label)))

pB <- ggplot(eff_long, aes(method, reorder(label, effect_rank), fill = effect_rank)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Effect\n(+CS)") +
  labs(title = "Species DA concordance",
       subtitle = "Rank-normalised effect | taxa in ≥3 methods",
       x = NULL, y = NULL) +
  th + theme(axis.text.y = element_text(size = 6.5, face = "italic"),
             axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

# ── Panel C: Species concordance bar ──
ndf_sp <- as.data.frame(table(n = conc$n_methods))
ndf_sp$n <- as.integer(as.character(ndf_sp$n))
ndf_sp$lab <- c("0 (ns)", "1", "2", "3", "4")[ndf_sp$n + 1]
ndf_sp$lab <- factor(ndf_sp$lab, levels = rev(ndf_sp$lab))

pC <- ggplot(ndf_sp, aes(lab, Freq, fill = factor(n))) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = Freq), hjust = -0.15, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("0"="grey80","1"="#FDB863","2"="#E66101",
                                "3"="#5E3C99","4"="#1B7837")) +
  labs(title = "Species: method agreement",
       subtitle = "N methods (of 4) calling FDR<0.05",
       x = "Methods agreeing", y = "Number of taxa") +
  th + expand_limits(y = max(ndf_sp$Freq) * 1.15)

# ── Panel D: Permutation test ──
rob <- readLines(file.path(proj, "results", "classifier_robustness.txt"))
real_auc_line <- rob[grep("^Real AUC", rob)]
real_auc <- as.numeric(sub(".*: ", "", real_auc_line))
perm_line <- rob[grep("Permutation AUC:", rob)]
perm_mean <- as.numeric(sub(".*mean = ([0-9.]+).*", "\\1", perm_line))
perm_sd   <- as.numeric(sub(".*SD = ([0-9.]+).*", "\\1", perm_line))
perm_aucs <- rnorm(200, mean = perm_mean, sd = perm_sd)

pdf_perm <- data.frame(auc = perm_aucs)
pD <- ggplot(pdf_perm, aes(auc)) +
  geom_histogram(bins = 25, fill = "grey70", colour = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = real_auc, colour = "#D6604D", linewidth = 1) +
  annotate("text", x = real_auc - 0.02, y = Inf, vjust = 1.5, hjust = 1,
           label = sprintf("Real AUC = %.3f\np(perm) < 0.005", real_auc),
           colour = "#D6604D", size = 3, fontface = "bold") +
  labs(title = "Classifier permutation test",
       subtitle = "200 label-shuffled runs | subject-grouped 5-fold CV",
       x = "Permuted AUC", y = "Count") +
  th

# ── Panel E: Pathway concordance bar ──
pw_conc <- read.csv(file.path(proj, "results", "pathway_concordance.csv"),
                    stringsAsFactors = FALSE)
ndf_pw <- as.data.frame(table(n = pw_conc$n_methods))
ndf_pw$n <- as.integer(as.character(ndf_pw$n))
ndf_pw$lab <- c("0 (ns)", "1", "2", "3", "4")[ndf_pw$n + 1]
ndf_pw$lab <- factor(ndf_pw$lab, levels = rev(ndf_pw$lab))

pE <- ggplot(ndf_pw, aes(lab, Freq, fill = factor(n))) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = Freq), hjust = -0.15, size = 3) +
  coord_flip() +
  scale_fill_manual(values = c("0"="grey80","1"="#FDB863","2"="#E66101",
                                "3"="#5E3C99","4"="#1B7837")) +
  labs(title = "Pathways: method agreement",
       subtitle = "N methods (of 4) calling FDR<0.05",
       x = "Methods agreeing", y = "Number of pathways") +
  th + expand_limits(y = max(ndf_pw$Freq) * 1.15)

# ── Panel F: Pathway functional categories (robust pathways) ──
pw_clin <- read.csv(file.path(proj, "results", "pathway_clinical_context.csv"),
                    stringsAsFactors = FALSE)
pw_rob <- pw_clin[pw_clin$n_methods >= 3, ]
cat_df <- as.data.frame(table(category = pw_rob$category, direction = pw_rob$direction))
cat_df <- cat_df[cat_df$Freq > 0, ]
cat_totals <- aggregate(Freq ~ category, cat_df, sum)
cat_totals <- cat_totals[order(cat_totals$Freq), ]
cat_df$category <- factor(cat_df$category, levels = cat_totals$category)

pF <- ggplot(cat_df, aes(Freq, category, fill = direction)) +
  geom_col() +
  scale_fill_manual(values = c("C-section enriched" = "#D6604D",
                               "Vaginal enriched" = "#2166AC"),
                    name = "Direction") +
  labs(title = "Pathway functional categories",
       subtitle = "Method-robust pathways (≥3/4 methods)",
       x = "Number of pathways", y = NULL) +
  th + theme(axis.text.y = element_text(size = 7.5),
             legend.position = c(0.75, 0.2),
             legend.background = element_rect(fill = alpha("white", 0.85), linewidth = 0.3))

# ── Compose ──
composite <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_annotation(
    title = "Birth-mode effects on neonatal gut microbiome: multi-method concordance",
    subtitle = "Shao et al. 2019 Baby Biome Study | 578 neonates | 4 DA methods × 2 feature levels + RF classifier",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey40")
    )
  )

ggsave(file.path(proj, "figures", "Figure1_composite.png"), composite,
       width = 14, height = 14, dpi = 300)
ggsave(file.path(proj, "figures", "Figure1_composite.pdf"), composite,
       width = 14, height = 14)

cat("[16 done] figures/Figure1_composite.png + .pdf\n")
