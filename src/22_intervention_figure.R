# ─────────────────────────────────────────────────────────────────────────────
# 22 — Multi-panel publication figure, REORGANISED BY INTERVENTION QUESTION
#
# Each panel answers one clinical decision question, not one analysis method.
# Layout (2 rows × 3 cols):
#   A  Is there ANY signal?                  → PCoA
#   B  Which species are real?               → DA-RF triangulation
#   C  Does the gap close by 3 months?       → Longitudinal trajectories
#   D  Is it CS or antibiotics?              → CS-effect vs ABx-effect
#   E  Elective vs emergency — same?         → CS subtype scatter
#   F  The actionable target list            → Stratification funnel
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ggplot2); library(vegan); library(patchwork); library(dplyr)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running

th <- theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey40"),
        legend.key.size = unit(0.35, "cm"),
        plot.margin = margin(4, 6, 4, 6))

# ── Panel A: PCoA (Is there ANY signal?) ──
d1 <- readRDS(file.path(proj, "results", "01_processed.rds"))
bc <- vegdist(t(d1$abundance), method = "bray")
pc <- cmdscale(bc, k = 2, eig = TRUE)
eig_pct <- round(100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0]), 1)
pcdf <- data.frame(PC1 = pc$points[, 1], PC2 = pc$points[, 2],
                   birth = d1$meta$born_method)
pA <- ggplot(pcdf, aes(PC1, PC2, colour = birth)) +
  geom_point(alpha = 0.45, size = 1.1) +
  stat_ellipse(level = 0.68, linewidth = 0.5) +
  scale_colour_manual(values = c(vaginal="#2166AC", c_section="#D6604D"),
                      labels = c("C-section","Vaginal"), name="") +
  labs(x = paste0("PCoA1 (", eig_pct[1], "%)"),
       y = paste0("PCoA2 (", eig_pct[2], "%)"),
       title = "A. Is there any birth-mode signal?",
       subtitle = "Bray-Curtis | PERMANOVA R²=0.057, p=0.001") +
  th + theme(legend.position = c(0.85, 0.12))

# ── Panel B: DA-RF triangulation (Which species are real?) ──
tri <- read.csv(file.path(proj, "results", "triangulation_table.csv"), stringsAsFactors=FALSE)
tri$rank_top15 <- pmin(tri$rf_rank, 16); tri$rank_top15[is.na(tri$rank_top15)] <- 16
tri$label <- gsub("_", " ", tri$taxon)
tri$label <- paste0(toupper(substr(tri$label,1,1)), substr(tri$label,2,nchar(tri$label)))
pB <- ggplot(tri[!is.na(tri$rf_imp), ],
             aes(da_methods, rank_top15, colour = direction)) +
  geom_hline(yintercept = 15.5, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 2.5,  linetype = "dashed", colour = "grey50") +
  geom_jitter(aes(size = rf_imp), width = 0.12, height = 0.2, alpha = 0.75) +
  scale_y_reverse(breaks = c(1,5,10,15,16), labels = c("1","5","10","15",">15")) +
  scale_x_continuous(breaks = 0:4) +
  scale_colour_manual(values = c(`CS-enriched`="#D6604D", `Vaginal-enriched`="#2166AC"), name="") +
  scale_size_continuous(range = c(0.6, 4), guide = "none") +
  labs(title = "B. Which species are real?",
       subtitle = "Top-right = supported by both DA concordance + RF",
       x = "DA methods agreeing (of 4)", y = "RF rank") +
  th + theme(legend.position = c(0.18, 0.12))

# ── Panel C: Longitudinal trajectories (Does the gap close?) ──
le <- read.csv(file.path(proj, "results", "longitudinal_effects.csv"), stringsAsFactors=FALSE)
le$window <- factor(le$window, levels = c("d0_7","d8_30","d91+"))
tr <- read.csv(file.path(proj, "results", "longitudinal_trajectories.csv"), stringsAsFactors=FALSE)
le$dir <- tr$baseline_direction[match(le$species, tr$species)]
pC <- ggplot(le, aes(window, gap, group = species, colour = dir)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_colour_manual(values = c(`CS-enriched`="#D6604D", `Vaginal-enriched`="#2166AC"), name="") +
  labs(title = "C. Does the d0-7 gap close by 3 months?",
       subtitle = "Most CS-enriched species attenuate; persistent deficits are Bacteroidota-centred",
       x = "Infant age window", y = "CLR(CS) − CLR(vaginal)") +
  th + theme(legend.position = "none")

# ── Panel D: Antibiotic vs CS effects (Which factor drives each?) ──
ax <- read.csv(file.path(proj, "results", "antibiotic_deconvolution.csv"), stringsAsFactors=FALSE)
ax$label <- gsub("species:|genus:", "", ax$species)
pD <- ggplot(ax, aes(cs_effect, abx_effect, colour = dominant)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = c(`CS-dominant`="#D6604D",
                                  `Antibiotic-dominant`="#9970AB",
                                  `Both contribute`="#5AAE61"), name="") +
  labs(title = "D. Is it CS or peripartum antibiotics?",
       subtitle = "5 species antibiotic-dominant → stewardship target, not probiotic",
       x = "CS effect (CLR)", y = "Antibiotic effect (CLR)") +
  th + theme(legend.position = c(0.78, 0.18),
             legend.text = element_text(size = 7))

# ── Panel E: Elective vs emergency CS ──
ct <- read.csv(file.path(proj, "results", "cs_type_stratification.csv"), stringsAsFactors=FALSE)
pE <- ggplot(ct, aes(elective_vs_vag, emergency_vs_vag, colour = pattern)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(size = 2, alpha = 0.8) +
  scale_colour_manual(values = c(Similar="grey55",
                                  Quantitative_difference="#E08214",
                                  Differential="#762A83"), name="") +
  labs(title = "E. Elective vs emergency CS — same disruption?",
       subtitle = "20/27 species sit on the diagonal → CS subtype does not matter",
       x = "Elective CS vs vaginal (CLR)",
       y = "Emergency CS vs vaginal (CLR)") +
  th + theme(legend.position = c(0.22, 0.85),
             legend.text = element_text(size = 7))

# ── Panel F: Stratification funnel ──
funnel <- data.frame(
  stage = factor(c("Single-method DA\n(any of 4)",
                    "Method-robust\n(≥3 of 4)",
                    "Persists / widens\nby d91+",
                    "CS-driven\n(not ABx)",
                    "Vaginal-enriched\n(restore target)"),
                  levels = c("Single-method DA\n(any of 4)",
                             "Method-robust\n(≥3 of 4)",
                             "Persists / widens\nby d91+",
                             "CS-driven\n(not ABx)",
                             "Vaginal-enriched\n(restore target)")),
  n = c(33, 27, 8, 6, 5))
pF <- ggplot(funnel, aes(n, stage, fill = stage)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), hjust = -0.25, size = 4, fontface = "bold") +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = c("grey80","#FDB863","#E08214","#D6604D","#762A83")) +
  labs(title = "F. The actionable target list",
       subtitle = "Successive stratifications: 33 → 5 species (4 Bacteroidota + C. aerofaciens)",
       x = "Number of species", y = NULL) +
  expand_limits(x = max(funnel$n) * 1.15) + th

# ── Compose ──
composite <- (pA | pB | pC) / (pD | pE | pF) +
  plot_annotation(
    title = "Which CS-born infants should get which microbiome intervention, and when?",
    subtitle = "Shao 2019 Baby Biome Study | 1470 samples, 579 infants, longitudinal d0-d428",
    theme = theme(plot.title = element_text(size = 14, face = "bold"),
                  plot.subtitle = element_text(size = 10, colour = "grey40"))
  )

ggsave(file.path(proj, "figures", "Figure1_intervention_composite.png"),
       composite, width = 16, height = 10, dpi = 200)
ggsave(file.path(proj, "figures", "Figure1_intervention_composite.pdf"),
       composite, width = 16, height = 10)

cat("[22 done] figures/Figure1_intervention_composite.{png,pdf}\n")
