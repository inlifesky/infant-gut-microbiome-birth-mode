# ─────────────────────────────────────────────────────────────────────────────
# 09 — Multi-method concordance analysis
# Compare 4 DA methods: CLR+Wilcoxon (03), ALDEx2 (06), ANCOM-BC2 (07), MaAsLin2 (08)
# This is the methodological core — which taxa are method-robust?
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(ggplot2) })
set.seed(42)
proj <- getwd()  # set working directory to repo root before running

# --- Load results from all 4 methods ---
r03 <- read.csv(file.path(proj, "results", "diff_abundance.csv"), stringsAsFactors = FALSE)
r06 <- read.csv(file.path(proj, "results", "DA_aldex2.csv"), stringsAsFactors = FALSE)
r07 <- read.csv(file.path(proj, "results", "DA_ancombc2.csv"), stringsAsFactors = FALSE)
r08 <- read.csv(file.path(proj, "results", "DA_maaslin2_unadjusted.csv"), stringsAsFactors = FALSE)

# Standardise taxon names: all methods to lowercase "genus_species" key
# 03/06/07 format: "species:Veillonella parvula" → "veillonella_parvula"
# 08 MaAsLin2 format: "species.Veillonella.parvula" → "veillonella_parvula"
clean_taxon <- function(x) {
  x <- sub("^(species|genus)[.:]", "", x)  # strip type prefix
  x <- gsub("[. ]+", "_", x)               # dots/spaces to underscores
  tolower(x)
}
r03$taxon_clean <- clean_taxon(r03$taxon)
r06$taxon_clean <- clean_taxon(r06$taxon)
r07$taxon_clean <- clean_taxon(r07$taxon)
r08$taxon_clean <- clean_taxon(r08$feature)

# Significant sets (FDR < 0.05)
sig03 <- r03$taxon_clean[r03$fdr < 0.05]
sig06 <- r06$taxon_clean[r06$wi_fdr < 0.05]
sig07 <- r07$taxon_clean[r07$fdr < 0.05]
sig08 <- r08$taxon_clean[r08$qval < 0.05]

all_taxa <- unique(c(r03$taxon_clean, r06$taxon_clean, r07$taxon_clean, r08$taxon_clean))

# Build membership matrix
mem <- data.frame(
  taxon     = all_taxa,
  CLR_Wilcoxon = as.integer(all_taxa %in% sig03),
  ALDEx2       = as.integer(all_taxa %in% sig06),
  ANCOM_BC2    = as.integer(all_taxa %in% sig07),
  MaAsLin2     = as.integer(all_taxa %in% sig08),
  stringsAsFactors = FALSE
)
mem$n_methods <- rowSums(mem[, 2:5])
mem <- mem[order(-mem$n_methods, mem$taxon), ]

# --- Summary stats ---
cat("=== Multi-method concordance summary ===\n")
cat("Total unique taxa tested:", length(all_taxa), "\n")
cat("Significant per method:\n")
cat("  CLR+Wilcoxon:", length(sig03), "\n")
cat("  ALDEx2:      ", length(sig06), "\n")
cat("  ANCOM-BC2:   ", length(sig07), "\n")
cat("  MaAsLin2:    ", length(sig08), "\n")
for (k in 4:1) {
  n <- sum(mem$n_methods >= k)
  cat(sprintf("  Significant in >= %d methods: %d taxa\n", k, n))
}
cat("\nMethod-robust taxa (all 4 methods agree):\n")
robust <- mem$taxon[mem$n_methods == 4]
cat(paste(" ", robust, collapse = "\n"), "\n")

# --- Effect direction concordance ---
# Harmonise signs so positive = C-section enriched for ALL methods.
# 03: clr_diff = csection - vaginal → +ve = CS-enriched (already correct)
# 06: ALDEx2 effect = group2(vaginal) - group1(c_section) → FLIP sign
# 07: ANCOM-BC2 lfc = vaginal as ref, CS as contrast → +ve = higher in CS?
#     Actually ANCOM-BC2 reports lfc relative to reference (vaginal), so +ve = CS-enriched
#     But empirically B.longum (vaginal-enriched) shows +1.64 → sign is FLIPPED
# 08: MaAsLin2 coef for born_methodc_section → -ve = lower in CS = vaginal-enriched (correct)
eff <- data.frame(taxon = all_taxa, stringsAsFactors = FALSE)
eff$dir_03 <- sign(r03$clr_diff[match(eff$taxon, r03$taxon_clean)])
eff$dir_06 <- sign(-r06$effect[match(eff$taxon, r06$taxon_clean)])   # FLIP
eff$dir_07 <- sign(-r07$lfc[match(eff$taxon, r07$taxon_clean)])      # FLIP
eff$dir_08 <- sign(r08$coef[match(eff$taxon, r08$taxon_clean)])

# Among robust taxa, check direction agreement
if (length(robust) > 0) {
  eff_robust <- eff[eff$taxon %in% robust, ]
  eff_robust$all_agree <- apply(eff_robust[, 2:5], 1, function(x) {
    x <- x[!is.na(x)]; length(unique(x)) == 1
  })
  cat("\nEffect direction agreement among robust taxa:\n")
  cat("  All agree:", sum(eff_robust$all_agree), "/", nrow(eff_robust), "\n")
  disagree <- eff_robust$taxon[!eff_robust$all_agree]
  if (length(disagree)) cat("  Disagree:", paste(disagree, collapse = ", "), "\n")
}

# --- Pairwise Jaccard similarity ---
methods <- list(CLR_Wilcoxon = sig03, ALDEx2 = sig06, ANCOM_BC2 = sig07, MaAsLin2 = sig08)
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
jmat <- matrix(NA, 4, 4, dimnames = list(names(methods), names(methods)))
for (i in seq_along(methods)) for (j in seq_along(methods))
  jmat[i, j] <- jaccard(methods[[i]], methods[[j]])
cat("\nPairwise Jaccard similarity (significant taxa):\n")
print(round(jmat, 3))

# --- Save full concordance table ---
write.csv(mem, file.path(proj, "results", "DA_concordance.csv"), row.names = FALSE)

# --- UpSet-style bar chart (how many taxa in each intersection) ---
mem$category <- apply(mem[, 2:5], 1, function(x) {
  nms <- c("CLR+Wil", "ALDEx2", "ANCOM", "MaAsLin2")
  if (all(x == 0)) return("none")
  paste(nms[x == 1], collapse = " ∩ ")
})
# Summarise the n_methods distribution
ndf <- as.data.frame(table(n_methods = mem$n_methods))
ndf$n_methods <- as.integer(as.character(ndf$n_methods))
ndf$label <- c("0 (not significant)", "1 method", "2 methods", "3 methods", "4 methods")[ndf$n_methods + 1]
ndf$label <- factor(ndf$label, levels = rev(ndf$label))

p1 <- ggplot(ndf, aes(label, Freq, fill = factor(n_methods))) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = Freq), hjust = -0.2, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("0"="grey80","1"="#FDB863","2"="#E66101",
                                "3"="#5E3C99","4"="#1B7837")) +
  labs(title = "Multi-method concordance: birth-mode DA",
       subtitle = "CLR+Wilcoxon × ALDEx2 × ANCOM-BC2 × MaAsLin2",
       x = NULL, y = "Number of taxa") +
  theme_bw(base_size = 12) +
  expand_limits(y = max(ndf$Freq) * 1.15)
ggsave(file.path(proj, "figures", "da_concordance.png"), p1, width = 7, height = 4, dpi = 200)

# --- Heatmap: effect size across methods for top taxa ---
# Rank-normalise effects to [-1, 1] for visual comparison
rank_norm <- function(x) { r <- rank(x, na.last = "keep"); 2 * (r - 1) / (sum(!is.na(r)) - 1) - 1 }
eff$rn_03 <- rank_norm(r03$clr_diff[match(eff$taxon, r03$taxon_clean)])
eff$rn_06 <- rank_norm(-r06$effect[match(eff$taxon, r06$taxon_clean)])   # FLIP
eff$rn_07 <- rank_norm(-r07$lfc[match(eff$taxon, r07$taxon_clean)])      # FLIP
eff$rn_08 <- rank_norm(r08$coef[match(eff$taxon, r08$taxon_clean)])

top_taxa <- mem$taxon[mem$n_methods >= 3]
if (length(top_taxa) > 25) top_taxa <- top_taxa[1:25]

eff_top <- eff[eff$taxon %in% top_taxa, ]
eff_long <- reshape(eff_top[, c("taxon","rn_03","rn_06","rn_07","rn_08")],
                    direction = "long",
                    varying = list(c("rn_03","rn_06","rn_07","rn_08")),
                    v.names = "effect_rank",
                    timevar = "method",
                    times = c("CLR+Wilcoxon","ALDEx2","ANCOM-BC2","MaAsLin2"))
eff_long$short <- gsub("s__", "", eff_long$taxon)

p2 <- ggplot(eff_long, aes(method, reorder(short, effect_rank), fill = effect_rank)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       name = "Rank-normalised\neffect (+C-sec)") +
  labs(title = "Effect-size concordance across DA methods",
       subtitle = "Taxa significant in ≥3 methods",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(proj, "figures", "da_effect_heatmap.png"), p2, width = 7, height = 8, dpi = 200)

cat("[09 done] results/DA_concordance.csv + figures/da_concordance.png + da_effect_heatmap.png\n")
