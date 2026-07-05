# ─────────────────────────────────────────────────────────────────────────────
# 15 — Pathway multi-method concordance
# Mirrors Step 09 (species concordance) for HUMAnN MetaCyc pathways.
# 4 methods: CLR+Wilcoxon (05), ALDEx2 (12), ANCOM-BC2 (13), MaAsLin2 (14)
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(ggplot2) })
set.seed(42)
proj <- getwd()  # set working directory to repo root before running

# --- Load results ---
r05 <- read.csv(file.path(proj, "results", "pathway_diff.csv"), stringsAsFactors = FALSE)
r12 <- read.csv(file.path(proj, "results", "pathway_aldex2.csv"), stringsAsFactors = FALSE)
r13 <- read.csv(file.path(proj, "results", "pathway_ancombc2.csv"), stringsAsFactors = FALSE)
r14 <- read.csv(file.path(proj, "results", "pathway_maaslin2_unadjusted.csv"), stringsAsFactors = FALSE)

# Standardise pathway names: extract the MetaCyc ID (e.g. "PWY0-1296")
# which is stable across all 4 outputs, then lowercase for matching.
# Format examples:
#   05/12/13: "PWY0-1296: purine ribonucleosides degradation"
#   14 MaAsLin2: "PWY0.1296..purine.ribonucleosides.degradation"
extract_pw_id <- function(x) {
  # Try colon-delimited first (05/12/13 format)
  id <- sub(":.*$", "", x)
  # For MaAsLin2 dot format: take up to first ".." (double dot = was ": ")
  id <- sub("\\.\\..+$", "", id)
  # Normalise remaining dots/hyphens to hyphens, lowercase
  id <- gsub("\\.", "-", id)
  tolower(trimws(id))
}
r05$pw_clean <- extract_pw_id(r05$pathway)
r12$pw_clean <- extract_pw_id(r12$pathway)
r13$pw_clean <- extract_pw_id(r13$pathway)
r14$pw_clean <- extract_pw_id(r14$feature)
# Keep a readable label from r05 for display
pw_labels <- setNames(r05$pathway, r05$pw_clean)

# Significant sets
sig05 <- r05$pw_clean[r05$fdr < 0.05]
sig12 <- r12$pw_clean[r12$wi_fdr < 0.05]
sig13 <- r13$pw_clean[r13$fdr < 0.05]
sig14 <- r14$pw_clean[r14$qval < 0.05]

all_pw <- unique(c(r05$pw_clean, r12$pw_clean, r13$pw_clean, r14$pw_clean))

mem <- data.frame(
  pathway      = all_pw,
  CLR_Wilcoxon = as.integer(all_pw %in% sig05),
  ALDEx2       = as.integer(all_pw %in% sig12),
  ANCOM_BC2    = as.integer(all_pw %in% sig13),
  MaAsLin2     = as.integer(all_pw %in% sig14),
  stringsAsFactors = FALSE
)
mem$n_methods <- rowSums(mem[, 2:5])
mem <- mem[order(-mem$n_methods, mem$pathway), ]

# --- Summary ---
cat("=== Pathway multi-method concordance ===\n")
cat("Total unique pathways tested:", length(all_pw), "\n")
cat("Significant per method:\n")
cat("  CLR+Wilcoxon:", length(sig05), "\n")
cat("  ALDEx2:      ", length(sig12), "\n")
cat("  ANCOM-BC2:   ", length(sig13), "\n")
cat("  MaAsLin2:    ", length(sig14), "\n")
for (k in 4:1) cat(sprintf("  >= %d methods: %d pathways\n", k, sum(mem$n_methods >= k)))

robust_pw <- mem$pathway[mem$n_methods == 4]
cat("\nMethod-robust pathways (all 4 agree):", length(robust_pw), "\n")
if (length(robust_pw) <= 30) cat(paste(" ", robust_pw, collapse = "\n"), "\n")

# --- Effect direction concordance ---
# Harmonise: +ve = C-section enriched
# Sign conventions differ from species layer because factor reference levels changed:
# 05 clr_diff = CS - vaginal → +ve = CS-enriched (correct)
# 12 ALDEx2 effect: group ordering puts c_section first → FLIP
# 13 ANCOM-BC2 lfc: reference = c_section (alphabetical) → +ve = vaginal-enriched → ALREADY
#    means -ve = CS-enriched; no flip needed (opposite of species where ref was vaginal)
# 14 MaAsLin2 coef: reports vaginal coefficient → +ve = vaginal-enriched → FLIP
eff <- data.frame(pathway = all_pw, stringsAsFactors = FALSE)
eff$dir_05 <- sign(r05$clr_diff[match(eff$pathway, r05$pw_clean)])
eff$dir_12 <- sign(-r12$effect[match(eff$pathway, r12$pw_clean)])
eff$dir_13 <- sign(r13$lfc[match(eff$pathway, r13$pw_clean)])      # NO flip
eff$dir_14 <- sign(-r14$coef[match(eff$pathway, r14$pw_clean)])    # FLIP

if (length(robust_pw) > 0) {
  eff_robust <- eff[eff$pathway %in% robust_pw, ]
  eff_robust$all_agree <- apply(eff_robust[, 2:5], 1, function(x) {
    x <- x[!is.na(x)]; length(unique(x)) == 1
  })
  cat("Direction agreement among robust pathways:",
      sum(eff_robust$all_agree), "/", nrow(eff_robust), "\n")
}

# --- Pairwise Jaccard ---
methods <- list(CLR_Wilcoxon = sig05, ALDEx2 = sig12, ANCOM_BC2 = sig13, MaAsLin2 = sig14)
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
jmat <- matrix(NA, 4, 4, dimnames = list(names(methods), names(methods)))
for (i in seq_along(methods)) for (j in seq_along(methods))
  jmat[i, j] <- jaccard(methods[[i]], methods[[j]])
cat("\nPairwise Jaccard (significant pathways):\n")
print(round(jmat, 3))

write.csv(mem, file.path(proj, "results", "pathway_concordance.csv"), row.names = FALSE)

# --- Bar chart ---
ndf <- as.data.frame(table(n_methods = mem$n_methods))
ndf$n_methods <- as.integer(as.character(ndf$n_methods))
ndf$label <- c("0 (ns)", "1 method", "2 methods", "3 methods", "4 methods")[ndf$n_methods + 1]
ndf$label <- factor(ndf$label, levels = rev(ndf$label))

p1 <- ggplot(ndf, aes(label, Freq, fill = factor(n_methods))) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = Freq), hjust = -0.2, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c("0"="grey80","1"="#FDB863","2"="#E66101",
                                "3"="#5E3C99","4"="#1B7837")) +
  labs(title = "Pathway multi-method concordance: birth-mode DA",
       subtitle = "CLR+Wilcoxon × ALDEx2 × ANCOM-BC2 × MaAsLin2",
       x = NULL, y = "Number of pathways") +
  theme_bw(base_size = 12) +
  expand_limits(y = max(ndf$Freq) * 1.15)
ggsave(file.path(proj, "figures", "pathway_concordance.png"), p1, width = 7, height = 4, dpi = 200)

# --- Effect heatmap for robust pathways ---
rank_norm <- function(x) { r <- rank(x, na.last = "keep"); 2*(r-1)/(sum(!is.na(r))-1)-1 }
eff$rn_05 <- rank_norm(r05$clr_diff[match(eff$pathway, r05$pw_clean)])
eff$rn_12 <- rank_norm(-r12$effect[match(eff$pathway, r12$pw_clean)])
eff$rn_13 <- rank_norm(r13$lfc[match(eff$pathway, r13$pw_clean)])     # NO flip
eff$rn_14 <- rank_norm(-r14$coef[match(eff$pathway, r14$pw_clean)])   # FLIP

top_pw <- mem$pathway[mem$n_methods >= 3]
if (length(top_pw) > 30) top_pw <- top_pw[1:30]

if (length(top_pw) >= 3) {
  eff_top <- eff[eff$pathway %in% top_pw, ]
  eff_long <- reshape(eff_top[, c("pathway","rn_05","rn_12","rn_13","rn_14")],
                      direction = "long",
                      varying = list(c("rn_05","rn_12","rn_13","rn_14")),
                      v.names = "effect_rank",
                      timevar = "method",
                      times = c("CLR+Wilcoxon","ALDEx2","ANCOM-BC2","MaAsLin2"))
  # Use readable labels from r05; shorten for display
  eff_long$full <- pw_labels[eff_long$pathway]
  eff_long$full[is.na(eff_long$full)] <- eff_long$pathway[is.na(eff_long$full)]
  eff_long$short <- ifelse(nchar(eff_long$full) > 55,
                           paste0(substr(eff_long$full, 1, 53), ".."),
                           eff_long$full)

  p2 <- ggplot(eff_long, aes(method, reorder(short, effect_rank), fill = effect_rank)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         name = "Rank-norm\neffect (+CS)") +
    labs(title = "Pathway effect concordance across DA methods",
         subtitle = "Pathways significant in ≥3 methods",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          axis.text.y = element_text(size = 7))
  ggsave(file.path(proj, "figures", "pathway_effect_heatmap.png"), p2,
         width = 8, height = max(5, length(top_pw) * 0.3), dpi = 200)
}

# --- Functional category summary for robust pathways ---
# Tag broad categories based on pathway name keywords
tag_category <- function(pw) {
  pw <- tolower(pw)
  dplyr::case_when(
    grepl("purine|pyrimidine|nucleotide|nucleoside", pw) ~ "Nucleotide metabolism",
    grepl("amino.?acid|lysine|methionine|threonine|isoleucine|valine|leucine|arginine|histidine|tryptophan|serine|glycine|alanine|proline|glutam|aspart|phenyl|tyrosine|cysteine", pw) ~ "Amino acid metabolism",
    grepl("sugar|glucose|galactose|fructose|mannose|starch|glycolysis|pentose|glucuronate|hexuronate|carbohydrate|sucrose|lactose|maltose|xylose|arabinose|rhamnose", pw) ~ "Carbohydrate metabolism",
    grepl("fatty.?acid|lipid|phospholipid|sterol|isoprenoid|mevalonate|terpenoid", pw) ~ "Lipid metabolism",
    grepl("cofactor|vitamin|folate|biotin|thiamin|riboflavin|cobalamin|pantothenate|nad|coenzyme|menaquinone|ubiquinone", pw) ~ "Cofactor/vitamin biosynthesis",
    grepl("ferment|acetyl.?coa|tca|citrate|pyruvate|acetate|butyrate|propionate|scfa|mixed.?acid", pw) ~ "Central/fermentation",
    grepl("cell.?wall|peptidoglycan|lipopolysaccharide|lps|teichoic", pw) ~ "Cell wall",
    grepl("sulfur|sulfate|thio", pw) ~ "Sulfur metabolism",
    TRUE ~ "Other"
  )
}

if (length(robust_pw) > 0) {
  robust_df <- mem[mem$n_methods == 4, ]
  robust_df$category <- tag_category(robust_df$pathway)
  robust_df$dir <- ifelse(
    eff$dir_05[match(robust_df$pathway, eff$pathway)] > 0,
    "C-section enriched", "Vaginal enriched"
  )
  cat("\n--- Functional categories of method-robust pathways ---\n")
  cat_tab <- table(robust_df$category, robust_df$dir)
  print(cat_tab)
  write.csv(robust_df[, c("pathway","n_methods","category","dir")],
            file.path(proj, "results", "pathway_robust_categories.csv"), row.names = FALSE)
}

cat("\n[15 done] results/pathway_concordance.csv + figures/pathway_concordance.png\n")
