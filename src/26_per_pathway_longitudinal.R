# ─────────────────────────────────────────────────────────────────────────────
# 26 — Per-pathway longitudinal trajectory (NOT category-aggregated)
#
# Step 24 aggregated pathways to category level, which diluted high-magnitude
# individual pathway signals (e.g. Bifidobacterium shunt -2.68 CLR got washed
# out in the Carbohydrate-metabolism category mean of -0.086). To honestly
# answer "do functional differences persist?", we need to track INDIVIDUAL
# method-robust pathways across age windows.
#
# Output:
#   - per-pathway trajectory table (gap at d0-7, d8-30, d91+; trajectory class)
#   - top persisting/widening pathways with their categories and direction
#   - cross-reference: do persistent pathways belong to functions of the
#     persistent Bacteroidota species?
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(curatedMetagenomicData); library(SummarizedExperiment)
  library(ggplot2); library(dplyr); library(tidyr)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running

# Load the same pathway TSE as step 24
pw_full <- curatedMetagenomicData("ShaoY_2019.pathway_abundance",
                                  dryrun = FALSE, counts = FALSE,
                                  rownames = "short")[[1]]

d <- readRDS(file.path(proj, "results", "01b_processed_longitudinal.rds"))
meta <- d$meta
samples_keep <- intersect(meta$sample_id, colnames(pw_full))
pw <- pw_full[, samples_keep]
meta <- meta[match(samples_keep, meta$sample_id), ]
ab <- as.matrix(assay(pw))

# Strip unintegrated / species-stratified rows; keep only community-level pathways
rn <- rownames(ab)
keep_rows <- !grepl("\\|", rn)         # remove species-stratified entries
ab <- ab[keep_rows, ]
cat(sprintf("Community-level pathways: %d × %d samples\n", nrow(ab), ncol(ab)))

# Match to method-robust pathway IDs (from step 11 output)
pcat <- read.csv(file.path(proj, "results", "pathway_clinical_context.csv"),
                 stringsAsFactors = FALSE)
# normalize both sides to MetaCyc ID lowercase
norm_pw <- function(x) {
  x <- sub(":.*$", "", x); tolower(gsub("[. _]+", "-", trimws(x)))
}
rn_id  <- norm_pw(rownames(ab))
pcat$id_norm <- tolower(gsub("[. _]+", "-", pcat$pathway))
matched_idx <- which(rn_id %in% pcat$id_norm)
ab_m <- ab[matched_idx, ]
rn_m_id <- rn_id[matched_idx]
cat(sprintf("Matched %d / %d method-robust pathways\n", length(matched_idx), nrow(pcat)))

# CLR
clr <- function(x) { x <- x + 1e-6; log(x) - mean(log(x)) }
clr_mat <- apply(ab_m, 2, clr)

# Age window
meta$age_bin <- cut(meta$infant_age, c(-1,7,30,90,Inf),
                    labels = c("d0_7","d8_30","d31_90","d91+"))
keep_idx <- which(meta$age_bin %in% c("d0_7","d8_30","d91+") &
                  !is.na(meta$born_method))
mk <- meta[keep_idx, ]; clr_k <- clr_mat[, keep_idx]

# Per-pathway per-window CS-vaginal gap
windows <- c("d0_7","d8_30","d91+")
res <- expand.grid(pw_id = rn_m_id, window = windows, stringsAsFactors = FALSE)
res$gap <- NA_real_; res$mean_cs <- NA_real_; res$mean_vag <- NA_real_
res$n_cs <- NA_integer_; res$n_vag <- NA_integer_; res$p <- NA_real_
pw_to_row <- setNames(seq_along(rn_m_id), rn_m_id)
for (i in seq_len(nrow(res))) {
  row_i <- pw_to_row[res$pw_id[i]]
  w <- res$window[i]
  mask <- mk$age_bin == w
  v <- clr_k[row_i, mask & mk$born_method == "vaginal"]
  c_ <- clr_k[row_i, mask & mk$born_method == "c_section"]
  res$n_vag[i] <- length(v); res$n_cs[i] <- length(c_)
  res$mean_vag[i] <- mean(v); res$mean_cs[i] <- mean(c_)
  res$gap[i]     <- mean(c_) - mean(v)
  res$p[i] <- tryCatch(wilcox.test(c_, v)$p.value, error = function(e) NA)
}
res$q <- p.adjust(res$p, method = "BH")

# Wide format: gap per window
wide <- res |>
  select(pw_id, window, gap) |>
  pivot_wider(names_from = window, values_from = gap, names_prefix = "gap_") |>
  rename(gap_d0_7 = `gap_d0_7`, gap_d8_30 = `gap_d8_30`, gap_d91p = `gap_d91+`)

# Attach category + direction
wide$category  <- pcat$category[match(wide$pw_id, pcat$id_norm)]
wide$direction <- pcat$direction[match(wide$pw_id, pcat$id_norm)]
wide$readable  <- pcat$readable[match(wide$pw_id, pcat$id_norm)]
wide$n_methods <- pcat$n_methods[match(wide$pw_id, pcat$id_norm)]

# Trajectory class (same logic as step 17)
wide$attenuation_pct <- 100 * (1 - abs(wide$gap_d91p) / pmax(abs(wide$gap_d0_7), 1e-9))
wide$class <- with(wide, ifelse(
  is.na(gap_d91p) | is.na(gap_d0_7), "insufficient_data",
  ifelse(sign(gap_d0_7) != sign(gap_d91p) & abs(gap_d91p) < 0.3 * abs(gap_d0_7),
         "closed_or_reversed",
  ifelse(attenuation_pct >= 50, "attenuated",
  ifelse(attenuation_pct >= 0,  "persists", "widened")))))

cat("\n=== Per-pathway trajectory class distribution ===\n")
print(table(wide$class))
cat("\n=== By baseline direction ===\n")
print(table(wide$class, wide$direction))

# Highlight persisting/widening with substantial magnitude at d91+
substantial <- wide |>
  filter(class %in% c("persists","widened"),
         abs(gap_d91p) >= 0.5)
cat(sprintf("\n=== Persisting/widening pathways with |d91+ gap| >= 0.5 CLR (n=%d) ===\n",
            nrow(substantial)))
print(substantial[order(-abs(substantial$gap_d91p)),
                  c("readable","direction","category","gap_d0_7","gap_d91p","class")])

# Save
write.csv(wide |> select(pw_id, readable, category, direction, n_methods,
                          gap_d0_7, gap_d8_30, gap_d91p, attenuation_pct, class),
          file.path(proj, "results", "per_pathway_trajectory.csv"),
          row.names = FALSE)
write.csv(substantial,
          file.path(proj, "results", "persistent_pathways_substantial.csv"),
          row.names = FALSE)

# ── Figure: trajectory by class + magnitude ──
plt <- res
plt$window <- factor(plt$window, levels = windows)
plt$direction <- pcat$direction[match(plt$pw_id, pcat$id_norm)]
plt$category  <- pcat$category[match(plt$pw_id, pcat$id_norm)]
# Class for line styling
plt$class <- wide$class[match(plt$pw_id, wide$pw_id)]

# Focus plot: only pathways with |gap_d0_7| >= 1 to reduce clutter
substantial_d07 <- wide$pw_id[abs(wide$gap_d0_7) >= 1]
plt_sub <- plt |> filter(pw_id %in% substantial_d07)

g <- ggplot(plt_sub, aes(window, gap, group = pw_id, colour = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  geom_point(alpha = 0.6, size = 1.3) +
  scale_colour_manual(values = c(`CS-section enriched`="#D6604D",
                                  `C-section enriched`="#D6604D",
                                  `Vaginal enriched`="#2166AC"), name="") +
  labs(title = sprintf("Per-pathway longitudinal trajectory (n=%d with |d0-7 gap|≥1 CLR)",
                       length(substantial_d07)),
       subtitle = "Individual method-robust pathways — NOT category-aggregated",
       x = "Infant age window", y = "CLR(CS) − CLR(vaginal)") +
  theme_bw(base_size = 11)
ggsave(file.path(proj, "figures", "per_pathway_trajectory.png"),
       g, width = 10, height = 7, dpi = 200)

# Class summary by direction (matches step 17 output style)
cls_tab <- as.data.frame(table(class = wide$class, direction = wide$direction))
cls_tab <- cls_tab[cls_tab$Freq > 0, ]
g2 <- ggplot(cls_tab, aes(class, Freq, fill = direction)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = Freq), position = position_dodge(0.9),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(`CS-section enriched`="#D6604D",
                                `C-section enriched`="#D6604D",
                                `Vaginal enriched`="#2166AC")) +
  labs(title = "Per-pathway trajectory class × baseline direction",
       subtitle = sprintf("All %d method-robust pathways (≥3/4 methods)", nrow(wide)),
       x = NULL, y = "Number of pathways") +
  theme_bw(base_size = 11)
ggsave(file.path(proj, "figures", "per_pathway_class_summary.png"),
       g2, width = 9, height = 5, dpi = 200)

cat("\n[26 done] per_pathway_trajectory.csv + 2 figures\n")
