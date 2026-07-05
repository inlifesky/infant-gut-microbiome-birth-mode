# ─────────────────────────────────────────────────────────────────────────────
# 24 — Pathway-level longitudinal + antibiotic stratification
#
# Apply the steps 17/18 logic to the 13 functional categories of method-robust
# pathways. Per-category mean CLR (across pathways in that category) as the
# unit of analysis — this is cleaner than per-pathway because (a) some
# categories have only 2 pathways, (b) clinically actionable inference is at
# category level, not individual MetaCyc ID level.
#
# Two questions per category:
#   - Does the d0-7 birth-mode gap persist, attenuate, or close by d91+?
#   - At d0-7, is the gap CS-driven or antibiotic-driven?
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(curatedMetagenomicData); library(SummarizedExperiment)
  library(ggplot2); library(dplyr); library(tidyr)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running

# Get the pathway category map from step 11 output
pcat <- read.csv(file.path(proj, "results", "pathway_clinical_context.csv"),
                 stringsAsFactors = FALSE)
cat(sprintf("Method-robust pathways: %d in %d categories\n",
            nrow(pcat), length(unique(pcat$category))))

# Pull longitudinal pathway abundance (HUMAnN3 from cMD)
cat("\nFetching HUMAnN3 pathway abundance for ShaoY_2019 (large download)...\n")
pw_full <- curatedMetagenomicData("ShaoY_2019.pathway_abundance",
                                  dryrun = FALSE, counts = FALSE,
                                  rownames = "short")[[1]]
cat("Pathway TSE dim:", dim(assay(pw_full)), "\n")

# Use the longitudinal sample metadata from step 01b
d <- readRDS(file.path(proj, "results", "01b_processed_longitudinal.rds"))
meta <- d$meta
samples_keep <- intersect(meta$sample_id, colnames(pw_full))
pw <- pw_full[, samples_keep]
meta <- meta[match(samples_keep, meta$sample_id), ]
ab <- as.matrix(assay(pw))
cat(sprintf("\nPathway abundance: %d pathways x %d samples\n", nrow(ab), ncol(ab)))

# Show example rownames to debug naming
cat("First 3 pathway rownames:\n"); print(head(rownames(ab), 3))
cat("First 3 method-robust IDs:\n"); print(head(pcat$pathway, 3))

# Normalize names: pathway rownames vs pcat$pathway IDs (lowercase, hyphenated)
norm_pw <- function(x) {
  x <- sub("\\|.*$", "", x)          # strip stratification
  x <- sub(":.*$", "", x)             # strip description
  tolower(gsub("[. _]+", "-", trimws(x)))
}
rn <- rownames(ab)
ab_id <- norm_pw(rn)
pcat$id_norm <- tolower(gsub("[. _]+", "-", pcat$pathway))

# Match
matched_idx <- which(ab_id %in% pcat$id_norm)
cat(sprintf("\nMatched %d / %d method-robust pathways to longitudinal data\n",
            length(matched_idx), nrow(pcat)))

ab_m <- ab[matched_idx, ]
ab_m_norm <- ab_id[matched_idx]
# Pathway → category map
cat_map <- setNames(pcat$category, pcat$id_norm)
row_cat <- cat_map[ab_m_norm]
# Pathway → direction map (baseline)
dir_map <- setNames(pcat$direction, pcat$id_norm)
row_dir <- dir_map[ab_m_norm]

cat("\nMatched pathways per category:\n")
print(table(row_cat, row_dir))

# CLR transform (pathway proportions are also relative)
clr <- function(x) { x <- x + 1e-6; log(x) - mean(log(x)) }
clr_mat <- apply(ab_m, 2, clr)

# Aggregate per-category mean CLR, per-sample
# (only over pathways with the baseline direction this row's category)
# We compute the CLR mean separately for vaginal-enriched and CS-enriched
# pathways within each category, since they tell different stories.
meta$age_bin <- cut(meta$infant_age, c(-1,7,30,90,Inf),
                    labels = c("d0_7","d8_30","d31_90","d91+"))

# Helper: per-sample category mean
build_cat_means <- function(category, direction) {
  rows <- which(row_cat == category & row_dir == direction)
  if (length(rows) == 0) return(NULL)
  colMeans(clr_mat[rows, , drop = FALSE])
}

categories <- unique(row_cat)
sub_long <- list()
for (ct in categories) for (dr in unique(row_dir)) {
  vals <- build_cat_means(ct, dr)
  if (is.null(vals)) next
  sub_long[[paste(ct, dr, sep = " | ")]] <-
    data.frame(category = ct, direction_baseline = dr,
               sample = names(vals), cat_clr = vals,
               stringsAsFactors = FALSE)
}
long <- do.call(rbind, sub_long)
long$born <- meta$born_method[match(long$sample, meta$sample_id)]
long$age_bin <- as.character(meta$age_bin[match(long$sample, meta$sample_id)])
long$abx <- meta$antibiotics_current_use[match(long$sample, meta$sample_id)]
long <- long[long$age_bin %in% c("d0_7","d8_30","d91+") & !is.na(long$born), ]

# ── Trajectory analysis ──
traj <- long |>
  group_by(category, direction_baseline, age_bin, born) |>
  summarise(mean_clr = mean(cat_clr), n = n(), .groups = "drop") |>
  pivot_wider(names_from = born, values_from = c(mean_clr, n)) |>
  mutate(gap = mean_clr_c_section - mean_clr_vaginal)

traj_summary <- traj |>
  filter(age_bin %in% c("d0_7", "d91+")) |>
  select(category, direction_baseline, age_bin, gap) |>
  pivot_wider(names_from = age_bin, values_from = gap) |>
  rename(gap_d0_7 = d0_7, gap_d91p = `d91+`) |>
  mutate(attenuation_pct = 100 * (1 - abs(gap_d91p)/pmax(abs(gap_d0_7), 1e-6)),
         class = case_when(
           is.na(gap_d91p) | is.na(gap_d0_7) ~ "insufficient_data",
           sign(gap_d0_7) != sign(gap_d91p) & abs(gap_d91p) < 0.3 * abs(gap_d0_7) ~ "closed_or_reversed",
           attenuation_pct >= 50 ~ "attenuated",
           attenuation_pct >= 0 ~ "persists",
           TRUE ~ "widened"))
cat("\n=== Pathway category trajectory classification ===\n")
print(as.data.frame(traj_summary))
write.csv(traj_summary,
          file.path(proj, "results", "pathway_category_trajectory.csv"),
          row.names = FALSE)

# ── Antibiotic deconvolution (d0-7 only) ──
abx_long <- long |> filter(age_bin == "d0_7", !is.na(abx))
deconv <- abx_long |>
  group_by(category, direction_baseline, born, abx) |>
  summarise(mean_clr = mean(cat_clr), n = n(), .groups = "drop") |>
  pivot_wider(names_from = c(born, abx), values_from = c(mean_clr, n)) |>
  mutate(cs_effect  = (mean_clr_c_section_no + mean_clr_c_section_yes)/2 -
                      (mean_clr_vaginal_no   + mean_clr_vaginal_yes)/2,
         abx_effect = (mean_clr_c_section_yes + mean_clr_vaginal_yes)/2 -
                      (mean_clr_c_section_no  + mean_clr_vaginal_no )/2,
         dominant = case_when(
           abs(cs_effect)  > 2 * abs(abx_effect) ~ "CS-dominant",
           abs(abx_effect) > 2 * abs(cs_effect)  ~ "Antibiotic-dominant",
           TRUE                                  ~ "Both contribute"))
cat("\n=== Pathway category antibiotic-vs-CS dominance ===\n")
print(as.data.frame(deconv[, c("category","direction_baseline","cs_effect","abx_effect","dominant")]))
write.csv(deconv,
          file.path(proj, "results", "pathway_category_antibiotic_deconv.csv"),
          row.names = FALSE)

# ── Plot: trajectory by category ──
plt <- traj |> filter(age_bin %in% c("d0_7","d8_30","d91+"))
plt$age_bin <- factor(plt$age_bin, levels = c("d0_7","d8_30","d91+"))
plt$strip <- paste0(plt$category, "\n(", plt$direction_baseline, ")")
g <- ggplot(plt, aes(age_bin, gap, group = strip)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(aes(colour = direction_baseline), linewidth = 0.7) +
  geom_point(aes(colour = direction_baseline), size = 1.7) +
  scale_colour_manual(values = c(`CS-section enriched`="#D6604D",
                                  `C-section enriched`="#D6604D",
                                  `Vaginal enriched`="#2166AC"), name="") +
  facet_wrap(~ strip, scales = "free_y", ncol = 4) +
  labs(title = "Pathway-category trajectory: CS-vaginal gap across age windows",
       subtitle = "Mean CLR per category × direction-of-baseline, per age window",
       x = "Infant age window", y = "Mean category CLR(CS) − CLR(vaginal)") +
  theme_bw(base_size = 10) + theme(strip.text = element_text(size = 7))
ggsave(file.path(proj, "figures", "pathway_category_trajectory.png"),
       g, width = 14, height = 10, dpi = 200)

cat("\n[24 done] pathway category trajectory + deconvolution + figure\n")
