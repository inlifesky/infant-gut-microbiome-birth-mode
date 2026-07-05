# ─────────────────────────────────────────────────────────────────────────────
# 17 — Longitudinal recovery: does the CS-vaginal gap close?
#
# Clinical question: for the 27 method-robust species (step 09 ≥3/4 methods),
# does the depletion/enrichment seen at day 4-7 persist, attenuate, or close
# by day 91+? This directly determines the INTERVENTION WINDOW — if Bifido
# recovers naturally by 3 months, probiotic supplementation should target the
# first 2 months. If it persists, longer intervention is warranted.
#
# Design:
#   - Three natural age windows: d0-7 (n=827), d8-30 (n=349), d91+ (n=291)
#   - For each robust species, compute CLR(CS) - CLR(vaginal) per window
#   - Plot trajectory; flag species that close (>50% attenuation by d91+)
#     vs persist vs widen
#
# Subjects contribute multiple samples but each sample × species pair is
# independent within a window — we report effects within window, not mixed-model
# subject-level slopes (BBS has only 2-3 samples per subject in different
# windows, so per-subject slopes are noisy). The within-window cross-sectional
# effect IS the clinically actionable quantity.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running

d <- readRDS(file.path(proj, "results", "01b_processed_longitudinal.rds"))
meta <- d$meta
ab   <- d$abundance  # species × sample (proportions, 830 × 1470)

# Age windows + sample filter
meta$age_bin <- cut(meta$infant_age, breaks = c(-1, 7, 30, 90, Inf),
                    labels = c("d0_7", "d8_30", "d31_90", "d91+"))
keep_idx <- which(meta$age_bin %in% c("d0_7", "d8_30", "d91+") &
                  !is.na(meta$born_method))
meta_k <- meta[keep_idx, ]
ab_k   <- ab[, keep_idx]

cat("Samples per (age window × birth mode):\n")
print(table(meta_k$age_bin, meta_k$born_method))

# CLR transform (add small pseudocount; proportions in [0,1])
clr <- function(x) {
  x <- x + 1e-6
  log(x) - mean(log(x))
}
clr_mat <- apply(ab_k, 2, clr)  # species × sample, CLR-transformed

# Robust species from step 09 (≥3/4 methods)
conc <- read.csv(file.path(proj, "results", "DA_concordance.csv"),
                 stringsAsFactors = FALSE)
robust <- conc[conc$n_methods >= 3, "taxon"]

# Match robust names to TSE rownames. Longitudinal: "species:Bifidobacterium longum".
# Step09: "bifidobacterium_longum". Normalise both to "bifidobacterium_longum".
sp_names <- rownames(ab_k)
norm <- function(x) {
  x <- sub("^species[: ]+", "", x)
  x <- sub("^genus[: ]+", "", x)
  x <- gsub("\\[|\\]", "", x)      # strip MetaPhlAn brackets like [Eubacterium]
  tolower(gsub("[. _ ]+", "_", trimws(x)))
}
sp_norm <- norm(sp_names)
rob_norm <- norm(robust)
matched <- sp_names[match(rob_norm, sp_norm)]
matched <- matched[!is.na(matched)]
cat(sprintf("\nMatched %d / %d robust species to longitudinal abundance table\n",
            length(matched), length(robust)))

# Per-window CS-vaginal effect, per species
windows <- c("d0_7", "d8_30", "d91+")
res <- expand.grid(species = matched, window = windows,
                   stringsAsFactors = FALSE)
res$clr_cs <- NA_real_; res$clr_vag <- NA_real_
res$gap <- NA_real_; res$se <- NA_real_
res$p <- NA_real_; res$n_cs <- NA_integer_; res$n_vag <- NA_integer_

for (i in seq_len(nrow(res))) {
  sp <- res$species[i]; w <- res$window[i]
  idx <- meta_k$age_bin == w
  vals_cs  <- clr_mat[sp, idx & meta_k$born_method == "c_section"]
  vals_vag <- clr_mat[sp, idx & meta_k$born_method == "vaginal"]
  res$n_cs[i]   <- length(vals_cs)
  res$n_vag[i]  <- length(vals_vag)
  res$clr_cs[i]  <- mean(vals_cs)
  res$clr_vag[i] <- mean(vals_vag)
  res$gap[i]     <- mean(vals_cs) - mean(vals_vag)  # +ve = CS-enriched
  pooled_sd <- sqrt(var(vals_cs)/length(vals_cs) + var(vals_vag)/length(vals_vag))
  res$se[i] <- pooled_sd
  res$p[i]  <- tryCatch(wilcox.test(vals_cs, vals_vag)$p.value,
                        error = function(e) NA_real_)
}
res$q <- p.adjust(res$p, method = "BH")

# Trajectory classification: per species, gap at d0_7 vs gap at d91+
traj <- res |>
  select(species, window, gap) |>
  pivot_wider(names_from = window, values_from = gap,
              names_prefix = "gap_") |>
  rename(gap_d0_7 = `gap_d0_7`, gap_d8_30 = `gap_d8_30`, gap_d91p = `gap_d91+`)
traj$attenuation_pct <- 100 * (1 - abs(traj$gap_d91p) / pmax(abs(traj$gap_d0_7), 1e-9))
traj$class <- with(traj, ifelse(
  is.na(gap_d91p) | is.na(gap_d0_7), "insufficient_data",
  ifelse(sign(gap_d0_7) != sign(gap_d91p) & abs(gap_d91p) < 0.3 * abs(gap_d0_7),
         "closed_or_reversed",
  ifelse(attenuation_pct >= 50, "attenuated",
  ifelse(attenuation_pct >= 0,  "persists",
                                "widened")))))

cat("\nTrajectory class counts:\n")
print(table(traj$class))

# Attach direction (from step 09 baseline)
r03 <- read.csv(file.path(proj, "results", "diff_abundance.csv"),
                stringsAsFactors = FALSE)
clean_taxon <- function(x) {
  x <- sub("^species[:.]?", "", x); x <- sub("^genus[:.]?", "", x)
  tolower(gsub("[. ]+", "_", x))
}
r03$tc <- clean_taxon(r03$taxon)
# norm() gives bifidobacterium_longum; clean_taxon() gives same
traj$baseline_direction <- ifelse(
  r03$clr_diff[match(norm(traj$species), r03$tc)] > 0,
  "CS-enriched", "Vaginal-enriched"
)
# Fallback for species not in r03 (use within-window gap)
miss <- is.na(traj$baseline_direction)
traj$baseline_direction[miss] <- ifelse(
  traj$gap_d0_7[miss] > 0, "CS-enriched", "Vaginal-enriched")

write.csv(res,  file.path(proj, "results", "longitudinal_effects.csv"), row.names = FALSE)
write.csv(traj, file.path(proj, "results", "longitudinal_trajectories.csv"), row.names = FALSE)

# ── Plot trajectory ──
plt <- res |>
  mutate(window = factor(window, levels = windows),
         label  = gsub("_", " ", species))
plt$label <- paste0(toupper(substr(plt$label,1,1)), substr(plt$label,2,nchar(plt$label)))
plt$species_dir <- traj$baseline_direction[match(plt$species, traj$species)]

g <- ggplot(plt, aes(window, gap, group = species, colour = species_dir)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(alpha = 0.4, linewidth = 0.5) +
  geom_point(aes(size = -log10(pmax(q, 1e-10))), alpha = 0.7) +
  scale_colour_manual(values = c("CS-enriched" = "#D6604D",
                                  "Vaginal-enriched" = "#2166AC"),
                      name = "Baseline direction") +
  scale_size_continuous(name = "-log10(q)", range = c(0.5, 4)) +
  labs(x = "Infant age window",
       y = "CLR(CS) - CLR(vaginal)",
       title = "Longitudinal trajectory of method-robust species",
       subtitle = sprintf("%d species across 3 age windows; +y = CS-enriched, -y = vaginal-enriched",
                          length(unique(plt$species)))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "right")

ggsave(file.path(proj, "figures", "longitudinal_trajectories.png"),
       g, width = 9, height = 6, dpi = 200)

# Per-class summary plot
cls_tab <- as.data.frame(table(class = traj$class,
                               direction = traj$baseline_direction))
cls_tab <- cls_tab[cls_tab$Freq > 0, ]
g2 <- ggplot(cls_tab, aes(class, Freq, fill = direction)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = Freq), position = position_dodge(0.9), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("CS-enriched" = "#D6604D",
                                "Vaginal-enriched" = "#2166AC")) +
  labs(title = "Trajectory class × baseline direction",
       subtitle = "Does the day 0-7 gap close, persist, or widen by day 91+?",
       x = NULL, y = "Number of species") +
  theme_bw(base_size = 11)
ggsave(file.path(proj, "figures", "longitudinal_class_summary.png"),
       g2, width = 8, height = 5, dpi = 200)

cat("\n[17 done] longitudinal_effects.csv + trajectories.csv + 2 figures\n")
