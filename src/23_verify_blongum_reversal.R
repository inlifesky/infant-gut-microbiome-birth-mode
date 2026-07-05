# ─────────────────────────────────────────────────────────────────────────────
# 23 — Verify: is the B. longum d91+ reversal a real natural recovery,
#       or an artifact of differential antibiotic exposure between birth modes?
#
# Literature warning (Mazzola 2016, Stearns 2017):
# IAP exposure → sustained B. longum reduction at 1 year. In our cohort,
# vaginal infants have HIGHER IAP exposure (27%) than CS infants (12%).
# If vaginal-IAP infants have very low B. longum, that could pull the vaginal
# group mean down at d91+ and create an APPARENT "CS catch-up" that is
# actually "vaginal-IAP suppression".
#
# Test: stratify B. longum trajectory by antibiotic exposure within each
# birth mode. If the reversal disappears in the antibiotic-naive subset,
# the reversal is an antibiotic artifact, not a CS-natural-recovery story.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
proj <- getwd()  # set working directory to repo root before running

d <- readRDS(file.path(proj, "results", "01b_processed_longitudinal.rds"))
meta <- d$meta; ab <- d$abundance
meta$age_bin <- cut(meta$infant_age, c(-1,7,30,90,Inf),
                    labels = c("d0_7","d8_30","d31_90","d91+"))

clr <- function(x) { x <- x + 1e-6; log(x) - mean(log(x)) }
clr_mat <- apply(ab, 2, clr)

species_of_interest <- c("species:Bifidobacterium longum",
                         "species:Escherichia coli")

# Build long-format per-sample data
build_df <- function(sp) {
  data.frame(
    sample_id = colnames(ab),
    clr = clr_mat[sp, ],
    born = meta$born_method,
    abx  = meta$antibiotics_current_use,
    age_bin = meta$age_bin
  )
}

results <- list()
for (sp in species_of_interest) {
  cat("\n══════════════════════════════════════════════════\n")
  cat("Species:", sp, "\n")
  cat("══════════════════════════════════════════════════\n")
  df <- build_df(sp)
  df <- df[df$age_bin %in% c("d0_7","d8_30","d91+") & !is.na(df$born) & !is.na(df$abx), ]
  df$stratum <- paste0(df$born, "_", ifelse(df$abx == "yes", "ABx", "noABx"))

  # Cell sizes
  cat("\nSample sizes per cell (stratum × age):\n")
  print(table(df$stratum, df$age_bin))

  # Means per cell
  cell_means <- df |>
    group_by(stratum, age_bin) |>
    summarise(mean_clr = mean(clr), n = n(), .groups = "drop")
  cat("\nMean CLR per cell:\n"); print(as.data.frame(cell_means))

  # The critical question:
  # Original (uncontrolled) gap at d91+:
  vag_d91 <- df$clr[df$born == "vaginal"   & df$age_bin == "d91+"]
  cs_d91  <- df$clr[df$born == "c_section" & df$age_bin == "d91+"]
  orig_gap <- mean(cs_d91) - mean(vag_d91)
  cat(sprintf("\nUncontrolled CS-vaginal gap at d91+: %+.3f CLR\n", orig_gap))

  # Antibiotic-naive subset only:
  vag_naive_d91 <- df$clr[df$born == "vaginal"   & df$abx == "no" & df$age_bin == "d91+"]
  cs_naive_d91  <- df$clr[df$born == "c_section" & df$abx == "no" & df$age_bin == "d91+"]
  naive_gap <- mean(cs_naive_d91) - mean(vag_naive_d91)
  cat(sprintf("ABx-naive only CS-vaginal gap at d91+: %+.3f CLR (n_vag=%d, n_cs=%d)\n",
              naive_gap, length(vag_naive_d91), length(cs_naive_d91)))

  # Antibiotic-exposed subset only:
  vag_abx_d91 <- df$clr[df$born == "vaginal"   & df$abx == "yes" & df$age_bin == "d91+"]
  cs_abx_d91  <- df$clr[df$born == "c_section" & df$abx == "yes" & df$age_bin == "d91+"]
  cat(sprintf("ABx-exposed only CS-vaginal gap at d91+: %+.3f CLR (n_vag=%d, n_cs=%d)\n",
              mean(cs_abx_d91) - mean(vag_abx_d91),
              length(vag_abx_d91), length(cs_abx_d91)))

  results[[sp]] <- list(uncontrolled = orig_gap, abx_naive = naive_gap,
                         cell_means = cell_means)

  # Plot trajectory by stratum
  plot_df <- df |>
    group_by(stratum, age_bin) |>
    summarise(mean_clr = mean(clr), se = sd(clr)/sqrt(n()), n = n(), .groups = "drop")
  plot_df$born <- ifelse(grepl("^c_section", plot_df$stratum), "CS", "Vaginal")
  plot_df$abx  <- ifelse(grepl("ABx$", plot_df$stratum), "+Antibiotics", "No antibiotics")

  g <- ggplot(plot_df, aes(age_bin, mean_clr, colour = born, linetype = abx, group = stratum)) +
    geom_line(linewidth = 0.9) +
    geom_point(aes(size = n)) +
    geom_errorbar(aes(ymin = mean_clr - se, ymax = mean_clr + se),
                  width = 0.1, alpha = 0.6) +
    scale_colour_manual(values = c(CS = "#D6604D", Vaginal = "#2166AC")) +
    scale_size_continuous(range = c(2, 5), name = "n") +
    labs(title = paste("Stratified trajectory:", sub("species:", "", sp)),
         subtitle = "Does the d91+ reversal survive antibiotic stratification?",
         x = "Infant age window", y = "Mean CLR") +
    theme_bw(base_size = 11)
  ggsave(file.path(proj, "figures",
                   sprintf("verify_%s_stratified.png",
                           gsub("[ :.]", "_", sub("species:", "", sp)))),
         g, width = 8, height = 5, dpi = 200)
}

# Save summary
out_summary <- do.call(rbind, lapply(names(results), function(sp) {
  data.frame(species = sub("species:", "", sp),
             uncontrolled_d91_gap = results[[sp]]$uncontrolled,
             abx_naive_d91_gap   = results[[sp]]$abx_naive)
}))
write.csv(out_summary, file.path(proj, "results", "verify_reversal_summary.csv"),
          row.names = FALSE)

cat("\n[23 done] verification + stratified figures saved\n")
