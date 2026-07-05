# ─────────────────────────────────────────────────────────────────────────────
# 08 — Differential abundance: MaAsLin2 (with gender as covariate)
# Generalised linear model with automatic normalisation + covariate adjustment.
# Unique value: can control for gender (255F / 269M / 54 NA in our data),
#   which the other 3 methods cannot do in their standard invocation.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(Maaslin2); library(ggplot2) })
set.seed(42)
proj <- getwd()  # set working directory to repo root before running
d <- readRDS(file.path(proj, "results", "01_processed.rds"))
ab <- d$abundance; meta <- d$meta

# MaAsLin2 wants samples-as-rows data frame + metadata data frame
input_data <- as.data.frame(t(ab))
input_meta <- meta

# Drop samples with missing gender for the covariate-adjusted model
has_gender <- !is.na(input_meta$gender)
cat("Samples with gender info:", sum(has_gender), "/ dropped:", sum(!has_gender), "\n")

# --- Model 1: birth_mode only (comparable to other methods) ---
outdir1 <- file.path(proj, "results", "maaslin2_birthmode")
fit1 <- Maaslin2(
  input_data      = input_data,
  input_metadata  = input_meta,
  output          = outdir1,
  fixed_effects   = "born_method",
  normalization   = "CLR",
  transform       = "NONE",
  analysis_method = "LM",
  correction      = "BH",
  min_prevalence  = 0.10,
  plot_heatmap    = FALSE,
  plot_scatter    = FALSE,
  cores           = 1
)

res1 <- fit1$results
res1 <- res1[order(res1$qval), ]
write.csv(res1, file.path(proj, "results", "DA_maaslin2_unadjusted.csv"), row.names = FALSE)
n1 <- sum(res1$qval < 0.05, na.rm = TRUE)
cat("MaAsLin2 (birth_mode only): significant =", n1, "of", length(unique(res1$feature)), "\n")

# --- Model 2: birth_mode + gender (covariate-adjusted) ---
outdir2 <- file.path(proj, "results", "maaslin2_adjusted")
fit2 <- Maaslin2(
  input_data      = input_data[has_gender, , drop = FALSE],
  input_metadata  = input_meta[has_gender, , drop = FALSE],
  output          = outdir2,
  fixed_effects   = c("born_method", "gender"),
  normalization   = "CLR",
  transform       = "NONE",
  analysis_method = "LM",
  correction      = "BH",
  min_prevalence  = 0.10,
  plot_heatmap    = FALSE,
  plot_scatter    = FALSE,
  cores           = 1
)

res2 <- fit2$results[fit2$results$metadata == "born_method", ]
res2 <- res2[order(res2$qval), ]
write.csv(res2, file.path(proj, "results", "DA_maaslin2_gender_adjusted.csv"), row.names = FALSE)
n2 <- sum(res2$qval < 0.05, na.rm = TRUE)
cat("MaAsLin2 (+ gender covariate): significant =", n2, "of", length(unique(res2$feature)), "\n")

# Gender effect summary
res_gender <- fit2$results[fit2$results$metadata == "gender", ]
n_gender <- sum(res_gender$qval < 0.05, na.rm = TRUE)
cat("Gender main effect: significant =", n_gender, "\n")
write.csv(res_gender[order(res_gender$qval), ],
          file.path(proj, "results", "DA_maaslin2_gender_effect.csv"), row.names = FALSE)

# --- Volcano plot (unadjusted model, for concordance comparison) ---
res1$sig <- ifelse(res1$qval < 0.05, "FDR<0.05", "ns")
p <- ggplot(res1, aes(coef, -log10(qval), colour = sig)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey70") +
  scale_colour_manual(values = c("FDR<0.05" = "#1B9E77", "ns" = "grey70")) +
  labs(title = "MaAsLin2 differential abundance by birth mode",
       subtitle = sprintf("CLR + LM (unadjusted) — %d / %d FDR<0.05", n1,
                          length(unique(res1$feature))),
       x = "Coefficient (C-section vs vaginal)", y = "-log10(FDR)", colour = NULL) +
  theme_bw(base_size = 12)
ggsave(file.path(proj, "figures", "da_maaslin2.png"), p, width = 8, height = 5, dpi = 200)

cat("[08 done] results/DA_maaslin2_*.csv + figures/da_maaslin2.png\n")
