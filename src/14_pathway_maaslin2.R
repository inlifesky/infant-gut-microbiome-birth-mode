# ─────────────────────────────────────────────────────────────────────────────
# 14 — Pathway differential abundance: MaAsLin2 (± gender)
# Uses pathway matrix saved by Step 12.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(Maaslin2); library(ggplot2) })
set.seed(42)
proj <- getwd()  # set working directory to repo root before running
d <- readRDS(file.path(proj, "results", "05_pathway_processed.rds"))
M <- d$abundance; meta <- d$meta

input_data <- as.data.frame(t(M))
input_meta <- as.data.frame(meta)
rownames(input_meta) <- input_meta$sample_id

has_gender <- !is.na(input_meta$gender)
cat("Pathway samples with gender:", sum(has_gender), "/ dropped:", sum(!has_gender), "\n")

# --- Model 1: birth_mode only ---
outdir1 <- file.path(proj, "results", "pathway_maaslin2_birthmode")
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
write.csv(res1, file.path(proj, "results", "pathway_maaslin2_unadjusted.csv"), row.names = FALSE)
n1 <- sum(res1$qval < 0.05, na.rm = TRUE)
cat("Pathway MaAsLin2 (birth_mode only): significant =", n1, "of",
    length(unique(res1$feature)), "\n")

# --- Model 2: birth_mode + gender ---
outdir2 <- file.path(proj, "results", "pathway_maaslin2_adjusted")
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
write.csv(res2, file.path(proj, "results", "pathway_maaslin2_gender_adjusted.csv"), row.names = FALSE)
n2 <- sum(res2$qval < 0.05, na.rm = TRUE)
cat("Pathway MaAsLin2 (+ gender): significant =", n2, "of",
    length(unique(res2$feature)), "\n")

res_gender <- fit2$results[fit2$results$metadata == "gender", ]
n_gender <- sum(res_gender$qval < 0.05, na.rm = TRUE)
cat("Gender main effect on pathways: significant =", n_gender, "\n")

# --- Volcano (unadjusted) ---
res1$sig <- ifelse(res1$qval < 0.05, "FDR<0.05", "ns")
p <- ggplot(res1, aes(coef, -log10(qval), colour = sig)) +
  geom_point(alpha = 0.5, size = 1.2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c("FDR<0.05" = "#1B9E77", "ns" = "grey70")) +
  labs(title = "MaAsLin2: metabolic pathway DA by birth mode",
       subtitle = sprintf("CLR + LM (unadjusted) — %d / %d FDR<0.05", n1,
                          length(unique(res1$feature))),
       x = "Coefficient (C-section vs vaginal)", y = "-log10(FDR)", colour = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(proj, "figures", "pathway_maaslin2.png"), p, width = 8, height = 5, dpi = 200)

cat("[14 done] results/pathway_maaslin2_*.csv + figures/pathway_maaslin2.png\n")
