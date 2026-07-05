# ─────────────────────────────────────────────────────────────────────────────
# 18 — Antibiotic vs CS effect deconvolution
#
# Clinical question: is the "CS effect" actually a peripartum/postnatal
# antibiotic effect? If yes, the intervention target is antibiotic stewardship,
# not delivery-mode-specific probiotics. If no, the CS effect is independent
# and probiotic interventions remain warranted regardless of antibiotic policy.
#
# Design: 2x2 stratification (CS × antibiotics), restricted to neonatal window
# (d0-7) where most clinical decisions are made. For each method-robust species,
# compute CLR mean by cell, then estimate independent CS and antibiotic effects.
#
# Note on coverage: vaginal infants actually have HIGHER antibiotic rate (27%)
# than CS (12%) in this cohort — peripartum antibiotic prophylaxis differs from
# what one might naively assume.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running

d <- readRDS(file.path(proj, "results", "01b_processed_longitudinal.rds"))
meta <- d$meta; ab <- d$abundance

# Restrict to neonatal window with both birth mode and antibiotic data
meta$age_bin <- cut(meta$infant_age, breaks = c(-1, 7, 30, 90, Inf),
                    labels = c("d0_7", "d8_30", "d31_90", "d91+"))
keep <- which(meta$age_bin == "d0_7" &
              !is.na(meta$born_method) &
              !is.na(meta$antibiotics_current_use))
mk <- meta[keep, ]; abk <- ab[, keep]

mk$cell <- paste0(ifelse(mk$born_method == "c_section", "CS", "VAG"), "_",
                  ifelse(mk$antibiotics_current_use == "yes", "ABx", "noABx"))
cat("Cell sizes (d0-7 only):\n")
print(table(mk$cell))

# CLR
clr <- function(x) { x <- x + 1e-6; log(x) - mean(log(x)) }
clr_mat <- apply(abk, 2, clr)

# Robust species (matched to longitudinal namespace)
conc <- read.csv(file.path(proj, "results", "DA_concordance.csv"),
                 stringsAsFactors = FALSE)
robust <- conc[conc$n_methods >= 3, "taxon"]
norm <- function(x) {
  x <- sub("^species[: ]+", "", x); x <- sub("^genus[: ]+", "", x)
  x <- gsub("\\[|\\]", "", x); tolower(gsub("[. _ ]+", "_", trimws(x)))
}
sp_names <- rownames(abk)
matched <- sp_names[match(norm(robust), norm(sp_names))]
matched <- matched[!is.na(matched)]
cat(sprintf("Matched %d / %d robust species\n", length(matched), length(robust)))

# Per-species 2x2 cell means + 2-way ANOVA-style effect decomposition
res <- data.frame(species = matched,
                  mean_CS_ABx = NA_real_, mean_CS_noABx = NA_real_,
                  mean_VAG_ABx = NA_real_, mean_VAG_noABx = NA_real_,
                  cs_effect = NA_real_, abx_effect = NA_real_,
                  interaction = NA_real_,
                  cs_p = NA_real_, abx_p = NA_real_,
                  stringsAsFactors = FALSE)

for (i in seq_along(matched)) {
  sp <- matched[i]
  y <- clr_mat[sp, ]
  df <- data.frame(y = y, cs = mk$born_method == "c_section",
                   abx = mk$antibiotics_current_use == "yes")
  res$mean_CS_ABx[i]    <- mean(y[df$cs & df$abx])
  res$mean_CS_noABx[i]  <- mean(y[df$cs & !df$abx])
  res$mean_VAG_ABx[i]   <- mean(y[!df$cs & df$abx])
  res$mean_VAG_noABx[i] <- mean(y[!df$cs & !df$abx])
  # Marginal effects (averaged over the other factor)
  res$cs_effect[i]  <- mean(c(res$mean_CS_ABx[i], res$mean_CS_noABx[i])) -
                       mean(c(res$mean_VAG_ABx[i], res$mean_VAG_noABx[i]))
  res$abx_effect[i] <- mean(c(res$mean_CS_ABx[i], res$mean_VAG_ABx[i])) -
                       mean(c(res$mean_CS_noABx[i], res$mean_VAG_noABx[i]))
  res$interaction[i] <- (res$mean_CS_ABx[i] - res$mean_CS_noABx[i]) -
                        (res$mean_VAG_ABx[i] - res$mean_VAG_noABx[i])
  # Two-way ANOVA p-values
  fit <- tryCatch(anova(lm(y ~ cs + abx, data = df)),
                  error = function(e) NULL)
  if (!is.null(fit)) {
    res$cs_p[i]  <- fit$`Pr(>F)`[1]
    res$abx_p[i] <- fit$`Pr(>F)`[2]
  }
}
res$cs_q  <- p.adjust(res$cs_p,  method = "BH")
res$abx_q <- p.adjust(res$abx_p, method = "BH")

# Classify each species by which factor dominates
res$dominant <- ifelse(abs(res$cs_effect) > 2 * abs(res$abx_effect),
                       "CS-dominant",
                ifelse(abs(res$abx_effect) > 2 * abs(res$cs_effect),
                       "Antibiotic-dominant",
                       "Both contribute"))
cat("\nDominant-factor classification:\n"); print(table(res$dominant))

write.csv(res, file.path(proj, "results", "antibiotic_deconvolution.csv"),
          row.names = FALSE)

# ── Plot: effect-effect scatter ──
plt <- res
plt$label <- gsub("species:|genus:", "", plt$species)
g <- ggplot(plt, aes(cs_effect, abx_effect, colour = dominant)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.5, alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 30,
                           segment.size = 0.2) +
  scale_colour_manual(values = c("CS-dominant" = "#D6604D",
                                  "Antibiotic-dominant" = "#9970AB",
                                  "Both contribute" = "#5AAE61")) +
  labs(title = "Birth-mode vs antibiotic effect, per robust species (d0-7)",
       subtitle = "+x = CS-enriched; +y = antibiotic-enriched",
       x = "CS effect (CLR units)", y = "Antibiotic effect (CLR units)") +
  theme_bw(base_size = 11) + theme(legend.position = "right")
ggsave(file.path(proj, "figures", "antibiotic_vs_cs_effects.png"),
       g, width = 9, height = 7, dpi = 200)

# Headline numbers for the report
cat("\nHeadline:\n")
cat(sprintf("  CS-dominant species:        %d\n", sum(res$dominant == "CS-dominant")))
cat(sprintf("  Antibiotic-dominant:        %d\n", sum(res$dominant == "Antibiotic-dominant")))
cat(sprintf("  Both contribute:            %d\n", sum(res$dominant == "Both contribute")))
cat(sprintf("  Median CS effect magnitude: %.2f CLR\n", median(abs(res$cs_effect))))
cat(sprintf("  Median ABx effect magnitude: %.2f CLR\n", median(abs(res$abx_effect))))

cat("\n[18 done] antibiotic_deconvolution.csv + figure\n")
