# ─────────────────────────────────────────────────────────────────────────────
# 19 — Elective vs emergency CS stratification
#
# Clinical question: do elective and emergency CS disrupt the microbiome
# equally? They differ in:
#   - Timing of antibiotic prophylaxis relative to delivery
#   - Stress hormones at delivery (emergency CS = labour started)
#   - Maternal vaginal microbe exposure (some emergency CS happen after labour
#     onset, giving some vaginal microbe contact)
#
# If elective CS disrupts MORE than emergency CS (or vice versa), intervention
# recommendations should differentiate.
#
# Design: within d0-7 window, compare CS_elective vs CS_emergency vs vaginal
# for each method-robust species. Report CLR means + pairwise tests.
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running

d <- readRDS(file.path(proj, "results", "01b_processed_longitudinal.rds"))
meta <- d$meta; ab <- d$abundance
meta$age_bin <- cut(meta$infant_age, c(-1,7,30,90,Inf),
                    labels = c("d0_7","d8_30","d31_90","d91+"))

# Three-level birth-mode variable
meta$mode3 <- ifelse(meta$born_method == "vaginal", "Vaginal",
              ifelse(meta$c_section_type == "Elective_CS", "CS_elective",
              ifelse(meta$c_section_type == "Emergency_CS", "CS_emergency", NA)))

keep <- which(meta$age_bin == "d0_7" & !is.na(meta$mode3))
mk <- meta[keep, ]; abk <- ab[, keep]
cat("Cell sizes (d0-7):\n"); print(table(mk$mode3))

clr <- function(x) { x <- x + 1e-6; log(x) - mean(log(x)) }
clr_mat <- apply(abk, 2, clr)

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

# Per-species: mean CLR per group + pairwise p
res <- data.frame(species = matched,
                  mean_vag = NA_real_, mean_elective = NA_real_, mean_emergency = NA_real_,
                  elective_vs_vag = NA_real_, emergency_vs_vag = NA_real_,
                  elective_vs_emergency = NA_real_,
                  p_elec_vag = NA_real_, p_emerg_vag = NA_real_, p_elec_emerg = NA_real_,
                  stringsAsFactors = FALSE)
for (i in seq_along(matched)) {
  sp <- matched[i]; y <- clr_mat[sp, ]
  v <- y[mk$mode3 == "Vaginal"]
  e <- y[mk$mode3 == "CS_elective"]
  m <- y[mk$mode3 == "CS_emergency"]
  res$mean_vag[i] <- mean(v); res$mean_elective[i] <- mean(e); res$mean_emergency[i] <- mean(m)
  res$elective_vs_vag[i]      <- mean(e) - mean(v)
  res$emergency_vs_vag[i]     <- mean(m) - mean(v)
  res$elective_vs_emergency[i] <- mean(e) - mean(m)
  res$p_elec_vag[i]   <- tryCatch(wilcox.test(e, v)$p.value, error=function(e) NA)
  res$p_emerg_vag[i]  <- tryCatch(wilcox.test(m, v)$p.value, error=function(e) NA)
  res$p_elec_emerg[i] <- tryCatch(wilcox.test(e, m)$p.value, error=function(e) NA)
}
res$q_elec_vag   <- p.adjust(res$p_elec_vag,   "BH")
res$q_emerg_vag  <- p.adjust(res$p_emerg_vag,  "BH")
res$q_elec_emerg <- p.adjust(res$p_elec_emerg, "BH")

# Classify
res$pattern <- with(res, ifelse(
  sign(elective_vs_vag) != sign(emergency_vs_vag) |
    (abs(elective_vs_vag) > 1 & abs(emergency_vs_vag) < 0.3 * abs(elective_vs_vag)) |
    (abs(emergency_vs_vag) > 1 & abs(elective_vs_vag) < 0.3 * abs(emergency_vs_vag)),
  "Differential",
  ifelse(abs(elective_vs_emergency) > 0.5 * pmax(abs(elective_vs_vag), abs(emergency_vs_vag)),
         "Quantitative_difference", "Similar")))
cat("\nElective vs emergency pattern:\n"); print(table(res$pattern))

write.csv(res, file.path(proj, "results", "cs_type_stratification.csv"), row.names = FALSE)

# Plot: per-species elective vs emergency CS effects
plt <- res
plt$label <- gsub("species:|genus:", "", plt$species)
g <- ggplot(plt, aes(elective_vs_vag, emergency_vs_vag, colour = pattern)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60") +
  geom_point(size = 2.5, alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 30,
                           segment.size = 0.2) +
  scale_colour_manual(values = c("Similar" = "grey55",
                                  "Quantitative_difference" = "#E08214",
                                  "Differential" = "#762A83")) +
  labs(title = "Elective vs emergency CS effect (d0-7, vs vaginal baseline)",
       subtitle = "Dashed line = equal disruption; points off-line differentiate CS subtypes",
       x = "Elective CS vs vaginal (CLR diff)",
       y = "Emergency CS vs vaginal (CLR diff)") +
  theme_bw(base_size = 11)
ggsave(file.path(proj, "figures", "cs_type_stratification.png"),
       g, width = 9, height = 7, dpi = 200)

cat("\n[19 done] cs_type_stratification.csv + figure\n")
