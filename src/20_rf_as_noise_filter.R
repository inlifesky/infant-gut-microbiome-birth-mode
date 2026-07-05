# ─────────────────────────────────────────────────────────────────────────────
# 20 — RF as a third noise-filter line (NOT as a classifier)
#
# The original RF analysis (step 04, 10) reported AUC, calibration, and
# permutation p-values — classifier metrics. But birth mode is ALREADY KNOWN
# from the medical record; the classifier has no diagnostic utility.
#
# What RF IS useful for: feature importance is a DATA-DRIVEN signal-strength
# measure that's independent of any DA model's distributional assumptions.
# A taxon that ranks high in RF importance AND survives 4-method DA concordance
# is supported by two independent analytical paradigms — stronger evidence of
# real biology than either alone.
#
# This script reframes the RF outputs as the third triangulation line:
#   - Line 1: 4-method DA concordance (steps 06-09)
#   - Line 2: 4-method pathway concordance (steps 12-15)
#   - Line 3: RF importance triangulation (THIS, replacing classifier framing)
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
proj <- getwd()  # set working directory to repo root before running

# Inputs
imp  <- read.csv(file.path(proj, "results", "classifier_importance_stability.csv"),
                 stringsAsFactors = FALSE)
conc <- read.csv(file.path(proj, "results", "DA_concordance.csv"),
                 stringsAsFactors = FALSE)
r03  <- read.csv(file.path(proj, "results", "diff_abundance.csv"),
                 stringsAsFactors = FALSE)
clean <- function(x) { x <- sub("^species[.:]", "", x); x <- sub("^genus[.:]", "", x)
                       tolower(gsub("[. ]+", "_", x)) }
imp$tc  <- clean(imp$taxon)
conc$tc <- conc$taxon
r03$tc  <- clean(r03$taxon)

# Build triangulation table
df <- data.frame(taxon = unique(c(imp$tc, conc$tc)), stringsAsFactors = FALSE)
df$rf_imp     <- imp$mean_imp[match(df$taxon, imp$tc)]
df$rf_rank    <- rank(-df$rf_imp, na.last = "keep", ties.method = "min")
df$da_methods <- conc$n_methods[match(df$taxon, conc$tc)]
df$da_methods[is.na(df$da_methods)] <- 0
df$clr_diff   <- r03$clr_diff[match(df$taxon, r03$tc)]
df$direction  <- ifelse(df$clr_diff > 0, "CS-enriched", "Vaginal-enriched")

# Triangulation classification
df$triangulated <- ifelse(df$da_methods >= 3 & !is.na(df$rf_rank) & df$rf_rank <= 15,
                          "Both lines support",
                   ifelse(df$da_methods >= 3, "DA only",
                   ifelse(!is.na(df$rf_rank) & df$rf_rank <= 15, "RF only",
                                                                 "Neither")))
tt <- table(df$triangulated)
cat("Triangulation table:\n"); print(tt)
cat(sprintf("\nAgreement rate (DA-robust taxa that are also RF-top-15): %.0f%%\n",
            100 * tt["Both lines support"] /
              (tt["Both lines support"] + tt["DA only"])))

df <- df[order(-df$da_methods, df$rf_rank), ]
write.csv(df, file.path(proj, "results", "triangulation_table.csv"), row.names = FALSE)

# Headline plot: DA-methods (x) vs RF importance rank (y), shaded by direction
plt <- df[!is.na(df$rf_imp), ]
plt$rank_top15 <- pmin(plt$rf_rank, 16)   # collapse beyond 15
plt$label <- gsub("_", " ", plt$taxon)
plt$label <- paste0(toupper(substr(plt$label,1,1)), substr(plt$label,2,nchar(plt$label)))

g <- ggplot(plt, aes(da_methods, rank_top15, colour = direction)) +
  geom_hline(yintercept = 15.5, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 2.5,  linetype = "dashed", colour = "grey50") +
  geom_jitter(aes(size = rf_imp), width = 0.15, height = 0.25, alpha = 0.75) +
  ggrepel::geom_text_repel(
    data = plt[plt$da_methods >= 3 & plt$rank_top15 <= 15, ],
    aes(label = label), size = 2.5, max.overlaps = 50, segment.size = 0.2) +
  scale_y_reverse(breaks = c(1, 5, 10, 15, 16),
                   labels = c("1","5","10","15",">15")) +
  scale_x_continuous(breaks = 0:4) +
  scale_colour_manual(values = c("CS-enriched"="#D6604D","Vaginal-enriched"="#2166AC")) +
  scale_size_continuous(range = c(1, 5), name = "RF importance") +
  labs(title = "Triangulation: DA methods agreement × RF importance rank",
       subtitle = "Top-right quadrant (≥3 DA methods + RF top-15) = strongest evidence",
       x = "Number of DA methods (of 4) calling FDR<0.05",
       y = "RF importance rank") +
  theme_bw(base_size = 11)
ggsave(file.path(proj, "figures", "triangulation_da_rf.png"),
       g, width = 9, height = 7, dpi = 200)

cat("\n[20 done] triangulation_table.csv + figure\n")
