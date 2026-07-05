# ─────────────────────────────────────────────────────────────────────────────
# 10 — Classifier robustness: permutation test, calibration, importance stability
# Extends Step 04 (subject-grouped RF) with 3 additional evaluation layers:
#   1. Permutation AUC distribution (is 0.853 better than chance?)
#   2. Calibration curve (are predicted probabilities reliable?)
#   3. Feature importance stability across CV folds
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(randomForest); library(ggplot2) })
set.seed(42)
proj <- getwd()  # set working directory to repo root before running
d <- readRDS(file.path(proj, "results", "01_processed.rds"))
ab <- d$abundance; meta <- d$meta

clr <- function(mat){ m <- mat + min(mat[mat>0])/2; lg <- log(m); sweep(lg,2,colMeans(lg),"-") }
X <- t(clr(ab)); y <- meta$born_method; grp <- meta$subject_id
ytrue <- as.integer(y == "c_section")

auc <- function(score, label){
  pos <- score[label==1]; neg <- score[label==0]
  if(!length(pos)||!length(neg)) return(NA)
  r <- rank(c(pos,neg)); (sum(r[seq_along(pos)]) - length(pos)*(length(pos)+1)/2)/(length(pos)*length(neg))
}

# Subject-grouped fold assignment (same as Step 04)
subs <- unique(grp); set.seed(42); subs <- sample(subs)
fold_of_sub <- setNames(rep(1:5, length.out = length(subs)), subs)
folds <- fold_of_sub[as.character(grp)]

# ─── Real model: out-of-fold predictions + per-fold importance ───
oof <- rep(NA_real_, nrow(X))
fold_auc <- numeric(5)
imp_list <- list()

for (k in 1:5) {
  tr <- which(folds != k); te <- which(folds == k)
  rf <- randomForest(x = X[tr,,drop=FALSE], y = y[tr], ntree = 800, importance = TRUE)
  pr <- predict(rf, X[te,,drop=FALSE], type = "prob")[, "c_section"]
  oof[te] <- pr
  fold_auc[k] <- auc(pr, ytrue[te])
  imp_list[[k]] <- importance(rf)[, "MeanDecreaseGini"]
}
real_auc <- auc(oof, ytrue)
cat(sprintf("Real model: mean fold AUC = %.3f, pooled AUC = %.3f\n", mean(fold_auc), real_auc))

# ─── 1. Permutation test: shuffle labels 200 times ───
cat("Running 200 permutations...\n")
n_perm <- 200
perm_auc <- numeric(n_perm)
for (i in seq_len(n_perm)) {
  set.seed(i + 1000)
  y_shuf <- sample(y)
  oof_p <- rep(NA_real_, nrow(X))
  for (k in 1:5) {
    tr <- which(folds != k); te <- which(folds == k)
    rf_p <- randomForest(x = X[tr,,drop=FALSE], y = y_shuf[tr], ntree = 400)
    oof_p[te] <- predict(rf_p, X[te,,drop=FALSE], type = "prob")[, "c_section"]
  }
  perm_auc[i] <- auc(oof_p, as.integer(y_shuf == "c_section"))
  if (i %% 50 == 0) cat(sprintf("  permutation %d/%d done\n", i, n_perm))
}
perm_p <- mean(perm_auc >= real_auc)
cat(sprintf("Permutation p-value: %.4f (real AUC %.3f vs perm mean %.3f)\n",
            perm_p, real_auc, mean(perm_auc)))

p_perm <- ggplot(data.frame(auc = perm_auc), aes(auc)) +
  geom_histogram(bins = 30, fill = "grey70", colour = "white") +
  geom_vline(xintercept = real_auc, colour = "#D7301F", linewidth = 1.2) +
  annotate("text", x = real_auc - 0.01, y = Inf, vjust = 2, hjust = 1,
           label = sprintf("Real AUC = %.3f\np(perm) = %.4f", real_auc, perm_p),
           colour = "#D7301F", size = 3.8) +
  labs(title = "Permutation test: birth-mode classifier",
       subtitle = sprintf("%d label-shuffled runs  |  subject-grouped 5-fold CV", n_perm),
       x = "Permuted AUC", y = "Count") +
  theme_bw(base_size = 12)
ggsave(file.path(proj, "figures", "classifier_permutation.png"), p_perm,
       width = 6, height = 4, dpi = 200)

# ─── 2. Calibration curve ───
# Bin predicted probabilities and compare to observed frequency
bins <- seq(0, 1, by = 0.1)
cal <- data.frame(
  bin_mid = (head(bins, -1) + tail(bins, -1)) / 2,
  pred_mean = tapply(oof, cut(oof, bins, include.lowest = TRUE), mean),
  obs_freq  = tapply(ytrue, cut(oof, bins, include.lowest = TRUE), mean),
  n         = tapply(ytrue, cut(oof, bins, include.lowest = TRUE), length)
)
cal <- cal[!is.na(cal$obs_freq), ]
brier <- mean((oof - ytrue)^2)

p_cal <- ggplot(cal, aes(pred_mean, obs_freq)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_point(aes(size = n), colour = "#1B7837") +
  geom_line(colour = "#1B7837") +
  scale_size_continuous(range = c(2, 8), name = "n samples") +
  labs(title = "Calibration curve: birth-mode classifier",
       subtitle = sprintf("Brier score = %.4f  |  perfect calibration = diagonal", brier),
       x = "Mean predicted P(C-section)", y = "Observed fraction C-section") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_bw(base_size = 12)
ggsave(file.path(proj, "figures", "classifier_calibration.png"), p_cal,
       width = 5.5, height = 5, dpi = 200)

# ─── 3. Feature importance stability across folds ───
imp_mat <- do.call(cbind, imp_list)
colnames(imp_mat) <- paste0("fold_", 1:5)
imp_df <- data.frame(
  taxon    = rownames(imp_mat),
  mean_imp = rowMeans(imp_mat),
  sd_imp   = apply(imp_mat, 1, sd),
  cv_imp   = apply(imp_mat, 1, sd) / rowMeans(imp_mat)
)
imp_df <- imp_df[order(-imp_df$mean_imp), ]

# Top 20 features
top20 <- head(imp_df, 20)
top20$short <- gsub("s__", "", top20$taxon)

p_imp <- ggplot(top20, aes(reorder(short, mean_imp), mean_imp)) +
  geom_col(fill = "#2C7FB8", alpha = 0.85) +
  geom_errorbar(aes(ymin = pmax(0, mean_imp - sd_imp), ymax = mean_imp + sd_imp),
                width = 0.3, colour = "grey30") +
  coord_flip() +
  labs(title = "Feature importance stability (top 20 taxa)",
       subtitle = "Mean ± SD Gini importance across 5 CV folds",
       x = NULL, y = "Mean Decrease Gini") +
  theme_bw(base_size = 10)
ggsave(file.path(proj, "figures", "classifier_importance_stability.png"), p_imp,
       width = 7.5, height = 5.5, dpi = 200)

# ─── 4. Sensitivity + Specificity + F1 at optimal threshold ───
thresholds <- sort(unique(oof))
metrics <- data.frame(t(sapply(thresholds, function(th) {
  pred <- as.integer(oof >= th)
  tp <- sum(pred == 1 & ytrue == 1); fp <- sum(pred == 1 & ytrue == 0)
  fn <- sum(pred == 0 & ytrue == 1); tn <- sum(pred == 0 & ytrue == 0)
  sens <- tp / (tp + fn); spec <- tn / (tn + fp)
  prec <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  f1   <- ifelse(prec + sens > 0, 2 * prec * sens / (prec + sens), 0)
  c(threshold = th, sensitivity = sens, specificity = spec,
    precision = prec, f1 = f1, youden = sens + spec - 1)
})))
best <- metrics[which.max(metrics$youden), ]
cat(sprintf("\nOptimal threshold (Youden): %.3f\n", best$threshold))
cat(sprintf("  Sensitivity: %.3f  Specificity: %.3f  F1: %.3f  Precision: %.3f\n",
            best$sensitivity, best$specificity, best$f1, best$precision))

# ─── Save everything ───
sink(file.path(proj, "results", "classifier_robustness.txt"))
cat("=== Classifier robustness evaluation ===\n\n")
cat("Model: Random Forest (800 trees), subject-grouped 5-fold CV\n")
cat("Features: CLR-transformed species profiles\n\n")
cat(sprintf("Pooled out-of-fold AUC: %.3f\n", real_auc))
cat(sprintf("Per-fold AUC: %s\n", paste(sprintf("%.3f", fold_auc), collapse = ", ")))
cat(sprintf("\n--- Permutation test (%d runs) ---\n", n_perm))
cat(sprintf("Permutation AUC: mean = %.3f, SD = %.3f, range = [%.3f, %.3f]\n",
            mean(perm_auc), sd(perm_auc), min(perm_auc), max(perm_auc)))
cat(sprintf("Empirical p-value: %.4f\n", perm_p))
cat(sprintf("\n--- Calibration ---\nBrier score: %.4f (0 = perfect)\n", brier))
cat(sprintf("\n--- Optimal operating point (Youden) ---\n"))
cat(sprintf("Threshold: %.3f\n", best$threshold))
cat(sprintf("Sensitivity: %.3f  Specificity: %.3f\n", best$sensitivity, best$specificity))
cat(sprintf("Precision: %.3f  F1: %.3f\n", best$precision, best$f1))
cat(sprintf("\n--- Top 10 most stable features (by mean Gini) ---\n"))
print(head(imp_df[, c("taxon","mean_imp","sd_imp","cv_imp")], 10), row.names = FALSE)
sink()

write.csv(imp_df, file.path(proj, "results", "classifier_importance_stability.csv"),
          row.names = FALSE)

cat("[10 done] results/classifier_robustness.txt + 3 figures\n")
