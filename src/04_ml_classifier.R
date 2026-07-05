# ─────────────────────────────────────────────────────────────────────────────
# 04 — ML classifier: predict birth mode from CLR species profile
# LEAKAGE-SAFE: subject-grouped 5-fold CV (no subject in both train & test).
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({ library(randomForest); library(ggplot2) })
set.seed(42)
proj <- getwd()  # set working directory to repo root before running
d <- readRDS(file.path(proj, "results", "01_processed.rds"))
ab <- d$abundance; meta <- d$meta

# CLR features (samples x features)
clr <- function(mat){ m <- mat + min(mat[mat>0])/2; lg <- log(m); sweep(lg,2,colMeans(lg),"-") }
X <- t(clr(ab)); y <- meta$born_method; grp <- meta$subject_id

# rank-based AUC (no extra package)
auc <- function(score, label){            # label: 1 = positive (c_section)
  pos <- score[label==1]; neg <- score[label==0]
  if(!length(pos)||!length(neg)) return(NA)
  r <- rank(c(pos,neg)); (sum(r[seq_along(pos)]) - length(pos)*(length(pos)+1)/2)/(length(pos)*length(neg))
}

# subject-grouped 5-fold assignment
subs <- unique(grp); set.seed(42); subs <- sample(subs)
fold_of_sub <- setNames(rep(1:5, length.out = length(subs)), subs)
folds <- fold_of_sub[as.character(grp)]

oof <- rep(NA_real_, nrow(X)); ytrue <- as.integer(y == "c_section")
fold_auc <- numeric(5)
for (k in 1:5) {
  tr <- which(folds != k); te <- which(folds == k)
  rf <- randomForest(x = X[tr,,drop=FALSE], y = y[tr], ntree = 800)
  pr <- predict(rf, X[te,,drop=FALSE], type = "prob")[, "c_section"]
  oof[te] <- pr
  fold_auc[k] <- auc(pr, ytrue[te])
  cat(sprintf("fold %d: n_test=%d  AUC=%.3f\n", k, length(te), fold_auc[k]))
}
overall_auc <- auc(oof, ytrue)

# importance from a full-data model (for reporting top taxa)
rf_full <- randomForest(x = X, y = y, ntree = 800, importance = TRUE)
imp <- importance(rf_full)[, "MeanDecreaseGini"]
imp <- sort(imp, decreasing = TRUE)
topimp <- data.frame(taxon = sub(".*s__","",names(head(imp,15))), gini = head(imp,15))

sink(file.path(proj, "results", "ml_performance.txt"))
cat("Random forest — predict birth mode from CLR species profile\n")
cat("Validation: subject-grouped 5-fold CV (leakage-safe)\n\n")
cat("per-fold AUC:", paste(sprintf("%.3f", fold_auc), collapse=", "), "\n")
cat(sprintf("mean fold AUC: %.3f (SD %.3f)\n", mean(fold_auc), sd(fold_auc)))
cat(sprintf("pooled out-of-fold AUC: %.3f\n", overall_auc))
cat("\nTop 15 taxa by Gini importance:\n"); print(topimp, row.names = FALSE)
sink()

# ROC from pooled out-of-fold predictions
thr <- sort(unique(oof), decreasing = TRUE)
roc <- data.frame(t(sapply(thr, function(t){
  pred <- as.integer(oof >= t)
  tpr <- sum(pred==1 & ytrue==1)/sum(ytrue==1)
  fpr <- sum(pred==1 & ytrue==0)/sum(ytrue==0); c(fpr=fpr, tpr=tpr)
})))
roc <- rbind(c(0,0), roc, c(1,1))
p <- ggplot(roc, aes(fpr, tpr)) +
  geom_abline(slope=1, intercept=0, linetype=2, colour="grey60") +
  geom_path(colour="#1B7837", linewidth=1) +
  annotate("text", x=.62, y=.18, label=sprintf("pooled AUC = %.3f", overall_auc), size=4.2) +
  labs(title="Birth-mode classifier (subject-grouped CV)",
       subtitle="Random forest on CLR species profile  |  Shao 2019",
       x="False positive rate", y="True positive rate") +
  coord_equal() + theme_bw(base_size=12)
ggsave(file.path(proj,"figures","roc.png"), p, width=4.8, height=4.8, dpi=200)
cat(sprintf("[04 done] mean fold AUC %.3f, pooled AUC %.3f -> results/ml_performance.txt + figures/roc.png\n",
            mean(fold_auc), overall_auc))
