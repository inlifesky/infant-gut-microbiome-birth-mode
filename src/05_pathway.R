# ─────────────────────────────────────────────────────────────────────────────
# 05 — Functional angle: metabolic pathway differences by birth mode
# Uses HUMAnN pathway abundances (MetaCyc) — connects taxonomic shifts to
# metabolic/functional potential (the project's "mitochondrial/metabolic" theme).
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(curatedMetagenomicData); library(SummarizedExperiment); library(ggplot2)
})
set.seed(42)
proj <- getwd()  # set working directory to repo root before running

pw <- curatedMetagenomicData("2021-03-31.ShaoY_2019.pathway_abundance",
                             dryrun = FALSE, rownames = "short")[[1]]
# use the SAME one-per-infant samples selected in 01 (no pseudoreplication / no mothers)
ids <- readRDS(file.path(proj, "results", "01_processed.rds"))$sample_ids
pw <- pw[, intersect(colnames(pw), ids)]
cat("pathway samples matched to step-01 selection:", ncol(pw), "\n")
grp <- factor(pw$born_method, levels = c("vaginal","c_section"))

M <- assay(pw)
# keep stratified-free, named pathways (drop UNMAPPED/UNINTEGRATED + species-stratified rows)
rn <- rownames(M)
M <- M[!grepl("\\|", rn) & !grepl("UNMAPPED|UNINTEGRATED", rn), , drop = FALSE]
if (max(colSums(M), na.rm=TRUE) > 1.5) M <- sweep(M, 2, colSums(M), "/")  # to proportions
prev <- rowMeans(M > 0); M <- M[prev >= 0.10, , drop = FALSE]
cat("pathways tested:", nrow(M), " samples:", ncol(M), "\n")

clr <- function(mat){ m <- mat + min(mat[mat>0])/2; lg <- log(m); sweep(lg,2,colMeans(lg),"-") }
cl <- clr(M)
res <- data.frame(pathway = rownames(cl),
                  clr_diff = rowMeans(cl[,grp=="c_section",drop=FALSE]) -
                             rowMeans(cl[,grp=="vaginal",drop=FALSE]),
                  p = apply(cl, 1, function(x) wilcox.test(x ~ grp)$p.value))
res$fdr <- p.adjust(res$p, "BH"); res <- res[order(res$fdr), ]
write.csv(res, file.path(proj,"results","pathway_diff.csv"), row.names = FALSE)
cat("significant pathways (FDR<0.05):", sum(res$fdr<0.05), "\n")

top <- head(res[res$fdr < 0.05, ], 15)
top$lab <- ifelse(nchar(top$pathway) > 48, paste0(substr(top$pathway,1,46),".."), top$pathway)
top$dir <- ifelse(top$clr_diff > 0, "C-section", "Vaginal")
p <- ggplot(top, aes(reorder(lab, clr_diff), clr_diff, fill = dir)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c("Vaginal"="#2C7FB8","C-section"="#D95F0E")) +
  labs(title="Top differential metabolic pathways by birth mode",
       subtitle="HUMAnN MetaCyc pathways, CLR + Wilcoxon (FDR<0.05)  |  Shao 2019",
       x=NULL, y="CLR difference (C-section - vaginal)", fill="Enriched in") +
  theme_bw(base_size = 10)
ggsave(file.path(proj,"figures","pathways.png"), p, width = 8, height = 5.5, dpi = 200)
cat("[05 done] results/pathway_diff.csv + figures/pathways.png\n")
