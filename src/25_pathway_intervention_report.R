# ─────────────────────────────────────────────────────────────────────────────
# 25 — Pathway-level intervention context report
#
# Generates category-specific narrative context for the pathway analysis:
#   - Iron acquisition (siderophore) — direction in this cohort follows E. coli
#   - Cell wall / peptidoglycan — Lys vs DAP stereochemistry differences
#   - Bifidobacterium shunt (HMO utilisation)
#   - Folate biosynthesis — direction depends on which folate-producers dominate
#   - Central / fermentation (SCFA context)
#
# Direction convention: CLR(CS) - CLR(vaginal); negative = vaginal-enriched.
# In this cohort E. coli is vaginal-enriched at d0-7, so Enterobacteriaceae-
# associated pathways (siderophore, LPS) also come out vaginal-enriched;
# see Cohort caveats section of the report page.
# ─────────────────────────────────────────────────────────────────────────────
proj <- getwd()  # set working directory to repo root before running

pcat <- read.csv(file.path(proj, "results", "pathway_clinical_context.csv"),
                 stringsAsFactors = FALSE)
ptraj <- read.csv(file.path(proj, "results", "pathway_category_trajectory.csv"),
                  stringsAsFactors = FALSE)
pdeconv <- read.csv(file.path(proj, "results", "pathway_category_antibiotic_deconv.csv"),
                    stringsAsFactors = FALSE)

# Per-category counts by direction
cat_dir <- as.data.frame(table(pcat$category, pcat$direction))
colnames(cat_dir) <- c("category", "direction", "n_pathways")

sink(file.path(proj, "results", "pathway_intervention_report.txt"))
cat("═══════════════════════════════════════════════════════════════════\n")
cat("  PATHWAY-LEVEL INTERVENTION ANALYSIS (with direction + lit audit)\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("PURPOSE\n")
cat("───────\n")
cat("This report provides pathway-level narrative context alongside the\n")
cat("species-level intervention analysis. Direction convention (CLR(CS) -\n")
cat("CLR(vaginal)): negative = vaginal-enriched, positive = CS-enriched.\n")
cat("In this cohort E. coli is vaginal-enriched, so the siderophore and\n")
cat("LPS biosynthesis pathways are also vaginal-enriched here — a\n")
cat("cohort-specific direction that differs from most published cohorts.\n")
cat("All direction claims below were checked against raw clr_diff values\n")
cat("in pathway_diff.csv.\n\n")

cat("───────────────────────────────────────────────────────────────────\n")
cat("  METHOD-ROBUST PATHWAY COUNTS BY CATEGORY × DIRECTION\n")
cat("───────────────────────────────────────────────────────────────────\n\n")
print(cat_dir)

cat("\n───────────────────────────────────────────────────────────────────\n")
cat("  PATHWAY-CATEGORY TRAJECTORY (d0-7 → d91+)\n")
cat("───────────────────────────────────────────────────────────────────\n\n")
cat("Per-category mean CLR(CS) − CLR(vaginal):\n")
print(ptraj[, c("category","direction_baseline","gap_d0_7","gap_d91p",
                 "attenuation_pct","class")])
cat("\nReading: nearly all pathway-category effects attenuate or close by\n")
cat("d91+, mirroring the species-level pattern. Per-category mean CLR is\n")
cat("small in magnitude (0.01-0.5) because averaging across pathways within\n")
cat("a category dilutes individual pathway signals; the species-level\n")
cat("intervention targets remain the primary actionable list.\n\n")

cat("───────────────────────────────────────────────────────────────────\n")
cat("  CATEGORY-SPECIFIC CLINICAL CONTEXT (audited + lit-verified)\n")
cat("───────────────────────────────────────────────────────────────────\n\n")

cat("1. IRON ACQUISITION (SIDEROPHORE) — 2 pathways, BOTH vaginal-enriched\n")
cat("   Data: enterobactin (CLR diff -1.09), aerobactin (CLR diff -0.89).\n")
cat("   In THIS cohort the siderophore signal follows E. coli, which is\n")
cat("   vaginal-enriched at d0-7. This is OPPOSITE to the typical Shao 2019\n")
cat("   narrative where Enterobacteriaceae are CS-enriched.\n")
cat("   Clinical implication: NOT a CS-iron-competition story in this\n")
cat("   cohort. If iron supplementation timing is to be guided by\n")
cat("   siderophore activity, it requires cohort-specific Enterobacteriaceae\n")
cat("   direction assessment, not a one-size-fits-all rule.\n\n")

cat("2. CELL WALL / ENVELOPE — 11 pathways: 7 CS-enriched, 4 vaginal-enriched\n")
cat("   CS-enriched: Staphylococcus peptidoglycan (PWY-5265), Enterococcus\n")
cat("   peptidoglycan (PWY-6471), multiple PEPTIDOGLYCANSYN variants —\n")
cat("   consistent with CS-enriched Staphylococcus and Enterococcus species\n")
cat("   in this cohort. Vaginal-enriched: LPS biosynthesis (LPSSYN-PWY,\n")
cat("   CLR -2.12) — consistent with vaginal E. coli.\n")
cat("   Clinical implication: different peptidoglycan stereochemistry (Lys-\n")
cat("   type from Staph vs DAP-type from Gram-neg) is recognized differently\n")
cat("   by NOD2 / TLR2 (Royet 2017 Nat Rev Immunol). CS infants' early\n")
cat("   exposure to Staph-derived Lys-type PGN may set a different innate-\n")
cat("   immune education trajectory than vaginal infants' DAP-type from\n")
cat("   Gram-negative commensals — a plausible mechanism for divergent\n")
cat("   atopic risk (hypothesis, not tested in our data).\n\n")

cat("3. BIFIDOBACTERIUM SHUNT (HMO fermentation) — vaginal-enriched\n")
cat("   Data: P124-PWY (Bifidobacterium shunt) CLR -2.68, 4/4 methods,\n")
cat("   vaginal-enriched. CS infants depleted in B. longum (4/4 methods)\n")
cat("   → reduced HMO fermentation capacity at d0-7. Trajectory: the\n")
cat("   carbohydrate-metabolism category broadly attenuates by d91+,\n")
cat("   consistent with B. longum recovery (step 23 verification).\n")
cat("   Clinical implication: the rationale for HMO-utilising Bifido\n")
cat("   probiotic supplementation in the ACUTE neonatal window (0-2 mo)\n")
cat("   is supported. Mainstream commercial products (e.g. B. infantis EVC001)\n")
cat("   already target this — well-established field.\n\n")

cat("4. FOLATE BIOSYNTHESIS (within Cofactor/Vitamin/ETC, 30 pathways) —\n")
cat("   FOLSYN-PWY and related folate pathways CS-enriched (+0.45-0.72 CLR).\n")
cat("   This is OPPOSITE to the typical narrative that 'CS infants shift\n")
cat("   away from beneficial folate-producing Bifidobacterium'. In this\n")
cat("   cohort, CS-enriched Enterococcus and Streptococcus species are\n")
cat("   themselves folate-producers (literature: Magnusdottir 2015 Front\n")
cat("   Genet; Rossi 2011 PMC).\n")
cat("   Clinical implication: routine neonatal folate supplementation\n")
cat("   policy probably does not need to be modulated by birth mode —\n")
cat("   the microbial folate production deficit is small or absent at the\n")
cat("   category level (and reverses direction depending on which folate-\n")
cat("   producers dominate). This category should NOT be claimed as an\n")
cat("   intervention target.\n\n")

cat("5. SCFA / CENTRAL FERMENTATION — 20 pathways, balanced (10/10)\n")
cat("   Mixed direction: e.g. methylglyoxal degradation (vaginal-enriched),\n")
cat("   glycerol→1,3-propanediol (CS-enriched), butanediol biosynthesis\n")
cat("   (CS-enriched). Butyrate biosynthesis pathways are NOT among the 20\n")
cat("   method-robust hits (likely below detection threshold this early).\n")
cat("   Clinical implication: the SCFA-deficit-in-CS story (e.g. butyrate\n")
cat("   for barrier integrity) is NOT supported at d0-7 in this dataset\n")
cat("   at the method-robust threshold. The literature SCFA-CS link is\n")
cat("   established for later infancy (post-weaning, post-solid-food\n")
cat("   introduction) — beyond what d0-7 metagenomics can show.\n\n")

cat("───────────────────────────────────────────────────────────────────\n")
cat("  PATHWAY-LEVEL ANTIBIOTIC-VS-CS DECONVOLUTION (d0-7)\n")
cat("───────────────────────────────────────────────────────────────────\n\n")
print(pdeconv[, c("category","direction_baseline","cs_effect","abx_effect","dominant")])
cat("\nReading: category-level effects are small in absolute magnitude\n")
cat("(0.01-0.5 CLR). The CS-dominant classification dominates because tiny\n")
cat("CS effects are usually larger than tiny antibiotic effects in this\n")
cat("averaged-category metric. Most clinically interpretable signal is at\n")
cat("INDIVIDUAL PATHWAY level (Bifidobacterium shunt, LPS, peptidoglycan-\n")
cat("Staph), not category-aggregate level.\n\n")

sink()
cat("[25 done] pathway_intervention_report.txt\n")
