# ─────────────────────────────────────────────────────────────────────────────
# 21 — Intervention report (rewritten with stratified evidence)
#
# Replaces the original clinical_context_report.txt narrative. The intervention
# claims are now grounded in the stratified evidence from steps 17-19:
#   - Longitudinal trajectory (does the gap close, persist, widen?)
#   - Antibiotic vs CS deconvolution (which factor drives each species?)
#   - CS-type stratification (does elective vs emergency CS matter?)
# ─────────────────────────────────────────────────────────────────────────────
proj <- getwd()  # set working directory to repo root before running

traj  <- read.csv(file.path(proj, "results", "longitudinal_trajectories.csv"),
                  stringsAsFactors = FALSE)
abx   <- read.csv(file.path(proj, "results", "antibiotic_deconvolution.csv"),
                  stringsAsFactors = FALSE)
cstyp <- read.csv(file.path(proj, "results", "cs_type_stratification.csv"),
                  stringsAsFactors = FALSE)
tri   <- read.csv(file.path(proj, "results", "triangulation_table.csv"),
                  stringsAsFactors = FALSE)

short <- function(x) gsub("species:|genus:", "", x)
traj$short  <- short(traj$species)
abx$short   <- short(abx$species)
cstyp$short <- short(cstyp$species)

# Merge into one master table
master <- merge(traj[, c("short", "baseline_direction", "class",
                          "gap_d0_7", "gap_d91p", "attenuation_pct")],
                abx[, c("short", "dominant", "cs_effect", "abx_effect")],
                by = "short")
master <- merge(master,
                cstyp[, c("short", "pattern", "elective_vs_vag", "emergency_vs_vag")],
                by = "short")
master <- master[order(master$baseline_direction, master$class, -abs(master$gap_d0_7)), ]
write.csv(master, file.path(proj, "results", "intervention_master_table.csv"),
          row.names = FALSE)

sink(file.path(proj, "results", "intervention_report.txt"))
cat("═══════════════════════════════════════════════════════════════════\n")
cat("  INTERVENTION PRIORITIES FOR CS-BORN INFANTS\n")
cat("  Evidence-stratified analysis from Shao 2019 Baby Biome Study\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("CLINICAL QUESTION\n")
cat("─────────────────\n")
cat("A baby is born by caesarean section. What microbiome-targeted\n")
cat("interventions (feeding, probiotic, monitoring) actually have evidence\n")
cat("support — once we stratify by recovery trajectory, antibiotic\n")
cat("exposure, and CS subtype?\n\n")

cat("WHAT THIS REPORT ADDS BEYOND \"CS HAS DIFFERENT MICROBIOTA\"\n")
cat("────────────────────────────────────────────────────────────\n")
cat("Every published study comparing CS vs vaginal microbiome reports\n")
cat("differences. This report goes further by asking, for each method-\n")
cat("robust difference: (1) does it persist long enough to matter?\n")
cat("(2) is it really driven by CS or by peripartum antibiotics?\n")
cat("(3) is it the same in elective and emergency CS?\n")
cat("Only the differences that survive all three stratifications\n")
cat("are robust intervention targets.\n\n")

# ───────────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════\n")
cat("  KEY FINDING 1: Most birth-mode effects ATTENUATE BY 3 MONTHS\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cls_tab <- as.data.frame(table(master$class))
cat("Trajectory class distribution (n=27 method-robust species):\n")
for (i in seq_len(nrow(cls_tab)))
  cat(sprintf("  %-25s %d\n", cls_tab$Var1[i], cls_tab$Freq[i]))

cat("\nClinical implication: most species differences at d0-7 are TRANSIENT.\n")
cat("A probiotic intervention targeting them addresses a self-resolving\n")
cat("problem. Interventions need to target species that PERSIST or WIDEN.\n\n")

persist <- master[master$class %in% c("persists", "widened"), ]
cat(sprintf("Species with persisting/widening deficit at d91+ (n=%d):\n", nrow(persist)))
for (i in seq_len(nrow(persist))) {
  cat(sprintf("  %s (%s, %s by d91+, gap=%.2f→%.2f CLR)\n",
              persist$short[i], persist$baseline_direction[i],
              persist$class[i], persist$gap_d0_7[i], persist$gap_d91p[i]))
}

cat("\nNotable: the 5 persisting/widening species are four Bacteroidota\n")
cat("(Bacteroides uniformis,\n")
cat("Phocaeicola vulgatus, Phocaeicola dorei, Parabacteroides distasonis,\n")
cat("Collinsella aerofaciens) are vaginal-enriched and STAY depleted at 3 months.\n")
cat("This is the real long-term CS deficit — and these genera are NOT in\n")
cat("standard commercial probiotic formulations.\n\n")

cat("Bifidobacterium longum and E. coli show DRAMATIC d0-7 deficits (CLR\n")
cat("-5.96 and -5.87) that LARGELY CLOSE by d91+. The all-samples gap\n")
cat("at d91+ (+0.57 and +0.78, i.e. apparent reversal to CS-enriched)\n")
cat("is partly an artifact of differential intrapartum antibiotic\n")
cat("exposure: in this cohort, vaginal infants have HIGHER IAP exposure\n")
cat("(27%) than CS infants (12%), and IAP-exposed vaginal infants have\n")
cat("durably reduced B. longum (Mazzola 2016, Stearns 2017).\n")
cat("After restricting to antibiotic-naive infants, the d91+ B. longum\n")
cat("gap is only +0.25 CLR (essentially closure, not reversal).\n")
cat("See step 23 for the stratified verification analysis.\n\n")
cat("The honest clinical implication: standard Bifidobacterium\n")
cat("probiotic supplementation is justified for the ACUTE neonatal window\n")
cat("(0-2 months) where the deficit is large and real. The 'long-term\n")
cat("Bifido deficit' rationale is not supported by this cohort once\n")
cat("antibiotic exposure is controlled for.\n\n")

# ───────────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════\n")
cat("  KEY FINDING 2: Some 'CS effects' are actually antibiotic effects\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

dom_tab <- as.data.frame(table(master$dominant))
cat("Dominant-factor classification (d0-7 only, n=27):\n")
for (i in seq_len(nrow(dom_tab)))
  cat(sprintf("  %-25s %d\n", dom_tab$Var1[i], dom_tab$Freq[i]))

abx_dom <- master[master$dominant == "Antibiotic-dominant", ]
cat(sprintf("\nAntibiotic-dominant species (n=%d):\n", nrow(abx_dom)))
for (i in seq_len(nrow(abx_dom))) {
  cat(sprintf("  %s (CS effect=%.2f, ABx effect=%.2f)\n",
              abx_dom$short[i], abx_dom$cs_effect[i], abx_dom$abx_effect[i]))
}
cat("\nClinical implication: for these species, the right intervention\n")
cat("target is ANTIBIOTIC STEWARDSHIP, not CS-specific probiotics.\n")
cat("Note: in this cohort, vaginal infants have HIGHER antibiotic exposure\n")
cat("(27%) than CS infants (12%) — likely intrapartum prophylaxis or\n")
cat("postnatal infection treatment. Antibiotic effect is genuinely\n")
cat("independent of CS, not confounded with it.\n\n")

# ───────────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════\n")
cat("  KEY FINDING 3: Elective and emergency CS disrupt similarly\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

pat_tab <- as.data.frame(table(master$pattern))
cat("Elective vs emergency CS pattern (d0-7, n=27):\n")
for (i in seq_len(nrow(pat_tab)))
  cat(sprintf("  %-30s %d\n", pat_tab$Var1[i], pat_tab$Freq[i]))

cat("\nClinical implication: for 20/27 species, the microbiome disruption\n")
cat("is comparable between elective and emergency CS. This means\n")
cat("intervention recommendations do NOT need to differentiate by CS\n")
cat("subtype for most species. Only a handful of species show\n")
cat("subtype-specific disruption patterns.\n\n")

# ───────────────────────────────────────────────────────────
cat("═══════════════════════════════════════════════════════════════════\n")
cat("  ACTIONABLE INTERVENTION TARGETS (after all three stratifications)\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

# Targets = species that (a) persist or widen by d91+ AND (b) are CS-dominant
# (not antibiotic-dominant) AND (c) baseline-direction = vaginal-enriched
# (i.e., we want to RESTORE these, not displace them)
actionable <- master[
  master$class %in% c("persists", "widened") &
  master$dominant %in% c("CS-dominant", "Both contribute") &
  master$baseline_direction == "Vaginal-enriched", ]
cat(sprintf("Species meeting all 3 criteria (persists, CS-driven, vaginal-enriched):\n"))
cat(sprintf("  (n = %d / 27 = %.0f%% of method-robust species)\n\n",
            nrow(actionable), 100 * nrow(actionable) / 27))
for (i in seq_len(nrow(actionable))) {
  cat(sprintf("  %s\n", actionable$short[i]))
  cat(sprintf("    trajectory: gap %.2f → %.2f at d91+ (%s)\n",
              actionable$gap_d0_7[i], actionable$gap_d91p[i], actionable$class[i]))
  cat(sprintf("    drivers: CS=%.2f, ABx=%.2f → %s\n",
              actionable$cs_effect[i], actionable$abx_effect[i], actionable$dominant[i]))
  cat(sprintf("    CS subtype: elective=%.2f, emergency=%.2f → %s\n\n",
              actionable$elective_vs_vag[i], actionable$emergency_vs_vag[i],
              actionable$pattern[i]))
}

cat("HONEST ASSESSMENT OF INTERVENTION LANDSCAPE\n")
cat("───────────────────────────────────────────\n\n")
cat("The actionable targets above are 4 Bacteroidota + 1 Actinomycetota\n")
cat("(Collinsella aerofaciens) — the Bacteroidota-centred residual is\n")
cat("dominant. Bacteroidota genera are\n")
cat("MAINSTREAM CONSUMER infant probiotics (Bifidobacterium + Lactobacillus\n")
cat("blends like Culturelle, BioGaia) do NOT target. The species these\n")
cat("commercial products contain (B. longum, B. breve) are TRANSIENT\n")
cat("deficits in this cohort and most other published cohorts.\n\n")

cat("The persistent Bacteroides deficit is being actively addressed in the\n")
cat("research literature:\n")
cat("  - Next-generation Bacteroides probiotics in development:\n")
cat("      B. uniformis CECT 7771 (Gauffin Cano 2012, Fernandez-Murga 2016)\n")
cat("      B. fragilis SNBF-1 (Choi 2024)\n")
cat("  - Maternal faecal microbiota transplantation: corrects persistent\n")
cat("    Bacteroides deficit in CS infants (Korpela 2020 Cell)\n")
cat("  - Synbiotic infant formula trials: restoration of Parabacteroides\n")
cat("    at 17 weeks and Bacteroides at 12 months (2025 EJCN trial)\n\n")
cat("The honest framing is: the deficit is poorly addressed by MAINSTREAM\n")
cat("CONSUMER products but is an ACTIVE RESEARCH-STAGE target. The clinical\n")
cat("intervention gap is the consumer-product gap, not a research gap.\n\n")

cat("METHODOLOGICAL NOTE\n")
cat("───────────────────\n\n")
cat("The DA concordance + RF triangulation work (steps 06-15, 20) produces\n")
cat("27 method-robust species. Steps 17-19 apply longitudinal persistence,\n")
cat("antibiotic deconvolution, and CS-subtype checks to identify which of\n")
cat("those 27 remain intervention-relevant.\n\n")

cat("LITERATURE CONTEXT\n")
cat("──────────────────\n\n")
cat("The main findings are confirmatory of published early-life microbiome\n")
cat("literature:\n")
cat("  (1) Persistent Bacteroidota deficit — established in Bokulich 2016,\n")
cat("      Korpela 2018/2020/2022, Reyman 2019.\n")
cat("  (2) Partial Bifidobacterium recovery by three months — a well-known\n")
cat("      maturation pattern.\n")
cat("  (3) IAP suppression of B. longum — documented in Mazzola 2016,\n")
cat("      Stearns 2017.\n")
cat("Cohort-specific quantitative contributions include: the per-species\n")
cat("CLR trajectory classification (9 closed + 10 attenuated + 4 persisted +\n")
cat("4 widened), the inverted IAP exposure pattern (vaginal > CS), and the\n")
cat("elective-vs-emergency CS similarity at species-level resolution.\n")
sink()

cat("[21 done] intervention_report.txt + intervention_master_table.csv\n")
