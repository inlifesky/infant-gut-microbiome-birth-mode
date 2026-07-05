# Caesarean-birth microbiome signals: a stratified re-analysis of the Baby Biome Study

A reproducible R workflow that re-analyses public shotgun-metagenomic profiles
from the Baby Biome Study (Shao et al., *Nature* 2019) to test which
caesarean-birth-associated infant gut microbiome signals remain after
(a) method-robustness across four differential-abundance approaches,
(b) longitudinal persistence to three months, and
(c) deconfounding from peripartum antibiotic exposure.

**Author:** Wei (Rita) Yuan

---

## Rendered reports

Two HTML report pages are included in this repository:

| Page | Purpose | Location |
|---|---|---|
| **Portfolio** | Short case-study version aimed at hiring panels and general readers | [`portfolio/index.html`](portfolio/index.html) |
| **Full scientific report** | Detailed methods, per-filter results, pathway view, glossary and references | [`web/index.html`](web/index.html) |

Both are static self-contained pages and can be viewed by opening the HTML
files directly or served locally (`python -m http.server` at the repo root).

---

## Repository layout

```
src/                      26 ordered R scripts + run_all.R
results/                  CSV output tables + session info
figures/                  PNG figures produced by the pipeline
web/                      Full scientific report (HTML, self-contained)
portfolio/                Portfolio case-study page (HTML, self-contained)
README.md                 This file
LICENSE                   MIT
.gitignore
```

---

## Data source

All input data are downloaded at runtime from the
[`curatedMetagenomicData`](https://doi.org/10.1038/nmeth.4468)
Bioconductor package (MetaPhlAn3 species and HUMAnN3 MetaCyc pathway profiles
of Shao et al. 2019). No raw sequencing reads are stored in the repository.

---

## Reproducing the analysis

### Requirements

- R ≥ 4.6 with Bioconductor 3.21
- CRAN packages: `vegan`, `randomForest`, `ggplot2`, `patchwork`, `dplyr`,
  `tidyr`, `ggrepel`
- Bioconductor packages: `curatedMetagenomicData`, `SummarizedExperiment`,
  `ANCOMBC`, `ALDEx2`, `mia`, `Maaslin2`, `microbiome`

### Run

```bash
git clone <repo-url>
cd <repo-name>
# From the repo root (this working directory matters — scripts use getwd()):
Rscript src/run_all.R
```

`run_all.R` executes the 26 analysis scripts in order and writes
`results/sessionInfo.txt`. The first execution downloads the cohort via
`curatedMetagenomicData` (~2 GB, cached by `ExperimentHub` on subsequent runs);
end-to-end runtime is approximately 45 minutes on a laptop.

### Individual scripts

Each script in `src/` is self-contained and can be run independently after
`01_load_clean.R` (species pipeline) or `01b_load_full_longitudinal.R`
(longitudinal pipeline) has been executed once. Scripts read intermediate
`.rds` files from `results/`. See headers within each script for details.

---

## Pipeline overview

The full pipeline is described in `web/index.html` (Analysis pipeline section)
and `portfolio/index.html` (What I did section). In summary:

1. **Profiles** — MetaPhlAn3 species and HUMAnN3 MetaCyc pathway abundances,
   accessed via `curatedMetagenomicData`
2. **Differential abundance** — four independent methods run in parallel on
   the same profile (CLR + Wilcoxon, ALDEx2, ANCOM-BC2, MaAsLin2)
3. **Filter 1** — species significant in ≥3 of 4 methods (FDR < 0.05)
4. **Filter 2** — longitudinal persistence to day 91+
5. **Filter 3** — antibiotic-versus-CS deconvolution (2×2 stratification)
6. **Filter 4** — direction filter for restore-target framing
7. **Pathway repetition** — the same four filters applied at the functional level

---

## Development workflow

This project was author-directed and developed with AI assistance.
Claude (Anthropic) was used for code drafting, exploratory literature lookup
and technical writing. The research question, filter design, statistical
validation, error verification and scientific interpretation were performed
by the author. Numerical outputs were verified against source data files,
and literature claims were checked against primary references before
inclusion in either report page.

---

## Acknowledgments

This project is a secondary analysis of publicly available shotgun metagenomic
profiles from the Baby Biome Study. The author gratefully acknowledges:

- **The mothers and infants** who consented to the Baby Biome Study, and the
  Shao et al. team at the Wellcome Sanger Institute and the Wellcome–MRC
  Institute of Metabolic Science, whose cohort makes this re-analysis possible.
- **The [`curatedMetagenomicData`](https://waldronlab.io/curatedMetagenomicData/)
  maintainers** (Pasolli et al.) for providing uniformly-processed
  MetaPhlAn3 / HUMAnN3 profiles as a Bioconductor package.
- **The methodology teams** whose tools this analysis relies on:
  MetaPhlAn3 / HUMAnN3 (Segata group), ALDEx2 (Gloor group), ANCOM-BC2
  (Peddada group), MaAsLin2 (Huttenhower group).

No re-identifiable participant data are included in this repository; all
downstream results are group-level aggregates.

---

## Licence

MIT — see `LICENSE`.

---

## Citation

If you use this workflow or its filter design for your own analysis, please
cite the source cohort:

> Shao Y, Forster SC, Tsaliki E, et al. Stunted microbiota and opportunistic
> pathogen colonization in caesarean-section birth. *Nature* 574, 117–121
> (2019). [doi:10.1038/s41586-019-1560-1](https://doi.org/10.1038/s41586-019-1560-1)
