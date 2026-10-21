# PMI estimation — BrainSeq DLPFC

Exploratory **methodological reproduction**, not an exact numerical replication, of Yang et al. (2026), *Forensic Science International* 389:113146, "Transcriptome-based postmortem interval prediction in the human dorsolateral prefrontal cortex using covariate-aware machine learning".

## Status (2026-10-21)

First 5-fold split only; folds 1–3 reached panel evaluation, fold 4 yielded only **2** BH-FDR-significant candidate genes (so the minimum 5-gene panel is unavailable), fold 5 was not run. **Do not report averages over folds 1–3 as complete 5-fold CV.**

| Fold | Train | Test | TMM reference | Candidate genes (BH-FDR <0.05 and abs(r)>0.2) |
|---|---:|---:|---|---:|
| 1 | 213 | 54 | R3374 | 1386 |
| 2 | 213 | 54 | R3991 | 1849 |
| 3 | 214 | 53 | R2999 | 917 |
| 4 | 214 | 53 | R4394 | 2 |
| 5 | 214 | 53 | not analyzed | not analyzed |

Fold 4: 16,968 genes after CPM filtering, 121 genes with abs(r)>0.2 and 2 genes with BH-FDR<0.05. PMI 90 h was held out in test, while fold-4 training max PMI was 68 h; whether this caused candidate collapse has **not** been established.

## Major limitation

Paper: 270 controls, GRCh38/HISAT2/featureCounts/Ensembl release 115. Here: **267 controls** in a separately released BrainSeq Phase 1 hg19/TopHat2/Ensembl v75 dataset. Expression estimates, included samples, and sometimes gene IDs differ.

The local preparation of `rse_final` (gene count `SummarizedExperiment`) and `meta_model` (267 rows, sample-aligned) **was not available in the recorded code** and is not fabricated. Supply the original preprocessing script or prepare these two objects locally before using this repository. The repository code is an **independent reconstruction from the documented method and prior session**, NOT a byte-for-byte export of the user's RStudio history.

## Files

- `scripts/01_validate_and_split.R`: checks prepared objects and generates the original exploratory outer split
- `scripts/02_fold_pipeline.R`: fold-specific CPM filter, TMM, covariate-adjusted association, bootstrap elastic-net ranking, and panel-size evaluation
- `scripts/03_run_and_save.R`: execute folds and save local, non-identifying summary tables
- `results/observed_fold_summary.csv`: observed screening counts from RStudio
- `results/observed_panel_metrics_partial.csv`: observed first ten panel MAE/RMSE for folds 1–3 (not complete)
- `.gitignore`: excludes individual-level data and large local objects

## Re-run

In RStudio (working directory = repo root), load your own `rse_final` and `meta_model` (with matching sample order, gene counts, PMI, Age, RIN, Sex, Manner_Of_Death, snpPC1–snpPC5). Then:

```r
source("scripts/01_validate_and_split.R")
source("scripts/02_fold_pipeline.R")
source("scripts/03_run_and_save.R")
```

IMPORTANT: script 03 is deliberately configured to **run zero folds by default**. Choose the folds to run in that file explicitly. Results can differ from the original RStudio run due to unresolved implementation details and random seeds.

## Implementational assumptions to verify against original scripts

- Train TMM reference: sample with 75th-percentile CPM closest to mean across training samples.
- Test normalization: pairwise TMM against held-out training reference, anchored to training reference norm factor (this is a reconstruction assumption).
- Correlation p-value: df = df.residual(PMI covariate model) - 1.
- Bootstrap elastic net: alpha=0.5, 100 resamples; 5-fold cv.glmnet and lambda.min used for selection (bootstrap lambda choice not explicitly documented in source article).
- Final gene-only panel: cv.glmnet(alpha=0.5, nfolds=5, lambda.min).
- Fold assignment: `set.seed(123)` followed by `sample(rep(1:5, length.out=267))`, not the paper authors' fold assignments.
- The article assessed **five repeats of five outer folds**. This repository captures only **one** outer split; no completion or equivalence is claimed.

**Never change FDR thresholds or reassign the PMI=90 h subject only to make a panel feasible.** A fixed-number ranking would be a separate, explicitly modified analysis, evaluated afresh on all folds.
