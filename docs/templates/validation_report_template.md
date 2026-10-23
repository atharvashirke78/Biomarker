# Independent validation of the locked diagnostic model (lock-diagnostic-v1): {{RUN_FOLDER}}

<!-- TEMPLATE. Placeholders are {{NAME}}. Sections marked FIXED TEXT are copied into every report unchanged (owner instruction 2026-10-16);
     they may only be changed with the owner's approval. Every number in the filled report is labelled VERIFIED (re-read from a named output
     file and recomputed by the second check script) or NOT VERIFIED. -->

**Technical validation on public data. Not clinical validation.** The locked model was applied **once** to the 15 pre-specified independent samples; nothing was tuned before or after. Run folder `{{RUN_FOLDER}}`; code commit `{{COMMIT}}`; date `{{DATE_UTC}}`.

## 1. What was tested (FIXED TEXT)

- Model: tag `lock-diagnostic-v1`, files in `lock/diagnostic_v1/` (9 genes: 8 up, 1 down; cutoff 1.352325902; call tumour if score > cutoff). Lock hashes checked before scoring: {{LOCK_HASH_CHECK}}.
- Pass rule (locked in `config/candidate_rules.yaml`, `classifier_test.pass_rule_15_samples`): at most 1 classification error per class (tumour, normal), and the same within each study.
- AUC is **reported only**; it is not the pass test. Its confidence interval is the DeLong interval produced by `scripts/apply_locked_model.R` (no other interval method is used).
- Rule 4 (candidate reproduction): number and percentage of the 1,390 candidates reproduced; no pass threshold.
- QC exclusion rule (locked): a sample with QC status FAIL is excluded automatically and listed with its reason; WARN samples are kept and listed; the cohort outlier flag is flag-only. No sample is excluded or kept by manual choice.

## 2. Samples (FIXED TEXT, counts filled from the files)

- Pre-specified list: `validation_run/samplesheets/validation_15.csv`, sha256 `14a0bfb9a6378c7fc91269aefc9eafc9b65e8a3df5e6b993cf92aa2402b47893` (git-tracked since commit `1e1f14b`, 2026-10-12, before the lock; identical to the version at tag `lock-diagnostic-v1`).
- Folder-based input: `{{INPUT_DIR}}` (symbolic links to the downloaded FASTQ), `metadata.csv` sha256 `{{METADATA_SHA256}}`, built by script from the list (not typed by hand).
- Counts from the list: SRP344867 4 normal / 5 tumour; SRP245232 3 normal / 3 tumour; total 7 normal / 8 tumour; 9 patients.
- Excluded by QC in this run: {{QC_EXCLUDED}}. Kept with WARN: {{QC_WARN}}. Scored: {{N_SCORED}}.

## 3. Result against the pass rule

| Scope | Tumour errors / n | Normal errors / n | Within rule (<= 1 per class) |
|--|--|--|--|
| All | {{ALL_T}} | {{ALL_N}} | {{ALL_OK}} |
| SRP344867 | {{S1_T}} | {{S1_N}} | {{S1_OK}} |
| SRP245232 | {{S2_T}} | {{S2_N}} | {{S2_OK}} |

**Locked pass rule: {{PASS_OR_FAIL}}.** A failure is reported as it is; nothing is adjusted.

Per-sample scores and calls: `{{SCORES_FILE}}`.

## 4. Reported only

- AUC (not the pass test): {{AUC}} (DeLong 95 % CI {{AUC_CI}}). See the limitations in section 7 on paired samples.
- Rule 4: {{RULE4_N}} of 1,390 candidates reproduced ({{RULE4_PCT}} %).

## 5. Reference genes (report only; FIXED TEXT, values filled)

The 10 frozen reference genes (`config/reference_genes.yaml`, status LOCKED) are **not** used in the score, the cutoff or the pass/fail. Their stability on the 15 samples is reported separately: {{REFGENE_SUMMARY}} (`{{REFGENE_FILE}}`). Using reference genes inside a score would be a new model version with its own lock and its own independent test, never a change to `lock/diagnostic_v1`.

## 6. Second check

{{SECOND_CHECK}} (independent script recomputing scores, calls, per-class and per-study errors and rule 4 from the counts and the lock files).

## 7. Independence and limitations (FIXED TEXT)

1. **Independence from the training samples: VERIFIED at record level.** The training set used studies SRP219837 and SRP301216. None of the 15 samples shares a run, GEO sample or project with them, and none of the 15 was used to choose genes, rules or the cutoff.
2. **Person-level overlap with the training samples: NOT VERIFIED.** Whether any person appears in both cannot be checked from public metadata; there is no sign of it.
3. **SRP344867 is part of the 15** (9 of the 15 samples); it was never used in training, so it is part of the independent test, not a separate check.
4. **Paired samples: 6 of the 9 patients contribute a tumour-normal pair** (SRP344867: 46676, 450585, 380321; SRP245232: CRC1, CRC2, CRC3). This does not affect the pass rule, which counts errors per sample, but the AUC confidence interval assumes more independence than the data have, so it is likely too narrow. The AUC stays report-only.
5. **Not independent of Hosseini 2023:** all 15 runs are listed in that paper's supplement. Comparisons with Hosseini's gene lists on these samples would be circular and are not made.
6. **TCGA overlap: NOT VERIFIED.** The metadata of the 15 does not mention TCGA; patient-level overlap cannot be checked.
7. The test set is small (15 samples); all intervals are wide. Results are technical, on public data, not clinical validation.
