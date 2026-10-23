# DRAFT: validation scoring ledger (not built, not approved; nothing is scored by this document)

Purpose: the 15 held-out validation samples are scored with the locked model **once**. This ledger makes a second scoring, or any use of these samples in a discovery run, impossible by accident. It protects against mistakes, not against someone who deliberately edits files.

## 1. Three files

| File | Content | Where |
|---|---|---|
| **Registry** `validation_samples.tsv` | the 15 validation samples: `sample_id, group, study, R1_bytes, R2_bytes, R1_md5, R2_md5` (from `validation_run/samplesheets/validation_15.csv` and the ENA file report) | in git (`validation/`), its sha256 goes into the freeze record |
| **Ledger** `scoring_ledger.tsv` | append-only event log, one row per event (section 3) | **decision needed:** in the shared install folder (`<install_dir>/validation_ledger/`, so every checkout and worktree sees the same file) or inside git. Recommendation: install folder as the working copy, and a copy committed into the freeze record after the scoring |
| **Frozen commit** `FROZEN_COMMIT.txt` | the commit that may score the validation samples | in git, written at freeze |

## 2. Rules (checked at the start of every run, and again right before the scoring task)

| # | Rule | Message when it stops |
|---|---|---|
| R1 | A sample is a validation sample if its ID is in the registry **or** its FASTQ sizes match a registry entry (so renaming does not hide it). A validation sample may only be used with `--mode apply --cohort validation` | "SRR... is one of the 15 validation samples; it must not be used in discover mode or with --cohort development" |
| R2 | `--cohort validation` needs **exactly** the 15 registry samples (no subset, no extra sample) | "validation scoring needs all 15 registered samples; found 9 of them. Scoring part of the set is not allowed" |
| R3 | If the ledger already holds `SCORING_STARTED` for this sample set, the run stops | "the validation samples were already scored on 2026-..-.. in /path (commit abc1234). A second scoring is not allowed. Owner override: --validation_override '<reason>' (recorded in the ledger)" |
| R4 | **Resume is allowed** only into the same `output_dir` as the ledger row (an interrupted scoring can finish; it is not a second scoring) | |
| R5 | `--cohort validation` needs a clean git state (0 uncommitted files) and HEAD equal to `FROZEN_COMMIT.txt` | "the code is not the frozen commit" |
| R6 | After quantification the FASTQ md5 in `manifest_quant.json` must equal the registry md5 | "FASTQ of SRR... differs from the registered file" |
| R7 | `--stop_after quant` is not a scoring: it writes no ledger row and can be repeated | |

## 3. Ledger format (tab-separated, append only)

```
event_utc            event               run_name   output_dir        commit   lock_sha256   sample_set_sha256   n   operator   note
2026-..-..T..:..:..Z SCORING_STARTED     <name>     /path/val15_run1  abc1234  37657d53...   9f3c...            15  <user>     written BEFORE any score exists
2026-..-..T..:..:..Z SCORING_COMPLETED   <name>     /path/val15_run1  abc1234  37657d53...   9f3c...            15  <user>     scores.tsv sha256 <hash>
2026-..-..T..:..:..Z OVERRIDE            <name>     /path/val15_run2  abc1234  ...                               15  <user>     reason text from --validation_override
```
- `SCORING_STARTED` is written **before** the scoring task computes anything. A crash after that still counts as "used": the only ways forward are R4 (resume into the same folder) or an explicit override.
- `sample_set_sha256` = sha256 of the sorted registry sample IDs; `lock_sha256` = sha256 of `lock/diagnostic_v1/sha256sums_final.txt`.
- The file is written atomically (temp file, then rename, under a file lock) and its sha256 is stored in `run_info.json` and in the freeze record.

## 4. Where it is built (when approved)

1. `nextflow/bin/ledger_guard.py` with three commands: `check` (read only, at the start), `start` (append `SCORING_STARTED`), `complete` (append `SCORING_COMPLETED` with the scores hash).
2. Groovy call in `folder_input.nf` at the start of `RUN_PIPELINE`; the `LOCKED_SCORE` process calls `start` just before the R script and `complete` after it.
3. Parameters: `--validation_registry`, `--validation_ledger` (defaults as in section 1), `--validation_override '<reason>'`.

## 5. Tests (no real scoring; fixtures with fake IDs in a temporary folder)

first scoring allowed · second scoring stops (R3) · resume into the same folder allowed (R4) · subset stops (R2) · a registry sample in discover mode or `--cohort development` stops (R1) · renamed ID with matching FASTQ sizes stops (R1) · dirty git or wrong commit stops (R5) · override is recorded · `--stop_after quant` writes no row (R7). Each gets a PASS line in the stub suite.

## 6. Decisions needed from the owner

1. Ledger location: install folder (recommended) or inside git.
2. Is "exactly the 15 samples" (R2) right, or may a QC-failed sample be missing from the scored set? (Proposal: the 15 must all be **present in the input**; QC-failed ones are excluded by the normal QC rule and listed as not evaluable.)
3. Who may use `--validation_override`, and is a second scoring ever allowed?
