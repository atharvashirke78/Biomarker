# Lock notes, diagnostic v1 (tag lock-diagnostic-v1, commit 10bd636)

Written after the lock, in plain words. These notes do not change the lock; they explain three things in this folder.

**(a) Two hash files.** `sha256sums.txt` in this folder is the PRE-LOCK draft record (made when the files were still in results_rc1prep_19/lock_proposal_v2/). It no longer matches the final `lock_statement.md` or the committed config files. `sha256sums_final.txt` is the FINAL record: the sha256 of every file in this folder and of `config/candidate_rules.yaml` and `config/reference_genes.yaml`, computed from git at the tag (`git show lock-diagnostic-v1:<path> | sha256sum`), so they match the tag exactly.

**(b) Words that still say DRAFT.** Line 6 of `lock_statement.md` ("not fitted in this draft") and `step12a_v2_summary.md` ("DRAFT, nothing committed") still use draft wording. The locked versions are the tagged commit and the committed config files (`config/candidate_rules.yaml`, `config/reference_genes.yaml`), not the wording of those two lines.

**(c) Rule 4 model term.** Rule 4 in `config/candidate_rules.yaml` says the model is "~ study + group". The code (`nextflow/bin/deseq2.R`) uses the column name `batch` in the sample sheet. For the current sample sheet batch equals study, so the two are the same model. A different sample sheet where batch and study differ would need this to be clarified.
