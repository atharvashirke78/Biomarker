# Install plan: glmnet and pROC (PLAN ONLY, nothing installed)

Package search (micromamba 2.9.0, conda-forge, linux-64), dry-run solve with `--dry-run` (no download, no env created):
- r-glmnet 5.0 (build r45ha36cffa_0) and 4.1_10 (r45ha36cffa_1) exist for R 4.5. PROPOSED: 4.1_10 is the older, long-used line; 5.0 is a new major version. Owner to choose; default proposal 5.0 as solved unless you prefer the 4.x line.
- r-proc 1.19.1 (r45hfa9885e_0); 1.19.0.1 also exists.
- Solver also pulls r-survival 3.8_12 (dependency of glmnet) plus r-foreach, r-iterators, r-shape, r-matrix, r-rcpp. So survival arrives as a side effect; ReactomePA and reactome.db are NOT included.

Proposed env: NEW env `r_classifier` at `/workspace/Biomarker_pipeline/env/prefix/r_classifier` (r-base=4.5.3, r-glmnet, r-proc). Reason: do not touch r_dge (frozen for the DESeq2 results) or r_pathway. 101 packages in the solve; 17 MB to download, the rest is already in the cache.
Disk: download 17 MB. Env size: not measured in a dry run; files are hard-linked from the cache on the same filesystem, so extra space is small (<~1 GB worst case if copied; r_dge for comparison is 4.2 GB with Bioconductor).
Space now: /workspace 379 TB free (network FS), root disk 4.7 GB free of 5.0 GB.
Package cache: `CONDA_PKGS_DIRS=/workspace/Biomarker_pipeline/env/mamba/pkgs` (already 19 GB in use there), must be set on every call because micromamba otherwise also lists /root/.mamba/pkgs. Envs dir is /workspace/envs by default; I use `-p` with an explicit prefix, so nothing goes to root.
Pinning: `env/r_classifier.yml` (exact versions, same style as env/r_dge.yml) plus `env/r_classifier.lock.txt` from `micromamba env export --explicit` (URLs and md5). Recorded in docs/code_sources.md and docs/tool_versions.tsv (version, channel, date, lock md5). Note: the other envs have only .yml files, no lock files.
Side effect already caused (read-only check): msigdbr wrote its gene-set cache to /root/.cache/R/msigdbr (85 MB, root disk). Proposed: set `R_USER_CACHE_DIR=/workspace/Biomarker_pipeline/env/cache` for later R runs and, with your OK, move that folder there.
