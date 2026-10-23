# Third-party software and data used by the pipeline

Status of this file: working inventory for the licence and ownership review (upgrades.md, H4f). **It is not legal advice. Counsel must confirm every entry before anything is hosted, shipped or sold.** No licence is declared here for the pipeline's own code (owner decision pending).

How to read the "Status" column:
- **verified (upstream)**: the licence text of the version-tagged upstream repository was read by the research step that supplied this table. Not re-read when this file was written.
- **installed metadata**: taken from the conda package metadata (`conda-meta/*.json`, field `license`) or the installed distribution of the package in `env/prefix/`. This is the packager's statement and can be wrong (example: the bioconda metadata of STAR says GPL, the upstream tag says MIT). It has not been compared with upstream.
- **not verified**: no source read. Needs checking.

Versions are those installed under `env/prefix/` (conda-meta) and recorded in `contract/manifest_quant.json` of the demo run.

## Programs

| Name | Installed version | Licence | Source of the statement | Status |
|---|---|---|---|---|
| STAR | 2.7.11b | MIT | upstream tag `2.7.11b`, LICENSE. Bioconda metadata says GPL-3.0-or-later, which is wrong for this tag | verified (upstream); the conda metadata disagrees |
| Salmon | 1.10.3 | GPL-3.0 | upstream tag `v1.10.3`, LICENSE. Later Salmon master is BSD-3-Clause: do not assume that for 1.10.3 | verified (upstream) |
| subread (featureCounts) | 2.1.1 | GPL-3.0 | upstream release; conda metadata GPL-3.0-only | verified (upstream); installed metadata agrees |
| fastp | 1.4.0 | MIT | conda metadata (MIT) | installed metadata (owner brief lists MIT as verified) |
| samtools | 1.24 | MIT | conda metadata (MIT) | installed metadata (owner brief lists MIT as verified) |
| Nextflow | 24.10.5 | Apache-2.0 | upstream | verified (upstream, per owner brief) |
| MultiQC | 1.35 | GPL-3.0 | `LICENSE` file inside the installed distribution (`multiqc-1.35.dist-info/licenses`) and the trove classifier "GPLv3" | verified (installed distribution), not checked against the upstream tag |
| FastQC (only if `run_fastqc` is on; off by default) | 0.13.0 | GPL >= 3 | conda metadata | installed metadata |
| micromamba (environment manager; not part of any result) | 2.9.0 | not verified | - | not verified |
| OpenJDK (Nextflow runtime) | 17.0.18 (nf env) | GPL-2.0 with Classpath exception (conda metadata for the JDK 25 package in the fastqc env; the JDK 17 package was not read) | conda metadata | not verified for the version in use |
| Python | 3.12.15 (pyqc, multiqc) | Python-2.0 (PSF) | conda metadata | installed metadata |
| numpy / pandas (pyqc) | 2.5.3 / 2.3.3 | BSD-3-Clause | conda metadata | installed metadata |
| R | 4.5.3 | GPL-2.0-or-later | conda metadata | installed metadata |

## R packages (names as in `contract/manifest_quant.json` and the scripts in `nextflow/bin/`)

| Name | Installed version | Licence | Source of the statement | Status |
|---|---|---|---|---|
| DESeq2 | 1.50.2 | LGPL (>= 3) | conda metadata, equals owner brief | verified (upstream, per owner brief) |
| apeglm | 1.32.0 | GPL-2 | conda metadata, equals owner brief | verified (upstream, per owner brief) |
| tximport | 1.38.2 | LGPL (>= 2) | conda metadata, equals owner brief | verified (upstream, per owner brief) |
| fgsea | 1.36.2 | MIT + file LICENCE | conda metadata; owner brief says MIT | verified (upstream, per owner brief) |
| clusterProfiler | 4.18.4 | Artistic-2.0 | conda metadata, equals owner brief | verified (upstream, per owner brief) |
| enrichplot | 1.30.5 | Artistic-2.0 | conda metadata | installed metadata |
| org.Hs.eg.db | 3.22.0 | Artistic-2.0 | conda metadata (annotation data package: the underlying NCBI/Ensembl data terms were not checked) | installed metadata; data terms not verified |
| GO.db | 3.22.0 | Artistic-2.0 | conda metadata (the Gene Ontology data itself is CC BY 4.0, see below) | installed metadata |
| BiocParallel | 1.44.0 | GPL-3.0-only | conda metadata | installed metadata |
| glmnet | 4.1-10 | GPL-2 | conda metadata (GPL-2.0-only), equals owner brief | verified (upstream, per owner brief) |
| pROC | 1.19.1 | GPL (>= 3) | conda metadata (GPL-3.0-or-later), equals owner brief | verified (upstream, per owner brief) |
| msigdbr | 26.1.1 | MIT (the package); the MSigDB data it ships has its own terms, see below | conda metadata, equals owner brief | verified (upstream, per owner brief) |
| ggplot2 | 4.0.3 | MIT | conda metadata | installed metadata |
| ggrepel | 0.9.8 | GPL-3.0-only | conda metadata | installed metadata |
| pheatmap | 1.0.13 | GPL-2.0-only | conda metadata | installed metadata |
| jsonlite | 2.0.0 | MIT | conda metadata | installed metadata |
| readr | 2.2.0 | MIT | conda metadata | installed metadata |
| dplyr | 1.2.1 | MIT | conda metadata | installed metadata |
| data.table | 1.18.6.1 | MPL-2.0 | conda metadata | installed metadata |
| yaml | 2.3.12 | BSD-3-Clause | conda metadata | installed metadata |
| ReactomePA | NOT INSTALLED | - | not used in the current run (recorded as "NOT INSTALLED") | - |

R packages that the scripts load but whose licence was not read in this inventory: none found beyond the table (checked with a search for `library(...)` in `nextflow/bin/*.R` and `scripts/*.R`); transitive dependencies of these packages (hundreds of conda packages) were NOT reviewed. A full software bill of materials from `conda-meta` is needed for a complete review.

## Data and reference files

| Name | Licence or terms | Source of the statement | Status |
|---|---|---|---|
| GENCODE v50 (annotation GTF, transcripts FASTA) | The GENCODE data page states no licence; EBI terms of use apply | owner brief | verified that none is stated on the data page (per owner brief); terms need counsel |
| GRCh38 primary assembly FASTA | not verified | - | not verified |
| MSigDB (Hallmark and other collections, via msigdbr) | CC BY 4.0 for some versions, with exceptions for some collections; the version in use must be checked | owner brief | verified (per owner brief) for the general statement; exact version and collection terms not checked |
| Gene Ontology (GO) | CC BY 4.0 | owner brief | verified (per owner brief) |
| KEGG | needs a commercial licence for use beyond academic use. The pipeline keeps KEGG off (`run_kegg = false`); do not turn it on in a hosted or commercial service | owner brief | verified (per owner brief) |
| Reactome | not used in the current run | - | - |
| Input data of users (FASTQ, metadata) | belong to their owners; no licence is granted by this pipeline | - | - |

## Open points for the owner and counsel

1. Our own code is MIT-licensed (`LICENSE`, 2026-10-23, owner decision). Ownership of AI-assisted code and the licences of third-party tools run as separate processes still need counsel's confirmation.
2. Salmon 1.10.3 and subread are GPL-3.0, MultiQC and several R packages are GPL. The pipeline calls them as separate programs (command line, `Rscript`). Whether and how that constrains distribution of a bundle or hosted service must be decided by counsel; this file makes no such statement.
3. The Salmon 2.x line was reported as BSD-3 on master and as "deterministic" in its release notes; if the pipeline moves to it, re-read the LICENSE of that tag and re-validate the numbers (upgrades.md H1b, H4g).
4. Check the exact MSigDB release behind msigdbr 26.1.1 and the terms of any collection used.
5. Complete the "not verified" rows and produce a full dependency list.

Nothing in this file is legal advice.
