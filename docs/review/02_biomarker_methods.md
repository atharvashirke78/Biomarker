# Domain 2 of 4: Methods and tools for building and honestly validating RNA-seq gene-signature classifiers

Review date: 2026-10-23. Reviewer: Claude agent (Sonnet 5.5). Status: protocol written BEFORE any search.

## 1. Protocol (fixed before searching)

Question: what methods and tools exist for building and honestly validating RNA-seq gene-signature classifiers (selection bias, nested CV, stability selection, cutoff selection, optimism correction, external validation, locked models), and how does our pipeline compare?

Our pipeline (comparator): DESeq2 candidate rules (padj, |shrunk LFC|>=1, baseMean, two-route concordance, leave-one-out >=90%) -> glmnet LASSO lambda.1se -> balanced up/down score (mean log2CPM up minus down) -> median-of-inner-fold Youden cutoff. Internal validation: repeated nested CV of the LASSO+cutoff stage only; Harrell bootstrap optimism; bootstrap cutoff interval; selection frequency. Locked 9-gene model applied to independent samples with a pre-registered rule.

Sources and planned exact search strings (date of each search recorded in section 2):
- PubMed (via PubMed MCP search_articles): 
  - Q1 `selection bias cross-validation gene expression classifier`
  - Q2 `nested cross-validation microarray OR RNA-seq classifier bias`
  - Q3 `MAQC-II OR SEQC2 predictive model development guidelines`
  - Q4 `optimal cutpoint estimation biomarker bootstrap`
  - Q5 `stability selection high-dimensional genomic`
  - Q6 `optimism bootstrap internal validation prediction model`
  - Q7 `RNA-seq classifier small sample size gene signature validation`
  - Q8 `multi-cohort gene signature validation meta-analysis`
- Europe PMC: web search restricted to europepmc.org with the same themes (no direct API tool available; use WebSearch/WebFetch of europepmc.org REST if reachable).
- bioRxiv: MCP tool has NO keyword search (category/date only); will use WebSearch restricted to biorxiv.org with the Q1-Q8 themes. Limitation declared.
- CRAN, Bioconductor, GitHub: WebSearch / WebFetch of package pages (cran.r-project.org/package=X, bioconductor.org/packages/X, github.com releases). Search strings: `CRAN nested cross-validation`, `CRAN stability selection`, `CRAN cutpoint optimal cutpoint`, `Bioconductor gene signature classifier`, `GitHub nested cross validation RNA-seq`.

Inclusion: methods / benchmarking / guideline papers, or maintained software (release within last 3 years, i.e. since 2023-10-23), on signature or classifier building from expression data (or general-purpose methods directly applicable to it). Foundational papers older than 3 years are included as methods references (the 3-year rule applies to software only).
Exclusion: single-disease claims without methods content; software without release in 3 years; duplicates; non-expression-data tools.
Screening target: at least 50 records screened; log identified / screened / included / excluded with reasons.
Extraction per item: aim, input, output, bias control, evidence of performance, licence, last release, URL.
Verification marking: each fact V (page or paper opened this session) or S (search snippet / metadata only). Prefer V. Pages that failed to open are listed.
Required topics: selection bias in CV (Ambroise & McLachlan; Varma & Simon), MAQC-II/SEQC2 model-building guidance, cutpoint estimation, stability selection, Harrell optimism bootstrap, multi-cohort validation, small-sample RNA-seq classifiers.

---

## 2. Execution log (all searches run on 2026-10-23)

| Source | Exact string / route actually used | Result |
|---|---|---|
| PubMed (MCP search_articles) | `selection bias cross-validation gene expression classifier`; `nested cross-validation microarray OR RNA-seq classifier bias`; `MAQC-II OR SEQC2 predictive model development guidelines`; `optimal cutpoint estimation biomarker bootstrap`; `stability selection high-dimensional genomic classifier`; `optimism bootstrap internal validation prediction model variable selection`; `RNA-seq classifier small sample size gene signature validation overfitting`; `multi-cohort gene signature validation meta-analysis framework` | 74 PMIDs returned (totals 24/204/24/3/24/226/0/6). Titles read for 33; the other 41 PMIDs were NOT read (the tool returns IDs only; metadata dump is large). PubMed web pages (pubmed.ncbi.nlm.nih.gov) failed to open (cookie wall). |
| Europe PMC (REST, via WebFetch) | 5 queries, sorted by citations, pageSize 25: (1) `"nested cross-validation" AND (RNA-seq OR "gene signature") AND (selection bias OR optimism)`; (2) `("optimal cutoff" OR "optimal cutpoint") AND (bias OR bootstrap) AND (biomarker OR expression)`; (3) `("stability selection" OR "selection stability") AND (RNA-seq OR transcriptom* OR "gene expression") AND (classifier OR signature)`; (4) `("multi-cohort" OR "multicohort") AND "gene expression" AND (signature OR classifier) AND (validation OR framework) AND (method* OR benchmark*)`; (5) `("small sample" OR "limited sample") AND RNA-seq AND (classifier OR "gene signature") AND (benchmark* OR comparison OR overfitting)` | 125 records listed, titles screened; precision very low (sorting by citations pulled unrelated highly cited papers); about 8 relevant. |
| bioRxiv | MCP tool has no keyword search (date/category only), so not usable for this. Substitute: WebSearch restricted to biorxiv.org/medrxiv.org: `bioRxiv nested cross-validation RNA-seq biomarker signature information leakage preprint` | 9 links, 5 relevant; 2 preprint pages failed to open (HTTP 429). |
| CRAN | WebFetch of cran.r-project.org/package=X for nestedcv, cutpointr, stabs, sharp, pminternal, glmnet, rms, mlr3, caret, rsample, pROC, OptimalCutpoints, MetaIntegrator. No CRAN full-text search tool available; candidates came from my prior knowledge plus web searches, which is a selection bias of its own. | 13 pages opened |
| Bioconductor | WebFetch of bioconductor.org/packages/X for ClassifyR, DaMiRseq, MLSeq, mixOmics, MetaIntegrator, multiclassPairs | 6 opened; MetaIntegrator and multiclassPairs URLs returned a "removed packages" index, not a package page |
| GitHub | WebSearch (allowed_domains github.com): `github nested cross-validation RNA-seq gene signature classifier lasso Youden bootstrap R package release`; WebFetch of github.com/Ghrieb/rnaSentry and github.com/myles-lewis/nestedcv | 9 links; mostly study-specific analysis repos (excluded); 2 package repos opened |
| General web (for foundational papers) | 17 WebSearch queries naming the paper (Ambroise; Varma; MAQC-II; SEQC; Meinshausen; Harrell/Steyerberg; Krstajic; Bernau; Altman/Leeflang; Sweeney/Khatri; TRIPOD+AI; RNA-seq small sample; frozen normalization; Tibshirani; Riley&Collins; etc.) | about 150 links |

### Record flow (approximate; I did not de-duplicate programmatically)
- Identified: about 350 (74 PubMed + 125 Europe PMC + about 150 web/software links).
- Screened (title/snippet actually read): about 215 (33 + 125 + about 60 distinct web titles). The 41 unread PubMed PMIDs are counted as identified, not screened.
- Included in the extraction tables: 41 (23 papers/guidelines/preprints, 18 software items).
- Excluded: about 175 after screening. Reasons: single-disease or single-signature claim without methods content (the large majority of PubMed hits, e.g. AD 3-gene classifier, HCC 10-gene signature, MM nomograms); off-topic (Europe PMC citation-sorted noise: CRISPR, guidelines, microbiome); software with no release in the last 3 years (OptimalCutpoints 1.1-5, 2021-10-07); software removed from repository (MetaIntegrator removed from CRAN 2024-10-14 per CRAN page; multiclassPairs apparently in the Bioconductor removed list, unconfirmed); study-specific analysis repos without a maintained package; duplicates.
- Screening-target check: >=50 screened, met.

Verification key: V = I opened the page/paper full text or the vendor page this session and the fact comes from it (the fetch tool summarizes pages with a small model, so "V" means page content, not that I read every word; two long pages were read only to character 100000). S = search snippet, title or metadata only.

Pages that failed to open: pubmed.ncbi.nlm.nih.gov/{31985777, 22205940, 27891226, 20676074, 18281551} (cookie wall); link.springer.com Genome Biology 2015 SEQC neuroblastoma (redirect to login); jitc.bmj.com Dobbin 2016 (403); jclinepi.com Steyerberg 2003 (403); biorxiv.org Parvandeh 2019, leakage preprint 2026.02.05.704016, LOL package preprint (HTTP 429); arxiv 2606.01566 (opened but body unreadable; only title, author, licence); Altman 1998 letter page opened but body not shown; multiclassPairs and MetaIntegrator Bioconductor pages (index of removed packages only).

## 3. Foundational methods and guidance (papers)

Columns condensed: Aim / Bias control / Evidence / Licence. "Input/Output" is omitted where it is simply "expression matrix -> error estimate".

| # | Item (URL) | Aim | Bias control and claims | Evidence | Licence | V/S |
|---|---|---|---|---|---|---|
| P1 | Ambroise & McLachlan, PNAS 2002;99:6562 (https://pmc.ncbi.nlm.nih.gov/articles/PMC124442) | Show gene-selection bias in CV error for expression classifiers | Gene selection must be repeated inside each training fold ("external" CV); prefers 10-fold over LOO; .632+ bootstrap external to selection | Colon and leukemia data; permuted labels still gave near-zero uncorrected error for 128 genes, corrected error about 0.40-0.45 | not stated on page opened | V |
| P2 | Varma & Simon, BMC Bioinformatics 2006;7:91 (https://pmc.ncbi.nlm.nih.gov/articles/PMC1397873) | Bias when CV error used both to tune and to report | Nested CV (tuning inside outer loop) nearly unbiased | Null data (true error 50%): tuned CV error 37.8% (shrunken centroids), 41.7% (SVM); nested 54.2% | CC BY 2.0 | V |
| P3 | MAQC Consortium (Shi et al.), Nat Biotechnol 2010;28:827 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3315840) | Model-building and validation practice, 36 teams, 6 datasets, 13 endpoints, >30,000 models | Feature selection and tuning inside each internal validation iteration; at least one external validation; do not reuse a validation set after modifying protocol; document the process; most teams used 10 repeats of 5-fold CV | Endpoint and team skill dominated; unbiased CV tracked external performance | not stated | V (first 100k of 117k characters) |
| P4 | SEQC/MAQC-III neuroblastoma, Genome Biol 2015;16:133 (https://link.springer.com/article/10.1186/s13059-015-0694-1) | RNA-seq vs microarray models for 6 clinical endpoints, 498 tumours, 360 models | Fixed train/validation split | Endpoint nature mattered most; platform/pipeline/feature level did not | not checked | S |
| P5 | Meinshausen & Buhlmann, JRSS-B 2010;72:417 (doi 10.1111/j.1467-9868.2010.00740.x) | Stability selection with error control | Subsample half the data, keep variables with selection frequency above threshold; PFER bound q^2/((2pi-1)p) under exchangeability assumption | theory; I did not see the paper's own proof | n/a | S |
| P6 | Shah & Samworth 2013 (complementary pairs stability selection) | Tighter bounds | via stabs/sharp docs | - | - | S |
| P7 | Steyerberg et al., J Clin Epidemiol 2003;56:441; Harrell, Lee & Mark 1996 | Bootstrap optimism correction | Bootstrap must repeat predictor selection each replicate, otherwise optimism is underestimated | simulation (per secondary summary) | n/a | S (primary papers not opened; Harrell 1996 not found in search) |
| P8 | Altman, Lausen, Sauerbrei, Schumacher, JNCI 1994;86:829; Altman letter BJC 1998 (https://pmc.ncbi.nlm.nih.gov/articles/PMC2063091) | Dangers of data-driven "optimal" cutpoints | min-p cutpoint gives about 40% false-positive rate vs nominal 5% (per secondary summary); cutpoints differ across studies | - | n/a | S for the numbers; V that the 1998 letter exists and cites Faraggi & Simon 1996 on CV for cutpoint choice |
| P9 | Leeflang, Moons, Reitsma, Zwinderman, Clin Chem 2008;54:729 | Bias in Se/Sp from data-driven optimal cutoffs | as title: bias mechanisms, magnitude, solutions | not seen | n/a | S |
| P10 | Thiele & Hirschfeld, cutpointr paper (https://ar5iv.labs.arxiv.org/html/2002.09209) | Estimate and validate optimal cutpoints | Out-of-bag bootstrap validation of the cutpoint; bootstrap/bagged cutpoint; 1000-2000 repeats for CIs | In-sample cutpoint optimistic (example in-bag 1.76 vs OOB 1.72); across 108 simulated scenarios bootstrapping best for small samples/small separation | not stated | V |
| P11 | Krstajic et al., J Cheminform 2014;6:10 (PMC3994246) | CV pitfalls in selection and assessment | Repeated grid-search CV; repeated nested CV (variance depends on split) | - | n/a | S |
| P12 | Tibshirani & Tibshirani, AoAS 2009;3:822; Bernau/Boulesteix weighted-mean correction (CMA package) | Cheap correction of minimum-CV-error bias | bias estimate from CV curves, no extra fits | later work: removes downward bias, sometimes over-corrects | n/a | S |
| P13 | Haury, Gestraud, Vert, PLoS One 2011 (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0028210) | Feature selection vs accuracy, stability, interpretability in 4 breast cancer sets | stability measured under hard perturbations; sample size main driver of instability | filters more stable than wrapper/embedded; t-test best accuracy+stability; ensembles did not beat it | CC BY | V |
| P14 | Parvandeh et al., Bioinformatics 2020 "Consensus features nested CV" (bioRxiv 10.1101/2019.12.31.891895) | Combine feature selection and classification inside nested CV | consensus of inner-fold features | not opened | - | S |
| P15 | Riley & Collins, Biom J 2023;65:2200302 (https://pmc.ncbi.nlm.nih.gov/articles/PMC10952221) | Stability of prediction models in small samples | Repeat whole model-building on >=200 (about 1000) bootstraps; instability, calibration and classification instability plots | LASSO does not remove instability in small samples; n=752 with 53 events and 27 predictors still unstable | CC BY 4.0 | V (first 100k of 125k characters) |
| P16 | Wu, Huang, Qin, JCO Precis Oncol 2021;5:PO.21.00103 (https://pmc.ncbi.nlm.nih.gov/articles/PMC8345919) | Effect of handling variation and normalization of test data on external validation of molecular classifiers | Frozen (training-derived) normalization of test data; uniform handling; report methods/code | miRNA arrays, 192 tumours, 100 runs/scenario; frozen median/quantile matched the no-handling-effect error, conventional test-set normalization often did not help; microarray data, not RNA-seq | CC BY 4.0 | V |
| P17 | TRIPOD+AI, BMJ 2024;385:e078378 (https://pmc.ncbi.nlm.nih.gov/articles/PMC11019967) | 27-item reporting guideline | Item 12c all model-building steps incl. tuning and internal validation; item 15 how classification thresholds were identified; item 22 full model details for third-party evaluation; item 10 study size | guideline | CC BY (page text inconsistent on commercial reuse) | V (first 100k of 169k) |
| P18 | Lewis et al. (nestedcv paper), Bioinformatics Advances 2023;3:vbad048 (https://pmc.ncbi.nlm.nih.gov/articles/PMC10125905) | Fast nested CV with filters for P>>n | Filters run inside outer loop on training data only | 50,000 noise predictors: filter-on-all-data gave inflated AUC, nested CV about 0.50; R4RA n=68: nested AUC 0.783 vs inflated inner pooled 0.863 | CC BY 4.0 | V |
| P19 | Sweeney et al., Sci Transl Med 2015;7:287ra71; Sweeney Lancet Respir Med 2016; Sci Transl Med 2016; Gupta et al., Lancet Respir Med 2020 (PMID 31958400) | Multi-cohort discovery with held-out cohorts | leave-cohort-out / independent cohort validation (titles and snippets only) | not opened | n/a | S |
| P20 | Dobbin et al., J Immunother Cancer 2016 (PMID 27891226), "Validation of biomarkers ... Volume II" | Biomarker validation guidance | not seen | - | - | S (title only; page 403) |
| P21 | Preprint: "Widespread data leakage inflates accuracy and corrupts biomarker discovery in cancer drug response prediction" (bioRxiv 10.64898/2026.02.05.704016) | Leakage in biomarker pipelines | leaky pipelines gave 18.1 stable features vs 2.2 under correct pipeline (per search snippet) | not peer reviewed | - | S |
| P22 | Preprints 2026: LOL (bioRxiv 10.64898/2026.09.26.754596; Python leakage-aware GEO phenotype prediction; nested macro-F1 0.737 vs naive 0.791 per snippet); medRxiv 10.64898/2026.07.11.26357823 (leakage-free nested CV + LASSO, blood lung cancer); arXiv 2606.01566 RobustModelMaker (bootstrap stability selection + leakage-safe nested CV; CC BY 4.0, rest unreadable) | recent applications | - | not peer reviewed | - | S (arXiv: title and licence V) |

Note on SEQC2: I searched for SEQC2 model-building guidance and found only the SEQC/MAQC-III 2015 neuroblastoma paper (P4). I did not find a distinct SEQC2 guideline on classifier building. This is a gap, not evidence of absence.

## 4. Software (maintained, release within 3 years unless stated)

All "last release" values are from the CRAN/Bioconductor page opened (V). Licence from the same page.

| Tool | Aim / input / output | Bias control | Evidence of performance | Licence | Last release | URL | V/S |
|---|---|---|---|---|---|---|---|
| nestedcv 0.9.0 | Nested k x l CV for glmnet/caret on expression matrices; returns outer-fold performance and final model | Filters (t-test, Wilcoxon, RF, ReliefF, Boruta) inside outer loop; alpha CV; balancing only on outer training folds | P18 simulations | MIT + file (CRAN); MIT (GitHub) | 2026-07-14 | https://CRAN.R-project.org/package=nestedcv | V |
| cutpointr 1.2.1 | Cutpoint estimation (Youden etc.) with bootstrap and OOB metrics | in-bag vs out-of-bag bootstrap; bootstrapped cutpoints | P10 | GPL-3 | 2025-06-13 | https://CRAN.R-project.org/package=cutpointr | V |
| stabs 0.7-1 | Stability selection with error control (standard and complementary pairs) around lasso/boosting or user selector | PFER-controlled selection frequency | theory (S) | GPL-2 | 2026-01-31 | https://CRAN.R-project.org/package=stabs | V |
| sharp 1.4.8 | Stability selection with calibration of penalty and threshold by maximising stability vs null; also consensus clustering | calibrated selection proportion | not checked | GPL (>=3) | 2025-07-23 | https://CRAN.R-project.org/package=sharp | V |
| pminternal 0.1.0 | Internal validation of binary prediction models: bootstrap optimism, .632, CV optimism; c-statistic, Brier, calibration; bootstrap instability | optimism correction. Whether the user model-fitting function is re-run in full per replicate (including upstream gene filtering) is NOT confirmed from the page | not checked | GPL-3 | 2025-03-18 | https://CRAN.R-project.org/package=pminternal | V |
| glmnet 5.1 | Lasso/elastic net paths, cv.glmnet | built-in CV only for tuning; not a validation tool | - | GPL-2 | 2026-09-24 | https://CRAN.R-project.org/package=glmnet | V |
| rms 8.1-1 | Regression modelling strategies (Harrell); description lists "validation" | bootstrap validate() is the Harrell optimism tool, but the CRAN page does not state this | - | GPL-2 or GPL-3 | 2026-02-18 | https://CRAN.R-project.org/package=rms | V for version/licence; S for validate() |
| mlr3 1.8.0 | ML framework with resampling objects | nested resampling not mentioned on the page opened | - | LGPL-3 | 2026-08-21 | https://CRAN.R-project.org/package=mlr3 | V |
| caret 7.0-1 | Training/plotting of models | - | - | GPL-2 or GPL-3 | 2024-12-10 | https://CRAN.R-project.org/package=caret | V |
| rsample 1.3.2 | Resampling objects (bootstrap, CV) | nested_cv not named on page | - | MIT + file | 2026-01-30 | https://CRAN.R-project.org/package=rsample | V |
| pROC 1.19.1 | ROC curves, CIs for AUC/ROC | page does not mention coords/Youden | - | GPL (>=3) | 2026-08-27 | https://CRAN.R-project.org/package=pROC | V |
| ClassifyR 3.16.0 | Bioconductor framework for cross-validated classification, feature selection as a pipeline stage | page does not mention nested CV or selection-bias safeguards | - | GPL-3 | BioC 3.23, no date shown | https://bioconductor.org/packages/ClassifyR | V |
| DaMiRseq 2.24.0 | RNA-seq counts -> normalization, feature selection, stacked ensemble classifier | page lists batch/confounder removal; selection-bias handling not stated on page | - | GPL (>=2) | BioC 3.23, no date shown | https://bioconductor.org/packages/DaMiRseq | V |
| MLSeq 2.30.0 | ML interface for RNA-seq (SVM, bagSVM, RF, CART) | not stated on page | - | GPL (>=2) | BioC 3.23, no date shown | https://bioconductor.org/packages/MLSeq | V |
| mixOmics 6.36.0 | sPLS-DA etc. | perf()/CV not mentioned on the page | - | GPL (>=2) | BioC 3.23, no date shown | https://bioconductor.org/packages/mixOmics | V |
| rnaSentry (Bioconductor devel 3.24; GitHub Ghrieb/rnaSentry) | Prognostic (survival) RNA-seq signature discovery with audit ledger | reproducibility lock against re-selection; events-per-parameter guardrail; README admits screening-internal CV is optimistic and only validate_external() is out-of-sample | README: label-permuted CV C 0.744-0.785 overlapped the real signature (0.78); external LUAD C 0.540, "transfer not demonstrated"; illustrative, not peer reviewed | MIT | no release/tag shown; 97 commits, 2 stars | https://github.com/Ghrieb/rnaSentry | V. Survival-focused, so only partly relevant; worth reading for the "locked" and "null comparison" ideas |
| CMA (Bernau) weighted_mcr | Weighted-mean tuning-bias correction | smooth alternative to nested CV | - | not checked | not checked (GitHub chbernau/CMA) | https://rdrr.io/github/chbernau/CMA/man/weighted_mcr.html | S |
| Excluded | OptimalCutpoints 1.1-5 (2021-10-07, older than 3 years); MetaIntegrator (CRAN removed/archived 2024-10-14, page says undeliverable maintainer email); multiclassPairs (Bioconductor page not obtained) | | | | | | V/S |

## 5. How our pipeline compares

Facts about our pipeline come from your brief only; I have not read the code. Items marked "check" depend on code I did not see.

| Concern | Literature standard | Ours | Verdict |
|---|---|---|---|
| Selection bias from the gene-filter stage | Every data-driven selection step, including univariate filters, goes inside the resampling loop (P1, P2, P3, P18). nestedcv simulation: filtering on all data with pure noise gave inflated AUC, full nesting gave about 0.50 | DESeq2 candidate rules (padj, shrunk LFC, baseMean, two-route concordance, LOO >=90%) are run once on all samples; repeated nested CV covers only LASSO + cutoff | Main weakness. The nested CV AUC/accuracy is conditional on a candidate list that has already seen every outer test sample. It is not an unbiased estimate of the whole procedure. It is still useful for the LASSO+cutoff stage, but it must be labelled that way |
| Bootstrap optimism | Bootstrap must replicate every selection step (P7, S) | Harrell bootstrap present; check whether DESeq2 + rules are re-run in each replicate | If DESeq2 is not re-run, optimism is underestimated (per P7, S). Check |
| Cutoff | Data-driven cutoff is optimistic and unstable (P8, P9, P10); estimate it inside training data and evaluate on held-out data; report interval | Median of inner-fold Youden, bootstrap interval, evaluated in nested CV | Sound design. Check that outer-fold sensitivity/specificity use the cutoff derived only from that outer fold's training data. Consider OOB bootstrap check as in cutpointr (V) |
| Stability | Selection frequency alone is descriptive; stability selection adds a threshold with an error bound (P5, S) and instability should be assessed by repeating the whole build (P15) | Selection frequency reported; lambda.1se; no error-controlled threshold | Partly covered. LASSO does not remove small-sample instability (P15, V). Frequencies come only from the LASSO stage unless DESeq2 is also resampled |
| Locked model and new-sample processing | Freeze normalization parameters learned in training and apply them to test data (P16, V on microarrays); report threshold identification and full model (TRIPOD+AI items 15, 22) | Locked 9-gene model, pre-registered rule | Strong. Open point (inference, not verified): log2CPM depends on the library size and gene universe of each new sample; the locked specification should state the exact CPM/prior-count definition, and ideally a minimum depth/QC rule |
| External validation | At least one independent validation; do not reuse the validation set after changing protocol (P3) | Independent samples, pre-registered rule | Matches the guidance. Precision depends on the number of independent samples; report CIs |
| Multi-cohort | Leave-cohort-out or multiple independent cohorts (P19, S) | Not stated in brief | Unknown |
| Null benchmark | Permutation of labels through the full pipeline (P1 used permuted labels; rnaSentry README shows null CV can match the real signature) | Not stated | Recommend |

## 6. Ranked recommendations

1. Put the DESeq2 candidate stage inside the outer folds (or run a full-pipeline repeated CV and bootstrap) and report that number as the headline internal estimate. Supports: P1, P2, P3, P18. If LOO >=90% is too slow per fold, run a documented cheaper version inside folds and report both, stating that the cheaper version is an approximation.
2. Make the Harrell bootstrap re-run the complete pipeline in every replicate (DESeq2 rules + LASSO + cutoff). Support is S-level (Steyerberg via secondary source). Compare with the LASSO-only optimism to quantify how much the unnested filter was hiding.
3. Run a label-permutation null through the full pipeline (for example 100-200 permutations) and show where the observed nested AUC and the 9-gene result fall. This detects residual leakage cheaply.
4. Add a bootstrap instability analysis of the whole procedure (P15): how often the same 9 genes (or the same score classification of each sample) appear. Optionally use stabs (V) or sharp (V) for an error-controlled frequency threshold, and say clearly that lambda.1se plus frequency is not a formal error bound.
5. Check that cutoff evaluation is strictly out-of-fold, and add OOB bootstrap cutpoint validation (cutpointr, V) as a cross-check on the median-of-inner-fold-Youden rule. Report cutoff interval alongside sensitivity/specificity at the locked cutoff, not just AUC.
6. Write the locked-model specification so a third party can reproduce it: gene list, direction, exact log2CPM definition, prior count, normalization applied to new samples, cutoff, QC exclusion rule. Use TRIPOD+AI (V) as the checklist, in particular items 12c, 15, 22.
7. If more than one independent set exists, report each separately (leave-cohort-out style) rather than pooling (P19, S).
8. Consider cross-checking against nestedcv (V) on the same data as an independent implementation of nested CV; it supports glmnet with filters inside the loop, though its filters are univariate tests, not DESeq2.
9. Report precision of external validation (CI on AUC/sens/spec) given sample size; avoid claims beyond that interval.

## 7. What remains unverified

- Original papers not opened: Harrell, Lee & Mark 1996; Steyerberg 2003 (and the exact claim that selection must be repeated per bootstrap, which is from a search summary); Meinshausen & Buhlmann 2010 and the PFER formula; Altman 1994 (the 40% false-positive figure is from a secondary summary); Leeflang 2008; Tibshirani & Tibshirani 2009; Bernau et al.; Krstajic 2014; Faraggi & Simon 1996.
- SEQC/MAQC-III Genome Biol 2015 page failed; I did not locate SEQC2 model-building guidance or the 2014 Nat Biotechnol SEQC paper.
- Sweeney/Khatri multicohort papers and Gupta 2020: titles and snippets only. Dobbin 2016: title only.
- Preprints from 2026 (leakage, LOL, RobustModelMaker, medRxiv lung): snippets only, not peer reviewed; some pages returned 429.
- Whether pminternal, rms::validate and caret/mlr3 can re-run an upstream DESeq2 step per resample: not verified; the CRAN pages I opened do not say. I did not read vignettes or code.
- CRAN/Bioconductor candidate list came from my own prior knowledge plus web hits, not from a systematic registry search; tools I did not think of are missing. Bioconductor pages did not show release dates.
- Europe PMC and PubMed searches were noisy; 41 PubMed hits were never read; recall of the 2023-2026 peer-reviewed benchmark literature on RNA-seq signature validation is probably incomplete. I did not find a 2023-2026 peer-reviewed benchmark that compares full-pipeline vs partial nested CV for RNA-seq signatures.
- Claims about our pipeline rest on your brief; code was not inspected. The log2CPM freezing concern in section 5 is my inference.
- The V flag means the page content was returned by the fetch tool's summarizer; I did not read raw HTML and two long pages were read only partially.
