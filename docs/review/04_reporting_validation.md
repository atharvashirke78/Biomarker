# Domain 4: Reporting guidelines, sample size, validation, regulatory (not legal advice)

Date of work: 2026-10-23. Marking: V = page/full text opened and read; S = search snippet or abstract only; M = from my own memory/calculation, not verified against a source.

## 1. Protocol (written before searching)
- Sources: EQUATOR Network; PubMed/PMC; journal sites (BMJ, Nature Biotechnology); FDA/CMS/EU public pages; MAQC/SEQC papers.
- Planned search strings (WebSearch unless noted; PubMed tool for PMIDs):
  1. "TRIPOD+AI statement updated guidance reporting clinical prediction models regression machine learning" (PubMed)
  2. "STARD 2015 an updated list of essential items for reporting diagnostic accuracy studies" (PubMed) and EQUATOR STARD page
  3. "REMARK reporting recommendations tumor marker prognostic studies McShane 2005 20-item checklist"
  4. "MAQC-II comprehensive study reliability predictive models Nature Biotechnology 2010"; "SEQC2 ... RNA-seq reproducibility"
  5. "Buderer 1996 sample size sensitivity specificity prevalence"; "Flahault 2005 sample size diagnostic test"; "Hajian-Tilaki 2014 sample size diagnostic test"; "Riley BMJ 2020 sample size clinical prediction model"; "Riley 2021 minimum sample size external validation binary"
  6. "Brown Cai DasGupta 2001 interval estimation binomial proportion"; "Schurch 2016 how many biological replicates RNA-seq"; "Conesa 2016 best practices RNA-seq"
  7. "FDA laboratory developed tests final rule vacated Eastern District of Texas 2025"; "FDA LDT rule rescinded Federal Register September 2025"; "laboratory developed tests 2026 update legislation"
  8. "EU IVDR 2017/746 in-house devices Article 5(5) Regulation 2024/1860"; "CMS CLIA 42 CFR 493.1253"
- Screening target: at least 40 records; log identified/screened/included/excluded below.

## 2. Screening log (approximate; I did not log counts programmatically)
- Identified: about 185 raw result lines across 20 WebSearch/PubMed queries; about 100 distinct URLs after removing duplicates and mirrors.
- Screened (title/domain): about 100 distinct. Excluded about 70: duplicate mirrors of the same paper (university repositories, Semantic-Scholar-style indexes), law-firm or vendor commentary repeating the same vacatur facts, wrong paper (see failed pages), commercial blogs, non-diagnostic topics.
- Included: 31 distinct sources cited here = 17 opened (V) + 14 snippet/abstract only (S).

Records screened and kept (named):
| # | Record | Status |
|---|---|---|
| 1 | TRIPOD+AI, BMJ 2024;385:e078378, PMC11019967 (PubMed full text) | V included |
| 2 | TRIPOD+AI erratum PMID 38636956 | S, not needed |
| 3 | TRIPOD+AI Korean translation PMID 40739973 | excluded (translation) |
| 4 | STARD 2015 E&E, Cohen 2016 BMJ Open, PMC5128957 (Table 2, 30 items) | V included |
| 5 | STARD 2015 BMJ h5527, PMC4623764 (abstract) | S included |
| 6 | EQUATOR STARD landing page | V included |
| 7 | EQUATOR STARD item 18 page | V included |
| 8 | REMARK E&E, PMC3362748 (Table 1, 20 items) | V included |
| 9 | REMARK 2005 BJC, PMC2361579 (text; table not rendered) | V included |
| 10 | MAQC-II, Nat Biotechnol 2010, PMC3315840 | V included |
| 11 | SEQC/MAQC-III RNA-seq, Nat Biotechnol 2014, PMC4321899 | V included |
| 12 | FDA LDT page (content current 09/19/2025) | V included |
| 13 | FDA RUO/IUO guidance page (Nov 2013) | V included (partial) |
| 14 | Arnold & Porter advisory, 9 Jul 2026 | V included |
| 15 | HPRA IVDR transitional provisions | V included |
| 16 | 42 CFR 493.1253 (Cornell LII) | V included |
| 17 | 21 CFR 809.10 (Cornell LII) | V included |
| 18 | epiR epi.ssdxsesp manual | V included |
| 19 | pmsampsize manual (rdrr.io) | V included |
| 20 | Foley, Venable, McDermott Hopkins, Mondaq, Exponent, GenomeWeb, Baker Donelson, Husch Blackwell, Hooper Lundy (LDT vacatur commentary) | S, one cited, rest excluded as duplicates |
| 21 | Federal Register public-inspection 2025-18239; AABB; FDA Law Blog (Sept 2025 reversion) | S included |
| 22 | AMP final judgment PDF | S, not opened |
| 23 | ASCP, CLP Mag, ARUP, Covington, CGM, ScreenMyGene (Enhancing CLIA Act 2026) | S included (status of H.R. 8890) |
| 24 | Hajian-Tilaki 2014 J Biomed Inform 48:193 | S (via epiR) |
| 25 | Buderer 1996 Acad Emerg Med 3:895 | S (via R package docs) |
| 26 | Flahault 2005 J Clin Epidemiol 58:859 | S (bibliographic only) |
| 27 | Riley et al. BMJ 2020;368:m441 | S (via pmsampsize) |
| 28 | Riley et al. Stat Med 2021;40:4230 | S |
| 29 | Brown, Cai, DasGupta 2001 Stat Sci 16:101 | S (abstract) |
| 30 | Schurch 2016 RNA 22:839 | S (abstract) |
| 31 | Conesa 2016 Genome Biol 17:13 | S (no replicate number retrieved) |
| 32-40 | Mirrors/duplicates of MAQC-II (Princeton, EPFL, UniTn, CNGB, ECNU, Cornell VIVO, Muni), Genohub blog, Frontiers plant paper, Thermo FAQ, Shimadzu whitepaper, HSE advice note, EAHP guidance, Danish Medicines Agency page, Springer Pathologie article | excluded: duplicate, vendor, or off-topic (Danish page cited as S for IVDR 5(5)(d)) |
| 41-45 | Bachmann/BMJ 2006 small-study survey; STARD 2015 PDF; GoodReports; ENCODE; Agresti-Coull 2001 | S or unreachable, not relied on |

### Failed or unusable pages
- bmj.com TRIPOD+AI article: HTTP 403 (obtained same text via PubMed full-text tool, PMC11019967).
- EQUATOR STARD-2015-checklist.pdf: binary PDF unreadable by fetcher (items obtained from the E&E Table 2 instead).
- pmc.ncbi.nlm.nih.gov via WebFetch: reCAPTCHA page (PMC2361579, PMC4321899); worked through the PubMed full-text tool.
- PMC3099432: a wrong paper (mitochondrial disorders); discarded. PMC7190074: wrong paper (insertion sequences); discarded.
- rnajournal.cshlp.org (Schurch): 403. fan.princeton.edu MAQC-II PDF: 403. projecteuclid (Brown 2001): empty. ecfr.gov: redirect to a block page (used Cornell LII instead).
- Federal Register 2025-18239 PDF: not opened (snippet only).
- Tables of the REMARK 2005 paper did not render; used the E&E Table 1 (same 20 items).

## 3. Important corrections to the question's premises
- SEQC2 is not an RNA-seq reproducibility project. SEQC2 (MAQC-IV) focuses on DNA sequencing (S, from search result summarising FDA MAQC/SEQC page). The relevant RNA-seq study is SEQC/MAQC-III (Nat Biotechnol 2014, V). Cite that, and cite SEQC2 only for DNA/NGS-assay-validation context.
- The only MAQC phase that addresses predictive models is MAQC-II (microarray, not RNA-seq). Its lessons transfer by analogy only.
- REMARK is written for prognostic markers, mainly single markers, with time-to-event outcomes (V). Its own text says it is not intended for models from very many candidate markers. Use it as a secondary checklist (items 1-11 on marker, assay, cutpoint, sample-size rationale), not the primary one.
- A tumour-vs-normal tissue classifier is not a prediction of a patient outcome. STARD fits better than TRIPOD+AI in spirit (index test vs histology reference standard), but a 9-gene score is a multivariable model, so TRIPOD+AI also applies. Recommended: report both, STARD as primary for the validation, TRIPOD+AI for model development/specification.

## 4. Reporting guidelines: what they require (V unless noted)

### 4.1 TRIPOD+AI (Collins et al., BMJ 2024;385:e078378; V)
URL: https://www.bmj.com/content/385/bmj-2023-078378 (403 to my fetcher); full text https://pmc.ncbi.nlm.nih.gov/articles/PMC11019967/ ; DOI 10.1136/bmj-2023-078378.
- 27 main items (52 sub-items), supersedes TRIPOD 2015 ("should no longer be used"). Items marked D = development, E = evaluation, D;E = both. A separate 13-item TRIPOD+AI for Abstracts exists (item 2). Authors are asked to submit the checklist with page/line numbers.
- Items (number, D/E, requirement, shortened) and our status. Status is inferred from the brief only (locked 9-gene score, 15 held-out samples, 7 normal/8 tumour, two studies, pre-registered at most 1 error per class, development performance labelled optimistic). "Unknown" means I cannot tell; confirm against the actual report.

| Item | D/E | Requirement (short) | Our status |
|---|---|---|---|
| 1 | D;E | Title identifies development/evaluation of a multivariable prediction model, target population, outcome | Unknown |
| 2 | D;E | Abstract per TRIPOD+AI for Abstracts | Unknown |
| 3a | D;E | Healthcare context (diagnostic/prognostic), rationale, existing models | Partial (need intended-use statement and comparison to existing colon signatures) |
| 3b | D;E | Target population, intended purpose in the care pathway, intended users | Missing (research use only; must state "no care-pathway use") |
| 3c | D;E | Known health inequalities between groups | Missing |
| 4 | D;E | Objectives; development vs validation | Present (pre-registered rule) |
| 5a | D;E | Data sources separately for development and evaluation, rationale, representativeness | Partial (two studies named; representativeness not assessed) |
| 5b | D;E | Dates of data collection | Missing/unknown (public datasets: give accession and deposit dates) |
| 6a | D;E | Setting, number and location of centres | Partial (two studies = at most two sites; give them) |
| 6b | D;E | Eligibility criteria | Missing/unknown |
| 6c | D;E | Treatments received and handling | Missing (pre-treatment status of tumours, neoadjuvant therapy) |
| 7 | D;E | Data pre-processing, quality checks, similar across sociodemographic groups | Partial (pipeline exists; QC thresholds and batch handling must be stated) |
| 8a | D;E | Outcome defined, time horizon, how assessed, consistent across groups | Partial (tumour vs normal by what reference standard: histology? database label?) |
| 8b | D;E | Qualifications of subjective outcome assessors | Missing (state if labels came from pathologists or from GEO metadata) |
| 8c | D;E | Blinding of outcome assessment | Missing |
| 9a | D | Choice of initial predictors, pre-selection | Partial (how the 9 genes were chosen from the candidate pool must be given) |
| 9b | D;E | Predictors defined, how/when measured, blinding | Partial |
| 9c | D;E | Assessor qualifications for subjective predictors | N/A (state: measured, not subjective) |
| 10 | D;E | How study size arrived at, separately for development and evaluation, any calculation | Missing (n=15 not justified; see section 6) |
| 11 | D;E | Missing-data handling | Unknown |
| 12a | D | How data were partitioned, considering sample size | Partial (held-out set exists; give partition logic and leakage controls) |
| 12b | D | Predictor handling (transformation, scaling) | Unknown |
| 12c | D | Model type, rationale, all model-building steps, hyperparameter tuning, internal validation method | Partial |
| 12d | D;E | Heterogeneity across clusters (studies/centres) | Partial (two studies: report per-study results) |
| 12e | D;E | All performance measures and plots, with rationale (discrimination, calibration, utility) | Partial (add calibration/score-distribution plot; threshold-based metrics only is incomplete) |
| 12f | E | Model updating after evaluation | Present if "locked" is true (state: no updating) |
| 12g | E | How predictions were calculated (formula, code) | Partial (publish weights, normalisation, threshold, code) |
| 13 | D;E | Class-imbalance methods | State "none used" |
| 14 | D;E | Fairness approaches | Missing (state sex/age/ancestry limits; public data may lack the metadata) |
| 15 | D | Model output; classification thresholds and how identified | Partial (state the threshold and that it was fixed before the held-out set) |
| 16 | D;E | Differences between development and evaluation data (setting, eligibility, outcome, predictors) | Missing (two-study heterogeneity table) |
| 17 | D;E | Ethics body and consent/waiver | Missing (public secondary data: state original approvals and that no new consent was needed) |
| 18a | D;E | Funding and role | Unknown |
| 18b | D;E | Conflicts of interest | Unknown |
| 18c | D;E | Where protocol accessible | Present if the pre-registration is public and timestamped |
| 18d | D;E | Registration details | Present/partial (give registry and number or say not registered) |
| 18e | D;E | Data availability | Unknown (accessions) |
| 18f | D;E | Code availability | Unknown |
| 19 | D;E | Patient and public involvement | Missing (state "none") |
| 20a | D;E | Participant flow, numbers with and without outcome | Partial (add flow diagram: raw samples, exclusions, split) |
| 20b | D;E | Characteristics per data source, key dates, predictors, sample size, events, missing data | Missing |
| 20c | E | Compare predictor distributions in evaluation vs development | Missing |
| 21 | D;E | Participants and events in each analysis | Partial |
| 22 | D | Full model specification to allow third-party predictions | Partial (locked score: publish coefficients, scaling, threshold) |
| 23a | D;E | Performance with confidence intervals, subgroups | Partial (CIs needed: exact intervals in 4.6) |
| 23b | D;E | Heterogeneity across clusters | Partial (per-study results) |
| 24 | E | Results of model updating | N/A if none |
| 25 | D;E | Interpretation including fairness | Partial |
| 26 | D;E | Limitations: non-representative sample, sample size, overfitting, missing data | Partial (dev-sample optimism is labelled, good; add n=15 imprecision) |
| 27a | D | How poor/unavailable input data should be handled at implementation | Missing (e.g., low RNA integrity, missing genes) |
| 27b | D | User interaction and expertise needed | Missing |
| 27c | D;E | Next steps, applicability, generalisability | Partial |

### 4.2 STARD 2015 (Bossuyt et al., BMJ 2015;351:h5527; item list from Cohen et al., BMJ Open 2016;6:e012799 Table 2; V)
URLs: https://www.equator-network.org/reporting-guidelines/stard/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC5128957/ ; PMC4623764. 30 items:
1 identify as diagnostic accuracy study (title/abstract); 2 structured abstract; 3 background, intended use and clinical role of index test; 4 objectives and hypotheses; 5 whether data collection was planned before index test and reference standard (prospective/retrospective); 6 eligibility; 7 how eligible participants identified; 8 where/when identified; 9 consecutive, random or convenience series; 10a index test in replicable detail; 10b reference standard in replicable detail; 11 rationale for reference standard; 12a/12b definition and rationale of positivity cut-offs (index/reference); 13a/13b blinding of index-test readers to reference/clinical info, and reverse; 14 methods for estimating/comparing accuracy; 15 indeterminate results; 16 missing data; 17 variability analyses (prespecified vs exploratory); 18 intended sample size and how determined; 19 flow diagram; 20 baseline characteristics; 21a severity distribution in target condition; 21b alternative diagnoses in non-diseased; 22 interval and interventions between tests; 23 cross-tabulation of index vs reference; 24 accuracy estimates with precision (95% CIs); 25 adverse events; 26 limitations (bias, uncertainty, generalisability); 27 implications, intended use and clinical role; 28 registration number; 29 where protocol can be accessed; 30 funding and role.
Our status: 3, 4, 12a (threshold locked), 29 present or partial; 5 (retrospective, public data: must say so), 9 (convenience series), 10a/10b, 13a/13b, 15, 16, 18, 19, 21a/21b (tumour stage, normal = adjacent vs healthy donor), 23 (2x2 table), 24 (CIs), 28, 30 missing or unknown. Item 18 wording (V): "Intended sample size and how it was determined"; page notes small studies give wide CIs and cites Flahault for methods.
Specific risk for us: normal samples that are tumour-adjacent normal vs healthy donor, and stage mix, change spectrum bias (STARD 21a/21b). State which.

### 4.3 REMARK (McShane et al. 2005; E&E Altman et al. 2012, BMC Med; V)
URLs: https://pmc.ncbi.nlm.nih.gov/articles/PMC2361579 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC3362748/ ; EQUATOR page https://www.equator-network.org/reporting-guidelines/reporting-recommendations-for-tumour-marker-prognostic-studies-remark/ . 20 items:
Introduction 1 marker, objectives, prespecified hypotheses. Methods: 2 patient characteristics, source, inclusion/exclusion; 3 treatments; 4 biological material and storage; 5 assay method, protocol, QC, reproducibility, blinding; 6 case selection, prospective/retrospective, time period; 7 endpoints; 8 candidate variables; 9 rationale for sample size, target power and effect size; 10 statistical methods, variable selection, missing data; 11 handling of marker values and cutpoint determination. Results: 12 flow of patients, numbers and events; 13 distributions of demographics and marker, missing values; 14 relation to standard prognostic variables; 15 univariable analyses; 16 multivariable effects with CIs; 17 marker plus standard variables regardless of significance; 18 further investigations incl. internal validation. Discussion: 19 interpretation vs prespecified hypotheses and limitations; 20 implications.
Relevance: items 1, 4, 5, 6, 8, 9, 10, 11, 12, 13, 18, 19 apply. Items 7, 14-17 are survival-oriented and mostly not applicable. The text says "not intended" for models from very many candidate markers (V, 2005 paper).

### 4.4 MAQC-II (Shi et al., Nat Biotechnol 2010;28:827; V) https://pmc.ncbi.nlm.nih.gov/articles/PMC3315840/
- 36 teams, 13 endpoints, more than 30,000 models; models frozen, then tested on a blinded validation set.
- Internal (cross-validation) performance was higher than external for most endpoints. Quote (V): "points to the need to subject every classifier to at least one external validation."
- Recommended steps: (1) design with "ratio of sample size to classifier complexity"; (2) pilot/internal validation by bootstrap or 10x repeated 5-fold CV, "feature selection and parameter optimization within each iteration"; (3) pivotal/external validation with complete separation of training and validation sets; do not reuse the validation set after modifying the protocol, because it then "turns the validation set into part of the design or training set."
- Endpoint predictability dominated performance (>65% of variance); negative-control (random label) and positive-control endpoints were included. Practical lesson: include a permutation/random-label negative control and a trivially predictable positive control (e.g., sex) in our pipeline.
- Authors say the guidelines "should not be construed as specific recommendations by the FDA."
- Applies to us: our held-out n=15 must be used once; any change to the 9-gene score or threshold after seeing the 15 invalidates them as a test set.

### 4.5 SEQC/MAQC-III RNA-seq (Su et al., Nat Biotechnol 2014;32:903; V) https://pmc.ncbi.nlm.nih.gov/articles/PMC4321899/
- Reference samples A (UHRR), B (brain), C, D (3:1 / 1:3 mixes) plus ERCC spike-ins; Illumina HiSeq 2000 sites with four library replicates per sample per site plus a vendor replicate.
- Findings: "low technical variance"; the authors "emphasize the value of biological replicates"; unfiltered differential expression gives high false positives; pipeline-dependent filters (p-value, fold-change, expression level) reduced empirical FDR; ERCC fold-change recovery needed about 90 million mapped fragments for log-concentration above 3; different analysis pipelines vary.
- Applies to us: the 9-gene score depends on normalisation and pipeline; report the pipeline version and show score stability across pipelines/sites.
- SEQC2 (DNA-focused) is not cited for RNA-seq claims (S).

### 4.6 Small-sample confidence intervals
- Brown, Cai, DasGupta 2001, Stat Sci 16:101 (S, abstract): recommend Wilson or equal-tailed Jeffreys intervals for small n, Agresti-Coull for larger n; Wald is erratic. Clopper-Pearson (exact) was not in that abstract; it is conservative (coverage at least nominal) and is the common regulatory default (M).
- Computed for our 15 samples (my calculation, formulas standard; see section 6):
  - Specificity 7/7: CP 0.590-1.000; Wilson 0.646-1.000. 6/7: CP 0.421-0.996; Wilson 0.487-0.974.
  - Sensitivity 8/8: CP 0.631-1.000; Wilson 0.676-1.000. 7/8: CP 0.473-0.997; Wilson 0.529-0.978.
  - Accuracy 15/15: CP 0.782-1.000; Wilson 0.796-1.000; 14/15: CP 0.681-0.998.
  - Meaning: even perfect results on 15 samples cannot exclude a true specificity near 0.6 or a sensitivity near 0.63. Report both CP and Wilson per class.
- Weakness of the pre-registered rule: with at most 1 error per class (7 normal, 8 tumour), probability of passing if true per-class accuracy is 0.95 is 0.90; 0.90 gives 0.69; 0.80 still gives 0.29; 0.70 gives 0.08 (my calculation, binomial). So passing is weak evidence; a failing result is more informative than a passing one.

## 5. Public regulatory pathways (not legal advice; US and EU only)
- FDA LDT final rule: published 6 May 2024 (89 Fed. Reg. 37286, amending 21 CFR 809.3(a)). Vacated and set aside in its entirety by E.D. Tex. on 31 March 2025 (V, FDA page; S, multiple law-firm summaries; ACLA/AMP plaintiffs). FDA did not appeal (V, Arnold & Porter 9 Jul 2026). On 19 September 2025 FDA issued a final rule reverting 809.3(a) to the pre-2024 text (V, FDA page "Content current as of 09/19/2025"; S, Federal Register public inspection 2025-18239). URL: https://www.fda.gov/medical-devices/in-vitro-diagnostics/laboratory-developed-tests . Discrepancy: the Arnold & Porter piece says the vacatur was in "May 2025"; FDA's page and the other sources say 31 March 2025; I use FDA's date.
- Current status as of this review: FDA is not enforcing the LDT rule; labs offering LDTs are governed by CLIA (CMS). A bill, H.R. 8890 "Enhancing Clinical Laboratory Innovation and Access Act of 2026" (introduced 19 May 2026), would put LDTs under an expanded CLIA framework; as of the 9 Jul 2026 advisory it had no action beyond introduction (V). A CMS Request for Information on CLIA regulations (CMS-3485) was reported as imminent (V, A&P). I could not confirm any event after about September 2026 or the current docket status of either item. Status after the FDA page date (19 Sep 2025) cannot be confirmed from FDA directly.
- CLIA analytical validation (V, 42 CFR 493.1253(b)(2), https://www.law.cornell.edu/cfr/text/42/493.1253): for modified or non-FDA-cleared tests a lab must, before reporting patient results, establish performance specifications: accuracy, precision (reproducibility), analytical sensitivity, analytical specificity (interfering substances), reportable range, reference intervals, and other needed characteristics. That is analytical, not clinical, validation.
- RUO/IUO labelling (V, 21 CFR 809.10(c)(2), https://www.law.cornell.edu/cfr/text/21/809.10 ; FDA guidance Nov 2013 page): RUO statement: "For Research Use Only. Not for use in diagnostic procedures." (laboratory research phase, not represented as effective IVD). IUO statement: "For Investigational Use Only. The performance characteristics of this product have not been established." The FDA RUO/IUO guidance page I opened did not give label text itself (the regulation does). Practical implication: a published research score is not a diagnostic; any patient-result use needs a CLIA-certified lab (and FDA clearance if sold as a kit/software, S via A&P statement that FDA still asserts jurisdiction over distributed kits, components and software).
- FDA premarket routes (510(k), De Novo, PMA): not opened; unverified here.
- EU IVDR 2017/746 (V, HPRA https://www.hpra.ie/regulation/medical-devices/manufacturers-and-authorised-representatives/ivdr-transitional-provisions ; date not shown on page): Regulation (EU) 2024/1860 in force 9 Jul 2024; legacy-device deadlines Class D 31 Dec 2027, C 31 Dec 2028, B and A-sterile 31 Dec 2029. In-house devices (health institutions, Art. 5(5)): Annex I applies from 26 May 2022; Art. 5(5)(b),(c),(e)-(i) from 26 May 2024; 5(5)(d) (justification that needs cannot be met by an equivalent marketed device) from 31 Dec 2030. A new risk-classified IVD (genomic cancer diagnostic likely class C, M) needs notified-body conformity assessment; research-use-only products are outside IVDR scope (M, not verified).

## 6. Sample size methods, formulas, numbers
Notation: z = 1.96 (95%), Se/Sp expected, d = absolute half-width, prev = tumour fraction.

### 6.1 Buderer / Hajian-Tilaki (precision of Se and Sp)
Formula (M, checked numerically; the epiR manual, V, states it follows Hajian-Tilaki 2014 and gives a worked example of 1255 subjects for Se 0.80, prevalence 0.10, absolute error 0.07; my formula gives 1254.35, rounds to 1255, so the formula is confirmed for the Se case):
- n_cases = z^2 Se(1-Se) / d^2 ; total N_Se = n_cases / prev
- n_controls = z^2 Sp(1-Sp) / d^2 ; total N_Sp = n_controls / (1-prev); N = max(N_Se, N_Sp).
Results per class (cases = tumour, controls = normal), 95% confidence:
| Expected Se=Sp | d | n per class | total at prev 0.5 | total at prev 0.4 |
|---|---|---|---|---|
| 0.95 | 0.05 | 73 | 146 | 183 |
| 0.95 | 0.10 | 19 | 38 | 46 |
| 0.90 | 0.05 | 139 | 278 | 346 |
| 0.90 | 0.10 | 35 | 70 | 87 |
| 0.85 | 0.10 | 49 | 98 | 123 |
| 0.80 | 0.10 | 62 | 124 | 154 |
Our 15 samples is below even the loosest row (19 per class). Flahault et al. 2005 (J Clin Epidemiol 58:859) is cited by STARD item 18 as a method source; I could not open it, so I do not state its formula (unverified).

### 6.2 Exact (Clopper-Pearson) zero/one-error designs (my calculation)
- With 0 errors in n per class, CP lower 95% (two-sided) bound = (0.025)^(1/n). To exclude accuracy below L: n >= ln(0.025)/ln(L).
  - L=0.80: n=17; L=0.85: 23; L=0.90: 36; L=0.95: 72 per class.
- Allowing 1 error in n per class and CP lower bound >= L: L=0.80: n=26; L=0.85: 35; L=0.90: 54.
- Hypothesis-test style (one-sided exact binomial, alpha <= 0.05, per class, pass if errors <= c):
  - H0 accuracy 0.80 vs H1 0.95: n=30 per class, c=2; alpha 0.044, power 0.81.
  - H0 0.80 vs H1 0.90: n=82, c=10; alpha 0.046, power 0.81.
  - H0 0.90 vs H1 0.99: n=46, c=1; alpha 0.048, power 0.92.
  - H0 0.85 vs H1 0.97: n=40, c=2; alpha 0.049, power 0.88.

### 6.3 AUC (Hajian-Tilaki)
Formula (M; the Hanley-McNeil-based variance as used by Hajian-Tilaki; not verified against the paper): a = sqrt(2) * Phi^-1(AUC); V = 0.0099 e^(-a^2/2) [(5a^2+8) + (a^2+8)/R]; n per group = z^2 V / d^2 for CI half-width d, R = n_normal/n_tumour.
- AUC 0.95: d=0.05 -> 50 per group; d=0.03 -> 137.
- AUC 0.90: d=0.05 -> 106. AUC 0.98: d=0.05 -> 15; d=0.03 -> 42.
- Current set (8 tumour/7 normal): Hanley-McNeil SE 0.060 at AUC 0.95 (half-width 0.118); 0.038 at AUC 0.98 (half-width 0.074).
Caution: AUC near 1 has unstable normal-approximation intervals; use bootstrap/exact methods and report the Mann-Whitney AUC with a CI, not an AUC of 1.000 alone.

### 6.4 Riley et al. for prediction models (S for method; formulas M)
- Development (BMJ 2020;368:m441, via pmsampsize manual V): binary outcome, three criteria: (i) expected shrinkage <= 10% (default 0.9), (ii) apparent vs adjusted Nagelkerke R^2 differ by <= 0.05, (iii) average risk estimated within +-0.05. Criterion (i): n = p / ((S-1) ln(1 - R2cs/S)).
  - Our 9-gene score is locked, so this is for refitting only. With p=9: R2cs 0.2 -> 359; 0.3 -> 222; 0.4 -> 154; 0.5 -> 111. Maximum R2cs at prevalence 0.5 is 0.75 (M). Criterion (iii) at prevalence 0.5, +-0.05: n=385; +-0.10: 97. If the 9 genes were chosen from thousands, effective p is much more than 9 and n must be larger; treat these numbers as lower bounds.
- External validation (Riley 2021, Stat Med 40:4230; S, abstract): sample size for precise calibration (O/E, calibration slope), C-statistic and net benefit; the abstract's example needed 9,835 participants (177 events) at event proportion 0.018, driven by calibration slope. Rule of thumb of at least 100 events and 100 non-events is commonly quoted (M; not verified here). For O/E with prevalence 0.5 and target SE(ln O/E)=0.10: n about 100 (approximation Var = (1-phi)/(n phi), M); SE 0.15: n about 45.
- Practical reading: a calibrated probability claim needs on the order of 100 per class; a threshold-based classification claim can be supported with 35-73 per class (6.1, 6.2).

### 6.5 RNA-seq replicate recommendations
- Schurch et al. 2016 (RNA 22:839; S, abstract): 48 biological replicates per condition in yeast; at least 6 replicates per condition, at least 12 when all fold changes must be detected; with 3 replicates most tools recovered only 20-40% of significant genes. This addresses differential-expression discovery, not a locked classifier.
- SEQC 2014 (V): four library replicates per sample per site; technical variance is low relative to biological variance; ERCC fold-change recovery at about 90 million mapped fragments.
- Conesa 2016 (S): survey of RNA-seq best practices; I did not retrieve its replicate number.
- For our locked score the relevant replicate question is score stability: technical (library) and batch/site replicates, not DE power. See validation plan stage C.

## 7. Gap table
| Item | Required by | Our status | Fix |
|---|---|---|---|
| Intended use and clinical role of the test | TRIPOD+AI 3a/3b; STARD 3, 27 | Missing | State "research use, no care-pathway use", name any future decision |
| Sample-size justification | TRIPOD+AI 10; STARD 18; REMARK 9 | Missing (n=15) | Add section 6 calculation; declare n=15 as exploratory/pilot |
| CIs on performance | TRIPOD+AI 23a; STARD 24 | Partial | Report CP and Wilson per class (4.6 numbers) |
| Reference standard definition | STARD 10b, 11; TRIPOD+AI 8a-c | Unknown | State label source (histology vs metadata), who assigned, blinding |
| Spectrum (stage, normal type) | STARD 21a/21b | Missing | Table of stage and normal source per study |
| Retrospective/convenience series | STARD 5, 9 | Missing | Declare both |
| Blinding of index test readers | STARD 13a/13b; MAQC-II | Unknown | State code ran on labels-hidden data; log timestamp |
| Threshold fixed before testing | STARD 12a; TRIPOD+AI 15 | Present (pre-registered) | Publish the registration hash/date |
| Flow diagram and per-study data | STARD 19, 20; TRIPOD+AI 20a/20b/16 | Missing | Add diagram and per-study table |
| Per-study heterogeneity | TRIPOD+AI 12d, 23b | Partial | Report each study separately |
| Full model specification and code | TRIPOD+AI 12g, 22, 18f | Partial | Publish gene IDs, coefficients, normalisation, threshold, code |
| Candidate-gene selection history | TRIPOD+AI 9a; REMARK 8, 10, 11 | Partial | Document selection from full gene set and the data used |
| Dev-set optimism | TRIPOD+AI 26; MAQC-II | Present (labelled optimistic) | Add repeated CV with selection inside the loop (MAQC-II) |
| Single-use of held-out set | MAQC-II | Unknown | Log that the 15 were used once; any re-lock requires new data |
| Negative/positive control endpoints | MAQC-II | Missing | Add permuted-label and sex-prediction controls |
| Calibration/score distribution | TRIPOD+AI 12e | Missing | Plot score by class and study |
| Pipeline/version reporting | TRIPOD+AI 7, 12b; REMARK 5; SEQC | Partial | Pin tool versions, reference annotation, QC thresholds |
| Poor-input handling | TRIPOD+AI 27a | Missing | Define RNA-quality and missing-gene rules |
| Fairness/ethics/funding/PPI/COI | TRIPOD+AI 3c, 14, 17, 18, 19 | Missing/unknown | One-line statements, including "not assessed" |
| Registration and protocol access | TRIPOD+AI 18c/d; STARD 28/29 | Partial | Public timestamped registry entry |
| Regulatory label | 21 CFR 809.10(c)(2) | Not applicable now | Add RUO statement to any distributed software/outputs |
| Analytical validation (accuracy, precision, specificity, range) | 42 CFR 493.1253(b)(2) | Missing | Only needed before patient reporting in a CLIA lab; plan stage C |

## 8. Validation plan with computed sample sizes
Stage A: lock and document (no data needed). Freeze genes, coefficients, normalisation, pipeline version, threshold, and the acceptance rule in a timestamped registration. Do not refit after any evaluation.
Stage B: independent external validation (new studies, not the 15 used).
1. Primary precision design: choose Se=Sp expected 0.90, d=0.10 -> 35 tumour + 35 normal (Buderer, 6.1). More robust: 0.95 expected, d=0.05 -> 73 + 73.
2. Primary decision rule (exact test): H0 accuracy 0.80 vs H1 0.95, per class: n=30, pass if at most 2 errors (alpha 0.044, power 0.81). Use for both Se and Sp separately; to claim both, require each to pass (family alpha about 0.09 unless one-sided 0.025 each; use alpha 0.025 per class, which needs a slightly larger n; recompute before registering).
3. Stricter claim: H0 0.90 vs H1 0.99: n=46 per class, pass if at most 1 error.
4. Recommended target: at least 50 per class from at least three independent studies/sites (covers 6.1 row 0.95/0.05 partly, AUC half-width 0.05 at AUC 0.95, and per-study heterogeneity). Report per-study results; no pooling that hides a failing study.
5. AUC: 50 per group gives CI half-width 0.05 at AUC 0.95 (6.3).
6. Calibration claims (if probabilities are reported): at least 100 per class (6.4); otherwise report the score and the binary call only.
7. Include controls: permuted-label negative control (expected AUC about 0.5) and a trivially predictable positive control.
8. Include spectrum: early and late stage, normal-adjacent and (if available) healthy-donor normal, both sexes, more than one age group; report each.
Stage C: technical robustness (needed for any move toward clinical use). Library replicates of at least 3 tumour and 3 normal samples in at least 2 batches/sites (SEQC design used 4 library replicates per site); report score CV, class-call agreement, pipeline-version invariance, input-quality limits. If clinical use is intended later, this stage maps to CLIA 493.1253(b)(2) performance specifications and to IVDR/FDA analytical requirements; plan a clinical-grade assay (RT-qPCR panel is typical, M) rather than whole RNA-seq.
Stage D (only if clinical use is pursued): prospective, consecutive, pre-specified intended-use cohort with STARD/TRIPOD+AI reporting, a regulatory pre-submission discussion, and appropriate approval/consent. Not required for research use.
Reporting for the current 15-sample result: describe as exploratory pilot, give CP and Wilson CIs, state the pass probability caveat from 4.6 (rule passes 29% of the time at true accuracy 0.80), state in-sample development performance as optimistic and the held-out result as one use of a small set.

## 9. What remains unverified
- Flahault 2005 formula and content (S, bibliographic only).
- Hajian-Tilaki 2014 AUC formula (M) and Buderer's original formula (checked only against the epiR 1255 example for the sensitivity case).
- Riley 2020/2021 formulas and the 100 events/100 non-events rule (M or S).
- Schurch 2016 abstract details (S); Conesa 2016 replicate guidance (not retrieved); Brown-Cai on Clopper-Pearson (not retrieved).
- STARD 2015 original checklist PDF (items taken from the E&E Table 2, V).
- REMARK 2005 checklist table (not rendered; used E&E Table 1, V).
- Current legal status of the FDA LDT rule after 19 Sep 2025 (FDA page date) and of H.R. 8890 or CMS-3485 after about September 2026; the court vacatur date differs in one source (May vs 31 Mar 2025); I cannot confirm there has been no later FDA action.
- FDA premarket pathways (510(k)/De Novo/PMA) and IVDR classification of this kind of test (M, not opened).
- SEQC2 scope (S); the claim that RUO products are outside IVDR scope (M).
- Our actual report contents: TRIPOD+AI and STARD statuses are inferred from the user's brief, not from reading the report.
