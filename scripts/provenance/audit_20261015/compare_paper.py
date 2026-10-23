#!/usr/bin/env python3
# Read-only comparison: Hosseini 2023 supplement (Table S3 DEG lists, S7 40 named DEGs, S8 genes) vs the NEW 19-sample fit (2026-10-14).
# Run from the project root: env/prefix/pyqc/bin/python3 -I results_rc1prep_19/audit_compare_20261015/compare_paper.py
# Symbol-level comparison (S3 gives symbols only, no fold change, no Ensembl IDs): one row per symbol, the tested row with the highest baseMean.
import os
import pandas as pd

O = 'results_rc1prep_19/audit_compare_20261015'
D = os.path.join(O, 'xlsx_dump/paper_MOESM2')
R = 'results_rc1prep_19'
log = []


def say(*a):
    m = ' '.join(str(x) for x in a)
    print(m)
    log.append(m)


# ---- paper lists (S3: row 1 title, row 2 Up/Down header, row 3 column names, data from row 4)
s3 = pd.read_csv(os.path.join(D, 'Table_S3.tsv'), sep='\t', header=None, skiprows=3, dtype=str, keep_default_na=False)
up_p = [g.strip() for g in s3[0] if g.strip()]
dn_p = [g.strip() for g in s3[4] if g.strip()]
say('paper S3 rows: up', len(up_p), 'down', len(dn_p), '| unique symbols: up', len(set(up_p)), 'down', len(set(dn_p)),
    '| in both lists:', len(set(up_p) & set(dn_p)))
pv = pd.to_numeric(pd.concat([s3[2], s3[6]]), errors='coerce').dropna()
ms = pd.to_numeric(pd.concat([s3[1], s3[5]]), errors='coerce').dropna()
say('paper S3 value range: pvalue max', pv.max(), '| metaSeq min', ms.min(), '(thresholds NOT stated in the supplement; range only)')
paper = pd.DataFrame({'symbol': list(set(up_p)) + list(set(dn_p)),
                      'paper_dir': ['up'] * len(set(up_p)) + ['down'] * len(set(dn_p))})
paper = paper[~paper.symbol.isin(set(up_p) & set(dn_p))]

fin = pd.read_csv(os.path.join(R, 'final_candidates/final_candidates_18of19.tsv'), sep='\t')
cand_sym = set(fin.symbol.dropna())
model = pd.read_csv('lock/diagnostic_v1/final_model_genes.tsv', sep='\t')
model_sym = set(model.symbol)


def ours(route_file):
    d = pd.read_csv(os.path.join(R, route_file), sep='\t')
    d = d[d.gene_name.notna() & (d.gene_name != '')]
    d = d.sort_values('baseMean', ascending=False).drop_duplicates('gene_name')
    sig_un = (d.padj < 0.05) & (d.log2FoldChange.abs() >= 1)
    sig_r1 = (d.padj < 0.05) & (d.log2FoldChange_shrunk.abs() >= 1)
    d = d.assign(dir_un=pd.Series(pd.NA, index=d.index, dtype='object'), dir_r1=pd.Series(pd.NA, index=d.index, dtype='object'))
    d.loc[sig_un, 'dir_un'] = d.loc[sig_un, 'log2FoldChange'].map(lambda x: 'up' if x > 0 else 'down')
    d.loc[sig_r1, 'dir_r1'] = d.loc[sig_r1, 'log2FoldChange_shrunk'].map(lambda x: 'up' if x > 0 else 'down')
    return d.set_index('gene_name')


rows = []
for route, f in [('featureCounts', 'deseq2_results_all.tsv'), ('Salmon', 'deseq2_salmon_results_all.tsv')]:
    d = ours(f)
    for rule, col in [('padj<0.05 & |unshrunk LFC|>=1 (pipeline deseq2.R)', 'dir_un'), ('rule 1: padj<0.05 & |shrunk LFC|>=1', 'dir_r1')]:
        for pdir in ['up', 'down']:
            P = paper[paper.paper_dir == pdir].symbol
            ours_dir = d[col]
            n_ours = int((ours_dir == pdir).sum())
            in_tab = P[P.isin(d.index)]
            matched = int((ours_dir.reindex(in_tab) == pdir).sum())
            opposite = int((ours_dir.reindex(in_tab) == ('down' if pdir == 'up' else 'up')).sum())
            rows.append(dict(route=route, rule=rule, direction=pdir, paper_n=len(P), ours_n=n_ours,
                             paper_genes_in_our_tested_table=len(in_tab), paper_genes_not_in_our_table=len(P) - len(in_tab),
                             overlap_same_direction=matched, overlap_opposite_direction=opposite,
                             recall_pct=round(100 * matched / len(P), 1), ours_also_in_paper_same_dir_pct=round(100 * matched / n_ours, 1) if n_ours else None,
                             direction_agreement_on_overlap_pct=round(100 * matched / (matched + opposite), 1) if matched + opposite else None))
summ = pd.DataFrame(rows)
summ.to_csv(os.path.join(O, 'paper_vs_ours_deg_counts.tsv'), sep='\t', index=False)
say(summ.to_string(index=False))

# ---- named genes: S7 (40 DEGs) + S8 (GUCA2A, COL3A1)
s7 = pd.read_csv(os.path.join(D, 'Table_S7_Co-Expression_matrix_.tsv'), sep='\t', header=None, skiprows=2, dtype=str, keep_default_na=False)
named = [g for g in s7[0] if g.strip()]
for g in ['GUCA2A', 'COL3A1']:
    if g not in named:
        named.append(g)
fc = ours('deseq2_results_all.tsv')
sal = ours('deseq2_salmon_results_all.tsv')
pdir = dict(zip(paper.symbol, paper.paper_dir))
out = []
for g in named:
    r = fc.loc[g] if g in fc.index else None
    s = sal.loc[g] if g in sal.index else None
    out.append(dict(symbol=g, source='S7' + ('+S8' if g in ('GUCA2A', 'COL3A1') else ''), paper_dir_S3=pdir.get(g, 'not in S3'),
                    in_our_fc_table=r is not None,
                    fc_padj=None if r is None else r.padj, fc_log2FC_unshrunk=None if r is None else r.log2FoldChange,
                    fc_log2FC_shrunk=None if r is None else r.log2FoldChange_shrunk,
                    fc_rule1_dir=None if r is None else r.dir_r1, salmon_rule1_dir=None if s is None else s.dir_r1,
                    in_1390_candidates=g in cand_sym, in_9gene_model=g in model_sym))
nm = pd.DataFrame(out)
nm.to_csv(os.path.join(O, 'paper_named_genes_status.tsv'), sep='\t', index=False)
say('named genes:', len(nm), '| in our fc table:', int(nm.in_our_fc_table.sum()),
    '| rule 1 fc same direction as S3:', int((nm.fc_rule1_dir == nm.paper_dir_S3).sum()),
    '| rule 1 fc opposite:', int(((nm.fc_rule1_dir.notna()) & (nm.paper_dir_S3.isin(['up', 'down'])) & (nm.fc_rule1_dir != nm.paper_dir_S3)).sum()),
    '| in 1,390:', int(nm.in_1390_candidates.sum()), sorted(nm.symbol[nm.in_1390_candidates]),
    '| in 9-gene model:', int(nm.in_9gene_model.sum()), sorted(nm.symbol[nm.in_9gene_model]))
# 9-gene model genes in the paper lists
m9 = pd.DataFrame({'symbol': sorted(model_sym)})
m9['paper_dir_S3'] = m9.symbol.map(lambda g: pdir.get(g, 'not in S3'))
m9['model_set'] = m9.symbol.map(dict(zip(model.symbol, model.set)))
m9.to_csv(os.path.join(O, 'model9_genes_in_paper.tsv'), sep='\t', index=False)
say(m9.to_string(index=False))
c = fin[['symbol', 'direction_fc']].dropna().drop_duplicates('symbol')
c['paper'] = c.symbol.map(lambda g: pdir.get(g, 'not in S3'))
say('1,390 candidates by paper S3 status (symbol level):', c.groupby(['direction_fc', 'paper']).size().to_dict())
open(os.path.join(O, 'compare_paper.log'), 'w').write('\n'.join(log) + '\n')
