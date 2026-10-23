#!/usr/bin/env python3
# NEW (refactor Step 2, design section 2.5). Writes manifest_quant.json / manifest_analysis.json. stdlib only.
# usage: write_manifest.py <stage quant|analysis> <run_meta.json> <params.json> <out.json> [--file role=path ...] [--fastq path ...] [--refs ref_dir] [--base base_dir] [--quant_manifest path] [--geneset_md5 file]
# Tool versions are read by running the tools from the project environments (absolute paths); a tool that cannot be run is recorded as "NOT AVAILABLE: <reason>", never skipped.
import hashlib, json, os, subprocess, sys, time
args = sys.argv[1:]
stage, meta_f, params_f, out_f = args[:4]; rest = args[4:]
files, fastq, refs, base, envbase, qman, gsmd5 = [], [], "", "", "", "", ""
i = 0
while i < len(rest):
    k = rest[i]
    if k == "--file": files.append(rest[i + 1]); i += 2
    elif k == "--fastq":
        i += 1
        while i < len(rest) and not rest[i].startswith("--"): fastq.append(rest[i]); i += 1
    elif k == "--refs": refs = rest[i + 1]; i += 2
    elif k == "--base": base = rest[i + 1]; i += 2
    elif k == "--envbase": envbase = rest[i + 1]; i += 2     # where env/prefix/<tool>/bin lives (the install dir); default = --base
    elif k == "--quant_manifest": qman = rest[i + 1]; i += 2
    elif k == "--geneset_md5": gsmd5 = rest[i + 1]; i += 2
    else: sys.exit("unknown argument " + k)

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""): h.update(b)
    return h.hexdigest()
def finfo(p, do_md5=True):
    p = os.path.realpath(p)
    if not os.path.exists(p): return {"path": p, "error": "NOT FOUND"}
    d = {"path": p, "size": os.path.getsize(p)}
    if do_md5: d["md5"] = md5(p)
    return d
def run(cmd, env_bin=None):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        t = (r.stdout + r.stderr).strip().splitlines()
        return t[0] if t else "NOT AVAILABLE: no output (exit %d)" % r.returncode
    except Exception as e: return "NOT AVAILABLE: %s" % e

meta = json.load(open(meta_f)); params = json.load(open(params_f))
m = {"manifest_schema": 1, "stage": stage, "written": time.strftime("%Y-%m-%d %H:%M:%S"), "run": meta, "parameters": params}
try: m["git_commit"] = subprocess.run(["git", "-C", base or ".", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip() or None
except Exception: m["git_commit"] = None

# ---- tool versions (absolute paths into the project environments)
tv = {}
if base or envbase:
    P = lambda env, b: os.path.join(envbase or base, "env", "prefix", env, "bin", b)
    tv["fastp"] = run([P("fastp", "fastp"), "--version"]); tv["STAR"] = run([P("star", "STAR"), "--version"])
    tv["samtools"] = run([P("samtools", "samtools"), "--version"]); tv["salmon"] = run([P("salmon", "salmon"), "--version"])
    tv["featureCounts"] = run([P("subread", "featureCounts"), "-v"]); tv["multiqc"] = run([P("multiqc", "multiqc"), "--version"])
    tv["python(pyqc)"] = run([P("pyqc", "python"), "--version"])
    for env, pk in (("r_dge", ["DESeq2", "tximport", "apeglm", "ggplot2", "jsonlite"]), ("r_pathway", ["clusterProfiler", "fgsea", "msigdbr", "org.Hs.eg.db", "GO.db", "ReactomePA", "enrichplot"])):
        rs = P(env, "Rscript")
        tv["R(%s)" % env] = run([rs, "--version"])
        code = "for (p in c(%s)) cat(p, tryCatch(as.character(packageVersion(p)), error=function(e) 'NOT INSTALLED'), '\\n')" % ",".join('"%s"' % x for x in pk)
        try:
            r = subprocess.run([rs, "-e", code], capture_output=True, text=True, timeout=300)
            for l in r.stdout.strip().splitlines():
                n, v = l.split(None, 1) if " " in l else (l, "?"); tv["R:%s (%s)" % (n, env)] = v.strip()
        except Exception as e: tv["R packages (%s)" % env] = "NOT AVAILABLE: %s" % e
m["tool_versions"] = {"nextflow": meta.get("nextflow"), "java": meta.get("java"), **tv}

# ---- inputs / outputs
roles = {}
for f in files:
    role, p = f.split("=", 1); roles[role] = finfo(p)
m["files"] = roles
if fastq: m["inputs_fastq"] = [finfo(p, do_md5=meta.get("fastq_md5", True)) for p in fastq]
if refs:
    R = lambda n: os.path.join(refs, n); rel = params.get("gencode_release", "v50"); ref = {}
    for k, n in (("genome_fasta", "GRCh38_primary_chr.fa"), ("gtf", "gencode.%s.primary_chr.annotation.gtf.gz" % rel), ("transcripts_fasta", "gencode.%s.transcripts.primary_chr.fa.gz" % rel), ("tx2gene", "tx2gene.tsv"), ("gene_id_to_symbol", "gene_id_to_symbol.tsv")):
        ref[k] = finfo(R(n))
    si = R("salmon_index"); ref["salmon_index"] = {d: finfo(os.path.join(si, d)) for d in ("info.json", "versionInfo.json", "pos.bin", "seq.bin", "complete_ref_lens.bin") if os.path.exists(os.path.join(si, d))}
    st = R("star_index")   # DECISION (design item 9): parameters file md5 + size of the big files, not a whole-index md5 (slow)
    ref["star_index"] = {d: finfo(os.path.join(st, d), do_md5=(d == "genomeParameters.txt")) for d in ("genomeParameters.txt", "SA", "SAindex", "Genome", "chrNameLength.txt") if os.path.exists(os.path.join(st, d))}
    ref["star_index_md5_policy"] = "md5 of genomeParameters.txt + file sizes of SA/SAindex/Genome (whole-index md5 not computed)"
    m["references"] = ref
if gsmd5 and os.path.exists(gsmd5):
    m["gene_sets"] = [{"md5": l.split(None, 1)[0], "path": l.split(None, 1)[1].strip()} for l in open(gsmd5) if l.strip()]
if qman and os.path.exists(qman) and os.path.basename(qman) != "NO_FILE":
    m["quant_manifest"] = {"file": finfo(qman), "content": json.load(open(qman))}
json.dump(m, open(out_f, "w"), indent=2, default=str)
print("wrote", out_f)
