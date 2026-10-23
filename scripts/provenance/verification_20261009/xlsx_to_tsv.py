# stdlib-only xlsx -> one TSV per sheet (cell grid, no interpretation). usage: xlsx_to_tsv.py in.xlsx outdir
import sys, zipfile, re, os, xml.etree.ElementTree as ET
ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main", "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
z = zipfile.ZipFile(sys.argv[1]); out = sys.argv[2]; os.makedirs(out, exist_ok=True)
ss = [("".join(t.text or "" for t in si.iter("{%s}t" % ns["m"]))) for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", ns)]
wb = ET.fromstring(z.read("xl/workbook.xml")); rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
def col(ref):
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(): n = n * 26 + ord(ch) - 64
    return n - 1
for sh in wb.find("m:sheets", ns):
    name = sh.get("name"); tgt = rels[sh.get("{%s}id" % ns["r"])]; tgt = tgt.lstrip("/"); tgt = tgt if tgt.startswith("xl/") else "xl/" + tgt
    rows = []
    for row in ET.fromstring(z.read(tgt)).iter("{%s}row" % ns["m"]):
        r = []
        for c in row.findall("m:c", ns):
            i = col(c.get("r")); v = c.find("m:v", ns); t = c.get("t")
            val = "" if v is None else (ss[int(v.text)] if t == "s" else v.text)
            if t == "inlineStr": val = "".join(x.text or "" for x in c.iter("{%s}t" % ns["m"]))
            r += [""] * (i + 1 - len(r)); r[i] = val.replace("\t", " ").replace("\n", " ")
        rows.append(r)
    with open(os.path.join(out, re.sub(r"\W+", "_", name) + ".tsv"), "w") as f:
        for r in rows: f.write("\t".join(r) + "\n")
    print(name, len(rows), "rows")
