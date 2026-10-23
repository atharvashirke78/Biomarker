#!/usr/bin/env python3
# Read-only: dump every sheet of an .xlsx to TSV (stdlib only; no openpyxl on this pod).
# usage: python3 -I xlsx_to_tsv.py <in.xlsx> <out_dir>
import sys, os, re, zipfile, xml.etree.ElementTree as ET

NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
      'r': 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
M = '{%s}' % NS['m']


def col_index(ref):
    letters = re.match(r'([A-Z]+)', ref).group(1)
    n = 0
    for ch in letters:
        n = n * 26 + ord(ch) - 64
    return n - 1


def main(path, out):
    z = zipfile.ZipFile(path)
    names = z.namelist()
    ss = []
    if 'xl/sharedStrings.xml' in names:
        for si in ET.fromstring(z.read('xl/sharedStrings.xml')):
            ss.append(''.join(t.text or '' for t in si.iter(M + 't')))
    wb = ET.fromstring(z.read('xl/workbook.xml'))
    rels = {r.get('Id'): r.get('Target') for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
    os.makedirs(out, exist_ok=True)
    for s in wb.find('m:sheets', NS):
        name = s.get('name')
        tgt = rels[s.get('{%s}id' % NS['r'])].lstrip('/')
        tgt = tgt if tgt.startswith('xl/') else 'xl/' + tgt
        rows = []
        for row in ET.fromstring(z.read(tgt)).iter(M + 'row'):
            vals = {}
            for c in row.findall('m:c', NS):
                v = c.find('m:v', NS)
                isv = c.find('m:is', NS)
                if c.get('t') == 's' and v is not None:
                    x = ss[int(v.text)]
                elif isv is not None:
                    x = ''.join(t.text or '' for t in isv.iter(M + 't'))
                elif v is not None:
                    x = v.text
                else:
                    x = ''
                vals[col_index(c.get('r'))] = x.replace('\t', ' ').replace('\n', ' ')
            if vals:
                rows.append([vals.get(i, '') for i in range(max(vals) + 1)])
        fn = re.sub(r'[^A-Za-z0-9_.-]+', '_', name) + '.tsv'
        with open(os.path.join(out, fn), 'w') as f:
            for r in rows:
                f.write('\t'.join(r) + '\n')
        print(f'{name}\t{len(rows)} rows\t{fn}')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
