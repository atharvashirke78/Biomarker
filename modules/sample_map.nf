process MAKE_SAMPLE_MAP {
    label 'process_low'
    executor 'local'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path samplesheet

    output:
    path "sample_map.tsv", emit: map

    script:
    """
    python3 - "${samplesheet}" <<'PY'
import csv
import sys

required = ["sample", "run", "study", "group"]
rows = []
seen_samples = set()
seen_runs = set()

with open(sys.argv[1], newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle)

    if not reader.fieldnames:
        raise SystemExit("ERROR: Empty sample sheet.")

    reader.fieldnames = [
        name.strip() if name else name
        for name in reader.fieldnames
    ]

    if len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise SystemExit("ERROR: Duplicate column names.")

    missing = [name for name in required if name not in reader.fieldnames]
    if missing:
        raise SystemExit("ERROR: Missing columns: " + ", ".join(missing))

    for line, row in enumerate(reader, start=2):
        if None in row:
            raise SystemExit(f"ERROR: Too many fields on line {line}.")

        values = {
            name: (row.get(name) or "").strip()
            for name in required
        }

        if not any(values.values()):
            continue

        empty = [name for name, value in values.items() if not value]
        if empty:
            raise SystemExit(
                f"ERROR: Empty fields on line {line}: "
                + ", ".join(empty)
            )

        if values["sample"] in seen_samples:
            raise SystemExit(f"ERROR: Duplicate sample on line {line}.")
        if values["run"] in seen_runs:
            raise SystemExit(f"ERROR: Duplicate run on line {line}.")

        seen_samples.add(values["sample"])
        seen_runs.add(values["run"])

        rows.append([
            values["run"],
            values["sample"],
            values["study"],
            values["group"]
        ])

if not rows:
    raise SystemExit("ERROR: No sample rows found.")

with open("sample_map.tsv", "w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out, delimiter="\\t", lineterminator="\\n")
    writer.writerow(["srr", "label", "study", "group"])
    writer.writerows(rows)

print(f"Created sample_map.tsv with {len(rows)} rows.")
PY
    """

    stub:
    """
    python3 - "${samplesheet}" <<'PY'
import csv
import sys

required = ["sample", "run", "study", "group"]
rows = []

with open(sys.argv[1], newline="", encoding="utf-8-sig") as handle:
    reader = csv.DictReader(handle)

    if not reader.fieldnames:
        raise SystemExit("ERROR: Empty sample sheet.")

    reader.fieldnames = [
        name.strip() if name else name
        for name in reader.fieldnames
    ]

    if len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise SystemExit("ERROR: Duplicate column names.")

    missing = [name for name in required if name not in reader.fieldnames]
    if missing:
        raise SystemExit("ERROR: Missing columns: " + ", ".join(missing))

    for line, row in enumerate(reader, start=2):
        if None in row:
            raise SystemExit(f"ERROR: Too many fields on line {line}.")

        values = {
            name: (row.get(name) or "").strip()
            for name in required
        }

        if not any(values.values()):
            continue

        empty = [name for name, value in values.items() if not value]
        if empty:
            raise SystemExit(
                f"ERROR: Empty fields on line {line}: "
                + ", ".join(empty)
            )

        if values["sample"] in {r[1] for r in rows}:
            raise SystemExit(f"ERROR:Duplicate sample on line {line}.")
        if values["run"] in {r[0] for r in rows}:
            raise SystemExit(f"ERROR: Duplicate run on line {line}.")

        rows.append([
            values["run"],
            values["sample"],
            values["study"],
            values["group"]
        ])

if not rows:
    raise SystemExit("ERROR: No sample rows found.")

with open("sample_map.tsv", "w", newline="", encoding="utf-8") as out:
    writer = csv.writer(out, delimiter="\\t", lineterminator="\\n")
    writer.writerow(["srr", "label", "study", "group"])
    writer.writerows(rows)
PY
    """
}
