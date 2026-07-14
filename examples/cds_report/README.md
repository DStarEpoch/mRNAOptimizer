# cdsopt report example

This example shows how to evaluate existing CDS sequences with the `cdsopt report` command.

The `report` command scans every input sequence and computes:

- CAI, tAI, CG content
- MFE / avg_MFE, AUP
- CPB (codon pair bias)
- **Forbidden motif counts** (new): restriction sites, polyT runs, polyA signals, homopolymers

## Files

| File | Description |
|------|-------------|
| `sequences.fa` | CDS sequences from the GL004 design4 vaccine construct |

## Run

```bash
python -m cdsopt report sequences.fa -o report.csv
```

With fixed 5'/3' flanking sequences (for full-length folding):

```bash
python -m cdsopt report sequences.fa \
  --prefix five_utr.txt \
  --suffix three_utr_polya.txt \
  -o report_full.csv
```

## Output

`report.csv` contains one row per sequence with a `motif` column like:

```text
AATAAT:1|BsaI:1|BspQ1:1|FseI:1|NdeI:1|SacI:1|SapI:1
```

Each entry is `motif_name:count`. The separator is `|` so the field stays CSV-safe.

## Real-world reference

For a full design report, see:
`D:\workspace\projects\mRNA\GL004\cds_project\delivery\design4_report_with_motif.csv`
