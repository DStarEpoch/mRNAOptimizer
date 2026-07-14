# cdsopt optimize example

This example shows how to run CDS optimization with the `cdsopt optimize` command.

## Files

| File | Description |
|------|-------------|
| `protein.fa` | Amino-acid sequence of the target protein (46 aa) |
| `init_cds.fa` | Optional initial CDS used to seed the population |
| `config.yaml` | Optimization parameters and forbidden-motif filters |

## Run

From this directory:

```bash
# Bash / Git Bash
python -m cdsopt optimize protein.fa -c config.yaml

# Windows CMD
python -m cdsopt optimize protein.fa -c config.yaml
```

Or pass options directly on the command line:

```bash
python -m cdsopt optimize protein.fa \
  --pop-size 50 --generations 30 \
  --target-cai 0.92 --target-cg 0.55 \
  --output-dir ./outputs
```

## Output

After running, the `./outputs` directory will contain:

- `pareto_front.fasta` – optimized CDS sequences
- `fitness.csv` – metrics for each sequence (CAI, CPB, avg_MFE, CG, AUP, motif, ...)
- `summary.json` – run configuration and final statistics
- `checkpoint.pkl` – checkpoint for resuming optimization

## Real-world reference

For a full mRNA vaccine design (e.g. GL004 epitope + MITD, 583 aa), see:
`D:\workspace\projects\mRNA\GL004\cds_project\epitope_mitd\config.yaml`
