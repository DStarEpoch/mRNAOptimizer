# cdsopt resume example

This example shows how to resume optimization from a checkpoint.

## Step 1: Run an initial short optimization

```bash
python -m cdsopt optimize protein.fa -c config.yaml
```

This creates `./outputs/checkpoint.pkl`.

## Step 2: Resume from the checkpoint

```bash
python -m cdsopt resume ./outputs/checkpoint.pkl --output-dir ./outputs_resumed
```

The resumed run continues from the last saved generation using the same configuration.

## Notes

- The checkpoint preserves population, fitness values, RNG state, and configuration.
- You can change `--processes` when resuming to use a different number of CPU cores.
- Forbidden motif constraints are restored from the checkpoint automatically.
