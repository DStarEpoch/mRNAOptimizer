#!/usr/bin/env bash
# Step 2: resume from checkpoint
set -e
python -m cdsopt resume ./outputs/checkpoint.pkl --output-dir ./outputs_resumed
