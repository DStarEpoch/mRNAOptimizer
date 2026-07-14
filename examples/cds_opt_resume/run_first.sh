#!/usr/bin/env bash
# Step 1: run initial optimization
set -e
python -m cdsopt optimize protein.fa -c config.yaml
