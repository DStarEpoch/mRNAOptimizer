#!/usr/bin/env bash
# Run cdsopt optimize example
set -e
python -m cdsopt optimize protein.fa -c config.yaml
