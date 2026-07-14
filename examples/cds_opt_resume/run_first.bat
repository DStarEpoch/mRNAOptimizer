@echo off
REM Step 1: run initial optimization
python -m cdsopt optimize protein.fa -c config.yaml
