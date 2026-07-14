@echo off
REM Step 2: resume from checkpoint
python -m cdsopt resume ./outputs/checkpoint.pkl --output-dir ./outputs_resumed
