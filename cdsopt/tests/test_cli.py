# -*- coding: utf-8 -*-
"""Integration tests for the CLI."""
from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
from click.testing import CliRunner

from cdsopt.main import app


@pytest.fixture
def runner() -> CliRunner:
    return CliRunner()


@pytest.fixture
def protein_txt(tmp_path: Path) -> Path:
    p = tmp_path / "protein.txt"
    p.write_text("MKW", encoding="utf-8")
    return p


@pytest.fixture
def protein_fasta(tmp_path: Path) -> Path:
    p = tmp_path / "protein.fa"
    p.write_text(
        ">sp|P12345|TEST\nMKW\n", encoding="utf-8"
    )
    return p


class TestOptimizeCommand:
    def test_runs_and_creates_output(self, runner: CliRunner, protein_txt: Path, tmp_path: Path):
        out_dir = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "optimize",
                str(protein_txt),
                "--pop-size", "4",
                "--generations", "2",
                "--processes", "1",
                "--seed", "42",
                "--output-dir", str(out_dir),
            ],
        )
        assert result.exit_code == 0, result.output
        assert (out_dir / "pareto_front.fasta").exists()
        assert (out_dir / "fitness.csv").exists()
        assert (out_dir / "summary.json").exists()
        assert (out_dir / "checkpoint.pkl").exists()

    def test_fasta_input(self, runner: CliRunner, protein_fasta: Path, tmp_path: Path):
        out_dir = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "optimize",
                str(protein_fasta),
                "--pop-size", "4",
                "--generations", "1",
                "--output-dir", str(out_dir),
            ],
        )
        assert result.exit_code == 0, result.output
        assert (out_dir / "pareto_front.fasta").exists()

    def test_yaml_config_overrides(self, runner: CliRunner, protein_txt: Path, tmp_path: Path):
        config = tmp_path / "config.yaml"
        config.write_text(
            "pop_size: 6\ngenerations: 3\nseed: 7\nenable_cpb: true\n",
            encoding="utf-8",
        )
        out_dir = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "optimize",
                str(protein_txt),
                "--config", str(config),
                "--output-dir", str(out_dir),
            ],
        )
        assert result.exit_code == 0, result.output
        assert (out_dir / "pareto_front.fasta").exists()
        # Verify that YAML values were used (pop_size=6 implies 6 sequences)
        fasta = (out_dir / "pareto_front.fasta").read_text(encoding="utf-8")
        assert fasta.count(">") == 6

    def test_cli_arg_overrides_yaml(self, runner: CliRunner, protein_txt: Path, tmp_path: Path):
        config = tmp_path / "config.yaml"
        config.write_text(
            "pop_size: 10\n",
            encoding="utf-8",
        )
        out_dir = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "optimize",
                str(protein_txt),
                "--config", str(config),
                "--pop-size", "4",
                "--generations", "1",
                "--output-dir", str(out_dir),
            ],
        )
        assert result.exit_code == 0, result.output
        fasta = (out_dir / "pareto_front.fasta").read_text(encoding="utf-8")
        assert fasta.count(">") == 4

    def test_yaml_not_mutate_idx_list(self, runner: CliRunner, protein_txt: Path, tmp_path: Path):
        config = tmp_path / "config.yaml"
        config.write_text(
            "pop_size: 4\ngenerations: 1\nnot_mutate_idx: [0, 1]\n",
            encoding="utf-8",
        )
        out_dir = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "optimize",
                str(protein_txt),
                "--config", str(config),
                "--output-dir", str(out_dir),
            ],
        )
        assert result.exit_code == 0, result.output

    def test_init_cds_accepted(self, runner: CliRunner, protein_txt: Path, tmp_path: Path):
        init_cds = tmp_path / "init.fa"
        init_cds.write_text(
            ">cds1\nAUGAAAAUGG\n>cds2\nAUGAAGUGG\n",
            encoding="utf-8",
        )
        out_dir = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "optimize",
                str(protein_txt),
                "--pop-size", "4",
                "--generations", "1",
                "--init-cds", str(init_cds),
                "--output-dir", str(out_dir),
            ],
        )
        assert result.exit_code == 0, result.output

    def test_init_cds_mismatch_raises(self, runner: CliRunner, protein_txt: Path, tmp_path: Path):
        init_cds = tmp_path / "init.fa"
        # AUGUUUAAA encodes MFK, not MKW
        init_cds.write_text(
            ">cds1\nAUGUUUAAA\n",
            encoding="utf-8",
        )
        out_dir = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "optimize",
                str(protein_txt),
                "--pop-size", "4",
                "--generations", "1",
                "--init-cds", str(init_cds),
                "--output-dir", str(out_dir),
            ],
        )
        assert result.exit_code != 0
        exc_msg = str(result.exception).lower() if result.exception else ""
        assert "match the protein" in exc_msg or "not a valid synonymous codon" in exc_msg

    def test_yaml_forbidden_motifs(self, runner: CliRunner, protein_txt: Path, tmp_path: Path):
        config = tmp_path / "config.yaml"
        config.write_text(
            "pop_size: 4\ngenerations: 1\nforbidden_motifs:\n  enzymes:\n    - BspQ1\n  polyt_min_len: 6\n",
            encoding="utf-8",
        )
        out_dir = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "optimize",
                str(protein_txt),
                "--config", str(config),
                "--output-dir", str(out_dir),
            ],
        )
        assert result.exit_code == 0, result.output
        assert (out_dir / "pareto_front.fasta").exists()


class TestResumeCommand:
    def test_resume_from_checkpoint(self, runner: CliRunner, protein_txt: Path, tmp_path: Path):
        out_dir = tmp_path / "out"
        # First run
        result = runner.invoke(
            app,
            [
                "optimize",
                str(protein_txt),
                "--pop-size", "4",
                "--generations", "2",
                "--seed", "42",
                "--output-dir", str(out_dir),
            ],
        )
        assert result.exit_code == 0, result.output
        checkpoint = out_dir / "checkpoint.pkl"
        assert checkpoint.exists()

        # Resume
        out_dir2 = tmp_path / "out2"
        result2 = runner.invoke(
            app,
            [
                "resume",
                str(checkpoint),
                "--output-dir", str(out_dir2),
            ],
        )
        assert result2.exit_code == 0, result2.output
        assert (out_dir2 / "pareto_front.fasta").exists()
