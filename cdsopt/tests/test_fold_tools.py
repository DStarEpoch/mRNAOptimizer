# -*- coding: utf-8 -*-
import shutil
import pytest
from cdsopt.utils.fold_tools import estimate_fold


class TestEstimateFold:
    def test_vienna_returns_required_keys(self):
        result = estimate_fold("GGGAAACCC", engine="vienna")
        assert "mfe" in result
        assert "structure" in result
        assert "aup" in result
        assert isinstance(result["mfe"], float)
        assert isinstance(result["structure"], str)
        assert isinstance(result["aup"], float)
        assert 0.0 <= result["aup"] <= 1.0

    def test_vienna_structure_length_matches(self):
        seq = "GGGAAACCC"
        result = estimate_fold(seq, engine="vienna")
        assert len(result["structure"]) == len(seq)

    def test_linearfold_not_installed_raises(self):
        if shutil.which("linearfold"):
            pytest.skip("linearfold is installed")
        with pytest.raises(FileNotFoundError):
            estimate_fold("GGGAAACCC", engine="linearfold")

    @pytest.mark.skipif(not shutil.which("linearfold"), reason="linearfold not installed")
    def test_linearfold_returns_required_keys(self):
        result = estimate_fold("GGGAAACCC", engine="linearfold")
        assert "mfe" in result
        assert "structure" in result
        assert "aup" in result

    def test_auto_fallback_to_vienna(self):
        if shutil.which("linearfold"):
            pytest.skip("linearfold installed, cannot test fallback")
        result = estimate_fold("GGGAAACCC", engine="auto")
        assert "mfe" in result
        assert "structure" in result
        assert "aup" in result
