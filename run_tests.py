# -*- coding: utf-8 -*-
"""
Unified test entry point for cdsopt.
Run with:  python run_tests.py
Or:        python run_tests.py -v
"""
import sys
import subprocess
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Run cdsopt test suite")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    parser.add_argument("-k", "--keyword", type=str, default="", help="Filter tests by keyword")
    parser.add_argument("--cov", action="store_true", help="Enable coverage report")
    parser.add_argument("--log-file", type=str, default="test_log.txt", help="Log file path")
    args = parser.parse_args()

    project_root = Path(__file__).parent
    test_dir = project_root / "cdsopt" / "tests"

    pytest_args = [
        sys.executable, "-m", "pytest",
        str(test_dir),
        "-v",
        "--tb=short",
        "-rA",
        f"--log-file={args.log_file}",
        "--log-file-level=INFO",
    ]

    if args.verbose:
        pytest_args.extend(["-s", "--log-cli-level=INFO"])
    if args.keyword:
        pytest_args.extend(["-k", args.keyword])
    if args.cov:
        pytest_args.extend(["--cov=cdsopt", "--cov-report=term-missing"])

    print("=" * 60)
    print("cdsopt Test Suite")
    print("=" * 60)
    print(f"Command: {' '.join(pytest_args)}")
    print(f"Test directory: {test_dir}")
    print("-" * 60)

    result = subprocess.run(pytest_args, cwd=str(project_root))

    print("-" * 60)
    if result.returncode == 0:
        print("All tests PASSED.")
    else:
        print(f"Tests FAILED with exit code {result.returncode}.")
    print("=" * 60)
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
