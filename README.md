# mRNAOptimizer

mRNA sequence optimizer based on multi-objective genetic algorithms.

## Installation

```bash
git clone https://github.com/DStarEpoch/mRNAOptimizer.git --recursive
pip install -e ".[dev]"
```

### Optional: LinearFold (for faster RNA folding)

This package supports **ViennaRNA** out of the box. If you want faster folding via **LinearFold**, build it from the bundled submodule and add it to your `PATH`:

#### 1. Clone / update the submodule

```bash
git submodule update --init --recursive
```

#### 2. Build LinearFold

**Linux / macOS**
```bash
cd submodules/LinearFold
make
```

**Windows (MinGW / MSYS2)**
> **Note:** `mingw32-make` may fail on the `chmod` step (which is not needed on Windows). If so, compile directly with `g++`:
```bash
cd submodules/LinearFold
> cd submodules/LinearFold
> g++ src/LinearFold.cpp -std=c++11 -O3 -Dlv -Dis_cube_pruning -Dis_candidate_list -o bin/linearfold_v.exe
> g++ src/LinearFold.cpp -std=c++11 -O3 -Dis_cube_pruning -Dis_candidate_list -o bin/linearfold_c.exe
> ```

The build produces two binaries:
- `submodules/LinearFold/bin/linearfold_v` (Vienna parameters)
- `submodules/LinearFold/bin/linearfold_c` (CONTRAfold parameters)

On Windows, the binaries have a `.exe` extension and a `linearfold.bat` wrapper is provided for `cdsopt` auto-detection.

#### 3. Add to PATH

**Linux / macOS**
```bash
export PATH="$(pwd)/submodules/LinearFold:$PATH"
```

**Windows (PowerShell)**
```powershell
$env:PATH = "$PWD\submodules\LinearFold;$env:PATH"
```

**Windows (CMD)**
```cmd
set PATH=%CD%\submodules\LinearFold;%PATH%
```

#### 4. Verify

```bash
linearfold --help
```

If `linearfold` is found, `cdsopt` will automatically prefer it over ViennaRNA. If not, ViennaRNA is used as the fallback.

## Usage

```bash
cdsopt optimize --help
```

## Running tests

```bash
python run_tests.py -v
```
