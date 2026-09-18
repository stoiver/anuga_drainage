# Installation

`anuga_drainage` assumes you already have a working **ANUGA** in your
environment — ANUGA is not installed by this package. On top of that:

```bash
pip install -e .            # the anuga_drainage package (needs numpy, pandas)
pip install pyswmm          # SWMM backend (standard PyPI release, >= 2.1)
pip install "pipedream-solver @ git+https://github.com/anuga-community/pipedream.git@anuga"
```

The package itself only depends on `numpy` and `pandas`; the two 1D backends are
optional extras, installed only for the backend(s) you use:

```bash
pip install -e .[swmm]            # pyswmm
pip install -e .[pipedream]       # pipedream (from the community fork, see below)
pip install -e .[test]            # pytest
```

## Backend notes

```{admonition} pipedream comes from the anuga-community fork, not PyPI
:class: warning
Neither the released `pipedream-solver` (0.2.2) nor upstream master works with
a current numpy/pandas: 0.2.2 uses `np.bool8`, removed in numpy 2.x, and
master's `SuperLink(...)` construction crashes under pandas 3. The
[anuga-community fork](https://github.com/anuga-community/pipedream)'s `anuga`
branch is upstream master plus both fixes (they are also open upstream as
mdbartos/pipedream#73 and #74). The `[pipedream]` extra installs it.
```

### SWMM / pyswmm 2.1 stepping constraints

Stock pyswmm 2.1 is **whole-second resolution**: the coupling stride must be an
integer number of seconds (`int(dt)`), and SWMM coupling therefore exchanges at
1-second granularity. The `Coupler` handles this for you. Sub-second coupling is
only available on the **pipedream** path (its step is pure Python).

## From-scratch conda environment

For a reproducible setup, the repository ships an `environment.yml` that builds
a conda environment with ANUGA from conda-forge plus this package and both
backends:

```bash
conda env create -f environment.yml
```

## Running the tests

The package's own physics is tested with `pytest`:

```bash
pip install -e .[test]
pytest
```

The pure-logic tests (geometry, the `.inp` parser, `calculate_Q` with an
explicit gravity) run without ANUGA; ANUGA / pyswmm / pipedream-dependent tests
skip automatically when those aren't installed.
