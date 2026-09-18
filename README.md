# anuga_drainage — coupling ANUGA to SWMM and pipedream

[![Documentation Status](https://readthedocs.org/projects/anuga-drainage/badge/?version=latest)](https://anuga-drainage.readthedocs.io/en/latest/?badge=latest)

Couples the 2D hydrodynamic shallow-water model **ANUGA** (overland / surface
flow) with a 1D stormwater / sewer network solver — **SWMM** (via `pyswmm`) or
**pipedream** (via `pipedream_solver`). The 2D model handles the surface, the 1D
model the underground pipe network, and the two exchange water at
inlets/manholes every timestep using weir/orifice equations.

📖 **Full documentation:** <https://anuga-drainage.readthedocs.io/>

## Installation

The quickest path is the bundled conda environment (ANUGA 3.3.6 from
conda-forge + this package + both backends), created from the repo root:

```bash
conda env create -f environment.yml
conda activate anuga_drainage
```

Or, on an **existing conda environment** that already provides ANUGA (it is
heavy and *not* pip-installed here):

```bash
pip install -e .              # the anuga_drainage package (numpy, pandas)
pip install -e .[swmm]        # + SWMM backend (standard pyswmm release)
pip install -e .[pipedream]   # + pipedream backend (from git; see note below)
```

> The `[pipedream]` extra installs pipedream from the
> [anuga-community fork](https://github.com/anuga-community/pipedream)'s
> `anuga` branch: upstream master plus the numpy 2 and pandas 3 fixes upstream
> has not merged. See [`CLAUDE.md`](CLAUDE.md) for this and other environment
> constraints (notably that stock pyswmm 2.1 couples in whole-second steps).

## Running an example

Each example runs from its own directory (scripts use relative paths):

```bash
cd examples/simple_culvert_example
python run_swmm_short.py     # ANUGA + SWMM, inlets hand-built
python run_pipedream.py      # ANUGA + pipedream, network hand-built
python run_boyd.py           # ANUGA-only reference (Boyd culvert operator)
python run_from_inp.py       # ANUGA + SWMM straight from the .inp (default)
python run_from_inp.py pipedream   # ...same .inp, pipedream backend
```

The `run_from_inp.py` scripts are the **living example of `couple_from_inp`**:
the sewer network *and* the ANUGA coupling inlets are built from a single SWMM
`.inp` for whichever backend you pick, so the only model-specific code is the
2D domain, the inflow, and the evolve loop — contrast the `run_swmm*` /
`run_pipedream.py` scripts, which hand-assemble the inlets (and, for pipedream,
the whole superjunction/superlink network). Both `simple_culvert_example/` and
`real_example/` ship one. With `log_hydrographs=True` they also write per-inlet
CSVs to `hydrographs_<backend>/` for the viewer.

ANUGA writes `<name>.sww`; scripts print a running mass-balance `loss` and may
save `Figure*.png`.

## The package

- `calculate_Q(...)` — the weir/orifice exchange-flux physics
  (Leandro & Martins, 2016). Positive Q = surface → pipe; negative = surcharge.
- `Coupler` with `SwmmBackend` / `PipedreamBackend` — drives the per-step
  exchange (read depths → `calculate_Q` → smooth → optional clamp → step the 1D
  model → feed the realised flow back to ANUGA). All the coupled `run_*` example
  scripts use it. Only `[JUNCTIONS]` couple to the surface; **outfalls are
  boundaries** — outfall water leaves the model at a free-draining boundary by
  default (tracked by `backend.outfall_volume()`, counted as `loss` in the
  audit), or set `outfall_inlet=<index>` to return it to the 2D surface at a
  chosen inlet. See the docs' *Where outfall water goes*.
- `initialize_inlets(...)` / `read_inp_coordinates(...)` — build ANUGA inlet
  operators from a SWMM `.inp` file.
- `InletSpec` / `INLET_LIBRARY` / `load_inlet_library(...)` — a named inlet
  **asset catalogue** (clear opening area + weir perimeter, with optional
  blockage that derates both). Assign one per junction via
  `couple_from_inp(inlet_specs={name: "Grate_600x600"}, blockage=...)` to drive
  that junction's exchange flux from a catalogued grate/lintel — **decoupled from
  the surface footprint** (the coupling region / pipedream storage still come
  from `manhole_area`/`inlet_polygons`). Catalogues can be loaded from TOML.
- `HydrographLogger` — optional per-inlet hydrograph logging. Pass
  `couple_from_inp(log_hydrographs=True)` (reachable at
  `coupling.coupler.logger`); records depth, 1D head, capture/surcharge and
  cumulative volumes, and `write_csv(dir)` emits one CSV per inlet — read by the
  viewer below.

## Hydrograph viewer

A Tkinter dashboard for inspecting per-inlet hydrograph CSVs (e.g. those written
by `HydrographLogger`):

```bash
pip install -e .[viewer]     # + matplotlib (Tkinter is stdlib)
anuga-drainage-viewer        # or: python -m anuga_drainage.viewer
```

Pick a folder of `hydrograph_*.csv`; per inlet it draws four diagnostic plots,
and a **View** menu adds a folder-combined hydrograph. The pure
`combine_hydrographs()` helper (used by the combined view) is importable and
tested on its own.

## Tests

```bash
pip install -e .[test]
pytest
```

The pure logic (flux physics, smoothing, geometry, `.inp` parsing) is tested
without requiring ANUGA; tests that need ANUGA's gravity constant skip when it
is absent.

## Examples

- `simple_culvert_example/` — canonical minimal channel + culvert; the clearest
  reference for the coupling pattern (SWMM, pipedream and Boyd variants).
- `culvert_example/`, `yshape/`, `real_example/` — further SWMM/pipedream cases.
- `pyswmm_with_openings/` — **legacy/reference only**: built against a patched
  pyswmm fork that is not in the standard release; not part of the active path.

See [`CLAUDE.md`](CLAUDE.md) for the architecture and coupling-loop details.
