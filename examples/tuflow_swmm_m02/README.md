# TUFLOW-SWMM tutorial Module 2 benchmark

An urban development drained by a 14-pit SWMM pipe network, under 2 hours of
rain on the grid (176.5 mm, peak 161 mm/h) plus an upstream creek inflow
(peak 92 m³/s), simulated for 3 hours. The inputs are TUFLOW's own
[TUFLOW-SWMM tutorial Module 2](https://wiki.tuflow.com/TUFLOW_SWMM_Tutorial_M02)
model, translated for ANUGA + `anuga_drainage`. The aim is a code-to-code
comparison against TUFLOW-SWMM.

```bash
python fetch_data.py                         # once: downloads ~100 MB, extracts Module 2
python run_from_inp.py [swmm|pipedream] [finaltime_s]   # default: swmm 10800
```

The data are free to download and run, but have no redistribution licence, so
`fetch_data.py` downloads them into `tuflow_data/` (git-ignored) instead of
committing them. The GIS translation needs `geopandas`, `rasterio`,
`shapely>=2` and `scipy`. These are needed only by this example, not by the package.

## What is translated (`tuflow_inputs.py`)

| TUFLOW input | ANUGA / anuga_drainage |
|---|---|
| `2d_code` polygon | mesh bounding polygon (25 m² triangles, 4 m² in the development) |
| `DEM.tif` → `Read GIS Z Shape` (road crests, MAX) → `Create TIN Zpts` (development) | same order, on the 0.5 m raster; crest lines also used as breaklines |
| `Read GIS Mat` + materials file | per-triangle Manning's n (deep-water value of TUFLOW's depth-varying n) |
| `Global Rainfall BC == RF_FC04` + material losses | one `Rate_operator` per loss class; TUFLOW's stepped mm-per-interval hyetograph, less initial (20 mm) and continuing (5 mm/h) loss on pervious materials |
| `QT` line, `FC01` | `Inlet_operator` just inside the QT line |
| `HQ` line | free-outflow (`Dirichlet` below the bed) boundary |
| `sw02_001.inp` | `sw02_anuga.inp`: adds `[OPTIONS]` (TUFLOW merges them from `sw01`), drops `[STREETS]`/`[INLETS]`, and sets each inlet pit's `Ymax` to its rim (inlet elevation − invert) instead of `0`, which SWMM reads as "up to the highest pipe crown" and would flood surcharged pits with the HGL still below ground |
| inlet usage layer (12 `Combo1` inlets) | the 12 pits with an inlet are coupled (`coupled_junctions`); `Pit3`/`Pit17` are sealed manholes |
| `SXZ` pit connections | one TUFLOW-cell (5 × 5 m) footprint per pit, lowered to the inlet elevation |
| `Combo1` (grate 0.6 × 0.3 P_BAR-50 + kerb 0.9 × 0.15) | `InletSpec`: orifice area 0.297 m², weir length 1.5 m |
| outfall `SX` line | the outfall flow is returned to the 2D along that line |

## Known differences from the TUFLOW run

These have to be resolved or quantified before the comparison is meaningful:

- The three SWMM culverts (`sw01_001.inp`) are left out. They use `[STORAGE]`
  nodes, which `read_inp` doesn't parse yet.
- Inlet capture uses `calculate_Q`'s weir/orifice equations, not SWMM's
  grate/kerb capture curves. The three on-grade inlets are treated as sag inlets.
- The outfall is held FIXED at 39.5 m. TUFLOW drives it with the 2D water level.
- The HQ rating boundary is replaced by a free outflow.
- There is no sub-grid sampling, and Manning's n is constant with depth.
- SWMM options: `MIN_SURFAREA` is SWMM's default 1.167 m² and
  `SKIP_STEADY_STATE` is `NO`. TUFLOW's value of 25 m² gave SWMM a 4–13 %
  flow-routing continuity error at these flows; see the comment in
  `tuflow_inputs.py`.

There are no TUFLOW reference results yet. The tutorial zip ships a Windows
TUFLOW executable that runs this model licence-free (`Tutorial Model == ON`),
so they have to be produced on Windows.
