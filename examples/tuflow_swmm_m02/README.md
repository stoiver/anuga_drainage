# TUFLOW-SWMM tutorial Module 2 benchmark

An urban development drained by a 14-pit SWMM pipe network, under 2 hours of
rain on the grid (176.5 mm, peak 161 mm/h) plus an upstream creek inflow
(peak 92 m³/s), simulated for 3 hours. The inputs are TUFLOW's own
[TUFLOW-SWMM tutorial Module 2](https://wiki.tuflow.com/TUFLOW_SWMM_Tutorial_M02)
model, translated for ANUGA + `anuga_drainage`. The aim is a code-to-code
comparison against TUFLOW-SWMM.

```bash
python fetch_data.py                         # once: downloads ~100 MB, extracts Module 2
python run_from_inp.py [swmm|pipedream] [finaltime_s] [hec22|none]   # default: swmm 10800 hec22
python compare_tuflow.py                     # against the TUFLOW reference run (below)
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
| `[INLETS]`/`[STREETS]` (SWMM 5.2 HEC-22 capture) | `InletCapture` (`anuga_drainage.hec22`, a port of SWMM's `inlet.c`) caps each pit's draw: sag inlets on the 2D depth, on-grade inlets on the approach flow measured in a ring around the pit |
| culverts `sw01_001.inp` + `2d_bc_SWMM_Culvert_Connections` | ANUGA `Boyd_box_operator`/`Boyd_pipe_operator` between the HX and SX lines, with the `.inp`'s barrels, section, inverts, n and entry+exit losses |

## Known differences from the TUFLOW run

These have to be resolved or quantified before the comparison is meaningful:

- The three culverts are ANUGA Boyd culverts, not SWMM conduits with
  `[STORAGE]` nodes. They carry 638,700 m³ against TUFLOW's 671,300 m³.
- The HEC-22 capture is a cap on `calculate_Q`'s exchange flow. TUFLOW gives
  the inlet curves to SWMM directly.
- The outfall is held FIXED at 39.5 m. TUFLOW drives it with the 2D water level.
- The HQ rating boundary is replaced by a free outflow.
- There is no sub-grid sampling, and Manning's n is constant with depth.
- SWMM options: `MIN_SURFAREA` is SWMM's default 1.167 m² and
  `SKIP_STEADY_STATE` is `NO`. TUFLOW's value of 25 m² gave SWMM a 4–13 %
  flow-routing continuity error at these flows; see the comment in
  `tuflow_inputs.py`.

## TUFLOW reference run (Linux, under Wine)

The tutorial zip ships the Windows TUFLOW 2026.0.0 build, which runs this model
licence-free (`Tutorial Model == ON`). On Linux it runs under Wine:

```bash
python fetch_data.py --exe      # also extracts the TUFLOW build to tuflow_data/exe/
./run_tuflow.sh                 # CPU, ~3 min; results in tuflow_data/TUFLOW/results/
```

This needs one-off setup in the Wine prefix: Microsoft's VC++ 2015–2022 x64
redistributable (`wine vc_redist.x64.exe /install /quiet /norestart`). The HPC
solver uses `VCOMP140.DLL`, and Wine's builtin version ignores the requested
OpenMP thread count. The other Wine workarounds are built into
`run_tuflow.sh` and commented there. CUDA doesn't reach Wine (no `nvcuda.dll`
in mainline Wine), so the run uses `Hardware == CPU`. Tested on Wine 10.0,
Ubuntu 26.04: a cumulative mass error of 0.05 %.

TUFLOW writes `results/TS02_5m_001_CPU_swmm.{rpt,out,inp}`: the SWMM side,
which includes the culverts. Its merged `.inp` confirms TUFLOW keeps
`Ymax = 0`, `MIN_SURFAREA 25` and `ROUTING_STEP 30`, and handles the
`[INLETS]` itself. The 2D side is in `results/grids/*_Max.tif` (peak depth,
level, velocity) and `results/plot/csv/` (plot-output series).

## First comparison (pipe network, SWMM backend)

The volume each pit takes over the 3 h, from TUFLOW's `swmm.rpt` Node Inflow
Summary compared with the `accepted` column of `run_from_inp.py swmm`:

| | TUFLOW-SWMM | anuga_drainage |
|---|---|---|
| total into the pipes | 10,101 m³ (no pit floods) | 27,278 m³ (5,900 m³ floods back) |
| peak pit inflow | 0.04–0.37 m³/s | 0.07–3.05 m³/s |
| peak pipe flow, Pipe16 to the outfall | 2.32 m³/s | 6.26 m³/s |

We capture about 2.7 times as much. The gap is worst at the on-grade inlet
`Pit15` (about 12 times) and at the large sag pits `Pit11`/`Pit13`. This
points at **inlet capture** as the main difference: SWMM's grate and kerb
curves (and on-grade bypass) limit a `Combo1` inlet to a few hundred L/s,
while `calculate_Q`'s orifice on 0.297 m² under deep ponding allows several
times more. The next step is to match the capture model. Pipe-network
peaks and flooding are secondary until that's done.

## With HEC-22 capture and the culverts

`python run_from_inp.py {swmm,pipedream} 10800 hec22`, then `python compare_tuflow.py`:

| | TUFLOW | SWMM | pipedream |
|---|---|---|---|
| total into the pipes | 10,101 m³ | 15,581 m³ (1.54×) | 13,967 m³ (1.38×) |
| peak level, mean ours − TUFLOW (wet in both) | | +0.04 m | +0.05 m |
| peak level, 90th percentile of the absolute difference | | 0.16 m | 0.16 m |

The culverts brought the 2D model into line: without them the creek ponded
behind the road embankments and peak levels were 0.22 m high on average. At
the pits, peak levels now agree within a few cm. The on-grade pits (7, 9, 15)
are within about 10 %, and the shallow sag pits are close too.

The remaining excess is at the deep low-point sag pits `Pit10`/`Pit11`/`Pit13`:
1.4 to 1.8 times TUFLOW's volume. Neither SWMM report shows surcharge, so the
pipes aren't the limit. Fed TUFLOW's own peak depths (`d_Max`), our HEC-22
equations reproduce TUFLOW's peak capture at the shallow sag pits but give
0.52–0.68 m³/s at these three, where TUFLOW took 0.25–0.37 m³/s. Something in
TUFLOW-SWMM limits the capture of deep sag pits.
