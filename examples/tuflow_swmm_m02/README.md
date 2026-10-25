# TUFLOW-SWMM tutorial Module 2 benchmark

An urban development drained by a 14-pit SWMM pipe network, under 2 hours of
rain on the grid (176.5 mm, peak 161 mm/h) plus an upstream creek inflow
(peak 92 m³/s), simulated for 3 hours. The inputs are TUFLOW's own
[TUFLOW-SWMM tutorial Module 2](https://wiki.tuflow.com/TUFLOW_SWMM_Tutorial_M02)
model, translated for ANUGA + `anuga_drainage`. The aim is a code-to-code
comparison against TUFLOW-SWMM.

```bash
python fetch_data.py                         # once: downloads ~100 MB, extracts Module 2
python run_from_inp.py [swmm|pipedream] [finaltime_s] [hec22|tuflow|none] [--gpu] [--pipedream-max-step=0.05]
                                             # default: swmm 10800 hec22, CPU
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

## TUFLOW keeps only the last `[INLETS]` line

The excess at the deep sag pits is not physics. TUFLOW-SWMM 2026.0 keeps
only the **last** `[INLETS]` line of a design. `Combo1` is defined as a GRATE
line then a CURB line, so TUFLOW runs it as a curb opening alone; EPA SWMM
5.2 combines the two lines into one combination inlet. With depth taken
above the inlet-usage `Elevation`, TUFLOW's deep sag peaks are exactly the
curb orifice 0.67·h·L·√(2g·d). The check: rerunning TUFLOW with the two lines
swapped makes its sag pits match grate-only HEC-22 instead, and its total
pipe intake goes from 10,101 to 13,785 m³.

`run_from_inp.py ... tuflow` reproduces this (`inlet_captures(last_line_wins=True)`):

| | TUFLOW | SWMM `tuflow` | pipedream `tuflow` |
|---|---|---|---|
| total into the pipes | 10,101 m³ | 9,940 m³ (0.98×) | 9,940 m³ (0.98×) |
| sag pits 2, 6, 10, 11, 13, 18, 19 | | 1.00–1.10× | 1.00–1.10× |
| shallow sag pits 4, 5 | | 1.36–1.49× | 1.36–1.49× |
| on-grade pits 7, 9, 15 | | 0.45–0.47× | 0.45–0.47× |
| peak level, mean ours − TUFLOW | | +0.06 m | +0.06 m |

The two backends agree because the pipes never limit capture. What remains
is the on-grade inlets, which capture about half of TUFLOW's. TUFLOW drives
them with a velocity from its 2D cells; we use the approach flow in a ring
around the pit. (Superseded: see the next section.)

## On-grade inlets: TUFLOW's connected cells

TUFLOW connects each pit to whole 5 m cells. A sag pit gets the cell it sits
in. An on-grade pit gets two: the cells whose centre cross-hairs its 2d_bc SX
line crosses, or (Pit9, `Conn_width -2`) the two nearest cells. TUFLOW
removes water from those cells without zeroing their velocity, and drives
the on-grade equations with the flow through them.
`tuflow_inputs.pit_cells()` reproduces the connections and matches TUFLOW's
check layer for every pit. Pit9 needs an explicit override
(`CELL_OVERRIDES`): its two nearest cell centres are 0.93 and 0.94 cells
away, and TUFLOW resolves the tie the other way. Every footprint is now its TUFLOW cells. The on-grade pits keep
their velocity (`inlet_zero_velocity`) and measure the approach flow over
those cells as mean h|u| × their width (`approach_regions`):

| pits | ratio to TUFLOW (`tuflow` capture, both backends) |
|---|---|
| sag 2, 6, 10, 11, 13, 18, 19 | 0.97–1.07 |
| on-grade 7, 9, 15 | 0.98, 0.91, 1.13 |
| shallow sag 4, 5 | 1.29–1.37 |
| **total** | **1.04** (10,472 against 10,101 m³) |

## GPU and run times

`--gpu` runs ANUGA in its unified compute mode, offloaded to the GPU. The
inlet, rain and culvert operators all run on the device, and the results are
identical to the CPU runs. For 3 h on a Ryzen 9 9955HX with an RTX 5070
Laptop GPU:

| | CPU, 12 threads | GPU |
|---|---|---|
| SWMM | 7:49 | 4:24 |
| pipedream, 0.02 s internal step | 8:42 | 5:03 |
| pipedream, 0.05 s internal step (default) | | 4:41 |
| TUFLOW HPC, 4 CPU threads, under Wine | 2:43 | |

The GPU is mostly waiting on the host: the serial Python exchange every 1 s
and the pipe solver. One pipedream exchange takes 3.5 ms at a 0.02 s internal
step, almost all of it Python overhead. BLAS threading makes no difference.
At 0.05 s the results are identical (pit volumes to 0.01 m³, heads within
1.4 mm), and pipedream's own mass loss falls from 113 to 86 m³. This network
never surcharges, so 0.05 s is not tested under surcharge.
