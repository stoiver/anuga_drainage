"""TUFLOW-SWMM tutorial Module 2 as an anuga_drainage benchmark (pipe network).

Rain on the grid over an urban development drained by a 14-pit SWMM pipe network,
with an upstream creek inflow. The inputs are TUFLOW's own (fetch them first with
``python fetch_data.py``), translated by ``tuflow_inputs.py``:

* 2D: the 2d_code area, DEM + road-crest Z shapes + development TIN, materials
  (Manning's n, rainfall losses), the FC01 inflow hydrograph (QT line) and an
  outflow boundary on the HQ line;
* rain: global rainfall RF_FC04 (176.5 mm, peak 161 mm/h) less each material's
  initial/continuing loss;
* 1D: sw02_001.inp. 12 pits carry an inlet (from the inlet-usage layer) and are
  coupled to the surface; Pit3/Pit17 are sealed manholes (pipe-only). The
  outfall discharges back onto the 2D creek bank along TUFLOW's SX line.

Inlet capture: calculate_Q (with an equivalent Combo1 InletSpec) sets the
head-driven exchange, and by default each pit's surface->pipe flow is capped at
its SWMM/HEC-22 inlet capacity (anuga_drainage.hec22, built from the .inp's
[INLETS]/[STREETS] and the inlet-usage layer), as TUFLOW-SWMM does: ON_SAG
inlets by ponded depth, ON_GRADE inlets by the approaching surface flow.
``hec22`` treats Combo1 as a grate + curb combination, as EPA SWMM does.
``tuflow`` reproduces TUFLOW-SWMM 2026.0, which keeps only the last [INLETS]
line of a design, so Combo1 captures as its curb opening alone. ``none``
drops the cap (calculate_Q alone).

Culverts: TUFLOW routes the three road culverts (sw01_001.inp) through SWMM;
here they are ANUGA Boyd culverts between TUFLOW's HX (upstream) and SX
(downstream) connection lines, with the .inp's barrels, sections, inverts,
roughness and losses. They matter: they carry the creek under the raised road
crests, and without them it ponds behind the embankments.

Not yet modelled (so this is not a like-for-like TUFLOW comparison yet): the
2D tailwater on the outfall (held FIXED at 39.5), the HQ rating (a
free-outflow boundary instead) and sub-grid sampling. Depth-varying Manning's n
is optional (--depth-varying-n): updated once per coupling step rather than
every timestep, and on the GPU it needs an anuga with
Domain.sync_friction_to_device (anuga-community/anuga_core#458).

    python run_from_inp.py [swmm|pipedream] [finaltime_s] [hec22|tuflow|none] [--gpu]
                           [--pipedream-max-step=0.05] [--depth-varying-n]
                                        # default: swmm 10800 hec22
"""
import sys
import warnings

import numpy as np
import anuga

from anuga_drainage import couple_from_inp, InletSpec
import tuflow_inputs as ti

warnings.filterwarnings("ignore", module="pyogrio")

flags = [a for a in sys.argv[1:] if a.startswith('--')]
args = [a for a in sys.argv[1:] if not a.startswith('--')]
gpu = '--gpu' in flags                 # ANUGA's unified compute mode, offloaded to the GPU
max_step = 0.05                        # pipedream internal step (s); see couple_from_inp
depth_n = '--depth-varying-n' in flags  # TUFLOW's depth-varying Manning's n, updated each dt
for f in flags:
    if f.startswith('--pipedream-max-step='):
        max_step = float(f.split('=', 1)[1])
    elif f not in ('--gpu', '--depth-varying-n'):
        sys.exit(f'unknown option {f!r}')
backend = args[0] if len(args) > 0 else 'swmm'
finaltime = float(args[1]) if len(args) > 1 else 3 * 3600.0   # .tcf End Time == 3 h
capture = args[2] if len(args) > 2 else 'hec22'               # inlet capture cap
if capture not in ('hec22', 'tuflow', 'none'):
    sys.exit(f'capture must be hec22, tuflow or none, got {capture!r}')
run = (f'{backend}_{capture}' + ('_gpu' if gpu else '') + ('_dn' if depth_n else '')
       + (f'_ms{max_step:g}' if backend == 'pipedream' and max_step != 0.05 else ''))
dt = 1.0                       # coupling/yield step (SWMM exchanges at whole seconds)
outputstep = 300.0             # .tcf Map Output Interval
base_area = 25.0               # max triangle area (m^2) — ~ TUFLOW's 5 m cell
dev_area = 4.0                 # finer in the development, where the pits are
pit_cell = ti.TUFLOW_CELL      # TUFLOW cell size: the width of each connected cell

# SWMM's Combo1 inlet (GRATE 0.6 x 0.3 m P_BAR-50 + CURB 0.9 x 0.15 m) as the
# area/perimeter calculate_Q needs. Orifice area: grate open area (P_BAR-50 is
# 90 % open in SWMM's grate table) + kerb opening. Weir length (HEC-22 sag
# combination): grate perimeter against the kerb (L + 2W) + kerb beyond the grate.
COMBO1 = InletSpec('Combo1', clear_area=0.6 * 0.3 * 0.90 + 0.9 * 0.15,
                   effective_perimeter=(0.6 + 2 * 0.3) + (0.9 - 0.6))

#------------------------------------------------------------------------------
# 1. ANUGA domain
#------------------------------------------------------------------------------
from shapely import contains_xy
from shapely.geometry import LineString, Polygon

polygon = ti.active_polygon()
bc = ti.boundary_lines()
_, qt_line = bc['QT']
_, hq_line = bc['HQ']
active = Polygon(polygon)

dev = ti._single_polygon(ti._layer(ti.TIN_LAYERS[0]).geometry.iloc[0]).intersection(
    active.buffer(-1.0))
interior_regions = [[[list(c[:2]) for c in dev.exterior.coords[:-1]], dev_area]]
# Road crests as breaklines so the mesh follows the raised Z shapes.
breaklines = []
for line in ti._lines(ti._layer(ti.ZSHAPE_LAYERS[0])):
    part = line.intersection(active.buffer(-1.0))
    for seg in getattr(part, 'geoms', [part]):
        if seg.length > 1.0:
            breaklines.append([list(c[:2]) for c in seg.coords])

domain = anuga.create_domain_from_regions(
    polygon, boundary_tags=ti.boundary_tags(polygon, hq_line, 'outflow'),
    maximum_triangle_area=base_area, interior_regions=interior_regions,
    breaklines=breaklines, use_cache=False, verbose=False)
domain.set_name(f'tuflow_m02_{run}')
if gpu:
    anuga.set_gpu_offload(True, strict=True)   # call before the domain goes unified
    domain.set_multiprocessor_mode(2)
domain.set_minimum_storable_height(0.001)
print(f'Mesh: {domain.number_of_triangles} triangles')

x, y, z = ti.elevation_grid()
cx, cy = domain.get_centroid_coordinates(absolute=True).T
zc = ti.elevation_function(x, y, z)(cx, cy)

# Coupling footprints: the TUFLOW cells each pit is connected to (one for a
# sag pit, two for an on-grade pit). The SX "Z" flag lowers them to the inlet level.
inlets = ti.pit_inlets()
cells = ti.pit_cells(inlets)
inlet_polygons = {}
for name, r in inlets.items():
    poly = ti.cells_polygon(cells[name])
    inlet_polygons[name] = poly
    mask = contains_xy(Polygon(poly), cx, cy)
    zc[mask] = np.minimum(zc[mask], r.Elevation)
domain.set_quantity('elevation', zc, location='centroids')
domain.set_quantity('stage', expression='elevation')

mats = ti.read_materials()
mat = ti.material_ids(cx, cy)
domain.set_quantity('friction', mats.loc[mat, 'n'].to_numpy(), location='centroids')
if depth_n:
    # Pasture, buildings and the vegetated creek have a depth-varying n. It is
    # updated once per coupling step from the centroid depth (TUFLOW does it
    # every timestep). Buildings matter most: n 0.02 for a thin film on a roof,
    # 0.4 once the water is deep enough to be blocked. On the GPU the solver
    # reads the device copy of friction, so each update is pushed there.
    sync_friction = getattr(domain, 'sync_friction_to_device', None)
    if gpu and sync_friction is None:
        sys.exit('--depth-varying-n with --gpu needs an anuga with '
                 'Domain.sync_friction_to_device (anuga_core#458): without it '
                 'the GPU silently keeps the initial n')
    dn_idx, n_of_depth = ti.depth_varying_n(mats, mat)
    friction_cv = domain.quantities['friction'].centroid_values
    stage_cv = domain.quantities['stage'].centroid_values
    bed_cv = domain.quantities['elevation'].centroid_values

    def update_friction():
        friction_cv[dn_idx] = n_of_depth(np.maximum(stage_cv[dn_idx] - bed_cv[dn_idx], 0.0))
        if sync_friction is not None:
            sync_friction()                  # no-op on the CPU
    update_friction()

# Rain on the grid: one Rate_operator per distinct loss pair, each applying the
# hyetograph less that material's initial/continuing loss.
rf_t, rf_mm = ti.rainfall('RF_FC04')
rain_ops = []
for (IL, CL), ids in mats.groupby(['IL', 'CL']).groups.items():
    idx = np.flatnonzero(np.isin(mat, list(ids)))
    if idx.size:
        rain_ops.append(anuga.Rate_operator(domain, rate=ti.net_rain_rate(rf_t, rf_mm, IL, CL),
                                            indices=idx, label=f'rain_IL{IL:g}_CL{CL:g}'))

# Upstream creek inflow FC01 along TUFLOW's QT line. That line lies on the
# domain edge, which no triangle crosses, so apply it just inside.
qt = LineString(qt_line)
qt = next(o for o in (qt.offset_curve(d) for d in (2.0, -2.0))
          if active.contains(o.interpolate(0.5, normalized=True)))
inflow_op = anuga.Inlet_operator(domain, [list(qt.coords[0]), list(qt.coords[-1])],
                                 Q=ti.hydrograph('FC01'))

# Road culverts (sw01_001.inp) as Boyd culverts between TUFLOW's HX/SX lines.
culvert_ops = []
for cv in ti.culverts():
    common = dict(losses=cv['losses'], barrels=cv['barrels'], manning=cv['n'],
                  exchange_lines=cv['exchange_lines'], invert_elevations=cv['inverts'],
                  label=cv['name'])
    if cv['shape'] == 'CIRCULAR':
        culvert_ops.append(anuga.Boyd_pipe_operator(domain, diameter=cv['height'], **common))
    else:
        culvert_ops.append(anuga.Boyd_box_operator(domain, width=cv['width'],
                                                   height=cv['height'], **common))

# HQ boundary approximated by a free outflow (stage below the bed); walls elsewhere.
z_out = zc[domain.tri_full_flag.astype(bool)].min()
domain.set_boundary({'outflow': anuga.Dirichlet_boundary([z_out - 1.0, 0.0, 0.0]),
                     'wall': anuga.Reflective_boundary(domain)})

#------------------------------------------------------------------------------
# 2. Coupling inlets + 1D network
#------------------------------------------------------------------------------
inp_path = ti.write_swmm_inp(f'sw02_anuga_{run}.inp', finaltime,   # per run: SWMM writes .rpt/.out beside it
                              rims={name: r.Elevation for name, r in inlets.items()})
# On-grade pits: like TUFLOW's SX cells, the footprint keeps its velocity, and
# the approach flow is the flow through the connected cells (h|u| times their
# width across the flow, one cell width each).
on_grade = [name for name, r in inlets.items() if r.Placement == 'ON_GRADE']
approach_regions = {name: (inlet_polygons[name], pit_cell * len(cells[name]))
                    for name in on_grade}
coupling = couple_from_inp(domain, inp_path, backend=backend,
                           coupled_junctions=list(inlets),
                           inlet_polygons=inlet_polygons,
                           inlet_specs={name: COMBO1 for name in inlets},
                           inlet_capture=None if capture == 'none' else
                           ti.inlet_captures(inlets, last_line_wins=capture == 'tuflow'),
                           approach_regions=approach_regions if capture != 'none' else None,
                           inlet_zero_velocity={name: False for name in on_grade},
                           time_average=10.0, clamp=True,
                           internal_links=4, pipedream_max_step=max_step,
                           log_hydrographs=True)
print(f'Coupled {len(coupling.inlets)} pits {sorted(coupling.inlets)}; '
      f'sealed: {sorted(set(coupling.inp.junctions.name) - set(inlets))}  '
      f'(backend={backend}, capture={capture})')

# The outfall discharges back onto the 2D creek bank along TUFLOW's SX line.
_, sx_line = next((r.Type, [list(c[:2]) for c in r.geometry.coords])
                  for _, r in ti._layer('2d_bc_SWMM_Pipe_Network_Connections_001_L').iterrows()
                  if r.Type == 'SX')
outfall_op = anuga.Inlet_operator(domain, [sx_line[0], sx_line[-1]], Q=0.0)

#------------------------------------------------------------------------------
# 3. VolumeBalance
#------------------------------------------------------------------------------
# The outfall flow leaves the pipe (outfall_volume) and re-enters ANUGA through
# outfall_op, so outfall_op is a source on the ANUGA side and the two cancel in
# `loss`. (It isn't a coupling inlet, so outfall_inlet stays None.)
coupling.add_volume_balance(inflow_operators=[inflow_op, outfall_op, *rain_ops])

#------------------------------------------------------------------------------
# 4. Evolve loop
#------------------------------------------------------------------------------
v_out = 0.0
for t in domain.evolve(yieldstep=dt, outputstep=outputstep, finaltime=finaltime):
    coupling.step(dt)
    if depth_n:
        update_friction()
    v_new = coupling.backend.outfall_volume()
    outfall_op.set_Q((v_new - v_out) / dt)      # this step's outfall flow, returned
    v_out = v_new
    if domain.yieldstep_counter % 300 == 0:
        domain.print_timestepping_statistics()

#------------------------------------------------------------------------------
# 5. Report
#------------------------------------------------------------------------------
print()
print(coupling.volume_balance.summary())
rain = sum(op.cumulative_influx for op in rain_ops)
for op in culvert_ops:
    print(f'Culvert {op.label}: {op.accumulated_flow:.0f} m^3')
print(f'Rain applied {rain:.0f} m^3, FC01 inflow {inflow_op.get_total_applied_volume():.0f} m^3, '
      f'pipe outfall {coupling.backend.outfall_volume():.0f} m^3')

# "requested" is calculate_Q's draw (what the logger records); "accepted" is
# the net volume that actually entered the pipe there (for SWMM, inflow minus
# flooding back out) — the difference is the pit rejecting/overflowing.
logger = coupling.coupler.logger
accepted = dict(zip(coupling.inlets, coupling.backend.coupling_inflow_volumes()))
rows = []
print(f'\n{"pit":6s} {"peak draw":>11s} {"requested":>11s} {"accepted":>11s}')
for name in sorted(coupling.inlets, key=lambda n: int(n[3:])):
    df = logger.to_dataframe(name)
    requested = df.Cum_Captured_m3.iloc[-1] - df.Cum_Surcharged_m3.iloc[-1]
    rows.append({'pit': name, 'peak_draw_cms': df.Captured_Q_cms.max(),
                 'requested_m3': requested, 'accepted_m3': accepted[name]})
    print(f'{name:6s} {df.Captured_Q_cms.max():6.3f} m3/s {requested:8.1f} m3 {accepted[name]:8.1f} m3')
# Per-pit summary for compare_tuflow.py.
import pandas as pd
pd.DataFrame(rows).to_csv(f'pits_{run}.csv', index=False)
paths = logger.write_csv(f'hydrographs_{run}')
print(f'Wrote {len(paths)} hydrograph CSV(s) to hydrographs_{run}/')

coupling.close()
