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

Not yet modelled (so this is not a like-for-like TUFLOW comparison yet): the
three SWMM culverts (sw01_001.inp), SWMM's inlet capture curves (we use
calculate_Q's weir/orifice with an equivalent Combo1 InletSpec; on-grade inlets
are treated as sag), the 2D tailwater on the outfall (held FIXED at 39.5), the
HQ rating (a free-outflow boundary instead), sub-grid sampling and
depth-varying Manning's n.

    python run_from_inp.py [swmm|pipedream] [finaltime_s]   # default: swmm 10800
"""
import sys
import warnings

import numpy as np
import anuga

from anuga_drainage import couple_from_inp, InletSpec
import tuflow_inputs as ti

warnings.filterwarnings("ignore", module="pyogrio")

backend = sys.argv[1] if len(sys.argv) > 1 else 'swmm'
finaltime = float(sys.argv[2]) if len(sys.argv) > 2 else 3 * 3600.0   # .tcf End Time == 3 h
dt = 1.0                       # coupling/yield step (SWMM exchanges at whole seconds)
outputstep = 300.0             # .tcf Map Output Interval
base_area = 25.0               # max triangle area (m^2) — ~ TUFLOW's 5 m cell
dev_area = 4.0                 # finer in the development, where the pits are
pit_cell = ti.TUFLOW_CELL      # coupling footprint: one TUFLOW cell per pit (SX)

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
domain.set_name(f'tuflow_m02_{backend}')
domain.set_minimum_storable_height(0.001)
print(f'Mesh: {domain.number_of_triangles} triangles')

x, y, z = ti.elevation_grid()
cx, cy = domain.get_centroid_coordinates(absolute=True).T
zc = ti.elevation_function(x, y, z)(cx, cy)

# Coupling footprints: one TUFLOW cell, aligned with TUFLOW's rotated grid,
# centred on each inlet pit. The SX "Z" flag lowers that cell to the inlet level.
inlets = ti.pit_inlets()
angle = np.arctan2(-200.0, 855.0)            # .tgc Origin -> Orientation
c, s = np.cos(angle), np.sin(angle)
square = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) * pit_cell / 2
square = square @ np.array([[c, s], [-s, c]])
inlet_polygons = {}
for name, r in inlets.items():
    poly = (square + [r.geometry.x, r.geometry.y]).tolist()
    inlet_polygons[name] = poly
    mask = contains_xy(Polygon(poly), cx, cy)
    zc[mask] = np.minimum(zc[mask], r.Elevation)
domain.set_quantity('elevation', zc, location='centroids')
domain.set_quantity('stage', expression='elevation')

mats = ti.read_materials()
mat = ti.material_ids(cx, cy)
domain.set_quantity('friction', mats.loc[mat, 'n'].to_numpy(), location='centroids')

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

# HQ boundary approximated by a free outflow (stage below the bed); walls elsewhere.
z_out = zc[domain.tri_full_flag.astype(bool)].min()
domain.set_boundary({'outflow': anuga.Dirichlet_boundary([z_out - 1.0, 0.0, 0.0]),
                     'wall': anuga.Reflective_boundary(domain)})

#------------------------------------------------------------------------------
# 2. Coupling inlets + 1D network
#------------------------------------------------------------------------------
inp_path = ti.write_swmm_inp('sw02_anuga.inp', finaltime,
                              rims={name: r.Elevation for name, r in inlets.items()})
coupling = couple_from_inp(domain, inp_path, backend=backend,
                           coupled_junctions=list(inlets),
                           inlet_polygons=inlet_polygons,
                           inlet_specs={name: COMBO1 for name in inlets},
                           time_average=10.0, clamp=True,
                           internal_links=4, pipedream_max_step=0.02,
                           log_hydrographs=True)
print(f'Coupled {len(coupling.inlets)} pits {sorted(coupling.inlets)}; '
      f'sealed: {sorted(set(coupling.inp.junctions.name) - set(inlets))}  (backend={backend})')

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
print(f'Rain applied {rain:.0f} m^3, FC01 inflow {inflow_op.get_total_applied_volume():.0f} m^3, '
      f'pipe outfall {coupling.backend.outfall_volume():.0f} m^3')

# "requested" is calculate_Q's draw (what the logger records); "accepted" is
# the net volume that actually entered the pipe there (for SWMM, inflow minus
# flooding back out) — the difference is the pit rejecting/overflowing.
logger = coupling.coupler.logger
accepted = dict(zip(coupling.inlets, coupling.backend.coupling_inflow_volumes()))
print(f'\n{"pit":6s} {"peak draw":>11s} {"requested":>11s} {"accepted":>11s}')
for name in sorted(coupling.inlets, key=lambda n: int(n[3:])):
    df = logger.to_dataframe(name)
    requested = df.Cum_Captured_m3.iloc[-1] - df.Cum_Surcharged_m3.iloc[-1]
    print(f'{name:6s} {df.Captured_Q_cms.max():6.3f} m3/s {requested:8.1f} m3 {accepted[name]:8.1f} m3')
paths = logger.write_csv(f'hydrographs_{backend}')
print(f'Wrote {len(paths)} hydrograph CSV(s) to hydrographs_{backend}/')

coupling.close()
