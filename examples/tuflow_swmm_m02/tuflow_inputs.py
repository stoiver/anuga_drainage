"""Translate the TUFLOW-SWMM Module 2 inputs into what an ANUGA run needs.

Each function reproduces one TUFLOW input. Where ANUGA can't match TUFLOW
exactly, the approximation is stated next to the code:

* ``elevation_grid``  — the 0.5 m DEM, then ``Read GIS Z Shape`` (road crests,
  MAX option) and ``Create TIN Zpts`` (the urban development), as in the .tgc.
  Order matters: TUFLOW applies them top to bottom, later layers win.
* ``material_ids`` / ``read_materials`` — ``Set Mat`` + ``Read GIS Mat`` layers,
  and the materials file (Manning's n, rainfall initial/continuing losses).
* ``net_rain_rate`` — a TUFLOW rainfall hyetograph (mm fallen between the
  previous and current time, applied as a STEPPED rate), minus the losses.
* ``hydrograph`` — a bc_dbase flow time series (linearly interpolated).
* ``pit_inlets`` — which SWMM junctions have an inlet (from the inlet-usage
  layer) and which are sealed manholes.
* ``write_swmm_inp`` — the pipe-network .inp made runnable on its own by pyswmm.

GIS dependencies (examples only, not the package): geopandas, rasterio,
shapely >= 2, scipy.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent / "tuflow_data" / "TUFLOW"
MODEL = ROOT / "model"
GPKG = MODEL / "gis" / "TS02_001.gpkg"
SWMM_INP = MODEL / "swmm" / "sw02_001.inp"            # the pipe network
INLET_USAGE = MODEL / "swmm" / "swmm_iu_TS02_001.gpkg"
MATERIALS = MODEL / "materials_TS02_001.csv"
BC_DBASE = ROOT / "bc_dbase"

# Layer names, in the order the .tgc / .tbc read them.
CODE_LAYER = "2d_code_M01_001_R"
ZSHAPE_LAYERS = ("2d_zsh_M02_rd_crest_001_L", "2d_zsh_M02_rd_crest_001_P")
TIN_LAYERS = ("2d_ztin_TS02_Dev_001_R", "2d_ztin_TS02_Dev_001_L", "2d_ztin_TS02_Dev_001_P")
MAT_LAYERS = ("2d_mat_M01_001_R", "2d_mat_TS02_Dev_Lots_001_R",
              "2d_mat_TS02_Dev_Buildings_001_R")
BC_LAYER = "2d_bc_M01_001_L"
TUFLOW_CELL = 5.0                                       # .tgc Cell Size


def _layer(name, path=GPKG):
    import geopandas as gpd
    return gpd.read_file(path, layer=name)


def _single_polygon(geom):
    if geom.geom_type == "MultiPolygon":
        if len(geom.geoms) != 1:
            raise ValueError("expected a single polygon")
        geom = geom.geoms[0]
    return geom


# --------------------------------------------------------------------------
# Domain extent and boundaries
# --------------------------------------------------------------------------
def active_polygon():
    """The 2d_code (Code 1) polygon: the active 2D area, as [[x, y], ...]."""
    code = _layer(CODE_LAYER)
    poly = _single_polygon(code.geometry.iloc[0])
    return [list(c[:2]) for c in poly.exterior.coords[:-1]]


def boundary_lines():
    """The 2d_bc lines: {'QT': (name, [[x, y], ...]), 'HQ': (slope, [[x, y], ...])}."""
    bc = _layer(BC_LAYER)
    out = {}
    for _, r in bc.iterrows():
        xy = [list(c[:2]) for c in r.geometry.coords]
        out[r.Type] = (r.Name if r.Type == "QT" else float(r.b), xy)
    return out


def boundary_tags(polygon, line, tag, other="wall", tol=2.0):
    """ANUGA boundary_tags for ``polygon``: segments whose midpoint lies within
    ``tol`` of ``line`` get ``tag``, the rest ``other``."""
    from shapely.geometry import LineString, Point
    ls = LineString(line)
    n = len(polygon)
    tags = {tag: [], other: []}
    for i in range(n):
        (x0, y0), (x1, y1) = polygon[i], polygon[(i + 1) % n]
        mid = Point((x0 + x1) / 2, (y0 + y1) / 2)
        tags[tag if ls.distance(mid) < tol else other].append(i)
    if not tags[tag]:
        raise ValueError(f"no polygon segment lies on the {tag!r} line")
    return tags


# --------------------------------------------------------------------------
# Elevation: DEM + Z-shape road crests + development TIN
# --------------------------------------------------------------------------
def _vertex_z(line, points):
    """Z at each vertex of ``line`` from the points sitting on it (TUFLOW puts
    the elevations on _P points at line vertices); vertices without a point
    are interpolated along the chainage."""
    import shapely
    xy = np.asarray(line.coords)[:, :2]
    pts = shapely.points(np.c_[points.geometry.x, points.geometry.y])
    z = np.full(len(xy), np.nan)
    for k, v in enumerate(shapely.points(xy)):
        d = shapely.distance(pts, v)
        j = int(np.argmin(d))
        if d[j] < 0.05 and np.isfinite(points.Z.iloc[j]):
            z[k] = points.Z.iloc[j]
    chain = np.r_[0.0, np.cumsum(np.hypot(*np.diff(xy, axis=0).T))]
    ok = np.isfinite(z)
    if ok.any():
        z = np.interp(chain, chain[ok], z[ok])
    return xy, chain, z


def _lines(gdf):
    for geom in gdf.geometry:
        yield from (geom.geoms if geom.geom_type.startswith("Multi") else [geom])


def _apply_zshape(z, X, Y):
    """Read GIS Z Shape with the MAX option: raise the ground along each line,
    over the line's Shape_Width (or one TUFLOW cell if unset), to the elevation
    interpolated between the vertex points."""
    import shapely
    lines, points = (_layer(n) for n in ZSHAPE_LAYERS)
    for (_, r), line in zip(lines.iterrows(), _lines(lines)):
        width = r.Shape_Widt if np.isfinite(r.Shape_Widt) else TUFLOW_CELL
        _, chain, vz = _vertex_z(line, points)
        x0, y0, x1, y1 = line.buffer(width / 2).bounds
        win = (X >= x0) & (X <= x1) & (Y >= y0) & (Y <= y1)
        pts = shapely.points(X[win], Y[win])
        near = shapely.distance(line, pts) <= width / 2
        s = shapely.line_locate_point(line, pts[near])
        idx = np.flatnonzero(win)[near]
        z.flat[idx] = np.maximum(z.flat[idx], np.interp(s, chain, vz))


def _apply_tin(z, X, Y, sample, spacing=1.0):
    """Create TIN Zpts: inside the TIN polygon, replace the ground with a TIN
    through the _P points and the _L breaklines (densified to ``spacing`` so a
    plain Delaunay honours them), tied to the existing ground along the polygon
    edge (the polygon has no points of its own, so TUFLOW samples the zpts).
    The layers' dZ attribute is not reproduced."""
    import shapely
    from scipy.interpolate import LinearNDInterpolator
    region, lines, points = (_layer(n) for n in TIN_LAYERS)
    poly = _single_polygon(region.geometry.iloc[0])
    nodes = [np.c_[points.geometry.x, points.geometry.y, points.Z]]
    for line in _lines(lines):
        xy, chain, vz = _vertex_z(line, points)
        if not np.isfinite(vz).any():
            continue
        s = np.linspace(0, chain[-1], max(2, int(chain[-1] / spacing) + 1))
        pts = shapely.line_interpolate_point(line, s)
        nodes.append(np.c_[shapely.get_x(pts), shapely.get_y(pts), np.interp(s, chain, vz)])
    ring = poly.exterior
    s = np.linspace(0, ring.length, int(ring.length / spacing) + 1)
    pts = shapely.line_interpolate_point(ring, s)
    ex, ey = shapely.get_x(pts), shapely.get_y(pts)
    nodes.append(np.c_[ex, ey, sample(ex, ey)])
    nodes = np.vstack(nodes)
    tin = LinearNDInterpolator(nodes[:, :2], nodes[:, 2])
    inside = shapely.contains_xy(poly, X, Y)
    zt = tin(X[inside], Y[inside])
    ok = np.isfinite(zt)
    z.flat[np.flatnonzero(inside)[ok]] = zt[ok]


def elevation_grid():
    """The modified 0.5 m elevation raster: returns (x_centres, y_centres
    ascending, z[ny, nx])."""
    import rasterio
    from scipy.interpolate import RegularGridInterpolator
    with rasterio.open(MODEL / "grid" / "DEM.tif") as r:
        z = r.read(1).astype(float)
        z[z == r.nodata] = np.nan
        t = r.transform
    ny, nx = z.shape
    x = t.c + t.a * (np.arange(nx) + 0.5)
    y = t.f + t.e * (np.arange(ny) + 0.5)          # descending (north-up raster)
    z, y = z[::-1], y[::-1]                         # make y ascending
    X, Y = np.meshgrid(x, y)
    dem = RegularGridInterpolator((y, x), z.copy(), bounds_error=False)
    _apply_zshape(z, X, Y)
    _apply_tin(z, X, Y, lambda px, py: dem(np.c_[py, px]))
    return x, y, z


def elevation_function(x, y, z):
    """ANUGA quantity function ``f(px, py)`` sampling the grid (bilinear).
    ``px``/``py`` are absolute coordinates."""
    from scipy.interpolate import RegularGridInterpolator
    interp = RegularGridInterpolator((y, x), z, bounds_error=False, fill_value=None)
    return lambda px, py: interp(np.c_[py, px])


# --------------------------------------------------------------------------
# Materials: Manning's n and rainfall losses
# --------------------------------------------------------------------------
def read_materials(path=MATERIALS):
    """Materials file -> DataFrame indexed by Material ID with columns
    ``n`` (Manning's n), ``IL`` (mm), ``CL`` (mm/h), ``description``.

    TUFLOW gives some materials a depth-varying n as ``"y1,n1,y2,n2"``: n1 up
    to depth y1, n2 beyond y2, linear between. ``n`` is the deep-water value
    n2; the ``y1``/``n1``/``y2``/``n2`` columns keep the whole curve (NaN for
    a constant n) for :func:`depth_varying_n`."""
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    rows = []
    for _, r in df.iterrows():
        n = [float(v) for v in str(r["Manning's n"]).split(",")]
        losses = [float(v) for v in str(r["Rainfall Losses"]).split(",")]
        curve = n if len(n) == 4 else [np.nan] * 4
        rows.append({"id": int(r["Material ID"]),
                     "n": n[3] if len(n) == 4 else n[0],
                     "y1": curve[0], "n1": curve[1], "y2": curve[2], "n2": curve[3],
                     "IL": losses[0], "CL": losses[1],
                     "description": str(r["! Description"]).lstrip("!")})
    return pd.DataFrame(rows).set_index("id")


def depth_varying_n(mats, mat):
    """For triangles with material ids ``mat``: ``(idx, n_of_depth)``, where
    ``idx`` are the triangles whose material has a depth-varying n and
    ``n_of_depth(h)`` gives their n at depths ``h`` (TUFLOW's y1,n1,y2,n2 rule)."""
    curve = mats.loc[mat, ["y1", "n1", "y2", "n2"]].to_numpy()
    idx = np.flatnonzero(np.isfinite(curve[:, 0]))
    y1, n1, y2, n2 = curve[idx].T

    def n_of_depth(h):
        f = np.clip((h - y1) / (y2 - y1), 0.0, 1.0)
        return n1 + f * (n2 - n1)
    return idx, n_of_depth


def material_ids(px, py):
    """Material ID at each point: ``Set Mat == 1`` then each ``Read GIS Mat``
    layer in .tgc order (later layers overwrite)."""
    import shapely
    mat = np.ones(len(px), dtype=int)
    for name in MAT_LAYERS:
        for _, r in _layer(name).iterrows():
            mat[shapely.contains_xy(r.geometry, px, py)] = int(r.Material)
    return mat


# --------------------------------------------------------------------------
# Boundary time series
# --------------------------------------------------------------------------
def _bc_series(name):
    db = pd.read_csv(BC_DBASE / "bc_dbase_TS02_001.csv")
    db.columns = [c.strip() for c in db.columns]
    row = db.set_index("Name").loc[name]
    ts = pd.read_csv(BC_DBASE / row["Source"]).dropna(how="all")
    return ts[row["Column 1"]].to_numpy(float), ts[row["Column 2"]].to_numpy(float)


def hydrograph(name):
    """A bc_dbase flow series as ``Q(t_seconds)`` [m^3/s], linearly
    interpolated (TUFLOW's default for flow boundaries), 0 outside."""
    t_h, q = _bc_series(name)
    return lambda t: float(np.interp(t / 3600.0, t_h, q, left=0.0, right=0.0))


def net_rain_rate(times_h, depths_mm, IL=0.0, CL=0.0):
    """Net rainfall rate ``r(t_seconds)`` [m/s] from a TUFLOW hyetograph.

    TUFLOW manual 8.4.3.1.1: each value is the depth (mm) that fell between the
    previous and the current time, applied as a constant (STEPPED) rate over
    that interval. Losses are taken off first: the initial loss ``IL`` (mm)
    absorbs the first rain, then the continuing loss ``CL`` (mm/h) is removed
    from each interval's rate (never below zero). Pure, so it is unit-testable.
    """
    times_h = np.asarray(times_h, float)
    depths = np.asarray(depths_mm, float)
    dt_h = np.diff(times_h)
    rain = depths[1:]                         # rain over (t[k-1], t[k]]
    il_left = IL
    net = np.zeros_like(rain)
    for k, p in enumerate(rain):
        taken = min(p, il_left)
        il_left -= taken
        p -= taken
        if il_left <= 0.0:
            p = max(0.0, p - CL * dt_h[k])
        net[k] = p
    rate = net / dt_h / 1000.0 / 3600.0       # mm per interval -> m/s
    edges = times_h * 3600.0

    def r(t):
        k = np.searchsorted(edges, t, side="right") - 1
        return float(rate[k]) if 0 <= k < len(rate) else 0.0
    return r


def rainfall(name):
    """The (times_h, depths_mm) hyetograph of a bc_dbase rainfall entry."""
    return _bc_series(name)


# --------------------------------------------------------------------------
# SWMM pipe network
# --------------------------------------------------------------------------
def pit_inlets():
    """Inlet-usage points -> {junction: row}, matched to the nearest junction.
    Junctions with no inlet are sealed manholes (pipe-only)."""
    import geopandas as gpd
    from anuga_drainage.inlet_initialization import read_inp_coordinates
    iu = gpd.read_file(INLET_USAGE)
    c = read_inp_coordinates(SWMM_INP)
    out = {}
    for _, r in iu.iterrows():
        d = np.hypot(c.X_Coord - r.geometry.x, c.Y_Coord - r.geometry.y)
        name = d.idxmin()
        if d.min() > 1.0:
            raise ValueError(f"inlet at {r.geometry} is {d.min():.1f} m from any junction")
        out[name] = r
    return out


# --------------------------------------------------------------------------
# TUFLOW 2D grid and the cells each pit is connected to
# --------------------------------------------------------------------------
GRID_ORIGIN = np.array([292725.0, 6177615.0])            # .tgc Origin
GRID_XAXIS = np.array([293580.0, 6177415.0]) - GRID_ORIGIN  # .tgc Orientation
INLET_LINKS = "2d_bc_SWMM_Inlet_Connections_001_L"
# Pit9's two nearest cell centres are 0.93 and 0.94 cells away: a near-tie
# that TUFLOW resolves the other way (from its check layer, 1d_to_2d_check_R).
CELL_OVERRIDES = {"Pit9": [(57, 111), (57, 112)]}


def _grid_axes():
    ex = GRID_XAXIS / np.hypot(*GRID_XAXIS)
    return ex, np.array([-ex[1], ex[0]])


def grid_cell(i, j):
    """Corner polygon [[x, y], ...] of TUFLOW cell (i, j)."""
    ex, ey = _grid_axes()
    d = TUFLOW_CELL
    return [list(GRID_ORIGIN + (i + a) * d * ex + (j + b) * d * ey)
            for a, b in ((0, 0), (1, 0), (1, 1), (0, 1))]


def _grid_ij(x, y):
    """Fractional (i, j) grid coordinates of points."""
    ex, ey = _grid_axes()
    r = np.c_[np.atleast_1d(x), np.atleast_1d(y)] - GRID_ORIGIN
    return r @ ex / TUFLOW_CELL, r @ ey / TUFLOW_CELL


def _crosshair_cells(fi, fj):
    """Cells whose centre cross-hairs a polyline (in grid coordinates)
    crosses: TUFLOW's cell selection for SX lines."""
    cells = set()
    for (i0, j0), (i1, j1) in zip(zip(fi[:-1], fj[:-1]), zip(fi[1:], fj[1:])):
        for a0, a1, b0, b1, swap in ((i0, i1, j0, j1, False), (j0, j1, i0, i1, True)):
            lo, hi = sorted((a0, a1))
            for k in range(int(np.ceil(lo - 0.5)), int(np.floor(hi - 0.5)) + 1):
                t = (k + 0.5 - a0) / (a1 - a0)        # where the segment crosses k + 0.5
                other = int(np.floor(b0 + t * (b1 - b0)))
                cells.add((other, k) if swap else (k, other))
    return sorted(cells)


def pit_cells(inlets=None):
    """{pit: [(i, j), ...]}: the TUFLOW cells each inlet pit is connected to,
    as TUFLOW does it. A pit with a 2d_bc SX line (reached from the pit by a
    CN line) gets the cells whose centre cross-hairs that line crosses;
    ``Conn_width`` -n gets the n cells with centres nearest the pit; otherwise
    the cell containing the pit (SXZ). With ``CELL_OVERRIDES`` for Pit9 (a
    near-tie), this matches TUFLOW's check layer for every pit."""
    import geopandas as gpd
    inlets = inlets if inlets is not None else pit_inlets()
    links = gpd.read_file(GPKG, layer=INLET_LINKS)
    sx = [g for t, g in zip(links.Type, links.geometry) if t == "SX"]
    cn = [g for t, g in zip(links.Type, links.geometry) if t == "CN"]
    out = {}
    for name, r in inlets.items():
        p = r.geometry
        line = None
        for c in cn:
            if c.boundary.distance(p) < 0.01:          # CN line starting at this pit
                line = min(sx, key=lambda g: g.distance(c))
        if line is not None:
            cells = _crosshair_cells(*_grid_ij(*np.array(line.coords).T))
        else:
            fi, fj = _grid_ij(p.x, p.y)
            n = int(-r.Conn_width) if np.isfinite(r.Conn_width) and r.Conn_width < 0 else 1
            i0, j0 = int(fi[0]), int(fj[0])
            cand = [(i, j) for i in range(i0 - 2, i0 + 3) for j in range(j0 - 2, j0 + 3)]
            centre = lambda ij: np.hypot(ij[0] + 0.5 - fi[0], ij[1] + 0.5 - fj[0])
            cells = sorted(sorted(cand, key=centre)[:n])
        out[name] = CELL_OVERRIDES.get(name, cells)
    return out


def cells_polygon(cells):
    """One polygon [[x, y], ...] covering adjacent TUFLOW cells."""
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    u = unary_union([Polygon(grid_cell(*ij)) for ij in cells]).simplify(1e-6)
    if u.geom_type != "Polygon":
        raise ValueError(f"cells {cells} are not one connected block")
    return [list(c) for c in u.exterior.coords[:-1]]


CULVERT_INP = MODEL / "swmm" / "sw01_001.inp"
CULVERT_LINKS = "2d_bc_SWMM_Culvert_Connections_001_L"


def culverts():
    """The three road culverts of sw01_001.inp as ANUGA Boyd-culvert parameters.

    Each is a SWMM conduit from a [STORAGE] node, which TUFLOW connects to the
    2D upstream with an HX line, to an [OUTFALLS] node, which it connects
    downstream with an SX line. Returns a list of dicts: ``name``, ``shape``
    (RECT_CLOSED/CIRCULAR), ``height``, ``width``, ``barrels``, ``length``,
    ``n``, ``losses`` (entry + exit), ``inverts`` (up, down, incl. offsets)
    and ``exchange_lines`` ([HX line, SX line])."""
    import shapely
    from anuga_drainage.inp import _read_sections
    sec = _read_sections(CULVERT_INP)
    elev = {r[0]: float(r[1]) for r in sec.get("STORAGE", []) + sec.get("OUTFALLS", [])}
    xy = {r[0]: (float(r[1]), float(r[2])) for r in sec["COORDINATES"]}
    xs = {r[0]: r for r in sec["XSECTIONS"]}
    loss = {r[0]: float(r[1]) + float(r[2]) for r in sec.get("LOSSES", [])}
    links = _layer(CULVERT_LINKS)
    lines = {t: [g for g in links.geometry[links.Type == t]] for t in ("HX", "SX")}

    def nearest(kind, node):
        p = shapely.Point(xy[node])
        g = min(lines[kind], key=p.distance)
        if p.distance(g) > 10.0:
            raise ValueError(f"no {kind} line near culvert node {node}")
        return [list(c[:2]) for c in g.coords]

    out = []
    for name, up, down, length, n, in_off, out_off, *_ in sec["CONDUITS"]:
        shape, g1, g2, *rest = xs[name][1:]
        out.append({"name": name, "shape": shape.upper(), "height": float(g1),
                    "width": float(g2) if shape.upper() == "RECT_CLOSED" else float(g1),
                    "barrels": int(rest[2]) if len(rest) > 2 else 1,
                    "length": float(length), "n": float(n), "losses": loss.get(name, 0.0),
                    "inverts": [elev[up] + float(in_off), elev[down] + float(out_off)],
                    "exchange_lines": [nearest("HX", up), nearest("SX", down)]})
    return out


def _last_inlet_lines(path=SWMM_INP):
    """{design: type of its last [INLETS] line} (GRATE, CURB, ...)."""
    last, section = {}, None
    for line in Path(path).read_text().splitlines():
        line = line.split(";", 1)[0].strip()
        if line.startswith("["):
            section = line.upper()
        elif line and section == "[INLETS]":
            name, kind = line.split()[:2]
            last[name] = kind.upper()
    return last


def inlet_captures(inlets=None, last_line_wins=False):
    """{junction: InletCapture} for each inlet pit: the SWMM/HEC-22 capture
    capacity TUFLOW-SWMM applies, built from sw02's [INLETS]/[STREETS] and the
    inlet-usage layer (placement, longitudinal slope, number, clogging, flow
    limit, local depression).

    By default a design given on several [INLETS] lines is a combination
    inlet, as in EPA SWMM. ``last_line_wins=True`` reproduces TUFLOW-SWMM
    2026.0, which keeps only a design's last line. Module 2's ``Combo1``
    (GRATE then CURB) then captures as a curb opening alone. This was
    verified by rerunning TUFLOW with the two lines swapped: the sag pits
    then match grate-only capture."""
    import dataclasses
    from anuga_drainage import InletCapture, read_inlet_designs
    designs, streets = read_inlet_designs(SWMM_INP)
    if last_line_wins:
        for name, kind in _last_inlet_lines().items():
            if kind == "CURB":
                designs[name] = dataclasses.replace(designs[name], grate=None)
            elif kind == "GRATE":
                designs[name] = dataclasses.replace(designs[name], curb=None)
    inlets = inlets if inlets is not None else pit_inlets()
    num = lambda v, default=0.0: float(v) if v is not None and np.isfinite(v) else default
    return {name: InletCapture(designs[r.Inlet], streets[r.StreetXSEC],
                               placement=r.Placement,
                               long_slope=num(r.SlopePct_Longitudinal) / 100.0,
                               number=int(num(r.Number, 1)),
                               clogged_pct=num(r.CloggedPct), q_max=num(r.Qmax),
                               local_depression=num(r.aLocal), local_width=num(r.wLocal))
            for name, r in inlets.items()}


# TUFLOW's options from the culvert .inp (sw01; TUFLOW merges [OPTIONS] across
# the .inp files) with the timing and ponding changed for pyswmm + our coupling:
# whole-second routing (pyswmm 2.1 strides in int seconds), and no ponding, so
# surcharge floods back to the 2D surface through the coupling. Two more
# changes cut SWMM's own flow-routing continuity error (standalone runs with
# 12 pits fed 0.01-0.5 m^3/s): TUFLOW's MIN_SURFAREA 25 gave 4-13 %, SWMM's
# default 1.167 with SKIP_STEADY_STATE NO gives ~1 % once flows build up.
_OPTIONS = """\
[TITLE]
TUFLOW-SWMM tutorial Module 2 pipe network (sw02_001.inp), made standalone for anuga_drainage

[OPTIONS]
FLOW_UNITS           CMS
FLOW_ROUTING         DYNWAVE
LINK_OFFSETS         DEPTH
FORCE_MAIN_EQUATION  H-W
IGNORE_RAINFALL      YES
IGNORE_SNOWMELT      YES
IGNORE_GROUNDWATER   YES
IGNORE_RDII          YES
IGNORE_QUALITY       YES
ALLOW_PONDING        NO
SKIP_STEADY_STATE    NO
START_DATE           01/01/2020
START_TIME           00:00:00
REPORT_START_DATE    01/01/2020
REPORT_START_TIME    00:00:00
END_DATE             01/01/2020
END_TIME             {end}
REPORT_STEP          00:01:00
WET_STEP             00:00:30
DRY_STEP             00:30:00
ROUTING_STEP         1
VARIABLE_STEP        0
MINIMUM_STEP         0.5
INERTIAL_DAMPING     PARTIAL
NORMAL_FLOW_LIMITED  BOTH
MIN_SURFAREA         1.167
MAX_TRIALS           20
HEAD_TOLERANCE       0.001
THREADS              1

"""

_DROP = {"[STREETS]", "[INLETS]"}      # SWMM's inlet model; our coupling replaces it


def write_swmm_inp(out_path, finaltime_s, rims=None):
    """Write a pyswmm-runnable copy of sw02_001.inp: add [OPTIONS] (it has
    none — TUFLOW takes them from sw01) and drop [STREETS]/[INLETS].

    ``rims`` ({junction: rim elevation}) sets those junctions' Ymax to rim -
    invert. The tutorial leaves Ymax = 0, which SWMM reads as "up to the highest
    pipe crown": a surcharged pipe then floods the pit (and the coupling puts
    that water on the surface) while the HGL is still metres below ground. With
    the rim as Ymax, a pit overflows only once its head reaches the surface."""
    rims = rims or {}
    secs = int(np.ceil(finaltime_s)) + 60
    end = f"{secs // 3600:02d}:{secs % 3600 // 60:02d}:{secs % 60:02d}"
    keep, section, lines = True, None, []
    for line in SWMM_INP.read_text().splitlines():
        s = line.strip()
        if s.startswith("["):
            section = s.upper()
            keep = section not in _DROP
        elif section == "[JUNCTIONS]" and s and not s.startswith(";"):
            f = s.split()
            if f[0] in rims:
                f[2] = f"{rims[f[0]] - float(f[1]):.3f}"
                line = "  ".join(f)
        if keep:
            lines.append(line)
    text = _OPTIONS.format(end=end) + "\n".join(lines) + "\n\n[REPORT]\nNODES ALL\nLINKS ALL\n"
    Path(out_path).write_text(text)
    return str(out_path)
