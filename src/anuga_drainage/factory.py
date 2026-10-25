"""One-call coupling setup from a single SWMM ``.inp``, for either backend.

``couple_from_inp`` parses a ``.inp``, builds the 1D backend (a SWMM
``Simulation`` or a pipedream ``SuperLink`` converted from the ``.inp``),
auto-creates the ANUGA inlet operators at each junction, and returns a ready
``Coupler``. So a coupled model becomes "write one ``.inp``, pick a backend,
run the evolve loop". The junctions are coupled to the surface; outfalls are
treated as boundaries (free drainage for pipedream; SWMM handles its own).
"""
from dataclasses import dataclass, field

import numpy as np

from .inp import read_inp, inp_to_pipedream
from .inlet_initialization import n_sided_inlet
from .coupler import Coupler, SwmmBackend, PipedreamBackend
from .hydrograph import HydrographLogger


@dataclass
class Coupling:
    """Result of :func:`couple_from_inp` — a ready coupling you drive directly.

    Drive the loop with :meth:`step`, optionally attach a volume audit with
    :meth:`add_volume_balance`, and release the backend with :meth:`close`::

        coupling = couple_from_inp(domain, 'net.inp', backend='pipedream')
        coupling.add_volume_balance(inflow_operators=[my_inflow_op])
        for t in domain.evolve(yieldstep=dt, finaltime=ft):
            coupling.step(dt)
        print(coupling.volume_balance.summary())
        coupling.close()

    The components are also exposed directly: ``coupler``, ``inlets``
    (name → ANUGA ``Inlet_operator``), ``backend``, ``handle`` (the pyswmm
    ``Simulation`` / pipedream ``SuperLink``), ``inp`` and ``domain``.
    """
    coupler: object
    inlets: dict          # junction name -> ANUGA Inlet_operator
    backend: object       # SwmmBackend / PipedreamBackend
    handle: object        # pyswmm Simulation (swmm) or pipedream SuperLink
    inp: object           # parsed InpNetwork
    domain: object        # the ANUGA domain
    volume_balance: object = None
    _prev_step: object = field(default=None, init=False, repr=False)

    def step(self, dt):
        """Run one coupled exchange step; if a VolumeBalance is attached, record
        it first (at the loop top, with the previous step, so the reads align)."""
        if self.volume_balance is not None:
            self.volume_balance.step(self.domain.get_time(), dt, self._prev_step)
        self._prev_step = self.coupler.step(dt)
        return self._prev_step

    def add_volume_balance(self, inflow_operators=(), outfall_inlet=None):
        """Attach a :class:`~anuga_drainage.VolumeBalance`; subsequent
        :meth:`step` calls update it. Returns the VolumeBalance."""
        from .volume_balance import VolumeBalance
        self.volume_balance = VolumeBalance(
            self.domain, list(self.inlets.values()), self.backend,
            inflow_operators=inflow_operators, outfall_inlet=outfall_inlet)
        return self.volume_balance

    def close(self):
        """Release backend resources (closes the SWMM simulation; no-op for
        pipedream)."""
        self.backend.close()


def _as_array(x, n):
    a = np.atleast_1d(np.asarray(x, dtype=float))
    return np.full(n, a[0]) if a.size == 1 else a


def _polygon_area(p):
    """Shoelace area of a polygon given as a list of ``[x, y]`` vertices."""
    a = np.asarray(p, dtype=float)
    x, y = a[:, 0], a[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def _polygon_perimeter(p):
    """Perimeter of a closed polygon given as a list of ``[x, y]`` vertices."""
    a = np.asarray(p, dtype=float)
    d = a - np.roll(a, 1, axis=0)
    return float(np.sum(np.hypot(d[:, 0], d[:, 1])))


def couple_from_inp(domain, inp_path, backend="swmm", *,
                    manhole_area=1.0, n_sides=6, rotation=0.0, inlet_polygons=None,
                    inlet_specs=None, library=None, blockage=0.0,
                    time_average=1.0, clamp=True, cw=0.67, co=0.67,
                    internal_links=20, pit_area=1.0, pipedream_max_step=None,
                    superlink_kwargs=None, log_hydrographs=False,
                    coupled_junctions=None):
    """Build a ready :class:`~anuga_drainage.Coupler` from a SWMM ``.inp``.

    Parameters
    ----------
    domain : the ANUGA domain (meshed, elevation set).
    inp_path : path to the SWMM ``.inp`` describing the sewer network.
    backend : ``"swmm"`` (pyswmm) or ``"pipedream"``.
    manhole_area : surface area of each inlet coupling region (scalar or one per
        junction); also used as the pipedream superjunction storage area.
    n_sides, rotation : geometry of the regular-polygon inlet regions (used for
        any junction not overridden by ``inlet_polygons``).
    inlet_polygons : optional ``{junction_name: [[x, y], ...]}`` to give specific
        inlets an explicit footprint instead of the auto regular polygon — needed
        when an inlet must **span a channel** (the ``.inp`` only carries a point,
        so the auto polygon can be too narrow and flow overtops past it). The
        coupling region's area and perimeter become that junction's
        ``manhole_area`` and ``weir_length`` (unless overridden by ``inlet_specs``).
    inlet_specs : optional ``{junction_name: spec}`` assigning a named inlet
        (catalogue key) or an :class:`~anuga_drainage.InletSpec` to a junction.
        The spec's ``operational_area`` / ``operational_perimeter`` then drive the
        weir/orifice flux (``area_manhole`` / ``length_weir`` in ``calculate_Q``)
        for that junction, **decoupled** from the surface coupling footprint:
        the ANUGA region (and pipedream storage) still come from ``manhole_area`` /
        ``inlet_polygons``, so a small grate opening doesn't shrink the footprint.
        Junctions without a spec keep the footprint-derived geometry (unchanged).
    library : optional ``{name: InletSpec}`` catalogue for resolving ``inlet_specs``
        keys (default: the built-in ``INLET_LIBRARY``).
    blockage : clogging fraction 0.0..1.0 applied to spec'd junctions; a scalar
        (all) or a ``{junction_name: fraction}`` dict. Derates the spec's area and
        perimeter. Ignored for junctions without an ``inlet_specs`` entry.
    time_average, clamp, cw, co : forwarded to the ``Coupler``.
    coupled_junctions : optional list of ``[JUNCTIONS]`` names to couple to the
        surface (default: all). The others stay in the 1D network but get no
        ANUGA inlet — e.g. sealed manholes that only join pipes. Per-junction
        ``manhole_area`` arrays are still given for *all* junctions, in
        ``[JUNCTIONS]`` order; ``inlet_polygons`` / ``inlet_specs`` /
        ``blockage`` keys must be coupled junctions.
    log_hydrographs : if True, attach a :class:`~anuga_drainage.HydrographLogger`
        that records a per-inlet hydrograph each step; access it via
        ``coupling.coupler.logger`` and dump CSVs with ``logger.write_csv(dir)``.
    internal_links, pit_area, superlink_kwargs : pipedream-only (discretisation,
        internal-junction storage, extra ``SuperLink`` kwargs).
    pipedream_max_step : pipedream-only cap on the solver's *internal* hydraulic
        timestep (s). The coupling ``dt`` (ANUGA yieldstep / exchange frequency)
        can stay coarse — e.g. 1 s — while each pipedream step is subdivided into
        ``ceil(dt/pipedream_max_step)`` sub-steps, which the semi-implicit solver
        needs for stability. ``None`` steps once at ``dt`` (was unstable at 1 s).
        The more ``internal_links``, the shorter each sub-conduit, so the smaller
        this must be (CFL): the default 20 links needs a finer step than the
        hand-built run_pipedream.py's 6 links @ 0.05 s.

    Returns
    -------
    Coupling
        A ready coupling (see :class:`Coupling`).
    """
    from anuga import Inlet_operator, Region   # lazy: pure callers don't need ANUGA
    from .inlet_catalogue import resolve_inlet_spec

    inp = read_inp(inp_path)
    all_jnames = list(inp.junctions["name"])
    if not all_jnames:
        raise ValueError(f"{inp_path}: no [JUNCTIONS] to couple")
    all_areas = _as_array(manhole_area, len(all_jnames))
    if coupled_junctions is None:
        jnames = all_jnames
    else:
        unknown = set(coupled_junctions) - set(all_jnames)
        if unknown:
            raise ValueError(f"coupled_junctions names not in [JUNCTIONS]: {sorted(unknown)}")
        jnames = [name for name in all_jnames if name in set(coupled_junctions)]
        if not jnames:
            raise ValueError("coupled_junctions is empty: nothing to couple")
    coords = inp.coordinates.set_index("node")
    areas_in = [all_areas[all_jnames.index(name)] for name in jnames]
    inlet_polygons = inlet_polygons or {}
    unknown = set(inlet_polygons) - set(jnames)
    if unknown:
        raise ValueError(f"inlet_polygons names not in [JUNCTIONS]: {sorted(unknown)}")

    # Resolve any inlet_specs to derated InletSpecs keyed by junction name.
    inlet_specs = inlet_specs or {}
    unknown = set(inlet_specs) - set(jnames)
    if unknown:
        raise ValueError(f"inlet_specs names not in [JUNCTIONS]: {sorted(unknown)}")
    specs = {name: resolve_inlet_spec(
                ref, library,
                blockage[name] if isinstance(blockage, dict) else blockage)
             for name, ref in inlet_specs.items()}

    # --- ANUGA inlet operators at each junction (backend-agnostic) ---
    # The polygon sets the surface coupling footprint (the ANUGA region, and the
    # pipedream storage area). The *hydraulic* area/perimeter fed to calculate_Q
    # come from an assigned inlet_spec if any, else the footprint geometry — so a
    # small grate opening drives the flux without shrinking the footprint.
    inlets, beds = [], []
    footprint_areas, hyd_weirs, hyd_areas = [], [], []
    for name, area in zip(jnames, areas_in):
        if name in inlet_polygons:
            vertices = [[float(x), float(y)] for x, y in inlet_polygons[name]]
            eff_area = _polygon_area(vertices)
            weir = _polygon_perimeter(vertices)
            # Honour the given footprint exactly: don't expand it across nearby
            # steep terrain (e.g. a channel bank), or returned surcharge gets
            # distributed onto those high cells and strands as a thin film.
            expand = False
        else:
            if name not in coords.index:
                raise ValueError(f"junction {name!r} has no [COORDINATES] entry")
            xy = [float(coords.loc[name, "x"]), float(coords.loc[name, "y"])]
            vertices, side = n_sided_inlet(n_sides, float(area), xy, rotation)
            eff_area = float(area)
            weir = n_sides * side
            expand = True   # a small auto polygon may not contain a cell centroid
        op = Inlet_operator(domain, Region(domain, polygon=vertices, expand_polygon=expand),
                            Q=0.0, zero_velocity=True)
        inlets.append(op)
        beds.append(op.inlet.get_average_elevation())
        footprint_areas.append(eff_area)
        spec = specs.get(name)
        if spec is not None:
            hyd_areas.append(spec.operational_area)
            hyd_weirs.append(spec.operational_perimeter)
        else:
            hyd_areas.append(eff_area)
            hyd_weirs.append(weir)
    beds = np.array(beds)
    footprint_areas = np.array(footprint_areas)
    hyd_weirs = np.array(hyd_weirs)
    hyd_areas = np.array(hyd_areas)

    # --- 1D backend, junctions ordered to match the inlets ---
    if backend == "swmm":
        from pyswmm import Simulation, Nodes
        sim = Simulation(inp_path)
        sim.start()
        nodes = Nodes(sim)
        be = SwmmBackend(sim, junctions=[nodes[name] for name in jnames])
        handle = sim
    elif backend == "pipedream":
        from pipedream_solver.hydraulics import SuperLink
        sj, sl = inp_to_pipedream(inp, manhole_area=float(footprint_areas[0]), pit_area=pit_area)
        superlink = SuperLink(sl, sj, internal_links=internal_links,
                              **(superlink_kwargs or {}))
        n_j = len(all_jnames)
        coupled = [all_jnames.index(name) for name in jnames]  # junctions are listed first
        outfalls = list(range(n_j, n_j + len(inp.outfalls)))  # outfalls follow them
        H_bc = superlink._z_inv_j.copy() if outfalls else None  # free-drain outfalls
        be = PipedreamBackend(superlink, coupled_indices=coupled, H_bc=H_bc,
                              outfall_indices=outfalls, max_step=pipedream_max_step)
        handle = superlink
    else:
        raise ValueError(f"backend must be 'swmm' or 'pipedream', got {backend!r}")

    logger = HydrographLogger(jnames) if log_hydrographs else None
    coupler = Coupler(inlets=inlets, beds=beds, weir_lengths=hyd_weirs,
                      manhole_areas=hyd_areas, backend=be,
                      time_average=time_average, clamp=clamp, cw=cw, co=co,
                      logger=logger)
    return Coupling(coupler=coupler, inlets=dict(zip(jnames, inlets)),
                    backend=be, handle=handle, inp=inp, domain=domain)
