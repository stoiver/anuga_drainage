"""Driver for the 2D (ANUGA) <-> 1D (SWMM / pipedream) coupling loop.

The per-step coupling sequence is the same across the example run scripts:

    1. read the 2D water depth at each inlet,
    2. read the 1D head at each inlet/junction,
    3. compute the exchange flux with calculate_Q (weir/orifice),
    4. smooth it in time to damp oscillations,
    5. optionally clamp the outflow to the water actually available in the 2D cell,
    6. advance the 1D model by dt using that flux,
    7. feed the *realised* 1D flow back into the 2D model via Inlet_operator.set_Q.

This module captures that sequence in `Coupler`, with the 1D-backend-specific
parts (head access, stepping, realised-flow sign/bookkeeping) behind a small
backend interface. The two backends differ in step 7: SWMM reports the flow it
actually accepted via node statistics, whereas pipedream takes the requested
flux as realised and feeds back -Q_in.

Volume-balance / logging / plotting stays in the run scripts (it varies per
example); the backends expose the underlying sim/superlink handles for that.

See CLAUDE.md for the coupling-loop background and the pyswmm 2.1.0 / numpy 2.x
environment constraints.
"""
from collections import namedtuple

import numpy as np

from .coupling import calculate_Q

CouplingStep = namedtuple("CouplingStep", ["Q_in", "anuga_flux"])


def smooth_Q(Q_new, Q_old, dt, time_average):
    """Time-average the coupling flux to damp oscillations.

    Q = ((time_average - dt) * Q_old + dt * Q_new) / time_average

    `time_average <= 0` disables smoothing (returns Q_new unchanged).
    """
    if time_average <= 0:
        return Q_new
    return ((time_average - dt) * Q_old + dt * Q_new) / time_average


def limit_outflow(Q_in, available_volume, dt, safety_factor=1.0):
    """Clamp positive (surface -> pipe) flux so a step cannot remove more water
    than is present in the 2D cell. Negative (surcharge) flux is untouched.
    """
    Q_limit = safety_factor * np.asarray(available_volume) / dt
    return np.where(Q_in > 0, np.minimum(Q_in, Q_limit), Q_in)


class SwmmBackend:
    """Coupling backend for the standard pyswmm release (>= 2.1).

    Heads come from the junction nodes; the model is advanced with whole-second
    strides (pyswmm 2.1.0 requires an int); the flow fed back to ANUGA is the
    flow SWMM actually accepted, derived from node statistics.
    """

    def __init__(self, sim, junctions=None, links=None, outfalls=None):
        from pyswmm import Nodes, Links

        self.sim = sim
        if junctions is None:
            junctions = [n for n in Nodes(sim) if n.is_junction()]
        self.junctions = list(junctions)
        self.links = list(links) if links is not None else list(Links(sim))
        if outfalls is None:
            outfalls = [n for n in Nodes(sim) if n.is_outfall()]
        self.outfalls = list(outfalls)
        # Every node that can hold water, coupled or not (e.g. a sealed manhole
        # that only joins pipes), so pipe_volume() doesn't depend on which
        # junctions happen to be coupled to the surface.
        self.storage_nodes = [n for n in Nodes(sim) if not n.is_outfall()]
        self._old_vol = np.array([self._inlet_vol(n) for n in self.junctions])
        self._outfall_vol = 0.0   # cumulative volume that left the network at outfalls

    @staticmethod
    def _inlet_vol(node):
        # Net volume that has left the 2D surface at this node: the generated
        # (lateral) inflow that entered the pipe, minus what flooded back out.
        s = node.statistics
        return -s["lateral_infow_vol"] + s["flooding_volume"]

    def get_heads(self):
        return np.array([n.head for n in self.junctions])

    # --- generic 1D-network state (same signature on both backends, so a run
    #     script can record/plot the pipe network without knowing the backend) ---
    def node_depths(self):
        """Water depth above the invert at each coupled node."""
        return np.array([n.depth for n in self.junctions])

    def conduit_names(self):
        """Name of each conduit (parallel to conduit_flows())."""
        return [link.linkid for link in self.links]

    def conduit_flows(self):
        """Flow in each conduit."""
        return np.array([link.flow for link in self.links])

    def step(self, Q_in, dt):
        for node, q in zip(self.junctions, Q_in):
            node.generated_inflow(q)
        self.sim.step_advance(int(dt))  # swmm_stride requires an int (whole seconds)
        next(self.sim)
        # Accumulate the volume leaving at outfalls (read post-step, matching the
        # outfall-return term the scripts add back to ANUGA).
        self._outfall_vol += sum(o.total_inflow for o in self.outfalls) * dt

    def anuga_flux(self, Q_in, dt):
        new = np.array([self._inlet_vol(n) for n in self.junctions])
        flux = (new - self._old_vol) / dt
        self._old_vol = new
        return flux

    def link_volume(self):
        return sum(link.volume for link in self.links)

    # --- independent pipe-side volume accounting (for VolumeBalance) ---
    def pipe_volume(self):
        """Water currently held in the network: conduits + node storage (all
        non-outfall nodes, including junctions not coupled to the surface)."""
        return self.link_volume() + sum(n.volume for n in self.storage_nodes)

    def coupling_inflow_volumes(self):
        """Per-junction cumulative net volume the surface injected (lateral
        inflow accepted minus what flooded back out), from SWMM's statistics."""
        return [n.statistics["lateral_infow_vol"] - n.statistics["flooding_volume"]
                for n in self.junctions]

    def coupling_inflow_volume(self):
        """Cumulative net volume the surface injected at the coupling junctions,
        measured independently of the ANUGA side."""
        return sum(self.coupling_inflow_volumes())

    def outfall_volume(self):
        """Cumulative volume that has left the network at outfalls."""
        return self._outfall_vol

    def close(self):
        """Close the SWMM simulation."""
        self.sim.close()


class PipedreamBackend:
    """Coupling backend for pipedream's SuperLink.

    Heads are the superjunction heads H_j; the requested flux is taken as
    realised, so the flow fed back to ANUGA is -Q_in.
    """

    def __init__(self, superlink, coupled_indices=None, H_bc=None, outfall_indices=None,
                 max_step=None):
        # coupled_indices: which superjunctions exchange with ANUGA (default all,
        #   matching the hand-built examples). couple_from_inp couples only the
        #   junctions and lists the outfalls as boundary (bc) superjunctions.
        # H_bc: fixed boundary heads for bc superjunctions (e.g. outfall inverts
        #   for free drainage); None keeps the previous no-bc behaviour exactly.
        # outfall_indices: bc superjunctions that shed water out of the system,
        #   for outfall_volume tracking.
        # max_step: cap on pipedream's *internal* hydraulic timestep. The coupling
        #   dt (ANUGA yieldstep / exchange frequency) can stay coarse (e.g. 1 s),
        #   but pipedream's semi-implicit solver is only stable at a small step, so
        #   each step(dt) is subdivided into ceil(dt/max_step) sub-steps of the
        #   same Q_in. None (default) steps once at dt, preserving prior behaviour.
        self.superlink = superlink
        self.max_step = max_step
        n = len(superlink.H_j)
        self.coupled = (np.arange(n) if coupled_indices is None
                        else np.asarray(coupled_indices, dtype=int))
        self.H_bc = None if H_bc is None else np.asarray(H_bc, dtype=float)
        self._injected = None   # per-coupled-junction cumulative injected volume
        # Superlink ends discharging at an outfall: Q_dk is the flow into a
        # superlink's downstream superjunction, Q_uk the flow out of its upstream
        # one — so outfall outflow = sum(Q_dk at outfall d/s ends) - sum(Q_uk at
        # outfall u/s ends), accumulated over time.
        outs = set() if outfall_indices is None else {int(i) for i in outfall_indices}
        self._outfall_dk = [k for k in range(len(superlink._J_dk))
                            if int(superlink._J_dk[k]) in outs]
        self._outfall_uk = [k for k in range(len(superlink._J_uk))
                            if int(superlink._J_uk[k]) in outs]
        self._outfall_vol = 0.0

    def get_heads(self):
        return self.superlink.H_j[self.coupled]

    # --- generic 1D-network state (mirrors SwmmBackend, so run scripts can
    #     record/plot the pipe network without knowing the backend) ---
    def node_depths(self):
        """Water depth above the invert at each coupled superjunction (>= 0)."""
        s = self.superlink
        return np.maximum(0.0, np.asarray(s.H_j[self.coupled] - s._z_inv_j[self.coupled]))

    def conduit_names(self):
        """Name of each conduit/superlink (parallel to conduit_flows())."""
        return list(self.superlink.superlinks["name"])

    def conduit_flows(self):
        """Flow in each conduit (the superlink upstream-end flow Q_uk)."""
        return np.asarray(self.superlink.Q_uk)

    def step(self, Q_in, dt):
        q = np.asarray(Q_in, dtype=float)
        full = np.zeros(len(self.superlink.H_j))
        full[self.coupled] = q
        self._injected = q * dt if self._injected is None else self._injected + q * dt
        # Refine pipedream's internal step without changing the exchange frequency:
        # hold Q_in fixed over dt and advance the solver in <= max_step sub-steps.
        nsub = 1 if not self.max_step else max(1, int(np.ceil(dt / self.max_step)))
        sub_dt = dt / nsub
        for _ in range(nsub):
            if self.H_bc is None:
                self.superlink.step(Q_in=full, dt=sub_dt)
            else:
                self.superlink.step(Q_in=full, H_bc=self.H_bc, dt=sub_dt)
            if self._outfall_dk or self._outfall_uk:
                s = self.superlink
                out = (sum(float(s.Q_dk[k]) for k in self._outfall_dk)
                       - sum(float(s.Q_uk[k]) for k in self._outfall_uk))
                self._outfall_vol += out * sub_dt

    def anuga_flux(self, Q_in, dt):
        return -np.asarray(Q_in)

    def link_volume(self):
        s = self.superlink
        return (s._A_ik * s._dx_ik).sum() + (s._A_SIk * s._h_Ik).sum()

    def node_volume(self):
        # Use pipedream's own storage-curve volume: it honours the functional/
        # tabular storage curves and clamps depth at min_depth, matching how
        # pipedream conserves mass. The linear A_sj*(H_j - z_inv_j) goes negative
        # when a superjunction's head drops below its invert (common upstream),
        # which otherwise shows up as a spurious R_pipe.
        return float(self.superlink.compute_storage_volumes().sum())

    def sewer_volume(self):
        return self.link_volume() + self.node_volume()

    # --- independent pipe-side volume accounting (for VolumeBalance) ---
    def pipe_volume(self):
        # Physical water volume (flow-area conduit + pit + superjunction). This
        # is exact, but pipedream's semi-implicit linearised continuity + the
        # Preissmann slot don't conserve *physical* volume to machine precision
        # during redistribution/surcharge, so VolumeBalance's R_pipe shows
        # pipedream's own ~0.1-1% numerical non-conservation (cf. SWMM's larger
        # finite-difference loss). That is the solver's behaviour, not a coupling
        # error or a measurement bug to "correct" away. See CLAUDE.md.
        return self.sewer_volume()

    def coupling_inflow_volumes(self):
        """Per-junction cumulative injected volume (pipedream takes the requested
        flux as realised, so this is the integral of Q_in per junction)."""
        return [] if self._injected is None else list(self._injected)

    def coupling_inflow_volume(self):
        """Cumulative volume injected at the coupling junctions."""
        return 0.0 if self._injected is None else float(self._injected.sum())

    def outfall_volume(self):
        """Cumulative volume shed at outfall (bc) superjunctions (0 if none)."""
        return self._outfall_vol

    def close(self):
        """No external resources to release for pipedream."""


class Coupler:
    """Drives the per-step 2D<->1D exchange for a set of inlets and a backend.

    `inlets` are ANUGA Inlet_operators, ordered to match the backend's heads
    (the SWMM junction order / pipedream superjunction order). `beds`,
    `weir_lengths` and `manhole_areas` are parallel arrays for calculate_Q.

    `capture` optionally caps each inlet's surface->pipe flow at its inlet
    capacity: one entry per inlet, each None (no cap) or a callable
    ``capture(depth, approach_Q) -> m^3/s`` such as
    :class:`~anuga_drainage.hec22.InletCapture`. calculate_Q still sets the
    head-driven exchange (and all surcharge); the cap only limits capture.

    `approach_indices` optionally gives, per inlet, the triangles to measure
    the approaching surface flow in (None = the inlet's own region). Inlet
    operators with ``zero_velocity`` zero the momentum inside their region,
    so an on-grade capture needs a ring of triangles just outside it.
    `approach_widths` optionally gives, per inlet, a flow width [m] (None =
    the default heuristic; see :meth:`approach_flows`).
    """

    def __init__(self, inlets, beds, weir_lengths, manhole_areas, backend,
                 time_average=0.0, clamp=False, safety_factor=1.0,
                 cw=0.67, co=0.67, g=None, logger=None, capture=None,
                 approach_indices=None, approach_widths=None):
        self.inlets = list(inlets)
        self.beds = np.asarray(beds, dtype=float)
        self.weir_lengths = np.asarray(weir_lengths, dtype=float)
        self.manhole_areas = np.asarray(manhole_areas, dtype=float)
        self.backend = backend
        self.time_average = time_average
        self.clamp = clamp
        self.safety_factor = safety_factor
        self.cw = cw
        self.co = co
        self.g = g  # gravity for calculate_Q; None -> ANUGA's value (see calculate_Q)
        self.logger = logger  # optional HydrographLogger; records each step if set
        if capture is not None and len(capture) != len(self.inlets):
            raise ValueError(f"capture has {len(capture)} entries for {len(self.inlets)} inlets")
        self.capture = list(capture) if capture is not None else None
        if approach_indices is not None and len(approach_indices) != len(self.inlets):
            raise ValueError(f"approach_indices has {len(approach_indices)} entries "
                             f"for {len(self.inlets)} inlets")
        self.approach_indices = (list(approach_indices) if approach_indices is not None
                                 else [None] * len(self.inlets))
        if approach_widths is not None and len(approach_widths) != len(self.inlets):
            raise ValueError(f"approach_widths has {len(approach_widths)} entries "
                             f"for {len(self.inlets)} inlets")
        self.approach_widths = (list(approach_widths) if approach_widths is not None
                                else [None] * len(self.inlets))
        self.Q_in = np.zeros(len(self.inlets))

    def depths(self):
        return np.array([op.inlet.get_average_depth() for op in self.inlets])

    def volumes(self):
        return np.array([op.inlet.get_total_water_volume() for op in self.inlets])

    def step(self, dt):
        """Advance the coupling by dt and return the (Q_in, anuga_flux) used."""
        depths = self.depths()
        heads = self.backend.get_heads()

        Q = calculate_Q(heads, depths, self.beds, self.weir_lengths,
                        self.manhole_areas, cw=self.cw, co=self.co, g=self.g)
        if self.capture is not None:
            Q = self._limit_capture(Q, depths)
        Q = smooth_Q(Q, self.Q_in, dt, self.time_average)
        if self.clamp:
            Q = limit_outflow(Q, self.volumes(), dt, self.safety_factor)
        self.Q_in = Q

        self.backend.step(Q, dt)

        flux = self.backend.anuga_flux(Q, dt)
        for op, f in zip(self.inlets, flux):
            op.set_Q(f)

        if self.logger is not None:
            self._log_step(dt, depths, heads, Q)

        return CouplingStep(Q_in=Q, anuga_flux=flux)

    def approach_flows(self):
        """Surface flow approaching each inlet [m^3/s].

        With an approach width W, this is the area-weighted mean of the
        specific-discharge magnitude h|u| over the approach triangles, times
        W: the flow through a strip of that width, as TUFLOW-SWMM sums it
        over an on-grade inlet's connected cells. Without one, it is the
        magnitude of the area-averaged momentum over the approach triangles
        (or the inlet's own region) times sqrt(inlet area), the same
        surface-side heuristic as the standalone capture model."""
        q = np.empty(len(self.inlets))
        for k, (op, idx, width) in enumerate(zip(self.inlets, self.approach_indices,
                                                 self.approach_widths)):
            if width is not None:
                dom = op.domain
                idx = op.inlet.triangle_indices if idx is None else idx
                a = dom.areas[idx]
                qmag = np.hypot(dom.quantities['xmomentum'].centroid_values[idx],
                                dom.quantities['ymomentum'].centroid_values[idx])
                q[k] = np.dot(qmag, a) / a.sum() * width
                continue
            if idx is None:
                uh, vh = op.inlet.get_average_xmom(), op.inlet.get_average_ymom()
            else:
                dom = op.domain
                a = dom.areas[idx]
                uh = np.dot(dom.quantities['xmomentum'].centroid_values[idx], a) / a.sum()
                vh = np.dot(dom.quantities['ymomentum'].centroid_values[idx], a) / a.sum()
            q[k] = np.hypot(uh, vh) * np.sqrt(op.inlet.get_area())
        return q

    def _limit_capture(self, Q, depths):
        """Cap positive (surface->pipe) Q at each inlet's capture capacity."""
        approach = self.approach_flows()
        cap = np.array([np.inf if c is None else c(d, a)
                        for c, d, a in zip(self.capture, depths, approach)])
        return np.where(Q > 0, np.minimum(Q, cap), Q)

    def _log_step(self, dt, depths, heads, Q):
        """Feed one row per inlet to the hydrograph logger."""
        time = self.inlets[0].domain.get_time() if self.inlets else 0.0
        self.logger.record(time, dt, depths, heads, self.approach_flows(), Q)
