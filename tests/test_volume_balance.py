"""Tests for the VolumeBalance residual algebra (with fakes — no ANUGA needed)."""
from collections import namedtuple

import pytest

from anuga_drainage.volume_balance import VolumeBalance

_Step = namedtuple("_Step", ["Q_in", "anuga_flux"])


class _FakeDomain:
    water_volume = 0.0
    boundary = 0.0

    def get_water_volume(self):
        return self.water_volume

    def get_boundary_flux_integral(self):
        return self.boundary


class _FakeOp:
    applied = 0.0

    def get_total_applied_volume(self):
        return self.applied


class _FakeBackend:
    pv = 0.0
    cin = 0.0
    out = 0.0

    def pipe_volume(self):
        return self.pv

    def coupling_inflow_volume(self):
        return self.cin

    def outfall_volume(self):
        return self.out


def _make():
    dom, inflow, inlet, be = _FakeDomain(), _FakeOp(), _FakeOp(), _FakeBackend()
    vb = VolumeBalance(dom, [inlet], be, inflow_operators=[inflow])
    vb.step(0.0)  # baseline (all zero)
    return vb, dom, inflow, inlet, be


def test_fully_consistent_system_has_zero_residuals():
    # Outfall water returns to ANUGA at inlet 0 (internal transfer).
    dom, inflow, inlet, be = _FakeDomain(), _FakeOp(), _FakeOp(), _FakeBackend()
    vb = VolumeBalance(dom, [inlet], be, inflow_operators=[inflow], outfall_inlet=0)
    vb.step(0.0)  # baseline (all zero)
    # 10 in via the source, 2 out via the boundary.
    inflow.applied = 10.0
    dom.boundary = -2.0
    # Coupling: 4 m^3 surface->pipe at the inlet; 1 m^3 pipe->surface at the outfall.
    inlet.applied = -4.0 + 1.0          # inlet exchange + outfall return (ANUGA's record)
    be.cin = 4.0                        # pipe received 4 at the inlet
    be.out = 1.0                        # pipe lost 1 at the outfall
    dom.water_volume = 10.0 - 2.0 - 3.0  # ANUGA closes -> 5
    be.pv = 4.0 - 1.0                    # pipe closes -> 3

    r = vb.step(1.0)
    assert r.R_anuga == pytest.approx(0.0, abs=1e-12)
    assert r.R_pipe == pytest.approx(0.0, abs=1e-12)
    assert r.R_couple == pytest.approx(0.0, abs=1e-12)
    assert r.loss == pytest.approx(0.0, abs=1e-12)


def test_outfall_that_leaves_the_system_is_a_loss_not_a_handoff():
    # No outfall_inlet: the 1 m^3 that leaves at the outfall exits the system,
    # so R_couple stays ~0 (inlet handoff consistent) and the leaving water
    # shows up in `loss`.
    dom, inflow, inlet, be = _FakeDomain(), _FakeOp(), _FakeOp(), _FakeBackend()
    vb = VolumeBalance(dom, [inlet], be, inflow_operators=[inflow])  # outfall leaves
    vb.step(0.0)
    inflow.applied = 10.0
    inlet.applied = -4.0       # 4 surface->pipe; nothing returns to ANUGA
    be.cin = 4.0               # pipe received 4
    be.out = 1.0               # pipe lost 1 at the outfall (exits the system)
    dom.water_volume = 10.0 - 4.0   # ANUGA closes -> 6
    be.pv = 4.0 - 1.0               # pipe closes -> 3
    r = vb.step(1.0)
    assert r.R_anuga == pytest.approx(0.0, abs=1e-12)
    assert r.R_pipe == pytest.approx(0.0, abs=1e-12)
    assert r.R_couple == pytest.approx(0.0, abs=1e-12)
    assert r.loss == pytest.approx(0.0, abs=1e-12)  # 9 - 9 + (1 - 0); the sink closes too


def test_loss_splits_into_the_three_residuals():
    vb, dom, inflow, inlet, be = _make()
    inflow.applied = 10.0
    inlet.applied = -4.0
    be.cin = 4.0
    # Introduce a pipe (SWMM-style) non-conservation: 0.3 m^3 unaccounted in the
    # pipe budget; ANUGA and the coupling stay consistent.
    dom.water_volume = 10.0 - 4.0       # ANUGA closes
    be.pv = 4.0 - 0.3                   # pipe is short by 0.3
    r = vb.step(1.0)
    assert r.R_anuga == pytest.approx(0.0, abs=1e-12)
    assert r.R_couple == pytest.approx(0.0, abs=1e-12)
    assert r.R_pipe == pytest.approx(-0.3, abs=1e-12)
    # the single loss is exactly the sum of the parts
    assert r.loss == pytest.approx(r.R_anuga + r.R_pipe + r.R_couple, abs=1e-12)


def test_coupling_mismatch_shows_in_R_couple():
    vb, dom, inflow, inlet, be = _make()
    # 4 left the surface but only 3 reached the pipe -> 1 lost at the handoff.
    inlet.applied = -4.0
    be.cin = 3.0
    dom.water_volume = -4.0   # ANUGA still closes (it removed 4)
    be.pv = 3.0               # pipe closes (it got 3)
    r = vb.step(1.0)
    assert r.R_anuga == pytest.approx(0.0, abs=1e-12)
    assert r.R_pipe == pytest.approx(0.0, abs=1e-12)
    assert r.R_couple == pytest.approx(-1.0, abs=1e-12)


class _VecBackend(_FakeBackend):
    """Backend with per-junction accepted volumes for the per-inlet breakdown."""
    vols = (0.0, 0.0)

    def coupling_inflow_volumes(self):
        return list(self.vols)

    def coupling_inflow_volume(self):
        return sum(self.vols)


def test_per_inlet_breakdown_catches_inlet_drying():
    dom = _FakeDomain()
    inlet0, inlet1 = _FakeOp(), _FakeOp()
    be = _VecBackend()
    vb = VolumeBalance(dom, [inlet0, inlet1], be)
    vb.step(0.0, dt=1.0, coupling_step=_Step(Q_in=[0.0, 0.0], anuga_flux=[0.0, 0.0]))

    # Request 2 into the pipe at inlet 0 and 3 out at inlet 1; the sewer takes it
    # all, but ANUGA only manages to remove 1.5 at inlet 0 (the cell dried out).
    be.vols = (2.0, -3.0)         # accepted by the sewer
    inlet0.applied = -1.5         # ANUGA removed only 1.5 (drying!)
    inlet1.applied = 3.0          # inlet 1 fully realised
    vb.step(1.0, dt=1.0, coupling_step=_Step(Q_in=[2.0, -3.0], anuga_flux=[0.0, 0.0]))

    p = vb.per_inlet[-1]
    assert list(p["requested"]) == pytest.approx([2.0, -3.0])
    assert list(p["accepted"]) == pytest.approx([2.0, -3.0])
    assert list(p["removed"]) == pytest.approx([-1.5, 3.0])
    assert p["drying"][0] == pytest.approx(0.5)     # inlet 0 over-drawn by 0.5
    assert p["drying"][1] == pytest.approx(0.0)     # inlet 1 consistent


def test_outfall_return_is_stripped_from_drying():
    # Inlet 1 receives a 1.0 m^3 outfall return on top of its weir exchange, so
    # its `removed` is inflated; with outfall_inlet=1 the `drying` column nets it
    # out and reads ~0 for a consistent handoff.
    dom = _FakeDomain()
    inlet0, inlet1 = _FakeOp(), _FakeOp()
    be = _VecBackend()
    vb = VolumeBalance(dom, [inlet0, inlet1], be, outfall_inlet=1)
    vb.step(0.0, dt=1.0, coupling_step=_Step(Q_in=[0.0, 0.0], anuga_flux=[0.0, 0.0]))

    be.vols = (2.0, -3.0)         # sewer accepted at the two junctions
    be.out = 1.0                  # 1.0 left the pipe at the outfall
    inlet0.applied = -2.0         # inlet 0: full weir exchange, no drying
    inlet1.applied = 3.0 + 1.0    # inlet 1: weir exchange (+3) plus outfall return (+1)
    vb.step(1.0, dt=1.0, coupling_step=_Step(Q_in=[2.0, -3.0], anuga_flux=[0.0, 0.0]))

    p = vb.per_inlet[-1]
    assert list(p["outfall_return"]) == pytest.approx([0.0, 1.0])
    assert p["drying"][0] == pytest.approx(0.0)   # 2 + (-2)
    assert p["drying"][1] == pytest.approx(0.0)   # -3 + 4 - 1 (outfall stripped)


def test_baseline_is_taken_on_first_step_not_construction():
    # Construct with a nonzero (e.g. unphysical) initial water volume; the
    # baseline should be the first step() reading, so residuals stay clean.
    dom, be = _FakeDomain(), _FakeBackend()
    dom.water_volume = -3555.0
    vb = VolumeBalance(dom, [], be)
    vb.step(0.0)
    assert vb.V_anuga0 == -3555.0
    dom.water_volume = -3555.0 + 7.0   # 7 m^3 added
    r = vb.step(1.0)
    assert r.R_anuga == pytest.approx(7.0, abs=1e-12)  # no inflow/inlets accounted -> shows as residual


class _FakeRateOp:
    """Mimics an ANUGA Rate_operator: no get_total_applied_volume, only the
    running cumulative_influx."""
    cumulative_influx = 0.0


def test_rate_operator_inflow_is_counted_via_cumulative_influx():
    # Rain on the grid (a Rate_operator) feeds the domain alongside an
    # Inlet_operator source; both must count as inflow for ANUGA to close.
    dom, inflow, inlet, be = _FakeDomain(), _FakeOp(), _FakeOp(), _FakeBackend()
    rain = _FakeRateOp()
    vb = VolumeBalance(dom, [inlet], be, inflow_operators=[inflow, rain])
    vb.step(0.0)
    inflow.applied = 10.0
    rain.cumulative_influx = 5.0
    dom.water_volume = 15.0
    r = vb.step(1.0)
    assert r.inflow == pytest.approx(15.0)
    assert r.R_anuga == pytest.approx(0.0, abs=1e-12)
    assert r.loss == pytest.approx(0.0, abs=1e-12)
