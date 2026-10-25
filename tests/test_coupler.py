"""Tests for the coupling driver (smooth_Q, limit_outflow, Coupler.step).

These exercise the pure flux mechanics with a fake backend and fake inlets, so
they need no ANUGA/SWMM/pipedream install. calculate_Q itself pulls anuga.g, so
Coupler.step tests that call it are skipped when ANUGA is unavailable.
"""
import numpy as np
import pytest

from anuga_drainage.coupler import smooth_Q, limit_outflow, Coupler


# --- pure helpers -----------------------------------------------------------

def test_smooth_Q_blends_old_and_new():
    Q_old = np.array([0.0, 0.0])
    Q_new = np.array([10.0, 10.0])
    # dt=1, time_average=10 -> 10% of the new value moves in.
    out = smooth_Q(Q_new, Q_old, dt=1.0, time_average=10.0)
    assert out == pytest.approx([1.0, 1.0])


def test_smooth_Q_disabled_returns_new():
    Q_new = np.array([3.0, -2.0])
    assert smooth_Q(Q_new, np.zeros(2), dt=1.0, time_average=0.0) is Q_new


def test_smooth_Q_dt_equals_time_average_is_full_step():
    Q_new = np.array([5.0])
    out = smooth_Q(Q_new, np.array([99.0]), dt=2.0, time_average=2.0)
    assert out == pytest.approx([5.0])


def test_limit_outflow_caps_positive_flux():
    Q = np.array([100.0, 0.5, -100.0])
    available = np.array([1.0, 1.0, 1.0])  # /dt=1 -> limit 1.0
    out = limit_outflow(Q, available, dt=1.0)
    # positive flux capped at available/dt; small positive and surcharge untouched.
    assert out == pytest.approx([1.0, 0.5, -100.0])


def test_limit_outflow_safety_factor():
    Q = np.array([100.0])
    out = limit_outflow(Q, np.array([2.0]), dt=2.0, safety_factor=0.5)
    # 0.5 * 2.0 / 2.0 = 0.5
    assert out == pytest.approx([0.5])


# --- Coupler.step orchestration with fakes ----------------------------------

class _FakeInlet:
    def __init__(self, depth, volume):
        self._depth = depth
        self._volume = volume
        self.Q_set = None

    class _I:
        pass

    @property
    def inlet(self):
        i = self._I()
        i.get_average_depth = lambda: self._depth
        i.get_total_water_volume = lambda: self._volume
        return i

    def set_Q(self, q):
        self.Q_set = q


class _FakeBackend:
    """Records the flux it was stepped with; feeds back -Q_in (pipedream-like)."""

    def __init__(self, heads):
        self._heads = np.asarray(heads, dtype=float)
        self.stepped_with = None

    def get_heads(self):
        return self._heads

    def step(self, Q_in, dt):
        self.stepped_with = np.array(Q_in)

    def anuga_flux(self, Q_in, dt):
        return -np.asarray(Q_in)


def test_coupler_step_feeds_realised_flux_to_inlets():
    inlets = [_FakeInlet(depth=1.0, volume=10.0), _FakeInlet(depth=1.0, volume=10.0)]
    backend = _FakeBackend(heads=[-1.0, -1.0])  # heads below bed -> weir inflow > 0
    coupler = Coupler(inlets, beds=[0.0, 0.0], weir_lengths=[2.0, 2.0],
                      manhole_areas=[1.0, 1.0], backend=backend, time_average=0.0,
                      g=9.81)

    step = coupler.step(dt=1.0)

    # Positive surface->pipe flux was passed to the backend ...
    assert np.all(step.Q_in > 0)
    assert backend.stepped_with == pytest.approx(step.Q_in)
    # ... and the realised flux (-Q_in here) was applied to each inlet.
    for op, f in zip(inlets, step.anuga_flux):
        assert op.Q_set == pytest.approx(f)
    assert np.all(step.anuga_flux < 0)


def test_coupler_smoothing_state_persists_across_steps():
    inlets = [_FakeInlet(depth=1.0, volume=10.0)]
    backend = _FakeBackend(heads=[-1.0])
    coupler = Coupler(inlets, beds=[0.0], weir_lengths=[2.0], manhole_areas=[1.0],
                      backend=backend, time_average=10.0, g=9.81)

    first = coupler.step(dt=1.0).Q_in.copy()
    second = coupler.step(dt=1.0).Q_in.copy()
    # With heads/depths constant, smoothing ramps the flux up toward the target,
    # so the second step exceeds the first.
    assert second[0] > first[0]


class _MovingInlet(_FakeInlet):
    """A fake inlet that also reports momentum and area (for approach flow)."""

    def __init__(self, depth, volume, uh=0.0, vh=0.0, area=4.0):
        super().__init__(depth, volume)
        self._uh, self._vh, self._area = uh, vh, area

    @property
    def inlet(self):
        i = super().inlet
        i.get_average_xmom = lambda: self._uh
        i.get_average_ymom = lambda: self._vh
        i.get_area = lambda: self._area
        return i


def test_coupler_capture_caps_surface_to_pipe_flow():
    # Inlet 0 is capped at 0.1 m^3/s; inlet 1 has no cap (None).
    inlets = [_MovingInlet(1.0, 10.0), _MovingInlet(1.0, 10.0)]
    backend = _FakeBackend(heads=[-1.0, -1.0])
    coupler = Coupler(inlets, beds=[0.0, 0.0], weir_lengths=[2.0, 2.0],
                      manhole_areas=[1.0, 1.0], backend=backend, g=9.81,
                      capture=[lambda d, q: 0.1, None])
    Q = coupler.step(dt=1.0).Q_in
    assert Q[0] == pytest.approx(0.1)
    assert Q[1] > 1.0                                  # uncapped weir/orifice flow


def test_coupler_capture_never_limits_surcharge():
    # Pipe head above the surface -> negative Q (surcharge); the cap is for
    # capture only, so it must leave surcharge untouched.
    inlets = [_MovingInlet(0.1, 10.0)]
    backend = _FakeBackend(heads=[2.0])
    free = Coupler(inlets, beds=[0.0], weir_lengths=[2.0], manhole_areas=[1.0],
                   backend=backend, g=9.81).step(dt=1.0).Q_in
    capped = Coupler(inlets, beds=[0.0], weir_lengths=[2.0], manhole_areas=[1.0],
                     backend=backend, g=9.81,
                     capture=[lambda d, q: 0.0]).step(dt=1.0).Q_in
    assert free[0] < 0
    assert capped == pytest.approx(free)


def test_coupler_capture_sees_depth_and_approach_flow():
    seen = []

    def cap(depth, approach):
        seen.append((depth, approach))
        return 1.0

    # |(uh, vh)| = 0.5 m^2/s over a width of sqrt(4) = 2 m -> 1.0 m^3/s approach.
    inlets = [_MovingInlet(0.3, 10.0, uh=0.3, vh=0.4, area=4.0)]
    Coupler(inlets, beds=[0.0], weir_lengths=[2.0], manhole_areas=[1.0],
            backend=_FakeBackend(heads=[-1.0]), g=9.81, capture=[cap]).step(dt=1.0)
    assert seen == [pytest.approx((0.3, 1.0))]


def test_coupler_capture_length_must_match_inlets():
    with pytest.raises(ValueError, match="capture"):
        Coupler([_MovingInlet(1.0, 1.0)], beds=[0.0], weir_lengths=[1.0],
                manhole_areas=[1.0], backend=_FakeBackend(heads=[0.0]),
                capture=[None, None])


def test_coupler_approach_flow_measured_in_given_triangles():
    # The inlet's own region has zero momentum (zero_velocity); the approach
    # triangles carry |(uh, vh)| = 0.5 m^2/s, so approach = 0.5 * sqrt(4) = 1.0.
    class _Q:
        def __init__(self, v):
            self.centroid_values = np.asarray(v, dtype=float)

    class _Dom:
        areas = np.array([1.0, 1.0, 2.0])
        quantities = {"xmomentum": _Q([0.0, 0.3, 0.3]), "ymomentum": _Q([0.0, 0.4, 0.4])}

    inlet = _MovingInlet(0.3, 10.0, uh=0.0, vh=0.0, area=4.0)
    inlet.domain = _Dom()
    seen = []
    Coupler([inlet], beds=[0.0], weir_lengths=[2.0], manhole_areas=[1.0],
            backend=_FakeBackend(heads=[-1.0]), g=9.81,
            capture=[lambda d, q: seen.append(q) or 1.0],
            approach_indices=[np.array([1, 2])]).step(dt=1.0)
    assert seen == [pytest.approx(1.0)]
