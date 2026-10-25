"""Tests for the PipedreamBackend coupled-subset support and couple_from_inp.

The PipedreamBackend test needs pipedream; the couple_from_inp test additionally
needs ANUGA. Both skip cleanly when those aren't installed.
"""
import numpy as np
import pytest

from anuga_drainage import read_inp, inp_to_pipedream, PipedreamBackend

_INP = """\
[JUNCTIONS]
;;Name  Elev  MaxD  InitD  Sur  Apond
J1      0.0   1.0   0      0    0
J2      0.0   1.0   0      0    0

[OUTFALLS]
O1      0.0   FREE         NO

[CONDUITS]
C1      J1    J2    8.0    0.013   0   0   0   0
C2      J2    O1    5.0    0.013   0   0   0   0

[XSECTIONS]
C1      CIRCULAR   0.5   0   0   0   1
C2      CIRCULAR   0.5   0   0   0   1

[COORDINATES]
J1      15.0   5.0
J2       8.0   5.0
O1       3.0   5.0
"""


@pytest.fixture
def inp_path(tmp_path):
    p = tmp_path / "net.inp"
    p.write_text(_INP)
    return str(p)


def test_pipedream_backend_couples_subset_and_holds_bc(inp_path):
    SuperLink = pytest.importorskip("pipedream_solver.hydraulics").SuperLink
    inp = read_inp(inp_path)
    sj, sl = inp_to_pipedream(inp)
    s = SuperLink(sl, sj, internal_links=4)
    njunc = len(inp.junctions)                     # 2 junctions, 1 outfall
    be = PipedreamBackend(s, coupled_indices=range(njunc), H_bc=s._z_inv_j.copy())

    assert len(be.get_heads()) == njunc            # only the junctions are coupled
    be.step(np.array([0.05, 0.0]), dt=1.0)         # inject at J1 only
    assert float(s.H_j[-1]) == pytest.approx(float(s._z_inv_j[-1]))  # outfall held at its invert
    assert s.H_j[0] > s._z_inv_j[0]                # J1 filled
    assert list(be.coupling_inflow_volumes()) == pytest.approx([0.05, 0.0])


def test_pipedream_backend_tracks_outfall_volume(inp_path):
    # outfall_volume should capture the flow shed at the bc outfall, so the pipe
    # budget closes: dV_pipe ~ injected - outfall (to pipedream's own residual).
    SuperLink = pytest.importorskip("pipedream_solver.hydraulics").SuperLink
    inp = read_inp(inp_path)
    sj, sl = inp_to_pipedream(inp)
    s = SuperLink(sl, sj, internal_links=4)
    njunc, noutfall = len(inp.junctions), len(inp.outfalls)
    be = PipedreamBackend(s, coupled_indices=range(njunc), H_bc=s._z_inv_j.copy(),
                          outfall_indices=range(njunc, njunc + noutfall))
    v0 = be.pipe_volume()
    injected = 0.0
    for _ in range(60):
        be.step(np.array([0.2, 0.0]), dt=1.0)
        injected += 0.2
    r_pipe = (be.pipe_volume() - v0) - (injected - be.outfall_volume())
    assert be.outfall_volume() > 5.0                 # water really left at the outfall
    assert abs(r_pipe) < 0.05 * injected             # pipe budget closes (pipedream residual)


def test_pipedream_backend_default_is_all_coupled_no_bc(inp_path):
    # Backwards compatibility: no coupled_indices/H_bc -> couple every
    # superjunction and never pass H_bc (matches the hand-built examples).
    SuperLink = pytest.importorskip("pipedream_solver.hydraulics").SuperLink
    sj, sl = inp_to_pipedream(read_inp(inp_path))
    s = SuperLink(sl, sj, internal_links=4)
    be = PipedreamBackend(s)
    assert len(be.get_heads()) == len(s.H_j)
    assert be.H_bc is None


def test_couple_from_inp_pipedream(inp_path):
    anuga = pytest.importorskip("anuga")
    pytest.importorskip("pipedream_solver.hydraulics")
    from anuga_drainage import couple_from_inp

    domain = anuga.rectangular_cross_domain(20, 10, len1=20.0, len2=10.0)
    domain.set_quantity("elevation", 0.0)
    domain.set_quantity("stage", 0.3)
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({"left": Br, "right": Br, "top": Br, "bottom": Br})

    c = couple_from_inp(domain, inp_path, backend="pipedream",
                        manhole_area=0.5, time_average=1.0, internal_links=4)

    # junctions are coupled, the outfall is not
    assert list(c.inlets) == ["J1", "J2"]
    assert len(c.coupler.inlets) == 2
    assert c.backend.H_bc is not None              # outfall is a free-drainage boundary

    c.add_volume_balance()                         # audit runs inside c.step()
    last = None
    for t in domain.evolve(yieldstep=1.0, finaltime=2.0):
        last = c.step(1.0)                          # one call: exchange + audit
    assert last.Q_in.shape == (2,)
    assert np.isfinite(last.Q_in).all()
    assert len(c.volume_balance.records) >= 1      # the audit recorded steps
    c.close()                                      # no-op for pipedream, but must not raise


def test_couple_from_inp_inlet_specs_decouple_hydraulics_from_footprint(inp_path):
    anuga = pytest.importorskip("anuga")
    pytest.importorskip("pipedream_solver.hydraulics")
    from anuga_drainage import couple_from_inp, INLET_LIBRARY

    domain = anuga.rectangular_cross_domain(20, 10, len1=20.0, len2=10.0)
    domain.set_quantity("elevation", 0.0)
    domain.set_quantity("stage", 0.3)
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({"left": Br, "right": Br, "top": Br, "bottom": Br})

    # J1 gets a catalogue grate (50% blocked); J2 keeps the footprint default.
    c = couple_from_inp(domain, inp_path, backend="pipedream",
                        manhole_area=0.5, internal_links=4,
                        inlet_specs={"J1": "Grate_600x600"}, blockage={"J1": 0.5})

    grate = INLET_LIBRARY["Grate_600x600"]
    # J1 (index 0): hydraulic area/perimeter come from the spec, derated 50%.
    assert c.coupler.manhole_areas[0] == pytest.approx(grate.clear_area * 0.5)
    assert c.coupler.weir_lengths[0] == pytest.approx(grate.effective_perimeter * 0.5)
    # J2 (index 1, no spec): keeps the footprint-derived hydraulic area (=manhole_area).
    assert c.coupler.manhole_areas[1] == pytest.approx(0.5)
    # Footprint is decoupled: J1's ANUGA region is far larger than the grate opening.
    assert c.coupler.inlets[0].inlet.get_area() > grate.operational_area
    c.close()


def test_couple_from_inp_logs_hydrographs(inp_path, tmp_path):
    anuga = pytest.importorskip("anuga")
    pytest.importorskip("pipedream_solver.hydraulics")
    from anuga_drainage import couple_from_inp

    domain = anuga.rectangular_cross_domain(20, 10, len1=20.0, len2=10.0)
    domain.set_quantity("elevation", 0.0)
    domain.set_quantity("stage", 0.3)
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({"left": Br, "right": Br, "top": Br, "bottom": Br})

    c = couple_from_inp(domain, inp_path, backend="pipedream",
                        manhole_area=0.5, internal_links=4, log_hydrographs=True)
    assert c.coupler.logger is not None

    for t in domain.evolve(yieldstep=1.0, finaltime=3.0):
        c.step(1.0)

    df = c.coupler.logger.to_dataframe("J1")
    assert not df.empty
    for col in ("Time_s", "Depth_m", "Head1D_m", "Captured_Q_cms",
                "Surcharge_Q_cms", "Cum_Captured_m3"):
        assert col in df.columns
    paths = c.coupler.logger.write_csv(directory=str(tmp_path))
    assert len(paths) == 2                          # one per junction (J1, J2)
    c.close()


def _flat_domain(anuga):
    domain = anuga.rectangular_cross_domain(20, 10, len1=20.0, len2=10.0)
    domain.set_quantity("elevation", 0.0)
    domain.set_quantity("stage", 0.3)
    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({"left": Br, "right": Br, "top": Br, "bottom": Br})
    return domain


def test_couple_from_inp_coupled_junctions_subset_pipedream(inp_path):
    # J2 is a sealed manhole: in the network but with no ANUGA inlet.
    anuga = pytest.importorskip("anuga")
    pytest.importorskip("pipedream_solver.hydraulics")
    from anuga_drainage import couple_from_inp

    domain = _flat_domain(anuga)
    c = couple_from_inp(domain, inp_path, backend="pipedream", manhole_area=[0.5, 0.7],
                        internal_links=4, coupled_junctions=["J1"])
    assert list(c.inlets) == ["J1"]
    assert list(c.backend.coupled) == [0]              # J1 only; J2 stays in the network
    assert len(c.backend.get_heads()) == 1
    assert c.inlets["J1"].inlet.get_area() >= 0.5       # J1's own manhole_area entry
    c.close()


def test_couple_from_inp_coupled_junctions_rejects_unknown_names(inp_path):
    anuga = pytest.importorskip("anuga")
    pytest.importorskip("pipedream_solver.hydraulics")
    from anuga_drainage import couple_from_inp

    with pytest.raises(ValueError, match="coupled_junctions"):
        couple_from_inp(_flat_domain(anuga), inp_path, backend="pipedream",
                        coupled_junctions=["J1", "nope"])


_SWMM_OPTIONS = """\
[OPTIONS]
FLOW_UNITS     CMS
FLOW_ROUTING   DYNWAVE
START_DATE     01/01/2020
START_TIME     00:00:00
END_DATE       01/01/2020
END_TIME       00:01:00
ROUTING_STEP   1

"""


def test_couple_from_inp_coupled_junctions_subset_swmm(tmp_path):
    # SWMM backend: only J1 is coupled, but pipe_volume must still count the
    # water stored at the uncoupled J2, or R_pipe would depend on the subset.
    anuga = pytest.importorskip("anuga")
    pytest.importorskip("pyswmm")
    from anuga_drainage import couple_from_inp

    inp_path = tmp_path / "net_swmm.inp"
    inp_path.write_text(_SWMM_OPTIONS + _INP)
    c = couple_from_inp(_flat_domain(anuga), str(inp_path), backend="swmm",
                        manhole_area=0.5, coupled_junctions=["J1"])
    assert list(c.inlets) == ["J1"]
    assert [n.nodeid for n in c.backend.junctions] == ["J1"]
    assert sorted(n.nodeid for n in c.backend.storage_nodes) == ["J1", "J2"]
    c.close()
