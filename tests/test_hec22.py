"""Tests for the HEC-22 inlet capture port (pure — no ANUGA/SWMM needed).

Expected values are worked by hand from the HEC-22 equations as SWMM 5.2's
inlet.c applies them (US units: 3.0 weir coefficient, 2g = 64.32 ft/s^2).
"""
import math

import pytest

from anuga_drainage.hec22 import (
    FT, CFS, Grate, Curb, InletDesign, Street, InletCapture,
    sag_capture, on_grade_capture, read_inlet_designs,
)

# A combination inlet and a street to test them with.
COMBO = InletDesign("Combo", grate=Grate(0.6, 0.3, "P_BAR-50"),
                    curb=Curb(0.9, 0.15, "VERTICAL"))
GRATE_ONLY = InletDesign("G", grate=Grate(0.6, 0.3, "P_BAR-50"))
STREET = Street(cross_slope=0.04, roughness=0.014, crown_width=7.0)


def test_sag_combo_weir_regime_matches_hand_calculation():
    # d = 0.1 m: grate weir on P = L + 2W over the depth averaged across the
    # cross slope, plus the curb opening beyond the grate (0.3 m) as a weir.
    d = 0.1 / FT
    Lg, Wg, Sx = 0.6 / FT, 0.3 / FT, 0.04
    di = d - (Wg / 2) * Sx
    q_grate = 3.0 * (Lg + 2 * Wg) * di ** 1.5
    q_sweep = 3.0 * (0.3 / FT) * d ** 1.5
    assert sag_capture(COMBO, 0.1, STREET) == pytest.approx((q_grate + q_sweep) * CFS)
    assert sag_capture(COMBO, 0.1, STREET) == pytest.approx(0.073, abs=0.001)


def test_sag_grate_orifice_regime():
    # Deep water: the grate switches to orifice flow on its open area (90 %).
    d = 0.5 / FT
    Lg, Wg, Sx = 0.6 / FT, 0.3 / FT, 0.04
    di = d - (Wg / 2) * Sx
    Ao = Lg * Wg * 0.90
    expected = 0.67 * Ao * math.sqrt(64.32 * di) * CFS
    assert sag_capture(GRATE_ONLY, 0.5, STREET) == pytest.approx(expected)


def test_sag_drop_grate_without_street_uses_full_perimeter():
    d = 0.05 / FT
    expected = 3.0 * 2 * (0.6 / FT + 0.3 / FT) * d ** 1.5 * CFS
    assert sag_capture(GRATE_ONLY, 0.05) == pytest.approx(expected)


def test_sag_capture_is_zero_when_dry_and_grows_with_depth():
    assert sag_capture(COMBO, 0.0, STREET) == 0.0
    q = [sag_capture(COMBO, h, STREET) for h in (0.02, 0.05, 0.1, 0.2, 0.5)]
    assert q == sorted(q)


def test_on_grade_captures_part_of_the_approach_flow():
    for q in (0.01, 0.1, 1.0):
        c = on_grade_capture(COMBO, q, STREET, long_slope=0.01)
        assert 0.0 < c < q
    # Efficiency falls as the gutter flow (and spread) grows.
    eff = [on_grade_capture(COMBO, q, STREET, 0.01) / q for q in (0.01, 0.1, 1.0)]
    assert eff == sorted(eff, reverse=True)


def test_on_grade_ignores_trickle_flow():
    assert on_grade_capture(COMBO, 1e-6, STREET, 0.01) == 0.0


def test_on_grade_inlets_in_series_capture_more():
    one = on_grade_capture(COMBO, 0.5, STREET, 0.01, number=1)
    two = on_grade_capture(COMBO, 0.5, STREET, 0.01, number=2)
    assert one < two < 0.5


def test_inlet_capture_usage_options():
    sag = InletCapture(COMBO, STREET, "ON_SAG")
    q = sag(0.2)
    assert InletCapture(COMBO, STREET, "ON_SAG", clogged_pct=50)(0.2) == pytest.approx(q / 2)
    assert InletCapture(COMBO, STREET, "ON_SAG", number=2)(0.2) == pytest.approx(2 * q)
    assert InletCapture(COMBO, STREET, "ON_SAG", q_max=0.05)(0.2) == pytest.approx(0.05)
    grade = InletCapture(COMBO, STREET, "ON_GRADE", long_slope=0.01)
    assert grade(0.2, 0.0) == 0.0                     # on grade: no approach flow, no capture
    auto = InletCapture(COMBO, STREET, "AUTOMATIC", long_slope=0.01)
    assert auto(0.2, 0.3) == pytest.approx(max(sag(0.2), grade(0.2, 0.3)))


def test_inlet_capture_rejects_bad_usage():
    with pytest.raises(ValueError, match="slope"):
        InletCapture(COMBO, STREET, "ON_GRADE")
    with pytest.raises(ValueError, match="placement"):
        InletCapture(COMBO, STREET, "SIDEWAYS")
    with pytest.raises(ValueError, match="grate type"):
        Grate(0.6, 0.3, "NOT_A_GRATE")


def test_read_inlet_designs(tmp_path):
    p = tmp_path / "net.inp"
    p.write_text("""\
[STREETS]
;;Name  Tcrown  Hcurb  Sx   nRoad  a    W    Sides  Tback  Sback  nBack
Xsec   7.0     0.2    4.0  0.014  0.0  0.0  1      5.0    2.0    0.016

[INLETS]
;;Name  Type   Param1  Param2  Param3
Combo  GRATE  0.6     0.3     P_BAR-50
Combo  CURB   0.9     0.15    VERTICAL
Drop1   GRATE  1.0     1.0     GENERIC  0.5  2.0
""")
    designs, streets = read_inlet_designs(str(p))
    assert designs["Combo"] == COMBO
    assert designs["Drop1"].grate == Grate(1.0, 1.0, "GENERIC", 0.5, 2.0)
    assert streets["Xsec"] == STREET
