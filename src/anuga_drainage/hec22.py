"""Street-inlet capture capacity from the HEC-22 equations, as implemented in SWMM 5.2.

A port of EPA SWMM 5.2's ``inlet.c`` (on-sag and on-grade grate, curb-opening
and combination inlets), which TUFLOW-SWMM also uses for its pit inlets. In
the coupling it supplies a per-inlet **capture limit**: ``calculate_Q`` still
sets the head-driven 2D<->1D exchange, but the flow an inlet can take off the
surface is capped at its HEC-22 capacity (see ``Coupler(capture=...)``).

* ``sag_capture(design, depth, street)`` — an inlet at a low point, where the
  capture depends only on the ponded depth (weir/orifice, HEC-22 Eqs 4-26..4-31).
* ``on_grade_capture(design, approach_Q, street, long_slope)`` — an inlet on a
  slope, which captures part of the gutter flow passing it and lets the rest
  bypass (flow spread, frontal/side efficiency, HEC-22 Eqs 4-2..4-24).
* ``InletCapture`` — one inlet usage (design + street + placement, number of
  inlets, clogging, flow limit) as a ``capture(depth, approach_Q)`` callable.
* ``read_inlet_designs(inp_path)`` — the ``[INLETS]`` and ``[STREETS]``
  sections of a SWMM ``.inp``.

Inputs and outputs are SI (m, m^3/s). Internally the functions work in SWMM's
US customary units with SWMM's coefficients, so each line can be checked
against ``inlet.c``. The module is pure (no ANUGA/SWMM needed), so it is
unit-testable standalone.
"""
from dataclasses import dataclass
import math

FT = 0.3048                    # m per ft
CFS = FT ** 3                  # m^3/s per cfs
_G2 = 2.0 * 32.16              # 2g in ft/s^2, as in inlet.c
_MIN_FLOW = 0.001              # cfs, SWMM's MIN_RUNOFF_FLOW

# Grate open-area ratios (HEC-22 Chart 9B) and splash-over velocity cubic
# coefficients (HEC-22 Chart 5B), copied from inlet.c.
GRATE_OPENING_RATIOS = {
    "P_BAR-50": 0.90, "P_BAR-50X100": 0.80, "P_BAR-30": 0.60,
    "CURVED_VANE": 0.35, "TILT_BAR-45": 0.17, "TILT_BAR-30": 0.34,
    "RETICULINE": 0.80, "GENERIC": 1.00,
}
_SPLASH_COEFFS = {
    "P_BAR-50": (2.22, 4.03, 0.65, 0.06), "P_BAR-50X100": (0.74, 2.44, 0.27, 0.02),
    "P_BAR-30": (1.76, 3.12, 0.45, 0.03), "CURVED_VANE": (0.30, 4.85, 1.31, 0.15),
    "TILT_BAR-45": (0.99, 2.64, 0.36, 0.03), "TILT_BAR-30": (0.51, 2.34, 0.2, 0.01),
    "RETICULINE": (0.28, 2.28, 0.18, 0.01),
}


@dataclass
class Grate:
    """A grate opening [m]. ``grate_type`` is a SWMM grate name; ``GENERIC``
    takes ``open_fraction`` and ``splash_velocity`` [m/s] explicitly."""
    length: float
    width: float
    grate_type: str = "P_BAR-50"
    open_fraction: float = None
    splash_velocity: float = None

    def __post_init__(self):
        self.grate_type = self.grate_type.upper()
        if self.grate_type not in GRATE_OPENING_RATIOS:
            raise ValueError(f"unknown grate type {self.grate_type!r}")
        if self.grate_type == "GENERIC" and self.open_fraction is None:
            raise ValueError("a GENERIC grate needs open_fraction")


@dataclass
class Curb:
    """A curb (kerb) opening [m]; ``throat`` is VERTICAL, INCLINED or HORIZONTAL."""
    length: float
    height: float
    throat: str = "VERTICAL"

    def __post_init__(self):
        self.throat = self.throat.upper()
        if self.throat not in ("VERTICAL", "INCLINED", "HORIZONTAL"):
            raise ValueError(f"unknown throat angle {self.throat!r}")


@dataclass
class InletDesign:
    """A SWMM inlet design: a grate, a curb opening, or both (combination)."""
    name: str
    grate: Grate = None
    curb: Curb = None


@dataclass
class Street:
    """A SWMM street cross-section (one side). ``cross_slope`` is m/m (the
    ``.inp`` gives %); ``crown_width`` is curb to crown [m]."""
    cross_slope: float
    roughness: float
    crown_width: float
    gutter_depression: float = 0.0
    gutter_width: float = 0.0
    sides: int = 1


class _Geom:
    """The street-geometry shared variables of inlet.c, in ft
    (getConduitGeometry), including an inlet's local depression."""

    def __init__(self, street, local_depression=0.0, local_width=0.0):
        self.Sx = street.cross_slope
        self.a = street.gutter_depression / FT
        self.W = street.gutter_width / FT
        self.n = street.roughness
        self.Tcrown = street.crown_width / FT
        self.sides = street.sides
        if local_depression * local_width > 0:
            self.a += local_depression / FT
            self.W = local_width / FT
        self.Sw = self.Sx + self.a / self.W if self.W * self.a > 0.0 else self.Sx


# --------------------------------------------------------------------------
# On sag
# --------------------------------------------------------------------------
def _sag_grate(grate, d, g):
    """findOnSagGrateFlows: (Qw, Qo) in cfs for depth d [ft]; g None = drop grate."""
    Lg, Wg = grate.length / FT, grate.width / FT
    if g is None:
        di, P = d, 2.0 * (Lg + Wg)
    else:
        if d <= Wg * g.Sw:                     # spread within the grate width
            Wg = d / g.Sw
        di = d - (Wg / 2.0) * g.Sw             # average depth over the grate
        P = Lg + 2.0 * Wg                      # effective perimeter
    ratio = grate.open_fraction if grate.grate_type == "GENERIC" else GRATE_OPENING_RATIOS[grate.grate_type]
    Ao = Lg * Wg * ratio
    if di <= 0.0:
        return 0.0, 0.0
    if d <= 1.79 * Ao / P:
        return 3.0 * P * di ** 1.5, 0.0                    # HEC-22 Eq (4-26)
    return 0.0, 0.67 * Ao * math.sqrt(_G2 * di)            # HEC-22 Eq (4-27)


def _curb_orifice(di, h, L, throat):
    """getCurbOrificeFlow [cfs]."""
    d = di
    if throat == "HORIZONTAL":
        d = di - h / 2.0
    elif throat == "INCLINED":
        d = di - (h / 2.0) * 0.7071
    return 0.67 * h * L * math.sqrt(_G2 * max(d, 0.0))     # HEC-22 Eq (4-31a)


def _sag_curb(curb, d, L, g):
    """findOnSagCurbFlows: (Qw, Qo) in cfs for an opening length L [ft]."""
    if L <= 0.0:
        return 0.0, 0.0
    h = curb.height / FT
    a = g.a if g is not None else 0.0
    W = g.W if g is not None else 0.0
    dorif = 1.4 * h
    if d > dorif:
        return 0.0, _curb_orifice(d, h, L, curb.throat)
    if a == 0.0 or L > 12.0:                   # uniform cross slope / long opening
        dweir = h
        if d < dweir:
            return 3.0 * L * d ** 1.5, 0.0                 # HEC-22 Eq (4-30)
        Qweir = 3.0 * L * dweir ** 1.5
    else:                                      # depressed gutter
        P = L + 1.8 * W
        dweir = h + a
        if d < dweir:
            return 2.3 * P * d ** 1.5, 0.0                 # HEC-22 Eq (4-28)
        Qweir = 2.3 * P * dweir ** 1.5
    # interpolate between Qweir at dweir and Qorif at dorif
    Qorif = _curb_orifice(dorif, h, L, curb.throat)
    r = (d - dweir) / (dorif - dweir)
    return (1.0 - r) * Qweir, r * Qorif


def _sag_cfs(design, d, g):
    """getOnSagInletCapture [cfs] for depth d [ft]."""
    Qgw = Qgo = Qsw = Qso = Qco = 0.0
    Lgrate = design.grate.length / FT if design.grate else 0.0
    if Lgrate > 0.0:
        Qgw, Qgo = _sag_grate(design.grate, d, g)
    if design.curb is not None:
        Lcurb = design.curb.length / FT
        Lsweep = Lcurb - Lgrate
        if Lsweep > 0.0:
            Qsw, Qso = _sag_curb(design.curb, d, Lsweep, g)
        if Qgo > 0.0:                          # curb alongside a grate in orifice flow
            _, Qco = _sag_curb(design.curb, d, Lgrate, g)
    return Qgw + Qgo + Qsw + Qso + Qco


def sag_capture(design, depth, street=None, local_depression=0.0, local_width=0.0):
    """Flow [m^3/s] one on-sag inlet captures at ponded ``depth`` [m].

    With a ``street`` the grate sits in the gutter (perimeter L + 2W, depth
    averaged over the cross slope); without one it is a drop inlet (perimeter
    2(L + W)), as in SWMM."""
    if depth <= 0.0:
        return 0.0
    g = _Geom(street, local_depression, local_width) if street is not None else None
    return _sag_cfs(design, depth / FT, g) * CFS


# --------------------------------------------------------------------------
# On grade
# --------------------------------------------------------------------------
def _eo(Sr, Ts, w):
    """getEo: gutter-to-total flow ratio, HEC-22 Eq (4-4)."""
    x = Sr / (Ts / w)
    x = (1.0 + x) ** 2.67 - 1.0
    return 1.0 / (1.0 + Sr / x)


def _spread(Q, f, g):
    """getFlowSpread [ft] for street flow Q [cfs] (Izzard, HEC-22 Eq 4-2)."""
    if g.a == 0.0:
        Ts1 = (Q / f) ** 0.375
    else:
        f1 = f * ((g.a / g.W) / g.Sx) ** 1.67
        Tw = (Q / f1) ** 0.375
        if Tw <= g.W:
            Ts1 = Tw
        else:
            Sr = (g.Sx + g.a / g.W) / g.Sx
            Ts1 = (Q / f) ** 0.375 - g.W
            if Ts1 <= 0:
                Ts1 = Tw - g.W
            Ts2 = Ts1
            for _ in range(10):
                Eo = _eo(Sr, Ts1, g.W)
                Ts2 = ((1.0 - Eo) * Q / f) ** 0.375        # HEC-22 Eq (4-6), (4-2)
                if abs(Ts2 - Ts1) < 0.01:
                    break
                Ts1 = Ts2
            Ts1 = Ts2 + g.W
    return min(Ts1, g.Tcrown)


def _gutter_flow_ratio(w, T, g):
    """getGutterFlowRatio."""
    if T <= w:
        return 1.0
    if g.a > 0.0:
        return _eo(g.Sw / g.Sx, T - w, w)
    return 1.0 - (1.0 - w / T) ** 2.67                     # HEC-22 Eq (4-16)


def _gutter_area_ratio(Wg, A, T, g):
    """getGutterAreaRatio."""
    if Wg >= g.W or T <= Wg:
        return 1.0
    if T <= g.W:
        return Wg / T
    As = 0.5 * (T - g.W) ** 2 * g.Sx
    Ag = Wg * ((T * g.Sx) + g.a - (Wg * g.Sw / 2.0))
    return Ag / (A - As)


def _grate_on_grade(grate, Q, T, f, g):
    """getGrateInletCapture [cfs] for a street gutter grate."""
    Lg, Wg = grate.length / FT, grate.width / FT
    Qo = Q
    if g.a == 0.0:                             # conventional gutter
        A = T * T * g.Sx / 2.0
        Eo = _gutter_flow_ratio(Wg, T, g)
        if T >= g.Tcrown:
            Qo = f * g.Tcrown ** 2.67
    else:                                      # composite gutter
        A = T * T * g.Sw / 2.0 if T <= g.W else (T * T * g.Sx + g.a * g.W) / 2.0
        Eo = _gutter_flow_ratio(g.W, T, g)
        if Eo < 1.0:
            if T >= g.Tcrown:
                Qo = f * g.Tcrown ** 2.67 / (1.0 - Eo)
            Eo = Eo * _gutter_area_ratio(Wg, A, T, g)      # HEC-22 Eq (4-20a)
    V = Qo / A
    if grate.grate_type == "GENERIC":
        Vo = (grate.splash_velocity or 0.0) / FT
    else:
        c0, c1, c2, c3 = _SPLASH_COEFFS[grate.grate_type]
        Vo = c0 + c1 * Lg - c2 * Lg * Lg + c3 * Lg ** 3
    Rf = 1.0 - 0.09 * (V - Vo) if V > Vo else 1.0          # HEC-22 Eq (4-18)
    Rs = 0.0
    if Eo < 1.0:
        Rs = 1.0 / (1.0 + (0.15 * V ** 1.8 / g.Sx / Lg ** 2.3))   # HEC-22 Eq (4-19)
    return Q * (Rf * Eo + Rs * (1.0 - Eo))                 # HEC-22 Eq (4-21)


def _curb_on_grade(Q, L, T, SL, g):
    """getCurbInletCapture [cfs] for an opening length L [ft]."""
    Se = g.Sx
    if g.a > 0.0:
        Eo = _eo(g.Sw / g.Sx, T - g.W, g.W)
        Se = g.Sx + (g.a / g.W) * Eo                       # HEC-22 Eq (4-24)
    Lt = 0.6 * Q ** 0.42 * SL ** 0.3 * (1.0 / (g.n * Se)) ** 0.6   # HEC-22 Eq (4-22a)
    E = 1.0
    if L < Lt:
        E = 1.0 - (1.0 - L / Lt) ** 1.8                    # HEC-22 Eq (4-23)
    return min(max(E, 0.0), 1.0) * Q


def _on_grade_cfs(design, Q, SL, f, g):
    """getOnGradeInletCapture [cfs] for one inlet seeing street flow Q [cfs]."""
    T = _spread(Q, f, g)
    Q1, Qc = Q, 0.0
    Lgrate = design.grate.length / FT if design.grate else 0.0
    if design.curb is not None:
        Lsweep = design.curb.length / FT - Lgrate
        if Lsweep > 0.0:
            Qc = _curb_on_grade(Q1, Lsweep, T, SL, g)
            Q1 -= Qc
    if Lgrate > 0.0 and Q1 > 0.0:
        if Q1 != Q:
            T = _spread(Q1, f, g)
        Qc += _grate_on_grade(design.grate, Q1, T, f, g)
    return Qc


def on_grade_capture(design, approach_Q, street, long_slope, number=1,
                     clog_factor=1.0, q_max=None, local_depression=0.0, local_width=0.0):
    """Flow [m^3/s] captured from ``approach_Q`` [m^3/s] of street flow by
    ``number`` on-grade inlets in series on a street of longitudinal slope
    ``long_slope`` (m/m); the remainder bypasses (getOnGradeCapturedFlow)."""
    q = approach_Q / CFS
    if q < _MIN_FLOW:
        return 0.0
    g = _Geom(street, local_depression, local_width)
    f = (0.56 / g.n) * long_slope ** 0.5 * g.Sx ** 1.67    # Izzard flow factor
    q_lim = q_max / CFS if q_max else math.inf
    q /= g.sides
    bypassed, captured = q, 0.0
    for _ in range(number):
        qc = min(_on_grade_cfs(design, bypassed, long_slope, f, g) * clog_factor, q_lim, bypassed)
        captured += qc
        bypassed -= qc
        if bypassed < _MIN_FLOW:
            break
    return captured * g.sides * CFS


# --------------------------------------------------------------------------
# One inlet usage, as the coupling's capture limit
# --------------------------------------------------------------------------
class InletCapture:
    """Capture capacity of one inlet usage, as ``capture(depth, approach_Q)``
    [m^3/s] — the form ``Coupler(capture=...)`` calls each step.

    ``placement`` is ``ON_SAG`` (depth-driven) or ``ON_GRADE`` (approach-flow
    driven; needs ``long_slope``). ``AUTOMATIC`` takes the larger of the two,
    as TUFLOW-SWMM does. ``clogged_pct`` and ``q_max`` follow SWMM's
    [INLET_USAGE] fields; ``q_max`` of 0/None means no limit."""

    def __init__(self, design, street, placement="ON_SAG", long_slope=None,
                 number=1, clogged_pct=0.0, q_max=None,
                 local_depression=0.0, local_width=0.0):
        placement = placement.upper()
        if placement not in ("ON_SAG", "ON_GRADE", "AUTOMATIC"):
            raise ValueError(f"unknown placement {placement!r}")
        if placement != "ON_SAG" and not long_slope:
            raise ValueError(f"{placement} needs a longitudinal slope")
        self.design, self.street, self.placement = design, street, placement
        self.long_slope, self.number = long_slope, int(number)
        self.clog_factor = 1.0 - (clogged_pct or 0.0) / 100.0
        self.q_max = q_max or None
        self._local = (local_depression, local_width)

    def sag(self, depth):
        q = sag_capture(self.design, depth, self.street, *self._local) * self.clog_factor
        if self.q_max:
            q = min(q, self.q_max)
        sides = self.street.sides if self.street is not None else 1
        return q * sides * self.number

    def on_grade(self, approach_Q):
        return on_grade_capture(self.design, approach_Q, self.street, self.long_slope,
                                self.number, self.clog_factor, self.q_max, *self._local)

    def __call__(self, depth, approach_Q=0.0):
        if self.placement == "ON_SAG":
            return self.sag(depth)
        if self.placement == "ON_GRADE":
            return self.on_grade(approach_Q)
        return max(self.sag(depth), self.on_grade(approach_Q))


# --------------------------------------------------------------------------
# .inp [INLETS] / [STREETS]
# --------------------------------------------------------------------------
def read_inlet_designs(inp_path):
    """Parse a SWMM ``.inp``'s ``[INLETS]`` and ``[STREETS]`` (SI units) into
    ``({name: InletDesign}, {name: Street})``. A design listed with both a
    GRATE and a CURB line is a combination inlet, as in SWMM."""
    from .inp import _read_sections
    sec = _read_sections(inp_path)
    designs = {}
    for row in sec.get("INLETS", []):
        name, kind = row[0], row[1].upper()
        d = designs.setdefault(name, InletDesign(name))
        if kind in ("GRATE", "DROP_GRATE"):
            gtype = row[4].upper()
            d.grate = Grate(float(row[2]), float(row[3]), gtype,
                            float(row[5]) if gtype == "GENERIC" else None,
                            float(row[6]) if gtype == "GENERIC" and len(row) > 6 else None)
        elif kind in ("CURB", "DROP_CURB"):
            d.curb = Curb(float(row[2]), float(row[3]), row[4] if len(row) > 4 else "VERTICAL")
        else:
            raise NotImplementedError(f"inlet {name!r}: {kind} inlets are not supported")
    streets = {}
    for row in sec.get("STREETS", []):
        # Name Tcrown Hcurb Sx nRoad a W Sides Tback Sback nBack
        streets[row[0]] = Street(cross_slope=float(row[3]) / 100.0, roughness=float(row[4]),
                                 crown_width=float(row[1]),
                                 gutter_depression=float(row[5]) if len(row) > 5 else 0.0,
                                 gutter_width=float(row[6]) if len(row) > 6 else 0.0,
                                 sides=int(row[7]) if len(row) > 7 else 1)
    return designs, streets
