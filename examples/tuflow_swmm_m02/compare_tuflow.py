"""Compare the anuga_drainage runs with the TUFLOW-SWMM reference, pit by pit.

Reads TUFLOW's SWMM report (Node Inflow Summary: the lateral inflow each pit's
inlet delivered to the pipes) and every ``pits_<backend>_<capture>.csv`` that
``run_from_inp.py`` has written, and prints the volume each pit takes over
the event, plus the peaks. Then compares peak water levels with TUFLOW's
h_Max grid, over the wet area and at each pit, for every ``tuflow_m02_*.sww``
(peak over the SWW output steps, so slightly under the true peak). Levels,
not depths: with sub-grid sampling TUFLOW's d_Max is measured from the lowest
point in a 5 m cell, which no ANUGA triangle sees.

    ./run_tuflow.sh                      # reference -> tuflow_data/TUFLOW/results/
    python run_from_inp.py swmm          # -> pits_swmm_hec22.csv
    python compare_tuflow.py
"""
import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "tuflow_data" / "TUFLOW" / "results"
RPT = RESULTS / "TS02_5m_001_CPU_swmm.rpt"
H_MAX = RESULTS / "grids" / "TS02_5m_001_CPU_h_Max.tif"
WET = 0.05  # [m] depth that counts as wet


def tuflow_pits(rpt=RPT):
    """Per-pit peak lateral inflow [m^3/s] and lateral volume [m^3] from a
    SWMM report's Node Inflow Summary (volumes are reported in 10^6 L)."""
    text = Path(rpt).read_text()
    sec = text[text.index("Node Inflow Summary"):text.index("Node Surcharge Summary")]
    rows = re.findall(r"^\s+(\S+)\s+JUNCTION\s+([\d.]+)\s+[\d.]+\s+\d+\s+[\d:]+\s+([\d.]+)",
                      sec, re.M)
    df = pd.DataFrame(rows, columns=["pit", "peak_cms", "volume_m3"]).set_index("pit").astype(float)
    df["volume_m3"] *= 1000.0
    return df


def sww_peaks(sww):
    """Centroid x, y, bed and peak stage over the stored time steps."""
    import netCDF4
    with netCDF4.Dataset(sww) as d:
        d.set_auto_mask(False)
        tri = d["volumes"][:]
        x = d["x"][:][tri].mean(axis=1) + d.xllcorner
        y = d["y"][:][tri].mean(axis=1) + d.yllcorner
        z = d["elevation_c"][:]
        z = z[0] if z.ndim == 2 else z
        return x, y, z, d["stage_c"][:].max(axis=0)


def sample_grid(tif, x, y):
    """Grid values at points; NaN where the grid has no data (dry)."""
    import rasterio
    with rasterio.open(tif) as g:
        v = np.array([s[0] for s in g.sample(np.c_[x, y])], float)
        v[v == g.nodata] = np.nan
    return v


def compare_levels():
    import tuflow_inputs as ti
    pits = ti.pit_inlets()
    px = np.array([r.geometry.x for r in pits.values()])
    py = np.array([r.geometry.y for r in pits.values()])
    tuf_pit = sample_grid(H_MAX, px, py)
    area, at_pits = {}, {"TUFLOW": tuf_pit}
    for sww in sorted(glob.glob(str(HERE / "tuflow_m02_*.sww"))):
        name = Path(sww).stem[len("tuflow_m02_"):]
        x, y, z, h = sww_peaks(sww)
        tuf = sample_grid(H_MAX, x, y)
        ours_wet, tuf_wet = h - z > WET, np.isfinite(tuf)
        both = ours_wet & tuf_wet
        dh = (h - tuf)[both]
        area[name] = {"wet_both": both.sum(), "wet_only_ours": (ours_wet & ~tuf_wet).sum(),
                      "wet_only_tuflow": (tuf_wet & ~ours_wet).sum(),
                      "mean_dh_m": dh.mean(), "median_dh_m": np.median(dh),
                      "p90_abs_dh_m": np.percentile(np.abs(dh), 90)}
        k = np.hypot(x[:, None] - px, y[:, None] - py).argmin(axis=0)
        at_pits[name] = h[k]
    print(f"\nPeak water level vs TUFLOW h_Max (triangle counts; dh = ours - TUFLOW where both wet)")
    print(pd.DataFrame(area).T.round(3).to_string())
    lv = pd.DataFrame(at_pits, index=list(pits))
    lv = lv.loc[sorted(lv.index, key=lambda n: int(n[3:]))]
    print("\nPeak water level at each pit [m AHD]  (difference from TUFLOW)")
    print(lv.assign(**{c: lv[c] - lv["TUFLOW"] for c in lv.columns if c != "TUFLOW"}).round(3).to_string())


def main():
    tuf = tuflow_pits()
    runs = {Path(p).stem[len("pits_"):]: pd.read_csv(p).set_index("pit")
            for p in sorted(glob.glob(str(HERE / "pits_*.csv")))}
    if not runs:
        raise SystemExit("no pits_*.csv yet — run run_from_inp.py first")
    pits = sorted(next(iter(runs.values())).index, key=lambda n: int(n[3:]))

    vol = pd.DataFrame({"TUFLOW": tuf.loc[pits, "volume_m3"]})
    peak = pd.DataFrame({"TUFLOW": tuf.loc[pits, "peak_cms"]})
    for name, df in runs.items():
        vol[name] = df.loc[pits, "accepted_m3"]
        peak[name] = df.loc[pits, "peak_draw_cms"]
    vol.loc["total"] = vol.sum()

    pd.set_option("display.width", 200)
    print("Volume into the pipes per pit [m^3]  (ours: accepted, net of flooding)")
    print(vol.round(0).to_string())
    print("\nRatio to TUFLOW")
    print((vol.div(vol["TUFLOW"], axis=0)).drop(columns="TUFLOW").round(2).to_string())
    print("\nPeak inlet flow per pit [m^3/s]  (ours: the requested draw, before any flooding)")
    print(peak.round(3).to_string())
    if H_MAX.exists():
        compare_levels()


if __name__ == "__main__":
    main()
