#!/usr/bin/env bash
# Run the TUFLOW-SWMM Module 2 model (TS02_5m_001.tcf) with the TUFLOW 2026.0.0
# Windows build bundled in the tutorial zip, under Wine. Licence-free: the .tcf
# sets "Tutorial Model == ON". Produces the reference results the anuga_drainage
# run (run_from_inp.py) is compared against, in tuflow_data/TUFLOW/results/.
#
#   python fetch_data.py --exe      # once: model + TUFLOW executable
#   ./run_tuflow.sh [cpu|gpu] [TUFLOW flags...]   # default: cpu
#
# The .tcf asks for "Hardware == GPU"; CUDA through Wine is hit-and-miss, so
# by default a copy of the .tcf with "Hardware == CPU" is run instead.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
runs="$here/tuflow_data/TUFLOW/runs"
exe="$here/tuflow_data/exe/2026.0.0/TUFLOW_iSP_w64.exe"
hardware="${1:-cpu}"
shift || true

[[ -f "$exe" ]] || { echo "missing $exe — run: python fetch_data.py --exe" >&2; exit 1; }
command -v wine >/dev/null || { echo "wine not found — install it: sudo apt install wine" >&2; exit 1; }

cd "$runs"
case "$hardware" in
  cpu) tcf=TS02_5m_001_CPU.tcf
       sed 's/^Hardware == GPU/Hardware == CPU/' TS02_5m_001.tcf > "$tcf" ;;
  gpu) tcf=TS02_5m_001.tcf ;;
  *)   echo "usage: $0 [cpu|gpu] [TUFLOW flags...]" >&2; exit 2 ;;
esac

# Wine workarounds (Wine 10, Ubuntu 26.04):
# - Intel's OpenMP (libiomp5md.dll) calls GetNumaNodeProcessorMaskEx, which Wine
#   doesn't implement; KMP_AFFINITY=disabled skips that topology detection.
# - The HPC solver uses Microsoft's OpenMP (VCOMP140.DLL, not shipped with
#   TUFLOW). Wine's builtin vcomp140 ignores the requested thread count, so
#   HPC fails with "Wanted number of threads = 4, Actual = <all CPUs>". Install
#   Microsoft's VC++ 2015-2022 x64 redistributable into the Wine prefix
#   (wine vc_redist.x64.exe /install /quiet /norestart) and prefer the native DLL.
# - TUFLOW's Fortran runtime writes to the Windows console (CONOUT$), which
#   only exists when Wine has a terminal: without one, run under `script`.
export KMP_AFFINITY=disabled
export WINEDLLOVERRIDES="vcomp140=n,b${WINEDLLOVERRIDES:+;$WINEDLLOVERRIDES}"
export WINEDEBUG=-all
# -b: batch mode (no "press any key" at the end), as in the tutorial's .bat.
cmd=(wine "$exe" -b "$@" "$tcf")
if [[ -t 1 ]]; then
  "${cmd[@]}"
else
  script -qefc "$(printf '%q ' "${cmd[@]}")" /dev/null
fi
