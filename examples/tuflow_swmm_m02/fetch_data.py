"""Download the TUFLOW-SWMM tutorial dataset and extract Module 2's complete model.

The TUFLOW tutorial models are free to download and run, but carry no licence for
redistributing the data, so they are fetched on demand rather than committed:

    python fetch_data.py            # -> tuflow_data/TUFLOW/{model,bc_dbase,runs}

The zip (~100 MB, every module plus a Windows TUFLOW executable) is cached in
``tuflow_data/`` so re-running only re-extracts. ``--exe`` also extracts that
TUFLOW 2026.0.0 build to ``tuflow_data/exe/`` (see run_tuflow.sh), which runs
this model licence-free (``Tutorial Model == ON``).
"""
import shutil
import sys
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath

URL = ("https://downloads.tuflow.com/TUFLOW/Wiki_Tute_Models/"
       "TUFLOW_SWMM_Tutorial_Models_QGIS_GPKG.zip")
ROOT = "TUFLOW_SWMM_Tutorial_Models_QGIS_GPKG/"
MODEL_PREFIX = ROOT + "TUFLOW_SWMM_Module_02/Complete_Model/"
EXE_PREFIX = ROOT + "exe/"

HERE = Path(__file__).resolve().parent
DATA = HERE / "tuflow_data"


def download(zip_path):
    if zip_path.exists():
        print(f"Using cached {zip_path.name}")
        return
    print(f"Downloading {URL}")
    tmp = zip_path.with_suffix(".part")
    with urllib.request.urlopen(URL) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f)
    tmp.rename(zip_path)


def extract(zip_path, prefix, dest):
    """Extract the members under ``prefix``, stripping it."""
    n = 0
    with zipfile.ZipFile(zip_path) as z:
        for info in z.infolist():
            if not info.filename.startswith(prefix) or info.is_dir():
                continue
            rel = PurePosixPath(info.filename[len(prefix):])
            if rel.is_absolute() or ".." in rel.parts:      # never write outside dest
                raise ValueError(f"unsafe path in archive: {info.filename}")
            target = dest.joinpath(*rel.parts)
            target.parent.mkdir(parents=True, exist_ok=True)
            with z.open(info) as src, open(target, "wb") as out:
                shutil.copyfileobj(src, out)
            n += 1
    if n == 0:
        raise RuntimeError(f"no files under {prefix!r} in {zip_path}")
    print(f"Extracted {n} files to {dest}")


def main():
    DATA.mkdir(exist_ok=True)
    zip_path = DATA / Path(URL).name
    download(zip_path)
    extract(zip_path, MODEL_PREFIX, DATA)
    if "--exe" in sys.argv[1:]:
        extract(zip_path, EXE_PREFIX, DATA / "exe")
    return 0


if __name__ == "__main__":
    sys.exit(main())
