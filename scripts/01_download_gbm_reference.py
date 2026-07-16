"""
01_download_gbmap_reference.py
------------------------------
Downloads the GBmap harmonized glioblastoma scRNA-seq reference from NCBI GEO
(GSE211376, Ruiz-Moreno et al. 2022) and saves it locally for use in the
human GBM extension pipeline.

Outputs (saved to REF_DIR):
    GBmap_raw_counts.tsv.gz   — raw count matrix (27,102 genes x 39,355 cells)
    GBmap_metadata.csv.gz     — cell metadata including predicted.high_hierarchy
                                (cell type labels used by cell2location)

Run this once before running any GBM extension scripts.

Usage:
    python scripts/01_download_gbmap_reference.py
    python scripts/01_download_gbmap_reference.py --ref-dir /path/to/custom/dir
"""

import os
import sys
import gzip
import hashlib
import argparse
import urllib.request

# ── Argument parsing ──────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Download GBmap reference from NCBI GEO")
parser.add_argument(
    "--ref-dir",
    default=os.path.join(os.path.dirname(__file__), "..", "data", "GBmap_reference"),
    help="Directory to save reference files (default: data/GBmap_reference/)"
)
parser.add_argument(
    "--force",
    action="store_true",
    help="Re-download even if files already exist"
)
args = parser.parse_args()

REF_DIR = os.path.abspath(args.ref_dir)
os.makedirs(REF_DIR, exist_ok=True)

# ── Download URLs ─────────────────────────────────────────────────────────────
# Direct links from NCBI GEO supplementary files for GSE211376
URLS = {
    "counts": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211376/suppl/"
        "GSE211376_raw_counts_Ruiz2022_all_samples_filtered_cells.tsv.gz"
    ),
    "metadata": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211376/suppl/"
        "GSE211376_metadata_Ruiz2022_all_samples_filtered_cells.csv.gz"
    ),
}

OUTPUT_FILES = {
    "counts":   os.path.join(REF_DIR, "GBmap_raw_counts.tsv.gz"),
    "metadata": os.path.join(REF_DIR, "GBmap_metadata.csv.gz"),
}


def md5(filepath):
    h = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            h.update(chunk)
    return h.hexdigest()


def download(name, url, output_path, force=False):
    if os.path.exists(output_path) and not force:
        size_mb = os.path.getsize(output_path) / 1024 / 1024
        print(f"  Already exists ({size_mb:.1f} MB) — skipping. Use --force to re-download.")
        return True

    print(f"  Downloading {name}...")
    print(f"  URL: {url}")
    try:
        urllib.request.urlretrieve(url, output_path)
        size_mb = os.path.getsize(output_path) / 1024 / 1024
        print(f"  Done: {size_mb:.1f} MB  |  MD5: {md5(output_path)}")
        return True
    except Exception as e:
        print(f"  FAILED: {e}")
        return False


def verify(name, filepath):
    if not os.path.exists(filepath):
        print(f"  MISSING: {filepath}")
        return False
    # Try reading first line to confirm file is not corrupted
    try:
        with gzip.open(filepath, "rt") as f:
            first = f.readline().strip()
        print(f"  OK ({os.path.getsize(filepath)/1024/1024:.1f} MB) — first line: {first[:80]}...")
        return True
    except Exception as e:
        print(f"  CORRUPTED or unreadable: {e}")
        return False


# ── Main ──────────────────────────────────────────────────────────────────────
print("=" * 60)
print("DOWNLOADING GBMAP REFERENCE (GSE211376)")
print(f"Saving to: {REF_DIR}")
print("=" * 60)

success = True
for name, url in URLS.items():
    print(f"\n[{name}]")
    ok = download(name, url, OUTPUT_FILES[name], force=args.force)
    if not ok:
        success = False

print("\n" + "=" * 60)
print("VERIFYING FILES")
print("=" * 60)

for name, filepath in OUTPUT_FILES.items():
    print(f"\n[{name}]")
    verify(name, filepath)

print("\n" + "=" * 60)
if success:
    print("Download complete. Run 02_prepare_gbmap_reference.py next.")
else:
    print("One or more downloads FAILED. Check your internet connection and retry.")
    sys.exit(1)
print("=" * 60)
