# DEM Differencing Toolkit

**Version:** 1.0.0

## What is this repo?

This repository provides generalized scripts that enable accurate and interpretable differencing of digital elevation models (DEMs) with:
- Differing grids and pixel sizes
- Different horizontal and vertical coordinate reference systems
- Varying extents

**Future** It includes options for local and supercomputer systems, plus working examples with data to help users gain confidence with the workflow. 

---

## How does it work?

This repository contains (or will contain) three major components:

### 1. Basic DEM Differencing Script
A standalone script for differencing two DEMs on a local machine.

### 2. Supercomputer-Adapted Script
A version optimized for HPC systems (like UAF's Chinook) with accompanying `.sh` batch files.

### 3. Example Implementation
Working examples using provided data to demonstrate both scripts.

---

## Why create this?

The goal is to establish a **standard pipeline for DEM differencing**. 

I developed this toolkit because I needed to difference ~40 DEMs of one region with vastly different:
- Coordinate reference systems
- Spatial extents  
- Resolutions

While xDEM provides valuable standardization, my workflow required additional functionality not covered by their established pipelines.

**For GEOS 694:** I've successfully adapted this code into a class format (Task #1). Next steps include adapting for Chinook (Task #2) and implementing effective parallelization (Secondary Task #1).

---

## Repository Contents

| File | Description |
|------|-------------|
| `dem_diff.py` | Main script for differencing DEMs with different grids, CRS, and formats. Requires users to specify horizontal/vertical CRS and NoData values. |
| `requirements.txt` | Python package dependencies for running the scripts. |
| `demDiff.yml` | Conda environment file for creating a local environment. |

---

## Installation

```bash
# Clone the repository
git clone https://github.com/miajac/diffDEM.git
cd diffDEM

# Create conda environment
conda env create -f demDiff.yml
conda activate demDiff

# Or install dependencies manually
pip install -r requirements.txt
```

---

## Usage

### Command Line

Edit the parameters at the bottom of `dem_diff.py` in the `if __name__ == "__main__":` block:

```python
differencer = DEMDifferencer(
    path_dem1="path/to/dem1.tif",
    nickname_dem1="DEM1",
    src_vcrs_dem1="NAVD88",  # or "Ellipsoid", "EGM96", etc. 
    src_hcrs_dem1="EPSG:3338",
    nodata_dem1=-9999,
    path_dem2="path/to/dem2.tif",
    nickname_dem2="DEM2",
    src_vcrs_dem2="Ellipsoid",
    src_hcrs_dem2="EPSG:32606",
    nodata_dem2=-9999,
    path_dest="path/to/output/directory",
    coregister=True,  # Optional: set True for Nuth & Kääb coregistration
    roi=None  # Optional: provide vector for region of interest
)
```

Then run:

```bash
python dem_diff.py
```

The differenced DEM will be saved to: `path_dest/nickname_dem2_nickname_dem1.tif`

---

### As a Python Module (Optional)

You can also import and use the class in your own scripts:

```python
from dem_diff import DEMDifferencer

differencer = DEMDifferencer(
    path_dem1="path/to/dem1.tif",
    # ... (same parameters as above)
)
differencer.run()
```

---

## Acknowledgments

This toolkit builds heavily on [xDEM](https://github.com/GlacioHack/xdem), an open-source Python package for analyzing digital elevation models. 

**xDEM citation:**
> xDEM contributors. (2024). xDEM (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.11492983
