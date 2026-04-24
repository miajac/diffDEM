# DEM Differencing Toolkit

**Version:** 2.0.0

Version 2.0.0 introduces `dem_diff_parallel.py`, a parallelized adaptation of the original `dem_diff.py` script optimized for larger datasets and multi-core machines.

## What is this repo?

This repository provides generalized scripts that enable accurate and interpretable differencing of digital elevation models (DEMs) with:
- Differing grids and pixel sizes
- Different horizontal and vertical coordinate reference systems
- Varying extents

It includes options for parallel and singular computing, plus a working example with data to help users gain confidence with the workflow.

---

## How does it work?

This repository contains three major components:

### 1. Basic DEM Differencing Script
A standalone script for differencing two DEMs on a local machine.

### 2. Parallel-Adapted Script
A version of the DEM differencing script optimized for larger datasets and 
multi-core machines. Key differences from the basic script include:

- **Sector-based differencing** — the DEM is split into a grid of sectors 
  that are differenced independently and then mosaicked back together,
  allowing the most computationally expensive step to run in parallel
- **Parallel DEM preparation** — CRS assignment, vertical datum conversion, 
  and reprojection of both DEMs run simultaneously rather than sequentially
- **Configurable workers and sectors** — the number of parallel workers and 
  sector grid size are controlled via command-line arguments 
  (`--workers`, `--num-sectors`)
- **Same config file format** — uses identical `.yml` parameter files as the 
  basic script, so no changes to your config are needed to switch between versions

Run with:
```bash
python dem_diff_parallel.py config_ifsar_lidar.yml --num-sectors 4 --workers auto
```

### 3. Example Implementation
Working examples using two cropped DEMs covering a small ~200 m x 200 m area 
of the Canwell Glacier, Alaska. The toy dataset includes:

- `toy_IFSAR_DTM_2010.tif` — a 2010 IFSAR DTM in NAD83 / Alaska Albers 
  (EPSG:3338) with NAVD88 vertical datum
- `toy_Lidar2025.tif` — a 2025 lidar DTM in UTM Zone 6N (EPSG:32606) with 
  ellipsoidal vertical datum

These two DEMs have different resolutions, horizontal CRS, and vertical datums, 
making them a realistic test of the full pipeline. Run the toy example with:

```bash
python dem_diff.py config_toy.yml
# or
python dem_diff_parallel.py config_toy.yml --num-sectors 4 --workers auto
```

---

## Why create this?

The goal is to pull together various python packages to create a **standard pipeline for DEM differencing**. 

I developed this toolkit because I needed to difference ~40 DEMs of one region with vastly different:
- Coordinate reference systems
- Spatial extents  
- Resolutions

While xDEM provides valuable standardization, my workflow required additional functionality not covered by their established pipelines.

**For GEOS 694:** I've successfully adapted this code into a class format (Task #1), have created parallelized script (Task #1), and have implemented a Parameter Input System (Task #2).

---

## Repository Contents

| File | Description |
|------|-------------|
| `dem_diff.py` | Main script for differencing two DEMs with different grids, CRS, and formats. |
| `dem_diff_parallel.py` | Parallelized version of `dem_diff.py`, optimized for local machines with multiple cores. |
| `config_template.yml` | Template config file with all parameters and descriptions. Copy and rename for each new run. |
| `config_toy.yml` | Pre-filled config file for running the toy example. |
| `demDiff.yml` | Conda environment file for creating a local environment. |
| `requirements.txt` | Python package dependencies for running the scripts. |
| `exampleData/` | Contains two cropped toy DEMs for testing (`toy_IFSAR_DTM_2010.tif`, `toy_Lidar2025.tif`). |

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

Parameters are passed via a YAML config file. Copy `config_template.yml`, 
rename it, and fill in your values:

```yaml
dem1:
  path: "/path/to/your/first_dem.tif"
  nickname: "ShortName1"
  src_vcrs: "EGM96"        # options: "Ellipsoid", "EGM96", "NAVD88"
  src_hcrs: "EPSG:4326"    # any valid EPSG code, see https://epsg.io/
  nodata: -9999

dem2:
  path: "/path/to/your/second_dem.tif"
  nickname: "ShortName2"
  src_vcrs: "Ellipsoid"
  src_hcrs: "EPSG:32606"
  nodata: -9999

options:
  path_dest: "/path/to/output/folder/"  # trailing slash required
  roi: null                              # path to vector file, or null
  coregister: false                      # true or false
```

Then run:

```bash
# Basic script
python dem_diff.py config_myrun.yml

# Parallel script
python dem_diff_parallel.py config_myrun.yml --num-sectors 4 --workers auto
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
