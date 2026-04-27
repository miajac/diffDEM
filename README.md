# DEM Differencing Toolkit

**Version:** 3.0.0

Version 3.0.0 introduces `dem_diff_batch.py`, an MPI-enabled batch script for differencing multiple DEM pairs across supercomputer nodes with SLURM job scheduling.

## What is this repo?

This repository provides generalized scripts that enable accurate and interpretable differencing of digital elevation models (DEMs) with:
- Differing grids and pixel sizes
- Different horizontal and vertical coordinate reference systems
- Varying extents

It includes options for singular, parallel, and supercomputer-based computing, plus a working example with data to help users gain confidence with the workflow.

---

## How does it work?

This repository contains four major components:

### 1. Basic DEM Differencing Script
A standalone script for differencing two DEMs on a local machine.

### 2. Parallel-Adapted Script
A version of the DEM differencing script optimized for larger datasets and 
multi-core machines. Key differences from the basic script include:

- **Sector-based differencing**: the DEM is split into a grid of sectors 
  that are differenced independently and then mosaicked back together,
  allowing the most computationally expensive step to run in parallel
- **Parallel DEM preparation**: CRS assignment, vertical datum conversion, 
  and reprojection of both DEMs run simultaneously rather than sequentially
- **Configurable workers and sectors**: the number of parallel workers and 
  sector grid size are controlled via command-line arguments 
  (`--workers`, `--num-sectors`)
- **Same config file format**: uses identical `.yml` parameter files as the 
  basic script, so no changes to your config are needed to switch between versions

Run with:
```bash
python dem_diff_parallel.py config_ifsar_lidar.yml --num-sectors 4 --workers auto
```

### 3. Supercomputer Batch Script
An MPI-enabled version of the DEM differencing script designed for 
supercomputer environments with SLURM job scheduling. Features include:

- **MPI task distribution**: each MPI task processes a different DEM pair 
  independently, so N pairs run simultaneously across N nodes
- **Flexible pair modes**: supports sequential, all-combinations, or 
  explicit user-defined pairs with optional per-pair coregistration overrides
- **SLURM integration**: submit via `sbatch`, with logging and error files 
  generated per job
- **Same config format**: uses a similar `.yml` parameter file as the other 
  scripts, extended to support multiple DEMs and pair definitions

Run with:
```bash
sbatch dem_diff_batch.sh
```

### 4. Example Implementation
Working examples using cropped DEMs covering a small area of the Canwell 
Glacier and Isabel Pass, Alaska. The toy dataset includes:

- `toy_IsabelIFSAR_2000.tif`: a 2000 IFSAR DSM in UTM Zone 6N (EPSG:32606) 
  with EGM96 vertical datum
- `toy_IFSAR_DTM_2010.tif`: a 2010 IFSAR DTM in NAD83 / Alaska Albers 
  (EPSG:3338) with NAVD88 vertical datum
- `toy_Lidar2025.tif`: a 2025 lidar DTM in UTM Zone 6N (EPSG:32606) with 
  ellipsoidal vertical datum

These DEMs have different resolutions, horizontal CRS, and vertical datums, 
making them a realistic test of the full pipeline. Run the toy example with:

```bash
python dem_diff.py config_toy.yml
# or
python dem_diff_parallel.py config_toy.yml --num-sectors 4 --workers auto
# or on a supercomputer
sbatch dem_diff_batch.sh
```

---

## Why create this?

The goal is to pull together various python packages to create a **standard pipeline for DEM differencing**. 

I developed this toolkit because I needed to difference ~40 DEMs of one region with vastly different:
- Coordinate reference systems
- Spatial extents  
- Resolutions

While xDEM provides valuable standardization, my workflow required additional functionality not covered by their established pipelines.

**For GEOS 694:** I've successfully adapted this code into a class format (Task #1), have created a parallelized script (Task #1), implemented a Parameter Input System (Task #2), and developed a supercomputer-compatible MPI batch script (Task #2).

---

## Repository Contents

| File | Description |
|------|-------------|
| `dem_diff.py` | Main script for differencing two DEMs with different grids, CRS, and formats. |
| `dem_diff_parallel.py` | Parallelized version of `dem_diff.py`, optimized for local machines with multiple cores. |
| `dem_diff_batch.py` | MPI-enabled batch script for differencing multiple DEM pairs across supercomputer nodes. |
| `dem_diff_batch.sh` | SLURM batch script for submitting `dem_diff_batch.py` to a supercomputer. |
| `config_template.yml` | Template config file with all parameters and descriptions. Copy and rename for each new run. |
| `config_toy.yml` | Pre-filled config file for running the two-DEM toy example. |
| `config_batch_template.yml` | Template config file for batch runs. Copy and rename for each new run. |
| `config_batch_toy.yml` | Pre-filled config file for running the batch toy example with three DEMs. |
| `demDiff.yml` | Conda environment file for creating a local environment. |
| `requirements.txt` | Python package dependencies for running the scripts. |
| `exampleData/` | Contains cropped toy DEMs for testing. |
| `customCanwell/` | Example config and parameter files for a real Canwell Glacier run. |

---

## Installation

### Local Machine

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

### Supercomputer

```bash
# Clone the repository
git clone https://github.com/miajac/diffDEM.git
cd diffDEM

# Set up conda environment
conda env create -f demDiff.yml
conda activate demDiff

# Install OpenMPI-compatible mpi4py
mamba install -c conda-forge mpi4py openmpi

# Copy the NAVD88 geoid grid to the pyproj data directory
# Download from https://cdn.proj.org/us_noaa_geoid09_ak.tif and place in:
python -c "import pyproj; print(pyproj.datadir.get_data_dir())"
```

---

## Usage

### Basic and Parallel Scripts

Parameters are passed via a YAML config file. Copy `config_template.yml`, 
rename it, and fill in your values.
Then run:

```bash
# Basic script
python dem_diff.py config_myrun.yml

# Parallel script
python dem_diff_parallel.py config_myrun.yml --num-sectors 4 --workers auto
```

### Batch Script (Supercomputer)

Copy `config_batch_template.yml`, rename it, and fill in your values.
Update `dem_diff_batch.sh` with the correct number of nodes (one per DEM pair), then submit:

```bash
sbatch dem_diff_batch.sh
```

Check job status with:
```bash
squeue -u yourusername
```
---

## Running the Examples

### Basic and Parallel Scripts (Two-DEM Toy Example)
Uses `toy_IFSAR_DTM_2010.tif` and `toy_Lidar2025.tif` (two DEMs of Canwell Glacier with different CRS and vertical datums).

```bash
# Basic script
python dem_diff.py config_toy.yml

# Parallel script
# add --workers auto at the end if you'd like to specify 
python dem_diff_parallel.py config_toy.yml --num-sectors 4 
```

### Batch Script (Three-DEM Toy Example)
Uses `toy_IsabelIFSAR_2000.tif`, `toy_IFSAR_DTM_2010.tif`, and `toy_Lidar2025.tif` (three DEMs of Canwell Glacier that run simultaneously across nodes). 

```bash
sbatch dem_diff_batch.sh
```

Output files will be saved to the directory specified in `path_dest` in your config file, named `nickname_dem2_nickname_dem1.tif`.

---

## Acknowledgments

This toolkit builds heavily on [xDEM](https://github.com/GlacioHack/xdem), an open-source Python package for analyzing digital elevation models. 

**xDEM citation:**
> xDEM contributors. (2024). xDEM (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.11492983
