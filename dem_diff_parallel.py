"""
dem_diff_parallel.py — Parallelized DEM differencing with sector-based processing.

This version parallelizes:
1. DEM preparation (loading, CRS assignment, vertical conversion, reprojection)
2. Differencing via sector-based processing (split, difference in parallel, mosaic)

Usage:
    $ python dem_diff_parallel.py config.yml [--num-sectors 10] [--workers auto]

Where config.yml is a parameter file (same format as dem_diff.py).
"""

import os
import sys
import yaml
import xdem
import pyproj
import geoutils as gu
import numpy as np
import tempfile
import shutil
from pyproj.transformer import TransformerGroup
from multiprocessing import Pool, cpu_count
from functools import partial
from pathlib import Path


class DEMDifferencerParallel:
    """
    Parallelized DEM differencer with sector-based processing.
    
    Parameters: (same as DEMDifferencer, plus additional parallelization params)
    num_sectors : int, optional
        Number of sectors along each axis (default 10 → 10x10=100 sectors).
    num_workers : int, optional
        Number of worker processes (default: cpu_count()).
    """

    TARGET_HCRS = "EPSG:32606"
    TARGET_VCRS = "EGM96"
    NAVD88_GRID = "us_noaa_geoid09_ak.tif"

    def __init__(
        self,
        path_dem1,
        path_dem2,
        path_dest,
        nickname_dem1,
        nickname_dem2,
        src_vcrs_dem1,
        src_hcrs_dem1,
        nodata_dem1,
        src_vcrs_dem2,
        src_hcrs_dem2,
        nodata_dem2,
        roi=None,
        coregister=False,
        num_sectors=10,
        num_workers=None,
    ):
        self.path_dem1 = path_dem1
        self.path_dem2 = path_dem2
        self.path_dest = path_dest
        self.nickname_dem1 = nickname_dem1
        self.nickname_dem2 = nickname_dem2
        self.src_vcrs_dem1 = src_vcrs_dem1
        self.src_hcrs_dem1 = src_hcrs_dem1
        self.nodata_dem1 = nodata_dem1
        self.src_vcrs_dem2 = src_vcrs_dem2
        self.src_hcrs_dem2 = src_hcrs_dem2
        self.nodata_dem2 = nodata_dem2
        self.roi = roi
        self.coregister = coregister
        self.num_sectors = num_sectors
        self.num_workers = num_workers or cpu_count()

        self.output_path = (
            f"{self.path_dest}{self.nickname_dem2}_{self.nickname_dem1}.tif"
        )

        # Temp directory for sector files
        self.temp_dir = None

        self.dem1 = None
        self.dem2 = None
        self.diff_dem = None

    def _check_grids(self):
        """Verify that required projection grid files are available."""
        needs_navd88 = any(
            vcrs in ("NAVD88", "EPSG:5703")
            for vcrs in [self.src_vcrs_dem1, self.src_vcrs_dem2]
        )

        if needs_navd88:
            data_dir = pyproj.datadir.get_data_dir()
            grid_path = os.path.join(data_dir, self.NAVD88_GRID)
            if os.path.exists(grid_path):
                print(
                    f"[grid check] {self.NAVD88_GRID} found locally at {grid_path}"
                )
            else:
                print(
                    f"[grid check] {self.NAVD88_GRID} not found locally, "
                    f"attempting download from cdn.proj.org..."
                )
                try:
                    tg = TransformerGroup("EPSG:4269", "EPSG:5703")
                    tg.download_grids(verbose=True)
                    if os.path.exists(grid_path):
                        print(f"[grid check] download successful: {grid_path}")
                    else:
                        raise FileNotFoundError(
                            f"Grid {self.NAVD88_GRID} could not be downloaded "
                            f"to {data_dir}. Download it manually from "
                            f"https://cdn.proj.org/{self.NAVD88_GRID} "
                            f"and place it in {data_dir}."
                        )
                except Exception as e:
                    raise RuntimeError(
                        f"Failed to download required PROJ grid "
                        f"{self.NAVD88_GRID}. Check your internet connection "
                        f"or download manually from https://cdn.proj.org/"
                        f"{self.NAVD88_GRID} and place in {data_dir}."
                    ) from e

    def load(self):
        """Load both DEMs from file."""
        self.dem1 = xdem.DEM(self.path_dem1, nodata=self.nodata_dem1)
        self.dem2 = xdem.DEM(self.path_dem2, nodata=self.nodata_dem2)

        print("Before corrections:")
        print(f"{self.nickname_dem1} source file: {self.path_dem1}")
        print(f"{self.nickname_dem1} assigned nodata: {self.nodata_dem1}")
        print(f"{self.nickname_dem1} source hCRS: {self.src_hcrs_dem1}")
        print(f"{self.nickname_dem1} source vCRS: {self.src_vcrs_dem1}")
        print(
            f"{self.nickname_dem1} -> will be converted to hCRS: "
            f"{self.TARGET_HCRS}, vCRS: {self.TARGET_VCRS}."
        )
        print()
        print(f"{self.nickname_dem2} source file: {self.path_dem2}")
        print(f"{self.nickname_dem2} assigned nodata: {self.nodata_dem2}")
        print(f"{self.nickname_dem2} source hCRS: {self.src_hcrs_dem2}")
        print(f"{self.nickname_dem2} source vCRS: {self.src_vcrs_dem2}")
        print(
            f"{self.nickname_dem2} -> will be converted to hCRS: "
            f"{self.TARGET_HCRS}, vCRS: {self.TARGET_VCRS}."
        )
        print()
        print(f"{self.nickname_dem1} raster info:")
        self.dem1.info()
        print(f"{self.nickname_dem2} raster info:")
        self.dem2.info()

    def _prepare_dem(self, dem, src_hcrs, src_vcrs, nodata, nickname, src_path):
        """Prepare a single DEM (CRS assignment, vertical conversion, reprojection)."""
        if dem.data.dtype in [np.uint8, np.uint16, np.uint32]:
            print(
                f"[{nickname}] converting dtype {dem.data.dtype} -> float32 to "
                f"support negative values."
            )
            dem = dem.astype(np.float32)

        dem.crs = src_hcrs
        print(f"[{nickname}] hCRS assigned: {src_hcrs}")

        if src_vcrs in ("EPSG:5703", "NAVD88"):
            print(f"[{nickname}] vCRS is NAVD88 — using two-step conversion:")
            print(
                f"[{nickname}] Step 1: NAVD88 -> Ellipsoid (inverse of "
                f"{self.NAVD88_GRID})"
            )
            dem.set_vcrs(self.NAVD88_GRID)
            dem.to_vcrs("Ellipsoid")
            print(f"[{nickname}] Step 2: Ellipsoid -> {self.TARGET_VCRS}")
            dem.set_vcrs("Ellipsoid")
            dem.to_vcrs(self.TARGET_VCRS)
        else:
            dem.set_vcrs(src_vcrs)
            print(f"[{nickname}] vCRS assigned: {src_vcrs}")
            dem.to_vcrs(self.TARGET_VCRS)
            print(
                f"[{nickname}] vCRS converted: {src_vcrs} -> {self.TARGET_VCRS}"
            )

        if dem.nodata != nodata:
            dem.nodata = nodata
            print(f"[{nickname}] nodata set: {nodata}")
        else:
            print(f"[{nickname}] nodata already correct: {nodata}, skipping")

        if dem.crs.to_epsg() != int(self.TARGET_HCRS.split(":")[1]):
            dem = dem.reproject(crs=self.TARGET_HCRS, resampling="bilinear")
            print(
                f"[{nickname}] hCRS reprojected: {src_hcrs} -> {self.TARGET_HCRS}"
            )
        else:
            print(
                f"[{nickname}] hCRS already {self.TARGET_HCRS}, skipping "
                f"reprojection"
            )

        dem.set_vcrs(self.TARGET_VCRS)
        print(f"[{nickname}] vCRS re-stamped after reproject: {self.TARGET_VCRS}")
        source_info = (
            dem.filename if dem.filename
            else f"in-memory (reprojected from {src_path})"
        )
        print(f"[{nickname}] source file: {source_info}")
        return dem

    def prepare_parallel(self):
        """Parallelize preparation of both DEMs using multiprocessing."""
        print("\nPreparing DEMs in parallel:")
        
        with Pool(processes=min(2, self.num_workers)) as pool:
            results = pool.starmap(
                self._prepare_dem,
                [
                    (
                        self.dem1,
                        self.src_hcrs_dem1,
                        self.src_vcrs_dem1,
                        self.nodata_dem1,
                        self.nickname_dem1,
                        self.path_dem1,
                    ),
                    (
                        self.dem2,
                        self.src_hcrs_dem2,
                        self.src_vcrs_dem2,
                        self.nodata_dem2,
                        self.nickname_dem2,
                        self.path_dem2,
                    ),
                ],
            )
            self.dem1, self.dem2 = results

        print(f"\nAfter vertical conversion and horizontal reprojection:")
        print(f"{self.nickname_dem1}:")
        self.dem1.info()
        print(f"{self.nickname_dem2}:")
        self.dem2.info()

    def align(self):
        """Align DEMs to a common grid."""
        if self.roi is not None:
            self.dem1 = self.dem1.crop(self.roi)
            self.dem2 = self.dem2.crop(self.roi)
            print("\nAfter ROI crop:")
            print(f"{self.nickname_dem1}:")
            self.dem1.info()
            print(f"{self.nickname_dem2}:")
            self.dem2.info()

        res1 = self.dem1.res[0]
        res2 = self.dem2.res[0]

        print(f"\n[{self.nickname_dem1}] pixel size: {res1}m")
        print(f"[{self.nickname_dem2}] pixel size: {res2}m")

        if res1 >= res2:
            print(
                f"{self.nickname_dem1} is coarser ({res1}m >= {res2}m) — "
                f"reprojecting {self.nickname_dem2} to match"
            )
            self.dem2 = self.dem2.reproject(self.dem1, resampling="bilinear")
            self.dem2.set_vcrs(self.TARGET_VCRS)
            self.dem1.set_vcrs(self.TARGET_VCRS)
            print(f"Reference grid: {self.nickname_dem1}")
        else:
            print(
                f"{self.nickname_dem2} is coarser ({res2}m > {res1}m) — "
                f"reprojecting {self.nickname_dem1} to match"
            )
            self.dem1 = self.dem1.reproject(self.dem2, resampling="bilinear")
            self.dem1.set_vcrs(self.TARGET_VCRS)
            self.dem2.set_vcrs(self.TARGET_VCRS)
            print(f"Reference grid: {self.nickname_dem2}")

        print("\nAfter full alignment:")
        print(f"{self.nickname_dem1}:")
        self.dem1.info()
        print(f"{self.nickname_dem2}:")
        self.dem2.info()

    def _coregister(self):
        """Apply Nuth & Kääb coregistration."""
        print("\nCoregistering DEMs (Nuth & Kääb)...")
        nuth_kaab = xdem.coreg.NuthKaab()
        nuth_kaab.fit(self.dem2, self.dem1)

        shift_x = nuth_kaab.meta["outputs"]["affine"]["shift_x"]
        shift_y = nuth_kaab.meta["outputs"]["affine"]["shift_y"]
        shift_z = nuth_kaab.meta["outputs"]["affine"]["shift_z"]

        print(f"Detected shift_x: {shift_x:.3f}m")
        print(f"Detected shift_y: {shift_y:.3f}m")
        print(f"Detected shift_z: {shift_z:.3f}m")

        if abs(shift_z) > 2:
            print(f"\n*** WARNING: shift_z of {shift_z:.2f}m suggests a ")
            print(f"*** residual vertical datum offset.")
            print(f"*** Coregistration is correcting for this automatically,")
            print(f"*** but the root cause should ideally be fixed by ")
            print(f"*** verifying the src_vcrs parameter for each DEM is ")
            print(f"*** correct. Common cause in Alaska: IFSAR data labeled ")
            print(f"*** as EGM96 but actually stored on the WGS84 ellipsoid ")
            print(f"*** ~14-16m offset.")

        self.dem1 = nuth_kaab.apply(self.dem1)
        self.dem1 = self.dem1.reproject(self.dem2, resampling="bilinear")
        self.dem1.set_vcrs(self.TARGET_VCRS)

        print(f"\nAfter coregistration:")
        print(f"{self.nickname_dem1}:")
        self.dem1.info()

    def _define_sectors(self):
        """
        Define sectors based on pixel indices.
        Returns list of tuples: (sector_id, y_start, y_end, x_start, x_end)
        """
        height, width = self.dem1.shape
        sector_height = height // self.num_sectors
        sector_width = width // self.num_sectors

        sectors = []
        sector_id = 0
        for i in range(self.num_sectors):
            for j in range(self.num_sectors):
                y_start = i * sector_height
                x_start = j * sector_width
                y_end = min((i + 1) * sector_height, height)
                x_end = min((j + 1) * sector_width, width)
                sectors.append((sector_id, y_start, y_end, x_start, x_end))
                sector_id += 1

        print(f"\nDefined {len(sectors)} sectors ({self.num_sectors}x{self.num_sectors})")
        print(f"Sector dimensions: ~{sector_height}x{sector_width} pixels")
        return sectors

    def _process_sector(self, sector_info):
        """
        Process a single sector: extract, difference, save.
        
        Parameters:
        sector_info : tuple
            (sector_id, y_start, y_end, x_start, x_end)
        
        Returns:
        str : path to temporary sector file
        """
        sector_id, y_start, y_end, x_start, x_end = sector_info

        try:
            # Crop sectors using spatial bounds instead of array slicing
            # so the result stays an xdem.DEM object
            transform = self.dem1.transform
            left  = transform.c + x_start * transform.a
            top   = transform.f + y_start * transform.e
            right = transform.c + x_end   * transform.a
            bottom= transform.f + y_end   * transform.e

            dem1_sector = self.dem1.crop(
                (left, bottom, right, top), inplace=False
            )
            dem2_sector = self.dem2.crop(
                (left, bottom, right, top), inplace=False
            )

            # Difference — result is now an xdem.DEM
            diff_sector = dem2_sector - dem1_sector

            # Preserve CRS and nodata metadata
            diff_sector.set_vcrs(self.TARGET_VCRS)
            diff_sector.nodata = self.dem2.nodata

            # Save to temp file with georeferencing
            temp_path = os.path.join(
                self.temp_dir, f"sector_{sector_id:05d}.tif"
            )
            diff_sector.to_file(temp_path)

            print(f"[Sector {sector_id:05d}] processed and saved to {temp_path}")
            return temp_path

        except Exception as e:
            print(f"[Sector {sector_id:05d}] ERROR: {e}")
            raise

        
    def difference_sectors_parallel(self):
        """
        Difference DEMs by splitting into sectors and processing in parallel.
        """
        print("\nDifferencing DEMs via parallel sector processing:")

        # Create temp directory
        self.temp_dir = tempfile.mkdtemp(prefix="dem_diff_sectors_")
        print(f"Created temp directory: {self.temp_dir}")

        # Define sectors
        sectors = self._define_sectors()

        # Process sectors in parallel
        print(f"Processing {len(sectors)} sectors with {self.num_workers} workers...")
        with Pool(processes=self.num_workers) as pool:
            sector_files = pool.map(self._process_sector, sectors)

        print(f"\nAll {len(sector_files)} sectors processed successfully")

        # Mosaic sector files back into single DEM
        print("\nMosaicking sectors into final DEM...")
        self.diff_dem = self._mosaic_sectors(sector_files)
        print("Mosaicking complete")

        # Clean up temp directory
        shutil.rmtree(self.temp_dir)
        print(f"Cleaned up temp directory: {self.temp_dir}")

    def _mosaic_sectors(self, sector_files):
        """
        Mosaic sector files into a single DEM.
        
        Uses geoutils to read and combine sectors while preserving
        georeferencing.
        """
        # Load all sector files
        sectors_data = [gu.Raster(f) for f in sector_files]

        # Stack and merge (geoutils handles georeferencing)
        # Since sectors are non-overlapping and adjacent, we can use simple
        # concatenation along spatial dimensions or use rasterio.merge equivalent
        mosaicked = self._merge_rasters(sectors_data)
        return mosaicked

    def _merge_rasters(self, raster_list):
        """
        Merge a list of non-overlapping geoutils Rasters into one.
        Assumes rasters are georeferenced and non-overlapping.
        """
        # Simple approach: use numpy concatenation with georeferencing
        # For more robustness, consider rasterio.merge.merge
        if len(raster_list) == 1:
            return raster_list[0]

        # Load data and build merged array
        # (This is simplified; production code should use rasterio for robustness)
        try:
            # Try using rasterio for proper merging
            import rasterio
            from rasterio.merge import merge

            with rasterio.open(raster_list[0].filename) as src:
                profile = src.profile

            # Merge all raster files
            merged_data, merged_transform = merge(
                [r.filename for r in raster_list]
            )

            # Create a new geoutils Raster from merged data
            merged_raster = gu.Raster(
                merged_data,
                transform=merged_transform,
                crs=raster_list[0].crs,
                nodata=raster_list[0].nodata,
            )
            return merged_raster

        except ImportError:
            print(
                "Warning: rasterio not available, using basic concatenation. "
                "Install rasterio for robust mosaicking."
            )
            # Fallback: basic stacking (assumes grid-aligned sectors)
            # This is a simplified fallback and may not handle edge cases
            all_data = np.concatenate([r.data for r in raster_list], axis=0)
            merged_raster = gu.Raster(
                all_data,
                transform=raster_list[0].transform,
                crs=raster_list[0].crs,
                nodata=raster_list[0].nodata,
            )
            return merged_raster

    def check_stable_terrain(self):
        """Check elevation differences on flat stable terrain."""
        print("\nStable terrain check:")
        slope = xdem.terrain.slope(self.diff_dem)
        flat_mask = (slope.data < 5) & (~np.ma.getmaskarray(self.diff_dem.data))
        diff_flat = self.diff_dem.data.data[flat_mask]

        mean_offset = np.mean(diff_flat)
        median_offset = np.median(diff_flat)
        std_offset = np.std(diff_flat)

        print(f"Flat terrain pixel count: {len(diff_flat)}")
        print(f"Mean offset: {mean_offset:.2f}m")
        print(f"Median offset: {median_offset:.2f}m")
        print(f"Std dev: {std_offset:.2f}m")

        if abs(median_offset) > 2:
            print(
                f"\n*** WARNING: median offset of {median_offset:.2f}m on flat "
                f"stable terrain."
            )
            print(f"*** suggests a residual vertical datum mismatch.")
            if 10 < abs(median_offset) < 20:
                print(f"*** Offset of ~14-15m in Alaska typically indicates ")
                print(f"*** the DEM is on the WGS84 ellipsoid rather than ")
                print(f"*** EGM96. Try changing src_vcrs from 'EGM96' to")
                print(f"*** 'Ellipsoid' and rerun.")
        else:
            print(
                f"Offset within acceptable range — vertical datums appear "
                f"consistent."
            )

    def save(self):
        """Save the differenced DEM to the output path."""
        self.diff_dem.to_file(self.output_path)
        print(f"\nDifferenced DEM saved to: {self.output_path}")

    def run(self):
        """Run the full parallelized pipeline."""
        self._check_grids()
        self.load()
        self.prepare_parallel()  # Parallel preparation
        self.align()
        if self.coregister:
            self._coregister()
        self.difference_sectors_parallel()  # Parallel sector differencing
        self.check_stable_terrain()
        self.save()


def load_config(config_path):
    """Load a YAML config file."""
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    required = {
        "dem1": ["path", "nickname", "src_vcrs", "src_hcrs", "nodata"],
        "dem2": ["path", "nickname", "src_vcrs", "src_hcrs", "nodata"],
        "options": ["path_dest"],
    }
    for section, fields in required.items():
        if section not in cfg:
            raise ValueError(f"Config missing required section: '{section}'")
        for field in fields:
            if field not in cfg[section]:
                raise ValueError(
                    f"Config section '{section}' missing required field: "
                    f"'{field}'"
                )
    # Expand ~ to full home directory path
    cfg["dem1"]["path"] = os.path.expanduser(cfg["dem1"]["path"])
    cfg["dem2"]["path"] = os.path.expanduser(cfg["dem2"]["path"])
    cfg["options"]["path_dest"] = os.path.expanduser(
        cfg["options"]["path_dest"]
    )

    return cfg


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Parallelized DEM differencing with sector-based processing"
    )
    parser.add_argument("config", help="Path to config YAML file")
    parser.add_argument(
        "--num-sectors",
        type=int,
        default=10,
        help="Number of sectors per axis (default: 10 → 10x10 grid)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of worker processes (default: auto-detect CPU count)",
    )

    args = parser.parse_args()

    config_path = args.config
    if not os.path.exists(config_path):
        print(f"Error: config file not found: {config_path}")
        sys.exit(1)

    print(f"Loading config: {config_path}")
    cfg = load_config(config_path)

    differencer = DEMDifferencerParallel(
        path_dem1=cfg["dem1"]["path"],
        nickname_dem1=cfg["dem1"]["nickname"],
        src_vcrs_dem1=cfg["dem1"]["src_vcrs"],
        src_hcrs_dem1=cfg["dem1"]["src_hcrs"],
        nodata_dem1=cfg["dem1"]["nodata"],
        path_dem2=cfg["dem2"]["path"],
        nickname_dem2=cfg["dem2"]["nickname"],
        src_vcrs_dem2=cfg["dem2"]["src_vcrs"],
        src_hcrs_dem2=cfg["dem2"]["src_hcrs"],
        nodata_dem2=cfg["dem2"]["nodata"],
        roi=cfg["options"].get("roi", None),
        coregister=cfg["options"].get("coregister", False),
        num_sectors=args.num_sectors,
        num_workers=args.workers,
        path_dest=cfg["options"]["path_dest"],
    )

    differencer.run()