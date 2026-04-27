"""
dem_diff_parallel.py
Parallelized DEM differencing with sector-based processing.

This version parallelizes:
1. DEM preparation (loading, CRS assignment, vertical conversion, reprojection)
2. Differencing via sector-based processing (split, difference in parallel, 
   mosaic)

Usage:
    $ python dem_diff_parallel.py config.yml --num-sectors 10 --workers auto

Where config.yml is a parameter file (same format as dem_diff.py).
"""

import os
import sys
import yaml
import xdem
import pyproj
import geoutils as gu
import rasterio
import numpy as np
import tempfile
import shutil
import argparse

        

from pyproj.transformer import TransformerGroup
from multiprocessing import Pool, cpu_count
from functools import partial
from pathlib import Path
from rasterio.merge import merge


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
                    f"Config section '{section}' missing required field: '{field}'"
                )
    # Expand ~ to full home directory path
    cfg["dem1"]["path"] = os.path.expanduser(cfg["dem1"]["path"])
    cfg["dem2"]["path"] = os.path.expanduser(cfg["dem2"]["path"])
    cfg["options"]["path_dest"] = os.path.expanduser(
        cfg["options"]["path_dest"]
    )

    return cfg


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
                    f"[grid check] {self.NAVD88_GRID} found locally at "
                    f"{grid_path}"
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

        print(f"\nDefined {len(sectors)} sectors" 
              f"({self.num_sectors}x{self.num_sectors})")
        print(f"Sector dimensions: ~{sector_height}x{sector_width} pixels")
        return sectors


    def _process_sector(self, sector_info):
        sector_id, y_start, y_end, x_start, x_end = sector_info
        try:
            transform = self.dem1.transform
            left = transform.c + x_start * transform.a
            top = transform.f + y_start * transform.e
            right = transform.c + x_end * transform.a
            bottom = transform.f + y_end * transform.e

            dem1_sector = self.dem1.crop((left, bottom, right, top), 
                                         inplace=False)
            dem2_sector = self.dem2.crop((left, bottom, right, top), 
                                         inplace=False)

            diff_sector = dem2_sector - dem1_sector

            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                diff_sector.set_vcrs(self.TARGET_VCRS)
                diff_sector.nodata = self.dem2.nodata

            temp_path = os.path.join(self.temp_dir, 
                                     f"sector_{sector_id:05d}.tif")
            diff_sector.to_file(temp_path)
            return temp_path

        except Exception as e:
            print(f"[Sector {sector_id:05d}] ERROR: {e}")
            raise
    
    def _merge_rasters(self, raster_list):
        """
        Merge a list of non-overlapping geoutils Rasters into one.
        Uses rasterio.merge on the temp sector files directly.
        """
        # Open all sector files
        src_files = [rasterio.open(r.filename) for r in raster_list]

        try:
            merged_data, merged_transform = merge(src_files)
            profile = src_files[0].profile.copy()
            profile.update({
                "height": merged_data.shape[1],
                "width": merged_data.shape[2],
                "transform": merged_transform,
                "tiled": False,
            })
            profile.pop("blockxsize", None)
            profile.pop("blockysize", None)

            # Write merged result to output directory so it survives cleanup
            os.makedirs(self.path_dest, exist_ok=True)
            merged_path = os.path.join(self.path_dest, "merged_temp.tif")
            with rasterio.open(merged_path, "w", **profile) as dst:
                dst.write(merged_data)

        finally:
            for src in src_files:
                src.close()

        # Load back as xdem.DEM to preserve type
        merged_dem = xdem.DEM(merged_path, nodata=self.dem2.nodata)
        merged_dem.load() # force data into memory before file deletion
        merged_dem.set_vcrs(self.TARGET_VCRS)
        os.remove(merged_path)              
        return merged_dem


    def _mosaic_sectors(self, sector_files):
        """
        Mosaic sector files into a single DEM.
        Uses geoutils to read and combine sectors while preserving
        georeferencing.
        """
        # Load all sector files
        sectors_data = [gu.Raster(f) for f in sector_files]

        # Stack and merge
        # Since sectors are non-overlapping and adjacent, we reference 
        # _merge_rasters (above)
        mosaicked = self._merge_rasters(sectors_data)
        return mosaicked
    

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
        print(f"Processing {len(sectors)} sectors with "
              f"{self.num_workers} workers...")
        with Pool(processes=self.num_workers) as pool:
            sector_files = pool.map(self._process_sector, sectors)

        print(f"\nAll {len(sector_files)} sectors processed successfully")

        # Mosaic sector files back into single DEM using _mosaic_sectors (above)
        print("\nMosaicking sectors into final DEM...")
        self.diff_dem = self._mosaic_sectors(sector_files)
        print("Mosaicking complete")

        # Clean up temp directory
        shutil.rmtree(self.temp_dir)
        print(f"Cleaned up temp directory: {self.temp_dir}")


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


if __name__ == "__main__":

    # Command-line argument parsing
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