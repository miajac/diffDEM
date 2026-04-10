"""
dem_diff.py was created to difference digital elevation models (DEMs) that have 
differing grids, pixel sizes, horizontal coordinate systems, and vertical 
coordinate systems. It will save a raster file of the differenced dem to your 
folder of choice.  

To call this file:

    $ python dem_diff.py

Must customize run parameters at the bottom of this script. 

Requirements: os, sys, geoutils, numpy, xdem, pyproj
"""

import os
import sys
import xdem
import pyproj
import geoutils as gu
import numpy as np

from pyproj.transformer import TransformerGroup

class DEMDifferencer:
    """
    A class to difference two DEMs that may have differing grids, pixel sizes,
    horizontal coordinate systems, and vertical coordinate systems.

    All DEMs are reprojected to EPSG:32606 (UTM Zone 6N) horizontally and 
    EGM96 vertically before differencing.

    Parameters: 
    path_dem1 : str
        File path to the first DEM.
    path_dem2 : str
        File path to the second DEM (used as the reference grid).
    path_dest : str
        File path for the final product (the differenced DEM).
    nickname_dem1 : str
        Shortened name for DEM 1, used in plot titles and output filename.
    nickname_dem2 : str
        Shortened name for DEM 2, used in plot titles and output filename.
    src_vcrs_dem1 : str
        Vertical coordinate reference system (CRS) of DEM 1 (ex. "Ellipsoid", "EGM96").
    src_hcrs_dem1 : str
        Horizontal CRS of DEM 1 (as an EPSG code; ex. "EPSG:4326" for WGS84). For more information on EPSG strings see https://epsg.io/.
    nodata_dem1 : float
        Nodata value for DEM 1 (ex. 0, -9999, -3.4028235e+038).
    src_vcrs_dem2 : str
        Vertical coordinate reference system (CRS) of DEM 2 (ex. "Ellipsoid", "EGM96").
    src_hcrs_dem2 : str
        Horizontal CRS of DEM 2 (as an EPSG code; ex. "EPSG:4326" for WGS84). For more information on EPSG strings see https://epsg.io/.
    nodata_dem2 : float
        Nodata value for DEM 2 (ex. 0, -9999, -3.4028235e+038).
    roi : optional
        Optional region of interest used to clip DEMs to a common extent before 
        differencing. Should be a vector.
    coregister : bool, optional
        Select True or False to clarify whether Nuth & Kääb coregistration
        should be applied after alignment to correct for residual horizontal and 
        vertical offsets. Default is False. Only use when CRS metadata is known 
        to be incorrect and cannot be fixed at source. When True, a warning is 
        raised if shift_z > 2m suggesting a vertical datum issue that should 
        ideally be fixed in src_vcrs instead.
    """

    # Final target CRS for the DEMs
    TARGET_HCRS = "EPSG:32606" # WGS 84 / UTM zone 6N
    TARGET_VCRS = "EGM96" # same as "EPSG:5773", geoid

    # NAVD88 GEOID09 grid file for converting between NAD83 (EPSG:4269) and 
    # NAVD88 (EPSG:5703) coordinate systems. Auto-downloaded from cdn.proj.org 
    # on first use
    NAVD88_GRID = "us_noaa_geoid09_ak.tif"

    # Initiate variables user defined, with some variables having default 
    # settings if users do not specify
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
        roi = None,
        coregister = False,
    ):
        self.path_dem1     = path_dem1
        self.path_dem2     = path_dem2
        self.path_dest     = path_dest
        self.nickname_dem1 = nickname_dem1
        self.nickname_dem2 = nickname_dem2
        self.src_vcrs_dem1 = src_vcrs_dem1
        self.src_hcrs_dem1 = src_hcrs_dem1
        self.nodata_dem1   = nodata_dem1
        self.src_vcrs_dem2 = src_vcrs_dem2
        self.src_hcrs_dem2 = src_hcrs_dem2
        self.nodata_dem2   = nodata_dem2
        self.roi           = roi
        self.coregister    = coregister

        # File path for differenced DEM to be saved to
        OUTPUT_DIR  = self.path_dest

        # Output path: <dem2_nickname>_<dem1_nickname>.tif
        self.output_path = (
            f"{self.OUTPUT_DIR}{self.nickname_dem2}_{self.nickname_dem1}.tif"
        )

        # Placeholders for processed DEMs and output
        self.dem1     = None
        self.dem2     = None
        self.diff_dem = None
    
    def _check_grids(self):
        """
        Verify that required projection grid files are available locally or 
        downloadable. Raises an error early if a grid is missing.
        """
        
        # is the vcrs in dem 1 or dem 2 NAVD88 or its EPSG equivalent? if yes, 
            # continue
        needs_navd88 = any(
            vcrs in ("NAVD88", "EPSG:5703") 
            for vcrs in [self.src_vcrs_dem1, self.src_vcrs_dem2]
        )

        # determine whether the grid to convert from NAVD88 is downloaded, if 
            # not, download it
        if needs_navd88:
            data_dir = pyproj.datadir.get_data_dir()
            grid_path = os.path.join(data_dir, self.NAVD88_GRID)
            if os.path.exists(grid_path):
                print(f"  [grid check] {self.NAVD88_GRID} found locally at {grid_path}")
            else:
                print(f"  [grid check] {self.NAVD88_GRID} not found locally, attempting download from cdn.proj.org...")
                try:
                    tg = TransformerGroup("EPSG:4269", "EPSG:5703")
                    tg.download_grids(verbose=True)
                    if os.path.exists(grid_path):
                        print(f"  [grid check] download successful: {grid_path}")
                    else:
                        raise FileNotFoundError(
                            f"Grid {self.NAVD88_GRID} could not be downloaded to {data_dir}. "
                            f"Download it manually from https://cdn.proj.org/{self.NAVD88_GRID} "
                            f"and place it in {data_dir}"
                        )
                except Exception as e:
                    raise RuntimeError(
                        f"Failed to download required PROJ grid {self.NAVD88_GRID}. "
                        f"Check your internet connection or download manually from "
                        f"https://cdn.proj.org/{self.NAVD88_GRID} and place in {data_dir}"
                    ) from e


    def load(self):

        """Load both DEMs from file and print info before any corrections."""

        # add DEMs and manually assert their no data values
        self.dem1 = xdem.DEM(self.path_dem1, nodata=self.nodata_dem1)
        self.dem2 = xdem.DEM(self.path_dem2, nodata=self.nodata_dem2)
        
        print("Before corrections:")
        print(f"{self.nickname_dem1} source file: {self.path_dem1}")
        print(f"{self.nickname_dem1} assigned nodata: {self.nodata_dem1}")
        print(f"{self.nickname_dem1} source hCRS: {self.src_hcrs_dem1}")
        print(f"{self.nickname_dem1} source vCRS: {self.src_vcrs_dem1}")
        print(f"{self.nickname_dem1} -> will be converted to hCRS: {self.TARGET_HCRS}, vCRS: {self.TARGET_VCRS}")
        print()
        print(f"  {self.nickname_dem2} source file : {self.path_dem2}")
        print(f"  {self.nickname_dem2} assigned nodata : {self.nodata_dem2}")
        print(f"  {self.nickname_dem2} source hCRS : {self.src_hcrs_dem2}")
        print(f"  {self.nickname_dem2} source vCRS : {self.src_vcrs_dem2}")
        print(f"  {self.nickname_dem2} -> will be converted to hCRS: {self.TARGET_HCRS}, vCRS: {self.TARGET_VCRS}")
        print()
        print(f"{self.nickname_dem1} raster info:")
        self.dem1.info()
        print(f"{self.nickname_dem2} raster info:")
        self.dem2.info()

    def _prepare_dem(self, dem, src_hcrs, src_vcrs, nodata, nickname, src_path):

        """
        Assign source CRS and nodata, convert vertical datum, then reproject
        to the target horizontal CRS. Vertical conversion is done first so the
        geoid conversion uses the original native pixel positions.

        Parameters:

        dem      : xdem.DEM
        src_hcrs : str
        src_vcrs : str
        nodata   : float
        nickname : str
        src_path : str

        Returns: 

        xdem.DEM
            DEM reprojected to TARGET_HCRS with TARGET_VCRS vertical datum.
        """

        # Convert uint to float32 to support negative difference values
        if dem.data.dtype in [np.uint8, np.uint16, np.uint32]:
            print(f"  [{nickname}] converting dtype {dem.data.dtype} -> float32 to support negative values")
            dem = dem.astype(np.float32)

        # Assign source horizontal CRS
        dem.crs = src_hcrs
        print(f"  [{nickname}] hCRS assigned: {src_hcrs}")

        # Handle vertical CRS conversion
        # If converting from NAVD88 a two-step conversion is required since xdem
        # cannot go directly from NAVD88 to EGM96:
        #   Step 1: NAVD88 to Ellipsoid using NOAA GEOID09 Alaska grid
        #   Step 2: Ellipsoid to EGM96 (handled by xdem)

        if src_vcrs in ("EPSG:5703", "NAVD88"):
            print(f"[{nickname}] vCRS is NAVD88 — using two-step conversion:")
            print(f"[{nickname}]   Step 1: NAVD88 -> Ellipsoid (inverse of {self.NAVD88_GRID})")
            print(f"[{nickname}]   Note: {self.NAVD88_GRID} will be auto-downloaded from cdn.proj.org if not cached")
            dem.set_vcrs(self.NAVD88_GRID)
            dem.to_vcrs("Ellipsoid")
            print(f"[{nickname}]   Step 2: Ellipsoid -> {self.TARGET_VCRS}")
            dem.set_vcrs("Ellipsoid")
            dem.to_vcrs(self.TARGET_VCRS)
        else:
            dem.set_vcrs(src_vcrs)
            print(f"[{nickname}] vCRS assigned: {src_vcrs}")
            dem.to_vcrs(self.TARGET_VCRS)
            print(f"[{nickname}] vCRS converted: {src_vcrs} -> {self.TARGET_VCRS}")

        # Only set nodata if it differs from what is already assigned to avoid 
            # masking legitimate pixel values that match the nodata value
        if dem.nodata != nodata:
            dem.nodata = nodata
            print(f"[{nickname}] nodata set: {nodata}")
        else:
            print(f"[{nickname}] nodata already correct: {nodata}, skipping")

        # Reproject to target horizontal CRS
        if dem.crs.to_epsg() != int(self.TARGET_HCRS.split(":")[1]):
            dem = dem.reproject(crs=self.TARGET_HCRS, resampling="bilinear")
            print(f"[{nickname}] hCRS reprojected: {src_hcrs} -> {self.TARGET_HCRS}")
        else:
            print(f"[{nickname}] hCRS already {self.TARGET_HCRS}, skipping reprojection")

        # Re-stamp vertical CRS since reproject drops it
        dem.set_vcrs(self.TARGET_VCRS)
        print(f"[{nickname}] vCRS re-stamped after reproject: {self.TARGET_VCRS}")
        print(f"[{nickname}] source file: {dem.filename if dem.filename else "in-memory (reprojected from " + src_path + ")"}")

        return dem

    def prepare(self):

        """Apply CRS assignment, vertical conversion, and horizontal
        reprojection to both DEMs."""

        print("\nPreparing DEMs:")
        self.dem1 = self._prepare_dem(
            self.dem1, self.src_hcrs_dem1, self.src_vcrs_dem1, 
            self.nodata_dem1, self.nickname_dem1, self.path_dem1
        )
        self.dem2 = self._prepare_dem(
            self.dem2, self.src_hcrs_dem2, self.src_vcrs_dem2, 
            self.nodata_dem2, self.nickname_dem2, self.path_dem2
        )
        print(f"\nAfter vertical conversion and horizontal reprojection:")
        print(f"{self.nickname_dem1}:")
        self.dem1.info()
        print(f"{self.nickname_dem2}:")
        self.dem2.info()

    def align(self):

        """
        Align DEMs to a common grid. The DEM with the larger pixel size is used
        as the reference to avoid resampling to a finer resolution than the
        coarser input (which would create false precision). The finer DEM is
        resampled to match the coarser one's grid.
        """

        if self.roi is not None:
            self.dem1 = self.dem1.crop(self.roi)
            self.dem2 = self.dem2.crop(self.roi)
            print("\nAfter ROI crop:")
            print(f"{self.nickname_dem1}:")
            self.dem1.info()
            print(f"{self.nickname_dem2}:")
            self.dem2.info()

        # Determine which DEM has the coarser resolution
        # res is a tuple of (x_pixel_size, y_pixel_size), use x as representative
        res1 = self.dem1.res[0]
        res2 = self.dem2.res[0]

        print(f"\n  [{self.nickname_dem1}] pixel size: {res1}m")
        print(f"  [{self.nickname_dem2}] pixel size: {res2}m")

        if res1 >= res2:
            # dem1 is coarser or equal — reproject dem2 onto dem1's grid
            print(f"  {self.nickname_dem1} is coarser ({res1}m >= {res2}m) — reprojecting {self.nickname_dem2} to match")
            self.dem2 = self.dem2.reproject(self.dem1, resampling="bilinear")
            self.dem2.set_vcrs(self.TARGET_VCRS)
            self.dem1.set_vcrs(self.TARGET_VCRS)
            print(f"  Reference grid: {self.nickname_dem1}")
        else:
            # dem2 is coarser — reproject dem1 onto dem2's grid
            print(f"  {self.nickname_dem2} is coarser ({res2}m > {res1}m) — reprojecting {self.nickname_dem1} to match")
            self.dem1 = self.dem1.reproject(self.dem2, resampling="bilinear")
            self.dem1.set_vcrs(self.TARGET_VCRS)
            self.dem2.set_vcrs(self.TARGET_VCRS)
            print(f"  Reference grid: {self.nickname_dem2}")

        print("\nAfter full alignment:")
        print(f"{self.nickname_dem1}:")
        self.dem1.info()
        print(f"{self.nickname_dem2}:")
        self.dem2.info()
    
    def _coregister(self):

        """
        Optional: Apply Nuth & Kääb coregistration to correct for residual 
        horizontal and vertical offsets between the two DEMs after alignment.
        Only use when CRS metadata is known to be incorrect and cannot be fixed
        at source.
        """

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
            print(f"\n*** WARNING: shift_z of {shift_z:.2f}m suggests a residual")
            print(f"*** vertical datum offset.")
            print(f"*** Coregistration is correcting for this automatically,")
            print(f"*** but the root cause should ideally be fixed by verifying")
            print(f"*** the src_vcrs parameter for each DEM is correct.")
            print(f"*** Common cause in Alaska: IFSAR data labeled as EGM96")
            print(f"*** but actually stored on the WGS84 ellipsoid (~14-16m offset).")

        # Apply correction and reproject back onto dem2's grid
        self.dem1 = nuth_kaab.apply(self.dem1)
        self.dem1 = self.dem1.reproject(self.dem2, resampling="bilinear")
        self.dem1.set_vcrs(self.TARGET_VCRS)

        print(f"\nAfter coregistration:")
        print(f"{self.nickname_dem1}:")
        self.dem1.info()

    def difference(self):

        """Difference the two aligned DEMs (DEM2 - DEM1) and print stats."""
        
        self.diff_dem = self.dem2 - self.dem1

        # Stamp vertical CRS on diff DEM since arithmetic drops it
        self.diff_dem.set_vcrs(self.TARGET_VCRS)
        print("\nDifference DEM stats:")
        self.diff_dem.info(stats=True)

    def check_stable_terrain(self):

        """
        Check elevation differences on flat stable terrain after differencing.
        A systematic offset on flat terrain indicates a vertical datum mismatch.
        Mean offset > 2m on flat terrain (slope < 5 degrees) triggers a warning.
        """

        print("\nStable terrain check:")
        slope = xdem.terrain.slope(self.diff_dem)
        flat_mask = (slope.data < 5) & (~np.ma.getmaskarray(self.diff_dem.data))
        diff_flat = self.diff_dem.data.data[flat_mask]

        mean_offset   = np.mean(diff_flat)
        median_offset = np.median(diff_flat)
        std_offset    = np.std(diff_flat)

        print(f"Flat terrain pixel count : {len(diff_flat)}")
        print(f"Mean offset              : {mean_offset:.2f}m")
        print(f"Median offset            : {median_offset:.2f}m")
        print(f"Std dev                  : {std_offset:.2f}m")

        if abs(median_offset) > 2:
            print(f"\n*** WARNING: median offset of {median_offset:.2f}m on flat stable terrain")
            print(f"*** suggests a residual vertical datum mismatch.")
            if 10 < abs(median_offset) < 20:
                print(f"*** Offset of ~14-15m in Alaska typically indicates the DEM is")
                print(f"*** on the WGS84 ellipsoid rather than EGM96.")
                print(f"*** Try changing src_vcrs from 'EGM96' to 'Ellipsoid' and rerun.")
        else:
            print(f" Offset within acceptable range — vertical datums appear consistent")

    def plot(self):
    
        """Plot both input DEMs and the differenced DEM."""
        
        self.dem1.plot(
            cmap="RdYlBu", vmin=0, vmax=2000,
            cbar_title=f"{self.nickname_dem1} (m)"
        )
        self.dem2.plot(
            cmap="RdYlBu", vmin=0, vmax=2000,
            cbar_title=f"{self.nickname_dem2} (m)"
        )
        self.diff_dem.plot(
            cmap="RdYlBu", vmin=-20, vmax=20,
            cbar_title=f"{self.nickname_dem2} - {self.nickname_dem1} (m)"
        )

    def save(self):

        """Save the differenced DEM to the output path."""

        self.diff_dem.to_file(self.output_path)
        print(f"\nDifferenced DEM saved to: {self.output_path}")

    def run(self):

        """Run the full pipeline: load, prepare, align, difference, plot, save."""
        
        self._check_grids()
        self.load()
        self.prepare()
        self.align()
        if self.coregister:
            self._coregister()
        self.difference()
        self.check_stable_terrain()
        # self.plot()
        self.save()

# Run
if __name__ == "__main__":

    differencer = DEMDifferencer(
        path_dem1     = "~/REPOS/diffDEMs/exampleData/IFSAR-Horz-AlbersConicalEqualArea/IFSAR_DTM_Summer_2010/IFSAR_DTM_Summer_2010.tif",
        nickname_dem1 = "IFSAR_DTM_2010",
        src_vcrs_dem1 = "NAVD88",    # stored as unknown in metadata, known to be NAVD88
                                        # two-step conversion: NAVD88 -> Ellipsoid -> EGM96, via us_noaa_geoid09_ak.tif (GEOID09 Alaska)
        src_hcrs_dem1 = "EPSG:3338", # NAD83 / Alaska Albers 
        nodata_dem1   = -9999,
        path_dem2     = "~/REPOS/diffDEMs/exampleData/Canwell_4Aug25_DTM.tif", # most recent DEM
        nickname_dem2 = "Lidar2025",
        src_vcrs_dem2 = "Ellipsoid",
        src_hcrs_dem2 = "EPSG:32606",
        coregister = True,
        nodata_dem2   = -9999,
        roi           = None,
        path_dest = "~/REPOS/diffDEMs/differencedDEMs" # file path for differenced DEMs
    )

    differencer.run()

## Example parameters for each dem:
## user must still specify whether each DEM is 1 or 2

## IFSAR DSM:
# path_dem     = "~/REPOS/diffDEMs/exampleData/IFSAR-Horz-AlbersConicalEqualArea/IFSAR_DSM_Summer_2010/IFSAR_DSM_Summer_2010.tif",
# nickname_dem = "IFSAR_DSM_2010",
# src_vcrs_dem = "NAVD88",    # stored as unknown in metadata, known to be NAVD88
#                                 # two-step conversion: NAVD88 -> Ellipsoid -> EGM96, via us_noaa_geoid09_ak.tif (GEOID09 Alaska)
# src_hcrs_dem = "EPSG:3338", # NAD83 / Alaska Albers 
# nodata_dem   = -9999,
# coregister   = True,

## IFSAR DTM:
# path_dem     = "~/REPOS/diffDEMs/exampleData/IFSAR-Horz-AlbersConicalEqualArea/IFSAR_DTM_Summer_2010/IFSAR_DTM_Summer_2010.tif",
# nickname_dem = "IFSAR_DTM_2010",
# src_vcrs_dem = "NAVD88",    # stored as unknown in metadata, known to be NAVD88
#                                 # two-step conversion: NAVD88 -> Ellipsoid -> EGM96, via us_noaa_geoid09_ak.tif (GEOID09 Alaska)
# src_hcrs_dem = "EPSG:3338", # NAD83 / Alaska Albers 
# nodata_dem   = -9999,
# coregister   = True,

## Drone:
# path_dem     = "~/REPOS/diffDEMs/exampleData/Drone_2025Aug13_DEM/Drone_2025Aug13_DEM.tif",
# nickname_dem = "CanwellDrone08132025",
# src_vcrs_dem = "Ellipsoid",
# src_hcrs_dem = "EPSG:32606",
# nodata_dem   = -9999,
# coregister   = False,

## Isabel pass:
# path_dem     = "~/REPOS/diffDEMs/exampleData/Isabel-Horz-Alaska-AEAC/IsabelPass_DSM_2000.tif",
# nickname_dem = "IsabelIFSAR2000",
# src_vcrs_dem = "EGM96",
# src_hcrs_dem = "EPSG:32606",
# nodata_dem   = 0,
# coregister   = True,

## Arctic DEM:
# path_dem     = "~/REPOS/diffDEMs/exampleData/ArcticDEM_20090606_EPSG_4326_most.tif",
# nickname_dem = "Arctic2009_06_06",
# src_vcrs_dem = "Ellipsoid",
# src_hcrs_dem = "EPSG:4326",
# nodata_dem   = -9999,
# coregister   = False,

## 2025 Lidar Survey:
# path_dem     = "~/REPOS/diffDEMs/exampleData/Canwell_4Aug25_DTM.tif",
# nickname_dem = "Lidar2025",
# src_vcrs_dem = "Ellipsoid",
# src_hcrs_dem = "EPSG:32606",
# nodata_dem   = -9999,                       