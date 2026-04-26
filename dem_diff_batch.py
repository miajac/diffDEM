"""
dem_diff_batch.py 
Version of dem_diff for batch processing of multiple DEM pairs with 
supercomputer access.

Each MPI task processes a different DEM pair independently, making this well
suited for supercomputer environments where many DEM pairs need to be 
differenced.

Usage:
    $ sbatch dem_diff_batch.sh

    Where dem_diff_batch.sh is the SLURM batch script that calls this script
    with mpirun. The batch script reads from a config YAML file that lists all
    DEMs and specifies how to pair them. See config_batch_template.yml for a
    full template.


Where config_batch.yml lists all DEMs and specifies how to pair them.
See config_batch_template.yml for a full template.

Requirements: demDiff environment and a supercomputer with MPI support.

Pair modes:
    sequential: differences consecutive DEMs (e.g. DEM1-DEM2, DEM2-DEM3)
    all: differences every possible combination
    explicit: user-defined list of pairs by nickname, with optional per-pair 
              coregister override
"""

import os
import sys
import yaml
import xdem
import pyproj
import geoutils as gu
import numpy as np
from itertools import combinations
from pyproj.transformer import TransformerGroup
from mpi4py import MPI


class DEMDifferencerBatch:
    """
    Batch DEM differencer that distributes pairs across tasks.

    Each task processes its assigned DEM pair using the full dem_diff
    pipeline (load, prepare, align, coregister, difference, save).
    
    Coregistration can be set globally via coregister_default, and overridden
    per pair in the config file.

    Parameters:
    dems : list of dict
        List of DEM parameter dictionaries, each containing:
        path, nickname, src_vcrs, src_hcrs, nodata.
    pairs : list of tuple
        List of (dem1_dict, dem2_dict, coregister) tuples to process.
        coregister is a bool that overrides the global default for that pair.
    path_dest : str
        Output directory for differenced DEMs.
    coregister_default : bool
        Global default for coregistration, used when not specified per pair.
    roi : str or None
        Path to region of interest vector file, or None.
    """

    TARGET_HCRS = "EPSG:32606"
    TARGET_VCRS = "EGM96"
    NAVD88_GRID = "us_noaa_geoid09_ak.tif"


    def __init__(
        self, dems, pairs, path_dest, coregister_default=False, roi=None
    ):
        self.dems = dems
        self.pairs = pairs
        self.path_dest = path_dest
        self.coregister_default = coregister_default
        self.roi = roi

        # MPI setup
        self.comm = MPI.COMM_WORLD
        self.task = self.comm.Get_rank()
        self.num_tasks = self.comm.Get_size()


    def _check_grids(self, vcrs_list):
        """Verify NAVD88 grid is available if needed."""
        needs_navd88 = any(
            vcrs in ("NAVD88", "EPSG:5703") for vcrs in vcrs_list
        )
        if needs_navd88:
            data_dir = pyproj.datadir.get_data_dir()
            grid_path = os.path.join(data_dir, self.NAVD88_GRID)
            if not os.path.exists(grid_path):
                print(
                    f"[Task {self.task}] {self.NAVD88_GRID} not found, "
                    f"attempting download..."
                )
                try:
                    tg = TransformerGroup("EPSG:4269", "EPSG:5703")
                    tg.download_grids(verbose=False)
                except Exception as e:
                    raise RuntimeError(
                        f"Failed to download {self.NAVD88_GRID}. "
                        f"Download manually from https://cdn.proj.org/"
                        f"{self.NAVD88_GRID} and place in {data_dir}."
                    ) from e
            else:
                print(
                    f"[Task {self.task}] {self.NAVD88_GRID} found at "
                    f"{grid_path}"
                )


    def _prepare_dem(self, dem, src_hcrs, src_vcrs, nodata, nickname):
        """Prepare a single DEM: dtype, CRS, vertical datum, reproject."""
        if dem.data.dtype in [np.uint8, np.uint16, np.uint32]:
            dem = dem.astype(np.float32)
            print(f"[Task {self.task}][{nickname}] converted to float32")

        dem.crs = src_hcrs

        if src_vcrs in ("EPSG:5703", "NAVD88"):
            dem.set_vcrs(self.NAVD88_GRID)
            dem.to_vcrs("Ellipsoid")
            dem.set_vcrs("Ellipsoid")
            dem.to_vcrs(self.TARGET_VCRS)
            print(
                f"[Task {self.task}][{nickname}] vCRS: NAVD88 -> "
                f"Ellipsoid -> {self.TARGET_VCRS}"
            )
        else:
            dem.set_vcrs(src_vcrs)
            dem.to_vcrs(self.TARGET_VCRS)
            print(
                f"[Task {self.task}][{nickname}] vCRS: {src_vcrs} -> "
                f"{self.TARGET_VCRS}"
            )

        if dem.nodata != nodata:
            dem.nodata = nodata

        if dem.crs.to_epsg() != int(self.TARGET_HCRS.split(":")[1]):
            dem = dem.reproject(crs=self.TARGET_HCRS, resampling="bilinear")
            print(
                f"[Task {self.task}][{nickname}] hCRS: {src_hcrs} -> "
                f"{self.TARGET_HCRS}"
            )

        dem.set_vcrs(self.TARGET_VCRS)
        return dem


    def _process_pair(self, dem1_cfg, dem2_cfg, coregister):
        """
        Run the full differencing pipeline for a single DEM pair.
        Called independently by each MPI task for its assigned pair.

        Parameters:
        dem1_cfg : dict
            DEM 1 config (older DEM).
        dem2_cfg : dict
            DEM 2 config (more recent DEM).
        coregister : bool
            Whether to apply Nuth & Kääb coregistration for this pair.
        """
        n1 = dem1_cfg["nickname"]
        n2 = dem2_cfg["nickname"]
        print(
            f"\n[Task {self.task}] Starting pair: {n2} - {n1} "
            f"(coregister={coregister})"
        )

        # Check grids
        self._check_grids([dem1_cfg["src_vcrs"], dem2_cfg["src_vcrs"]])

        # Load
        dem1 = xdem.DEM(dem1_cfg["path"], nodata=dem1_cfg["nodata"])
        dem2 = xdem.DEM(dem2_cfg["path"], nodata=dem2_cfg["nodata"])
        print(f"[Task {self.task}] Loaded {n1} and {n2}")

        # Prepare
        dem1 = self._prepare_dem(
            dem1, dem1_cfg["src_hcrs"], dem1_cfg["src_vcrs"],
            dem1_cfg["nodata"], n1
        )
        dem2 = self._prepare_dem(
            dem2, dem2_cfg["src_hcrs"], dem2_cfg["src_vcrs"],
            dem2_cfg["nodata"], n2
        )

        # Align to coarser grid
        res1 = dem1.res[0]
        res2 = dem2.res[0]
        if res1 >= res2:
            dem2 = dem2.reproject(dem1, resampling="bilinear")
            dem2.set_vcrs(self.TARGET_VCRS)
            dem1.set_vcrs(self.TARGET_VCRS)
            print(f"[Task {self.task}] Aligned to {n1} grid ({res1:.1f}m)")
        else:
            dem1 = dem1.reproject(dem2, resampling="bilinear")
            dem1.set_vcrs(self.TARGET_VCRS)
            dem2.set_vcrs(self.TARGET_VCRS)
            print(f"[Task {self.task}] Aligned to {n2} grid ({res2:.1f}m)")

        # Crop to ROI if provided
        if self.roi is not None:
            dem1 = dem1.crop(self.roi)
            dem2 = dem2.crop(self.roi)

        # Coregister — uses per-pair setting
        if coregister:
            print(f"[Task {self.task}] Coregistering (Nuth & Kääb)...")
            nuth_kaab = xdem.coreg.NuthKaab()
            nuth_kaab.fit(dem2, dem1)
            shift_z = nuth_kaab.meta["outputs"]["affine"]["shift_z"]
            print(f"[Task {self.task}] shift_z: {shift_z:.3f}m")
            if abs(shift_z) > 2:
                print(
                    f"[Task {self.task}] *** WARNING: shift_z of "
                    f"{shift_z:.2f}m suggests a residual vertical datum "
                    f"offset."
                )
            dem1 = nuth_kaab.apply(dem1)
            dem1 = dem1.reproject(dem2, resampling="bilinear")
            dem1.set_vcrs(self.TARGET_VCRS)

        # Difference
        diff_dem = dem2 - dem1
        diff_dem.set_vcrs(self.TARGET_VCRS)
        print(f"[Task {self.task}] Differenced {n2} - {n1}")

        # Stable terrain check
        slope = xdem.terrain.slope(diff_dem)
        flat_mask = (slope.data < 5) & (~np.ma.getmaskarray(diff_dem.data))
        diff_flat = diff_dem.data.data[flat_mask]
        median_offset = np.median(diff_flat)
        print(
            f"[Task {self.task}] Stable terrain median offset: "
            f"{median_offset:.2f}m"
        )
        if abs(median_offset) > 2:
            print(
                f"[Task {self.task}] *** WARNING: median offset of "
                f"{median_offset:.2f}m on flat terrain suggests vertical "
                f"datum mismatch."
            )

        # Save
        os.makedirs(self.path_dest, exist_ok=True)
        output_path = os.path.join(self.path_dest, f"{n2}_{n1}.tif")
        diff_dem.to_file(output_path)
        print(f"[Task {self.task}] Saved: {output_path}")
        return output_path


    def run(self):
        """
        Distribute pairs across MPI tasks and process each independently.
        Task 0 reports a summary when all tasks are done.
        """
        total_pairs = len(self.pairs)

        if self.task == 0:
            print(
                f"\nBatch DEM differencing: {total_pairs} pairs, "
                f"{self.num_tasks} MPI tasks"
            )
            for i, (d1, d2, coreg) in enumerate(self.pairs):
                print(
                    f"  Pair {i:03d}: {d2['nickname']} - {d1['nickname']} "
                    f"(coregister={coreg})"
                )

        # Distribute pairs across tasks — round-robin
        my_pairs = [
            (i, self.pairs[i])
            for i in range(total_pairs)
            if i % self.num_tasks == self.task
        ]

        print(
            f"\n[Task {self.task}] assigned {len(my_pairs)} pair(s): "
            f"{[self.pairs[i][1]['nickname'] + '-' + 
                self.pairs[i][0]['nickname'] for i, _ in my_pairs]}"
        )

        # Process assigned pairs
        results = []
        for pair_idx, (dem1_cfg, dem2_cfg, coregister) in my_pairs:
            try:
                output = self._process_pair(dem1_cfg, dem2_cfg, coregister)
                results.append((pair_idx, "SUCCESS", output))
            except Exception as e:
                print(f"[Task {self.task}] ERROR on pair {pair_idx}: {e}")
                results.append((pair_idx, "FAILED", str(e)))

        # Gather results on task 0 and print summary
        if MPI_AVAILABLE:
            all_results = self.comm.gather(results, root=0)
        else:
            all_results = [results]

        if self.task == 0:
            print("\n" + "="*60)
            print("BATCH COMPLETE — Summary:")
            flat_results = [
                r for task_results in all_results for r in task_results
            ]
            flat_results.sort(key=lambda x: x[0])
            for pair_idx, status, info in flat_results:
                d1 = self.pairs[pair_idx][0]["nickname"]
                d2 = self.pairs[pair_idx][1]["nickname"]
                print(
                    f"  Pair {pair_idx:03d} ({d2}-{d1}): {status} — {info}"
                )
            print("="*60)


def build_pairs(dems, pair_mode, coregister_default, explicit_pairs=None):
    """
    Build list of (dem1, dem2, coregister) tuples from DEM list.

    Parameters:
    dems : list of dict
        List of DEM configs.
    pair_mode : str
        "sequential", "all", or "explicit".
    coregister_default : bool
        Global default coregister value, used for sequential and all modes,
        and as fallback for explicit pairs that don't specify it.
    explicit_pairs : list of list, optional
        For "explicit" mode: list of [nickname1, nickname2] or
        [nickname1, nickname2, coregister] entries.

    Returns:
    list of (dem1_dict, dem2_dict, coregister)
    """
    dem_by_nickname = {d["nickname"]: d for d in dems}

    if pair_mode == "sequential":
        return [
            (dems[i], dems[i + 1], coregister_default)
            for i in range(len(dems) - 1)
        ]

    elif pair_mode == "all":
        return [
            (d1, d2, coregister_default)
            for d1, d2 in combinations(dems, 2)
        ]

    elif pair_mode == "explicit":
        if not explicit_pairs:
            raise ValueError(
                "pair_mode is 'explicit' but no explicit_pairs list provided "
                "in config."
            )
        pairs = []
        for p in explicit_pairs:
            n1, n2 = p[0], p[1]
            # Use per-pair coregister if provided, otherwise use global default
            coreg = p[2] if len(p) > 2 else coregister_default
            if n1 not in dem_by_nickname:
                raise ValueError(f"Nickname '{n1}' not found in dems list.")
            if n2 not in dem_by_nickname:
                raise ValueError(f"Nickname '{n2}' not found in dems list.")
            pairs.append((dem_by_nickname[n1], dem_by_nickname[n2], coreg))
        return pairs

    else:
        raise ValueError(
            f"Unknown pair_mode '{pair_mode}'. "
            f"Choose from: 'sequential', 'all', 'explicit'."
        )


def load_config(config_path):
    """Load and validate a batch YAML config file."""
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    if "dems" not in cfg or not cfg["dems"]:
        raise ValueError("Config missing required 'dems' list.")
    if "options" not in cfg:
        raise ValueError("Config missing required 'options' section.")
    if "path_dest" not in cfg["options"]:
        raise ValueError(
            "Config 'options' missing required field: 'path_dest'."
        )

    required_dem_fields = ["path", "nickname", "src_vcrs", "src_hcrs", "nodata"]
    for i, dem in enumerate(cfg["dems"]):
        for field in required_dem_fields:
            if field not in dem:
                raise ValueError(
                    f"DEM entry {i} missing required field: '{field}'"
                )
        dem["path"] = os.path.expanduser(dem["path"])

    cfg["options"]["path_dest"] = os.path.expanduser(
        cfg["options"]["path_dest"]
    )

    return cfg


if __name__ == "__main__":
    import argparse

    # Parse command-line argument for config file path
    parser = argparse.ArgumentParser(
        description="MPI batch DEM differencing across multiple pairs"
    )
    parser.add_argument("config", help="Path to batch config YAML file")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        print(f"Error: config file not found: {args.config}")
        sys.exit(1)

    cfg = load_config(args.config)

    pair_mode = cfg["options"].get("pairs", "sequential")
    coregister_default = cfg["options"].get("coregister", False)
    explicit_pairs = cfg["options"].get("explicit_pairs", None)

    pairs = build_pairs(
        cfg["dems"], pair_mode, coregister_default, explicit_pairs
    )

    differencer = DEMDifferencerBatch(
        dems=cfg["dems"],
        pairs=pairs,
        path_dest=cfg["options"]["path_dest"],
        coregister_default=coregister_default,
        roi=cfg["options"].get("roi", None),
    )

    differencer.run()