#!/usr/bin/env python3
"""
Offline test for the AWI-ESM-v3.4.1-CAV-ICB iceberg plugin.

This test executes the iceberg generation plugin using provided test data
without requiring a full esm_tools simulation environment.

Hardcoded paths for DARS2cav mesh and test data.

Usage:
    source ~/.bashrc
    conda activate esm_tools
    python test_iceberg_plugin.py
"""

import os
import sys
import shutil
import tempfile

# Add parent directory to path to import the plugin
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from icb_apply_distribution_functions import IcebergCalving


# Hardcoded paths
MESH_DIR = "/work/ab0995/a270186/model_inputs/fesom2/mesh/DARS2cav/"
FW_FILE = "/home/a/a270186/esm_tools/plugins/AWI-ESM-v3.4.x-CAV-ICB/test/data/fw.fesom.1600.nc"
CALVING_FILE = "/home/a/a270186/esm_tools/plugins/AWI-ESM-v3.4.x-CAV-ICB/test/data/calving_AA.fesom.1600.nc"
RESTART_FILE = "/home/a/a270186/esm_tools/plugins/AWI-ESM-v3.4.x-CAV-ICB/test/output/iceberg.restart.ISM"
BASIN_FILE = "/home/a/a270186/esm_tools/plugins/basins/basins_antarctica.nc"
SEED_YEAR = 1600
OUTPUT_DIR = "/home/a/a270186/esm_tools/plugins/AWI-ESM-v3.4.x-CAV-ICB/test/output/"


def verify_paths():
    """Verify that all hardcoded paths exist."""
    errors = []
    
    if not os.path.exists(MESH_DIR):
        errors.append(f"Mesh directory not found: {MESH_DIR}")
    elif not os.path.exists(os.path.join(MESH_DIR, "nod2d.out")):
        errors.append(f"nod2d.out not found in mesh directory: {MESH_DIR}")
    
    if not os.path.exists(FW_FILE):
        errors.append(f"Freshwater file not found: {FW_FILE}")
    
    if not os.path.exists(CALVING_FILE):
        errors.append(f"Calving file not found: {CALVING_FILE}")
    
    if not os.path.exists(BASIN_FILE):
        errors.append(f"Basin file not found: {BASIN_FILE}")
    
    if errors:
        print("ERROR: Missing required files:")
        for err in errors:
            print(f"  - {err}")
        sys.exit(1)


def run_iceberg_test(output_dir):
    """
    Run the iceberg generation test.
    """
    print("\n" + "="*60)
    print("ICEBERG PLUGIN OFFLINE TEST")
    print("="*60)
    
    print(f"\n* Mesh directory: {MESH_DIR}")
    print(f"* Freshwater file: {FW_FILE}")
    print(f"* Calving file: {CALVING_FILE}")
    print(f"* Basin file: {BASIN_FILE}")
    print(f"* Output directory: {output_dir}")
    print(f"* Seed year: {SEED_YEAR}")
    
    # Initialize the IcebergCalving class
    print("\n* Initializing IcebergCalving...")
    ib = IcebergCalving(
        mesh_path=MESH_DIR,
        icb_path=output_dir,
        latest_restart_file=RESTART_FILE,
        abg=[50, 15, -90],
        scaling_factor=[100, 50, 25, 1, 1, 1],
        seed=1000,
        bcavities=True,
        ibareamax=100, #km2
        domain="SH",
        basin_file=BASIN_FILE,
        fw_file=FW_FILE,
        calving_file=CALVING_FILE
    )
    
    # Create dataframe
    print("\n* Creating dataframe...")
    ib.create_dataframe()
    
    # Generate icebergs
    print("\n* Generating icebergs...")
    ib._icb_generator(fmode="w")
    
    # Check output files
    print("\n* Checking output files...")
    expected_files = [
        "icb_longitude.dat",
        "icb_latitude.dat",
        "icb_length.dat",
        "icb_height.dat",
        "icb_scaling.dat",
        "icb_felem.dat"
    ]
    
    generated_files = []
    total_icebergs = 0
    for fname in expected_files:
        fpath = os.path.join(output_dir, fname)
        if os.path.exists(fpath):
            size = os.path.getsize(fpath)
            n_lines = sum(1 for _ in open(fpath))
            print(f"  ✓ {fname}: {size} bytes, {n_lines} entries")
            generated_files.append(fname)
            if fname == "icb_latitude.dat":
                total_icebergs = n_lines
        else:
            print(f"  ✗ {fname}: NOT FOUND")
    
    print(f"\n* Successfully generated {len(generated_files)}/{len(expected_files)} iceberg files")
    print(f"* Total icebergs generated: {total_icebergs}")
    
    return len(generated_files) == len(expected_files)


def main():
    """Main test function."""
    print("AWI-ESM-v3.4.1-CAV-ICB Iceberg Plugin Offline Test")
    print("="*60)
    
    # Verify all hardcoded paths exist
    print("* Verifying input paths...")
    verify_paths()
    print("* All input paths verified.")
    
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"* Output directory: {OUTPUT_DIR}")
    
    # Run the test
    success = False
    try:
        success = run_iceberg_test(OUTPUT_DIR)
    except Exception as e:
        print(f"\nERROR during test execution: {e}")
        import traceback
        traceback.print_exc()
    
    # Exit with appropriate code
    if success:
        print(f"\n* Output files written to: {OUTPUT_DIR}")
        print("\n" + "="*60)
        print("TEST PASSED: Iceberg files generated successfully!")
        print("="*60)
        sys.exit(0)
    else:
        print("\n" + "="*60)
        print("TEST FAILED: Could not generate all iceberg files")
        print("="*60)
        sys.exit(1)


if __name__ == "__main__":
    main()
