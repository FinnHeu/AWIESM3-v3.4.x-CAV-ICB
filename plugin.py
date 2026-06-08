import f90nml
import os
import shutil
import glob

from icb_apply_distribution_functions import IcebergCalving
from esm_runscripts.namelists import Namelist

def prep_icebergs(config):
    """
    Prepares iceberg files for FESOM simulation.
    
    Inputs:
    - config: esm_tools configuration dictionary
    """

    print("\n=== PREPARE ICEBERGS ===\n")

    # get config values for iceberg usage
    with_icb   = config["fesom"].get("with_icb", False)
    run_number = config["general"].get("run_number", 0)
    # if iceberg usage is enabled...
    if with_icb:
        # if first year of simulation...
        if run_number == 1:
            print("* first year of simulation")
            icb_ini_dir = config["fesom"]["ini_iceberg_dir"]
            # copy iceberg files from initial directory to experiment work dir.
            icb_files = glob.glob(os.path.join(icb_ini_dir, "icb*.dat"))
            for f in icb_files:
                shutil.copy(f, config["general"]["thisrun_work_dir"])
            print(f"* using iceberg initial files from: {icb_ini_dir}")
            # also copy the iceberg.restart.ISM file from icb_ini_dir if it exists for a hot restart
            icb_restart_file = os.path.join(icb_ini_dir, "iceberg.restart.ISM")
            if os.path.exists(icb_restart_file):
                shutil.copy(icb_restart_file, config["general"]["thisrun_work_dir"])
                print(f"* using iceberg restart file from: {icb_ini_dir}")
            else:
                print(f"* no iceberg restart file found in: {icb_ini_dir}")

        else:
            print("* not first year of simulation")    
            print(f"* updating icebergs for run number: {run_number}")
            # run update icebergs, which creates the icb*.dat files
            config = update_icebergs(config)

        return config

def update_icebergs(config):
    """
    Updates iceberg files for FESOM simulation.
    Computes the residual of cavity melt and available Antarctic solid runoff and converts it to iceberg calving.
    
    Inputs:
    - config: esm_tools configuration dictionary
    """
    print(" * start updating icebergs")
    
    import os
    
    mesh_dir = config["fesom"]["mesh_dir"]
    basin_file = config["fesom"].get("basin_file", "")
    icb_restart_file = config["fesom"]["restart_in_sources"].get("icb_restart_ISM", "")
    scaling_factor = config["fesom"].get("scaling_factor", [1, 1, 1, 1, 1, 1])
    bcavities = config["fesom"].get("use_cav", True)
    icb_path = config["general"]["thisrun_work_dir"]
    ibareamax = config["fesom"].get("ibareamax", 400)
    
    # Get previous year for reading previous year's output files
    prev_date = config["general"]["prev_date"]
    prev_year = prev_date.year if hasattr(prev_date, 'year') else int(str(prev_date)[:4])
    
    # Read fw and calving files directly from outdata/fesom/ (previous year's output)
    experiment_dir = config["general"]["experiment_dir"]
    outdata_fesom_dir = os.path.join(experiment_dir, "outdata", "fesom")
    
    fw_file = os.path.join(outdata_fesom_dir, f"fw.fesom.{prev_year}.nc")
    calving_file = os.path.join(outdata_fesom_dir, f"calving_AA.fesom.{prev_year}.nc")
    
    # print all inputs for sanity check
    print(f"*---> Inputs:")
    print(f" * mesh_dir: {mesh_dir}")
    print(f" * basin_file: {basin_file}")
    print(f" * icb_restart_file: {icb_restart_file}")
    print(f" * scaling_factor: {scaling_factor}")
    print(f" * bcavities: {bcavities}")
    print(f" * icb_path: {icb_path}")
    print(f" * ibareamax: {ibareamax}")
    print(f" * prev_year: {prev_year}")
    print(f" * fw_file: {fw_file}")
    print(f" * calving_file: {calving_file}")
    
    # Check if files exist
    if not os.path.exists(fw_file):
        raise FileNotFoundError(f"Freshwater flux file not found: {fw_file}")
    if not os.path.exists(calving_file):
        raise FileNotFoundError(f"Calving file not found: {calving_file}")
    
    # Use current year as seed for reproducible iceberg generation
    current_date = config["general"]["current_date"]
    seed_year = current_date.year if hasattr(current_date, 'year') else int(str(current_date)[:4])
    
    ib = IcebergCalving(
                        mesh_path=mesh_dir, 
                        icb_path=icb_path, 
                        latest_restart_file=icb_restart_file,
                        abg=[50,15,-90],
                        scaling_factor=scaling_factor,
                        seed=seed_year, 
                        bcavities=bcavities, 
                        ibareamax=ibareamax,
                        domain="SH",
                        basin_file=basin_file,               # basin file ?needed?
                        #calving_basin_id_file=fesom_basin_id_file,       # file containing the basin IDs for solid runoff/enthalpy on fesom grid
                        fw_file=fw_file,                  # fesom freshwater flux file
                        calving_file=calving_file                 # 
    )
    ib.create_dataframe()
    ib._icb_generator(fmode="w")
    
    return config

def apply_iceberg_calving_to_namelists(config):
    """
    Updates the namelist.config file with the correct number of icebergs.

    Inputs:
    - config: esm_tools configuration dictionary
    """

    # get config values for iceberg usage
    with_icb   = config["fesom"].get("with_icb", False)
    run_number = config["general"].get("run_number", 0)
    icb_path   = config["general"]["thisrun_work_dir"]

    if with_icb:
        print("\n*---> Counting new icebergs...")
        # Count the lines with data in icb_latitude file
        with open(os.path.join(icb_path, "icb_latitude.dat"), "r") as f:
            num_new_icebergs = sum(1 for line in f)
            print(f"* Number of new icebergs: {num_new_icebergs}")

        if os.path.isfile(os.path.join(icb_path, "iceberg.restart.ISM")):
            # Count the lines with data in iceberg.restart.ISM file
            print("\n*---> Counting old icebergs...")
            with open(os.path.join(icb_path, "iceberg.restart.ISM"), "r") as f:
                num_old_icebergs = sum(1 for line in f)
            print(f"* Number of old icebergs: {num_old_icebergs}")
        else:
            print("\n*---> No iceberg.restart.ISM file found, assuming no old icebergs")
            num_old_icebergs = 0

        with open(os.path.join(icb_path, "num_non_melted_icb_file"), 'w') as f:
            f.write(str(num_old_icebergs))
            print(f"* Created num_non_melted_icb_file with {num_old_icebergs} non-melted icebergs")

    # Open the namelist.config in the current work directory...
    nml_path = os.path.join(icb_path, "namelist.config")
    # And replace the ib_num entry with the sum of new and old icebergs
    print(f"\n*---> Updating namelist.config with ib_num = {num_new_icebergs + num_old_icebergs}")
    with open(nml_path, "r") as f:
        nml_content = f.read()
    nml_content = nml_content.replace("ib_num = 1", f"ib_num = {num_new_icebergs + num_old_icebergs}")
    with open(nml_path, "w") as f:
        f.write(nml_content)

    print("\n*---> Updating icebergs done!")

    return config