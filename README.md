# AWI-ESM-v3.4.x-CAV-ICB

Iceberg Calving Plugin for AWIESM3-v3.4.x with ice shelf cavity support.

## Overview

This plugin generates iceberg input files for FESOM2 by computing the residual between sub-shelf cavity melt and available Antarctic solid runoff (accumulated snow on ANtarctica). The difference represents the mass flux that should be released as icebergs into the ocean.

The plugin distributes this mass flux across the calving fronts of the Antarctic cavities and generates icebergs following an observational power-law size distribution.

## Features

- Reads FESOM freshwater flux (`fw.fesom.<year>.nc`) and calving output (`calving_AA.fesom.<year>.nc`)
- Computes annual mean cavity melt and calving fluxes
- Distributes iceberg calving across Antarctic drainage basins
- Generates icebergs with realistic size distributions (power-law scaling)
- Supports reproducible iceberg generation via year-based seeding
- Integrates seamlessly with `esm_tools` workflow

## Installation

Install the plugin in the environment used by esm tools.

### From source (editable mode)

```bash
cd /path/to/esm_tools/plugins/AWI-ESM-v3.4.x-CAV-ICB
pip install -e .
```

### Verify installation

```bash
pip show AWI-ESM-v3.4.x-CAV-ICB
```

## Dependencies

- `numpy`
- `pandas`
- `xarray`
- `netCDF4`
- `pyfesom2`
- `powerlaw`
- `numexpr`
- `tqdm`
- `scipy`

## Configuration

### Runscript settings

Enable the plugin in your AWI-ESM runscript:

```yaml
fesom:
    with_icb: true
    use_cav: true
    basin_file: "/path/to/basins_antarctica.nc"
    ini_iceberg_dir: "/path/to/initial/iceberg/files"
    scaling_factor: [1, 1, 1, 1, 1, 1]  # Optional: per-basin scaling
```

### esm_tools integration

The plugin registers two entry points that are automatically called during the `prepcompute` phase when `fesom.with_icb: true`:

1. **`prep_icebergs`**: Prepares iceberg files
   - First run: Copies initial iceberg files from `ini_iceberg_dir`
   - Subsequent runs: Generates new icebergs from previous year's output

2. **`apply_iceberg_calving_to_namelists`**: Updates FESOM namelist
   - Counts new and existing icebergs
   - Sets `ib_num` parameter in `namelist.config`

AWIESM3-v3.4.x get their prepcompute recipe from the oifs.yaml file.
Add the following to the oifs.yaml file:

```yaml
prepcompute_recipe:
       - "compile_model"
       - "_show_simulation_info"
       - "create_new_files"
       - "create_empty_folders"
       - "prepare_coupler_files"
       - "assemble"
       - "log_used_files"
       - "_write_finalized_config"
       - "wait_for_iterative_coupling"
       - "copy_files_to_thisrun"
       - "write_env"
       - "preprocess"
       - "modify_namelists"
       - "append_to_namelist"
       - "modify_files"
       - "copy_files_to_work"
       - "report_missing_files"
       - "compute_and_log_file_checksums"
       # see https://github.com/esm-tools/esm_tools/discussions/774
       # - "add_vcs_info"
       #- "check_vcs_info_against_last_run"
       - "database_entry"

    # Conditionally add iceberg plugin steps when with_icb is enabled
    choose_fesom.with_icb:
        True:
            add_prepcompute_recipe:
                - "prep_icebergs"                       # Plugin
                - "apply_iceberg_calving_to_namelists"  # Plugin
    # this does not work as expected, solution is work in progress
    #   - "postprocess"

## Authors

Finn Ole Heukamp, 2026

## License

MIT License - see LICENSE file
