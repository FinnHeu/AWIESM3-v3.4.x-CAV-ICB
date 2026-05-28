# AWI-ESM-v3.4.x-CAV-ICB

Iceberg Calving Plugin for AWIESM3-v3.4.x with ice shelf cavity support.

## Overview

This plugin generates iceberg input files for FESOM2 by computing the residual between sub-shelf cavity melt and available Antarctic solid runoff (accumulated snow on ANtarctica). The difference represents the mass flux that should be released as icebergs into the ocean.

The plugin distributes this mass flux across the calving fronts of the Antarctic cavities and generates icebergs following an observational power-law size distribution.

## Features

- Reads FESOM freshwater flux (`fw.fesom.<year>.nc`) and Antarctic calving (accumulated snow) output (`calving_AA.fesom.<year>.nc`)
- Computes the residual between cavity melt and available Antarctic solid runoff
- Distributes resulting iceberg calving across Antarctic drainage basins by calving front length in the respective basin
- Generates icebergs with realistic size distributions (power-law scaling)
- Assigns calving days with austral summer intensification for smaller icebergs
- Supports reproducible iceberg generation via year-based seeding
- Integrates seamlessly with `esm_tools` workflow

## Calving Day Seasonal Cycle

The plugin assigns a calving day (1-365) to each iceberg, controlling when FESOM starts processing it during the simulation year. The distribution depends on iceberg size class:

### Large icebergs (scaling factor = 1)
- **Uniform random distribution** across all days (1-365)
- Applies to the 3 largest size classes
- Represents year-round calving of large tabular icebergs

### Small icebergs (scaling factor > 1)
- **Austral summer-weighted distribution** using von Mises (circular normal) distribution
- Centered on day 15 (mid-January, peak austral summer) with concentration parameter κ=1.5
- Results in ~60% of small icebergs calving during austral summer (Nov-Feb)
- Represents the observed seasonal intensification of smaller calving events

The von Mises distribution ensures smooth wrap-around at year boundaries, making it ideal for cyclic day-of-year data.

## Installation

Install the plugin in the environment used by esm tools.

### From source (editable mode)

```bash
cd /path/to/esm_tools/plugins/AWI-ESM-v3.4.x-CAV-ICB
pip install -e .
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

AWIESM3-v3.4.x get their prepcompute recipe from the `oifs.yaml` file.
Add the following lines to the end of the `esm_tools/configs/components/oifs/oifs.yaml` file:

```
# Conditionally add iceberg plugin steps when with_icb is enabled
choose_fesom.with_icb:
   True:
      add_prepcompute_recipe:
         - "prep_icebergs"                       # Plugin
         - "apply_iceberg_calving_to_namelists"  # Plugin
```
## Authors

Finn Ole Heukamp, 2026

## License

MIT License - see LICENSE file
