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

## Authors

Finn Ole Heukamp, 2026

## License

MIT License - see LICENSE file
