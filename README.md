# AWI-ESM-v3.4.1-CAV-ICB

Iceberg Calving Plugin for AWI-ESM v3.4.1 with ice shelf cavity support.

## Overview

This plugin generates iceberg input files for FESOM2 simulations by computing the residual between sub-shelf cavity melt and available Antarctic solid runoff (calving). The difference represents the mass flux that should be released as icebergs into the ocean.

### Physical Background

In coupled ice-ocean simulations with ice shelf cavities:
- **Cavity melt**: Warm ocean water melts ice shelves from below, adding freshwater to the ocean
- **Solid runoff (calving)**: Ice sheet models provide calving fluxes at the ice front
- **Residual**: The difference between these fluxes determines the iceberg mass to be released

The plugin distributes this mass flux across Antarctic drainage basins and generates icebergs following an observational power-law size distribution.

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
cd /path/to/esm_tools/plugins/AWI-ESM-v3.4.1-CAV-ICB
pip install -e .
```

### Verify installation

```bash
pip show AWI-ESM-v3.4.1-CAV-ICB
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

## Plugin Functions

### `prep_icebergs(config)`

Main entry point that determines whether to use initial iceberg files or generate new ones.

**Workflow:**
- Run 1: Copy `icb_*.dat` files from `ini_iceberg_dir` to work directory
- Run 2+: Call `update_icebergs()` to generate new icebergs

### `update_icebergs(config)`

Generates iceberg files from FESOM output.

**Input files (from previous year's `outdata/fesom/`):**
- `fw.fesom.<prev_year>.nc` - Freshwater flux
- `calving_AA.fesom.<prev_year>.nc` - Antarctic calving

**Output files (written to work directory):**
| File | Description |
|------|-------------|
| `icb_latitude.dat` | Iceberg latitude positions |
| `icb_longitude.dat` | Iceberg longitude positions |
| `icb_length.dat` | Iceberg horizontal dimensions (m) |
| `icb_height.dat` | Iceberg thickness (m) |
| `icb_scaling.dat` | Scaling factors (represents multiple icebergs) |
| `icb_felem.dat` | FESOM mesh element indices |

### `apply_iceberg_calving_to_namelists(config)`

Updates FESOM namelist with the total iceberg count.

**Counts:**
- New icebergs from `icb_latitude.dat`
- Existing icebergs from `iceberg.restart.ISM` (if present)

**Updates:**
- `ib_num` in `namelist.config`

## IcebergCalving Class

The core `IcebergCalving` class handles the physics:

```python
IcebergCalving(
    mesh_path,              # FESOM mesh directory
    icb_path,               # Output directory for iceberg files
    latest_restart_file,    # Previous iceberg restart file
    abg=[50, 15, -90],      # Euler rotation angles
    scaling_factor,         # Per-basin scaling factors
    seed,                   # Random seed (year-based)
    bcavities,              # Use cavity melt (True/False)
    ibareamax,              # Maximum iceberg area (km²)
    domain,                 # Domain ("SH" for Southern Hemisphere)
    basin_file,             # Antarctic drainage basins NetCDF
    fw_file,                # FESOM freshwater flux file
    calving_file            # FESOM calving file
)
```

### Processing steps:

1. Read basin definitions from `basins_antarctica.nc`
2. Read and integrate FESOM freshwater flux (cavity melt)
3. Read and integrate FESOM calving flux
4. Compute residual: `iceberg_flux = calving - cavity_melt`
5. Distribute flux across basins weighted by calving front length
6. Generate icebergs following power-law size distribution
7. Assign icebergs to FESOM mesh elements
8. Write output `.dat` files

## File Structure

```
AWI-ESM-v3.4.1-CAV-ICB/
├── plugin.py                          # esm_tools plugin entry points
├── icb_apply_distribution_functions.py # IcebergCalving class
├── pyproject.toml                     # Package configuration
├── __init__.py
├── LICENSE
├── README.md
└── test/
    ├── test_iceberg_plugin.py         # Offline test script
    ├── run_test.sh                    # Test runner
    ├── plot_icb.py                    # Visualization script
    ├── data/                          # Test input data
    └── output/                        # Test output (gitignored)
```

## Testing

Run the offline test to verify plugin functionality:

```bash
cd test/
./run_test.sh
```

See `test/README.md` for detailed testing instructions.

## Troubleshooting

### FileNotFoundError for fw/calving files

Ensure the previous year's FESOM output exists in `outdata/fesom/`:
- `fw.fesom.<year>.nc`
- `calving_AA.fesom.<year>.nc`

### Plugin not executing

1. Verify installation: `pip show AWI-ESM-v3.4.1-CAV-ICB`
2. Check `with_icb: true` in runscript
3. Verify plugin steps in `oifs.yaml` `prepcompute_recipe`

### KeyError in config

Ensure all required FESOM config keys are set:
- `mesh_dir`
- `basin_file`
- `ini_iceberg_dir` (for first run)

## Authors

AWI-ESM Development Team

## License

MIT License - see LICENSE file
