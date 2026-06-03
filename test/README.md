# AWI-ESM-v3.4.1-CAV-ICB Iceberg Plugin Test

This directory contains an offline test for the iceberg calving plugin that generates iceberg files from FESOM output.

## Files

- `test_iceberg_plugin.py` - Main test script (Python)
- `run_test.sh` - Shell script to run the test with proper conda environment setup
- `data/` - Test data directory containing:
  - `fw.fesom.1600.nc` - FESOM freshwater flux output
  - `calving_AA.fesom.1600.nc` - FESOM Antarctic calving output

## Requirements

- `esm_tools` conda environment with:
  - numpy
  - pandas
  - xarray
  - netcdf4
  - pyfesom2
  - powerlaw
  - numexpr
  - tqdm
  - scipy

## Running the Test

### Option 1: Using the shell script (recommended)

```bash
cd /home/a/a270186/esm_tools/plugins/AWI-ESM-v3.4.1-CAV-ICB/test
./run_test.sh
```

This automatically:
- Sources your bashrc
- Initializes conda
- Activates the `esm_tools` environment
- Runs the test

### Option 2: Manual execution

If your conda environment is already activated:

```bash
conda activate esm_tools
python test_iceberg_plugin.py
```

## What the Test Does

1. **Locates input files**: Finds the FESOM output files (`fw.fesom.*.nc`, `calving_AA.fesom.*.nc`) in `data/` and the basins file in `../basins/`

2. **Mesh handling**:
   - First tries to find the DARS2cav mesh at `/work/ab0995/a270186/model_inputs/fesom2/mesh/DARS2cav/`
   - If not available, creates a minimal synthetic mesh based on the input data

3. **Runs iceberg generation**:
   - Initializes `IcebergCalving` class with test data
   - Creates discharge dataframe from cavity melt and calving data
   - Generates iceberg distribution using power-law scaling
   - Writes output files to a temporary directory

4. **Verifies outputs**: Checks that all 7 iceberg output files are created:
   - `icb_longitude.dat` - Longitude positions
   - `icb_latitude.dat` - Latitude positions  
   - `icb_length.dat` - Iceberg lengths (m)
   - `icb_height.dat` - Iceberg heights/depths (m)
   - `icb_scaling.dat` - Scaling factors (number of icebergs per entry)
   - `icb_felem.dat` - FESOM element indices
   - `icb_calving_day.dat` - Calving day of year (1-365)

## Output

On successful completion, the test preserves the output directory and prints its location:

```
* Output files preserved in: /tmp/icb_test_output_XXXXXX
```

You can inspect the generated iceberg files there.

## Calving Day Seasonal Cycle

The plugin generates calving days for each iceberg with different distributions based on iceberg size:

### Large icebergs (scaling factor = 1)
- **Uniform random distribution** across all days (1-365)
- Represents the 3 largest size classes that calve year-round

### Small icebergs (scaling factor > 1)
- **Austral summer-weighted distribution** using von Mises (circular normal) distribution
- Centered on day 15 (mid-January) with κ=1.5
- Results in ~60% of small icebergs calving during austral summer (Nov-Feb)
- Represents smaller, more seasonal calving events

The von Mises distribution is ideal for cyclic data like day-of-year, ensuring smooth wrap-around at year boundaries.

## Troubleshooting

### Missing conda environment
If `esm_tools` environment is not found, check available environments:
```bash
conda env list
```

### Missing mesh files
If the real DARS2cav mesh is not available, the test automatically creates a minimal synthetic mesh. This may produce different iceberg distributions than the real mesh but validates the code execution.

### Import errors
Ensure all dependencies are installed in the `esm_tools` environment:
```bash
conda activate esm_tools
pip install pyfesom2 powerlaw numexpr tqdm
```
