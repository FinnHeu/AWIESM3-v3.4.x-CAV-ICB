#!/bin/bash
# Run the iceberg plugin offline test
# 
# Usage:
#   ./run_test.sh
#
# This script:
# 1. Sources bashrc
# 2. Initializes conda
# 3. Activates the esm_tools conda environment
# 4. Runs the test

set -e  # Exit on error

echo "========================================"
echo "AWI-ESM-v3.4.1-CAV-ICB Test Runner"
echo "========================================"

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source bashrc if it exists
echo "* Sourcing bashrc..."
if [ -f "$HOME/.bashrc" ]; then
    source "$HOME/.bashrc"
elif [ -f "$HOME/.bash_profile" ]; then
    source "$HOME/.bash_profile"
fi

# Initialize conda for bash shell
echo "* Initializing conda..."
if command -v conda &> /dev/null; then
    # Conda is in PATH, ensure it's initialized
    eval "$(conda shell.bash hook)"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "/sw/spack-levante/miniconda3-24.3.0-0/etc/profile.d/conda.sh" ]; then
    # Levante cluster specific
    source "/sw/spack-levante/miniconda3-24.3.0-0/etc/profile.d/conda.sh"
else
    echo "WARNING: Could not find conda initialization script"
    echo "Attempting to proceed anyway..."
fi

# Activate esm_tools conda environment
echo "* Activating esm_tools conda environment..."
conda activate esm_tools

# Verify activation
echo "* Python version: $(python --version)"
echo "* Conda environment: $CONDA_DEFAULT_ENV"

# Change to test directory
cd "$SCRIPT_DIR"

# Run the test
echo ""
echo "========================================"
echo "Running iceberg plugin test..."
echo "========================================"
python test_iceberg_plugin.py

# Capture exit code
TEST_EXIT_CODE=$?

if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "========================================"
    echo "TEST COMPLETED SUCCESSFULLY"
    echo "========================================"
else
    echo ""
    echo "========================================"
    echo "TEST FAILED (exit code: $TEST_EXIT_CODE)"
    echo "========================================"
fi

exit $TEST_EXIT_CODE
