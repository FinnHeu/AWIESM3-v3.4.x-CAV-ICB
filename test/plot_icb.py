"""
Plot initial iceberg locations near Antarctica.

Run this in the so_ase environment (NOT esm_tools):
    conda activate so_ase
    python plot_icb.py
"""
import sys
sys.path.insert(0, '/home/a/a270186/python_modules/SO-ASE')

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import so_ase as so

# Hardcoded paths (same as test_iceberg_plugin.py)
MESH_DIR = "/work/ab0995/a270186/model_inputs/fesom2/mesh/DARS2cav/"
OUTPUT_DIR = "/home/a/a270186/esm_tools/plugins/AWI-ESM-v3.4.1-CAV-ICB/test/output/"


def load_iceberg_data(output_dir):
    """Load iceberg location data from output files."""
    lon = np.loadtxt(f"{output_dir}/icb_longitude.dat")
    lat = np.loadtxt(f"{output_dir}/icb_latitude.dat")
    length = np.loadtxt(f"{output_dir}/icb_length.dat")
    return lon, lat, length


def load_mesh_and_cavity(mesh_dir):
    """Load FESOM mesh and cavity information."""
    mesh_diag_file = f"{mesh_dir}/fesom.mesh.diag.nc"
    cavity_file = f"{mesh_dir}/cavity_elvls.out"

    # Load mesh diagnostics
    mesh = xr.open_dataset(mesh_diag_file)

    # Load cavity flags
    cavity = np.loadtxt(cavity_file)
    cavity_mask = cavity > 1

    return mesh, cavity_mask


def plot_iceberg_locations(mesh_dir=MESH_DIR, output_dir=OUTPUT_DIR, output_file="iceberg_locations.png"):
    """
    Create a South Polar Stereo plot of iceberg locations and cavity elements.

    Parameters:
    -----------
    mesh_dir : str
        Path to FESOM mesh directory
    output_dir : str
        Path to iceberg output directory
    output_file : str
        Output filename for the plot
    """
    print("Loading iceberg data...")
    icb_lon, icb_lat, icb_length = load_iceberg_data(output_dir)
    print(f"  Loaded {len(icb_lon)} icebergs")

    print("Loading mesh and cavity data...")
    mesh, cavity_mask = load_mesh_and_cavity(mesh_dir)
    print(f"  Mesh: {mesh.sizes['nod2']} nodes, {mesh.sizes['elem']} elements")
    print(f"  Cavity nodes: {np.sum(cavity_mask)}")

    # Get coordinates
    lon = mesh['lon'].values
    lat = mesh['lat'].values

    # Get element connectivity (face_nodes gives 3 nodes per element)
    # face_nodes shape is (n3, elem) - need to transpose to (elem, 3)
    face_nodes = mesh['face_nodes'].values.T - 1  # Convert to 0-indexed

    # Create figure with South Polar Stereo projection
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo())

    # Set map extent (-90N to -60N)
    ax.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())

    # Add map features using so_ase functions if available, otherwise cartopy
    try:
        # Try so_ase land feature
        so.add_land(ax, color='0.8')
    except (AttributeError, TypeError):
        # Fallback to cartopy
        ax.add_feature(cfeature.LAND, facecolor='0.8', edgecolor='k', zorder=1)
        ax.add_feature(cfeature.COASTLINE, edgecolor='k', linewidth=0.5, zorder=2)

    # Plot cavity elements as tripcolor background
    print("Plotting cavity elements...")

    # Create a mask for cavity elements (elements where all 3 nodes are cavity nodes)
    elem_cavity = cavity_mask[face_nodes[:, 0]] & cavity_mask[face_nodes[:, 1]] & cavity_mask[face_nodes[:, 2]]

    # Plot all cavity elements with a light blue color
    tri = ax.tripcolor(lon, lat, face_nodes[elem_cavity],
                       np.ones(np.sum(elem_cavity)),
                       cmap='Blues', alpha=0.3, transform=ccrs.PlateCarree(),
                       vmin=0, vmax=1)

    # Scatter iceberg locations
    print("Plotting iceberg locations...")

    # Scale marker size by iceberg length
    sizes = np.clip(icb_length / 1000, 1, 100)  # Scale factor for visibility

    sc = ax.scatter(icb_lon, icb_lat, s=sizes, c=icb_length, cmap='viridis',
                    alpha=0.6, edgecolors='none', transform=ccrs.PlateCarree(),
                    zorder=5)

    # Add colorbar for iceberg size
    cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', pad=0.05, shrink=0.6)
    cbar.set_label('Iceberg Length (m)', fontsize=10)

    # Add gridlines
    gl = ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False,
                      linewidth=0.5, color='gray', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False

    # Title
    ax.set_title(f'Initial Iceberg Locations near Antarctica\n({len(icb_lon)} icebergs)', fontsize=12)

    # Save figure
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")

    # Close mesh dataset
    mesh.close()

    return fig, ax


if __name__ == "__main__":
    plot_iceberg_locations()