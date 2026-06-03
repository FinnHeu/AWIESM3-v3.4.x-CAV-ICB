import os
import sys
import numpy.random as random
import math
import time
import numpy as np
import pandas as pd
import xarray as xr
import numexpr as ne
import pyfesom2 as pf
import powerlaw
import warnings
warnings.filterwarnings("ignore")

from tqdm import tqdm
from scipy.spatial import cKDTree
from scipy.stats import vonmises
from datetime import datetime

class IcebergCalving:
    def __init__(self,
                mesh_path,
                icb_path,
                latest_restart_file="",
                abg=[50,15,-90],
                scaling_factor=[1, 1, 1, 1, 1, 1],
                seed=0,
                bcavities=True,
                ibareamax=400,
                domain="SH",
                basin_file=None,               # basins file (basins_antarctica.nc)
                fw_file=None,                  # fesom freshwater flux file to compute cavity subshelf melt
                calving_file=None              # fesom Antarctic calving file (calving_AA.<year>.nc)
                ):
        # set seed for random number generation
        random.seed(seed)
       
        print(" * seed = ", seed)
        
        # Store input type for later processing
        self.icb_path = icb_path
        self.basin_file = basin_file
        self.fw_file = fw_file
        self.calving_file = calving_file
        self.mesh_path = mesh_path
        self.mesh_diag_file = f"{self.mesh_path}fesom.mesh.diag.nc"
        self.mesh_diag = xr.open_dataset(self.mesh_diag_file)
        self.node_area = self.mesh_diag['nod_area'].max(dim='nz')
        self.nod2d_file = os.path.join(mesh_path, "nod2d.out")
        self.elem2d_file = os.path.join(mesh_path, "elem2d.out")
        self.cavity_elvls_file = os.path.join(mesh_path, "cavity_elvls.out")
        self.latest_restart_file = latest_restart_file
        self.abg = abg
        self.rho_ice = 850
        self.rho_water = 1000
        self.bins = [0.1, 1, 10, 100, 1000]
        self.weights_area = [0.0005, 0.008, 0.025, 0.074, 0.893]
        self.weights_dist = [0.75, 0.175, 0.05, 0.02, 0.005]
        self.area_mean = [0.01, 0.1, 1, 10, 100]
        self.area_min = 0.01             #[km2]
        self.area_max = ibareamax           #[km2]
        self.min_disch_in_cell = 0.0   #[kg m-2 year-1]
        self.scaling_factor = np.array(scaling_factor)
        self.thick = np.array([0.25, 0.25, 0.25, 0.25, 0.25, 0.25])
        self.thick_max = 0.25
        self.depth = self.thick * 7/8
        self.height = self.thick - self.depth
        self.seed = seed


        self.domain    = domain
        self.bcavities = bcavities
        
    
        # Read basins file first (needed for weight distribution)
        if self.basin_file:
            self._read_basins_file()
        
        # FESOM freshwater flux section - new processing
        self._read_fesom_fw_file()
        self._apply_cavity_mask_and_integrate()
        self._read_fesom_calving_file()
        self._integrate_calving()
        self._compute_diff()
        self._convert_annual_mean_to_annual()
        self._convert_water_volume_to_ice_volume()
        self._distribute_by_basin_weights()
        self._get_fesom_coords()

        # FESOM mesh section
        self._read_mesh()
        self._read_nod2d_file()
        self._read_elem2d_file()
        
        if self.bcavities and os.path.isfile(self.cavity_elvls_file):
            self._read_cavity_elvls_file()

        # Check for restart file
        if not self.latest_restart_file=="" and os.path.exists(self.latest_restart_file):
            print(" * restart file ", self.latest_restart_file)
            self._get_full_cells()
        else:
            print("* no restart file found, continue anyway...")
            self.full_elems = []

    def create_dataframe(self):
        """
        Create a pandas DataFrame with the discharge data.
        """
        self._find_basins()
        
        if self.bcavities:
            # New approach: weights based on calving front element counts (proxy for calving front length)
            print(f" * total ice volume flux: {self.total_calving_flux.values * 1e-9:.2f} km3/year (ice)")
            
            # Create empty df_agg first
            self.df_agg = pd.DataFrame(index=self.basin_ids)
            
            # Find calving front elements (ocean elements adjacent to cavities) per basin
            self._find_calving_front_elements()
            
            # Calculate weights based on number of calving front elements per basin
            elem_counts = np.array([len(elems) for elems in self.df_agg["elems"]])
            total_elems = elem_counts.sum()
            
            if total_elems > 0:
                self.basin_weights = elem_counts / total_elems
            else:
                # Fallback to equal weights if no elements found
                self.basin_weights = np.ones(self.n_basins) / self.n_basins
            
            print(f"\n * Basin weights (based on calving front element counts):")
            for i, basin_id in enumerate(self.basin_ids):
                print(f"   Basin {basin_id}: {elem_counts[i]} elements, weight = {self.basin_weights[i]:.4f}")
            
            # Calculate per-basin discharge based on weights
            self.basin_discharge = {}
            for i, basin_id in enumerate(self.basin_ids):
                self.basin_discharge[basin_id] = self.total_calving_flux.values * self.basin_weights[i]
            
            # Update df_agg with discharge values
            self.df_agg["disch"] = [self.basin_discharge[bid] for bid in self.basin_ids]
            
            print(f"\n * Per-basin discharge (km3/y): {[f'{v * 1e-9:.2f}' for v in self.basin_discharge.values()]}")
        else:
            # Original approach: melt-weighted distribution
            self._find_FESOM_elem()
        
            disch_m3_year = self.data
            print(f" * total ice volume flux: {np.sum(disch_m3_year * 1e-9):.2f} km3/year (ice)")
            
            self.df = pd.DataFrame({
                    "disch": disch_m3_year, #[m3/year] 
                    "elems": self.indices1D,
                    "basin": self.basins1D,
                    })
            self.df.dropna(inplace=True)
            self.df_agg = self.df[["disch", "basin"]].groupby("basin").sum()

            elem_tmp = []
            neigh_tmp = []
            
            # iterate over basins
            for b in self.df.groupby("basin"):

                # get all FESOM elements nearest to discharge location within this basin 
                # every element shall occur only once
                elem_tmp.append(b[1]["elems"].unique())
                n=[]
                for x in elem_tmp[-1]:
                    
                    # get all neighbouring elements for each FESOM element
                    # every above found FESOM element is associated with a list of neighbouring elements
                    n.append(self._get_FESOM_neighbours(x))
                    
                    ################################################   
                    #tmp = self._get_FESOM_neighbours(x)
                    #m=[]
                    #for y in tmp:
                    #    m = m + self._get_FESOM_neighbours(y)
                    #n.append(np.unique(m))
                    ################################################
                neigh_tmp.append(n)
            
            self.df_agg["elems"] = elem_tmp
            self.df_agg["neigh."] = neigh_tmp
    
    def _read_mesh(self):
        self.mesh = pf.load_mesh(self.mesh_path, self.abg, usepickle=False)

    ###----------------------------------------------------------------------------###
    # Methods                                                                        #
    ###----------------------------------------------------------------------------###
    
    def _read_fesom_fw_file(self):
        """
        Read FESOM freshwater flux file (fw variable in m/s).
        Handles all number of timesteps, e.g. daily (365 timesteps), monthly (12 timesteps) and annual (1 timestep) means.
        Returns: fw_field (xarray.DataArray) as annual mean.
        """
        print(" \n*---> Reading FESOM freshwater flux file and computing annual mean")
        print(f" * Opening FESOM freshwater flux file: {self.fw_file}")
        self.ds_fw = xr.open_dataset(self.fw_file)
        
        # Get the fw variable
        self.fw_field = self.ds_fw['fw']
        
        # Check number of timesteps and compute annual mean
        n_time = len(self.ds_fw.time)
        print(f" * Number of timesteps in file: {n_time}")
        
        if n_time == 12:
            # Monthly means - compute annual mean with weights for seconds of the month
            print(f" * Detected {n_time} timesteps, computing annual mean based on weighted seconds per months")
            spm = seconds_per_month(self.seed)
            self.ds_fw['spm'] = (('time'), spm)
            self.fw_annual_mean = self.fw_field.weighted(self.ds_fw['spm']).mean(dim='time')
        elif n_time == 1:
            # Annual mean - use directly
            print(" * Detected annual data, using directly")
            self.fw_annual_mean = self.fw_field.squeeze(dim='time')
        elif (n_time == 365) | (n_time == 366):
            # Daily mean - use directly
            print(" * Detected daily data, using directly")
            self.fw_annual_mean = self.fw_field.squeeze(dim='time')
        else:
            print(f" * WARNING: unexpected number of timesteps ({n_time}), using mean over all")
            self.fw_annual_mean = self.fw_field.mean(dim='time')
        
    def _apply_cavity_mask_and_integrate(self):
        """
        Apply cavity mask to select only cavity gridcells.
        Uses fesom.mesh.diag.nc.
        """
        print(" \n*---> Applying cavity mask to select only cavity gridcells")
        print(" * Building cavity mask for nodes")
        
        # Build cavity mask for nodes
        self.cavity_mask = self.mesh_diag['nod_area'].isel(nz=0).values == 0
        n_cavity = np.sum(self.cavity_mask)
        print(f" * Found {n_cavity} cavity nodes out of {len(self.cavity_mask)} total nodes")
        
        # Apply mask to fw field - keep only cavity nodes
        print(" * Applying mask to freshwater field and node areas")
        self.fw_annual_mean_cavity = self.fw_annual_mean.values[self.cavity_mask]
        
        # Integrate over cavity nodes --> m/s per cell to total m3/s
        self.total_fw_cavity = -np.sum(self.fw_annual_mean_cavity * self.node_area[self.cavity_mask]) # m3/s
        print(f" * Total freshwater flux from cavities: {self.total_fw_cavity.values} m3/s (water)")

    def _read_fesom_calving_file(self):
        """
        Read FESOM calving file (calving variable in kg/(m2 * s)), which contains Antarctic runoff only.
        Handles all number of timesteps, e.g. daily (365 timesteps), monthly (12 timesteps) and annual (1 timestep) means.
        Returns: calving flux (xarray.DataArray) as annual mean.
        """
        print(" \n*---> Reading FESOM Antarctic calving flux file and computing annual mean")
        print(f" * Opening FESOM calving flux file: {self.calving_file}")
        self.ds_calving = xr.open_dataset(self.calving_file)
        
        # Get the runoff_solid variable
        self.calving_field = self.ds_calving['calving_AA']
        
        # Check number of timesteps and compute annual mean
        n_time = len(self.ds_calving.time)
        print(f" * Number of timesteps in file: {n_time}")
        
        if n_time > 1:
            # Monthly/Daily/... means - compute annual mean
            print(f" * Detected {n_time} timesteps, computing annual mean")
            self.calving_annual_mean = self.calving_field.mean(dim='time')
        elif n_time == 1:
            # Annual mean - use directly
            print(" * Detected annual data, using directly")
            self.calving_annual_mean = self.calving_field.squeeze(dim='time')
        else:
            print(f" * WARNING: unexpected number of timesteps ({n_time}), using mean over all")
            self.calving_annual_mean = self.calving_field.mean(dim='time')

    def _integrate_calving(self):
        """
        Integrate (Antarctica) calving.
        """
        print(" \n*---> Integrate Antarctic calving")
        
        # Integrate
        self.total_calving_basin_66 = np.sum(self.calving_annual_mean * self.node_area) # m3/s
        print(f" * Total calving flux in basin 66: {self.total_calving_basin_66.values} m3/s (water)")

    def _compute_diff(self):
        """
        Compute the difference between the total solid runoff in basin 66 and the total freshwater flux in cavities.
        """
        print(" \n*---> Computing residual between the total solid runoff in basin 66 and the total freshwater flux in cavities")
        self.residual = np.abs(self.total_calving_basin_66) - np.abs(self.total_fw_cavity)
        print(f" * Residual: {self.residual.values} m3/s (water)")

        if self.residual < 0:
            print(" ***********************************************")
            print(" * WARNING:Residual is negative, setting to 0! *")
            print(" ***********************************************")
            self.residual = 0

    def _convert_annual_mean_to_annual(self):
        """
        Convert annual mean to total annual value.
        """
        
        print(" \n*---> Converting annual mean to annual")
        year = self.ds_calving.time.dt.year[0].values
        if (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0):
            days_in_year = 366
        else:
            days_in_year = 365
        self.residual_annual = self.residual * days_in_year * 86400

        print(f" * Total annual calving flux: {self.total_calving_basin_66.values * days_in_year * 86400 * 1e-9} Gt/y")
        print(f" * Total annual subshelf melt flux: {self.total_fw_cavity.values * days_in_year * 86400 * 1e-9} Gt/y")
        print(f" * Total annual residual iceberg flux: {self.residual_annual.values} m3/y (water) or {self.residual_annual.values * 1e-9} Gt/y")


    def _convert_water_volume_to_ice_volume(self):
        """
        Convert water volume to ice volume.
        """
        print(" \n*---> Converting water volume to ice volume")
        self.total_calving_flux = self.residual_annual * (self.rho_water / self.rho_ice)
        print(f" * Total iceberg flux <before> conversion to ice volume: {self.residual_annual.values} m3/y (water) or {self.residual_annual.values * 1e-9} Gt/y")
        print(f" * Total iceberg flux <after> conversion to ice volume: {self.total_calving_flux.values} m3/y (ice)")

    def _distribute_by_basin_weights(self):
        """
        Identify unique basins from the basins file.
        Weight calculation is deferred to create_dataframe() after calving front elements are found.
        """
        print(" \n*---> Identifying basins for flux distribution")
        
        # Find all unique basins from the basins file
        unique_basins = np.unique(self.basins.values)
        unique_basins = unique_basins[~np.isnan(unique_basins)]  # Remove NaN
        unique_basins = unique_basins[unique_basins > 0]  # Keep only valid basins (> 0)
        
        self.n_basins = len(unique_basins)
        self.basin_ids = unique_basins.astype(int)
        
        print(f"   Number of basins: {self.n_basins}")
        print(f"   Basin IDs: {self.basin_ids}")
        print(f" * Total iceberg flux: {self.total_calving_flux.values * 1e-9:.2f} km3/y (ice)") 
        
    def _get_fesom_coords(self):
            """
            Get lon/lat coordinates for calving nodes from mesh diagnostics.
            """
            print(" \n*---> Extracting lon/lat coordinates for cavity nodes")
            # Get lon/lat from mesh diagnostics
            print(" * Reading coordinates from mesh diagnostics")
            all_lons = self.mesh_diag['lon'].values
            all_lats = self.mesh_diag['lat'].values
            
            # Select coordinates for calving nodes
            print(" * Selecting coordinates for cavity node indices")
            self.lons = list(all_lons[self.cavity_mask])
            self.lats = list(all_lats[self.cavity_mask])
            
            # Normalize longitudes to -180 to 180 range
            print(" * Normalizing longitudes to -180 to 180 range")
            self.lons = [lon if lon < 180 else lon - 360 for lon in self.lons]
            
            print(f" * Extracted coordinates for {len(self.lons)} cavity locations")

    def _read_basins_file(self):
        fl = xr.open_dataset(self.basin_file)
        if "basin" in fl:
            self.basins = fl.squeeze().basin
        elif "basins" in fl:
            self.basins = fl.squeeze().basins
        else:
            print("No basins in basin file")
            return -1
    
    def _get_nearest_lon_lat(self, ds, lon, lat):
        #https://stackoverflow.com/questions/58758480/xarray-select-nearest-lat-lon-with-multi-dimension-coordinates
        abslat = np.abs(ds.lat-lat)
        abslon = np.abs(ds.lon-lon)
        c = np.maximum(abslon, abslat)
    
        # Find minimum - may have multiple matches, take first one
        min_indices = np.where(c == np.min(c))
        yloc = min_indices[0][0]
        xloc = min_indices[1][0]
        point_ds = ds.isel(x=xloc, y=yloc)
        return point_ds

    def _find_basins(self):
        """
        Find basins for all cavity nodes using fast KDTree nearest neighbor search.
        """
        # Build KDTree from basin coordinates
        basin_lons = self.basins.lon.values.flatten()
        basin_lats = self.basins.lat.values.flatten()
        basin_coords = np.column_stack([basin_lons, basin_lats])

        # Create query points from cavity node coordinates
        node_coords = np.column_stack([self.lons, self.lats])

        # Build tree and query all points at once
        tree = cKDTree(basin_coords)
        _, indices = tree.query(node_coords, k=1)

        # Get basin values for nearest neighbors
        basin_flat = self.basins.values.flatten()
        self.basins1D = [int(basin_flat[i]) for i in indices]

    def _read_nod2d_file(self):
        self.nod2d = pd.read_csv(self.nod2d_file, header=0, names=["lon", "lat", "coastal"], sep='\s+', index_col=0)
        
    def _read_elem2d_file(self):
        self.elem2d = pd.read_csv(self.elem2d_file, header=0, names=["nod1", "nod2", "nod3"], sep='\s+')
        self._load_elem_adjacency_from_mesh()

    def _load_elem_adjacency_from_mesh(self):
        """
        Build element adjacency from node connectivity.
        Two elements are neighbors if they share exactly 2 nodes (an edge).
        Caches result to mesh folder for faster subsequent runs.
        """
        cache_file = os.path.join(self.mesh_path, "elem_adjacency.npz")
        
        # Try to load from cache
        if os.path.exists(cache_file):
            print(f" * Loading element adjacency from cache: {cache_file}")
            data = np.load(cache_file, allow_pickle=True)
            self._elem_neighbors = data['elem_neighbors'].item()
            print(f" * Loaded adjacency for {len(self._elem_neighbors)} elements")
            return
        
        print(" * Building element adjacency from node connectivity...")

        nelems = len(self.elem2d)
        n1 = self.elem2d['nod1'].values
        n2 = self.elem2d['nod2'].values
        n3 = self.elem2d['nod3'].values

        # Build node-to-elements mapping
        node_to_elems = {}
        for elem_idx in range(nelems):
            for node in (n1[elem_idx], n2[elem_idx], n3[elem_idx]):
                if node not in node_to_elems:
                    node_to_elems[node] = []
                node_to_elems[node].append(elem_idx)

        # Build adjacency: elements sharing 2 nodes are neighbors
        self._elem_neighbors = {}
        for elem_idx in range(nelems):
            # Get all elements connected to any of this element's nodes
            candidate_elems = set()
            for node in (n1[elem_idx], n2[elem_idx], n3[elem_idx]):
                candidate_elems.update(node_to_elems[node])
            candidate_elems.discard(elem_idx)  # Remove self

            # Keep only those sharing exactly 2 nodes (share an edge)
            neighbors = []
            elem_nodes = {n1[elem_idx], n2[elem_idx], n3[elem_idx]}
            for other_idx in candidate_elems:
                other_nodes = {n1[other_idx], n2[other_idx], n3[other_idx]}
                shared_nodes = elem_nodes & other_nodes
                if len(shared_nodes) == 2:  # Share an edge
                    neighbors.append(other_idx)

            self._elem_neighbors[elem_idx] = neighbors

        print(f" * Adjacency built for {len(self._elem_neighbors)} elements")
        
        # Save to cache
        print(f" * Saving element adjacency to cache: {cache_file}")
        np.savez(cache_file, elem_neighbors=self._elem_neighbors)

    def _read_cavity_elvls_file(self):
        tmp = pd.read_csv(self.cavity_elvls_file, names=["cavity"], sep='\s+')
        self.cavity_flags = tmp[tmp>1].notna()

    def _get_full_cells(self):
        """
        Identify mesh elements that already host an iceberg from the restart file.
        
        FESOM only allows one iceberg per mesh element (cell_saturation=4), so any element containing
        an existing iceberg must be excluded from seeding new icebergs.
        
        Column 18 in the restart file contains the FESOM element ID (1-indexed).
        """
        df = pd.read_csv(self.latest_restart_file, header=None, delim_whitespace=True)
        
        # Get unique element IDs from column 18 (1-indexed in restart file)
        occupied_elems = df[18].unique()
        
        # Filter valid element IDs and convert to 0-indexed
        full_elems_tmp = []
        for felem in occupied_elems:
            felem_int = int(felem)
            # Skip invalid element IDs (0 or out of range)
            if felem_int == 0 or felem_int > len(self.mesh.voltri):
                continue
            # Convert to 0-indexed
            full_elems_tmp.append(felem_int - 1)
        
        print(f" * Found {len(full_elems_tmp)} occupied elements (excluded from seeding)")
        self.full_elems = full_elems_tmp

    def _remove_cavities(self):
        print(" *--> remove cavity elements")
        #self.elem2d = self.elem2d[~self.cavity_flags]
        #self.elem2d = self.elem2d[~self.cavity_flags.values].reset_index()[["nod1", "nod2", "nod3"]]
        self.elem2d = self.elem2d[~self.cavity_flags.values][["nod1", "nod2", "nod3"]]
        print(" * list of good elements: ", self.elem2d.index.values[:])

    def _find_FESOM_elem(self):
        """
        Find nearest FESOM element for each cavity node using KDTree acceleration.
        """
        print(" * Building KDTree for element centroids...")

        # Get element vertex coordinates
        lon1 = self.nod2d.lon[self.elem2d.nod1].values
        lat1 = self.nod2d.lat[self.elem2d.nod1].values
        lon2 = self.nod2d.lon[self.elem2d.nod2].values
        lat2 = self.nod2d.lat[self.elem2d.nod2].values
        lon3 = self.nod2d.lon[self.elem2d.nod3].values
        lat3 = self.nod2d.lat[self.elem2d.nod3].values

        # Compute triangle centroids
        elem_centroids_lon = (lon1 + lon2 + lon3) / 3.0
        elem_centroids_lat = (lat1 + lat2 + lat3) / 3.0
        centroid_coords = np.column_stack([elem_centroids_lon, elem_centroids_lat])

        # Build KDTree on centroids
        tree = cKDTree(centroid_coords)

        # Query points (cavity node locations)
        node_coords = np.column_stack([self.lons, self.lats])

        # Find K nearest centroid candidates for each query point
        K = 20  # Check top 20 closest centroids
        print(f" * Querying KDTree for {len(node_coords)} points (K={K})...")
        _, candidate_indices = tree.query(node_coords, k=min(K, len(centroid_coords)))

        # For each query point, find the actual closest triangle among candidates
        points = []
        indices = []

        with tqdm(total=len(self.lons), file=sys.stdout, desc='find FESOM elements') as pbar:
            for i, (lon, lat) in enumerate(zip(self.lons, self.lats)):
                # Get candidate elements for this query point
                cands = candidate_indices[i] if K > 1 else [candidate_indices[i]]

                # Extract coordinates for candidate elements
                c_lon1, c_lat1 = lon1[cands], lat1[cands]
                c_lon2, c_lat2 = lon2[cands], lat2[cands]
                c_lon3, c_lat3 = lon3[cands], lat3[cands]

                # Compute summed distance to all 3 vertices for each candidate
                d1 = (c_lon1 - lon)**2 + (c_lat1 - lat)**2
                d2 = (c_lon2 - lon)**2 + (c_lat2 - lat)**2
                d3 = (c_lon3 - lon)**2 + (c_lat3 - lat)**2
                total_dist = d1 + d2 + d3

                # Find closest candidate
                min_idx = np.argmin(total_dist)
                closest_elem_global_idx = cands[min_idx]

                # Store result
                p1 = point(lon1[closest_elem_global_idx], lat1[closest_elem_global_idx])
                p2 = point(lon2[closest_elem_global_idx], lat2[closest_elem_global_idx])
                p3 = point(lon3[closest_elem_global_idx], lat3[closest_elem_global_idx])
                points.append([p1, p2, p3])
                indices.append(self.elem2d.index[closest_elem_global_idx])

                pbar.update(1)

        self.points1D = points
        self.indices1D = indices

    def _get_FESOM_neighbours(self, ind):
        """
        Fast lookup of neighboring elements using pre-built adjacency map.
        """
        return self._elem_neighbors.get(ind, [])

    def _find_calving_front_elements(self):
        """
        Find ocean elements suitable for iceberg seeding near the calving front.
        
        Strategy:
        1. Find calving front elements (ocean elements adjacent to cavity elements)
        2. Find all elements sharing 1-2 nodes with calving front elements
        3. Keep only those that are pure ocean elements (not cavity)
        
        This gives us "second row" elements that are safely in open ocean,
        avoiding FESOM's rejection of elements with any cavity_depth != 0 nodes.
        
        Ocean elements: cavity_elvls == 1
        Cavity elements: cavity_elvls > 1
        
        Modifies self.df_agg to contain the seeding elements per basin.
        """
        # Use new cache file name since logic changed
        cache_file = os.path.join(self.mesh_path, "elem_icb_seeding.npz")
        
        # Read raw cavity levels (1 = ocean, > 1 = cavity)
        lev_cav = pd.read_csv(self.cavity_elvls_file, names=["cavity"], sep='\s+')['cavity'].values
        is_ocean = lev_cav == 1
        is_cavity = lev_cav > 1
        
        nelems = len(self.elem2d)
        
        # Get element node connectivity (1-indexed to 0-indexed)
        n1 = self.elem2d['nod1'].values - 1
        n2 = self.elem2d['nod2'].values - 1
        n3 = self.elem2d['nod3'].values - 1
        
        # Try to load from cache
        if os.path.exists(cache_file):
            print(f" * Loading iceberg seeding elements from cache: {cache_file}")
            data = np.load(cache_file, allow_pickle=True)
            all_calving_front = data['all_calving_front']
            all_seeding_elems = data['all_seeding_elems']
            elems_of_node = data['elems_of_node']
            print(f"   Total elements: {nelems}, Ocean: {is_ocean.sum()}, Cavity: {is_cavity.sum()}")
            print(f"   Calving front elements: {len(all_calving_front)}")
            print(f"   Iceberg seeding elements: {len(all_seeding_elems)}")
        else:
            print(" * Finding calving front and iceberg seeding elements...")
            print(f"   Total elements: {nelems}, Ocean: {is_ocean.sum()}, Cavity: {is_cavity.sum()}")
            
            # Build elems_of_node: for each node, list of elements containing that node
            num_nodes = len(self.nod2d)
            elems_of_node = [[] for _ in range(num_nodes)]
            for eidx in range(nelems):
                elems_of_node[n1[eidx]].append(eidx)
                elems_of_node[n2[eidx]].append(eidx)
                elems_of_node[n3[eidx]].append(eidx)
            
            # Step 1: Find calving front elements (ocean elements adjacent to cavity elements)
            ocean_near_cavity = np.zeros(nelems, dtype=bool)
            for eidx in range(nelems):
                if is_ocean[eidx]:
                    elem_nodes = (n1[eidx], n2[eidx], n3[eidx])
                    for node in elem_nodes:
                        if any(is_cavity[e] for e in elems_of_node[node]):
                            ocean_near_cavity[eidx] = True
                            break
            
            all_calving_front = np.where(ocean_near_cavity)[0]
            calving_front_set = set(all_calving_front)
            print(f"   Calving front elements: {len(all_calving_front)}")
            
            # Step 2: Find elements sharing 1-2 nodes with calving front elements
            # These are potential seeding elements
            potential_seeding = set()
            for cf_eidx in all_calving_front:
                cf_nodes = (n1[cf_eidx], n2[cf_eidx], n3[cf_eidx])
                for node in cf_nodes:
                    # All elements containing this node
                    for neighbor_eidx in elems_of_node[node]:
                        # Skip if it's a calving front element itself
                        if neighbor_eidx not in calving_front_set:
                            potential_seeding.add(neighbor_eidx)
            
            # Step 3: Keep only pure ocean elements (not cavity)
            all_seeding_elems = np.array([e for e in potential_seeding if is_ocean[e]], dtype=np.int64)
            print(f"   Iceberg seeding elements (ocean, not calving front): {len(all_seeding_elems)}")
            
            # Save to cache
            print(f" * Saving iceberg seeding elements to cache: {cache_file}")
            np.savez(cache_file, 
                     all_calving_front=all_calving_front, 
                     all_seeding_elems=all_seeding_elems,
                     elems_of_node=np.array(elems_of_node, dtype=object))
        
        # Build node-to-basin mapping for fast lookup
        # self.basins1D was created by _find_basins() for each cavity node
        cavity_nodes_all = np.where(self.cavity_mask)[0]  # Global node indices of cavity nodes
        node_to_basin = {}
        for i, node_idx in enumerate(cavity_nodes_all):
            node_to_basin[node_idx] = self.basins1D[i]
        
        # Build calving front to basin mapping first
        # A calving front element belongs to basin B if it's adjacent to a cavity element of basin B
        calving_front_to_basins = {eidx: set() for eidx in all_calving_front}
        for eidx in all_calving_front:
            elem_nodes = (n1[eidx], n2[eidx], n3[eidx])
            for node in elem_nodes:
                for cav_eidx in elems_of_node[node]:
                    if is_cavity[cav_eidx]:
                        cav_nodes = (n1[cav_eidx], n2[cav_eidx], n3[cav_eidx])
                        for cav_node in cav_nodes:
                            if cav_node in node_to_basin:
                                calving_front_to_basins[eidx].add(node_to_basin[cav_node])
        
        # Assign seeding elements to basins based on which calving front elements they neighbor
        seeding_by_basin = {basin_id: set() for basin_id in self.basin_ids}
        calving_front_set = set(all_calving_front)
        
        for seed_eidx in all_seeding_elems:
            seed_nodes = (n1[seed_eidx], n2[seed_eidx], n3[seed_eidx])
            for node in seed_nodes:
                for neighbor_eidx in elems_of_node[node]:
                    if neighbor_eidx in calving_front_to_basins:
                        # This seeding element neighbors a calving front element
                        for basin_id in calving_front_to_basins[neighbor_eidx]:
                            if basin_id in seeding_by_basin:
                                seeding_by_basin[basin_id].add(seed_eidx)
        
        # For each basin, collect seeding elements
        elem_tmp = []
        neigh_tmp = []
        
        for basin_idx in self.basin_ids:
            seeding_list = sorted(list(seeding_by_basin[basin_idx]))
            elem_tmp.append(seeding_list)
            neigh_tmp.append([])  # No neighbors needed
            
            print(f"   Basin {basin_idx}: {len(seeding_list)} iceberg seeding elements")
        
        self.df_agg["elems"] = elem_tmp
        self.df_agg["neigh."] = neigh_tmp

    def _generate_calving_day(self, scaling_factor):
        """
        Generate a calving day (1-365) for an iceberg based on its scaling factor.
        
        For large icebergs (scaling_factor == 1): 
            Uniform random distribution across all days.
            
        For smaller icebergs (scaling_factor > 1):
            Austral summer-weighted distribution using von Mises distribution
            centered on mid-January (day ~15), with higher probability during
            December-February (austral summer).
        
        Parameters
        ----------
        scaling_factor : int
            The scaling factor for the iceberg (1 = large, >1 = smaller)
            
        Returns
        -------
        int
            Calving day of year (1-365)
        """
        if scaling_factor == 1:
            # Large icebergs: uniform random distribution
            return random.randint(1, 366)  # 1-365 inclusive
        else:
            # Smaller icebergs: austral summer intensification
            # Use von Mises distribution centered on day 15 (mid-January)
            # kappa controls concentration (higher = more peaked)
            # kappa=1.5 gives moderate summer preference while still having some winter calving
            
            # Convert day of year to angle (0 = Jan 1, 2*pi = Dec 31)
            # Center on day 15 (mid-January) = peak austral summer
            mu_day = 15  # Peak calving around mid-January
            mu_angle = (mu_day / 365.0) * 2 * np.pi
            
            # kappa controls how concentrated the distribution is
            # Higher kappa = more concentrated around summer
            # kappa=1.5 gives ~60% of calving in Nov-Feb, ~40% rest of year
            kappa = 1.5
            
            # Draw from von Mises distribution
            angle = vonmises.rvs(kappa, loc=mu_angle)
            
            # Convert angle back to day of year
            day = int((angle % (2 * np.pi)) / (2 * np.pi) * 365) + 1
            
            # Ensure day is in valid range
            day = max(1, min(365, day))
            
            return day

    def _create_icebergs_within_basin(self, df, idx):
    ######################################
    # input:    data frame for one basin: discharge and FESOM cell corners (p1, p2, p3)
    # output:   iceberg volume array
    ######################################
        # maximal time in seconds to wait for iceberg generation to finish (for one basin)
        dtime_MAX = 30
        
        # mu and sigma for lognormal distribution after Tournadre et al. (2011)
        mu, sigma = 12.3, 1.55**0.5
        xmin = self.area_min

        ############################################################
        if self.domain.lower() == "sh":
            # Southern Hemisphere: alpha for powerlaw after Tournadre et al. (2015) 
            a = 1.52
        elif self.domain.lower() == "greenland":
            # Greenland: alpha for powerlaw after Shiggins et al. (2023)
            if idx == 3:    # SKJI and UI
                a = 2.02
            elif idx == 5:  # KNS
                a = 2.38
            else:
                a = 2.2 
        else:
            a, xmin = 1.52, 0.01
        print(" * domain is {}, using alpha={} for basin {}".format(self.domain, a, idx))
        ############################################################


        median = 2**(1/(a-1))*xmin

        # get values within basin
        vals = abs(df.disch)
    
        # get total discharge within basin in [km3 year-1]
        disch_tot = vals / 1e9
        print(f" * total ice discharge within basin: {np.round(disch_tot, 2)} km3 year-1")
    
        # get total iceberg area within basin in [km2 year-1]
        # assuming constant iceberg height
        area_tot = disch_tot / self.thick_max
    
        # create iceberg areas according to Tournadre et al. (2015)
        # divide icebergs into classes of different area sizes (0.1-1, 1-10, 10-100, ... [km2])
        # and draw from powerlaw distribution with alpha=1.52 except for the icebergs from
        # smallest class. Get total number of icebergs with share of smallest class (WEIGHTS_N)
        # and mean size within smalles class (SMEAN_1). Get number of icebergs of each other class
        # with corresponding share. 
        N = int(area_tot / median)

        # generates random variates of power law distribution
        vrs = powerlaw.Power_Law(xmin=xmin, xmax=self.area_max, parameters=[a]).generate_random(N)

        x = vrs
        corr = area_tot / sum(x)
        
        x = x * corr

        # correction with respect to iceberg volume and not iceberg area
        thick = x**(1/2)
        thick[thick>self.thick_max] = self.thick_max
        vol = x * thick
        corr = disch_tot / sum(vol) 
        #corr = area_tot / sum(x)
        
        x = x * corr
        
        # correction with respect to iceberg volume and not iceberg area
        thick = x**(1/2)
        thick[thick>self.thick_max] = self.thick_max
        vol = x * thick
        vol_sum_0 = sum(vol)
        #x_sum_0 = sum(x)
        
        x = x[x>=xmin]
        x = x[x<=self.area_max]
        if sum(x) == 0:
            print(" * no icebegs")
            return pd.DataFrame()
        
        # correction with respect to iceberg volume and not iceberg area
        thick = x**(1/2)
        thick[thick>self.thick_max] = self.thick_max
        vol = x * thick
        vol_sum_1 = sum(vol)
        #x_sum_1 = sum(x)
        corr = vol_sum_0 / vol_sum_1
        #corr = x_sum_0 / x_sum_1
        x = x * corr

        x_tot = x
        tstart = time.time()
        while N > 0:
            N = int((vol_sum_0 - vol_sum_1) / vol_sum_0 * N)
            print(" N = ", str(N))
            if N==0:
                break
            vrs = powerlaw.Power_Law(xmin=xmin, xmax=self.area_max, parameters=[a]).generate_random(N)
            x = vrs
            x_too_small = x[x<xmin]
            x_too_large = x[x>self.area_max]
            x = x[x>=xmin]
            x = x[x<=self.area_max]
            x_tot = np.concatenate([x_tot, x])
            if sum(x_tot) == 0:
                print(" * no icebegs")
                return pd.DataFrame()
        
            # correction with respect to iceberg volume and not iceberg area
            thick = x_tot**(1/2)
            thick[thick>self.thick_max] = self.thick_max
            vol = x_tot * thick
            corr = disch_tot / sum(vol)
            #corr = area_tot / sum(x_tot)
            x_tot = x_tot * corr
            tend = time.time()
            dtime = tend - tstart
            if dtime >= dtime_MAX:
                print("elapsed time = ", str(tend - tstart))
                print("start iceberg generation again for this basin")
                return -1
            
        # correction with respect to iceberg volume and not iceberg area
        vol = x_tot * thick

        area = x_tot
        bins = np.digitize(area, self.bins, right=True)

        #for a, b in zip(area, bins):
        #    vol = np.concatenate([vol, [a * self.thick[b]]])

        # create data frame with iceberg elements: area, volume, bin
        ib_elems = pd.DataFrame({"area": area, 
                                "volume": vol,
                                "bin": np.digitize(area, self.bins, right=True)})
    
        ib_elems_ = ib_elems.where(ib_elems.area >= self.area_min).dropna()
    
#        print("*** Check for validity:")
#        print("***      assumed iceberg thickness [km]:         ", self.thick)
#        print("***      total discharge [km3 year-1]:        ", disch_tot)
#        print("***      summed iceberg volume [km3]:         ", sum(ib_elems.volume))
#        print("***      summed iceberg volume [km3]:         ", sum(ib_elems_.volume))
#        print("***      total iceberg area [km2 year-1]:     ", area_tot)
#        print("***      summed area (of generated ib) [km2]: ", sum(ib_elems.area)) 
#        print("***      summed area (of generated ib) [km2]: ", sum(ib_elems_.area)) 
#        print("***      total number of icebergs:            ", len(ib_elems))
#        print("***      total number of icebergs:            ", len(ib_elems_))
#        print(ib_elems_)
        return ib_elems_

    def _scale_icebergs(self, df):
    ######################################
    # input:    data frame: area, volume, bin
    # output:   data frame: length, scaling, depth
    ######################################
        # loop over all bins
        with tqdm(total=len(self.scaling_factor), file=sys.stdout, desc='go through all bins') as pbar:
            for i, (s, d) in enumerate(zip(self.scaling_factor, self.thick)):
                
                # get icb elements of particular size class
                ib_bin = df.where(df.bin==i).dropna()
    
                if not ib_bin.empty:
                    # split iceberg array of size class into chunks with length s
                    chunks = np.array_split(ib_bin, math.ceil(len(ib_bin)/s))
                    # get mean of each chunk
                    chunks_mean_area = np.array([chunk.area.mean(axis=0) for chunk in chunks])
                    chunks_mean_volume = np.array([chunk.volume.mean(axis=0) for chunk in chunks])
    
                    # check if arrays are initialized
                    if not 'length' in locals():
                        # get mean length of icebergs for each chunk
                        length = ne.evaluate('chunks_mean_area**(1/2)')
                        # get scaling factor (length of each chunk)
                        scaling = np.array([len(chunk) for chunk in chunks])
                        # get mean height of icebergs for each chunk
                        depth = ne.evaluate('chunks_mean_volume / chunks_mean_area')
                        ## get depth
                        #depth = np.array([d] * len(chunks))
    
                    else:
                        # get mean length of icebergs for each chunk
                        length = np.append(length, ne.evaluate('chunks_mean_area**(1/2)'))
                        # get scaling factor (length of each chunk)
                        scaling = np.append(scaling, np.array([len(chunk) for chunk in chunks]))
                        # get mean height of icebergs for each chunk
                        depth = np.append(depth, ne.evaluate('chunks_mean_volume / chunks_mean_area'))
                        ## get depth
                        #depth = np.append(depth, np.array([d] * len(chunks)))
                    pbar.update(1)
                else:
                    print("*** bin is empty")
                    pbar.update(1)
        
        # create data frame with scaled iceberg elements: length, scaling, depth
        df_out = pd.DataFrame({"length": length,
                                "scaling": scaling,
                                "depth": depth})
        
#        print("*** Check for validity:")
#        print("***      BEFORE SCALING:")
#        print("***      total iceberg area [km2]:   ", df.sum(axis=0).area)
#        print("***      total iceberg volume [km3]: ", df.sum(axis=0).volume)
#        print("***      total amount of icebergs:   ", len(df))
#        print("***      AFTER SCALING:")
#        print("***      total iceberg area [km2]:   ", np.sum(df_out.length * df_out.length * df_out.scaling))
#        print("***      total iceberg volume [km3]: ", np.sum(df_out.length * df_out.length * df_out.scaling * df_out.depth))
#        print("***      total amount of icebergs:   ", np.sum(df_out.scaling))
#        print("***      total am. of sim. icebergs: ", len(df_out))
        return df_out
    
    #generate icebergs
    def _icb_generator(self, fmode="w"):
        ###############################
        # bisher verwendet!
        mu, sigma = 12.3, 1.55**0.5     #Tournadre et al. 2011
        a = 1.52
    
        ib_elems_scaled = pd.DataFrame()
        ib_elems_loc = pd.DataFrame()
    
        points = []
        height = [] #height=depth*8/7=length*8/7*2/3=length*16/21
        
        # Track iceberg counts per basin for summary
        iceberg_counts_per_basin = {}
    
        with tqdm(total=len(self.df_agg), file=sys.stdout, desc='go through basins') as pbar:
            for basin_idx in self.df_agg.index:
                # inner loop to enable redo if generation of icebergs takes too long
                # https://stackoverflow.com/questions/36573486/redo-for-loop-iteration-in-python
                while True:
                    b = self.df_agg.loc[basin_idx]
                    print("*****************************")
                    print("*** BASIN = ", basin_idx)

                    # create icebergs for basin [m3]
                    ib_tmp= self._create_icebergs_within_basin(b, idx=basin_idx)
                    if isinstance(ib_tmp, int):
                        if ib_tmp == -1:
                            continue            # equivalent to redo
                    elif ib_tmp.empty:
                        iceberg_counts_per_basin[basin_idx] = 0
                        break                   # now equivalent to continue
                    ib_elems = self._scale_icebergs(ib_tmp)

                    # make list of fesom elements and it's neighbours
                    felems = list(b.elems)
                    if not self.bcavities:
                        for n in b["neigh."]:
                            felems = felems + list(n)
                    felems = list(set(felems))

                    ##############################################################
                    # exclude coastal nodes (and full cells)
                    elems_to_drop = list(self.full_elems)  # Copy to avoid modifying original
                
                    for felem in felems:
                        nodes = self.elem2d.loc[felem].values
                        coastal = False
                        for node in nodes:
                            lon, lat, tmp = self.nod2d.loc[node]
                            if (tmp == 1 or coastal == 1):
                                coastal = True
                    
                        if coastal == 1:
                            elems_to_drop.append(felem)

                    #print(" * drop these element indices: ", elems_to_drop) 
                    new_felems = [elem for elem in felems if elem not in elems_to_drop]
                    felems = new_felems
                    ##############################################################

                    if len(felems) != 0:
                        # Shuffle elements to spread icebergs across calving front
                        np.random.shuffle(felems)
                        tmp = felems * int(len(ib_elems) / len(felems)) + felems[:len(ib_elems)%len(felems)]
                        felems = tmp
                        
                        # Track iceberg count for this basin
                        iceberg_counts_per_basin[basin_idx] = len(ib_elems)

                        with tqdm(total=len(self.df_agg), file=sys.stdout, desc='initialize icebergs') as pbar:
                            for felem, ib_idx in zip(felems, ib_elems.index):
                                ib_elem = ib_elems.loc[ib_idx]
                                
                                nod1, nod2, nod3 = self.elem2d.loc[felem].values
                                lon1, lat1, tmp = self.nod2d.loc[nod1].values
                                lon2, lat2, tmp = self.nod2d.loc[nod2].values
                                lon3, lat3, tmp = self.nod2d.loc[nod3].values

                                r1 = random.rand()
                                r2 = random.rand()
                                
                                lower_bound = 0.25
                                upper_bound = 0.75
    
                                r1 = r1 * (upper_bound - lower_bound) + lower_bound
                                r2 = r2 * (upper_bound - lower_bound) + lower_bound
                                #https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle
                                try:
                                    lon = (1-np.sqrt(r1))*lon1 + (np.sqrt(r1)*(1-r2))*lon2 + (r2*np.sqrt(r1))*lon3
                                    lat = (1-np.sqrt(r1))*lat1 + (np.sqrt(r1)*(1-r2))*lat2 + (r2*np.sqrt(r1))*lat3
                                except:
                                    continue
                            
                                # Generate calving day based on scaling factor
                                calving_day = self._generate_calving_day(ib_elem.scaling)
                                
                                if ib_elems_loc.empty:
                                    ib_elems_loc = pd.DataFrame({"length": [ib_elem.length], 
                                                                "depth": [ib_elem.depth],
                                                                "scaling": [ib_elem.scaling],
                                                                "lon": [lon], "lat": [lat],
                                                                "felem": [felem],
                                                                "calving_day": [calving_day]})
                                else:
                                    ib_elems_loc = pd.concat([ib_elems_loc, pd.DataFrame({"length": [ib_elem.length], 
                                                                                        "depth": [ib_elem.depth],
                                                                                        "scaling": [ib_elem.scaling],
                                                                                        "lon": [lon], "lat": [lat],
                                                                                        "felem": [felem],
                                                                                        "calving_day": [calving_day]})])
                                pbar.update(1)
                    pbar.update(1)
                    break

        if not ib_elems_loc.empty:
            with open(os.path.join(self.icb_path, "icb_longitude.dat"), fmode) as f:
                np.savetxt(f, ib_elems_loc.lon.values)
                f.close()
            with open(os.path.join(self.icb_path, "icb_latitude.dat"), fmode) as f:
                np.savetxt(f, ib_elems_loc.lat.values)
                f.close()
            with open(os.path.join(self.icb_path, "icb_length.dat"), fmode) as f:
                np.savetxt(f, ib_elems_loc.length.values * 1e3)
                f.close()
            with open(os.path.join(self.icb_path, "icb_height.dat"), fmode) as f:
                np.savetxt(f, ib_elems_loc.depth.values * 1e3)
                f.close()
            with open(os.path.join(self.icb_path, "icb_scaling.dat"), fmode) as f:
                np.savetxt(f, ib_elems_loc.scaling.values, fmt='%d')
                f.close()
            with open(os.path.join(self.icb_path, "icb_felem.dat"), fmode) as f:
                np.savetxt(f, ib_elems_loc.felem.values + 1, fmt='%d')
                f.close()
            with open(os.path.join(self.icb_path, "icb_calving_day.dat"), fmode) as f:
                np.savetxt(f, ib_elems_loc.calving_day.values, fmt='%d')
                f.close()
            
            # Verify no new iceberg is seeded in an already occupied element
            if hasattr(self, 'full_elems') and len(self.full_elems) > 0:
                new_elems = set(ib_elems_loc.felem.values)
                occupied_elems = set(self.full_elems)
                conflicts = new_elems.intersection(occupied_elems)
                if len(conflicts) == 0:
                    print(f" * VERIFICATION PASSED: No new icebergs seeded in occupied elements")
                    print(f"   (checked {len(new_elems)} new elements against {len(occupied_elems)} occupied elements)")
                else:
                    # Convert to 1-indexed for reporting (FESOM convention)
                    conflicts_1idx = [e + 1 for e in conflicts]
                    print(f" * WARNING: {len(conflicts)} new icebergs seeded in occupied elements!")
                    print(f"   Conflicting elements (1-indexed): {conflicts_1idx[:10]}{'...' if len(conflicts) > 10 else ''}")
            
            # Print final summary
            self._print_iceberg_summary(iceberg_counts_per_basin)
    
    def _print_iceberg_summary(self, iceberg_counts_per_basin):
        """
        Print a summary of iceberg generation results.
        """
        print("\n")
        print("=" * 80)
        print("                        ICEBERG GENERATION SUMMARY")
        print("=" * 80)
        
        # Read generated files to verify
        try:
            lon_file = os.path.join(self.icb_path, "icb_longitude.dat")
            total_from_files = len(np.loadtxt(lon_file))
        except:
            total_from_files = "N/A"
        
        # Read calving days and scaling for statistics
        try:
            calving_days = np.loadtxt(os.path.join(self.icb_path, "icb_calving_day.dat")).astype(int)
            scaling = np.loadtxt(os.path.join(self.icb_path, "icb_scaling.dat")).astype(int)
            has_calving_stats = True
        except:
            has_calving_stats = False
        
        # Total flux
        total_flux_km3 = self.total_calving_flux.values * 1e-9
        print(f"\n  TOTAL ICEBERG FLUX:")
        print(f"    From calving calculation:  {total_flux_km3:.2f} km³/year (ice)")
        print(f"    Total icebergs generated:  {sum(iceberg_counts_per_basin.values())}")
        print(f"    Total icebergs in files:   {total_from_files}")
        
        # Per-basin summary
        print(f"\n  PER-BASIN BREAKDOWN:")
        print(f"  {'-' * 76}")
        print(f"  {'Basin':>6} | {'Weight':>8} | {'Flux (km³/y)':>14} | {'Icebergs':>10} | {'CF Elements':>12}")
        print(f"  {'-' * 76}")
        
        for i, basin_id in enumerate(self.basin_ids):
            weight = self.basin_weights[i]
            flux = self.basin_discharge[basin_id] * 1e-9
            n_icebergs = iceberg_counts_per_basin.get(basin_id, 0)
            n_elems = len(self.df_agg.loc[basin_id, "elems"])
            print(f"  {basin_id:>6} | {weight:>8.4f} | {flux:>14.2f} | {n_icebergs:>10} | {n_elems:>12}")
        
        print(f"  {'-' * 76}")
        print(f"  {'TOTAL':>6} | {sum(self.basin_weights):>8.4f} | {total_flux_km3:>14.2f} | {sum(iceberg_counts_per_basin.values()):>10} | {sum(len(self.df_agg.loc[bid, 'elems']) for bid in self.basin_ids):>12}")
        print(f"  {'-' * 76}")
        
        # Calving day statistics
        if has_calving_stats and len(calving_days) > 0:
            print(f"\n  CALVING DAY DISTRIBUTION:")
            print(f"  {'-' * 76}")
            
            # Define austral summer as Nov-Feb (days 305-365 and 1-59)
            is_summer = ((calving_days >= 305) | (calving_days <= 59))
            
            # Large icebergs (scaling == 1)
            large_mask = scaling == 1
            if large_mask.sum() > 0:
                large_summer_pct = 100 * is_summer[large_mask].sum() / large_mask.sum()
                print(f"    Large icebergs (scaling=1):   {large_mask.sum():>6} total, {large_summer_pct:>5.1f}% in austral summer (Nov-Feb)")
            
            # Small icebergs (scaling > 1)
            small_mask = scaling > 1
            if small_mask.sum() > 0:
                small_summer_pct = 100 * is_summer[small_mask].sum() / small_mask.sum()
                print(f"    Small icebergs (scaling>1):   {small_mask.sum():>6} total, {small_summer_pct:>5.1f}% in austral summer (Nov-Feb)")
            
            print(f"  {'-' * 76}")
        
        print(f"\n  Output directory: {self.icb_path}")
        print("=" * 80)
        print("\n")


class point:
        def __init__(self, x, y):
            self.x = x
            self.y = y
        def to_dict(self):
            return {'x': self.x, 'y': self.y}

# https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle
def sign(p1, p2, p3):
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)

def PointInTriangle(pt, v1, v2, v3):
    d1 = sign(pt, v1, v2)
    d2 = sign(pt, v2, v3)
    d3 = sign(pt, v3, v1)

    has_neg = (d1 < 0) or (d2 < 0) or (d3 < 0)
    has_pos = (d1 > 0) or (d2 > 0) or (d3 > 0)

    return not (has_neg and has_pos)

def PointTriangle_distance(lon0, lat0, lon1, lat1, lon2, lat2, lon3, lat3):
    d1 = ne.evaluate('(lon1 - lon0)**2 + (lat1 - lat0)**2')
    d2 = ne.evaluate('(lon2 - lon0)**2 + (lat2 - lat0)**2')
    d3 = ne.evaluate('(lon3 - lon0)**2 + (lat3 - lat0)**2')
    
    dis = d1+d2+d3
    ind = np.where(dis == np.amin(dis))

    p1 = point(lon1[ind], lat1[ind])
    p2 = point(lon2[ind], lat2[ind])
    p3 = point(lon3[ind], lat3[ind])

    return [p1, p2, p3], ind[0][0]

def seconds_per_month(years):
    """
    Return the number of seconds in each month for given year(s).

    Parameters
    ----------
    years : int or iterable of int
        Year or list of years (e.g., 1980 or [1980, 1981, 1982])

    Returns
    -------
    numpy array
        Array of seconds in each month
    """
    
    # Accept single year or iterable
    if isinstance(years, int):
        years = [years]
    
    secs = []
    
    for y in years:
        for m in range(1, 13):
            # start of month
            start = datetime(y, m, 1)
            # start of next month (handle December→January rollover)
            if m == 12:
                end = datetime(y + 1, 1, 1)
            else:
                end = datetime(y, m + 1, 1)
            secs.append((end - start).total_seconds())
    
    return np.array(secs)










