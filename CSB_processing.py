# -*- coding: utf-8 -*-
"""
Created on Tue May 21 12:49:05 2024
Updated on Fri Jun 26 2025

@author: Anthony.R.Klemm

This script now supports two options for reference bathymetry:
  - Automated BlueTopo download.
  - User-provided local BAG or GeoTIFF file.

This script integrates multiple stages of CSB data processing:
1.  Initial processing and tide correction.
2.  (Optional) Post-processing analysis including histograms, offset application,
    and outlier detection for individual transits.
3.  (Optional) A final gridding stage to produce averaged GeoTIFFs and
    filtered GeoPackages, either as a single file or tiled using a
    tessellation scheme.
4.  (Optional) Organization of final rasters by CRS and VRT creation.

"""

import tkinter as tk
from tkinter import filedialog, ttk
import geopandas as gpd
import pandas as pd
import numpy as np
import requests
import os
from osgeo import gdal, osr
import rasterio
from scipy.interpolate import interp1d
from rasterio.features import shapes, rasterize
from shapely.geometry import shape, LineString, Point
from shapely.validation import make_valid
from rasterio.transform import from_origin
from rasterio.warp import calculate_default_transform, reproject, Resampling
from datetime import datetime, timedelta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import threading
from rasterio.merge import merge
from rasterio.enums import Resampling as rioResampling
from shapely.ops import unary_union
from skimage.morphology import binary_dilation, binary_erosion
import shutil
import glob
import subprocess
import duckdb
import gc
import math

# --- IMPORTS for post-processing analysis ---
import seaborn as sns
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from scipy.ndimage import uniform_filter1d


# Global variables
title = ""
csb = ""
BAG_filepath = ""
fp_zones = ""
output_dir = ""
resolution = 20  # resolution of quick-look geotiff raster
MASTER_OFFSET_FILE = os.path.join(output_dir, "master_offsets.csv")

pd.set_option('display.max_columns', None)


def loadCSB():
    print('*****Reading CSB input csv file*****')
    
    df = pd.read_csv(csb)
    
    # Convert 'depth' to numeric and 'time' to datetime; invalid parsing results in NaN
    df['depth'] = pd.to_numeric(df['depth'], errors='coerce')
    df['time'] = pd.to_datetime(df['time'], errors='coerce')
    
    # Drop rows where 'time' or 'depth' conversion failed
    df = df.dropna(subset=['time', 'depth'])
    
    # Filter rows based on depth criteria
    df = df[(df['depth'] > 0.5) & (df['depth'] < 1000)]
    
    # Define time bounds and filter by time
    lower_bound = pd.to_datetime("2014")
    upper_bound = pd.to_datetime(str(datetime.now().year + 1))
    df = df[(df['time'] > lower_bound) & (df['time'] < upper_bound)]
    
    # Drop duplicate rows based on key columns
    df = df.drop_duplicates(subset=['lon', 'lat', 'depth', 'time', 'unique_id'])
    
    # Create a GeoDataFrame using lon/lat columns
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs="EPSG:4326")
    return gdf

def insert_into_duckdb(gdf, duckdb_path):
    """
    Inserts (or appends) data, automatically upgrading the table schema if necessary.
    -- FINAL, ROBUST VERSION --
    """
    # Ensure the directory for the DuckDB file exists.
    duckdb_dir = os.path.dirname(duckdb_path)
    if not os.path.exists(duckdb_dir):
        os.makedirs(duckdb_dir, exist_ok=True)

    df = gdf.copy()
    if 'geometry' in df.columns:
        df['wkb_geom'] = df['geometry'].apply(lambda geom: geom.wkb if geom is not None else None)
        df = df.drop(columns='geometry')
    else:
        df['wkb_geom'] = None

    # This master list defines the complete, final schema.
    master_columns_with_types = {
        "ControlStn": "VARCHAR", "Raster_Value": "DOUBLE", "Uncertainty_Value": "DOUBLE",
        "accuracy_score": "DOUBLE", "date_range": "VARCHAR", "depth_new": "DOUBLE",
        "depth_old": "DOUBLE", "depthfinal": "DOUBLE", "lat": "FLOAT", "lon": "FLOAT",
        "offset_value": "DOUBLE", "platform_name_x": "VARCHAR", "provider": "VARCHAR",
        "std_dev": "DOUBLE", "tile_name": "VARCHAR", "time": "VARCHAR", "unique_id": "VARCHAR",
        "wkb_geom": "BLOB", "diff": "DOUBLE", "depth_mod": "DOUBLE", "uncertainty_vert": "DOUBLE",
        "uncertainty_hori": "DOUBLE", "Outlier": "BOOLEAN", "transit_id": "VARCHAR",
        "vessel_speed_smoothed": "DOUBLE"
    }
    master_columns = list(master_columns_with_types.keys())

    # Ensure the incoming DataFrame has all columns, filling missing with None
    for col in master_columns:
        if col not in df.columns:
            df[col] = None
    df = df[master_columns] # Ensure consistent order

    try:
        with duckdb.connect(database=duckdb_path, read_only=False) as con:
            # Check if the 'csb' table exists
            table_exists = con.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='csb'").fetchone()

            if table_exists:
                # --- SCHEMA MIGRATION LOGIC ---
                # Table exists, so check and add any missing columns.
                existing_columns_df = con.execute("DESCRIBE csb").fetchdf()
                existing_columns = existing_columns_df['column_name'].tolist()
                
                for col_name, col_type in master_columns_with_types.items():
                    if col_name not in existing_columns:
                        print(f"Schema mismatch detected. Adding missing column '{col_name}' to the 'csb' table.")
                        con.execute(f"ALTER TABLE csb ADD COLUMN {col_name} {col_type};")
            else:
                # Table does not exist, create it with the full schema.
                columns_for_create = ", ".join([f"{name} {dtype}" for name, dtype in master_columns_with_types.items()])
                con.execute(f"CREATE TABLE csb ({columns_for_create})")

            # Now, the table schema is guaranteed to match our DataFrame.
            # Use an explicit column list in the INSERT statement for maximum safety.
            col_list_for_insert = ", ".join(master_columns)
            con.register("temp_df", df)
            con.execute(f"INSERT INTO csb ({col_list_for_insert}) SELECT * FROM temp_df")
            
        print(f"Data successfully appended to DuckDB table at {duckdb_path}.")
    except Exception as e:
        print(f"CRITICAL ERROR inserting into DuckDB: {e}")
        raise e

            
def draft_corr(master_offsets_df): # MODIFIED: Accept DataFrame
    csb_corr = derive_draft(master_offsets_df) # MODIFIED: Pass DataFrame

    # Merge the CSB data with the master offsets based on unique vessel ID
    # This will now include any newly derived offsets from the step above
    master_offsets_updated = read_master_offsets() # Read the potentially updated file
    csb_corr1 = csb_corr.merge(master_offsets_updated, on='unique_id', how='left')

    # Apply the offset correction
    # Fill missing offsets with 0 so the calculation doesn't fail
    csb_corr1['offset_value'] = csb_corr1['offset_value'].fillna(0)
    csb_corr1['depthfinal'] = csb_corr1['depth_new'] - csb_corr1['offset_value']
    csb_corr1['depthfinal'] = csb_corr1['depthfinal'] * -1  

    # try to drop some unneeded columns
    csb_corr1 = csb_corr1.drop(columns=['s', 'f', 'q', 'DataProv', 'ControlS_2', 'ControlS_1', 'row_id', 'platform_name_y'], errors='ignore')
    
    print('*****Processed CSB data ready*****')
    # Instead of immediately exporting to geopackage, return the processed GeoDataFrame.
    return csb_corr1

def rasterize_CSB(master_offsets_df): # MODIFIED: Accept DataFrame
    csb_corr1 = draft_corr(master_offsets_df) # MODIFIED: Pass DataFrame
    
    #Based on GUI options, insert processed data into DuckDB.
    if duckdb_option_var.get():
        duckdb_path = os.path.join(output_dir, "csb.duckdb")
        insert_into_duckdb(csb_corr1, duckdb_path)
    
    # Optionally export as geopackage if the checkbox is selected.
    if export_gp_var.get():
        gpkg_path = os.path.join(output_dir, 'csb_processed_'+ title +'.gpkg')
        print('*****Exporting processed CSB data to geopackage*****')
        csb_corr1.to_file(gpkg_path, driver='GPKG', layer='csb')
        print(f"Geopackage exported to {gpkg_path}")
    
    return csb_corr1

def reproject_tiff(input_dir, output_dir, target_epsg='EPSG:3395'):
    """
    Searches the input directory for files ending in .tif or .tiff,
    reprojects each to the target EPSG code, and writes the output to output_dir.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.lower().endswith((".tif", ".tiff")):
                input_path = os.path.join(root, file)
                relative_path = os.path.relpath(root, input_dir)
                base_filename, file_extension = os.path.splitext(file)
                new_filename = f"{base_filename}_reprojected{file_extension}"
                output_path = os.path.join(output_dir, relative_path, new_filename)
                output_folder = os.path.dirname(output_path)
                if not os.path.exists(output_folder):
                    os.makedirs(output_folder)
                command = [
                    'gdalwarp', '-t_srs', target_epsg,
                    input_path, output_path
                ]
                try:
                    subprocess.run(command, check=True)
                    # print(f"Reprojected {input_path} to {output_path}")
                except Exception as e:
                    print(f"[DEBUG] Error reprojecting {input_path}: {e}")

def mosaic_tiles(tiles_dir):
    """
    Searches the given tiles directory for reprojected TIFF files (with the "_reprojected" suffix).
    Then it merges the available TIFF files using rasterio.merge.merge, 
    writes the merged raster to a file, and returns the path of the merged raster.
    """
    mosaic_raster_path = os.path.join(tiles_dir, 'merged_tiles.tif')
    
    # Call the reproject_tiff function (make sure it also works with both .tif and .tiff files)
    reproject_tiff(tiles_dir, tiles_dir)
    
    # List to store the opened rasters
    raster_list = []
    print("[DEBUG] Searching for reprojected TIFF files in:", tiles_dir)
    for root, dirs, files in os.walk(tiles_dir):
        for file in files:
            # Use lower-case for a case-insensitive check
            if file.lower().endswith('_reprojected.tif') or file.lower().endswith('_reprojected.tiff'):
                raster_path = os.path.join(root, file)
                print(f"[DEBUG] Found reprojected file: {raster_path}")
                try:
                    src = rasterio.open(raster_path)
                    raster_list.append(src)
                except Exception as e:
                    print(f"[DEBUG] Error opening {raster_path}: {e}")
    
    if not raster_list:
        raise RuntimeError("No input dataset specified. No reprojected TIFF files were found in " + tiles_dir)
    
    print(f"[DEBUG] Total input datasets found: {len(raster_list)}")
    try:
        merged_raster, out_transform = merge(raster_list)
    except Exception as e:
        raise RuntimeError(f"Error during merging: {e}")
    
    # Use metadata from the last opened raster in the list
    out_meta = raster_list[-1].meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": merged_raster.shape[1],
        "width": merged_raster.shape[2],
        "transform": out_transform,
        "count": raster_list[-1].count
    })
    
    try:
        with rasterio.open(mosaic_raster_path, "w", **out_meta) as dest:
            for i in range(1, out_meta["count"] + 1):
                dest.write(merged_raster[i - 1, :, :], i)
        print(f"[DEBUG] Merged raster written to {mosaic_raster_path}")
    except Exception as e:
        raise RuntimeError(f"Error writing merged raster: {e}")
    
    # Close all opened raster files
    for src in raster_list:
        src.close()
    
    return mosaic_raster_path


def create_convex_hull_and_download_tiles(csb_data_path, output_dir, use_bluetopo=True):
    """
    Loads CSB data, builds a convex hull, and writes it to a shapefile.
    If use_bluetopo is True, it creates a dedicated Modeling folder,
    downloads BlueTopo tiles, copies them to a separate archive folder,
    and then builds a VRT from all GeoTIFF files found recursively in that folder.
    If use_bluetopo is False, it returns the user‐provided BAG_filepath.
    """

    # Load CSB data and create convex hull
    csb_data = pd.read_csv(csb_data_path)
    gdf = gpd.GeoDataFrame(csb_data, geometry=gpd.points_from_xy(csb_data.lon, csb_data.lat))
    gdf = gdf.set_crs(4326, allow_override=True)
    convex_hull_polygon = gdf.unary_union.convex_hull
    convex_hull_gdf = gpd.GeoDataFrame(geometry=[convex_hull_polygon], crs=gdf.crs)
    convex_hull_shapefile = os.path.join(output_dir, "convex_hull_polygon.shp")
    convex_hull_gdf.to_file(convex_hull_shapefile)
    print(f"Convex hull shapefile written to: {convex_hull_shapefile}")

    if use_bluetopo:
        # Create the 'Modeling' folder for downloading tiles
        bluetopo_tiles_dir = os.path.join(output_dir, "Modeling")
        os.makedirs(bluetopo_tiles_dir, exist_ok=True)

        # Download BlueTopo tiles
        from nbs.bluetopo import fetch_tiles
        fetch_tiles(bluetopo_tiles_dir, convex_hull_shapefile, data_source='modeling')

        # Use the CSV file's title (assumed to be stored in the global variable 'title')
        bluetopo_tiles_copy = os.path.join(output_dir, f"BlueTopo_Tiles_{title}")
        if os.path.exists(bluetopo_tiles_copy):
            shutil.rmtree(bluetopo_tiles_copy)
        shutil.copytree(bluetopo_tiles_dir, bluetopo_tiles_copy)
        print(f"Copied BlueTopo tiles to {bluetopo_tiles_copy}")
        
        # Build a VRT from the copied tiles using a unique folder name that includes the title.
        vrt_dir = os.path.join(output_dir, f"BlueTopo_VRT_{title}")
        os.makedirs(vrt_dir, exist_ok=True)
        
        # Create glob patterns to find both .tif and .tiff files in the unique tiles folder.
        tif_pattern = os.path.join(bluetopo_tiles_copy, '**', '*.tif')
        tiff_pattern = os.path.join(bluetopo_tiles_copy, '**', '*.tiff')
        print(f"[DEBUG] Glob pattern for .tif: {tif_pattern}")
        print(f"[DEBUG] Glob pattern for .tiff: {tiff_pattern}")
        tile_files = glob.glob(tif_pattern, recursive=True) + glob.glob(tiff_pattern, recursive=True)
        print("[DEBUG] Found the following tile files for VRT building:")
        for f in tile_files:
            print("  ", f)
        
        if not tile_files:
            raise RuntimeError("No BlueTopo GeoTIFF files were found in " + bluetopo_tiles_copy)
        
        # Save the VRT file with the title appended to the file name.
        vrt_path = os.path.join(vrt_dir, f"merged_tiles_{title}.vrt")
        vrt = gdal.BuildVRT(vrt_path, tile_files)
        vrt = None  # Close the VRT dataset
        print(f"Created VRT at {vrt_path}")
        return vrt_path
    else:
        # If not using automated download, assume BAG_filepath is provided by the user.
        return BAG_filepath

    
def reproject_raster(input_raster, output_raster, dst_crs):
    with rasterio.open(input_raster) as src:
        transform, width, height = calculate_default_transform(
            src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height,
            'driver': 'GTiff',
            'count': 1,
            'dtype': src.dtypes[0]
        })

        with rasterio.open(output_raster, 'w', **kwargs) as dst:
            reproject(
                source=rasterio.band(src, 1),  # Reproject only the first band
                destination=rasterio.band(dst, 1),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=transform,
                dst_crs=dst_crs,
                resampling=Resampling.nearest)


def fetch_tide_data(station_id, start_date, end_date, product, interval=None, attempt_great_lakes=False):
    base_url = "https://api.tidesandcurrents.noaa.gov/api/prod/datagetter"
    params = {
        "begin_date": start_date,
        "end_date": end_date,
        "station": station_id,
        "datum": "MLLW" if not attempt_great_lakes else "LWD",
        "time_zone": "gmt",
        "units": "metric",
        "format": "json",
        "product": product
    }

    if interval:
        params["interval"] = interval

    request_url = requests.Request('GET', base_url, params=params).prepare().url
    print(f"Requesting URL: {request_url}")

    response = requests.get(base_url, params=params)
    data = response.json()

    if 'predictions' in data:
        df = pd.json_normalize(data['predictions'])
        data_type = "predicted data"
    elif 'data' in data:
        df = pd.json_normalize(data['data'])
        data_type = "observed data"
    else:
        print(f"No data returned for URL: {request_url}")
        return pd.DataFrame()

    df['t'] = pd.to_datetime(df['t'])
    # Convert 'v' to numeric, coercing errors to NaN
    df['v'] = pd.to_numeric(df['v'], errors='coerce')
    if df['v'].isna().any():
        print("Warning: Some tide values could not be converted to numeric and will be dropped.")
        df = df.dropna(subset=['v'])

    print(f"Pulled {data_type} for station {station_id} from {start_date} to {end_date}")
    return df


def check_for_gaps(dataframe, max_gap_duration='1h'):
    gaps = dataframe['t'].diff() > pd.Timedelta(max_gap_duration)
    return gaps.any()

def cosine_interpolation(df, start_date, end_date):
    df['time_num'] = (df['t'] - pd.Timestamp("1970-01-01")) // pd.Timedelta('1s')
    df = df.sort_values('time_num')
    interp_func = interp1d(df['time_num'], df['v'], kind='cubic')
    time_num_grid = np.linspace(df['time_num'].min(), df['time_num'].max(), num=3000)
    v_grid = interp_func(time_num_grid)
    time_grid = pd.to_datetime(time_num_grid, unit='s')
    
    # Trimming
    trim_start = pd.to_datetime(start_date) + pd.Timedelta(hours=12)
    trim_end = pd.to_datetime(end_date) - pd.Timedelta(hours=12)
    trimmed_df = pd.DataFrame({'t': time_grid, 'v': v_grid})
    trimmed_df = trimmed_df[(trimmed_df['t'] >= trim_start) & (trimmed_df['t'] <= trim_end)]
    #print(trimmed_df)
    return trimmed_df

def create_survey_outline(raster_path, output_dir, title, desired_resolution=8, dilation_iterations=3, erosion_iterations=2):
    print("starting create_survey_outline() function")
    with rasterio.open(raster_path) as raster:
        # Resample the raster
        data = raster.read(
            1,  # Reading only the first band
            out_shape=(
                raster.height // desired_resolution,
                raster.width // desired_resolution
            ),
            resampling=Resampling.bilinear
        )

        # Create a binary mask
        nodata = raster.nodatavals[0] or 1000000
        binary_mask = (data != nodata).astype(np.uint8)

        # Apply dilation and erosion
        for _ in range(dilation_iterations):
            binary_mask = binary_dilation(binary_mask)
        for _ in range(erosion_iterations):
            binary_mask = binary_erosion(binary_mask)

        # Ensure binary_mask is of type uint8
        binary_mask = binary_mask.astype(np.uint8)

        # Generate shapes from the binary mask
        transform = raster.transform * raster.transform.scale(
            (raster.width / data.shape[-1]),
            (raster.height / data.shape[-2])
        )

        #polygons = [shape(geom) for geom, val in shapes(binary_mask, mask=binary_mask, transform=transform) if val == 1]

        # Perform unary union
        #unified_geometry = unary_union(polygons)
        # Generate polygons from the binary mask and make them valid
        polygons = [make_valid(shape(geom)) for geom, val in shapes(binary_mask, mask=binary_mask, transform=transform) if val == 1]

        # Simplify polygons to reduce complexity
        simplified_polygons = [polygon.simplify(tolerance=0.001, preserve_topology=True) for polygon in polygons]

        # Perform unary union on simplified, valid polygons
        unified_geometry = unary_union(simplified_polygons)
        # Reproject unified geometry to WGS84 before simplification
        geo_df = gpd.GeoDataFrame(geometry=[unified_geometry], crs=raster.crs)
        geo_df = geo_df.to_crs(epsg=4326)

        # Simplify the geometry
        geo_df['geometry'] = geo_df.geometry.simplify(tolerance=0.001)

        # Check and fix bad topology if necessary
        geo_df['geometry'] = geo_df.geometry.apply(lambda geom: geom.buffer(0) if not geom.is_valid else geom)

        # Save to a shapefile
        bathy_polygon_shp = f"{output_dir}/{title}_bathy_polygon.shp"
        geo_df.to_file(bathy_polygon_shp, driver='ESRI Shapefile')

        print('Bathymetry polygon shapefile created.')
        return bathy_polygon_shp

def tides():
    gdf = loadCSB()
    print('CSB data from csv file loaded. Starting tide correction')

    zones = gpd.read_file(fp_zones)
    join = gpd.sjoin(gdf, zones, how='inner', predicate='within')
    join = join.astype({'time': 'datetime64[ns]'})
    join = join.sort_values('time')

    def generate_date_ranges(dates):
        dates.sort()
        date_ranges = []
        for date in dates:
            if not date_ranges or date - pd.Timedelta(days=1) > date_ranges[-1][1]:
                date_ranges.append([date, date])
            else:
                date_ranges[-1][1] = date
        date_ranges = [(start_date - pd.Timedelta(days=1), end_date + pd.Timedelta(days=1)) for start_date, end_date in date_ranges]
        return [(start_date.strftime('%Y%m%d'), end_date.strftime('%Y%m%d')) for start_date, end_date in date_ranges]

    tdf = []
    known_subordinate_stations = set()
    known_great_lakes_stations = set()

    for station_id in join['ControlStn'].unique():
        station_dates = join[join['ControlStn'] == station_id]['time'].dt.floor('d').unique()
        date_ranges = generate_date_ranges(list(station_dates))

        for start_date, end_date in date_ranges:
            # Try to fetch observed data first
            verified_data = fetch_tide_data(station_id, start_date, end_date, product="water_level")
            if not verified_data.empty and not check_for_gaps(verified_data):
                tdf.append(verified_data)
            else:
                # If there are gaps, try fetching 6-minute predicted data
                predicted_data = fetch_tide_data(station_id, start_date, end_date, product="predictions")
                if not predicted_data.empty and not check_for_gaps(predicted_data):
                    tdf.append(predicted_data)
                else:
                    # If that doesn't work, fallback to predicted hilo data
                    hilo_predictions = fetch_tide_data(station_id, start_date, end_date, product="predictions", interval='hilo')
                    if not hilo_predictions.empty:
                        interpolated_hilo = cosine_interpolation(hilo_predictions, start_date, end_date)
                        tdf.append(interpolated_hilo)
                        known_subordinate_stations.add(station_id)
                    else:
                        great_lakes_data = fetch_tide_data(station_id, start_date, end_date, product="water_level", attempt_great_lakes=True)
                        if not great_lakes_data.empty:
                            tdf.append(great_lakes_data)
                            known_great_lakes_stations.add(station_id)
                        else:
                            print(f"No water level data available for station {station_id}.")

    if tdf:
        tdf = pd.concat(tdf)
        print("Concatenated tdf shape:", tdf.shape)

        tdf = tdf.sort_values('t')
        jtdf = pd.merge_asof(join, tdf, left_on='time', right_on='t')
        print("jtdf shape before column drop:", jtdf.shape)

        columns_to_drop = ['Shape__Are', 'Shape__Len', 'Input_FID', 'id', 'name', 'state', 'affil', 
                           'latitude', 'longitude', 'data', 'metaapi', 'dataapi', 'Shape_Le_2']
        jtdf.drop(columns=columns_to_drop, inplace=True, errors='ignore')
        print("jtdf shape after column drop:", jtdf.shape)

        jtdf = jtdf.dropna(subset=['depth', 'time', 'geometry'])
        print("jtdf shape after dropna:", jtdf.shape)

        jtdf['t_corr'] = jtdf['t'] + pd.to_timedelta(jtdf['ATCorr'], unit='m')

        newdf = jtdf[['t_corr', 'v']].copy()
        print("newdf shape before dropna:", newdf.shape)

        newdf = newdf.rename(columns={'v': 'v_new', 't_corr': 't_new'})
        newdf = newdf.sort_values('t_new').dropna()
        print("newdf shape after dropna:", newdf.shape)

        csb_corr = pd.merge_asof(jtdf, newdf, left_on='time', right_on='t_new', direction='nearest')
        print("csb_corr shape before dropna:", csb_corr.shape)

        print("csb_corr shape after dropna:", csb_corr.shape)

        csb_corr['depth_new'] = csb_corr['depth'] - (csb_corr['RR'] * csb_corr['v_new'])
        print("csb_corr shape after applying tide corrections:", csb_corr.shape)

        csb_corr = gpd.GeoDataFrame(csb_corr, geometry='geometry', crs='EPSG:4326')
        csb_corr['time'] = csb_corr['time'].dt.strftime("%Y%m%d %H:%M:%S")

        csb_corr = csb_corr[(csb_corr['depth'] > 1.5) & (csb_corr['depth'] < 1000)]
        csb_corr = csb_corr.rename(columns={'depth': 'depth_old'}).drop(columns=['index_right', 'ATCorr', 'RR', 'ATCorr2', 'RR2', 'Shape_Leng', 'Shape_Area', 'Shape_Le_1', 't', 'v', 't_corr', 't_new', 'v_new'])
        return csb_corr
    else:
        print("No tide data available for the specified period.")
        return pd.DataFrame()

def BAGextract():
    print("starting BAGextract() function")
    global BAG_filepath
    BAG_filepath = os.path.abspath(BAG_filepath)  # Ensure it's an absolute path
    print("DEBUG - BAG_filepath:", BAG_filepath)
    print('*****Starting to import BAG bathy and aggregate to 8m geotiff*****')

    # Translate options and creation of multi-band GeoTIFF
    translate_options = gdal.TranslateOptions(bandList=[1, 2], creationOptions=['COMPRESS=LZW'])
    intermediate_raster_path = output_dir + '/' + title + '_intermediate.tif'

    dataset = gdal.Open(BAG_filepath, gdal.GA_ReadOnly)
    if dataset is None:
        print("Error: Unable to open the file. Check the file path.")
        return None, None

    gdal.Translate(intermediate_raster_path, dataset, options=translate_options)

    # Check the CRS
    crs = osr.SpatialReference(wkt=dataset.GetProjection())
    dataset = None  # Close the dataset

    if crs.IsGeographic():
        print("Reprojecting GeoTIFF to a projected CRS...")
        reprojected_raster_path = output_dir + '/' + title + '_reprojected.tif'
        gdal.Warp(reprojected_raster_path, intermediate_raster_path, dstSRS='EPSG:3395', creationOptions=['COMPRESS=LZW'])
        intermediate_raster_path = reprojected_raster_path

    # Open the intermediate raster to check its resolution
    with rasterio.open(intermediate_raster_path) as src:
        res_x, res_y = src.res

    # Check if the resolution is coarser than 8m
    if max(res_x, res_y) > 8:
        output_raster = intermediate_raster_path
    else:
        # Resample to desired resolution
        output_raster_resampled = output_dir + '/' + title + '_5m_MLLW.tif'
        gdal.Warp(output_raster_resampled, intermediate_raster_path, xRes=8, yRes=8, creationOptions=['COMPRESS=LZW'])
        output_raster = output_raster_resampled
        # Clean up intermediate file
        os.remove(intermediate_raster_path)

    # Replace NaN values and update nodata value in the raster
    with rasterio.open(output_raster) as src:
        data = src.read()
        meta = src.meta

    data = np.where(np.isnan(data), 1000000, data)
    meta.update(nodata=1000000)

    with rasterio.open(output_raster, 'w', **meta) as dst:
        dst.write(data)

    # Reproject the raster to WGS84
    output_raster_wgs84 = output_dir + '/' + title + '_wgs84.tif'
    gdal.Warp(output_raster_wgs84, output_raster, dstSRS='EPSG:4326', creationOptions=['COMPRESS=LZW'])

    # Call create_survey_outline to generate the bathymetry polygon shapefile
    bathy_polygon_shp = create_survey_outline(output_raster_wgs84, output_dir, title)

    # Clean up intermediate files
    try:
        os.remove(intermediate_raster_path)
        os.remove(output_raster_resampled)
        os.remove(output_raster_wgs84)
    except Exception:
        pass

    return output_raster_wgs84, bathy_polygon_shp

def read_master_offsets():
    """Reads the master offsets from a CSV file."""
    MASTER_OFFSET_FILE = os.path.join(output_dir, "master_offsets.csv")

    if os.path.exists(MASTER_OFFSET_FILE):
        return pd.read_csv(MASTER_OFFSET_FILE)
    else:
        return pd.DataFrame(columns=['unique_id', 'platform_name', 'offset_value', 'std_dev', 'accuracy_score', 'date_range', 'tile_name'])

def update_master_offsets(unique_id, platform_name, new_offset, std_dev, date_range, tile_name):
    global MASTER_OFFSET_FILE
    MASTER_OFFSET_FILE = os.path.join(output_dir, "master_offsets.csv")
    master_offsets = read_master_offsets()


    accuracy_score = 1 / std_dev if std_dev != 0 else 0

    #print('checking for existing offset by unique_id and platform_name')
    existing_index = master_offsets[(master_offsets['unique_id'] == unique_id)].index

    new_row = pd.DataFrame([{
        'unique_id': unique_id,
        'platform_name': platform_name,
        'offset_value': new_offset,
        'std_dev': std_dev,
        'accuracy_score': accuracy_score,
        'date_range': date_range,
        'tile_name': tile_name
    }])

    # Exclude empty or all-NA entries before concatenation
    new_row = new_row.dropna(how='all')

    if existing_index.empty:
        master_offsets = pd.concat([master_offsets, new_row], ignore_index=True)
    else:
        if master_offsets.loc[existing_index[0], 'accuracy_score'] <= accuracy_score:
            master_offsets.loc[existing_index[0], list(new_row.columns)] = new_row.iloc[0]

    try:
        master_offsets.to_csv(MASTER_OFFSET_FILE, index=False)
        #print(f"Master offsets updated successfully in {MASTER_OFFSET_FILE}.")
    except Exception as e:
        print(f"Failed to update master offsets: {e}")


def get_raster_values_vectorized(coords, raster_path, batch_size=100000):
    """
    Given a list of (x, y) coordinate pairs and a raster file path,
    returns a list of lists, each containing the pixel values from
    all bands at that coordinate.
    
    If the number of coordinates exceeds batch_size, processing is done in batches.
    Pixels with nodata values are replaced with np.nan.
    """
    all_samples = []
    
    with rasterio.open(raster_path) as src:
        nodata = src.nodatavals  # Tuple of nodata values for each band
        
        # Process in batches if necessary to manage memory usage
        if len(coords) > batch_size:
            for i in range(0, len(coords), batch_size):
                batch_coords = coords[i:i+batch_size]
                batch_samples = list(src.sample(batch_coords))
                all_samples.extend(batch_samples)
        else:
            all_samples = list(src.sample(coords))
    
    # Process the samples to replace nodata values with np.nan
    processed_samples = []
    for sample in all_samples:
        processed = []
        for i, val in enumerate(sample):
            if nodata[i] is not None and val == nodata[i]:
                processed.append(np.nan)
            else:
                processed.append(val)
        processed_samples.append(processed)
        
    return processed_samples
        

def derive_draft(master_offsets_df): # MODIFIED: Accept DataFrame
    output_raster, raster_boundary_shp = BAGextract()
    csb_corr = tides()

    # --- NEW: Get a list of vessels that already have an offset ---
    vessels_with_offsets = master_offsets_df['unique_id'].unique().tolist()
    print(f"Found {len(vessels_with_offsets)} vessels with pre-existing offsets.")

    # Read and fix geometries in the raster boundary and CSB data
    raster_boundary = gpd.read_file(raster_boundary_shp)
    raster_boundary['geometry'] = raster_boundary['geometry'].apply(
        lambda geom: geom if geom.is_valid else geom.buffer(0)
    )
    csb_corr['geometry'] = csb_corr['geometry'].apply(
        lambda geom: geom if geom.is_valid else geom.buffer(0)
    )
    csb_corr['row_id'] = range(len(csb_corr))
    
    # Compute the unary union of the raster boundary geometries once
    boundary_union = raster_boundary.geometry.unary_union
    
    # Get the bounding box of the boundary_union: (minx, miny, maxx, maxy)
    bbox = boundary_union.bounds
    
    # Use the spatial index to find candidate points that intersect the bbox
    possible_matches_index = csb_corr.sindex.query(boundary_union, predicate="intersects")
    possible_matches = csb_corr.iloc[possible_matches_index]
    
    # Now apply the precise .within() filter on this subset
    csb_corr_subset = possible_matches[possible_matches.geometry.within(boundary_union)]
    
    # --- MODIFIED: Filter out vessels that already have an offset ---
    csb_for_offset_derivation = csb_corr_subset[~csb_corr_subset['unique_id'].isin(vessels_with_offsets)]
    
    if csb_for_offset_derivation.empty:
        print("All vessels within reference data already have offsets. Skipping new offset calculation.")
        # We still need to sample the raster for all points for post-processing diff calculation
        # So we continue, but the offset derivation part will be skipped.
    else:
        print(f"Found {csb_for_offset_derivation['unique_id'].nunique()} new vessels requiring offset derivation.")

    # Instead of sampling 1000 points per vessel, use all points.
    sampled_csb_corr_subset = csb_corr_subset.copy()
    
    # Compute date ranges for each vessel (unique_id) using all available points.
    date_ranges = {}
    for name, group in csb_corr_subset.groupby('unique_id'):
        # Ensure time is in datetime format.
        group['time'] = pd.to_datetime(group['time'])
        min_timestamp = group['time'].min()
        max_timestamp = group['time'].max()
        if pd.notnull(min_timestamp) and pd.notnull(max_timestamp):
            date_ranges[name] = (min_timestamp.strftime('%Y%m%d'), max_timestamp.strftime('%Y%m%d'))
        else:
            date_ranges[name] = ("19700101", "19700101")

    # Use the new vectorized raster sampling function over all points.
    try:
        # Extract (x, y) coordinates from the GeoDataFrame.
        coords = [(geom.x, geom.y) for geom in sampled_csb_corr_subset.geometry]
        # Get raster values (supports batching if needed).
        raster_samples = get_raster_values_vectorized(coords, output_raster)
        # Assuming the raster has two bands: assign first as Raster_Value and second as Uncertainty_Value.
        sampled_csb_corr_subset['Raster_Value'] = [vals[0] for vals in raster_samples]
        sampled_csb_corr_subset['Uncertainty_Value'] = [
            vals[1] if len(vals) > 1 else np.nan for vals in raster_samples
        ]
    except KeyError as e:
        print(f"KeyError encountered during selection of Raster_Value from reference bathy: {e}")
    except Exception as e:
        print(f"Unexpected error encountered during selection of Raster_Value from reference bathy: {e}")

    # Merge the new raster values back into the main CSB dataframe.
    try:
        csb_corr = csb_corr.merge(
            sampled_csb_corr_subset[['row_id', 'Raster_Value', 'Uncertainty_Value']],
            on='row_id', how='left'
        )
    except KeyError as e:
        print(f"KeyError encountered during merging: {e}")
    except Exception as e:
        print(f"Unexpected error encountered during merging: {e}")

    # --- MODIFIED: Use the filtered DataFrame for offset derivation ---
    # Filter based on uncertainty threshold.
    filtered_csb_corr = csb_for_offset_derivation.merge(
        csb_corr[['row_id', 'Raster_Value', 'Uncertainty_Value']], on='row_id', how='left'
    )
    filtered_csb_corr = filtered_csb_corr[filtered_csb_corr['Uncertainty_Value'] < 4]

    # Calculate the tide correction difference.
    try:
        if 'Raster_Value' in filtered_csb_corr.columns and 'depth_new' in filtered_csb_corr.columns:
            filtered_csb_corr['diff'] = filtered_csb_corr['depth_new'] - (filtered_csb_corr['Raster_Value'] * -1)
    except KeyError as e:
        print(f"KeyError encountered calculating diff: {e}")
    except Exception as e:
        print(f"Unexpected error encountered calculating diff: {e}")

    # Aggregate differences by vessel (unique_id).
    if not filtered_csb_corr.empty and 'diff' in filtered_csb_corr.columns:
        try:
            out = filtered_csb_corr.groupby('unique_id')['diff'].agg(['mean', 'std', 'count']).reset_index()
            
            out.loc[:, 'mean'] = out['mean'].fillna(0)
            out.loc[:, 'std'] = out['std'].fillna(999)
            out.loc[:, 'count'] = out['count'].fillna(0)
            out.loc[(out['mean'] > 3) | (out['mean'] < -11), ['mean', 'std', 'count']] = [0, 999, 0]
            out.loc[(out['std'] > 7), ['mean', 'std', 'count']] = [0, 999, 0]
            out.to_csv(output_dir + '/VESSEL_OFFSETS_csb_corr_' + title + '.csv', mode='a')

            platform_mapping = filtered_csb_corr[['unique_id', 'platform_name']].drop_duplicates()
            out_with_platform = out.merge(platform_mapping, on='unique_id', how='left')

            for index, row in out_with_platform.iterrows():
                unique_id = row['unique_id']
                platform_name = row['platform_name']
                new_offset = row['mean']
                std_dev = row['std']
                date_range = date_ranges.get(unique_id, ("19700101", "19700101"))
                update_master_offsets(unique_id, platform_name, new_offset, std_dev, date_range, title)
        except KeyError as e:
            print(f"KeyError encountered creating aggregation dataframe: {e}")
        except Exception as e:
            print(f"Unexpected error encountered creating aggregation dataframe: {e}")
    else:
        print("No new offsets to calculate in this run.")

    return csb_corr


def reproject_to_mercator(geodataframe):
    # Reproject to EPSG:3395 (World Mercator)
    return geodataframe.to_crs(epsg=3395)

def rasterize_with_rasterio(geodataframe, output_path, resolution=50, nodatavalue=1000000):
    # Reproject the geodataframe
    geodataframe = reproject_to_mercator(geodataframe)

    # Debugging steps: Check if the GeoDataFrame is empty
    if geodataframe.empty:
        print("Error: The GeoDataFrame is empty.")
        return
    else:
        print(f"Number of features in the GeoDataFrame: {len(geodataframe)}")
        print(geodataframe.head())

    # Define the bounds of your raster
    minx, miny, maxx, maxy = geodataframe.total_bounds
    print("Spatial extent:", minx, miny, maxx, maxy)

    # Check if bounds are reasonable
    if (maxx-minx) <= 0 or (maxy-miny) <= 0:
        print("Error: Invalid spatial extent.")
        return

    # Calculate the dimensions of the raster
    x_res = int((maxx - minx) / resolution)
    y_res = int((maxy - miny) / resolution)

    # Debugging step: Check if resolution is leading to zero dimensions
    if x_res <= 0 or y_res <= 0:
        print("Error: Resolution too high or invalid spatial extent, leading to zero dimensions.")
        return

    # Define the transform
    transform = from_origin(minx, maxy, resolution, resolution)

    # Define the output raster dimensions and CRS
    out_meta = {
        'driver': 'GTiff',
        'height': y_res,
        'width': x_res,
        'count': 1,
        'dtype': 'float32',
        'crs': geodataframe.crs.to_string(),
        'transform': transform,
        'nodata': nodatavalue
    }

    # Rasterize the geometries
    with rasterio.open(output_path, 'w', **out_meta) as out_raster:
        out_raster.write(rasterize(
            ((geom, value) for geom, value in zip(geodataframe.geometry, geodataframe['depthfinal'])),
            out_shape=(y_res, x_res),
            transform=transform,
            fill=nodatavalue
        ), 1)

def get_utm_zone_nad83(lat, lon):
    """
    Calculates the NAD83 UTM zone EPSG code for a given latitude and longitude.
    """
    zone_number = int((lon + 180) // 6) + 1
    if lat < 0:
        raise ValueError("NAD83 UTM zones are generally for northern hemisphere data only.")
    return 26900 + zone_number

def plot_vessel_tracks(gdf):
    # Check which platform_name variant exists and use it
    platform_col = None
    for potential_name in ['platform_name', 'platform_name_x', 'platform_name_y']:
        if potential_name in gdf.columns:
            platform_col = potential_name
            break
    if platform_col is None:
        raise ValueError("No platform_name column found in DataFrame.")
    
    # Ensure 'gdf' is in WGS84
    gdf = gdf.to_crs(epsg=4326)
    
    # Ensure 'time' column is in datetime format
    gdf['time'] = pd.to_datetime(gdf['time'])

    # Initialize an empty dictionary for LineStrings with platform names as keys
    lines_by_platform = {}

    # Group by the platform name column found
    for platform_name, platform_group in gdf.groupby(platform_col):
        # Further group each platform's data by 'unique_id'
        for unique_id, group in platform_group.groupby('unique_id'):
            group = group.sort_values('time')
            current_line = []

            for i, row in group.iterrows():
                if len(current_line) > 0:
                    time_diff = row['time'] - current_line[-1][1]
                    if time_diff.total_seconds() > 600:  # 10 minutes
                        if len(current_line) > 1:
                            line = LineString([point for point, _ in current_line])
                            lines_by_platform.setdefault(platform_name, []).append(line)
                        current_line = []
                current_line.append((row.geometry, row['time']))
            
            if len(current_line) > 1:
                line = LineString([point for point, _ in current_line])
                lines_by_platform.setdefault(platform_name, []).append(line)

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot each platform's lines with a unique color and add to the legend
    for platform_name, lines in lines_by_platform.items():
        for line in lines:
            ax.plot(*line.xy, label=platform_name)
    
    # Handling duplicate labels in legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='center left', bbox_to_anchor=(1, 0.5))

    # Adjust the visualization
    ax.set_title(title + " CSB Vessel Tracks")
    plt.tight_layout()
    
    # Save the plot
    plt.savefig(output_dir + '/csb_plot_' + title + '.png', dpi=300)
    plt.close()

# --- START: POST-PROCESSING ANALYSIS FUNCTIONS ---

def run_histograms_calibration_points(db_path, hist_export_dir):
    print("Starting Post-Processing Step 1: Histograms and Calibration Points")
    os.makedirs(hist_export_dir, exist_ok=True)
    
    with duckdb.connect(database=db_path, read_only=False) as con:
        columns_query = "DESCRIBE csb"
        columns_df = con.execute(columns_query).fetchdf()
        if 'diff' in columns_df['column_name'].values:
            print("Dropping incorrect 'diff' column in DuckDB...")
            con.execute("ALTER TABLE csb DROP COLUMN diff")

        columns_df = con.execute(columns_query).fetchdf()
        if 'diff' not in columns_df['column_name'].values:
            print("Creating 'diff' column in DuckDB with correct calculation...")
            con.execute("ALTER TABLE csb ADD COLUMN diff DOUBLE DEFAULT NULL")
            con.execute("UPDATE csb SET diff = (depth_new *-1 - Raster_Value)*-1 WHERE diff IS NULL")

        unique_ids_query = "SELECT DISTINCT unique_id FROM csb"
        unique_ids = con.execute(unique_ids_query).fetchdf()['unique_id']

        for unique_id in unique_ids:
            print(f'Calculating offset histogram for {unique_id}')
            data_query = f"""
            SELECT unique_id, platform_name_x, diff, lon, lat
            FROM csb
            WHERE unique_id = '{unique_id}' 
              AND diff > -12 AND diff < 12 
              AND Raster_Value > -20 
              AND Uncertainty_Value < 3
            """
            data_df = con.execute(data_query).fetchdf()

            if not data_df.empty:
                platform_name = data_df['platform_name_x'].iloc[0]
                output_csv_path = os.path.join(hist_export_dir, f"{unique_id}_csb_offset_analysis.csv")
                data_df.to_csv(output_csv_path)

                plt.figure(figsize=(10, 6))
                sns.histplot(data_df['diff'], bins=30, kde=True, color="skyblue", label='Histogram')
                plt.axvline(data_df['diff'].mean(), color='green', linestyle='--', label=f'Mean: {data_df["diff"].mean():.2f}')
                plt.axvline(data_df['diff'].mean() - data_df['diff'].std(), color='purple', linestyle='--', 
                            label=f'-1 Std Dev: {(data_df["diff"].mean() - data_df["diff"].std()):.2f}')
                plt.axvline(data_df['diff'].mean() + data_df['diff'].std(), color='purple', linestyle='--', 
                            label=f'+1 Std Dev: {(data_df["diff"].mean() + data_df["diff"].std()):.2f}')
                plt.text(data_df['diff'].mean() - data_df['diff'].std(), plt.ylim()[1] * 0.95, 
                         f'-1 SD: {data_df["diff"].std():.2f}', horizontalalignment='right', color='purple')
                plt.text(data_df['diff'].mean() + data_df['diff'].std(), plt.ylim()[1] * 0.95, 
                         f'+1 SD: {data_df["diff"].std():.2f}', horizontalalignment='left', color='purple')
                plt.title(f'Distribution of Diff Values for {unique_id} ({platform_name})')
                plt.xlabel('Diff')
                plt.ylabel('Frequency')
                plt.legend()
                plt.savefig(os.path.join(hist_export_dir, f"{unique_id}_histogram.png"))
                plt.close()

    print("Completed Step 1.")

def run_apply_best_offsets(db_path):
    print("Starting Post-Processing Step 2: Apply Best Offsets")
    with duckdb.connect(database=db_path, read_only=False) as con:
        # Check if columns exist; if not, add them. This is safer for iterative runs.
        columns_df = con.execute("DESCRIBE csb").fetchdf()
        if 'depth_mod' not in columns_df['column_name'].values:
            con.execute("ALTER TABLE csb ADD COLUMN depth_mod DOUBLE;")
        if 'uncertainty_vert' not in columns_df['column_name'].values:
            con.execute("ALTER TABLE csb ADD COLUMN uncertainty_vert DOUBLE;")
        if 'uncertainty_hori' not in columns_df['column_name'].values:
            con.execute("ALTER TABLE csb ADD COLUMN uncertainty_hori DOUBLE;")

        # --- MODIFIED: Apply offsets only to new rows where depth_mod IS NULL ---
        print("Applying best offsets to new data...")
        con.execute("""
        UPDATE csb
        SET depth_mod = (depth_new - sub.average_diff) * -1
        FROM (
            SELECT unique_id, AVG(diff) AS average_diff
            FROM csb
            WHERE diff > -12 AND diff < 12 AND Raster_Value > -20 AND Uncertainty_Value < 3
            GROUP BY unique_id
        ) AS sub
        WHERE csb.unique_id = sub.unique_id AND csb.depth_mod IS NULL;
        """)
        
        # Fallback for points that couldn't be calibrated against reference data
        update_query = """
        UPDATE csb
        SET depth_mod = depthfinal
        WHERE depth_mod IS NULL; 
        """
        con.execute(update_query)
        print("Updated remaining depth_mod values with initial depthfinal.")

        # --- MODIFIED: Calculate uncertainty only for new rows ---
        print("Calculating CATZOC uncertainty for new data...")
        uncert_vert_query = "UPDATE csb SET uncertainty_vert = (2 + (depth_mod * -0.05)) WHERE uncertainty_vert IS NULL"
        uncert_hori_query = "UPDATE csb SET uncertainty_hori = 10 WHERE uncertainty_hori IS NULL"
        con.execute(uncert_vert_query)
        con.execute(uncert_hori_query)
        print("Uncertainty values calculated.")

    print("Completed Step 2.")

def run_export_transits(db_path, exports_folder):
    print("Starting Post-Processing Step 3: Outlier Detection & Transit ID Assignment")
    
    # --- Helper Functions ---
    def haversine(lat1, lon1, lat2, lon2):
        R = 6371000
        phi1, phi2 = np.radians(lat1), np.radians(lat2)
        dphi = np.radians(lat2 - lat1)
        dlambda = np.radians(lon2 - lon1)
        a = np.sin(dphi / 2)**2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2)**2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        return R * c

    def calculate_vessel_speed(group):
        group = group.sort_values('time').copy()
        group['lat'] = pd.to_numeric(group['lat'], errors='coerce')
        group['lon'] = pd.to_numeric(group['lon'], errors='coerce')
        group['time_diff'] = group['time'].diff().dt.total_seconds()
        group['distance'] = haversine(
            group['lat'].shift(), group['lon'].shift(),
            group['lat'], group['lon']
        )
        group['vessel_speed'] = group['distance'] / group['time_diff']
        group['vessel_speed'] = group['vessel_speed'].replace([np.inf, -np.inf], 0)
        group['vessel_speed'] = group['vessel_speed'].fillna(0)
        group['vessel_speed_smoothed'] = uniform_filter1d(group['vessel_speed'], size=5)
        group['vessel_speed_smoothed'] = group['vessel_speed_smoothed'].replace([np.inf, -np.inf], 0)
        return group

    def detect_outliers(data, scaler, threshold_percentile, original_gdf, return_smoothed=False):
        data_scaled = scaler.fit_transform(data)
        imputer = IterativeImputer(
            estimator=LinearRegression(),
            max_iter=15,
            random_state=42,
            sample_posterior=False
        )
        data_imputed = imputer.fit_transform(data_scaled)
        smoothed_depth = uniform_filter1d(data_imputed[:, 2], size=50)
        residuals = np.abs(data_scaled[:, 2] - smoothed_depth)
        threshold = np.percentile(residuals, threshold_percentile)
        outliers = residuals > threshold
        smoothed_depth_denorm = smoothed_depth * scaler.scale_[2] + scaler.mean_[2]
        original_gdf.loc[data.index[outliers], 'Outlier'] = True
        outlier_count = np.sum(outliers)
        if return_smoothed:
            full_smoothed_depth = pd.Series(index=original_gdf.index, dtype=float)
            full_smoothed_depth[data.index] = smoothed_depth_denorm
            return full_smoothed_depth, outlier_count
        return data[~outliers], outlier_count

    def create_geotiff(gdf, filename, resolution=8):
        try:
            bounds = gdf.total_bounds
            x_min, y_min, x_max, y_max = bounds
            if x_max == x_min or y_max == y_min:
                raise ValueError("Invalid geographic bounds. All points may be identical or too close.")
            x_res = int((x_max - x_min) / resolution)
            y_res = int((y_max - y_min) / resolution)
            transform = from_origin(x_min, y_max, resolution, resolution)
            out_meta = {
                'driver': 'GTiff',
                'height': y_res,
                'width': x_res,
                'count': 2,
                'dtype': 'float32',
                'crs': gdf.crs.to_string(),
                'transform': transform,
                'nodata': 1000000,
                'compress': 'lzw',
                'interleave': 'band'
            }
            with rasterio.open(filename, "w", **out_meta) as dest:
                for idx, col in enumerate(['depth', 'uncertainty'], start=1):
                    array = np.full((y_res, x_res), out_meta['nodata'], dtype='float32')
                    for point, value in zip(gdf.geometry, gdf[col]):
                        col_idx = int((point.x - x_min) / resolution)
                        row_idx = int((y_max - point.y) / resolution)
                        if 0 <= col_idx < x_res and 0 <= row_idx < y_res:
                            array[row_idx, col_idx] = value
                    dest.write(array, idx)
        except Exception as e:
            print(f"Failed to create GeoTIFF for {filename}. Error: {str(e)}")

    def create_transit_ids(df, max_hours_gap, max_days_duration):
        df = df.sort_values(by='time')
        current_transit_id = None
        last_time = None
        current_start_time = None
        transit_ids = []
        for index, row in df.iterrows():
            if last_time is None or (row['time'] - last_time > timedelta(hours=max_hours_gap)) \
               or ((row['time'] - current_start_time) > timedelta(days=max_days_duration)):
                current_transit_id = f"{row['unique_id']}_{row['time'].strftime('%Y-%m-%d_%H-%M-%S')}"
                current_start_time = row['time']
            transit_ids.append(current_transit_id)
            last_time = row['time']
        df['transit_id'] = transit_ids
        return df

    os.makedirs(exports_folder, exist_ok=True)
    MAX_HOURS_GAP = 4
    MAX_DAYS_DURATION = 7
    with duckdb.connect(database=db_path, read_only=False) as con:
        con.install_extension('spatial')
        con.load_extension('spatial')

        # Ensure columns exist
        columns_df = con.execute("DESCRIBE csb").fetchdf()
        if 'Outlier' not in columns_df['column_name'].values:
            con.execute("ALTER TABLE csb ADD COLUMN Outlier BOOLEAN DEFAULT FALSE;")
        if 'transit_id' not in columns_df['column_name'].values:
            con.execute("ALTER TABLE csb ADD COLUMN transit_id VARCHAR;")
        if 'vessel_speed_smoothed' not in columns_df['column_name'].values:
            con.execute("ALTER TABLE csb ADD COLUMN vessel_speed_smoothed DOUBLE;")

        # --- MODIFIED: Get unique_ids ONLY from unprocessed data ---
        unprocessed_ids_query = "SELECT DISTINCT unique_id FROM csb WHERE transit_id IS NULL"
        unique_ids = con.execute(unprocessed_ids_query).fetchall()
        unique_ids = [x[0] for x in unique_ids]
        
        if not unique_ids:
            print("No new data to process for outlier detection. Skipping.")
            print("Completed Step 3.")
            return
            
        print(f"Found {len(unique_ids)} unique_id(s) with new data to process.")

        for i, unique_id in enumerate(unique_ids, start=1):
            print(f"\nProcessing new data for unique_id {i}/{len(unique_ids)}: {unique_id}")
            # --- MODIFIED: Query ONLY unprocessed data for the current vessel ---
            query = f"""
            SELECT
                rowid, unique_id, platform_name_x AS platform_name, time,
                depth_mod AS depth, uncertainty_vert AS uncertainty, uncertainty_hori,
                lat, lon, Raster_Value
            FROM csb
            WHERE unique_id = '{unique_id}'
              AND depth_mod IS NOT NULL
              AND transit_id IS NULL  -- The key change!
            ORDER BY time;
            """
            df = con.execute(query).df()
            if df.empty:
                print(f"  No new data with depth_mod for {unique_id}, skipping.")
                continue
            
            df['time'] = pd.to_datetime(df['time'], format='%Y%m%d %H:%M:%S')
            df = create_transit_ids(df, MAX_HOURS_GAP, MAX_DAYS_DURATION)
            for transit_id, group in df.groupby('transit_id'):
                print(f"  Processing transit: {transit_id}")
                group = group.sort_values('time').copy()
                group['lat'] = pd.to_numeric(group['lat'], errors='coerce')
                group['lon'] = pd.to_numeric(group['lon'], errors='coerce')
                group['Outlier'] = False
                group = calculate_vessel_speed(group)
                excess_speed_points = group[group['vessel_speed_smoothed'] > 10.3]
                print(f"    {len(excess_speed_points)} points exceed 20 knots (not filtered out).")
                data_for_outlier = group[['lat', 'lon', 'depth']].copy()
                scaler = StandardScaler()
                print("    First Pass (99th percentile):")
                filtered_data_1, outlier_count_1 = detect_outliers(data_for_outlier.copy(), scaler, 99, original_gdf=group)
                print(f"    Outliers detected in Pass 1: {outlier_count_1}")
                print("    Second Pass (98th percentile):")
                filtered_data_2, outlier_count_2 = detect_outliers(filtered_data_1.copy(), scaler, 98, original_gdf=group)
                print(f"    Outliers detected in Pass 2: {outlier_count_2}")
                print("    Third Pass (98th percentile, strict):")
                final_smoothed_depth, outlier_count_3 = detect_outliers(filtered_data_2.copy(), scaler, 98, original_gdf=group, return_smoothed=True)
                print(f"    Outliers detected in Pass 3: {outlier_count_3}")
                group['Final_Smoothed_Depth'] = final_smoothed_depth

                plot_filename = os.path.join(exports_folder, f"{unique_id}_{transit_id}_outlier_plot.png")
                fig, ax = plt.subplots(figsize=(12, 6))
                mask_valid = group['Outlier'] == False
                mask_outlier = group['Outlier'] == True
                ax.scatter(group.index[mask_valid], group.loc[mask_valid, 'depth'],
                           color='blue', s=1, label="Valid")
                ax.scatter(group.index[mask_outlier], group.loc[mask_outlier, 'depth'],
                           color='red', s=10, label="Outlier")
                ax.set_title(f"Outlier Detection for Transit {transit_id} (unique_id: {unique_id})")
                ax.set_xlabel("Record Index")
                ax.set_ylabel("Depth")
                ax.legend()
                plt.savefig(plot_filename)
                plt.close()
                print(f"    Saved outlier plot to {plot_filename}")

                updates = []
                for idx, row in group.iterrows():
                    updates.append((row['rowid'], row['transit_id'], str(row['Outlier']).upper()))
                if updates:
                    updates_df = pd.DataFrame(updates, columns=['rowid', 'transit_id', 'outlier'])
                    con.register('updates_df', updates_df)
                    con.execute("""
                        UPDATE csb
                        SET transit_id = updates_df.transit_id,
                            Outlier = CAST(updates_df.outlier AS BOOLEAN)
                        FROM updates_df
                        WHERE csb.rowid = updates_df.rowid;
                    """)
                    con.unregister('updates_df')
                    print(f"    ...updating {len(updates)} rows for Outlier and transit_id.")
                del updates
                gc.collect()

                speed_updates = []
                for idx, row in group.iterrows():
                    speed_updates.append((row['rowid'], row['vessel_speed_smoothed']))
                if speed_updates:
                    speed_updates_df = pd.DataFrame(speed_updates, columns=['rowid', 'speed'])
                    con.register('speed_updates_df', speed_updates_df)
                    con.execute("""
                        UPDATE csb
                        SET vessel_speed_smoothed = speed_updates_df.speed
                        FROM speed_updates_df
                        WHERE csb.rowid = speed_updates_df.rowid;
                    """)
                    con.unregister('speed_updates_df')
                    print(f"    ...updating {len(speed_updates)} rows for vessel_speed_smoothed.")
                del speed_updates
                gc.collect()

                # --- Transit export is optional ---
                if export_transits_var.get():
                    gdf = gpd.GeoDataFrame(group, geometry=gpd.points_from_xy(group.lon, group.lat))
                    gdf.set_crs(epsg=4326, inplace=True)
                    if gdf.empty:
                        print(f"No valid points left in transit {transit_id}. Skipping export.")
                        continue
                    
                    non_outlier_gdf = gdf[gdf['Outlier'] == False]
                    if non_outlier_gdf.empty:
                        print(f"No non-outlier points left in transit {transit_id}. Skipping export.")
                        continue

                    avg_lat = non_outlier_gdf['lat'].mean()
                    avg_lon = non_outlier_gdf['lon'].mean()
                    try:
                        epsg_zone = get_utm_zone_nad83(avg_lat, avg_lon)
                    except ValueError as ve:
                        print(f"Could not determine NAD83 UTM zone for lat={avg_lat}, lon={avg_lon}: {ve}")
                        continue

                    zone_folder = os.path.join(exports_folder, f"zone_{epsg_zone}")
                    os.makedirs(zone_folder, exist_ok=True)
                    non_outlier_gdf.to_crs(epsg=epsg_zone, inplace=True)

                    start_date = group['time'].min().strftime('%Y%m%d%H%M%S')
                    end_date = group['time'].max().strftime('%Y%m%d%H%M%S')
                    gpkg_filename = f"{unique_id}_{transit_id}_{start_date}_{end_date}.gpkg"
                    gpkg_path = os.path.join(zone_folder, gpkg_filename)

                    non_outlier_gdf.to_file(gpkg_path, driver='GPKG')
                    print(f"Exported GeoPackage {gpkg_path}")

                    tiff_filename = gpkg_filename.replace('.gpkg', '.tif')
                    tiff_path = os.path.join(zone_folder, tiff_filename)
                    create_geotiff(non_outlier_gdf, tiff_path)
                    print(f"Exported GeoTIFF {tiff_path}")
                    del gdf
                
                del group
                gc.collect()

    print("Completed Step 3.")

# --- END: POST-PROCESSING ANALYSIS FUNCTIONS ---

# --- START: FINAL GRIDDING AND EXPORT FUNCTIONS ---

def points_to_raster_average(gdf, out_raster_path, value_col='depth', nodata=1000000):
    """
    Creates a GeoTIFF by averaging point values within each grid cell.
    """
    resolution = float(grid_resolution_var.get())
    if gdf.crs is None:
        raise ValueError("GeoDataFrame has no CRS. Please set or reproject first.")

    x_min, y_min, x_max, y_max = gdf.total_bounds
    width = int(np.ceil((x_max - x_min) / resolution))
    height = int(np.ceil((y_max - y_min) / resolution))

    if width <= 0 or height <= 0:
        print(f"Warning: Raster dimensions are zero or negative for {os.path.basename(out_raster_path)}. Skipping.")
        return

    transform = from_origin(x_min, y_max, resolution, resolution)

    sum_array   = np.zeros((height, width), dtype=np.float32)
    count_array = np.zeros((height, width), dtype=np.float32)

    for geom, value in zip(gdf.geometry, gdf[value_col]):
        if geom is None or pd.isna(value):
            continue
        col = int((geom.x - x_min) // resolution)
        row = int((y_max - geom.y) // resolution)
        if 0 <= col < width and 0 <= row < height:
            sum_array[row, col]   += value
            count_array[row, col] += 1

    avg_array = np.where(count_array > 0, sum_array / count_array, nodata)

    with rasterio.open(
            out_raster_path, 'w', driver='GTiff',
            height=height, width=width, count=1,
            dtype=np.float32, crs=gdf.crs.to_string(),
            transform=transform, nodata=nodata,
            compress='lzw'
    ) as dst:
        dst.write(avg_array, 1)

    print(f"Final gridded GeoTIFF created at {out_raster_path}")

def organize_by_epsg(input_dir):
    """
    Moves TIFF files into subdirectories named by their EPSG code.
    """
    print("\nOrganizing final GeoTIFFs by EPSG code...")
    for fn in os.listdir(input_dir):
        if not fn.lower().endswith(('.tif', '.tiff')):
            continue

        src_path = os.path.join(input_dir, fn)
        if not os.path.isfile(src_path): continue

        try:
            with rasterio.open(src_path) as src:
                epsg = src.crs.to_epsg()
        except Exception as e:
            print(f"[ERROR] could not open {fn}: {e}")
            continue

        if epsg is None:
            print(f"[WARN] {fn} has no recognized EPSG code, skipping.")
            continue

        out_folder = os.path.join(input_dir, f"EPSG_{epsg}")
        os.makedirs(out_folder, exist_ok=True)

        dest_path = os.path.join(out_folder, fn)
        shutil.move(src_path, dest_path)
        print(f"Moved {fn} → {out_folder}")

def create_vrts_for_epsg_folders(base_dir):
    """
    Scans for 'EPSG_' subfolders and builds a VRT for the TIFFs in each.
    """
    print("\nBuilding VRTs for each EPSG folder...")
    for folder_name in os.listdir(base_dir):
        if folder_name.startswith("EPSG_") and os.path.isdir(os.path.join(base_dir, folder_name)):
            directory = os.path.join(base_dir, folder_name)
            pattern = os.path.join(directory, "*.tif")
            files = glob.glob(pattern)

            if not files:
                print(f"No GeoTIFFs found in {directory}, skipping VRT creation.")
                continue

            print(f"Found {len(files)} files in {directory}. Building VRT...")
            vrt_filename = os.path.join(directory, f"mosaic_{folder_name}.vrt")
            
            gdal.BuildVRT(vrt_filename, files)
            print("VRT created:", vrt_filename)

            # Build overviews
            print("Building overviews for VRT...")
            ds = gdal.Open(vrt_filename)
            if ds:
                gdal.SetConfigOption('COMPRESS_OVERVIEW', 'LZW')
                ds.BuildOverviews("AVERAGE", [2,4,8,16,32,64])
                ds = None
                print(f"Overviews built for {vrt_filename}")


def run_final_gridding_and_export():
    """
    Main function for the final gridding and export stage.
    """
    print("\n***** Starting Final Gridding & Export Stage *****")
    db_path = os.path.join(output_dir, "csb.duckdb")
    tessellation_shp = tessellation_shp_var.get()
    output_folder = os.path.join(output_dir, "final_products")
    os.makedirs(output_folder, exist_ok=True)
    with duckdb.connect(database=db_path, read_only=False) as con:
        if tessellation_shp and os.path.exists(tessellation_shp):
            print(f"Using tessellation scheme: {tessellation_shp}")
            polygons_gdf = gpd.read_file(tessellation_shp)
            if polygons_gdf.crs.to_epsg() != 4326:
                polygons_gdf = polygons_gdf.to_crs(epsg=4326)

            for idx, poly_row in polygons_gdf.iterrows():
                poly_geom = poly_row.geometry
                polygon_id = str(poly_row.get('GRID_ID', idx))
                print(f"\n--- Processing Tile: {polygon_id} ---")

                minx, miny, maxx, maxy = poly_geom.bounds
                query = f"""
                    SELECT lat, lon, "Outlier" AS outlier, depth_mod AS depth, unique_id,
                    platform_name_x as platform_name, time, uncertainty_vert
                    FROM csb
                    WHERE (lat BETWEEN {miny} AND {maxy} AND lon BETWEEN {minx} AND {maxx})
                    AND (Raster_Value IS NULL OR ABS(Raster_Value - depth_mod) <= (uncertainty_vert * 3.5))
                """
                df_points = con.execute(query).df()
                

                print(f"  Found {len(df_points)} points within the bounding box that passed the quality filter.")

                if df_points.empty:
                    continue

                points_gdf_4326 = gpd.GeoDataFrame(
                    df_points,
                    geometry=[Point(xy) for xy in zip(df_points.lon, df_points.lat)],
                    crs="EPSG:4326"
                )

                points_gdf_4326 = points_gdf_4326[points_gdf_4326.geometry.within(poly_geom)]
                

                print(f"  {len(points_gdf_4326)} points remain after precise clipping to the polygon.")

                if points_gdf_4326.empty:
                    continue
                
                if export_final_gpkg_var.get():
                    gpkg_path = os.path.join(output_folder, f"{polygon_id}_points.gpkg")
                    print(f"  Saving {len(points_gdf_4326)} points to GeoPackage...")
                    points_gdf_4326.to_file(gpkg_path, driver="GPKG")
                    print(f"  Saved points GeoPackage (EPSG:4326): {gpkg_path}")

                lat_c, lon_c = poly_geom.centroid.y, poly_geom.centroid.x
                try:
                    epsg_zone = get_utm_zone_nad83(lat_c, lon_c)
                    points_gdf_utm = points_gdf_4326.to_crs(epsg=epsg_zone)
                    
                    points_for_raster = points_gdf_utm[points_gdf_utm['outlier'] == False]

                    print(f"  Found {len(points_for_raster)} non-outlier points to create raster from.")

                    if not points_for_raster.empty:
                        tif_path = os.path.join(output_folder, f"{polygon_id}_gridded.tif")
                        points_to_raster_average(points_for_raster, tif_path, value_col='depth')
                except ValueError as e:
                    print(f"  Skipping raster for {polygon_id}: {e}")
                    continue

        else: # This is the case for a single file output
            print("No tessellation scheme provided. Processing all data into a single file.")
            query = """
                SELECT lat, lon, "Outlier" AS outlier, depth_mod AS depth, unique_id,
                platform_name_x as platform_name, time, uncertainty_vert 
                FROM csb 
                WHERE 
                    Raster_Value IS NULL 
                    OR 
                    ABS(Raster_Value - depth_mod) <= (uncertainty_vert * 3.5)
            """
            df_points = con.execute(query).df()


            print(f"Initial query returned {len(df_points)} points from the database.")

            if df_points.empty:
                print("No points passed the quality filter. Aborting final gridding.")
                print("This can happen if all points intersected reference data but failed the quality check.")
                return
                
            points_gdf_4326 = gpd.GeoDataFrame(
                df_points,
                geometry=[Point(xy) for xy in zip(df_points.lon, df_points.lat)],
                crs="EPSG:4326"
            )
            
            if export_final_gpkg_var.get():
                gpkg_path = os.path.join(output_folder, "csb_final_points.gpkg")
                print(f"Saving {len(points_gdf_4326)} points to GeoPackage...")
                points_gdf_4326.to_file(gpkg_path, driver="GPKG")
                print(f"Saved final points GeoPackage (EPSG:4326): {gpkg_path}")
            
            try:
                world_centroid = points_gdf_4326.unary_union.centroid
                epsg_zone = get_utm_zone_nad83(world_centroid.y, world_centroid.x)
                print(f"Determined overall EPSG zone for raster as: {epsg_zone}")
                points_gdf_utm = points_gdf_4326.to_crs(epsg=epsg_zone)
                
                points_for_raster = points_gdf_utm[points_gdf_utm['outlier'] == False]

                print(f"Found {len(points_for_raster)} non-outlier points to create raster from.")

                if not points_for_raster.empty:
                    tif_path = os.path.join(output_folder, "csb_final_gridded.tif")
                    points_to_raster_average(points_for_raster, tif_path, value_col='depth')
                else:
                    print("No valid non-outlier points to create final raster.")
                
            except ValueError as e:
                print(f"Could not process single file output: {e}")

    if organize_vrt_var.get():
        organize_by_epsg(output_folder)
        create_vrts_for_epsg_folders(output_folder)

    print("***** Final Gridding & Export Stage Complete *****")


def cleanup_interim_files(output_dir, title):
    """
    Safely removes all temporary and intermediate files and folders created during processing.
    """
    print(f"\n--- Cleaning up interim files for {title} ---")
    
    # List of folder paths to remove
    folders_to_remove = [
        os.path.join(output_dir, "Modeling")
    ]
    
    for folder in folders_to_remove:
        try:
            if os.path.exists(folder):
                shutil.rmtree(folder)
                print(f"Removed folder: {folder}")
        except Exception as e:
            print(f"Error removing folder {folder}: {e}")

    # List of file patterns to remove
    file_patterns_to_remove = [
        os.path.join(output_dir, f"{title}_5m_MLLW.tif"),
        os.path.join(output_dir, f"{title}_wgs84.tif"),
        os.path.join(output_dir, f"{title}_intermediate.tif"),
        os.path.join(output_dir, "convex_hull_polygon.*"),
        os.path.join(output_dir, f"{title}_bathy_polygon.*")
    ]
    
    for pattern in file_patterns_to_remove:
        files = glob.glob(pattern)
        for f in files:
            try:
                if os.path.exists(f):
                    os.remove(f)
                    print(f"Removed file: {f}")
            except Exception as e:
                print(f"Error removing file {f}: {e}")

# --- END: FINAL GRIDDING AND EXPORT FUNCTIONS ---
    
def open_file_dialog(var, file_types):
    filename = filedialog.askopenfilename(filetypes=file_types)
    var.set(filename)

def open_folder_dialog(var):
    foldername = filedialog.askdirectory()
    var.set(foldername)

def process_csb_threaded():
    csb_directory = csb_var.get()
    if not os.path.isdir(csb_directory):
        print("Selected path is not a directory.")
        return
    processing_thread = threading.Thread(target=process_csb_for_directory, args=(csb_directory, root))
    processing_thread.start()

def process_csb_for_directory(csb_directory, root_window):
    start_time = time.time()
    global csb, title, fp_zones, use_bluetopo, output_dir, BAG_filepath
    csb = csb_var.get()
    fp_zones = fp_zones_var.get()
    output_dir = output_dir_var.get()
    csb_directory = os.path.abspath(csb_directory)
    fp_zones = os.path.abspath(fp_zones)
    output_dir = os.path.abspath(output_dir)
    print(f"output_dir is: {output_dir}")

    # --- NEW: Load master offsets once at the start ---
    MASTER_OFFSET_FILE = os.path.join(output_dir, "master_offsets.csv")
    if os.path.exists(MASTER_OFFSET_FILE):
        print(f"Found existing master offsets file at: {MASTER_OFFSET_FILE}")
        master_offsets_df = pd.read_csv(MASTER_OFFSET_FILE)
    else:
        print("No master_offsets.csv found. Will create a new one.")
        master_offsets_df = pd.DataFrame(columns=['unique_id', 'platform_name', 'offset_value', 'std_dev', 'accuracy_score', 'date_range', 'tile_name'])

    csv_files = [os.path.join(csb_directory, f) for f in os.listdir(csb_directory) if f.endswith('.csv')]
    for csb_file in csv_files:
        csb = csb_file
        title = os.path.splitext(os.path.basename(csb_file))[0]
        # We check for a final product to determine if we should skip
        final_product_check = os.path.join(output_dir, "final_products", "csb_final_gridded.tif")
        if os.path.exists(final_product_check):
            print(f"Skipping already processed file based on existing final products: {title}")
            continue

        print(f"Processing {csb} with title: {title}")
        try:
            modeling_dir_path = os.path.join(output_dir, "Modeling")
            try:
                shutil.rmtree(modeling_dir_path)
            except FileNotFoundError:
                pass # It's ok if it doesn't exist
            except Exception as e:
                print(f"Error deleting old Modeling folder: {e}")
            
            if bluetopo_var.get():
                BAG_filepath = create_convex_hull_and_download_tiles(csb, output_dir, use_bluetopo=True)
            else:
                BAG_filepath = BAG_filepath_var.get()
            
            # --- MODIFIED: Pass the master_offsets_df to rasterize_CSB ---
            rasterize_CSB(master_offsets_df)
            
        except Exception as e:
            print(f"An error occurred during initial processing of {csb}: {e}")
        finally:
            cleanup_interim_files(output_dir, title)

    if run_analysis_var.get():
        print("\n***** Starting Post-Processing Analysis *****")
        db_path = os.path.join(output_dir, "csb.duckdb")
        hist_export_dir = os.path.join(output_dir, "histograms")
        exports_folder = os.path.join(output_dir, "transit_exports")

        if not os.path.exists(db_path):
            print(f"Error: DuckDB file not found at {db_path}. Cannot run analysis.")
        else:
            run_histograms_calibration_points(db_path, hist_export_dir)
            run_apply_best_offsets(db_path)
            run_export_transits(db_path, exports_folder) # This now checks internally if it should run
            
    if run_final_grid_var.get():
        run_final_gridding_and_export()
            
    end_time = time.time()
    duration = end_time - start_time
    minutes, seconds = divmod(duration, 60)
    print(f"***** ALL STAGES DONE! Total processing time: {int(minutes)} minutes and {seconds:.1f} seconds")
    print("processing complete. closing application window.")
    root_window.destroy()


# --- GUI SETUP ---
root = tk.Tk()
root.title("CSB Processing Pipeline")
root.geometry("900x650")

# --- Main Frame ---
main_frame = ttk.Frame(root, padding="10")
main_frame.pack(fill=tk.BOTH, expand=True)

# --- Variable Definitions ---
csb_var = tk.StringVar()
BAG_filepath_var = tk.StringVar()
fp_zones_var = tk.StringVar()
output_dir_var = tk.StringVar()
tessellation_shp_var = tk.StringVar()
grid_resolution_var = tk.StringVar(value="10")
export_final_gpkg_var=tk.BooleanVar(value=True)


# BooleanVars for options
bluetopo_var = tk.BooleanVar()
duckdb_option_var = tk.BooleanVar(value=True)
export_gp_var = tk.BooleanVar(value=False)
run_analysis_var = tk.BooleanVar(value=True)
export_transits_var = tk.BooleanVar(value=False)
run_final_grid_var = tk.BooleanVar(value=True)
organize_vrt_var = tk.BooleanVar(value=True)

# --- Input Files Section ---
input_frame = ttk.LabelFrame(main_frame, text="1. Input Files & Folders", padding="10")
input_frame.pack(fill=tk.X, expand=True, pady=5)

ttk.Label(input_frame, text='Raw CSB CSV Directory').grid(row=0, column=0, sticky='w', padx=5, pady=2)
csb_dir_entry = ttk.Entry(input_frame, textvariable=csb_var, width=60)
csb_dir_entry.grid(row=0, column=1)
ttk.Button(input_frame, text='Browse', command=lambda: open_folder_dialog(csb_var)).grid(row=0, column=2, padx=5)

ttk.Label(input_frame, text='Tide Zone Shapefile').grid(row=1, column=0, sticky='w', padx=5, pady=2)
fp_zones_entry = ttk.Entry(input_frame, textvariable=fp_zones_var, width=60)
fp_zones_entry.grid(row=1, column=1)
ttk.Button(input_frame, text='Browse', command=lambda: open_file_dialog(fp_zones_var, [("Shapefile", "*.shp")])).grid(row=1, column=2, padx=5)

ttk.Label(input_frame, text='Output Directory').grid(row=2, column=0, sticky='w', padx=5, pady=2)
output_dir_entry = ttk.Entry(input_frame, textvariable=output_dir_var, width=60)
output_dir_entry.grid(row=2, column=1)
ttk.Button(input_frame, text='Browse', command=lambda: open_folder_dialog(output_dir_var)).grid(row=2, column=2, padx=5)

# --- Reference Bathymetry Section ---
ref_bathy_frame = ttk.LabelFrame(main_frame, text="2. Reference Bathymetry", padding="10")
ref_bathy_frame.pack(fill=tk.X, expand=True, pady=5)

bluetopo_checkbox = ttk.Checkbutton(ref_bathy_frame, text="Use Automated BlueTopo Download", variable=bluetopo_var)
bluetopo_checkbox.grid(row=0, column=0, columnspan=3, sticky='w', padx=5, pady=2)

ttk.Label(ref_bathy_frame, text='Local BAG/GeoTiff File').grid(row=1, column=0, sticky='w', padx=5, pady=2)
BAG_filepath_entry = ttk.Entry(ref_bathy_frame, textvariable=BAG_filepath_var, width=60)
BAG_filepath_entry.grid(row=1, column=1)
ttk.Button(ref_bathy_frame, text='Browse', command=lambda: open_file_dialog(BAG_filepath_var, [("BAG or GeoTIFF file", "*.bag;*.tif;*.tiff")])).grid(row=1, column=2, padx=5)

# --- Processing Options Section ---
options_frame = ttk.LabelFrame(main_frame, text="3. Processing & Export Options", padding="10")
options_frame.pack(fill=tk.X, expand=True, pady=5)

ttk.Checkbutton(options_frame, text="Insert into DuckDB (Required for all post-processing)", variable=duckdb_option_var).grid(row=0, column=0, columnspan=2, sticky='w', padx=5)
ttk.Checkbutton(options_frame, text="Export Initial Processed Geopackage (per input file)", variable=export_gp_var).grid(row=1, column=0, columnspan=2, sticky='w', padx=5)
ttk.Separator(options_frame, orient='horizontal').grid(row=2, columnspan=3, sticky='ew', pady=5)

# Post-Processing Options
ttk.Checkbutton(options_frame, text="Run Post-Processing (Outlier Flagging, etc.)", variable=run_analysis_var).grid(row=3, column=0, columnspan=2, sticky='w', padx=5)
export_transits_checkbox = ttk.Checkbutton(options_frame, text="Export Individual Transit Files (GPKG & GeoTIFF)", variable=export_transits_var)
export_transits_checkbox.grid(row=4, column=0, sticky='w', padx=25)
ttk.Separator(options_frame, orient='horizontal').grid(row=5, columnspan=3, sticky='ew', pady=5)

# Final Gridding Options
final_grid_checkbox = ttk.Checkbutton(options_frame, text="Run Final Gridding & Export", variable=run_final_grid_var)
final_grid_checkbox.grid(row=6, column=0, sticky='w', padx=5)

gpkg_export_checkbox=ttk.Checkbutton(options_frame, text="Export Final Points GeoPackage (in epsg:4326)", variable=export_final_gpkg_var)
gpkg_export_checkbox.grid(row=7,column=0,columnspan=2,sticky='w',padx=25)

ttk.Label(options_frame, text='Optional Tessellation Shapefile').grid(row=8, column=0, sticky='w', padx=25)
tess_entry = ttk.Entry(options_frame, textvariable=tessellation_shp_var, width=45)
tess_entry.grid(row=8, column=1, sticky='w')
tess_button = ttk.Button(options_frame, text='Browse', command=lambda: open_file_dialog(tessellation_shp_var, [("Shapefile", "*.shp")]))
tess_button.grid(row=8, column=2, padx=5)

ttk.Label(options_frame, text='Grid Resolution (meters)').grid(row=9, column=0, sticky='w', padx=25)
res_entry = ttk.Entry(options_frame, textvariable=grid_resolution_var, width=10)
res_entry.grid(row=9, column=1, sticky='w')

organize_vrt_checkbox = ttk.Checkbutton(options_frame, text="Organize GeoTIFFs by EPSG and Create VRTs", variable=organize_vrt_var)
organize_vrt_checkbox.grid(row=10, column=0, columnspan=2, sticky='w', padx=25)

# --- Process Button ---
process_button = ttk.Button(main_frame, text='Start Processing', command=process_csb_threaded)
process_button.pack(pady=15)

# --- GUI Logic ---
def on_bluetopo_check():
    if bluetopo_var.get():
        BAG_filepath_entry.config(state='disabled')
        BAG_filepath_var.set("")
    else:
        BAG_filepath_entry.config(state='normal')

def toggle_analysis_options():
    state = 'normal' if run_analysis_var.get() else 'disabled'
    export_transits_checkbox.config(state=state)

def toggle_final_grid_options():
    state = 'normal' if run_final_grid_var.get() else 'disabled'
    tess_entry.config(state=state)
    tess_button.config(state=state)
    res_entry.config(state=state)
    organize_vrt_checkbox.config(state=state)
    if state == 'disabled':
        tessellation_shp_var.set("")

bluetopo_checkbox.config(command=on_bluetopo_check)
run_analysis_var.trace_add('write', lambda *args: toggle_analysis_options())
run_final_grid_var.trace_add('write', lambda *args: toggle_final_grid_options())

# Initial state setup
on_bluetopo_check()
toggle_analysis_options()
toggle_final_grid_options()

if __name__ == "__main__":
    root.mainloop()