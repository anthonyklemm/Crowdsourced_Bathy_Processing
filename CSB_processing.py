# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 14:51:41 2023

@author: anthonyklemm
"""

import tkinter as tk
from tkinter import filedialog
import geopandas as gpd
import pandas as pd
import numpy as np
import requests
import os
from osgeo import gdal, osr
import subprocess
import rasterio
from scipy.interpolate import interp1d
from rasterio.features import shapes, rasterize
from shapely.geometry import shape
from rasterio.transform import from_origin
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rasterio.plot import show
import time
import threading
from rasterio.merge import merge
from rasterio.enums import Resampling
from shapely.ops import unary_union
from skimage.morphology import binary_dilation, binary_erosion



# Global variables
title = ""
#output_crs = ""
csb = ""
BAG_filepath = ""
fp_zones = ""
output_dir = ""
resolution = 20 # resolution of quick-look geotiff raster

pd.set_option('display.max_columns', None)

def loadCSB():

    print('*****Reading CSB input csv file***** ')
    df1=gpd.read_file(csb) #read in CSB data in CSV
    #df1 = df1.head(30000)  # Limit to top 0000 records for debugging
    
    df1 = df1.astype({'depth' : 'float'}) #turn depth field numeric (float64)
    df2 = df1[df1['depth'] > .5] #remove depth values less that 0.5m (filtering out probable noise)
    df2 = df2[df2['depth'] < 1000] #remove depths values greater than possibleNYC depths in the area - this value can change based on area)
    # Data cleaning - remove records with erroneous datetimes
    df2 = df2.astype({'time':'datetime64'}, errors='ignore')
    df2 = df2.dropna(subset=['time'])
    df2 = df2[df2['time'] > '2014']
    df2 = df2[df2['time'] < '2025']
    
    df2 = df2.drop_duplicates(subset=['lon', 'lat', 'depth', 'time', 'unique_id'])
    gdf = gpd.GeoDataFrame(df2, geometry=gpd.points_from_xy(df2.lon, df2.lat)) #create vector point file from lat and lon
    gdf = gdf.set_crs(4326, allow_override=True)
    return gdf

def create_convex_hull_and_download_tiles(csb_data_path, output_dir, use_bluetopo=True):
    # Load CSB data
    csb_data = pd.read_csv(csb_data_path)
    gdf = gpd.GeoDataFrame(csb_data, geometry=gpd.points_from_xy(csb_data.lon, csb_data.lat))
    gdf = gdf.set_crs(4326, allow_override=True)


    # Create convex hull
    convex_hull_polygon = gdf.unary_union.convex_hull
    convex_hull_gdf = gpd.GeoDataFrame(geometry=[convex_hull_polygon], crs=gdf.crs)

    # Export to shapefile
    convex_hull_shapefile = os.path.join(output_dir, "convex_hull_polygon.shp")
    convex_hull_gdf.to_file(convex_hull_shapefile)

    if use_bluetopo:
        # Define a dedicated directory for BlueTopo tiles with the project title
        bluetopo_tiles_dir = os.path.join(output_dir, title, "Modeling")

        # Ensure the directory exists
        os.makedirs(bluetopo_tiles_dir, exist_ok=True)

        # Download BlueTopo tiles to the dedicated directory
        from nbs.bluetopo import fetch_tiles
        fetch_tiles(bluetopo_tiles_dir, convex_hull_shapefile, data_source='modeling')

        # Mosaic downloaded tiles using the new directory
        return mosaic_tiles(bluetopo_tiles_dir) 

def mosaic_tiles(tiles_dir):
    # Directory where the tiles are stored is now passed as a parameter
    mosaic_raster_path = os.path.join(tiles_dir, 'merged_tiles.tif')

    # List to store the opened rasters
    raster_list = []

    # Recursively search for geotiff files within the tiles_dir and open them
    for root, dirs, files in os.walk(tiles_dir):
        for file in files:
            if file.endswith('.tiff'):
                raster_path = os.path.join(root, file)
                src = rasterio.open(raster_path)
                raster_list.append(src)

    # Merge the rasters using only the first band
    merged_raster, out_transform = merge(raster_list)

    # Write the merged raster
    out_meta = src.meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": merged_raster.shape[1],
        "width": merged_raster.shape[2],
        "transform": out_transform,
        "count": 1  # Only one band
    })

    with rasterio.open(mosaic_raster_path, "w", **out_meta) as dest:
        dest.write(merged_raster[0, :, :], 1)  # Write only the first band

    # Close the opened rasters
    for src in raster_list:
        src.close()

    return mosaic_raster_path



def fetch_tide_data(station_id, start_date, end_date, interval=None, attempt_great_lakes=False):
    base_url = "https://tidesandcurrents.noaa.gov/api/datagetter"
    params = {
        "begin_date": start_date,
        "end_date": end_date,
        "station": station_id,
        "product": "predictions",
        "datum": "MLLW",
        "time_zone": "gmt",
        "units": "metric",
        "application": "NOAA_Coast_Survey",
        "format": "json"
    }

    if interval == 'hilo':
        params["interval"] = interval
    elif attempt_great_lakes:
        # Parameters for Great Lakes water level data
        params.update({
            "product": "water_level",
            "datum": "LWD"
        })

    # Construct the URL for debugging
    request_url = requests.Request('GET', base_url, params=params).prepare().url
    print(f"Requesting URL: {request_url}")  # Print the URL

    # Make the request
    response = requests.get(request_url)
    data = response.json()

    if 'predictions' in data:
        df = pd.json_normalize(data['predictions'])
        if 't' in df.columns and 'v' in df.columns:
            df['v'] = df['v'].astype(float)
            df['t'] = pd.to_datetime(df['t'])
            return df
    elif 'data' in data:
        df = pd.json_normalize(data['data'])
        if 't' in df.columns and 'v' in df.columns:
            df['v'] = pd.to_numeric(df['v'], errors='coerce')
            df['t'] = pd.to_datetime(df['t'], errors='coerce')
            return df
        else:
            print(f"Warning: Missing 't' or 'v' column in response for station {station_id}.")
            return pd.DataFrame()
    else:
        print(f"Warning: No 'predictions' in response for station {station_id}.")
        return pd.DataFrame()


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

        polygons = [shape(geom) for geom, val in shapes(binary_mask, mask=binary_mask, transform=transform) if val == 1]

        # Perform unary union
        unified_geometry = unary_union(polygons)

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
    join = join.astype({'time':'datetime64'})

    # Sort the join DataFrame by 'time' column
    join = join.sort_values('time')

    ts = join.groupby('ControlStn').agg(['min','max'])
    ts = ts['time']
    ts = ts.reset_index(drop=False)
    
    # Filter out non-numeric ControlStn values
    ts = ts[ts['ControlStn'].str.match(r'^\d+$')]  # Keep only rows where ControlStn is numeric

    # Format timestamps for NOAA CO-OPS Tidal Data API
    ts['min'] = ts['min'].dt.strftime("%Y%m%d %H:%M")
    ts['max'] = ts['max'].dt.strftime("%Y%m%d %H:%M")
    ts = ts.astype({'min': 'datetime64'})
    ts = ts.astype({'max': 'datetime64'})
    ts = ts.astype({'ControlStn':'int32'})

    ts.rename(columns = {'min':'StartDate', 'max':'EndDate'}, inplace=True)
    print(ts)
    ts = ts.set_index('StartDate')
    new_df = pd.DataFrame()
    for i, data in ts.iterrows():
        data = data.to_frame().transpose()
        data = data.reindex(pd.date_range(start=data.index[0], end=data.EndDate[0])).fillna(method='ffill').reset_index().rename(columns={'index': 'StartDate'})
        data = data.set_index('StartDate')
        data = data.resample('MS').bfill()
        data = data.reset_index()
        data['EndDate'] = data['StartDate']+pd.Timedelta(days=31)   
        new_df = pd.concat([new_df, data])
    
    new_df = new_df[['ControlStn','StartDate', 'EndDate']]

    new_df.rename(columns = {'StartDate':'min','EndDate': 'max'}, inplace=True)
    new_df['min'] = new_df['min'].dt.strftime("%Y%m%d %H:%M")
    new_df['max'] = new_df['max'].dt.strftime("%Y%m%d %H:%M")
    new_df = new_df.reset_index()
    new_df = new_df.astype({'ControlStn':'str'})
    print('*****retrieving tide data from NOAA COOPS API*****')
    
    tdf = []
    known_subordinate_stations = set()
    known_great_lakes_stations = set()

    for ind in new_df.index:
        station_id = new_df['ControlStn'][ind]

        # Directly fetch water level data for known Great Lakes stations
        if station_id in known_great_lakes_stations:
            waterlevel_data = fetch_tide_data(station_id, new_df['min'][ind], new_df['max'][ind], attempt_great_lakes=True)
            if not waterlevel_data.empty:
                tdf.append(waterlevel_data)
            else:
                print(f"No water level data available for Great Lakes station {station_id}.")
            continue  # Skip to next station

        # Attempt to fetch 6-minute predictions for non-subordinate stations
        if station_id not in known_subordinate_stations:
            min_predictions = fetch_tide_data(station_id, new_df['min'][ind], new_df['max'][ind])
            if not min_predictions.empty:
                tdf.append(min_predictions)
                continue

        # Fetch high-low data and interpolate for subordinate stations
        hilo_predictions = fetch_tide_data(station_id, new_df['min'][ind], new_df['max'][ind], interval='hilo')
        if not hilo_predictions.empty:
            interpolated_hilo = cosine_interpolation(hilo_predictions, new_df['min'][ind], new_df['max'][ind])
            tdf.append(interpolated_hilo)
            known_subordinate_stations.add(station_id)
            continue

        # Attempt to fetch Great Lakes water level data if other data not available
        waterlevel_data = fetch_tide_data(station_id, new_df['min'][ind], new_df['max'][ind], attempt_great_lakes=True)
        if not waterlevel_data.empty:
            tdf.append(waterlevel_data)
            known_great_lakes_stations.add(station_id)
        else:
            print(f"No water level data available for station {station_id}.")
 
    if tdf:
        tdf = pd.concat(tdf)
        print("Concatenated tdf shape:", tdf.shape)
        
        tdf = tdf.sort_values('t')
        jtdf = pd.merge_asof(join, tdf, left_on='time', right_on='t')
        print("jtdf shape before column drop:", jtdf.shape)
    
        # Dropping unneeded columns
        columns_to_drop = ['Shape__Are', 'Shape__Len', 'Input_FID', 'id', 'name', 'state', 'affil', 
                           'latitude', 'longitude', 'data', 'metaapi', 'dataapi', 'Shape_Le_2']
        jtdf.drop(columns=columns_to_drop, inplace=True, errors='ignore')
        print("jtdf shape after column drop:", jtdf.shape)
    
        # Check NaN values after dropping columns
        #print("NaN values in jtdf:", jtdf.isna().sum())
    
        # Drop rows based on NaNs in critical columns
        jtdf = jtdf.dropna(subset=['depth', 'time', 'geometry'])  # Assuming these are critical
        print("jtdf shape after dropna:", jtdf.shape)
        
        jtdf['t_corr'] = jtdf['t'] + pd.to_timedelta(jtdf['ATCorr'], unit='m')

        newdf = jtdf[['t_corr','v']].copy()
        print("newdf shape before dropna:", newdf.shape)

        newdf = newdf.rename(columns={'v':'v_new', 't_corr':'t_new'})
        newdf = newdf.sort_values('t_new').dropna()
        print("newdf shape after dropna:", newdf.shape)

        csb_corr = pd.merge_asof(jtdf, newdf, left_on='time', right_on='t_new', direction='nearest')
        print("csb_corr shape before dropna:", csb_corr.shape)

        #csb_corr = csb_corr.dropna()
        print("csb_corr shape after dropna:", csb_corr.shape)

        csb_corr['depth_new'] = csb_corr['depth'] - (csb_corr['RR'] * csb_corr['v_new'])   
        print("csb_corr shape after applying tide corrections:", csb_corr.shape)
             
        csb_corr = gpd.GeoDataFrame(csb_corr, geometry='geometry', crs='EPSG:4326')
        csb_corr['time'] = csb_corr['time'].dt.strftime("%Y%m%d %H:%M:%S")
        
        # Filter out depths
        csb_corr = csb_corr[(csb_corr['depth'] > 1.5) & (csb_corr['depth'] < 1000)]
        csb_corr = csb_corr.rename(columns={'depth':'depth_old'}).drop(columns=['index_right', 'ATCorr', 'RR', 'ATCorr2', 'RR2', 'Shape_Leng', 'Shape_Area', 'Shape_Le_1', 't', 'v', 't_corr', 't_new', 'v_new'])
        #print(csb_corr)
        return csb_corr
    else:
        print("No tide data available for the specified period.")
        return pd.DataFrame()


def BAGextract():
    print("DEBUG - BAG_filepath:", BAG_filepath)
    print('*****Starting to import BAG bathy and aggregate to 8m geotiff*****')
    #file_extension = os.path.splitext(BAG_filepath)[1].lower()

    # Translate options and creation of single-band GeoTIFF
    translate_options = gdal.TranslateOptions(bandList=[1], creationOptions=['COMPRESS=LZW'])
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
        gdal.Warp(reprojected_raster_path, intermediate_raster_path, dstSRS='EPSG:3395')
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
        gdal.Warp(output_raster_resampled, intermediate_raster_path, xRes=8, yRes=8)
        output_raster = output_raster_resampled
        # Clean up intermediate file
        os.remove(intermediate_raster_path)

    # Replace NaN values and update nodata value in the raster
    with rasterio.open(output_raster) as src:
        data = src.read(1)
        meta = src.meta
    
    data = np.where(np.isnan(data), 1000000, data)
    meta.update(nodata=1000000)
    
    with rasterio.open(output_raster, 'w', **meta) as dst:
        dst.write(data, 1)
    
    # Reproject the raster to WGS84
    output_raster_wgs84 = output_dir + '/' + title + '_wgs84.tif'
    gdal.Warp(output_raster_wgs84, output_raster, dstSRS='EPSG:4326')
    
    # Call create_survey_outline to generate the bathymetry polygon shapefile
    bathy_polygon_shp = create_survey_outline(output_raster_wgs84, output_dir, title)
    
    # Clean up intermediate files
    try:
        os.remove(intermediate_raster_path)
        intermediate_raster_path = output_dir + '/' + title + '_intermediate.tif'
        os.remove(intermediate_raster_path)
        os.remove(output_raster_resampled)
        os.remove(output_raster_wgs84)
    except Exception:
        pass
    
    return output_raster_wgs84, bathy_polygon_shp



def get_raster_value(x, y, raster):
    """ Get the value of the raster at the specified coordinates """
    ds = gdal.Open(raster)
    gt = ds.GetGeoTransform()
    rb = ds.GetRasterBand(1)
    nodata = rb.GetNoDataValue()

    px = int((x - gt[0]) / gt[1])
    py = int((y - gt[3]) / gt[5])

    if px < 0 or py < 0 or px >= ds.RasterXSize or py >= ds.RasterYSize:
        return np.nan  # Point is outside the raster

    value = rb.ReadAsArray(px, py, 1, 1)[0][0]

    if value == nodata:
        return np.nan  # Handle nodata values
    #print(value)
    return value

def derive_draft():
    output_raster, raster_boundary_shp = BAGextract()
    csb_corr = tides()
    
    print('*****Starting derive_draft function*****')
    
    # Load the raster boundary as a GeoDataFrame
    raster_boundary = gpd.read_file(raster_boundary_shp)

    # Ensure all geometries are valid
    raster_boundary['geometry'] = raster_boundary['geometry'].apply(lambda geom: geom if geom.is_valid else geom.buffer(0))
    csb_corr['geometry'] = csb_corr['geometry'].apply(lambda geom: geom if geom.is_valid else geom.buffer(0))

    # Assign a unique row identifier to the original dataframe
    csb_corr['row_id'] = range(len(csb_corr))
    
    
    # Perform spatial join to create a subset
    csb_corr_subset = csb_corr[csb_corr.geometry.within(raster_boundary.geometry.unary_union)]
    print("csb_corr_subset is...")
    #print(csb_corr_subset)

    # Sample from the subset
    max_sample_size = 300  # Maximum number of samples per unique_id
    sampled_csb_corr_subset = pd.DataFrame()
    for name, group in csb_corr_subset.groupby('unique_id'):
        if len(group) > max_sample_size:
            sampled_group = group.sample(n=max_sample_size)
        else:
            sampled_group = group
        sampled_csb_corr_subset = pd.concat([sampled_csb_corr_subset, sampled_group])
    print("sampled_csb_corr_subset shape after sampling:", sampled_csb_corr_subset.shape)
    print('Performing the raster depth extraction with the sampled subset')
    #print("sampled_csb_corr_subset is...")
    #print(sampled_csb_corr_subset)

    # Extract raster values for each point in the sampled subset
    sampled_csb_corr_subset['Raster_Value'] = sampled_csb_corr_subset.apply(
        lambda row: get_raster_value(row.geometry.x, row.geometry.y, output_raster), axis=1
    )

    # Merge the results back into the original dataframe
    csb_corr = csb_corr.merge(sampled_csb_corr_subset[['row_id', 'Raster_Value']], on='row_id', how='left')


    # Merge the results back into the original dataframe
    #csb_corr = csb_corr.merge(csb_corr_subset[['row_id', 'Raster_Value']], on='row_id', how='left')

    # Calculate the depth difference with the merged data
    print('*****Raster depth extraction complete. Calculating depth difference statistics.*****')
    csb_corr['diff'] = csb_corr['depth_new'] - (csb_corr['Raster_Value'] * -1)

    # Group and calculate statistics efficiently
    out = csb_corr.groupby('unique_id')['diff'].agg(['mean', 'std', 'count']).reset_index()

    print(out)
    out.to_csv(output_dir + '/VESSEL_OFFSETS_csb_corr_' + title + '.csv', mode='a')

    return csb_corr

def draft_corr():    
    #join vessel static offset table with csb data (1m std threshold for defining vertical transducer offset)
    #need to add in code to create and look for draft values in a master data-derived transducer offset table
    csb_corr = derive_draft()
    df = pd.read_csv(output_dir + "/VESSEL_OFFSETS_csb_corr_"+ title +'.csv')
    #df_master = pd.read_csv('master_vessel_offset_file.csv')
    df['std'] = pd.to_numeric(df['std'], errors='coerce').astype('float')
    df['mean'] = pd.to_numeric(df['mean'], errors='coerce').astype('float')
    csb_corr1 = pd.merge(csb_corr, df[df['std'] < 5.0], on='unique_id')
    csb_corr1['depth_final'] = ((csb_corr1['depth_new'] - csb_corr1['mean'])*-1) #correct csb bathy with transducer offset
    
    print('*****exporting tide and vessel offset corrected csb to shapefile*****')
    csb_corr1.to_file(output_dir + '/csb_OFFSETS_APPLIED_'+ title +'.shp', driver='ESRI Shapefile')
    return(csb_corr1)

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
            ((geom, value) for geom, value in zip(geodataframe.geometry, geodataframe['depth_final'])),
            out_shape=(y_res, x_res),
            transform=transform,
            fill=nodatavalue
        ), 1)
        
def rasterize_CSB():
    csb_corr1 = draft_corr()
    print('*****Rasterizing CSB data and exporting geotiff*****')
    output_raster_path = output_dir + '/csb_OFFSETS_APPLIED_' + title + '.tif'
    rasterize_with_rasterio(csb_corr1, output_raster_path, resolution)    
	
	# Define the path for the PNG plot
    plot_png_path = output_dir + '/csb_plot_' + title + '.png'

    # Create and save the plot as a PNG file
    with rasterio.open(output_raster_path) as src:
        fig, ax = plt.subplots(figsize=(12, 12))  # 1200x1200 pixels
        show(src, ax=ax, title= title + 'Quick-Look CSB Raster')
        plt.savefig(plot_png_path, dpi=200)  # Save plot as PNG
        plt.close(fig)  # Close the figure

    # Open the PNG file using the default image viewer
    if os.name == 'nt':  # Windows
        os.system(f'start {plot_png_path}')
    elif os.name == 'posix':  # macOS, Linux
        os.system(f'xdg-open {plot_png_path}')

    # Open explorer window of CSB processing results
    mod_output_dir = output_dir.replace("/", "\\")
    subprocess.Popen(r'explorer "'+ mod_output_dir)

def open_file_dialog(var, file_types):
    filename = filedialog.askopenfilename(filetypes=file_types)
    var.set(filename)

def open_folder_dialog(var):
    foldername = filedialog.askdirectory()
    var.set(foldername)
	
def process_csb_threaded():
    # Create a thread for process_csb function
    processing_thread = threading.Thread(target=process_csb)
    processing_thread.start()
	
def process_csb():
    start_time = time.time()  # Start the timer
    print("Processing...")
    global title, csb, BAG_filepath, fp_zones, output_dir

    # Update global variables
    title = title_var.get()
    csb = csb_var.get()
    fp_zones = fp_zones_var.get()
    output_dir = output_dir_var.get()
    
    # Update BAG_filepath with the selected file if BlueTopo is not used
    if not bluetopo_var.get():
        BAG_filepath = BAG_filepath_var.get()

    if bluetopo_var.get():
        # Update BAG_filepath with the path of the mosaicked raster
        BAG_filepath = create_convex_hull_and_download_tiles(csb, output_dir)
        rasterize_CSB()
    else:
        # If BlueTopo is not selected, BAG_filepath remains as set by the user
        rasterize_CSB()
    os.remove(output_dir + '/' + title + '_wgs84.tif')
    end_time = time.time()  # End the timer
    duration = end_time - start_time
    minutes, seconds = divmod(duration, 60)
    print(f"***** DONE! Thanks for all the CSB! *****\nTotal processing time: {int(minutes)} minutes and {seconds:.1f} seconds")


# Tkinter GUI setup
root = tk.Tk()
root.title("CSB Processing Input Files")

# Global variables as StringVar
title_var = tk.StringVar()
#output_crs_var = tk.StringVar()
csb_var = tk.StringVar()
BAG_filepath_var = tk.StringVar()
fp_zones_var = tk.StringVar()
output_dir_var = tk.StringVar()

# Layout and widgets
tk.Label(root, text='1. Title of CSB Data (do not use spaces between words)').grid(row=0, column=0, sticky='w')
title_entry = tk.Entry(root, textvariable=title_var)
title_entry.grid(row=0, column=1)

title_entry = tk.Entry(root, textvariable=title_var)

#button if you want to specify EPSG code for output geotiff instead of using hardcoded 3395
#tk.Label(root, text='EPSG Code for projected output geotiff CRS').grid(row=1, column=0, sticky='w')
#output_crs_entry = tk.Entry(root, textvariable = output_crs_var)
#output_crs_entry.grid(row=1, column=1)

tk.Label(root, text='2. Raw CSB data in *.csv format').grid(row=1, column=0, sticky='w')
csb_entry = tk.Entry(root, textvariable= csb_var)
csb_entry.grid(row=1, column=1)
tk.Button(root, text='Browse', command=lambda: open_file_dialog(csb_var, [("CSV file", "*.csv")])).grid(row=1, column=2)

tk.Label(root, text='3b. Input BAG or GeoTiff file for comparison bathymetry').grid(row=3, column=0, sticky='w')
BAG_filepath_entry = tk.Entry(root, textvariable= BAG_filepath_var)
BAG_filepath_entry.grid(row=3, column=1)
tk.Button(root, text='Browse', command=lambda: open_file_dialog(BAG_filepath_var, [("BAG or GeoTIFF file", "*.bag;*.tif;*.tiff")])).grid(row=3, column=2)

tk.Label(root, text='4. Tide Zone file in *.shp format').grid(row=4, column=0, sticky='w')
fp_zones_entry = tk.Entry(root, textvariable= fp_zones_var)
fp_zones_entry.grid(row=4, column=1)
tk.Button(root, text='Browse', command=lambda: open_file_dialog(fp_zones_var, [("Shapefile", "*.shp")])).grid(row=4, column=2)

tk.Label(root, text='5. Specify output folder').grid(row=5, column=0, sticky='w')
output_dir_entry = tk.Entry(root, textvariable=output_dir_var)
output_dir_entry.grid(row=5, column=1)
tk.Button(root, text='Browse', command=lambda: open_folder_dialog(output_dir_var)).grid(row=5, column=2)


tk.Button(root, text='Process', command=process_csb_threaded).grid(row=6, column=1, sticky='e')

bluetopo_var = tk.BooleanVar()
bluetopo_checkbox = tk.Checkbutton(root, text="3a. Use Automated BlueTopo Download... or,", variable=bluetopo_var)
bluetopo_checkbox.grid(row=2, column=0, columnspan=2, sticky='w')

def on_bluetopo_check():
    if bluetopo_var.get():
        BAG_filepath_entry.config(state='disabled')
    else:
        BAG_filepath_entry.config(state='normal')

bluetopo_checkbox.config(command=on_bluetopo_check)

if __name__ == "__main__":
    root.mainloop()
