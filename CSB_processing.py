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
from osgeo import gdal
import subprocess
import rasterio
from rasterio.features import shapes, rasterize
from shapely.geometry import shape
from scipy.ndimage import binary_dilation
from rasterio.transform import from_origin
import matplotlib.pyplot as plt
from rasterio.plot import show
import time



# Global variables
title = ""
#output_crs = ""
csb = ""
BAG_filepath = ""
fp_zones = ""
output_dir = ""
resolution = 50 # resolution of quick-look geotiff raster

pd.set_option('display.max_columns', None)

def loadCSB():

    print('*****Reading CSB input csv file***** ')
    df1=gpd.read_file(csb) #read in CSB data in CSV
    #df1=df[df['platform_name'] != "Anonymous"] #filter out "Anonymous" vessels (can't do a vessel by vessel offset comparison against valid data)
    #Data cleaning = make sure dates are correct
    
    df1 = df1.astype({'depth' : 'float'}) #turn depth field numeric (float64)
    df2 = df1[df1['depth'] > .5] #remove depth values less that 1.2m (filtering out probable noise)
    df2 = df2[df2['depth'] < 1000] #remove depths values greater than possibleNYC depths in the area - this value can change based on area)
    #data cleaning - remove records with erroneous datetimes
    df2 = df2.astype({'time':'datetime64'}, errors='ignore')
    df2 = df2.dropna(subset=['time'])
    df2 = df2[df2['time'] > '2014']
    df2 = df2[df2['time'] < '2025']
    
    df2 = df2.drop_duplicates(subset=['lon', 'lat', 'depth', 'time', 'unique_id'])
    gdf = gpd.GeoDataFrame(df2, geometry=gpd.points_from_xy(df2.lon, df2.lat)) #create vector point file from lat and lon
    gdf = gdf.set_crs(4326, allow_override=True)
    return gdf
    
def tides():      
    gdf = loadCSB()
    print('CSB data from csv file loaded. Starting tide correction')

    zones=gpd.read_file(fp_zones)
    join = gpd.sjoin(gdf, zones, how='inner', predicate='within')
    join = join.astype({'time':'datetime64'})
    #print(join)
    
    ts = join.groupby('ControlStn').agg(['min','max'])
    ts = ts['time']
    ts = ts.reset_index(drop=False)
    
    
    #format timestamps for NOAA CO-OPS Tidal Data API
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
    for ind in new_df.index:
        try:
        
            URL_API = 'https://tidesandcurrents.noaa.gov/api/datagetter?begin_date='+new_df['min'][ind]+'&end_date='+new_df['max'][ind]+'&station='+new_df['ControlStn'][ind]+'&product=predictions&datum=mllw&units=metric&time_zone=gmt&application=NOAA_Coast_Survey&format=json'
            print(URL_API)
            response = requests.get(URL_API)
            json_dict = response.json()

        except Exception:
                continue
        try:
            data = pd.json_normalize(json_dict["predictions"])
            data[['v']] = data[['v']].apply(pd.to_numeric)
            data = data.astype({'t':'datetime64'})
            tdf.append(data)
        except Exception:
            continue
    try:
        tdf = pd.concat(tdf)
    except Exception:
        pass
    tdf = tdf.sort_values('t')
    join = join.sort_values('time')
          
    jtdf = pd.merge_asof(join, tdf, left_on='time', right_on='t')
    jtdf = jtdf.drop(columns=['ControlS_1'])
    jtdf = jtdf.drop(columns=['ControlS_2'])
    jtdf = jtdf.drop(columns=['DataProv'])
    jtdf = jtdf.dropna()
    jtdf['t_corr'] = jtdf['t'] + pd.to_timedelta(jtdf['ATCorr'], unit='m')
    
    
    
    newdf = jtdf[['t_corr','v']].copy()
    newdf = newdf.rename(columns={'v':'v_new', 't_corr':'t_new'})
    newdf = newdf.sort_values('t_new')
    newdf = newdf.dropna() 
    
    csb_corr = pd.merge_asof(jtdf, newdf, left_on='time', right_on='t_new', direction='nearest')
    csb_corr = csb_corr.dropna()
    
    csb_corr['depth_new'] = csb_corr['depth'] - (csb_corr['RR'] * csb_corr['v_new'])                
    csb_corr = gpd.GeoDataFrame(csb_corr, geometry='geometry', crs='EPSG:4326')
    csb_corr['time'] = csb_corr['time'].dt.strftime("%Y%m%d %H:%M:%S")
    
    #filter out depths less than 1.5m and greater than 1000m
    csb_corr = csb_corr[csb_corr['depth'] > 1.5]
    csb_corr = csb_corr[csb_corr['depth'] < 1000]
    csb_corr = csb_corr.rename(columns={'depth':'depth_old'})
    csb_corr = csb_corr.drop(columns=['index_right','ATCorr','RR','ATCorr2','RR2','Shape_Leng','Shape_Area','Shape_Le_1','t','v','t_corr','t_new','v_new'])
    return csb_corr

def BAGextract():
    print("DEBUG - BAG_filepath:", BAG_filepath)  # Debug print
    print('*****Starting to import BAG bathy and aggregate to 8m geotiff*****')
    dataset = gdal.Open(BAG_filepath, gdal.GA_ReadOnly)
    if dataset is None:
        print("Error: Unable to open the BAG file. Check the file path.")
        return None, None  # Return early to avoid further processing

    
    # Specify the options for gdal.Translate
    translate_options = gdal.TranslateOptions(bandList=[1],  # Use only the first band
                                              creationOptions=['COMPRESS=LZW'])  # Apply LZW compression

    # Create a single-band GeoTIFF with LZW compression
    gdal.Translate(output_dir + '/' + title + '.tif', dataset, options=translate_options)

    dataset = None

    dsReprj = gdal.Warp(output_dir + '/' + title + '_5m_MLLW.tif', output_dir + '/' + title + '.tif', xRes=8, yRes=8)  # resamples to chosen raster resolution
    dsReprj = None

    str1 = "gdal_translate -b 1 -of AAIGrid " + output_dir + '/' + title + '_5m_MLLW.tif ' + output_dir + '/' + title + '_MB_5m_elevation.tif'
    os.system(str1)  # extracts the elevation band and saves as new geotiff

    raster_file = output_dir + '/' + title + '_MB_5m_elevation.tif'
    input_raster = gdal.Open(raster_file)
    output_raster = output_dir + '/' + title + '_MB_5m_elevation_wgs84.tif'
    warp = gdal.Warp(output_raster, input_raster, dstSRS='EPSG:4326')
    raster_file = None
    input_raster = None
    warp = None

    os.remove(output_dir + '/' + title + '_5m_MLLW.tif')
    os.remove(output_dir + '/' + title + '.tif')
    os.remove(output_dir + '/' + title + '_MB_5m_elevation.tif')

    output_raster = output_dir + '/' + title + '_MB_5m_elevation_wgs84.tif'

    with rasterio.open(output_raster) as src:
        array = src.read(1)
        transform = src.transform
        nodata = src.nodatavals[0]

        # Create a binary array: 1 where there's data, 0 where there's no data
        binary_array = np.where(array == nodata, 0, 1).astype('int16')

        # Dilate the binary array and convert to a compatible dtype
        dilated_array = binary_dilation(binary_array, iterations=25)
        dilated_array = dilated_array.astype('uint8')  # Convert to uint8

    # Generate shapes from the dilated binary array
    geometry = []
    for shp, val in shapes(dilated_array, transform=transform):
        if val == 1:  # Only consider shapes with value 1 (dilated areas)
            geom = shape(shp).simplify(tolerance=0.001)  # Apply simplification
            if geom.is_valid:
                geometry.append(geom)

    # Create a GeoDataFrame
    geo_df = gpd.GeoDataFrame(geometry=geometry, crs=src.crs)
    

    # Save the bathymetry polygon shapefile
    bathy_polygon_shp = output_dir + '/' + title + '_bathy_polygon.shp'
    geo_df.to_file(bathy_polygon_shp, driver='ESRI Shapefile')

    print('Bathymetry polygon shapefile created.')
    return output_raster, bathy_polygon_shp

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

    print('Performing the raster depth extraction with the subset')

    # Extract raster values for each point in the subset
    csb_corr_subset['Raster_Value'] = csb_corr_subset.apply(
        lambda row: get_raster_value(row.geometry.x, row.geometry.y, output_raster), axis=1
    )

    # Merge the results back into the original dataframe
    csb_corr = csb_corr.merge(csb_corr_subset[['row_id', 'Raster_Value']], on='row_id', how='left')

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

    with rasterio.open(output_raster_path) as src:
        fig, ax = plt.subplots()
        show(src, ax=ax, title='CSB Raster')
        plt.show()
        
    # Open explorer window of CSB processing results
    mod_output_dir = output_dir.replace("/", "\\")
    subprocess.Popen(r'explorer "'+ mod_output_dir)

def open_file_dialog(var, file_types):
    filename = filedialog.askopenfilename(filetypes=file_types)
    var.set(filename)

def open_folder_dialog(var):
    foldername = filedialog.askdirectory()
    var.set(foldername)

def process_csb():
    start_time = time.time()  # Start the timer
    print("Processing...")
    global title, output_crs, csb, BAG_filepath, fp_zones, output_dir

    # Update global variables
    title = title_var.get()
    #output_crs = output_crs_var.get()
    csb = csb_var.get()
    BAG_filepath = BAG_filepath_var.get()
    fp_zones = fp_zones_var.get()
    output_dir = output_dir_var.get()

    # Call the main processing function
    rasterize_CSB()
    
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
tk.Label(root, text='Title of CSB Data').grid(row=0, column=0, sticky='w')
title_entry = tk.Entry(root, textvariable=title_var)
title_entry.grid(row=0, column=1)

title_entry = tk.Entry(root, textvariable=title_var)

#button if you want to specify EPSG code for output geotiff instead of using hardcoded 3395
#tk.Label(root, text='EPSG Code for projected output geotiff CRS').grid(row=1, column=0, sticky='w')
#output_crs_entry = tk.Entry(root, textvariable = output_crs_var)
#output_crs_entry.grid(row=1, column=1)

tk.Label(root, text='Raw CSB data in *.csv format').grid(row=1, column=0, sticky='w')
csb_entry = tk.Entry(root, textvariable= csb_var)
csb_entry.grid(row=1, column=1)
tk.Button(root, text='Browse', command=lambda: open_file_dialog(csb_var, [("CSV file", "*.csv")])).grid(row=1, column=2)

tk.Label(root, text='Input BAG file for comparison bathymetry').grid(row=2, column=0, sticky='w')
BAG_filepath_entry = tk.Entry(root, textvariable= BAG_filepath_var)
BAG_filepath_entry.grid(row=2, column=1)
tk.Button(root, text='Browse', command=lambda: open_file_dialog(BAG_filepath_var, [("BAG file", "*.bag")])).grid(row=2, column=2)

tk.Label(root, text='Tide Zone file in *.shp format').grid(row=3, column=0, sticky='w')
fp_zones_entry = tk.Entry(root, textvariable= fp_zones_var)
fp_zones_entry.grid(row=3, column=1)
tk.Button(root, text='Browse', command=lambda: open_file_dialog(fp_zones_var, [("Shapefile", "*.shp")])).grid(row=3, column=2)

tk.Label(root, text='Specify output folder').grid(row=4, column=0, sticky='w')
output_dir_entry = tk.Entry(root, textvariable=output_dir_var)
output_dir_entry.grid(row=4, column=1)
tk.Button(root, text='Browse', command=lambda: open_folder_dialog(output_dir_var)).grid(row=4, column=2)


tk.Button(root, text='Process', command=process_csb).grid(row=5, column=1, sticky='e')

if __name__ == "__main__":
    root.mainloop()
