# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 13:57:44 2020
Updated 11/13/2023

This script automates the tide correction of crowdsourced bathymetry (CSB) files (in CSV format)
downloaded from the International Hydrographic Organization's (IHO) Data Centre for Digital Bathymetry 
Crowdsourced Bathymetry Database. 

The raw CSB data is tide corrected using a discrete zoned tides model, and then compared against known 
bathymetry (in BAG format) in a common high-traffic area.
The mean difference and standard deviation is tabulated for each contributor vessel, and any mean difference
that has a standard deviation less than 2m is used as a vertical static transducer offset value for that vessel.
Those vertical offsets are applied to the tide-corrected data to create a CSB bathymetry solution.

The output is a shapefile and a basic geotiff without sophisticated interpolation (recommend changing to IDW). 
The next development step is to compile the script as a standalone executable. 

@author: Anthony Klemm
"""

import PySimpleGUI as sg
import geopandas as gpd
import pandas as pd
import numpy as np
import requests
import sys
import os
from osgeo import gdal, ogr, osr
from rasterstats import point_query
from geocube.api.core import make_geocube
import subprocess
import time
import rasterio
from rasterio.features import shapes
from shapely.geometry import shape
from shapely.ops import unary_union
from scipy.ndimage import binary_dilation



sg.theme('DarkTeal2')
if len(sys.argv) == 1:
    event, values = sg.Window('CSB Processing Input Files',
                    [[sg.Text('Title of CSB Data', size=(40,1)), sg.InputText()],
                    [sg.Text('EPSG Code for projected output geotiff CRS', size=(40,1)), sg.InputText()],
                    [sg.Text('Raw CSB data in *.csv format')],
                    [sg.In(), sg.FileBrowse(file_types=(('CSV file', '*.csv'),))],
                    [sg.Text('Input BAG file for comparison bathymetry')],
                    [sg.In(), sg.FileBrowse(file_types=(('BAG file', '*.bag'),))],
                    [sg.Text('Tide Zone file in *.shp format')],
                    [sg.In(), sg.FileBrowse(file_types=(('Shapefile', '*.shp'),))],
                    [sg.Text('Specify output folder')],
                    [sg.In(), sg.FolderBrowse()],
                    [sg.Open(), sg.Cancel()]]).read(close=True)
    
    title = values[0]
    output_crs = "epsg:"+values[1]
    csb = values[2]
    BAG_filepath = values[3]
    fp_zones = values[4]
    output_dir = values[5]
else:
    csb = sys.argv[1]

if not csb:
    sg.popup("Cancel", "No filename for csb raw file supplied")
    raise SystemExit("Cancelling: no filename supplied")
else:
    sg.popup('The CSB filename you chose was', csb)
    start_time = time.time()


pd.set_option('display.max_columns', None)
resolution = 50 #desired resolution of output raster in the units of your output_crs (meters or decimal degrees, most likely)

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
    df2 = df2[df2['time'] < '2024']
    
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
            #print(json_dict)
            #export out as individual json files of the reference tide station data
            #out_file = open("E:/csb/jsondump1/"+'8770733'+'.json', 'w')
            #json.dump(json_dict, out_file, indent = 6)
            #out_file.close()
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
    print('*****Starting to import BAG bathy and aggregate to 8m geotiff*****')
    # Converting BAG to geotiff and extracting the elevation layer
    dataset = gdal.Open(BAG_filepath, gdal.GA_ReadOnly)
    gdal.Translate(output_dir + '/' + title + '.tif', dataset)  # converts BAG to dual band geotiff
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

def derive_draft():
    output_raster, raster_boundary_shp = BAGextract()
    csb_corr = tides()

    print('*****Starting derive_draft function*****')
    
    
    # Load the raster boundary as a GeoDataFrame
    raster_boundary = gpd.read_file(raster_boundary_shp)

    # Ensure all geometries are valid
    raster_boundary['geometry'] = raster_boundary['geometry'].apply(lambda geom: geom if geom.is_valid else geom.buffer(0))
    csb_corr['geometry'] = csb_corr['geometry'].apply(lambda geom: geom if geom.is_valid else geom.buffer(0))

    # Remove any null geometries
    raster_boundary = raster_boundary[raster_boundary['geometry'].notnull()]
    #csb_corr = csb_corr[csb_corr['geometry'].notnull()]

    print('Performing the spatial join to create a subset')
    #csb_corr_subset = csb_corr[csb_corr.geometry.within(raster_boundary.geometry.unary_union)]


    # Assign a unique row identifier to the original dataframe
    csb_corr['row_id'] = range(len(csb_corr))

    # Perform spatial join to create a subset
    csb_corr_subset = csb_corr[csb_corr.geometry.within(raster_boundary.geometry.unary_union)]

    print('Performing the raster depth extraction with the subset')
    pts = point_query(csb_corr_subset, output_raster)

    # Add the results to the subset
    csb_corr_subset["Raster_Value"] = pts

    # Merge the results back into the original dataframe
    # We will use a left merge to ensure all original records in csb_corr are preserved
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
    
def rasterize_CSB():    
    #Rasterize the CSB results and export geotiff depth grid
    csb_corr1=draft_corr()
    print('*****rasterizing CSB data and exporting geotiff*****')

    geo_grid = make_geocube(vector_data=
        csb_corr1,
        measurements=["depth_final"],
        output_crs=output_crs,
        resolution=((resolution * -1), resolution),
        )
        
    geo_grid['depth_final'].plot()
    geo_grid["depth_final"].rio.to_raster(output_dir + '/csb_OFFSETS_APPLIED_'+ title + '.tif')
    mod_output_dir = output_dir.replace("/", "\\")
    subprocess.Popen(r'explorer "'+ mod_output_dir)
    end_time = time.time()
    execution_time = end_time - start_time
    


    print('***** DONE! Thanks for all the CSB! *****')
    print(f"Total execution time: {execution_time} seconds")

if __name__ == "__main__":
    rasterize_CSB()
    