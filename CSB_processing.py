# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 13:57:44 2020
Updated 2/6/2023

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
import requests
import sys
import os
from osgeo import gdal
from rasterstats import point_query
from geocube.api.core import make_geocube
import subprocess


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

pd.set_option('display.max_columns', None)
resolution = 50 #desired resolution of output raster in the units of your output_crs (meters or decimal degrees, most likely)

def loadCSB():
    print('*****Reading CSB input csv file***** ')
    df=gpd.read_file(csb) #read in CSB data in CSV
    df1=df[df['platform'] != "Anonymous"] #filter out "Anonymous" vessels (can't do a vessel by vessel offset comparison against valid data)
    #Data cleaning = make sure dates are correct
    
    df1 = df1.astype({'depth' : 'float'}) #turn depth field numeric (float64)
    df2 = df1[df1['depth'] > .5] #remove depth values less that 1.2m (filtering out probable noise)
    df2 = df2[df2['depth'] < 1000] #remove depths values greater than possible depths in the area - this value can change based on area)
    #data cleaning - remove records with erroneous datetimes
    df2 = df2.astype({'time':'datetime64'})
    df2 = df2[df2['time'] > '2014']
    df2 = df2[df2['time'] < '2023']
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
    csb_corr = csb_corr.drop(columns=['index_right','ATCorr','RR','ATCorr2','RR2','DataProv','Shape_Leng','Shape_Area','Shape_Le_1','t','v','t_corr','t_new','v_new'])
    return csb_corr

def BAGextract():
    print('*****Starting to import BAG bathy and aggregate to 5m geotiff*****')
    #how to convert BAG to geotiff and extract just the elevation layer (and resample if wanted)
    
    #should try and not create so many actual files (keep rasters in memory instead)
    dataset = gdal.Open(BAG_filepath, gdal.GA_ReadOnly)
    gdal.Translate(output_dir +'/'+ title+ '.tif', dataset) #converts BAG to dual band geotiff
    dataset = None
    dsReprj = gdal.Warp(output_dir + '/' + title + '_5m_MLLW.tif', output_dir + '/' + title + '.tif', xRes = 8, yRes=8) #resamples the input geotiff to chosen raster resolution
    str1 = "gdal_translate -b 1 -of AAIGrid " + output_dir + '/' + title +'_5m_MLLW.tif '+ output_dir + '/' + title + '_MB_5m_elevation.tif'
    #str1 = "gdal_translate -b 1 -of AAIGrid " + 'E:/CSB_13SEPT2021/' + title +'_5m_MLLW.tif E:/CSB_13SEPT2021/' + title + '_MB_5m_elevation.tif'
    dsReprj = None
    os.system(str1) #extracts the elevation band (band1) and saves as new geotiff
    raster_file = output_dir + '/' + title + '_MB_5m_elevation.tif'
    input_raster = gdal.Open(raster_file)
    output_raster = output_dir + '/' + title + '_MB_5m_elevation_wgs84.tif'
    warp = gdal.Warp(output_raster, input_raster, dstSRS='EPSG:4326')
    raster_file = None
    input_raster = None
    warp = None
    os.remove(output_dir + '/' + title + '_5m_MLLW.tif')
    os.remove(output_dir +'/'+ title+ '.tif')
    os.remove(output_dir + '/' + title + '_MB_5m_elevation.tif')
    print('****bathy file processed, 5m geotiff creation complete****')
    return output_raster

def derive_draft():
    output_raster = BAGextract()
    csb_corr = tides()
    print('*****Starting raster depth extraction at CSB point positions. This step can take a long time depending on how many CSB records are in your file*****')

    pts = point_query(csb_corr, output_raster) #this seems to take a long time
    print('*****raster depth extraction complete. Calculating and aggregating depth difference statistics. Aggregated results will be stored in the pandas dataframe named: out*****')
    #next, take resulting pts1 list and add as new column to pts gdf. 
    csb_corr["Raster_Value"] = pts
    csb_corr['diff'] = csb_corr['depth_new'] - (csb_corr['Raster_Value'] * -1) #create new column of diff between raster value and observed csb depth
    out=csb_corr.groupby('platform')['diff'].agg(['mean','std','count']).reset_index() #create new dataframe of vessels showing summary stats for diff
    print(out)
    #should append Master Vessel Offset Table instead of creating individual tables per area
    #how to update fields (mean diff and stdev) only where count is greater than existing record???
    out.to_csv(output_dir +'/VESSEL_OFFSETS_csb_corr_'+ title +'.csv', mode = 'a')
    return(csb_corr)

def draft_corr():    
    #join vessel static offset table with csb data (1m std threshold for defining vertical transducer offset)
    #need to add in code to create and look for draft values in a master data-derived transducer offset table
    csb_corr = derive_draft()
    df = pd.read_csv(output_dir + "/VESSEL_OFFSETS_csb_corr_"+ title +'.csv')
    #df_master = pd.read_csv('master_vessel_offset_file.csv')
    df['std'] = pd.to_numeric(df['std'], errors='coerce').astype('float')
    df['mean'] = pd.to_numeric(df['mean'], errors='coerce').astype('float')
    csb_corr1 = pd.merge(csb_corr, df[df['std'] < 5.0], on='platform')
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
    
    print('***** DONE! Thanks for all the CSB! *****')

if __name__ == "__main__":
    rasterize_CSB()
    