# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 13:57:44 2020
Updated 9/29/2022

This script automates the tide correction of crowdsourced bathymetry (CSB) files (in CSV format)
downloaded from the Internation Hydrographic Organization's (IHO) Data Centre for Digital Bathymetry 
Crowdsourced Bathymetry Database. 

The raw CSB data is compared to tide corrected, and then compared against known bathymetry (in BAG format).
The mean difference and standard deviation is tabulated for each contributor vessel, and any mean difference
that has a standard deviation less than 2m is used as a vertical static transducer offset value for that vessel.
Those vertical offsets are applied to the tide-corrected data to create a CSB bathymetry solution.

The output is a shapefile.


@author: Anthony Klemm
"""


import geopandas as gpd
import pandas as pd
import requests
import os
from osgeo import gdal
from rasterstats import point_query


pd.set_option('display.max_columns', None)
title = 'False_Pass_CSB'

fp_zones = r"F:\csb\tide zone polygons\tide_zone_polygons_new_WGS84.shp"
csb = r"E:\CSB_13SEPT2021\False_Pass.csv"
BAG_filepath = "E:/CSB_13SEPT2021/H12632_MBVB_4m_MLLW_Combined.bag"

print('*****Reading CSB input csv file***** ')
df=gpd.read_file(csb) #read in CSB data in CSV
df1=df[df['platform'] != "Anonymous"] #filter out "Anonymous" vessels (can't do a vessel by vessel offset comparison against valid data)
df1 = df1.astype({'depth' : 'float'}) #turn depth field numeric (float64)
df2=df1[df1['depth'] > .5] #remove depth values less that 1.2m (filtering out probable noise)
df2=df2[df2['depth'] < 1000] #remove depths values greater than possible depths in the area - this value can change based on area)
gdf = gpd.GeoDataFrame(df2, geometry=gpd.points_from_xy(df2.lon, df2.lat)) #create vector point file from lat and lon
gdf = gdf.set_crs(4326, allow_override=True)
#print(gdf)
print('CSV data csv file loaded. Starting tide correction')


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

ts.rename(columns = {'min':'StartDate', 'max':'EndDate'}, inplace=True)
print(ts)
ts = ts.set_index('StartDate')
new_df = pd.DataFrame()
for i, data in ts.iterrows():
    data = data.to_frame().transpose()
    data = data.reindex(pd.date_range(start=data.index[0], end=data.EndDate[0])).fillna(method='ffill').reset_index().rename(columns={'index': 'StartDate'})
    new_df = pd.concat([new_df, data])

new_df = new_df[['ControlStn','StartDate', 'EndDate']]
new_df['EndDate'] = new_df['StartDate']+pd.Timedelta(days=1)
new_df.rename(columns = {'StartDate':'min','EndDate': 'max'}, inplace=True)
new_df['min'] = new_df['min'].dt.strftime("%Y%m%d %H:%M")
new_df['max'] = new_df['max'].dt.strftime("%Y%m%d %H:%M")
new_df = new_df.reset_index()
print(new_df)
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
print(tdf)
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
print(csb_corr)


print('*****Starting to import BAG bathy and aggregate to 5m geotiff*****')
#how to convert BAG to geotiff and extract just the elevation layer (and resample if wanted)


dataset = gdal.Open(BAG_filepath, gdal.GA_ReadOnly)
gdal.Translate('E:/CSB_13SEPT2021/' + title+ '.tif', dataset) #converts BAG to dual band geotiff
dsReprj = gdal.Warp('E:/CSB_13SEPT2021/' + title+'_5m_MLLW.tif', 'E:/CSB_13SEPT2021/' + title + '.tif', xRes = 5, yRes=5) #resamples the input geotiff to 10m raster

str1 = "gdal_translate -b 1 -of AAIGrid " + 'E:/CSB_13SEPT2021/' + title +'_5m_MLLW.tif E:/CSB_13SEPT2021/' + title + '_MB_5m_elevation.tif'
os.system(str1) #extracts the elevation band (band1) and saves as new geotiff
raster_file = 'E:/CSB_13SEPT2021/' + title + '_MB_5m_elevation.tif'
input_raster = gdal.Open(raster_file)
output_raster = 'E:/CSB_13SEPT2021/' + title + '_MB_5m_elevation_wgs84.tif'
warp = gdal.Warp(output_raster, input_raster, dstSRS='EPSG:4326')
warp = None

print('*****bathy file processed, 5m geotiff creation complete. Starting raster depth extraction at CSB point positions. This step can take a long time depending on how many CSB records are in your file*****')

#how to extract raster values that intersect with the point locations and add new column to geodataframe
#need to make sure bag and points are in the same spatial/coordinate reference system/projection


pts = point_query(csb_corr, output_raster) #this seems to take a long time
print('*****raster depth extraction complete. Calculating and aggregating depth difference statistics. Aggregated results will be stored in the pandas dataframe named: out*****')
#next, take resulting pts1 list and add as new column to pts gdf. 
csb_corr["Raster_Value"] = pts
csb_corr['diff'] = csb_corr['depth_new'] - (csb_corr['Raster_Value'] * -1) #create new column of diff between raster value and observed csb depth
out=csb_corr.groupby('platform')['diff'].agg(['mean','std','count']).reset_index() #create new dataframe of vessels showing summary stats for diff
print(out)
out.to_csv('E:/CSB_13SEPT2021/VESSEL_OFFSETS_csb_corr_'+ title +'.csv')

#join vessel static offset table with csb data (1m std threshold for defining vertical transducer offset)
csb_corr1 = pd.merge(csb_corr, out[out['std'] < 2.0], on='platform')
csb_corr1['depth_final'] = csb_corr1['depth_new'] - csb_corr1['mean'] #correct csb bathy with transducer offset

print('*****exporting tide and vessel offset corrected csb to shapefile*****')
csb_corr1.to_file('F:/csb/results/csb_OFFSETS_APPLIED_'+ title +'.shp', driver='ESRI Shapefile')

print('***** DONE! Thanks for all the CSB! *****')
