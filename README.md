# Crowdsourced_Bathy_Processing
Processing crowdsourced bathymetry from International Hydrographic Office's Data Centre for Digital Bathymetry (CSV data processing). Corrects for water levels/tides and data-driven/postulated vessel transducer offset values. 

Basic steps are: 
1. Read in CSV file and perform basic outlier filtering (also filters out Anonymous vessels)
2. Spatially join CSB points with tide zone polygons
3. Download waterlevel/tide data from tide control stations for time period of CSB
4. Perform waterlevel correction to MLLW
5. Extract 5m geotiff from best available NOAA BAG bathymetry in the area
6. Extract raster values to CSB points, perform aggregated summary stats of difference between CSB depth and BAG bathymetry.              
7. Use this table choose unapplied static vertical transducer offsets of each vessel (any mean diff that has a std dev < 1.0m)
8. Apply the data-derived static offsets to the CSB data by performing an inner join with the offset table and CSB points
9. Export CSB points to shapefile for further bathy processing in ArcGIS Pro (IDW gridding, etc...) -***Would like to integrate gridding into python script***

Images: Comparison of CSB grid (IDW algorithm) vs Hydrographic Survey H13387 in Houston, TX Harbor

![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_8-46-32.png?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_8-47-30.png?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_9-10-00.png?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_9-20-11.png?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_9-13-34.png?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_16-18-02.png?raw=true)
