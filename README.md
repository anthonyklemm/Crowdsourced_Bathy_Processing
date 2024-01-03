# Crowdsourced Bathymtery Processing
This tool processes crowdsourced bathymetry from International Hydrographic Organization's Data Centre for Digital Bathymetry (CSV data from their CSB point store). It corrects for water levels/tides
 and data-derived vessel transducer offset values. 

***Note: This is the first iteration of this tool, and is intended to provide users reconnaissance information from existing CSB data.***

***Contact Anthony Klemm anthony.r.klemm@noaa.gov for any questions regarding this tool*** 

To use this tool, you must download csb data from the IHO Data Centre for Digital Bathymetry https://www.ncei.noaa.gov/maps/iho_dcdb/.
You must also download a BAG file (or geotiff) for some comparison data in the same area where your CSB tracklines intersect. This is important for the data-derived transducer offset step.
The script accomodates comparison bathymetry BlueTopo tiles in geotiff format. Use the downloader tool manually in Pydro if interested, OR you can click the option button to have the script auto-download all the bluetopo model tiles automatically (THIS IS THE PREFERRED METHOD in areas where the bluetopo model is available).

The tide zone polygons are from NOAA COOPS and a copy of the polygons are located here: "C:\Pydro22\NOAA\site-packages\Python38\svn_repo\HSTB\CSB_processing\tide_zone_polygons.shp" 
You may find an area of interest with missing tide zone polygons. You can edit this shapefile to include new custom polygons (just understand that this may not be as accurate as the published COOPS zones). 
If there is a tide reference station nearby, you just need to update the ControlStn attribute of your custom polygon with the correlating tide reference station control number (found on the COOPS website),
as well as the ATCorr to 0 (which means there is a time correction of zero since it's nearby the tide station), and a RR of 1 (the magnitude coefficient). Also, a beta version of a tie zone polygon dataset is also available with auto-generated zones around subordinate tide prediction stations outside the official tide zone definitions provided by NOAA CO-OPS. **USE WITH CAUTION** They were generated using Voronoi diagram/Thiessen Plogyons algorithm from the point locations of all subordinate tide prediction stations, but were not edited to respect shoreline so there are errors. The script pulls the predicted highs and lows and interpolates between the values using cosine interpolation (which approximates tide changes but adds uncertainty to the estimated values). 

Basic steps of the script are: 

1. Read in CSV file and perform basic outlier filtering
2. Spatially join CSB points with tide zone polygons
3. Download waterlevel/tide data from tide control stations for time period of CSB from NOAA COOPS datagetter web API
4. Perform waterlevel correction to MLLW
5. Extract geotiff from best available NOAA BAG bathymetry in the area
6. Extract raster values to CSB points, perform aggregated summary stats of difference between CSB depth and BAG bathymetry.              
7. Use this table to choose unapplied static vertical transducer offsets of each vessel 
8. Apply the data-derived static offsets to the CSB data by performing an inner join with the offset table and CSB points
9. Export CSB points to shapefile for further bathy processing in ArcGIS Pro (IDW gridding, etc...) - The attribute with the final corrected depths is called "depth_fina" and is in WGS84 crs (unprojected).
10. Export out a quick-look geotiff (not using any interpolation) at EPSG: 3395 (World Mercator) and 50m resolution

Images: Comparison of CSB grid (IDW algorithm) vs Hydrographic Survey H13387 in Houston, TX Harbor

![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/csb_vs_BAG.gif?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_9-10-00.png?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/False_Pass.gif?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-10-05_22-33-19.png?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_8-46-32.png?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_8-47-30.png?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_16-13-34.png?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_9-20-11.png?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_9-13-34.png?raw=true)
![Screenshot](https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/images/2022-09-29_16-18-02.png?raw=true)
