# Crowdsourced_Bathy_Processing
Processing crowdsourced bathymetry from International Hydrographic Office's Data Centre for Digital Bathymetry (CSV data processing). Corrects for water levels/tides and data-driven/postulated vessel transducer offset values. 

Basic steps are: 
1. Read in CSV file and perform basic outlier filtering (also filters out Anonymous vessels)
2. Spatially join CSB points with tide zone polygons
3. Download waterlevel/tide data from tide control stations for time period of CSB
4. Perform waterlevel correction to MLLW
5. Extract 5m geotiff from best available NOAA BAG bathymetry in the area
6. Extract raster values to CSB points, perform aggregated summary stats of difference between CSB depth and BAG bathymetry. (see example below)
                                platform mean_diff       std
                  0   CAPT JERRY J.WILTZ -2.698360  0.828494
                  1            CAPT PETE -2.364760  0.367794
                  2               CIBOLO -2.909903  0.470642
                  3                COOKE -2.467871  0.395678
                  4              GRAYSON -2.212950  0.466788
                  5      Genesis Patriot -4.604306  0.286585
                  6             JOE PYNE  4.094317  4.411822
                  7       JOHN T MCMAHAN -2.105020  0.517740
                  8             Joe Pyne  4.334097  4.221233
                  9       Kathryn Louise -2.618789  0.290383
                  10        MISS CYNTHIA -2.242290  0.566003
                  11       POINT MALLARD -2.349711  0.429325
                  12          REX DOBSON -2.086876  0.384315
                  13         ROY MICHAEL -2.499080  0.638724
                  14                STBL -0.666195  5.328128
                  15        THREE RIVERS -2.419527  0.329376
                  16            joe pyne  4.457161  4.869262
7. Use this table choose unapplied static vertical transducer offsets of each vessel (any mean diff that has a std dev < 1.0m)
8. Apply the data-derived static offsets to the CSB data by performing an inner join with the offset table and CSB points
9. Export CSB points to shapefile for further bathy processing in ArcGIS Pro (IDW gridding, etc...) -***Would like to integrate gridding into python script***
