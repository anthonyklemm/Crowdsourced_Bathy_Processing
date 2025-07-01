# Crowdsourced Bathymetry Processing
This tool processes crowdsourced bathymetry from International Hydrographic Organization's Data Centre for Digital Bathymetry (CSV data from their CSB point store). It corrects for water levels/tides
 and data-derived vessel transducer offset values. 
 
An SOP on how to use the tool can be found here (copy and paste this URL in your web browser):  https://github.com/anthonyklemm/Crowdsourced_Bathy_Processing/blob/main/CSB%20Processing%20SOP%20v1.pdf


***Note: This is the first iteration of this tool, and is intended to provide users reconnaissance information from existing CSB data.***

***Contact Anthony Klemm anthony.r.klemm@noaa.gov for any questions regarding this tool*** 

To use this tool, you must download csb data from the IHO Data Centre for Digital Bathymetry https://www.ncei.noaa.gov/maps/iho_dcdb/.
# Crowdsourced Bathymetry (CSB) Processing Pipeline

This application provides a full end-to-end workflow for processing raw Crowdsourced Bathymetry (CSB) data into analysis-ready geospatial products. It features a graphical user interface (GUI) to manage inputs, processing steps, and outputs, making complex geospatial operations accessible.

The pipeline ingests raw CSV files, applies tidal and vertical offset corrections, performs statistical outlier detection, and generates final products such as a consolidated points GeoPackage and gridded GeoTIFF rasters.

## Features

- **Tide Correction:** Applies tidal adjustments to depth soundings using data from the NOAA Tides and Currents API and a user-provided tide zone shapefile. location of supplied zoned tide model polygons shapefile:
C:\Pydro24\NOAA\site-packages\Python3\svn_repo\HSTB\CSB_processing\BETA_subordinate_tide_zones\tide_zone_polygons_new_WGS84_merge.shp
- **Vertical Offset Calibration:** Compares CSB data against a reference bathymetric surface to calculate and apply a vertical offset for each vessel, normalizing data from different sources.
- **Flexible Reference Surfaces:**
    - **Automated BlueTopo Download:** Automatically fetches the latest NOAA BlueTopo bathymetric tiles for the specific area covered by the input CSB data.
    - **Local File Support:** Allows the user to provide their own local reference bathymetry in BAG or GeoTIFF format.
- **Advanced Outlier Detection:** Employs a multi-pass statistical routine using `scikit-learn`'s `IterativeImputer` to flag potential outliers in each vessel's transit.
- **Database Backend:** Utilizes DuckDB for efficient intermediate data storage, allowing for complex and fast queries on the processed data.
- **Incremental Processing & Append-Only Workflow:**
    - The pipeline is designed for continuous data integration. If an existing csb.duckdb database and master_offsets.csv file are found in the output directory, the script will use them as a baseline.
    - New CSB data is processed and appended to the existing database.
    - It leverages existing vessel offsets from master_offsets.csv to avoid re-calculating them, improving consistency and speed.
    - Post-processing steps are applied only to the newly added data.
- **Flexible Final Products:**
    - **Intermediate Products (Optional):** Per-transit GeoPackages and plots.
    - **Final Products (Optional):**
        - A consolidated points GeoPackage (`csb_final_points.gpkg`) in EPSG:4326.
        - A gridded GeoTIFF (`csb_final_gridded.tif`) created by averaging non-outlier points.
- **Tiled Processing:** Supports providing a tessellation shapefile to process large areas into smaller, tiled GeoTIFFs and GeoPackages.
- **VRT Creation (Optional):** Automatically organizes tiled GeoTIFFs by their coordinate system and builds a VRT (Virtual Raster) with overviews for each group, making them easy to use in GIS software.
- **GUI:** A user-friendly graphical interface built with Tkinter to manage all inputs and processing options.
- **Cleanup:** Automatically removes temporary files after processing.

## How to Use

1.  **Prepare Your Input Data:**
    - **Raw CSB Directory:** A folder containing one or more CSB data files in `.csv` format.
    - **Tide Zone Shapefile:** A shapefile that defines the geographic zones for tide correction. It must contain the columns `ControlStn`, `ATCorr`, and `RR`.
    - **Output Directory:** An empty folder where all output products will be saved.
    - **(Optional) Reference Bathymetry:** A local BAG or GeoTIFF file if you are not using the automated BlueTopo download.
    - **(Optional) Tessellation Shapefile:** A polygon shapefile to use for tiled processing.

2.  **Launch the Application:**
    Run the script from your terminal:
    ```
    python csb_processing.py
    ```

3.  **Configure the Processing Run via the GUI:**

    - **Section 1: Input Files & Folders:** Fill in the paths to your CSB data, tide zone file, and output directory.
    - **Section 2: Reference Bathymetry:**
        - Check "Use Automated BlueTopo Download" to have the script fetch reference data automatically.
        - OR, uncheck it and provide the path to your own local BAG or GeoTIFF file.
    - **Section 3: Processing & Export Options:**
        - **Insert into DuckDB:** (Recommended: ON) This is required for all post-processing and final gridding steps.
        - **Run Post-Processing:** (Recommended: ON) Enables the vertical offset calculation, uncertainty assignment, and outlier flagging. The "Final Gridding" option depends on this step.
        - **Export Individual Transit Files:** (Optional) If you want to inspect the results of the outlier flagging for every single transit, check this box. Be aware that this can create a very large number of files.
        - **Run Final Gridding & Export:** (Recommended: ON) Enables the final, primary output stage.
        - **Export Final Points GeoPackage:** (Recommended: ON) Creates the final `csb_final_points.gpkg` file.
        - **Optional Tessellation Shapefile:** Provide a shapefile here to enable tiled processing mode. If left blank, the script will produce a single set of final products.
        - **Grid Resolution (meters):** Set the cell size for the final output GeoTIFF.
        - **Organize GeoTIFFs ... and Create VRTs:** (Recommended: ON for tiled mode) Automatically sorts tiled outputs and builds virtual rasters.

4.  **Start Processing:** Click the "Start Processing" button. The application window will close automatically when the entire process is complete. Monitor the console for progress updates.

## Output Products

All outputs will be located in the specified **Output Directory**.

- `csb.duckdb`: The intermediate DuckDB database file containing all processed data.
- `master_offsets.csv`: A CSV file containing the calculated vertical offsets for each vessel.
- `/histograms/`: A folder containing PNG histograms showing the distribution of offsets for each vessel.
- `/transit_exports/`: (If enabled) A folder containing plots and GeoPackages for each individual vessel transit.
- `/final_products/`: This folder contains the primary deliverables.
    - `csb_final_points.gpkg`: (If enabled) A GeoPackage containing all processed points in EPSG:4326, with an `outlier` field indicating their status.
    - `csb_final_gridded.tif` or `<id>_gridded.tif`: The final gridded raster(s) in the appropriate UTM zone(s).
    - `/EPSG_.../`: (If VRT creation is enabled) Subfolders where the GeoTIFFs are organized by their coordinate system. Each folder will contain a `mosaic_EPSG_....vrt` file, which is the best file to load into GIS software.


***Contact Anthony Klemm anthony.r.klemm@noaa.gov for any questions regarding this tool*** 



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
