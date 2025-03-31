# IMPACTS_PHD README
README.md

Jon Seibert, Ph.D. Research
IMPACTS Assimilation Project Code Documentation
Last updated: 27 Mar 2025

This respository contains all scripts utilized in the creation of analysis plots and tables for my Ph.D. at PSU, entitled "Assimilating NASA IMPACTS Data for Feb 7 2020 Storm to Improve Winter Cyclone Prediction", submitted in Spring 2025. 

For details on each script, see internal code comments.
Below are descriptions of the source [code] for each figure in the dissertation, followed by descriptions of the data sources, and then of the primary mechanisms of the important Matlab scripts.

--------------------------------------------------------------------------------------------------------------------------------------------------------------
FIGURES

Fig. 1: 3 panels from WPC Sfc Analysis Archives, modified  
Fig. 2: 6 panels via ptp_analysis_ens_full.m (d01 500mb GPH, with combo_gph_full = true)  
Fig. 3: 1 panel from Tropical Tidbits.com, modified  
Fig. 4: 1 panel from NOAA National Gridded Snowfall Analysis, modified  
Fig. 5: 4 panels via plot_gis_2020_v3.m  
Fig. 6: 1 panel via p3_plotting.m  
Fig. 7: 1 panel via refl_struct_analysis_ex.m  
Fig. 8: Provided by Yunji Zhang, Ph.D.  
Fig. 9: Created in PowerPoint  
Fig. 10: 1 panel via ptp_analysis_ens_full.m, modified  
Fig. 11: 6 panels via ptp_analysis_ens_full.m  
Fig. 12: 3 panels via ptp_analysis_ens_full.m, modified  
Fig. 13: 15 panels via np_analysis.m  
Fig. 14: 2 panels via plot_error_tables.m  
Fig. 15: 2 panels via plot_error_tables.m, modified  
Fig. 16: 2 panels via plot_error_tables.m, modified  
Fig. 17: 1 panel via refl_struct_analysis_mb.m, modified  
Fig. 18: 3 panels via refl_struct_analysis_hgt.m  
Fig. 19: 1 panel via wrf_3d_ref_to_base.m  
Fig. 20: 2 panels via ptp_analysis_ens_full.m, modified  
Fig. 21: 6 panels via ptp_analysis_ens_full.m, modified  
Fig. 22: 6 panels via ptp_analysis_ens_full.m, modified  
Fig. 23: 6 panels via ptp_analysis_ens_full.m, modified  
Fig. 24: 4 panels via refl_struct_analysis_mb.m  
Fig. 25: 4 panels via refl_struct_analysis_mb.m  
Fig. 26: 1 panel from NOAA HYSPLIT Model Back Trajectory  
Fig. 27: 9 panels via ptp_analysis_ens_full.m, modified  
Fig. 28: 20 panels via ptp_analysis_ens_full.m, modified  
Fig. 29: 10 panels via ptp_analysis_ens_full.m, modified  
Fig. 30: 11 panels via ptp_analysis_ens_full.m, modified  
Fig. 31: 3 panels (b-d) via ptp_analysis_ens_full.m, 1 panel (a) made in Matlab command line from NOAA National Snowfall Gridded Analysis data  
Fig. 32: 4 panels via ptp_analysis_ens_full.m, modified  
Fig. 33: 1 panel (a) via ptp_analysis_ens_full.m, 1 panel (b) via snowband_tracking.m [Previous version of snowband_tracking_obs.m]  
Fig. 34: 5 panels via snowband_tracking_obs.m  
Fig. 35: 6 panels via snowband_tracking_obs.m  
Fig. 36: 6 panels via snowband_tracking_obs.m  
Fig. 37: 2 panels via snowband_tracking_obs.m  

Table 1: manually assembled from output tables of plot_error_tables.m  
Table 2: manually assembled from output tables of plot_error_tables.m  
Table 3: manually assembled from output tables of plot_error_tables.m  
Table 4: significance_test.m  
Table 5: snowband_dif.m  
Table 6: snowband_dif.m  
TAble 7: condense_snowband_tables.m  
Table 8: condense_snowband_tables.m  

-------------------------------------------------------------------------------------------------------------------------------------------------------------
DATA

Conventional Obs  
	Provided by Yunji Zhang via NCEP MADIS Site (https://madis.ncep.noaa.gov/)  

IMPACTS Obs  
	Aircraft obs provided as NetCDF and Excel files by Steven Guimond of NASA and Song Zhang  
	Soundings via NASA Earthdata (https://search.earthdata.nasa.gov/)  

WRF Outputs  
	PSU WRF-EnKF nested 9km-3km domain model run from 0 UTC Feb 7 - 0 UTC Feb 8, 40 members (Details in dissertation text)  
	NetCDF files in format "wrf_enkf_output_d01_002.nc" (analysis) and "wrfinput_d01_2020-02-07_14-00-00_002.nc" (background)  
	d01/d02 refers to outer domain (9km resolution) or inner (3km)  

NEXRAD Obs  
	IEM WSR-88D 0.5-deg radar reflectivity mosaics for 0Z Feb 7 to 0Z Feb 9 downloaded from IEM website.  
	Imported as images into matlab, converted from RGB to dBZ values.  
	(https://mesonet.agron.iastate.edu/docs/nexrad_mosaic/)  

Snowfall totals  
	NOAA National Gridded Snowfall analysis: provided figure 4, as well as data used to create figure 31  
	(https://www.nohrsc.noaa.gov/snowfall/)  

ERA-5 MSLP  
	(https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview)  

-------------------------------------------------------------------------------------------------------------------------------------------------------------

NASA IMPACTS Assimilation:  
	[ptp_analysis_ens_full.m]  
	NetCDF files for background and analysis downloaded from Roar Collab via Globus  
	Imported by variable from NetCDF files into Matlab  

1) Loop over available domains to find largest NEXRAD domain that is fully contained by the WRF domain, save indices and lat/lon values  
For each timestamp:  
	2) Trim NEXRAD data and lat/lon grid according to above  
	For each ensemble member:  
		3) Loop over all trimmed NEXRAD grid points:  
			- Find 4 nearest neighbors to target point in WRF grid (harder than it sounds, given irregular matching)  
			- Use bilinear interpolation to compute the value at the NEXRAD point  
			- Any points that could not find nearest neighbors are set to a filler value of -35 dBZ  
			- Replicate any NaNs present in the NEXRAD data  
		4) Apply radar mask to remove any data from the WRF fields that could not bee seen by the NEXRAD network  
		5) With newly regridded WRF data, compute direct difference and error statistics between WRF and NEXRAD reflectivity  
		6) Plot all raw inputs and differences  
     END  
	7) Output error table  
END  

-------------------------------------------------------------------------------------------------------------------------------------------------------------

To create artificial base reflectivity:  
[wrf_3d_ref_to_base.m]  
	
1) Import coordinates and elevations of all WSR-88D radar stations.  
2) Using an equation derived from the curvature of the Earth, compute the height of the beam for each grid point and each radar within its maximum range.  
For each timestamp:  
	For each ensemble member:  
		3) [Linear] Interpolate the simulated 3D reflectivity values to the radar beam heights.  
		4) Take the maximum value in each column from among those interpolated.  
		5) Plot  
	END  
	6) Save as .mat file  
END  

-------------------------------------------------------------------------------------------------------------------------------------------------------------

To create error tables and line plots:  
[plot_error_tables.m]  

For each experiment (AIR, CONV, etc):  
	For each variable (error value):  
		For each timestamp:  
			1) Read in error values from ptp_analysis output tables  
			2) Compute std of each variable from the ensemble for error bars  
			3) Store and plot each section of data  
		END  
	END  
END  
For each variable:  
	4) Plot stored error data across all times and experiments 
END 
5) Save large summary table of error values  

-------------------------------------------------------------------------------------------------------------------------------------------------------------

To create SLINK-clustered "snowband" objects and Demons-adjusted versions thereof:  
[snowband_tracking_obs.m]  

Must be run once separately for each timestep, due to limits on available computational power  

For each ensemble member:  
	1) Read in WRF output data and OBS data  
	2) Regrid OBS data to match WRF output (v.v. causes an overlarge matrix and memory issues)  
	3) Apply dBZ threshold to data  
	Twice:  
		4) On first iteration, choose original data. On second iteration, compute Demons image registration of WRF data towards OBS and use that.  
		5) Compute distance matrix between each above-threshold point in data  
		Until the distance between clusters exceeds the distance threshold:  
			6) Find point/cluster pairing with smallest remaining distance  
			7) Combine into 1 cluster, delete matrix entry of smaller one, and update distances  
		END  
		For each remaining cluster:  
			8) Check minimum size requirements  
			9) Check length:width ratio and size  
			10) Compute centroid, bounding box, and major/minor axes of cluster object  
		END  
		11) Save clustering data as .mat file  
		12) Plot points color coded by cluster with bounding boxes  
	END  
END  

-------------------------------------------------------------------------------------------------------------------------------------------------------------

To compute differences between "snowband" objects between two sets:  
[snowband_dif.m]  

For each timestamp:  
	For each ensemble member:  
		1) Load in SLINK clustering data from .mat files (WRF-D and OBS)  
		2) Remove clusters below minimum size req and outside region of interest  
		3) For each cluster in each dataset, find closest (Euclidean) cluster in opposing set  
		4) Use Demons displacement matrix to match WRF-D clusters to closest (Euclidean) WRF clusters  
		5) Compute distance between pre-Demons WRF clusters to newly matched OBS, weight according to cluster size  
		6) Save as table and boxplot  
	END  
END  
7) Make one giant boxplot of all distance error values  

-------------------------------------------------------------------------------------------------------------------------------------------------------------
