Downscaling code for GCM, developed by Mamad Tamamadin, Department of Meteorology, Institut Teknologi Bandung. The method is adopted from Lima et al (2021) with paper title "A Bayesian Kriging model applied for spatial downscaling of daily rainfall from GCMs".  

Guide for executing this python code:
1. Install python 
2. Install the libraries: netCDF4, scikit-gstat, pykrige, matplotlib, numpy, pandas, scipy, pyshp, pyproj, pathlib, sridentify
3. Create the home directory, for example: "climate_change", and create some directories under this directory, as follows:
4. Create the directory "script" and store the this code file in this directory. 
5. Create the directory "csv" and store the observed rainfall data in each csv files in this directory.  
6. Create the directory "subdas" and store the shapefile of subdas in this directory.
7. Create the directory "output". You can see the output result in this directory.
8. Create the directory "original". You can see the original csv file of GCM after you finished the downscaling.
9. Create the directory "GCM" and store the GCM files in this directory.
10. Command for executing this code: "python downscaling_GCM_bayesian_kriging.py subbasin_name gcm_file start_year end_year resolution_in_km"
11. For example: "python downscaling_GCM_bayesian_kriging.py indonesia 01b-pr_day_CMCC-CM2-SR5_ssp585_r1i1p1f1_gn_20400101-20641231 2020 2030 20"
