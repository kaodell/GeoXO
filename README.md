# README
This repository contains scripts and python3 environments used to produce the analysis and figures in “Public health benefits from improved identification of severe air pollution events with geostationary satellite data” by O’Dell et al. soon to be submitted to PNAS. If you have any questions about these codes, want to use them, or find any errors, please contact Kate O’Dell at the contact email provided with this github page. 
## Instructions
To complete the analysis and create figures in the paper, the python scripts in this folder are run in the following order. Note that the local and remote codes use different python environments. Both environments are provided in this folder and are entitled py_env_local.yml and py_env_remote.yml, for the local and remote environments, respectively. 

### Step 1) HIA
#### Prep_pop_data.py
This script loads the NASA SECAD GPWv4.11 and regrids the population density using a nearest neighbors interpolation. Run locally.

Inputs:
- 0.01 x 0.01 degree grid used for lat/lon for satellite-based PM2.5 datasets created by Hai Zhang
- Population density from NASA SEDAC GPW v4.11, available: https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11, accessed 27 Oct 2021
- 2018 TIGER/Line national shapefile, available: https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html, accessed 27 Oct 2021
- 2020 TIGER/Line states shapefile, avilable: https://www.census.gov/geographies/mapping-files/2020/geo/tiger-line-file.html, accessed 27 Oct 2021

Outputs:
- netCDF file containing regridded population, population density, grid area, and a US mask and california mask for the grid

#### Grid_baseline_mort.py
This script loads county level baseline rates from the CDC wonder database and grids them to the 0.01 degree grid. Run remote.

Inputs:
- County-level and state-level baseline mortality rates for 2015-2019 from the CDC WONDER database, https://wonder.cdc.gov/ucd-icd10.html, accessed 20 Sept 2022
- Grid for the satellite-based PM2.5 datasets
- 2019 TIGER/Line shapefiles for counties, available: https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html, accessed 14 Sept 2022

Outputs:
- Gridded baseline mortality rates at the county level

#### calc_annual_alert_HIA.py (checked, running now)
This script calculates the annual averages, observation completeness, and alert days for the three satellite based PM2.5 datasets. It also runs a health impact assessment for the ABI-daytime and ABI-1pm datasets. The script needs to be run separately for each dataset and is run on a remote, high memory server.

Inputs:
- Gridded baseline mortality rates output from grid_baseline_mort.py
- Gridded population output from prep_pop_data.py
- Gridded surface pm2.5 estimates based on the three satellite-based AOD datasets from Zhang et al., 2022 https://doi.org/10.1175/WAF-D-22-0114.1 regridded to a 0.01 x 0.01 degree grid by co-author Zigang Wei.

Outputs:
- netCDF containing health impact assessments, counts of days with observations, alert days, area covered each day, and population-weighted PM

#### Calc_national_ervs.py
This script calculates a HIA of asthma emergency room visits (ERVs) for daily reductions in pm exposure for the ABI-Daytime and ABI-1pm datasets for the year 2020. Run remotely.

Inputs:
- Gridded population output from prep_pop_data.py
- Gridded surface pm2.5 estimates based on the two ABI-based AOD datasets from Zhang et al., 2022 https://doi.org/10.1175/WAF-D-22-0114.1 regridded to a 0.01 x 0.01 degree grid by co-author Zigang Wei.

Outputs:
- netCDF file containing annual total pm-attributable asthma ERVs for each dataset baseline case and alert-behavior modification exposure reduction case 

### Step 2) Process AirNow Data
#### Airnow_2020_dailyareameanmax.py (written by co-author Barron H. Henderson w/small modifications by Katelyn O'Dell) - not included here, access by contacting manuscript authors.
This script calculates the daily mean and max pm2.5 for each reporting area using the airnow monitor data. The script was written by Barron Henderson with some minor modifications from myself (Kate) and is run locally. The edits I made to this script also draw from earlier code written by Barron, AirNowEquivalent.pdf

Inputs:
- TIGER/Line ZCTA shapefile accessed from IPUMS NHGIS for the 2010 census ZCTAs, access date 10 May 2023
- State shapefile from https://www.census.gov/geographies/mapping-files/2018/geo/carto-boundary-file.html, accessed 28 Sept 2022

Outputs:
- .dat file containing the mean and max pm2.5 by reporting area for each day in 2020
- Csv files for each day in 2020 with the airnow-reporting area pm2.5 concentrations assigned by ZCTAs. 

#### count_AQS_coverage.py
This script assigns the individual airnow monitors to populations based on census tracts and counts number of observations, alerts, and alerted populations from monitors if only assigned to local census tracts. Run locally. Takes about 1.5 hours to run.

Inputs:
- TIGER/Line census tract shapefile and population estimates accessed from IPUMS NHGIS for the 2020 census, standardized to 2010 geometries, accessed 10 May 2023

Outputs:
- Csv of count of observations, alert days, and population by census tract

#### Count_airnowalert_coverage.py
This script uses the assignment of monitors to reporting areas to count observations, alerts, and alerted populations for the airnow dataset by ZCTA. Run locally.

Inputs:
- TIGER/Line ZCTAs shapefile and population estimates accessed from IPUMS NHGIS for the 2020 census, standardized to 2010 geometries, accessed 10 May 2023
- ZCTAs and county geoid crosswalk file from US Census Bureau, available at: https://www.census.gov/geographies/reference-files/time-series/geo/relationship-files.2010.html#zcta accessed 9 May 2023
- Csv files output from Airnow_2020_dailyareameanmax.py
- Csv of state fips codes and US postal codes

Outputs:
- Csv with count of observations, alert days, and population by ZCTA as defined for the 2010 census

### Step 3) Make figures
#### Plot_national_results.py (checked - final adjustment needed for final annual HIA data)
This script takes inputs of counts, alerts, populations, and health impact assessments for the different pm2.5 datasets. conducts statistical analysis, and creates all figures in the paper. Run locally. Needs to be run separately for sensitivity analysis to the percent PM reduction due to behavior modification on alert days.

Inputs:
- Abi, viirs, and abi1pm HIA files output from calc_annual_alert_HIA.py
- HIA ERV totals output from calc_national_ervs.py
- Regridded population output from prep_pop_data.py
- Grided baseline mortality rates output from grid_baseline_mort.py
- County-level baseline mortality rates for 2015-2019 from the CDC WONDER database (for checking gridded results)
- Csv of monitor results assigned to census tracts output from count_AQS_coverage.py
- Csv of monitor results assigned to reporting areas output from count_airnowalert_coverage.py
- TIGER/Line ZCTAs shapefile accessed from IPUMS NHGIS for the 2020 census, standardized to 2010 geometries, accessed 10 May 2023
- 2019 TIGER/Line shapefiles for counties, available: https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html, accessed 14 Sept 2022
- State shapefile from https://www.census.gov/geographies/mapping-files/2018/geo/carto-boundary-file.html, accessed 28 Sept 2022

Outputs:
- All figures included in the manuscript

