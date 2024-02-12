# atm_bathymetry_toolkit
This repository contains MATLAB® functions for deriving supraglacial lake bathymetry from ATM laser altimetry data products

# Airborne Topographic Mapper (ATM) Bathymetry Toolkit

The Airborne Topographic Mapper (ATM) was a scanning lidar developed and used by NASA for observing the Earth’s topography for several scientific applications, foremost of which was the measurement of changing Arctic and Antarctic ice sheets, glaciers and sea ice. ATM measured topography to an accuracy of better than 5 centimeters by incorporating measurements from GPS (global positioning system) receivers and inertial navigation system (INS) attitude sensors.

This collection of MATLAB® functions is for working with ATM laser altimetry data products in HDF5 waveform format as well as ATM’s Continuous Airborne Mapping by Optical Translator (CAMBOT) three-channel, natural color, red, green, and blue (RGB) digital images. The functions provide tools to geolocate laser footprints, identify supraglacial hydrological features such as lakes and channels and estimate their water depth from lidar data.

The data products are freely available at the National Snow and Ice Data Center (NSIDC) at https://nsidc.org/data/icebridge and can also be downloaded from the NASA Earthdata portal at https://earthdata.nasa.gov/

The ILATMW1B waveform data is available at NSIDC:  
https://nsidc.org/data/ILNSAW1B/versions/1 (narrow swath)  
https://nsidc.org/data/ILATMW1B/versions/1 (wide swath)

The CAMBOT Level 0 raw data and the much larger Level 1B geolocated and orthorectified CAMBOT data product are available here:  
https://nsidc.org/data/IOCAM0/  
https://nsidc.org/data/iocam1b  

The methods for estimating water depths of supraglacial lakes on the Greenland Ice Sheet are described in a publication by Studinger et al., 2022 (https://doi.org/10.5194/tc-16-3649-2022). The code release for this publication is archived on Zenodo at: https://zenodo.org/records/6341230

__See also__: 

* User guide for NASA's Airborne Topographic Mapper HDF5 Waveform Data: Products: https://doi.org/10.5281/zenodo.7246097
* Collection of MATLAB® functions for working with ATM (Airborne Topographic Mapper, laser altimetry data products in HDF5 waveform format: https://github.com/mstudinger/ATM-waveform-tools
* NASA's Airborne Topographic Mapper (ATM) ground calibration data for waveform data products: https://doi.org/10.5281/zenodo.7225936
* An __archived version of the code__ that was used in the publication by Studinger et al., 2022 (https://doi.org/10.5194/tc-16-3649-2022) is available on Zenodo: https://zenodo.org/records/6341230
