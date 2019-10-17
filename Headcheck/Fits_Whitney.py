#!/usr/bin/env python
# coding: utf-8

# # Head Check Analysis
# *Eric G. Suchanek, Ph.D. v.1, 3/9/19*
# 

# In[103]:


#
# Analyze fit files from the Axiom at Chews Ridge as part of Whitney's 'headcheck' procedure.
# This code reads .fits frames from a specified path,extracts the obs-date, exposure time, 
# airmass, time start, object, RA, DEC.
# The program computes the time mid, time end, HA, Airmass, 
# in UTC and LST by using the coordinates of Chews Ridge, This is used to compute the object's
# HA, altitude and airmass. These results are displayed in tabular form. The raw image and a stretched
# version are shown side by side. The stretched version is normalized over the image intensities
# code automatically adjusts the timezone offset by checking to see if DST is in effect.
# 
# -egs- 3/9/19 

from glob import glob
import os

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Latitude, Longitude
from astropy.io import fits
from astropy.time import Time
import time
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize,PercentileInterval)
import datetime

#
# DST calculation for a given date
# expects date to be .unix format
def is_dst(date):
    return bool(time.localtime(date).tm_isdst)

# OOS coordinates taken from Google Earth
sitelong = '-121:34:1'
sitelat = '36:18:20'
height = 1516*u.m
site_location = EarthLocation(sitelong, sitelat, height=height)

# Base constant UTC offset for US/Pacific. This can change to -7 hour when DST is in effect.
# the local var utcoffset reflects the ACTUAL utcoffset for the given date, taking DST into account
UTC_offset = -8*u.hour

# grab only the light frames. this command produces a list of files matching
# the string L*.fits
path = './PHOT030727/C/'
files = glob(path + 'L*.fits')

start_time = time.time()
count = 0

# Loop over all files
for fits_image_file in files:
    data, hdr = fits.getdata(fits_image_file, header=True, ext=0)
    # collect parameters from the hdr
    date_obs = hdr['date-obs']
    exptime = hdr['exptime']
    obj_name = hdr['object']
    exposure_time = exptime*u.second
    airmass = hdr['airmass']
    obj_ra = hdr['ra']
    obj_dec = hdr['dec']
    obj_detector = hdr['detector']
    obj_filter = hdr['filter']
    obj_darktime = hdr['darktime']
    
    # convert the strings for ra and dec into actual angle objects for later use
    ra = Longitude(obj_ra, unit=u.hourangle)
    dec = Latitude(obj_dec, unit=u.deg)

    # Do the time calculations - based on camera time from the hdr
  
    UTC_time_observation = Time(date_obs,format='fits',location=site_location)
    UTC_mid_exposure_time = UTC_time_observation + exposure_time / 2
    UTC_end_exposure_time = UTC_time_observation + exposure_time
    # check if dst is in effect
    if is_dst(UTC_time_observation.unix):
        utcoffset = UTC_offset + 1*u.hour
        dst_string = 'PDT'
    else:
        dst_string = 'PST'
            
    local_time_observation = Time(UTC_time_observation + utcoffset, location=site_location)
    local_time_observation_mid = local_time_observation + exposure_time / 2
    local_time_observation_end = local_time_observation + exposure_time
    
    sidereal_time_type = 'apparent' # or 'mean'
    sidereal_time_local_start = local_time_observation.sidereal_time(sidereal_time_type) 
    sidereal_time_local_mid = local_time_observation_mid.sidereal_time(sidereal_time_type) 
    sidereal_time_local_end = local_time_observation_end.sidereal_time(sidereal_time_type) 

    sidereal_time_UTC_start = UTC_time_observation.sidereal_time(sidereal_time_type) 
    sidereal_time_UTC_mid = UTC_mid_exposure_time.sidereal_time(sidereal_time_type) 
    sidereal_time_UTC_end = UTC_end_exposure_time.sidereal_time(sidereal_time_type) 

    # camera's coordinates and time as input:  UTC
    object_Sky = SkyCoord(ra,dec)
    object_altaz_camera_start = object_Sky.transform_to(AltAz(obstime=UTC_time_observation,
                                                              location=site_location))
    computed_airmass_start = object_altaz_camera_start.secz
    
    object_altaz_camera_mid = object_Sky.transform_to(AltAz(obstime=(UTC_mid_exposure_time),
                                                            location=site_location))
    computed_airmass_mid = object_altaz_camera_mid.secz
    hourangle_mid = sidereal_time_local_mid - object_observed_camera_mid.ra
    
    print('File:      ', fits_image_file)
    print('Object:    ', obj_name)
    print('Exposure:  ', exposure_time)
    print("Local Time: {} {}".format(local_time_observation.iso, dst_string))
    
    norm = ImageNormalize(data, interval=PercentileInterval(.1))
    fig = plt.figure(figsize=(8,3),dpi=120)
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.imshow(data, origin='lower')
    ax2 = fig.add_subplot(1,2,2)
    im2 = ax2.imshow(data, origin='lower',norm=norm)
    fig.colorbar(im2)
    plt.show()
    
    print('RA:                        ', ra)
    print('Dec:                       ', dec)
    print('HA mid:                    ', hourangle_mid)
    print("Altitude (mid exp):         {0.alt:.6}".format(object_altaz_camera_mid))
    print('Exposure:                  ', exposure_time)
    print('Airmass (hdr):             ', airmass)
    print('Airmass mid,(computed):     {0:.4}'.format(computed_airmass_mid))
    print('JD UTC start exp, (hdr):   ', UTC_time_observation.jd)
    print('JD UTC mid exp:            ', UTC_mid_exposure_time.jd)
    print('JD UTC end exp:            ', UTC_end_exposure_time.jd)
    print('UTC start exp:             ', UTC_time_observation.iso)
    print('UTC mid exp:               ', UTC_mid_exposure_time.iso)
    print('UTC end exp:               ', UTC_end_exposure_time.iso)
    print('Local Time start exp:      ', local_time_observation.iso)
    print('Local Time mid exp:        ', local_time_observation_mid.iso)
    print('Local Time end exp:        ', local_time_observation_end.iso)
    print('LST start:                 ', sidereal_time_local_start)
    print('LST mid:                   ', sidereal_time_local_mid)
    print('LST end:                   ', sidereal_time_local_end)
    print('GST start:                 ', sidereal_time_UTC_start)
    print('GST mid:                   ', sidereal_time_UTC_mid)
    print('GST end:                   ', sidereal_time_UTC_end)
    print('---------------------------------------------------')
    count += 1
    # fits.writeto(fits_image_file, data, hdr, overwrite=True)
end_time = time.time()
elapsed_time = (end_time - start_time)*u.s

print('Processed', count, 'files in:', elapsed_time)


# In[ ]:




