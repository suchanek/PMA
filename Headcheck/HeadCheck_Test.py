#!/usr/bin/env python
# coding: utf-8

# In[1]:


# 
# Create a program to read input from the user per the Headcheck manual, load the .fits file and then calculate the 
# various times and positions needed. Interactively take the data and append into lists for bulk processing
#

from glob import glob
import os
import time
import datetime
import sys

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Latitude, Longitude
from astropy.io import fits
from astropy.time import Time

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from astropy.visualization import (SqrtStretch,ImageNormalize,PercentileInterval)

#
# DST calculation for a given date
# expects date to be .unix format

def is_dst(date):
    return bool(time.localtime(date).tm_isdst)

#
# Main code to input and analyze .fits files per the Headcheck protocol
# Set up our Observatory location.
# OOS coordinates taken from Google Earth
# convert to geodetic coords via the EarthLocation class. We need this for sidereal calculations

sitelong = '-121:34:00'
sitelat = '36:18:20'
height = 1525*u.m
observatory_location = EarthLocation(sitelong, sitelat, height=height)

# Base constant UTC offset for US/Pacific. This can change to -7 hour when DST is in effect.
# the local var utcoffset reflects the ACTUAL utcoffset for the given date, taking DST into account
UTC_offset = -8*u.hour

# default path is local dir

# global vars
_path = './'
_ra_str = ''
_dec_str = ''
_airmass = 0.0
_time_obs = ''
_file = ''
_exposure = 0.0
_utc_time_obs = ''
_file_list = [] # hold all files entered for processing
_file_count = 0

# toggle image plotting
Plot_images = 'False'

getting_input = 'True'

path_prompt = 'Enter a path for the data files relative to the path of this script: '
fname_prompt = 'Enter a file name (no path): '
ra_prompt = 'Enter RA (hh:mm:SS): '
dec_prompt = 'Enter Dec (deg:mm:SS): '
airmass_prompt = 'Enter Airmass: '
exp_prompt = 'Enter Exposure (sec): '

cmd_prompt = 'Enter a command, q to exit and process: '
cmd_list = "\n******************* HeadCheck Input *******************\np --- Enter path\nf --- Enter a filename to be processed \nl --- List files\no --- Enter Observation Date and time \nr --- Enter RA \nd --- Enter Dec \nm --- Enter Airmass \nv --- Review Input \nx --- Enter Exposure \nz --- Reset all \nq --- Quit and Process"

path = '.'
from pathlib import Path

def exit():
    getting_input = False
    return

def get_path():
    path = input(path_prompt)
    my_path = Path(path)
    if my_path.exists():
        ("")
    else:
        print('Path: ', path, 'cannot be found!')
        input('Press Return to continue')
        return ''
    return path

# list only .fits files
def list_files(path):
    i = 0
    files = glob(path + '*.fits')
    if len(files) == 0:
        print("No files found.")
    else:
        print('\nFiles in path:', path)
        for file in files:
            print(file)
            i+=1
    return(i)

def enter_ra():
    _r = input(ra_prompt)
    return _r

def enter_dec():
    _d = input(dec_prompt)
    return _d

def enter_airmass():
    _am = input(airmass_prompt)
    return _am

def enter_exposure():
    _exp = input(exp_prompt)
    return _exp

def add_file():
    _fname = input(fname_prompt)
    _file = _path + _fname      
    my_file = Path(_file)
    if my_file.is_file():
        ("")
    else:
        print('File', _file, 'cannot be opened. Check path and filename!')
        input('Press Return to continue')
        return ''
    return _file

def reset_all():
    _airmass = 0.0
    _ra_str = ''
    _dec_str = ''
    _time_obs = ''
    _file_list = []
    _file = ''
    _path = './'
    _exposure = 0.0
    _utc_time_obs = ''
    
    return

def enter_time_obs():
    _timeobs = input("Enter date/time of the observation: eg 2019-3-09T11:34:22")
    try:
        _tst = Time(_timeobs,format='fits',location=observatory_location)
    except:
        print('Invalid time format. Try again')
        return('')
    return _timeobs

def review_input():
    print('*********** Current Perameters *************')
    print(' Path: ', _path)
    print(' File: ', _file)
    print(' Ra: ', _ra_str)
    print(' Dec: ', _dec_str)
    print(' Observation date/time: ', _utc_time_obs)
    print(' Exposure: ', _exposure)
    print(' Airmass: ', _airmass)
    print('************************')
    return 


def str_to_functions_to_strings(argument):
    switcher = {
        'f': add_file(_path),
        'l': list_files(_path),
        'o': enter_time_obs,
        'r': enter_ra,
        'd': enter_dec,
        'p': get_path,
        'm': enter_airmass,
        'v': review_input,
        'x': enter_exposure,
        'v': review_input,
        'z': reset_all,
        'q': exit
    }
    # Get the function from switcher dictionary
    func = switcher.get(argument, lambda: "nothing")
    # Execute the function
    return func()

getting_input = True

while getting_input is True:
    # main input loop
    os.system('clear')  # on linux / os x
    print(cmd_list)
    cmd = input(cmd_prompt)
    
    if cmd is 'q':
        getting_input = False
    elif cmd is 'f':
        _file = add_file()
        if _file is not '':
            _file_list.append(_file)
    elif cmd is 'l':
        file_count = list_files(_path)
        input('Press return to continue')
    elif cmd is 'r':
        _ra_str = enter_ra()
    elif cmd is 'd':
        _dec_str = enter_dec()
    elif cmd is 'p':
        _path = get_path()
    elif cmd is 'o':
        _utc_time_obs = enter_time_obs()
    elif cmd is 'm':
        _airmass = enter_airmass()
    elif cmd is 'x':
        _exposure = enter_exposure()
    elif cmd is 'v':
        review_input()
        input('Press return to continue')
    elif cmd is 'z':
        reset_all()
    else:
        print('Invalid Command!\n')
            
    #print(getting_input)

print('Done')


# In[ ]:




