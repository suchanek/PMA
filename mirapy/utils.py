"""mirapy module for astronomical data analysis, Eric G. Suchanek, Ph.D., MIRA"""
# @Author: Eric G. Suchanek, Ph.D.
# @License: GNU
# @Version: __version__
# Microsoft Azure PAT: u3orvr6bwdm7fpmulnhw7bvv63c3njc2s2wqt3etnidb5joxdsga
#
# get rid of the pylint warnings about invalid cases for variable names, line length
# pylint: disable=C0103
# pylint: disable=c0303
# opylint: disable=C0111

import os
from os.path import expanduser
from time import time
import datetime
import socket
import numpy as np

from astropy.coordinates import spherical_to_cartesian, SkyCoord
import astropy.units as u

from astroquery.simbad import Simbad

# make a dict of system names and system information. Easily extended with more hosts.

host_systems = {"Atlas": "PC, 3.0 GHz i5, 16 GB DDR4 2166 MHz, Ubuntu 18.04",
                "KFlux.Home": "iMac, 3.3 GHz i5, 32 GB DDR3 1867 MHz, OSX 10.15",
                "parallels": "PC, 3.0 GHz i5', 16 GB DDR3 1867 MHz, Windows 10, Parallels",
                "TBD": "Generic PC",
                "newton": "Macbook Pro, 2.5 i7, 16GB 1066 MHz DDR3, OSX 10.15",
                "spider": "Dell laptop, 2.5 GHz, i5, 8 GB, Ubuntu 18.04",
                "MacbookBootcamp": "Macbook Pro, 2.5GHz i7, 16GB 1066 MHz DDR3, Windows 10",
                "TBD2": "iMac, 3.3 GHz i5, 32 GB DDR3 1867 MHz, Ubuntu 18.04"
               }


# returns the extents of each input vector shunk by sh (decimal percent)
def myshrink(X, Y, Z, sh):
    """shrinks the input vector extents by the input sh factor

    Parameters:
    ----------
    X, Y, Z : arrays to be scaled
    sh      : scaling percentage

    Returns:
    ----------
    xmin, xmax, ymin, ymax, zmin, zmax : shrunken values
    """
    xmin = np.min(X)
    xmax = np.max(X)
    ymin = np.min(Y)
    ymax = np.max(Y)
    zmin = np.min(Z)
    zmax = np.max(Z)
    zmin = zmin + sh * zmin
    zmax = zmax - sh * zmax
    xmin = xmin + sh * xmin
    xmax = xmax - sh * xmax
    ymin = ymin + sh * ymax
    ymax = ymax - sh * ymax
    return xmin, xmax, ymin, ymax, zmin, zmax

def get_home_dir():
    """Return the directory for the current user. Works correctly for Windows, Linux, OSX
    Parameters:
    ----------
    null

    Returns:
    ----------
    home directory
    """
    return expanduser("~")

def ensure_dir_path(file_path):
    """check to see if the file_path exists. Create the path if needed
    Parameters:
    ----------
    file_path : full path

    Returns:
    ----------
    null, creates path if needed
    """
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)

# return a pretty string with formatted time and machine architecture

def pprint_delta(secs):
    """Print a pretty string with formatted time and machine architecture

    Parameters:
    ----------
    elapsed : elapsed time in seconds

    Returns:
    -------
    null
    
    Output:
    --------
    pprinted string
    """
    host = str(socket.gethostname())
    
    specs = host_systems.get(host)
    if specs:
        output2 = "Specs: " + host_systems.get(host)
    else:
        output2 = "Undefined system!"
    
    output = "Elapsed: " + str(datetime.timedelta(seconds=secs)) + " on: "
    print(output)
    print(output2)

# return a pretty string with formatted time and machine architecture
def pprint_elapsed(start, verbose=False):
    """Print a pretty string with formatted time and machine architecture

    Parameters:
    ----------
    start: starting time, eg start = time.Time()

    Returns:
    ----------
    null

    Output:
    ----------
    pprinted string
    """
    host = ""
    host = str(socket.gethostname())
    if verbose:
        print('Host is: ' + host)
    
    specs = host_systems.get(host)
    if specs:
        output2 = " a " + host_systems.get(host)
    else:
        output2 = "an undefined system!"
    
    end = time()
    output = "Elapsed:" + str(datetime.timedelta(seconds=end-start)) + \
    " on " + host + output2
    
    print(output)
 
def make_output_filename(project_dir, file_prefix, ra, dec, ra_width, dec_width,
                         sig=0.0, etype='none'):
    """construct a fully-qualified file name representing a Gaia PM & velocity file

    Parameters:
    ----------

    project_dir: (str) path to the directory where the project data are stored

    file_prefix: (str) the prefix for the data filenames

    ra: ra coordinates in string form for the query center

    dec: dec coordinates in string form the for the query center
    
    ra_width: width (deg) for ra extent
    
    dec_width: width (deg) for dec extent

    sig: significance value, if > 0 then the filename has the significance figure added

    etype: signficance type, either 'para' for parallax, 'both' for pm and parallax,
            or 'none' for none, (default)

    Returns:
    -------
    Full path for the file. (str)
    """
    home_dir = get_home_dir()
    project_name = home_dir + project_dir + file_prefix

    if sig > 0.0:
        # format the filenames using the format: projectname_ra1_dec1_rawidthxdecwidth_deg.csv
        # raw output filename for the GAIA query
        # output filename for the data file with appended 3d velocity data and PM angle
        vel_filename = str.format("{}{}_{}_{}x{}_{}sig_{}_vel.csv", project_name, ra, dec,
                                  round(ra_width, 2), round(dec_width, 2), sig, etype)
    else:
        # output filename for the data file with appended 3d velocity data and PM angle
        vel_filename = str.format("{}{}_{}_{}x{}_all_vel.csv", project_name, ra, dec,
                                  round(ra_width, 2), round(dec_width, 2))
    return vel_filename

def normalized(a, axis=-1, order=2):
    """normalize a ndim vector
    Parameters:
    ----------
    a : vector to be normalized

    axis: x

    order: normalization order, 2 is euclidean

    Returns:
    -------
    normalized vector
    """
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2 == 0] = 1
    return a / np.expand_dims(l2, axis)

# take a df of stars, extract ra, dec, distance, convert and return x, y, z
def extract_cartesian_coords(stars):
    """Return cartesian coordinates for a DataFrame containing ra, dec, and distance.
    
    Parameters:
    ------
        stars: DataFrame 
            DataFrame containing ra, dec and distance values
        
    Returns:
    -------
        x: array
            The first cartesian coordinate
            
        y: array
            The second cartesian coordinate
            
        z: array
            The third cartesian coordinate
    """
    
    ra = np.asarray(stars['ra'])
    dec = np.asarray(stars['dec'])
    # print(ra[:10])
    # NB: the conversion function requires all units to be the same, hence
    # the use of parallax rather than 1/parallax here.
    dist = np.asarray(stars['parallax'])
    x, y, z = spherical_to_cartesian(dist, dec, ra)
    return x, y, z
   
# end of file

def query_Simbad_obj_SC(obj_name):
    """return a SkyCoord object representing the center point of the input obj_name

    Parameters:
    ----------
    obj_name - object name for the query, e.g. 'NGC5000'

    Returns
    ----------
    SkyCoord object initialized with the object coordinates
    """
    result_table = Simbad.query_object(obj_name)
    df = result_table.to_pandas()
    df['RA'] = df['RA'].astype('str')
    df['DEC'] = df['DEC'].astype('str')
    df['RA'] = df['RA'].str.replace(' ', ':')
    df['DEC'] = df['DEC'].str.replace(' ', ':')
    ra = df['RA'].iloc[0]
    dec = df['DEC'].iloc[0]
    sc = SkyCoord(ra, dec, unit=(u.hourangle, u.degree))
    return sc

def query_Simbad_obj_pos(obj_name):
    """return ra, dec as formatted strings for an input object

    Parameters:
    ----------
    obj_name - object name for the query, e.g. 'NGC5000'

    Returns
    ----------
    ra - string representing ra of object center
    dec - string representing dec of object center
    """
    result_table = Simbad.query_object(obj_name)
    # convert table object to pandas DataFrame
    df = result_table.to_pandas()
    # the coordinates are dtype(object), convert to string
    df['RA'] = df['RA'].astype('str')
    df['DEC'] = df['DEC'].astype('str')
    # replace the spaces with : to make a valid positional string
    df['RA'] = df['RA'].str.replace(' ', ':')
    df['DEC'] = df['DEC'].str.replace(' ', ':')
    return df['RA'].iloc[0], df['DEC'].iloc[0]

#end of file
