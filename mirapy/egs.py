#"""Subroutines for astronomical data analysis, authored by Eric G. Suchanek, Ph.D."""
# Author: Eric G. Suchanek, Ph.D.
# License: GNU
# These routines utilize a number of astronomical packages to perform data analysis
# on sets of stars. This particular set is tuned towards data extracted from the Gaia database.
# The nature of the datastructures used is important. There are currently two.
# 1) pandas dataframes - these are used as the primary way of storing the stellar data imported
# from Gaia, since the .csv import naturally leads to one. It is easy to read and manipulate
# the individual rows and columns within the dataframe (just like in a spreadsheet).
#
# Microsoft Azure PAT: u3orvr6bwdm7fpmulnhw7bvv63c3njc2s2wqt3etnidb5joxdsga
#
# get rid of the pylint warning about invalid cases for variable names and a few others
# pylint: disable=C0103
# pylint: disable=c0303
# pylint: disable=C0111

# pylint: disable=E1101
# pylint: disable=C0302

from time import time
import os
import numpy as np
import pandas as pd
from pandas import DataFrame
from pandas import read_csv
import progressbar

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import Distance
from astropy.coordinates import Latitude, Longitude, SkyCoord, ICRS

from astroquery.gaia import Gaia

import mirapy.utils as ut

# SQL query strings for various GAIA astronomical searches
#
# The significance calculations are taken from:
# https://www.gaia.ac.uk/data/gaia-data-release-1/adql-cookbook
#
# extract stars whose pmra and pmdec significance is X sigma, non-null parallax
# box query that only pulls non-null PMRA, PMDEC and PARALLAX
sel_str_box_pm_sig = "SELECT DISTANCE(POINT('ICRS',ra,dec), POINT('ICRS',{},{})) AS dist, * \
FROM gaiadr2.gaia_source  \
WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS',{},{},{},{}))=1 \
AND pmra*pmra + pmdec*pmdec > {}*SQRT(pmra*pmra*pmra_error*pmra_error + pmdec*pmdec*pmdec_error*pmdec_error) \
AND parallax > 0 \
AND pmra IS NOT NULL \
AND pmdec IS NOT NULL; \
order by dist asc"

# pull stars whose parallax and pmra, pmdec significance is X sigma
sel_str_box_both_sig = "SELECT DISTANCE(POINT('ICRS',ra,dec), POINT('ICRS',{},{})) AS dist, * \
FROM gaiadr2.gaia_source  \
WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS',{},{},{},{}))=1 \
AND pmra*pmra + pmdec*pmdec > {}*SQRT(pmra*pmra*pmra_error*pmra_error + pmdec*pmdec*pmdec_error*pmdec_error) \
AND parallax*parallax >{}*SQRT(parallax*parallax*parallax_error*parallax_error) \
AND parallax > 0 \
AND pmra IS NOT NULL \
AND pmdec IS NOT NULL; \
order by dist asc"

# pull stars whose parallax significance is X sigma, non-null pmra, pmdec
sel_str_box_para_sig = "SELECT DISTANCE(POINT('ICRS',ra,dec), POINT('ICRS',{},{})) AS dist, * \
FROM gaiadr2.gaia_source  \
WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),BOX('ICRS',{},{},{},{}))=1 \
AND parallax > 0 \
AND parallax*parallax >{}*SQRT(parallax*parallax*parallax_error*parallax_error) \
AND pmra IS NOT NULL \
AND pmdec IS NOT NULL; \
order by dist asc"

#
# https://www.gaia.ac.uk/data/gaia-data-release-1/adql-cookbook
# box query that only pulls non-null PMRA, PMDEC and PARALLAX

# pull non-null pm data. Can easily add additional columns
# calculate and add a distance column
sel_str_box_filt = "SELECT DISTANCE(POINT('ICRS',ra,dec), POINT('ICRS',{},{})) AS dist, \
source_id, ra, dec, ra_error, dec_error, parallax, parallax_over_error, pmra,\
pmra_error, pmdec, pmdec_error, phot_g_mean_mag FROM gaiadr2.gaia_source \
WHERE CONTAINS(POINT('ICRS', ra, dec), BOX('ICRS',{},{},{},{}))=1 \
AND parallax > 0 \
AND pmra IS NOT NULL AND pmdec IS NOT NULL; \
order by dist asc"

sel_str_box_all = "SELECT DISTANCE(POINT('ICRS',ra,dec), POINT('ICRS',{},{})) AS dist, * \
FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS', ra, dec), BOX('ICRS',{},{},{},{}))=1 \
AND parallax IS NOT NULL \
AND pmra IS NOT NULL AND pmdec IS NOT NULL; \
order by dist asc"

# astronomy calculations begin
# compute proper motion angle, return in degrees, correct for quadrant
# bruce's formula forced into 0-360


def bruce(x, y):
    """compute the PM angle between created by the x, y, cartesian coordinates

    Parameters
    ----------
    x, y position

    Returns
    ----------
    PM angle between 0-360 degrees
    """
    # rho = np.sqrt(x**2 + y**2)
    ra = x * np.pi / 180.0
    dec = y * np.pi / 180.0

    phi = np.arctan(ra/dec)
    phi *= 180.0 / np.pi  # to degrees
    if dec < 0.0:
        phi += 180.0
    return phi % 360


def bruce2(x, y):
    """compute the PM angle between created by the x, y, cartesian coordinates

    Parameters
    ----------
    x, y position arrays

    Returns
    ---------
    PM angle between 0-360
    """
    rho = np.sqrt(x**2 + y**2)
    ra = x * np.pi / 180.0
    dec = y * np.pi / 180.0
    phi = np.arctan(ra/dec)
    phi *= 180.0 / np.pi  # to degrees
    if dec < 0.0:
        phi += 180.0
    return rho, phi % 360

# the 90 - phi statement rotates into dec as 0 frame

def pm2ang(x, y):
    """compute the PM angle between created by the x, y, cartesian coordinates

    Parameters:
    ----------
    x, y position (arrays)

    Returns
    ----------
    PM angle between 0-360
    """
    # rho = np.sqrt(x**2 + y**2)
    ra = x * np.pi / 180.0
    dec = y * np.pi / 180.0
    phi = np.arctan2(ra, dec)
    phi *= 180.0 / np.pi  # to degrees
    phi = 90 - phi
    return phi % 360

# pandas dataframe based routines
# compute proper motion angle, return in degrees, correct for quadrant


def compute_pm_angles_pandas(stars):
    """compute the PM angle between created by the x, y, PM values as an array
    Parameters: 
    ------------
    pandas Dataframe containing pmra and pmdec values

    Returns: 
    ------------
    PM angle between 0-360
    """
    tot = len(stars)
    ang = [0.0] * tot
    for index, row in stars.iterrows():
        pmra = row['pmra']
        pmdec = row['pmdec']
        ang[index] = bruce(pmra, pmdec)
    return ang
#
# regular methods
# works with Gaia job result list


def compute_pm_angles(stars):
    """compute the PM angle between pmra, pmdec  values as an array
    Parameters: 
    -----------
    array containing pmra and pmdec values

    Returns: 
    -----------
    PM angle between 0-360
    """
    tot = len(stars)
    ang = [0.0] * tot
    for i in range(tot):
        ang[i] = bruce(stars[i]['pmra'], stars[i]['pmdec'])
    return ang


def compute_pm_angles_and_rho(stars):
    """compute the PM angle between created by the PMRA, PMDec vectors, and the PM vector length.
    PM and rho values as arrays

    Parameters
    ----------
    array containing pmra and pmdec values

    Returns
    ----------
    PM angle array between 0-360 degrees
    """
    tot = len(stars)
    ang = [0.0] * tot
    rho = [0.0] * tot
    for i in range(tot):
        rho[i], ang[i] = bruce2(stars[i]['pmra'], stars[i]['pmdec'])
    return rho, ang


def compute_pm_angles_coords(pmra, pmdec):
    """compute the PM angle between created by the PMRA, PMDec vectors, and the PM vector length.

    Parameters
    ----------
    pmra, pmdec: arrays containing pmra and pmdec values

    Returns
    ----------
    PM angle array between 0-360 degress.
    """
    tot = len(pmra)
    ang = []
    i = 0
    for i in range(tot):
        _pmra = pmra[i]
        _pmdec = pmdec[i]
        ang.append(bruce(_pmra, _pmdec))
    return ang

# convert to polar coordinates from cartesion


def pol2cart(rho, phi):
    """convert polar to cartesian coordinates

    Parameters: array containing pmra and pmdec values

    Returns: PM angle array between 0-360, vector length array
    """
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return x, y
#
# SkyCoord setup and manipulation for Proper motion and Velocity calculations
#
# return a list of SkyCoords representing the individual stars from the dataframe, and
# containing only the positional and proper motion parameters


def make_skycoord_list(df):
    """make a list of SkyCoords representing the individual
    stars from the dataframe, containing only the positional and proper motion parameters

    Parameters
    ----------
    df: Gaia DataFrame object

    Returns
    ----------
     array containing individual SkyCoord objects, initialized with Proper Motion params.
    """
    star_coords = []
    star_coords = [SkyCoord(row['ra']*u.deg, row['dec']*u.deg, pm_ra_cosdec=row['pmra']*u.mas/u.yr,
                            pm_dec=row['pmdec']*u.mas/u.yr,
                            distance=row['parallax']*u.mas) for index, row in df.iterrows()]
    return np.asarray(star_coords)

# create a new DataFrame from the input DataFrame with a new column 'SkyCoord'
# containing a position, velocity and distance initialized SkyCoord object.
# Return the new dataframe.


def add_skycoords(df):
    """create a new DataFrame from the input DataFrame with a new column 'SkyCoord'
    containing a position, velocity and distance initialized SkyCoord object.

    Parameters
    ----------
    df: Gaia DataFrame object

    Returns
    ----------
    array containing individual SkyCoord objects, initialized with Proper Motion params.
    """
    star_coords = []
    df2 = df.copy()
    star_coords = make_skycoord_list(df)
    df2['SkyCoord'] = star_coords
    return df2

# create a SkyCoord object containing initialized position and PM data and return
# SINGLE SkyCoord object # using the input DataFrame object. This is used in the
# 3D separation calculations.

def make_skycoords_pm_single(df):
    """create a single SkyCoord object containing all of the Proper Motion
     parameters initialized from the input DataFrame (normal Gaia query result).

    Parameters:
    -----------
      df: normal Gaia DataFrame object

    Returns: 
    --------
    sc: a single SkyCoord object, containing the elements initialized with
    Proper Motion params.
    """
    ra = np.asarray(df['ra'])*u.deg
    dec = np.asarray(df['dec'])*u.deg
    pmra = np.asarray(df['pmra'])*u.mas/u.yr
    pmdec = np.asarray(df['pmdec'])*u.mas/u.yr
    # rad_velocity = np.asarray(df['radial_velocity'])*u.km/u.s
    dist = Distance(np.asarray(df['parallax'])*u.mas, allow_negative=True)
    sc = SkyCoord(ra, dec, distance=dist, pm_ra_cosdec=pmra, pm_dec=pmdec)
    return sc
#
# return a SkyCoord object with position and Proper Motion variables initialized with
# correct unit. Input is a stellar dataframe, index is the row index (0-based) for the
# star within the dataframe
#


def make_single_skycoord(stars, index):
    """create a SkyCoord object containing all of the Proper Motion parameters
    initialized from the input DataFrame (normal Gaia query result).

    Parameters:
        stars: normal Gaia DataFrame object

        index: the numeric index into the array that you want to retrieve and use

    Returns: a single SkyCoord object, representing the star of interest with
    Proper Motion params initialized.
    """
    pos = stars.iloc[index]
    ra = pos['ra']*u.deg
    dec = pos['dec']*u.deg
    pmra = pos['pmra']*u.mas/u.yr
    pmdec = pos['pmdec']*u.mas/u.yr
    distance = Distance(pos['parallax']*u.mas, allow_negative=True)
    sc = SkyCoord(ra, dec, pm_ra_cosdec=pmra, pm_dec=pmdec,
                  distance=distance, frame='icrs')
    return sc
#
# find a list of SkyCoords containing all stars within an input 3D distance cutoff (pc)


def extract_stars_by_separation(star, neighborhood, cutoff):
    """create a list of SkyCoords containing all stars within an input 3D distance cutoff (pc)

    Parameters:
        star: a SkyCoord object with Proper Motion data initialized, usually by
            make_single_skycoord().

        neighborhood: a normal Pandas dataframe with stellar data

        cutoff: distance cutoff (pc)
    Returns: a SkyCoord list containing the neighborhood of stars within a cutoff
    """
    # make a single SkyCoord object containing all stars in the neighborhood
    # we need this for the 3d search
    scoords = make_skycoords_pm_single(neighborhood)
    # make a SkyCoord array for each star in the neighborhood. we need this to
    sc_list = make_skycoord_list(neighborhood)
    neighbors = []
    sep = star.separation_3d(scoords)
    # build an array containing the SkyCoords of the stars based on the
    neighbors = [sc_list[index]
                 for index in range(0, len(sep)) if sep[index] < cutoff]
    return len(neighbors), np.asarray(neighbors)

#
# return a list of 3d separations from the input SkyCoord object and the neighborhood.


def compute_3d_separations(star, neighborhood):
    """create a list of 3d separations from the input SkyCoord object and the
    neighborhood of surrounding stars. This is a single object whose coordinates
    are bundled by make_single_skycoord.

    Parameters:
    ----------
       star: a SkyCoord object with Proper Motion data initialized, usually by
        make_single_skycoord().

       neighborhood: Pandas dataframe with stellar data as processed by
        read_gaia_stars(). This processing adds VX, VY, VZ, Distance and
        pm_angles columns to the raw Gaia query result file.

    Returns:
    --------
       list of 3D (pc) separations from the input star's position to each star in the
       neighborhood.
    """
    # make a single SkyCoord object containing the search neighborhood with Proper
    # Motion fields initialized
    scoords = make_skycoords_pm_single(neighborhood)
    sep = star.separation_3d(scoords)
    return sep

# return a list containg VX, VY, VZ from the input stellar dataframein the
# form: [[vx,vy,vz],[vx2,vy2,vz2]...] the pandas dataframe has no units.


def make_velocity_array(df):
    """create a list containing 3d velocity vectors

    Input:
    ------------
       df: a normal DataFrame with stellar data as processed by
        read_gaia_stars(). This processing adds VX, VY, VZ, Distance and
        pm_angles columns to the raw Gaia query result file.
    Returns:
    ------------
       a unitless array in the form: [[vx,vy,vz],[vx2,vy2,vz2]...]
    """
    vel_list = []
    vel_list = [[row['VX'], row['VY'], row['VZ']]
                for index, row in df.iterrows()]
    return vel_list


def make_velocity_array2d(df):
    """create a list containing 3d velocity vectors

    Input:
    ------------
       df: a normal DataFrame with stellar data as processed by
          read_gaia_stars(). This processing adds VX, VY, VZ, Distance and
          pm_angles columns to the raw Gaia query result file.
    Returns:
    ------------
       an unitless array in the form: [[vx,vy],[vx2,vy2]...]
    """
    vel_list = []
    vel_list = [[row['VX'], row['VY']]
                for index, row in df.iterrows()]
    return vel_list

# return a 3xn array of unitless values representing the X,Y,Z coordinates from an input
# SkyCoord list whose Proper motion values have been initialized, typically through the
# egs.make_skycoord_list function call.


def make_cartesian_list_vals(sc_list):
    """create a 3xn array of unitless values representing the X,Y,Z coordinates from
    an input SkyCoord list whose Proper motion values have been initialized, typically through the
    make_skycoord_list function call.

    Input
    ----------
       sc_list: list containing initialized SkyCoord objects
    Returns
    ----------
       an unitless array in the form: [[vx,vy,vz],[vx2,vy2,vz2]...]
    """
    cartpos_list = []
    cartpos_list = [[sc.cartesian.x.value, sc.cartesian.y.value, sc.cartesian.z.value]
                    for sc in sc_list]
    return cartpos_list

info = pd.DataFrame()

def make_cartesian_list(sc_list):
    """create a 3xn array of unitless values representing the X,Y,Z coordinates from
    an input SkyCoord list whose Proper motion values have been initialized, typically
    through the make_skycoord_list function call.

    Parameters:
    ---------
       sc_list: list containing initialized SkyCoord objects

    Returns:
    ----------
       an array in the form: [[x,y,z],[x2,y2,z2]...] without units
    """
    cartpos_list = []
    cartpos_list = [[sc.cartesian.x, sc.cartesian.y, sc.cartesian.z]
                    for sc in sc_list]
    return cartpos_list

# take an input SkyCoord list of positions and return X, Y, and Z vectors


def extract_cartesian_coords_vals(sc_list):
    """create 3 arrays of unitless values representing the X,Y,Z coordinates from
    an input SkyCoord list whose Proper motion values have been initialized, typically
    through the make_skycoord_list function call.

    Parameters:
    ------------
       sc_list: list containing initialized SkyCoord objects

    Returns:
    ------------
       X, Y, Z arrays without units
    """
    xpos_list = []
    ypos_list = []
    zpos_list = []

    for sc in sc_list:
        xpos_list.append(sc.cartesian.x.value)
        ypos_list.append(sc.cartesian.y.value)
        zpos_list.append(sc.cartesian.z.value)

    return xpos_list, ypos_list, zpos_list


def return_constraint(df, constraint, low, high):
    """returns dataframe containing rows from the input between the constraints
    low and high. The constraint is given in the string form 'Distance', 'VX'
    etc.
    Parameters:
    ------------
        df: star dataframe, low cutoff, high cutoff   
    Returns:
    ------------
    pmra, pmdec, ra, dec lists containing stars (rows) between the low and high 
    Parameters, inclusive
    """
    dfout = [row for index, row in df.iterrows()
             if row[constraint] > low and row[constraint] < high]
    return DataFrame(dfout)

# normalize the velocity vectors VX, VY, VZ contained in the dataframe
# return a new DataFrame with these velocities


def normalize_velocity(df):
    """returns a dataframe with vx, vy, vz normalized by their magnitudes

        Parameter:
        ----------
        df: star list as pandas DF, containing vx, vy, vz data

        Returns:
        ---------
        output: dataframe with normalized velocities
    """
    #length = []
    length = np.sqrt(df['VX']**2 + df['VY']**2 + df['VZ']**2)
    df_out = df.copy()
    df_out['VX'] = df['VX']/length
    df_out['VY'] = df['VY']/length
    df_out['VZ'] = df['VZ']/length
    return df_out

# the following 'slice' functions are obsolete. Use the return_constraint function instead.
# left in so as to not break old code. -egs-


def slice_range_coords_pandas(stars, low, high):
    """returns 4 lists containing rows from a pandas df satisfying the low and 
    high angle constraints.

    Parameters:
    ----------

    stars: pandas dataframe containing stellar data

    low: low angular cutoff

    high: angular cutoff

    Returns:
    ------------
    pmra, pmdec, ra, dec: lists containing stars (rows) satisfying the constraints, inclusive
    """
    pmra_list = []
    pmdec_list = []
    ra_list = []
    dec_list = []
    cnt = 0

    for _, row in stars.iterrows():
        pmra = row['pmra']
        pmdec = row['pmdec']
        ra = row['ra']
        dec = row['dec']
        if (bruce(pmra, pmdec) >= low) and (bruce(pmra, pmdec) <= high):
            pmra_list.append(pmra)
            pmdec_list.append(pmdec)
            ra_list.append(ra)
            dec_list.append(dec)
            cnt += 1
    return pmra_list, pmdec_list, ra_list, dec_list, cnt

#
# return ra, dec, distance lists from a pandas dataframe between low and high distances


def slice_distance_stars_pandas(stars, low, high):
    """create 4 lists containing rows from a pandas df satisfying the low and high
    angle constraints.
        input: star list, low angular cutoff, high angular cutoff

        output: pmra, pmdec, ra, dec lists containing stars (rows) satisfying the
        constraints, inclusive
    """
    dist_list = []
    ra_list = []
    dec_list = []
    cnt = 0

    for _, row in stars.iterrows():
        dist = row['Distance']
        ra = row['ra']
        dec = row['dec']
        if ((dist >= low) and (dist < high)):
            ra_list.append(ra)
            dec_list.append(dec)
            dist_list.append(dist)
            cnt += 1
    return np.asarray(ra_list), np.asarray(dec_list), np.asarray(dist_list), cnt

# return a new DataFrame containing stars between low and high parallax cutoff


# return a list containing radial velocities converted to GSR coordinates
# works with Gaia job result only
def get_radial_vals(star_list):
    """create a list containing non-zero Radial Velocities.

    The input star_list is a list object returned by job.results()

    Parameters
    ----------
    star_list : list returned by job.results()

    Returns
    -------
    v_gsr : list of coordinates containing non-zero radial_velicities
    """
    i = 0
    count = len(star_list)
    x = []

    # make a fancy progressbar to monitor progress since this this iterates over
    # typically long star lists
    widgets = [
        ' [', progressbar.Timer(), '] ', progressbar.Bar(
        ), ' (', progressbar.ETA(), ') ',
    ]

    with progressbar.ProgressBar(max_value=count, widgets=widgets) as bar2:
        for star in star_list:
            ra = star['ra'] * u.deg
            dec = star['dec'] * u.deg
            radial = star['radial_velocity'] * u.km/u.s
            star_icrs = ICRS(ra=ra, dec=dec, radial_velocity=radial)
            # convert to GSR from barycentric
            x.append(rv_to_gsr(star_icrs).value)
            count = count + 1
            i = i + 1
            bar2.update(i)
    return x

# Function converts the barycentric radial velocity to GSR coordinate frame


def rv_to_gsr(c, v_sun=None):
    """Transform a barycentric radial velocity to the Galactic Standard of Rest
    (GSR).

    The input radial velocity must be passed in as a ICRS object

    Parameters
    ----------
    c : `~astropy.coordinates.BaseCoordinateFrame` subclass instance
        The radial velocity, associated with a sky coordinate, to be
        transformed.
    v_sun : `~astropy.units.Quantity` (optional)
        The 3D velocity of the solar system barycenter in the GSR frame.
        Defaults to the same solar motion as in the
        `~astropy.coordinates.Galactocentric` frame.

    Returns
    -------
    v_gsr : `~astropy.units.Quantity`
        The input radial velocity transformed to a GSR frame.

    """
    if v_sun is None:
        v_sun = coord.Galactocentric.galcen_v_sun.to_cartesian()

    gal = c.transform_to(coord.Galactic)
    cart_data = gal.data.to_cartesian()
    unit_vector = cart_data / cart_data.norm()

    v_proj = v_sun.dot(unit_vector)
    return c.radial_velocity + v_proj

# compute 3d Velocities for a pandas dataframe and add 'VX', 'VY', 'VZ' to the dataframe
# from: http://www.astronexus.com/a-a/motions-long-term (km/sec). Ignores vR since we
# generally don't have radial velocities


def compute_3d_velocities_pandas(df):
    """compute 3D Velocities for a pandas DataFrame object and
    add 'VX', 'VY', 'VZ' to the dataframe (km/sec)
    from: http://www.astronexus.com/a-a/motions-long-term

    Parameters
    ----------
    df : normal Gaia pandas stellar DataFrame

    Returns
    -------
    stars : new DataFrame containing the velocity columns
    """
    ra = df['ra'] * np.pi / 180.0
    dec = df['dec'] * np.pi / 180.0

    # input is in mas, convert to arcseconds
    pmra = df['pmra'] / 1000.0
    pmdec = df['pmdec'] / 1000.0
    parallax = df['parallax'] / 1000.0

    # in general we don't have radial velocities for the stars
    vR = 0.0
    dist = 1.0 / parallax
    # make a copy of the input so we can safely add to it
    df2 = df.copy()
    # transverse distances in km/sec
    vtA = dist * pmra * 4.74
    vtD = dist * pmdec * 4.74
    # add VX, VY, VZ columns to the dataframe
    df2['VX'] = (vR * np.cos(dec) * np.cos(ra)) - \
        (vtA * np.sin(ra)) - (vtD * np.sin(dec) * np.cos(ra))
    df2['VY'] = (vR * np.cos(dec) * np.sin(ra)) + \
        (vtA * np.cos(ra)) - (vtD * np.sin(dec) * np.sin(ra))
    df2['VZ'] = vR * np.sin(dec) + vtD * np.cos(dec)
    return df2

# convert from km/sec to pc/year, ignores radial velocities


def compute_3d_velocities_pandas_pc(df):
    """compute 3D Velocities for a pandas DataFrame object and
    add 'VX', 'VY', 'VZ' to the dataframe (pc/yr)
    from: http://www.astronexus.com/a-a/motions-long-term

    Parameters
    ----------
    df : normal Gaia pandas stellar DataFrame

    Returns
    -------
    stars : new DataFrame containing the velocity columns
    """

    df2 = compute_3d_velocities_pandas(df)
    df2['VX'] = df2['VX'] / 977780.0
    df2['VY'] = df2['VY'] / 977780.0
    df2['VZ'] = df2['VZ'] / 977780.0
    return df2

# compute 3D velocities with input vectors for ra, dec, pmra, pmdec and parallax


def compute_3d_velocities(ra, dec, pmra, pmdec, parallax, vR=0.0):
    """compute 3D Velocities for a list of coordinates and proper motion parameters (km/sec)
    from: http://www.astronexus.com/a-a/motions-long-term

    Parameters
    ----------
    ra          : Right Ascension
    dec         : Declination
    pmra        : Proper Motion in RA
    pmdec       : Proper Motion in Dec
    parallax    : parallax
    vR          : radial velocity

    Returns
    -------
    vx, vy, vz  : components of velocity vectors in x, y, z
    """
    dist = 1.0 / parallax
    # transverse distances in km/sec
    vtA = dist * pmra * 4.74
    vtD = dist * pmdec * 4.74

    vx = (vR * np.cos(dec.rad) * np.cos(ra.rad)) - \
        (vtA * np.sin(ra.rad)) - (vtD * np.sin(dec.rad) * np.cos(ra.rad))
    vy = (vR * np.cos(dec.rad) * np.sin(ra.rad)) + \
        (vtA * np.cos(ra.rad)) - (vtD * np.sin(dec.rad) * np.sin(ra.rad))
    vz = vR * np.sin(dec.rad) + vtD * np.cos(dec.rad)
    return vx, vy, vz

#
# returns in pc/yr


def compute_3d_velocities_pc(ra, dec, pmra, pmdec, parallax, vR=0.0):
    """compute 3D Velocities for a list of coordinates and proper motion parameters (pc/yr)
    from: http://www.astronexus.com/a-a/motions-long-term

    Parameters
    ----------
    ra          : Right Ascension
    dec         : Declination
    pmra        : Proper Motion in RA
    pmdec       : Proper Motion in Dec
    parallax    : parallax
    vR          : radial velocity

    Returns
    -------
    vx, vy, vz  : components of velocity vectors in x, y, z (pc/yr)
    """
    vx, vy, vz = compute_3d_velocities(ra, dec, pmra, pmdec, parallax, vR)
    vx /= 977780.0
    vy /= 977780.0
    vz /= 977780.0

    return vx, vy, vz

# return a nx3 matrix of cartesian [x,y,z] coordinates from a DataFrame stellar data objects


def make_cartesian_coords_pandas(stars):
    """create a cnt x 3 list containing x, y, z unit sphere coords [[x,y,z]...]
    from a stellar DataFrame object

    Parameters
    ----------
    stars : DataFrame returned by job.results()

    Returns
    -------
     cartc : list of shape cnt,3 [[vx,vy,vz]...]
    """
    coord_list = []

    for _, row in stars.iterrows():
        x, y, z = ra_dec_to_xyz(row['ra'], row['dec'])
        vec = [x, y, z]
        coord_list.append(vec)
    return np.asarray(coord_list)


def ra_dec_to_xyz(ra, dec):
    """Convert ra & dec to Euclidean points on a unit sphere

    Parameters
    ----------
    ra, dec : ndarrays

    Returns
    x, y, z : ndarrays
    """
    sin_ra = np.sin(ra * np.pi / 180.)
    cos_ra = np.cos(ra * np.pi / 180.)

    sin_dec = np.sin(np.pi / 2 - dec * np.pi / 180.)
    cos_dec = np.cos(np.pi / 2 - dec * np.pi / 180.)

    xyz = [cos_ra * sin_dec, sin_ra * sin_dec, cos_dec]

    return np.asarray(xyz)

# function to do the rectangluar Gaia extraction


def extract_Gaia_stars(_RA, _DEC, width, height, project_dir, file_prefix, sig=-1.0, etype='',
                       verbose=False, Filtered=True):
    """Routine to extract stars from GAIA DR2 using a rectangular query centered at the ra, dec
    points, width x height, significance and an 'etype' which is either 'para', 'both' or
    'pm'. This refers to parameter the significance calculation is applied to, (either Parallax
    alone, both PMRA and PMDEC AND Parallax, or just PMRA and PMDEC). Do not specify this parameter
    at all to ignore. The pm angles and 3D velocities are calculated and written to a file whose
    path and prefix is specified. The significance calculation is discussed in:
        https://arxiv.org/pdf/0711.3593v1.pdf and described in the ADQL cookbook:

        https://www.gaia.ac.uk/data/gaia-data-release-1/adql-cookbook

    Parameters
    ----------
    _RA, _DEC : strings for the query center point RA and Dec, eg: "20h10m20s", "37deg10m,39s"

    width, height: strings for total search width, eg: "30m"

    sig: significance filter level - 0 = no filter

    project_dir: string representing the path to the directory for the output files, 
                 eg: "/MIRA/data/PM/"

    file_prefix: string for filename prefix, eg: "pelican_pm_"

    etype: string for significance screen type, eg: 'para' for parallax screening, 'both' for
    PM and Parallax screening, omit for just PM screening (the default)

    verbose: if True, print the query string

    Filtered: if False, no filtering for null parameters

    Returns
    ---------
     starcount: stars read

     seconds: elapsed time for query in seconds

     vel_filename: output velocity filename written
    """
    home_dir = ut.get_home_dir()
    full_path = home_dir + project_dir
    project_name = full_path + file_prefix

    # create directory path if neeeded
    ut.ensure_dir_path(full_path)

   # convert to astropy.Coordinate objects for easier manipulation
    ra = Longitude(_RA, unit=u.hourangle)
    dec = Latitude(_DEC, unit=u.degree)

    ra_width = Longitude(width, unit=u.degree)
    dec_width = Latitude(height, unit=u.degree)

    # make the overall query strings by printing into the sel_str_box templates
    # select significance for both PM and Parallax, Parallax alone, PM alone, respectively
    _ra = ra.deg
    _dec = dec.deg

    q_both = str.format(sel_str_box_both_sig, _ra, _dec, _ra,
                        _dec, ra_width.deg, dec_width.deg, sig, sig)
    q_para = str.format(sel_str_box_para_sig, _ra, _dec, _ra,
                        _dec, ra_width.deg, dec_width.deg, sig)
    q_pm = str.format(sel_str_box_pm_sig, _ra, _dec, _ra,
                      _dec, ra_width.deg, dec_width.deg, sig)
    q_all = str.format(sel_str_box_all, _ra, _dec, _ra,
                       _dec, ra_width.deg, dec_width.deg)
    q_filt = str.format(sel_str_box_filt, _ra, _dec, _ra,
                        _dec, ra_width.deg, dec_width.deg)

    _sig = sig
    _etype = etype

    # i'm turning this off for now due to the changing queries 10/5/19 -egs-
    _sig = -1.0
    _etype = 'none'
    
    if _sig > 0.0:
        # format the filenames using the format: projectname_ra1_dec1_rawidthxdecwidth_deg.csv
        # raw output filename for the GAIA query
        # expand create_vel_filename to include these cases
        data_filename = str.format("{}{}_{}_{}x{}_{}sig_{}.csv", project_name, _RA, _DEC,
                                   round(ra_width.deg, 2), round(dec_width.deg, 2), sig, etype)
        if etype == 'para':
            q = q_para
        elif etype == 'both':
            q = q_both
        else:
            q = q_pm
    else:
        data_filename = str.format("{}{}_{}_{}x{}_all.csv", project_name, _RA, _DEC,
                                   round(ra_width.deg, 2), round(dec_width.deg, 2))

    vel_filename = ut.make_output_filename(project_dir, file_prefix, _RA, _DEC,
                                           ra_width.deg, dec_width.deg, _sig, _etype)
    if Filtered:
        q = q_filt
    else:
        q = q_all
    start = time()

    if verbose:
        print("Query is:", q)
        print("Velocity filename is: ", vel_filename)
        print("Data filename is: ", data_filename)

    job = Gaia.launch_job_async(query=q, dump_to_file=True, output_file=data_filename,
                                output_format='csv')
    stars = job.get_results()
    star_count = len(stars)

    # PM and parallax are in mas
    pmra = stars['pmra']
    pmdec = stars['pmdec']
    # Proper Motion length vector
    pmvec = np.sqrt(pmra**2 + pmdec**2)

    # convert to arcsecs to do the 3d Veclocity and distance calculations
    pmra = stars['pmra'] / 1000.0
    pmdec = stars['pmdec'] / 1000.0
    parallax = stars['parallax'] / 1000.0

    # convert to long and lat to get the units right
    ra = Longitude(stars['ra'], unit=u.deg)
    dec = Latitude(stars['dec'], unit=u.deg)

    # compute the velocity vectors
    # 3d velocity components initialized to empty lists for the calculation
    vx = []
    vy = []
    vz = []
    vx, vy, vz = compute_3d_velocities_pc(ra, dec, pmra, pmdec, parallax)

    # compute the pm angles
    angles = []
    angles = compute_pm_angles(stars)

    # read the data back into a dataframe
    data = read_csv(data_filename)
    # add new columns to the data for the computed 3D velocities
    data['VX'] = vx
    data['VY'] = vy
    data['VZ'] = vz
    data['PM_Angles'] = angles
    data['PM_Vec'] = pmvec
    distances = []
    distances = 1.0 / parallax
    data['Distance'] = distances
    # write the file with the additional columns to the vel_filename output file
    data.to_csv(vel_filename)
    end = time()
    return star_count, end-start, vel_filename

# read in the .csv file represented by the input filename
# fail if it doesn't exist. This does no processing, merely reading.

def read_Gaia_stars(filename, verbose=False):
    """Read a .csv file containing stellar data and create a DataFrame object \
        indexed by source_id.
        
        Parameter: 
        ----------
         filename: string representing full filename path
         verbose: True if one wants to see the full file path on error
         
        Returns: 
        --------
         df: DataFrame containing the stellar data with index set to source_id
    """
    if os.path.isfile(filename):
        df = pd.read_csv(filename)
        df.set_index('source_id')
        return df
    else:
        if verbose:
            print("File:", filename, " does not exist! Can't return DataFrame()")
        return False
# end of file
