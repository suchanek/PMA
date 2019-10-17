# To add a new cell, type '#%%'
# To add a new markdown cell, type '#%% [markdown]'
#%% Change working directory from the workspace root to the ipynb file location. Turn this addition off with the DataScience.changeDirOnImportExport setting
# ms-python.python added
import os
try:
	os.chdir(os.path.join(os.getcwd(), 'NAmerican'))
	print(os.getcwd())
except:
	pass
#%%
from IPython import get_ipython

#%% [markdown]
# # NGC 7000 Cluster Analysis
# 
# *Eric G. Suchanek, Ph.D. 9/26/19*
# 
# This notebook shows extracts from the Gaia DR2 repository that were created via my iPython notebook: <a href="./ExtractStars_batch_all.ipynb">ExtractStars_batch_all.ipynb</a>. The code blocks utilize
# the python-based <a href="https://docs.astropy.org/en/stable/index.html">*astropy*</a> object and method library and its associated <a href="https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html">*SkyCoord* class</a> to compute 3D stellar separations from positional and parallax data extracted from Gaia. See: <a href=https://het.as.utexas.edu/HET/Software/Astropy-1.0/coordinates/matchsep.html> Astropy Separations </a> for more information on the technique used. Clustering techniques are described in the <a href="https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fclusterdata.html#scipy.cluster.hierarchy.fclusterdata"> SciPy documentation</a>.
# 
# <i>NB: the code relies on the routines stored in the file: <a href="./egs.py">egs.py</a> which should be located in the same directory as this notebook.</i>
#%% [markdown]
# ## Preamble and globals. Always run the following cell first.

#%%
# Author: Eric G. Suchanek, Ph.D.
# Written in association with the Monterey Institute for Research in Astronomy
# License: BSD

# Important libraries and global variables that the following cells rely on. Run this first. -egs-
# Last modification: 9/5/19
# Suppress warnings. Comment this out if you wish to see the warning messages

import warnings
warnings.filterwarnings('ignore')

import numpy as np

import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

# plotting and visualization libraries
import plotly.offline
import plotly_express as px

# libraries for interactive plots
import ipywidgets as widgets 
from ipywidgets import Layout, Box, widgets

# pandas routines
from pandas import DataFrame

# astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table

# Aladin
import ipyaladin as ipyal
from ipyaladin import Aladin

# misc
from time import time
import random

from scipy.cluster.hierarchy import fclusterdata

# the mirapy module must be in the python library path!
from mirapy import *
from mirapy.utils import pprint_elapsed, get_home_dir, create_vel_filename
import mirapy.egs as egs
from mirapy.egs import make_velocity_array2d


# Standard plotly imports
#import chart_studio.plotly as py
import plotly.graph_objs as go
from plotly.offline import iplot, init_notebook_mode

# Using plotly + cufflinks in offline mode
# We go offline with cufflinks first, since this is all local
#
init_notebook_mode(connected=False)


#%%
# Global variables used by the following cells. If you export the cells to straight Python remember to
# include the appropriate libraries and variables from this area of the notebook!

# specify the absolute path prefix for the location of the extracted stars
# will need to use this when reading as well!
# these files are stored in the following sub directory under my home tree
# you have to use escape on Windows due to the assinine path formats

project_dir = "/MIRA/data/PM_subset/"
file_prefix = "NGC7000PM_"
# Center of NGC 7000
ra_center = "20h58m47s" 
dec_center = "44d19m48s"

# A string name for the object we're analyzing. Used in the graphs
obj_name = "NGC7000"

# extent of the current analysis
_extent = .75

# one of linux, imac, parallels, pc, laptop, bootcamp
_sys = imac_ubuntu

# make the filename for the given extent. These were created by extract_batch_stars_all
# with no statistical constraints. The only constraints were for non-null PMRA and 
# PMDec and Parallax > 0.

vel_filename = create_vel_filename(project_dir, file_prefix, ra_center, 
                                   dec_center, _extent, _extent)
#vel_filename = "/home/suchanek/repos/Mira1/NAmerican/ngc7000_30sec_circle.csv"

start = time()

# read the data into a data frame
# this will be our working global dataframe for the following cells.

STARS = egs.read_Gaia_stars(vel_filename)
if STARS is False:
    print("Can't open file: " + vel_filename)
    Star_count = -1
else:
    print('Read in:',len(STARS),'stars', 'for extent', _extent, 'degrees.')
    t = Table.read(vel_filename,format="ascii.csv")
    Star_count = len(STARS)

pprint_elapsed(start,_sys)

# it's simple to display the target using the current vel files FOV
pov = Aladin(target=obj_name, fov=_extent)
pov


#%%
pov.add_table(t)


#%%
# make a plotly express 3d scatterplot and histogram from the stellar data stored in STARS
# plot 3D scatterplot of the normally computed PM angle parameter:

title_str = str.format("PM Angle for {} stars at {}, {} Extent: {}x{} deg", 
                       len(STARS), ra_center, dec_center, _extent, _extent)


fig = px.scatter_3d(STARS, x="ra", y="dec", z="PM_Angles", title=title_str, 
                    color="PM_Angles", size="PM_Angles",
                    template='plotly_dark')
fig.show()

# plot a histogram of PM angles:
fig2 = px.histogram(STARS, x='PM_Angles', title=title_str, marginal='violin',
                    template='plotly_dark')
fig2.show()


#%%


#%% [markdown]
# ## Clustering trials using various metrics and algorithms. 
#%% [markdown]
# ### Attempt to cluster by PM Magnitude 
# $$mag = \sqrt{pmra^2 + pmdec^2}$$
# *The following cells operate on the slice of stars whose PM Angle is in the peak, shown above*

#%%
# slice out the stars in the peak range by pm angle
# calculate 3d separations
# cluster by 3d separations
low = 211.0
high = 222.5
narrow_pm = egs.return_constraint(STARS, 'PM_Angles', low, high)
print("Slicing between", low, "to", high, "resulted in:", len(narrow_pm), "stars")


#%%
#
# Compute number of clusters by PM Magnitude
# Iterate over a number of cluster distance cutoffs and determine the number of clusters found,
# plot them on a scatterplot
#
# make a list for the cutoff we want to analyze
start = time()
cutoffs = [i for i in np.arange(.1, 31, .25)]

# make a dataframe by defining a dict for our two columns
pm = DataFrame()
pm['PM Vec'] = narrow_pm['PM_Vec']

# create a list of the maximum number of clusters by iterating over the distance cutoffs
# using list comprehension:
max_clusters = [fclusterdata(pm, cut,'distance').max() for cut in cutoffs]

# make a dataframe for plotting
data = DataFrame()
data['Cutoff (PM Mag)'] = cutoffs
data['Clusters'] = max_clusters
title_str = str.format("Number of clusters by PM Mag for {} stars at {}, {} Extent: {}x{} deg.", 
                       len(narrow_pm), ra_center, dec_center, _extent, _extent)

# make a plotly express scatterplot from the dataframe
fig = px.scatter(data, x="Cutoff (PM Mag)", y="Clusters", template='plotly_dark', title=title_str)
fig.show()

# plot a histogram of PM angles:
fig2 = px.histogram(narrow_pm, x='PM_Vec', title=title_str, marginal='violin',
                    template='plotly_dark')
fig2.show()
pprint_elapsed(start, _sys)

#%% [markdown]
# 
#%% [markdown]
# ## Pick a cutoff and plot cluster membership

#%%
#
# plot cluster membership by Position for a given cutoff

pm = DataFrame()
pm['PM Vec'] = narrow_pm['PM_Vec']

# pick a cutoff based on the above result
cut = .02
clust = fclusterdata(pm, cut,'distance')

# maximum number of clusters
max_clust = clust.max()
data = DataFrame()
data['ra'] = narrow_pm['ra']
data['dec'] = narrow_pm['dec']
data['clust'] = clust

title_str = str.format("{} clusters by PM Mag for {} stars at {}, {} Extent: {}x{} deg, cutoff {}", 
                       max_clust, len(narrow_pm), ra_center, dec_center, _extent, _extent, cut)

# make a plotly express scatterplot from the dataframe
fig = px.scatter(data, x="ra", y="dec",  title=title_str, color='clust',
                 color_continuous_scale=px.colors.sequential.Viridis,
                template='plotly_dark')
fig.show()

data2 = DataFrame()
data2['Cluster'] = clust
fig2 = px.histogram(data2, x='Cluster', title=title_str, marginal='violin',
                    template='plotly_dark')
fig2.show()


#%%


#%% [markdown]
# ## Distributions

#%%



#%%
# Look at timing for building 3d separation distance matrices
# !!!
lower = 300.0
start = time()
star_count_list = []

time_list = []
cutoff_list = [400, 750, 1000, 1250, 1500, 2000, 2500, 3000]

for high in cutoff_list:
    start2 = time()
    close = egs.return_constraint(STARS,'Distance', lower, high)
    star_count = len(close)
    sep_list = [egs.compute_3d_separations(egs.make_single_skycoord(close,i), close) for i in range(star_count)]
    print('-- Returned:', star_count, 'stars')
    pprint_elapsed(start2, arch=_sys)
    print(' ')
    star_count_list.append(star_count)
    time_list.append(time() - start2)

print("Overall elapsed time:")
pprint_elapsed(start, _sys)

df = DataFrame()
df['Starcount'] = star_count_list
df['Elapsed'] = time_list
df['Cutoffs'] = cutoff_list

title_str = "Star count by upper distance cutoff with lower distance: " + str(lower) + " Extent:" + str(_extent)
title_str2 = "Elapsed time to compute 3d seps for lower distance: " + str(lower) + " Extent:" + str(_extent)

# make a plotly express scatterplot from the dataframe
fig = px.scatter(df, x="Cutoffs", y="Starcount", template='plotly_dark', title=title_str)
fig2 = px.scatter(df, x="Cutoffs", y="Elapsed", template='plotly_dark', title=title_str2)

fig.show()
fig2.show()


#%%


#%% [markdown]
# ## Create a cone plot where the directions are the 3D velocity vectors

#%%
import plotly.graph_objects as go
import pandas as pd
from mirapy.egs import make_velocity_array2d

low = 500
high = 1000

close = egs.return_constraint(STARS,'Distance', low, high)
star_count = len(close)

vels = make_velocity_array2d(close)
norms = np.sqrt(close['VX']**2 + close['VY']**2)

titl_str = str.format("3D Velocities, {} stars at {}, {}, {}x{} deg,\r Distance between: {} and {}pc", 
                      star_count, ra_center, dec_center, _extent, _extent, low, high)

Data = [go.Cone(x=close['ra'], y=close['dec'], z=close['Distance'], u=close['VX'], v=close['VY'],w=norms,
                colorscale='Blues', sizemode="scaled", sizeref=1)]
    
fig = go.Figure(data = Data,
                layout_title_text = titl_str,
                layout_template = 'plotly_dark'
               )

#fig = go.Figure(data=data, layout=layout)
fig.update_layout(title_text=titl_str, xaxis_title="RA", yaxis_title="Dec", title_font_size=18, width=900, height=900)
fig.show()


#%%
size_list = [str(size)+"m" for size in range(5,125,5)]
print(size_list)


#%%



#%%
from astroquery.sdss import SDSS
from astropy import coordinates as coords
pos = coords.SkyCoord(ra_center, dec_center, frame='icrs')
xid = SDSS.query_region(pos, spectro=True)
print(xid)


