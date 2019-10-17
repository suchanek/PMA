# Introduction

Jupyter notebooks for Astronomical analysis. Written by Eric G. Suchanek, Ph.D.,
Monterey Institute for Research in Astronomy (MIRA).

## Getting Started

These notebooks have been developed using *Anaconda* and has been tested under OSX 10.15,
Ubuntu 18.04, Windows 7 Professional and Windows 10. I have used Microsoft's
*Visual Studio Code* IDE for development and have been in general quite pleased with it.
It plays well with Anaconda and can correctly use the virtual python environment described below.

*NB:* As with any Python project there are a number of libraries that must be installed. Anaconda/conda greatly facilitates this process. This will be described more fully below:

## Virtual Environment Installation/Creation

1. Install Anaconda (<http://anaconda.org>)
    - Create a new environment using python 3.7
    - Activate the environment

2. Build the environment (manually):
    - Install the following libraries from within the environment created above:
        - astropy
        - astroquery
        - astroml
        - ipyvol
        - pandas
        - plotly
        - plotly_express
        - scikit-learn
        - pip

3. Use pip to install:
    - ipyaladin

## Latest release

The current master branch has the latest release

## Running the notebooks in Binder

Click below to launch the Jupyter Notebook browser in Binder.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/suchanek/PMA.git/master)
