# Introduction

Jupyter notebooks for Astronomical analysis. Written by Eric G. Suchanek, Ph.D.,
Monterey Institute for Research in Astronomy (MIRA).

## Getting Started

These notebooks have been developed using *Anaconda* and has been tested under OSX 10.15,
Ubuntu 18.04, Windows 7 Professional and Windows 10. I have used Microsoft's
*Visual Studio Code* IDE for development and have been in general quite pleased with it.
It plays well with Anaconda and can correctly use the virtual python environment described below.

*NB:* As with any Python project there are a number of libraries that must be installed. Anaconda/conda greatly facilitates this process. This will be described more fully below:

## Base Environment Installation

1. Install Anaconda (<http://anaconda.org>)
    - create a new environment using python 3.7
    - activate the environment

2. Install VisualStudio Code (<https://code.visualstudio.com/download>)

3. Build the environment:
    - Install the folloowing libraries from within the environment created above:
        - astropy
        - astroquery
        - astroml
        - ipyvol
        - plotly
        - plotly_express
        - scikit-learn
        - pip

4. Use pip to install:
    - ipyaladin

5. Latest release
    - the current master branch has the latest working release

## Running the notebooks in Binder

Click below to launch the Jupyter Notebook browser in *Binder*

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/suchanek/PMA.git/master)
