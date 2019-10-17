# Introduction

Jupyter notebooks for Astronomical analysis. Written by Eric G. Suchanek, Ph.D.,
Monterey Institute for Research in Astronomy (MIRA).

## Getting Started

These notebooks have been developed and using *Anaconda* and have been tested under OSX 10.15,
Ubuntu 18.04, Windows 7 Professional and Windows 10. I have used Microsoft's
*Visual Studio Code* IDE for development and Jupyter Notebook for interactive visualization.

*NB:* As with any Python project there are a number of libraries that must be installed. Anaconda/conda greatly facilitates this process. This will be described more fully below:

## Virtual Environment Installation/Creation

1. Install Anaconda (<http://anaconda.org>)
    - Create a new environment using python 3.7
    - Activate the environment

2. Build the environment
    - Manually
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
    - Using the *environment.yml* spec file (alternative, from a shell prompt)
     - clone this repo into your working repository directory
     - cd into the resulting directory
     - execute: *conda env create --file=environment.yaml*

3. Use pip to install:
    - ipyaladin

## Latest release

The current master branch has the latest release

## Running the notebooks in Binder

Click below to launch the Jupyter Notebook browser in Binder.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/suchanek/PMA.git/master)
