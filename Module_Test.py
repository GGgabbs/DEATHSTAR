#
#         # # # # #
#       #           #
#     #           #   #
#   #             #   # #
#   #               #   #
#   #                   #
#   # # # # # # # # # # #
#   #                   #
#     #               #
#       #           #
#         # # # # #
#
# Filename:     Module_Test.py
# Created by:   Gabrielle Ross
# Last updated: 2/10/2024
# Github:       https://github.com/GGgabbs/DEATHSTAR/tree/main
#
# DEATHSTAR:    A system for confirming planets and identifying
#               false positive signals in TESS data using 
#               ground-based time domain surveys
#
# Purpose:      Runs automatically along with setup.sh to ensure
#               all dependencies are properly installed within
#               the conda virtual environment DEATHSTAR_environment


# Modules
# Arrays and handling data
import numpy as np
import pandas as pd
import math

# Astronomical files
from astropy.io import fits
from astroquery.mast import Mast
from astropy.table import Table

# Plotting data
import matplotlib
from matplotlib import pyplot as plot
from mpl_toolkits.axes_grid1 import make_axes_locatable # Color log
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes # Zoom
from mpl_toolkits.axes_grid1.inset_locator import mark_inset # Zoom

# Pulling data
import requests
import json
from urllib.parse import quote as urlencode
import urllib.request
from ztfquery import query

# Accessing local files
import os
import sys
import pickle

# Creating apertures
from photutils.aperture import EllipticalAperture
from photutils.aperture import EllipticalAnnulus
from photutils.aperture import ApertureStats
from photutils.aperture import aperture_photometry
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
from astropy.wcs import WCS # Converting between pixels and WCS
from astropy.coordinates import Angle # Converting from degrees for aperture
import astropy.coordinates
from mpfit import mpfit # Making a gaussian fit for the aperture

from IPython.display import clear_output # Clearing output to not lag notebook

import pdb # For viewing and stopping the program

# Fixing time offsets
from astropy.time import Time
import datetime
import time
import calendar

from PyAstronomy import pyasl


print("All modules are correctly installed!\nPlease press ENTER key to close this window and end program.")
