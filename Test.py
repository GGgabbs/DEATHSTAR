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
# Filename:     DEATHSTAR_Example.py
# Created by:   Gabrielle Ross
# Last updated: 2/20/2024
# Github:       https://github.com/GGgabbs/DEATHSTAR/tree/main
# 
# Please cite our paper if you use this code!
# ADS link: https://ui.adsabs.harvard.edu/abs/2023MNRAS.tmp.3722R/abstract
#
# DEATHSTAR:    A system for confirming planets and identifying
#               false positive signals in TESS data using 
#               ground-based time domain surveys
#
# Purpose:      Runs a test target on already confirmed TOI 4148.01
#               / TIC 137157546 as it has a fast download time
#               to test image download and package installation



# Modules
from DEATHSTAR import setup # Importing DEATHSTAR and all its functionality


setup(137157546, signal_tic_ID = 137157545) # Example of running full program all the way through