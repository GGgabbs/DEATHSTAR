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
# Filename:     DEATHSTAR.py
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
# Purpose:      Runs automatically along with setup.sh to ensure
#               all dependencies are properly installed within
#               the conda virtual environment DEATHSTAR_environment



# Modules
import os.path
from Extracting_Lightcurves import setup as get_lightcurves
from Plotting import setup as get_plots



# Private functions
def setup(
    tic_ID # Only required argument
    , is_ATLAS = False, is_overwrite = False # Arguments for both programs
    , is_testing = False, is_processing = True, sigma_coefficient = 1, new_fitter = 0 # Additional arguments for Extracting_Lightcurves.py
    , comparison_removal = [], signal_tic_ID = 0, revised_period = 0, size = 100, xlim = "zoomed", frame_number = 1, is_plotting = True, is_showing_index = True, is_saving = True, is_plotting_MAD_vs_MAG = True, is_lcbin = False, is_period_revision = False # Additional arguments for Plotting.py
    ):
    
    if os.path.isfile("TICs/" + str(tic_ID) + "/" + "All_" + str(tic_ID) + "_Star_Data.csv"): # If the star data exists
        if is_overwrite:
            get_lightcurves(tic_ID, is_testing, is_processing, is_ATLAS, sigma_coefficient = sigma_coefficient, new_fitter = new_fitter) # TIC ID, is testing, is processing, and is ATLAS. Can add the sigma coefficient (sigma_coefficient) and the TIC ID (new_fitter) of new fitter
        else:
            print("The extracted light curves .csv file already exists. Please turn on is_overwrite = True in order to rerun this target.")
    else:
        get_lightcurves(tic_ID, is_testing, is_processing, is_ATLAS, sigma_coefficient = sigma_coefficient, new_fitter = new_fitter) # TIC ID, is testing, is processing, and is ATLAS (vs. ZTF). Can add sigma coefficient (sigma_coefficient), and TIC ID of new fitter (new_fitter)
    
    if not is_testing:
        get_plots(tic_ID, is_ATLAS, comparison_removal = comparison_removal, signal_tic_ID = signal_tic_ID, revised_period = revised_period, size = size, xlim = xlim, frame_number = frame_number, is_plotting = is_plotting, is_showing_index = is_showing_index, is_saving = is_saving, is_plotting_MAD_vs_MAG = is_plotting_MAD_vs_MAG, is_lcbin = is_lcbin, is_period_revision = is_period_revision, is_done = False) # TIC ID and is ATLAS (vs. ZTF). Can add comparisons to remove (comparison_removal), the TIC ID of signal (signal_tic_ID), revised period if found (revised_period), size of radius for stars to plot (size), x-limits of the period plots (xlim), the frame number of the reference image to display (frame_number), is plotting (is_plotting), is showing the index or the number of working frames for each star within the reference plot (is_showing_index), if saving the images (is_saving), if plotting noise diagnostic (is_plotting_MAD_vs_MAG), if plotting binned light curves (is_lcbin), is preforming and plotting period revision (is_period_revision), and if the function will rerun or is its last run (is_done)
