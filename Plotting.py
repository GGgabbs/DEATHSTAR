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
# Filename:     Plotting.py
# Created by:   Gabrielle Ross
# Last updated: 2/16/2024
# Github:       https://github.com/GGgabbs/DEATHSTAR/tree/main
# 
# Please cite our paper if you use this code!
# ADS link: https://ui.adsabs.harvard.edu/abs/2023MNRAS.tmp.3722R/abstract
#
# DEATHSTAR:    A system for confirming planets and identifying
#               false positive signals in TESS data using 
#               ground-based time domain surveys
#
# Purpose:      Accesses the All_[TIC ID]_Star_Data.csv,
#               plots the reference image and light curve sheets,
#               can plot binned and period revsion comparisons,
#               and saves images in preferred format and full report



# Modules
# Arrays and handling data
import numpy as np
import pandas as pd

# Plotting data
import matplotlib
from matplotlib import pyplot as plot
from mpl_toolkits.axes_grid1 import make_axes_locatable # Color log
from astropy.io import fits
from astropy.wcs import WCS
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar # Scalebar
from matplotlib.transforms import Bbox # Scalebar position
import matplotlib.font_manager as fm # Font properties

# Calculations
import math
from scipy import stats
from scipy.interpolate import make_interp_spline # For smoothing the errorbar line data

# Pulling data
import requests
import urllib.request
import urllib3
import json

# Converting time
from astropy.time import Time
import datetime
import time
import calendar

# Local files
import sys # Terminating program
import os
from pathlib import Path

from fpdf import FPDF # Saving data

from IPython.display import clear_output # Clearing output to not lag notebook

import lcbin as lcb # Binning data

from astropy.timeseries import BoxLeastSquares # Period revision



# Custom look and feel variables
font_size = 20
target_star_color = "#D81B60" # Pink
comparison_star_color = "#1E88E5" # Light blue
saturated_star_color = "#e8b007" # Orange
dim_star_color = "#00725c" # Teal green
too_few_points_color = "white" # White
transit_boundary_color = "m" # Magenta
plots_per_page = 12



# Private functions
def get_TOI_info(mini_name_identifier): # Getting the info from tev
    url = "https://exofop.ipac.caltech.edu/tess/target.php?id=" + mini_name_identifier + "&json"
    response = urllib.request.urlopen(url)
    data = json.loads(response.read())
    return {"ra": str(data["coordinates"]["ra"]), "dec": str(data["coordinates"]["dec"]), "period": float(data["planet_parameters"][1]["per"]), "period_uncertainty": float(data["planet_parameters"][1]["per_e"]), "depth": float(data["planet_parameters"][1]["dep_p"]), "epoch": float(data["planet_parameters"][1]["epoch"]), "epoch_uncertainty": float(data["planet_parameters"][1]["epoch_e"]), "duration": float(data["planet_parameters"][1]["dur"]), "Tmag": float(data["magnitudes"][0]["value"])} # ra, dec, period (days), period_uncertainty, transit_depth (ppm / parts per million), transit_epoch (BJD time of the center of the transit), transit_epoch_uncertainty, transit_duration (hours), target_magnitude (Tmag / TESS magnitude)

# Time conversions
def display_period_uncertainty(mini_all_star_dataframe, toi_info): # Period uncertainty
    uncertainty_on_time_of_future_transits = []
    for i in (min, max): # Getting the first and last image date
        reference_image_date = float(i(list(mini_all_star_dataframe["time"]))) # Tt
        number_of_transits = abs(toi_info["epoch"] - reference_image_date) / toi_info["period"] # n
        uncertainty_on_time_of_future_transit = math.sqrt(((number_of_transits * toi_info["period_uncertainty"]) ** 2) + ((toi_info["epoch_uncertainty"]) ** 2))
        uncertainty_on_time_of_future_transits.append(uncertainty_on_time_of_future_transit)
    
    print("TOI period duration: " + str(toi_info["duration"] / 24)) # Convert from hours to days
    
    if toi_info["old_period"] == toi_info["period"]: # Display normal 5-year propogated uncertainty
        print("Uncertainty on time of future transit: " + str(max(uncertainty_on_time_of_future_transits)))
    else: # Display uncertainty of the edited period
        toi_info["period_uncertainty"] = 12 * toi_info["period"] * toi_info["duration"] / abs(float(max(list(mini_all_star_dataframe["time"]))) - float(min(list(mini_all_star_dataframe["time"]))))
        print("Uncertainty on time of future transit: " + str(toi_info["period_uncertainty"]))
    
    return toi_info["period_uncertainty"]

def convert_BJD_and_JD(days, toi_info, BJD_to_JD): # Converting between BJD and JD time
    URL = "https://astroutils.astronomy.osu.edu/time/convert.php?JDS=" + str(days) + "&RA=" + str(toi_info["ra"]) + "&DEC=" + str(toi_info["dec"]) + "&FUNCTION=utc2bjd"
    if BJD_to_JD: # Converting from BJD to JD-UTC
        URL = URL.replace("utc2bjd", "bjd2utc")
    
    response = requests.get(URL, verify = False)
    content = response.content.decode("utf-8") # Returned data
    
    return content.replace("\n", "")

def get_transit_BJD_times_offset(toi_info): # Getting the start and end times of transit from the epoch time and the transit duration
    julian_days = convert_BJD_and_JD(toi_info["epoch"], toi_info, True)
    return float(float(julian_days) - (toi_info["duration"] / (2 * 24))), float(float(julian_days) + (toi_info["duration"] / (2 * 24)))

def display_rounded_uncertainties(number, uncertainty): # Rounding the period and epoch. Take the uncertainty, get 2 actual numbers after all the 0s (if <1), then count all the decimal places with those 2 rounded and put that number of decimal places in the transit epoch and period. So if uncertainty is 0.000076543, then rounded is 0.000077, and use 6 decimal places for the actual number
    number_display = number
    uncertainty_display = uncertainty
    
    if str(uncertainty_display)[0] == "0":
        number_array = []
        for i in range(len(str(uncertainty_display)[2:])):
            if sum(x != "0" for x in number_array) == 2: #if sum(x is not "0" for x in number_array) == 2:
                break
            number_array.append(str(uncertainty_display)[2:][i])
        
        if len(number_array) < len(str(uncertainty_display)[2:]): # Need to round
            format_string = "{0:." + str(len(number_array)) + "f}"
            uncertainty_display = float(format_string.format(uncertainty_display, len(number_array) + 1))
        else:
            uncertainty_display = float("0." + ''.join(number_array))
        
        format_string = "{0:." + str(len(number_array)) + "f}"
        number_display = float(format_string.format(number, len(number_array)))
    else:
        number_display = float("{:.3f}".format(number))
        uncertainty_display = int(uncertainty_display)
    
    return number_display, uncertainty_display

def make_bjd_time_mod_period(times, mini_period, mini_transit_epoch): # Mod the times by the period
  tmodn = (times % mini_period) - (mini_transit_epoch % mini_period)
  tmodn = tmodn + (mini_period * (tmodn < -0.5 * mini_period)) - (mini_period * (tmodn > 0.5 * mini_period))
  return tmodn


# Helpers
def factor_int(n): # Getting the size of the plot grid. This function takes in a number and outputs 2 numbers that form the most square-like grid possible
    val = math.ceil(math.sqrt(n))
    val2 = int(n / val)
    while (val2 * val) != float(n):
        val -= 1
        val2 = int(n / val)
    return max(val, val2), min(val, val2)

def get_comparison_sums(dataframe_column, ics): # Combining all of the brightnesses of the comparison stars
    sum = 0
    for i in ics:
        sum += np.array(dataframe_column[i])
    
    return sum


# Populating the sheet subplot
def plot_ZTF_scatter(mini_ax, column_name, immediate_comparison_stars, horizontal, markersize, toi_info, cone, clean_dataframe, lists):
    # Separating each filter
    ZTF_g_dataframe = clean_dataframe[(clean_dataframe["image_filter"] == "ZTF_g") | (clean_dataframe["image_filter"] == "ZTF g")]
    ZTF_r_dataframe = clean_dataframe[(clean_dataframe["image_filter"] == "ZTF_r") | (clean_dataframe["image_filter"] == "ZTF r")]
    ZTF_i_dataframe = clean_dataframe[(clean_dataframe["image_filter"] == "ZTF_i") | (clean_dataframe["image_filter"] == "ZTF i")]
    
    # Getting the y limits
    star_g_with_comparison = (np.array(ZTF_g_dataframe[column_name]) / get_comparison_sums(ZTF_g_dataframe[immediate_comparison_stars].copy(), immediate_comparison_stars)) / np.median((np.array(ZTF_g_dataframe[column_name]) / get_comparison_sums(ZTF_g_dataframe[immediate_comparison_stars].copy(), immediate_comparison_stars)))
    median_absolute_deviation_g = stats.median_abs_deviation(star_g_with_comparison)
    star_r_with_comparison = (np.array(ZTF_r_dataframe[column_name]) / get_comparison_sums(ZTF_r_dataframe[immediate_comparison_stars].copy(), immediate_comparison_stars)) / np.median((np.array(ZTF_r_dataframe[column_name]) / get_comparison_sums(ZTF_r_dataframe[immediate_comparison_stars].copy(), immediate_comparison_stars)))
    median_absolute_deviation_r = stats.median_abs_deviation(star_r_with_comparison)
    star_i_with_comparison = (np.array(ZTF_i_dataframe[column_name]) / get_comparison_sums(ZTF_i_dataframe[immediate_comparison_stars].copy(), immediate_comparison_stars)) / np.median((np.array(ZTF_i_dataframe[column_name]) / get_comparison_sums(ZTF_i_dataframe[immediate_comparison_stars].copy(), immediate_comparison_stars)))
    median_absolute_deviation_i = stats.median_abs_deviation(star_i_with_comparison)
    
    lists["MAD"]["g"].append(median_absolute_deviation_g / np.median(star_g_with_comparison))
    lists["MAD"]["r"].append(median_absolute_deviation_r / np.median(star_r_with_comparison))
    lists["MAD"]["i"].append(median_absolute_deviation_i / np.median(star_i_with_comparison))
    lists["MAG"].append(cone.at[list(cone["ID"]).index(int(column_name.split("_")[0])), "Tmag"])
    if column_name in lists["comparison_stars"]:
        lists["MAD_comparison"]["g"].append(median_absolute_deviation_g / np.median(star_g_with_comparison))
        lists["MAD_comparison"]["r"].append(median_absolute_deviation_r / np.median(star_r_with_comparison))
        lists["MAD_comparison"]["i"].append(median_absolute_deviation_i / np.median(star_i_with_comparison))
        lists["MAG_comparison"].append(cone.at[list(cone["ID"]).index(int(column_name.split("_")[0])), "Tmag"])
        lists["MAD_vs_MAG_comparison_names"].append(column_name)
    
    max_median_absolute_deviation = max((1 / (median_absolute_deviation_g ** 2)), (1 / (median_absolute_deviation_r ** 2)), (1 / (median_absolute_deviation_i ** 2)))
    if max_median_absolute_deviation == (1 / (median_absolute_deviation_i ** 2)):
        mini_ax.set_ylim([min([(horizontal - 0.15), (((np.median(star_i_with_comparison) - (5 * median_absolute_deviation_i)) / np.median(star_i_with_comparison)) - 0.15)]), (((np.median(star_i_with_comparison) + (5 * median_absolute_deviation_i)) / np.median(star_i_with_comparison)) + 0.1)])
    elif max_median_absolute_deviation == (1 / (median_absolute_deviation_r ** 2)):
        mini_ax.set_ylim([min([(horizontal - 0.15), (((np.median(star_r_with_comparison) - (5 * median_absolute_deviation_r)) / np.median(star_r_with_comparison))) - 0.15]), (((np.median(star_r_with_comparison) + (5 * median_absolute_deviation_r)) / np.median(star_r_with_comparison)) + 0.1)])
    else:
        mini_ax.set_ylim([min([(horizontal - 0.15), (((np.median(star_g_with_comparison) - (5 * median_absolute_deviation_g)) / np.median(star_g_with_comparison))) - 0.15]), (((np.median(star_g_with_comparison) + (5 * median_absolute_deviation_g)) / np.median(star_g_with_comparison)) + 0.1)])
    
    # Polygon uncertainty on transit time
    plot.fill(-toi_info["period_uncertainty"], mini_ax.get_ylim()[0], toi_info["period_uncertainty"], mini_ax.get_ylim()[1], color = "gray", alpha = 0.5, zorder = 5)
    
    max_median_absolute_deviation = max((1 / (median_absolute_deviation_g ** 2)), (1 / (median_absolute_deviation_r ** 2)))
    # Actual plotting code
    if not ZTF_g_dataframe.empty:
        g_alpha = (1 / (median_absolute_deviation_g ** 2)) / max_median_absolute_deviation
        if g_alpha < 0.15:
            g_alpha = 0.15
        mini_ax.plot(make_bjd_time_mod_period(np.array(ZTF_g_dataframe["time"]), toi_info["period"], toi_info["epoch"]), star_g_with_comparison, ".", color = "b", alpha = g_alpha, markersize = markersize)
    if not ZTF_r_dataframe.empty:
        r_alpha = (1 / (median_absolute_deviation_r ** 2)) / max_median_absolute_deviation
        if r_alpha < 0.15:
            r_alpha = 0.15
        mini_ax.plot(make_bjd_time_mod_period(np.array(ZTF_r_dataframe["time"]), toi_info["period"], toi_info["epoch"]), star_r_with_comparison, ".", color = "r", alpha = r_alpha, markersize = markersize)
    if not ZTF_i_dataframe.empty:
        mini_ax.plot(make_bjd_time_mod_period(np.array(ZTF_i_dataframe["time"]), toi_info["period"], toi_info["epoch"]), star_i_with_comparison, ".", color = "brown", alpha = 1, markersize = markersize)
    
    if column_name == str(toi_info["signal_tic_ID"]) + "_brightness":
        star_with_comparison = list(star_g_with_comparison) + list(star_r_with_comparison) + list(star_i_with_comparison)
        all_times = list(np.array(ZTF_g_dataframe["time"])) + list(np.array(ZTF_r_dataframe["time"])) + list(np.array(ZTF_i_dataframe["time"]))
        signal_brightnesses = {"flux": star_with_comparison, "time": all_times, "filter": (["g"] * len(star_g_with_comparison) + ["r"] * len(star_r_with_comparison) + ["i"] * len(star_i_with_comparison))}
        pd.DataFrame(signal_brightnesses).to_csv("TICs/" + str(toi_info["tic_ID"]) + "/" + str(toi_info["signal_tic_ID"]) + "_Signal_Brightnesses.csv")
        signal_y_bottom, signal_y_top = mini_ax.get_ylim()
        toi_info["signal_y_bottom"] = signal_y_bottom
        toi_info["signal_y_top"] = signal_y_top

def plot_ATLAS_scatter(mini_ax, column_name, immediate_comparison_stars, horizontal, markersize, toi_info, clean_dataframe):
    # Separating each filter
    ZTF_o_dataframe = clean_dataframe[(clean_dataframe["image_filter"] == "o")]
    ZTF_c_dataframe = clean_dataframe[(clean_dataframe["image_filter"] == "c")]
    
    # Getting the y limits
    star_o_with_comparison = (np.array(ZTF_o_dataframe[column_name]) / get_comparison_sums(ZTF_o_dataframe[immediate_comparison_stars].copy(), immediate_comparison_stars)) / np.median((np.array(ZTF_o_dataframe[column_name]) / get_comparison_sums(ZTF_o_dataframe[immediate_comparison_stars].copy(), immediate_comparison_stars)))
    median_absolute_deviation_o = stats.median_abs_deviation(star_o_with_comparison)
    star_c_with_comparison = (np.array(ZTF_c_dataframe[column_name]) / get_comparison_sums(ZTF_c_dataframe[immediate_comparison_stars].copy(), immediate_comparison_stars)) / np.median((np.array(ZTF_c_dataframe[column_name]) / get_comparison_sums(ZTF_c_dataframe[immediate_comparison_stars].copy(), immediate_comparison_stars)))
    median_absolute_deviation_c = stats.median_abs_deviation(star_c_with_comparison)
    
    MAD_list["o"].append(median_absolute_deviation_o / np.median(star_o_with_comparison))
    MAD_list["c"].append(median_absolute_deviation_c / np.median(star_c_with_comparison))
    MAG_list.append(cone.at[list(cone["ID"]).index(int(column_name.split("_")[0])), "Tmag"])
    if column_name in comparison_stars_list:
        MAD_comparison_list["o"].append(median_absolute_deviation_o / np.median(star_o_with_comparison))
        MAD_comparison_list["c"].append(median_absolute_deviation_c / np.median(star_c_with_comparison))
        MAG_comparison_list.append(cone.at[list(cone["ID"]).index(int(column_name.split("_")[0])), "Tmag"])
        MAG_vs_MAD_comparison_name_list.append(column_name)
    
    max_median_absolute_deviation = max((1 / (median_absolute_deviation_c ** 2)), (1 / (median_absolute_deviation_o ** 2)))
    if max_median_absolute_deviation == (1 / (median_absolute_deviation_o ** 2)):
        mini_ax.set_ylim([min([(horizontal - 0.15), (((np.median(star_o_with_comparison) - (5 * median_absolute_deviation_o)) / np.median(star_o_with_comparison)) - 0.15)]), (((np.median(star_o_with_comparison) + (5 * median_absolute_deviation_o)) / np.median(star_o_with_comparison)) + 0.1)])
    else:
        mini_ax.set_ylim([min([(horizontal - 0.15), (((np.median(star_c_with_comparison) - (5 * median_absolute_deviation_c)) / np.median(star_c_with_comparison))) - 0.15]), (((np.median(star_c_with_comparison) + (5 * median_absolute_deviation_c)) / np.median(star_c_with_comparison)) + 0.1)])
    
    # Actual plotting code
    if not ZTF_o_dataframe.empty:
        o_alpha = (1 / (median_absolute_deviation_o ** 2)) / max_median_absolute_deviation
        if o_alpha < 0.15:
            o_alpha = 0.15
        mini_ax.plot(make_bjd_time_mod_period(np.array(ZTF_o_dataframe["time"]), period, transit_epoch), star_o_with_comparison, ".", color = "r", alpha = o_alpha, markersize = markersize)
    if not ZTF_c_dataframe.empty:
        c_alpha = (1 / (median_absolute_deviation_c ** 2)) / max_median_absolute_deviation
        if c_alpha < 0.15:
            c_alpha = 0.15
        mini_ax.plot(make_bjd_time_mod_period(np.array(ZTF_c_dataframe["time"]), period, transit_epoch), star_c_with_comparison, ".", color = "b", alpha = c_alpha, markersize = markersize)
    
    if column_name == str(signal_tic_ID.split()[1]) + "_brightness":
        star_with_comparison = list(star_o_with_comparison) + list(star_c_with_comparison)
        all_times = list(np.array(ZTF_o_dataframe["time"])) + list(np.array(ZTF_c_dataframe["time"]))
        signal_brightnesses = {"flux": star_with_comparison, "time": all_times}
        pd.DataFrame(signal_brightnesses).to_csv("TICs/" + name_identifier + "/" + str(signal_tic_ID.split()[1]) + "_Signal_Brightnesses.csv")

def change_ax_colors(mini_ax, color): # Highlighting certain subplots (like target, saturated, dim, etc.)
    mini_ax.xaxis.label.set_color(color)
    for side in ["top", "right", "left", "bottom"]:
        mini_ax.spines[side].set_color(color)
    mini_ax.tick_params(axis = "x", colors = color)
    mini_ax.tick_params(axis = "y", colors = color)

def periodigram_ax(figure, mini_ax, toi_info, markersize = 6):
    signal_brightnesses = pd.read_csv("TICs/" + str(toi_info["tic_ID"]) + "/" + str(toi_info["signal_tic_ID"]) + "_Signal_Brightnesses.csv")
    top = 1.03
    bottom_int = int(toi_info["h"] * 100)
    
    all_periods = []
    all_powers = []
    '''
    for i in range(bottom_int - 5, bottom_int + 5):
        bottom = i / 100
        if bottom < 1:
            good = np.argwhere(np.logical_and(np.array(signal_brightnesses["flux"]) < top, np.array(signal_brightnesses["flux"] > bottom)))
            
            # Getting the periodigram
            periods = np.linspace(((-10 * toi_info["period_uncertainty"]) + toi_info["period"]), ((10 * toi_info["period_uncertainty"]) + toi_info["period"]), 1000)
            model = BoxLeastSquares(np.ravel(np.array(signal_brightnesses["time"])[good]), np.ravel(np.array(signal_brightnesses["flux"])[good]))
            bls_power = model.power(periods, 0.2)
            
            # Plotting the periodigram
            all_periods += list(periods)
            all_powers += list(bls_power.power)
    '''
    
    bottom = bottom_int / 100
    good = np.argwhere(np.logical_and(np.array(signal_brightnesses["flux"]) < top, np.array(signal_brightnesses["flux"] > bottom)))
    
    # Getting the periodigram
    periods = np.linspace(((-30 * toi_info["period_uncertainty"]) + toi_info["period"]), ((30 * toi_info["period_uncertainty"]) + toi_info["period"]), 1000)
    model = BoxLeastSquares(np.ravel(np.array(signal_brightnesses["time"])[good]), np.ravel(np.array(signal_brightnesses["flux"])[good]))
    bls_power = model.power(periods, 0.2)
    
    # Plotting the periodigram
    all_periods += list(periods)
    all_powers += list(bls_power.power)
    
    # Plotting the periodigram
    mini_ax.plot(all_periods, all_powers, "-", color = "blue", alpha = 0.5, markersize = markersize)
    mini_ax.tick_params(axis = "x", labelsize = font_size)
    mini_ax.tick_params(axis = "y", labelsize = font_size)
    
    # Updating my variables
    toi_info["old_period"] = toi_info["period"]
    toi_info["period"] = all_periods[np.argmax(all_powers)]
    print("Revised period: " + str(toi_info["period"]) + " days")

def ax_setup(figure, mini_ax, column_index, toi_info, all_star_dataframe, cone, clean_dataframe, lists, mini_is_ATLAS, markersize = 6):
    figure.text(0.5, 0.01, "Time from mid-transit\nFolded on period " + str(toi_info["period_display"]) + " " + u"\u00B1" + " " + str(toi_info["period_uncertainty_display"]) + " Epoch " + str(toi_info["epoch_display"]) + " " + u"\u00B1" + " " + str(toi_info["epoch_uncertainty_display"]) + " in BJD - 2457000", ha = "center", fontsize = font_size, fontweight = "normal") # Tess Julian Days = Julian Days - 2457000
    figure.text(0.01, 0.5, "Relative Flux", va = "center", rotation = "vertical", fontsize = font_size, fontweight = "normal")
    figure.subplots_adjust(bottom = 0.11, top = 0.9, left = 0.06, right = 0.975, wspace = 0.15, hspace = 0.6)
    
    if toi_info["xlim"] == "zoomed":
        mini_ax.set_xlim([- 4 * (abs((toi_info["transit_start"] % toi_info["period"]) - (toi_info["transit_end"] % toi_info["period"]))), (4 * abs((toi_info["transit_start"] % toi_info["period"]) - (toi_info["transit_end"] % toi_info["period"])))]) # 4x the duration of the transit
    elif toi_info["xlim"] != "all": # If xlim == "all" then let auto-zoom
        mini_ax.set_xlim([-toi_info["xlim"], toi_info["xlim"]])
    
    column_name = str(list(clean_dataframe.columns[2:])[column_index])
    immediate_comparison_stars = lists["comparison_stars"].copy()
    if column_name in lists["comparison_stars"]: # Removing the comparison star if it is the 'target' for the subplot so it is not dividing by itself
        if column_name in immediate_comparison_stars:
            immediate_comparison_stars.remove(column_name)
    
    # Getting a horizontal line where the depth of the transit should be
    faint_star_magnitude = cone.at[int(np.where(cone["ID"] == int(str((list(clean_dataframe.columns[2:])[column_index]).split("_")[0])))[0]), "Tmag"] # Where signal is actually coming from
    faint_star_flux = 2.512 ** (toi_info["Tmag"] - faint_star_magnitude) # Converting from magnitude to flux
    horizontal = ((1 - toi_info["depth"]) * (faint_star_flux + 1) - 1) / faint_star_flux # Assuming the normal flux is 1
    if column_index == toi_info["target_index"]: # If the transit is on-target
        horizontal = 1 - toi_info["depth"] # The one above is a buffer and is deeper than intended but this is the actual depth
    mini_ax.axhline(horizontal, linestyle = "--", color = transit_boundary_color, zorder = 10)
    if column_index == toi_info["signal_index"]:
        toi_info["h"] = horizontal
    
    # Subplot title
    subplot_x_axis_title = "T" + str(column_index) + ", TIC " + str((list(clean_dataframe.columns[2:])[column_index]).split("_")[0])
    subplot_x_axis_title += r", $\rho$=" + str(round(np.mean(list(all_star_dataframe[str(str((list(clean_dataframe.columns[2:])[column_index]).split("_")[0]) + "_distance_to_real_target")])), 1)) + "''" # Convert from degrees to arcseconds
    subplot_x_axis_title += ",\n" + r"$PA=" + str(round(float(str(stats.mode(list(all_star_dataframe[str(str((list(clean_dataframe.columns[2:])[column_index]).split("_")[0]) + "_position_angle")]), keepdims = True).mode[0]).replace("rad", "")) * 180 / np.pi, 1)) + r"^\circ$" # Convert from decimal radians to degrees
    subplot_x_axis_title += ", $\Delta$$T=" + str(round(toi_info["Tmag"] - faint_star_magnitude, 1)) + "$"
    subplot_x_axis_title += ", $T=" + str(round(faint_star_magnitude, 1)) + "$"
    
    if mini_is_ATLAS:
        plot_ATLAS_scatter(mini_ax, column_name, immediate_comparison_stars, horizontal, markersize, toi_info, cone, clean_dataframe, lists)
    else:
        plot_ZTF_scatter(mini_ax, column_name, immediate_comparison_stars, horizontal, markersize, toi_info, cone, clean_dataframe, lists)
    
    if column_index == toi_info["target_index"]: # Highlight target
        change_ax_colors(mini_ax, target_star_color)
    elif list(clean_dataframe.columns[2:])[column_index] in lists["saturated_stars"]: # Highlight saturated stars
        change_ax_colors(mini_ax, saturated_star_color)
    elif (not list(clean_dataframe.columns[2:])[column_index] in lists["comparison_stars"]) and (not list(clean_dataframe.columns[2:])[column_index] in lists["comparison_stars"]): # Highlight dim stars
        change_ax_colors(mini_ax, dim_star_color)
    else:
        change_ax_colors(mini_ax, comparison_star_color)
    
    # Adding the start, middle, and end transit lines
    mini_ax.axvline(make_bjd_time_mod_period(toi_info["transit_start"], toi_info["period"], toi_info["epoch"]), linestyle = "--", color = transit_boundary_color, zorder = 10)
    mini_ax.axvline(make_bjd_time_mod_period(toi_info["transit_end"], toi_info["period"], toi_info["epoch"]), linestyle = "--", color = transit_boundary_color, zorder = 10)
    
    mini_ax.xaxis.set_label_position("top")
    mini_ax.set_xlabel(subplot_x_axis_title, fontsize = font_size)
    mini_ax.tick_params(axis = "x", labelsize = font_size)
    mini_ax.tick_params(axis = "y", labelsize = font_size)
    
    return column_index + 1

def make_subplots(mini_starting_index, total_number, toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list):
    plot_y_grid, plot_x_grid = factor_int(total_number)
    
    column_index = mini_starting_index
    if plot_x_grid == 1:
        figure, ax = plot.subplots(plot_y_grid, sharex = True)
        if plot_y_grid == 1:
            column_index = ax_setup(figure, ax, column_index, toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list["ATLAS"])
        else:
            for y_ax in range(plot_y_grid):
                column_index = ax_setup(figure, ax[y_ax], column_index, toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list["ATLAS"])
        
        if is_list["saving"]:
            save_plot(("TICs/" + str(toi_info["tic_ID"]) + "/TIC_" + str(toi_info["tic_ID"]) + "_Lightcurve" + str(int(math.ceil(mini_starting_index / plots_per_page)) + 1)), figure, False, False)
    else:
        figure, ax = plot.subplots(plot_y_grid, plot_x_grid, sharex = True)
        for y_ax in range(plot_y_grid): 
            for x_ax in range(plot_x_grid):
                column_index = ax_setup(figure, ax[y_ax, x_ax], column_index, toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list["ATLAS"])
        
        if is_list["saving"]:
            save_plot(("TICs/" + str(toi_info["tic_ID"]) + "/TIC_" + str(toi_info["tic_ID"]) + "_Lightcurve" + str(int(math.ceil(mini_starting_index / plots_per_page)) + 1)), figure, False, False)
    
    make_plot_fullscreen()

def plot_MAD_vs_MAG(toi_info, lists, is_list, size, frame_number, comparison_removal, is_done):
    figure, ax = plot.subplots()
    if is_list["ATLAS"]:
        ax.plot(lists["MAG"], lists["MAD"]["o"], ".", color = "orange", label = "Orange Filter")
        ax.plot(lists["MAG"], lists["MAD"]["c"], ".", color = "cyan", label = "Cyan Filter")
        ax.plot(lists["MAG_comparison"], lists["MAD_comparison"]["o"], "*", color = "orange", label = "Orange Filter", markersize = 5)
        ax.plot(lists["MAG_comparison"], lists["MAD_comparison"]["c"], "*", color = "cyan", label = "Cyan Filter", markersize = 5)
        for i in range(len(lists["MAD_vs_MAG_comparison_names"])):
            ax.text(lists["MAG_comparison"][i], lists["MAD_comparison"]["c"][i] + 0.05, lists["MAD_vs_MAG_comparison_names"][i].split("_")[0])
            if lists["MAD_comparison"]["c"][i] > 0.1:
                print('comparison_stars_list.remove("' + str(lists["MAD_vs_MAG_comparison_names"][i]) + '")')
    else:
        ax.plot(lists["MAG"], lists["MAD"]["g"], ".", color = "b", label = "Green Filter")
        ax.plot(lists["MAG"], lists["MAD"]["r"], ".", color = "r", label = "Red Filter")
        ax.plot(lists["MAG"], lists["MAD"]["i"], ".", color = "brown", label = "IR Filter")
        ax.plot(lists["MAG_comparison"], lists["MAD_comparison"]["g"], "*", color = "b", label = "Green Filter", markersize = 5)
        ax.plot(lists["MAG_comparison"], lists["MAD_comparison"]["r"], "*", color = "r", label = "Red Filter", markersize = 5)
        ax.plot(lists["MAG_comparison"], lists["MAD_comparison"]["i"], "*", color = "brown", label = "IR Filter", markersize = 5)
        
        comparison_removal_values = []
        for i in range(len(lists["MAD_vs_MAG_comparison_names"])):
            ax.text(lists["MAG_comparison"][i], lists["MAD_comparison"]["g"][i] + 0.05, lists["MAD_vs_MAG_comparison_names"][i].split("_")[0])
            if (lists["MAD_comparison"]["g"][i] > 0.1) or (lists["MAD_comparison"]["r"][i] > 0.1) or (lists["MAD_comparison"]["i"][i] > 0.1):
                comparison_removal_values.append(max(lists["MAD_comparison"]["g"][i], lists["MAD_comparison"]["r"][i], lists["MAD_comparison"]["i"][i])) # Adding the highest of the 3 filters
            else:
                comparison_removal_values.append(0) # Just to keep the indicies correct
        
        if not is_done:
            ax.set_ylim([0, 1.5])
            ax.axhline(0.1, linestyle = "--", color = transit_boundary_color, zorder = 10)
            make_plot_fullscreen()
            
            if comparison_removal_values == [0] * len(comparison_removal_values): # Re-run without the worst comparison star
                make_plot_fullscreen()
                
                if is_list["saving"]:
                    save_plot(("TICs/" + str(toi_info["tic_ID"]) + "/TIC_" + str(toi_info["tic_ID"]) + "_MAD_vs_MAG_FINAL"), figure, False, False)
                
                setup(toi_info["tic_ID"], is_list["ATLAS"], comparison_removal = comparison_removal, signal_tic_ID = toi_info["signal_tic_ID"], revised_period = toi_info["period"], size = size, xlim = toi_info["xlim"], frame_number = frame_number, is_plotting = is_list["plotting"], is_showing_index = is_list["showing_index"], is_saving = is_list["saving"], is_plotting_MAD_vs_MAG = False, is_lcbin = is_list["lcbin"], is_period_revision = is_list["period_revision"], is_done = True)
            else:
                new_comparison_removal_index = comparison_removal_values.index(max(comparison_removal_values)) # Get the highest filter value of any of the comparison stars and remove just that star
                new_comparison_removal = int(lists["MAD_vs_MAG_comparison_names"][new_comparison_removal_index].split("_")[0])
                print("Removing " + str(new_comparison_removal))
                comparison_removal.append(new_comparison_removal)
                
                if is_list["saving"]:
                    save_plot(("TICs/" + str(toi_info["tic_ID"]) + "/TIC_" + str(toi_info["tic_ID"]) + "_MAD_vs_MAG_" + str(comparison_removal).replace("[", "").replace("]", "")), figure, False, False)
                    
                setup(toi_info["tic_ID"], is_list["ATLAS"], comparison_removal = comparison_removal, signal_tic_ID = toi_info["signal_tic_ID"], revised_period = toi_info["period"], size = size, xlim = toi_info["xlim"], frame_number = frame_number, is_plotting = is_list["plotting"], is_showing_index = is_list["showing_index"], is_saving = is_list["saving"], is_plotting_MAD_vs_MAG = is_list["plotting_MAD_vs_MAG"], is_lcbin = is_list["lcbin"], is_period_revision = is_list["period_revision"])

def plot_signal_vs_target(toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list):
    if toi_info["tic_ID"] == toi_info["signal_tic_ID"]:
        figure, ax = plot.subplots(1, 1) # 1 row, 1 columns (1 cell total)
        _ = ax_setup(figure, ax, toi_info["target_index"], toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list["ATLAS"], 12)
        
        # Changing the colors and spacing that are different just for this half-page plot
        change_ax_colors(ax, "black")
        
        figure.subplots_adjust(bottom = 0.15, top = 0.9, left = 0.115, right = 0.975, wspace = 0.1, hspace = 0.35)
        
        if is_list["saving"]:
            save_plot(("TICs/" + str(toi_info["tic_ID"]) + "/Target_Lightcurve"), figure, False, True)
    
        figure.set_size_inches((9.6, 6.8), forward = True) # If True, overwrite existing window
    else:
        figure, ax = plot.subplots(1, 2) # 1 row, 2 columns (2 cells total)
        _ = ax_setup(figure, ax[0], toi_info["target_index"], toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list["ATLAS"], 12)
        _ = ax_setup(figure, ax[1], toi_info["signal_index"], toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list["ATLAS"], 12)
        
        # Changing the colors and spacing that are different just for this half-page plot
        change_ax_colors(ax[0], "black")
        change_ax_colors(ax[1], "black")
        
        figure.subplots_adjust(bottom = 0.15, top = 0.9, left = 0.075, right = 0.975, wspace = 0.1, hspace = 0.35)
        
        if is_list["saving"]:
            save_plot(("TICs/" + str(toi_info["tic_ID"]) + "/Target_vs_Signal_Lightcurve"), figure, False, True)
    
        figure.set_size_inches((19.2, 6.8), forward = True) # If True, overwrite existing window

def plot_lightcurves(toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list, size, frame_number, comparison_removal, is_done): # Plotting the scatter lightcurves in multiple plots for better legibility
    starting_index = 0
    total_indexes = int(len(clean_dataframe.columns) - 2)
    mod_length = int(total_indexes % plots_per_page)

    # Plotting each 5x4 grid
    if not (total_indexes - mod_length) == 0:
        for i in range(int((total_indexes - mod_length) / plots_per_page)):
            make_subplots(starting_index, plots_per_page, toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list)
            starting_index += plots_per_page

    # Plotting leftover from ____ % plots_per_page
    if not mod_length == 0:
        closest_square = math.floor(math.sqrt(mod_length)) ** 2
        if not closest_square == 1:
            make_subplots(starting_index, closest_square, toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list)
            if not (mod_length - closest_square) == 0:
                make_subplots((starting_index + closest_square), (mod_length - closest_square), toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list)
        else:
            make_subplots(starting_index, mod_length, toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list)
    
    if is_list["plotting"]:
        plot_signal_vs_target(toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list)
        
        if is_list["plotting_MAD_vs_MAG"]:
            plot_MAD_vs_MAG(toi_info, lists, is_list, size, frame_number, comparison_removal, is_done)
        
        plot.show()

def lcbin(toi_info, is_list):
    print("Binning lightcurves...")
    signal_brightnesses = pd.read_csv("TICs/" + str(toi_info["tic_ID"]) + "/" + str(toi_info["signal_tic_ID"]) + "_Signal_Brightnesses.csv") # Can add as a ZTF .csv and an ATLAS one for plotting ATLAS vs. ZTF as it is the same plot
    
    figure, ax = plot.subplots(1, 2) # 1 row, 2 columns (2 cells total)
    
    g = signal_brightnesses[(signal_brightnesses["filter"] == "g")]
    r = signal_brightnesses[(signal_brightnesses["filter"] == "r")]
    i = signal_brightnesses[(signal_brightnesses["filter"] == "i")]
    t = np.array(make_bjd_time_mod_period(np.array(signal_brightnesses["time"]), toi_info["period"], toi_info["epoch"]))
    f = np.array(signal_brightnesses["flux"])
    tbin, fbin, ebin, mask = lcb.lcbin(t, f, int(3 * toi_info["period"] / (toi_info["duration"] / 24))) # The greater number, the noisier the points. 5 is how many points within the transit I want BUT IT HAS TO BE ODD and hope that there isn't a large gap between -pi/2 and pi/2
    
    point_size = 20
    
    ax[0].plot(make_bjd_time_mod_period(np.array(g["time"]), toi_info["period"], toi_info["epoch"]), g["flux"], '.', alpha = 0.4, markersize = point_size, color = "b")
    ax[0].plot(make_bjd_time_mod_period(np.array(r["time"]), toi_info["period"], toi_info["epoch"]), r["flux"], '.', alpha = 0.4, markersize = point_size, color = "r")
    ax[0].plot(make_bjd_time_mod_period(np.array(i["time"]), toi_info["period"], toi_info["epoch"]), i["flux"], '.', alpha = 0.4, markersize = point_size, color = "brown")
    ax[1].plot(t, f, '.', alpha = 0.2, markersize = 20, color = "grey")
    ax[1].errorbar(tbin, fbin, yerr = ebin, fmt = "o", color = "r") # Scatter of averaged bins
    
    # Line connecting error dots
    ax[1].plot(tbin, fbin, color = "r")
    
    if toi_info["xlim"] == "zoomed":
        ax[0].set_xlim([- 4 * (abs((toi_info["transit_start"] % toi_info["period"]) - (toi_info["transit_end"] % toi_info["period"]))), (4 * abs((toi_info["transit_start"] % toi_info["period"]) - (toi_info["transit_end"] % toi_info["period"])))]) # 4x the duration of the transit
        ax[1].set_xlim([- 4 * (abs((toi_info["transit_start"] % toi_info["period"]) - (toi_info["transit_end"] % toi_info["period"]))), (4 * abs((toi_info["transit_start"] % toi_info["period"]) - (toi_info["transit_end"] % toi_info["period"])))]) # 4x the duration of the transit
        ax[0].hlines(toi_info["h"], xmin = - 4 * (abs((toi_info["transit_start"] % toi_info["period"]) - (toi_info["transit_end"] % toi_info["period"]))), xmax = (4 * abs((toi_info["transit_start"] % toi_info["period"]) - (toi_info["transit_end"] % toi_info["period"]))), linestyle = "--", color = "magenta", zorder = 10) #77
        ax[1].hlines(toi_info["h"], xmin = - 4 * (abs((toi_info["transit_start"] % toi_info["period"]) - (toi_info["transit_end"] % toi_info["period"]))), xmax = (4 * abs((toi_info["transit_start"] % toi_info["period"]) - (toi_info["transit_end"] % toi_info["period"]))), linestyle = "--", color = "magenta", zorder = 10) #77
    elif toi_info["xlim"] != "all": # If xlim == "all" then let auto-zoom
        ax[0].set_xlim([-toi_info["xlim"], toi_info["xlim"]])
        ax[1].set_xlim([-toi_info["xlim"], toi_info["xlim"]])
        ax[0].hlines(toi_info["h"], xmin = -toi_info["xlim"], xmax = toi_info["xlim"], linestyle = "--", color = "magenta", zorder = 10) #77
        ax[1].hlines(toi_info["h"], xmin = -toi_info["xlim"], xmax = toi_info["xlim"], linestyle = "--", color = "magenta", zorder = 10) #77
    else:
        signal_x_bottom, signal_x_top = ax[0].get_xlim()
        ax[0].hlines(toi_info["h"], xmin = signal_x_bottom, xmax = signal_x_top, linestyle = "--", color = "magenta", zorder = 10) #77
        ax[1].hlines(toi_info["h"], xmin = signal_x_bottom, xmax = signal_x_top, linestyle = "--", color = "magenta", zorder = 10) #77
        ax[1].hlines(1, xmin = signal_x_bottom, xmax = signal_x_top, linestyle = "--", color = "magenta", zorder = 10) #77
    
    ax[0].set_ylim([toi_info["signal_y_bottom"], toi_info["signal_y_top"]])
    ax[1].set_ylim([toi_info["signal_y_bottom"], toi_info["signal_y_top"]])
    ax[0].vlines(toi_info["duration"] / 24 / 2, ymin = toi_info["signal_y_bottom"], ymax = toi_info["signal_y_top"], linestyle = "--", color = "magenta", zorder = 10)
    ax[0].vlines(-toi_info["duration"] / 24 / 2, ymin = toi_info["signal_y_bottom"], ymax = toi_info["signal_y_top"], linestyle = "--", color = "magenta", zorder = 10)
    ax[1].vlines(toi_info["duration"] / 24 / 2, ymin = toi_info["signal_y_bottom"], ymax = toi_info["signal_y_top"], linestyle = "--", color = "magenta", zorder = 10)
    ax[1].vlines(-toi_info["duration"] / 24 / 2, ymin = toi_info["signal_y_bottom"], ymax = toi_info["signal_y_top"], linestyle = "--", color = "magenta", zorder = 10)
    
    ax[0].tick_params(axis = "x", labelsize = font_size)
    ax[0].tick_params(axis = "y", labelsize = font_size)
    ax[1].tick_params(axis = "x", labelsize = font_size)
    ax[1].tick_params(axis = "y", labelsize = font_size)
    
    figure.text(0.28, 0.94, "ZTF Light Curve", ha = "center", fontsize = font_size, fontweight = "normal")
    figure.text(0.76, 0.94, "Averaging of ZTF Light Curve", ha = "center", fontsize = font_size, fontweight = "normal") # Or 'ATLAS Light Curve'
    figure.text(0.5, 0.01, "Time from mid-transit\nFolded on period " + str(toi_info["period_display"]) + " " + u"\u00B1" + " " + str(toi_info["period_uncertainty_display"]) + " Epoch " + str(toi_info["epoch_display"]) + " " + u"\u00B1" + " " + str(toi_info["epoch_uncertainty_display"]) + " in BJD - 2457000", ha = "center", fontsize = font_size, fontweight = "normal") # Tess Julian Days = Julian Days - 2457000
    figure.text(0.01, 0.5, "Relative Flux", va = "center", rotation = "vertical", fontsize = font_size, fontweight = "normal")
    
    figure.set_size_inches((19.2, 6.8), forward = True) # If True, overwrite existing window
    figure.subplots_adjust(bottom = 0.15, top = 0.9, left = 0.075, right = 0.975, wspace = 0.1, hspace = 0.35)
    if is_list["plotting"]:
        plot.show()
    
    figure.set_size_inches((19.2, 6.8), forward = False)
    figure.savefig("TICs/" + str(toi_info["tic_ID"]) + "/Binned.png", dpi = 100)
    figure.savefig("TICs/" + str(toi_info["tic_ID"]) + "/Binned.pdf", dpi = 100)
    figure.savefig("TICs/" + str(toi_info["tic_ID"]) + "/Binned.svg", dpi = 100)

def plot_period_revision(toi_info, mini_all_star_dataframe, mini_cone, mini_clean_dataframe, lists, is_list):
    print("Revising period...")
    figure, ax = plot.subplots(1, 3) # 1 row, 3 columns (3 cells total)
    ax_setup(figure, ax[0], toi_info["signal_index"], toi_info, mini_all_star_dataframe, mini_cone, mini_clean_dataframe, lists, is_list["ATLAS"])
    periodigram_ax(figure, ax[1], toi_info, 12)
    ax_setup(figure, ax[2], toi_info["signal_index"], toi_info, mini_all_star_dataframe, mini_cone, mini_clean_dataframe, lists, is_list["ATLAS"])
    
    # Remove the axis labels and colors that have happened due to ax_setup()
    for text in figure.texts:
        text.set_visible(False)
    ax[0].set(xlabel = None)
    ax[2].set(xlabel = None)
    change_ax_colors(ax[0], "#000000")
    change_ax_colors(ax[2], "#000000")
    
    # Bottom titles
    figure.text(0.2, 0.01, "\nTime from Mid-transit (days)", ha = "center", fontsize = font_size, fontweight = "normal")
    figure.text(0.5, 0.01, "\nPeriod (days)", ha = "center", fontsize = font_size, fontweight = "normal")
    figure.text(0.85, 0.01, "\nTime from Mid-transit (days)", ha = "center", fontsize = font_size, fontweight = "normal")
    
    # Top titles
    figure.text(0.2, 0.93, "Light Curve Folded on TOI Period", ha = "center", fontsize = font_size, fontweight = "normal")
    figure.text(0.5, 0.93, "BLS Periodogram", ha = "center", fontsize = font_size, fontweight = "normal")
    figure.text(0.85, 0.93, "Light Curve Folded on Revised Period", ha = "center", fontsize = font_size, fontweight = "normal")
    
    figure.subplots_adjust(bottom = 0.11, top = 0.9, left = 0.06, right = 0.975, wspace = 0.23, hspace = 0.6) # Spacing
    
    # y-axis
    figure.text(0.005, 0.5, "Relative Flux", va = "center", rotation = "vertical", fontsize = font_size, fontweight = "normal")
    figure.text(0.655, 0.5, "Relative Flux", va = "center", rotation = "vertical", fontsize = font_size, fontweight = "normal")
    figure.text(0.33, 0.5, "Power", va = "center", rotation = "vertical", fontsize = font_size, fontweight = "normal")
    
    if is_list["saving"]:
        save_plot(("TICs/" + str(toi_info["tic_ID"]) + "/Period_Revision"), figure, False, True)
    
    figure.set_size_inches((19.2, 6.8), forward = True) # If True, overwrite existing window
    if is_list["plotting"]:
        plot.show()


# Reference field image rendering
def add_reference_image_patch(ax, mini_all_star_dataframe, mini_clean_dataframe, is_list, image_index, mini_comparison_star_index, mini_comparison_star_patch_x, mini_comparison_star_patch_y, mini_patch_color): # Adding the aperture to the reference field image
    edge_alpha = 0.8
    comparison_star_patch_x = mini_comparison_star_patch_x
    comparison_star_patch_y = mini_comparison_star_patch_y
    comparison_star_patch = matplotlib.patches.Ellipse((comparison_star_patch_x - 0.5, comparison_star_patch_y - 0.5), float(list(mini_all_star_dataframe["sigma_x"])[image_index]) * 2 * 3, float(list(mini_all_star_dataframe["sigma_y"])[image_index]) * 2 * 3, angle = float(list(mini_all_star_dataframe["theta"])[image_index]), fill = False, color = mini_patch_color, alpha = edge_alpha) # Each Position is - 0.5 because the image pixels are defined by the corners and the centroid is difined by the center so it is half a pixel off in each direction. Angle in Degrees
    
    ax.add_patch(comparison_star_patch)
    
    if is_list["ATLAS"]: # Displaying the annulus
        inner_x = (float(list(mini_all_star_dataframe["sigma_x"])[image_index]) * 2 * 3) + (float(list(mini_all_star_dataframe["sigma_x"])[image_index]) * 3)
        inner_y = (float(list(mini_all_star_dataframe["sigma_y"])[image_index]) * 2 * 3) + (float(list(mini_all_star_dataframe["sigma_y"])[image_index]) * 3)
        annulus_star_patch1 = matplotlib.patches.Ellipse((comparison_star_patch_x - 0.5, comparison_star_patch_y - 0.5), inner_x, inner_y, float(list(mini_all_star_dataframe["theta"])[image_index]), fill = False, color = mini_patch_color, alpha = 0.5)
        ax.add_patch(annulus_star_patch1)
        
        outer_x = (float(list(mini_all_star_dataframe["sigma_x"])[image_index]) * 2 * 3) + (float(list(mini_all_star_dataframe["sigma_x"])[image_index]) * 3) + (float(list(mini_all_star_dataframe["sigma_x"])[image_index]) * 2 * 3 * 1.5)
        outer_y = (float(list(mini_all_star_dataframe["sigma_y"])[image_index]) * 2 * 3) + (float(list(mini_all_star_dataframe["sigma_y"])[image_index]) * 3) + (float(list(mini_all_star_dataframe["sigma_y"])[image_index]) * 2 * 3 * 1.5)
        annulus_star_patch2 = matplotlib.patches.Ellipse((comparison_star_patch_x - 0.5, comparison_star_patch_y - 0.5), outer_x, outer_y, float(list(mini_all_star_dataframe["theta"])[image_index]), fill = False, color = mini_patch_color, alpha = 0.5)
        ax.add_patch(annulus_star_patch2)
    
    comparison_star_patch_label = 0
    if is_list["showing_index"]: # Whether to label based on how many frames from the dataframe or the index (shown in the plots)
        comparison_star_patch_label = ax.annotate(("T" + str(mini_comparison_star_index)), xy = (comparison_star_patch_x, comparison_star_patch_y + 7.5), ha = "center", color = mini_patch_color, fontsize = font_size - 4, alpha = edge_alpha)
    else:
        comparison_star_patch_label = ax.annotate(str(mini_all_star_dataframe[str(list(mini_clean_dataframe.columns[2:])[mini_comparison_star_index])].isnull().sum()), xy = (comparison_star_patch_x, comparison_star_patch_y + 7.5), ha = "center", color = mini_patch_color, alpha = edge_alpha)
    
    comparison_star_patch_label.set_color(mini_patch_color)
    
    comparison_star_patch.set_edgecolor(mini_patch_color)
    return comparison_star_patch

def plot_reference_image(image_index, toi_info, mini_all_star_dataframe, mini_clean_dataframe, mini_lists, mini_size, is_list): # Displaying the reference image    
    # Getting the image
    fits_file = fits.open(mini_all_star_dataframe["local_image_link"][image_index]) 
    fits_header = fits_file[0].header
    fits_image = fits_file[0].data
    
    # Fixing the WCS header
    if is_list["ATLAS"]:
        toss = ["RP_SCHIN", "RADECSYS"] 
        toss += [x for x in fits_header if x.startswith("PV")]
        toss += [x for x in fits_header if x.startswith("CNPIX")]
        for kk in toss:
            fits_header.remove(kk, remove_all = True, ignore_missing = True)
        
        fits_header["CTYPE1"] = "RA---TAN-SIP"
        fits_header["CTYPE2"] = "DEC--TAN-SIP"
        
        fits_header.update()
        
        fits_image = fits_image - np.median(fits_image.min(axis = 0))
        fits_image[fits_image <= 0] = 1
    
    wcs = WCS(fits_header)
    figure = plot.figure()
    ax = plot.subplot(projection = wcs)
    
    
    # Arrows
    pixel_x1 = 250
    pixel_y1 = 250
    if is_list["ATLAS"]: # Put the star in the corner more because the field is larger than ZTF
        pixel_x1 = 350
        pixel_y1 = 350
    sky = wcs.pixel_to_world(pixel_x1, pixel_y1)
    north_ra = sky.ra.degree
    north_dec = sky.dec.degree + (30 / 3600) # 1 degree = 3600 arcseconds
    if is_list["ATLAS"]: # Make the arrow look like the same size in ZTF (field in ATLAS is larger, so the arrow as the same size looks smaller)
        north_dec = sky.dec.degree + (50 / 3600) # 1 degree = 3600 arcseconds
    pixel_x2, pixel_y2 = wcs.wcs_world2pix(north_ra, north_dec, 1)
    north_arrow_width = pixel_x2 - pixel_x1
    north_arrow_length = pixel_y2 - pixel_y1
    ax.arrow(pixel_x1, pixel_y1, north_arrow_width, north_arrow_length, width = 5, head_width = 15, head_length = 10, color = "b")
    if pixel_y2 - pixel_y1 > 0:
        ax.annotate("N", xy = (pixel_x2, pixel_y2 + 20), ha = "center", color = "b", fontsize = font_size)
    else:
        ax.annotate("N", xy = (pixel_x2, pixel_y2 - 20), ha = "center", color = "b", fontsize = font_size)
    
    theta = np.arctan(north_arrow_length / -north_arrow_width) # dy / dx
    if pixel_y2 - pixel_y1 > 0: # East arrow pointing right
        ax.arrow(pixel_x1 + (2.75), pixel_y1 - (2.5 * np.sin(theta)) - 0.15, north_arrow_length, -north_arrow_width, width = 5, head_width = 15, head_length = 10, color = "r")
        ax.annotate("E", xy = (pixel_x1 + 2.75 + north_arrow_length + 10, pixel_y1 - (2.5 * np.sin(theta)) - 0.15) + north_arrow_width + 5, ha = "center", color = "r", fontsize = font_size)
    else: # East arrow pointing left
        ax.arrow(pixel_x1 - (2.75), pixel_y1 + (2.5 * np.sin(theta)) - 0.15, north_arrow_length, -north_arrow_width, width = 5, head_width = 15, head_length = 10, color = "r")
        ax.annotate("E", xy = (pixel_x1 - (2.75) + north_arrow_length - 10, pixel_y1 + (2.5 * np.sin(theta)) - 0.15) - north_arrow_width - 5, ha = "center", color = "r", fontsize = font_size)
    
    
    img = ax.imshow(np.log(fits_image - np.min(fits_image) + 1), cmap = "Greys") # Log color scale to more easily see all the potential stars
    
    # Adding the patches
    comparison_star_patches = {"target_patch": 0, "comparison_patch": 0, "saturated_patch": 0, "dim_patch": 0, "missing_patch": 0}
    for comparison_star_index in range(int(len(mini_clean_dataframe.columns) - 2)):
        if comparison_star_index == toi_info["target_index"]:
            comparison_star_patch = add_reference_image_patch(ax, mini_all_star_dataframe, mini_clean_dataframe, is_list, image_index, comparison_star_index, (float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_px")])[image_index]) + float(list(mini_all_star_dataframe["x0"])[image_index]) - (float(list(mini_all_star_dataframe["subarray_length"])[image_index]) / 2)), (float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_py")])[image_index]) + float(list(mini_all_star_dataframe["y0"])[image_index]) - (float(list(mini_all_star_dataframe["subarray_length"])[image_index]) / 2)), target_star_color)
            if is_list["ATLAS"]:
                comparison_star_patch = add_reference_image_patch(ax, mini_all_star_dataframe, mini_clean_dataframe, is_list, image_index, comparison_star_index, float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_px")])[image_index]), float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_py")])[image_index]), target_star_color)
            comparison_star_patches["target_patch"] = comparison_star_patch
        elif list(mini_clean_dataframe.columns[2:])[comparison_star_index] in mini_lists["comparison_stars"]: # Comparisons have no highlight
            comparison_star_patch = add_reference_image_patch(ax, mini_all_star_dataframe, mini_clean_dataframe, is_list, image_index, comparison_star_index, (float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_px")])[image_index]) + float(list(mini_all_star_dataframe["x0"])[image_index]) - (float(list(mini_all_star_dataframe["subarray_length"])[image_index]) / 2)), (float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_py")])[image_index]) + float(list(mini_all_star_dataframe["y0"])[image_index]) - (float(list(mini_all_star_dataframe["subarray_length"])[image_index]) / 2)), comparison_star_color)
            if is_list["ATLAS"]:
                comparison_star_patch = add_reference_image_patch(ax,mini_all_star_dataframe,  mini_clean_dataframe, is_list, image_index, comparison_star_index, float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_px")])[image_index]), float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_py")])[image_index]), comparison_star_color)
            comparison_star_patches["comparison_patch"] = comparison_star_patch
        elif list(mini_clean_dataframe.columns[2:])[comparison_star_index] in mini_lists["saturated_stars"]:
            comparison_star_patch = add_reference_image_patch(ax, mini_all_star_dataframe, mini_clean_dataframe, is_list, image_index, comparison_star_index, (float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_px")])[image_index]) + float(list(mini_all_star_dataframe["x0"])[image_index]) - (float(list(mini_all_star_dataframe["subarray_length"])[image_index]) / 2)), (float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_py")])[image_index]) + float(list(mini_all_star_dataframe["y0"])[image_index]) - (float(list(mini_all_star_dataframe["subarray_length"])[image_index]) / 2)), saturated_star_color)
            if is_list["ATLAS"]:
                comparison_star_patch = add_reference_image_patch(ax, mini_all_star_dataframe, mini_clean_dataframe, is_list, image_index, comparison_star_index, float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_px")])[image_index]), float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_py")])[image_index]), saturated_star_color)
            comparison_star_patches["saturated_patch"] = comparison_star_patch
        else:
            if is_list["ATLAS"]:
                if list(mini_clean_dataframe.columns[2:])[comparison_star_index] == str(toi_info["signal_tic_ID"]) + "_brightness":
                    comparison_star_patch = add_reference_image_patch(ax, mini_all_star_dataframe, mini_clean_dataframe, is_list, image_index, comparison_star_index, float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_px")])[image_index]), float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_py")])[image_index]), dim_star_color)
                    comparison_star_patches["dim_patch"] = comparison_star_patch
            else:
                comparison_star_patch = add_reference_image_patch(ax, mini_all_star_dataframe, mini_clean_dataframe, is_list, image_index, comparison_star_index, (float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_px")])[image_index]) + float(list(mini_all_star_dataframe["x0"])[image_index]) - (float(list(mini_all_star_dataframe["subarray_length"])[image_index]) / 2)), (float(list(mini_all_star_dataframe[str(str((list(mini_clean_dataframe.columns[2:])[comparison_star_index]).split("_")[0]) + "_py")])[image_index]) + float(list(mini_all_star_dataframe["y0"])[image_index]) - (float(list(mini_all_star_dataframe["subarray_length"])[image_index]) / 2)), dim_star_color)
                comparison_star_patches["dim_patch"] = comparison_star_patch 
    
    for comparison_star_index in range(len(list(mini_all_star_dataframe.columns[14::6]))):
        name_of_thing = str(str(list(mini_all_star_dataframe.columns[14::6])[comparison_star_index]).split("_")[0])
        if np.median(mini_all_star_dataframe[name_of_thing + "_distance_to_real_target"]) < mini_size:
            if not list(mini_all_star_dataframe.columns[14::6])[comparison_star_index] in mini_clean_dataframe.columns:
                comparison_star_patch = add_reference_image_patch(ax, mini_all_star_dataframe, mini_clean_dataframe, is_list, image_index, False, (float(list(mini_all_star_dataframe[name_of_thing + "_px"])[image_index]) + float(list(mini_all_star_dataframe["x0"])[image_index]) - (float(list(mini_all_star_dataframe["subarray_length"])[image_index]) / 2)), (float(list(mini_all_star_dataframe[name_of_thing + "_py"])[image_index]) + float(list(mini_all_star_dataframe["y0"])[image_index]) - (float(list(mini_all_star_dataframe["subarray_length"])[image_index]) / 2)), too_few_points_color)
                comparison_star_patches["missing_patch"] = comparison_star_patch
    
    # Legend
    comparison_star_patches_dictionary = {"patch": [], "label": []}
    for patch in comparison_star_patches: # patch is a key, not a value
        if not isinstance(comparison_star_patches[patch], int): # If it is an int it was not something that appeared in the plot (the dictionary as empty was an int originally)
            comparison_star_patches_dictionary["patch"].append(comparison_star_patches[patch])
            if patch == "target_patch":
                comparison_star_patches_dictionary["label"].append("Target")
            elif patch == "comparison_patch":
                comparison_star_patches_dictionary["label"].append("Comparison")
            elif patch == "saturated_patch":
                comparison_star_patches_dictionary["label"].append("Too Saturated for Comparison")
            elif patch == "dim_patch":
                comparison_star_patches_dictionary["label"].append("Too Faint for Comparison")
            else:
                comparison_star_patches_dictionary["label"].append("Too Few Data Points")
    
    ax.legend(comparison_star_patches_dictionary["patch"], comparison_star_patches_dictionary["label"], bbox_to_anchor = (0.6, 0.15))
    ax.set_axis_off()
    fontprops = fm.FontProperties(size = font_size)
    scalebar = AnchoredSizeBar(ax.transData, 100, "100 arc-seconds", "lower center", color = "black", frameon = False, size_vertical = 1, fontproperties = fontprops, bbox_to_anchor = Bbox.from_bounds(0, 0.05, 0.7, -1), bbox_transform = ax.figure.transFigure) # 1 pixel = 1 arcsecond
    if is_list["ATLAS"]:
        scalebar = AnchoredSizeBar(ax.transData, 54, "100 arc-seconds", "lower center", color = "black", frameon = False, size_vertical = 1, fontproperties = fontprops, bbox_to_anchor = Bbox.from_bounds(0, 0.05, 0.7, -1), bbox_transform = ax.figure.transFigure) # 1 pixel = 1.86 arcseconds
    ax.add_artist(scalebar)
    
    make_plot_fullscreen()
    
    if is_list["saving"]:
        save_plot(("TICs/" + str(toi_info["tic_ID"]) + "/TIC_" + str(toi_info["tic_ID"]) + "_" + mini_all_star_dataframe["image_filter"][image_index] + "_Filter_Reference_Image"), figure, True, False)
    
    if is_list["plotting"]:
        plot.show()


# Outputs
def make_plot_fullscreen():
    manager = plot.get_current_fig_manager()
    #manager.window.showMaximized() # Display plot fullscreen to the user

def save_plot(mini_image_name, figure, mini_is_reference, mini_is_ExoFOP): # Saving the plot as specific image types
    # Save to an landscape 8.5 by 11 image for normal images
    width = 19.2
    height = 10.8
    
    if mini_is_ExoFOP: # For the 2 subplot one
        height = 6.8
    
    if mini_is_reference: # For the reference image
        width = 10.8
        
    figure.set_size_inches((width, height), forward = False)
    figure.savefig(mini_image_name + ".png", dpi = 100)
    figure.savefig(mini_image_name + ".pdf", dpi = 100)
    figure.savefig(mini_image_name + ".svg", dpi = 100)

def create_report(tic_ID): # PDF of all the saved images
    print("Creating report...")
    pdf = FPDF()
    filenames = []
    for filename in list(os.listdir("TICs/" + str(tic_ID))):
        filename = "TICs/" + str(tic_ID) + "/" + filename
        if filename.endswith(".png"):
            if "Reference_Image" in filename:
                filenames.insert(0, filename)
            else:
                filenames.append(filename)
    
    for filename in filenames:
        pdf.add_page()
        pdf.rotate(270)
        if "Reference_Image" in filename:
            pdf.image(filename, 10, -200, 225)
        else:
            pdf.image(filename, 10, -160, 275)
    
    pdf.output(("TICs/" + str(tic_ID) + "/" + str(tic_ID) + "_Lightcurve_Report.pdf"), "F")
    print("Done!")



def setup(tic_ID, is_ATLAS, comparison_removal = [], signal_tic_ID = 0, revised_period = 0, size = 100, xlim = "zoomed", frame_number = 1, is_plotting = True, is_showing_index = True, is_saving = True, is_plotting_MAD_vs_MAG = True, is_lcbin = False, is_period_revision = False, is_done = False): # Running everything together
    '''
    Arguments:
        tic_ID -- TIC ID of the target to run
        is_ATLAS -- If is_ATLAS, data from ATLAS will be plotted (otherwise ZTF data will). ATLAS data will automatically run if the declination is less than -28 degrees as that is the current bound for ZTF's telescope imaging. Default for is_ATLAS is set to False
        comparison_removal -- array of TIC IDs of noisy comparison stars that need to be removed. This array is automatically added to while the program is running based on the median-absolute deviation (MAD) of the star; greater MADs show noisier light curves. If a noisy star is used as a comparison but not automatically removed, add the TIC ID of the star to this list. In order to figure out if a comparison star needs to be removed, consult the MAD versus TMag plot as an output of this program. Default for comparison_removal is set to an empty array and is populated automatically
        signal_tic_ID -- tic ID of the signal star. If the signal_tic_ID is not 0 and is a valid TIC ID, the target versus signal plot will be created. Default for signal_tic_ID is set to 0
        revised_period -- found revised period in days. When revised_period is not 0, light curves will be folded on the value of revised_period. Default for revised_period is set to 0
        size -- radius size of field to chose comparison stars and plot light curves in arcseconds. The radius must be large enough to incorporate 2 valid comparison stars for light curve creation
        xlim -- x-limits of the light curve plots. When xlim = "zoomed", the x-limits of the light curves are auotmatically set to 8 times the length of the transit duration. When xlim = "all", the entirety of the datapoints folded along the period are displayed. xlim can also be manually set to any positive number for a custom display (based on the times folded on the period). Default for xlim is set to "zoomed"
        frame_number -- image index to display as the reference photo. Chosing indexes depends on the order of the All_[TIC ID]_Star_Data.csv. It is really only necessary to change the frame_number if the auto-chosen frame looks bad. Default for frame_number is set to 1
        is_plotting -- if plots are displayed or just saved. Default for is_plotting is set to True
        is_showing_index -- if the reference field image is displaying the corresponding index of the star (respective to its position within the returned light curve sheets) when set to True, or if each star displays the number of frames that that star has successfully processed and incorporated into light curve creation. This argument when set to False is mearly a diagnostic for checking as some ZTF images along its limits are smaller than requested image size resulting in stars not being processed. Default for is_showing_index is set to True
        is_saving -- is the outputted plots are saved. When True, images of the reference field and each sheet will be saved as .pdf, .png, and .eps form for any later purpose. Also, a .pdf report of all the .pngs combined will be rendered into 1 singular .pdf for more easy sharing. Default for is_saving is set to True
        is_plotting_MAD_vs_MAG -- if is plotting noise diagnostic plot (max median absolute deviation versus TESS magnitude). This diagnostic helps identify noisy comparison stars if they need to be removed and gives a visual how clean the light curves of all the stars in the field are. Default for is_plotting_MAD_vs_MAG is set to True
        is_lcbin -- if is plotting the binned light curves for easier viewing of shallower or messier transit light curves. The averaged light curves will be displayed on the tic_ID (or signal_tic_ID if provided). Default for is_lcbin is set to False
        is_period_revision -- if is preforming a Box-Least Squares period revison function. If set to True, a period will be tested on the light curve tic_ID (or the signal_tic_ID if provided) and display the transit folded on the original versus edited periods. Default for is_period_revision is set to False
        is_done -- interior argument used exclusively for running this program. DO NOT CHANGE. If a compsrison star has a large max median absolute deviation, it will be automatically removed and is_done will be updated to rerun this program without that comparison star. This process will repeat until all comparison stars have a max median absolute deviation less than 0.2. Default for is_done is set to False
        '''
    
    
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning) # Removing HTTPS warning
    
    
    clear_output(wait = True) # Clearing the output so the jupyter doesn't crash
    plot.close("all") # Closing all the plots so we don't run into any memory issues
    print(tic_ID)
    
    if signal_tic_ID == 0:
        signal_tic_ID = tic_ID
    
    
    # Opening the .csvs
    toi_info = get_TOI_info(str(tic_ID))
    toi_info["tic_ID"] = tic_ID
    toi_info["signal_tic_ID"] = signal_tic_ID
    toi_info["old_period"] = toi_info["period"]
    toi_info["h"] = 1
    if (not revised_period == toi_info["period"]) and (not revised_period == 0):
        toi_info["period"] = revised_period
    
    toi_info["xlim"] = xlim # zoomed (4x the duration), all (auto-generated), or a float (manually setting the xlim)
    all_star_dataframe = pd.read_csv("TICs/" + str(toi_info["tic_ID"]) + "/All_" + str(toi_info["tic_ID"]) + "_Star_Data.csv")
    if not "time" in list(all_star_dataframe.columns):
        print("Run program not finished or failed. Renaming columns...")
        all_star_dataframe.rename(columns = {"jd_time": "time"}, inplace = True) # In case program crashed
    all_star_dataframe = all_star_dataframe.dropna(subset = list(all_star_dataframe.columns[14:]), how = "all")
    cone = pd.read_csv("TICs/" + str(toi_info["tic_ID"]) + "/" + str(toi_info["tic_ID"]) + "_Cone.csv")
    
    
    # Converting and getting the display values
    toi_info["period_uncertainty"] = display_period_uncertainty(all_star_dataframe, toi_info)
    toi_info["depth"] = toi_info["depth"] / 1000000 # Converting transit depth from ppm to a fraction
    toi_info["transit_start"], toi_info["transit_end"] = get_transit_BJD_times_offset(toi_info)
    toi_info["epoch_display"], toi_info["epoch_uncertainty_display"] = display_rounded_uncertainties(toi_info["epoch"] - 2457000, toi_info["epoch_uncertainty"])
    toi_info["period_display"], toi_info["period_uncertainty_display"] = display_rounded_uncertainties(toi_info["period"], toi_info["period_uncertainty"])
    
    
    # Getting the usable data and collecting the correct comparison stars
    delta_m = (np.log(1 / toi_info["depth"]) / np.log(2.512)) + 0.5 # Positive number
    
    old_cone = cone
    cone = cone[(cone["Tmag"] < (cone.at[int(np.where(cone["ID"] == toi_info["tic_ID"])[0]), "Tmag"] + delta_m)) | (cone["ID"] == toi_info["tic_ID"])]
    cone = cone.reset_index(drop = True)
    
    lists = {"comparison_stars": [], "saturated_stars": [], "MAD": {"g": [], "r": [], "i": [], "o": [], "c": []}, "MAD_comparison": {"g": [], "r": [], "i": [], "o": [], "c": []}, "MAG": [], "MAG_comparison": [], "MAD_vs_MAG_comparison_names": []}
    is_list = {"ATLAS": is_ATLAS, "plotting": is_plotting, "showing_index": is_showing_index, "saving": is_saving, "plotting_MAD_vs_MAG": is_plotting_MAD_vs_MAG, "lcbin": is_lcbin, "period_revision": is_period_revision}
    
    # Checking if I need to force use ATLAS
    if float(toi_info["dec"]) < -28: # ZTF limits are -28 degrees dec
        print("dec is past -28 (" + str(dec) + "), which is too far South for ZTF processing. Automatically downloading ATLAS images...")
        is_list["ATLAS"] = True
    
    
    clean_dataframe = all_star_dataframe[["time", "image_filter"]].copy()
    if is_list["ATLAS"]:
        minimum_magnitude_for_comparison = 12 # ZTF lower bound good comparison TMag
        maximum_magnitude_for_comparison = 15.5 # ZTF upper bound good comparison TMag
    else:
        minimum_magnitude_for_comparison = 14 # ATLAS lower bound good comparison TMag
        maximum_magnitude_for_comparison = 17 # ATLAS upper bound good comparison TMag
    for column in list(cone["ID"]):
        if all_star_dataframe[str(column) + "_brightness"].isnull().sum() < (len(all_star_dataframe[str(toi_info["tic_ID"]) + "_brightness"]) / 4):
            if cone.at[int(np.where(cone["ID"] == column)[0]), "Tmag"] > minimum_magnitude_for_comparison:
                if np.nanpercentile(list(all_star_dataframe[str(column) + "_max_pixel"]), 95) < 65000: # Making sure it is not saturated
                    if cone.at[int(np.where(cone["ID"] == column)[0]), "Tmag"] < maximum_magnitude_for_comparison: # Good comparison star
                        if np.median(all_star_dataframe[str(column) + "_distance_to_real_target"]) < size:
                            if not str(cone.at[int(np.where(cone["ID"] == column)[0]), "ID"]) == str(toi_info["tic_ID"]):
                                clean_dataframe[str(column) + "_brightness"] = all_star_dataframe[str(column) + "_brightness"]
                                lists["comparison_stars"].append(str(column) + "_brightness")
                    else: # Really dim star
                        if np.median(all_star_dataframe[str(column) + "_distance_to_real_target"]) < size:
                            if not str(cone.at[int(np.where(cone["ID"] == column)[0]), "ID"]) == str(toi_info["tic_ID"]):
                                clean_dataframe[str(column) + "_brightness"] = all_star_dataframe[str(column) + "_brightness"]
                else: # Saturated star
                    if np.median(all_star_dataframe[str(column) + "_distance_to_real_target"]) < size:
                        if not str(cone.at[int(np.where(cone["ID"] == column)[0]), "ID"]) == str(toi_info["tic_ID"]):
                            clean_dataframe[str(column) + "_brightness"] = all_star_dataframe[str(column) + "_brightness"]
                            lists["saturated_stars"].append(str(column) + "_brightness")
                            print("line 167    " + str(cone.at[int(np.where(cone["ID"] == column)[0]), "Tmag"]) + "   " + str(max(list(all_star_dataframe[str(column) + "_max_pixel"]))) + "     " + str(column))
            else: # Saturated star
                if np.median(all_star_dataframe[str(column) + "_distance_to_real_target"]) < size:
                    if not str(cone.at[int(np.where(cone["ID"] == column)[0]), "ID"]) == str(toi_info["tic_ID"]):
                        clean_dataframe[str(column) + "_brightness"] = all_star_dataframe[str(column) + "_brightness"]
                        lists["comparison_stars"].append(str(column) + "_brightness")
    
    
    # Remove stars with really large MADs
    if (str(toi_info["signal_tic_ID"]) + "_brightness") in lists["comparison_stars"]:
        lists["comparison_stars"].remove(str(toi_info["signal_tic_ID"]) + "_brightness")
    
    for bad_comparison_tic_ID in comparison_removal:
        if not str(bad_comparison_tic_ID) == str(signal_tic_ID):
            lists["comparison_stars"].remove(str(bad_comparison_tic_ID) + "_brightness")
    
    clean_dataframe[str(toi_info["signal_tic_ID"]) + "_brightness"] = all_star_dataframe[str(toi_info["signal_tic_ID"]) + "_brightness"] # Adding the signal
    clean_dataframe[str(toi_info["tic_ID"]) + "_brightness"] = all_star_dataframe[str(toi_info["tic_ID"]) + "_brightness"] # Adding the target
    
    # Clearing the nan rows so all of the rows have data
    pd.set_option("display.max_columns", None)
    clean_dataframe = clean_dataframe.dropna()
    if is_list["ATLAS"]: # Needs an extra cleaning up because some values are 0
        clean_dataframe.drop(clean_dataframe[(clean_dataframe.values == 0).any(axis = 1)].index.values.tolist(), axis = 0, inplace = True)
    
    print(str(len(clean_dataframe.columns[2:])) + " stars are usable out of " + str(len(old_cone)))
    print(str(len(old_cone) - len(clean_dataframe.columns[2:])) + " stars that could be bright enough did not have enough points ")
    print(str(len(lists["comparison_stars"])) + " comparison stars")
    if len(lists["comparison_stars"]) < 2:
        print("There are not enough comparison stars to plot.")
        sys.exit()
    
    
    # Getting the target
    toi_info["target_index"] = list(clean_dataframe.columns).index(str(toi_info["tic_ID"]) + "_brightness") - 2 # Need to exclude this column from the comparisons
    toi_info["signal_index"] = list(clean_dataframe.columns).index(str(toi_info["signal_tic_ID"]) + "_brightness") - 2 # Need to exclude this column from the comparisons
    print("Target star (TIC " + str(cone.at[int(np.where(cone["ID"] == toi_info["tic_ID"])[0]), "ID"]) + ") is at index " + str(toi_info["target_index"]))
    
    
    # Creating plots
    plot_reference_image(frame_number, toi_info, all_star_dataframe, clean_dataframe, lists, size, is_list)
    plot_lightcurves(toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list, size, frame_number, comparison_removal, is_done)
    
    if is_period_revision:
        plot_period_revision(toi_info, all_star_dataframe, cone, clean_dataframe, lists, is_list)
    
    if is_lcbin:
        lcbin(toi_info, is_list)
    
    if is_done:
        create_report(tic_ID)
