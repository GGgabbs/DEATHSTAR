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
# Filename:     Extracting_Lightcurves.py
# Created by:   Gabrielle Ross
# Last updated: 2/21/2024
# Github:       https://github.com/GGgabbs/DEATHSTAR/tree/main
# 
# Please cite our paper if you use this code!
# ADS link: https://ui.adsabs.harvard.edu/abs/2023MNRAS.tmp.3722R/abstract
#
# DEATHSTAR:    A system for confirming planets and identifying
#               false positive signals in TESS data using 
#               ground-based time domain surveys
#
# Purpose:      Collects all stars within the field,
#               downloads either ZTF or ATLAS ground data,
#               extracts brightnesses of each star in the field,
#               and saves brightnesses in All_[TIC ID]_Star_Data.csv



# Modules
# Arrays and handling data
import numpy as np
import pandas as pd
import math

# Astronomical files
from astropy.io import fits
from astropy.table import Table
from astroquery.mast import Mast
from astroquery.vizier import Vizier

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
import urllib3
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
from astropy import units as u
from mpfit import mpfit # Making a gaussian fit for the aperture

from IPython.display import clear_output # Clearing output to not lag notebook

import pdb # For viewing and stopping the program

# Fixing time offsets
from astropy.time import Time
import datetime
import time
import calendar

from PyAstronomy import pyasl



# Private functions
def vizier_query(mini_ra, mini_dec, mini_radius): # Query using Vizier to get the cone and all of the stars inside the cone
    all_catalogs = Vizier(columns = ["all"], row_limit = 200).query_region(astropy.coordinates.SkyCoord(ra = mini_ra, dec = mini_dec, unit = (u.deg, u.deg)), radius = Angle(str(mini_radius) + "d"))
    
    if "IV/39/tic82" in list(all_catalogs.keys()): # TIC
        print("Using Vizier's TIC catalog...")
        return all_catalogs["IV/39/tic82"]
    elif "IV/39/gaiadr3" in list(all_catalogs.keys()): # GAIA DR3
        print("Using Vizier's GAIA DR3 catalog...")
        return all_catalogs["I/355/gaiadr3"]
    else:
        print("Cone retrieval from Vizier TIC and GAIA DR3 catalog has also failed.")
        sys.exit()
        

# Query using MAST to get the cone and all the stars inside the cone
def mast_query(request):
    request_url = "https://mast.stsci.edu/api/v0/invoke" # Base URL
    
    version = ".".join(map(str, sys.version_info[:3])) # Grab Python version
    
    headers = { # Create HTTP header variables
        "Content-type": "application/x-www-form-urlencoded",
        "Accept": "text/plain",
        "User-agent": "python-requests/" + str(version)
    }
    
    # Encoding the request as a JSON string
    req_string = json.dumps(request)
    req_string = urlencode(req_string)
    
    response = requests.post(str(request_url), data = "request=" + str(req_string), headers = headers) # Perform the HTTP request
    
    # Pull out the headers and response content
    head = response.headers # Content where head is the response HTTP headers
    content = response.content.decode("utf-8") # Returned data
    return head, content

def tic_advanced_search_position_rows(mini_ra, mini_dec, mini_radius):
    request = {"service": "Mast.Catalogs.Filtered.Tic.Position.Rows",
               "format": "json",
               "params": {
                   "columns": "*",
                   "filters": [],
                   "ra": mini_ra,
                   "dec": mini_dec,
                   "radius": mini_radius
               }
    }
    
    headers, out_string = mast_query(request)
    out_data = json.loads(out_string)
    #out_data["status"] = "fail" # For testing Vizier
    
    print("MAST retrieval: " + str(out_data["status"]))
    if str(out_data["status"]) == "COMPLETE":
        dataframe = pd.DataFrame(out_data["data"])
        
        # Replace with zero will make sure there is no time offset but it won't crash because not all the pmRA and pmDEC values are in the cone. Not needed from Vizier because every value is filled in
        dataframe["pmRA"] = dataframe["pmRA"].fillna(0)
        dataframe["pmDEC"] = dataframe["pmDEC"].fillna(0)
    else: # Needs to query Vizier
        print("Cone retrieval from MAST has failed. Querying star data from Vizier...")
        dataframe = vizier_query(mini_ra, mini_dec, mini_radius)
        
        dataframe = pd.DataFrame(np.array(dataframe))
        dataframe = pd.DataFrame({"ID": list(dataframe["TIC"]), "ra": list(dataframe["RAJ2000"]), "dec": list(dataframe["DEJ2000"]), "Tmag": list(dataframe["Tmag"]), "pmRA": list(dataframe["pmRA"]), "pmDEC": list(dataframe["pmDE"])})
    
    return dataframe


def get_TOI_info(mini_tic_ID): # Getting the info from tev
    url = "https://exofop.ipac.caltech.edu/tess/target.php?id=" + str(mini_tic_ID) + "&json"
    response = urllib.request.urlopen(url)
    data = json.loads(response.read())
    return float(data["coordinates"]["ra"]), float(data["coordinates"]["dec"]), float(data["planet_parameters"][1]["dep_p"]) # ra, dec, transit_depth (ppm / parts per million)


def get_distance_from_target(ra_degrees1, dec_degrees1, ra_degrees2, dec_degrees2): # Getting the distance between 2 points of a great circle
    return np.array(pyasl.getAngDist(np.array(ra_degrees1), np.array(dec_degrees1), np.array(ra_degrees2), np.array(dec_degrees2))) * 3600


# Time conversions
def calculate_time_offset(julian_days, mini_ra, mini_dec, pm_ra, pm_dec): # Changing for the 2000 time offset
    times = [julian_days]
    t = Time(times, format = "jd", scale = "utc")
    date = datetime.datetime.strptime(t.iso[0][:10], "%Y-%m-%d")
    if calendar.isleap(int(t.iso[0][:4])):
        days = 366
    else:
        days = 365
    
    seconds_in_months = []
    for month in range(0, int(t.iso[0][5:7]) - 1): # Iterating through one less month than the date's month because those days can just be tacked on
        seconds_in_months.append((calendar.monthrange(int(t.iso[0][:4]), month + 1)[1]) * 24 * 60 * 60) # Seconds in month
    seconds_in_months.append(int(t.iso[0][8:10]) * 24 * 60 * 60)
    
    time_string = time.strptime(t.iso[0][11:],"%H:%M:%S.%f")
    seconds_in_time = datetime.timedelta(hours = time_string.tm_hour, minutes = time_string.tm_min, seconds = time_string.tm_sec).total_seconds() # Seconds of my timestamp from the TIME ONLY
    
    datetime_decimal = int(t.iso[0][:4]) + ((sum(seconds_in_months) + seconds_in_time) / (days * 24 * 60 * 60))
    
    pm_ra = pm_ra / 3600000 # 1 degree = 3600000 milliarcseconds. pmRA and pmDEC need to be converted from (milliarc seconds)/year to degrees/year
    pm_dec = pm_dec / 3600000
    delta_ra = ((datetime_decimal - 2000.0) * pm_ra) / (np.cos(np.deg2rad(mini_dec)))
    corrected_ra = mini_ra + delta_ra
    delta_dec = (datetime_decimal - 2000.0) * pm_dec
    corrected_dec = mini_dec + delta_dec
    
    return corrected_ra, corrected_dec

def convert_BJD_and_JD(days, mini_ra, mini_dec, BJD_to_JD):
    URL = "https://astroutils.astronomy.osu.edu/time/convert.php?JDS=" + str(days) + "&RA=" + str(mini_ra) + "&DEC=" + str(mini_dec) + "&FUNCTION=utc2bjd"
    if BJD_to_JD: # Converting from BJD to JD-UTC
        URL = URL.replace("utc2bjd", "bjd2utc")
    
    response = requests.get(URL, verify = False)
    content = response.content.decode("utf-8") # Returned data
    
    return list(map(float, content.split("\n")[:-1]))


# mpyfit aperture brightness extraction
def create_2D_gaussian(x_data_tuple, parameters): # Creating the 2D gaussian
    (x, y) = x_data_tuple 
    
    amplitude, x0, sigma_x, y0, sigma_y, theta, offset = parameters
    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    
    g = offset + amplitude*np.exp( - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0)+ c*((y-y0)**2)))
    return g.ravel()

def myfunct(parameters, fjac = None, x = None, y = None, err = None, brightness = None): # If fjac == None, then partial derivatives should not be computed. It will always be None if MPFIT is called with default flag
    model = create_2D_gaussian((x, y), parameters)
    status = 0 # Non-negative status value means MPFIT should continue, negative means stop the calculation
    return [status, ((brightness - model.reshape(int(math.sqrt(model.shape[0])), int(math.sqrt(model.shape[0]))))).reshape(model.shape[0])]

def find_ellipse_from_2D_gaussian_fit(image_max, image_median, subarray):
    max_y, max_x = np.unravel_index(np.argmax(subarray), subarray.shape)
    parameters = [image_max - image_median, max_x, 2, max_y, 2, 0, image_median] # A (amplitude -- like the brightest pixel - the noise), x0 (actual x center of star), y0 (actual y center of star), sigma_x (half the width of the star in the x direction), sigma_y (half the width of the star in the y direction), c (vertical shift -- noise)
    
    x = np.arange(subarray.shape[0]) # Columns of subarray -- NOT BRIGHTNESSES
    y = np.arange(subarray.shape[1]) # Rows of subarray -- NOT BRIGHTNESSES
    x, y = np.meshgrid(x, y)
    
    data = create_2D_gaussian((x, y), parameters)
    
    start_parameters = parameters # Choose starting values for parameter you are trying to fit. Roughly like what they should be
    parameters_info = [ # Step sizes for parameters. If fixed == True, then it will stay constant and not fit for that parameter, so we want all of them to be False
        {"fixed": False, "step": 1e1, "mpprint": 0} # Amplitude
        , {"fixed": False, "step": 1e-2, "mpprint": 0} # x0
        , {"fixed": False, "step": 1e-3, "mpprint": 0} # sigma_x
        , {"fixed": False, "step": 1e-2, "mpprint": 0} # y0
        , {"fixed": False, "step": 1e-3, "mpprint": 0} # sigma_y
        , {"fixed": False, "step": 1e-2, "mpprint": 0} # Theta
        , {"fixed": False, "step": 1e0, "mpprint": 0} # Offset
    ]
    
    fa = {"x": x, "y": y, "brightness": subarray}
    model = mpfit(myfunct, start_parameters, functkw = fa, parinfo = parameters_info, maxiter = 100, nprint = 0) #, quiet = "quiet")
    
    
    return model.params[1], model.params[2], model.params[3], model.params[4], model.params[5] # x0 (center), sigma_x (semimajor), y0 (center), sigma_y (semiminor), theta

def get_brightness(fits_image, subarray_length, px, py, x0, sigma_x, y0, sigma_y, theta, mini_is_ATLAS):
    number_of_sigma = 3
    major_axis = max(number_of_sigma * sigma_x, number_of_sigma * sigma_y) # a is always bigger
    minor_axis = min(number_of_sigma * sigma_x, number_of_sigma * sigma_y) # b is always smaller
    
    positions_dictionary = {"px": np.array(px) + x0 - (subarray_length / 2), "py": np.array(py) + y0 - (subarray_length / 2)}
    if mini_is_ATLAS:
        positions_dictionary = {"px": np.array(px), "py": np.array(py)}
    
    positions_dataframe = pd.DataFrame(positions_dictionary)
    positions_dataframe["point"] = list(zip(positions_dataframe["px"], positions_dataframe["py"]))
    #print(list(positions_dataframe["px"]))
    #hi
    #print("1")
    aperture = EllipticalAperture(list(positions_dataframe["point"]), major_axis, minor_axis, theta = Angle(str(theta) + "d")) # Either takes Angle() or radians but Angle() was more accurate (takes in degrees rather than conversion)
    #print("1.1")
    aperture_stats = ApertureStats(fits_image, aperture)
    #print("1.2")
    aperture_area = aperture.area_overlap(fits_image)
    #print("1.3")
    total_background = np.median(fits_image) * np.array(aperture_area)
    #print("2")
    if mini_is_ATLAS:
        sigma_clip = SigmaClip(sigma = 3, maxiters = 10)
        annulus = EllipticalAnnulus(list(positions_dataframe["point"]), (major_axis + 0.5 * major_axis), (major_axis + 0.5 * major_axis + major_axis * 1.5), (minor_axis + 0.5 * minor_axis), theta = Angle(str(theta) + "d"))
        annulus_stats = ApertureStats(fits_image, annulus, sigma_clip = sigma_clip)
        total_background = annulus_stats.median * aperture_stats.sum_aper_area.value
    #print("3")
    phot_table = aperture_photometry(fits_image, aperture)
    phot_background_subtraction = phot_table["aperture_sum"] - total_background
    #print("4")
    
    return np.array(phot_background_subtraction), np.array(aperture_stats.max)

def get_offset(fits_image, px, py, mini_sigma_coefficient):
    _, median, _ = sigma_clipped_stats(fits_image, sigma = 5)
    fits_image_without_median = fits_image - median # Subtract background from the data
    
    subarray = fits_image_without_median
    for subtraction_index in range(10, 1, -1):
        if ((int(px) - subtraction_index) > 0) and ((int(py) - subtraction_index) > 0) and ((int(px) + subtraction_index) < len(fits_image[0])) and ((int(py) + subtraction_index) < len(fits_image)):
            subarray = fits_image_without_median[int(py) - subtraction_index : int(py) + subtraction_index, int(px) - subtraction_index : int(px) + subtraction_index]
            break
        #else:
        #    print((int(px) - subtraction_index), (int(px) + subtraction_index), (int(py) - subtraction_index), (int(py) + subtraction_index))
    
    print("Subarray size (pixels): (" + str(subarray.shape[0]) + " by " + str(int(subarray.shape[1])) + ")")
    x0, sigma_x, y0, sigma_y, theta = find_ellipse_from_2D_gaussian_fit(np.max(fits_image), median, subarray) # Degrees
    sigma_x = sigma_x * mini_sigma_coefficient # 2 will get more of the fallout of brightness, whereas 0.9 will be a tighter aperture
    sigma_y = sigma_y * mini_sigma_coefficient
    
    
    return fits_image_without_median, int(subarray.shape[0]), x0, sigma_x, y0, sigma_y, theta


def add_new_frame(jd_time, local_image_link, image_median, image_filter, exposure_time, airmass, x0, sigma_x, y0, sigma_y, theta, subarray_length, mini_all_star_data): # Updating .csv file
    mini_all_star_data.append({"jd_time": jd_time, "local_image_link": local_image_link, "image_median": image_median, "image_filter": image_filter, "exposure": exposure_time, "airmass": airmass, "x0": x0, "sigma_x": sigma_x, "y0": y0, "sigma_y": sigma_y, "theta": theta, "processing_time": 0, "subarray_length": subarray_length})


def open_fits_file(mini_tic_ID, mini_frame_index, fits_URL, mini_cone, mini_target_star_index, mini_ra, mini_dec, mini_is_processing, mini_is_ATLAS, mini_sigma_coefficient, all_star_data): # Extracting light curves from the .fits file
    start_frame_time = time.time() # Starting the timer
    
    fits_file = fits.open(fits_URL) # fits can open URLS directly without needing to save the data as a .fits file
    fits_image = np.array(fits_file[0].data)
    fits_header = fits_file[0].header
    
    # Fixing the WCS header
    if mini_is_ATLAS:
        print("Editing WCS headers...")
        toss = ["RP_SCHIN", "RADECSYS"] 
        toss += [x for x in fits_header if x.startswith("PV")]
        toss += [x for x in fits_header if x.startswith("CNPIX")]
        for kk in toss:
            fits_header.remove(kk, remove_all = True, ignore_missing = True)
        
        fits_header["CTYPE1"] = "RA---TAN-SIP"
        fits_header["CTYPE2"] = "DEC--TAN-SIP"
        
        fits_header.update()
        
        # Changing the contrast
        fits_image = fits_image - np.median(fits_image.min(axis = 0))
        fits_image[fits_image <= 0] = 1 # Making all the things less than or equal to zero the new second minimum so they don't crash when plotting
    
    target_w = WCS(fits_header)
    print("Picture size (pixels): " + str(len(fits_image[0])) + " by " + str(len(fits_image)))
    
    if mini_is_processing:
        # Fixing the WCS header
        time_fits_header = 0
        exposure_fits_header = 0
        if mini_is_ATLAS:
            time_fits_header = float(fits_header["MJD-OBS"]) + 2400000.5
            exposure_fits_header = fits_header["EXPTIME"]
        else:
            time_fits_header = float(fits_header["OBSJD"])
            exposure_fits_header = fits_header["exposure"]
        
        # Getting the pixel locations from ra and dec
        all_star_ra, all_star_dec = calculate_time_offset(float(time_fits_header), np.array(mini_cone["ra"]), np.array(mini_cone["dec"]), np.array(mini_cone["pmRA"]), np.array(mini_cone["pmDEC"]))
        print("WCS (ra and dec in degrees): (" + str(all_star_ra[int(np.where(mini_cone["ID"] == mini_tic_ID)[0][0])]) + ", " + str(all_star_dec[int(np.where(mini_cone["ID"] == mini_tic_ID)[0][0])]) + ") ")
        all_star_px, all_star_py = target_w.wcs_world2pix(all_star_ra, all_star_dec, 1)
        
        if mini_is_ATLAS: # Need to correct for offset due to fits WCS irregularity. The actual center of the image is where the target is
            try:
                all_star_px, all_star_py = target_w.all_world2pix(all_star_ra, all_star_dec, 1)
            except:
                difference_in_ra = mini_ra - all_star_ra # Degrees
                difference_in_dec = mini_dec - all_star_dec # Degrees
                all_star_px = 200 - (difference_in_ra * 3600 / 1.86 * np.cos(all_star_dec)) # Conversion from degree differences to pixels. 3600 arcseconds in a degree. Each pixel is 1.86"
                all_star_py = 200 - (difference_in_dec * 3600 / 1.86 / np.sin(all_star_ra)) # Conversion from degree differences to pixels. 3600 arcseconds in a degree. Each pixel is 1.86"
        
        print("Target pixels (ra and dec converted to pixels): (" + str(all_star_px[mini_target_star_index]) + ", " + str(all_star_py[mini_target_star_index]) + ")")
        
        # Calculating the offsets with the gaussian fit
        print("Calculating guassian fit")
        #print(all_star_px[mini_target_star_index], all_star_py[mini_target_star_index])
        fits_image_without_median, subarray_length, x0, sigma_x, y0, sigma_y, theta = get_offset(fits_image, float(all_star_px[mini_target_star_index]), float(all_star_py[mini_target_star_index]), mini_sigma_coefficient)
        print("Saving offsets")
        add_new_frame(time_fits_header, fits_URL, np.median(fits_image), fits_header["filter"], exposure_fits_header, fits_header["airmass"], x0, sigma_x, y0, sigma_y, theta, subarray_length, all_star_data)
        
        # Getting all the brightnesses
        print("Getting other stars' brightnesses")
        all_star_from_cone_brightness, all_star_from_cone_max_pixel = get_brightness(fits_image, subarray_length, all_star_px, all_star_py, x0, sigma_x, y0, sigma_y, theta, mini_is_ATLAS)
        print("Saving other stars")
        for i in range(len(all_star_from_cone_brightness)):
            all_star_data[-1][str(mini_cone.iloc[i]["ID"]) + "_brightness"] = all_star_from_cone_brightness[i]
            all_star_data[-1][str(mini_cone.iloc[i]["ID"]) + "_px"] = all_star_px[i]
            all_star_data[-1][str(mini_cone.iloc[i]["ID"]) + "_py"] = all_star_py[i]
            all_star_data[-1][str(mini_cone.iloc[i]["ID"]) + "_max_pixel"] = all_star_from_cone_max_pixel[i]
            all_star_data[-1][str(mini_cone.iloc[i]["ID"]) + "_distance_to_real_target"] = get_distance_from_target(float(all_star_ra[int(np.where(mini_cone["ID"] == mini_tic_ID)[0][0])]), float(all_star_dec[int(np.where(mini_cone["ID"] == mini_tic_ID)[0][0])]), float(mini_cone.iloc[i]["ra"]), float(mini_cone.iloc[i]["dec"]))
            all_star_data[-1][str(mini_cone.iloc[i]["ID"]) + "_position_angle"] = float(str(astropy.coordinates.position_angle(float(all_star_ra[int(np.where(mini_cone["ID"] == mini_tic_ID)[0][0])]) * float(np.pi) / float(180), float(all_star_dec[int(np.where(mini_cone["ID"] == mini_tic_ID)[0][0])]) * float(np.pi) / float(180), float(mini_cone.iloc[i]["ra"]) * float(np.pi) / float(180), float(mini_cone.iloc[i]["dec"]) * float(np.pi) / float(180))).replace("rad", ""))
        
        # Saving the dataframe to a .csv file
        all_star_data[-1]["processing_time"] = time.time() - start_frame_time
        all_star_dataframe = pd.DataFrame(all_star_data)
        all_star_dataframe.to_csv("TICs/" + str(mini_tic_ID) + "/" + "All_" + str(mini_tic_ID) + "_Star_Data.csv")
    else:
        figure, ax = plot.subplots()
        figure.set_figheight(20)
        figure.set_figwidth(20)
        img = ax.imshow(np.log(fits_image - np.min(fits_image) + 1), cmap = "Greys")
        
        # Patches of the new and old targets
        print("Original target (" + str(mini_tic_ID) + ") location (degrees): (" + str(mini_ra) + ", " + str(mini_dec) + ")")
        print("New target (" + str(mini_cone.at[mini_target_star_index, "ID"]) + ") location (degrees): (" + str(float(mini_cone.at[mini_target_star_index, "ra"])) + ", " + str(float(mini_cone.at[mini_target_star_index, "dec"])) + ")")
        
        old_px, old_py = target_w.wcs_world2pix(float(mini_ra), float(mini_dec), 1)
        new_px, new_py = target_w.wcs_world2pix(float(mini_cone.at[mini_target_star_index, "ra"]), float(mini_cone.at[mini_target_star_index, "dec"]), 1)
        standard_cutout_patch = matplotlib.patches.Circle((len(fits_image[0]) / 2, len(fits_image) / 2), 150, fill = False, color = "white")
        if mini_is_ATLAS: # Need to correct for offset due to fits WCS irregularity. The actual center of the image is where the target is
            old_px, old_py = target_w.all_world2pix(float(mini_ra), float(mini_dec), 1)
            new_px, new_py = target_w.all_world2pix(float(mini_cone.at[mini_target_star_index, "ra"]), float(mini_cone.at[mini_target_star_index, "dec"]), 1)
            standard_cutout_patch = matplotlib.patches.Circle((len(fits_image[0]) / 2, len(fits_image) / 2), 150 / 1.86, fill = False, color = "white")
        
        old_target_patch = matplotlib.patches.Circle((old_px, old_py), 10, fill = False, color = "orange")
        print("Original target (" + str(mini_tic_ID) + ") location (pixels): (" + str(old_px) + ", " + str(old_py) + ")")
        ax.add_patch(old_target_patch)
        new_target_patch = matplotlib.patches.Circle((new_px, new_py), 10, fill = False, color = "magenta")
        print("New target (" + str(mini_cone.at[mini_target_star_index, "ID"]) + ") location (pixels): (" + str(new_px) + ", " + str(new_py) + ")")
        ax.add_patch(new_target_patch)
        ax.add_patch(standard_cutout_patch)
        
        ax.legend([old_target_patch, new_target_patch, standard_cutout_patch], [("Old Target: " + str(mini_tic_ID)), ("New Target: " + str(mini_cone.at[mini_target_star_index, "ID"])), "2.5' Cutout"])
        
        
        # Patches for rest of stars that are bright enough to produce signal
        if len(mini_cone) < 30:
            for comparison_star_index in range(len(mini_cone["ID"])):
                comparison_cutout_patch_x, comparison_cutout_patch_y = target_w.wcs_world2pix(float(mini_cone["ra"][comparison_star_index]), float(mini_cone["dec"][comparison_star_index]), 1)
                if mini_is_ATLAS:
                    comparison_cutout_patch_x, comparison_cutout_patch_y = target_w.all_world2pix(float(mini_cone["ra"][comparison_star_index]), float(mini_cone["dec"][comparison_star_index]), 1)
            
                comparison_cutout_patch = matplotlib.patches.Circle((comparison_cutout_patch_x, comparison_cutout_patch_y), 5, fill = False, color = "white")
                ax.add_patch(comparison_cutout_patch)
        
        
        # Colorbar
        fits_image = np.log(fits_image - np.min(fits_image) + 1)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size = "5%", pad = 0.05)
        colorbar = figure.colorbar(img, cax = cax, orientation = "vertical")
        colorbar.set_ticks([np.mean(fits_image), np.median(fits_image), np.max(fits_image), np.min(fits_image)])
        colorbar.set_ticklabels(["Mean: " + str(np.mean(fits_image)), "Median: " + str(np.median(fits_image)), "Max: " + str(np.max(fits_image)), "Min: " + str(np.min(fits_image))])
        
        plot.show()
        
        pdb.set_trace()
    
    fits_file.close()
    print("Success!\n\n")


def download_ZTF_images(mini_tic_ID, mini_ra, mini_dec, mini_cutout_size, mini_is_testing_images): # Download the ZTF image from ztfquery
    # Get the metadata
    print("Getting the ZTF metadata for a " + mini_cutout_size + " arcsecond cutout")
    zquery = query.ZTFQuery()
    zquery.load_metadata(kind = "sci", radec = [mini_ra, mini_dec], size = (int(mini_cutout_size) / 3600)) # Degrees
    zquery.metatable.to_csv("TICs/" + str(mini_tic_ID) + "/" + str(mini_tic_ID) + "_Metadata.csv")
    metatable = pd.read_csv("TICs/" + str(mini_tic_ID) + "/" + str(mini_tic_ID) + "_Metadata.csv", index_col = 0)
    zquery_reloaded  = query.ZTFQuery(metatable = metatable, kind = "sci")
    
    # Download the images
    print("Downloading the ZTF Images")
    if mini_is_testing_images:
        zquery.download_data("sciimg.fits", nprocess = 4, show_progress = False, cutouts = True, cutout_size = (int(mini_cutout_size) * 2), radec = [mini_ra, mini_dec], indexes = np.arange(0, 10, 1))
    else:
        zquery.download_data("sciimg.fits", nprocess = 4, show_progress = False, cutouts = True, cutout_size = (int(mini_cutout_size) * 2), radec = [mini_ra, mini_dec])
    
    return zquery.get_local_data("sciimg.fits") # The local downloaded data



def setup(tic_ID, is_testing, is_processing, is_ATLAS, sigma_coefficient = 1, new_fitter = 0): # Running everything together
    '''
    Arguments:
        tic_ID -- TIC ID of the target to run.
        is_testing -- If is_testing, only the first 10 images will run (otherwise the  entirety of dataset's image availibility will run). Default for is_testing is set to False
        is_processing -- If is_processing, lightcurve brightness data points will be extracted (otherwise image data will be displayed for visual help). Default for is_processing is set to True
        is_ATLAS -- If is_ATLAS, data from ATLAS will be downloaded and analyzed (otherwise ZTF data will be analyzed). ATLAS data will automatically run if the declination is less than -28 degrees as that is the current bound for ZTF's telescope imaging. Default for is_ATLAS is set to False
        sigma_coefficient -- ratio sizing of the aperture as sometimes a smaller or larger aperture than the mpyfit fit ellipse is preferred. Default for sigma_coefficient is set to 1
        new_fitter -- if the automaticly chosen star for mpyfit aperture fitting is saturated or too close to another star. When new_fitter is not 0, if it is a valid other TIC ID within the star field, it will fit using that chosen star
    '''
    
    
    start_time = time.time() # Starting the timer
    
    
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning) # Removing HTTPS warning
    
    # Getting the cutout size
    cutout_size = "150" # 150 arcseconds = 2.5 arcminutes
    radius = str(float(cutout_size) * float(math.sqrt(2)) / float(3600)) # Search radius. Cutout_size * root 2 / 2 for each side. Needs to be in decimal degrees 1/120 degrees. Root 2 is for the corners of the image and adds a slight buffer
    #radius = float(cutout_size) / float(3600) # Search radius (not in a square)
    ra, dec, transit_depth = get_TOI_info(tic_ID)
    
    
    # Creating the folders
    if not os.path.exists("TICs/" + str(tic_ID)):
        os.makedirs("TICs/" + str(tic_ID))
    if not os.path.exists("Data/Atlas/" + str(tic_ID)):
        os.makedirs("Data/Atlas/" + str(tic_ID))
    
    
    # Checking if I need to force use ATLAS
    if dec < -28: # ZTF limits are -28 degrees dec
        print("dec is past -28 (" + str(dec) + "), which is too far South for ZTF processing. Automatically downloading ATLAS images...")
        is_ATLAS = True
        sys.exit()
    
    
    # Cone
    cone = tic_advanced_search_position_rows(ra, dec, radius) # get_cone(ra, dec, radius)
    old_cone_length = len(cone)
    transit_depth = transit_depth / 1000000 # Converting transit depth from ppm to a fraction
    delta_m = (np.log(1 / transit_depth) / np.log(2.512)) + 0.5 # Positive number
    cone = cone[(cone["Tmag"] < (cone.at[int(np.where(cone["ID"] == tic_ID)[0]), "Tmag"] + delta_m)) | (cone["ID"] == tic_ID)]
    print(str(len(cone)) + " stars are usable out of " + str(old_cone_length))
    cone = cone.reset_index(drop = True)
    if not os.path.exists("TICs"):
        os.makedirs("TICs")
    if not os.path.exists("TICs/" + str(tic_ID)):
        os.makedirs("TICs/" + str(tic_ID))
    cone.to_csv("TICs/" + str(tic_ID) + "/" + str(tic_ID) + "_Cone.csv")
    
    
    all_star_data = []
    
    
    # Finding the target star index from the cone
    target_star_index = int(np.where(cone["ID"] == tic_ID)[0]) # Where the target is
    print("Target star (TIC " + str(cone.at[int(np.where(cone["ID"] == tic_ID)[0]), "ID"]) + ") is at index " + str(target_star_index))
    
    
    # Getting frame length and fitter star
    frames = 0
    frames_length = 10
    if is_ATLAS:
        sigma_coefficient = 0.5
        frames = [("Data/Atlas/" + str(tic_ID) + "/" + filename) for filename in os.listdir("Data/Atlas/" + str(tic_ID)) if filename.endswith(".fits")]
    else: # ZTF data
        if float(cone.at[target_star_index, "Tmag"]) < 14: # If the target is saturated
            potential_replacement_targets = cone[(cone["Tmag"] < (cone.at[int(np.where(cone["ID"] == tic_ID)[0]), "Tmag"] + delta_m)) & (cone["Tmag"] > 14) & (cone["Tmag"] < 15)].copy() # Good magnitude for targets
            potential_replacement_targets["distances_to_target"] = get_distance_from_target(([ra] * len(potential_replacement_targets["ra"])), ([dec] * len(potential_replacement_targets["dec"])), potential_replacement_targets["ra"], potential_replacement_targets["dec"])
            
            found_new_fitter = False
            while found_new_fitter == False:
                target_star_index = int(potential_replacement_targets[potential_replacement_targets["distances_to_target"] == min(list(potential_replacement_targets["distances_to_target"]))].index.tolist()[0]) # Choose the closest fitter as a potential new target
                
                # Get all the stars within 10" of the potential new fitter
                nearby_stars_cone = cone[(cone["Tmag"] < (cone.at[int(np.where(cone["ID"] == tic_ID)[0]), "Tmag"] + delta_m))] # Good magnitude for targets
                nearby_stars_cone["distances_to_target"] = get_distance_from_target(([cone.at[target_star_index, "ra"]] * len(cone["ra"])), ([cone.at[target_star_index, "dec"]] * len(cone["dec"])), cone["ra"], cone["dec"])
                nearby_stars_cone = nearby_stars_cone[(nearby_stars_cone["distances_to_target"] < 20)]
                count = 0
                for nearby_star in nearby_stars_cone["ID"]:
                    nearby_star_index = int(np.where(cone["ID"] == nearby_star)[0])
                    if cone.at[nearby_star_index, "Tmag"] < 14:
                        count += 1
                
                if count > 0:
                    print("Fitter " + str(cone.at[target_star_index, "ID"]) + " may be contaminated. Finding another fitter")
                    found_new_fitter = False
                    potential_replacement_targets = potential_replacement_targets.drop([target_star_index])
                else:
                    found_new_fitter = True
            
            print("New target star (TIC " + str(cone.at[target_star_index, "ID"]) + ") is at index " + str(target_star_index))
        
        frames = download_ZTF_images(tic_ID, ra, dec, cutout_size, is_testing) # Getting the links
    
    if not new_fitter == 0:
        if new_fitter in list(cone["ID"]):
            target_star_index = int(np.where(cone["ID"] == new_fitter)[0]) # Just assigning it to a chosen star from MAST
            print("Chosen target star (TIC " + str(cone.at[target_star_index, "ID"]) + ") is at index " + str(target_star_index))
        else:
            print("New fitter not in cone.")
            sys.exit()
    
    if not is_testing:
        frames_length = len(frames)
    
    for frame_index in range(frames_length):
        if not is_testing:
            clear_output(wait = True) # Clearing the output so the jupyter doesn't crash
            print("[" + ("=" * int(frame_index / frames_length * 10)) + (" " * (10 - int(frame_index / frames_length * 10))) + "] " + str(int(frame_index / frames_length * 10) * 10) + "%")
        
        print("Frame " + str(frame_index) + " of " + str(frames_length - 1) + " frames")
        fits_URL = frames[frame_index]
        
        # Opening and analyzing the .fits file
        try:
            open_fits_file(tic_ID, frame_index, fits_URL, cone, target_star_index, ra, dec, is_processing, is_ATLAS, sigma_coefficient, all_star_data)
        except ValueError as e: # Just needs to be except to pass correctly
            print("Failed and skipped: " + str(e))
            print("..................\n\n")
            pass
    
    
    if is_processing:
        print("Converting JD to BJD")
        all_star_dataframe = pd.read_csv("TICs/" + str(tic_ID) + "/" + "All_" + str(tic_ID) + "_Star_Data.csv", index_col = 0)
        non_converted_times = list(all_star_dataframe["jd_time"])[:frames_length] # Getting the JD times
        
        # Converting JD to BJD times (the PHP request can only handle 400 max times at once)
        starting_index = 0
        total_indexes = len(non_converted_times)
        mod_length = int(total_indexes % 400)
        full_converted_times = []
        
        if not (total_indexes - mod_length) == 0:
            for i in range(int((total_indexes - mod_length) / 400)):
                full_converted_times += convert_BJD_and_JD(str(non_converted_times[starting_index : starting_index + 400]).replace("[", "").replace("]", ""), float(ra), float(dec), False)
                starting_index += 400
        
        if not mod_length == 0: # Less than 400
            full_converted_times += convert_BJD_and_JD(str(non_converted_times[starting_index:]).replace("[", "").replace("]", ""), float(ra), float(dec), False)
        
        # Saving the fully converted list into the all star dataframe
        all_star_dataframe = all_star_dataframe.drop("jd_time", axis = 1) # 1 means column, 0 means row
        all_star_dataframe.insert(0, "time", full_converted_times)
        all_star_dataframe.to_csv("TICs/" + str(tic_ID) + "/" + "All_" + str(tic_ID) + "_Star_Data.csv")
    
    
    print("\n\n\nDone!")
    print("Time elapsed: " + str(datetime.timedelta(seconds = (time.time() - start_time)))) # Printing the change in time and formatting it in a more readable format than just seconds
