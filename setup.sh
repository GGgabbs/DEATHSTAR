#
#        # # # # #
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
# Filename:     setup.sh
# Created by:   Gabrielle Ross
# Last updated: 2/10/2024
# Github:       https://github.com/GGgabbs/DEATHSTAR/tree/main
#
# DEATHSTAR:    A system for confirming planets and identifying
#               false positive signals in TESS data using 
#               ground-based time domain surveys
#
# Purpose:      Creates conda virtual environment 
#               DEATHSTAR_environment and installs all 
#               dependencies for DEATHSTAR


#!/bin/bash


if [[ $(conda env list) == *"DEATHSTAR_environment"* ]]; then
	echo "DEATHSTAR environment already exists!";
	conda activate DEATHSTAR_environment;
else
	conda create --name DEATHSTAR_environment python=3.9.15 numpy matplotlib;
	
	conda activate DEATHSTAR_environment;
	
	python3 -m pip install pandas;
	python3 -m pip install astropy;
	python3 -m pip install astroquery;
	
	python3 -m pip install ztfquery;
	
	python3 -m pip install photutils;
	python3 -m pip install scipy;
	
	python3 -m pip install PyAstronomy;
fi


python3 Module_Test.py;
read