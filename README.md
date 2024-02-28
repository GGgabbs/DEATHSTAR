# DEATHSTAR
## **D**etecting and **E**valuating **A** **T**ransit: finding its **H**idden **S**ource in **T**ime-domain **A**rchival **R**ecords
### A system for confirming planets and identifying false positive signals in TESS data using ground-based time domain surveys

![](README_Assets/DEATHSTAR_gif.gif)

### Created by: **Gabrielle Ross**
### Last updated: **2/19/2024**

Please see this **Google Doc** for the most up-to-date documentation: **[https://docs.google.com/document/d/1XhSLHx4Errv8sN3Wgqgwl7IM7kJosu0pbtIqBg6fUhQ/edit](https://docs.google.com/document/d/1XhSLHx4Errv8sN3Wgqgwl7IM7kJosu0pbtIqBg6fUhQ/edit)**

Please **cite our paper** if you use this code: **[https://ui.adsabs.harvard.edu/abs/2023MNRAS.tmp.3722R/abstract](https://ui.adsabs.harvard.edu/abs/2023MNRAS.tmp.3722R/abstract)**


---

## How DEATHSTAR Kills Planets:

1. 


---

## Downloading DEATHSTAR:

1. **Install [Anaconda](https://www.anaconda.com/download)** on your computer as this package and its dependencies will be installed inside of a conda environment
2. DEATHSTAR is a **complete pipeline** available on **[Github](https://github.com/GGgabbs/DEATHSTAR/tree/main)**

      Go to the GitHub page and download the code as a **.zip file**
3. **Unzip DEATHSTAR** and its contents:

     Unzip and extract the files wherever your project’s code is on your computer. This means that they need to share the same directory when running your own code (or change the path)!
     You will either be able to run DEATHSTAR in **your own .py file** as shown in the example **jupyter notebook**.
4. Creating the **conda environment** and **installing dependencies**:

     Open your **Anaconda Prompt** (for Windows) or **Terminal** (for Mac). This is important because this is what has conda installed
     Type in and run `conda create -n DEATHSTAR python=3.9.15 numpy matplotlib scipy` in your Anaconda Prompt/ Terminal to create the environment
     Activate the DEATHSTAR environment using `conda activate DEATHSTAR`
     **Install dependencies** using `pip install ztfquery` and then subsequently `pip3 install pandas astropy astroquery photutils pyastronomy fpdf ipython notebook`
     If your program uses **additional dependencies**, use `pip3 install [PACKAGE NAME]` to install them
5. Logins and accounts for datasets:

     **[Zwicky Transient Facility (ZTF)](https://irsa.ipac.caltech.edu/frontpage/)** login: in order to retrieve ZTF data, you need a login on their website (you only need to input your username and password the first time you run the code)
6. Opening DEATHSTAR:

     Navigate to the DEATHSTAR project folder (wherever you have extracted it) within your anaconda prompt using `cd [FOLDER NAME]` and replace `[FOLDER NAME]` with your own directory name
     Test package installation with the download program `Test.py` using `python Test.py` OR `python3 Test.py`
   
     **Note:** If you get an error saying `there is no fpdf module`, deactivate the DEATHSTAR conda environment using `conda deactivate`, install fpdf using `pip install fpdf` (which installs in general conda base environment), then reactivate DEATHSTAR using `conda activate DEATHSTAR` and then rerun `python Test.py` OR `python3 Test.py`
   
     This will prompt you to fill in your ZTF login information you just created. Then the program will go through the full extracting and plotting process in 1 go as a complete pipeline example
     In order to view example outputs, open a new Anaconda Prompt (or normal Terminal on Mac) via navigating to the project folder using `cd [FOLDER NAME]` and then activating the DEATHSTAR conda environment using `conda activate DEATHSTAR`
   
     Open the DEATHSTAR_Example.ipynb using the following command `jupyter notebook`
   
     Go to the browser where the Jupyter notebook has opened and open the .ipynb file
   
     **Note:** Jupyter will open in the browser window that you last used!

### Now this battle station is fully operational!


---

## Planet Murder with DEATHSTAR:

1. 


---
