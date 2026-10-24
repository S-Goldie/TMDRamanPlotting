# TMDRamanPlotting #
*Analysis of Raman/Photoluminescence spectra of transition metal dichalcogenides (currently only sulfides)*

Python script for fitting Raman spectra from map data used to gain statistical understanding, and wide-range spectra showing Raman peaks and PL response.

Included in the repository are the following:
* Instructions for use
* Main program file
* Functions required for the primary program

## Dependencies ##

The current version is verified to work with Python >=3.9.
A standard Anaconda installation includes most required packages [click here for installation](https://docs.anaconda.com/free/anaconda/install/); however, lmfit and tqdm will need to be installed.

#### Required Packages ####

* numpy
* matplotlib
* scipy
* lmfit (not included in Anaconda)
* tqdm (not included in Anaconda)

## Installation ##

#### Python Code ####

* Clone the repository: `git clone https://github.com/S-Goldie/TMDRamanPlotting.git`
* Alternatively, copy `FolderAnalysis_Mk2.py` and `FunctionBank.py` to a local folder
* Run `FolderAnalysis_Mk2.py` in your chosen Python environment

## Usage ##

#### Input Data ####
Place all measurement files in a single directory. This directory will also be used to save the output files.
The laser wavelength is defined within `FolderAnalysis_Mk2.py` as `lam = 532 #nm`. If other laser wavelengths are used, this should be changed.

Raman File:
```text
SampleName_map.txt
```
*If multiple maps are present, the largest one will be used.*
*Map data must currently be in Renishaw's standard format, with columns for x coordinate, y coordinate, Raman shift, and intensity. Other file formats may be supported in the future.*

Photoluminescence Files:
```text
SampleName_PL-#.txt
```
*Multiple spectra can be analysed simultaneously for comparison; use '#' to number the spectra.*
*Save each PL spectrum as an individual two-column x,y file, with the Raman shift in the first column.*

#### Workflow ####
1. Launch the software.
2. Select a folder containing Raman and/or PL datasets.
3. The software uses file name conventions to identify files.
4. Spectra are fitted and statistical checks are run.
5. Results and figures are written to the output directory.

#### Output Files
- map_fitting.csv = all fitted peak parameters from each spectrum in the Raman map
- map_spectra.png = overlay plot showing all Raman map spectra individually. The modal spectrum is highlighted in darker colour and with a thicker line width.
- normalised_spectra.csv = all PL spectra combined and normalised to the $\mathrm{A}_{1g}$ intensity.
- PL.png = combined plot of all normalised PL spectra, rescaled to energy.
- PL_data.csv = all PL spectra combined and normalised to the $\mathrm{A}_{1g}$ intensity, plotted against photon energy rather than Raman shift wavenumber.
- PL_fitting.csv = all fitted peak parameters from each PL spectrum
- representative_spectrum.csv = a single Raman spectrum most representative of the Raman map. If the normality test was past, this is closest to the mean peak spacing. If the data is not normally distributed, the spectrum closest to the modal peak spacing is used.

### Features - Raman Analysis ###
Automatically fits Raman spectra using:

- Linear background correction
- $\mathrm{E}_{2g}^{1}$ mode
- $\mathrm{A}_{1g}$ mode
- 2LA(M) mode
- Silicon reference peak (when appropriate)

The exact model used is defined dynamically; only when additional peaks improve the fit above random noise is a more complex model used.
Statistical analysis of the fitted map data is displayed in the Python interpreter window. A Shapiro-Wilk normality test is automatically applied, and a warning is given if the data set is non-normal. Summary statistics such as the mean and standard deviation are returned, but the user should use their discretion when reporting such average values from non-normally distributed data sets.

### Features - Photoluminescence (PL) Analysis
Automatically fits PL spectra using:

- Baseline offset
- A exciton (XA)
- B exciton (XB)
- Trion (XT⁻, optional)

Optional fit components are only included when they significantly improve the quality of the fit.

__Disclaimer:__ This software was developed to automate the data analysis of Raman spectra in line with current literature best practice as understood by the authors at the time of publication. Access is provided free of charge, but the authors cannot guarantee functionality and accuracy across different data sets.
