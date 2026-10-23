# TMDRamanPlotting #
*Analysis of Raman/Photoluminescence spectra of transition metal dichalcogenides (currently only sulfides)*

Python script for fitting Raman spectra from map data used to gain statistical understanding, and wide range spectra showing Raman peaks and PL response.

Included in the repository are the following:
* Instructions for use
* Main program file
* Functions required for primary program

## Dependencies ##

The current version is verified to work on Python >=3.9.
A standard anaconda installation includes most required packages [click here for installation](https://docs.anaconda.com/free/anaconda/install/), however lmfit and tqdm will need to be installed.

#### Required Packages ####

* numpy
* matplotlib
* scipy
* lmfit (not included in Anaconda)
* tqdm (not included in Anaconda) 

## Installation ##

#### Python Code ####

* Clone the repository `git clone https://github.com/SGoldie4/A2DfromUV/program.git`
* Alternatively, copy `FolderAnalysis_Mk2.py` and `FunctionBank.py` to a local folder
* Run `FolderAnalysis_Mk2.py` in your chosen python environment

## Usage ##

#### Input Data ####
Place all measurement files in a single directory.

Raman Files:
```text
SampleName_map.txt
```

Photoluminescence Files:
```text
SampleName_PL.txt
```

#### Workflow ####
1. Launch the software.
2. Select a folder containing Raman and/or PL datasets.
3. The software used file name conventions to identify files.
4. Spectra are fitted and statistical checks run.
5. Results and figures are written to the output directory.

#### Output Structure



### Features - Raman Analysis ###
Automatically fits Raman spectra using:

- Linear background correction
- $\mathrm{E}_{2g}^{1}$ mode
- $\mathrm{A}_{1g}$ mode
- 2LA(M) mode
- Silicon reference peak (when appropriate)

The software will include a silicon calibration peak when doing so improves fit quality, decided automatically with internal fit comparisons.

### Features - Photoluminescence (PL) Analysis
Automatically fits PL spectra using:

- Baseline offset
- A exciton (XA)
- B exciton (XB)
- Trion (XT⁻, optional)


Optional fit components are only included when they significantly improve the quality of the fit.

__Disclaimer:__ this software was developed to automate the data analysis of Raman spectra in line with current literature best-practice as understood by the authors at time of publication. Access is provided free of charge, but the authors cannot guarantee functionality and accuracy across different data sets.
