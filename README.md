# HopPyBar

## Overview
A python program to import, analyze, and export Split-Hopkinson Pressure Bar (Kolsky Bar) data.

This program serves as a developmental platform to "whitebox" the data analysis process (in constrast to similar "blackbox" programs).

Data streams can be captured (in-situ) to enable advanced or unconventional analyses, statistics, and comparisons. General workflow will import SHPB data in one of a number of formats, identify relevant portions of data signals, and convert to stress-strain-strainrate to show material behavior as a function of dynamic testing. Dispersion correction functionality is built in. Several export options exist, with an option to add additional metadata for easier import to Granta (or similar) database. 

## Required Distribution/Packages
HopPyBar should work on any modern python distribution, and development currently uses Python 3.11. Older versions of Python3 or required packages should work in most cases (and may throw a deprecation warning in the command line output). Otherwise, common packages are used throughout to make installation and use easier with conda/macports/homebrew/etc. Generally, most of the heavy lifting is done by matplotlib, numpy, pandas, and PyQt5. Occasionally, there are non-standard (but still common) required packages (e.g. openpyxl for output to excel spreadsheets). 

## Usage
To use the program, just run the HopPyBar.py main file and follow the prompts. A window will pop up to allow for selection of input file and other options. Sometimes, plots will pop up which require action (prompting a user to select a window of data). For example, to select the incident pulse, a plot will appear with the imported data on the left and a zoomed view on the right. Click and drag to select the incident pulse, at which point the zoomed view will update. If changes are desired at this point, a user can click and drag on the left again, until they are happy with the selection. To continue, simply close the plot window and the program will continue. Not every pop-up requires an action... Some are informational (like the plot of 1-, 2-, and 3-wave results). Some are optional (e.g. the background selection, transmitted pulse selection, or pulse timing adjustment plots can all accept user input, but if the user closes the plot without making a selection a default behavior is used instead). Once a user reaches the end of the program, it will terminate on it's own. To run a new file (or re-analyze a previous dataset), just run the main program again.

## Major Contributions
Created by Ben Morrow, Los Alamos National Laboratory. Originally based loosely on a previous analysis routine of Chris Meredith (ARL). Major contributions from Ginny Euser (LANL). Dispersion correction based on work of Anatolii Bragov (NNSTU). 

## Other Info
Full license info can be found [here](LICENSE)

© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are.
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare.
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit.
others to do so.

Approved for Public Release: C# O4628