#!/usr/bin/python
# Import, analyze, export SHPB data

# © 2023. Triad National Security, LLC. All rights reserved.
# This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), 
# which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. 
# All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
# Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, 
# irrevocable worldwide license in this material to reproduce, prepare derivative works, distribute copies to the public, 
# perform publicly and display publicly, and to permit others to do so.

# For full list of previous versions, see VERSION_HISTORY.md
version = '20230628'    # First open-source release!!! Added banners and updated licensing info.
                        # Added dispersion correction method select (including "None"); 
                        # Organized GRANTA options in separate folder, no more hard-coding dropbox values;
                        # Possible BUG: FLAG might have variable timestep, which would screw with things if .raw file were used.

#%%###########################################

import unitconversion  # separate homebrew python module for unit conversions
import output_plots
import output_excel
import data_funcs   # various smoothing functions to be added later
import os, sys  # for filehandling and whatnot. Removed shutil, csv because we weren't using them anymore
from glob import glob   # for filehandling again, but differently
from datetime import datetime  # need to convert unix timestamps to date & time
import numpy as np
import pandas as pd # also required openpyxl to be installed for verbose features, imported below if needed
import matplotlib
matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector, Slider, Button
# dispersion correction import happens later if needed

#%%##########################################
#         INPUT CONSTANTS                ####
#############################################
# Standard Units:
# Young's Modulus (E): GPa
# Sound/Wave Speed (c): m/s
# Density (rho): kg/m^3
# Lenth/Dimensions: mm
# Excitation: V

# Unit conversion functions available in module: unitconversion.py

now = datetime.now()
date = now.strftime('%Y%m%d')
cwd = os.getcwd()   # fetch the working directory. Eventually we want to save the dispersion correction calcs for each bar in there.

#%%###########################################
### Open GUI from separate file (input_ui.py) to choose dataset to analyze (and set test platform for import style). This is read in later down the program.
### GUI checks if analysis folder already exists, handles overwrite options if applicable
##############################################
import input_ui
app = input_ui.QApplication(sys.argv)
import_ui = input_ui.ImportData()
import_ui.show()
app.exec()

platforms = {   # variable name (short): Fullname
    'mini':'Mini-Kolsky',
    'LANL':'LANL HE/Main Bar',
    'FLAG':'LANL FLAG Output',
    'APS':'APS/DCS'}

#fileid = input_ui.fileid
fileid, platform, overwrite, dispersion_select, suppress_output, suppress_granta , eng_out, data_dir = import_ui.get_output()
file_path, filename = os.path.split(fileid)
#platform = input_ui.platform
#overwrite = input_ui.overwrite
if dispersion_select.currentIndex() == 0:
    dispersion_select = 'Bragov2022'
else:
    dispersion_select = 'None'
    print('WARNING: DISPERSION CORRECTION SET TO NONE! NO CORRECTIONS WILL BE CALCULATED!')

if suppress_output == 0:
    suppress_output = False
    import openpyxl # need this for verbose mode
else:
    suppress_output = True

if suppress_granta == 0:
    suppress_granta = False
    granta_dir = os.path.join('.','GRANTA') + os.sep
else:
    suppress_granta = True

if eng_out == 0:
    eng_out = False
else:
    eng_out = True

print(f'Selected Platform = {platforms[platform]}')
print(f'Dispersion Correction Method = {dispersion_select}')
print(f'Suppress Output = {suppress_output}')
print(f'Overwrite: {str(overwrite)}')
print(f'Engineering Output = {eng_out}')
print(f'Suppress Granta = {suppress_granta}')

#%%##############################################
###  RETREIVE FILES, ORGANIZE, IMPORT FROM CSV
#################################################
if platform =='mini':
    files = sorted(glob(os.path.join(file_path, '*Ch*')), key=os.path.basename)        # Couldn't get the file finding popups to work for some reason, some wacky stuff seems to happen with the widgets in Spyder, VKE - 05/2020
    # print(files)    # DEBUG ONLY
    # For standard scope (single channel per file) data. This for loop sets the pathname to account
    # for the correct channel number and loops through to grab all channels in the folder.
    df = pd.read_csv(files[0], usecols=[3], names = ['time'])               # Reads time data in from first file (ch1) as dataframe (assuming all time data is the same across channels, seems to be valid for present data), VKE - 05/2020
    df = df * 1e6                                                           # Convert time from seconds to microseconds, VKE - 05/2020
    chs = [pd.read_csv(f, usecols=[4], header = None) for f in files]       # Get channel columns from files and put into df dataframe, VKE - 05/2020
    for i in range(0,4):
        print(f'Importing data from Channel {str(i+1)} of 4.')
        df['ch' + str(i+1)] = chs[i]
    df.loc[:, 'ch1':'ch4'] = df.loc[:, 'ch1':'ch4'] * 1000                  # V to mV conversion, VKE - 05/2020
    del chs  # delete this variable just to spruce things up, since we don't need this anymore
    inputid = glob(os.path.join(file_path, '*Input*'))
    test_params = pd.read_csv(inputid[0], delimiter = '\t', usecols=[0,1], encoding = 'cp1252', header=None, skip_blank_lines=True, names=['property', 'data']).dropna()    #encoding = 'cp1252'
    timestep = float((df.time.iloc[1]-df.time.iloc[0])/1e6)
#
elif platform =='LANL':
    # Expecting 1 file (.raw)
    try:
        df = pd.read_csv(fileid, delimiter = '\t', usecols=[2,3], header = 1, names = ['ch1', 'ch2'])   # this captures one extra row, which happens to be the timestep, then jettisons that row after reading the value
    except:
        df = pd.read_csv(fileid, delimiter = '\t', encoding = 'cp1252', usecols=[2,3], header = 1, names = ['ch1', 'ch2'])

    timestep = float(df.ch1.iloc[0])
    df = df.drop(0).reset_index(drop=True)  #BUG: this wasn't overwriting before. Should be fixed now.
    df.loc[:, 'ch1':'ch2'] = df.loc[:, 'ch1':'ch2'] * 1000                  # V to mV conversion, VKE - 05/2020
    # we need a time column created from scratch because we decided we didn't need discrete time for whatever stupid reason.

    time = []
    #timestep = 1E-7
    for i in range(len(df)):
        time = np.append(time, i*timestep)
    #
    df['time'] = time * 1e6     # convert to microseconds here.
    #
    # So this is incredibly stupid, but just so we don't have to do anything to the rest of the code,
    # we will duplicate ch2 to ch3 & 4, and set ch2 to be a duplicate of ch1. This all collapses back
    # when we average the channels together to get ch_i and ch_t. Will fix someday.
    df['ch3'] = -1*df['ch2']    # the factor of -1 is so it matches the sign of the minibar data, unless the analysis gets crazy
    df['ch4'] = -1*df['ch2']    # see above
    df['ch2'] = -1*df['ch1']    # see above above
    df['ch1'] = -1*df['ch1']
    # For LANL files we can pull all the necessary parameters automatically and repopulate the variables.
    # The raw file uses Western encoding "cp1252" (Windows 1252) and we want to drop any rows with NaN.
    inputid = fileid
    test_params = pd.read_csv(inputid, delimiter = '\t', encoding = 'cp1252', usecols=[0,1], header=None,
    skip_blank_lines=False, keep_default_na=False, verbose = True, names=['property', 'data']).dropna()
#
elif platform == 'FLAG':
    #First read test parameters from the "headers" file. This is copy-pasted from above LANL/mini versions, so we should streamline later.
    inputid = glob(os.path.join(file_path, '*Input*'))
    test_params = pd.read_csv(inputid[0], delimiter = '\t', usecols=[0,1], encoding = 'cp1252', header=None, skip_blank_lines=False, names=['property', 'data']).dropna()    #encoding = 'cp1252'
    # NOTE: Typical value of c should be something on the order of 5 mm/us
    c_bi = float(test_params.data[31])  # wave speed (rod, not longitudinal or shear) in the incident bar in m/s;
    c_bt = float(test_params.data[48])  # wave speed (rod) in the transmitted bar in m/s;
    #c_bi = float(test_params.data[32])  # longitudinal wave speed in the incident bar in m/s;
    #c_bt = float(test_params.data[49])  # longitudinal wave speed in the incident bar in m/s;
    print('WARNING: FLAG data will use rod sound speed assuming 2D cylindrical geometry. Uncomment code to use longitudinal sound speed for 2D cartesian.')
    excitation_i = float(test_params.data[37])  # excitation voltage (V)
    excitation_t = float(test_params.data[54])  # excitation voltage (V)
    gain_i = float(test_params.data[16])  # calibrated gain; 1000 microstrain = 10mV, so 1mV=100 microstrain
    gain_t = float(test_params.data[17])  # calibrated gain
    # # So hopefully 10 * 100 = 1000 microstrain/V
    gauge_factor_i = float(test_params.data[24])  # Incident Gauge Correction Factor (if needed, otherwise 1); gauge factors are in ??? units (typically 2, otherwise requires Correction)
    gauge_factor_t = float(test_params.data[43])  # Transmitted Gauge Correction Factor (if needed, otherwise 1); gauge factors are in ??? units (typically 2, otherwise requires Correction)
    gcf_i = 2 / gauge_factor_i  # If gauge_factor isn't 2, it will scale everything so that gcf is correct.
    gcf_t = 2 / gauge_factor_t  # If gauge_factor isn't 2, it will scale everything so that gcf is correct.
    # Now import data itself. FLAG output creates a false strain gauge, so some factors above (like gain, voltage, etc) aren't actually used to calculate strain.
    #files = sorted(glob(os.path.dirname("".join(file_path)) + '/*dat.*'), key=os.path.basename)
    files = sorted(glob(os.path.join(file_path, '*dat.*')), key=os.path.basename)
    dataframes = [pd.read_csv(f, delim_whitespace = True, skiprows = [0], header = None,
                              names = ['cycle', 'time', 'position1', 'position2', 'vel1', 'vel2', 'pressure']) for f in files]
    #
    distance_striker_ibar = 0.2 #Distance between striker and incident bar (cm)
    striker_vel = 0.0015 #striker velocity (cm/microsec)
    time_strikercontact = (distance_striker_ibar/striker_vel)
    #
    i = dataframes[0]   #left ibar tracer
    i2 = dataframes[1]  #center ibar tracer
    i3 = dataframes[2]  #right ibar tracer
    t = dataframes[3]   #left tbar tracer
    t2 = dataframes[4]  #center tbar tracer
    t3 = dataframes[5]  #right tbar tracer
    #Time, pressure, velocity
    timestep = float((i2.time.iloc[1]-i2.time.iloc[0])/1e6)
    flg_time_i = i2['time'] - time_strikercontact   #Subtract time it takes for striker to make contact with incident bar
    flg_p_i = i2['pressure'] * 100000               #Convert from MBar to MPa
    flg_vel_i = i2['vel2'] * 10                     #Convert from cm to mm
    #
    flg_time_t = t2['time'] - time_strikercontact   #Subtract time it takes for striker to make contact with incident bar
    flg_p_t = t2['pressure'] * 100000               #Convert from MBar to MPa
    flg_vel_t = t2['vel2'] * 10                     #Convert from cm to mm
    #Position data
    flg_pos_i = i['position1'] * 10                 #Convert from cm to mm
    flg_pos_i3 = i3['position1'] * 10               #Convert from cm to mm
    flg_pos_t = t['position1'] * 10                 #Convert from cm to mm
    flg_pos_t3 = t3['position1'] * 10               #Convert from cm to mm
    #Initial "gage" length
    lo_i = flg_pos_i3[0] - flg_pos_i[0]
    lo_t = flg_pos_t3[0] - flg_pos_t[0]
    #Strain from change in gage length
    flg_e_i = (((flg_pos_i3 - flg_pos_i)-lo_i)/lo_i)*-1
    #flg_e_i = flg_e_i[flg_time_i >= 0]
    #flg_time_i = flg_time_i[flg_time_i >= 0]
    #
    flg_e_t = (((flg_pos_t3 - flg_pos_t)-lo_t)/lo_t)*-1
    #flg_e_t = flg_e_t[flg_time_t >= 0]
    #flg_time_t = flg_time_t[flg_time_t >= 0]
    #Strain from velocity
    time_i = i2['time'] - time_strikercontact
    flg_e2_i = (flg_vel_i/(c_bi))*-1
    # just so we don't have to break things later, we will convert this to voltage from gauge strain, so it can be converted back later. *sigh*
    flg_V_i = flg_e_i * (gain_i * excitation_i*1000) / (-1 * gcf_i)
    flg_V_t = flg_e_t * (gain_t * excitation_t*1000) / (-1 * gcf_t)
    #
    #df = pd.DataFrame({'time': time_i, 'ch1': -1*flg_e_i, 'ch2': -1*flg_e_i, 'ch3': -1*flg_e_t, 'ch4': -1*flg_e_t})
    df = pd.DataFrame({'time': time_i, 'ch1': flg_V_i, 'ch2': flg_V_i, 'ch3': flg_V_t, 'ch4': flg_V_t})
#
elif platform == 'APS':
    # Expecting: 1 file (so we don't need to sort files), 5 columns (time + 4 channels), with 21 row header
    # This represents collecting all the data on one scope. Some of the channels are probably bonded together (averaged) and one channel will show x-ray times.
    df = pd.read_csv(fileid, delimiter=',', skiprows=21, names = ['time', 'ch1', 'ch2', 'ch3', 'ch4'])    # here we could only skip 20 rows and let pandas read in the headers as column names, but they'd be all caps.
    df['xray']=df['ch4']    # I think this is the x-ray column.
    df['ch4']=df['ch3']     # transmitted should be physically bonded together on ch3, so we'll duplicate and it will collapse again when we average.
    timestep = float((df.time.iloc[1]-df.time.iloc[0]))
#
else:
    print('Unrecognized test platform. Exiting...')
    sys.exit()

#%%#############################################################################
### # Assign input data to variables. This only works because we have a quasi-standard input format.
###  (Either the first 2 columns of the LANL .dat file or a standalone file that approximates this).
################################################################################
E_bi = float(test_params.data[29])  # Young's Modulus of incident bar (GPa).
E_bt = float(test_params.data[46])  # Young's Modulus of transmitted bar (GPa)
# NOTE: Typical value of c should be something on the order of 5 mm/us
if platform == 'FLAG':  # For some reason FLAG uses the longitudinal sound speed instead of rod. Someday I'll figure out what's going on there.
    c_bi = float(test_params.data[32])  # longitudinal wave speed in the incident bar in m/s;
    c_bt = float(test_params.data[49])  # longitudinal wave speed in the incident bar in m/s;
else:
    c_bi = float(test_params.data[31])  # wave speed (rod, not longitudinal or shear) in the incident bar in m/s;
    c_bt = float(test_params.data[48])  # wave speed (rod) in the transmitted bar in m/s;
#
poisson_bi = float(test_params.data[30])  # assuming both bars are the same for now
poisson_bt = float(test_params.data[47])  # assuming both bars are the same for now
rho_bi = float(test_params.data[34])*1000           # bar density in kg/m^3; Used only for dispersion correction, currently.
rho_bt = float(test_params.data[51])*1000           # transmitted bar density in kg/m^3; Used only for dispersion correction.
### Specimen Dimentions
thick = float(test_params.data[6])  # specimen thickness (mm)
d = float(test_params.data[5])  # specimen dimension for square/rectangular (mm)
d_bar_i = float(test_params.data[35])  # diameter of the incident bar (mm)
d_bar_t = float(test_params.data[52])  # diameter of the transmitted bar (mm)
d_bar_i_equiv = d_bar_i
d_bar_t_equiv = d_bar_t
A_bar_i = (np.pi/4)*d_bar_i**2     # Area of round bars
A_bar_t = (np.pi/4)*d_bar_t**2     # Area of round bars
A_sample = (np.pi/4)*d**2 # Area of a round sample, if applicable
A_ratio_i = A_bar_i / A_sample  # Ratio of incident bar to sample areas
A_ratio_t = A_bar_t / A_sample  # Ratio of transmitted bar to sample areas
excitation_i = float(test_params.data[37])  # excitation voltage (V)
excitation_t = float(test_params.data[54])  # excitation voltage (V)
gain_i = float(test_params.data[16])  # calibrated gain; 1000 microstrain = 10mV, so 1mV=100 microstrain
gain_t = float(test_params.data[17])  # calibrated gain
# So hopefully 10 * 100 = 1000 microstrain/V
gauge_factor_i = float(test_params.data[24])  # Incident Gauge Correction Factor (if needed, otherwise 1); gauge factors are in ??? units (typically 2, otherwise requires Correction)
gauge_factor_t = float(test_params.data[43])  # Transmitted Gauge Correction Factor (if needed, otherwise 1); gauge factors are in ??? units (typically 2, otherwise requires Correction)
gcf_i = 2 / gauge_factor_i  # The strain equations have a 2/gauge_factor term, which is shifted here and called "gage correction factor"
gcf_t = 2 / gauge_factor_t
igauge2sample = float(test_params.data[36])  # Distance from incident gauge to sample (mm)
tgauge2sample = float(test_params.data[53])
filedate = datetime.utcfromtimestamp(os.path.getmtime(fileid)).strftime('%m/%d/%Y')  # gets the Unix timestamp of modification of the data file, and converts to readable string
filetime = datetime.utcfromtimestamp(os.path.getmtime(fileid)).strftime('%I:%M:%S %p')  # gets the Unix timestamp of mod. of the data file, converts to 12-hour time. %H instead of %I would be 24-hr time.
sysid = test_params.data[2]  # System ID (HopPyBar for this program, or HE, MB, KB for HE, Main Bar or Kolsky Bar)
test_loading = test_params.data[3]  # Tension or compression (or something else)
Sample_Description = test_params.data[4]  # Sample description
temperature = float(test_params.data[7])  # Test temperature in degrees C
breech_pressure = float(test_params.data[8])  # breech pressure in psi
len_striker = float(test_params.data[9])  # Striker length in mm
tip_matl = test_params.data[10]  # tip material
sample_lube = test_params.data[11]  # sample lubrication
sample_matl = test_params.data[12]  # sample material
calibration_id_i = test_params.data[13]  # incident bar calibration id
calibration_id_t = test_params.data[14]  # transmitted bar calibration id
c_bar_i_refl_time = float(test_params.data[26]) # Incident Bar Reflected Time Sound Speed(m/s)
#file_name = test_params.data[15]  # reads the file name (name.extension) from the file. Will be different if name is changed after.
file_name = os.path.split(fileid)[1]  # reads the file name (name.extension) from the file name itself
v_striker_desired = float(test_params.data[18])  # desired striker velocity (m/s)
# v_striker = 0                                           # actual striker velocity (m/s) (This was commented out because we set it above near the appearance of gauge strains)
gap = float(test_params.data[19])  # gap between bars?
bar_i_id = test_params.data[22]
bar_t_id = test_params.data[41]
bar_i_description = test_params.data[23]
bar_i_refl_time = float(test_params.data[25])
bar_t_description = test_params.data[42]
bar_i_gage_desc = test_params.data[28]
bar_t_gage_desc = test_params.data[45]
bar_i_calibration_date = test_params.data[27]
bar_t_calibration_date = test_params.data[44]
calibrated_gain_i = float(test_params.data[38])
calibrated_gain_t = float(test_params.data[55])
c_bar_i_long = float(test_params.data[32])#'Incident Bar Longitudinal Sound Speed(m/s)'
c_bar_i_shear = float(test_params.data[33])#'Incident Bar Shear Sound Speed(m/s)'
c_bar_t_long = float(test_params.data[49])#'Transmitted Bar Longitudinal Sound Speed(m/s)'
c_bar_t_shear = float(test_params.data[50])#'Transmitted Bar Shear Sound Speed(m/s)'
test_date = test_params.data[0]
test_time = test_params.data[1]

#%%#############################################################################
### # Data input check for unit errors. If suspect data shown, kick warning.
################################################################################
# basically trying to catch input errors from unit mismatch. Generally the number will be a few orders of magnitude different from expectations.
if rho_bi < 500:
    print('Data input warning: Input bar density less than 0.5 g/cm^3 (500 kg/m^3). Check input units.')

if rho_bt < 500:
    print('Data input warning: Transmitted bar density less than 0.5 g/cm^3 (500 kg/m^3). Check input units.')

if c_bi < 500:
    print('Data input warning: Input bar rod sound speed less than 500 m/s (0.5 mm/us). Typical values ~5000 m/s. Check input units.')

if c_bt < 500:
    print('Data input warning: Transmitted bar rod sound speed less than 500 m/s (0.5 mm/us). Typical values ~5000 m/s. Check input units.')

if E_bi > 1000:
    print('Data input warning: Incident bar modulus > 1000 GPa. Check input units.')

if E_bt > 1000:
    print('Data input warning: Transmitted bar modulus > 1000 GPa. Check input units.')

#%%#############################################################################
# Check if Analysis folder exists and create it if it doesn't. Moved here to be able to make testname a little more robust instead of using folder name.
################################################################################
testname = os.path.basename(file_path)  # set testname to folder name
# attempts to set the testname to some permutation of the date or filename
#testdatetime = ' '.join([test_date, test_time])
#testdatestamp = datetime.strptime(testdatetime, "%m/%d/%y %I:%M:%S %p")
#testname = testdatestamp.strftime("%Y%m%d_%H%M%S")+sysid
#testname = os.path.splitext(filename)[0]
print(f'Test name = {testname}')
data_dir = data_dir + os.sep   # note: os.sep should just be the OS native slash (/ for mac, \ for win)
#save_dir = data_dir + '{}_Analysis_'.format(testname) + str(date) + '/'
if suppress_output == False:
    save_dir = data_dir + 'Analysis_{}'.format(str(date)) + os.sep
    #if not os.path.exists(data_dir + 'Analysis_{}'.format(str(date))):
    #    os.makedirs(data_dir + 'Analysis_{}'.format(str(date)))
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    if overwrite == 'yes':
        print(f'Analysis data in {save_dir} will be overwritten.')
        # empty the directory before continuing. Some of the filenames depend on input, so they aren't always overwritten naturally.
        for root, dirs, files in os.walk(save_dir):
            for file in files:
                os.remove(os.path.join(root, file))
    elif overwrite =='no':
        print('Overwrite of previous data declined. Aborting data analysis.')
        sys.exit()  # this closes the whole program, stopping analysis.
    else:
        print(f'Analysis data will be written to {save_dir}')
else:
    print('Suppressing Plots/Output: Only final .dat file will be output.')

if suppress_granta == True:
    print('Suppressing GRANTA Output. No GRANTA files will be output.')
else:
    print('GRANTA output will be written to base folder with .dat file.')
    granta_wfm = pd.DataFrame({'TimeRaw(s)': df.time*1e-6, 'Incident/Reflected Waveform':(df.ch1+df.ch2)/2000, 'Transmitted Waveform':(df.ch3+df.ch4)/2000})  # have to convert back to seconds & Volts

pulse_width = (len_striker*2)/unitconversion.m_s2mm_us(c_bi)  # really it should be the sound speed in the striker, but that's not something we typically report, so we will assume the bar and striker are the same, which should get us close.
print(f'Estimated Incident Pulse Width (\u03BCs): {pulse_width:0.1f}')

# We need this at the end, but we want to capture the raw signals before we edit them.
if platform != "LANL":
    if df['ch1'].mean() > 0:    # sometimes the scope collects as opposite sign of what is expected, so check and flip if necessary
        raw_wfm = pd.DataFrame({'Incident/Reflected Waveform':(df.ch1+df.ch2)/2000, 'Transmitted Waveform':(df.ch3+df.ch4)/2000})   # grabbing a copy to reconstruct a raw file if needed
    else:
        raw_wfm = pd.DataFrame({'Incident/Reflected Waveform':-1*(df.ch1+df.ch2)/2000, 'Transmitted Waveform':-1*(df.ch3+df.ch4)/2000})
#%%###########################
### BASELINE SUBTRACTION ###
### Plots the data and allows a user to select a snippet of baseline voltage, which is then removed from all the data channels to zero them out initially
##############################
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(121, facecolor='#FFFFCC')  # 1 row, 2 columns of subfigures, subplot 1
plt.title('Baseline Sub')
plt.xlabel('Time (microsecond)')
plt.ylabel('Strain (mm/mm)')

ax.plot(df.time, df.ch1, 'g-', label='Ch1')
ax.plot(df.time, df.ch2, 'r-', label='Ch2')
ax.plot(df.time, df.ch3, 'k-', label='Ch3')
ax.plot(df.time, df.ch4, 'c-', label='Ch4')
ax.set_title('Click and Drag to Select Baseline (pre-test)')
plt.legend()
plt.grid()

ax2 = fig.add_subplot(122, facecolor='#FFFFCC')
ax2.set_title('Zoomed View')
line2, = ax2.plot(df.time, df.ch1, 'g-')

class Baseline_Select():
    """Container for baseline voltage selection"""
    def __init__(self):
        self.coords = {}
    def __call__(self, xmin, xmax):
        indmin, indmax = np.searchsorted(df.time, (xmin, xmax))
        indmax = min(len(df.time)-1, indmax)
        thisx = df.time[indmin:indmax]
        thisy = df.ch1[indmin:indmax]
        thisy2 = df.ch2[indmin:indmax]
        thisy3 = df.ch3[indmin:indmax]
        thisy4 = df.ch4[indmin:indmax]
        self.coords['time'] = thisx
        self.coords['ch1'] = thisy
        self.coords['ch2'] = thisy2
        self.coords['ch3'] = thisy3
        self.coords['ch4'] = thisy4
        self.basemax = df.time[indmax]
        self.basemin = df.time[indmin]
        self.avg1 = np.mean(thisy)
        self.avg2 = np.mean(thisy2)
        self.avg3 = np.mean(thisy3)
        self.avg4 = np.mean(thisy4)
        line2.set_data(thisx, thisy)
        ax2.set_xlim(thisx.min(), thisx.max())
        ax2.set_ylim(thisy.min(), thisy.max())
        fig.canvas.draw()

baseselect = Baseline_Select()

try:
    span = SpanSelector(ax, baseselect, 'horizontal', useblit=True, props=dict(alpha=0.5, facecolor='red') )
except:
    span = SpanSelector(ax, baseselect, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red') )
    print('Time to update matplotlib.')

plt.show(block=True)

# Rezero to remove baseline voltage
try:                                    # standard error handling, just to make the code a bit more robust.
    avg_ch1 = baseselect.avg1
    avg_ch2 = baseselect.avg2
    avg_ch3 = baseselect.avg3
    avg_ch4 = baseselect.avg4
except:                                 # If we wanted to carve out specific exceptions we could, but we'll leave it generic to just ignore the baseline voltage offset in case of error (usually because no baseline was selected)
    print('No baseline voltage selected, setting baseline offset to zero.')
    avg_ch1 = 0
    avg_ch2 = 0
    avg_ch3 = 0
    avg_ch4 = 0
    baseselect.basemax = df.time[0]          # and also set the last point in the "baseline" section to the minumum so the rest doesn't break.
    baseselect.basemin = df.time[0]          # and also set the last point in the "baseline" section to the minumum so the rest doesn't break.

basemax = baseselect.basemax
basemin = baseselect.basemin

#%%#################################################################
### THIS IS WHERE WE START MODIFYING CHANNELS, SPECIFIC TO EACH TEST
####################################################################

#Remove the baseline portions from subsequent analysis (this is mostly for convenience, just so we don't have a bunch of extra garbage in front). Ginny's version doesn't remove only the baseline snippet, just starts after end of baseline.
df = df[df.time > basemax]

# Set variable names & subtract baseline (background) voltage. Conversion from V to mV done at import.
ch1 = (df.ch1 - avg_ch1)
ch2 = (df.ch2 - avg_ch2)
ch3 = (df.ch3 - avg_ch3)
ch4 = (df.ch4 - avg_ch4)
time = df.time

# Average appropriate channels to get incident and transmitted signals
ch_i = (ch1 + ch2) / 2  # avg(channels 1 & 2) for input bar
ch_t = (ch3 + ch4) / 2  # avg(channels 3 & 4) for transmitted bar
# Note: Strain gauge factor correction (if needed) is done later during conversion to strain.

#%%######################################
### Select Incident Pulse
### This plots the incident data, allows a user to select the incident pulse, then it writes that out as dyni
### Note: reflected pulse will be calculated by this selection plus the offset time based on sound speed. There is still a chance to tweak things afterwards.
#########################################
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(121, facecolor='#FFFFCC') # 1 row, 2 columns of subfigures, subplot 1

plt.title('Incident Pulse Selection')
plt.xlabel('Time (microsecond)')
plt.ylabel('Voltage (mV)')

x = time
y = ch_i

ax.plot(x, y, 'k-', label='Incident')
ax.set_title('Click and Drag to Select Incident Pulse')
#ax.set_xlim(np.min(x), np.max(x)/2)            # just so we can zoom in a little on the useful data
plt.legend()
plt.grid()

ax2 = fig.add_subplot(122, facecolor='#FFFFCC')
ax2.set_title('Zoomed View')
line2, = ax2.plot(x, y, 'b-')
#
class Data_Range_Select():
    """Container for selection of ranges from displayed data."""
    def __init__(self):
        self.coords = {}
    def __call__(self, xmin, xmax):
        indmin, indmax = np.searchsorted(x, (xmin, xmax))
        indmax = min(len(x)-1, indmax)
        thisx = x.iloc[indmin:indmax]
        thisy = y.iloc[indmin:indmax]
        self.coords['x'] = thisx
        self.coords['y'] = thisy
        line2.set_data(thisx, thisy)
        ax2.set_xlim(thisx.min(), thisx.max())  # Used to be ax2.set_xlim(thisx[0], thisx[-1]), but this wasn't working for some reason so changed to match y, VKE - 05/2020
        ax2.set_ylim(thisy.min(), thisy.max())
        fig.canvas.draw()

dyni_select = Data_Range_Select()

try:
    span = SpanSelector(ax, dyni_select, 'horizontal', useblit=True, props=dict(alpha=0.5, facecolor='red') )
except:
    span = SpanSelector(ax, dyni_select, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red') )
    print('Time to update matplotlib.')

plt.show(block=True)                            # For interactive python, have to include block=True so that program stops until you've made plot selection and exited out of plot, VKE - 05/2020

dyni = pd.DataFrame({'time_i': dyni_select.coords['x'], 'v_i': dyni_select.coords['y']})
dyni_start = min(dyni.time_i)
dyni_end = max(dyni.time_i)

print(f'Incident Pulse Time [\u03BCs]: {dyni_start:0.3f}')

#%%######################################
### Select Transmitted Pulse
### This plots the transmitted data, allows a user to select the transmitted pulse, then it writes that out as dynt
#########################################
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(121, facecolor='#FFFFCC')  # 1 row, 2 columns of subfigures, subplot 1

plt.title('Transmitted Pulse Selection')
plt.xlabel('Time (microsecond)')
plt.ylabel('Voltage (mV)')

x = time    # not really necessary; this was already set
y = ch_t

ax.plot(x, y, 'k-', label='Transmitted')
ax.set_title('Click and Drag to Select Transmitted Pulse\n or Close Window to Use Calculated Value')
#ax.set_xlim(np.min(x), np.max(x) / 2)  # just so we can zoom in a little on the useful data
plt.legend()
plt.grid()

ax2 = fig.add_subplot(122, facecolor='#FFFFCC')  # 1 row, 2 columns of subfigures, subplot 2
ax2.set_title('Zoomed View')
line2, = ax2.plot(x, y, 'b-')

dynt_select = Data_Range_Select()
try:
    span = SpanSelector(ax, dynt_select, 'horizontal', useblit=True, props=dict(alpha=0.5, facecolor='red'))
except:
    span = SpanSelector(ax, dynt_select, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red'))
plt.show(block=True)

tt_dr_shift = (igauge2sample / unitconversion.m_s2mm_us(c_bi)) + (tgauge2sample / unitconversion.m_s2mm_us(c_bt))        # Time shift for transmitted pulse, here c is in mm/us to match gauge position unit (mm) and time unit (us)
    # Need sound speed through sample to numerically determine dynt_start - VKE, 05/2020
print(f'tt_dr_shift: {tt_dr_shift:0.3f}')

# This bit will use the calculated (dead-reckoning) transmitted pulse time if no window is selected. This is helpful for when things are really sluggish going through the material and it may be hard to deterimine timing otherwise
if not dynt_select.coords:  # Check for the dictionary of transmitted window coordinates. Empty if no window was selected by user.
    print('No transmitted selection made, using dead calculated values instead.')
    dynt_start = dyni_start + tt_dr_shift   # calculate the transmitted wave times by adding the time shift value to what was selected for incident.
    dynt_end = dyni_end + tt_dr_shift
    dynt_time = time[time>=dynt_start]  # this and the next step shave down the
    dynt_time = dynt_time[dynt_time<=dynt_end]
    dynt_v = ch_t.loc[min(dynt_time.index):max(dynt_time.index)]    # section the voltage values for just the times from the window
    dynt = pd.DataFrame({'time_t': dynt_time, 'v_t': dynt_v})
else:   # if a window is selected, use this, as usual
    dynt = pd.DataFrame({'time_t': dynt_select.coords['x'], 'v_t': dynt_select.coords['y']})
    dynt_start = min(dynt_select.coords['x'])

print(f'Transmitted Pulse Time (\u03BCs): {dynt_start:0.3f}')

#%%######################################
### Select Reflected Pulse
### This is based off of the
#########################################
## Let's try to dead-reckon (dr) the timing from the sound speed given above:
tr_dr_shift = (2 * igauge2sample) / unitconversion.m_s2mm_us(c_bi)  # time for reflected pulse to arrive at gauge (micoseconds?)
print(f'tr_dr_shift: {tr_dr_shift:0.3f}')

dynr_start = dyni_start + tr_dr_shift
print(f'Reflected Pulse Time [\u03BCs]: {dynr_start:0.3f}')

#%%######################################
### Resample data so it's all the same size
#########################################
window = max(len(dyni), len(dynt))  # Make the window as wide as the largest array. We don't care about the length of dynr, because it will be the same as dyni by definition.
print(f'Reshaping array to max length: {window}')

dyni_ind = np.searchsorted(time, dyni_start)
dynr_ind = np.searchsorted(time, dynr_start)
dynt_ind = np.searchsorted(time, dynt_start)

label_dict = {'time':'Time (microsecond)', 'v':'Voltage (mV)'} # make a label dictionary to convert between long and short names for asthetics.
dyni = pd.DataFrame({'time': time.iloc[dyni_ind:dyni_ind+window], 'v': ch_i.iloc[dyni_ind:dyni_ind+window]}) # Changed output to be a dataframe
dynr = pd.DataFrame({'time': time.iloc[dynr_ind:dynr_ind+window],'v': ch_i.iloc[dynr_ind:dynr_ind+window]})
dynt = pd.DataFrame({'time': time.iloc[dynt_ind:dynt_ind+window], 'v': ch_t.iloc[dynt_ind:dynt_ind+window]})

dyni = dyni.reset_index(drop=True) # Have to reset indices of dataframes, otherwise indices no longer start at 0
dynr = dynr.reset_index(drop=True)
dynt = dynt.reset_index(drop=True)

# Output raw data for each pulse. Rename the columns to something sane using the dictionary, then export.
if suppress_output == False:
    dyni.rename(label_dict, axis=1).to_csv(save_dir + 'dyni.txt', index=None, sep='\t')
    dynr.rename(label_dict, axis=1).to_csv(save_dir + 'dynr.txt', index=None, sep='\t')
    dynt.rename(label_dict, axis=1).to_csv(save_dir + 'dynt.txt', index=None, sep='\t')

# Plot the pulses together, at their original arrival times.
if suppress_output == False:
    output_plots.plot_pulses(dyni, dynr, dynt, save_dir)

### Output shifted strain time data.
tr_dr = dynr.time - tr_dr_shift        # these shift values were calculated earlier, in the reflected pulse code segment
tt_dr = dynt.time - tt_dr_shift
#
if suppress_output == False:
    output_plots.plot_pulses_DRtimeshift(dyni, dynr, dynt, tt_dr, tr_dr, save_dir)

#%%#######################################################
### TAKE DELTA TIME AND CONVERT FROM VOLTAGE TO STRAIN ###
# Polarity of voltages from the LANL HE bar was switched w.r.t. the mini bar to make stress positive in compression
# NOTE: This is gage strain in this section. It gets converted to more traditional strain later.
xti = dyni.time - dyni.time[0]
si = -1 * (dyni.v * gcf_i) / (gain_i * excitation_i * 1000)  # The factor of 1000 is to convert excitation from V to mV
incident = pd.DataFrame({'time': xti, 'strain': si})

xtt = dynt.time - dynt.time[0]
st = -1 * (dynt.v * gcf_t) / (gain_t * excitation_t * 1000)
transmitted = pd.DataFrame({'time': xtt, 'strain': st})

xtr = dynr.time - dynr.time[0]
sr = (dynr.v * gcf_i) / (gain_i * excitation_i * 1000)
reflected = pd.DataFrame({'time': xtr, 'strain': sr})


# Let's try to pick off the striker velocity here by detecting the rising and falling edges of incident pulse and determining the striker velocity via particle velocity.
trigger_val = 0.85 * np.max(incident.strain)  # set trigger value (like a scope), since we're doing it after the fact, we can safely use a percentage of the maximum value. It's negative and we use min instead of max because the incident signal is always negative.
### Ginny's method... sorts the dataframe for data above the trigger, then uses the first and last true values as the limit. Seems more robust against large oscillations.
xval = incident.time[incident['strain'] > trigger_val]
xmin, xmax = xval.iloc[[0, -1]]

si_plateau = incident.strain[incident.time.between(xmin, xmax)]
plateau_avg = np.mean(si_plateau)
#print(f'Plateau: {plateau_avg}')   # DEBUG ONLY

v_striker = 2 * plateau_avg * c_bi                     # striker velocity is 2 * particle velocity (= 2 * gaugestrain_incident * sound speed of incident bar). units: m/s.
print(f'Calculated Striker Velocity [m/s]: {v_striker:0.2f}')

try:
    print(f'Sytem: {sysid}')
    print(f'Striker Length (mm) {len_striker}')
    print(f'Breech Pressure (psi): {breech_pressure}')
except:
    print('System ID, Breech Pressure, or Striker Length Unknown!')

#%%#################################################
### Export gauge voltages and gauge strains to Excel. This is mostly if we want to compare directly with a model.
####################################################
if suppress_output == False:
    ## TO DO: I'm sure we don't need to make a new dataframe for this. We can probably just clean up.
    df_volt = pd.DataFrame({'Incident Time (microsecond)': dyni.time, 'Incident Gauge Voltage (mV)': dyni.v,
                            'Reflected Time (microsecond)': dynr.time, 'Reflected Gauge Voltage (mV)': dynr.v,
                            'Transmitted Time (microsecond)': dynt.time, 'Transmitted Gauge Voltage (mV)': dynt.v})
    df_gauge = pd.DataFrame({'Incident Time (microsecond)': incident.time, 'Incident Gauge Strain (mm/mm)': incident.strain,
                             'Reflected Time (microsecond)': reflected.time, 'Reflected Gauge Strain (mm/mm)': reflected.strain,
                             'Transmitted Time (microsecond)': transmitted.time, 'Transmitted Gauge Strain (mm/mm)': transmitted.strain})
    with pd.ExcelWriter(save_dir + testname + '_gaugevoltage&strain.xlsx') as book_gauge:
        df_volt.to_excel(book_gauge, sheet_name='gauge voltage', index=False)
        df_gauge.to_excel(book_gauge, sheet_name='gauge strain', index=False)

### PLOT RAW DATA
fig, ax = plt.subplots()
ax.plot(transmitted.time, transmitted.strain, 'g-', label='Transmitted')
ax.plot(incident.time, incident.strain, 'r-', label='Incident')
ax.plot(reflected.time, reflected.strain, 'b-', label='Reflected')
# ax.plot(xtt, stn_d,'k-', label='Incident-Reflected')
plt.title(testname)
plt.xlabel('Time (microsecond)')
plt.ylabel('Gauge Strain (mm/mm)')
plt.legend()
plt.grid()
if suppress_output == False:
    plt.savefig(save_dir + 'Strain-Time_noshift.png')
#plt.show(block=True)     # # DEBUG ONLY; Plotted again when we do the timeshift so this is redundant
plt.close()


#%%############################################################################
###                                                                         ###
###                 DISPERSION CORRECTION STARTS HERE                       ###
###                                                                         ###
###############################################################################
if dispersion_select == 'Bragov2022':
    from dispersionLib_Bragov2022 import barWithDispersion     # library for dispersion correction via Bragov method
    ### Bragov et al. method:
    # b1 = barWithDispersion(E = 200e9, rho = 7850, nu = 0.29, r0 = 10e-3)
    ## Here I assume E is in Pa, rho in kg/m^3, nu is Poisson ratio and r0 is radius of the bar in m.
    # disp_inc = barWithDispersion(E = E_bi(MPa)*1e6, rho = rho_bi, nu = poisson_bi, r0 = (d_bar_i_equiv/2)/1000)
    disp_inc = barWithDispersion(E = E_bi*1e9, rho = rho_bi, nu = poisson_bi, r0 = (d_bar_i_equiv/2)/1000)

    if (os.sep in calibration_id_i) or (os.sep in calibration_id_t):
        print(f'Path separator detected in calibration ID. Incident = {calibration_id_i}; Transmitted = {calibration_id_t}. Replacing with hyphen...')
        calibration_id_i = calibration_id_i.replace(os.sep, '-')   # if the calibration id has a slash, it won't be able to find the dispersion correction file later
        calibration_id_t = calibration_id_t.replace(os.sep, '-')   # if the calibration id has a slash, it won't find the dispersion correction file later
        print(f'New calibration IDs: Incident = {calibration_id_i}; Transmitted = {calibration_id_t}.')

    if not os.path.exists(os.path.join(cwd, 'dispersion_calcs')):
        os.makedirs(os.path.join(cwd, 'dispersion_calcs'))

    if not os.path.exists(os.path.join(cwd, 'dispersion_calcs', calibration_id_i + '_disp_inc.txt')):
        print('Calculating input bar dispersion corrections...')
        disp_inc.solve()
        print('Input bar dispersion correction calculations complete!')
        disp_inc.save(os.path.join(cwd, 'dispersion_calcs', calibration_id_i + '_disp_inc.txt'))
    else:
        print(f'Dispersion calculations found in {cwd}{os.sep}dispersion_calcs. Delete {calibration_id_i}_disp_inc.txt and run program again to force recalculation.')

    disp_inc.load(os.path.join(cwd, 'dispersion_calcs', calibration_id_i + '_disp_inc.txt'))

    # We will do incident and transmitted separately in case they're different.
    disp_trans = barWithDispersion(E = E_bt*1e9, rho = rho_bt, nu = poisson_bt, r0 = (d_bar_t_equiv/2)/1000)

    #if not os.path.exists(save_dir + 'disp_trans.txt'):
    if not os.path.exists(os.path.join(cwd, 'dispersion_calcs', calibration_id_t + '_disp_trans.txt')):
        print('Calculating transmitted bar dispersion corrections...')
        disp_trans.solve()
        print('Transmitted bar dispersion correction calculations complete!')
        #disp_trans.save(save_dir + 'disp_trans.txt')
        disp_trans.save(os.path.join(cwd, 'dispersion_calcs', calibration_id_t + '_disp_trans.txt'))
    else:
        print(f'Dispersion calculations found in {cwd}{os.sep}dispersion_calcs. Delete {calibration_id_t}_disp_trans.txt and run program again to force recalculation.')

    disp_trans.load(os.path.join(cwd, 'dispersion_calcs', calibration_id_t + '_disp_trans.txt'))


    #t = np.linspace(0, 0.001, 1000)     # time step
    #p = [pulse(tt, rise_time = 50e-6, fall_time = 1e-6, plateau_time = 1e-6) for tt in t]
    ## If we want original time and gage strain, we can use dyni.time with incident.strain (and similar for other pulses). The backshift could then be performed automatically, I think.
    plt.plot(incident.time, incident.strain, 'b')
    plt.title('Incident Pulse')
    plt.xlabel('time, s')
    plt.ylabel('strain')
    #plt.show(block=True)   # DEBUG ONLY
    plt.close()

    ## When we do dispersion, we want to shift things to be coincident at the sample, so we need forward dispersion of
    ## incident pulse over the length of half the incident bar, then reverse the dispersion of the reflected and transmitted pulses by
    ## their respective distances (half incident and half transmitted bar, respectively, assuming gauges are in center). We already have
    ## gauge positions igauge2sample and tgauge2sample so we can use those directly.

    # here time is in seconds, so we divide by 1e6 to convert from microsec.
    #s_inc = disp_inc.disp_shift_wave(t = (incident.time/1e6).tolist(), y = incident.strain.tolist(), dz= (igauge2sample/1000), backshift=True)  # dz is a shift factor. From Bragov2019, this looks like a DISTANCE. Here we assume the gage is halfway down the bar and this is in meters. Backshift just does a simple shift to get it back to the right time after dispersion correction. - Ben
    incident['disp_corr_strain'] = disp_inc.disp_shift_wave(t = (incident.time/1e6).tolist(), y = incident.strain.tolist(), dz= (igauge2sample/1000), backshift=True)  # dz is a shift factor. From Bragov2019, this looks like a DISTANCE. Here we assume the gage is halfway down the bar and this is in meters. Backshift just does a simple shift to get it back to the right time after dispersion correction. - Ben
    #s_inc_back = disp_inc.disp_shift_wave(t = (incident.time/1e6).tolist(), y = s_inc, dz= -1*(igauge2sample/1000), backshift=True) # let's try shifting this back to see if it looks the same as what we started with.
    #s3 = b1.simple_shift_wave(t = (incident.time/1e6).tolist(), y = incident.strain.tolist(), tshift = (igauge2sample/1000)/b1.cb)   #, 2/b1.cb is timeshift; I think this literally just shifts everything by an offset without doing anything else.
    plt.plot(incident.time, incident.strain, 'b', label='original data')
    plt.plot(incident.time, incident.disp_corr_strain, 'g', label='wave shift')
    #plt.plot(incident.time, s_inc_back, 'r', label='reverse dispersion to original')
    #plt.plot(incident.time, s3, 'r', label='simple shift')
    plt.xlabel('Time (\u03BCs)')
    plt.ylabel('Strain')
    plt.title('Incident Pulse Dispersion Correction')
    plt.legend()
    #plt.show(block=True)   # DEBUG ONLY
    plt.close()

    ## If we want original time and gage strain, we can use dyni.time with incident.strain (and similar for other pulses). The backshift could then be performed automatically, I think.
    plt.plot(dynr.time, reflected.strain, 'b')  # somehow when I ran this from code I picked up an extra row in reflected, so there was a mismatch in dimension. Note: this was below 1-wave calc before, so double check there.
    plt.xlabel('time, s')
    plt.ylabel('strain')
    #plt.show(block=True)   # DEBUG ONLY
    ## NOTE BELOW: Changed from dynr.time to reflected.time due to length mismatch between dataframes.
    reflected['disp_corr_strain'] = disp_inc.disp_shift_wave(t = (reflected.time/1e6).tolist(), y = reflected.strain.tolist(), dz= -1*(igauge2sample/1000), backshift=True)  # dz is a shift factor. From Bragov2019, this looks like a DISTANCE. Here we assume the gage is halfway down the bar and this is in meters. Backshift just does a simple shift to get it back to the right time after dispersion correction. - Ben
    #s_refl = disp_inc.disp_shift_wave(t = (reflected.time/1e6).tolist(), y = reflected.strain.tolist(), dz= -1*(igauge2sample/1000), backshift=True)  # dz is a shift factor. From Bragov2019, this looks like a DISTANCE. Here we assume the gage is halfway down the bar and this is in meters. Backshift just does a simple shift to get it back to the right time after dispersion correction. - Ben
    #s_refl_back = disp_inc.disp_shift_wave(t = (reflected.time/1e6).tolist(), y = s_refl, dz= (igauge2sample/1000), backshift=True)
    #s3 = b1.simple_shift_wave(t = (incident.time/1e6).tolist(), y = incident.strain.tolist(), tshift = (igauge2sample/1000)/b1.cb)   #, 2/b1.cb is timeshift; I think this literally just shifts everything by an offset without doing anything else.
    plt.plot(reflected.time, reflected.strain, 'b', label='original data')
    plt.plot(reflected.time, reflected.disp_corr_strain, 'g', label='wave shift')
    #plt.plot(reflected.time, s_refl_back, 'r', label='reverse to original')
    #plt.plot(incident.time, s3, 'r', label='simple shift')
    plt.xlabel('Time (\u03BCs)')    # time in microseconds
    plt.ylabel('strain')
    plt.title('Reflected Pulse Dispersion Correction')
    plt.legend()
    #plt.show(block=True)   # DEBUG ONLY
    plt.close()

    ## For the transmitted signal, reverse dispersion from sample to transmitted gauge.
    ## NOTE BELOW: As with reflected, changed dynt.time to transmitted.time so lengths matched
    transmitted['disp_corr_strain'] = disp_trans.disp_shift_wave(t = (transmitted.time/1e6).tolist(), y = transmitted.strain.tolist(), dz= -1*(tgauge2sample/1000), backshift=True)
    #t_refl = disp_trans.disp_shift_wave(t = (transmitted.time/1e6).tolist(), y = transmitted.strain.tolist(), dz= -1*(tgauge2sample/1000), backshift=True)
    #t_refl_back = disp_trans.disp_shift_wave(t = (transmitted.time/1e6).tolist(), y = t_refl, dz= (tgauge2sample/1000), backshift=True) # dispersion correct back to original to check they're the same.
    plt.plot(transmitted.time, transmitted.strain, 'b*', label='original data')
    plt.plot(transmitted.time, transmitted.disp_corr_strain, 'g', label='wave shift')
    #plt.plot(transmitted.time, t_refl_back, 'r', label='reverse dispersion to original')
    #plt.plot(incident.time, s3, 'r', label='simple shift')
    plt.xlabel('Time (\u03BCs)')    # time in microseconds
    plt.ylabel('Strain')
    plt.title('Transmitted Pulse Dispersion Correction')
    plt.legend()
    #plt.show(block=True)   # DEBUG ONLY
    plt.close()

else: # No Dispersion Correction! Just pipes everything through to a new column without changes.
    incident['disp_corr_strain'] = incident['strain']
    reflected['disp_corr_strain'] = reflected['strain']
    transmitted['disp_corr_strain'] = transmitted['strain']

#%%######################################
### FINE-TUNE ALIGN TRANSMITTED & REFLECTED PULSES
### This plots the transmitted data atop the incident, allows the user to align the rising portions of each, then resets the time portion of the transmitted pulse
#########################################
var_i = pd.DataFrame({'time': incident.time, 'strain':incident.disp_corr_strain})
var_r = pd.DataFrame({'time': reflected.time, 'strain':reflected.disp_corr_strain})
var_t = pd.DataFrame({'time': transmitted.time, 'strain':transmitted.disp_corr_strain})
var = pd.DataFrame({'time': incident.time, 'variance':(var_i.strain - (var_r.strain + var_t.strain))})

fig = plt.figure(figsize=(12, 6))
fig.subplots_adjust(bottom=0.35)

plt.title('Pulse Alignment')
plt.xlabel('Time (microsecond)')
plt.ylabel('Strain (mm/mm)')

ttshift = 0.00  # Initial default for time shift on the slider (below).
trshift = 0.00
multishift = 2.00   # start the multiplier on 2, and it also helps to be able to see the slider
transmitted.time += ttshift  # then shift the data with this timeshift
reflected.time += trshift
deltatt = 0.00
deltatr = 0.00

#ax.plot(incident.time, incident.disp_corr_strain, 'r-', label='Incident')
ax = fig.add_subplot(121, facecolor='#FFFFCC')
k, = ax.plot(incident.time, incident.disp_corr_strain, 'r-', label='Incident')
l, = ax.plot(transmitted.time, transmitted.disp_corr_strain, 'g-', label='Transmitted')
m, = ax.plot(reflected.time, reflected.disp_corr_strain, 'b-', label='Reflected')
# plt.axis([0, 1, -10, 10])

ax2 = fig.add_subplot(122, facecolor='#FFFFCC')
ax2.set_title('Variance')
line2, = ax2.plot(var.time, var.variance, 'k-')

ax.legend()
ax.grid()
ax2.grid()
axcolor = 'lightgoldenrodyellow'
ax_tt = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
ax_tr = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
ax_multi = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)
slide_trans = Slider(ax_tt, 'Transmitted Offset', -30.0, 30.0, valinit=ttshift, color='green')  # changed min value from -30 to 0, because shifting time was way easier without resampling
slide_refl = Slider(ax_tr, 'Reflected Offset', -30.0, 30.0, valinit=trshift, color='blue')
slide_multi = Slider(ax_multi, 'Transmitted Multiplier', 1, 25, valinit=multishift, color='red')  # multiplier to make it easier to align transmitted wave. Does not get applied to data, only to the plot.

#BUG (FIXED): The time shift was previously set to SUBTRACT delta time instead of ADD. I think this was backwards. -BMM 20210330
# def update(val, deltatt, deltatr):
def update(val):
    l.set_xdata(transmitted.time + slide_trans.val)
    l.set_ydata(transmitted.disp_corr_strain * slide_multi.val)
    m.set_xdata(reflected.time + slide_refl.val)
    var_i = pd.DataFrame({'time': incident.time, 'strain':incident.disp_corr_strain})
    var_r = pd.DataFrame({'time': reflected.time+slide_refl.val, 'strain':reflected.disp_corr_strain})
    var_t = pd.DataFrame({'time': transmitted.time+slide_trans.val, 'strain':transmitted.disp_corr_strain})
    # Here's where things get clunky, but let's trim to just the time values that overlap, then we can calculate variance.
    var_min_time = max(min(var_i.time), min(var_r.time), min(var_t.time))
    var_max_time = min(max(var_i.time), max(var_r.time), max(var_t.time))
    var_i = var_i[var_i.time > var_min_time]
    var_i = var_i[var_i.time < var_max_time]
    var_r = var_r[var_r.time > var_min_time]
    var_r = var_r[var_r.time < var_max_time]
    var_t = var_t[var_t.time > var_min_time]
    var_t = var_t[var_t.time < var_max_time]
    x1 = var_i.to_numpy()
    x2 = var_r.to_numpy()
    x3 = var_t.to_numpy()

    try:    #BUG: it's possible for the two arrays below (new_time and new_var) to be different lengths, which throws a non-fatal (but annoyingly verbose) error
        new_time = (x1[:,0]+x2[:,0]+x3[:,0])/3  # new time is the average of all three
        #print('Min/max new time: {}, {}'.format(min(new_time), max(new_time)))
        new_var = x1[:,1]-(x2[:,1]+x3[:,1])
        #print('Mean variance: {}'.format(np.mean(new_var)))
    except:
        return

    line2.set_xdata(new_time)
    line2.set_ydata(new_var)
    fig.canvas.draw_idle()
#
slide_trans.on_changed(update)
slide_refl.on_changed(update)
slide_multi.on_changed(update)
#
resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', hovercolor='0.975')
#
def reset(event):
    slide_trans.reset()
    slide_refl.reset()
    slide_multi.reset()
#
button.on_clicked(reset)
#
plt.show(block=True)

print(f'Transmitted Pulse Timeshift (\u03BCs): {slide_trans.val:0.2f}')
print(f'Reflected Pulse Timeshift (\u03BCs): {slide_refl.val:0.2f}')

# Now that we know timeshift values, let's resample the transmitted and reflected data to capture all the important bits.
transmitted.time = transmitted.time + slide_trans.val
reflected.time = reflected.time + slide_refl.val

# Roughly calculate sound speed in bars based on time between arrival of pulses and distance from gauge to sample and back.
c_bi_xp = (2 * igauge2sample) / ((dynr.time[0] + slide_refl.val) - dyni.time[0])
print(f'Incident Bar Sound Speed (input) (m/s): {c_bi:0.1f}')
print(f'Incident Bar Sound Speed (exp. calc.) (m/s): {unitconversion.mm_us2m_s(c_bi_xp):0.1f}')
if abs((unitconversion.mm_us2m_s(c_bi_xp)-c_bi)/c_bi) > 0.01:
    print('WARNING: Sound speed mismatch (input vs. experimental) greater than 1%. Please check data.')

transmitted = transmitted[transmitted.time >= 0].reset_index(drop=True)
reflected = reflected[reflected.time >= 0].reset_index(drop=True)
incident = incident.reset_index(drop=True)

### Output shifted strain time data.
if suppress_output == False:
    output_plots.plot_disp_corr_strain_time(incident, reflected, transmitted, save_dir)

#%%######################################
### Set cutoff value so we only process the important stuff
#########################################
# window = 1001   # one more than an even number to account for the wonky way for loops are handled later (increment from index 1 instead of 0 so we can subtract the previous values)
window = min(len(incident.time), len(reflected.time), len(transmitted.time))  # We'll try this (set to the minumum array length) instead of an arbitrary constant

incident = incident[0:window]
reflected = reflected[0:window]
transmitted = transmitted[0:window]
print(f'window: {window}')

# Plot and export the variance, since the incident pulse should have the same magnitude as the transmitted and reflected combined.
variance = pd.DataFrame({'deltastrain':(incident.strain - (reflected.strain + transmitted.strain)), 'time':transmitted.time})

if suppress_output == False:
    output_plots.plot_variance(variance, save_dir)
    variance.to_csv(save_dir + testname + '_variance.txt', sep='\t', index=False)   # export variance

#%%########################################################
###########################################################
#####     ONE-WAVE ANALYSIS
###########################################################
###########################################################

# x = np.array([])           # initialize numpy array, put x = np.array([0]) if you want a zero as the first element
engstrain_1wave = np.array([0])  # initialize numpy array for engineering strain

for i in range(1, len(reflected)):  # shifting the indices by one because MATLAB is stupid.
    engstrain_1wave = np.append(engstrain_1wave, (engstrain_1wave[i - 1] + ((2.0 * unitconversion.m_s2mm_us(c_bi) / thick) * (reflected.disp_corr_strain[i] * (reflected.time[i] - reflected.time[i - 1])))))  # Should be integral over time of Gray Eq 12; units(mm/mm)
    #engstrain_1wave = np.append(engstrain_1wave, (engstrain_1wave[i - 1] + 2.0 * unitconversion.m_s2mm_us(c_bi) / thick * reflected.strain[i] * (reflected.time[i] - reflected.time[i - 1])))  # Should be integral over time of Gray Eq 12; units(mm/mm)

#import scipy.integrate as integrate
#engstrain_1wave_new = integrate.cumtrapz(reflected.disp_corr_strain, reflected.time, initial = 0)
#engstrain_1wave_new = engstrain_1wave_new*3

onewave = pd.DataFrame({'time':reflected.time-reflected.time[0],
    'eng_stress':(A_ratio_t * (E_bt*1000) * transmitted.disp_corr_strain),   # This is basically Eq 13 from Gray2000; units(MPa)
    'eng_rate': ((2 * unitconversion.m_s2mm_us(c_bi) * reflected.disp_corr_strain) / thick) * 1e6,  # This is basically Eq 12 from Gray2000; units(/s)
#    'eng_stress':(A_ratio_t * (E_bt*1000) * transmitted.strain),   # This is basically Eq 13 from Gray2000; units(MPa)
#    'eng_rate': ((2 * unitconversion.m_s2mm_us(c_bi) * reflected.strain) / thick) * 1e6,  # This is basically Eq 12 from Gray2000; units(/s)
#    'eng_strain': engstrain_1wave   # Should be integral over time of Gray Eq 12; units(mm/mm)
    'eng_strain': engstrain_1wave
    })

fromrate = onewave.eng_rate * 1e-6 * 0.1
engstrain_1wave_fromrate = np.cumsum(fromrate)

# plots comparison of integration methods. Even if not suppressing output, this is of limited use.
if suppress_output == False:
    output_plots.plot_1W_integral_compare(reflected, engstrain_1wave, engstrain_1wave_fromrate, save_dir)

onewave['true_stress'] = onewave.eng_stress * (1. - onewave.eng_strain)
onewave['true_strain'] = -1 * np.log(1 - onewave.eng_strain)

repsj1 = np.array([0])  # initialize numpy array for true strain rate

for i in range(1, len(onewave)):  # shifting the indices by one because MATLAB is stupid.
    repsj1 = np.append(repsj1, (onewave.true_strain[i] - onewave.true_strain[i - 1]) / (onewave.time[i] - onewave.time[i - 1]) * 1e6)

onewave['true_rate'] = repsj1

if suppress_output == False:
    output_plots.plot_1W(onewave, testname, save_dir)   # this is actually a bunch of plots
    onewave.to_csv(save_dir + testname + '_1wave.txt', sep='\t', index=False)

# truerate2 = (engrate_1wave/(1-engstrain_1wave))*1e12   # from Bill Blumenthal, gives same result as other method
# print(truerate2)

#%%########################################################
###########################################################
#####     TWO-WAVE ANALYSIS
###########################################################
###########################################################

engstrain_2wave = np.array([0])  # initialize numpy array for engineering strain

for i in range(1, len(transmitted)):  # shifting the indices by one because MATLAB is stupid.
    engstrain_2wave = np.append(engstrain_2wave, (engstrain_2wave[i-1] - ((2.0 / thick) * (
            unitconversion.m_s2mm_us(c_bt) * (transmitted.disp_corr_strain[i] * (transmitted.time[i] - transmitted.time[i-1])) - (unitconversion.m_s2mm_us(c_bi) * incident.disp_corr_strain[i] * (incident.time[i] - incident.time[i-1]))))))  # From Blumenthal's notes
    #engstrain_2wave = np.append(engstrain_2wave, (engstrain_2wave[i-1] - (2.0 / thick) * (
    #        unitconversion.m_s2mm_us(c_bt) * (transmitted.strain[i] * (transmitted.time[i] - transmitted.time[i-1])) - (unitconversion.m_s2mm_us(c_bi) * incident.strain[i] * (incident.time[i] - incident.time[i-1])))))  # From Blumenthal's notes


twowave = pd.DataFrame({'time':reflected.time-reflected.time[0],
    'eng_stress':(A_ratio_i * (E_bi*1000) * (incident.disp_corr_strain - reflected.disp_corr_strain)),   # This is basically Eq 13 from Gray2000; units(MPa)
    'eng_rate': ((2 * ((unitconversion.m_s2mm_us(c_bt) * transmitted.disp_corr_strain) - (unitconversion.m_s2mm_us(c_bi) * incident.disp_corr_strain)) * 1e6) / thick),  # This is basically Eq 12 from Gray2000; units(/s) # Derived from Chris Meredith's telecon 201712 slide 6? New notes from Blumenthal
#    'eng_stress':(A_ratio_i * (E_bi*1000) * (incident.strain - reflected.strain)),   # This is basically Eq 13 from Gray2000; units(MPa)
#    'eng_rate': ((2 * ((unitconversion.m_s2mm_us(c_bt) * transmitted.strain) - (unitconversion.m_s2mm_us(c_bi) * incident.strain)) * 1e6) / thick),  # This is basically Eq 12 from Gray2000; units(/s) # Derived from Chris Meredith's telecon 201712 slide 6? New notes from Blumenthal
    'eng_strain': engstrain_2wave,  # Should be integral over time of Gray Eq 12; units(mm/mm)
    })

twowave['true_stress'] = twowave.eng_stress * (1. - twowave.eng_strain)
twowave['true_strain'] = -1 * np.log(1 - twowave.eng_strain)

repsj2 = np.array([0])  # initialize numpy array for true strain rate
for i in range(1, len(twowave)):  # shifting the indices by one because MATLAB is stupid.
    repsj2 = np.append(repsj2, (twowave.true_strain[i] - twowave.true_strain[i - 1]) / (twowave.time[i] - twowave.time[i - 1]) * 1e6)

twowave['true_rate'] = repsj2

if suppress_output == False:
    # Compare forces at each bar-sample interface; Only needed for verbose output option, so we have it here so it's not caluclated unless needed.
    barforce = pd.DataFrame({'time':incident.time,
        'incident': (A_bar_i * (E_bi*1000) * (incident.strain - reflected.strain)),   # ch2; Based on Gray2000 Eq9, but the equation in the paper assumes reflected strain is negative. Units = MPa
        'transmitted': (A_bar_t * (E_bt*1000) * transmitted.strain),           # ch2; Based on Gray2000 Eq10; units: MPa
        })
    barforce['difference'] = barforce.incident-barforce.transmitted
#
    output_plots.plot_2W(twowave, testname, save_dir)   # several plots available in this fuction
    twowave.to_csv(save_dir + testname + '_2wave.txt', sep='\t', index=False)
    output_plots.plot_barforce(incident, transmitted, barforce, save_dir)
    barforce.to_csv(save_dir+'barforce.txt', index=None, sep='\t', lineterminator='\r\n')

#%%########################################################
###########################################################
#####     THREE-WAVE ANALYSIS
###########################################################
###########################################################
engstrain_3wave = np.array([0])
for i in range(1, len(transmitted)):
    #      engstress_3wave = np.append(engstress_3wave, ((A_bar * E_bi(MPa) * (si[i]-sr[i]))+(A_bar * E_bt(MPa) * st[i]))/(2*A_sample))     # Gray Eqs. 14 & 15
    engstrain_3wave = np.append(engstrain_3wave, engstrain_3wave[i - 1] + unitconversion.m_s2mm_us(c_bi) / thick * (incident.strain[i] + reflected.strain[i] - transmitted.strain[i]) * (reflected.time[i] - reflected.time[i - 1]))  # ch1

threewave = pd.DataFrame({
    'time': incident.time,
    'eng_stress': (0.5 * (((E_bi*1000) * (incident.disp_corr_strain - reflected.disp_corr_strain) * A_ratio_i) + ((E_bt*1000) * transmitted.disp_corr_strain * A_ratio_t))),    # Gray Eqs. 14 & 15; made the first thing si-sr to account for the opposite sign - BMM, units = MPa
    #'eng_stress': (0.5 * (((E_bi*1000) * (incident.strain - reflected.strain) * A_ratio_i) + ((E_bt*1000) * transmitted.strain * A_ratio_t))),    # Gray Eqs. 14 & 15; made the first thing si-sr to account for the opposite sign - BMM, units = MPa
    'eng_strain': engstrain_3wave
    })

threewave['true_stress'] = threewave.eng_stress * (1. - threewave.eng_strain)
threewave['true_strain'] = -np.log(1 - threewave.eng_strain)

repsj3 = np.array([0])
for i in range(1, len(threewave)):
    repsj3 = np.append(repsj3, ((threewave.true_strain[i] - threewave.true_strain[i - 1]) / (threewave.time[i] - threewave.time[i - 1]) * 1e6))

threewave['true_rate'] = repsj3

if suppress_output == False:
    output_plots.plot_3W(threewave, testname, save_dir)
    threewave.to_csv(save_dir + testname + '_3wave.txt', sep='\t', index=False)

fig, ax = plt.subplots()
if eng_out == False:
    ax.plot(onewave.true_strain, onewave.true_stress, 'b-', label='1-Wave')
    ax.plot(twowave.true_strain, twowave.true_stress, 'g-', label='2-Wave')
    ax.plot(threewave.true_strain, threewave.true_stress, 'r-', label='3-wave')
    plt.xlabel('True Strain (mm/mm)')
    plt.ylabel('True Stress (MPa)')
else:
    ax.plot(onewave.eng_strain, onewave.eng_stress, 'b-', label='1-Wave')
    ax.plot(twowave.eng_strain, twowave.eng_stress, 'g-', label='2-Wave')
    ax.plot(threewave.eng_strain, threewave.eng_stress, 'r-', label='3-wave')
    plt.xlabel('Eng. Strain (mm/mm)')
    plt.ylabel('Eng. Stress (MPa)')

plt.title(testname + ' (1&2&3-Wave Analysis Comparison)')
plt.legend()
plt.grid()
if suppress_output == False:
    plt.savefig(save_dir + '1and2and3-Wave_Compare_True_Stress-Strain.png')

plt.show()
plt.close()

print('Data analysis complete.')
#%%########################################################
###########################################################
#####     EXPORT DATA AS EXCEL
###########################################################
###########################################################
data_crop_threshold = 0.6   # the DA program crops the unloading portion at some threshold, reproduced here
# note: we should still use the true values to find cutoffs, otherwise necking may cause weirdness with engineering stress/strain
xmax = onewave.loc[onewave.true_stress > (max(onewave.true_stress)*data_crop_threshold)].index[-1]

header1 = 'Time(s), 1W-True_Strain, True_Stress_MPa, Rate_1/s, 2W-True_Strain, True_Stress_MPa, Rate_1/s'
#ugh. the text importer is setup with regex, and someone named multiple columns the same so we have to manually write the column headers, then suppress them below (because they're wrong in the pandas so it doesn't overwrite the same column)
granta_1W_True = pd.DataFrame(
    {'Time(s)': onewave.time*1e-6, '1W-True_Strain': onewave.true_strain, 'True_Stress_MPa': onewave.true_stress, 'Rate_1/s': onewave.true_rate})
granta_2W_True = pd.DataFrame(
    {'2W-True_Strain': twowave.true_strain, '2W-True_Stress_MPa':twowave.true_stress, '2W-Rate_1/s':twowave.true_rate})
granta_1W_True = granta_1W_True.head(xmax)
granta_2W_True = granta_2W_True.head(xmax)

# only bother with this if we need it.
if eng_out == True:
    header1 = 'Time(s), 1W-Eng_Strain, Eng_Stress_MPa, Rate_1/s, 2W-Eng_Strain, Eng_Stress_MPa, Rate_1/s'   # re-write this header to reflect engineering values if needed.
    granta_1W_Eng = pd.DataFrame(
        {'Time(s)': onewave.time*1e-6, '1W-Eng_Strain': onewave.eng_strain, 'Eng_Stress_MPa': onewave.eng_stress, 'Rate_1/s': onewave.eng_rate})
    granta_2W_Eng = pd.DataFrame(
        {'2W-Eng_Strain': twowave.eng_strain, '2W-Eng_Stress_MPa':twowave.eng_stress, '2W-Rate_1/s':twowave.eng_rate})
    granta_1W_Eng = granta_1W_Eng.head(xmax)
    granta_2W_Eng = granta_2W_Eng.head(xmax)

##############################################
### Export to a dat file, for better exchange with Sierra Peaks Kolsky Program
##############################################
if eng_out == False:
    strain_stream = granta_1W_True['1W-True_Strain']
    rate_stream = granta_1W_True['Rate_1/s']
    label_stream = 'True Strain Rate'
else:
    strain_stream = granta_1W_Eng['1W-Eng_Strain']
    rate_stream = granta_1W_Eng['Rate_1/s']
    label_stream = 'Eng. Strain Rate'

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(121, facecolor='#FFFFCC')  # 1 row, 2 columns of subfigures, subplot 1
plt.title('Avg. Strain Rate Selector')
plt.xlabel('Strain (mm/mm)')
plt.ylabel('Strain Rate (/s)')
#
ax.plot(strain_stream, rate_stream, 'g-', label=label_stream)
ax.set_title('Click and Drag to Select Strain Rate Range')
plt.legend()
plt.grid()
#
ax2 = fig.add_subplot(122, facecolor='#FFFFCC')
ax2.set_title('Zoomed View')
line2, = ax2.plot(strain_stream, rate_stream, 'g-')
#
class Rate_Select():
    """Container for strain rate selection"""
    def __init__(self):
        self.coords = {}
    def __call__(self, xmin, xmax):
        indmin, indmax = np.searchsorted(strain_stream, (xmin, xmax))
        indmax = min(len(strain_stream)-1, indmax)
        thisx = strain_stream[indmin:indmax]
        thisy = rate_stream[indmin:indmax]
        self.coords['strain'] = thisx
        self.coords['rate'] = thisy
        self.smax = strain_stream[indmax]
        self.smin = strain_stream[indmin]
        self.avg1 = np.mean(thisy)
        line2.set_data(thisx, thisy)
        ax2.set_xlim(thisx.min(), thisx.max())
        ax2.set_ylim(thisy.min(), thisy.max())
        fig.canvas.draw()
#
rateselect = Rate_Select()
#
try:
    span = SpanSelector(ax, rateselect, 'horizontal', useblit=True, props=dict(alpha=0.5, facecolor='red') )
except:
    span = SpanSelector(ax, rateselect, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='red') )
    print('Time to update matplotlib.')
plt.show(block=True)

if eng_out == False:
    try:
        final_strain = int(granta_1W_True['1W-True_Strain'].iloc[-1]*100)
    except:
        final_strain = 0
else:
    final_strain = int(granta_1W_Eng['1W-Eng_Strain'].iloc[-1]*100)
print(f'Final Strain = {final_strain}%')

try:
    avg_rate = int(round(rateselect.avg1, -2))  # should round to the nearest 100 /s
    #range_rate = [int(round(rateselect.smin * 100)), int(round(rateselect.smax * 100))]    # rounds to integer of strain
    range_rate = [round(rateselect.smin * 100, 1), round(rateselect.smax * 100, 1)]         # rounds to 1 decimal of % strain
    print(f'Strain Rate = {avg_rate} /s @ {range_rate[0]}-{range_rate[1]}% strain')
    sample_cond = f'Analyzed {date}; timing {dyni_start:.1f}/{dyni_end:.1f}/{dynt_start+slide_trans.val:.1f}; Strain Rate = {avg_rate} /s @ {range_rate[0]}-{range_rate[1]}% strain; Final strain = {final_strain}%'
except:
    avg_rate = 'XXXX'
    range_rate = ['err', 'ERR']

# Here we will make a dataframe that can be exported with the rest of the output below.
dfin = pd.DataFrame([
    ['Test Date', test_date],
    ['Test Time', test_time],
    ['System ID', sysid],
    ['Test Loading', test_loading],
    ['Sample Description', Sample_Description],
    ['Initial Sample Diameter(mm)', d],
    ['Initial Sample Length(mm)', thick],
    ['Test Temperature(C)', temperature],
    ['Breech Pressure(psi)', breech_pressure],
    ['Striker Length(mm)', len_striker],
    ['Tip Material', tip_matl],
    ['Sample Lubrication', sample_lube],
    ['Sample Material', sample_matl],
    ['Incident Calibration ID', calibration_id_i],
    ['Transmitted Calibration ID', calibration_id_t],
    ['Raw Data File Name', file_name],
    ['Incident/Reflected Bar Gain', int(gain_i)],
    ['Transmitted Bar Gain', int(gain_t)],
    ['Desired Striker Velocity(m/s)', v_striker_desired],
    ['Gap(mm)', gap],
    ['', '']
])
# mu μ = \u03BC in unicode, but we're putting the symbol here because we're outputting in western encoding
df_bi = pd.DataFrame([  # Create Data Frame for incident bar parameters (mostly calibration data)
    ['Incident Bar Properties', ''],
    ['Incident Bar ID', bar_i_id],
    ['Incident Bar Description', bar_i_description],
    ['Incident Bar Calibrated Gage Factor', gauge_factor_i],
    ['Incident Bar Reflected Time(μs)', bar_i_refl_time],
    ['Incident Bar Reflected Time Sound Speed(m/s)', c_bar_i_refl_time],
    ['Incident Bar Calibration Date', bar_i_calibration_date],
    ['Incident Bar Gage Description', bar_i_gage_desc],
    ['Incident Bar Youngs Modulus(GPa)', E_bi],
    ['Incident Bar Poisson Ratio', poisson_bi],
    ['Incident Bar Rod Sound Speed(m/s)', c_bi],
    ['Incident Bar Longitudinal Sound Speed(m/s)', c_bar_i_long],
    ['Incident Bar Shear Sound Speed(m/s)', c_bar_i_shear],
    ['Incident Bar Density(g/cm^3)', rho_bi/1000],
    ['Incident Bar Diameter(mm)', d_bar_i_equiv],
    ['Incident Bar Strain Gage Location(mm)', igauge2sample],
    ['Incident Bar Excitation Voltage(V)', excitation_i],
    ['Incident Bar Calibrated Gain', int(calibrated_gain_i)],
    ['', '']
])
df_bt = pd.DataFrame([  # Create Data Frame for transmitted bar parameters (mostly calibration data)
    ['Transmitted Bar Properties', ''],
    ['Transmitted Bar ID', bar_t_id],
    ['Transmitted Bar Description', bar_t_description],
    ['Transmitted Bar Calibrated Gage Factor', gauge_factor_t],
    ['Transmitted Calibration Date', bar_t_calibration_date],
    ['Transmitted Bar Gage Description', bar_t_gage_desc],
    ['Transmitted Bar Youngs Modulus(GPa)', E_bt],
    ['Transmitted Bar Poisson Ratio', poisson_bt],
    ['Transmitted Bar Rod Sound Speed(m/s)', c_bt],
    ['Transmitted Bar Longitudinal Sound Speed(m/s)', c_bar_t_long],
    ['Transmitted Bar Shear Sound Speed(m/s)', c_bar_t_shear],
    ['Transmitted Bar Density(g/cm^3)', rho_bt/1000],
    ['Transmitted Bar Diameter(mm)', d_bar_t_equiv],
    ['Transmitted Bar Strain Gage Location(mm)', tgauge2sample],
    ['Transmitted Bar Excitation Voltage(V)', excitation_t],
    ['Transmitted Bar Calibrated Gain', int(calibrated_gain_t)],
    ['', '']
])
df_finalvalues = pd.DataFrame([  # Create Data Frame for final sample values (dimensions, striker velocity, etc.))
    ['Final Values', ''],
    ['Incident Wave Initial Time(μs)', dyni_start],
    ['Incident Wave Final Time(μs)', dyni_end],
    ['Reflected Wave Initial Time(μs)', dynr_start + slide_refl.val],
    ['Transmitted Wave Initial Time(μs)', dynt_start+slide_trans.val],
    ['Incident/Reflected Wave Baseline Voltage(V)', (avg_ch1 + avg_ch2) / 2000],  # pulled from baseline section above, and /1000 to convert back to V
    ['Incident/Reflected Baseline V Initial Time(μs)', basemin],
    ['Incident/Reflected Baseline V Final Time(μs)', basemax],
    ['Transmitted Wave Baseline Voltage(V)', (avg_ch3 + avg_ch4) / 2000],  # pulled from baseline section above, and /1000 to convert back to V
    ['Transmitted Baseline V Initial Time(μs)', basemin],
    ['Transmitted Baseline V Final Time(μs)', basemax],
    ['Final Sample Diameter(mm)', ''],
    ['Final Sample Length(mm)', ''],
    ['Final Sample Condition', sample_cond],
    ['Actual Striker Velocity(m/s)', v_striker],
    ['Desired Striker Velocity(m/s)', v_striker_desired]
])
temp = dfin.iloc[0]
#df_input = dfin.drop(0).append([df_bi, df_bt, df_finalvalues], ignore_index=True) #CLEANUP
# BUG: The above line works, but kicks a deprecation warning in pandas, so the following line was subbed.
df_input = pd.concat([dfin.drop(0), df_bi, df_bt, df_finalvalues], axis=0, ignore_index = True)
df_input.columns = [temp[0],temp[1]]

df1 = pd.DataFrame({'Time (microsecond)': onewave.time, 'True Strain (mm/mm) (1-wave)': onewave.true_strain,
                'True Strain Rate (/s) (1-wave)': onewave.true_rate, 'True Stress (MPa) (1-wave)': onewave.true_stress})
df2 = pd.DataFrame({'Time (microsecond)': twowave.time, 'True Strain (mm/mm) (2-wave)': twowave.true_strain,
                'True Strain Rate (/s) (2-wave)': twowave.true_rate, 'True Stress (MPa) (2-wave)': twowave.true_stress})
df3 = pd.DataFrame({'Time (microsecond)': threewave.time, 'True Strain (mm/mm) (3-wave)': threewave.true_strain,
                'True Strain Rate (/s) (3-wave)': threewave.true_rate, 'True Stress (MPa) (3-wave)': threewave.true_stress})
df4 = pd.DataFrame({'Time (microsecond)': onewave.time, 'Eng. Strain (mm/mm) (1-wave)': onewave.eng_strain,
                 'Eng. Strain Rate (/s) (1-wave)': onewave.eng_rate, 'Eng. Stress (MPa) (1-wave)': onewave.eng_stress})
df5 = pd.DataFrame({'Time (microsecond)': twowave.time, 'Eng. Strain (mm/mm) (2-wave)': twowave.eng_strain,
                 'Eng. Strain Rate (/s) (2-wave)': twowave.eng_rate, 'Eng. Stress (MPa) (2-wave)': twowave.eng_stress, })
df6 = pd.DataFrame({'Time (microsecond)': threewave.time, 'Eng. Strain (mm/mm) (3-wave)': threewave.eng_strain,
                 'Eng. Stress (MPa) (3-wave)': threewave.eng_stress})  # , 'True Strain Rate (/s) (3-wave)':repsj3 })

if suppress_output == False:
    output_excel.excel_final(df_input, df1, df2, df3, df4, df5, df6, testname, save_dir)

# Reconstruct .raw file if input is mini-bar, especially for Granta upload ease-of-use
if platform != 'LANL':
    raw_file = data_dir + file_name[:-8] + '.raw'
    raw_head = df_input.head(56)  # everything but the final values
    # there are a couple lines in the raw files that aren't actual data that need to be tacked on to the top, so we create some with neg. index, then sort to put them at the top
    raw_wfm.loc[-2] = [-5.00E-05, -5.00E-05]
    raw_wfm.loc[-1] = [timestep, timestep]
    raw_wfm.index = raw_wfm.index + 2  # shifting index
    raw_wfm = raw_wfm.sort_index().reset_index(drop=True)  # sorting by index
    raw_out = pd.concat([raw_head, raw_wfm], axis=1)
    #granta_wfm.to_csv(raw_file, index=None, sep='\t', lineterminator='\r\n', mode = 'a')
    raw_out.to_csv(raw_file, index = None, sep = '\t', lineterminator='\r\n')

# Write the output file in LANL format (.dat)
if eng_out == False:
    dat_df = pd.concat([df_input, df1.drop(columns='Time (microsecond)'), df2.drop(columns='Time (microsecond)')], axis=1)
else:
    dat_df = pd.concat([df_input, df4.drop(columns='Time (microsecond)'), df5.drop(columns='Time (microsecond)')], axis=1)

dat_df.to_csv(data_dir + file_name[:-4] + '_' + str(avg_rate) +'s.dat', index=None, sep='\t', lineterminator='\r\n')#, encoding = 'cp1252')

#%%########################################################
###########################################################
#####     EXPORT DATA AS GRANTA FORMAT
###########################################################
###########################################################
encode_out = 'utf-8'    # 'cp1252' for western windows/ANSI
if suppress_granta == False:
    ### Clean file first... _clean.txt: Expected output, 2-column csv, headers "True_Strain" and "True_Stress_MPa"
    clean_file = data_dir + file_name[:-4] + '_' + str(avg_rate) +'s_clean.txt'
    if eng_out == False:
        clean_data = pd.DataFrame({'True_Strain':granta_1W_True['1W-True_Strain'], 'True_Stress_MPa':granta_1W_True['True_Stress_MPa']})
    else:
        clean_data = pd.DataFrame({'Eng_Strain':granta_1W_Eng['1W-Eng_Strain'], 'Eng_Stress_MPa':granta_1W_Eng['Eng_Stress_MPa']})
    clean_data.to_csv(clean_file, sep=',', encoding = encode_out, lineterminator='\r\n', index = False, header = True)    # Must specify lineterminator or it crashes DA on import. Might need to specify encoding = 'cp1252', but it's okay for me without.
## Collect GRANTA input from GUI
    from GRANTA import granta_input_ui
    app = granta_input_ui.QApplication(sys.argv)
    main_window = granta_input_ui.InputGrantaData()
    main_window.show()
    app.exec()
#
    dataSensitivity, dataOriginator, testContact, testOperator, testLab, testGroup, testLoc, specimenOrient = main_window.get_output()
    # dataSensitivity = str(granta_input_ui.dataSensitivity.currentText())
    # dataOriginator = str(granta_input_ui.dataOriginator.currentText())
    # testContact = str(granta_input_ui.testContact.currentText())
    # testOperator = str(granta_input_ui.testOperator.currentText())
    # testLab = str(granta_input_ui.testLab.currentText())
    # testGroup = str(granta_input_ui.testGroup.currentText())
    # testLoc = str(granta_input_ui.testLoc.currentText())
    # specimenOrient = str(granta_input_ui.specimenOrient.currentText())
    # also possible: granta_input_ui.dataSensitivity.itemText(0) where zero is the index of the combobox option, granta_input_ui.testContact.currentIndex()
#
    granta_file = data_dir + file_name[:-4] + '_' + str(avg_rate) +'s_Granta.txt'
    granta_in = pd.read_csv(granta_dir + 'granta_defaults.txt', delimiter = ',', names = ['Field','Value','Units'])  #encoding = 'cp1252',
    granta_in['Value'].loc[granta_in['Field']=='Specimen ID']= file_name[:-4] + '_' + str(avg_rate) +'s'
    granta_in['Value'].loc[granta_in['Field']=='Data Sensitivity']= dataSensitivity
    granta_in['Value'].loc[granta_in['Field']=='Data Originator']= dataOriginator
    granta_in['Value'].loc[granta_in['Field']=='Test Lab']= testLab
    granta_in['Value'].loc[granta_in['Field']=='Test Group']= testGroup
    granta_in['Value'].loc[granta_in['Field']=='Test Location']= testLoc
    granta_in['Value'].loc[granta_in['Field']=='Testing Contact']= testContact
    granta_in['Value'].loc[granta_in['Field']=='Operator']= testOperator
    granta_in['Value'].loc[granta_in['Field']=='Specimen Orientation']= specimenOrient
    granta_in['Value'].loc[granta_in['Field']=='Strain Rate']= str(avg_rate)
    granta_in['Value'].loc[granta_in['Field']=='Analysis Notes']= sample_cond
    separator = '*********************************************************' # a bunch of *'s which go between sections in the Granta file
#
    granta_in.to_csv(granta_file, header = None, index=None, sep=',', encoding = 'cp1252', lineterminator='\r\n')   # seems like this encoding needs to be 1252 to deal with degrees sign (like 45deg IP), otherwise it chokes on import via Granta Toolbox
#
    with open(granta_file,'a', newline='', encoding = encode_out) as f: f.write(separator+'\r\n')
#
    df_input.to_csv(granta_file, index=None, sep=',', encoding = encode_out, lineterminator='\r\n', mode = 'a')
    # Note: with old text importer we needed to save with a unique separator and replace with the DA ' , ' separator (appending line by line)
#
    with open(granta_file,'a', newline='', encoding = encode_out) as f: f.write(separator+'\r\n')
    # Data crop and 1W & 2W packaging happens above so we can use it in the strain rate determination to crop the end data a bit before measuring avg strain rate
    if eng_out == False:
        granta_1W2W = pd.concat([granta_1W_True,granta_2W_True], axis=1)
    else:
        granta_1W2W = pd.concat([granta_1W_Eng,granta_2W_Eng], axis=1)
#
    with open(granta_file,'a', newline='', encoding = encode_out) as f: f.write(header1+'\r\n')
#
    granta_1W2W.to_csv(granta_file, header = None, index=None, sep=',', encoding = encode_out, lineterminator='\r\n', mode = 'a')   # make sure to suppress the header because it's wrong and we already wrote the correct one.
#
    header2 = 'Wave, 1'
    with open(granta_file,'a', newline='', encoding = encode_out) as f: f.write(header2+'\r\n')
#
    if eng_out == False:
        granta_1W_True['Wave'] = 1  # fill in all 1's for the wave number because why not?
        granta_1W_True.drop(columns='Time(s)', inplace=True)
        granta_1W_True.to_csv(granta_file, header = None, index=None, sep=',', encoding = encode_out, lineterminator='\r\n', mode = 'a')   # make sure to suppress the header because it's wrong and we already wrote the correct one.
    else:
        granta_1W_Eng['Wave'] = 1  # fill in all 1's for the wave number because why not?
        granta_1W_Eng.drop(columns='Time(s)', inplace=True)
        granta_1W_Eng.to_csv(granta_file, header = None, index=None, sep=',', encoding = encode_out, lineterminator='\r\n', mode = 'a')   # make sure to suppress the header because it's wrong and we already wrote the correct one.
#
    header3 = 'Wave, 2'
    with open(granta_file,'a', newline='', encoding = encode_out) as f: f.write(header3+'\r\n')    #, encoding = 'cp1252'
#
    if eng_out == False:
        granta_2W_True['Wave'] = 2  # fill in all 2's for the wave number because why not?
        granta_2W_True.to_csv(granta_file, header = None, index=None, sep=',', encoding = encode_out, lineterminator='\r\n', mode = 'a')   # make sure to suppress the header because it's wrong and we already wrote the correct one.
    else:
        granta_2W_Eng['Wave'] = 2  # fill in all 2's for the wave number because why not?
        granta_2W_Eng.to_csv(granta_file, header = None, index=None, sep=',', encoding = encode_out, lineterminator='\r\n', mode = 'a')   # make sure to suppress the header because it's wrong and we already wrote the correct one.
    #header4 = '\n\rTimeRaw(s), Incident/Reflected Waveform, Transmitted Waveform'
    #with open(granta_file,'a', newline='') as f: f.write(header4+'\r\n')
    with open(granta_file,'a', newline='', encoding = encode_out) as f: f.write('\r\n')    # just a blank space this time, no separator
#
    granta_wfm.to_csv(granta_file, index=None, sep=',', encoding = encode_out, lineterminator='\r\n', mode = 'a')   # this was defined above (after initial data import) with correct headers, hopefully.

print('Final data export complete.')
