############################################
####      Functions to Edit Data        ####
############################################
# WARNING: RESETTING THE INDEX MIGHT CAUSE PROBLEMS LATER, ESPECIALLY IF (FOR EXAMPLE) THE TIME IS SET WITH INDEX*TIMESTEP.
# BE CAREFUL!
import pandas as pd # also required openpyxl to be installed for verbose features
#import matplotlib
#matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt
import numpy as np

def window_average(data_in, window_size):
    """Smooths data using rolling window."""
    data_out = data_in.rolling(window_size).mean().dropna().reset_index(drop=True)
    show_smooth(data_in, data_out)
    return data_out

def bin_data(data_in, bin_size):
    """Smooths data by binning."""
    data_out = data_in.groupby(np.arange(len(data_in.index))//bin_size).mean()  # This one automatically resets the index.
    show_smooth(data_in, data_out)
    return data_out

def use_every_nth(data_in, n):
    """Reduces data by discarding all but every n^th point."""
    data_out = data_in[::n].reset_index(drop=True)
    show_smooth(data_in, data_out)
    return data_out

def show_smooth(data_in, data_out):
    """Show a comparison plot to identify problems with smoothing."""
    fig, ax = plt.subplots()
    ax.plot(data_in, 'g-', label='Original')
    ax.plot(data_out, 'r-', label='Cropped')
    plt.title('Summary of smoothing results')
    plt.xlabel('x-axis')
    plt.ylabel('y-axis')
    plt.legend()
    plt.grid()
    plt.show()
    plt.close()