############################################
####      Functions to Create Plots     ####
############################################
import pandas as pd # also required openpyxl to be installed for verbose features
import matplotlib
matplotlib.use('QT5Agg')
import matplotlib.pyplot as plt
import numpy as np

def plot_pulses(dyni, dynr, dynt, save_dir):              # Convert from mm to in for some reason.
    """Plots the pulses together, at their original arrival times."""
    fig, ax = plt.subplots()
    ax.plot(dyni.time, dyni.v, 'r-', label='Incident')
    ax.plot(dynt.time, dynt.v, 'g-', label='Transmitted')
    ax.plot(dynr.time, dynr.v, 'b-', label='Reflected')
    plt.title('Original Arrival Times')
    plt.xlabel('Time (microsecond)')
    plt.ylabel('Voltage (mV)')
    plt.legend()
    plt.grid()
    plt.savefig(save_dir + 'CroppedData.png')
    #plt.show(block=True)         # DEBUG ONLY
    plt.close()
    return

def plot_pulses_DRtimeshift(dyni, dynr, dynt, tt_dr, tr_dr, save_dir):
    """Plots the timeshifted pulses where timeshift is calculated via dead-reckoning (sound speed & length)."""
    fig, ax = plt.subplots()
    ax.plot(tt_dr, dynt.v,'g-', label='Transmitted')
    ax.plot(dyni.time, dyni.v,'r-', label='Incident')
    ax.plot(tr_dr, dynr.v,'b-', label='Reflected')
    plt.title('Time-Shifted Signals')
    plt.xlabel('Time (microsecond)')
    plt.ylabel('Voltage (mV)')
    plt.legend()
    plt.grid()
    plt.savefig(save_dir + 'DR_Calculated_Voltage-Time_shift.png')
    #plt.show(block=True)   # DEBUG ONLY
    plt.close()
    return

def plot_disp_corr_strain_time(incident, reflected, transmitted, save_dir):
    """Plots dispersion corrected strain vs. time."""
    fig, ax = plt.subplots()
    ax.plot(transmitted.time, transmitted.disp_corr_strain, 'g-', label='Transmitted')
    ax.plot(incident.time, incident.disp_corr_strain, 'r-', label='Incident')
    ax.plot(reflected.time, reflected.disp_corr_strain, 'b-', label='Reflected')
    plt.title('Dispersion Corrected Gauge Strains')
    plt.xlabel('Time (microsecond)')
    plt.ylabel('Strain (mm/mm)')
    plt.legend()
    plt.grid()
    plt.savefig(save_dir + 'Strain-Time_shift.png')
    # plt.show()     # This is plotted again when we do the timeshift so this is redundant
    plt.close()
    return

def plot_variance(variance, save_dir):
    """Plots variance between incident pulse and transmitted plus reflected pulses."""
    fig, ax = plt.subplots()
    ax.plot(variance.time, variance.deltastrain, 'k-', label='Variance')
    plt.title('Variance Between Gauges')
    plt.xlabel('Time (microsecond)')
    plt.ylabel('Strain Variance Between Pulses (mm/mm)')
    plt.legend()
    plt.grid()
    plt.savefig(save_dir + 'Strain-Time_shift_var.png')
    # plt.show(block=True)     # This is plotted again when we do the timeshift so this is redundant
    plt.close()
    return

def plot_1W_integral_compare(reflected, engstrain_1wave, engstrain_1wave_fromrate, save_dir):
    """Plots comparison of 1-wave strain from various integration methods. Generally unnecessary to have this plot."""
    fig, ax = plt.subplots()
    ax.plot(reflected.time, engstrain_1wave, 'g-', label='Original')
    #ax.plot(reflected.time, engstrain_1wave_new, 'r-', label='New Integration')
    ax.plot(reflected.time, engstrain_1wave_fromrate, 'b--', label='Strain From Rate')
    plt.title('1-wave Integration Comparison')
    plt.xlabel('Time (microsecond)')
    plt.ylabel('Eng. Strain')
    plt.legend()
    plt.grid()
    plt.savefig(save_dir + '1-wave_Integral_Comparison.png')
    #plt.show(block=True)   # DEBUG ONLY
    plt.close()
    return

def plot_1W(onewave, testname, save_dir):
    """Creates various plots for 1-wave analysis."""
# eng stress-strain
    fig, ax = plt.subplots()
    ax.plot(onewave.eng_strain, onewave.eng_stress, label='engstrain_1wave')
    ax.plot(onewave.eng_strain[onewave.eng_stress.argmax(axis=0)], np.max(onewave.eng_stress), 'ro')
    ax.annotate('max=%.1f MPa, %.1f pct' % (np.max(onewave.eng_stress), onewave.eng_strain[onewave.eng_stress.argmax(axis=0)] * 100),
                xy=(onewave.eng_strain[onewave.eng_stress.argmax(axis=0)], np.max(onewave.eng_stress)), textcoords='data')
    plt.title(testname + ' (1-Wave Analysis)')
    plt.xlabel('Eng. Strain (mm/mm)')
    plt.ylabel('Eng. Stress (MPa)')
    plt.legend()
    plt.grid()
    plt.savefig(save_dir + '1-Wave_Eng_Stress-Strain.png')
    #plt.show(block=True)     # # DEBUG ONLY
    plt.close()
# true stress-strain
    fig, ax = plt.subplots()
    ax.plot(onewave.true_strain, onewave.true_stress)
    ax.plot(onewave.true_strain[onewave.true_stress.argmax(axis=0)], np.max(onewave.true_stress), 'ro')
    ax.annotate('max=%.1f MPa, %.1f pct' % (np.max(onewave.true_stress), onewave.true_strain[onewave.true_stress.argmax(axis=0)] * 100),
                xy=(onewave.true_strain[onewave.true_stress.argmax(axis=0)], np.max(onewave.true_stress)), textcoords='data')
    plt.title(testname + ' (1-Wave Analysis)')
    plt.xlabel('True Strain (mm/mm)')
    plt.ylabel('True Stress (MPa)')
    plt.grid()
    plt.savefig(save_dir + '1-Wave_True_Stress-Strain.png')
    #plt.show(block=True)   # DEBUG ONLY
    plt.close()
# true strain-strain rate
    fig, ax = plt.subplots()
    ax.plot(onewave.true_strain, onewave.true_rate)
    plt.title(testname)
    plt.xlabel('True Strain (mm/mm)')
    plt.ylabel('Strain Rate (/s)')
    plt.grid()
    plt.savefig(save_dir + '1-Wave_True_Strain-Rate.png')
    #plt.show() # DEBUG ONLY
    plt.close()  # cleanup
# True strain-time
    fig, ax = plt.subplots()
    ax.plot(onewave.time, onewave.true_strain)
    plt.title(testname)
    plt.xlabel('Time (microsecond)')
    plt.ylabel('True Strain (mm/mm)')
    plt.grid()
    plt.savefig(save_dir + '1-Wave_True_Strain-Time.png')
    # plt.show()
    plt.close()  # cleanup
    return

def plot_2W(twowave, testname, save_dir):
    """Creates various plots for 2-wave analysis."""
# Engineering stress-strain
    fig, ax = plt.subplots()
    # ax.plot(engstrain_1wave, engstress_1wave, label='engstrain_1wave')
    ax.plot(twowave.eng_strain, twowave.eng_stress, label='engstrain_2wave')
    ax.plot(twowave.eng_strain[twowave.eng_stress.argmax(axis=0)], np.max(twowave.eng_stress), 'ro')
    ax.annotate('max=%.1f MPa, %.1f pct' % (np.max(twowave.eng_stress), twowave.eng_strain[twowave.eng_stress.argmax(axis=0)] * 100),
                xy=(twowave.eng_strain[twowave.eng_stress.argmax(axis=0)], np.max(twowave.eng_stress)), textcoords='data')
    plt.title(testname + ' (2-Wave Analysis)')
    plt.xlabel('Eng. Strain (mm/mm)')
    plt.ylabel('Eng. Stress (MPa)')
    plt.legend()
    plt.grid()
    plt.savefig(save_dir + '2-Wave_Eng_Stress-Strain.png')
    # plt.show(block=True)     # For brittle samples there won't be much difference, so we'll just skip showing this plot
    plt.close()
# True stress-strain
    fig, ax = plt.subplots()
    ax.plot(twowave.true_strain, twowave.true_stress)
    ax.plot(twowave.true_strain[twowave.true_stress.argmax(axis=0)], np.max(twowave.true_stress), 'ro')
    ax.annotate('max=%.1f MPa, %.1f pct' % (np.max(twowave.true_stress), twowave.true_strain[twowave.true_stress.argmax(axis=0)] * 100),
                xy=(twowave.true_strain[twowave.true_stress.argmax(axis=0)], np.max(twowave.true_stress)), textcoords='data')
    plt.title(testname + ' (2-Wave Analysis)')
    plt.xlabel('True Strain (mm/mm)')
    plt.ylabel('True Stress (MPa)')
    plt.grid()
    plt.savefig(save_dir + '2-Wave_True_Stress-Strain.png')
    #plt.show(block=True)   # DEBUG ONLY
    plt.close()
# True stress vs. strain rate; Generally unnecessary
    # fig4, ax = plt.subplots()
    # ax.plot(twowave.true_strain, twowave.true_rate)
    ##ax.plot(twowave.true_strain.argmax(axis=0)], np.max(twowave.true_rate), 'ro')
    ##ax.annotate('max=%.1f MPa, %.1f pct' % (np.max(twowave.true_rate), twowave.true_strain[twowave.true_rate.argmax(axis=0)]*100), xy=(twowave.true_strain[twowave.true_rate.argmax(axis=0)],np.max(twowave.true_rate)), textcoords='data')
    # plt.title(testname)
    # plt.xlabel('True Strain (mm/mm)')
    # plt.ylabel('Strain Rate (/s)')
    # plt.grid()
    # plt.savefig(save_dir + '2-Wave_True_Strain-Rate.png')
    # plt.show()    # DEBUG ONLY
    # plt.close() #cleanup
# True strain vs. time; Generally unnecessary
    # fig5, ax = plt.subplots()
    # ax.plot(twowave.time, twowave.true_strain)
    # plt.title(testname)
    # plt.xlabel('Time (microsecond)')
    # plt.ylabel('True Strain (mm/mm)')
    # plt.grid()
    # plt.savefig(save_dir + '2-Wave_True_Strain-Time.png')
    # plt.show()    # DEBUG ONLY
    # plt.close() #cleanup

def plot_barforce(incident, transmitted, barforce, save_dir):
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(121, facecolor='#FFFFFF')
    ax.plot(incident.time, barforce.incident, 'b.', label='xti2_sigi')
    ax.plot(transmitted.time, barforce.transmitted, 'r.', label='xtt2_sigt')
    ax2 = fig.add_subplot(122, facecolor='#FFFFFF')
    ax2.plot(incident.time, barforce.difference, 'k', label='F_inc - F_trans')
    plt.title(' Force Balance Between Bars')
    ax.set_xlabel('Time (\u03BCs)')
    ax.set_ylabel('Bar Force (N)')
    ax2.set_xlabel('Time (\u03BCs)')
    ax2.set_ylabel('Difference in Bar Force (N)')
    plt.legend()
    plt.grid()
    plt.savefig(save_dir + 'StressesBetweenBars.png')
    #plt.show() # DEBUG ONLY
    plt.close()

def plot_3W(threewave, testname, save_dir):
    """Creates various plots for 3-wave analysis."""
    fig = plt.plot(threewave.true_strain, threewave.true_stress)
    plt.title(testname + ' (3-Wave Analysis)')
    plt.xlabel('True Strain (mm/mm)')
    plt.ylabel('True Stress (MPa)')
    plt.grid()
    plt.savefig(save_dir + '3-Wave_True_Stress-Strain.png')
    #plt.show() # DEBUG ONLY
    plt.close()
