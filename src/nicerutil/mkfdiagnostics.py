"""
Module to create some diagnostics plots for NICER - mainly developed to look
at the light-leak issue, but could also be useful to study one's own data to make
sensibile decisions on cleaning
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import argparse

from nicerutil.nicermkf import MkfFileOps, readmkffile, define_nicerdetloc

sys.dont_write_bytecode = True


#
def mkf_diagnostics(mkfFile, sunshine=2, sunAngLR=45, sunAngUR=180, moonAngLR=0, moonAngUR=180, brearthLR=0,
                    brearthUR=180, sunAzLR=-180, sunAzUR=180, timepostleak=True, writetocsv=False):
    # read mkf
    mkf_table = readmkffile(mkfFile)

    # Tracking filter
    trackingfiltermkf_table = MkfFileOps(mkf_table).trackingfiltermkf()
    print('Size of MKF file after standard SAA, Pointing, and On-Target tracking filtering is {}'.format(
        np.shape(trackingfiltermkf_table)[0]))

    # Writing initial table after tracking filtering to CSV file
    if writetocsv:
        print('Writing it to CSV file...')
        MkfFileOps(trackingfiltermkf_table).write_mkf_to_csv(mkfFile.split(".")[0] + '_mkf')

    # Time filtering
    if timepostleak:
        timefiltered_mkf = MkfFileOps(trackingfiltermkf_table).timefiltermkf(
            gtilist=np.array([[296229602.000, 596229602.000]]))
        print('Size of MKF file post-leak {}'.format(np.shape(timefiltered_mkf)[0]))
    else:
        timefiltered_mkf = trackingfiltermkf_table

    # Sunshine filtering
    sunshinefiltered_mkf = MkfFileOps(timefiltered_mkf).sunshinefiltermkf(sunshine=sunshine)
    print('Size of MKF file after sunshine filtering {}'.format(np.shape(sunshinefiltered_mkf)[0]))

    # Sun angle filtering
    sunanglefiltered_mkf = MkfFileOps(sunshinefiltered_mkf).sunanglefiltermkf(sunang_ll=sunAngLR, sunang_ul=sunAngUR)
    print('Size of MKF file for sun angle filtering [{}, {}] degrees = {}'.format(str(sunAngLR), str(sunAngUR),
                                                                                  np.shape(sunanglefiltered_mkf)[0]))

    # Moon angle filtering
    moonanglefiltered_mkf = MkfFileOps(sunanglefiltered_mkf).moonanglefiltermkf(moonang_ll=moonAngLR,
                                                                                moonang_ul=moonAngUR)
    print('Size of MKF file for moon angle filtering [{}, {}] degrees = {}'.format(str(moonAngLR), str(moonAngUR),
                                                                                   np.shape(sunanglefiltered_mkf)[0]))

    # Bright earth filtering
    brightearthfiltered_mkf = (MkfFileOps(moonanglefiltered_mkf).brightearthanglefiltermkf(brightearth_ll=brearthLR,
                                                                                           brightearth_ul=brearthUR))
    print('Size of MKF file for bright earth angle filtering [{}, {}] degrees = {}'.format(str(brearthLR),
                                                                                           str(brearthUR),
                                                                                           np.shape(
                                                                                               brightearthfiltered_mkf)[
                                                                                               0]))

    # Sun Azimuth (clocking) filtering
    sunazfiltered_mkf = MkfFileOps(brightearthfiltered_mkf).sunazfiltermkf(sunaz_ll=sunAzLR, sunaz_ul=sunAzUR)
    print('Size of MKF file for sun azimuth filtering [{}, {}] degrees = {}'.format(str(sunAzLR), str(sunAzUR),
                                                                                    np.shape(sunazfiltered_mkf)[0]))

    # Checking if the dataframe after filtering is empty or not
    if sunazfiltered_mkf.empty:
        print('DataFrame after all filtering is empty - Exiting')
        return

    # Spread of time interval that defines our filtering criteria
    timespread = (np.max(np.sort(sunazfiltered_mkf['tNICERmkf'])) -
                  np.mean(np.sort(sunazfiltered_mkf['tNICERmkf']))) / 86400
    print('Spread in time after filtering is {} days'.format(timespread))
    exposure = (len(sunazfiltered_mkf['tNICERmkf']))
    print('Exposure after filtering is {} seconds'.format(exposure))
    avg_moonang = np.mean(sunazfiltered_mkf['MOON_ANGLE'])
    print('Average Moon angle is {} degrees'.format(avg_moonang))
    avg_brightearth = np.mean(sunazfiltered_mkf['brightEarth'])
    print('Average Bright earth angle is {} degrees'.format(avg_brightearth))
    avg_elevation = np.mean(sunazfiltered_mkf['elevation'])
    print('Average elevation angle is {} degrees'.format(avg_elevation))

    average_ancilliary_info = {'timespread': timespread, 'avg_moon': avg_moonang,
                               'avg_brightearth': avg_brightearth, 'avg_elevation': avg_elevation}

    # Defining NICER detectors in geographical location
    # We will be plotting things per NICER detector
    nicDET_geograph = define_nicerdetloc()

    # Undershoot parameters
    under_perFPM_mean = np.zeros(len(nicDET_geograph))
    under_perFPM_stdv = np.zeros(len(nicDET_geograph))
    under_perFPM_median = np.zeros(len(nicDET_geograph))
    for ll, det_num in enumerate(nicDET_geograph):
        under_perFPM_mean[ll] = sunazfiltered_mkf['FPM_under' + det_num].mean()
        under_perFPM_stdv[ll] = sunazfiltered_mkf['FPM_under' + det_num].std()
        under_perFPM_median[ll] = sunazfiltered_mkf['FPM_under' + det_num].median()

    average_undershoot_perFPM = pd.DataFrame(np.vstack((under_perFPM_mean, under_perFPM_stdv, under_perFPM_median)).T,
                                             columns=["average", "stdv", "median"], index=nicDET_geograph)

    return sunazfiltered_mkf, average_undershoot_perFPM, average_ancilliary_info


def plot_under_sunAz_sunAngle(sunazfiltered_mkf, nicDET_geograph, outputfile):
    # Deriving maxundershoot (for plotting purposes)
    all_undershoots = []
    for ll, det_num in enumerate(nicDET_geograph):
        all_undershoots = np.append(all_undershoots, sunazfiltered_mkf['FPM_under' + det_num])
    maxundershoot = np.nanmax(all_undershoots)

    SUN_ANGLE = sunazfiltered_mkf['SUN_ANGLE']
    AZ_SUN = sunazfiltered_mkf['AZ_SUN']

    # Difining the plot axes
    fig, axs = plt.subplots(7, 8, figsize=(50, 50), dpi=100, facecolor='w', edgecolor='k', sharex=True,
                            gridspec_kw={'width_ratios': [1, 1, 1, 1, 1, 1, 1, 1],
                                         'height_ratios': [1, 1, 1, 1, 1, 1, 1]})
    axs = axs.ravel()

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.05, right=0.9, left=0.05, bottom=0.05)

    # Plot each detector separately
    for ll, det_num in enumerate(nicDET_geograph):

        axs[ll].tick_params(axis='both', labelsize=40)
        axs[ll].xaxis.offsetText.set_fontsize(40)
        axs[ll].ticklabel_format(style='plain', axis='y', scilimits=(0, 0), useMathText=True)
        axs[ll].xaxis.offsetText.set_fontsize(40)
        axs[ll].yaxis.offsetText.set_fontsize(40)

        # Undershoot per FPM
        FPM_under_perFPM = sunazfiltered_mkf['FPM_under' + det_num]

        axs[ll].scatter(AZ_SUN, FPM_under_perFPM, c=SUN_ANGLE, cmap='copper', label='FPM' + det_num)

        for axis in ['top', 'bottom', 'left', 'right']:
            axs[ll].spines[axis].set_linewidth(2)
            axs[ll].tick_params(width=2)

        axs[ll].legend(loc='upper right', fontsize=40, frameon=False, markerscale=0)

        axs[ll].set_ylim(0, maxundershoot)

        if det_num in ["07", "16", "17", "27", "37", "47", "57",
                       "15", "25", "26", "35", "36", "46", "56",
                       "14", "24", "34", "44", "45", "54", "55",
                       "13", "23", "33", "43", "53", "66", "67",
                       "12", "22", "32", "42", "52", "64", "65",
                       "11", "21", "31", "41", "51", "62", "63",
                       "10", "20", "30", "40", "50", "60", "61"]:
            axs[ll].set_yticklabels([])

    # axes labels
    axs[52].set_xlabel('Sun clocking/Azimuth angle (degrees)', fontsize=40)
    axs[24].set_ylabel('Under_count per FPM (counts)', fontsize=40)

    # Creating map for color bar
    map1 = axs[55].imshow(np.stack([SUN_ANGLE, SUN_ANGLE]), cmap='copper', aspect='auto')

    # position for the colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(map1, cax=cbaxes)
    cbar.ax.tick_params(labelsize=40)
    cbar.set_label("Sun Angle", fontsize=40)

    # Saving figure
    plotName = outputfile + '_under_sunAz_sunAngle.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def plot_under_sunAz_moonAngle(sunazfiltered_mkf, nicDET_geograph, outputfile):
    # Deriving maxundershoot (for plotting purposes)
    all_undershoots = []
    for ll, det_num in enumerate(nicDET_geograph):
        all_undershoots = np.append(all_undershoots, sunazfiltered_mkf['FPM_under' + det_num])
    maxundershoot = np.nanmax(all_undershoots)

    MOON_ANGLE = sunazfiltered_mkf['MOON_ANGLE']
    AZ_SUN = sunazfiltered_mkf['AZ_SUN']

    # Difining the plot axes
    fig, axs = plt.subplots(7, 8, figsize=(50, 50), dpi=100, facecolor='w', edgecolor='k', sharex=True,
                            gridspec_kw={'width_ratios': [1, 1, 1, 1, 1, 1, 1, 1],
                                         'height_ratios': [1, 1, 1, 1, 1, 1, 1]})
    axs = axs.ravel()

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.05, right=0.9, left=0.05, bottom=0.05)

    # Plot each detector separately
    for ll, det_num in enumerate(nicDET_geograph):

        axs[ll].tick_params(axis='both', labelsize=40)
        axs[ll].xaxis.offsetText.set_fontsize(40)
        axs[ll].ticklabel_format(style='plain', axis='y', scilimits=(0, 0), useMathText=True)
        axs[ll].xaxis.offsetText.set_fontsize(40)
        axs[ll].yaxis.offsetText.set_fontsize(40)

        # Undershoot per FPM
        FPM_under_perFPM = sunazfiltered_mkf['FPM_under' + det_num]

        axs[ll].scatter(AZ_SUN, FPM_under_perFPM, c=MOON_ANGLE, cmap='copper', label='FPM' + det_num)

        for axis in ['top', 'bottom', 'left', 'right']:
            axs[ll].spines[axis].set_linewidth(2)
            axs[ll].tick_params(width=2)

        axs[ll].legend(loc='upper right', fontsize=40, frameon=False, markerscale=0)

        axs[ll].set_ylim(0, maxundershoot)

        if det_num in ["07", "16", "17", "27", "37", "47", "57",
                       "15", "25", "26", "35", "36", "46", "56",
                       "14", "24", "34", "44", "45", "54", "55",
                       "13", "23", "33", "43", "53", "66", "67",
                       "12", "22", "32", "42", "52", "64", "65",
                       "11", "21", "31", "41", "51", "62", "63",
                       "10", "20", "30", "40", "50", "60", "61"]:
            axs[ll].set_yticklabels([])

    # axes labels
    axs[52].set_xlabel('Sun clocking/Azimuth angle (degrees)', fontsize=40)
    axs[24].set_ylabel('Under_count per FPM (counts)', fontsize=40)

    # Creating map for color bar
    map1 = axs[55].imshow(np.stack([MOON_ANGLE, MOON_ANGLE]), cmap='copper', aspect='auto')

    # position for the colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(map1, cax=cbaxes)
    cbar.ax.tick_params(labelsize=40)
    cbar.set_label("Moon Angle", fontsize=40)

    # Saving figure
    plotName = outputfile + '_under_sunAz_moonAngle.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def plot_under_sunAz_brightearth(sunazfiltered_mkf, nicDET_geograph, outputfile):
    # Deriving maxundershoot (for plotting purposes)
    all_undershoots = []
    for ll, det_num in enumerate(nicDET_geograph):
        all_undershoots = np.append(all_undershoots, sunazfiltered_mkf['FPM_under' + det_num])
    maxundershoot = np.nanmax(all_undershoots)

    brightEarth = sunazfiltered_mkf['brightEarth']
    AZ_SUN = sunazfiltered_mkf['AZ_SUN']

    # Difining the plot axes
    fig, axs = plt.subplots(7, 8, figsize=(50, 50), dpi=100, facecolor='w', edgecolor='k', sharex=True,
                            gridspec_kw={'width_ratios': [1, 1, 1, 1, 1, 1, 1, 1],
                                         'height_ratios': [1, 1, 1, 1, 1, 1, 1]})
    axs = axs.ravel()

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.05, right=0.9, left=0.05, bottom=0.05)

    # Plot each detector separately
    for ll, det_num in enumerate(nicDET_geograph):

        axs[ll].tick_params(axis='both', labelsize=40)
        axs[ll].xaxis.offsetText.set_fontsize(40)
        axs[ll].ticklabel_format(style='plain', axis='y', scilimits=(0, 0), useMathText=True)
        axs[ll].xaxis.offsetText.set_fontsize(40)
        axs[ll].yaxis.offsetText.set_fontsize(40)

        # Undershoot per FPM
        FPM_under_perFPM = sunazfiltered_mkf['FPM_under' + det_num]

        axs[ll].scatter(AZ_SUN, FPM_under_perFPM, c=brightEarth, cmap='copper', label='FPM' + det_num)

        for axis in ['top', 'bottom', 'left', 'right']:
            axs[ll].spines[axis].set_linewidth(2)
            axs[ll].tick_params(width=2)

        axs[ll].legend(loc='upper right', fontsize=40, frameon=False, markerscale=0)

        axs[ll].set_ylim(0, maxundershoot)

        if det_num in ["07", "16", "17", "27", "37", "47", "57",
                       "15", "25", "26", "35", "36", "46", "56",
                       "14", "24", "34", "44", "45", "54", "55",
                       "13", "23", "33", "43", "53", "66", "67",
                       "12", "22", "32", "42", "52", "64", "65",
                       "11", "21", "31", "41", "51", "62", "63",
                       "10", "20", "30", "40", "50", "60", "61"]:
            axs[ll].set_yticklabels([])

    # axes labels
    axs[52].set_xlabel('Sun clocking/Azimuth angle (degrees)', fontsize=40)
    axs[24].set_ylabel('Under_count per FPM (counts)', fontsize=40)

    # Creating map for color bar
    map1 = axs[55].imshow(np.stack([brightEarth, brightEarth]), cmap='copper', aspect='auto')

    # position for the colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(map1, cax=cbaxes)
    cbar.ax.tick_params(labelsize=40)
    cbar.set_label("Bright Earth", fontsize=40)

    # Saving figure
    plotName = outputfile + '_under_sunAz_brightEarth.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def plot_under_sunAz_elevation(sunazfiltered_mkf, nicDET_geograph, outputfile):
    # Deriving maxundershoot (for plotting purposes)
    all_undershoots = []
    for ll, det_num in enumerate(nicDET_geograph):
        all_undershoots = np.append(all_undershoots, sunazfiltered_mkf['FPM_under' + det_num])
    maxundershoot = np.nanmax(all_undershoots)

    elevation = sunazfiltered_mkf['elevation']
    AZ_SUN = sunazfiltered_mkf['AZ_SUN']

    # Difining the plot axes
    fig, axs = plt.subplots(7, 8, figsize=(50, 50), dpi=100, facecolor='w', edgecolor='k', sharex=True,
                            gridspec_kw={'width_ratios': [1, 1, 1, 1, 1, 1, 1, 1],
                                         'height_ratios': [1, 1, 1, 1, 1, 1, 1]})
    axs = axs.ravel()

    fig.tight_layout()
    fig.subplots_adjust(wspace=0.05, hspace=0.05, right=0.9, left=0.05, bottom=0.05)

    # Plot each detector separately
    for ll, det_num in enumerate(nicDET_geograph):

        axs[ll].tick_params(axis='both', labelsize=40)
        axs[ll].xaxis.offsetText.set_fontsize(40)
        axs[ll].ticklabel_format(style='plain', axis='y', scilimits=(0, 0), useMathText=True)
        axs[ll].xaxis.offsetText.set_fontsize(40)
        axs[ll].yaxis.offsetText.set_fontsize(40)

        # Undershoot per FPM
        FPM_under_perFPM = sunazfiltered_mkf['FPM_under' + det_num]

        axs[ll].scatter(AZ_SUN, FPM_under_perFPM, c=elevation, cmap='copper', label='FPM' + det_num)

        for axis in ['top', 'bottom', 'left', 'right']:
            axs[ll].spines[axis].set_linewidth(2)
            axs[ll].tick_params(width=2)

        axs[ll].legend(loc='upper right', fontsize=40, frameon=False, markerscale=0)

        axs[ll].set_ylim(0, maxundershoot)

        if det_num in ["07", "16", "17", "27", "37", "47", "57",
                       "15", "25", "26", "35", "36", "46", "56",
                       "14", "24", "34", "44", "45", "54", "55",
                       "13", "23", "33", "43", "53", "66", "67",
                       "12", "22", "32", "42", "52", "64", "65",
                       "11", "21", "31", "41", "51", "62", "63",
                       "10", "20", "30", "40", "50", "60", "61"]:
            axs[ll].set_yticklabels([])

    # axes labels
    axs[52].set_xlabel('Sun clocking/Azimuth angle (degrees)', fontsize=40)
    axs[24].set_ylabel('Under_count per FPM (counts)', fontsize=40)

    # Creating map for color bar
    map1 = axs[55].imshow(np.stack([elevation, elevation]), cmap='copper', aspect='auto')

    # position for the colorbar
    cbaxes = fig.add_axes([0.91, 0.1, 0.02, 0.8])
    cbar = fig.colorbar(map1, cax=cbaxes)
    cbar.ax.tick_params(labelsize=40)
    cbar.set_label("Elevation", fontsize=40)

    # Saving figure
    plotName = outputfile + '_under_sunAz_elevation.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def plot_sunAz_under(sunazfiltered_mkf, nicDET_geograph, outputfile):
    AZ_SUN = sunazfiltered_mkf['AZ_SUN']

    # Define the axes
    fig, ax1 = plt.subplots(1, figsize=(6, 4), dpi=80, facecolor='w', edgecolor='k')
    ax1.tick_params(axis='both', labelsize=8)
    ax1.set_xlabel('Sun clocking/Azimuth angle (degrees)', fontsize=8)
    ax1.set_ylabel('Under_count per FPM (counts)', fontsize=8)
    ax1.xaxis.offsetText.set_fontsize(8)
    colCycle = plt.cm.brg(np.linspace(0, 1, len(nicDET_geograph)))

    # Plot each detector separately
    for ll, det_num in enumerate(nicDET_geograph):
        FPM_under_perFPM = sunazfiltered_mkf['FPM_under' + det_num]
        ax1.plot(AZ_SUN, FPM_under_perFPM, color=colCycle[ll], marker='.', markersize='2', label='FPM' + det_num,
                 alpha=0.5, ls='')

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.tick_params(width=2)

    # ax1.legend(loc='upper right', fontsize=40, frameon=False, markerscale=0)
    ax1Leg = ax1.legend(loc='lower right', fontsize=5, bbox_to_anchor=(1.3, 0), framealpha=None, ncol=2)
    ax1Leg.get_frame().set_linewidth(1)
    ax1Leg.get_frame().set_edgecolor('k')

    # Saving figure
    fig.tight_layout()
    plotName = outputfile + '_sunAz_under.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def plot_sunang_under(sunazfiltered_mkf, nicDET_geograph, outputfile):
    SUN_ANGLE = sunazfiltered_mkf['SUN_ANGLE']

    # Define the axes
    fig, ax1 = plt.subplots(1, figsize=(6, 4), dpi=80, facecolor='w', edgecolor='k')
    ax1.tick_params(axis='both', labelsize=8)
    ax1.set_xlabel('Sun angle (degrees)', fontsize=8)
    ax1.set_ylabel('Under_count per FPM (counts)', fontsize=8)
    ax1.xaxis.offsetText.set_fontsize(8)
    colCycle = plt.cm.brg(np.linspace(0, 1, len(nicDET_geograph)))

    # Plot each detector separately
    for ll, det_num in enumerate(nicDET_geograph):
        FPM_under_perFPM = sunazfiltered_mkf['FPM_under' + det_num]
        ax1.plot(SUN_ANGLE, FPM_under_perFPM, color=colCycle[ll], marker='.', markersize='2',
                 label='FPM' + det_num, alpha=0.5, ls='')

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.tick_params(width=2)

    # ax1.legend(loc='upper right', fontsize=40, frameon=False, markerscale=0)
    ax1Leg = ax1.legend(loc='lower right', fontsize=5, bbox_to_anchor=(1.3, 0), framealpha=None, ncol=2)
    ax1Leg.get_frame().set_linewidth(1)
    ax1Leg.get_frame().set_edgecolor('k')

    # Saving figure
    fig.tight_layout()
    plotName = outputfile + '_sunangle_under.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def plot_moonangle_under(sunazfiltered_mkf, nicDET_geograph, outputfile):
    MOON_ANGLE = sunazfiltered_mkf['MOON_ANGLE']

    # Define the axes
    fig, ax1 = plt.subplots(1, figsize=(6, 4), dpi=80, facecolor='w', edgecolor='k')
    ax1.tick_params(axis='both', labelsize=8)
    ax1.set_xlabel('Moon angle (degrees)', fontsize=8)
    ax1.set_ylabel('Under_count per FPM (counts)', fontsize=8)
    ax1.xaxis.offsetText.set_fontsize(8)
    colCycle = plt.cm.brg(np.linspace(0, 1, len(nicDET_geograph)))

    # Plot each detector separately
    for ll, det_num in enumerate(nicDET_geograph):
        FPM_under_perFPM = sunazfiltered_mkf['FPM_under' + det_num]
        ax1.plot(MOON_ANGLE, FPM_under_perFPM, color=colCycle[ll], marker='.', markersize='2',
                 label='FPM' + det_num, alpha=0.5, ls='')

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.tick_params(width=2)

    # ax1.legend(loc='upper right', fontsize=40, frameon=False, markerscale=0)
    ax1Leg = ax1.legend(loc='lower right', fontsize=5, bbox_to_anchor=(1.3, 0), framealpha=None, ncol=2)
    ax1Leg.get_frame().set_linewidth(1)
    ax1Leg.get_frame().set_edgecolor('k')

    # Saving figure
    fig.tight_layout()
    plotName = outputfile + '_moonangle_under.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def plot_elevation_under(sunazfiltered_mkf, nicDET_geograph, outputfile):
    elevation = sunazfiltered_mkf['elevation']

    # Define the axes
    fig, ax1 = plt.subplots(1, figsize=(6, 4), dpi=80, facecolor='w', edgecolor='k')
    ax1.tick_params(axis='both', labelsize=8)
    ax1.set_xlabel('Elevation (degrees)', fontsize=8)
    ax1.set_ylabel('Under_count per FPM (counts)', fontsize=8)
    ax1.xaxis.offsetText.set_fontsize(8)
    colCycle = plt.cm.brg(np.linspace(0, 1, len(nicDET_geograph)))

    # Plot each detector separately
    for ll, det_num in enumerate(nicDET_geograph):
        FPM_under_perFPM = sunazfiltered_mkf['FPM_under' + det_num]
        ax1.plot(elevation, FPM_under_perFPM, color=colCycle[ll], marker='.', markersize='2',
                 label='FPM' + det_num, alpha=0.5, ls='')

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.tick_params(width=2)

    # ax1.legend(loc='upper right', fontsize=40, frameon=False, markerscale=0)
    ax1Leg = ax1.legend(loc='lower right', fontsize=5, bbox_to_anchor=(1.3, 0), framealpha=None, ncol=2)
    ax1Leg.get_frame().set_linewidth(1)
    ax1Leg.get_frame().set_edgecolor('k')

    # Saving figure
    fig.tight_layout()
    plotName = outputfile + '_elevation_under.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def plot_brightearth_under(sunazfiltered_mkf, nicDET_geograph, outputfile):
    brightEarth = sunazfiltered_mkf['brightEarth']

    # Define the axes
    fig, ax1 = plt.subplots(1, figsize=(6, 4), dpi=80, facecolor='w', edgecolor='k')
    ax1.tick_params(axis='both', labelsize=8)
    ax1.set_xlabel('Bright Earth angle (degrees)', fontsize=8)
    ax1.set_ylabel('Under_count per FPM (counts)', fontsize=8)
    ax1.xaxis.offsetText.set_fontsize(8)
    colCycle = plt.cm.brg(np.linspace(0, 1, len(nicDET_geograph)))

    # Plot each detector separately
    for ll, det_num in enumerate(nicDET_geograph):
        FPM_under_perFPM = sunazfiltered_mkf['FPM_under' + det_num]
        ax1.plot(brightEarth, FPM_under_perFPM, color=colCycle[ll], marker='.', markersize='2',
                 label='FPM' + det_num, alpha=0.5, ls='')

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.tick_params(width=2)

    # ax1.legend(loc='upper right', fontsize=40, frameon=False, markerscale=0)
    ax1Leg = ax1.legend(loc='lower right', fontsize=5, bbox_to_anchor=(1.3, 0), framealpha=None, ncol=2)
    ax1Leg.get_frame().set_linewidth(1)
    ax1Leg.get_frame().set_edgecolor('k')

    # Saving figure
    fig.tight_layout()
    plotName = outputfile + '_brightearth_under.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def plot_averageunder_perfpm(average_undershoot_perFPM, nicDET_geograph, outputfile):
    # Defining plot axes
    fig, ax1 = plt.subplots(1, figsize=(12, 4), dpi=80, facecolor='w', edgecolor='k')
    ax1.tick_params(axis='both', labelsize=8)
    ax1.set_xlabel('FPM', fontsize=8)
    ax1.set_ylabel('Under_count per FPM (counts)', fontsize=8)
    ax1.xaxis.offsetText.set_fontsize(8)

    under_perFPM_mean = average_undershoot_perFPM["average"]
    under_perFPM_stdv = average_undershoot_perFPM["stdv"]

    unique, rev = np.unique(nicDET_geograph, return_inverse=True)

    ax1.set_xticks(range(len(unique)))
    ax1.set_xticklabels(unique)

    ax1.errorbar(rev, under_perFPM_mean, yerr=under_perFPM_stdv, color='k', fmt='s', zorder=10, markersize=5)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.tick_params(width=2)

    # ax1.set_ylim(0,maxundershoot)

    # Saving figure
    fig.tight_layout()
    plotName = outputfile + '_averageunder_perfpm.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def plot_medianunder_perfpm(average_undershoot_perFPM, nicDET_geograph, outputfile):
    # Defining plot axes
    fig, ax1 = plt.subplots(1, figsize=(12, 4), dpi=80, facecolor='w', edgecolor='k')
    ax1.tick_params(axis='both', labelsize=8)
    ax1.set_xlabel('FPM', fontsize=8)
    ax1.set_ylabel('Under_count per FPM (counts)', fontsize=8)
    ax1.xaxis.offsetText.set_fontsize(8)

    under_perFPM_median = average_undershoot_perFPM['median']

    unique, rev = np.unique(nicDET_geograph, return_inverse=True)

    ax1.set_xticks(range(len(unique)))
    ax1.set_xticklabels(unique)

    ax1.plot(rev, under_perFPM_median, marker='s', color='k', zorder=10, markersize=5, ls='')

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.tick_params(width=2)

    # Saving figure
    fig.tight_layout()
    plotName = outputfile + '_medianunder_perfpm.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def plot_moonang_dist(sunazfiltered_mkf, outputfile):
    MOON_ANGLE = sunazfiltered_mkf['MOON_ANGLE']
    # Define the axes
    fig, ax1 = plt.subplots(1, figsize=(6, 4), dpi=80, facecolor='w', edgecolor='k')
    ax1.tick_params(axis='both', labelsize=8)
    ax1.set_xlabel('Moon angle (degrees)', fontsize=8)
    ax1.set_ylabel('Count', fontsize=8)
    ax1.xaxis.offsetText.set_fontsize(8)

    # Plot each detector separately
    ax1.hist(MOON_ANGLE, bins=int(len(MOON_ANGLE) / 25), density=False, color='k', label='Moon angle')

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.tick_params(width=2)

    ax1.legend(loc='upper right', fontsize=8, frameon=False, markerscale=0)

    # Saving figure
    fig.tight_layout()
    plotName = outputfile + '_moon_dist.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def plot_brightearth_dist(sunazfiltered_mkf, outputfile):
    brightEarth = sunazfiltered_mkf['brightEarth']
    # Define the axes
    fig, ax1 = plt.subplots(1, figsize=(6, 4), dpi=80, facecolor='w', edgecolor='k')
    ax1.tick_params(axis='both', labelsize=8)
    ax1.set_xlabel('Bright Earth angle (degrees)', fontsize=8)
    ax1.set_ylabel('Count', fontsize=8)
    ax1.xaxis.offsetText.set_fontsize(8)

    # Plot each detector separately
    ax1.hist(brightEarth, bins=int(len(brightEarth) / 25), density=False, color='k', label='Bright earth')

    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
        ax1.tick_params(width=2)

    ax1.legend(loc='upper right', fontsize=8, frameon=False, markerscale=0)

    # Saving figure
    fig.tight_layout()
    plotName = outputfile + '_brightEarth_dist.png'
    fig.savefig(plotName, format='png', dpi=200)
    plt.close()

    return


def createalldiagnosticsplots(mkffile, average_undershoot_perFPM, nicDET_geograph, outputfile):
    # Create under_sunAz_sunAngle plot
    plot_under_sunAz_sunAngle(mkffile, nicDET_geograph, outputfile=outputfile)

    # Create under_sunAz_moonAngle plot
    plot_under_sunAz_moonAngle(mkffile, nicDET_geograph, outputfile=outputfile)

    # Create under_sunAz_brightearth plot
    plot_under_sunAz_brightearth(mkffile, nicDET_geograph, outputfile=outputfile)

    # Create under_sunAz_elevation plot
    plot_under_sunAz_elevation(mkffile, nicDET_geograph, outputfile=outputfile)

    # Create sunAz_under plot
    plot_sunAz_under(mkffile, nicDET_geograph, outputfile=outputfile)

    # Create averageunder_perfpm plot
    plot_averageunder_perfpm(average_undershoot_perFPM, nicDET_geograph, outputfile=outputfile)

    # Create medianunder_perfpm plot
    plot_medianunder_perfpm(average_undershoot_perFPM, nicDET_geograph, outputfile=outputfile)

    # Create Sun angle perfpm plot
    plot_sunang_under(mkffile, nicDET_geograph, outputfile=outputfile)

    # Create Moon angle perfpm plot
    plot_moonangle_under(mkffile, nicDET_geograph, outputfile=outputfile)

    # Create Elevation angle perfpm plot
    plot_elevation_under(mkffile, nicDET_geograph, outputfile=outputfile)

    # Create Bright Earth angle perfpm plot
    plot_brightearth_under(mkffile, nicDET_geograph, outputfile=outputfile)

    # Create moon angle distribution plot
    plot_moonang_dist(mkffile, outputfile=outputfile)

    # Create bright earth distribution plot
    plot_brightearth_dist(mkffile, outputfile=outputfile)


def main():
    parser = argparse.ArgumentParser(description="Diagnose using MKF file")
    parser.add_argument("mkfFile", help="A NICER MKF file", type=str)
    parser.add_argument("-ss", "--sunshine", help="Filtering for sunshine, 0 for night, 1 for day, "
                                                  "and 2 for no filtering (default=1)", type=int, default=2)
    parser.add_argument("-sl", "--sunAngLR", help="Filtering for sun angle, lower-range", type=float,
                        default=45)
    parser.add_argument("-su", "--sunAngUR", help="Filtering for sun angle, upper-range", type=float,
                        default=180)

    parser.add_argument("-ml", "--moonAngLR", help="Filtering for moon angle, lower-range", type=float,
                        default=0)
    parser.add_argument("-mu", "--moonAngUR", help="Filtering for moon angle, upper-range", type=float,
                        default=180)
    parser.add_argument("-bel", "--brearthLR", help="Filtering for bright Earth angle, lower-range",
                        type=float, default=0)
    parser.add_argument("-beu", "--brearthUR", help="Filtering for bright Earth angle, upper-range",
                        type=float, default=180)
    parser.add_argument("-al", "--sunAzLR", help="Filtering for sun azimuth (clocking) angle, "
                                                 "lower-range", type=float,
                        default=-180)
    parser.add_argument("-au", "--sunAzUR", help="Filtering for sun azimuth (clocking) angle, "
                                                 "upper-range", type=float,
                        default=180)
    parser.add_argument("-tf", "--timepostleak", help="Accept only time post repair (default=True)",
                        type=bool, default=True, action=argparse.BooleanOptionalAction)
    parser.add_argument("-wc", "--writetocsv", help="Bool to write mkf to .txt file", type=bool,
                        default=False, action=argparse.BooleanOptionalAction)
    parser.add_argument("-of", "--outputfile", help="name of output .png light curve showing flagged "
                                                    "bins", type=str, default='mkf_diagnostics')
    args = parser.parse_args()

    sunazfiltered_mkf, average_undershoot_perFPM, average_ancilliary_info = mkf_diagnostics(args.mkfFile, args.sunshine,
                                                                                            args.sunAngLR,
                                                                                            args.sunAngUR,
                                                                                            args.moonAngLR,
                                                                                            args.moonAngUR,
                                                                                            args.brearthLR,
                                                                                            args.brearthUR,
                                                                                            args.sunAzLR, args.sunAzUR,
                                                                                            args.timepostleak,
                                                                                            args.writetocsv)
    # Checking if the dataframe after filtering is empty or not
    if sunazfiltered_mkf.empty:
        print('DataFrame after all filtering is empty - Exiting')
        return

    nicDET_geograph = define_nicerdetloc()

    # Creating the diagnostics plots
    createalldiagnosticsplots(sunazfiltered_mkf, average_undershoot_perFPM, nicDET_geograph, args.outputfile)

    return


if __name__ == '__main__':
    main()
