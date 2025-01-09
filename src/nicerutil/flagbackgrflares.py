"""
Module to help clean NICER data from high energy flares (known as PREL, see
https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/flares/). Given
a NICER event file, this module will create a light curve in the energy range 
12-15 keV with a time-bin of 5 seconds (both of which optional inputs), where 
NICER effective area is practically 0. It will then search for bins with a number
of counts that is too large to be considered random Poisson fluctuation and flag
it as a potential background flare. If desired the user could request an xselect,
and NICERDAS-compatible GTI that corresponds to the flagged events. A plot of the
flagged events is also produced.

|
To do:
- Comparison to in-band light curve: create light curve in a sensible energy range
while highlighting eliminated "flare" bins
- derive proxy signal-to-noise ratio for in-band rate to background rate with and
without burst elimination
- Plot the mkf FPM_OVERONLY_COUNT and COR_SAX
- Create an overarching plot of all cleaning
- Summarizing the lost GTIs in seconds and fraction of total exposure

|
Warning: make sure you have initialized heasoft before running this module
"""
import matplotlib.pyplot as plt
import sys
import argparse
import os
import logging
import numpy as np

# Custom modules
from crimp.eventfile import EvtFileOps
from nicerutil.lightcurve import lightcurve
from nicerutil.bursts import burstsearch, mergesamebursts
from nicerutil.correctrateforfpmsel import correctrateforfpmsel

sys.dont_write_bytecode = True

# Log config
############
logFormatter = logging.Formatter('[%(asctime)s] %(levelname)8s %(message)s ' +
                                 '(%(filename)s:%(lineno)s)', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('nicerutil_log')
logger.setLevel(logging.DEBUG)

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
consoleHandler.setLevel(logging.WARNING)
logger.addHandler(consoleHandler)


def flagbackgrflares(eventfile, eneLow_back=12, eneHigh_back=15, timebin=5., lcthresh=0.1, probLim=0.01,
                     outputFile='flagged_flares', creategti=None, eneLow_src=1, eneHigh_src=5):
    """
    flagbackgrflares creates a light curves, and flags bins that are considered to be outliers
    according to a Poisson mass function
    :param eventfile: name of the fits event file
    :type eventfile: str
    :param eneLow_back: low energy cutoff
    :type eneLow_back: float
    :param eneHigh_back: high energy cutoff
    :type eneHigh_back: float
    :param timebin: time binsize of light curve
    :type timebin: float
    :param lcthresh: exclude time bins with exposure less than lcthresh*timebin
    :type lcthresh: float
    :param probLim: probability mass function below hich to consider a bin part of a burst
    :type probLim: float
    :param outputFile: Name of file to save the burst search (.txt and .pdf) (default = "flagged_flares")
    :type outputFile: str
    :param creategti: If given, name of .txt gti file that is compatible with xselect functionality maketime
    :type creategti: str
    :param eneLow_src: low energy cutoff for source (for plotting purposes only default = 1 keV)
    :type eneLow_src: float
    :param eneHigh_src: high energy cutoff for source (for plotting purposes only default = 5 keV)
    :type eneHigh_src: float
    :return: flares, distinct_flares (all falgged bins and a merged list as well if bins are consecutive)
    :rtype: list
    """
    # Reading data
    EF = EvtFileOps(eventfile)
    evtFileKeyWords, gtiList = EF.readGTI()

    # Dealing with the background
    #############################
    # Energy filtering, reading TIME column as numpy array
    dataTP_eneFlt = EF.filtenergy(eneLow=eneLow_back, eneHigh=eneHigh_back)
    TIME_back = dataTP_eneFlt['TIME'].to_numpy()

    # Create light curve
    binnedLC, _ = lightcurve(TIME_back, gtiList, timebin=timebin, lcthresh=lcthresh)
    binnedLC_corr = correctrateforfpmsel(eventfile, binnedLC)

    # FLare search
    burstsearch_result = burstsearch(binnedLC_corr, probLim=probLim, outputfile=outputFile)
    flares = burstsearch_result['bin_bursts']

    # Dealing with the source
    #############################
    # Creating a light curve within source energy range
    dataTP_eneFlt_src = EF.filtenergy(eneLow=eneLow_src, eneHigh=eneHigh_src)
    TIME_src = dataTP_eneFlt_src['TIME'].to_numpy()
    binnedLC_src, _ = lightcurve(TIME_src, gtiList, timebin=timebin, lcthresh=lcthresh)
    binnedLC_src_corr = correctrateforfpmsel(eventfile, binnedLC_src)
    
    # Create a simple diagnostic plot of the flares
    #############################
    plot_flare_diagnostics(binnedLC_corr, flares, binnedLC_src_corr, outputFile=outputFile)

    # Create .txt and if desired .fits GTI files
    #############################
    if not flares.empty:
        # Merge flares separated by timebin
        distinct_flares = mergesamebursts(flares, tsameburst=1.5 * timebin)

        # Creating a GTI file
        if creategti is not None:
            # Creating a GTI excluding these time intervals for later use with nicerl2
            f = open(creategti + ".txt", "w+")
            # Have to think about how heasoft likes the times written
            for jj, (burst_nbr, burst_df) in enumerate(distinct_flares.items()):
                # If one flare
                if (jj == 0) and (len(distinct_flares) == 1):
                    burststart_tmp = ((burst_df['lcBins'].loc[0]) - (burst_df['lcBinsRange'].loc[0] / 2))
                    burstend_tmp = ((burst_df['lcBins'].loc[burst_df.index[-1]]) +
                                    (burst_df['lcBinsRange'].loc[burst_df.index[-1]] / 2))
                    f.write('(TIME<{}).or.(TIME>{})'.format(burststart_tmp, burstend_tmp))
                # If first flare but multiple flares
                elif jj == 0 and (len(distinct_flares) != 1):
                    burststart_tmp = ((burst_df['lcBins'].loc[0]) - (burst_df['lcBinsRange'].loc[0] / 2))
                    burstend_tmp = ((burst_df['lcBins'].loc[burst_df.index[-1]]) +
                                    (burst_df['lcBinsRange'].loc[burst_df.index[-1]] / 2))
                    f.write('(TIME<{}).or.((TIME>{})'.format(burststart_tmp, burstend_tmp))
                # If last flare
                elif jj == (len(distinct_flares) - 1):
                    burststart_tmp = ((burst_df['lcBins'].loc[0]) - (burst_df['lcBinsRange'].loc[0] / 2))
                    burstend_tmp = ((burst_df['lcBins'].loc[burst_df.index[-1]]) +
                                    (burst_df['lcBinsRange'].loc[burst_df.index[-1]] / 2))
                    f.write('.and.(TIME<{})).or.(TIME>{})'.format(burststart_tmp, burstend_tmp))
                # If middle flares
                else:
                    burststart_tmp = ((burst_df['lcBins'].loc[0]) - (burst_df['lcBinsRange'].loc[0] / 2))
                    burstend_tmp = ((burst_df['lcBins'].loc[burst_df.index[-1]]) +
                                    (burst_df['lcBinsRange'].loc[burst_df.index[-1]] / 2))
                    f.write('.and.(TIME<{})).or.((TIME>{})'.format(burststart_tmp, burstend_tmp))

            f.close()

            command = 'maketime ' + eventfile + ' ' + creategti + '.fits @' + creategti + '.txt anything anything TIME no clobber=yes'
            os.system(command)

        # Adding some info to logging file about the cleaning process
        logger.info('\n Cleaning statistics:')
        logger.info('\n Average rate of background light curve before any cleaning : {}'.format(
            burstsearch_result['average_rate']))
        logger.info('\n Average rate of background light curve after cleaning : {}'.format(
            burstsearch_result['average_rate_nobursts']))

    else:
        distinct_flares = None
        logger.info('\n No flares found.\n')

    return flares, distinct_flares


def plot_flare_diagnostics(back_lightcurve, flare_lightcurve, src_lightcurve, outputFile='flare_diagnostics'):

    fig, axs = plt.subplots(2, 1, figsize=(10, 8), dpi=200, facecolor='w', edgecolor='k', sharex=True,
                            sharey=False, gridspec_kw={'width_ratios': [1], 'height_ratios': [1, 1]})
    axs = axs.ravel()

    # Plotting the "background" light curve,
    # i.e., out-of-bounds events, and those flagged as flares
    #########################################################
    axs[0].tick_params(axis='both', labelsize=12)
    axs[0].set_ylabel(r'$\,\mathrm{Rate\,(counts/\,s)}$', fontsize=12)
    axs[0].ticklabel_format(style='plain', axis='y', scilimits=(0, 0))
    axs[0].xaxis.offsetText.set_fontsize(12)
    axs[0].yaxis.offsetText.set_fontsize(12)
    # Plotting all out-of-band events
    axs[0].errorbar(back_lightcurve['lcBins'], back_lightcurve['ctrate'], xerr=back_lightcurve['lcBinsRange'] / 2,
                    yerr=back_lightcurve['ctrateErr'], fmt='ok', label='Good GTIs')
    # Plotting the flare interval
    axs[0].errorbar(flare_lightcurve['lcBins'], flare_lightcurve['ctrate'], xerr=flare_lightcurve['lcBinsRange'] / 2,
                    yerr=flare_lightcurve['ctrateErr'], fmt='or', label='Flare GTIs')
    axs[0].set_title(r'$\,\mathrm{Background\,energy\,range}$', fontsize=12)

    # Plotting the "source" light curve,
    # i.e., in-bound events, and highlight the bins that correspond to the flares
    #############################################################################
    axs[1].tick_params(axis='both', labelsize=12)
    axs[1].set_xlabel(r'$\,\mathrm{Time\,(MET,\,s)}$', fontsize=12)
    axs[1].set_ylabel(r'$\,\mathrm{Rate\,(counts/\,s)}$', fontsize=12)
    axs[1].ticklabel_format(style='plain', axis='y', scilimits=(0, 0))
    axs[1].xaxis.offsetText.set_fontsize(12)
    axs[1].yaxis.offsetText.set_fontsize(12)
    # Plotting all in-band events
    axs[1].errorbar(src_lightcurve['lcBins'], src_lightcurve['ctrate'], xerr=src_lightcurve['lcBinsRange'] / 2,
                    yerr=src_lightcurve['ctrateErr'], fmt='ok')

    # Plotting the flare interval atop the source light curve
    mask = np.where(src_lightcurve.lcBins.isin(flare_lightcurve.lcBins), True, False)
    flare_src_lightcurve = src_lightcurve[mask]

    axs[1].errorbar(flare_src_lightcurve['lcBins'], flare_src_lightcurve['ctrate'],
                    xerr=flare_src_lightcurve['lcBinsRange'] / 2, yerr=flare_src_lightcurve['ctrateErr'], fmt='or')
    axs[1].set_title(r'$\,\mathrm{Source\,energy\,range}$', fontsize=12)

    for axis in ['top', 'bottom', 'left', 'right']:
        axs[0].spines[axis].set_linewidth(1.5)
        axs[0].tick_params(width=1.5)
        axs[1].spines[axis].set_linewidth(1.5)
        axs[1].tick_params(width=1.5)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.1)

    outPlot = outputFile + '.pdf'
    fig.savefig(outPlot, format='pdf', dpi=200)
    plt.close(fig)

    return


def main():
    parser = argparse.ArgumentParser(
        description="Creating time intervals that are likely associated with high-energy background flares - built for "
                    "NICER")
    parser.add_argument("evtFile", help="Fits event file", type=str)
    parser.add_argument("-elb", "--eneLow_back", help="Low energy filter in event file, default=12",
                        type=float, default=12)
    parser.add_argument("-ehb", "--eneHigh_back", help="High energy filter in event file, default=15",
                        type=float, default=15)
    parser.add_argument("-tb", "--timebin", help="Time resolution of light curve in seconds",
                        type=float, default=5.)
    parser.add_argument("-lt", "--lcthresh", help="exclude time bins with exposure less than lcthresh*timebin",
                        type=float, default=0.1)
    parser.add_argument("-pb", "--problim", help="Probability mass function limit for burst search, "
                                                 "default = 0.01", type=float, default=0.01)
    parser.add_argument("-of", "--outputFile",
                        help="name of output .pdf light curve showing flagged bins", type=str, default='flagged_flares')
    parser.add_argument("-cg", "--creategti", help="Name of gti fits file compatible with "
                                                   "NICERDAS, default = None", type=str, default=None)
    parser.add_argument("-els", "--eneLow_src", help="low energy cutoff for source (for plotting "
                                                     "purposes only), default=1", type=float, default=1)
    parser.add_argument("-ehs", "--eneHigh_src", help="high energy cutoff for source (for plotting "
                                                      "purposes only), default=5", type=float, default=5)
    args = parser.parse_args()

    flagbackgrflares(args.evtFile, args.eneLow_back, args.eneHigh_back, args.timebin, args.lcthresh, args.problim,
                     args.outputFile, args.creategti, args.eneLow_src, args.eneHigh_src)


if __name__ == '__main__':
    main()
