"""
Module to help clean NICER data from high energy flares (known as PREL, see
https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/flares/). Given
a NICER event file and a mkf file, this module will create a light curve in
the energy range 12-15 keV with a time-bin of 5 seconds (both of which optional
inputs), where NICER effective area is practically 0. It will then search for bins
with a number of counts that is too large to be considered random Poisson fluctuation
and flagit as a potential background flare. It will also produce a xselect and
NICERDAS-compatible GTI that corresponds to the flagged events. A plot of the
flagged events is also produced. Other optional parameters are allowed.

Warning: make sure you have initialized heasoft before running this module
"""
import matplotlib.pyplot as plt
import sys
import argparse
import os
import numpy as np
from astropy.io import fits
import pandas as pd

# Custom modules
from nicerutil.eventfile import EvtFileOps
from nicerutil.lightcurve import lightcurve, plotlightcurve
from nicerutil.bursts import burstsearch, mergesamebursts
from nicerutil.correctrateforfpmsel import correctrateforfpmsel
from nicerutil.nicermkf import MkfFileOps, readmkffile
from nicerutil.nicerutil_logging import get_logger
from nicerutil.build_gti import write_gti

sys.dont_write_bytecode = True

# Log config
############
logger = get_logger(__name__)


def flagbackgrflares(eventfile, mkffile, eneLow_back=12, eneHigh_back=15, timebin=5., lcthresh=0.1, probLim=0.01,
                     eneLow_src=1, eneHigh_src=5, outputFile='flagged_flares'):
    """
    flagbackgrflares creates a light curves, and flags bins that are considered to be outliers
    according to a Poisson mass function
    :param eventfile: name of the fits event file
    :type eventfile: str
    :param mkffile: name of the nicer MKF file
    :type mkffile: str
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
    :param eneLow_src: low energy cutoff for source (for plotting purposes only default = 1 keV)
    :type eneLow_src: float
    :param eneHigh_src: high energy cutoff for source (for plotting purposes only default = 5 keV)
    :type eneHigh_src: float
    :param outputFile: Name of file to save diagnostic .pdf plot, gti _gti.txt and _gti.fits files (default = "flagged_flares")
    :type outputFile: str
    :return: flares, distinct_flares (all falgged bins and a merged list if bins are consecutive)
    :rtype: list
    """

    logger.info('\n Running flagbackgrflares module with input parameters :'
                '\n eventfile : ' + str(eventfile) +
                '\n mkffile : ' + str(mkffile) +
                '\n eneLow_back (keV) : ' + str(eneLow_back) +
                '\n eneHigh_back (keV) : ' + str(eneHigh_back) +
                '\n timebin (s) : ' + str(timebin) +
                '\n lcthresh : ' + str(lcthresh) +
                '\n probLim : ' + str(probLim) +
                '\n eneLow_src (keV) : ' + str(eneLow_src) +
                '\n eneHigh_src (keV) : ' + str(eneHigh_src) +
                '\n outputFile : ' + str(outputFile) + '\n')

    # Reading data
    EF = EvtFileOps(eventfile)
    evtFileKeyWords, gtiList = EF.readGTI()
    full_exposure = evtFileKeyWords['ONTIME']

    # Guard against non-existent GTI in event file (but with a valid ONTIME keyword!)
    # GTI with the same START and STOP times (but also with a valid ONTIME keyword!)
    # 0 exposure event files
    # These instances become relevant for orbit day data after optical light leak
    if gtiList.size == 0:
        logger.warning('No valid GTI in event file {} - skipping'.format(eventfile))
        return None, None
    elif np.sum(gtiList[:, -1] - gtiList[:, 0]) == 0:
        logger.warning('GTI sum to 0 in event file {} - skipping'.format(eventfile))
        return None, None
    elif full_exposure == 0:
        logger.warning('Exposure of event file {} is 0 - skipping'.format(eventfile))
        return None, None

    # Dealing with the background
    #############################
    # Energy filtering, reading TIME column as numpy array
    dataTP_eneFlt = EF.filtenergy(eneLow=eneLow_back, eneHigh=eneHigh_back)
    TIME_back = dataTP_eneFlt['TIME'].to_numpy()

    # Create light curve
    binnedLC, _ = lightcurve(TIME_back, gtiList, timebin=timebin, lcthresh=lcthresh)
    binnedLC_corr = correctrateforfpmsel(eventfile, binnedLC)
    plotlightcurve(binnedLC_corr, outputFile=f"{outputFile}_bg_lc")

    # FLare search
    burstsearch_result = burstsearch(binnedLC_corr, probLim=probLim, outputfile=outputFile)
    flares = burstsearch_result['bin_bursts']

    # Dealing with the source
    #########################
    # Creating a light curve within source energy range
    dataTP_eneFlt_src = EF.filtenergy(eneLow=eneLow_src, eneHigh=eneHigh_src)
    TIME_src = dataTP_eneFlt_src['TIME'].to_numpy()
    binnedLC_src, _ = lightcurve(TIME_src, gtiList, timebin=timebin, lcthresh=lcthresh)
    binnedLC_src_corr = correctrateforfpmsel(eventfile, binnedLC_src)
    plotlightcurve(binnedLC_src_corr, outputFile=f"{outputFile}_src_lc")

    # Create .txt and .fits GTI files
    #################################
    if not flares.empty:
        # Merge flares separated by timebin
        distinct_flares = mergesamebursts(flares, tsameburst=1.5 * timebin)
        flare_start_stop_df = flares_dict_to_intervals(distinct_flares)
        write_gti(flare_start_stop_df, mkffile, outputFile)

        # Create a simple diagnostic plot of the flares
        ###############################################
        # Reading mkf data for full GTI of observation
        mkf_table = readmkffile(mkffile, over='MPU_OVER_COUNT')
        timefiltered_mkf = MkfFileOps(mkf_table).timefiltermkf(gtiList)
        mkf_over_cor = timefiltered_mkf[['tNICERmkf', 'OVER_ONLY_COUNT', 'corSax']].copy()
        # Reading mkf data of clean GTI (flare excluded)
        hdulist_gti = fits.open(outputFile + '_gti.fits')
        gtiList_clean = (np.vstack((hdulist_gti["STDGTI"].data.field("START"),
                                    hdulist_gti["STDGTI"].data.field("STOP")))).T
        timefiltered_mkf_clean = MkfFileOps(mkf_table).timefiltermkf(gtiList_clean)
        mkf_over_cor_clean = timefiltered_mkf_clean[['tNICERmkf', 'OVER_ONLY_COUNT', 'corSax']].copy()
        # Creating the plot
        plot_flare_diagnostics(binnedLC_corr, flares, binnedLC_src_corr, mkf_over_cor, mkf_over_cor_clean,
                               outputFile=outputFile)

        # Calculate exposures
        clean_exposure = np.sum(burstsearch_result['bins_nobursts']['lcBinsRange'])
        fractional_loss = (full_exposure - clean_exposure) / full_exposure

        # Adding some info to logging file about the cleaning process
        logger.info('\n Cleaning statistics:')
        logger.info('\n Average rate of background light curve before any cleaning : {}'.format(
            burstsearch_result['average_rate']))
        logger.info('\n Average rate of background light curve after cleaning : {}'.format(
            burstsearch_result['average_rate_nobursts']))
        logger.info('\n Full exposure : {}\n Clean exposure : {}\n Fractional loss : {}'.format(
            full_exposure, clean_exposure, fractional_loss))

    else:
        distinct_flares = None
        logger.info('\n No flares found.')

    return flares, distinct_flares


def flares_dict_to_intervals(distinct_flares: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """
    Convert a dict of light curve DataFrames into a (tstart, tstop) DataFrame.
    Each DF must have ['lcBins','lcBinsRange'] and be ordered by time - as in the result of
    the lightcurve() function in lightcurve.py module
    This is really only used for filtering out flaring time intervals
    :param distinct_flares: dict of flare light-curve dataframes
    :type distinct_flares: dict[str, pd]
    :return: DataFrame of tstart and tstop times of the flare intervals
    :rtype: pandas.DataFrame
    """
    rows = []
    for key, burst_df in distinct_flares.items():
        start = burst_df['lcBins'].iloc[0] - (burst_df['lcBinsRange'].iloc[0] / 2.0)
        stop = burst_df['lcBins'].iloc[-1] + (burst_df['lcBinsRange'].iloc[-1] / 2.0)
        rows.append((start, stop, key))
    return pd.DataFrame(rows, columns=['tstart', 'tstop', 'label'])[['tstart', 'tstop']]


def plot_flare_diagnostics(back_lightcurve, flare_lightcurve, src_lightcurve, mkf_over_cor, mkf_over_cor_clean,
                           outputFile='flare_diagnostics'):
    """
   Creating a plot of the flagged elevated, background events. Mainly called by flagbackgrflares function and
   not meant to be a standalone function.
   """

    fig, axs = plt.subplots(4, 1, figsize=(10, 14), dpi=200, facecolor='w', edgecolor='k', sharex=True,
                            sharey=False, gridspec_kw={'width_ratios': [1], 'height_ratios': [1, 1, 1, 1]})
    axs = axs.ravel()

    # Plotting the "background" light curve,
    # i.e., out-of-bounds events, and those flagged as flares
    #########################################################
    axs[0].tick_params(axis='both', labelsize=14)
    axs[0].set_ylabel(r'$\,\mathrm{Rate\,(counts/\,s)}$', fontsize=14)
    axs[0].ticklabel_format(style='plain', axis='y', scilimits=(0, 0))
    axs[0].xaxis.offsetText.set_fontsize(14)
    axs[0].yaxis.offsetText.set_fontsize(14)
    # Plotting all out-of-band events
    axs[0].errorbar(back_lightcurve['lcBins'], back_lightcurve['ctrate'], xerr=back_lightcurve['lcBinsRange'] / 2,
                    yerr=back_lightcurve['ctrateErr'], fmt='ok', label='Good GTIs')
    # Plotting the flare interval
    axs[0].errorbar(flare_lightcurve['lcBins'], flare_lightcurve['ctrate'], xerr=flare_lightcurve['lcBinsRange'] / 2,
                    yerr=flare_lightcurve['ctrateErr'], fmt='or', label='Flare GTIs')
    axs[0].set_title(r'$\,\mathrm{Background\,energy\,range}$', fontsize=14)

    # Plotting the "source" light curve,
    # i.e., in-bound events, and highlight the bins that correspond to the flares
    #############################################################################
    axs[1].tick_params(axis='both', labelsize=14)
    axs[1].set_ylabel(r'$\,\mathrm{Rate\,(counts/\,s)}$', fontsize=14)
    axs[1].ticklabel_format(style='plain', axis='y', scilimits=(0, 0))
    axs[1].xaxis.offsetText.set_fontsize(14)
    axs[1].yaxis.offsetText.set_fontsize(14)
    # Plotting all in-band events
    axs[1].errorbar(src_lightcurve['lcBins'], src_lightcurve['ctrate'], xerr=src_lightcurve['lcBinsRange'] / 2,
                    yerr=src_lightcurve['ctrateErr'], fmt='ok')

    # Plotting the flare interval atop the source light curve
    mask = np.where(src_lightcurve.lcBins.isin(flare_lightcurve.lcBins), True, False)
    flare_src_lightcurve = src_lightcurve[mask]

    axs[1].errorbar(flare_src_lightcurve['lcBins'], flare_src_lightcurve['ctrate'],
                    xerr=flare_src_lightcurve['lcBinsRange'] / 2, yerr=flare_src_lightcurve['ctrateErr'], fmt='or')
    axs[1].set_title(r'$\,\mathrm{Source\,energy\,range}$', fontsize=14)

    # Plotting the MKF relevant columns
    #############################################################################
    # TOT_OVER_COUNT
    axs[2].tick_params(axis='both', labelsize=14)
    axs[2].set_ylabel(r'$\,\mathrm{Rate\,(counts/\,s)}$', fontsize=14)
    axs[2].ticklabel_format(style='plain', axis='y', scilimits=(0, 0))
    axs[2].xaxis.offsetText.set_fontsize(14)
    axs[2].yaxis.offsetText.set_fontsize(14)
    # Plotting all TOT_OVER_COUNT
    axs[2].errorbar(mkf_over_cor['tNICERmkf'], mkf_over_cor['OVER_ONLY_COUNT'], xerr=1 / 2, fmt='ok')
    # Plotting the flare interval atop
    mask = np.where(~mkf_over_cor.tNICERmkf.isin(mkf_over_cor_clean.tNICERmkf), True, False)
    mkf_over_cor_flare = mkf_over_cor[mask]
    axs[2].errorbar(mkf_over_cor_flare['tNICERmkf'], mkf_over_cor_flare['OVER_ONLY_COUNT'], xerr=1 / 2, fmt='or')
    axs[2].set_title(r'$\,\mathrm{TOT\_OVER\_COUNT}$', fontsize=14)

    # COR_SAX
    axs[3].tick_params(axis='both', labelsize=14)
    axs[3].set_xlabel(r'$\,\mathrm{Time\,(MET,\,s)}$', fontsize=14)
    axs[3].set_ylabel(r'$\,\mathrm{GeV/\,c}$', fontsize=14)
    axs[3].ticklabel_format(style='plain', axis='y', scilimits=(0, 0))
    axs[3].xaxis.offsetText.set_fontsize(14)
    axs[3].yaxis.offsetText.set_fontsize(14)
    # Plotting all COR_SAX
    axs[3].errorbar(mkf_over_cor['tNICERmkf'], mkf_over_cor['corSax'], xerr=1 / 2, fmt='ok')
    # Plotting the flare interval atop
    axs[3].errorbar(mkf_over_cor_flare['tNICERmkf'], mkf_over_cor_flare['corSax'], xerr=1 / 2, fmt='or')
    axs[3].hlines(y=1.5, xmin=mkf_over_cor['tNICERmkf'].min(), xmax=mkf_over_cor['tNICERmkf'].max(),
                  linewidth=2, color='r', linestyles='dashed')
    axs[3].set_title(r'$\,\mathrm{COR\_SAX}$', fontsize=14)

    #############################################################################
    for axis in ['top', 'bottom', 'left', 'right']:
        axs[0].spines[axis].set_linewidth(1.5)
        axs[0].tick_params(width=1.5)
        axs[1].spines[axis].set_linewidth(1.5)
        axs[1].tick_params(width=1.5)
        axs[2].spines[axis].set_linewidth(1.5)
        axs[2].tick_params(width=1.5)
        axs[3].spines[axis].set_linewidth(1.5)
        axs[3].tick_params(width=1.5)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.15)

    outPlot = outputFile + '.pdf'
    fig.savefig(outPlot, format='pdf', dpi=200)
    plt.close(fig)

    return


def main():
    parser = argparse.ArgumentParser(
        description="Creating time intervals that are likely associated with high-energy background flares - built for "
                    "NICER")
    parser.add_argument("evtFile", help="Fits event file", type=str)
    parser.add_argument("mkffile", help="Fits nicer MKF file", type=str)
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
    parser.add_argument("-els", "--eneLow_src", help="low energy cutoff for source (for plotting "
                                                     "purposes only), default=1", type=float, default=1)
    parser.add_argument("-ehs", "--eneHigh_src", help="high energy cutoff for source (for plotting "
                                                      "purposes only), default=5", type=float, default=5)
    parser.add_argument("-of", "--outputFile",
                        help="name of output .pdf light curve showing flagged bins", type=str, default='flagged_flares')
    args = parser.parse_args()

    flagbackgrflares(args.evtFile, args.mkffile, args.eneLow_back, args.eneHigh_back, args.timebin, args.lcthresh,
                     args.problim, args.eneLow_src, args.eneHigh_src, args.outputFile)


if __name__ == '__main__':
    main()
