"""
# correctfpmsel.py
"""

import numpy as np
import pandas as pd

# Custom modules
from nicerutil.eventfile import EvtFileOps

import sys

sys.dont_write_bytecode = True


def correctrateforfpmsel(eventfile, binnedlightcurve):
    # Reading the selected FPM per 1-second interval from event file
    EF = EvtFileOps(eventfile)
    _, FPMSEL_table_condensed = EF.read_fpmsel()

    lcBins = binnedlightcurve["lcBins"].to_numpy()
    lcBinsRange = binnedlightcurve["lcBinsRange"].to_numpy()
    ctrate = binnedlightcurve["ctrate"].to_numpy()
    ctrateErr = binnedlightcurve["ctrateErr"].to_numpy()
    ctsbin = binnedlightcurve["ctsbin"].to_numpy()

    for ii, timeofbin in enumerate(lcBins):
        bin_start = timeofbin - lcBinsRange[ii] / 2
        bin_end = timeofbin + lcBinsRange[ii] / 2

        # filtering FPM selection table to within each ToA
        # The -1 is to ensure we register the number of detectors at the start of each timebin
        mkf_toa_filtered = FPMSEL_table_condensed.loc[
            ((FPMSEL_table_condensed['TIME'] > (bin_start - 1)) & (FPMSEL_table_condensed['TIME'] <= bin_end))]
        mkf_toa_filtered = mkf_toa_filtered.reset_index(drop=True)

        # Here we measure the number of selected detectors, naively it is the sum of FPM_SEL within each time bin
        # However, the first and last bin may have covered a fraction of 1 second, hence we normalize accordingly
        nbr_sel_det = ((((1 - (bin_start -
                               mkf_toa_filtered['TIME'].head(1).to_numpy())) *
                         mkf_toa_filtered['TOTFPMSEL'].head(1).to_numpy()) +
                        ((bin_end - mkf_toa_filtered['TIME'].tail(1).to_numpy()) *
                         mkf_toa_filtered['TOTFPMSEL'].tail(1).to_numpy())) +
                       np.sum(mkf_toa_filtered['TOTFPMSEL'].loc[1:len(mkf_toa_filtered) - 2]))

        # If number of selected detectors is 0 during a certain supposed good time interval, move on
        # This is okay here, but keep in mind it does mess up your good-time-intervals
        # if nbr_sel_det < 1:
        #    continue

        # Here we measure the number of expected detectors, naively it is 52 * number of time bins
        # However, the first and last bins may have covered a fraction of 1 second, hence we normalize accordingly
        exp_nbr_det = ((1 - (bin_start - mkf_toa_filtered['TIME'].head(1).to_numpy())) * 52 +
                       (bin_end - mkf_toa_filtered['TIME'].tail(1).to_numpy()) * 52 +
                       52 * len(mkf_toa_filtered.loc[1:len(mkf_toa_filtered) - 2]))

        # Summing all detectors - old way and wrong at edges
        # nbr_sel_det = np.sum(mkf_toa_filtered['TOTFPMSEL'])
        # total number of detectors if 52 were operating
        # exp_nbr_det = 52 * len(mkf_toa_filtered['TIME'])

        # Measuring correction factor
        correction_factor = exp_nbr_det / nbr_sel_det
        # Changing rate values with corrected ones
        ctrate[ii] *= correction_factor
        ctrateErr[ii] *= correction_factor
        ctsbin[ii] = (ctsbin[ii] * correction_factor).astype(int)

    corrbinnedLC = {'lcBins': lcBins, 'lcBinsRange': lcBinsRange, 'ctrate': ctrate, 'ctrateErr': ctrateErr,
                    'ctsbin': ctsbin}
    corrbinnedLC = pd.DataFrame.from_dict(corrbinnedLC)

    return corrbinnedLC
