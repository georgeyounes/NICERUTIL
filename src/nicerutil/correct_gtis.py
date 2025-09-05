"""
correct_gtis.py is a module that filters out any GTIs in a level 2 event file of minimum length with 0 events.
For reasons which I do not know, there exists some instances when a valid GTI of considerable length is saved,
yet the GTI does not seem to have any registered events. These are removed. Interestingly, if the user runs nicerl3-lc
or nicerl3-spect, these GTIs are indeed filtered out. Craig Markwardt, the NICERDAS author, is aware of this bug
and a fix should be implemented at a later version of the software.
"""

import argparse
import numpy as np
from astropy.io import fits
from astropy.table import Table

import sys
from nicerutil.nicerutil_logging import get_logger
from nicerutil.eventfile import EvtFileOps
sys.dont_write_bytecode = True

# Log config
############
logger = get_logger(__name__)


def correct_gtis(eventfile, min_length=5):
    """
    Function that filters out a row in a GTI extension in a nicer event file if the exposure is of length
    min_len (seconds) within which no events are registered
    :param eventfile: name of the level 2 fits event file
    :type eventfile: str
    :param min_length: duration in seconds of GTI row (default = 5 seconds)
    :type min_length: float
    """
    EF = EvtFileOps(eventfile)
    evtFileKeyWords, gti = EF.readGTI()

    # Get full list of events
    df_events_T_PI = EF.filtenergy(eneLow=0.1, eneHigh=20)  # We consider the full energy range of NICER

    bad_rows_indices = []
    for idx, row in enumerate(gti):
        if row[1] - row[0] < min_length:
            logger.info(f"Length of GTI row {idx} in observation ID {evtFileKeyWords['OBS_ID']} is less than desired "
                        f"{min_length} seconds. Keeping this GTI row without checking EVENTS.")
            continue
        else:
            df_events_T_PI_pergti = df_events_T_PI[(df_events_T_PI['TIME'] >= row[0]) &
                                                   (df_events_T_PI['TIME'] <= row[1])]
            if df_events_T_PI_pergti.empty:
                bad_rows_indices.append(idx)
                logger.info(f"No events within row {idx} of observation ID {evtFileKeyWords['OBS_ID']}\n."
                            f"GTI row START, STOP is {row[0]}, {row[1]}; total exposure is {row[1] - row[0]} s.")

    bad_rows_indices = np.array(bad_rows_indices, dtype=int)

    # Let's read the GTI header
    hdulist_gti = fits.open(eventfile)
    gti_hdr = hdulist_gti['GTI'].header
    hdulist_gti.close()

    # Filtering out the GTI table - this is done through TABLE class of astropy.table
    gti_table = Table.read(eventfile, format='fits', hdu='GTI')
    mask = np.ones(len(gti_table), dtype=bool)
    mask[bad_rows_indices] = False
    gti_table_flt = gti_table[mask]

    # Changing types to match those of the original
    gti_table_flt['START'].info.dtype = '>f8'
    gti_table_flt['STOP'].info.dtype = '>f8'

    # Updating event file
    newhdulEFPH = fits.BinTableHDU(data=gti_table_flt, header=gti_hdr, name='GTI')
    fits.update(eventfile, newhdulEFPH.data, newhdulEFPH.header, 'GTI')

    return eventfile


def main():
    """
    Run correct_gtis.py from command line
    """
    parser = argparse.ArgumentParser(description="Filter GTI table to exclude rows with no registered events")
    parser.add_argument("eventfile", help="Name of NICER event fits file", type=str)
    parser.add_argument("-ml", "--min_length", help="Minimum duration of GTI row for filtering",
                        type=float, default=5)
    args = parser.parse_args()

    correct_gtis(args.eventfile, args.min_length)


if __name__ == '__main__':
    main()
