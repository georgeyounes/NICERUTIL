"""
correct_gtis.py is a module that filters out any GTIs in a level 2 event file of minimum length with 0 events.
For reasons which I do not know, there exists some instances when a valid GTI of considerable length is saved,
yet the GTI does not seem to have any registered events. These are removed. Interestingly, if the user runs nicerl3-lc
or nicerl3-spect, these GTIs are indeed filtered out. Craig Markwardt, the NICERDAS author, is aware of this bug
and a fix should be implemented at a later version of the software.

Update (11-10-2025): while analyzing observation ID 6533042408 (thanks to Joanna Berteau for pointing it out), we
noticed that, in a level 2 event file (e.g., the default output of nicerl2), several ks separates the START of the
first GTI and the first event in the EVENTS table (when considering orbit-day data). Likewise for the separation
between the last recorded EVENT and the END of the last GTI. This is an edge case of the above and may cause several
problems, one of it is rate calculation, and the other is barycentric correction, in the case where the GTI times
land outside an orbit file timestamp - which is the case in this specific observation.
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


def correct_gtis(eventfile, min_length=5, new_eventfile=None):
    """
    Function that filters out a row in a GTI extension in a nicer event file if the exposure is of length
    min_len (seconds) within which no events are registered
    :param eventfile: name of the level 2 fits event file
    :type eventfile: str
    :param min_length: duration in seconds of GTI row (default = 5 seconds)
    :type min_length: float
    :param new_eventfile: Name of new event file (default = None, ie, update in-place)
    :type new_eventfile: str | None
    """
    EF = EvtFileOps(eventfile)
    evtFileKeyWords, gti = EF.readGTI()

    # Get full list of events
    df_events_T_PI = EF.filtenergy(eneLow=0.1, eneHigh=20)  # We consider the full energy range of NICER
    event_times = df_events_T_PI['TIME']
    if event_times.empty:
        logger.info("Event file has 0 events - returning the original input event file")
        return eventfile

    # Initilize bad GTIs
    bad_rows_indices = []

    # let's take care of GTIs that are larger than 5 seconds in length, and yet do not have any registered events
    # they will be removed
    for idx, row in enumerate(gti):
        if row[1] - row[0] < min_length:
            logger.info(f"Length of GTI row {idx} in observation ID {evtFileKeyWords['OBS_ID']} is less than desired "
                        f"{min_length} seconds. Keeping this GTI row without checking nbr of EVENTS.")
            continue
        else:
            df_events_T_PI_pergti = df_events_T_PI[(df_events_T_PI['TIME'] >= row[0]) &
                                                   (df_events_T_PI['TIME'] <= row[1])]
            if df_events_T_PI_pergti.empty:
                bad_rows_indices.append(idx)
                logger.info(f"No events within row {idx} of observation ID {evtFileKeyWords['OBS_ID']}\n."
                            f"GTI row START, STOP is {row[0]}, {row[1]}; total exposure is {row[1] - row[0]} s."
                            f"This GTI will be eliminated.")

    # let's also take care of the case where GTI START and/or GTI END is more than 5-seconds from the first and/or last
    # event in the EVENTS table, respectively - remove those rows
    # Taking care of start
    event_start = event_times.iloc[0] - 5  # we are considering a 5-second buffer
    near_start_idx = np.where(gti[:, 1] < (event_start - 5))[0]  # Flag GTIs whose end is < (event_start - 5)
    bad_rows_indices.append(near_start_idx)
    # Taking care of end
    event_end = event_times.iloc[-1] + 5  # we are considering a 5-second buffer
    after_end_idx = np.where(gti[:, 0] > (event_end + 5))[0]  # Flag GTIs whose start is > (event_end + 5)
    bad_rows_indices.append(after_end_idx)

    # Merge, clean, and sort bad-row indices
    bad_lists = []
    for x in bad_rows_indices:
        x = np.atleast_1d(x)
        if x.size:
            bad_lists.append(x.astype(int))
    if bad_lists:
        bad_rows_indices = np.sort(np.unique(np.concatenate(bad_lists)))
    else:
        logger.info(f"No bad GTIs exist - returning the original input event file")
        return eventfile

    # Build filtered GTI table from file
    # Original gti table
    gti_table = Table.read(eventfile, format='fits', hdu='GTI')
    # filter the GTI table according to the "bad" GTI row indices
    mask = np.ones(len(gti_table), dtype=bool)
    mask[bad_rows_indices] = False
    gti_table_flt = gti_table[mask]

    # If everything was filtered, bail with a clear error
    if len(gti_table_flt) == 0:
        raise RuntimeError(
            f"All GTI rows were filtered for {evtFileKeyWords.get('OBS_ID', '<unknown>')} — no valid GTIs remain."
        )

    # Update the START and STOP columns of the GTI table
    gti_table_flt['START'] = gti_table_flt['START'].astype('float64')
    gti_table_flt['STOP'] = gti_table_flt['STOP'].astype('float64')

    # New start and end times from filtered GTI
    cleaned_gti_start = float(np.min(gti_table_flt['START']))
    cleaned_gti_end = float(np.max(gti_table_flt['STOP']))

    # (Re)build GTI HDU and set its header keywords properly
    gti_hdu = fits.table_to_hdu(gti_table_flt)
    gti_hdu.name = 'GTI'
    gti_hdu.header['TSTART'] = cleaned_gti_start
    gti_hdu.header['TSTOP'] = cleaned_gti_end

    # Write a NEW event file with updated GTI, and sync headers in EVENTS/FPM_SEL
    if new_eventfile is None:
        # Update in place
        with fits.open(eventfile, mode='update') as hdul:
            if 'GTI' in hdul:
                orig_extver = hdul['GTI'].header.get('EXTVER')
                hdul['GTI'] = gti_hdu
                if orig_extver is not None:
                    hdul['GTI'].header['EXTVER'] = orig_extver
            else:
                hdul.append(gti_hdu)

            # Sync TSTART/TSTOP in EVENTS and FPM_SEL
            for hdu_name in ('EVENTS', 'FPM_SEL'):
                if hdu_name in hdul:
                    hdr = hdul[hdu_name].header
                    hdr['TSTART'] = cleaned_gti_start
                    hdr['TSTOP'] = cleaned_gti_end

            hdul.flush()
        outname = eventfile

    else:
        # Write to a new file
        with fits.open(eventfile) as hdul:
            if 'GTI' in hdul:
                orig_extver = hdul['GTI'].header.get('EXTVER')
                hdul['GTI'] = gti_hdu
                if orig_extver is not None:
                    hdul['GTI'].header['EXTVER'] = orig_extver
            else:
                hdul.append(gti_hdu)

            for hdu_name in ('EVENTS', 'FPM_SEL'):
                if hdu_name in hdul:
                    hdr = hdul[hdu_name].header
                    hdr['TSTART'] = cleaned_gti_start
                    hdr['TSTOP'] = cleaned_gti_end

            hdul.writeto(new_eventfile, overwrite=True)
        outname = new_eventfile

    return outname


def update_gti(GTI, tstart, tend):
    # Restructuring GTI if tstart and tend are provided
    if tstart is None and tend is None:
        logger.info('\n No tstart or tend provided. Using full GTI table instead\n')
        return GTI

    elif tstart is not None and tend is None:
        logger.info('\n Fixing GTI start time to match user defined tstart')
        includegti_idx = (GTI[:, 1] > tstart)
        GTI = GTI[includegti_idx]
        if tstart > GTI[0, 0]:
            GTI[0, 0] = tstart
        return GTI

    elif tstart is None and tend is not None:
        logger.info('\n Fixing GTI end time to match user defined tend')
        includegti_idx = (GTI[:, 0] < tend)
        GTI = GTI[includegti_idx]
        if tend < GTI[-1, -1]:
            GTI[-1, -1] = tend
        return GTI

    elif tstart is not None and tend is not None:
        logger.info('\n Fixing GTI start and end time to match user defined tstart and tend')
        # Fixing tstart
        includegti_idx = (GTI[:, 1] > tstart)
        GTI = GTI[includegti_idx]
        if tstart > GTI[0, 0]:
            GTI[0, 0] = tstart
        # Fixing tend
        includegti_idx = (GTI[:, 0] < tend)
        GTI = GTI[includegti_idx]
        if tend < GTI[-1, -1]:
            GTI[-1, -1] = tend
        return GTI


def main():
    """
    Run correct_gtis.py from command line
    """
    parser = argparse.ArgumentParser(description="Filter GTI table to exclude rows with no registered events")
    parser.add_argument("eventfile", help="Name of NICER event fits file", type=str)
    parser.add_argument("-ml", "--min_length", help="Minimum duration of GTI row for filtering",
                        type=float, default=5)
    parser.add_argument("-op", "--output",
                        help="Name of output FITS file (default = None, i.e., overwrite input)",
                        default=None)
    args = parser.parse_args()

    correct_gtis(args.eventfile, args.min_length, args.output)


if __name__ == '__main__':
    main()
