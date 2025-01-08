"""
correctfpmsel.py is a module that filters the FPM_SEL extension in a nicer event
file to match its GTI start and stop times. This is done for proper barycentering
correction of the FMP_SEL table, which on few occasions had time stamps that were
outside the orbit file. Filtering within GTI times ensures that all FPM_SEL time
stamps fall within the observation orbit file.
"""

import argparse
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column

import sys

sys.dont_write_bytecode = True


def correctfpmsel(eventfile):
    """
    Function that filters the FPM_SEL extension in a nicer event file to match its GTI start and stop times
    :param eventfile: name of the fits event file
    :type eventfile: str
    """
    hdulist = fits.open(eventfile)

    # Getting TSTART and TSTOP from EVENTS table
    TSTART = hdulist['EVENTS'].header['TSTART']
    TSTOP = hdulist['EVENTS'].header['TSTOP']
    # #xtracting GTI table
    ST_GTI = hdulist["GTI"].data['START']
    ET_GTI = hdulist["GTI"].data['STOP']
    gtiTable = (np.vstack((ST_GTI, ET_GTI))).T

    # Reading FPM_SEL extension
    hdulist_fpmsel = fits.open(eventfile)
    fpmsel_ext = hdulist_fpmsel['FPM_SEL']
    fpmsel_hdr = hdulist_fpmsel['FPM_SEL'].header

    # Changing TSTART and TSTOP in FPMSEL table to match EVENTS table
    fpmsel_hdr['TSTART'] = TSTART
    fpmsel_hdr['TSTOP'] = TSTOP

    # Getting indices of good FPM_SEL times, i.e., match the cleaned event file gtis
    indices_goodfpmsel = []
    for jj in range(len(gtiTable[:, 0])):
        indices_goodfpmsel_pergti = np.where((fpmsel_ext.data['TIME'] >= np.floor(gtiTable[jj, 0])) & (
                    fpmsel_ext.data['TIME'] <= np.ceil(gtiTable[jj, 1])))
        indices_goodfpmsel = np.append(indices_goodfpmsel, indices_goodfpmsel_pergti)
    indices_goodfpmsel = indices_goodfpmsel.astype(int)

    # Filtering out the FPM_SEL table - this is done through TABLE class of astropy.table
    fpmsel_table = Table.read(eventfile, format='fits', hdu='FPM_SEL')
    fpmsel_table_flt = fpmsel_table[indices_goodfpmsel]

    # Changing types to match those of the original
    fpmsel_table_flt['TIME'].info.dtype = '>f8'
    fpmsel_table_flt['FPM_ON'] = fpmsel_table_flt['FPM_ON'].astype(
        'B')  # I have a vague idea why the above does not work on FPM_ON and FPM_SEL but this works
    fpmsel_table_flt['FPM_SEL'] = fpmsel_table_flt['FPM_SEL'].astype(np.int16)

    # Updating event file
    newhdulEFPH = fits.BinTableHDU(data=fpmsel_table_flt, header=fpmsel_hdr, name='FPM_SEL')
    fits.update(eventfile, newhdulEFPH.data, newhdulEFPH.header, 'FPM_SEL')

    return


def main():
    """
    Run correctfpmsel.py from command line
    """
    parser = argparse.ArgumentParser(description="Time-filter FPM_SEL table in a NICER event file to match GTI")
    parser.add_argument("eventfile", help="Name of NICER event fits file", type=str)
    args = parser.parse_args()

    correctfpmsel(args.eventfile)


if __name__ == '__main__':
    main()
