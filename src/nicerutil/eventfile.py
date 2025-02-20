"""
eventfile.py is a module to perform simple operations o nicer event file,
reading useful header words, filter for energy etc.
This is a simpler version of the eventfile.py package in CRIMP
"""

import sys
import argparse
import numpy as np
import pandas as pd
import logging

from astropy.table import Table, Column
from astropy.io import fits

sys.dont_write_bytecode = True

# Log config
############
logger = get_logger(__name__)

class EvtFileOps:
    """
        A class to operate on a fits event file

        Attributes
        ----------
        evtFile : str
            name of the fits event file

        Methods
        -------
        readEF(): reads essential keywords from an event file
        read_fpmsel(): reads an FPM_SEL table from an event file
        readGTI(): reads GTI table from an event
        filtEneEF(): filters the event list accoring to energy (in keV)
        """

    def __init__(self, evtFile: str):
        """
        Constructs the necessary attribute

        :param evtFile: name of the fits event file
        :type evtFile: str
        """
        self.evtFile = evtFile

    #################################################################
    def readEF(self):  # Reading fits event file from X-ray satellites
        """
        Reads essential keywords from an event file
        :return: evtFileKeyWords - dictionary of essential keywords
        :rtype: dict
        """
        # Opening the fits file
        hdulist = fits.open(self.evtFile)

        # Reading some essential keywords
        TELESCOPE = hdulist['EVENTS'].header['TELESCOP']
        INSTRUME = hdulist['EVENTS'].header['INSTRUME']
        TSTART = hdulist['EVENTS'].header['TSTART']
        TSTOP = hdulist['EVENTS'].header['TSTOP']
        TIMESYS = hdulist['EVENTS'].header['TIMESYS']
        DATEOBS = hdulist['EVENTS'].header['DATE-OBS']
        TIMEZERO = hdulist['EVENTS'].header['TIMEZERO']
        OBS_ID = hdulist['EVENTS'].header['OBS_ID']
        LIVETIME = hdulist['EVENTS'].header['LIVETIME']
        MJDREF = hdulist['EVENTS'].header['MJDREFI'] + hdulist['EVENTS'].header['MJDREFF']

        evtFileKeyWords = {'TELESCOPE': TELESCOPE, 'INSTRUME': INSTRUME, 'OBS_ID': OBS_ID, 'TSTART': TSTART,
                           'TSTOP': TSTOP, 'LIVETIME': LIVETIME, 'TIMESYS': TIMESYS,
                           'MJDREF': MJDREF, 'TIMEZERO': TIMEZERO, 'DATEOBS': DATEOBS}

        return evtFileKeyWords

    #################################################################
    def read_fpmsel(self):
        """
        Reads FPM_SEL extension from a NICER event file
        :return: FPMSEL_table_condensed - same as FPM_SEL but with total number of detectors per time stamp
        :rtype: pandas.DataFrame
        """
        evtFileKeyWords = self.readEF()
        TELESCOPE = evtFileKeyWords["TELESCOPE"]


        hdulist = fits.open(self.evtFile)

        # Reading the table
        FPMSEL_table = hdulist["FPM_SEL"].data
        TIME = FPMSEL_table['TIME']
        FPMsel = FPMSEL_table['FPM_SEL']
        FPMon = FPMSEL_table['FPM_ON']
        # Adding up all selected FPMs
        totfpmsel = np.zeros(np.size(TIME))
        for k in range(len(TIME)):
            totfpmsel[k] = (np.sum(FPMsel[k]))
        # Adding up all FPMs that were on
        totfpmon = np.zeros(np.size(TIME))
        for k in range(len(TIME)):
            totfpmon[k] = (np.sum(FPMon[k]))

        FPMSEL_table_condensed = pd.DataFrame(np.vstack((TIME, totfpmsel, totfpmon)).T, columns=['TIME', 'TOTFPMSEL', 'TOTFPMON'])

        return FPMSEL_table, FPMSEL_table_condensed

    #################################################################
    def readGTI(self):  # Reading fits event file GTI lists
        """
        Reads GTI table from an event
        :return: gtiList
        :rtype: numpy.ndarray
        """
        # Reading EF for some necessary keywords
        evtFileKeyWords = self.readEF()
        TELESCOPE = evtFileKeyWords["TELESCOPE"]

        hdulist = fits.open(self.evtFile)
        GTIdata = hdulist["GTI"].data
        ST_GTI = GTIdata.field("START")
        ET_GTI = GTIdata.field("STOP")
        gtiList = (np.vstack((ST_GTI, ET_GTI))).T

        return evtFileKeyWords, gtiList

    ################################################################################
    def filtenergy(self, eneLow: float, eneHigh: float):  # Filtering event file according to energy
        """
        Filters the event list according to energy (in keV)
        :param eneLow: low energy cutoff
        :type eneLow: float
        :param eneHigh: high energy cutoff
        :type eneHigh: float
        :return: dataTP_eneFlt, pandas dataframe of TIME and PI, filtered for energy
        :rtype: pandas.DataFrame
        """
        # Reading columns TIME and PI (pulse-invariant - proxy for photon energy) from binary table
        hdulist = fits.open(self.evtFile)
        tbdata = hdulist['EVENTS'].data

        piLow = eneLow / 0.01
        piHigh = eneHigh / 0.01
        dataTP = pd.DataFrame(np.vstack((tbdata.field('TIME'), tbdata.field('PI'))).T, columns=['TIME', 'PI'])
        dataTP_eneFlt = dataTP.loc[((dataTP['PI'] >= piLow) & (dataTP['PI'] <= piHigh))]

        return dataTP_eneFlt
