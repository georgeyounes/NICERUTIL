"""
nicermkf.py is a module to perform simple operations on nicer MKF files
Likely to substantially change with time 
"""

import numpy as np
import pandas as pd

from astropy.io import fits

import sys
import argparse

from crimp.eventfile import EvtFileOps

sys.dont_write_bytecode = True


############################################################
# Class that performs simple operations on nicer MKF files #
############################################################

class MkfFileOps:
    """
        A class to operate on a nicer MKF files

        Attributes
        ----------
        mkffile : str
            name of the fits MKF file

        Methods
        -------
        readmkf(): reads essential column from PREFILTER table in MKF file
        filtertime(): Filters table according to a GTI nx2 array
        filtercolumn(): filters table according to a range in a given column
        """

    def __init__(self, mkffile: str):
        """
        Constructs the necessary attribute for the Phases object.

        :param mkffile: name of the fits nicer MKF file
        :type mkffile: str
        """
        self.mkffile = mkffile

    #################################################################
    def readmkf(self):
        """
        Reads essential column from PREFILTER table in nicer MKF file
        :return: mkfData_df - dataframe of read-in columns
        :rtype: pandas.DataFrame
        """
        hdulist = fits.open(self.mkffile)

        # Reading full PREFILTER table data from MKF file
        tbdata = hdulist['PREFILTER'].data

        # MKF time stamps
        tNICERmkf = tbdata.field('TIME')
        MJDREF = hdulist['PREFILTER'].header['MJDREFF'] + hdulist['PREFILTER'].header['MJDREFI']
        tNICERmkf_mjd = tNICERmkf / 86400 + MJDREF

        # Sun related
        RA_sun = tbdata.field('SUN_RA')
        DEC_sun = tbdata.field('SUN_DEC')
        SUN_ANGLE = tbdata.field('SUN_ANGLE')
        SUNSHINE = tbdata.field('SUNSHINE')
        AZ_SUN = tbdata.field('SUN_BODY_AZIMUTH')

        # Moon related
        moonAng = tbdata.field('MOON_ANGLE')

        # Pointing related
        RA_pointing = tbdata.field('RA')
        DEC_pointing = tbdata.field('DEC')
        ROLL = tbdata.field('ROLL')
        ang_dist = tbdata.field('ang_dist')
        elevation = tbdata.field('ELV')
        brightEarth = tbdata.field('BR_EARTH')
        starTrackerValid = tbdata.field('ST_VALID')

        # Tracking related
        ATT_MODE = tbdata.field('ATT_MODE')
        ATT_SUBMODE_AZ = tbdata.field('ATT_SUBMODE_AZ')
        ATT_SUBMODE_EL = tbdata.field('ATT_SUBMODE_EL')

        # SAA related
        inSAA = tbdata.field('SAA')

        # Instrument related
        corSax = tbdata.field('COR_SAX')
        TOToverCount = tbdata.field('TOT_OVER_COUNT')
        TOTunderCount = tbdata.field('TOT_UNDER_COUNT')
        MPUoverCount = tbdata.field('MPU_OVERONLY_COUNT')
        MPUunderCount = tbdata.field('MPU_UNDERONLY_COUNT')

        # ATTITUDE AND POINTING COLUMNS
        ATT_ANG_AZ = tbdata.field('ATT_ANG_AZ')
        ATT_ANG_EL = tbdata.field('ATT_ANG_EL')

        # This one does not have the per FPM data - a pain to implement in PANDAS
        mkfData_tmp = np.vstack(
            (tNICERmkf, tNICERmkf_mjd, RA_sun, DEC_sun, SUN_ANGLE, SUNSHINE, AZ_SUN, moonAng, RA_pointing,
             DEC_pointing, ROLL, ang_dist, elevation, brightEarth, starTrackerValid, ATT_MODE, ATT_SUBMODE_AZ,
             ATT_SUBMODE_EL, inSAA, corSax, ATT_ANG_AZ, ATT_ANG_EL, TOToverCount,
             TOTunderCount))

        # Reading OVER_COUNT per FPM
        # Falttening out the per MPU table, which includes FPM information and include them as single column in the overall pandas data frame
        overCountALLFPMs = np.zeros((np.size(tNICERmkf)))

        for ii in range(7):  # MPUs
            for jj in range(8):  # FPM per MPU
                overCountALLFPMs = np.vstack((overCountALLFPMs, MPUoverCount[:, ii, jj]))

        overCountALLFPMs = overCountALLFPMs[1:]  # Removing the 0 preallocation - ugly for now

        # Reading UNDER_COUNT per FPM
        # Falttening out the per MPU table, which includes FPM information and include them as single column in the overall pandas data frame
        underCountALLFPMs = np.zeros((np.size(tNICERmkf)))

        for ii in range(7):  # MPUs
            for jj in range(8):  # FPM per MPU
                underCountALLFPMs = np.vstack((underCountALLFPMs, MPUunderCount[:, ii, jj]))

        underCountALLFPMs = underCountALLFPMs[1:]  # Removing the 0 preallocation - ugly for now

        # Merging all information so far into a single mkfData
        mkfData = np.vstack((mkfData_tmp, overCountALLFPMs, underCountALLFPMs)).T

        # converting the above to a dataframe
        mkfData_df = pd.DataFrame(mkfData,
                                  columns=['tNICERmkf', 'tNICERmkf_mjd', 'RA_sun', 'DEC_sun', 'SUN_ANGLE', 'SUNSHINE',
                                           'AZ_SUN', 'moonAng', 'RA_pointing', 'DEC_pointing', 'ROLL', 'ang_dist',
                                           'elevation', 'brightEarth',
                                           'starTrackerValid', 'ATT_MODE', 'ATT_SUBMODE_AZ', 'ATT_SUBMODE_EL', 'inSAA',
                                           'corSax',
                                           'ATT_ANG_AZ', 'ATT_ANG_EL', 'TOToverCount', 'TOTunderCount', 'FPM_over00',
                                           'FPM_over01', 'FPM_over02', 'FPM_over03', 'FPM_over04', 'FPM_over05',
                                           'FPM_over06', 'FPM_over07', 'FPM_over10', 'FPM_over11', 'FPM_over12',
                                           'FPM_over13', 'FPM_over14', 'FPM_over15', 'FPM_over16', 'FPM_over17',
                                           'FPM_over20', 'FPM_over21', 'FPM_over22', 'FPM_over23', 'FPM_over24',
                                           'FPM_over25', 'FPM_over26', 'FPM_over27', 'FPM_over30', 'FPM_over31',
                                           'FPM_over32', 'FPM_over33', 'FPM_over34', 'FPM_over35', 'FPM_over36',
                                           'FPM_over37', 'FPM_over40', 'FPM_over41', 'FPM_over42', 'FPM_over43',
                                           'FPM_over44', 'FPM_over45', 'FPM_over46', 'FPM_over47', 'FPM_over50',
                                           'FPM_over51', 'FPM_over52', 'FPM_over53', 'FPM_over54', 'FPM_over55',
                                           'FPM_over56', 'FPM_over57', 'FPM_over60', 'FPM_over61', 'FPM_over62',
                                           'FPM_over63', 'FPM_over64', 'FPM_over65', 'FPM_over66', 'FPM_over67',
                                           'FPM_under00', 'FPM_under01', 'FPM_under02', 'FPM_under03', 'FPM_under04',
                                           'FPM_under05', 'FPM_under06', 'FPM_under07', 'FPM_under10', 'FPM_under11',
                                           'FPM_under12', 'FPM_under13', 'FPM_under14', 'FPM_under15', 'FPM_under16',
                                           'FPM_under17', 'FPM_under20', 'FPM_under21', 'FPM_under22', 'FPM_under23',
                                           'FPM_under24', 'FPM_under25', 'FPM_under26', 'FPM_under27', 'FPM_under30',
                                           'FPM_under31', 'FPM_under32', 'FPM_under33', 'FPM_under34', 'FPM_under35',
                                           'FPM_under36', 'FPM_under37', 'FPM_under40', 'FPM_under41', 'FPM_under42',
                                           'FPM_under43', 'FPM_under44', 'FPM_under45', 'FPM_under46', 'FPM_under47',
                                           'FPM_under50', 'FPM_under51', 'FPM_under52', 'FPM_under53', 'FPM_under54',
                                           'FPM_under55', 'FPM_under56', 'FPM_under57', 'FPM_under60', 'FPM_under61',
                                           'FPM_under62', 'FPM_under63', 'FPM_under64', 'FPM_under65', 'FPM_under66',
                                           'FPM_under67'])
        return mkfData_df

    #################################################################
    def timefiltermkf(self, gtilist):
        """
        Filters the mkf file according to a GTI table in nx2 numpy array format
        :param gtilist: GTI table
        :type gtilist: numpy.ndarray
        :return: timeflt_mkfData, pandas dataframe of TIME and PI, filtered for energy
        :rtype: pandas.DataFrame
        """
        fullmkf = self.readmkf()
        timefiltered_mkf = pd.DataFrame()
        for kk, startstop in enumerate(gtilist):
            timefiltered_mkf_tmp = fullmkf[
                (fullmkf['tNICERmkf'] > startstop[0]) & (fullmkf['tNICERmkf'] < startstop[1])]
            timefiltered_mkf = pd.concat([timefiltered_mkf, timefiltered_mkf_tmp], ignore_index=True)

        return timefiltered_mkf


if __name__ == '__main__':
    """
    Main function for nicermkf.py
    This runs the method timefiltermkf from class MkfFileOps
    """
    parser = argparse.ArgumentParser(description="For testing purposes only")
    parser.add_argument("mkfFile", help="A NICER MKF file", type=str)
    parser.add_argument("evtFile", help="A NICER event file", type=str)
    args = parser.parse_args()

    # GTI from event file
    _, GTI = EvtFileOps(args.evtFile).readGTI()
    # Perform the time filtering
    timefiltered_mkf_test = MkfFileOps(args.mkfFile).timefiltermkf(GTI)
    print(timefiltered_mkf_test)
