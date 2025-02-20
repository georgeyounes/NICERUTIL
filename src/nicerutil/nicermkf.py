"""
nicermkf.py is a module to perform simple operations on nicer MKF files
Likely to substantially change with time 
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from astroplan import moon_illumination

from nicerutil.nicerutil_logging import get_logger

import sys

sys.dont_write_bytecode = True

# Log config
############
logger = get_logger(__name__)


class MkfFileOps:
    """
        A class to operate on a nicer MKF table, as read from nicermkf.readmkf function

        Attributes
        ----------
        mkftable : pandas.DataFrame
            name of the MKF table

        Methods
        -------
        timefiltermkf(): Filters table according to a GTI
        trackingfiltermkf(): Filters table according to a standard tracking
        sunshinefiltermkf(): filters table according to orbit (night=0, day=1, both=2)
        sunanglefiltermkf(): filters table according to sun angle
        """

    def __init__(self, mkftable: str):
        """
        Operates on nicer mkf file
        :param mkftable: name of the mkf table as read from nicermkf.readmkf function
        :type mkftable: pandas.DataFrame
        """
        self.mkftable = mkftable

    #################################################################
    def timefiltermkf(self, gtilist):
        """
        Filters the mkf file according to a GTI table in nx2 numpy array format
        :param gtilist: GTI table
        :type gtilist: numpy.ndarray
        :return: timeflt_mkfData
        :rtype: pandas.DataFrame
        """
        timefiltered_mkf = pd.DataFrame()
        for kk, startstop in enumerate(gtilist):
            timefiltered_mkf_tmp = self.mkftable[
                (self.mkftable['tNICERmkf'] > startstop[0]) & (self.mkftable['tNICERmkf'] < startstop[1])]
            timefiltered_mkf = pd.concat([timefiltered_mkf, timefiltered_mkf_tmp], ignore_index=True)

        return timefiltered_mkf

    def trackingfiltermkf(self):
        """
        Filters the mkf file according to standard tracking filters
        From nimaketime:
        1- "inSAA == 0"
        2- (ATT_MODE==1 && ATT_SUBMODE_AZ==2 && ATT_SUBMODE_EL==2)"
        3- "ANG_DIST < DIST"      DIST = 0.015
        4- "st_valid = YES"
        5- "ELV > MINELV"         MINELV = 15 degrees
        6- "BR_EARTH > MIN_BR_EARTH"         MIN_BR_EARTH = 30 degrees
        :return: trackingfiltered_mkf
        :rtype: pandas.DataFrame
        """
        trackingfiltered_mkf = self.mkftable.loc[
            (self.mkftable['inSAA'] == 0) & (self.mkftable['starTrackerValid'] == 1)
            & (self.mkftable['ATT_MODE'] == 1) & (self.mkftable['ATT_SUBMODE_AZ'] == 2) & (
                    self.mkftable['ATT_SUBMODE_EL'] == 2)
            & (self.mkftable['ang_dist'] < 0.015) & (self.mkftable['elevation'] > 15) & (
                    self.mkftable['brightEarth'] > 20)
            ]

        return trackingfiltered_mkf

    def sunshinefiltermkf(self, sunshine=1):
        """
        Filters the mkf file according to day or night orbit
        :param sunshine: sunshine keyword
        :type sunshine: int
        :return: sunshinefiltered_mkf
        :rtype: pandas.DataFrame
        """
        valid_under = {1, 2}
        if sunshine not in valid_under:
            raise ValueError("sunshine: must be one of %r." % valid_under)

        sunshinefiltered_mkf = self.mkftable.loc[(self.mkftable['SUNSHINE'] == sunshine)]

        return sunshinefiltered_mkf

    def sunanglefiltermkf(self, sunang_ll=0, sunang_ul=180):
        """
        Filters the mkf file according to sun angle
        :param sunang_ll: sunangle lower-limit
        :type sunang_ll: float
        :param sunang_ul: sunangle upper-limit
        :type sunang_ul: float
        :return: sunanglefiltered_mkf
        :rtype: pandas.DataFrame
        """
        sunanglefiltered_mkf = self.mkftable.loc[((self.mkftable['SUN_ANGLE'] >= sunang_ll) &
                                                  (self.mkftable['SUN_ANGLE'] <= sunang_ul))]

        return sunanglefiltered_mkf

    def moonanglefiltermkf(self, moonang_ll=0, moonang_ul=180):
        """
        Filters the mkf file according to moon angle
        :param moonang_ll: Moon angle lower-limit
        :type moonang_ll: float
        :param moonang_ul: Moon angle upper-limit
        :type moonang_ul: float
        :return: moonanglefiltered_mkf
        :rtype: pandas.DataFrame
        """
        moonanglefiltered_mkf = self.mkftable.loc[((self.mkftable['MOON_ANGLE'] >= moonang_ll) &
                                                   (self.mkftable['MOON_ANGLE'] <= moonang_ul))]

        return moonanglefiltered_mkf

    def moonphasefiltermkf(self, moonphase_ll=0, moonphase_ul=100):
        """
        Filters the mkf file according to moon phase
        :param moonphase_ll: moon phase lower-limit
        :type moonphase_ll: float
        :param moonphase_ul: moon phase upper-limit
        :type moonphase_ul: float
        :return: moonphasefiltered_mkf
        :rtype: pandas.DataFrame
        """
        moonphasefiltered_mkf = self.mkftable.loc[((self.mkftable['MOONFRACTION'] >= moonphase_ll) &
                                                   (self.mkftable['MOONFRACTION'] <= moonphase_ul))]

        return moonphasefiltered_mkf

    def brightearthanglefiltermkf(self, brightearth_ll=0, brightearth_ul=180):
        """
        Filters the mkf file according to bright earth angle
        :param brightearth_ll: bright Earth lower-limit
        :type brightearth_ll: float
        :param brightearth_ul: bright Earth upper-limit
        :type brightearth_ul: float
        :return: sunanglefiltered_mkf
        :rtype: pandas.DataFrame
        """
        brightearthanglefiltered_mkf = self.mkftable.loc[((self.mkftable['brightEarth'] >= brightearth_ll) &
                                                          (self.mkftable['brightEarth'] <= brightearth_ul))]

        return brightearthanglefiltered_mkf

    def sunazfiltermkf(self, sunaz_ll=0, sunaz_ul=180):
        """
        Filters the mkf file according to sun azimuth angle
        :param sunaz_ll: sun azimuth (clocking) lower-limit
        :type sunaz_ll: float
        :param sunaz_ul: sun azimuth (clocking) upper-limit
        :type sunaz_ul: float
        :return: sunanglefiltered_mkf
        :rtype: pandas.DataFrame
        """
        sunazfiltered_mkf = self.mkftable.loc[((self.mkftable['AZ_SUN'] >= sunaz_ll) &
                                               (self.mkftable['AZ_SUN'] <= sunaz_ul))]

        return sunazfiltered_mkf

    def write_mkf_to_csv(self, output_file, saveindex=True):
        """
        Write mkf table to csv file output file .mkf.txt
        :param output_file: name of output CSV file
        :type output_file: str
        :param saveindex: whether to save index to output file
        :type saveindex: bool
        """
        self.mkftable.to_csv(output_file + '.txt', sep=' ', index=saveindex)

        return

    def merge_mkf_tables(self, mkftable_2):
        """
        Merge two mkf files
        :param mkftable_2: name of second mkf table
        :type mkftable_2: pandas.DataFrame
        """
        merged_mkfs = pd.concat([self.mkftable, mkftable_2], ignore_index=True)

        return merged_mkfs

    def addmoonfraction(self):
        """
        Add moon illuminated fraction to mkftable
        :return: moonadded_mkf
        :rtype: pandas.DataFrame
        """
        moonfraction = moonilluminatedfraction(self.mkftable['tNICERmkf_mjd'])
        mkftable_moonphase = self.mkftable.copy()
        mkftable_moonphase.loc[:, 'MOONFRACTION'] = moonfraction

        return mkftable_moonphase


def readmkffile(mkffile, under='MPU_UNDERONLY_COUNT'):

    valid_under = {'MPU_UNDERONLY_COUNT', 'MPU_UNDER_COUNT'}
    if under not in valid_under:
        raise ValueError("readmkffile: under input parameter must be one of %r." % valid_under)

    # open mkf file
    hdulist = fits.open(mkffile)

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
    KP_index = tbdata.field('KP')
    SUN_BETA = tbdata.field('BETA_ANGLE')

    if 'SUN_BODY_AZIMUTH' in hdulist['PREFILTER'].columns.names:
        # Check if sun clocking angle (SUN_BODY_AZIMUTH) in PREFILTER table as calculated with nicerl2
        AZ_SUN = tbdata.field('SUN_BODY_AZIMUTH')
    elif 'SUN_BODY_AZIMUTH' in hdulist['ORIG_PREFILTER'].columns.names:
        # Check if sun clocking angle (SUN_BODY_AZIMUTH) in ORIG_PREFILTER table as calculated with nicerl2
        AZ_SUN = hdulist['ORIG_PREFILTER'].data.field('SUN_BODY_AZIMUTH')
    elif 'AZ_SUN' in hdulist['PREFILTER'].columns.names:
        # Check if sun clocking angle (AZ_SUN) in PREFILTER table as calculated with nicerutil
        AZ_SUN = tbdata.field('AZ_SUN').T
    elif 'SUN_AZ' in hdulist['PREFILTER'].columns.names:
        # Check if sun clocking angle (SUN_AZ) in PREFILTER table as calculated with old version of nicerutil
        AZ_SUN = tbdata.field('SUN_AZ').T
    else:
        logger.warning('Cannot find Sun clocking angle (SUN_AZ) - if the purpose of your work is diagnostics, '
                       'nicerutil may break')
        AZ_SUN = np.array([None] * len(RA_sun))

    # Moon related
    MOON_ANGLE = tbdata.field('MOON_ANGLE')

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
    MPUoverCount = tbdata.field('MPU_OVERONLY_COUNT')

    if under == 'MPU_UNDERONLY_COUNT':
        MPUunderCount = tbdata.field('MPU_UNDERONLY_COUNT')
    else:
        MPUunderCount = tbdata.field('MPU_UNDER_COUNT')

    # ATTITUDE AND POINTING COLUMNS
    ATT_ANG_AZ = tbdata.field('ATT_ANG_AZ')
    ATT_ANG_EL = tbdata.field('ATT_ANG_EL')

    # This one does not have the per FPM data - a pain to implement in PANDAS
    mkfData_tmp = np.vstack(
        (tNICERmkf, tNICERmkf_mjd, RA_sun, DEC_sun, SUN_ANGLE, SUNSHINE, KP_index, SUN_BETA, AZ_SUN, MOON_ANGLE,
         RA_pointing, DEC_pointing, ROLL, ang_dist, elevation, brightEarth, starTrackerValid, ATT_MODE, ATT_SUBMODE_AZ,
         ATT_SUBMODE_EL, inSAA, corSax, ATT_ANG_AZ, ATT_ANG_EL))

    # Reading OVER_COUNT per FPM
    # Falttening out the per MPU table, which includes FPM information and include them as single column in the overall pandas data frame
    overCountALLFPMs = np.zeros((np.size(tNICERmkf)))

    for ii in range(7):  # MPUs
        for jj in range(8):  # FPM per MPU
            overCountALLFPMs = np.vstack((overCountALLFPMs, MPUoverCount[:, ii, jj]))

    overCountALLFPMs = overCountALLFPMs[1:]  # Removing the 0 preallocation - ugly for now
    overCountALLFPMs = np.where(overCountALLFPMs == 0, np.nan, overCountALLFPMs)  # Replace 0 with NAN

    # Reading UNDER_COUNT per FPM
    # Falttening out the per MPU table, which includes FPM information and include them as single column in the overall pandas data frame
    underCountALLFPMs = np.zeros((np.size(tNICERmkf)))

    for ii in range(7):  # MPUs
        for jj in range(8):  # FPM per MPU
            underCountALLFPMs = np.vstack((underCountALLFPMs, MPUunderCount[:, ii, jj]))

    underCountALLFPMs = underCountALLFPMs[1:]  # Removing the 0 preallocation - ugly for now
    underCountALLFPMs = np.where(underCountALLFPMs == 0, np.nan, underCountALLFPMs)

    # Merging all information so far into a single mkfData
    mkfData = np.vstack((mkfData_tmp, overCountALLFPMs, underCountALLFPMs)).T

    # converting the above to a dataframe
    mkftable_df = pd.DataFrame(mkfData,
                               columns=['tNICERmkf', 'tNICERmkf_mjd', 'RA_sun', 'DEC_sun', 'SUN_ANGLE', 'SUNSHINE',
                                        'KP_index', 'SUN_BETA', 'AZ_SUN', 'MOON_ANGLE', 'RA_pointing', 'DEC_pointing',
                                        'ROLL', 'ang_dist', 'elevation', 'brightEarth',
                                        'starTrackerValid', 'ATT_MODE', 'ATT_SUBMODE_AZ', 'ATT_SUBMODE_EL', 'inSAA',
                                        'corSax', 'ATT_ANG_AZ', 'ATT_ANG_EL', 'FPM_over00',
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
    return mkftable_df


def moonilluminatedfraction(eventtimes):
    """
    Calculate Moon illumination fraction
    :param eventtimes: TIME column from a mkf (or event) file
    :type eventtimes: numpy.ndarray
    :return: moonfraction
    :rtype: numpy.ndarray
    """
    eventtimes = Time(eventtimes, format='mjd')
    moonfraction = moon_illumination(eventtimes, ephemeris='jpl')

    return moonfraction
