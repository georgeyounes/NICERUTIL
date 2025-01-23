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


# Class that performs simple operations on nicer MKF files #


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
        fullmkf = self.mkftable
        timefiltered_mkf = pd.DataFrame()
        for kk, startstop in enumerate(gtilist):
            timefiltered_mkf_tmp = fullmkf[
                (fullmkf['tNICERmkf'] > startstop[0]) & (fullmkf['tNICERmkf'] < startstop[1])]
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
        fullmkf = self.mkftable
        trackingfiltered_mkf = fullmkf.loc[(fullmkf['inSAA'] == 0) & (fullmkf['starTrackerValid'] == 1)
                                           & (fullmkf['ATT_MODE'] == 1) & (fullmkf['ATT_SUBMODE_AZ'] == 2) & (
                                                   fullmkf['ATT_SUBMODE_EL'] == 2)
                                           & (fullmkf['ang_dist'] < 0.015) & (fullmkf['elevation'] > 15) & (
                                                   fullmkf['brightEarth'] > 20)
                                           ]

        return trackingfiltered_mkf

    def sunshinefiltermkf(self, sunshine=2):
        """
        Filters the mkf file according to day or night orbit
        :param sunshine: sunshine keyword
        :type sunshine: int
        :return: sunshinefiltered_mkf
        :rtype: pandas.DataFrame
        """
        fullmkf = self.mkftable
        if sunshine == 0:
            sunshinefiltered_mkf = fullmkf.loc[(fullmkf['SUNSHINE'] == sunshine)]
        elif sunshine == 1:
            sunshinefiltered_mkf = fullmkf.loc[(fullmkf['SUNSHINE'] == sunshine)]
        else:
            raise Exception("Orbit filtering: 0=night, 1=day, 2=everything")

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
        fullmkf = self.mkftable
        sunanglefiltered_mkf = fullmkf.loc[((fullmkf['SUN_ANGLE'] >= sunang_ll) &
                                            (fullmkf['SUN_ANGLE'] <= sunang_ul))]

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
        fullmkf = self.mkftable
        moonanglefiltered_mkf = fullmkf.loc[((fullmkf['MOON_ANGLE'] >= moonang_ll) &
                                            (fullmkf['MOON_ANGLE'] <= moonang_ul))]

        return moonanglefiltered_mkf

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
        fullmkf = self.mkftable
        brightearthanglefiltered_mkf = fullmkf.loc[((fullmkf['brightEarth'] >= brightearth_ll) &
                                                    (fullmkf['brightEarth'] <= brightearth_ul))]

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
        fullmkf = self.mkftable
        sunazfiltered_mkf = fullmkf.loc[((fullmkf['AZ_SUN'] >= sunaz_ll) &
                                         (fullmkf['AZ_SUN'] <= sunaz_ul))]

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


def readmkffile(mkffile):
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
        raise Exception('Sorry cannot find Sun clocking angle')

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
    TOToverCount = tbdata.field('TOT_OVER_COUNT')
    TOTunderCount = tbdata.field('TOT_UNDER_COUNT')
    MPUoverCount = tbdata.field('MPU_OVERONLY_COUNT')
    MPUunderCount = tbdata.field('MPU_UNDERONLY_COUNT')

    # ATTITUDE AND POINTING COLUMNS
    ATT_ANG_AZ = tbdata.field('ATT_ANG_AZ')
    ATT_ANG_EL = tbdata.field('ATT_ANG_EL')

    # This one does not have the per FPM data - a pain to implement in PANDAS
    mkfData_tmp = np.vstack(
        (tNICERmkf, tNICERmkf_mjd, RA_sun, DEC_sun, SUN_ANGLE, SUNSHINE, AZ_SUN, MOON_ANGLE, RA_pointing,
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
                                        'AZ_SUN', 'MOON_ANGLE', 'RA_pointing', 'DEC_pointing', 'ROLL', 'ang_dist',
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
    return mkftable_df


def define_nicerdetloc():
    """
    Define the NICER detectors in their geographical coordinates
    @return: nicerdetloc_geograph
    @rtype: numpy
    :return: nicerdetloc_geograph
    :rtype: numpy array
    """
    # Defining NICER detectors in geographical location
    # We will be plotting things per NICER detector
    nicerdetloc_geograph = (["06", "07", "16", "17", "27", "37", "47", "57",
                             "05", "15", "25", "26", "35", "36", "46", "56",
                             "04", "14", "24", "34", "44", "45", "54", "55",
                             "03", "13", "23", "33", "43", "53", "66", "67",
                             "02", "12", "22", "32", "42", "52", "64", "65",
                             "01", "11", "21", "31", "41", "51", "62", "63",
                             "00", "10", "20", "30", "40", "50", "60", "61"])

    return nicerdetloc_geograph


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
