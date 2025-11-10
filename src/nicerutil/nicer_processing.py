import sys

from nicerutil.get_nicer_aws import get_data
from nicerutil.convenience import *
import yaml
import argparse

from nicerutil.correctfpmsel import correctfpmsel
from nicerutil.correct_gtis import correct_gtis
from nicerutil.flagbackgrflares import flagbackgrflares
from pathlib import Path
from nicerutil.nicerutil_logging import get_logger
from nicerutil.eventfile import EvtFileOps
from nicerutil.nicermkf import pick_mkf

sys.dont_write_bytecode = True

# Log config
############
logger = get_logger(__name__)

# Setting up environment variables ti query heasoft from python
os.environ["HEADASNOQUERY"] = ""
os.environ["HEADASPROMPT"] = "/dev/null"


def run_nicerl2(obsID, indir='.', evtsuffix='NONE', tasks='all',
                threshfilter='night', gtifiles='NONE', clobber='no'):
    """
    Run a simple version of nicerl2
    """
    command_nicerl2 = (f"nicerl2 indir={indir}/{obsID} evtsuffix={evtsuffix} tasks={tasks} "
                       f"threshfilter={threshfilter} gtifiles={gtifiles} clobber={clobber}")
    execute(command_nicerl2)

    if evtsuffix == 'NONE':
        level2_evtfile = f"{indir}/{obsID}/xti/event_cl/ni{obsID}_0mpu7_cl.evt"
    else:
        level2_evtfile = f"{indir}/{obsID}/xti/event_cl/ni{obsID}_0mpu7_cl_{evtsuffix}.evt"

    return level2_evtfile


def run_barycentering(eventfile, outfile, orbitfile, ra, dec, clobber):
    """
    Run barycentering correction
    """
    command_barycorr = (f"barycorr infile={eventfile} outfile={outfile} ra={ra} dec={dec} "
                        f"orbitfiles={orbitfile} clobber={clobber} refframe=ICRS ephem=JPLEPH.430 chatter=5")
    execute(command_barycorr)

    return outfile


def run_ftmerge(event_list, outputfile_merged):
    """
    ftmerge: merge a list of eventfile as it appears in a text file
    """
    execute(f"ftmerge @{event_list} {outputfile_merged}")

    return outputfile_merged


def load_yaml_source_param(config_path):
    """
    Load a YAML source parameter file for processing as a Python dictionary.
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def process_obsid(obsID, indir, radec, evtsuffix='', tasks='all', threshfilter='night',
                  gtifiles='NONE', clobber='no', eneLow_back=12, eneHigh_back=15,
                  timebin=5, lcthresh=0.1, probLim=0.01, eneLow_src=1,
                  eneHigh_src=5, outputFile="flagged_flares"):
    if not folder_exists(obsID):
        raise Exception(f"Folder {obsID} does not exist in current working directory.")

    # 1- run nicerl2
    ################
    os.system('nigeodown')
    run_nicerl2(obsID, indir=indir, evtsuffix=evtsuffix, tasks=tasks,
                threshfilter=threshfilter, gtifiles=gtifiles, clobber=clobber)

    cd(f"{obsID}/xti/event_cl/")
    if (evtsuffix.lower() == 'none') or (evtsuffix == ''):
        level2_evtfile = f"ni{obsID}_0mpu7_cl.evt"
    else:
        level2_evtfile = f"ni{obsID}_0mpu7_cl_{evtsuffix}.evt"

    # 2- run correct_gtis
    #####################
    correct_gtis(level2_evtfile)

    # 3- run correctfpmsel
    ######################
    correctfpmsel(level2_evtfile)

    # 4- run flagbackgrflares, and then nicerl2 again if necessary
    ##############################################################
    auxil_folder = Path("../../auxil/")  # finding the mkffile
    mkffile = str(pick_mkf(auxil_folder, obsID))

    _, distinct_flares = flagbackgrflares(level2_evtfile, mkffile, eneLow_back=eneLow_back,
                                          eneHigh_back=eneHigh_back, timebin=timebin, lcthresh=lcthresh,
                                          probLim=probLim, eneLow_src=eneLow_src, eneHigh_src=eneHigh_src,
                                          outputFile=outputFile)

    if (distinct_flares is not None) and (gtifiles.lower() != 'none'):
        flagged_gtis = f"{outputFile}_gti.fits"  # naming convention of flare GTI
        # Note that STDGTI is name of HDU table when gti is built with heasoft tools
        execute(f"ftmgtime {flagged_gtis},{gtifiles} flaregti_and_usergti.fits AND clobber=YES")
        # These are not necessary with ftmgtime
        #execute(f"ftmergesort 'flaregti_and_usergti.fits[STDGTI]' flaregti_and_usergti_sorted.fits START")
        #execute(f"mv -f flaregti_and_usergti_sorted.fits flaregti_and_usergti.fits")
        # rerun nicerl2
        cd("../../../")  # You must be in the base directory where obsID resides to run nicerl2
        run_nicerl2(obsID, indir='.', evtsuffix=evtsuffix, tasks='SCREEN', threshfilter=threshfilter,
                    gtifiles=f"{obsID}/xti/event_cl/flaregti_and_usergti.fits", clobber='YES')
        cd(f"{obsID}/xti/event_cl/")

    elif (distinct_flares is not None) and (gtifiles.lower() == 'none'):
        flagged_gtis = f"{outputFile}_gti.fits"  # naming convention of flare GTI
        cd("../../../")  # You must be in the base directory where obsID resides to run nicerl2
        run_nicerl2(obsID, indir='.', evtsuffix=evtsuffix, tasks='SCREEN', threshfilter=threshfilter,
                    gtifiles=f"{obsID}/xti/event_cl/{flagged_gtis}", clobber='YES')
        cd(f"{obsID}/xti/event_cl/")

    else:
        logger.info('\n No flares found.')

    # 5- barycenter correct data - only if data is available
    ########################################################
    exposure_afterprocessing = EvtFileOps(level2_evtfile).readEF()['ONTIME']
    if exposure_afterprocessing == 0:
        logger.info('\n Exposure total 0 after all filtering - skipping barycentering.')
    else:
        ra, dec = radec
        file = next(auxil_folder.glob(f"ni{obsID}.orb*"), None)
        orbfile = str(file)
        bary_outfile = f"{level2_evtfile.split('.')[0]}_bc.evt"
        run_barycentering(level2_evtfile, bary_outfile, orbfile, ra, dec, clobber='YES')

        logger.info('\n Barycentering correction complete.')

    # 6- returning to base directory
    ################################
    cd("../../../")

    return None


def process_source(nicer_process_config_file="./nicer_process_config.yaml"):
    """
    Full processing of nicer data for a given source for observation ids in a given time interval
    """
    # Get yaml parameters as dictionary
    nicer_cf = load_yaml_source_param(nicer_process_config_file)

    # Create output directory if it does not exist
    outdir_download = os.path.join(nicer_cf['outdir'], 'nicer')
    mkdir(outdir_download)

    # Downloading all requested obsids
    # oids: observation id full path
    oids = get_data(nicer_cf['radius'], nicer_cf['start'], nicer_cf['end'],
                    radec=nicer_cf['radec'], outdir=outdir_download)

    # let's first cd to the directory where all observations (oids) reside
    cd(outdir_download)

    for oid_dir in oids:
        obsID = oid_dir.split('/')[-1]

        # Process each obs ID according to threshfilter
        if nicer_cf['threshfilter'] == 'night':
            process_obsid(obsID, './', nicer_cf['radec'], evtsuffix=nicer_cf['evtsuffix_night'],
                          tasks=nicer_cf['tasks'], threshfilter=nicer_cf['threshfilter'],
                          gtifiles=nicer_cf['gtifiles'], clobber=nicer_cf['clobber'],
                          eneLow_back=nicer_cf['eneLow_back'], eneHigh_back=nicer_cf['eneHigh_back'],
                          timebin=nicer_cf['timebin'], lcthresh=nicer_cf['lcthresh'], probLim=nicer_cf['probLim'],
                          eneLow_src=nicer_cf['eneLow_src'], eneHigh_src=nicer_cf['eneHigh_src'],
                          outputFile=f"ni{obsID}_ff_night")

        elif nicer_cf['threshfilter'] == 'day':
            process_obsid(obsID, './', nicer_cf['radec'], evtsuffix=nicer_cf['evtsuffix_day'],
                          tasks=nicer_cf['tasks'], threshfilter=nicer_cf['threshfilter'],
                          gtifiles=nicer_cf['gtifiles'], clobber=nicer_cf['clobber'],
                          eneLow_back=nicer_cf['eneLow_back'], eneHigh_back=nicer_cf['eneHigh_back'],
                          timebin=nicer_cf['timebin'], lcthresh=nicer_cf['lcthresh'], probLim=nicer_cf['probLim'],
                          eneLow_src=nicer_cf['eneLow_src'], eneHigh_src=nicer_cf['eneHigh_src'],
                          outputFile=f"ni{obsID}_ff_day")

        elif nicer_cf['threshfilter'] == 'all':
            # First we extract night data
            process_obsid(obsID, './', nicer_cf['radec'], evtsuffix=nicer_cf['evtsuffix_night'],
                          tasks=nicer_cf['tasks'], threshfilter='night',
                          gtifiles=nicer_cf['gtifiles'], clobber=nicer_cf['clobber'],
                          eneLow_back=nicer_cf['eneLow_back'], eneHigh_back=nicer_cf['eneHigh_back'],
                          timebin=nicer_cf['timebin'], lcthresh=nicer_cf['lcthresh'], probLim=nicer_cf['probLim'],
                          eneLow_src=nicer_cf['eneLow_src'], eneHigh_src=nicer_cf['eneHigh_src'],
                          outputFile=f"ni{obsID}_ff_night")

            # Second we extract day data - tasks = SCREEN no need to run everything from scratch
            process_obsid(obsID, './', nicer_cf['radec'], evtsuffix=nicer_cf['evtsuffix_day'],
                          tasks='SCREEN', threshfilter='day',
                          gtifiles=nicer_cf['gtifiles'], clobber=nicer_cf['clobber'],
                          eneLow_back=nicer_cf['eneLow_back'], eneHigh_back=nicer_cf['eneHigh_back'],
                          timebin=nicer_cf['timebin'], lcthresh=nicer_cf['lcthresh'], probLim=nicer_cf['probLim'],
                          eneLow_src=nicer_cf['eneLow_src'], eneHigh_src=nicer_cf['eneHigh_src'],
                          outputFile=f"ni{obsID}_ff_day")

        else:
            raise ValueError("threshfilter must be one of either 'day', 'night', or 'all'.")

    return oids


def main():
    parser = argparse.ArgumentParser(description="Processing nicer observations of a given target")
    parser.add_argument("yamlfile", help="NICER processing details written in yaml format", type=str)
    args = parser.parse_args()

    process_source(nicer_process_config_file=args.yamlfile)


if __name__ == '__main__':
    main()
