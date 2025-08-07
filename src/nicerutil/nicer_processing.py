import os
from nicerutil.get_nicer_aws import get_data
from convenience import * 


def get_nicer_data(datadir,nm,radius,ra,dec,start,end,make_toas,verbose=False):

    mkdir(os.path.join(datadir,'nicer'))
    oids = get_data('', radius, (ra, dec), start, end, os.path.join(datadir,'nicer'))

    for oid_dir in oids:
        oid = oid_dir.split('/')[-1]
        if not os.path.isfile(f"{oid_dir}/xti/event_cl/ni{oid}_0mpu7_cl.evt"):
            execute(f"nicerl2 {oid_dir}")
        subdir = f"{oid_dir}/xti/event_cl"
        execute(f"correctfpmsel {subdir}/ni{oid}_0mpu7_cl.evt")
        execute(f"flaghighenergyflares {subdir}/ni{oid}_0mpu7_cl.evt ni{oid}.mkf")
        if os.path.isfile(f"{oid_dir}/auxil/ni{oid}.orb"):
            orbfil = f"{oid_dir}/auxil/ni{oid}.orb"
        elif os.path.isfile(f"{oid_dir}/auxil/ni{oid}.orb.gz"):
            orbfil = f"{oid_dir}/auxil/ni{oid}.orb.gz"
        if os.path.isfile(f"{subdir}/ni{oid}_0mpu7_cl_barycorr.evt"):
            execute(f"rm {subdir}/ni{oid}_0mpu7_cl_barycorr.evt")
        execute(f"barycorr infile={subdir}/ni{oid}_0mpu7_cl.evt " + \
                f"outfile={subdir}/ni{oid}_0mpu7_cl_barycorr.evt " + \
                f"orbitfiles={orbfil} ra={ra} dec={dec} " + \
                f"barytime=no")
        execute(f"ls {subdir}/ni{oid}_0mpu7_cl_barycorr.evt | awk '" + \
                "{print $1" + '"[EVENTS]"}' + f"' >> evt_fils.txt")
    execute(f"ftmerge @evt_fils.txt merged_evtsonly.evt")
    execute(f"sed 's/EVENTS/GTI/g' evt_fils.txt > tmp; mv tmp nicer_evt_fils.txt")
    execute(f"ftmerge @nicer_evt_fils.txt merged_gtisonly.evt")
    execute(f"ftmerge 'merged_evtsonly.evt[GTI],merged_gtisonly.evt[GTI]' {nm}_nicer_merged.evt")

