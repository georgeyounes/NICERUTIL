import argparse
import numpy as np
import os
from astropy.coordinates import SkyCoord
from glob import glob
from nicerutil.get_nicer_aws import get_data
from subprocess import call
from xmm_processing import get_xmm_data
from nicer_processing import get_nicer_data
from convenience import *
 


def process_magnetar(model, srcname, datadir, 
                     radius, start, end, numiter,
                     numcomp, do_nicer,
                     do_xmm, elow=0.1,
                     ehigh=10., extent=10):

    with open(model,'r') as tmp:
        lns = tmp.readlines()
    for ln in lns:
        if ln.split()[0] == "RAJ":
            ra = ln.split()[1].rstrip('\n')
        elif ln.split()[0] == "DECJ":
            dec = ln.split()[1].rstrip('\n')   
        elif ln.split()[0] in ["PSR", "PSRJ"]:
            nm = ln.split()[1].rstrip('\n') 
    if srcname is not None:
        nm = srcname
    else:
        srcname = nm
    pos = SkyCoord(ra,dec,unit=('hourangle','deg'))
    ra, dec = pos.ra.value, pos.dec.value

    if datadir is None:
        datadir = srcname

    build_directories(datadir)
    headdir = os.getcwd()
    if do_nicer:
        get_nicer_data(datadir,nm,radius,ra,dec,start,end,verbose=True)
    if do_xmm:
        get_xmm_data(nm, pos, os.path.join(headdir,datadir), start, end, 
                     radius, elow=elow, ehigh=ehigh, extent=extent,
                     verbose=True)
    cd(headdir)
    cleanup(nm)
                


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download NICER data for a provided source " + \
                                                 "and conduct standard TOA generation and timing")
    parser.add_argument("-m", type=str, required=True,
                        help="timing model (.par) file")
    parser.add_argument("-n", type=str, default=None,
                        help="name of source (default read from par file)")
    parser.add_argument("-o", type=str, default=None,
                        help="output directory to store data (default [src name])")
    parser.add_argument("-r", type=float, default=0.05,
                        help="radius of search for observations [deg], default 0.05")
    parser.add_argument("-s", type=str, default="2017-01-01",
                        help="date of earliest observation to include (YYYY-MM-DD)")
    parser.add_argument("-e", type=str, default="2032-12-31",
                        help="date of latest observation to include (YYYY-MM-DD)")
    parser.add_argument("-i", type=int, default=3,
                        help="number of iterations for fitting template")
    parser.add_argument("-c", type=int, default=3,
                        help="number of components in profile template")
    parser.add_argument("-el", type=float, default=0.1,
                        help="low energy limit [keV]")
    parser.add_argument("-eh", type=float, default=10.,
                        help="high energy limit [keV]")
    parser.add_argument("-x", type=bool, default=False,
                        help="Download and process XMM data",
                        action=argparse.BooleanOptionalAction)
    parser.add_argument("-N", type=bool, default=False,
                        help="Download and process NICER data",
                        action=argparse.BooleanOptionalAction)
    args = parser.parse_args()


    process_magnetar(args.m, args.n, args.o, 
                     args.r, args.s, args.e, 
                     args.i, args.c, 
                     args.N, args.x)


