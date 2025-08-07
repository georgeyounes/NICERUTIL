from glob import glob
from subprocess import call
import os
import shutil


def build_directories(datadir,srcname):
    for fld in [datadir,'results']:
        mkdir(fld)
        mkdir(f"{fld}/{srcname}")

def mkdir(newdir,verbose=False):
    if not os.path.exists(newdir):
        os.mkdir(newdir)
        if verbose:
            print(f"{newdir} created")
    else:
        if verbose:
            print(f"{newdir} alredy exists")

def cd(newdir, createdir=False, verbose=False):
    if not os.path.isdir(newdir) and createdir:
        if verbose:
            print(f"Creating {newdir}")
        mkdir(newdir,verbose)
    
    os.chdir(newdir)
    if verbose: 
        print(f"moved to {os.getcwd()}")

def execute(cmd,verbose=False):
    if verbose:
        print(cmd)
    call(cmd,shell=True)

def mv(fil,newloc,verbose=False):
    if os.path.exists(os.path.join(newloc,fil)):
        os.remove(os.path.join(newloc,fil))
    shutil.move(fil,os.path.join(newloc,fil))
    if verbose:
        print(f"moved {fil} to {newloc}")
    return newloc+'/'+fil

def cleanup(sourcename):
    mvstr = ""
    for glb in ["flagged_flares.txt", "nicerutil.log", f"{sourcename}*", 
                "ToA*", "merged*", "evt_fils.txt","*fit"]:
        if 'par' not in glb and len(glob(glb))>0:
            mvstr += f" {glb}"
    if mvstr != "":
        execute(f"mv {mvstr} results/{sourcename}")

