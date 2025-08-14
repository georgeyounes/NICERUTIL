import os
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from astropy import units as u
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
from astroquery.heasarc import Heasarc
from subprocess import call
from nicerutil.convenience import * 
import shutil


def find_xmm_beams(position, start, end, heasarc=None, radius=0.2):
    if not heasarc:
        heasarc = Heasarc()
    query = heasarc.query_region(position, catalog="xmmssc", radius=radius * u.deg)

    start = Time(start,format='iso',scale='utc').mjd
    end = Time(end,format='iso',scale='utc').mjd

    date_cnd = np.array(query['time'] > start*u.day) & np.array(query['time'] < end*u.day)
    query = query[date_cnd]

    beams = np.unique([src[1:11] for src in query['srcid']]).astype(str)
    print(f"{len(beams)} XMM observations to process")
    return beams


def getXMMData(obsID,srcName,datadir,verbose=False):
    #############################################################
    ## This code downloads XMM data using an observation ID in the current path
    ## If the data already exists in the name "files.tar" the code skips the download
    ## It then untars it, and creates the proper setup for the observation
    ## by running cifbuild, odfingest, epchain
    ##
    ##
    ## Input:
    ## 1- exposureID :  10 digit XMM observation ID
    ## 2- srcName : Source name or field name
    ##
    ## output:
    ##
    ##
    ## Written by George Younes 2017 June 7
    ##
    ## Future work:
    ## 1- Needs more testing
    ##
    #############################################################

    """
    Changes made by Alex McEwen, 2024
    - added verbose flag
    - added calls to make_dir(obsID)
    - added obsID to 'files.tar' name
    - removed 'srcName' input, as it isn't used
    - removed emchain
    - added step to check for the existence of the final products of spectral step,
      and if they exist, move to the next step in the pipeline
    """

    srd = os.getcwd()
    cd(datadir,verbose=verbose,createdir=True)
    cd('xmm',verbose=verbose,createdir=True)
    cd(obsID,verbose=verbose,createdir=True)
    cd('odf',verbose=verbose,createdir=True)
    sasglob = glob("*.SAS")
    
    if len(sasglob) > 0:
        os.putenv("SAS_ODF", os.path.join(os.getcwd(),sasglob[0].split('/')[-1]))
        os.putenv("SAS_CCF", os.path.join(os.getcwd(),"ccf.cif"))
        return 0
        

    # If 'files_obsID.tar' does not exist, download and untar
    if not os.path.isfile(f'files_{obsID}.tar'):
        if verbose:
            print(' ------------------------------\n' + \
                 f' Downloading XMM obs ID {obsID}\n' + \
                  ' ------------------------------')
        
        path = f'"https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno={obsID}&level=ODF&instname=PN"'
        if verbose:
            command = f'curl -o files_{obsID}.tar {path}'
        else:
            command = f'curl -s -o files_{obsID}.tar {path}'    
        execute(command,verbose=verbose)


    if verbose:
        print(' -------------------------------\n' + \
             f' Untarring XMM obs ID {obsID}\n' + \
              ' -------------------------------')
    if verbose:
        command = f'tar -xvf files_{obsID}.tar >> {srcName}.log 2>&1 '
    else: 
        command = f'tar -xf files_{obsID}.tar  >> {srcName}.log 2>&1 '

    execute(command,verbose=verbose)
    if len(glob(f"????_{obsID}.TAR")) == 0:
        print(f"WARNING: XMM {obsID} download seems to have failed.")
        return 1

    # Removing unwanted files.tar
    command= f'rm -rf files_{obsID}.tar'
    execute(command,verbose=verbose)
    
    # Setting up SAS for observation
    if verbose:
        print(' ---------------------------------------------------------------------- \n' + \
              ' Running cifbuild, odfingest, epchain, and emchain. Creating setup file \n' + \
              ' ----------------------------------------------------------------------')

    TARSCIFILE = glob('*.TAR')
    if verbose:
        command = f'tar -xvf {TARSCIFILE[0]} >> {srcName}.log 2>&1'
    else:
        command = f'tar -xf {TARSCIFILE[0]} >> {srcName}.log 2>&1'
    execute(command,verbose=verbose)
    

    # Defining the ODF env variable
    os.putenv("SAS_ODF", os.getcwd())

    # Running cifbuild
    command = f' cifbuild  >> {srcName}.log 2>&1'
    execute(command,verbose=verbose)

    # Defining the cifbuild env variables
    os.putenv("SAS_CCF", os.path.join(os.getcwd(),"ccf.cif"))

    # Running odfingest
    command = f' odfingest  >> {srcName}.log 2>&1'
    execute(command,verbose=verbose)

    # Getting *.SAS filename
    odfSAS = glob("*.SAS")

    # Defining the SAS_ODF env variables
    os.putenv("SAS_ODF", os.path.join(os.getcwd(),odfSAS[0]))

    # Running epchain
    command = f' epchain  >> {srcName}.log 2>&1'
    execute(command,verbose=verbose)
    cd(srd,verbose=verbose)
    if verbose:
        print(' ------------------------------------ \n' + \
                ' End of download and setup XMM script \n' + \
                ' ------------------------------------ ')
    return 0

###############################################
## Function to correct PN flaring background ##
###############################################
def corrFlBackPN(evtFilePN,obsID,verbose=False):
    if verbose:
        print(' ----------------------------------------------------- \n' + \
                ' Correcting EPIC-PN event files for flaring background \n' + \
                ' ----------------------------------------------------- ')

    pnEvtFile = 'evtFilePN_'+ obsID + '_flBackCorr.fit'
    
    command = f" evselect  table={evtFilePN}:EVENTS withrateset=yes " + \
              f"rateset=fullFOVPNLC100s_{obsID}_above10keV.fit " + \
               "maketimecolumn=yes timecolumn=TIME timebinsize=100 makeratecolumn=yes withfilteredset=yes " + \
               "expression='(PATTERN == 0)&&#XMMEA_EP&&(FLAG == 0)&&(PI in [10000:12000])' " + \
               "filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes >> {srcName}.log 2>&1"
    execute(command,verbose=verbose)

    if verbose:
        print(' ----------------------- \n' + \
                ' Filtering PN event file \n' + \
                ' ----------------------- ')

    command = f" tabgtigen  table=fullFOVPNLC100s_{obsID}_above10keV.fit "+\
              f"expression='RATE<=0.4' gtiset=PNgti.fit >> {srcName}.log 2>&1"
    execute(command,verbose=verbose)

    
    command = f" evselect table={evtFilePN} withfilteredset=yes filteredset={pnEvtFile} destruct=Y " + \
               "keepfilteroutput=T expression='gti(PNgti.fit,TIME)' >> {srcName}.log 2>&1"
    execute(command,verbose=verbose)

    command = 'rm -f filtered.fits'
    execute(command,verbose=verbose)

    if verbose:
        print(' ----------------------------------------------------------------------------------- \n' + \
                ' Correcting EPIC-PN event files for flaring background done. Please check for errors \n' + \
                ' -----------------------------------------------------------------------------------')

    return pnEvtFile


###############################################
## Function to create images of the FoV ##
###############################################
def crtImPN(evtFilePN,obsID,elow,ehigh,srcName,verbose=False,tag=None):

    if verbose:
        print(' -------------------------------------------------------- \n' + \
                f' Creating EPIC-PN band images for event file {evtFilePN.split("/")[-1]}\n' + \
                 ' --------------------------------------------------------')

    print(os.getcwd())
    if tag is None:
        tag = f"_{elow:.1f}-{ehigh:.1f}.fits" 
    imFil = obsID + tag
    command = f" evselect  table={evtFilePN}:EVENTS withimageset=yes imageset={imFil} " + \
               "xcolumn=X ycolumn=Y imagebinning=imageSize ximagesize=600 yimagesize=600 " + \
               "withfilteredset=yes expression='(PATTERN <= 12)&&#XMMEA_EP&&(FLAG == 0)" + \
              f"&&(PI in [{elow}:{ehigh}])' filtertype=expression "+\
              f"keepfilteroutput=yes updateexposure=yes filterexposure=yes >> {srcName}.log 2>&1"
    execute(command,verbose=verbose)

    if verbose:
        print(' ------------------------------------------------ \n' + \
              ' Finished creating image. Please check for errors \n' + \
              ' ------------------------------------------------ ')

    command = 'rm -f filtered.fits'
    execute(command,verbose=verbose)

    return imFil

#######################################################
## Function to remove sources provided in input file ##
#######################################################

def extPtSrcsPN_fromfile(evtFilePN,imagePN,attFile,obsID,srcsToExtract,srcName,verbose=False):

    if verbose:
        print(' ------------------------ \n' + \
                ' Extracting input sources \n' + \
                ' ------------------------ ')

    outfile_bn = f"evtFilePN_{obsID}_allSrcsFlt"
    command = f" evselect  table={evtFilePN}:EVENTS withimageset=no withfilteredset=yes " + \
              f"expression='((PATTERN <= 4)&&#XMMEA_EP&&(FLAG == 0)&&(PI in [200:12000])"
    for src in srcsToExtract:

        ratmp  = src['RA']
        dectmp = src['DEC']
        exttmp = float([src['EP_EXTENT'] if src['EP_EXTENT']>0 else 20][0])*u.arcsec.to('deg')
        command += f"&&!((RA,DEC) in CIRCLE({ratmp},{dectmp},{exttmp}))"

    command += ")' filtertype=expression keepfilteroutput=yes updateexposure=yes " + \
              f"filterexposure=yes filteredset={outfile_bn}.fit >> {srcName}.log 2>&1"
    execute(command,verbose=verbose)

    return f"{outfile_bn}.fit"

########################################
## Function to create PN region files ##
########################################
def crtRegPN(image,pnEvtFile,attFile,obsID,srcName,ra,dec,extent,oid,cat,
             verbose=False):

    if verbose:
        print(' ------------------------------------- \n' + \
                ' Determining optimal extraction radius \n' + \
                ' ------------------------------------- ')

    srcRegFile, backRegFile = None, None
    band_suffix = image.split('/')[-1].split('_')[-1].split('.fits')[0]
    band = np.array(band_suffix.split('-')).astype(float)

    ################################################
    ################################################
    ################################################

    if verbose:
        print(' ----------------------------------------- \n' + \
                f' Creating radial profile for {srcName} \n' + \
                 ' ----------------------------------------- ')

    print('image: '+image,band_suffix)
    psfenergy = np.mean(band)

    othersources_cnd = np.array(cat['OBS_ID'] == obsID) & \
                       np.array(cat['IAUNAME'] != srcName.replace('_',' ').replace('p','+'))
    pnEvtFile_swisscheese = extPtSrcsPN_fromfile(pnEvtFile,
                                                 image,
                                                 attFile,
                                                 obsID,
                                                 cat['RA','DEC','EP_EXTENT'][othersources_cnd],
                                                 srcName,
                                                 verbose=verbose,
                                                 )
    swiss_image = crtImPN(pnEvtFile_swisscheese,
                          obsID,
                          band[0],
                          band[1],
                          srcName,
                          verbose=verbose,
                          tag=f'_{srcName}_swisscheese.im',
                          )
    swiss_image = mv(swiss_image,os.path.join(oid,srcName,swiss_image),verbose=verbose)
    command = f" eradial imageset={swiss_image} " + \
              f"srcexp='(RA,DEC) in circle({ra},{dec},{20*u.arcsec.to('deg')})'" + \
              f" psfenergy={psfenergy} centroid=yes >> {srcName}.log 2>&1"

    execute(command,verbose=verbose)

    # save profile fit, read in values
    prof_fit = os.path.join(oid,srcName,f"radProf_{srcName}_{obsID}_{band_suffix}.fit")
    command = f"mv radprof.ds {prof_fit}"
    execute(command,verbose=verbose)

    with fits.open(prof_fit) as hdulist:
        data = hdulist[1].data


    plt.xscale('log')
    plt.yscale('log')
    plt.grid()

    plt.ylabel('Radial Profile [cts/arcsec^2]')
    plt.xlabel('Inner edge of Radial bin [arcsec]')
    plt.title(prof_fit.split('/')[-1])

    outfile = os.path.join(oid,srcName,f"radProf_{srcName}_{obsID}_{band_suffix}.png")
    x = data['RAD_LO']
    y = data['RPROF']
    yerr = data['RPROF_ERR']
    cnd = np.array(x>0)&np.array(x<150)&np.array(y>0)
    x = x[cnd]
    y = y[cnd]
    yerr = yerr[cnd]

    try:
        xfine = np.linspace(x.min(),x.max(),1000)
    except:
        raise ValueError(f"Problem with radial profile; check {prof_fit}")

    fit_params = create_params(amp   = {'value':2,   'vary':True},
                               alpha = {'value':8,   'vary':True},
                               r0    = {'value':100, 'vary':True},
                               dc    = {'value':2,   'vary':True}
                              )

    try:
        result = minimize(king_fit,
                          fit_params,
                          args=(x,y),
                          method='nelder',
                         )
    except:
        raise ValueError(f"Problem with radial profile fit; check {prof_fit}")


    amp,alpha,r0,dc = [result.params[key].value for key in result.params.keys()]
    k = king_dc(xfine,amp,alpha,r0,dc)
    norm_k = (k-dc)/np.trapz((k-dc),xfine)

    csum = np.cumsum(norm_k)*np.diff(xfine)[0]
    try:
        r_opt = xfine[csum > 0.99][0]
    except:
        raise ValueError(f"Problem with radial profile fit/normalization; check {prof_fit}")

    plt.plot(xfine,k,ls='--',color='black',lw=1,label='King+BG Fit')
    plt.axvline(r_opt,ls='--',color='lime',label=r"R$_{99}$"+f" = {r_opt:2.3f} arcsec")

    plt.errorbar(x,y,yerr=yerr,marker='x',lw=0,elinewidth=1,markersize=2,label='Profile')
    plt.axhline(dc,ls='--',color='red',label=f"{dc:2.2e}"+r' ct/as$^2$')


    plt.legend(loc='lower left')
    plt.xlim([1,200])
    plt.savefig(outfile)
    plt.clf()


    ################################################
    ################################################
    ################################################

    if verbose:
        print(' -------------------------------\n' + \
                ' Creating source DS9 regionfile \n' + \
                ' -------------------------------')

    srcRegFile = f"src_{srcName}_{obsID}_{band_suffix}_pn.reg"

    circPos = f"circle({ra},{dec},{r_opt*u.arcsec.to('deg')})"
    f = open(srcRegFile,'w+')
    f.write('# Region file format: DS9 version 4.1\n' + \
            'global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" ' + \
            'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n' + \
            'physical\n' + \
           f'j2000; {circPos} # text=' + "{" + srcName.split('_')[-1] + '}\n')
    f.close()

    all_src_reg = os.path.join(oid,'all_srcs.reg')
    if os.path.exists(all_src_reg):
        f = open(all_src_reg,'a')
    else:
        f = open(all_src_reg,'w')
        f.write('# Region file format: DS9 version 4.1\n' + \
                'global color=green dashlist=8 3 width=2 font="helvetica 3 normal roman" ' + \
                'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n' + \
                'physical\n')
    f.write(f"j2000; {circPos} # text=" + "{" + srcName.split('_')[-1] + "}\n")
    f.close()

    ################################################
    ################################################
    ################################################

    if verbose:
        print(' -----------------------------------\n' + \
                ' Creating background DS9 regionfile \n' + \
                ' -----------------------------------')

    backRegFile = f"back_{srcName}_{obsID}_{band_suffix}_pn.reg"

    f = open(backRegFile,'w+')
    f.write('# Region file format: DS9 version 4.1\n' + \
            'global color=red dashlist=8 3 width=2 font="helvetica 10 normal roman" ' + \
            'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n' + \
            'physical\n' + \
           f'j2000; annulus({ra},{dec},0.04861,0.06944)\n')
    f.close()



    ################################################
    ################################################
    ################################################
    if verbose:
        print(' --------------------------------- \n' + \
                ' Running crtRegFiles.py -- Success \n' + \
                ' --------------------------------- ')

    return srcRegFile, backRegFile, pnEvtFile_swisscheese


#######################################################################
## Function to create spectral files for a source and its background ##
#######################################################################
def crtSpecPN(evtFilePN,
              evtFilePN_allSrcFlt,
              obsID,
              srcName,
              srcRegFile,
              backRegFile,
              prefix,
              verbose=False,
              ):

    if verbose:
        print(' ------------------------------------------------------------- \n' + \
               f' Extracting PN spectral files for source {srcName}\n' + \
                ' -------------------------------------------------------------')

    
    ############################
    # Reading source regionfiles
    ############################
    data_file = open(srcRegFile,'r')
    circPos = data_file.readlines()[-1].rstrip('\n').split(';')[-1].split('#')[0]
    data_file.close()
            
    ################################    
    # Reading background regionfiles
    ################################
    data_file = open(backRegFile,'r')
    annPos = data_file.readlines()[-1].rstrip('\n').split(';')[-1]
    data_file.close()
            
    ########################################
    # Spectral name definition for later use
    ########################################
    srcPhaFile    = f"{prefix}{srcName}_{obsID}_pnSrc.pha"
    backPhaFile   = f"{prefix}{srcName}_{obsID}_pnBack.pha"
    respFILE      = f"{prefix}{srcName}_{obsID}_pnSrc.rmf"
    arfFile       = f"{prefix}{srcName}_{obsID}_pnSrc.arf"
    grpFile5      = f"{prefix}{srcName}_{obsID}_pnSrc_grp5.pha"
    
    ############################
    # Extracting source spectrum
    ############################
    command = f" evselect  table={evtFilePN}:EVENTS expression='(PATTERN <= 4)&&#XMMEA_EP&&(FLAG == 0)" + \
              f"&&(PI in [200:10000])&&((RA,DEC) IN {circPos})' filtertype=expression keepfilteroutput=yes " + \
               "updateexposure=yes filterexposure=yes withfilteredset=yes withspectrumset=yes " + \
              f"spectrumset={srcPhaFile} spectralbinsize=5 withspecranges=yes specchannelmin=0 " + \
               "specchannelmax=20479 energycolumn=PI"
    execute(command,verbose=verbose)

    
    ################################
    # Extracting background spectrum
    ################################
    command = f" evselect  table={evtFilePN_allSrcFlt}:EVENTS expression='(PATTERN <= 4)&&#XMMEA_EP&&(FLAG == 0)" + \
              f"&&(PI in [200:10000])&&((RA,DEC) IN {annPos})' filtertype=expression keepfilteroutput=yes " + \
               "updateexposure=yes filterexposure=yes withfilteredset=yes withspectrumset=yes " + \
              f"spectrumset={backPhaFile} spectralbinsize=5 withspecranges=yes specchannelmin=0 " + \
               "specchannelmax=20479 energycolumn=PI"
    execute(command,verbose=verbose)


    #####################################################################################
    # Scaling region sizes, extracting response and ancillary files, and grouping spectra
    #####################################################################################
    command = f" backscale  spectrumset={srcPhaFile} badpixlocation={evtFilePN}"
    execute(command,verbose=verbose)
    command = f" backscale  spectrumset={backPhaFile} badpixlocation={evtFilePN_allSrcFlt}"
    execute(command,verbose=verbose)

    # Extracting response file
    command = f" rmfgen  spectrumset={srcPhaFile} rmfset={respFILE}"
    execute(command,verbose=verbose)

    # Extracting ancillary file
    command = f" arfgen  arfset={arfFile} spectrumset={srcPhaFile} withrmfset=yes " + \
              f"rmfset={respFILE} badpixlocation={evtFilePN}"
    execute(command,verbose=verbose)

    # Grouping spectra
    command = f" specgroup  spectrumset={srcPhaFile} addfilenames=yes backgndset={backPhaFile} " + \
              f"rmfset={respFILE} arfset={arfFile} mincounts=5 groupedset={grpFile5}"
    execute(command,verbose=verbose)

    if verbose:
        print(' --------------------------------------------------- \n' + \
             ' Extracting PN spectra done. Please check for errors \n' + \
             ' --------------------------------------------------- ')

    command = 'rm -f filtered.fits'
    execute(command,verbose=verbose)


########################################################################
## Function to create PN timing files for a source and its background ##
########################################################################
def crtTimPN(evtFilePN,
             evtFilePN_allSrcFlt,
             obsID,
             srcName,
             srcRegFile,
             backRegFile,
             ra,
             dec,
             band,
             oid,
             verbose=False,
             overwrite=False,
             ):

    if verbose:
        print(' ----------------------------------------------------------- \n' + \
             f' Extracting PN Timing files for source {srcName} \n' + \
              ' ----------------------------------------------------------- ')


    ############################
    # Reading source regionfiles
    ############################
    data_file = open(srcRegFile,'r')
    circPos = data_file.readlines()[-1].rstrip('\n').split(';')[-1].split('#')[0]
    data_file.close()

    ################################
    # Reading background regionfiles
    ################################
    data_file = open(backRegFile,'r')
    annPos = data_file.readlines()[-1].rstrip('\n').split(';')[-1]
    data_file.close()

    ###############################################################################################
    # Looping through bands to search and creating source/background eventfiles, then barycentering
    ###############################################################################################
    filescreated = False
    band_str_lb = [str(band[0]).replace('.','p') if str(band[0]).split('.')[-1] != '0' else int(band[0])][0]
    band_str_ub = [str(band[1]).replace('.','p') if str(band[1]).split('.')[-1] != '0' else int(band[1])][0]
    srcTimFile  =  os.path.join(oid,
                                srcName,
                                f"{srcName}_{obsID}_pn_tim{band_str_lb}-{band_str_ub}keV.fit")
    backTimFile =  os.path.join(oid,
                                srcName,
                                f"{srcName}_{obsID}_pnBack_tim{band_str_lb}-{band_str_ub}keV.fit")
    if not (os.path.isfile(srcTimFile) and os.path.isfile(backTimFile)) or overwrite:
        command = f" evselect  table={evtFilePN}:EVENTS expression='(PATTERN <= 12)&&#XMMEA_EP&&(FLAG == 0)" + \
                  f"&&(PI in [{int(1000*band[0])}:{int(1000*band[1])}])&&((RA,DEC) IN {circPos})' " + \
                   "filtertype=expression keepfilteroutput=yes updateexposure=yes filterexposure=yes " + \
                  f"withfilteredset=yes filteredset={srcTimFile}"
        execute(command,verbose=verbose)

        command = f" evselect table={evtFilePN_allSrcFlt}:EVENTS expression='(PATTERN <= 12)" + \
                  f"&&#XMMEA_EP&&(FLAG == 0)&&(PI in [{int(1000*band[0])}:{int(1000*band[1])}])" + \
                  f"&&((RA,DEC) IN {annPos})' filtertype=expression keepfilteroutput=yes updateexposure=yes "+\
                  f"filterexposure=yes withfilteredset=yes filteredset={backTimFile}"
        execute(command,verbose=verbose)
        command = f' barycen table={srcTimFile}:EVENTS withsrccoordinates=yes srcra={ra} ' + \
                  f'srcdec={dec} ephemeris=DE405'
        execute(command,verbose=verbose)

        filescreated = True
    else:
        r_str = f'Timing files already exist in {os.path.join(oid,srcName,"timingFiles")}' + \
                 '; moving to next source'

    if verbose:
        print(' -------------------------------------------------------- \n' + \
                ' Extracting PN Timing files done. Please check for errors \n' + \
                ' -------------------------------------------------------- ')

    command = 'rm -f filtered.fits'
    execute(command,verbose=verbose)
    return srcTimFile

def get_xmm_data(srcname, position, datadir, start, 
                 end, radius, elow=0.1, ehigh=10, extent=10, verbose=True):

    cat = Table.read("4XMM_DR13cat_v1.0_full.fits")
    heasarc = Heasarc()
    beams = find_xmm_beams(position,start,end,heasarc,radius)
    for beam in beams:
        retcode = getXMMData(beam,srcname,datadir,verbose=verbose)
        if retcode == 1:
            print(f"WARNING: {beam} has failed")
            continue
        odfdir = os.path.join(datadir,'xmm',beam,'odf')
        oiddir = os.path.join(datadir,'xmm',beam)
        pnEvtFile = os.path.join(odfdir,f"evtFilePN_{beam}_flBackCorr.fit")
        if not os.path.exists(pnEvtFile):

            try:
                pnEvtFileRaw = glob(odfdir+'/P*PN*EV*')[0]
            except:
                if verbose:
                    print(f"ERROR: Raw event file is missing - download may"+\
                          f" have failed. Skipping {beam}.")   
                continue
            if verbose:
                print(f"making event file {pnEvtFile}")
            pnEvtFile = corrFlBackPN(pnEvtFileRaw,beam,verbose=verbose)
            attFile = odfdir + '/P' + beam + 'OBX000ATTTSR0000.FIT'
            pnSrcReg  = os.path.join(oiddir,f"src_{srcname}_{beam}_pn.reg") 
            pnBackReg = os.path.join(oiddir,f"back_{srcname}_{beam}_pn.reg")
            img = crtImPN(pnEvtFile,beam,elow,ehigh,srcname,verbose=verbose)
            img = mv(img,oiddir,verbose=verbose)
            pnSrcReg, pnBackReg, pnSwiss = crtRegPN(img, pnEvtFile, attFile,
                                                    beam, srcname,
                                                    position.ra.value,
                                                    position.dec.value,
                                                    extent*u.arcsecond.to('deg'),
                                                    oiddir,
                                                    cat,
                                                    verbose=verbose)
            crtTimPN(pnEvtFile, pnSwiss, beam, srcname, pnSrcReg, 
                      pnBackReg, position.ra.value, position.dec.value,
                      [elow,ehigh], oiddir, verbose=verbose)
            mv(pnEvtFile,"results/{srcname}/")



