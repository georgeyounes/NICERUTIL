# nicerutil - Scripts to aid in [NICER](https://heasarc.gsfc.nasa.gov/docs/nicer/) data analysis

## Description

Several small modules and scripts that aid in the analysis of NICER data. Mainly, there is a 
background filtering routine that has been beneficial to my research, which flags what is likely 
background flares caused by high-energy particles. This component in NICER is known as PREL 
(precipitating electrons). For more information, I highly recommend reading 
[this page](https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/flares/).

Another script called **correctfpmsel** will filter the FPM_SEL table in an event file so that its 
time stamps (TIME column) fall within the GTIs. This is not always the case and, when some of the 
time stamps in FPM_SEL aren't covered by the orbit file of the observation, it causes the HEASoft 
barycentering routine ```BARYCORR``` to fail when attempting to barycenter the FPM_SEL TIME column. For 
this reason, ```BARYCORR```, in HEASoft versions 6.31 to 6.35, skips over the barycentering of the TIME column in 
the FPM_SEL table of an event file. Below, I will summarize how to regain the ```BARYCORR``` ability to 
apply barycenteric correction to the FPM_SEL table, and the use of the
script **correctfpmsel**. Note, that HEASoft 6.36 now corrects for this, and **correctfpmsel** can be safely skipped.  

There is a module that allows the user to download nicer data from AWS (get_nicer_aws.py). It 
can be accessed through the CLI **getniceraws**. The usual **getniceraws -h** will provide a list 
of arguments. Example usage:

```bash
getniceraws -sn "1E 1547.0-5408" -od ./1e1547/nicer/ -s 2022-12-29 -e 2023-07-25
```

where -sn (--srcname) and -od (--outdir) are required arguments, and represent the name of the source and the output
directory where the observation ids will be downloaded. The -sn should match a known source name; it will be parsed by
astropy SkyCoord for RA and DEC coordinates, which will be passed to Heasarc().query_region (from the package 
astroquery). Otherwise, the user could provide the RA and DEC directly through the argument -rd (--radec; 
e.g., -rd 85.046667 -69.331722) instead of -sn. The -od is where the observation ids will be downloaded. Note that the 
script does not create any directories by default (e.g. srcname), except the observation ids themselves. If in the -od 
an observation_id directory exists and is not empty, the script will skip the download of that observation. Finally, 
the -s (--start) and -e (--end) are optional arguments and allow the user to specify time range for the download, 
rather than the full nicer list of observations (default). Note that you will need to download the [AWS CLI 
tool](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html), which can be downloaded here. 



Just in case you are not aware of it, there is a much more extensive
nicer data analysis tools developed by other members of the NICER
team, primarly Paul Ray, which I highly recommend you to check
out. This can be accessed here [https://github.com/paulray/NICERsoft](https://github.com/paulray/NICERsoft)

## Acknowledgements

I am grateful to the discussions I have had with Craig Markwardt, Jeremy Hare, and the rest of the 
NICER team at Goddard Space Flight Center on NICER data, best analysis practices, and pitfalls.

## Installation

At the moment, nicerutil can be installed locally after cloning the directory or simply downloading 
and untarring it. Then from the NICERUTIL root directory:

```bash
  python -m pip install .
```

You may add the '-e' option to pip to install as an editable package.

The code has been tested on Python 3.9.16, and it requires matplotlib, pandas, scipy, astropy.

## Quick example usage

Upon installation, these command line scripts will be available to you: **createlightcurve**, 
**correctfpmsel**, and **flaghighenergyflares**. You can get their help messages and the list 
of required and optional arguments with the usual '-h' option. Here we shall cover the latter two. 

### Running correctfpmsel

**correctfpmsel** is really only useful if you wish to barycenter-correct the TIME column in the 
FPM_SEL table of a NICER event file, and only in the cases when the said TIME column has some time 
stamps outside the valid orbit file for the given observation, which causes ```BARYCORR``` to fail.
From [HEASoft version 6.31](https://heasarc.gsfc.nasa.gov/FTP/software/ftools/release/archive/Release_Notes_6.31)) 
to version 6.35,  ```BARYCORR``` skips over barycentering the TIME column of the FPM_SEL table. Hence, re-enable 
barycentering of the TIME column in the FPM_SEL table is required.

Starting in [HEASoft version 6.36](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/release_notes/RelNotes_6.36.0.html)
```BARYCORR``` again barycenters the FPM_SEL table, hence, you do not need to run **correctfpmsel** anymore. Running it
should not cause any harm, however.

#### Enabling FPM_SEL barycentering through BARYCORR

If you are running HEASoft version 6.31 to 6.35, and you wish to allow the heasoft tool ```BARYCORR``` 
to perform barycentering on the FPM_SEL table of a NICER event file, you can follow these steps: 
(1) in axBary.c, which can be found under heasoft-version/heagen/barycorr/, remove the block that states 
that "always ignore NICER FPM_SEL extension" (code on line 577), 
(2) navigate to heasoft-version/heagen/your-build/BUILD_DIR, and (3) run **./hmake** and then **./hmake install**.

In the above, heasoft-version directory refers to the version number, e.g., heasoft-6.33.1. Also, 
the your-build directory refers to your machine's architecture, e.g., aarch64-apple-darwin23.4.0. 
Your ```BARYCORR``` should now barycenter-correct the TIME column in the FPM_SEL table.

For the most part, ```BARYCORR``` will run without any errors. However, there are few instances when 
that is not the case, and ```BARYCORR``` exits with an error (something like TIME outside of orbit file).
The tool **correctfpmsel** will fix this issue by filtering out these time stamps.


To run **correctfpmsel**, all that you require is an event file

```bash
>> correctfpmsel ni6533062201_0mpu7_cl.evt
```

This will eliminate any 1-second TIME stamps in the FPM_SEL table of the event file that are not 
covered by the GTI table. You may now run ```BARYCORR``` on that event file, and it shall run through 
the barycentering of the TIME column in FPM_SEL without any errors.

### Running **flaghighenergyflares** 

In its simplest form, **flaghighenergyflares** can be run as follows: 

```bash
>> flaghighenergyflares ni7020500115_0mpu7_cl.evt ni7020500115.mkf
```

This script will produce a light curve in the 12-15 keV range (default parameters), 
where the NICER effective area is practically 0. The default time-bin used is 5 seconds. 
Any flare-like structure in the light curve is to be considered background. Under this 
assumption, the script will search for bins with a number of counts that is too large 
compared to the full light curve mean to be considered random Poisson fluctuation and 
flag it as a potential background flare. A xselect and NICERDAS-compatible GTI .txt and  
.fits file that correspond to the flagged events will be produced (outputFile"_gti.txt" and 
outputFile"_gti.fits" where outputFile is "flagged_flares" by default). A plot of the flagged 
events is also produced (outputFile".pdf"). A simple log file (nicerutil.log) is output with
essential information and any warning/error messages. **The user must initialize HEASoft 
before running this script for the GTI .fits file to be created!**

There are multiple flags that can help diagnose the extent of the flare contamination 
to your data. The -elb and -ehb control the energy range in keV within which to search 
for the flares. The -els and -ehs control the energy range in keV where source events 
supposedly lie (default is 1 and 5, respectively), this is used for plotting purposes 
only. The -tb is the time bin of the light curve in seconds (default is 5). The 
aggressiveness with which you wish to filter the background can be adjusted with the 
optional parameter --problim (or -pb=0.01 by default); the lower the number the less 
aggressive the cleaning is. The -of flag (--outputFile) defines the name of the output 
_gti.fits and _gti.txt GTI files, a .txt file that define the flagged bins, and a .pdf 
plot that shows a 4-row panels, the elb-ehb keV light curve, the els-ehs keV light curve, 
the FPM_OVERONLY_COUNT, and COR_SAX. The flare intervals are flagged in red in all four 
panels. 

The output of this particular example run, and the files required to run it,
can be found in the folder data. The output diagnosis plot

[flagged_flares.pdf](data%2Fflagged_flares.pdf)

shows the four aforementioned panels. The red bins are the ones flagged as flares and 
ultimately removed from the output GTI .fits file, **flagged_flares_gti.fits**. Note that 
the flare-like bins that are present in the source energy range light curve (second panel)
but absent from the out-of-range light curve (upper-panel) are intrinsic to the source; 
these are magnetar bursts (the source in question is 1E 1841-045 and the observation took 
place on August 21, at the height of its 2024 active period). 

It is worth noting that it really depends on your science whether, and how much, 
particle-background cleaning you want to apply. I personally work with sources on the 
faint end of things (~0.1 to few counts/s), and I care about observing the pulsed emission, 
which is sensitive to background light especially if the pulsed flux is low. Hence, I tend 
to be more conservative in my approach. Yet, for one of the sources with a count rate 
approaching the 60 counts/s, I tend to be more loose, since it is a relatively bright 
target. So, think of your science case and the source brightness when applying this 
filtering. If you decide to utilize this script, I encourage you to experiment with the 
input parameters to get a feel for the filtering done.

## Disclaimer

This code is distributed in the hope that it will be useful, but WITHOUT ANY EXPRESSED OR 
IMPLIED WARRANTY. Use this code solely if you understand what it is doing. Feel free to 
reach out with any questions.

nicerutil is unrelated to [NICERDAS](https://heasarc.gsfc.nasa.gov/docs/nicer/nicer_analysis.html), the 
official data analysis software package that is part of the HEASoft distribution. Nothing 
here should be construed as formally recommended by the NICER instrument team.

## License

[MIT](https://choosealicense.com/licenses/mit/)
