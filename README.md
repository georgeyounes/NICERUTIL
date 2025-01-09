# nicerutil - Scripts to aid in [NICER](https://heasarc.gsfc.nasa.gov/docs/nicer/) data analysis

## Description

Several small modules and scripts that aid in the analysis of NICER data. Mainly, there is a 
background filtering routine that has been beneficial to my research, which flags what is likely 
background flares caused by high-energy particles. This component in NICER is known as PREL 
(precipitating electrons). For more information, I highly recommend reading 
[this page](https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/flares/).


Another script called **correctfpmsel** will filter the FPM_SEL table in an event file so that its time stamps (TIME column) fall within the GTIs. This is not always the case and, when some of the time stamps in FPM_SEL aren't covered by the orbit file of the observation, it causes the HEASoft barycentering routine "BARYCORR" to fail when attempting to barycenter the FPM_SEL TIME column. For this reason, "BARYCORR", in its current version, skips over the barycentering of the TIME column in the FPM_SEL table of an event file. Below, I will summarize how to regain the "BARYCORR" ability to apply barycenteric correction to the FPM_SEL table, and the use of the script **correctfpmsel**.

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

The code has been tested on Python 3.9.16, and it requires [CRIMP](https://github.com/georgeyounes/CRIMP), 
matplotlib, pandas, scipy, astropy.

## Quick example usage

Upon installation, these command line scripts will be available to you: **createlightcurve**, 
**correctfpmsel**, and **flaghighenergyflares**. You can get their help messages and the list 
of required and optional arguments with the usual '-h' option. Here we shall cover the latter two. 

### Running correctfpmsel

**correctfpmsel** is really only useful if you wish to barycenter-correct the TIME column in the FPM_SEL table of a NICER event file, and only in the cases when the said TIME column has some time stamps outside of the valid orbit file for the given observation, which causes **barycorr** to fail. Note that in its current form, **barycorr** skips over barycentering the TIME column of the FPM_SEL table (since [HEASoft version 6.31](https://heasarc.gsfc.nasa.gov/FTP/software/ftools/release/archive/Release_Notes_6.31)). Hence, we first need to re-enable barycentering of the TIME column in the FPM_SEL table.

#### Enabling FPM_SEL barycentering through BARYCORR

If you wish to allow the heasoft tool **barycorr** to perform barycentering on the FPM_SEL table of a NICER event file, you can follow these steps: (1) in axBary.c, which can be found under heasoft-version/heagen/barycorr/, remove the block that states that NICER FPM_SEL should stay untouched, (2) navigate to heasoft-version/heagen/your-build/BUILD_DIR, and (3) run **./hmake** and then **./hmake install**.

In the above, heasoft-version directory refers to the version number, e.g., heasoft-6.33.1. Also, the your-build directory refers to your machine's architecture, e.g., aarch64-apple-darwin23.4.0. Your **barycorr** should now barycenter-correct the TIME column in the FPM_SEL table.

For the most part, **barycorr** will run without any errors. However, there are few instances when that is not the case, and **barycorr** exits with an error (something like TIME outside of orbit file). The tool **correctfpmsel** will fix this issue by filtering out these time stamps.


To run **correctfpmsel**, all that you require is an event file

```bash
>> correctfpmsel ni6533062201_0mpu7_cl.evt
```

This will eliminate any 1-second TIME stamps in the FPM_SEL table of the event file that are not covered by the GTI table. You may now run **barycorr** on that event file, and it shall run through the barycentering of the the TIME column in FPM_SEL without any errors.

### Running **flaghighenergyflares** 

In its simplest form, **flaghighenergyflares** can be run as follows: 

```bash
>> flaghighenergyflares ni6533062201_0mpu7_cl.evt
```

This script will produce a light curve in the 12-15 keV range, where the NICER effective area is practically 0. The default time-bin used is 5 seconds. Any flare-like structure in the light curve is to be considered background. The script will search for bins with a number of counts that is too large compared to the mean to be considered random Poisson fluctuation and flag it as a potential background flare. If desired the user could request an xselect, and NICERDAS-compatible GTI fits file that corresponds to the flagged events. A plot of the flagged events is also produced.

As an example, I have added to the data folder an event file with such high-energy particle background (ni6020040101_0mpu7_cl.evt). Running 

```bash
>> flaghighenergyflares ni6020040101_0mpu7_cl.evt -elb 12 -ehb 15 -tb 5 -of ni6020040101_flares -cg ni6020040101_flares
```
The -elb and -ehb control the energy range in keV within which to search for the flares. The -els and -ehs control the energy range in keV where source events supposedly lies, this is used for plotting purposes only. The -tb is the time bin of the light curve. The aggressiveness with which you wish to filter the background can be adjusted with the optional parameter --problim (or -pb); the lower the number the less aggressive the cleaning is. The -cg argument (--creategti) will create a NICERDAS compatible GTI fits file that can be used with nicerl2 (note that HEASOFT should be initialized for this to properly run).

The output of this **flaghighenergyflares** run will produce three output files. ni6020040101_flares.pdf is a plot highlighting the eliminated flare intervals, ni6020040101_flares.fits is the corresponding GTI fits file that can be used with nicerl2 to filter out those flare intervals, and a log file nicerutil_logfile.log which holds some useful information. 

## Disclaimer

This code is distributed in the hope that it will be useful, but WITHOUT ANY EXPRESSED OR IMPLIED 
WARRANTY. Use this code solely if you understand what it is doing. Feel free to reach out with any questions.

nicerutil is unrelated to [NICERDAS](https://heasarc.gsfc.nasa.gov/docs/nicer/nicer_analysis.html), the 
official data analysis software package that is part of the HEASoft distribution. Nothing here should 
be construed as formally recommended by the NICER instrument team.

## License

[MIT](https://choosealicense.com/licenses/mit/)