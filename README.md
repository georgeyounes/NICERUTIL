# nicerutil - Scripts to aid in [NICER](https://heasarc.gsfc.nasa.gov/docs/nicer/) data analysis

## Description

Several small modules and scripts that aid in the analysis of NICER data. Mainly, there is a 
background filtering routine that has been beneficial to my research, which flags what is likely 
background flares caused by high-energy particles. This component in NICER is known as PREL 
(precipitating electrons). For more information, I highly recommend reading 
[this page](https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/flares/).


Another script called **correctfpmsel** will filter the FPM_SEL table in an event file so that its 
time stamps match those of the EVENTS table. This is not always the case and has caused the HEASoft 
barycentering routine "BARYCORR" to fail. For this reason, "BARYCORR", in its current version, skips 
over the barycentering of the TIME column and keywords in the FPM_SEL table of an event file. Below, 
I will summarize how to regain the "BARYCORR" ability to apply barycenteric correction to the FPM_SEL 
table, and when to use the script **correctfpmsel**.

## Acknowledgements

I am grateful to the discussions I have had with Craig Markwardt, Jeremy Hare, Keith Gendreau, and 
the rest of the NICER team at Goddard Space Flight Center on NICER data, best analysis practices, and pitfalls.


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

To run **correctfpmsel**, all that you require is an event file

```bash
>> correctfpmsel ni6533062201_0mpu7_cl.evt
```

This will eliminate any 1-second TIME stamps in the FPM_SEL table of the event file that is not 
covered by the GTI table.

Running **flaghighenergyflares** in its simplest form is as follows: 

```bash
>> flaghighenergyflares ni6533062201_0mpu7_cl.evt
```

This script will produce a light curve in the 12-15 keV range, where the NICER effective area is 
practically 0. The default time-bin used is 5 seconds. Any flare-like structure in the light curve 
is to be considered background. The script will search for bins with a number of counts that 
is too large compared to the mean to be considered random Poisson fluctuation and flag it as a 
potential background flare. If desired the user could request an xselect, and NICERDAS-compatible 
GTI fits file that corresponds to the flagged events. A plot of the flagged events is also produced.

## Enabling FPM_SEL barycentering through BARYCORR

If you wish to allow the heasoft tool **barycorr** to perform barycentering on the FPM_SEL table 
of a NICER event file, you can follow these steps: (1) in axBary.c, which can be found under 
heasoft-version/heagen/barycorr/, remove the block that states that NICER FPM_SEL should stay 
untouched, (2) navigate to heasoft-version/heagen/your-build/BUILD_DIR, and (3) 
run **./hmake** and then **./hmake install**.

In the above, heasoft-version directory refers to the version number, e.g., heasoft-6.33.1. Also, 
the your-build directory refers to your machine's architecture, e.g., aarch64-apple-darwin23.4.0.

For the most part, the **barycorr** tool will run without any errors. However, there are few 
instances when the FPM_SEL has some time stamps that do not match the orbit file of the given 
observation, which causes barycorr to fail. The tool **correctfpmsel** will fix this issue by 
filtering out these time stamps.

## Disclaimer

This code is distributed in the hope that it will be useful, but WITHOUT ANY EXPRESSED OR IMPLIED 
WARRANTY. Use this code solely if you understand what it is doing. Feel free to reach out with any questions.

nicerutil is unrelated to [NICERDAS](https://heasarc.gsfc.nasa.gov/docs/nicer/nicer_analysis.html), the 
official data analysis software package that is part of the HEASoft distribution. Nothing here should 
be construed as formally recommended by the NICER instrument team.

## License

[MIT](https://choosealicense.com/licenses/mit/)