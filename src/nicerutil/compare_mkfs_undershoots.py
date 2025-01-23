"""
Compare the undershoots in N sets of MKF files
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
# import glob

from nicerutil.nicermkf import MkfFileOps, readmkffile, define_nicerdetloc, define_nicerdetloc
from nicerutil.mkfdiagnostics import mkf_diagnostics, createalldiagnosticsplots

import sys
import argparse
import os

sys.dont_write_bytecode = True


def comparemkfundershoots(listofmkffiles, sunshine=2, sunAngLR=45, sunAngUR=180, moonAngLR=0, moonAngUR=180,
                          brearthLR=0, brearthUR=180, sunAzLR=-180, sunAzUR=180, timepostleak=True,
                          number_largestfpms_to_flag=12, diagnosticsplots=False, outputdir='compare_mkfs'):
    # Create folder in working directory to drop plots and files in
    if not os.path.exists(outputdir):
        # if the parent directory is not present, create it
        os.makedirs("outputdir")

    # Reading mkfs
    for kk, mkffile in enumerate(listofmkffiles):
        filteredmkf, average_undershoot_perFPM, _ = mkf_diagnostics(mkffile, sunshine=sunshine, sunAngLR=sunAngLR,
                                                                    sunAngUR=sunAngUR, moonAngLR=moonAngLR,
                                                                    moonAngUR=moonAngUR, brearthLR=brearthLR,
                                                                    brearthUR=brearthUR, sunAzLR=sunAzLR,
                                                                    sunAzUR=sunAzUR, timepostleak=timepostleak)

        # Move on if mkf is empty after filtering
        if filteredmkf.empty:
            print('{} file: DataFrame  after all filtering is empty - moving on'.format(mkffile))
            return

        # Suffix for several file names
        suffix = mkffile.split(".")[0]

        # Create working directory for specific cuts
        try:
            directory_specific_cuts = ('ss' + sunshine +'_sl' + sunAngLR + '_su' + sunAngUR + '_ml' + moonAngLR + '_mu'
                                       + moonAngUR + '_bel' + brearthLR + '_beu' + brearthUR + '_al' + sunAzLR + '_au'
                                       + sunAzUR)
            os.mkdir(directory_specific_cuts)
            command = 'mv ' + directory_specific_cuts + ' ' + outputdir
            os.system(command)
        except FileExistsError:
            print(f"Directory '{directory_specific_cuts}' already exists.")

        # Create all diagnostics plots and move them to their own directory 'outputdir'
        if diagnosticsplots:
            try:
                os.mkdir('diagnostics_plots')
                command = 'mv diagnostics_plots ' + outputdir + '/' + directory_specific_cuts
                os.system(command)
            except FileExistsError:
                print(f"Directory diagnostics_plots already exists.")
            nicDET_geograph = define_nicerdetloc()
            createalldiagnosticsplots(filteredmkf, average_undershoot_perFPM, nicDET_geograph, suffix)
            command = 'mv *.png ./' + outputdir + '/' + directory_specific_cuts + '/diagnostics_plots/'
            os.system(command)

        # Merging the mkf files
        if kk == 0:
            merged_average_undershoot_perFPM = average_undershoot_perFPM
            merged_average_undershoot_perFPM = merged_average_undershoot_perFPM.rename(
                columns={'average': 'average_' + str(kk), 'stdv': 'stdv_' + str(kk), 'median': 'median_' + str(kk)})
        else:
            merged_average_undershoot_perFPM = pd.merge(merged_average_undershoot_perFPM, average_undershoot_perFPM,
                                                        how='outer', left_index=True, right_index=True)
            merged_average_undershoot_perFPM = merged_average_undershoot_perFPM.rename(
                columns={'average': 'average_' + str(kk), 'stdv': 'stdv_' + str(kk), 'median': 'median_' + str(kk)})

    # Drop detector 63 which is consistently larger - should be same column for all
    merged_average_undershoot_perFPM = merged_average_undershoot_perFPM.drop(pd.Index(['63']))

    # remove rows of the 9 detectors with largest undershoot in each column
    # Count those detectors and we hope they are close to 9 rather than close to 63
    max_indices_all = merged_average_undershoot_perFPM.apply(lambda col: col.nlargest(
        number_largestfpms_to_flag).index.tolist())
    max_indices_all_flat_unique = np.unique(max_indices_all.to_numpy().flatten())
    max_indices = pd.Index(max_indices_all_flat_unique)
    merged_average_undershoot_perFPM_maxremoved = merged_average_undershoot_perFPM.drop(max_indices)

    # Remove all nans (many will be removed in the steps above)
    nan_indices = merged_average_undershoot_perFPM_maxremoved.loc[
        (merged_average_undershoot_perFPM_maxremoved['median'].isna())].index
    merged_average_undershoot_perFPM_maxnanremoved = merged_average_undershoot_perFPM_maxremoved.drop(
        nan_indices)

    # Separate merged dataframe into average and median (and stdv)
    merged_average_undershoot_perFPM_maxnanremoved_dropstdv = merged_average_undershoot_perFPM_maxnanremoved[
        merged_average_undershoot_perFPM_maxnanremoved.columns.drop(
            list(merged_average_undershoot_perFPM_maxnanremoved.filter(regex='stdv')))]
    # Get median out
    merged_median_undershoot_perFPM_allclean = merged_average_undershoot_perFPM_maxnanremoved_dropstdv[
        merged_average_undershoot_perFPM_maxnanremoved_dropstdv.columns.drop(
            list(merged_average_undershoot_perFPM_maxnanremoved_dropstdv.filter(regex='average')))]
    # Get average out
    merged_average_undershoot_perFPM_allclean = merged_average_undershoot_perFPM_maxnanremoved_dropstdv[
        merged_average_undershoot_perFPM_maxnanremoved_dropstdv.columns.drop(
            list(merged_average_undershoot_perFPM_maxnanremoved_dropstdv.filter(regex='median')))]

    # Create correlation matrix and pair plot for each of the above
    corr_matrix_median = plotcorrandpair(merged_median_undershoot_perFPM_allclean, outputfile=suffix)
    corr_matrix_average = plotcorrandpair(merged_average_undershoot_perFPM_allclean, outputfile=suffix)

    # Writing final pandas dataframe to .csv file
    MkfFileOps(merged_average_undershoot_perFPM_maxnanremoved).write_mkf_to_csv(suffix+'_clean_averageinfo')

    # Get indices (detectors) and save them to csv file
    indices_good_det = merged_average_undershoot_perFPM_allclean.DataFrame(
        merged_average_undershoot_perFPM_allclean.index, columns=['Index'])
    # Save the index DataFrame to a CSV file
    MkfFileOps(indices_good_det).write_mkf_to_csv(suffix+'_accepteddetectors', saveindex=False)

    # Moving files to proper directory
    command = 'mv *.png *.txt ./' + outputdir + '/' + directory_specific_cuts
    os.system(command)

    # Visualize and create a list of on/off detectors
    # to-do

    return (merged_median_undershoot_perFPM_allclean, merged_average_undershoot_perFPM_allclean,
            corr_matrix_median, corr_matrix_average)


def readmkfsfromtextfile(mkfsintextfile):
    with open(mkfsintextfile) as file:
        listofmkffiles = [line.rstrip() for line in file]
    return listofmkffiles


def plotcorrandpair(mkf_averages, outputfile):
    corr_matrix = mkf_averages.corr()

    # Plot the correlation matrix using seaborn
    plt.figure(figsize=(18, 16))
    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(corr_matrix, mask=mask, cmap='coolwarm', vmin=0, vmax=1, annot=True, square=True,
                linewidths=.5, cbar_kws={"shrink": .5})
    plt.title('Correlation Matrix')
    # Save the plot
    plt.tight_layout()
    plotname = outputfile + '_corrmatrix.png'
    plt.savefig(plotname, format='png', dpi=200)
    plt.close()

    # Plotting the pair plots
    sns.pairplot(mkf_averages, corner=True, diag_kind='kde')
    # Save the plot
    plt.tight_layout()
    plotname = outputfile + '_pairplot.png'
    plt.savefig(plotname, format='png', dpi=200)
    plt.close()

    return corr_matrix


def main():
    parser = argparse.ArgumentParser(description="Compare mkf files to a baseline mkf file")
    parser.add_argument("mkffiles", help=".txt list of NICER MKF files", type=str)
    parser.add_argument("-ss", "--sunshine", help="Filtering for sunshine, 0 for night, 1 for day, "
                                                  "and 2 for no filtering (default=1)", type=int, default=2)
    parser.add_argument("-sl", "--sunAngLR", help="Filtering for sun angle, lower-range", type=float,
                        default=45)
    parser.add_argument("-su", "--sunAngUR", help="Filtering for sun angle, upper-range", type=float,
                        default=180)

    parser.add_argument("-ml", "--moonAngLR", help="Filtering for moon angle, lower-range", type=float,
                        default=0)
    parser.add_argument("-mu", "--moonAngUR", help="Filtering for moon angle, upper-range", type=float,
                        default=180)
    parser.add_argument("-bel", "--brearthLR", help="Filtering for bright Earth angle, lower-range",
                        type=float, default=0)
    parser.add_argument("-beu", "--brearthUR", help="Filtering for bright Earth angle, upper-range",
                        type=float, default=180)
    parser.add_argument("-al", "--sunAzLR", help="Filtering for sun azimuth (clocking) angle, "
                                                 "lower-range", type=float,
                        default=-180)
    parser.add_argument("-au", "--sunAzUR", help="Filtering for sun azimuth (clocking) angle, "
                                                 "upper-range", type=float,
                        default=180)
    parser.add_argument("-tf", "--timepostleak", help="Accept only time post repair (default=True)",
                        type=bool, default=True, action=argparse.BooleanOptionalAction)
    parser.add_argument("-of", "--outputfile", help="name of output plotfiles", type=str,
                        default='mkf_diagnostics')
    args = parser.parse_args()

    # Reading mkf .txt file
    listofmkffiles = readmkfsfromtextfile(args.mkffiles)
    # Running primary function
    comparemkfundershoots(listofmkffiles, args.sunshine, args.sunAngLR, args.sunAngUR, args.moonAngLR, args.moonAngUR,
                          args.brearthLR, args.brearthUR, args.sunAzLR, args.sunAzUR, args.timepostleak,
                          args.outputfile)


if __name__ == '__main__':
    main()
