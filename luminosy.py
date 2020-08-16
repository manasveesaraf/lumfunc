"""Spreadsheet Column Printer

This script allows the user to print to the console all columns in the
spreadsheet. It is assumed that the first row of the spreadsheet is the
location of the columns.

This tool accepts comma separated value files (.csv) as well as excel
(.xls, .xlsx) files.

This script requires that `pandas` be installed within the Python
environment you are running this script in.

This file can also be imported as a module and contains the following
functions:

    * get_maggy - returns the column headers of the file
    * get_maggy_inv_var - the main function of the script
"""
# -----------------------
# Package Imports
# -----------------------
from scipy.optimize import curve_fit

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from kmeans_radec import KMeans, kmeans_sample
import kmeans_radec
import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)

# -----------------------
# Methods
# -----------------------
def get_maggy(app_mag_list: np.ndarray) -> np.ndarray:
    """
    Gets and prints the spreadsheet's header columns

    Parameters
    ----------
    app_mag_list : np.ndarray, list
        numpy array of all apparent magnitude values

    Returns
    -------
    np.ndarray, list
        numpy array of all correspoding maggy values

    """
    maggies_list = 10**(app_mag_list / (-2.5))
    return maggies_list

def get_maggy_inv_var(maggies_list: np.ndarray, app_mag_err_list: np.ndarray) -> np.ndarray:
    """
    Gets and prints the spreadsheet's header columns

    Parameters
    ----------
    maggies_list : np.ndarray, list
        numpy array of all maggy values
    app_mag_err_list : np.ndarray, list
        numpy array of all correspoding errors on apparent magnitude values

    Returns
    -------
    np.ndarray, list
        numpy array of all correspoding maggy inverse variance values 

    """
    inv_var_list = (0.4 * np.log(10) * maggies_list * app_mag_err_list)**(-2)
    return inv_var_list

def get_rest_mag(redshift_list: np.ndarray, app_mag_list: np.ndarray, maggy_ratio_list: np.ndarray) -> np.ndarray:
    """
    Gets and prints the spreadsheet's header columns

    Parameters
    ----------
    redshift_list : np.ndarray, list
        numpy array of all redshifts
    app_mag_list : np.ndarray, list
        numpy array of all corresponding apparent magnitudes
    maggy_ratio_list : np.ndarray, list
        numpy array of all corresponding maggy ratios

    Returns
    -------
    np.ndarray, list
        numpy array of all corresponding rest-frame magnitudes

    """
    # calculate luminosity distance
    lum_dist_list = cosmo.luminosity_distance(redshift_list).value
    print('\tLuminosity distance calculated.')

    # calculate abs mag
    abs_mag_list = app_mag_list - (5 * np.log10(lum_dist_list)) - 25
    print('\tAbsolute magnitude calculated.')

    # calculate K corrections
    Kcorr_list = -2.5 * np.log10(maggy_ratio_list)
    print('\tK-corrections calculated.')

    # calculate rest mag
    rest_mag_list = abs_mag_list - Kcorr_list
    print('\tRest-frame magnitude calculated.')

    return rest_mag_list

def get_volume(survey_area: float, redshift_list: np.ndarray) -> np.ndarray:
    """
    Gets and prints the spreadsheet's header columns

    Parameters
    ----------
    survey_area : float
        float value of survey area in sq. deg. 
    redshift_list : np.ndarray, list
        numpy array of all redshifts

    Returns
    -------
    np.ndarray, list
        numpy array of all corresponding comoving volumes

    """
    # calculate comoving distance
    com_dist_list = cosmo.comoving_distance(redshift_list).value
    print('\tComoving distance calculated.')

    # convert survey area to steradian
    survey_steradian = survey_area * ((np.pi / 180.)**2)
    print('\tSurvey area converted.')

    # calculate comoving volume
    vol_list = (com_dist_list**3) * (survey_steradian / 3)
    print('\tComoving volume calculated.')

    return vol_list

def get_binned_phi(rest_mag_list: np.ndarray, Vmax_list: np.ndarray, n_mag_bins: int) -> np.ndarray:
    """
    Gets and prints the spreadsheet's header columns

    Parameters
    ----------
    rest_mag_list : np.ndarray
        numpy array of all rest-frame magnitudes
    Vmax_list : np.ndarray
        numpy array of all corresponding maximum volumes
    n_mag_bins: int
        integer value of number of magnitude bins required

    Returns
    -------
    np.ndarray
        numpy array of mid-magnitude (i.e. x) value of each bin
    np.ndarray
        numpy array of magnitude-width/2 (i.e. x-error) value of each bin
    np.ndarray
        numpy array of phi (i.e. y) value of each bin (with h = 0.7)

    """
    # get bin_edges for diving the rest_mags in n_bins
    counts, bin_edges = np.histogram(rest_mag_list, bins=n_mag_bins)

    # sort rest_mag and Vmax lists per increasing mag
    sorted_index = np.argsort(rest_mag_list)
    sorted_Vmax_list = np.array(Vmax_list)[sorted_index]
    sorted_rest_mag_list = np.sort(rest_mag_list)

    # ceate empty lists for mid_M, phi and M_err
    mid_M_list = np.empty(n_mag_bins)
    M_err_list = np.empty(n_mag_bins)
    phi_list = np.empty(n_mag_bins)

    # loop over each bin
    for i in range(n_mag_bins):

        # find min and max M of bin
        max_M = bin_edges[i + 1]
        min_M = bin_edges[i]

        # add mid_M to list
        mid_M_list[i] = (min_M + max_M) / 2

        # add M_err to list
        M_err_list[i] = (abs(min_M) - abs(max_M)) / 2

        # find indicies upto the max_M
        up_lim_indices = np.where(sorted_rest_mag_list <= max_M)[0]

        # limit M and Vmax corresponding to max_M
        up_lim_rest_mag_list = sorted_rest_mag_list[up_lim_indices]
        up_lim_Vmax_list = sorted_Vmax_list[up_lim_indices]

        # find indicies from min_M to max_M value of bin
        if i != 0:
            lim_indices = np.where(up_lim_rest_mag_list > min_M)[0]
        else:
            lim_indices = np.where(up_lim_rest_mag_list >= min_M)[0]

        # limit Vmax corresponding from min_M to max_M
        Vmax_values = up_lim_Vmax_list[lim_indices]

        # calculate 1/Vmax
        phi_values = np.reciprocal(Vmax_values)
        # sum 1/Vmax all in this bin
        phi = sum(phi_values)

        # convert 1/Vmax to phi and add to list
        h = 0.7
        phi_list[i] = phi * ((h)**3) / M_err_list[i]

    return mid_M_list, M_err_list, phi_list

def get_patches(RA_list: np.ndarray,
                DEC_list: np.ndarray,
                n_patches: int,
                center_guesses: np.ndarray,
                survey='kids': str,
                numba_installed=True: bool,
                plot_savename='none': str) -> np.ndarray:
    """
    Gets and prints the spreadsheet's header columns

    Parameters
    ----------
    RA_list : np.ndarray
        numpy array of all RA values
    DEC_list : np.ndarray
        numpy array of all corresponding Dec values
    n_patches: int
        integer value of number of patches required
    center_guesses: np.ndarray
        (n_patches, 2) numpy array of patch center guesses [RA,Dec]
    survey: str
        string with survey name
    numba_installed: bool
        boolean - mark as False if numba is not installed
    plot_savename: str
        string with name and extension to save plot as

    Returns
    -------
    np.ndarray
        numpy array of patch assignment labels for each RA entry

    """
    # MAKE SURE ALL PATCHES ARE SITCHED ON SKY
    # works for most surveys - GAMA, KiDS - check rest
    if survey == 'kids':
        corrected_RA_list = np.where(RA_list > 320., RA_list - 360., RA_list)
    # use if a survey patch covers/connects over 320 degrees RA
    # and does not connect over 360 to 0 degree RA
    if survey != 'kids':
        corrected_RA_list = RA_list

    # STACK RA AND DEC AS X
    X = np.column_stack((corrected_RA_list, DEC_list))

    # FIND LABELS TO DIVIDE X INTO EQUAL n_patches
    if numba_installed:
        km = KMeans(center_guesses, method='fast')
    else:
        km = KMeans(center_guesses)
    labels = km.find_nearest(X)

    # VISUALISE ON PLOT
    if plot_savename != 'none':
        
        colors = cm.tab20(np.linspace(0, 1, n_patches))
        plt.figure(figsize=(10, 10))
        plt.suptitle("Galaxy Patches", fontsize=20)

        # get patch counts on histogram
        plt.subplot(211)
        plt.grid(True)
        N, b, p = plt.hist(labels, bins=n_patches)
        for n in range(n_patches):
            p[n].set_facecolor(colors[n])
        plt.xlabel("Label", fontsize=20)
        plt.ylabel("Count", fontsize=20)

        # get patches on sky
        plt.subplot(212)
        plt.grid(True)
        for n in range(n_patches):
            subset_indices = np.where(labels == n)
            plt.scatter(corrected_RA_list[subset_indices],
                        DEC_list[subset_indices],
                        color=colors[n],
                        s=1)
        # if 'gama' in datasetname:
            # plt.xlim(120, 240)
            # plt.ylim(-10, 10)
        # if 'kids' in datasetname:
            # plt.xlim(-50, 250)
            # plt.ylim(-40, 10)
        plt.xlabel("RA(J2000)/ deg", fontsize=20)
        plt.ylabel("Dec(J2000)/ deg", fontsize=20)
        
        plt.savefig(plot_savename, dpi=300)
        plt.show()

    return labels




