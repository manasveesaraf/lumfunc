"""Luminosity Function Constructor and Modeller

This script allows the user to construct and model Galaxian Luminosity Functions using the 1/Vmax estimator and Schechter function. 

Rest-frame magnitudes and spatial variance on the counts can be obtained. 
Plotting function for easy visualisation are included.

This file can also be imported as a module and contains the following
functions:

    * get_maggy - converts magnitudes into maggies
    * get_maggy_inv_var - returns inverse variances on maggies
    * get_obs_maggies_file - saves file of calculated maggies and inverse variances
    * get_rec_maggies_files - saves file of reconstructed maggies at input redshift
    * get_rest_maggy_ratio_file - saves file of calculated rest-frame maggy ratios
    * get_rest_mag - converts apparent magnitudes into rest-frame magnitudes
    * get_maggy_ratio_file - saves file of calculated reconstructed maggy ratios
    * get_all_maggy_ratios_file - consolidates files of calculated maggy ratios
    * get_volume - returns comoving volume of input survey area and redshift
    * get_binned_phi - bins and weights galaxy counts per magnitude by 1/Vmax
    * get_patch_centers - saves file of centers of equal patches over survey area
    * get_patch_labels - divides survey into equal patches
    * get_binned_phi_error - returns spatial variance of the luminosity function 
    * get_plot - plots magnitude-binned and 1/Vmax weighted luminosity function
    * filter_plot_by_colour - plots luminosity functions filtered by galaxy colour
    * SchechterMagModel - single Schechter function in terms of magnitude
    * DoubleSchechterMagModel - double Schechter function in terms of magnitude 
    * get_gof - returns reduced chi squared estimate of goodness of fit
    * get_schechter_phi - best fits single Schechter function on data 
    * get_double_schechter_phi - best fits double Schechter function on data 

"""
# -----------------------
# Package Imports
# -----------------------
# import kcorrect

import numpy as np

from typing import Tuple

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

import kmeans_radec
from kmeans_radec import KMeans, kmeans_sample

from astropy.io import ascii

from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)

# -----------------------
# Methods
# -----------------------
def get_maggy(app_mag_list: np.ndarray) -> np.ndarray:
    """
    Converts magnitudes into maggies.
    
    Parameters
    ----------
    app_mag_list : np.ndarray
        apparent magnitude of each data point (galaxy)

    Returns
    -------
    np.ndarray
        all corresponding maggy values

    """
    maggies_list = 10**(app_mag_list / (-2.5))
    return maggies_list

def get_maggy_inv_var(maggies_list: np.ndarray, 
                      app_mag_err_list: np.ndarray) -> np.ndarray:
    """
    Returns inverse variances on maggies using maggies and magnitude errors.
    
    Parameters
    ----------
    maggies_list : np.ndarray
        maggy value of each data point (galaxy)
    app_mag_err_list : np.ndarray
        all corresponding errors on apparent magnitude values

    Returns
    -------
    np.ndarray
        all correspoding maggy inverse variance values 

    """
    
    inv_var_list = (0.4 * np.log(10) * maggies_list * app_mag_err_list)**(-2)
    return inv_var_list

def get_obs_maggies_file(obs_maggies_outfile_name: str,
                         bands: str,
                         redshift_list: np.ndarray,
                         u_app_mag_list: np.ndarray,
                         g_app_mag_list: np.ndarray,
                         r_app_mag_list: np.ndarray,
                         i_app_mag_list: np.ndarray,
                         Z_app_mag_list: np.ndarray,
                         Y_app_mag_list=np.empty(0),
                         J_app_mag_list=np.empty(0),
                         H_app_mag_list=np.empty(0),
                         Ks_app_mag_list=np.empty(0),
                         u_app_mag_err_list=np.empty(0),
                         g_app_mag_err_list=np.empty(0),
                         r_app_mag_err_list=np.empty(0),
                         i_app_mag_err_list=np.empty(0),
                         Z_app_mag_err_list=np.empty(0),
                         Y_app_mag_err_list=np.empty(0),
                         J_app_mag_err_list=np.empty(0),
                         H_app_mag_err_list=np.empty(0),
                         Ks_app_mag_err_list=np.empty(0)):
    '''
    Calculates maggy and inverse variance values from apparent magnitude and their error values
    and saves the values in a space delimited csv file with columns (without headers):
        
        redshift u_maggy g_maggy r_maggy... u_inv_var g_inv_var r_inv_var...
    
    File is required to be used with the get_rec_maggies_files function 
    or other kcorrect_python functions that best-fit SED coefficients.
    WARNING: pre-existing file with same name will be over-written.
    
    Parameters
    ----------
    obs_maggies_outfile_name : str
        name/path of file with '.csv' extention to save maggies and respective inverse variance values in
    bands : str
        'ugriz' or 'ugriZYJHKs' - refer source code if using other bands
    redshift_list : np.ndarray
        redshift of each data point (galaxy)
    u_app_mag_list : np.ndarray
        all corresponding apparent magnitudes in u-band
    g_app_mag_list : np.ndarray
        all corresponding apparent magnitudes in g-band
    r_app_mag_list : np.ndarray
        all corresponding apparent magnitudes in r-band
    i_app_mag_list : np.ndarray
        all corresponding apparent magnitudes in i-band
    Z_app_mag_list : np.ndarray
        all corresponding apparent magnitudes in Z-band
    Y_app_mag_list : np.ndarray
        all corresponding apparent magnitudes in Y-band
    J_app_mag_list : np.ndarray
        all corresponding apparent magnitudes in J-band
    H_app_mag_list : np.ndarray
        all corresponding apparent magnitudes in H-band
    Ks_app_mag_list : np.ndarray
        all corresponding apparent magnitudes in Ks-band
    u_app_mag_err_list : np.ndarray
        all corresponding errors on apparent magnitudes in u-band
    g_app_mag_err_list : np.ndarray
        all corresponding errors on apparent magnitudes in g-band
    r_app_mag_err_list : np.ndarray
        all corresponding errors on apparent magnitudes in r-band
    i_app_mag_err_list : np.ndarray
        all corresponding errors on apparent magnitudes in i-band
    Z_app_mag_err_list : np.ndarray
        all corresponding errors on apparent magnitudes in Z-band
    Y_app_mag_err_list : np.ndarray
        all corresponding errors on apparent magnitudes in Y-band
    J_app_mag_err_list : np.ndarray
        all corresponding errors on apparent magnitudes in J-band
    H_app_mag_err_list : np.ndarray
        all corresponding errors on apparent magnitudes in H-band
    Ks_app_mag_err_list : np.ndarray
        all corresponding errors on apparent magnitudes in Ks-band
    '''

    if bands == 'ugriz':
        maggy_inv_var_table = np.column_stack(
            (redshift_list, get_maggy(u_app_mag_list),
             get_maggy(g_app_mag_list), get_maggy(r_app_mag_list),
             get_maggy(i_app_mag_list), get_maggy(Z_app_mag_list),
             get_maggy_inv_var(get_maggy(u_app_mag_list), u_app_mag_err_list),
             get_maggy_inv_var(get_maggy(g_app_mag_list), g_app_mag_err_list),
             get_maggy_inv_var(get_maggy(r_app_mag_list), r_app_mag_err_list),
             get_maggy_inv_var(get_maggy(i_app_mag_list), i_app_mag_err_list),
             get_maggy_inv_var(get_maggy(Z_app_mag_list), Z_app_mag_err_list)))
        ascii.write(maggy_inv_var_table,
                    obs_maggies_outfile_name,
                    overwrite=True,
                    format='no_header',
                    names=[
                        'redshift', 'u_maggy', 'g_maggy', 'r_maggy', 'i_maggy',
                        'z_maggy', 'u_inv_var', 'g_inv_var', 'r_inv_var',
                        'i_inv_var', 'z_inv_var'
                    ])
        print(
            '\tRedshifts, and ' + bands +
            ' maggies and their inverse variances calculated, stacked and saved in '
            + obs_maggies_outfile_name + '.')

    elif bands == 'ugriZYJHKs':
        maggy_inv_var_table = np.column_stack(
            (redshift_list, get_maggy(u_app_mag_list),
             get_maggy(g_app_mag_list), get_maggy(r_app_mag_list),
             get_maggy(i_app_mag_list), get_maggy(Z_app_mag_list),
             get_maggy(Y_app_mag_list), get_maggy(J_app_mag_list),
             get_maggy(H_app_mag_list), get_maggy(Ks_app_mag_list),
             get_maggy_inv_var(get_maggy(u_app_mag_list), u_app_mag_err_list),
             get_maggy_inv_var(get_maggy(g_app_mag_list), g_app_mag_err_list),
             get_maggy_inv_var(get_maggy(r_app_mag_list), r_app_mag_err_list),
             get_maggy_inv_var(get_maggy(i_app_mag_list), i_app_mag_err_list),
             get_maggy_inv_var(get_maggy(Z_app_mag_list), Z_app_mag_err_list),
             get_maggy_inv_var(get_maggy(Y_app_mag_list), Y_app_mag_err_list),
             get_maggy_inv_var(get_maggy(J_app_mag_list), J_app_mag_err_list),
             get_maggy_inv_var(get_maggy(H_app_mag_list), H_app_mag_err_list),
             get_maggy_inv_var(get_maggy(Ks_app_mag_list), Ks_app_mag_err_list)))
        ascii.write(maggy_inv_var_table,
                    obs_maggies_outfile_name,
                    overwrite=True,
                    format='no_header',
                    names=[
                        'redshift', 'u_maggy', 'g_maggy', 'r_maggy', 'i_maggy',
                        'Z_maggy', 'Y_maggy', 'J_maggy', 'H_maggy', 'Ks_maggy',
                        'u_inv_var', 'g_inv_var', 'r_inv_var', 'i_inv_var',
                        'Z_inv_var', 'Y_inv_var', 'J_inv_var', 'H_inv_var',
                        'Ks_inv_var'
                    ])
        print(
            '\tRedshifts, and ' + bands +
            ' maggies and their inverse variances calculated, stacked and saved in '
            + obs_maggies_outfile_name + '.')

    else:
        print('\tOnly valid for bands ugriz or ugriZYJHKs.')
        print(
            '\tCheck the source code for basic structure of this function that creates the required file if using other bands.'
        )

# def get_rec_maggies_files(obs_maggies_file_path: str,
#                           n_bands: int,
#                           rec_z_list: np.ndarray,
#                           rec_maggies_outfile_affix='',
#                           survey='sdss',
#                           band_z_shift=0.0,
#                           template_vmatrix_file_path='vmatrix.default.dat',
#                           template_lambda_file_path='lambda.default.dat',
#                           filters_list_file_path='sdss_filters.dat'):
#     '''
#     Reconstructs the observed maggy values at required redshift values
#     by best-fitting galaxy SEDs on data using templates and filter transmission curves,
#     and saves the reconstructed maggy values in a space delimited csv file with columns (without headers):
        
#         redshift rec_u_maggy rec_g_maggy rec_r_maggy...
    
#     File is required to be used with the get_maggy_ratio_file or get_rest_maggy_ratio_file functions.
#     WARNING: pre-existing file with same name will be over-written.
    
#     Parameters
#     ----------
#     obs_maggies_file_path : str
#         path of '.csv' file with the observed maggies and respective inverse variance values. File can be obtained from the get_obs_maggies_file function
#     n_bands : int
#         number of bands used in the survey (and present in the obs_maggies_file)
#     rec_z_list : np.ndarray
#         redshift values required to reconstruct maggies at
#     rec_maggies_outfile_affix : str
#         output file identifier - reconstructed maggies will be saved in 'maggies_at_z[redshift-value]_[identifier].csv'
#     survey : str
#         name of survey being used. Set as 'sdss' by default - do not change if sdss-ugriz are being used
#     band_z_shift : float
#         redshift value to shift the bandpasses/filters by, default is set at 0.0 i.e. no shift
#     template_vmatrix_file_path : str
#         path of '.dat' file with vmatrix of SED templates - must change if survey parameter is not 'sdss'
#     template_lambda_file_path : str
#         path of '.dat' file with lambda of SED templates - must change if survey parameter is not 'sdss'
#     filters_list_file_path : str
#         path of '.dat' file with the list of '.dat' files corresponding to each band and containing its filter transmission curve - must change if survey parameter is not 'sdss'
#     '''

#     if survey == 'sdss':
#         kcorrect.load_templates()
#         print('\tTemplates loaded.')
#         kcorrect.load_filters(band_shift=band_z_shift)
#         print('\tFilters loaded.')
#     else:
#         kcorrect.load_templates(v=template_vmatrix_file_path,
#                                 l=template_lambda_file_path)
#         print('\tTemplates loaded.')
#         kcorrect.load_filters(filters_list_file_path, band_shift=band_z_shift)
#         print('\tFilters loaded.')

#     maggy_inv_var_table = np.genfromtxt(obs_maggies_file_path, delimiter=' ')
#     print('\tRead ' + obs_maggies_file_path + '.')

#     for rec_z in rec_z_list:
#         rec_maggies_outfile_name = 'maggies_at_z' + str(rec_z) + '_' + rec_maggies_outfile_affix + '.csv'
#         rec_maggies_stack = []
#         for i in range(len(maggy_inv_var_table[:, 0])):
#             redshift = maggy_inv_var_table[i, 0]
#             maggies = maggy_inv_var_table[i, 1:(n_bands + 1)]
#             maggies_inv_var = maggy_inv_var_table[i, (n_bands + 1):((2 * n_bands) + 1)]
#             coeffs = kcorrect.fit_nonneg(redshift, maggies, maggies_inv_var)
#             rec_maggies_row = kcorrect.reconstruct_maggies(coeffs, redshift=rec_z)
#             rec_maggies_stack.append(rec_maggies_row)
#         rec_maggies_table = np.array(rec_maggies_stack)
#         ascii.write(rec_maggies_table,
#                     rec_maggies_outfile_name,
#                     overwrite=True,
#                     format='no_header')
#         print('\t' + rec_maggies_outfile_name + ' saved.')
#     print('\tMaggies reconstructed at all redshifts in input array rec_z_list.')

def get_rest_maggy_ratio_file(ID_list: np.ndarray,
                              obs_maggies_file_path: str,
                              rest_maggies_file_path: str,
                              band_index: int,
                              rest_maggy_ratio_outfile_affix=''):
    '''
    Calculates rest-frame maggy ratios i.e. (obs_maggy/rest_maggy),
    and saves the maggy ratio values in a csv file with 3 space delimited columns, of headers:
        
        ID rest_z maggy_ratio
    
    File can be unpacked and used with get_rest_mag function to calculate rest-frame magnitudes.
    WARNING: pre-existing file with same name will be over-written.
    
    Parameters
    ----------
    ID_list: np.ndarray
        ID of each data point (galaxy)
    obs_maggies_file_path : str
        path of '.csv' file with the observed maggies and respective inverse variance values. File can be obtained from the get_obs_maggies_file function
    rest_maggies_file_path : str
        path of '.csv' file with the reconstructed maggies at redshift zero. File can be obtained from the get_rec_maggies_files function by setting rec_z_list to np.array([0.0])
    band_index : int
        band number of required maggy ratio (e.g. 3 for r maggy in ugriz bands)
    rest_maggy_ratio_outfile_affix : str
        output file identifier - rest-frame maggy ratios will be saved in 'rest_maggy_ratios_[identifier].csv'
    '''

    obs_maggies_table = np.genfromtxt(obs_maggies_file_path, delimiter=' ')
    rest_maggies_table = np.genfromtxt(rest_maggies_file_path, delimiter=' ')

    rest_z_list = rest_maggies_table[:, 0]

    obs_maggies_list = obs_maggies_table[:, band_index]
    rest_maggies_list = rest_maggies_table[:, band_index]
    rest_maggy_ratios_list = obs_maggies_list / rest_maggies_list

    rest_maggy_ratio_outfile_name = 'rest_maggy_ratios_' + rest_maggy_ratio_outfile_affix + '.csv'
    rest_maggy_ratios_table = np.column_stack(
        (ID_list, rest_z_list, rest_maggy_ratios_list))
    ascii.write(rest_maggy_ratios_table,
                rest_maggy_ratio_outfile_name,
                overwrite=True,
                names=['ID', 'rest_z', 'maggy_ratio'])
    print('\t' + rest_maggy_ratio_outfile_name + ' created.')

def get_rest_mag(redshift_list: np.ndarray, 
                 app_mag_list: np.ndarray, 
                 maggy_ratio_list: np.ndarray) -> np.ndarray:
    """
    Converts apparent magnitudes into rest-frame magnitudes.
    It uses the apparent magnitudes, redshifts and maggy ratios.
    
    Parameters
    ----------
    redshift_list : np.ndarray
        redshift of each data point (galaxy)
    app_mag_list : np.ndarray
        all corresponding apparent magnitudes
    maggy_ratio_list : np.ndarray
        all corresponding maggy ratios

    Returns
    -------
    np.ndarray
        all corresponding rest-frame magnitudes

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

def get_maggy_ratio_file(ID_list: np.ndarray,
                         rest_maggies_file_path: str,
                         rec_maggies_file_path: str,
                         rec_z: float,
                         band_index: int,
                         maggy_ratio_outfile_affix=''):
    '''
    Calculates reconstructed maggy ratios i.e. (rec_maggy/rest_maggy),
    and saves the maggy ratio values in a csv file with 3 space delimited columns, of headers:
        
        ID rec_z maggy_ratio
    
    WARNING: pre-existing file with same name will be over-written.
    
    Parameters
    ----------
    ID_list: np.ndarray
        ID of each data point (galaxy)
    rest_maggies_file_path : str
        path of '.csv' file with the reconstructed maggies at redshift zero. File can be obtained from the get_rec_maggies_files function by setting rec_z_list to np.array([0.0])
    rec_maggies_file_path : str
        path of '.csv' file with the reconstructed maggies at required reconstruction redshift (rec_z). File can be obtained from the get_rec_maggies_files function by setting rec_z_list to np.array([rec_z])
    rec_z : float
        redshift value where maggies have been reconstruct at
    band_index : int
        band number of required maggy ratio (e.g. 3 for r maggy in ugriz bands)
    rest_maggy_ratio_outfile_affix : str
        output file identifier - maggy ratios will be saved in 'maggy_ratios_at_z[redshift-value]_[identifier].csv'
    '''

    rec_maggies_table = np.genfromtxt(rec_maggies_file_path, delimiter=' ')
    rest_maggies_table = np.genfromtxt(rest_maggies_file_path, delimiter=' ')

    rec_z_list = rec_maggies_table[:, 0]

    rec_maggies_list = rec_maggies_table[:, band_index]
    rest_maggies_list = rest_maggies_table[:, band_index]
    maggy_ratios_list = rec_maggies_list / rest_maggies_list

    maggy_ratio_outfile_name = 'maggy_ratios_at_z' + str(rec_z) + '_' + maggy_ratio_outfile_affix + '.csv'
    maggy_ratios_table = np.column_stack(
        (ID_list, rec_z_list, maggy_ratios_list))
    ascii.write(maggy_ratios_table,
                maggy_ratio_outfile_name,
                overwrite=True,
                names=['ID', 'rec_z', 'maggy_ratio'])
    print('\t' + maggy_ratio_outfile_name + ' saved.')

def get_all_maggy_ratios_file(rec_z_list: np.ndarray, 
                              ID_list: np.ndarray, 
                              band_index: int, 
                              maggies_and_out_files_affix=''):
    '''
    Calculates reconstructed maggy ratios i.e. (rec_maggy/rest_maggy)
    and saves the maggy ratio values at each redshift value in rec_z_list
    in a separate csv file with 3 space delimited columns, of headers:
        
        ID rec_z maggy_ratio
    
    Finally, consolidates all maggy ratios by joining the above files in the order of rec_z_list
    in a single csv file with 3 space delimited columns, of headers:
        
        ID rec_z maggy_ratio

    File with all maggy ratios can be used to calculate z-max.
    WARNING: pre-existing file with same name will be over-written.
    
    Parameters
    ----------
    rec_z_list : np.ndarray
        redshift values where maggies have been reconstruct at - array must have 0.0 redshift value at index 0
    ID_list : np.ndarray
        ID of each data point (galaxy)
    band_index : int
        band number of required maggy ratio (e.g. 3 for r maggy in ugriz bands)
    maggies_and_out_files_affix : str
        output file identifier - values will be saved in 'maggy_ratios_at_z[redshift-value]_[identifier].csv' and 'all_maggy_ratios_[identifier].csv' - must be the same string as rec_maggies_outfile_affix parameter used in get_rec_maggies_files function
    '''

    rest_maggies_file_name = 'maggies_at_z' + str(rec_z_list[0]) + '_' + maggies_and_out_files_affix + '.csv'

    for rec_z in rec_z_list:
        rec_maggies_file_name = 'maggies_at_z' + str(rec_z) + '_' + maggies_and_out_files_affix + '.csv'
        get_maggy_ratio_file(ID_list,
                             rest_maggies_file_name,
                             rec_maggies_file_name,
                             rec_z,
                             band_index,
                             maggy_ratio_outfile_affix=maggies_and_out_files_affix)
    print('\tMaggy ratios calculated at all redshifts in input array rec_z_list.')

    all_maggy_ratios_outfile_name = 'all_maggy_ratios_' + maggies_and_out_files_affix + '.csv'

    rest_maggy_ratio_file_name = 'maggy_ratios_at_z' + str(rec_z_list[0]) + '_' + maggies_and_out_files_affix + '.csv'

    all_maggy_ratios_file = open(all_maggy_ratios_outfile_name, 'w')
    # first file:
    for line in open(rest_maggy_ratio_file_name):
        all_maggy_ratios_file.write(line)
    # now the rest:
    for i in range(len(rec_z_list) - 1):
        maggy_ratio_file_name = 'maggy_ratios_at_z' + str(rec_z_list[i + 1]) + '_' + maggies_and_out_files_affix + '.csv'
        maggy_ratio_file = open(maggy_ratio_file_name)
        maggy_ratio_lines = maggy_ratio_file.readlines()[1:]  # skip the header
        for line in maggy_ratio_lines:
            all_maggy_ratios_file.write(line)
        maggy_ratio_file.close()
    all_maggy_ratios_file.close()

    print('\tAll maggy ratios consolidated in file ' + all_maggy_ratios_outfile_name + '.')
    
def get_volume(survey_area: float, 
               redshift_list: np.ndarray) -> np.ndarray:
    """
    Returns comoving volume of input survey area and redshift.
    
    Parameters
    ----------
    survey_area : float
        survey area in sq. deg. 
    redshift_list : np.ndarray
        redshift of each data point (galaxy)

    Returns
    -------
    np.ndarray
        all corresponding comoving volumes

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

def get_binned_phi(rest_mag_list: np.ndarray, 
                   Vmax_list: np.ndarray, 
                   n_mag_bins: int) -> np.ndarray:
    """
    Bins and weighs galaxy counts per magnitude implementing the 1/Vmax estimator. 
    Returns phi using rest-frame magnitude, maximum observed volume and the number of bins.
    
    Parameters
    ----------
    rest_mag_list : np.ndarray
        rest-frame magnitude of each data point (galaxy)
    Vmax_list : np.ndarray
        all corresponding maximum volumes
    n_mag_bins: int
        number of magnitude bins required
        
    Returns
    -------
    np.ndarray
        mid-magnitude (i.e. x) value of each bin
    np.ndarray
        magnitude-width/2 (i.e. x-error) value of each bin
    np.ndarray
        phi (i.e. y) value of each bin (with h = 0.7)

    """

    # get bin_edges for diving the rest_mags in n_bins
    counts, bin_edges = np.histogram(rest_mag_list, bins=n_mag_bins)

    # sort rest_mag and Vmax lists per increasing mag
    sorted_index = np.argsort(rest_mag_list)
    sorted_Vmax_list = np.array(Vmax_list)[sorted_index]
    sorted_rest_mag_list = np.sort(rest_mag_list)

    # create empty lists for mid_M, phi and M_err
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

def get_patch_centers(uniform_random_RA_list: np.ndarray,
                      uniform_random_DEC_list: np.ndarray,
                      n_patches: int,
                      survey='kids',
                      max_iterations=int(100),
                      tolerance=1.0e-5,
                      patch_centers_outfile_affix=''):
    """
    Divides the input uniform random survey into equally distributed and equally sized patches. 
    Calculates n_patches centers [RA,Dec] from RA, Dec and number of patches and saves in a csv file 
    with 2 space delimited columns (without headers):
        
        RA Dec

    Function does not overwrite any existing file with the same name. File need not be updated with every run.

    Parameters
    ----------
    uniform_random_RA_list : np.ndarray
        RA values of each data point (galaxy) in a uniform random catalogue
    uniform_random_DEC_list : np.ndarray
        all corresponding Dec values in the uniform random catalogue
    n_patches : int
        number of equal survey area patches required
    survey : str, optional
        survey name - only change if survey area covers/connects over 320 degrees RA and does not connect over 360 to 0 degrees RA 
    max_iterations : int, optional
        maximum number of iterations to run
    tolerance : float, optional
        relative change in the average distance to centers, signifies convergence
    patch_centers_outfile_affix : str
        output file identifier - values will be saved in 'patch_centers_tol[tolerance]_[identifier].csv'
    """

    # MAKE SURE ALL PATCHES ARE SITCHED ON SKY
    # works for most surveys - GAMA, KiDS - check rest
    if survey == 'kids':
        corrected_uniform_random_RA_list = np.where(
            uniform_random_RA_list > 320., uniform_random_RA_list - 360.,
            uniform_random_RA_list)
    # use if a survey patch covers/connects over 320 degrees RA
    # and does not connect over 360 to 0 degree RA
    if survey != 'kids':
        corrected_uniform_random_RA_list = uniform_random_RA_list

    # STACK RA AND DEC AS uniform_random_X
    uniform_random_X = np.column_stack(
        (corrected_uniform_random_RA_list, uniform_random_DEC_list))

    # DIVIDE uniform_random_X INTO EQUAL n_patches
    uniform_random_km = kmeans_sample(uniform_random_X,
                                      n_patches,
                                      maxiter=max_iterations,
                                      tol=tolerance)
    center_guesses = uniform_random_km.centers
    ra_guesses = center_guesses[:, 0]
    dec_guesses = center_guesses[:, 1]
    centers_table = np.column_stack((ra_guesses, dec_guesses))

    patch_centers_outfile_name = 'patch_centers_tol' + str(tolerance) + '_' + patch_centers_outfile_affix + '.csv'
    ascii.write(centers_table,
                patch_centers_outfile_name,
                overwrite=False,
                format='no_header')
    print('Patch center guesses saved in '+ patch_centers_outfile_name)

def get_patch_labels(RA_list: np.ndarray,
                     DEC_list: np.ndarray,
                     n_patches: int,
                     patch_centers_file_path: str,
                     survey='kids',
                     numba_installed=True,
                     plot_savename='none') -> np.ndarray:
    """
    Divides survey into equally distributed and equally sized patches. Returns labels for patches from RA, Dec, number of patches and patch center guesses file.
    WARNING: does not display plot, must specify plot_savename parameter to save plot

    Parameters
    ----------
    RA_list : np.ndarray
        RA values of each data point (galaxy)
    DEC_list : np.ndarray
        all corresponding Dec values
    n_patches : int
        number of equal survey area patches required
    patch_centers_file_path : str
        path of '.csv' file with (n_patches x 2) patch center guesses (RA, Dec). File can be obtained from the get_patch_centers function
    survey : str, optional
        survey name - only change if survey area covers/connects over 320 degrees RA and does not connect over 360 to 0 degrees RA 
    numba_installed : bool, optional
        mark as False if numba is not installed
    plot_savename : str, optional
        name and extension to save plot as, plot will not be saved if not changed
        
    Returns
    -------
    np.ndarray
        array of patch assignment label for each data point

    """

    # MAKE SURE ALL PATCHES ARE STITCHED ON SKY
    # works for most surveys - GAMA, KiDS - check rest
    if survey == 'kids':
        corrected_RA_list = np.where(RA_list > 320., RA_list - 360., RA_list)
    # use if a survey patch covers/connects over 320 degrees RA
    # and does not connect over 360 to 0 degree RA
    if survey != 'kids':
        corrected_RA_list = RA_list

    # STACK RA AND DEC AS X
    X = np.column_stack((corrected_RA_list, DEC_list))

    #UNPACK PATCH CENTER GUESSES
    centers_table = np.genfromtxt(patch_centers_file_path, delimiter=' ')
    ra_guesses = centers_table[ : , 0]
    dec_guesses = centers_table[ : , 1]
    center_guesses = np.column_stack((ra_guesses, dec_guesses))

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

    return labels

def get_binned_phi_error(rest_mag_list: np.ndarray, 
                         Vmax_list: np.ndarray, 
                         labels: np.ndarray, 
                         n_patches: int, 
                         n_mag_bins: int) -> np.ndarray:
    """
    Spatial variance on galaxy number density per magnitude. 
    Returns error on phi from rest-frame magnitude, maximum observed volume, labels, number of patches and number of bins.
    
    Parameters
    ----------
    rest_mag_list : np.ndarray
        rest-frame magnitude of each data point (galaxy)
    Vmax_list : np.ndarray
        all corresponding maximum volumes
    labels : np.ndarray
        all corresponding survey patch assignment labels
    n_patches : int
        number of equal survey area patches required
    n_mag_bins : int
        number of magnitude bins required
    
    Returns
    -------
    np.ndarray
        phi error (i.e. y-error) value of each bin
        
    """

    # GET PHI VALUES USING ONLY VALUES IN EACH PATCH
    patch_phis = []
    for n in range(n_patches):
        patch_indices = np.where(labels == n)
        patch_M = rest_mag_list[patch_indices]
        patch_Vmax = Vmax_list[patch_indices] / n_patches
        mid_M_list, M_err_list, phi_list = get_binned_phi(
            patch_M, patch_Vmax, n_mag_bins)
        patch_phis.append(phi_list)

    # STANDARD ERRORS ON PHI VALUES BETWEEN EACH PATCH
    phi_err_list = np.std(patch_phis, axis=0)

    return phi_err_list

def get_plot(rest_mag_list: np.ndarray,
             Vmax_list: np.ndarray,
             n_mag_bins: int,
             RA_list: np.ndarray,
             DEC_list: np.ndarray,
             n_patches: int,
             patch_centers_file_path: str,
             survey='kids',
             numba_installed=True, 
             plot_savename='none') -> np.ndarray:
    """
    Plots the 1/Vmax weighted luminosity function from data, binned by magnitude.
    WARNING: does not display plot, must specify plot_savename parameter to save plot
    
    Parameters
    ----------
    rest_mag_list : np.ndarray
        rest-frame magnitude of each data point (galaxy)
    Vmax_list : np.ndarray
        all corresponding maximum volumes
    n_mag_bins : int
        number of magnitude bins required
    RA_list : np.ndarray
        all corresponding RA values
    DEC_list : np.ndarray
        all corresponding Dec values
    n_patches : int
        number of equal survey area patches required
    patch_centers_file_path : str
        path of '.csv' file with (n_patches x 2) patch center guesses (RA, Dec). File can be obtained from the get_patch_centers function
    survey : str, optional
        survey name - only change if survey area covers/connects over 320 degrees RA and does not connect over 360 to 0 degrees RA
    numba_installed : bool, optional
        mark as False if numba is not installed
    plot_savename : str, optional
        name and extension to save plot as, plot will not be saved if not changed
    
    Returns
    -------
    np.ndarray
        mid-magnitude (i.e. x) value of each bin
    np.ndarray
        magnitude-width/2 (i.e. x-error) value of each bin
    np.ndarray
        phi (i.e. y) value of each bin (with h = 0.7)
    np.ndarray
        phi error (i.e. y-error) value of each bin
        
    """

    # phi
    M_list, M_err_list, phi_list = get_binned_phi(rest_mag_list, Vmax_list, n_mag_bins)
    # patches
    labels = get_patch_labels(RA_list, DEC_list, n_patches, patch_centers_file_path, survey, numba_installed)
    # phi errors
    phi_err_list = get_binned_phi_error(rest_mag_list, Vmax_list, labels, n_patches, n_mag_bins)

    if plot_savename != 'none':
        plt.figure(figsize=(10, 10))

        # plot data
        plt.errorbar(M_list,
                     phi_list,
                     xerr=M_err_list,
                     yerr=phi_err_list,
                     fmt='gx',
                     mec='k',
                     label='galaxies:' + str(len(rest_mag_list)))

        plt.yscale('log')
        # plt.xlim(-26,-12)
        # plt.ylim(1e-8,0.9)
        plt.xlabel("rest-frame magnitude/ $(M_{r})_{cal}$/ mag", fontsize=20)
        plt.ylabel(
            "number density / $\Phi (M_{r})/ h_{70}^{3}Mpc^{-3}mag^{-1}$",
            fontsize=20)
        # plt.title(plot_savename, fontsize=20)

        plt.grid(True)
        plt.legend(loc='upper left')

        plt.savefig(plot_savename, dpi=300)

    return M_list, M_err_list, phi_list, phi_err_list

def filter_plot_by_colour(dichotomy_slope: float,
                          dichotomy_intercept: float,
                          rest_mag_list: np.ndarray,
                          higher_band_rest_mag_list: np.ndarray,
                          Vmax_list: np.ndarray,
                          n_mag_bins: int,
                          RA_list: np.ndarray,
                          DEC_list: np.ndarray,
                          n_patches: int,
                          patch_centers_file_path: str,
                          survey='kids',
                          numba_installed=True,
                          plot_savename='none') -> np.ndarray:
    """
    Plots the 1/Vmax weighted luminosity function from data, binned by magnitude and filtered by galaxy colours. The galaxy colours are filtered by red and blue with the help of the input colour dichotomy line parameters. The colour dichotomy line parameters can be inferred from a CMD plot.
    WARNING: does not display plot, must specify plot_savename parameter to save plot

    Parameters
    ----------
    dichotomy_slope : float
        slope of the colour dichotomy line
    dichotomy_intercept : float
        intercept of the colour dichotomy line
    rest_mag_list : np.ndarray
        rest-frame magnitude of each data point (galaxy)
    higher_band_rest_mag_list : np.ndarray
        rest-frame magnitudes of each data point (galaxy) from a higher wavelength band
    Vmax_list : np.ndarray
        all corresponding maximum volumes
    n_mag_bins : int
        number of magnitude bins required
    RA_list : np.ndarray
        all coressponding RA values
    DEC_list : np.ndarray
        all corresponding Dec values
    n_patches : int
        number of patches required
    patch_centers_file_path : str
        path of '.csv' file with (n_patches x 2) patch center guesses (RA, Dec). File can be obtained from the get_patch_centers function
    survey : str, optional
        survey name - only change if survey area covers/connects over 320 degrees RA and does not connect over 360 to 0 degrees RA
    numba_installed : bool, optional
        mark as False if numba is not installed
    plot_savename : str, optional
        name and extension to save plot as, plot will not be saved if not changed

    Returns
    -------
    np.ndarray
        all galaxies' LF's mid-magnitude (i.e. x) value of each bin
    np.ndarray
        all galaxies' LF's magnitude-width/2 (i.e. x-error) value of each bin
    np.ndarray
        all galaxies' LF's phi (i.e. y) value of each bin (with h = 0.7)
    np.ndarray
        all galaxies' LF's phi error (i.e. y-error) value of each bin
    np.ndarray
        red galaxies' LF's mid-magnitude (i.e. x) value of each bin
    np.ndarray
        red galaxies' LF's magnitude-width/2 (i.e. x-error) value of each bin
    np.ndarray
        red galaxies' LF's phi (i.e. y) value of each bin (with h = 0.7)
    np.ndarray
        red galaxies' LF's phi error (i.e. y-error) value of each bin
    np.ndarray
        blue galaxies' LF's mid-magnitude (i.e. x) value of each bin
    np.ndarray
        blue galaxies' LF's magnitude-width/2 (i.e. x-error) value of each bin
    np.ndarray
        blue galaxies' LF's phi (i.e. y) value of each bin (with h = 0.7)
    np.ndarray
        blue galaxies' LF's phi error (i.e. y-error) value of each bin
        
    """

    colour_mag_list = higher_band_rest_mag_list - rest_mag_list
    dichotomy_line = dichotomy_slope * rest_mag_list + dichotomy_intercept
    red_index = np.where(colour_mag_list >= dichotomy_line)[0]
    blue_index = np.where(colour_mag_list < dichotomy_line)[0]

    # all
    M_list, M_err_list, phi_list, phi_err_list = get_plot(
        rest_mag_list, Vmax_list, n_mag_bins, RA_list, DEC_list, n_patches,
        patch_centers_file_path, survey, numba_installed)

    # red
    red_M_list, red_M_err_list, red_phi_list, red_phi_err_list = get_plot(
        rest_mag_list[red_index], Vmax_list[red_index], n_mag_bins,
        RA_list[red_index], DEC_list[red_index], n_patches, patch_centers_file_path,
        survey, numba_installed)

    # blue
    blue_M_list, blue_M_err_list, blue_phi_list, blue_phi_err_list = get_plot(
        rest_mag_list[blue_index], Vmax_list[blue_index], n_mag_bins,
        RA_list[blue_index], DEC_list[blue_index], n_patches, patch_centers_file_path,
        survey, numba_installed)

    if plot_savename != 'none':
        plt.figure(figsize=(10, 10))

        # plot all data
        plt.errorbar(M_list,
                     phi_list,
                     xerr=M_err_list,
                     yerr=phi_err_list,
                     fmt='gx',
                     mec='k',
                     label='all:' + str(len(rest_mag_list)))

        # plot red data
        plt.errorbar(red_M_list,
                     red_phi_list,
                     xerr=red_M_err_list,
                     yerr=red_phi_err_list,
                     fmt='rx',
                     mec='k',
                     label='red:' + str(len(rest_mag_list[red_index])))

        # plot blue data
        plt.errorbar(blue_M_list,
                     blue_phi_list,
                     xerr=blue_M_err_list,
                     yerr=blue_phi_err_list,
                     fmt='bx',
                     mec='k',
                     label='blue:' + str(len(rest_mag_list[blue_index])))

        plt.yscale('log')
        # plt.xlim(-26,-12)
        # plt.ylim(1e-8,0.9)
        plt.xlabel("rest-frame r-magnitude/ $(M_{r})_{cal}$/ mag", fontsize=20)
        plt.ylabel(
            "number density / $\Phi (M_{r})/ h_{70}^{3}Mpc^{-3}mag^{-1}$",
            fontsize=20)
        # plt.title(plot_savename, fontsize=20)

        plt.grid(True)
        plt.legend(loc='upper left')

        plt.savefig(plot_savename, dpi=300)

    return M_list, M_err_list, phi_list, phi_err_list, red_M_list, red_M_err_list, red_phi_list, red_phi_err_list, blue_M_list, blue_M_err_list, blue_phi_list, blue_phi_err_list

def SchechterMagModel(M_list: np.ndarray, 
                      M_star: float, 
                      phi_star: float, 
                      alpha: float) -> np.ndarray:
    """
    Single Schechter luminosity function in terms of magnitude from 3 free parameters of the model.
    
    Parameters
    ----------
    M_list : np.ndarray
        array of magnitudes (i.e. x)
    M_star : float
        model parameter M_star
    phi_star : float 
        model parameter phi_star
    alpha : float
        model parameter alpha

    Returns
    -------
    np.ndarray
        array of Schechter modelled phi (i.e. y)
        
    """

    # FACTOR
    factor = (2 / 5) * np.log(10)

    # POWER
    Mstar_Mlist = M_star - M_list
    power = (2 / 5) * Mstar_Mlist

    # PART 1
    power1 = -10**(power)
    part1 = np.exp(power1)

    # PART 2
    index = alpha + 1
    power2 = power * index
    part2 = phi_star * 10**(power2)

    # PHI(M)
    phi_list = factor * part1 * part2

    return phi_list

def DoubleSchechterMagModel(M_list: np.ndarray, 
                            M_star: float, 
                            phi_star1: float, 
                            alpha1: float, 
                            phi_star2: float, 
                            alpha2: float) -> np.ndarray:
    """
    Double Schechter luminosity function in terms of magnitude from 5 free parameters of the model.
    
    Parameters
    ----------
    M_list : np.ndarray
        array of magnitudes (i.e. x)
    M_star : float
        model parameter M_star
    phi_star1 : float 
        model parameter phi_star1
    alpha1 : float
        model parameter alpha1
    phi_star2 : float 
        model parameter phi_star2
    alpha2 : float
        model parameter alpha2

    Returns
    -------
    np.ndarray
        array of Double Schechter modelled phi (i.e. y)
        
    """

    # FACTOR
    factor = (2 / 5) * np.log(10)

    # POWER
    Mstar_Mlist = M_star - M_list
    power = (2 / 5) * Mstar_Mlist

    # PART 1
    power1 = -10**(power)
    part1 = np.exp(power1)

    # PART 2
    index1 = alpha1 + 1
    power2 = power * index1
    part2 = phi_star1 * 10**(power2)

    # PART 3
    index2 = alpha2 + 1
    power3 = power * index2
    part3 = phi_star2 * 10**(power3)

    # PHI(M)
    phi_list = factor * part1 * (part2 + part3)

    return phi_list

def get_gof(obs: np.ndarray, 
            err: np.ndarray, 
            exp: np.ndarray, 
            m: int) -> float:
    """
    Returns reduced chi squared estimate of goodness of fit.
    
    Parameters
    ----------
    obs : np.ndarray
        observed values (e.g. phi from survey data)
    err : np.ndarray
        errors on observed values
    exp : np.ndarray
        expected values (e.g. phi from the Schechter function)
    m : int
        number of parameters used to calculate the expected values
    
    Returns
    -------
    np.ndarray
        reduced chi square
        
    """

    residuals = obs - exp
    rBYerr = residuals / err
    rBYerr_sq = rBYerr**2
    chi_sq = np.sum(rBYerr_sq)

    dof = len(obs) - m

    red_chi_sq = chi_sq / dof

    return red_chi_sq

def get_schechter_phi(M_list: np.ndarray, 
                      M_err_list: np.ndarray, 
                      phi_list: np.ndarray, 
                      phi_err_list: np.ndarray, 
                      guesses: np.ndarray, 
                      plot_savename='none') -> Tuple[np.ndarray, float, float, float, float, float, float, float]:
    """
    Least square fits single Schechter function model on data.
    Returns best fit phi, reduced chi squared estimate and the 3 Schechter parameters with their errors.
    WARNING: does not display plot, must specify plot_savename parameter to save plot
    
    Parameters
    ----------
    M_list : np.ndarray
        mid magnitude (i.e. x) value of each bin
    M_err_list : np.ndarray
        magnitudes error (i.e. x-error) value of each bin
    phi_list : np.ndarray
        phi (i.e. y) value of each bin
    phi_err_list : np.ndarray
        phi error (i.e. y-error) value of each bin
    guesses : np.ndarray
        array of Schechter parameter guesses in order [M_star, phi_star, aplha]
    plot_savename : str, optional
        name and extension to save plot as, plot will not be saved if not changed

    Returns
    -------
    np.ndarray
        Schechter modelled phi (i.e. y) of each bin
    float
        reduced chi square of the fit
    float
        fit parameter M_star
    float
        error on fit parameter M_star
    float
        fit parameter phi_star
    float
        error on fit parameter phi_star
    float
        fit parameter alpha
    float
        error on fit parameter alpha
        
    """

    popt, pcov = curve_fit(SchechterMagModel,
                           M_list,
                           phi_list,
                           p0=guesses,
                           sigma=phi_err_list)

    perr = np.sqrt(np.diag(pcov))

    M_star = popt[0]
    M_star_err = perr[0]
    phi_star = popt[1]
    phi_star_err = perr[1]
    alpha = popt[2]
    alpha_err = perr[2]

    model_phi_list = SchechterMagModel(M_list, M_star, phi_star, alpha)

    m = 3
    red_chi_sq = get_gof(phi_list, phi_err_list, model_phi_list, m)

    if plot_savename != 'none':

        plt.figure(figsize=(10, 10))

        # plot data
        plt.errorbar(M_list,
                     phi_list,
                     xerr=M_err_list,
                     yerr=phi_err_list,
                     fmt='yx',
                     mec='k',
                     label='Survey data')

        # plot model
        plt.plot(
            M_list,
            model_phi_list,
            'g--',
            label='Schechter, alpha: {0:.4f} $\pm$ {1:.4f}, $\chi^{2}$: {3:.4f}'
            .format(alpha, alpha_err, 2, red_chi_sq))

        # plot turning point
        plt.errorbar(
            M_star,
            phi_star,
            xerr=M_star_err,
            yerr=phi_star_err,
            fmt='c*',
            mec='b',
            label=
            '$M^{0}$: {1:.4f} $\pm$ {2:.4f}, $log\Phi^{3}$: {4:.4f} $\pm$ {5:.4f}'
            .format('*', M_star, M_star_err, '*', np.log10(phi_star),
                    np.log10(phi_star_err)))

        plt.yscale('log')
        # plt.xlim(-26, -12)
        # plt.ylim(1e-8, 0.9)
        plt.xlabel("rest-frame magnitude/ $(M)_{cal}$/ mag", fontsize=20)
        plt.ylabel("number density / $\Phi (M)/ h_{70}^{3}Mpc^{-3}mag^{-1}$",
                   fontsize=20)
        # plt.title(plot_savename, fontsize=20)
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.savefig(plot_savename, dpi=300)

    return model_phi_list, red_chi_sq, M_star, M_star_err, phi_star, phi_star_err, alpha, alpha_err

def get_double_schechter_phi(M_list: np.ndarray, 
                             M_err_list: np.ndarray, 
                             phi_list: np.ndarray, 
                             phi_err_list: np.ndarray, 
                             guesses: np.ndarray, 
                             plot_savename='none') -> Tuple[np.ndarray, float, float, float, float, float, float, float, float, float, float, float]:
    """
    Least square fits double Schechter function model on data.
    Returns best fit phi, reduced chi squared estimate and the 5 Schechter parameters with their errors.
    WARNING: does not display plot, must specify plot_savename parameter to save plot
    
    Parameters
    ----------
    M_list : np.ndarray
        mid magnitude (i.e. x) value of each bin
    M_err_list : np.ndarray
        magnitudes error (i.e. x-error) value of each bin
    phi_list : np.ndarray
        phi (i.e. y) value of each bin
    phi_err_list : np.ndarray
        phi error (i.e. y-error) value of each bin
    guesses : np.ndarray
        array of Schechter parameter guesses in order [M_star, phi_star, aplha]
    plot_savename : str, optional
        name and extension to save plot as, plot will not be saved if not changed

    Returns
    -------
    np.ndarray
        Schechter modelled phi (i.e. y) of each bin
    float
        reduced chi square of the fit
    float
        fit parameter M_star 
    float
        error on fit parameter M_star
    float
        fit parameter phi_star_1
    float
        error on fit parameter phi_star_1
    float
        fit parameter alpha_1
    float
        error on fit parameter alpha_1
    float
        fit parameter phi_star_2
    float
        error on fit parameter phi_star_2
    float
        fit parameter alpha_2
    float
        error on fit parameter alpha_2
        
    """

    popt, pcov = curve_fit(DoubleSchechterMagModel,
                           M_list,
                           phi_list,
                           p0=guesses,
                           sigma=phi_err_list)

    perr = np.sqrt(np.diag(pcov))

    M_star = popt[0]
    M_star_err = perr[0]
    phi_star_1 = popt[1]
    phi_star_err_1 = perr[1]
    alpha_1 = popt[2]
    alpha_err_1 = perr[2]
    phi_star_2 = popt[3]
    phi_star_err_2 = perr[3]
    alpha_2 = popt[4]
    alpha_err_2 = perr[4]
    model_phi_list = DoubleSchechterMagModel(M_list, M_star, phi_star_1,
                                             alpha_1, phi_star_2, alpha_2)

    m = 3
    red_chi_sq = get_gof(phi_list, phi_err_list, model_phi_list, m)

    if plot_savename != 'none':

        plt.figure(figsize=(10, 10))

        # plot data
        plt.errorbar(M_list,
                     phi_list,
                     xerr=M_err_list,
                     yerr=phi_err_list,
                     fmt='yx',
                     mec='k',
                     label='Survey data')

        # plot model
        plt.plot(M_list,
                 model_phi_list,
                 'g--',
                 label='Double Schechter, $\chi^{0}$: {1:.4f}'.format(
                     2, red_chi_sq))

        # plot turning point 1
        plt.errorbar(
            M_star,
            phi_star_1,
            xerr=M_star_err,
            yerr=phi_star_err_1,
            fmt='m*',
            mec='r',
            label=
            '$M^{0}$: {1:.2f} $\pm$ {2:.2f}, $log\Phi_{3}^{4}$: {5:.2f} $\pm$ {6:.2f}, alpha$_{7}$: {8:.2f} $\pm$ {9:.2f}'
            .format('*', M_star, M_star_err, 1, '*', np.log10(phi_star_1),
                    np.log10(phi_star_err_1), 1, alpha_1, alpha_err_1))

        # plot turning point 2
        plt.errorbar(
            M_star,
            phi_star_2,
            xerr=M_star_err,
            yerr=phi_star_err_2,
            fmt='c*',
            mec='b',
            label=
            '$M^{0}$: {1:.2f} $\pm$ {2:.2f}, $log\Phi_{3}^{4}$: {5:.2f} $\pm$ {6:.2f}, alpha$_{7}$: {8:.2f} $\pm$ {9:.2f}'
            .format('*', M_star, M_star_err, 2, '*', np.log10(phi_star_2),
                    np.log10(phi_star_err_2), 2, alpha_2, alpha_err_2))

        plt.yscale('log')
        # plt.xlim(-26, -12)
        # plt.ylim(1e-8, 0.9)
        plt.xlabel("rest-frame magnitude/ $(M)_{cal}$/ mag", fontsize=20)
        plt.ylabel("number density / $\Phi (M)/ h_{70}^{3}Mpc^{-3}mag^{-1}$",
                   fontsize=20)
        # plt.title(plot_savename, fontsize=20)
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.savefig(plot_savename, dpi=300)

    return model_phi_list, red_chi_sq, M_star, M_star_err, phi_star_1, phi_star_err_1, alpha_1, alpha_err_1, phi_star_2, phi_star_err_2, alpha_2, alpha_err_2
