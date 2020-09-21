"""This file tests the following functions:

    * get_plot
    * filter_plot_by_colour
    
"""

# -----------------------
# Package Imports
# -----------------------
import math 
import numpy as np
import pandas as pd
# import cv2
# import pytest
import sys
sys.path.insert(1, '../lumfunc/')
import lumfunc as lf

# -----------------------
# Unpack Test Data
# -----------------------
# test data (photometric galaxian survey)
data_table = pd.read_csv('__tests__/catalogue_test.csv')

ID_list = np.array(data_table['ID'])

RA_list = np.array(data_table['RA'])
Dec_list = np.array(data_table['Dec'])

u_app_mag_list = np.array(data_table['u_mag'])
u_app_mag_err_list = np.array(data_table['u_mag_err'])
g_app_mag_list = np.array(data_table['g_mag'])
g_app_mag_err_list = np.array(data_table['g_mag_err'])
r_app_mag_list = np.array(data_table['r_mag'])
r_app_mag_err_list = np.array(data_table['r_mag_err'])
i_app_mag_list = np.array(data_table['i_mag'])
i_app_mag_err_list = np.array(data_table['i_mag_err'])
Z_app_mag_list = np.array(data_table['Z_mag'])
Z_app_mag_err_list = np.array(data_table['Z_mag_err'])
Y_app_mag_list = np.array(data_table['Y_mag'])
Y_app_mag_err_list = np.array(data_table['Y_mag_err'])
J_app_mag_list = np.array(data_table['J_mag'])
J_app_mag_err_list = np.array(data_table['J_mag_err'])
H_app_mag_list = np.array(data_table['H_mag'])
H_app_mag_err_list = np.array(data_table['H_mag_err'])
K_app_mag_list = np.array(data_table['K_mag'])
K_app_mag_err_list = np.array(data_table['K_mag_err'])

z_photo_list = np.array(data_table['z_photo'])
z_spec_list = np.array(data_table['z_spec'])


# -----------------------
# Test Methods
# -----------------------

# def test_get_plot( ):
#     r_maggy_ratios_table = pd.read_csv('__tests__/rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
#     r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
#     r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)

#     zmax_table = pd.read_csv('__tests__/zmax_test.csv', delimiter=' ')
#     z_max_list = np.array(zmax_table['zmax'])
#     survey_area = 2.5 #sq. degrees
#     Vmax_list = lf.get_volume(survey_area, z_max_list)
    
#     n_bins = 10
    
#     n_patches = 10
#     ugriz_test_patch_centers_file_path = '__tests__/patch_centers_tol0.01_ugriz_test.csv'
#     centers_table = np.genfromtxt(ugriz_test_patch_centers_file_path, delimiter=' ')
#     ra_guesses = centers_table[ : , 0]
#     dec_guesses = centers_table[ : , 1]
#     ugriz_test_patch_centers_guesses = np.column_stack((ra_guesses, dec_guesses))

#     plot_M_list, plot_M_err_list, plot_phi_list, plot_phi_err_list = lf.get_plot(
#         r_rest_mag_list,
#         Vmax_list,
#         n_bins,
#         RA_list,
#         Dec_list,
#         n_patches,
#         ugriz_test_patch_centers_guesses,
#         survey='kids',
#         numba_installed=True,
#         plot_savename='pytest_LF.png')

#     assert list(plot_M_list) == [-24.62894308598833, -23.404512811278693, -22.18008253656906, -20.955652261859427, -19.731221987149794, 
#         -18.506791712440158, -17.282361437730522, -16.05793116302089, -14.833500888311256, -13.609070613601622]
#     assert list(plot_M_err_list) == [0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548164, 0.6122151373548181, 
#         0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548173, 0.6122151373548164]
#     assert list(plot_phi_list) == [290.491672773112, 265.7977860159212, 9.55747321361586e-05, 0.0002549444468315998, 0.0006247531885987247, 
#         0.001075916507027022, 0.0019105283915257355, 0.005624556120539741, 0.0038603784240661696, 0.06417684969053933]
#     assert list(plot_phi_err_list) == [810.9397645708841, 607.8170004077618, 4.36417468536423e-05, 0.00019704012356457505, 
#         0.0005486180652616492, 0.0004654318611417406, 0.0005773328573042754, 0.00459036071996224, 0.00221037276729238, 0.164362438153455]
    
#     # test_result = cv2.imread('__tests__/test_LF.png')
#     # pytest_result = cv2.imread('__tests__/pytest_LF.png')
#     # difference = cv2.subtract(test_result, pytest_result)
#     # b, g, r = cv2.split(difference)
#     # assert test_result.shape == pytest_result.shape
#     # assert cv2.countNonZero(b) == 0
#     # assert cv2.countNonZero(g) == 0
#     # assert cv2.countNonZero(r) == 0

def test_get_plot_no_plot( ):
    r_maggy_ratios_table = pd.read_csv('__tests__/rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)

    zmax_table = pd.read_csv('__tests__/zmax_test.csv', delimiter=' ')
    z_max_list = np.array(zmax_table['zmax'])
    survey_area = 2.5 #sq. degrees
    Vmax_list = lf.get_volume(survey_area, z_max_list)
    
    n_bins = 10
    
    n_patches = 10
    ugriz_test_patch_centers_file_path = '__tests__/patch_centers_tol0.01_ugriz_test.csv'
    centers_table = np.genfromtxt(ugriz_test_patch_centers_file_path, delimiter=' ')
    ra_guesses = centers_table[ : , 0]
    dec_guesses = centers_table[ : , 1]
    ugriz_test_patch_centers_guesses = np.column_stack((ra_guesses, dec_guesses))

    plot_M_list, plot_M_err_list, plot_phi_list, plot_phi_err_list = lf.get_plot(
        r_rest_mag_list,
        Vmax_list,
        n_bins,
        RA_list,
        Dec_list,
        n_patches,
        ugriz_test_patch_centers_guesses,
        survey='kids',
        numba_installed=True)

    assert list(plot_M_list) == [-24.62894308598833, -23.404512811278693, -22.18008253656906, -20.955652261859427, -19.731221987149794, 
        -18.506791712440158, -17.282361437730522, -16.05793116302089, -14.833500888311256, -13.609070613601622]
    assert list(plot_M_err_list) == [0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548164, 0.6122151373548181, 
        0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548173, 0.6122151373548164]
    assert list(plot_phi_list) == [290.491672773112, 265.7977860159212, 9.55747321361586e-05, 0.0002549444468315998, 0.0006247531885987247, 
        0.001075916507027022, 0.0019105283915257355, 0.005624556120539741, 0.0038603784240661696, 0.06417684969053933]
    assert list(plot_phi_err_list) == [810.9397645708841, 607.8170004077618, 4.36417468536423e-05, 0.00019704012356457505, 
        0.0005486180652616492, 0.0004654318611417406, 0.0005773328573042754, 0.00459036071996224, 0.00221037276729238, 0.164362438153455]
# -----------------------

# def test_filter_plot_by_colour( ):
#     r_maggy_ratios_table = pd.read_csv('__tests__/rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
#     r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
#     r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)

#     zmax_table = pd.read_csv('__tests__/zmax_test.csv', delimiter=' ')
#     z_max_list = np.array(zmax_table['zmax'])
#     survey_area = 2.5 #sq. degrees
#     Vmax_list = lf.get_volume(survey_area, z_max_list)
    
#     n_bins = 10
    
#     n_patches = 10
#     ugriz_test_patch_centers_file_path = '__tests__/patch_centers_tol0.01_ugriz_test.csv'
#     centers_table = np.genfromtxt(ugriz_test_patch_centers_file_path, delimiter=' ')
#     ra_guesses = centers_table[ : , 0]
#     dec_guesses = centers_table[ : , 1]
#     ugriz_test_patch_centers_guesses = np.column_stack((ra_guesses, dec_guesses))

#     g_maggy_ratios_table = pd.read_csv('__tests__/rest_maggy_ratios_g_ugriz_test.csv', delimiter=' ')
#     g_maggy_ratio_list = np.array(g_maggy_ratios_table['maggy_ratio'])

#     g_rest_mag_list = lf.get_rest_mag(z_photo_list, g_app_mag_list, g_maggy_ratio_list)

#     colour_cut_slope = 0.0
#     colour_cut_intercept = 0.65
#     all_M_list, all_M_err_list, all_phi_list, all_phi_err_list, red_M_list, red_M_err_list, red_phi_list, red_phi_err_list, blue_M_list, blue_M_err_list, blue_phi_list, blue_phi_err_list = lf.filter_plot_by_colour(
#         colour_cut_slope,
#         colour_cut_intercept,
#         r_rest_mag_list,
#         g_rest_mag_list,
#         Vmax_list,
#         n_bins,
#         RA_list,
#         Dec_list,
#         n_patches,
#         ugriz_test_patch_centers_guesses,
#         survey='kids',
#         numba_installed=True,
#         plot_savename='pytest_LF_colour.png')
    
#     assert list(all_M_list) == [-24.62894308598833, -23.404512811278693, -22.18008253656906, -20.955652261859427, -19.731221987149794, 
#         -18.506791712440158, -17.282361437730522, -16.05793116302089, -14.833500888311256, -13.609070613601622]
#     assert list(all_M_err_list) == [0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548164, 0.6122151373548181,     
#         0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548173, 0.6122151373548164]
#     assert list(all_phi_list) == [290.491672773112, 265.7977860159212, 9.55747321361586e-05, 0.0002549444468315998, 0.0006247531885987247, 
#         0.001075916507027022, 0.0019105283915257355, 0.005624556120539741, 0.0038603784240661696, 0.06417684969053933]
#     assert list(all_phi_err_list) == [810.9397645708841, 607.8170004077618, 4.36417468536423e-05, 0.00019704012356457505, 
#         0.0005486180652616492, 0.0004654318611417406, 0.0005773328573042754, 0.00459036071996224, 0.00221037276729238, 0.164362438153455]
#     assert list(red_M_list) == [-23.749705414732908, -23.22313054306747, -22.696555671402038, -22.1699807997366, -21.643405928071164, 
#         -21.11683105640573, -20.590256184740294, -20.06368131307486, -19.53710644140942, -19.01053156974399]
#     assert list(red_M_err_list) == [0.2632874358327175, 0.2632874358327175, 0.2632874358327175, 0.2632874358327175, 0.2632874358327193, 
#         0.26328743583271574, 0.2632874358327193, 0.2632874358327175, 0.2632874358327175, 0.2632874358327175]
#     assert list(red_phi_list) == [5.262220153129362e-06, 1.146322895875979e-05, 2.2815766050216087e-05, 3.063244885429312e-05, 
#         3.784760366971176e-05, 6.955865011705814e-05, 4.601876296333914e-05, 4.2237548701635644e-05, 0.0001626682945797204, 5.198919360306014e-05]
#     assert list(red_phi_err_list) == [1.6816801547847714e-05, 1.8848825092512034e-05, 2.1415807046810993e-05, 3.71536660000855e-05, 
#         5.644501842774335e-05, 3.6815620551692435e-05, 5.656805584781544e-05, 7.501902486045965e-05, 0.00015384519165375355, 0.00021127915259708341]
#     assert list(blue_M_list) == [-24.62894308598833, -23.404512811278693, -22.18008253656906, -20.955652261859427, -19.731221987149794, 
#         -18.506791712440158, -17.282361437730522, -16.05793116302089, -14.833500888311256, -13.609070613601622]
#     assert list(blue_M_err_list) == [0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548164, 0.6122151373548181, 
#         0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548173, 0.6122151373548164]
#     assert list(blue_phi_list) == [290.491672773112, 265.7977762939234, 6.601873861989515e-05, 0.00019806224916273766, 0.0005366319861440035, 
#         0.0010535581878029952, 0.0019105283915257355, 0.005624556120539741, 0.0038603784240661696, 0.06417684969053933]
#     assert list(blue_phi_err_list) == [810.939765667249, 607.8170014696902, 3.096420475680061e-05, 0.0001368281769262127, 0.0005482677763962965, 
#         0.00044155205816259833, 0.0005126210576546734, 0.004590031416295528, 0.002209831427396876, 0.16436038459664015]
    
# #     test_result = cv2.imread('__tests__/test_LF_colour.png')
# #     pytest_result = cv2.imread('__tests__/pytest_LF_colour.png')
# #     difference = cv2.subtract(test_result, pytest_result)
# #     b, g, r = cv2.split(difference)
# #     assert test_result.shape == pytest_result.shape
# #     assert cv2.countNonZero(b) == 0
# #     assert cv2.countNonZero(g) == 0
# #     assert cv2.countNonZero(r) == 0

def test_filter_plot_by_colour_no_plot( ):
    r_maggy_ratios_table = pd.read_csv('__tests__/rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)

    zmax_table = pd.read_csv('__tests__/zmax_test.csv', delimiter=' ')
    z_max_list = np.array(zmax_table['zmax'])
    survey_area = 2.5 #sq. degrees
    Vmax_list = lf.get_volume(survey_area, z_max_list)
    
    n_bins = 10
    
    n_patches = 10
    ugriz_test_patch_centers_file_path = '__tests__/patch_centers_tol0.01_ugriz_test.csv'
    centers_table = np.genfromtxt(ugriz_test_patch_centers_file_path, delimiter=' ')
    ra_guesses = centers_table[ : , 0]
    dec_guesses = centers_table[ : , 1]
    ugriz_test_patch_centers_guesses = np.column_stack((ra_guesses, dec_guesses))

    g_maggy_ratios_table = pd.read_csv('__tests__/rest_maggy_ratios_g_ugriz_test.csv', delimiter=' ')
    g_maggy_ratio_list = np.array(g_maggy_ratios_table['maggy_ratio'])

    g_rest_mag_list = lf.get_rest_mag(z_photo_list, g_app_mag_list, g_maggy_ratio_list)

    colour_cut_slope = 0.0
    colour_cut_intercept = 0.65
    all_M_list, all_M_err_list, all_phi_list, all_phi_err_list, red_M_list, red_M_err_list, red_phi_list, red_phi_err_list, blue_M_list, blue_M_err_list, blue_phi_list, blue_phi_err_list = lf.filter_plot_by_colour(
        colour_cut_slope,
        colour_cut_intercept,
        r_rest_mag_list,
        g_rest_mag_list,
        Vmax_list,
        n_bins,
        RA_list,
        Dec_list,
        n_patches,
        ugriz_test_patch_centers_guesses,
        survey='kids',
        numba_installed=True)
    
    assert list(all_M_list) == [-24.62894308598833, -23.404512811278693, -22.18008253656906, -20.955652261859427, -19.731221987149794, 
        -18.506791712440158, -17.282361437730522, -16.05793116302089, -14.833500888311256, -13.609070613601622]
    assert list(all_M_err_list) == [0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548164, 0.6122151373548181,     
        0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548173, 0.6122151373548164]
    assert list(all_phi_list) == [290.491672773112, 265.7977860159212, 9.55747321361586e-05, 0.0002549444468315998, 0.0006247531885987247, 
        0.001075916507027022, 0.0019105283915257355, 0.005624556120539741, 0.0038603784240661696, 0.06417684969053933]
    assert list(all_phi_err_list) == [810.9397645708841, 607.8170004077618, 4.36417468536423e-05, 0.00019704012356457505, 
        0.0005486180652616492, 0.0004654318611417406, 0.0005773328573042754, 0.00459036071996224, 0.00221037276729238, 0.164362438153455]
    assert list(red_M_list) == [-23.749705414732908, -23.22313054306747, -22.696555671402038, -22.1699807997366, -21.643405928071164, 
        -21.11683105640573, -20.590256184740294, -20.06368131307486, -19.53710644140942, -19.01053156974399]
    assert list(red_M_err_list) == [0.2632874358327175, 0.2632874358327175, 0.2632874358327175, 0.2632874358327175, 0.2632874358327193, 
        0.26328743583271574, 0.2632874358327193, 0.2632874358327175, 0.2632874358327175, 0.2632874358327175]
    assert list(red_phi_list) == [5.262220153129362e-06, 1.146322895875979e-05, 2.2815766050216087e-05, 3.063244885429312e-05, 
        3.784760366971176e-05, 6.955865011705814e-05, 4.601876296333914e-05, 4.2237548701635644e-05, 0.0001626682945797204, 5.198919360306014e-05]
    assert list(red_phi_err_list) == [1.6816801547847714e-05, 1.8848825092512034e-05, 2.1415807046810993e-05, 3.71536660000855e-05, 
        5.644501842774335e-05, 3.6815620551692435e-05, 5.656805584781544e-05, 7.501902486045965e-05, 0.00015384519165375355, 0.00021127915259708341]
    assert list(blue_M_list) == [-24.62894308598833, -23.404512811278693, -22.18008253656906, -20.955652261859427, -19.731221987149794, 
        -18.506791712440158, -17.282361437730522, -16.05793116302089, -14.833500888311256, -13.609070613601622]
    assert list(blue_M_err_list) == [0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548164, 0.6122151373548181, 
        0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548173, 0.6122151373548164]
    assert list(blue_phi_list) == [290.491672773112, 265.7977762939234, 6.601873861989515e-05, 0.00019806224916273766, 0.0005366319861440035, 
        0.0010535581878029952, 0.0019105283915257355, 0.005624556120539741, 0.0038603784240661696, 0.06417684969053933]
    assert list(blue_phi_err_list) == [810.939765667249, 607.8170014696902, 3.096420475680061e-05, 0.0001368281769262127, 0.0005482677763962965, 
        0.00044155205816259833, 0.0005126210576546734, 0.004590031416295528, 0.002209831427396876, 0.16436038459664015]
# -----------------------