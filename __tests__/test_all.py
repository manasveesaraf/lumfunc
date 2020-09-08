"""This file tests the following functions:

    * get_maggy - yes
    * get_maggy_inv_var - yes
    * get_rest_mag - yes
    * get_volume - yes
    * get_binned_phi - yes
    * get_patch_centers - no - ??
    * get_patch_labels - not yet (check how to assert the image)
    * get_binned_phi_error - yes
    * get_plot - not yet (check how to assert the image)
    * filter_plot_by_colour - not yet (check how to assert the image)
    * SchechterMagModel - yes
    * DoubleSchechterMagModel - yes
    * get_gof - not yet
    * get_schechter_phi - not yet (check how to assert the image)
    * get_double_schechter_phi - not yet (check how to assert the image)
    
"""

# -----------------------
# Package Imports
# -----------------------
import numpy as np
import pandas as pd
# import pytest
# import sys
# sys.path.insert(1, '../lumfunc/')
import lumfunc as lf

# -----------------------
# Unpack Test Data
# -----------------------
# test data (photometric galaxian survey)
# data_table = pd.read_csv('__tests__/test_catalogue.csv')
data_table = pd.read_csv('test_catalogue.csv')

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

def test_get_maggy( ):    
    r_maggies_list = lf.get_maggy(r_app_mag_list)
    assert list(lf.get_maggy(r_app_mag_list)[0:4]) == [2.17126084e-08, 
        1.88972757e-08, 9.39864400e-09, 3.74726494e-08]

def test_get_maggy_rudimentary( ):
    assert list(lf.get_maggy(np.array([19.15822, 19.309002, 20.067337, 18.565714]))) == [2.17126084e-08, 
        1.88972757e-08, 9.39864400e-09, 3.74726494e-08]
# -----------------------

def test_get_maggy_inv_var( ):
    r_maggies_list = lf.get_maggy(r_app_mag_list)
    r_maggy_inv_var_list = lf.get_maggy_inv_var(r_maggies_list, r_app_mag_err_list)
    assert list(r_maggy_inv_var_list[0:4]) == [2.61353653e+20, 2.21539925e+20, 2.63295704e+20, 1.52030876e+20]

def test_get_maggy_inv_var_rudimentary( ):
    result = lf.get_maggy_inv_var(
        np.array([2.17126084e-08, 1.88972757e-08, 9.39864400e-09, 3.74726494e-08]), 
        np.array([0.00309313, 0.0038601, 0.0071193, 0.00234987]))
    assert list(result) == [2.61353484e+20, 2.21540499e+20, 2.63295631e+20, 1.52031005e+20]
# -----------------------

def test_get_rest_mag( ):
    r_maggy_ratios_table = pd.read_csv('rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)
    assert list(r_rest_mag_list[0:4]) == [-22.51871096, -20.36706085, -23.67084707, -23.68118244]

def test_get_rest_mag_rudimentary( ):
    result = lf.get_rest_mag(np.array([0.34, 0.17, 0.61, 0.41]),
                np.array([19.15822, 19.309002, 20.067337, 18.565714]),
                np.array([0.69938735, 0.90226577, 0.43780755, 0.59193305]))
    assert list(result) == [-22.50048221, -20.3671756,  -23.61190369, -23.75133512]
# -----------------------

def test_get_volume( ):
    zmax_table = pd.read_csv('test_zmax.csv', delimiter=' ')
    z_max_list = np.array(zmax_table['zmax'])

    survey_area = 2.5 #sq. degrees
    Vmax_list = lf.get_volume(survey_area, z_max_list)
    assert list(Vmax_list[:4]) == [1756716.17055371, 178625.22629838, 2447025.53293128, 2287569.94863823]

def test_get_volume_rudimentary( ):
    result = lf.get_volume(2.5, np.array([0.50523681, 0.21884399, 0.57489149, 0.55985663]))
    assert list(result) == [1756716.14012229, 178625.22858948, 2447025.55434235, 2287569.98290078]
# -----------------------

def test_get_binned_phi( ):
    r_maggy_ratios_table = pd.read_csv('rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)

    zmax_table = pd.read_csv('test_zmax.csv', delimiter=' ')
    z_max_list = np.array(zmax_table['zmax'])
    survey_area = 2.5 #sq. degrees
    Vmax_list = lf.get_volume(survey_area, z_max_list)

    n_bins = 10
    M_list, M_err_list, phi_list = lf.get_binned_phi(r_rest_mag_list, Vmax_list, n_bins)
    assert list(M_list) == [-24.62894309, -23.40451281, -22.18008254, -20.95565226, -19.73122199,
        -18.50679171, -17.28236144, -16.05793116, -14.83350089, -13.60907061]
    assert list(M_err_list) == [0.61221514, 0.61221514, 0.61221514, 0.61221514, 0.61221514, 0.61221514,
          0.61221514, 0.61221514, 0.61221514, 0.61221514]
    assert list(phi_list) == [2.90491673e+02, 2.65797786e+02, 9.55747321e-05, 2.54944447e-04,
          6.24753189e-04, 1.07591651e-03, 1.91052839e-03, 5.62455612e-03, 3.86037842e-03, 6.41768497e-02]

def test_get_binned_phi_rudimentary( ):
    M_result, M_err_result, phi_result = lf.get_binned_phi(
        np.array([-23, -21, -19, -22, -23, -23, -22, -23, -22, -22, -19, -21]),
        np.array([
            8e+08, 2e+08, 2e+07, 3e+08, 6e+08, 6e+08, 4e+08, 7e+08, 5e+08, 6e+08,
            7e+06, 1e+08
        ]), 4)
    assert list(M_result) == [-22.5, -21.5, -20.5, -19.5]
    assert list(M_err_result) == [0.5, 0.5, 0.5, 0.5]
    assert list(phi_result) == [1.06411667e-08, 1.02900000e-08, 0.00000000e+00, 1.32300000e-07]
# -----------------------

def test_get_binned_phi_error( ):
    r_maggy_ratios_table = pd.read_csv('rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)

    zmax_table = pd.read_csv('test_zmax.csv', delimiter=' ')
    z_max_list = np.array(zmax_table['zmax'])
    survey_area = 2.5 #sq. degrees
    Vmax_list = lf.get_volume(survey_area, z_max_list)
    
    n_bins = 10
    
    n_patches = 10
    ugriz_test_patch_centers_file_path = 'patch_centers_tol0.01_ugriz_test.csv'
    centers_table = np.genfromtxt(ugriz_test_patch_centers_file_path, delimiter=' ')
    ra_guesses = centers_table[ : , 0]
    dec_guesses = centers_table[ : , 1]
    ugriz_test_patch_centers_guesses = np.column_stack((ra_guesses, dec_guesses))
    
    labels = lf.get_patch_labels(RA_list,
                                Dec_list,
                                n_patches,
                                ugriz_test_patch_centers_guesses,
                                survey='kids',
                                numba_installed=True,
                                plot_savename='test_patches.png')

    phi_err_list = lf.get_binned_phi_error(r_rest_mag_list, Vmax_list, labels, n_patches, n_bins)
    assert list(phi_err_list) == [8.10939765e+02, 6.07817000e+02, 4.36417469e-05, 1.97040124e-04,
          5.48618065e-04, 4.65431861e-04, 5.77332857e-04, 4.59036072e-03, 2.21037277e-03, 1.64362438e-01]

def test_get_binned_phi_error_rudimentary( ):
    result = lf.get_binned_phi_error(
    np.array([-23, -21, -19, -22, -23, -23, -22, -23, -22, -22, -19, -21]),
    np.array([
        8e+08, 2e+08, 2e+07, 3e+08, 6e+08, 6e+08, 4e+08, 7e+08, 5e+08, 6e+08,
        7e+06, 1e+08
    ]), np.array([1, 1, 2, 2, 3, 0, 1, 1, 2, 2, 3, 3]), 4, 4)
    assert list(result) == [9.86494122e-09, 9.90155712e-09, 0.00000000e+00, 1.55859031e-07]
# -----------------------

def test_SchechterMagModel( ):
    r_maggy_ratios_table = pd.read_csv('rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)

    zmax_table = pd.read_csv('test_zmax.csv', delimiter=' ')
    z_max_list = np.array(zmax_table['zmax'])
    survey_area = 2.5 #sq. degrees
    Vmax_list = lf.get_volume(survey_area, z_max_list)

    n_bins = 10
    M_list, M_err_list, phi_list = lf.get_binned_phi(r_rest_mag_list, Vmax_list, n_bins)

    M_star_guess = -20.7
    phi_star_guess = 9.5e-3
    alpha_guess = -1.3
    sch1_model_phi_list = lf.SchechterMagModel(M_list, M_star_guess, phi_star_guess, alpha_guess)
    assert list(sch1_model_phi_list) == [1.88907752e-19, 2.36778419e-08, 1.16643327e-04, 2.29997398e-03,
          7.59124212e-03, 1.40466857e-02, 2.15508182e-02, 3.11177839e-02, 4.40579218e-02, 6.19837431e-02]

def test_SchechterMagModel_rudimentary( ):
    result = lf.SchechterMagModel(
        np.array([
            -25.1487769, -23.86987184, -22.59096677, -21.31206171, -20.03315665,
            -18.75425159, -17.47534652, -16.19644146, -14.9175364, -13.63863134
        ]), -20.7, 9.5e-3, -1.3)
    assert list(result) == [1.85685828e-29, 3.25671116e-11, 1.72458835e-05, 1.27468679e-03,
        6.12395219e-03, 1.26803535e-02, 2.02617665e-02, 2.98927403e-02, 4.30310959e-02, 6.14770529e-02]
# -----------------------

def test_DoubleSchechterMagModel( ):
    r_maggy_ratios_table = pd.read_csv('rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)

    zmax_table = pd.read_csv('test_zmax.csv', delimiter=' ')
    z_max_list = np.array(zmax_table['zmax'])
    survey_area = 2.5 #sq. degrees
    Vmax_list = lf.get_volume(survey_area, z_max_list)

    n_bins = 10
    M_list, M_err_list, phi_list = lf.get_binned_phi(r_rest_mag_list, Vmax_list, n_bins)

    M_star_guess = -20.7
    phi_star_1_guess = 6.16e-3
    alpha_1_guess = -0.79
    phi_star_2_guess = 6.16e-3
    alpha_2_guess = -0.79
    sch2_model_phi_list = lf.DoubleSchechterMagModel(M_list, 
                                                M_star_guess,
                                                phi_star_1_guess,
                                                alpha_1_guess,
                                                phi_star_2_guess,
                                                alpha_2_guess)
    assert list(sch2_model_phi_list) == [1.55110526e-18, 1.09383000e-07, 3.03168335e-04, 3.36328048e-03,
        6.24552903e-03, 6.50199270e-03, 5.61245148e-03, 4.55946326e-03, 3.63199542e-03, 2.87485077e-03]

def test_DoubleSchechterMagModel_rudimentary( ):
    result = lf.DoubleSchechterMagModel(
        np.array([
            -25.1487769, -23.86987184, -22.59096677, -21.31206171, -20.03315665,
            -18.75425159, -17.47534652, -16.19644146, -14.9175364, -13.63863134
        ]), -20.7, 6.16e-3, -0.79, 6.16e-3, -0.79)
    assert list(result) == [1.94632943e-28, 1.87206188e-10, 5.43662993e-05, 2.20369343e-03,
        5.80607779e-03, 6.59304119e-03, 5.77743541e-03, 4.67441094e-03, 3.69017477e-03, 2.89121865e-03]
# -----------------------