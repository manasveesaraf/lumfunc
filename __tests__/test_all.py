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
    assert list(lf.get_maggy(r_app_mag_list)[0:4]) == [2.1712608416407457e-08, 
        1.8897275734216393e-08, 9.398644004513803e-09, 3.747264941992267e-08]

def test_get_maggy_rudimentary( ):
    assert list(lf.get_maggy(np.array([19.15822, 19.309002, 20.067337, 18.565714]))) == [2.1712608416407457e-08, 
        1.8897275734216393e-08, 9.398644004513803e-09, 3.747264941992267e-08]
# -----------------------

def test_get_maggy_inv_var( ):
    r_maggies_list = lf.get_maggy(r_app_mag_list)
    r_maggy_inv_var_list = lf.get_maggy_inv_var(r_maggies_list, r_app_mag_err_list)
    assert list(r_maggy_inv_var_list[0:4]) == [2.613536528041307e+20, 
        2.2153992456790612e+20, 2.6329570445628214e+20, 1.520308755766036e+20]

def test_get_maggy_inv_var_rudimentary( ):
    result = lf.get_maggy_inv_var(
        np.array([2.17126084e-08, 1.88972757e-08, 9.39864400e-09, 3.74726494e-08]), 
        np.array([0.00309313, 0.0038601, 0.0071193, 0.00234987]))
    assert list(result) == [2.613534842093864e+20, 
        2.2154049929330334e+20, 2.632956307424491e+20, 1.520310051334256e+20]
# -----------------------

def test_get_rest_mag( ):
    r_maggy_ratios_table = pd.read_csv('rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)
    assert list(r_rest_mag_list[0:4]) == [-22.518710955165023, 
        -20.36706085176511, -23.670847071363468, -23.681182442586575]

def test_get_rest_mag_rudimentary( ):
    result = lf.get_rest_mag(np.array([0.34, 0.17, 0.61, 0.41]),
                np.array([19.15822, 19.309002, 20.067337, 18.565714]),
                np.array([0.69938735, 0.90226577, 0.43780755, 0.59193305]))
    assert list(result) == [-22.50048221280549, 
        -20.367175595236922, -23.611903685248716, -23.751335116325944]
# -----------------------

def test_get_volume( ):
    zmax_table = pd.read_csv('test_zmax.csv', delimiter=' ')
    z_max_list = np.array(zmax_table['zmax'])

    survey_area = 2.5 #sq. degrees
    Vmax_list = lf.get_volume(survey_area, z_max_list)
    assert list(Vmax_list[:4]) == [1756716.170553711, 
        178625.22629838384, 2447025.5329312785, 2287569.9486382306]

def test_get_volume_rudimentary( ):
    result = lf.get_volume(2.5, np.array([0.50523681, 0.21884399, 0.57489149, 0.55985663]))
    assert list(result) == [1756716.140122295, 
        178625.2285894822, 2447025.5543423546, 2287569.9829007764]
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
    assert list(M_list) == [-24.62894308598833, -23.404512811278693, -22.18008253656906, -20.955652261859427, 
        -19.731221987149794, -18.506791712440158, -17.282361437730522, -16.05793116302089, -14.833500888311256, 
        -13.609070613601622]
    assert list(M_err_list) == [0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548164, 
    0.6122151373548181, 0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548173, 
    0.6122151373548164]
    assert list(phi_list) == [290.491672773112, 265.7977860159212, 9.55747321361586e-05, 0.0002549444468315998, 
    0.0006247531885987247, 0.001075916507027022, 0.0019105283915257355, 0.005624556120539741, 0.0038603784240661696, 
    0.06417684969053933]

def test_get_binned_phi_rudimentary( ):
    M_result, M_err_result, phi_result = lf.get_binned_phi(
        np.array([-23, -21, -19, -22, -23, -23, -22, -23, -22, -22, -19, -21]),
        np.array([
            8e+08, 2e+08, 2e+07, 3e+08, 6e+08, 6e+08, 4e+08, 7e+08, 5e+08, 6e+08,
            7e+06, 1e+08
        ]), 4)
    assert list(M_result) == [-22.5, -21.5, -20.5, -19.5]
    assert list(M_err_result) == [0.5, 0.5, 0.5, 0.5]
    assert list(phi_result) == [1.0641166666666664e-08, 1.0289999999999998e-08, 0.0, 1.3229999999999996e-07]
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
    assert list(phi_err_list) == [810.9397645708841, 607.8170004077618, 4.36417468536423e-05, 
        0.00019704012356457505, 0.0005486180652616492, 0.0004654318611417406, 0.0005773328573042754, 
        0.00459036071996224, 0.00221037276729238, 0.164362438153455]

def test_get_binned_phi_error_rudimentary( ):
    result = lf.get_binned_phi_error(
    np.array([-23, -21, -19, -22, -23, -23, -22, -23, -22, -22, -19, -21]),
    np.array([
        8e+08, 2e+08, 2e+07, 3e+08, 6e+08, 6e+08, 4e+08, 7e+08, 5e+08, 6e+08,
        7e+06, 1e+08
    ]), np.array([1, 1, 2, 2, 3, 0, 1, 1, 2, 2, 3, 3]), 4, 4)
    assert list(result) == [9.864941217372873e-09, 9.901557116602078e-09, 0.0, 1.5585903075108177e-07]
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
    assert list(sch1_model_phi_list) == [1.889077519830292e-19, 2.3677841853441523e-08, 0.00011664332685035941, 
        0.002299973977228317, 0.007591242124153264, 0.014046685725374422, 0.021550818203935722, 0.03111778390239495, 
        0.044057921833823524, 0.061983743062742416]

def test_SchechterMagModel_rudimentary( ):
    result = lf.SchechterMagModel(
        np.array([
            -25.1487769, -23.86987184, -22.59096677, -21.31206171, -20.03315665,
            -18.75425159, -17.47534652, -16.19644146, -14.9175364, -13.63863134
        ]), -20.7, 9.5e-3, -1.3)
    assert list(result) == [1.8568582821959972e-29, 3.256711157578086e-11, 1.7245883492815743e-05, 
        0.0012746867894433624, 0.006123952187528593, 0.012680353540368563, 0.0202617665019816, 
        0.02989274030134303, 0.04303109591644922, 0.061477052920961686]
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
    assert list(sch2_model_phi_list) == [1.5511052649970585e-18, 1.0938300029533067e-07, 0.0003031683345133834, 
        0.0033632804762819004, 0.006245529031705679, 0.006501992701830366, 0.005612451475589847, 
        0.0045594632577110345, 0.00363199541697254, 0.0028748507739080474]

def test_DoubleSchechterMagModel_rudimentary( ):
    result = lf.DoubleSchechterMagModel(
        np.array([
            -25.1487769, -23.86987184, -22.59096677, -21.31206171, -20.03315665,
            -18.75425159, -17.47534652, -16.19644146, -14.9175364, -13.63863134
        ]), -20.7, 6.16e-3, -0.79, 6.16e-3, -0.79)
    assert list(result) == [1.9463294285254365e-28, 1.8720618798131103e-10, 5.436629926965087e-05, 
        0.0022036934269993955, 0.005806077792944941, 0.006593041192333871, 0.00577743540908892, 
        0.004674410937670327, 0.0036901747679265527, 0.002891218646017407]
# -----------------------