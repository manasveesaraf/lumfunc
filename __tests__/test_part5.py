"""This file tests the following functions:

    * SchechterMagModel
    * DoubleSchechterMagModel
    * get_gof
    * get_schechter_phi
    * get_double_schechter_phi
    
"""

# -----------------------
# Package Imports
# -----------------------
import math 
import numpy as np
import pandas as pd
# import cv2
from pytest import approx
import sys
sys.path.insert(1, 'lumfunc/')
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

def test_SchechterMagModel( ):
    r_maggy_ratios_table = pd.read_csv('__tests__/rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)

    zmax_table = pd.read_csv('__tests__/zmax_test.csv', delimiter=' ')
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
    r_maggy_ratios_table = pd.read_csv('__tests__/rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)

    zmax_table = pd.read_csv('__tests__/zmax_test.csv', delimiter=' ')
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

def test_get_gof( ):
    r_maggy_ratios_table = pd.read_csv('__tests__/rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)

    zmax_table = pd.read_csv('__tests__/zmax_test.csv', delimiter=' ')
    z_max_list = np.array(zmax_table['zmax'])
    survey_area = 2.5 #sq. degrees
    Vmax_list = lf.get_volume(survey_area, z_max_list)

    n_bins = 10
    M_list, M_err_list, phi_list = lf.get_binned_phi(r_rest_mag_list, Vmax_list, n_bins)

    M_star_guess = -20.7
    phi_star_guess = 9.5e-3
    alpha_guess = -1.3
    sch1_model_phi_list = lf.SchechterMagModel(M_list, M_star_guess, phi_star_guess, alpha_guess)

    n_patches = 10
    ugriz_test_patch_centers_file_path = '__tests__/patch_centers_tol0.01_ugriz_test.csv'
    centers_table = np.genfromtxt(ugriz_test_patch_centers_file_path, delimiter=' ')
    ra_guesses = centers_table[ : , 0]
    dec_guesses = centers_table[ : , 1]
    ugriz_test_patch_centers_guesses = np.column_stack((ra_guesses, dec_guesses))
    
    labels = lf.get_patch_labels(RA_list,
                                Dec_list,
                                n_patches,
                                ugriz_test_patch_centers_guesses,
                                survey='kids')

    phi_err_list = lf.get_binned_phi_error(r_rest_mag_list, Vmax_list, labels, n_patches, n_bins)

    m = 3
    gof = lf.get_gof(phi_list, phi_err_list, sch1_model_phi_list, m)
    assert gof == 366.43103358282144

def test_get_gof_rudimentary( ):
    gof_result = lf.get_gof(
        np.array([
            2.78118218e+02, 2.54476157e+02, 6.57347457e-05, 1.98257155e-04,
            4.84943102e-04, 1.02149157e-03, 1.49165665e-03, 4.54012724e-03,
            5.08195775e-03, 6.14432455e-02
        ]),
        np.array([
            6.31512459e+02, 5.32152268e+02, 4.31666309e-05, 2.22841109e-04,
            4.81148550e-04, 3.16386417e-04, 6.52443936e-04, 4.68698737e-03,
            2.05929233e-03, 1.60744165e-01
        ]),
        np.array([
            1.94632963e-28, 1.87206201e-10, 5.43662983e-05, 2.20369342e-03,
            5.80607779e-03, 6.59304119e-03, 5.77743541e-03, 4.67441094e-03,
            3.69017477e-03, 2.89121864e-03
        ]), 3)
    assert gof_result == 79.66254082924551
# -----------------------

# def test_get_schechter_phi_all( ):
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
#         survey='kids')
    
#     M_star_guess = -20.7
#     phi_star_guess = 9.5e-3
#     alpha_guess = -1.3

#     all_sch1_model_phi_list, all_chi_sq, all_M_star, all_M_star_err, all_phi_star, all_phi_star_err, all_alpha_star, all_alpha_star_err = lf.get_schechter_phi(
#         all_M_list,
#         all_M_err_list,
#         all_phi_list,
#         all_phi_err_list,
#         np.array([M_star_guess, phi_star_guess, alpha_guess]),
#         plot_savename='pytest_all_Sch.png')
    
    # assert list(all_sch1_model_phi_list) == [2.8325898598722236e-09, 5.743891439174078e-06, 9.259399753521627e-05, 0.00031201159158978584, 
    #     0.000633385107084906, 0.0010912081614801337, 0.0017826941154962066, 0.0028627054171209425, 0.0045714941754318224, 0.007287133378452643]
    # assert round(all_chi_sq,5) == round(0.14910742282850892,5)
    # assert round(all_M_star,5) == round(-22.068531742285295,5)
    # assert round(all_M_star_err,5) == round(0.3557347014819093,5)
    # assert round(all_phi_star,5) == round(0.0003176940137059405,5)
    # assert round(all_phi_star_err,5) == round(0.0001288373384458377,5)
    # assert round(all_alpha_star,5) == round(-1.4126892538229192,5)
    # assert round(all_alpha_star_err,5) == round(0.06081125190828317,5)

#     # test_result = cv2.imread('__tests__/test_all_Sch.png')
#     # pytest_result = cv2.imread('__tests__/pytest_all_Sch.png')
#     # difference = cv2.subtract(test_result, pytest_result)
#     # b, g, r = cv2.split(difference)
#     # assert test_result.shape == pytest_result.shape
#     # assert cv2.countNonZero(b) == 0
#     # assert cv2.countNonZero(g) == 0
#     # assert cv2.countNonZero(r) == 0

def test_get_schechter_phi_all_no_plot( ):
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
        survey='kids')
    
    M_star_guess = -20.7
    phi_star_guess = 9.5e-3
    alpha_guess = -1.3

    all_sch1_model_phi_list, all_chi_sq, all_M_star, all_M_star_err, all_phi_star, all_phi_star_err, all_alpha_star, all_alpha_star_err = lf.get_schechter_phi(
        all_M_list,
        all_M_err_list,
        all_phi_list,
        all_phi_err_list,
        np.array([M_star_guess, phi_star_guess, alpha_guess]))
    
    assert list(all_sch1_model_phi_list) == [2.8325898598722236e-09, 5.743891439174078e-06, 9.259399753521627e-05, 0.00031201159158978584, 
        0.000633385107084906, 0.0010912081614801337, 0.0017826941154962066, 0.0028627054171209425, 0.0045714941754318224, 0.007287133378452643]
    assert round(all_chi_sq,5) == round(0.14910742282850892,5)
    assert round(all_M_star,5) == round(-22.068531742285295,5)
    assert round(all_M_star_err,5) == round(0.3557347014819093,5)
    assert round(all_phi_star,5) == round(0.0003176940137059405,5)
    assert round(all_phi_star_err,5) == round(0.0001288373384458377,5)
    assert round(all_alpha_star,5) == round(-1.4126892538229192,5)
    assert round(all_alpha_star_err,5) == round(0.06081125190828317,5)

# def test_get_schechter_phi_blue( ):
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
#         survey='kids')

#     M_star_guess = -20.7
#     phi_star_guess = 9.5e-3
#     alpha_guess = -1.3

#     blue_sch1_model_phi_list, blue_chi_sq, blue_M_star, blue_M_star_err, blue_phi_star, blue_phi_star_err, blue_alpha_star, blue_alpha_star_err = lf.get_schechter_phi(
#         blue_M_list,
#         blue_M_err_list,
#         blue_phi_list,
#         blue_phi_err_list,
#         np.array([M_star_guess, phi_star_guess, alpha_guess]),
#         plot_savename='pytest_blue_Sch.png')
    
    # assert list(blue_sch1_model_phi_list) == approx([1.9842034834819953e-10, 2.1809348172275195e-06, 6.210185218825129e-05, 0.00025711149229255473, 
    #     0.0005701699443248206, 0.0010330083316408773, 0.0017530015118458727, 0.0029124539166523883, 0.004805700054962961, 0.007912062557208274])
    # assert blue_chi_sq == approx(0.18163420708695324)
    # assert blue_M_star == approx(-21.842075975175316)
    # assert blue_M_star_err == approx(0.31845816378631797)
    # assert blue_phi_star == approx(0.0003029586014597913)
    # assert blue_phi_star_err == approx(0.00012126827264875354)
    # assert blue_alpha_star == approx(-1.4411669183679228)
    # assert blue_alpha_star_err == approx(0.06358938020533868)

#     # test_result = cv2.imread('__tests__/test_blue_Sch.png')
#     # pytest_result = cv2.imread('__tests__/pytest_blue_Sch.png')
#     # difference = cv2.subtract(test_result, pytest_result)
#     # b, g, r = cv2.split(difference)
#     # assert test_result.shape == pytest_result.shape
#     # assert cv2.countNonZero(b) == 0
#     # assert cv2.countNonZero(g) == 0
#     # assert cv2.countNonZero(r) == 0

def test_get_schechter_phi_blue_no_plot( ):
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
        survey='kids')

    M_star_guess = -20.7
    phi_star_guess = 9.5e-3
    alpha_guess = -1.3

    blue_sch1_model_phi_list, blue_chi_sq, blue_M_star, blue_M_star_err, blue_phi_star, blue_phi_star_err, blue_alpha_star, blue_alpha_star_err = lf.get_schechter_phi(
        blue_M_list,
        blue_M_err_list,
        blue_phi_list,
        blue_phi_err_list,
        np.array([M_star_guess, phi_star_guess, alpha_guess]))
    
    assert list(blue_sch1_model_phi_list) == approx([1.9842034834819953e-10, 2.1809348172275195e-06, 6.210185218825129e-05, 0.00025711149229255473, 
        0.0005701699443248206, 0.0010330083316408773, 0.0017530015118458727, 0.0029124539166523883, 0.004805700054962961, 0.007912062557208274])
    assert blue_chi_sq == approx(0.18163420708695324)
    assert blue_M_star == approx(-21.842075975175316)
    assert blue_M_star_err == approx(0.31845816378631797)
    assert blue_phi_star == approx(0.0003029586014597913)
    assert blue_phi_star_err == approx(0.00012126827264875354)
    assert blue_alpha_star == approx(-1.4411669183679228)
    assert blue_alpha_star_err == approx(0.06358938020533868)

# def test_get_schechter_phi_rudimentarily( ):
#     all_sch1_model_phi_result, chi_sq_sch1_result, M_star_result, M_star_err_result, phi_star_result, phi_star_err_result, alpha_star_result, alpha_star_err_result = lf.get_schechter_phi(np.array([
#         -24.7, -24.1, -23.5, -22.9, -22.3, -21.7, -21.1, -20.5, -19.9, -19.3,
#         -18.7, -18.1, -17.5, -16.9
#     ]),
#         np.ones(14) * 0.3,
#         np.array([
#             8.1e-07, 3.9e-06, 3.7e-05, 1.9e-04, 4.1e-04, 6.7e-04,
#             9.1e-04, 1.1e-03, 1.5e-03, 2.5e-03, 3.5e-03, 3.7e-03,
#             5.1e-03, 7.6e-03
#         ]),
#         np.array([
#             2.6e-07, 1.3e-06, 9.5e-06, 3.9e-05, 8.5e-05, 1.4e-04,
#             2.1e-04, 2.6e-04, 3.9e-04, 6.9e-04, 1.1e-03, 1.2e-03,
#             1.4e-03, 2.2e-03
#         ]),
#         np.array([-20.71, 9.5e-3, -1.3]),
#         plot_savename='pytest_rud_Sch.png')
#     # displays plot
    # assert list(all_sch1_model_phi_result) == [1.7105353558263635e-07, 5.1640448032129975e-06, 4.053707109679924e-05, 0.0001465834405778961, 
    #     0.00033931824553655666, 0.0006076603569602913, 0.0009387957310397152, 0.001332203535517021, 0.0018002389419955371, 0.002365197612525618, 
    #     0.003057535570029159, 0.003915871425980178, 0.004988344284986229, 0.006334965713301948]
    # assert round(chi_sq_sch1_result,5) == round(1.0209802688993401,5)
    # assert round(M_star_result,5) == round(-22.51627500778435,5)
    # assert round(M_star_err_result,5) == round(0.0964342301982251,5)
    # assert round(phi_star_result,5) == round(0.0007681235644217974,5)
    # assert round(phi_star_err_result,5) == round(0.00015735301981608952,5)
    # assert round(alpha_star_result,5) == round(-1.4248810024852225,5)
    # assert round(alpha_star_err_result,5) == round(0.06007607488402875,5)

#     # test_result = cv2.imread('__tests__/test_rud_Sch.png')
#     # pytest_result = cv2.imread('__tests__/pytest_rud_Sch.png')
#     # difference = cv2.subtract(test_result, pytest_result)
#     # b, g, r = cv2.split(difference)
#     # assert test_result.shape == pytest_result.shape
#     # assert cv2.countNonZero(b) == 0
#     # assert cv2.countNonZero(g) == 0
#     # assert cv2.countNonZero(r) == 0

def test_get_schechter_phi_rudimentarily_no_plot( ):
    all_sch1_model_phi_result, chi_sq_sch1_result, M_star_result, M_star_err_result, phi_star_result, phi_star_err_result, alpha_star_result, alpha_star_err_result = lf.get_schechter_phi(np.array([
        -24.7, -24.1, -23.5, -22.9, -22.3, -21.7, -21.1, -20.5, -19.9, -19.3,
        -18.7, -18.1, -17.5, -16.9
    ]),
        np.ones(14) * 0.3,
        np.array([
            8.1e-07, 3.9e-06, 3.7e-05, 1.9e-04, 4.1e-04, 6.7e-04,
            9.1e-04, 1.1e-03, 1.5e-03, 2.5e-03, 3.5e-03, 3.7e-03,
            5.1e-03, 7.6e-03
        ]),
        np.array([
            2.6e-07, 1.3e-06, 9.5e-06, 3.9e-05, 8.5e-05, 1.4e-04,
            2.1e-04, 2.6e-04, 3.9e-04, 6.9e-04, 1.1e-03, 1.2e-03,
            1.4e-03, 2.2e-03
        ]),
        np.array([-20.71, 9.5e-3, -1.3]))
    # displays plot
    assert list(all_sch1_model_phi_result) == [1.7105353558263635e-07, 5.1640448032129975e-06, 4.053707109679924e-05, 0.0001465834405778961, 
        0.00033931824553655666, 0.0006076603569602913, 0.0009387957310397152, 0.001332203535517021, 0.0018002389419955371, 0.002365197612525618, 
        0.003057535570029159, 0.003915871425980178, 0.004988344284986229, 0.006334965713301948]
    assert round(chi_sq_sch1_result,5) == round(1.0209802688993401,5)
    assert round(M_star_result,5) == round(-22.51627500778435,5)
    assert round(M_star_err_result,5) == round(0.0964342301982251,5)
    assert round(phi_star_result,5) == round(0.0007681235644217974,5)
    assert round(phi_star_err_result,5) == round(0.00015735301981608952,5)
    assert round(alpha_star_result,5) == round(-1.4248810024852225,5)
    assert round(alpha_star_err_result,5) == round(0.06007607488402875,5)
# -----------------------

# def test_get_double_schechter_phi_red( ):
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
#         survey='kids')

#     M_star_guess = -20.7
#     phi_star_1_guess = 6.16e-3
#     alpha_1_guess = -0.79
#     phi_star_2_guess = 6.16e-3
#     alpha_2_guess = -0.79

#     red_sch2_model_phi_list, red_chi_sq, red_M_star, red_M_star_err, red_phi_star_1, red_phi_star_err_1, red_phi_star_2, red_phi_star_err_2, red_alpha_star_1, red_alpha_star_err_1, red_alpha_star_2, red_alpha_star_err_2 = lf.get_double_schechter_phi(
#         red_M_list,
#         red_M_err_list,
#         red_phi_list,
#         red_phi_err_list,
#         np.array([M_star_guess, phi_star_1_guess, alpha_1_guess, phi_star_2_guess, alpha_2_guess]),
#         plot_savename='pytest_red_dSch.png')
#     # displays plot
    # assert list(red_sch2_model_phi_list) == approx([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.6375293296239083e-229, 2.3525320265182555e-141, 
    #     2.752009553823618e-87])
    # assert round(red_chi_sq,2) == round(1.2084645603920292,2)
    # assert round(red_M_star,2) == round(-13.256375067114597,2)
    # assert red_M_star_err == math.inf
    # assert round(red_phi_star_1,2) == round(-0.005143924152379018,2)
    # assert red_phi_star_err_1 == math.inf
    # assert round(red_phi_star_2,2) == round(-1.872910729853815,2)
    # assert red_phi_star_err_2 == math.inf
    # assert round(red_alpha_star_1,2) == round(0.012183946742584995,2)
    # assert red_alpha_star_err_1 == math.inf
    # assert round(red_alpha_star_2,2) == round(0.025603076393042268,2)
    # assert red_alpha_star_err_2 == math.inf

#     # test_result = cv2.imread('__tests__/test_red_dSch.png')
#     # pytest_result = cv2.imread('__tests__/pytest_red_dSch.png')
#     # difference = cv2.subtract(test_result, pytest_result)
#     # b, g, r = cv2.split(difference)
#     # assert test_result.shape == pytest_result.shape
#     # assert cv2.countNonZero(b) == 0
#     # assert cv2.countNonZero(g) == 0
#     # assert cv2.countNonZero(r) == 0

def test_get_double_schechter_phi_red_no_plot( ):
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
        survey='kids')

    M_star_guess = -20.7
    phi_star_1_guess = 6.16e-3
    alpha_1_guess = -0.79
    phi_star_2_guess = 6.16e-3
    alpha_2_guess = -0.79

    red_sch2_model_phi_list, red_chi_sq, red_M_star, red_M_star_err, red_phi_star_1, red_phi_star_err_1, red_phi_star_2, red_phi_star_err_2, red_alpha_star_1, red_alpha_star_err_1, red_alpha_star_2, red_alpha_star_err_2 = lf.get_double_schechter_phi(
        red_M_list,
        red_M_err_list,
        red_phi_list,
        red_phi_err_list,
        np.array([M_star_guess, phi_star_1_guess, alpha_1_guess, phi_star_2_guess, alpha_2_guess]))
    # displays plot
    assert list(red_sch2_model_phi_list) == approx([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.6375293296239083e-229, 2.3525320265182555e-141, 
        2.752009553823618e-87])
    assert round(red_chi_sq,2) == round(1.2084645603920292,2)
    assert round(red_M_star,2) == round(-13.256375067114597,2)
    assert red_M_star_err == math.inf
    assert round(red_phi_star_1,2) == round(-0.005143924152379018,2)
    assert red_phi_star_err_1 == math.inf
    assert round(red_phi_star_2,2) == round(-1.872910729853815,2)
    assert red_phi_star_err_2 == math.inf
    assert round(red_alpha_star_1,2) == round(0.012183946742584995,2)
    assert red_alpha_star_err_1 == math.inf
    assert round(red_alpha_star_2,2) == round(0.025603076393042268,2)
    assert red_alpha_star_err_2 == math.inf


# def test_get_double_schechter_phi_rudimentarily( ):
#     all_sch2_model_phi_result, chi_sq_sch2_result, M_star_result, M_star_err_result, phi_star_1_result, phi_star_err_1_result, phi_star_2_result, phi_star_err_2_result, alpha_star_1_result, alpha_star_err_1_result, alpha_star_2_result, alpha_star_err_2_result = lf.get_double_schechter_phi(np.array([
#         -24.7, -24.1, -23.5, -22.9, -22.3, -21.7, -21.1, -20.5, -19.9, -19.3,
#         -18.7, -18.1, -17.5, -16.9
#     ]),
#         np.ones(14) * 0.3,
#         np.array([
#             8.05e-07, 3.88e-06, 3.69e-05, 1.89e-04,
#             4.05e-04, 6.72e-04, 9.09e-04, 1.11e-03,
#             1.48e-03, 2.49e-03, 3.51e-03, 3.72e-03,
#             5.01e-03, 7.55e-03
#         ]),
#         np.array([
#             2.61e-07, 1.25e-06, 9.52e-06, 3.89e-05,
#             8.49e-05, 1.39e-04, 2.00e-04, 2.57e-04,
#             3.95e-04, 6.88e-04, 1.10e-03, 1.17e-03,
#             1.39e-03, 2.17e-03
#         ]),
#         np.array([-20.7, 6.16e-3, -0.79, 6.16e-3, -0.79]),
#         plot_savename='pytest_rud_dSch.png')
#     # displays plot
    # # round all assertions to 2 decimal places
    # all_result = list(all_sch2_model_phi_result)
    # all_result = [ round(x,2) for x in all_result ]
    # all_test = [8.521602535554413e-08, 4.304795096021841e-06, 4.252947712862992e-05, 
    #     0.00016513644802319284, 0.00037724853104172785, 0.0006409589905341704, 0.0009291455434703172, 0.001246599413378984, 
    #     0.0016250833276945204, 0.0021183671618024385, 0.002805526837713822, 0.003802654108449027, 0.0052833317077602675, 
    #     0.007510562710100609]
    # all_test = [ round(x,2) for x in all_test ]
    # assert all_real == all_test
    # assert round(chi_sq_sch2_result,2) == round(0.8888283543610924,2)
    # assert round(M_star_result,2) == round(-22.303908380116704,2)
    # assert round(M_star_err_result,2) == round(0.26464127945271887,2)
    # assert round(phi_star_1_result,2) == round(0.0009668887609189701,2)
    # assert round(phi_star_err_1_result,2) == round(0.000640187578339006,2)
    # assert round(phi_star_2_result,2) == round(-1.0900241221219484,2)
    # assert round(phi_star_err_2_result,2) == round(0.7987986322969173,2)
    # assert round(alpha_star_1_result,2) == round(0.0001418318772494868,2)
    # assert round(alpha_star_err_1_result,2) == round(0.0008399596540331241,2)
    # assert round(alpha_star_2_result,2) == round(-1.774506451062984,2)
    # assert round(alpha_star_err_2_result,2) == round(0.9946532141625982,2)

#     # test_result = cv2.imread('__tests__/test_rud_dSch.png')
#     # pytest_result = cv2.imread('__tests__/pytest_rud_dSch.png')
#     # difference = cv2.subtract(test_result, pytest_result)
#     # b, g, r = cv2.split(difference)
#     # assert test_result.shape == pytest_result.shape
#     # assert cv2.countNonZero(b) == 0
#     # assert cv2.countNonZero(g) == 0
#     # assert cv2.countNonZero(r) == 0

def test_get_double_schechter_phi_rudimentarily_no_plot( ):
    all_sch2_model_phi_result, chi_sq_sch2_result, M_star_result, M_star_err_result, phi_star_1_result, phi_star_err_1_result, phi_star_2_result, phi_star_err_2_result, alpha_star_1_result, alpha_star_err_1_result, alpha_star_2_result, alpha_star_err_2_result = lf.get_double_schechter_phi(np.array([
        -24.7, -24.1, -23.5, -22.9, -22.3, -21.7, -21.1, -20.5, -19.9, -19.3,
        -18.7, -18.1, -17.5, -16.9
    ]),
        np.ones(14) * 0.3,
        np.array([
            8.05e-07, 3.88e-06, 3.69e-05, 1.89e-04,
            4.05e-04, 6.72e-04, 9.09e-04, 1.11e-03,
            1.48e-03, 2.49e-03, 3.51e-03, 3.72e-03,
            5.01e-03, 7.55e-03
        ]),
        np.array([
            2.61e-07, 1.25e-06, 9.52e-06, 3.89e-05,
            8.49e-05, 1.39e-04, 2.00e-04, 2.57e-04,
            3.95e-04, 6.88e-04, 1.10e-03, 1.17e-03,
            1.39e-03, 2.17e-03
        ]),
        np.array([-20.7, 6.16e-3, -0.79, 6.16e-3, -0.79]))
    
    # round all assertions to 2 decimal places
    phi_result = list(all_sch2_model_phi_result)
    phi_result = [ round(x,2) for x in phi_result ]
    phi_test = [8.521602535554413e-08, 4.304795096021841e-06, 4.252947712862992e-05, 
        0.00016513644802319284, 0.00037724853104172785, 0.0006409589905341704, 0.0009291455434703172, 0.001246599413378984, 
        0.0016250833276945204, 0.0021183671618024385, 0.002805526837713822, 0.003802654108449027, 0.0052833317077602675, 
        0.007510562710100609]
    phi_test = [ round(x,2) for x in phi_test ]
    assert phi_result == phi_test
    assert round(chi_sq_sch2_result,2) == round(0.8888283543610924,2)
    assert round(M_star_result,2) == round(-22.303908380116704,2)
    assert round(M_star_err_result,2) == round(0.26464127945271887,2)
    assert round(phi_star_1_result,2) == round(0.0009668887609189701,2)
    assert round(phi_star_err_1_result,2) == round(0.000640187578339006,2)
    assert round(phi_star_2_result,2) == round(-1.0900241221219484,2)
    assert round(phi_star_err_2_result,2) == round(0.7987986322969173,2)
    assert round(alpha_star_1_result,2) == round(0.0001418318772494868,2)
    assert round(alpha_star_err_1_result,2) == round(0.0008399596540331241,2)
    assert round(alpha_star_2_result,2) == round(-1.774506451062984,2)
    assert round(alpha_star_err_2_result,3) == round(0.9946532141625982,3)
# -----------------------