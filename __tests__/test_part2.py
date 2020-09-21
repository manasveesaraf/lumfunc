"""This file tests the following functions:

    * get_volume
    * get_binned_phi
    
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

def test_get_volume( ):
    zmax_table = pd.read_csv('__tests__/zmax_test.csv', delimiter=' ')
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
    r_maggy_ratios_table = pd.read_csv('__tests__/rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)

    zmax_table = pd.read_csv('__tests__/zmax_test.csv', delimiter=' ')
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