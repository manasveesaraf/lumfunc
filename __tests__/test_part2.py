"""This file tests the following functions:

    * get_patch_centers
    * get_patch_labels
    * get_binned_phi_error
    
"""

# -----------------------
# Package Imports
# -----------------------
import math 
import numpy as np
import pandas as pd
import cv2
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

def test_get_patch_centers( ):
    uniform_data_table = pd.read_csv('test_uniform_catalogue.csv')
    uniform_RA_list = np.array(uniform_data_table['uniform_RA'])
    uniform_Dec_list = np.array(uniform_data_table['uniform_Dec'])

    n_patches = 10
    centers_guesses = lf.get_patch_centers(uniform_RA_list,
                        uniform_Dec_list,
                        n_patches,
                        survey='kids',
                        max_iterations=int(100),
                        tolerance=1.0e-2)
    assert len(centers_guesses[:, 0]) == n_patches
    assert len(centers_guesses[:, 1]) == n_patches
# -----------------------

def test_get_patch_labels( ):
    ugriz_test_patch_centers_file_path = 'patch_centers_tol0.01_ugriz_test.csv'
    centers_table = np.genfromtxt(ugriz_test_patch_centers_file_path, delimiter=' ')
    ra_guesses = centers_table[ : , 0]
    dec_guesses = centers_table[ : , 1]
    ugriz_test_patch_centers_guesses = np.column_stack((ra_guesses, dec_guesses))

    n_patches = 10
    labels = lf.get_patch_labels(RA_list,
                                Dec_list,
                                n_patches,
                                ugriz_test_patch_centers_guesses,
                                survey='kids',
                                numba_installed=True,
                                plot_savename='pytest_patches.png')
    assert list(labels[0:4]) == [1, 3, 5, 6]
    
    test_result = cv2.imread('test_patches.png')
    pytest_result = cv2.imread('pytest_patches.png')
    difference = cv2.subtract(test_result, pytest_result)
    b, g, r = cv2.split(difference)
    assert test_result.shape == pytest_result.shape
    assert cv2.countNonZero(b) == 0
    assert cv2.countNonZero(g) == 0
    assert cv2.countNonZero(r) == 0
    
def test_get_patch_labels_no_numba( ):
    ugriz_test_patch_centers_file_path = 'patch_centers_tol0.01_ugriz_test.csv'
    centers_table = np.genfromtxt(ugriz_test_patch_centers_file_path, delimiter=' ')
    ra_guesses = centers_table[ : , 0]
    dec_guesses = centers_table[ : , 1]
    ugriz_test_patch_centers_guesses = np.column_stack((ra_guesses, dec_guesses))

    n_patches = 10
    labels = lf.get_patch_labels(RA_list,
                                Dec_list,
                                n_patches,
                                ugriz_test_patch_centers_guesses,
                                survey='kids',
                                numba_installed=False,
                                plot_savename='pytest_patches_no_numba.png')
    assert list(labels[0:4]) == [1, 3, 5, 6]
    
    test_result = cv2.imread('test_patches.png')
    pytest_result = cv2.imread('pytest_patches_no_numba.png')
    difference = cv2.subtract(test_result, pytest_result)
    b, g, r = cv2.split(difference)
    assert test_result.shape == pytest_result.shape
    assert cv2.countNonZero(b) == 0
    assert cv2.countNonZero(g) == 0
    assert cv2.countNonZero(r) == 0
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
                                survey='kids')

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