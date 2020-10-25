"""This file tests the following functions:

    * get_maggy
    * get_maggy_inv_var
    * get_obs_maggies_file
    * get_rec_maggies_files
    * get_rest_maggy_ratio_file
    * get_rest_mag
    
"""

# -----------------------
# Package Imports
# -----------------------
import math 
import numpy as np
import pandas as pd
import os
from pytest import approx
import sys
sys.path.insert(1, 'lumfunc/')
import lumfunc as lf

# -----------------------
# Unpack Test Data
# -----------------------
# test data (photometric galaxian survey)
data_table = pd.read_csv('test/catalogue_test.csv')

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
    assert list(r_maggies_list[0:4]) == approx([2.1712608416407457e-08, 1.8897275734216393e-08, 9.398644004513803e-09, 
        3.747264941992267e-08], rel=1e-6, abs=9e-4)

def test_get_maggy_rudimentary( ):
    assert list(lf.get_maggy(np.array([19.15822, 19.309002, 20.067337, 18.565714]))) == approx([2.1712608416407457e-08, 
        1.8897275734216393e-08, 9.398644004513803e-09, 3.747264941992267e-08], rel=1e-6, abs=9e-4)
# -----------------------

def test_get_maggy_inv_var( ):
    r_maggies_list = lf.get_maggy(r_app_mag_list)
    r_maggy_inv_var_list = lf.get_maggy_inv_var(r_maggies_list, r_app_mag_err_list)
    assert list(r_maggy_inv_var_list[0:4]) == approx([2.613536528041307e+20, 2.2153992456790612e+20, 2.6329570445628214e+20, 
        1.520308755766036e+20], rel=1e-6, abs=9e-4)

def test_get_maggy_inv_var_rudimentary( ):
    result = lf.get_maggy_inv_var(
        np.array([2.17126084e-08, 1.88972757e-08, 9.39864400e-09, 3.74726494e-08]), 
        np.array([0.00309313, 0.0038601, 0.0071193, 0.00234987]))
    assert list(result) == approx([2.613534842093864e+20, 2.2154049929330334e+20, 2.632956307424491e+20, 
        1.520310051334256e+20], rel=1e-6, abs=9e-4)
# -----------------------

def test_get_obs_maggies_file_ugriz( ):

    os.chdir('__tests__')

    ugriz_test_obs_maggies_file_path = 'obs_maggies_ugriz_pytest.csv'
    ugriz_test_bands = 'ugriz'
    lf.get_obs_maggies_file(ugriz_test_obs_maggies_file_path,
                            ugriz_test_bands,
                            z_photo_list,
                            u_app_mag_list,
                            g_app_mag_list,
                            r_app_mag_list,
                            i_app_mag_list,
                            Z_app_mag_list,
                            u_app_mag_err_list = u_app_mag_err_list,
                            g_app_mag_err_list = g_app_mag_err_list,
                            r_app_mag_err_list = r_app_mag_err_list,
                            i_app_mag_err_list = i_app_mag_err_list,
                            Z_app_mag_err_list = Z_app_mag_err_list)

    test_result_table = np.genfromtxt('../test/obs_maggies_ugriz_test.csv', delimiter=' ')
    test_result_z_photo_list = np.array(test_result_table[ : , 0])
    test_result_u_maggies_list = np.array(test_result_table[ : , 1])
    test_result_u_maggy_err_list = np.array(test_result_table[ : , 6])
    test_result_g_maggies_list = np.array(test_result_table[ : , 2])
    test_result_g_maggy_err_list = np.array(test_result_table[ : , 7])
    test_result_r_maggies_list = np.array(test_result_table[ : , 3])
    test_result_r_maggy_err_list = np.array(test_result_table[ : , 8])
    test_result_i_maggies_list = np.array(test_result_table[ : , 4])
    test_result_i_maggy_err_list = np.array(test_result_table[ : , 9])
    test_result_Z_maggies_list = np.array(test_result_table[ : , 5])
    test_result_Z_maggy_err_list = np.array(test_result_table[ : , 10])
    pytest_result_table = np.genfromtxt('obs_maggies_ugriz_pytest.csv', delimiter=' ')
    pytest_result_z_photo_list = np.array(pytest_result_table[ : , 0])
    pytest_result_u_maggies_list = np.array(pytest_result_table[ : , 1])
    pytest_result_u_maggy_err_list = np.array(pytest_result_table[ : , 6])
    pytest_result_g_maggies_list = np.array(pytest_result_table[ : , 2])
    pytest_result_g_maggy_err_list = np.array(pytest_result_table[ : , 7])
    pytest_result_r_maggies_list = np.array(pytest_result_table[ : , 3])
    pytest_result_r_maggy_err_list = np.array(pytest_result_table[ : , 8])
    pytest_result_i_maggies_list = np.array(pytest_result_table[ : , 4])
    pytest_result_i_maggy_err_list = np.array(pytest_result_table[ : , 9])
    pytest_result_Z_maggies_list = np.array(pytest_result_table[ : , 5])
    pytest_result_Z_maggy_err_list = np.array(pytest_result_table[ : , 10])

    assert list(pytest_result_z_photo_list) == list(test_result_z_photo_list)
    assert list(pytest_result_u_maggies_list) == list(test_result_u_maggies_list)
    assert list(pytest_result_u_maggy_err_list) == list(test_result_u_maggy_err_list)
    assert list(pytest_result_g_maggies_list) == list(test_result_g_maggies_list)
    assert list(pytest_result_g_maggy_err_list) == list(test_result_g_maggy_err_list)
    assert list(pytest_result_r_maggies_list) == list(test_result_r_maggies_list)
    assert list(pytest_result_r_maggy_err_list) == list(test_result_r_maggy_err_list)
    assert list(pytest_result_i_maggies_list) == list(test_result_i_maggies_list)
    assert list(pytest_result_i_maggy_err_list) == list(test_result_i_maggy_err_list)
    assert list(pytest_result_Z_maggies_list) == list(test_result_Z_maggies_list)
    assert list(pytest_result_Z_maggy_err_list) == list(test_result_Z_maggy_err_list)

def test_get_obs_maggies_file_ugriZYJHKs( ):

    os.chdir('../__tests__')

    ugriZYJHKs_test_obs_maggies_file_path = 'obs_maggies_ugriZYJHKs_pytest.csv'
    ugriZYJHKs_test_bands = 'ugriZYJHKs'
    lf.get_obs_maggies_file(ugriZYJHKs_test_obs_maggies_file_path,
                            ugriZYJHKs_test_bands,
                            z_photo_list,
                            u_app_mag_list,
                            g_app_mag_list,
                            r_app_mag_list,
                            i_app_mag_list,
                            Z_app_mag_list,
                            Y_app_mag_list = Y_app_mag_list,
                            J_app_mag_list = J_app_mag_list,
                            H_app_mag_list = H_app_mag_list,
                            Ks_app_mag_list = K_app_mag_list,
                            u_app_mag_err_list = u_app_mag_err_list,
                            g_app_mag_err_list = g_app_mag_err_list,
                            r_app_mag_err_list = r_app_mag_err_list,
                            i_app_mag_err_list = i_app_mag_err_list,
                            Z_app_mag_err_list = Z_app_mag_err_list,
                            Y_app_mag_err_list = Y_app_mag_err_list,
                            J_app_mag_err_list = J_app_mag_err_list,
                            H_app_mag_err_list = H_app_mag_err_list,
                            Ks_app_mag_err_list = K_app_mag_err_list)

    test_result_table = np.genfromtxt('../test/obs_maggies_ugriZYJHKs_test.csv', delimiter=' ')
    test_result_z_photo_list = np.array(test_result_table[ : , 0])
    test_result_u_maggies_list = np.array(test_result_table[ : , 1])
    test_result_u_maggy_err_list = np.array(test_result_table[ : , 10])
    test_result_g_maggies_list = np.array(test_result_table[ : , 2])
    test_result_g_maggy_err_list = np.array(test_result_table[ : , 11])
    test_result_r_maggies_list = np.array(test_result_table[ : , 3])
    test_result_r_maggy_err_list = np.array(test_result_table[ : , 12])
    test_result_i_maggies_list = np.array(test_result_table[ : , 4])
    test_result_i_maggy_err_list = np.array(test_result_table[ : , 13])
    test_result_Z_maggies_list = np.array(test_result_table[ : , 5])
    test_result_Z_maggy_err_list = np.array(test_result_table[ : , 14])
    test_result_Y_maggies_list = np.array(test_result_table[ : , 6])
    test_result_Y_maggy_err_list = np.array(test_result_table[ : , 15])
    test_result_J_maggies_list = np.array(test_result_table[ : , 7])
    test_result_J_maggy_err_list = np.array(test_result_table[ : , 16])
    test_result_H_maggies_list = np.array(test_result_table[ : , 8])
    test_result_H_maggy_err_list = np.array(test_result_table[ : , 17])
    test_result_Ks_maggies_list = np.array(test_result_table[ : , 9])
    test_result_Ks_maggy_err_list = np.array(test_result_table[ : , 18])
    pytest_result_table = np.genfromtxt('obs_maggies_ugriZYJHKs_pytest.csv', delimiter=' ')
    pytest_result_z_photo_list = np.array(pytest_result_table[ : , 0])
    pytest_result_u_maggies_list = np.array(pytest_result_table[ : , 1])
    pytest_result_u_maggy_err_list = np.array(pytest_result_table[ : , 10])
    pytest_result_g_maggies_list = np.array(pytest_result_table[ : , 2])
    pytest_result_g_maggy_err_list = np.array(pytest_result_table[ : , 11])
    pytest_result_r_maggies_list = np.array(pytest_result_table[ : , 3])
    pytest_result_r_maggy_err_list = np.array(pytest_result_table[ : , 12])
    pytest_result_i_maggies_list = np.array(pytest_result_table[ : , 4])
    pytest_result_i_maggy_err_list = np.array(pytest_result_table[ : , 13])
    pytest_result_Z_maggies_list = np.array(pytest_result_table[ : , 5])
    pytest_result_Z_maggy_err_list = np.array(pytest_result_table[ : , 14])
    pytest_result_Y_maggies_list = np.array(pytest_result_table[ : , 6])
    pytest_result_Y_maggy_err_list = np.array(pytest_result_table[ : , 15])
    pytest_result_J_maggies_list = np.array(pytest_result_table[ : , 7])
    pytest_result_J_maggy_err_list = np.array(pytest_result_table[ : , 16])
    pytest_result_H_maggies_list = np.array(pytest_result_table[ : , 8])
    pytest_result_H_maggy_err_list = np.array(pytest_result_table[ : , 17])
    pytest_result_Ks_maggies_list = np.array(pytest_result_table[ : , 9])
    pytest_result_Ks_maggy_err_list = np.array(pytest_result_table[ : , 18])

    assert list(pytest_result_z_photo_list) == list(test_result_z_photo_list)
    assert list(pytest_result_u_maggies_list) == list(test_result_u_maggies_list)
    assert list(pytest_result_u_maggy_err_list) == list(test_result_u_maggy_err_list)
    assert list(pytest_result_g_maggies_list) == list(test_result_g_maggies_list)
    assert list(pytest_result_g_maggy_err_list) == list(test_result_g_maggy_err_list)
    assert list(pytest_result_r_maggies_list) == list(test_result_r_maggies_list)
    assert list(pytest_result_r_maggy_err_list) == list(test_result_r_maggy_err_list)
    assert list(pytest_result_i_maggies_list) == list(test_result_i_maggies_list)
    assert list(pytest_result_i_maggy_err_list) == list(test_result_i_maggy_err_list)
    assert list(pytest_result_Z_maggies_list) == list(test_result_Z_maggies_list)
    assert list(pytest_result_Z_maggy_err_list) == list(test_result_Z_maggy_err_list)
    assert list(pytest_result_Y_maggies_list) == list(test_result_Y_maggies_list)
    assert list(pytest_result_Y_maggy_err_list) == list(test_result_Y_maggy_err_list)
    assert list(pytest_result_J_maggies_list) == list(test_result_J_maggies_list)
    assert list(pytest_result_J_maggy_err_list) == list(test_result_J_maggy_err_list)
    assert list(pytest_result_H_maggies_list) == list(test_result_H_maggies_list)
    assert list(pytest_result_H_maggy_err_list) == list(test_result_H_maggy_err_list)
    assert list(pytest_result_Ks_maggies_list) == list(test_result_Ks_maggies_list)
    assert list(pytest_result_Ks_maggy_err_list) == list(test_result_Ks_maggy_err_list)

# -----------------------

# def test_get_rec_maggies_files_5bands( ):

#     os.chdir('../__tests__')

#     ugriz_test_obs_maggies_file_path = 'obs_maggies_ugriz_pytest.csv'

#     z_values = np.arange(0.00, 1.00, 0.01)
#     rec_z_list = np.around(z_values, decimals=2)

#     ugriz_test_n_bands = 5
#     lf.get_rec_maggies_files(ugriz_test_obs_maggies_file_path,
#                              ugriz_test_n_bands,
#                              rec_z_list,
#                              rec_maggies_outfile_affix='ugriz_pytest',
#                              survey='sdss',
#                              band_z_shift=0.0,
#                              template_vmatrix_file_path='vmatrix.default.dat',
#                              template_lambda_file_path='lambda.default.dat',
#                              filters_list_file_path='sdss_filters.dat')
    
#     for rec_z in rec_z_list:
#         test_result = '../test/maggies_at_z' + str(rec_z) + '_ugriz_test.csv'
#         test_result_table = np.genfromtxt(test_result, delimiter=' ')
#         test_result_z_photo_list = np.array(test_result_table[ : , 0])
#         test_result_u_maggies_list = np.array(test_result_table[ : , 1])
#         test_result_g_maggies_list = np.array(test_result_table[ : , 2])
#         test_result_r_maggies_list = np.array(test_result_table[ : , 3])
#         test_result_i_maggies_list = np.array(test_result_table[ : , 4])
#         test_result_Z_maggies_list = np.array(test_result_table[ : , 5])
#         pytest_result = 'maggies_at_z' + str(rec_z) + '_ugriz_pytest.csv'
#         pytest_result_table = np.genfromtxt(pytest_result, delimiter=' ')
#         pytest_result_z_photo_list = np.array(pytest_result_table[ : , 0])
#         pytest_result_u_maggies_list = np.array(pytest_result_table[ : , 1])
#         pytest_result_g_maggies_list = np.array(pytest_result_table[ : , 2])
#         pytest_result_r_maggies_list = np.array(pytest_result_table[ : , 3])
#         pytest_result_i_maggies_list = np.array(pytest_result_table[ : , 4])
#         pytest_result_Z_maggies_list = np.array(pytest_result_table[ : , 5])
#         assert list(pytest_result_z_photo_list) == list(test_result_z_photo_list)
#         assert list(pytest_result_u_maggies_list) == list(test_result_u_maggies_list)
#         assert list(pytest_result_g_maggies_list) == list(test_result_g_maggies_list)
#         assert list(pytest_result_r_maggies_list) == list(test_result_r_maggies_list)
#         assert list(pytest_result_i_maggies_list) == list(test_result_i_maggies_list)
#         assert list(pytest_result_Z_maggies_list) == list(test_result_Z_maggies_list)

# def test_get_rec_maggies_files_9bands( ):

#     os.chdir('../__tests__')

#     ugriZYJHKs_test_obs_maggies_file_path = 'obs_maggies_ugriZYJHKs_pytest.csv'

#     z_values = np.arange(0.00, 1.00, 0.01)
#     rec_z_list = np.around(z_values, decimals=2)

#     # need test template files to run
#     ugriZYJHKs_test_n_bands = 9
#     lf.get_rec_maggies_files(ugriZYJHKs_test_obs_maggies_file_path,
#                              ugriZYJHKs_test_n_bands,
#                              rec_z_list,
#                              rec_maggies_outfile_affix='ugriZYJHKs_pytest',
#                              survey='test',
#                              band_z_shift=0.0,
#                              template_vmatrix_file_path='vmatrix.test.dat',
#                              template_lambda_file_path='lambda.test.dat',
#                              filters_list_file_path='test_filters.dat')

#     for rec_z in rec_z_list:
#         test_result = '../test/maggies_at_z' + str(rec_z) + '_ugriZYJHKs_test.csv'
#         test_result_table = np.genfromtxt(test_result, delimiter=' ')
#         test_result_z_photo_list = np.array(test_result_table[ : , 0])
#         test_result_u_maggies_list = np.array(test_result_table[ : , 1])
#         test_result_g_maggies_list = np.array(test_result_table[ : , 2])
#         test_result_r_maggies_list = np.array(test_result_table[ : , 3])
#         test_result_i_maggies_list = np.array(test_result_table[ : , 4])
#         test_result_Z_maggies_list = np.array(test_result_table[ : , 5])
#         test_result_Y_maggies_list = np.array(test_result_table[ : , 6])
#         test_result_J_maggies_list = np.array(test_result_table[ : , 7])
#         test_result_H_maggies_list = np.array(test_result_table[ : , 8])
#         test_result_Ks_maggies_list = np.array(test_result_table[ : , 9])
#         pytest_result = 'maggies_at_z' + str(rec_z) + '_ugriZYJHKs_pytest.csv'
#         pytest_result_table = np.genfromtxt(pytest_result, delimiter=' ')
#         pytest_result_z_photo_list = np.array(pytest_result_table[ : , 0])
#         pytest_result_u_maggies_list = np.array(pytest_result_table[ : , 1])
#         pytest_result_g_maggies_list = np.array(pytest_result_table[ : , 2])
#         pytest_result_r_maggies_list = np.array(pytest_result_table[ : , 3])
#         pytest_result_i_maggies_list = np.array(pytest_result_table[ : , 4])
#         pytest_result_Z_maggies_list = np.array(pytest_result_table[ : , 5])
#         pytest_result_Y_maggies_list = np.array(pytest_result_table[ : , 6])
#         pytest_result_J_maggies_list = np.array(pytest_result_table[ : , 7])
#         pytest_result_H_maggies_list = np.array(pytest_result_table[ : , 8])
#         pytest_result_Ks_maggies_list = np.array(pytest_result_table[ : , 9])
#         assert list(pytest_result_z_photo_list) == list(test_result_z_photo_list)
#         assert list(pytest_result_u_maggies_list) == list(test_result_u_maggies_list)
#         assert list(pytest_result_g_maggies_list) == list(test_result_g_maggies_list)
#         assert list(pytest_result_r_maggies_list) == list(test_result_r_maggies_list)
#         assert list(pytest_result_i_maggies_list) == list(test_result_i_maggies_list)
#         assert list(pytest_result_Z_maggies_list) == list(test_result_Z_maggies_list)
#         assert list(pytest_result_Y_maggies_list) == list(test_result_Y_maggies_list)
#         assert list(pytest_result_J_maggies_list) == list(test_result_J_maggies_list)
#         assert list(pytest_result_H_maggies_list) == list(test_result_H_maggies_list)
#         assert list(pytest_result_Ks_maggies_list) == list(test_result_Ks_maggies_list)
# -----------------------

def test_get_rest_maggy_ratio_file_r( ):

    os.chdir('../__tests__')

    ugriz_test_obs_maggies_file_path = 'obs_maggies_ugriz_pytest.csv'

    r_test_band_index = 3
    ugriz_test_rest_maggies_file_path = 'maggies_at_z0.0_ugriz_pytest.csv'
    lf.get_rest_maggy_ratio_file(ID_list,
                                ugriz_test_obs_maggies_file_path,
                                ugriz_test_rest_maggies_file_path,
                                r_test_band_index,
                                rest_maggy_ratio_outfile_affix='r_ugriz_pytest')
    
    test_result_table = pd.read_csv('../test/rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    test_result_ID_list = np.array(test_result_table['ID'])
    test_result_rec_z_list = np.array(test_result_table['rest_z'])
    test_result_maggy_ratio_list = np.array(test_result_table['maggy_ratio'])
    pytest_result_table = pd.read_csv('rest_maggy_ratios_r_ugriz_pytest.csv', delimiter=' ')
    pytest_result_ID_list = np.array(pytest_result_table['ID'])
    pytest_result_rec_z_list = np.array(pytest_result_table['rest_z'])
    pytest_result_maggy_ratio_list = np.array(pytest_result_table['maggy_ratio'])

    assert list(pytest_result_ID_list) == list(test_result_ID_list)
    assert list(pytest_result_rec_z_list) == list(test_result_rec_z_list)
    assert list(pytest_result_maggy_ratio_list) == list(test_result_maggy_ratio_list)
    
def test_get_rest_maggy_ratio_file_g( ):

    os.chdir('../__tests__')

    ugriz_test_obs_maggies_file_path = 'obs_maggies_ugriz_pytest.csv'

    g_test_band_index = 2
    ugriz_test_rest_maggies_file_path = 'maggies_at_z0.0_ugriz_pytest.csv'
    lf.get_rest_maggy_ratio_file(ID_list,
                                ugriz_test_obs_maggies_file_path,
                                ugriz_test_rest_maggies_file_path,
                                g_test_band_index,
                                rest_maggy_ratio_outfile_affix='g_ugriz_pytest')
    
    test_result_table = pd.read_csv('../test/rest_maggy_ratios_g_ugriz_test.csv', delimiter=' ')
    test_result_ID_list = np.array(test_result_table['ID'])
    test_result_rec_z_list = np.array(test_result_table['rest_z'])
    test_result_maggy_ratio_list = np.array(test_result_table['maggy_ratio'])
    pytest_result_table = pd.read_csv('rest_maggy_ratios_g_ugriz_pytest.csv', delimiter=' ')
    pytest_result_ID_list = np.array(pytest_result_table['ID'])
    pytest_result_rec_z_list = np.array(pytest_result_table['rest_z'])
    pytest_result_maggy_ratio_list = np.array(pytest_result_table['maggy_ratio'])

    assert list(pytest_result_ID_list) == list(test_result_ID_list)
    assert list(pytest_result_rec_z_list) == list(test_result_rec_z_list)
    assert list(pytest_result_maggy_ratio_list) == list(test_result_maggy_ratio_list)
# -----------------------

def test_get_rest_mag( ):
    os.chdir('..')
    r_maggy_ratios_table = pd.read_csv('test/rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)
    assert list(r_rest_mag_list[0:4]) == approx([-22.518710955165023, -20.36706085176511, -23.670847071363468, 
        -23.681182442586575], rel=1e-6, abs=9e-4)

def test_get_rest_mag_rudimentary( ):
    result = lf.get_rest_mag(np.array([0.34, 0.17, 0.61, 0.41]),
                np.array([19.15822, 19.309002, 20.067337, 18.565714]),
                np.array([0.69938735, 0.90226577, 0.43780755, 0.59193305]))
    assert list(result) == approx([-22.50048221280549, -20.367175595236922, -23.611903685248716, 
        -23.751335116325944], rel=1e-6, abs=9e-4)
# -----------------------