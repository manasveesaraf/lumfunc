# import pytest
import numpy as np
import pandas as pd

import sys
sys.path.insert(1, '../lumfunc/')
import lumfunc as lf


data_table = pd.read_csv('__tests__/test_catalogue.csv')
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


def test_get_maggy_a():    
    r_maggies_list = lf.get_maggy(r_app_mag_list)
    assert list(lf.get_maggy(r_app_mag_list)[0:4]) == [2.1712608416407457e-08, 
        1.8897275734216393e-08, 9.398644004513803e-09, 3.747264941992267e-08]

def test_get_maggy_b( ):
    assert list(lf.get_maggy(np.array([19.15822, 19.309002, 20.067337, 18.565714]))) == [2.1712608416407457e-08,
        1.8897275734216393e-08, 9.398644004513803e-09, 3.747264941992267e-08]

def test_get_maggy_inv_var_a( ):
    r_maggies_list = lf.get_maggy(r_app_mag_list)
    r_maggy_inv_var_list = lf.get_maggy_inv_var(r_maggies_list, r_app_mag_err_list)
    assert list(r_maggy_inv_var_list[0:4]) == [2.613536528041307e+20, 
        2.2153992456790612e+20, 2.6329570445628214e+20, 1.520308755766036e+20]

def test_get_maggy_inv_var_b( ):
    result = lf.get_maggy_inv_var(
        np.array([2.17126084e-08, 1.88972757e-08, 9.39864400e-09, 3.74726494e-08]),
        np.array([0.00309313, 0.0038601, 0.0071193, 0.00234987]))
    assert list(result) == [2.613534842093864e+20, 2.2154049929330334e+20, 2.632956307424491e+20, 1.520310051334256e+20]

def test_get_rest_mag( ):
    maggy_ratios_table = pd.read_csv('__tests__/test_maggy_ratios.csv', delimiter=' ')
    r_maggy_ratio_list = np.array(maggy_ratios_table['maggy_ratio'])
    r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)
    assert list(r_rest_mag_list[0:4]) == [-22.50048221746272, 
        -20.367175598846956, -23.611903680288837, -23.751335108989057]

