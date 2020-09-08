# -----------------------
# Package Imports
# -----------------------
import lumfunc as lf
import numpy as np
import pandas as pd

# -----------------------
# Unpack Test Data
# -----------------------
# test data (photometric galaxian survey)
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
r_maggies_list = lf.get_maggy(r_app_mag_list) 
print(r_maggies_list[0:4])
# returns 
# [2.17126084e-08 1.88972757e-08 9.39864400e-09 3.74726494e-08]

# rudimentarily:
lf.get_maggy(np.array([19.15822, 19.309002, 20.067337, 18.565714]))
# returns
# array([2.17126084e-08, 1.88972757e-08, 9.39864400e-09, 3.74726494e-08])
# -----------------------


r_maggy_inv_var_list = lf.get_maggy_inv_var(r_maggies_list, r_app_mag_err_list)
print(r_maggy_inv_var_list[0:4])
# returns 
# [2.61353653e+20 2.21539925e+20 2.63295704e+20 1.52030876e+20]

# rudimentarily:
lf.get_maggy_inv_var(
    np.array([2.17126084e-08, 1.88972757e-08, 9.39864400e-09, 3.74726494e-08]),
    np.array([0.00309313, 0.0038601, 0.0071193, 0.00234987]))
# returns
# array([2.61353484e+20, 2.21540499e+20, 2.63295631e+20, 1.52031005e+20])
# -----------------------


ugriz_test_obs_maggies_file_path = 'obs_maggies_ugriz_test.csv'
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

ugriZYJHKs_test_obs_maggies_file_path = 'obs_maggies_ugriZYJHKs_test.csv'
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
# -----------------------


z_values = np.arange(0.00, 1.00, 0.01)
rec_z_list = np.around(z_values, decimals=2)

ugriz_test_n_bands = 5
lf.get_rec_maggies_files(ugriz_test_obs_maggies_file_path,
                         ugriz_test_n_bands,
                         rec_z_list,
                         rec_maggies_outfile_affix='ugriz_test',
                         survey='sdss',
                         band_z_shift=0.0,
                         template_vmatrix_file_path='vmatrix.default.dat',
                         template_lambda_file_path='lambda.default.dat',
                         filters_list_file_path='sdss_filters.dat')

# # need test template files to run
# ugriZYJHKs_test_n_bands = 9
# lf.get_rec_maggies_files(ugriZYJHKs_test_obs_maggies_file_path,
#                          ugriZYJHKs_test_n_bands,
#                          rec_z_list,
#                          rec_maggies_outfile_affix='ugriZYJHKs_test',
#                          survey='test',
#                          band_z_shift=0.0,
#                          template_vmatrix_file_path='vmatrix.test.dat',
#                          template_lambda_file_path='lambda.test.dat',
#                          filters_list_file_path='test_filters.dat')
# -----------------------


r_test_band_index = 3
ugriz_test_rest_maggies_file_path = 'maggies_at_z0.0_ugriz_test.csv'
lf.get_rest_maggy_ratio_file(ID_list,
                             ugriz_test_obs_maggies_file_path,
                             ugriz_test_rest_maggies_file_path,
                             r_test_band_index,
                             rest_maggy_ratio_outfile_affix='r_ugriz_test')

g_test_band_index = 2
ugriz_test_rest_maggies_file_path = 'maggies_at_z0.0_ugriz_test.csv'
get_rest_maggy_ratio_file(ID_list,
                             ugriz_test_obs_maggies_file_path,
                             ugriz_test_rest_maggies_file_path,
                             g_test_band_index,
                             rest_maggy_ratio_outfile_affix='g_ugriz_test')
# -----------------------


r_maggy_ratios_table = pd.read_csv('rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])

r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)
print(r_rest_mag_list[0:4])
# returns 
# [-22.51871096 -20.36706085 -23.67084707 -23.68118244]

# rudimentarily:
lf.get_rest_mag(np.array([0.34, 0.17, 0.61, 0.41]),
                np.array([19.15822, 19.309002, 20.067337, 18.565714]),
                np.array([0.69938735, 0.90226577, 0.43780755, 0.59193305]))
# returns
# array([-22.50048221, -20.3671756 , -23.61190369, -23.75133512])
# -----------------------


ugriz_test_rec_maggies_file_path = 'maggies_at_z0.01_ugriz_test.csv'
lf.get_maggy_ratio_file(ID_list,
                        ugriz_test_rest_maggies_file_path,
                        ugriz_test_rec_maggies_file_path,
                        0.01,
                        r_test_band_index,
                        maggy_ratio_outfile_affix='r_ugriz_test')
# -----------------------


lf.get_all_maggy_ratios_file(rec_z_list,
                             ID_list,
                             r_test_band_index,
                             maggies_and_out_files_affix='ugriz_test')
# -----------------------


zmax_table = pd.read_csv('test_zmax.csv', delimiter=' ')
z_max_list = np.array(zmax_table['zmax'])

survey_area = 2.5 #sq. degrees
Vmax_list = lf.get_volume(survey_area, z_max_list)
print(Vmax_list[:4])
# returns 
# [1756716.17055371  178625.22629838 2447025.53293128 2287569.94863823]

# rudimentarily:
lf.get_volume(2.5, np.array([0.50523681, 0.21884399, 0.57489149, 0.55985663]))
# returns
# array([1756716.14012229, 178625.22858948, 2447025.55434235, 2287569.98290078])
# -----------------------


n_bins = 10
M_list, M_err_list, phi_list = lf.get_binned_phi(r_rest_mag_list, Vmax_list, n_bins)
print(M_list)
# returns
# [-24.62894309 -23.40451281 -22.18008254 -20.95565226 -19.73122199
#  -18.50679171 -17.28236144 -16.05793116 -14.83350089 -13.60907061]
print(M_err_list)
# returns
# [0.61221514 0.61221514 0.61221514 0.61221514 0.61221514 
#  0.61221514 0.61221514 0.61221514 0.61221514 0.61221514]
print(phi_list)
# returns 
# [2.90491673e+02 2.65797786e+02 9.55747321e-05 2.54944447e-04 6.24753189e-04
#  1.07591651e-03 1.91052839e-03 5.62455612e-03 3.86037842e-03 6.41768497e-02]

# OR a rudimentarily example:
lf.get_binned_phi(
    np.array([-23, -21, -19, -22, -23, -23, -22, -23, -22, -22, -19, -21]),
    np.array([
        8e+08, 2e+08, 2e+07, 3e+08, 6e+08, 6e+08, 4e+08, 7e+08, 5e+08, 6e+08,
        7e+06, 1e+08
    ]), 4)
# returns 
# (array([-22.5, -21.5, -20.5, -19.5]),
#  array([0.5, 0.5, 0.5, 0.5]),
#  array([1.06411667e-08, 1.02900000e-08, 0.00000000e+00, 1.32300000e-07]))
# -----------------------


uniform_data_table = pd.read_csv('test_uniform_catalogue.csv')
uniform_RA_list = np.array(uniform_data_table['uniform_RA'])
uniform_Dec_list = np.array(uniform_data_table['uniform_Dec'])

n_patches = 10
lf.get_patch_centers(uniform_RA_list,
                     uniform_Dec_list,
                     n_patches,
                     survey='kids',
                     max_iterations=int(100),
                     tolerance=1.0e-2,
                     patch_centers_outfile_affix='ugriz_test')
# saves file patch_centers_tol0.01_ugriz_test.csv
# -----------------------


ugriz_test_patch_centers_file_path = 'patch_centers_tol0.01_ugriz_test.csv'
labels = lf.get_patch_labels(RA_list,
                             Dec_list,
                             n_patches,
                             ugriz_test_patch_centers_file_path,
                             survey='kids',
                             numba_installed=True,
                             plot_savename='test_patches.png')
# displays plot
# -----------------------


phi_err_list = lf.get_binned_phi_error(r_rest_mag_list, Vmax_list, labels, n_patches, n_bins)
print(phi_err_list)
# returns
# [8.10939765e+02 6.07817000e+02 4.36417469e-05 1.97040124e-04 5.48618065e-04
#  4.65431861e-04 5.77332857e-04 4.59036072e-03 2.21037277e-03 1.64362438e-01]

lf.get_binned_phi_error(
    np.array([-23, -21, -19, -22, -23, -23, -22, -23, -22, -22, -19, -21]),
    np.array([
        8e+08, 2e+08, 2e+07, 3e+08, 6e+08, 6e+08, 4e+08, 7e+08, 5e+08, 6e+08,
        7e+06, 1e+08
    ]), np.array([1, 1, 2, 2, 3, 0, 1, 1, 2, 2, 3, 3]), 4, 4)
# -----------------------


M_list, M_err_list, phi_list, phi_err_list = lf.get_plot(
    r_rest_mag_list,
    Vmax_list,
    n_bins,
    RA_list,
    Dec_list,
    n_patches,
    ugriz_test_patch_centers_file_path,
    survey='kids',
    numba_installed=True,
    plot_savename='test_LF.png')
# displays plot
# -----------------------


g_maggy_ratios_table = pd.read_csv('rest_maggy_ratios_g_ugriz_test.csv', delimiter=' ')
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
    ugriz_test_patch_centers_file_path,
    survey='kids',
    numba_installed=True,
    plot_savename='test_LF_colour.png')
# displays plot
# -----------------------


M_star_guess = -20.7
phi_star_guess = 9.5e-3
alpha_guess = -1.3
sch1_model_phi_list = lf.SchechterMagModel(M_list, M_star_guess, phi_star_guess, alpha_guess)
print(sch1_model_phi_list)
# returns
# [1.88907752e-19 2.36778419e-08 1.16643327e-04 2.29997398e-03 7.59124212e-03 
#  1.40466857e-02 2.15508182e-02 3.11177839e-02 4.40579218e-02 6.19837431e-02]

# OR a rudimentarily example:
lf.SchechterMagModel(
    np.array([
        -25.1487769, -23.86987184, -22.59096677, -21.31206171, -20.03315665,
        -18.75425159, -17.47534652, -16.19644146, -14.9175364, -13.63863134
    ]), -20.7, 9.5e-3, -1.3)
# -----------------------


M_star_guess = -20.7
phi_star_1_guess = 6.16e-3
alpha_1_guess = -0.79
phi_star_2_guess = 6.16e-3
alpha_2_guess = -0.79
sch2_model_phi_list = DoubleSchechterMagModel(M_list, 
                                              M_star_guess,
                                              phi_star_1_guess,
                                              alpha_1_guess,
                                              phi_star_2_guess,
                                              alpha_2_guess)
print(sch2_model_phi_list)
# returns
# [1.55110526e-18 1.09383000e-07 3.03168335e-04 3.36328048e-03 6.24552903e-03
#  6.50199270e-03 5.61245148e-03 4.55946326e-03 3.63199542e-03 2.87485077e-03]

# OR a rudimentarily example:
lf.DoubleSchechterMagModel(
    np.array([
        -25.1487769, -23.86987184, -22.59096677, -21.31206171, -20.03315665,
        -18.75425159, -17.47534652, -16.19644146, -14.9175364, -13.63863134
    ]), -20.7, 6.16e-3, -0.79, 6.16e-3, -0.79)
# -----------------------


m = 3
gof = lf.get_gof(phi_list, phi_err_list, sch1_model_phi_list, m)
print(gof)
# returns
# 366.43103358282144

# OR a rudimentarily example:
lf.get_gof(
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
# -----------------------


all_sch1_model_phi_list, all_chi_sq_1, all_M_star, all_M_star_err, all_phi_star, all_phi_star_err, all_alpha_star, all_alpha_star_err = lf.get_schechter_phi(
    all_M_list,
    all_M_err_list,
    all_phi_list,
    all_phi_err_list,
    np.array([M_star_guess, phi_star_guess, alpha_guess]),
    plot_savename='test_all_Sch.png')
# displays plot

blue_sch1_model_phi_list, blue_chi_sq_1, blue_M_star, blue_M_star_err, blue_phi_star, blue_phi_star_err, blue_alpha_star, blue_alpha_star_err = lf.get_schechter_phi(
    blue_M_list,
    blue_M_err_list,
    blue_phi_list,
    blue_phi_err_list,
    np.array([M_star_guess, phi_star_guess, alpha_guess]),
    plot_savename='test_blue_Sch.png')
# displays plot

# OR a rudimentarily example:
lf.get_schechter_phi(np.array([
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
    np.array([-20.71, 9.5e-3, -1.3]),
    plot_savename='rud_test_Sch.png')
# -----------------------


red_sch2_model_phi_list, red_chi_sq_1, red_M_star, red_M_star_err, red_phi_star_1, red_phi_star_err_1, red_phi_star_2, red_phi_star_err_2, red_alpha_star_1, red_alpha_star_err_1, red_alpha_star_2, red_alpha_star_err_2 = lf.get_double_schechter_phi(
    red_M_list,
    red_M_err_list,
    red_phi_list,
    red_phi_err_list,
    np.array([M_star_guess, phi_star_1_guess, alpha_1_guess, phi_star_2_guess, alpha_2_guess]),
    plot_savename='test_red_dSch.png')
# displays plot

# OR a rudimentarily example:
lf.get_double_schechter_phi(np.array([
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
    np.array([-20.7, 6.16e-3, -0.79, 6.16e-3, -0.79]),
    plot_savename='rud_test_dSch.png')
# -----------------------