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
data_table = pd.read_csv('catalogue_test.csv')

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
print(list(r_maggies_list[0:4]))
# returns 
# [2.17126084e-08 1.88972757e-08 9.39864400e-09 3.74726494e-08]

# rudimentarily:
r_maggies_result = lf.get_maggy(np.array([19.15822, 19.309002, 20.067337, 18.565714]))
print(list(r_maggies_result[0:4]))
# returns
# [2.17126084e-08 1.88972757e-08 9.39864400e-09 3.74726494e-08]
# -----------------------


r_maggy_inv_var_list = lf.get_maggy_inv_var(r_maggies_list, r_app_mag_err_list)
print(list(r_maggy_inv_var_list[0:4]))
# returns 
# [2.61353653e+20 2.21539925e+20 2.63295704e+20 1.52030876e+20]

# rudimentarily:
r_maggy_inv_var_result = lf.get_maggy_inv_var(
    np.array([2.17126084e-08, 1.88972757e-08, 9.39864400e-09, 3.74726494e-08]),
    np.array([0.00309313, 0.0038601, 0.0071193, 0.00234987]))
print(list(r_maggy_inv_var_result[0:4]))
# returns
# [2.61353484e+20 2.21540499e+20 2.63295631e+20 1.52031005e+20]
# -----------------------


r_maggy_ratios_table = pd.read_csv('rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])

r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)
print(list(r_rest_mag_list[0:4]))
# returns 
# [-22.51871096 -20.36706085 -23.67084707 -23.68118244]

# rudimentarily:
r_rest_mag_result = lf.get_rest_mag(np.array([0.34, 0.17, 0.61, 0.41]),
                                  np.array([19.15822, 19.309002, 20.067337, 18.565714]),
                                  np.array([0.69938735, 0.90226577, 0.43780755, 0.59193305]))
print(list(r_rest_mag_result[0:4]))
# returns
# [-22.50048221 -20.3671756  -23.61190369 -23.75133512]
# -----------------------


zmax_table = pd.read_csv('zmax_test.csv', delimiter=' ')
z_max_list = np.array(zmax_table['zmax'])

survey_area = 2.5 #sq. degrees
Vmax_list = lf.get_volume(survey_area, z_max_list)
print(list(Vmax_list[:4]))
# returns 
# [1756716.17055371  178625.22629838 2447025.53293128 2287569.94863823]

# rudimentarily:
Vmax_result = lf.get_volume(2.5, np.array([0.50523681, 0.21884399, 0.57489149, 0.55985663]))
print(list(Vmax_result[:4]))
# returns
# [1756716.14012229  178625.22858948 2447025.55434235 2287569.98290078]
# -----------------------


n_bins = 10
M_list, M_err_list, phi_list = lf.get_binned_phi(r_rest_mag_list, Vmax_list, n_bins)
print(list(M_list))
# returns
# [-24.62894309 -23.40451281 -22.18008254 -20.95565226 -19.73122199
#  -18.50679171 -17.28236144 -16.05793116 -14.83350089 -13.60907061]
print(list(M_err_list))
# returns
# [0.61221514 0.61221514 0.61221514 0.61221514 0.61221514 0.61221514
#  0.61221514 0.61221514 0.61221514 0.61221514]
print(list(phi_list))
# returns 
# [2.90491673e+02 2.65797786e+02 9.55747321e-05 2.54944447e-04
#  6.24753189e-04 1.07591651e-03 1.91052839e-03 5.62455612e-03
#  3.86037842e-03 6.41768497e-02]

# OR a rudimentarily example:
M_result, M_err_result, phi_result = lf.get_binned_phi(
    np.array([-23, -21, -19, -22, -23, -23, -22, -23, -22, -22, -19, -21]),
    np.array([
        8e+08, 2e+08, 2e+07, 3e+08, 6e+08, 6e+08, 4e+08, 7e+08, 5e+08, 6e+08,
        7e+06, 1e+08
    ]), 4)
print(list(M_result))
# returns 
# [-22.5 -21.5 -20.5 -19.5]
print(list(M_err_result))
# returns 
# [0.5 0.5 0.5 0.5]
print(list(phi_result))
# returns
# [1.06411667e-08 1.02900000e-08 0.00000000e+00 1.32300000e-07]
# -----------------------


uniform_data_table = pd.read_csv('uniform_catalogue_test.csv')
uniform_RA_list = np.array(uniform_data_table['uniform_RA'])
uniform_Dec_list = np.array(uniform_data_table['uniform_Dec'])

n_patches = 10
centers_guesses = lf.get_patch_centers(uniform_RA_list,
                     uniform_Dec_list,
                     n_patches,
                     survey='kids',
                     max_iterations=int(100),
                     tolerance=1.0e-2)
print(list(centers_guesses))
# returns
# [[ 2.23334587e+02  2.05711043e+00]
#  [ 2.23351178e+02  1.21001503e+00]
#  [ 2.23242242e+02 -6.92140179e-01]
#  [ 2.23163230e+02  1.42282290e+00]
#  [ 2.23238201e+02  7.33355022e-01]
#  [ 2.23291130e+02 -1.81978780e+00]
#  [ 2.23252168e+02  6.51768885e-02]
#  [ 2.23211561e+02  1.78274686e+00]
#  [ 2.23212695e+02  2.38909725e+00]
#  [ 2.23249748e+02  2.79812147e+00]]
# -----------------------


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
# displays plot
print(list(labels[0:4]))
# returns
# [1 3 5 6]
# -----------------------


phi_err_list = lf.get_binned_phi_error(r_rest_mag_list, Vmax_list, labels, n_patches, n_bins)
print(list(phi_err_list))
# returns
# [8.10939765e+02 6.07817000e+02 4.36417469e-05 1.97040124e-04
#  5.48618065e-04 4.65431861e-04 5.77332857e-04 4.59036072e-03
#  2.21037277e-03 1.64362438e-01]

phi_err_result = lf.get_binned_phi_error(
    np.array([-23, -21, -19, -22, -23, -23, -22, -23, -22, -22, -19, -21]),
    np.array([
        8e+08, 2e+08, 2e+07, 3e+08, 6e+08, 6e+08, 4e+08, 7e+08, 5e+08, 6e+08,
        7e+06, 1e+08
    ]), np.array([1, 1, 2, 2, 3, 0, 1, 1, 2, 2, 3, 3]), 4, 4)
print(list(phi_err_result))
# returns
# [9.86494122e-09 9.90155712e-09 0.00000000e+00 1.55859031e-07]
# -----------------------


plot_M_list, plot_M_err_list, plot_phi_list, plot_phi_err_list = lf.get_plot(
    r_rest_mag_list,
    Vmax_list,
    n_bins,
    RA_list,
    Dec_list,
    n_patches,
    ugriz_test_patch_centers_guesses,
    survey='kids',
    numba_installed=True,
    plot_savename='test_LF.png')
# displays plot
print(list(plot_M_list))
# returns
# [-24.62894309 -23.40451281 -22.18008254 -20.95565226 -19.73122199
#  -18.50679171 -17.28236144 -16.05793116 -14.83350089 -13.60907061]
print(list(plot_M_err_list))
# returns
# [0.61221514 0.61221514 0.61221514 0.61221514 0.61221514 0.61221514
#  0.61221514 0.61221514 0.61221514 0.61221514]
print(list(plot_phi_list))
# returns 
# [2.90491673e+02 2.65797786e+02 9.55747321e-05 2.54944447e-04
#  6.24753189e-04 1.07591651e-03 1.91052839e-03 5.62455612e-03
#  3.86037842e-03 6.41768497e-02]
print(list(plot_phi_err_list))
# returns
# [8.10939765e+02 6.07817000e+02 4.36417469e-05 1.97040124e-04
#  5.48618065e-04 4.65431861e-04 5.77332857e-04 4.59036072e-03
#  2.21037277e-03 1.64362438e-01]
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
    ugriz_test_patch_centers_guesses,
    survey='kids',
    numba_installed=True,
    plot_savename='test_LF_colour.png')
# displays plot
print(list(all_M_list))
# returns
# [-24.62894309 -23.40451281 -22.18008254 -20.95565226 -19.73122199
#  -18.50679171 -17.28236144 -16.05793116 -14.83350089 -13.60907061]
print(list(all_M_err_list))
# returns
# [0.61221514 0.61221514 0.61221514 0.61221514 0.61221514 0.61221514
#  0.61221514 0.61221514 0.61221514 0.61221514]
print(list(all_phi_list))
# returns 
# [2.90491673e+02 2.65797786e+02 9.55747321e-05 2.54944447e-04
#  6.24753189e-04 1.07591651e-03 1.91052839e-03 5.62455612e-03
#  3.86037842e-03 6.41768497e-02]
print(list(all_phi_err_list))
# returns
# [8.10939765e+02 6.07817000e+02 4.36417469e-05 1.97040124e-04
#  5.48618065e-04 4.65431861e-04 5.77332857e-04 4.59036072e-03
#  2.21037277e-03 1.64362438e-01]
print(list(red_M_list))
# returns
# [-23.74970541 -23.22313054 -22.69655567 -22.1699808  -21.64340593
#  -21.11683106 -20.59025618 -20.06368131 -19.53710644 -19.01053157]
print(list(red_M_err_list))
# returns
# [0.26328744 0.26328744 0.26328744 0.26328744 0.26328744 0.26328744
#  0.26328744 0.26328744 0.26328744 0.26328744]
print(list(red_phi_list))
# returns 
# [5.26222015e-06 1.14632290e-05 2.28157661e-05 3.06324489e-05
#  3.78476037e-05 6.95586501e-05 4.60187630e-05 4.22375487e-05
#  1.62668295e-04 5.19891936e-05]
print(list(red_phi_err_list))
# returns
# [1.68168015e-05 1.88488251e-05 2.14158070e-05 3.71536660e-05
#  5.64450184e-05 3.68156206e-05 5.65680558e-05 7.50190249e-05
#  1.53845192e-04 2.11279153e-04]
print(list(blue_M_list))
# returns
# [-24.62894309 -23.40451281 -22.18008254 -20.95565226 -19.73122199
#  -18.50679171 -17.28236144 -16.05793116 -14.83350089 -13.60907061]
print(list(blue_M_err_list))
# returns
# [0.61221514 0.61221514 0.61221514 0.61221514 0.61221514 0.61221514
#  0.61221514 0.61221514 0.61221514 0.61221514]
print(list(blue_phi_list))
# returns 
# [2.90491673e+02 2.65797776e+02 6.60187386e-05 1.98062249e-04
#  5.36631986e-04 1.05355819e-03 1.91052839e-03 5.62455612e-03
#  3.86037842e-03 6.41768497e-02]
print(list(blue_phi_err_list))
# returns
# [8.10939766e+02 6.07817001e+02 3.09642048e-05 1.36828177e-04
#  5.48267776e-04 4.41552058e-04 5.12621058e-04 4.59003142e-03
#  2.20983143e-03 1.64360385e-01]
# -----------------------


M_star_guess = -20.7
phi_star_guess = 9.5e-3
alpha_guess = -1.3
sch1_model_phi_list = lf.SchechterMagModel(M_list, M_star_guess, phi_star_guess, alpha_guess)
print(list(sch1_model_phi_list))
# returns
# [1.88907752e-19 2.36778419e-08 1.16643327e-04 2.29997398e-03
#  7.59124212e-03 1.40466857e-02 2.15508182e-02 3.11177839e-02
#  4.40579218e-02 6.19837431e-02]

# OR a rudimentarily example:
sch1_model_phi_result = lf.SchechterMagModel(
    np.array([
        -25.1487769, -23.86987184, -22.59096677, -21.31206171, -20.03315665,
        -18.75425159, -17.47534652, -16.19644146, -14.9175364, -13.63863134
    ]), -20.7, 9.5e-3, -1.3)
print(list(sch1_model_phi_result))
# returns
# [1.85685828e-29 3.25671116e-11 1.72458835e-05 1.27468679e-03
#  6.12395219e-03 1.26803535e-02 2.02617665e-02 2.98927403e-02
#  4.30310959e-02 6.14770529e-02]
# -----------------------


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
print(list(sch2_model_phi_list))
# returns
# [1.55110526e-18 1.09383000e-07 3.03168335e-04 3.36328048e-03
#  6.24552903e-03 6.50199270e-03 5.61245148e-03 4.55946326e-03
#  3.63199542e-03 2.87485077e-03]

# OR a rudimentarily example:
sch2_model_phi_result = lf.DoubleSchechterMagModel(
    np.array([
        -25.1487769, -23.86987184, -22.59096677, -21.31206171, -20.03315665,
        -18.75425159, -17.47534652, -16.19644146, -14.9175364, -13.63863134
    ]), -20.7, 6.16e-3, -0.79, 6.16e-3, -0.79)
print(list(sch2_model_phi_result))
# returns
# [1.94632943e-28 1.87206188e-10 5.43662993e-05 2.20369343e-03
#  5.80607779e-03 6.59304119e-03 5.77743541e-03 4.67441094e-03
#  3.69017477e-03 2.89121865e-03]
# -----------------------


m = 3
gof = lf.get_gof(phi_list, phi_err_list, sch1_model_phi_list, m)
print(gof)
# returns
# 366.43103358282144

# OR a rudimentarily example:
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
print(gof_result)
# returns
# 79.66254082924551
# -----------------------


all_sch1_model_phi_list, all_chi_sq, all_M_star, all_M_star_err, all_phi_star, all_phi_star_err, all_alpha_star, all_alpha_star_err = lf.get_schechter_phi(
    all_M_list,
    all_M_err_list,
    all_phi_list,
    all_phi_err_list,
    np.array([M_star_guess, phi_star_guess, alpha_guess]),
    plot_savename='test_all_Sch.png')
# displays plot
print(list(all_sch1_model_phi_list))
# returns
# [2.83258986e-09 5.74389144e-06 9.25939975e-05 3.12011592e-04
#  6.33385107e-04 1.09120816e-03 1.78269412e-03 2.86270542e-03
#  4.57149418e-03 7.28713338e-03]
print(all_chi_sq)
# returns
# 0.14910742282850892
print(all_M_star)
# returns
# -22.068531742285295
print(all_M_star_err)
# returns
# 0.35573470148190917
print(all_phi_star)
# returns
# 0.0003176940137059405
print(all_phi_star_err)
# returns
# 0.0001288373384458377
print(all_alpha_star)
# returns
# -1.4126892538229192
print(all_alpha_star_err)
# returns
# 0.06081125190828317


blue_sch1_model_phi_list, blue_chi_sq, blue_M_star, blue_M_star_err, blue_phi_star, blue_phi_star_err, blue_alpha_star, blue_alpha_star_err = lf.get_schechter_phi(
    blue_M_list,
    blue_M_err_list,
    blue_phi_list,
    blue_phi_err_list,
    np.array([M_star_guess, phi_star_guess, alpha_guess]),
    plot_savename='test_blue_Sch.png')
# displays plot
print(list(blue_sch1_model_phi_list))
# returns
# [1.98420348e-10 2.18093482e-06 6.21018522e-05 2.57111492e-04
#  5.70169944e-04 1.03300833e-03 1.75300151e-03 2.91245392e-03
#  4.80570005e-03 7.91206256e-03]
print(blue_chi_sq)
# returns
# 0.18163420708695324
print(blue_M_star)
# returns
# -21.842075975175316
print(blue_M_star_err)
# returns
# 0.31845816378631797
print(blue_phi_star)
# returns
# 0.0003029586014597913
print(blue_phi_star_err)
# returns
# 0.00012126827264875354
print(blue_alpha_star)
# returns
# -1.4411669183679228
print(blue_alpha_star_err)
# returns
# 0.06358938020533868

# OR a rudimentarily example:
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
    np.array([-20.71, 9.5e-3, -1.3]),
    plot_savename='test_rud_Sch.png')
# displays plot
print(list(all_sch1_model_phi_result))
# returns
# [1.71053536e-07 5.16404480e-06 4.05370711e-05 1.46583441e-04
#  3.39318246e-04 6.07660357e-04 9.38795731e-04 1.33220354e-03
#  1.80023894e-03 2.36519761e-03 3.05753557e-03 3.91587143e-03
#  4.98834428e-03 6.33496571e-03]
print(chi_sq_sch1_result)
# returns
# 1.0209802688993401
print(M_star_result)
# returns
# -22.51627500778435
print(M_star_err_result)
# returns
# 0.09643423019822513
print(phi_star_result)
# returns
# 0.0007681235644217974
print(phi_star_err_result)
# returns
# 0.00015735301981608952
print(alpha_star_result)
# returns
# -1.4248810024852225
print(alpha_star_err_result)
# returns
# 0.06007607488402875
# -----------------------


red_sch2_model_phi_list, red_chi_sq, red_M_star, red_M_star_err, red_phi_star_1, red_phi_star_err_1, red_phi_star_2, red_phi_star_err_2, red_alpha_star_1, red_alpha_star_err_1, red_alpha_star_2, red_alpha_star_err_2 = lf.get_double_schechter_phi(
    red_M_list,
    red_M_err_list,
    red_phi_list,
    red_phi_err_list,
    np.array([M_star_guess, phi_star_1_guess, alpha_1_guess, phi_star_2_guess, alpha_2_guess]),
    plot_savename='test_red_dSch.png')
# displays plot
print(list(red_sch2_model_phi_list))
# returns
# [0.00000000e+000 0.00000000e+000 0.00000000e+000 0.00000000e+000
#  0.00000000e+000 0.00000000e+000 0.00000000e+000 2.63752933e-229
#  2.35253203e-141 2.75200955e-087]
print(red_chi_sq)
# returns
# 1.2084645603920292
print(red_M_star)
# returns
# -13.256557144101373
print(red_M_star_err)
# returns
# inf
print(red_phi_star_1) 
# returns
# -0.005143924152379018
print(red_phi_star_err_1) 
# returns
# inf
print(red_phi_star_2)
# returns
# -1.8735454729853815
print(red_phi_star_err_2) 
# returns
# inf
print(red_alpha_star_1) 
# returns
# 0.012183946742584995
print(red_alpha_star_err_1) 
# returns
# inf
print(red_alpha_star_2)
# returns
# 0.025603076393042268
print(red_alpha_star_err_2)
# returns
# inf

# OR a rudimentarily example:
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
    np.array([-20.7, 6.16e-3, -0.79, 6.16e-3, -0.79]),
    plot_savename='test_rud_dSch.png')
# displays plot
print(list(all_sch2_model_phi_result))
# returns
# [8.52160254e-08 4.30479510e-06 4.25294771e-05 1.65136448e-04
#  3.77248531e-04 6.40958991e-04 9.29145543e-04 1.24659941e-03
#  1.62508333e-03 2.11836716e-03 2.80552684e-03 3.80265411e-03
#  5.28333171e-03 7.51056271e-03]
print(chi_sq_sch2_result)
# returns
# 0.8888283543610924
print(M_star_result)
# returns
# -22.303878380116704
print(M_star_err_result)
# returns
# 0.2646412794527086
print(phi_star_1_result) 
# returns
# 0.0009668887609189701
print(phi_star_err_1_result) 
# returns
# 0.000640187578339006
print(phi_star_2_result)
# returns
# -1.0900241221219484
print(phi_star_err_2_result) 
# returns
# 0.7987986322969173
print(alpha_star_1_result) 
# returns
# 0.0001418318772494868
print(alpha_star_err_1_result) 
# returns
# 0.0008399596540331241
print(alpha_star_2_result)
# returns
# -1.774506451062984
print(alpha_star_err_2_result)
# returns
# 0.9946532141625982
# -----------------------