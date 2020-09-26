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
# [2.1712608416407457e-08, 1.8897275734216393e-08, 9.398644004513803e-09, 3.747264941992267e-08]

# rudimentarily:
r_maggies_result = lf.get_maggy(np.array([19.15822, 19.309002, 20.067337, 18.565714]))
print(list(r_maggies_result[0:4]))
# returns
# [2.1712608416407457e-08, 1.8897275734216393e-08, 9.398644004513803e-09, 3.747264941992267e-08]
# -----------------------


r_maggy_inv_var_list = lf.get_maggy_inv_var(r_maggies_list, r_app_mag_err_list)
print(list(r_maggy_inv_var_list[0:4]))
# returns 
# [2.613536528041307e+20, 2.2153992456790612e+20, 2.6329570445628214e+20, 1.520308755766036e+20]

# rudimentarily:
r_maggy_inv_var_result = lf.get_maggy_inv_var(
    np.array([2.17126084e-08, 1.88972757e-08, 9.39864400e-09, 3.74726494e-08]),
    np.array([0.00309313, 0.0038601, 0.0071193, 0.00234987]))
print(list(r_maggy_inv_var_result[0:4]))
# returns
# [2.613534842093864e+20, 2.2154049929330334e+20, 2.632956307424491e+20, 1.520310051334256e+20]
# -----------------------


r_maggy_ratios_table = pd.read_csv('rest_maggy_ratios_r_ugriz_test.csv', delimiter=' ')
r_maggy_ratio_list = np.array(r_maggy_ratios_table['maggy_ratio'])

r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)
print(list(r_rest_mag_list[0:4]))
# returns 
# [-22.518710955165023, -20.36706085176511, -23.670847071363468, -23.681182442586575]

# rudimentarily:
r_rest_mag_result = lf.get_rest_mag(np.array([0.34, 0.17, 0.61, 0.41]),
                                  np.array([19.15822, 19.309002, 20.067337, 18.565714]),
                                  np.array([0.69938735, 0.90226577, 0.43780755, 0.59193305]))
print(list(r_rest_mag_result[0:4]))
# returns
# [-22.50048221280549, -20.367175595236922, -23.611903685248716, -23.751335116325944]
# -----------------------


zmax_table = pd.read_csv('zmax_test.csv', delimiter=' ')
z_max_list = np.array(zmax_table['zmax'])

survey_area = 2.5 #sq. degrees
Vmax_list = lf.get_volume(survey_area, z_max_list)
print(list(Vmax_list[:4]))
# returns 
# [1756716.170553711, 178625.22629838384, 2447025.5329312785, 2287569.9486382306]

# rudimentarily:
Vmax_result = lf.get_volume(2.5, np.array([0.50523681, 0.21884399, 0.57489149, 0.55985663]))
print(list(Vmax_result[:4]))
# returns
# [1756716.140122295, 178625.2285894822, 2447025.5543423546, 2287569.9829007764]
# -----------------------


n_bins = 10
M_list, M_err_list, phi_list = lf.get_binned_phi(r_rest_mag_list, Vmax_list, n_bins)
print(list(M_list))
# returns
# [-24.62894308598833, -23.404512811278693, -22.18008253656906, -20.955652261859427, -19.731221987149794, 
# -18.506791712440158, -17.282361437730522, -16.05793116302089, -14.833500888311256, -13.609070613601622]
print(list(M_err_list))
# returns
# [0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548164, 0.6122151373548181, 
# 0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548173, 0.6122151373548164]
print(list(phi_list))
# returns 
# [290.491672773112, 265.7977860159212, 9.55747321361586e-05, 0.0002549444468315998, 0.0006247531885987247, 
# 0.001075916507027022, 0.0019105283915257355, 0.005624556120539741, 0.0038603784240661696, 0.06417684969053933]

# OR a rudimentarily example:
M_result, M_err_result, phi_result = lf.get_binned_phi(
    np.array([-23, -21, -19, -22, -23, -23, -22, -23, -22, -22, -19, -21]),
    np.array([
        8e+08, 2e+08, 2e+07, 3e+08, 6e+08, 6e+08, 4e+08, 7e+08, 5e+08, 6e+08,
        7e+06, 1e+08
    ]), 4)
print(list(M_result))
# returns 
# [-22.5, -21.5, -20.5, -19.5]
print(list(M_err_result))
# returns 
# [0.5, 0.5, 0.5, 0.5]
print(list(phi_result))
# returns
# [1.0641166666666664e-08, 1.0289999999999998e-08, 0.0, 1.3229999999999996e-07]
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
# [array([223.25956531,   1.24525691]),
#  array([223.15175072,   2.24830099]),
#  array([223.28727953,  -1.83027433]),
#  array([223.25181969,   2.73045313]),
#  array([223.23336142,   1.68595717]),
#  array([223.36725045,   2.07319784]),
#  array([223.25230466,  -0.82353628]),
#  array([2.23246689e+02, 2.17238520e-01]),
#  array([223.25119842,   0.7459858 ]),
#  array([223.23769478,  -0.28535647])]
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
# [1, 3, 5, 6]
# -----------------------


phi_err_list = lf.get_binned_phi_error(r_rest_mag_list, Vmax_list, labels, n_patches, n_bins)
print(list(phi_err_list))
# returns
# [810.9397645708841, 607.8170004077618, 4.36417468536423e-05, 0.00019704012356457505, 0.0005486180652616492, 
# 0.0004654318611417406, 0.0005773328573042754, 0.00459036071996224, 0.00221037276729238, 0.164362438153455]

phi_err_result = lf.get_binned_phi_error(
    np.array([-23, -21, -19, -22, -23, -23, -22, -23, -22, -22, -19, -21]),
    np.array([
        8e+08, 2e+08, 2e+07, 3e+08, 6e+08, 6e+08, 4e+08, 7e+08, 5e+08, 6e+08,
        7e+06, 1e+08
    ]), np.array([1, 1, 2, 2, 3, 0, 1, 1, 2, 2, 3, 3]), 4, 4)
print(list(phi_err_result))
# returns
# [9.864941217372873e-09, 9.901557116602078e-09, 0.0, 1.5585903075108177e-07]
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
# [-24.62894308598833, -23.404512811278693, -22.18008253656906, -20.955652261859427, -19.731221987149794, 
# -18.506791712440158, -17.282361437730522, -16.05793116302089, -14.833500888311256, -13.609070613601622]
print(list(plot_M_err_list))
# returns
# [0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548164, 0.6122151373548181, 
# 0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548173, 0.6122151373548164]
print(list(plot_phi_list))
# returns 
# [290.491672773112, 265.7977860159212, 9.55747321361586e-05, 0.0002549444468315998, 0.0006247531885987247, 
# 0.001075916507027022, 0.0019105283915257355, 0.005624556120539741, 0.0038603784240661696, 0.06417684969053933]
print(list(plot_phi_err_list))
# returns
# [810.9397645708841, 607.8170004077618, 4.36417468536423e-05, 0.00019704012356457505, 0.0005486180652616492, 
# 0.0004654318611417406, 0.0005773328573042754, 0.00459036071996224, 0.00221037276729238, 0.164362438153455]
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
# [-24.62894308598833, -23.404512811278693, -22.18008253656906, -20.955652261859427, -19.731221987149794, 
# -18.506791712440158, -17.282361437730522, -16.05793116302089, -14.833500888311256, -13.609070613601622]
print(list(all_M_err_list))
# returns
# [0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548164, 0.6122151373548181, 
# 0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548173, 0.6122151373548164]
print(list(all_phi_list))
# returns 
# [290.491672773112, 265.7977860159212, 9.55747321361586e-05, 0.0002549444468315998, 0.0006247531885987247, 
# 0.001075916507027022, 0.0019105283915257355, 0.005624556120539741, 0.0038603784240661696, 0.06417684969053933]
print(list(all_phi_err_list))
# returns
# [810.9397645708841, 607.8170004077618, 4.36417468536423e-05, 0.00019704012356457505, 0.0005486180652616492, 
# 0.0004654318611417406, 0.0005773328573042754, 0.00459036071996224, 0.00221037276729238, 0.164362438153455]
print(list(red_M_list))
# returns
# [-23.749705414732908, -23.22313054306747, -22.696555671402038, -22.1699807997366, -21.643405928071164, 
# -21.11683105640573, -20.590256184740294, -20.06368131307486, -19.53710644140942, -19.01053156974399]
print(list(red_M_err_list))
# returns
# [0.2632874358327175, 0.2632874358327175, 0.2632874358327175, 0.2632874358327175, 0.2632874358327193, 
# 0.26328743583271574, 0.2632874358327193, 0.2632874358327175, 0.2632874358327175, 0.2632874358327175]
print(list(red_phi_list))
# returns 
# [5.262220153129362e-06, 1.146322895875979e-05, 2.2815766050216087e-05, 3.063244885429312e-05, 3.784760366971176e-05, 
# 6.955865011705814e-05, 4.601876296333914e-05, 4.2237548701635644e-05, 0.0001626682945797204, 5.198919360306014e-05]
print(list(red_phi_err_list))
# returns
# [1.6816801547847714e-05, 1.8848825092512034e-05, 2.1415807046810993e-05, 3.71536660000855e-05, 5.644501842774335e-05, 
# 3.6815620551692435e-05, 5.656805584781544e-05, 7.501902486045965e-05, 0.00015384519165375355, 0.00021127915259708341]
print(list(blue_M_list))
# returns
# [-24.62894308598833, -23.404512811278693, -22.18008253656906, -20.955652261859427, -19.731221987149794, 
# -18.506791712440158, -17.282361437730522, -16.05793116302089, -14.833500888311256, -13.609070613601622]
print(list(blue_M_err_list))
# returns
# [0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548164, 0.6122151373548181, 
# 0.6122151373548164, 0.6122151373548181, 0.6122151373548164, 0.6122151373548173, 0.6122151373548164]
print(list(blue_phi_list))
# returns 
# [290.491672773112, 265.7977762939234, 6.601873861989515e-05, 0.00019806224916273766, 0.0005366319861440035, 
# 0.0010535581878029952, 0.0019105283915257355, 0.005624556120539741, 0.0038603784240661696, 0.06417684969053933]
print(list(blue_phi_err_list))
# returns
# [810.939765667249, 607.8170014696902, 3.096420475680061e-05, 0.0001368281769262127, 0.0005482677763962965, 
# 0.00044155205816259833, 0.0005126210576546734, 0.004590031416295528, 0.002209831427396876, 0.16436038459664015]
# -----------------------


M_star_guess = -20.7
phi_star_guess = 9.5e-3
alpha_guess = -1.3
sch1_model_phi_list = lf.SchechterMagModel(M_list, M_star_guess, phi_star_guess, alpha_guess)
print(list(sch1_model_phi_list))
# returns
# [1.889077519830292e-19, 2.3677841853441523e-08, 0.00011664332685035941, 0.002299973977228317, 0.007591242124153264, 
# 0.014046685725374422, 0.021550818203935722, 0.03111778390239495, 0.044057921833823524, 0.061983743062742416]

# OR a rudimentarily example:
sch1_model_phi_result = lf.SchechterMagModel(
    np.array([
        -25.1487769, -23.86987184, -22.59096677, -21.31206171, -20.03315665,
        -18.75425159, -17.47534652, -16.19644146, -14.9175364, -13.63863134
    ]), -20.7, 9.5e-3, -1.3)
print(list(sch1_model_phi_result))
# returns
# [1.8568582821959972e-29, 3.256711157578086e-11, 1.7245883492815743e-05, 0.0012746867894433624, 0.006123952187528593, 
# 0.012680353540368563, 0.0202617665019816, 0.02989274030134303, 0.04303109591644922, 0.061477052920961686]
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
# [1.5511052649970585e-18, 1.0938300029533067e-07, 0.0003031683345133834, 0.0033632804762819004, 0.006245529031705679, 
# 0.006501992701830366, 0.005612451475589847, 0.0045594632577110345, 0.00363199541697254, 0.0028748507739080474]

# OR a rudimentarily example:
sch2_model_phi_result = lf.DoubleSchechterMagModel(
    np.array([
        -25.1487769, -23.86987184, -22.59096677, -21.31206171, -20.03315665,
        -18.75425159, -17.47534652, -16.19644146, -14.9175364, -13.63863134
    ]), -20.7, 6.16e-3, -0.79, 6.16e-3, -0.79)
print(list(sch2_model_phi_result))
# returns
# [1.9463294285254365e-28, 1.8720618798131103e-10, 5.436629926965087e-05, 0.0022036934269993955, 0.005806077792944941, 
# 0.006593041192333871, 0.00577743540908892, 0.004674410937670327, 0.0036901747679265527, 0.002891218646017407]
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
# [2.8325898598722236e-09, 5.743891439174078e-06, 9.259399753521627e-05, 0.00031201159158978584, 0.000633385107084906, 
# 0.0010912081614801337, 0.0017826941154962066, 0.0028627054171209425, 0.0045714941754318224, 0.007287133378452643]
print(all_chi_sq)
# returns
# 0.14910742282850892
print(all_M_star)
# returns
# -22.068531742285295
print(all_M_star_err)
# returns
# 0.3557347014819093
print(all_phi_star)
# returns
# 0.0003176940137059405
print(all_phi_star_err)
# returns
# 0.00012883733844583773
print(all_alpha_star)
# returns
# -1.4126892538229192
print(all_alpha_star_err)
# returns
# 0.06081125190828318


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
# [1.9842034834819953e-10, 2.1809348172275195e-06, 6.210185218825129e-05, 0.00025711149229255473, 0.0005701699443248206, 
# 0.0010330083316408773, 0.0017530015118458727, 0.0029124539166523883, 0.004805700054962961, 0.007912062557208274]
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
# [1.7105353558263635e-07, 5.1640448032129975e-06, 4.053707109679924e-05, 0.0001465834405778961, 0.00033931824553655666, 
# 0.0006076603569602913, 0.0009387957310397152, 0.001332203535517021, 0.0018002389419955371, 0.002365197612525618, 
# 0.003057535570029159, 0.003915871425980178, 0.004988344284986229, 0.006334965713301948]
print(chi_sq_sch1_result)
# returns
# 1.0209802688993401
print(M_star_result)
# returns
# -22.51627500778435
print(M_star_err_result)
# returns
# 0.0964342301982251
print(phi_star_result)
# returns
# 0.0007681235644217974
print(phi_star_err_result)
# returns
# 0.00015735301981608947
print(alpha_star_result)
# returns
# -1.4248810024852225
print(alpha_star_err_result)
# returns
# 0.06007607488402874
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
# [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.6375293296239083e-229, 2.3525320265182555e-141, 2.752009553823618e-87]
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
# [8.521602535554413e-08, 4.304795096021841e-06, 4.252947712862992e-05, 0.00016513644802319284, 0.00037724853104172785, 
# 0.0006409589905341704, 0.0009291455434703172, 0.001246599413378984, 0.0016250833276945204, 0.0021183671618024385, 
# 0.002805526837713822, 0.003802654108449027, 0.0052833317077602675, 0.007510562710100609]
print(chi_sq_sch2_result)
# returns
# 0.8888283543610924
print(M_star_result)
# returns
# -22.303878380116704
print(M_star_err_result)
# returns
# 0.26464127945271887
print(phi_star_1_result) 
# returns
# 0.0009668887609189701
print(phi_star_err_1_result) 
# returns
# 0.000640187578339018
print(phi_star_2_result)
# returns
# -1.0900241221219484
print(phi_star_err_2_result) 
# returns
# 0.7987986322969486
print(alpha_star_1_result) 
# returns
# 0.0001418318772494868
print(alpha_star_err_1_result) 
# returns
# 0.0008399596540331543
print(alpha_star_2_result)
# returns
# -1.774506451062984
print(alpha_star_err_2_result)
# returns
# 0.9946532141626322
# -----------------------