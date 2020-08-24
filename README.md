# Luminosity Function Constructor and Modeller

This package allows the user to construct and model Galaxian Luminosity Functions using the ![1/Vmax](https://render.githubusercontent.com/render/math?math=\frac{1}{V_{max}} ) estimator and Schechter function. 

![PyPI](https://img.shields.io/pypi/v/lumfunc?color=sucess)    ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/lumfunc)    ![PyPI - Downloads](https://img.shields.io/pypi/dm/lumfunc?color=blue&label=downloads%20%E2%AC%87)    [![GitHub issues](https://img.shields.io/github/issues/manasveesaraf/lumfunc)](https://github.com/manasveesaraf/lumfunc/issues)    [![GitHub stars](https://img.shields.io/github/stars/manasveesaraf/lumfunc)](https://github.com/manasveesaraf/lumfunc/stargazers)    [![GitHub forks](https://img.shields.io/github/forks/manasveesaraf/lumfunc)](https://github.com/manasveesaraf/lumfunc/network)    [![GitHub license](https://img.shields.io/github/license/manasveesaraf/lumfunc)](https://github.com/manasveesaraf/lumfunc/blob/master/LICENSE)

## Installation

Use the package manager [pip](https://pypi.org/project/lumfunc/) to install lumfunc.

```bash
pip install lumfunc
```
Keep the package up to date to access all commands. 

```bash
pip install --upgrade lumfunc
```

## Usage

Import the module in your Python code.

```python
import lumfunc as lf
```
Load the catalogued data from survey. Usually stored in .fits or .csv files.

```python
import numpy as np
import pandas as pd

# test data (photometric galaxian survey)
data_table = pd.read_csv('test_catalogue.csv')
RA_list = np.array(data_table['RA'])
Dec_list = np.array(data_table['Dec'])
g_app_mag_list = np.array(data_table['g_mag'])
r_app_mag_list = np.array(data_table['r_mag'])
r_app_mag_err_list = np.array(data_table['r_mag_err'])
z_photo_list = np.array(data_table['z_photo'])
```


### 1. K-correction and Malmquist bias reduction:

<details><summary><b>get_maggy( )</b>: Convert the measurements of flux in magnitudes to maggies for use with <a href="https://github.com/nirinA/kcorrect_python">kcorrect_python</a></summary>
<p>

Return maggies from magnitudes.
$$ f = 10^{\frac{m}{-2.5}}$$

```python
r_maggies_list = lf.get_maggy(r_app_mag_list) 
print(r_maggies_list[0:4])
# returns 
# [2.17126084e-08 1.88972757e-08 9.39864400e-09 3.74726494e-08]

# rudimentarily:
lf.get_maggy(np.array([19.15822, 19.309002, 20.067337, 18.565714]))
# returns
# array([12.17126084e-08, 1.88972757e-08, 9.39864400e-09, 3.74726494e-08])
```

</p>
</details>

<details><summary><b>get_maggy_inv_var( )</b>: Convert the magnitude errors to maggy inverse variances for use with <a href="https://github.com/nirinA/kcorrect_python">kcorrect_python</a></summary>
<p>

Return maggy inverse variances from maggies and magnitude errors.

```python
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
```

</p>
</details>

<details><summary><b>get_rest_mag( )</b>: Convert the measured apparent magnitudes into rest-frame magnitudes using the catalogue data and output from <a href="https://github.com/nirinA/kcorrect_python">kcorrect_python</a> functions</summary>
<p>
    
Load maggy ratios output file from <a href="https://github.com/nirinA/kcorrect_python">kcorrect_python</a>.

```python
maggy_ratios_table = pd.read_csv('test_maggy_ratios.csv', delimiter=' ')
r_maggy_ratio_list = np.array(maggy_ratios_table['maggy_ratio'])
```    
Return rest-frame magnitudes from the apparent magnitudes, redshifts and maggy ratios.

```python
r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)
print(r_rest_mag_list[0:4])
# returns 
# [-22.50048222 -20.3671756  -23.61190368 -23.75133511]

# rudimentarily:
lf.get_rest_mag(np.array([0.34, 0.17, 0.61, 0.41]),
                np.array([19.15822, 19.309002, 20.067337, 18.565714]),
                np.array([0.69938735, 0.90226577, 0.43780755, 0.59193305]))
# returns
# array([-22.50048222, -20.3671756 , -23.61190369, -23.75133512])
```

</p>
</details>

<details><summary><b>get_volume( )</b>: Convert the survey area in square degrees and respective redshift of each data point into comoving volumes. So, estimate <img src="https://render.githubusercontent.com/render/math?math={V_{max}}" alt="Vmax" /> from <img src = "https://render.githubusercontent.com/render/math?math={z_{max}}" alt="Zmax" /> values</summary>
<p>

Load the ![zmax](https://render.githubusercontent.com/render/math?math=z_{max} ) file.

```python
zmax_table = pd.read_csv('test_zmax.csv', delimiter=' ')
z_max_list = np.array(zmax_table['zmax'])
```

Return comoving volume from the survey area and redshifts.

```python
survey_area = 2.5 #sq. degrees
Vmax_list = lf.get_volume(survey_area, z_max_list)
print(Vmax_list[:4])
# returns 
# [1756716.17902236  178625.22666027 2447025.54638078 2287569.96087901]

# rudimentarily:
lf.get_volume(2.5, np.array([0.50523681, 0.21884399, 0.57489149, 0.55985663]))
# returns
# array([1756716.14859094,  178625.22895137, 2447025.56779186, 2287569.99514156])
```

</p>
</details>

<details><summary><b>get_binned_phi( )</b>: Bin and weigh galaxy counts per magnitude by <img src="https://render.githubusercontent.com/render/math?math=\frac{1}{V_{max}}" alt=:"1/Vmax"></summary>
<p>

Return M, M errors and phi from the rest-frame magnitudes,  ![Vmax](https://render.githubusercontent.com/render/math?math=V_{max} ) values and number of bins.
    
```python
n_bins = 10
M_list, M_err_list, phi_list = lf.get_binned_phi(r_rest_mag_list, Vmax_list, n_bins)
print(M_list)
# returns
# [-25.1487769  -23.86987184 -22.59096677 -21.31206171 -20.03315665
#  -18.75425159 -17.47534652 -16.19644146 -14.9175364  -13.63863134]
print(M_err_list)
# returns
# [0.63945253 0.63945253 0.63945253 0.63945253 0.63945253 
#  0.63945253 0.63945253 0.63945253 0.63945253 0.63945253]
print(phi_list)
# returns 
# [2.78118218e+02 2.54476157e+02 6.57347457e-05 1.98257155e-04 4.84943102e-04 
#  1.02149157e-03 1.49165665e-03 4.54012724e-03 5.08195775e-03 6.14432455e-02]

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
```

</p>
</details>


### 2. Spatial variances of the phi, <img src="https://render.githubusercontent.com/render/math?math=\phi" alt="phi">, values:

<details><summary><b>get_patch_centers( )</b>: First, divide uniformly and randomly simulated data points over the survey area into equally distributed and equally sized patches</summary>
<p>
Load RA and Dec from uniformly distributed catalogue.
    
Return patch centers as (RA, Dec) from the uniform RA, Dec and number of patches.

```python
n_patches = 10
centers_array = lf.get_patch_centers(uniform_RA_list,
                                     uniform_Dec_list,
                                     n_patches,
                                     survey='kids',
                                     max_iterations=int(100),
                                     tolerance=1.0e-1)
print(centers_array)
# returns
# [[223.14598923   1.41920976]
#  [223.20391197   0.51178337]
#  [223.25991937  -0.79624346]
#  [223.28288439   1.0680993 ]
#  [223.21604089  -1.45187853]
#  [223.25022217  -0.2376229 ]
#  [223.22077723   2.47724817]
#  [223.22337304   1.86789658]
#  [223.21732614  -1.81036826]
#  [223.12400092  -1.13952779]]
```

</p>
</details>

<details><summary><b>get_patch_labels( )</b>: Then, use the patch centers to label the survey data points by equally distributed and equally sized patches</summary>
<p>

Return patch labels for each data point from RA, Dec, number of patches and patch center guesses.

```python
labels = lf.get_patch_labels(RA_list,
                             Dec_list,
                             n_patches,
                             centers_array,
                             survey='kids',
                             numba_installed=True,
                             plot_savename='test_patches.png')
# displays plot
```

![get_patch_labels](https://raw.githubusercontent.com/manasveesaraf/lumfunc/master/test/test_patches.png)

</p>
</details>

<details><summary><b>get_binned_phi_error( )</b>: Finally, use the patch labels to compute the spatial variances of <img src="https://render.githubusercontent.com/render/math?math=\phi" alt="phi"> </summary>
<p>

Return error on phi from rest-frame magnitude, maximum observed volume, labels, number of patches and number of bins.

```python
phi_err_list = lf.get_binned_phi_error(r_rest_mag_list, Vmax_list, labels, 10, 10)
print(phi_err_list)
# returns
# [6.51069007e+02 5.73814814e+02 5.14270184e-05 1.05945659e-04 4.14049337e-04 
#  6.03358074e-04 6.94527484e-04 4.64177159e-03 3.28583313e-03 1.60789859e-01]
```

</p>
</details>


### 3. Visualisation:

<details><summary><b>get_plot( )</b>: Perform <code>get_binned_phi()</code> , <code>get_patch_labels()</code> and <code>get_binned_phi_error()</code> functions using only one composite function and visualise the luminsoity function</summary>
<p>

Plot the ![1/Vmax](https://render.githubusercontent.com/render/math?math=\frac{1}{V_{max}} ) weighted luminosity function, binned by magnitude.

```python
M_list, M_err_list, phi_list, phi_err_list = lf.get_plot(
    r_rest_mag_list,
    Vmax_list,
    n_bins,
    RA_list,
    Dec_list,
    n_patches,
    centers_array,
    survey='kids',
    numba_installed=True,
    plot_savename='test_LF.png')

# displays plot
```

![get_plot](https://raw.githubusercontent.com/manasveesaraf/lumfunc/master/test/test_LF.png)

</p>
</details>



<details><summary><b>filter_plot_by_colour( )</b>: Study the luminosity function by colour properties by specifying the colour dichotomy</summary>
<p>

Plot the ![1/Vmax](https://render.githubusercontent.com/render/math?math=\frac{1}{V_{max}} ) weighted luminosity function from data, binned by magnitude and filtered by galaxy colours. The galaxy colours are filtered by red and blue with the help of the input colour dichotomy line parameters. The colour dichotomy line parameters must be inferred first from a CMD plot.

```python
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
    centers_array,
    survey='kids',
    numba_installed=True,
    plot_savename='test_LF_colour.png')

# displays plot
```

![filter_plot_by_colour](https://raw.githubusercontent.com/manasveesaraf/lumfunc/master/test/test_LF_colour.png)

</p>
</details>


### 4. Modelling with Schechter functions:

<details><summary><b>SchechterMagModel( )</b></summary>
<p>

Return single Schechter luminosity function in terms of magnitude from 3 free parameters of the model.

```python
M_star_guess = -20.7
phi_star_guess = 9.5e-3
alpha_guess = -1.3
sch1_model_phi_list = lf.SchechterMagModel(M_list, M_star_guess, phi_star_guess, alpha_guess)
print(sch1_model_phi_list)
# returns
# [1.85685848e-29 3.25671139e-11 1.72458831e-05 1.27468679e-03 6.12395219e-03 
#  1.26803536e-02 2.02617665e-02 2.98927403e-02 4.30310959e-02 6.14770530e-02]
```

</p>
</details>

<details><summary><b>DoubleSchechterMagModel( )</b></summary>
<p>

Return double Schechter luminosity function in terms of magnitude from 5 free parameters of the model.

```python
M_star_guess = -20.7
phi_star_1_guess = 6.16e-3
alpha_1_guess = -0.79
phi_star_2_guess = 6.16e-3
alpha_2_guess = -0.79
sch2_model_phi_list = lf.DoubleSchechterMagModel(M_list, M_star_guess,
                                                 phi_star_1_guess,
                                                 alpha_1_guess,
                                                 phi_star_2_guess,
                                                 alpha_2_guess)
print(sch2_model_phi_list)
# returns
# [1.94632963e-28 1.87206201e-10 5.43662983e-05 2.20369342e-03 5.80607779e-03 
#  6.59304119e-03 5.77743541e-03 4.67441094e-03 3.69017477e-03 2.89121864e-03]
```

</p>
</details>

<details><summary><b>get_gof( )</b>: Estimate the goodness of the fit by the reduced chi square, <img src="https://render.githubusercontent.com/render/math?math=\chi_{\nu}^{2}" alt="redchisq"></summary>
<p>

Returns reduced chi squared estimate of goodness of fit from observed values, modelled values, errors and number of free parameters used in model.

```python
m = 3
gof = lf.get_gof(phi_list, phi_err_list, sch1_model_phi_list, m)
print(gof)
# returns
# 92.50772762457441
```

</p>
</details>

<details><summary><b>get_schechter_phi( )</b>: Least square fit single Schechter function on data and plot</summary>
<p>

Returns least square fit of phi with single Schechter function, reduced chi squared estimate and the 3 Schechter parameters with their errors.

```python
sch1_model_phi_list, chi_sq_1, M_star, M_star_err, phi_star, phi_star_err, alpha_star, alpha_star_err = lf.get_schechter_phi(
    M_list,
    M_err_list,
    all_phi_list,
    all_phi_err_list,
    np.array([M_star_guess, phi_star_guess, alpha_guess]),
    plot_savename='test_Sch.png')

# displays plot
```

![get_schechter_phi](https://raw.githubusercontent.com/manasveesaraf/lumfunc/master/test/test_Sch.png)

</p>
</details>

<details><summary><b>get_double_schechter_phi( )</b>: Least square fit double Schechter function on data and plot</summary>
<p>
    
Returns least square fit of phi with double Schechter function, reduced chi squared estimate and the 5 Schechter parameters with their errors.    

```python
sch2_model_phi_list, chi_sq_1, M_star, M_star_err, phi_star_1, phi_star_err_1, phi_star_2, phi_star_err_2, alpha_star_1, alpha_star_err_1, alpha_star_2, alpha_star_err_2 = lf.get_double_schechter_phi(
    M_list,
    M_err_list,
    all_phi_list,
    all_phi_err_list,
    np.array([M_star_guess, phi_star_1_guess, alpha_1_guess, phi_star_2_guess, alpha_2_guess]),
    plot_savename='test_dSch.png')

# displays plot
```

![get_double_schechter_phi](https://raw.githubusercontent.com/manasveesaraf/lumfunc/master/test/test_dSch.png)

</p>
</details>

## Dependencies
![PyPI](https://img.shields.io/pypi/v/astropy?label=astropy)    ![PyPI](https://img.shields.io/pypi/v/numpy?label=numpy)    ![PyPI](https://img.shields.io/pypi/v/scipy?label=scipy)    ![PyPI](https://img.shields.io/pypi/v/matplotlib?label=matplotlib)  ![https://github.com/nirinA/kcorrect_python](https://img.shields.io/github/v/tag/nirinA/kcorrect_python?label=kcorrect_python)

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://github.com/manasveesaraf/LuminosityFunction/blob/master/LICENSE)
