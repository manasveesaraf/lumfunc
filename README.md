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
r_app_mag_list = np.array(data_table['r_mag'])
r_app_mag_err_list = np.array(data_table['r_mag_err'])
z_photo_list = np.array(data_table['z_photo'])
```

Convert the measurements of flux in magnitudes to maggies for use with [kcorrect_python](https://github.com/nirinA/kcorrect_python):

<details><summary><b>get_maggy( )</b></summary>
<p>

Return maggies from magnitudes.

```python
r_maggies_list = lf.get_maggy(r_app_mag_list) 
print(r_maggies_list[0:4])
# returns 
# [1.83315843e-08 2.27614539e-08 1.33659552e-08 1.13031632e-07]

# rudimentarily:
lf.get_maggy(np.array([19.342, 19.107, 19.685, 17.367]))
# returns
# array([1.83315843e-08, 2.27614539e-08, 1.33659552e-08, 1.13031632e-07])
```

</p>
</details>

Convert the magnitude errors to maggy inverse variances for use with [kcorrect_python](https://github.com/nirinA/kcorrect_python):

<details><summary><b>get_maggy_inv_var( )</b></summary>
<p>

Return maggy inverse variances from maggies and magnitude errors.

```python
r_maggy_inv_var_list = lf.get_maggy_inv_var(r_maggies_list, r_app_mag_err_list)
print(r_maggy_inv_var_list[0:4])
# returns 
# [2.19244475e+20 5.68838063e+20 4.12409497e+20 9.22674759e+19]

# rudimentarily:
lf.get_maggy_inv_var(np.array([1.83315843e-08, 2.27614539e-08, 1.33659552e-08, 1.13031632e-07]),
                     np.array([0.004, 0.002, 0.004, 0.001]))
# returns
# array([2.19244474e+20, 5.68838064e+20, 4.12409494e+20, 9.22674766e+19])
```

</p>
</details>

Convert the measured apparent magnitudes into rest-frame magnitudes using the catalogue data and output from [kcorrect_python](https://github.com/nirinA/kcorrect_python) functions:

<details><summary><b>get_rest_mag( )</b></summary>
<p>
    
Load maggy ratios output file from [kcorrect_python](https://github.com/nirinA/kcorrect_python).

```python
maggy_ratios_table = pd.read_csv('test_maggy_ratios.csv', delimiter=' ')
r_maggy_ratio_list = np.array(maggy_ratios_table['maggy_ratio'])
```    
Return rest-frame magnitudes from the apparent magnitudes, redshifts and maggy ratios.

```python
r_rest_mag_list = lf.get_rest_mag(z_photo_list, r_app_mag_list, r_maggy_ratio_list)
print(r_rest_mag_list[0:4])
# returns 
# [-22.89979359 -21.51881811 -23.02717126 -20.79614551]

# rudimentarily:
lf.get_rest_mag(np.array([0.42, 0.24, 0.46, 0.09]),
                np.array([19.342, 19.107, 19.685, 17.367]),
                np.array([0.67165941, 0.81335927, 0.54066526, 0.91925443]))
# returns
# array([-22.8997936 , -21.51881811, -23.02717126, -20.79614551])
```

</p>
</details>

Convert the survey area in square degrees and respective redshift of each data point into comoving volumes. So, estimate ![Vmax](https://render.githubusercontent.com/render/math?math={V_{max}} ) from ![zmax](https://render.githubusercontent.com/render/math?math={z_{max}} ) values:

<details><summary><b>get_volume( )</b></summary>
<p>

Return comoving volume from the survey area and redshifts.

```python
survey_area = 100.0 #sq. degrees
V_list = lf.get_volume(survey_area, z_photo_list)
print(V_list[:4])
# returns 
# [43208407.50293904 9274338.02683353 54988309.45363603 546254.32632565]

# rudimentarily:
lf.get_volume(100.0, np.array([0.42, 0.24, 0.46, 0.09]))
# returns
# array([43208407.50293904, 9274338.02683353, 54988309.45363603, 546254.32632565])
```

</p>
</details>

Bin and weigh galaxy counts per magnitude by ![1/Vmax](https://render.githubusercontent.com/render/math?math=\frac{1}{V_{max}} ):

<details><summary><b>get_binned_phi( )</b></summary>
<p>

Return M, M errors and phi from the rest-frame magnitudes, ![Vmax](https://render.githubusercontent.com/render/math?math={V_{max}} ) values and number of bins.
    
```python
n_bins = 10
M_list, M_err_list, phi_list = lf.get_binned_phi(r_rest_mag_list, V_list, n_bins)
print(M_list)
# returns
# [-27.75116273 -26.26581137 -24.78046    -23.29510864 -21.80975727
   -20.32440591 -18.83905454 -17.35370318 -15.86835182 -14.38300045]
print(M_err_list)
# returns
# [0.74267568 0.74267568 0.74267568 0.74267568 0.74267568 
   0.74267568 0.74267568 0.74267568 0.74267568 0.74267568]
print(phi_list)
# returns 
# [5.12016808e-10 0.00000000e+00 6.87358202e-08 3.55674570e-06 1.18791217e-05 
   2.44735150e-05 5.43431411e-05 1.30067824e-04 1.04554476e-04 1.74886746e-03]

# OR a rudimentarily example:
lf.get_binned_phi(
    np.array([-23, -21, -19, -22, -23, -23, -22, -23, -22, -22, -19, -21]),
    np.array([
        8e+08, 2e+08, 2e+07, 3e+08, 6e+08, 6e+08, 4e+08, 7e+08, 5e+08, 6e+08,
        7e+06, 1e+08
    ]), 4)
# returns 
# (array([-22.5, -21.5, -20.5, -19.5]),
   array([0.5, 0.5, 0.5, 0.5]),
   array([1.06411667e-08, 1.02900000e-08, 0.00000000e+00, 1.32300000e-07]))
```

</p>
</details>

To get spatial variances of the phi (![phi](https://render.githubusercontent.com/render/math?math=\phi )) values, first 
divide uniformly and randomly simulated data points over the survey area into equally distributed and equally sized patches: 

<details><summary><b>get_patch_centers( )</b></summary>
<p>

Return patch centers as (RA, Dec) from the RA, Dec and number of patches.

```python
n_patches = 10
centers_array = lf.get_patch_centers(RA_list,
                                     Dec_list,
                                     n_patches,
                                     survey='kids',
                                     max_iterations=int(100),
                                     tolerance=1.0e-1)
print(centers_array)
# returns
# [[ 1.38832190e+02 -1.00733144e+00]
#  [ 2.17105380e+02  1.08365630e+00]
#  [ 1.80666296e+02 -2.73070692e-01]
#  [ 1.34335764e+02  1.31532218e-01]
#  [ 1.38831715e+02  2.15292944e+00]
#  [ 1.29005160e+02  1.01211250e+00]
#  [ 2.13883209e+02 -1.52070351e-02]
#  [ 1.32326750e+02  2.01815821e+00]
#  [ 2.21141020e+02  4.73369162e-01]
#  [ 1.38831187e+02  5.23810834e-01]]
```

</p>
</details>

Then use the patch centers to label the survey data points by equally distributed and equally sized patches: 

<details><summary><b>get_patch_labels( )</b></summary>
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

![get_patches](https://raw.githubusercontent.com/manasveesaraf/lumfunc/master/test/test_patches.png)

</p>
</details>

Using the patch labels, lastly compute the spatial variances of ![phi](https://render.githubusercontent.com/render/math?math=\phi ):

<details><summary><b>get_binned_phi_error( )</b></summary>
<p>

Return error on phi from rest-frame magnitude, maximum observed volume, labels, number of patches and number of bins.

```python
phi_err_list = lf.get_binned_phi_error(r_rest_mag_list, V_list, labels, 10, 10)
print(phi_err_list)
# returns
# [3.03839559e-06 7.40731159e-06 9.37491641e-06 1.52090965e-05
#  3.56343615e-05 5.44297508e-05 4.18036097e-05 1.39310857e-04
#  2.08627224e-04 3.58080092e-03]
```

</p>
</details>

Instead, perform get_binned_phi(), get_patches_centers(), get_patches() and get_binned_phi_error() functions using only one function and visualise the luminsoity function:

<details><summary><b>get_plot( )</b></summary>
<p>

Plot the ![1/Vmax](https://render.githubusercontent.com/render/math?math=\frac{1}{V_{max}} ) weighted luminosity function, binned by magnitude.

```python
M_list, M_err_list, phi_list, phi_err_list = lf.get_plot(
    r_rest_mag_list,
    V_list,
    10,
    RA_list,
    Dec_list,
    10,
    centers_array,
    survey='kids',
    numba_installed=True,
    plot_savename='test_LF.png')

# displays plot
```

![plot_LF](https://raw.githubusercontent.com/manasveesaraf/lumfunc/master/test/test_LF.png)

</p>
</details>
<details><summary><b>filter_plot_by_colour( )</b></summary>
<p>

Plots the 1/Vmax weighted luminosity function from data, binned by magnitude and filtered by galaxy colours. The galaxy colours are filtered by red and blue with the help of the input colour dichotomy line parameters. The colour dichotomy line parameters can be inferred from a CMD plot.

</p>
</details>

<details><summary><b>SchechterMagModel( )</b></summary>
<p>

Single Schechter luminosity function in terms of magnitude from 3 free parameters of the model.

</p>
</details>

<details><summary><b>DoubleSchechterMagModel( )</b></summary>
<p>

Double Schechter luminosity function in terms of magnitude from 5 free parameters of the model.

</p>
</details>

<details><summary><b>get_gof( )</b></summary>
<p>

Returns reduced chi squared estimate of goodness of fit.

</p>
</details>

<details><summary><b>get_schechter_phi( )</b></summary>
<p>

Least square fits single Schechter function model on data.
Returns best fit phi, reduced chi squared estimate and the 3 Schechter parameters with their errors.

</p>
</details>

<details><summary><b>get_double_schechter_phi( )</b></summary>
<p>
    
Least square fits double Schechter function model on data.
Returns best fit phi, reduced chi squared estimate and the 5 Schechter parameters with their errors.    

</p>
</details>

## Dependencies
![PyPI](https://img.shields.io/pypi/v/astropy?label=astropy)    ![PyPI](https://img.shields.io/pypi/v/numpy?label=numpy)    ![PyPI](https://img.shields.io/pypi/v/scipy?label=scipy)    ![PyPI](https://img.shields.io/pypi/v/matplotlib?label=matplotlib)

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://github.com/manasveesaraf/LuminosityFunction/blob/master/LICENSE)