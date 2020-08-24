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

Convert the measurements of flux in magnitudes to maggies for use with [kcorrect_python](https://github.com/nirinA/kcorrect_python).

<details><summary><b>get_maggy( )</b></summary>
<p>

Converts magnitudes into maggies.

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

Convert the magnitude errors to maggy inverse variances for use with [kcorrect_python](https://github.com/nirinA/kcorrect_python).

<details><summary><b>get_maggy_inv_var( )</b></summary>
<p>
    
Returns inverse variances on maggies using maggies and magnitude errors.

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

Convert the measured apparent magnitudes into rest-frame magnitudes using the catalogue data and output from [kcorrect_python](https://github.com/nirinA/kcorrect_python) functions.

```python
maggy_ratios_table = pd.read_csv('test_maggy_ratios.csv', delimiter=' ')
r_maggy_ratio_list = np.array(maggy_ratios_table['maggy_ratio'])
```

<details><summary><b>get_rest_mag( )</b></summary>
<p>
    
Converts apparent magnitudes into rest-frame magnitudes.
It uses the apparent magnitudes, redshifts and maggy ratios.

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

Convert the survey area in square degrees and respective redshift of each data point into comoving volumes. Use to estimate ![Vmax](https://render.githubusercontent.com/render/math?math={V_{max}} ) from ![zmax](https://render.githubusercontent.com/render/math?math={z_{max}} ) values.

<details><summary><b>get_volume( )</b></summary>
<p>
    
Returns comoving volume of input survey area and redshift.

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

<details><summary><b>get_binned_phi( )</b></summary>
<p>

Bins and weighs galaxy counts per magnitude implementing the 1/Vmax estimator. 
Returns phi using rest-frame magnitude, maximum observed volume and the number of bins.

</p>
</details>

<details><summary><b>get_patches_centers( )</b></summary>
<p>

Divides the input uniform random survey into equally distributed and equally sized patches. 
Returns n_patches centers as (RA,Dec) from RA, Dec and number of patches.

</p>
</details>

<details><summary><b>get_patches( )</b></summary>
<p>

Divides survey into equally distributed and equally sized patches. 
Returns labels for patches from RA, Dec, number of patches and patch center guesses.

```python
lf.get_patches(np.array([20, 21, 22, 20, 21, 22]),
            np.array([20, 21, 22, 20, 21, 22]),
            3,
            np.array([[20, 21], [22, 20], [21, 22], [20, 21], [22, 20],
                      [21, 22]]),
            survey='kids',
            numba_installed=True,
            plot_savename='test_patches.png')
# Displays the plot
```

![get_patches](https://raw.githubusercontent.com/manasveesaraf/lumfunc/master/docs/test_patches.png)

</p>
</details>

<details><summary><b>get_binned_phi_error( )</b></summary>
<p>

Spatial variance on galaxy number density per magnitude. 
Returns error on phi from rest-frame magnitude, maximum observed volume, labels, number of patches and number of bins.

</p>
</details>

<details><summary><b>plot_LF( )</b></summary>
<p>

Plots the 1/Vmax weighted luminosity function from data, binned by magnitude.

</p>
</details>
<details><summary><b>analyse_LF_by_colour( )</b></summary>
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