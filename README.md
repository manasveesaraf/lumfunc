# Luminosity Function Constructor and Modeller

This package allows the user to construct and model Galaxian Luminosity Functions using the ![1/Vmax](https://render.githubusercontent.com/render/math?math=\frac{1}{V_{max}} ) estimator and Schechter function.  

![PyPI](https://img.shields.io/pypi/v/lumfunc?color=sucess)    ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/lumfunc)    ![PyPI - Downloads](https://img.shields.io/pypi/dm/lumfunc?color=blue&label=downloads%20%E2%AC%87)    [![GitHub issues](https://img.shields.io/github/issues/manasveesaraf/lumfunc)](https://github.com/manasveesaraf/lumfunc/issues)    [![GitHub stars](https://img.shields.io/github/stars/manasveesaraf/lumfunc)](https://github.com/manasveesaraf/lumfunc/stargazers)    [![GitHub forks](https://img.shields.io/github/forks/manasveesaraf/lumfunc)](https://github.com/manasveesaraf/lumfunc/network)    [![GitHub license](https://img.shields.io/github/license/manasveesaraf/lumfunc)](https://github.com/manasveesaraf/lumfunc/blob/master/LICENSE)

## Installation

Use the package manager [pip](https://pypi.org/project/lumfunc/) to install lumfunc.

```bash
pip install lumfunc
```
To upgrade to the latest version: 

```bash
pip install --upgrade lumfunc
```

## Usage


```python
import lumfunc as lf
import numpy as np

lf.get_maggy(np.array([10, 100, 20])) # returns maggy values
```

<details><summary><b>get_maggy( )</b></summary>
<p>

Converts magnitudes into maggies.

```python
lf.get_maggy(np.array([10, 100, 20])) 
# returns array([1.e-04, 1.e-40, 1.e-08])
```

</p>
</details>

<details><summary><b>get_maggy_inv_var( )</b></summary>
<p>
    
Returns inverse variances on maggies using maggies and magnitude errors.

</p>
</details>

<details><summary><b>get_rest_mag( )</b></summary>
<p>
    
Converts apparent magnitudes into rest-frame magnitudes.
It uses the apparent magnitudes, redshifts and maggy ratios.

</p>
</details>

<details><summary><b>get_volume( )</b></summary>
<p>
    
Returns comoving volume of input survey area and redshift.

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