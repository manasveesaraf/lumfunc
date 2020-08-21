# Luminosity Function Constructor and Modeller

This packag allows the user to construct and model Galaxian Luminosity Functions using the $\frac{1}{V_{max}}$ estimator and Schechter function.  

![PyPI](https://img.shields.io/pypi/v/lumfunc?color=sucess) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/lumfunc)    ![PyPI - Downloads](https://img.shields.io/pypi/dm/lumfunc?color=blue&label=downloads%20%E2%AC%87)  [![GitHub issues](https://img.shields.io/github/issues/manasveesaraf/lumfunc)](https://github.com/manasveesaraf/lumfunc/issues) [![GitHub stars](https://img.shields.io/github/stars/manasveesaraf/lumfunc)](https://github.com/manasveesaraf/lumfunc/stargazers)   [![GitHub forks](https://img.shields.io/github/forks/manasveesaraf/lumfunc)](https://github.com/manasveesaraf/lumfunc/network)  [![GitHub license](https://img.shields.io/github/license/manasveesaraf/lumfunc)](https://github.com/manasveesaraf/lumfunc/blob/master/LICENSE)

## Installation

Use the package manager [pip](https://pypi.org/project/lumfunc/) to install lumfunc.

```bash
pip install lumfunc
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

<details><summary><b>get_patches( )</b></summary>
<p>

Divides survey into equally distributed and equally sized patches. Returns labels for patches from RA, Dec, number of patches and patch center guesses.

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

![get_patches](docs/test_patches.png)

</p>
</details>

## Dependencies
![PyPI](https://img.shields.io/pypi/v/astropy?label=astropy)    ![PyPI](https://img.shields.io/pypi/v/numpy?label=numpy)    ![PyPI](https://img.shields.io/pypi/v/scipy?label=scipy)    ![PyPI](https://img.shields.io/pypi/v/matplotlib?label=matplotlib)

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://github.com/manasveesaraf/LuminosityFunction/blob/master/LICENSE)