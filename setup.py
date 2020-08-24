# python setup.py sdist bdist_wheel
# twine check dist/*
# twine upload dist/*
# pip install --upgrade lumfunc

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='lumfunc',  
     version='0.2.7',
     author="Manasvee Saraf",
     author_email="saraf.manasvee@gmail.com",
     description="Galaxian Luminosity Function Constructor package using the 1/Vmax estimator and Schechter model.",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/manasveesaraf/LuminosityFunction",
    #  packages=['lumfunc'],
     py_modules=["lumfunc","kmeans_radec"],
     package_dir={'':'lumfunc'},
     install_requires=['scipy', 'matplotlib', 'numpy', 'astropy' ],
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
 )
