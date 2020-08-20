import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='luminosity',  
     version='0.1',
     scripts=['luminosity'] ,
     author="Manasvee Saraf",
     author_email="saraf.manasvee@gmail.com",
     description="Galaxian Luminosity Constructor package using the 1/Vmax estimator and Schechter model.",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/manasveesaraf/LuminosityFunction",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
 )