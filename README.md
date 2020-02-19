![My image](https://github.com/MariaTsantaki/FASMA-synthesis/blob/master/img/running_icon.png)


# FASMA 2.0
Stellar spectral analysis package.

FASMA delivers the atmospheric stellar parameters (effective temperature, surface gravity, metallicity, microturbulence, macroturbulence, and rotational velocity) based on the spectral synthesis technique. The principle of this technique relies on the comparison of synthetic spectra with observations to yield the best-fit parameters under a &chi;<sup>2</sup> minimization process. FASMA 2.0 is now updated to deliver also chemical abundances of 12 elements (Na, Mg, Al, Si, Ca, Sc, Ti, V, Cr, Mn, and Ni).


The python code runs with python 3 and is wrapped around the spectral synthesis package in fortran: [MOOG version 2019](http://www.as.utexas.edu/~chris/moog.html) (Sneden et al. 1973).

The spectral synthesis technique requires model atmospheres for the calculation of the synthetic spectra by MOOG. FASMA includes two grids of models: Kurucz and marcs.

Contact here for bugs: tsantaki@arcetri.astro.it

# Installation
Installing FASMA requires a few simple steps. MOOG will also be installed to create the synthetic spectra.

Run the installation script `./install_fasma.sh`. The user will be asked about system requirements:
"rh64" for 64-bit linux system, "rh" for 32-bit linux system, "maclap" for mac laptop, "macdesk" for mac desktop.

A list of basic python packages will be installed automatically with pip (see [requirements](https://github.com/MariaTsantaki/FASMA-synthesis/blob/master/requirements.txt)).

More information about MOOG is [here](http://www.as.utexas.edu/~chris/moog.html).

If you do not have sudo access, replace the pip command with:

```
pip install --user .
```

To uninstall FASMA:

```
pip uninstall FASMA
```

# Usage
FASMA is so easy. You can run FASMA from the terminal by configuring the options in the `config.yml` file and then run in the working directory:

```
fasma
```

A small tutorial is given [here](https://github.com/MariaTsantaki/FASMA-synthesis/blob/master/manual/Manual_fasma.pdf)
FASMA includes a log file `fasma.log` every time it is used which catches errors in order to inform the user.

For large lists of stars, it is preferable to use the terminal version of FASMA.
The configuration options are added from a dictionary.

```
from FASMA import FASMA

options = {'observations': '/home/FASMA-synthesis/FASMA/spectra/Sun_HARPS.fits',
           'minimize': True,
           'plot':True}
result = FASMA(**options)
```

The output is a dictionary with the final parameters, and can be saved to a file (appended to previous results if needed).

Note that the input spectra must have a specific format and should be corrected for radial velocity shifts. This is an example on how to correct for such [RV shifts](https://github.com/MariaTsantaki/RV_correction).

## Configuration file

The complete list of options of the configuration file has the following form with the correct paths for the input files (see [manual](https://github.com/MariaTsantaki/FASMA-synthesis/blob/master/manual/Manual_fasma.pdf) for examples):

```
star:
  teff: 5777
  logg: 4.44
  feh: 0.0
  vmac: 3.21
  vsini: 1.9
  vt: 1.0
  MOOGv: 2014
  damping: 1
  element: Na
  fix_feh: False
  fix_logg: False
  fix_teff: False
  fix_vmac: False
  fix_vsini: False
  fix_vt: False
  intervals_file: /home/user/FASMA-synthesis/FASMA/rawLinelist/intervals_elements.lst
  limb: 0.6
  linelist: /home/user/FASMA-synthesis/FASMA/rawLinelist/elements.lst
  minimize: False
  model: apogee_kurucz
  observations: /home/user/FASMA-synthesis/FASMA/spectra/Sun_HARPS.fits
  plot: False
  plot_res: False
  refine: False
  resolution: null
  save: False
  snr: null
  step_flux: 3.0
  step_wave: 0.01
```

The configuration file is space sensitive and it can include only the values which
are other than the defaults.

FASMA includes the line list for the derivation of stellar parameters (`linelist.lst`), the list of the spectral regions for the spectral synthesis (`intervals.lst`) as tested in [Tsantaki et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.473.5066T/abstract). For the chemical abundance determination, FASMA incorporates the line list of 12 elements (`elements.lst`) and the intervals (`intervals_elements.lst`) from [Adibekyan et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015A%26A...583A..94A/abstract). The above lists are in the `rawLinelist` folder. The correct paths should be provided in the `config.yml` file.

# Outputs

The delivered stellar parameters are saved in the `synthresults.dat` and for the chemical abundances in the `synthresults_elements.dat`.

For a set of parameters (Teff, logg, [M/H]), a model atmosphere is obtained by interpolating from the model grid. The model atmosphere is created in the `out.atm` file. It is a tabulated form of the basic physical properties (temperature, pressure, etc.) at each layer of the atmosphere.

The raw synthetic spectrum created by MOOG is saved in `summary.out`.

The `linelist.moog` contains the line list in the format of MOOG and `result.out` is a summary of all the parameters used for the synthesis of the spectrum.

The MOOG configuration file is `batch.par` and its options are also set from the `config.yml` file.

# AUTHORS

   * [M. Tsantaki](https://github.com/MariaTsantaki)
   * [D. Andreasen](https://github.com/DanielAndreasen)
   * [G. Teixeira](https://github.com/gdcteixeira)
