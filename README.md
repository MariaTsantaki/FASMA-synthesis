![My image](https://github.com/MariaTsantaki/FASMA-synthesis/blob/master/img/running_icon.png)


# FASMA 2.0
Stellar spectral analysis package. FASMA delivers the atmospheric stellar parameters (effective temperature, surface gravity, metallicity, microturbulence, macroturbulence, rotational velocity) and chemical abundances of 12 elements.

The python code is wrapped around the spectral synthesis package in fortran: [MOOG](http://www.as.utexas.edu/~chris/moog.html).

Contact here for bugs: tsantaki@arcetri.astro.it

# Installation
Installing FASMA requires a few simple steps.

Run the makefile file with
```
make install
```

A list of basic packages will be installed automatically with pip (see requirements.txt). It requires sudo privilege otherwise correct the makefile with pip install --user instead of sudo pip install.

FASMA requires MOOG which creates the synthetic spectra and it is installed separately from FASMA but it is provided in
the `FASMA/MOOG` folder. In case MOOG is not installed, edit line 29 of the `Moogsilent.f` in the `FASMA/MOOG` folder.
Then depending on the system, compile:

```
make -f Makefile.xxx           
```

where xxx = "rh64silent", "rhsilent", "maclapselent", "macdesksilent".

More instructions are [here](http://www.as.utexas.edu/~chris/moog.html). FASMA runs with python 3.

# Usage
FASMA is so easy. You can run FASMA from the terminal by configuring the options in the `StarMe_synth.cfg` file
and then run in the working directory:

```
fasma
```

A small tutorial is given [here](https://github.com/MariaTsantaki/FASMA-synthesis/blob/master/manual/Manual_fasma.pdf)
FASMA includes a log file `captain.log` every time it is used which catches errors in order to inform the user.

## Configuration file

A standard setting of the configuration file has this form with the correct paths (see [manual](https://github.com/MariaTsantaki/FASMA-synthesis/blob/master/manual/Manual_fasma.pdf) for more information):

`linelist teff logg [M/H] vt vmac vsini options`

```
\home\star\rawLinelist\linelist.lst 5777 4.44 0.0 1.0 3.21 1.9 observations:\home\star\spectra\Sun_HARPS.fits,resolution:115000,inter_file:\home\star\rawLinelist\intervals.lst,minimize,refine
```

The line list and the initial stellar parameters (teff logg [M/H] vt vmac vsini) are space separated. The options are comma separated (e.g. minimize,damping:1,plot).
The default options of FASMA can be changed in the configuration file `StarMe_synth.cfg`.

```
spt:          False
model:        marcs
MOOGv:        2017
save:         False
fix_teff:     False
fix_logg:     False
fix_feh:      False
fix_vt:       False
fix_vmac:     False
fix_vsini:    False
plot:         False
plot_res:     False
damping:      1
step_wave:    0.01
step_flux:    3.0
minimize:     False
element:      False
refine:       False
observations: False
inter_file:   intervals.lst
snr:          None
resolution:   None
limb:         0.6
```

FASMA includes the line list for the derivation of stellar parameters (`giraffe_sun_arcturus_calib.lst`), the list of the spectral regions for the spectral synthesis (`intervals.lst`) as tested in [Tsantaki et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.473.5066T/abstract). For the chemical abundance determination, FASMA incorporates the line list of 12 elements (`elements.lst`) and the intervals (`intervals_elements.lst`) from [Adibekyan et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015A%26A...583A..94A/abstract) . The above lists are in the `rawLinelist` folder. The correct paths should be provided in the `StarMe_synth.cfg` file.

# Outputs

Assuming FASMA works well, the delivered stellar parameters are saved in the `synthresults.dat` and for the
chemical abundances in the `synthresults_elements.dat`. The model atmosphere produced by the interpolation from
the model grid is created in the `out.atm` file. The output synthetic spectrum of MOOG is in `summary.out`. The
`linelist.moog` contains the line list in the format of MOOG and `result.out` is the summary of the parameters used
for the synthetic spectrum. The MOOG configuration file is `batch.par` and its options are also set from the  `StarMe_synth.cfg` file.

# AUTHORS

   * [M. Tsantaki](https://github.com/MariaTsantaki)
   * [D. Andreasen](https://github.com/DanielAndreasen)
   * [G. Teixeira](https://github.com/gdcteixeira)
