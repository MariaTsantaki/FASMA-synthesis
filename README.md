![My image](https://github.com/MariaTsantaki/FASMA-synthesis/blob/master/img/running_icon.png)


# FASMA 2.0
Stellar spectral analysis package. 

The python code is wrapped around the spectral synthesis package: [MOOG](http://www.as.utexas.edu/~chris/moog.html).

To install, simply run the makefile file with `make install`.

**WARNING**: Contact here for bugs: tsantaki@arcetri.astro.it

# Installation
Installing `FASMA-synthesis` requires a few simple steps. We highly recommend
using [Anaconda](https://www.continuum.io/) to manage your python packages.
If you do not have Anaconda installed, you need to change [line 13 - 14](https://github.com/MariaTsantaki/FASMA-synthesis/blob/master/makefile#L13-L14) in the `makefile`. Otherwise, just run `make`
to install all dependencies.
Now FASMA runs with python 3.

# Usage
FASMA is so easy. You can run FASMA 1) either with the GUI or 2) from the terminal. 

1) Run to open the GUI control:

```
python FASMA.py
```

2) Run the CLI version.

Add the options in the `StarMe_synth.cfg` file manually and then:

```
python synthDriver.py
```


A small tutorial is given [here](https://github.com/MariaTsantaki/FASMA-synthesis/blob/master/manual/Manual_fasma.pdf)

## Configuration file

A standard setting of the configuration file has this form:

`linelist teff logg [M/H] vt vmac vsini options`

```
linelist.lst 5777 4.44 0.0 1.0 3.21 1.9 observations:Sun_HARPS.fits,resolution:115000,minimize,refine
```

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
flag_vt:      False
flag_vmac:    False
plot:         False
plot_res:     False
damping:      1
step_wave:    0.01
step_flux:    3.0
minimize:     False
refine:       False
errors:       False
observations: False
inter_file:   intervals_hr10_15n.lst
snr:          None
resolution:   None
limb:         0.6
```

# AUTHORS

   * [M. Tsantaki](https://github.com/MariaTsantaki)
   * [D. Andreasen](https://github.com/DanielAndreasen)
   * [G. Teixeira](https://github.com/gdcteixeira)
