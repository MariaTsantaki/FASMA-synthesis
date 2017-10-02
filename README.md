![My image](https://github.com/MariaTsantaki/FASMA-synthesis/blob/master/img/running_icon.png)


# FASMA-synthesis
Python code for stellar spectral analysis.

Spectral synthesis around [MOOG](http://www.as.utexas.edu/~chris/moog.html).

To install, simply run the makefile file with `make install`.

**WARNING**: Contact here for bugs: m.tsantaki@crya.unam.mx

# Installation
Installing `FASMA-synthesis` requires a few simple steps. We highly recommend
using [Anaconda](https://www.continuum.io/) to manage your python packages.
If you do not have Anaconda installed, you need to change [line 13 - 14](https://github.com/MariaTsantaki/FASMA-synthesis/blob/master/makefile#L13-L14) in the `makefile`. Otherwise, just run `make`
to install all dependencies.

# Usage
Place your spectra in the `spectra` folder and run `python FASMA.py`
to open the GUI control. It is possible to create the `StarMe_synth.cfg`
manually and run the CLI version of FASMA with `python synthDriver.py`.

## Configuration file

A standard setting of the configuration file has this form:

`linelist teff logg [M/H] vt vmac vsini options`

```
giraffe_sun_arcturus_calib.lst 5777 4.44 0.0 1.0 3.21 1.9 observations:sun.fits,resolution:115000,minimize,refine
```

The default options of FASMA can be changed in the configuration file `StarMe_synth.cfg`.

```
'spt':          False
'model':        'kurucz95'
'MOOGv':        2014
'plotpars':     0
'save':         False
'fix_teff':     False
'fix_logg':     False
'fix_feh':      False
'fix_vt':       False
'fix_vmac':     False
'fix_vsini':    False
'flag_vt':      False
'flag_vmac':    False
'plot':         False
'plot_res':     False
'damping':      1
'step_wave':    0.01
'step_flux':    3.0
'minimize':     False
'refine':       False
'errors':       False
'observations': False
'inter_file':   'intervals_hr10_15n.lst'
'snr':          None
'resolution':   None
'limb':         0.6
```

# AUTHORS

   * [M. Tsantaki](https://github.com/MariaTsantaki)
   * [D. Andreasen](https://github.com/DanielAndreasen)
   * [G. Teixeira](https://github.com/gdcteixeira)

# LICENCE

FASMA uses the MIT licence.
Copyright © 2015 Maria Tsantaki, Daniel Andreasen, and Guilherme Teixeira.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

