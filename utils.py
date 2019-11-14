#!/usr/bin/env python
# -*- coding: utf8 -*-

from __future__ import division
import os
from itertools import islice
import numpy as np
from synthetic import broadening, _read_raw_moog

kurucz95 = {'teff': (3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000,
                     6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000, 8250, 8500,
                     8750, 9000, 9250, 9500, 9750, 10000, 10250, 10500, 10750,
                     11000, 11250, 11500, 11750, 12000, 12250, 12500, 12750, 13000,
                     14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000,
                     23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000,
                     32000, 33000, 34000, 35000, 36000, 37000, 38000, 39000),
             'feh': (-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, 0.0,
                      0.1, 0.2, 0.3, 0.5, 1.0),
             'logg': (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)}

apogee_kurucz = {'teff': (3500, 3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000,
                          6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000, 8250, 8500,
                          8750, 9000, 9250, 9500, 9750, 10000, 10250, 10500, 10750,
                          11000, 11250, 11500, 11750, 12000, 12250, 12500, 12750, 13000,
                          14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000,
                          23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000),
                 'feh': (-5.0, -4.5, -4.0, -3.5, -3.0, -2.75, -2.5, -2.25, -2.0, -1.75,
                         -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.5),
                 'logg': (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)}

marcs =  {'teff': (2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400,
                   3500, 3600, 3700, 3800, 3900, 4000, 4250, 4500, 4750, 5000,
                   5250, 5500, 5750, 6000, 6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000),
          'feh': (-5.0, -4.0, -3.0, -2.5, -2.0, -1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0),
          'logg': (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)}


class GetModels:
    '''
    Find the names of the closest grid points for a given effective
    temperature, surface gravity, and iron abundance (proxy for metallicity).

    Inputs
    ------
    teff : int
      The effective temperature (K) for the model atmosphere
    logg : float
      The surface gravity (logarithmic in cgs) for the model atmosphere
    feh : float
      The metallicity for the model atmosphere
    atmtype : str
      The type of atmosphere models to use. Currently only Kurucz from '95.
    '''

    def __init__(self, teff, logg, feh, atmtype):
        self.teff = teff
        self.logg = logg
        self.feh  = feh
        self.atmtype = atmtype

        atmmodels = {'kurucz95': [kurucz95, 'kurucz95'], 'apogee_kurucz': [apogee_kurucz, 'apogee_kurucz'], 'marcs': [marcs, 'marcs']}
        if atmtype in atmmodels.keys():
            self.grid = atmmodels[atmtype][0]
        else:
            raise NotImplementedError('You request for atmospheric models: %s is not available' % atmtype)
        self.grid['teff'] = np.asarray(self.grid['teff'])
        self.grid['logg'] = np.asarray(self.grid['logg'])
        self.grid['feh']  = np.asarray(self.grid['feh'])

        # Checking for bounds in Teff, logg, and [Fe/H]
        if (self.teff < self.grid['teff'][0]) or (self.teff > self.grid['teff'][-1]):
            raise ValueError('Teff out of bounds: %s' % self.teff)
        if (self.logg < self.grid['logg'][0]) or (self.logg > self.grid['logg'][-1]):
            raise ValueError('logg out of bounds: %s' % self.logg)
        if (self.feh < self.grid['feh'][0]) or (self.feh > self.grid['feh'][-1]):
            raise ValueError('[Fe/H] out of bounds: %s' % self.feh)

    def _model_path(self, teff_model, logg_model, feh_model):
        '''Create the path for atmosphere models given Teff, logg, and [Fe/H]

        Inputs
        ------
        teff_model : int
          The Teff from the model grid
        logg_model : float
          The logg from the model grid
        feh_model : float
          The [Fe/H] from the model grid

        Output
        ------
        name : str
          The path to the atmosphere model
        '''
        name = '/home/paranoia/Software/models/%s/' % self.atmtype
        name = 'models/%s/' % self.atmtype
        if feh_model < 0:
            name += 'm%s/' % str(abs(feh_model)).replace('.', '')
        else:
            name += 'p%s/' % str(abs(feh_model)).replace('.', '')
        name += '%ig%s.' % (teff_model, str(logg_model).replace('.', ''))
        if feh_model < 0:
            name += 'm%s.gz' % str(abs(feh_model)).replace('.', '')
        else:
            name += 'p%s.gz' % str(abs(feh_model)).replace('.', '')
        return name

    def _model_exists(self, teff_model, logg_model, feh_model, upper=True):
        '''Check if a model exists. If not lower/raise Teff

        Inputs
        ------
        teff_model : int
          The Teff from the model grid
        logg_model : float
          The logg from the model grid
        feh_model : float
          The [Fe/H] from the model grid
        upper : bool
          If True, then search for Teff higher than previous. False otherwise. (Default: True)

        Outputs
        -------
        fname : str
          Path for the model
        teff_model : int
          The new Teff. Same Teff is returned if the model exists at the right place
        '''

        fname = self._model_path(teff_model, logg_model, feh_model)
        if os.path.isfile(fname):
            return fname, teff_model, logg_model
        else:
            print('Models do not exist.')
            return False
        # Change the Teff (up or down) to compensate for the gap
        teff_model0 = teff_model
        idx = np.where(teff_model == self.grid['teff'])[0][0]
        while True:
            idx = idx+1 if upper else idx-1
            try:
                teff_model = self.grid['teff'][idx]
            except IndexError:
                teff_model = teff_model0
                break
            fname = self._model_path(teff_model, logg_model, feh_model)
            if os.path.isfile(fname):
                return fname, teff_model, logg_model

        # Change logg to compensate for missing values
        idx = np.where(logg_model == self.grid['logg'])[0][0]
        while True:
            idx += 1
            logg_model = self.grid['logg'][idx]
            fname = self._model_path(teff_model, logg_model, feh_model)
            if os.path.isfile(fname):
                return fname, teff_model, logg_model

    def _looping_models(self, teff_model, logg_model, feh_model):
        models = []
        for i, teff_m in enumerate(teff_model):
            for j, logg_m in enumerate(logg_model):
                for feh_m in feh_model:
                    upper = True if self.teff < teff_m else False
                    if self._model_exists(teff_m, logg_m, feh_m, upper) is False:
                        return False
                    fname, Te, ge = self._model_exists(teff_m, logg_m, feh_m, upper)
                    teff_model[i] = Te
                    logg_model[j] = ge
                    models.append(fname)
        return models, teff_model, logg_model

    def neighbour(self, arr, val, k=2):
        '''Return the K surrounding neighbours of an array, given a certain value.

        Inputs
        ------
        arr : array_like
          The array from which some neighbours should be found (assumed sorted).
        val : float
          The value to found neighbours around in arr
        k : int
          The number of neighbours to find

        Output
        ------
        array : list
          A list with the k surrounding neighbours
        '''
        for idx, (l1, l2) in enumerate(zip(arr, islice(arr, 1, None))):
            if l1 <= val <= l2:
                break
        if k == 2:
            return [ai for ai in arr[idx:idx+2]]
        elif k == 4:
            return [ai for ai in arr[idx-1:idx+3]]

    def getmodels(self):
        '''Get the atmosphere models surrounding the requested atmospheric
        parameters. This function should be called before using the interpolation
        code within FASMA.

        output
        ------
        models : list
          List with path to 8 models the two closest in each parameter space (4x2x2)
        teff_model : list
          The four closest effective temperatures in the grid
        logg_model : list
          The two closest surface gravities in the grid
        feh_model : list
          The two closest metallicities in the grid

        The last three return values are used for the interpolation to do some
        mapping. If only the paths to the models are needed, do not pay attention
        to them.
        '''
        # Get list of parameter values
        teff_model = self.neighbour(self.grid['teff'], self.teff, k=4)
        if len(teff_model) < 4:  # Means we are close to an edge
            teff_model = self.neighbour(self.grid['teff'], self.teff, k=2)
        logg_model = self.neighbour(self.grid['logg'], self.logg, k=2)
        feh_model = self.neighbour(self.grid['feh'], self.feh, k=2)

        ratio = 1 - (logg_model[1]-self.logg)/(logg_model[1]-logg_model[0])

        teff0 = tuple(teff_model)
        logg0 = tuple(logg_model)

        if self._looping_models(teff_model, logg_model, feh_model) is False:
            return False

        models, teff_model, logg_model = self._looping_models(teff_model, logg_model, feh_model)

        if logg_model != logg0:
            if len(np.unique(logg_model)) == 1:  # We have the same values, change one
                idx = np.where(logg_model[1] == self.grid['logg'])[0][0]
                logg_model[1] = self.grid['logg'][idx+1]
            models = []
            teff_model = list(teff0)
            models, teff_model, logg_model = self._looping_models(teff_model, logg_model, feh_model)

        self.logg = logg_model[1] - (1-ratio)*(logg_model[1]-logg_model[0])
        return {'models': models, 'teff': (self.teff, teff_model),
                'logg': (self.logg, logg_model), 'feh': (self.feh, feh_model)}

def _update_par_synth(start_wave, end_wave, **kwargs):
    '''Update the parameter file (batch.par) with new linelists, atmosphere
    models, or others.

    Inputs
    -----
    atmosphere_model    :   Location of your model atmosphere file
    line_list           :   Location of your line list

    Additional keyword arguments
    ----------------------------
    These additional keyword arguments allow the user to have full control
    over what is put into the MOOG input file. The default values are:

    terminal        'x11'
    atmosphere      1
    molecules       2
    trudamp         1
    lines           1
    flux/int        1
    damping         1
    units           0
    iraf            0
    plot            0
    obspectrum      1       Unless obspectrum is provided to the function.
    opacit          0
    freeform        0
    strong          0       Unless a strong lines list is provided.
    plotpars        0       0.75 Gaussian smoothing by default. Show full
                            synthesized spectral range with y:[0, 1.2]vsini
    histogram       0
    synlimits               Defaults to the wavelength range provided and
                            the given wavelength step size, and the delta
                            defaults to the wavelength step size.

    Outputs
    -------
    And updated parameter file
    '''

    default_kwargs = {
        'plotpars':      0,
        'plot':          0,
        'atmosphere':    1,
        'molecules':     1,
        'lines':         1,
        'trudamp':       1,
        'strong':        0,
        'units':         0,
        'opacit':        0,
        'terminal':      'x11',
        'flux/int':      0,
        'obspectrum':    0,
        'model_in':     'out.atm',
        'lines_in':     'linelist.moog',
        'smoothed_out': 'smooth.out',
        'summary':      'summary.out'}
    # Fill the keyword arguments with the defaults if they don't exist already
    moog_contents = "synth\n"\
                    "terminal           %s\n"\
                    "model_in          '%s'\n"\
                    "summary_out       '%s'\n"\
                    "smoothed_out      '%s'\n"\
                    "standard_out      'result.out'\n"\
                    "lines_in          'linelist.moog'\n"\
                    "synlimits\n"\
                    "      %s      %s       %s      %s\n"\
                    "damping        %s\n" % (default_kwargs['terminal'], default_kwargs['model_in'],
                                             default_kwargs['summary'], default_kwargs['smoothed_out'],
                                             start_wave, end_wave, kwargs['options']['step_wave'],
                                             kwargs['options']['step_flux'], kwargs['options']['damping'])

    # Fill the keyword arguments with the defaults if they don't exist already
    for key, value in default_kwargs.items():
        if key not in kwargs.keys():
            kwargs[key] = value

    settings = 'plotpars,plot,atmosphere,molecules,trudamp,lines,strong,flux/int,damping,'\
               'units,iraf,opacit,freeform,observed_in,obspectrum,histogram,'\
               'synlimits'.split(',')

    # plot and plotpar values are the same
    if 'plotpars' in kwargs:
        if kwargs['plotpars'] != 0:
            settings.append('plot')
            settings.append('plotpars')

    for setting in settings:
        if setting in kwargs:
            moog_contents += "%s %s\n" % (setting + ' ' * (14 - len(setting)), kwargs[setting])

    with open('batch.par', 'w') as moog:
        moog.writelines(moog_contents)

def _run_moog(par='batch.par', driver='synth'):
    '''Run MOOGSILENT with the given parameter file

    Inputs
    ------
    par : str
      The input file for MOOG (default: batch.par)
    driver : str
      Which driver to use MOOG in. Choices are: 'abfind', 'synth'. (Default: 'synth')

    Output
    ------
      Run MOOG once in silent mode
    '''
    if driver == 'abfind':
        os.system('MOOGSILENT > /dev/null')
    elif driver == 'synth':
        with open('stupid.tmp', 'w') as f:
            f.writelines('batch.par\nq')
        os.system('MOOGSILENT < stupid.tmp > /dev/null')
        os.remove('stupid.tmp')

def fun_moog_synth(x, atmtype, abund=0.0, par='batch.par', ranges=None, results='summary.out',
                   driver='synth', version=2014, **options):
    '''Run MOOG and create synthetic spectrum for the synth driver.
    :x: A tuple/list with values (teff, logg, [Fe/H], vt, vmic, vmac)
    :atmtype: model atmosphere
    :abund: abundance of specific species
    :ranges: array of intervals for the synthesis
    :par: The parameter file (batch.par)
    :results: The summary file
    :driver: Which driver to use when running MOOG
    :version: The version of MOOG
    :returns: w, f : wavelength and flux
    '''

    from interpolation import interpolator

    fnames = ['summary.out', 'result.out']
    for i in fnames:
        try:
            os.remove(i)
        except OSError:
            pass

    # Create an atmosphere model from input parameters
    teff, logg, feh, _, vmac, vsini = x
    teff  = int(teff)
    logg  = round(logg, 3)
    feh   = round(feh, 3)
    vsini = round(vsini, 3)
    vmac  = round(vmac, 3)

    interpolator(x[0:4], atmtype=atmtype, abund=abund, elem=options['element'])

    # Create synthetic spectrum for each of the intervals
    spec = []
    for i, ri in enumerate(ranges):
        _update_par_synth(ri[0], ri[1], options=options)
        _run_moog(driver='synth')
        x_synth, y_synth = _read_raw_moog('summary.out')

        if len(x_synth) < 1:
            print('This region does not exits in synthetic spectrum: ', ri)
            print('Or the region is to wide for synthesis.')

        spec.append(broadening(x_synth, y_synth, vsini, vmac, resolution=options['resolution'], epsilon=options['limb']))

    # Gather all individual spectra to one
    w = np.column_stack(spec)[0]
    f = np.column_stack(spec)[1]
    return w, f
