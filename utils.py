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

marcs = {'teff': (2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400,
                  3500, 3600, 3700, 3800, 3900, 4000, 4250, 4500, 4750, 5000,
                  5250, 5500, 5750, 6000, 6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000),
         'feh': (-5.0, -4.0, -3.0, -2.5, -2.0, -1.5, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0),
         'logg': (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)}

kurucz08 = {'teff': (3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000,
                     6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000, 8250, 8500,
                     8750, 9000, 9250, 9500, 9750, 10000, 10250, 10500, 10750,
                     11000, 11250, 11500, 11750, 12000, 12250, 12500, 12750, 13000,
                     14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000,
                     23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000,
                     32000, 33000, 34000, 35000, 3500, 36000, 37000, 38000, 39000),
            'feh': (-4.0, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, 0.0,
                    0.1, 0.2, 0.3, 0.5, 1.0),
            'logg': (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)}


class GetModels:
    '''
    Find the names of the closest grid points for a given effective
    temperature, surface gravity, and iron abundance (proxy for metallicity).

    Inputs
    ------
    teff : int
      The effective temperature(K) for the model atmosphere
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
        self.feh = feh
        self.atmtype = atmtype

        atmmodels = {'kurucz95': [kurucz95, 'kurucz95'], 'apogee_kurucz': [apogee_kurucz, 'apogee_kurucz'], 'marcs': [marcs, 'marcs'], 'kurucz08': [kurucz08, 'kurucz08']}
        if atmtype in atmmodels.keys():
            self.grid = atmmodels[atmtype][0]
        else:
            raise NotImplementedError('You request for atmospheric models: %s is not available' % atmtype)
        self.grid['teff'] = np.asarray(self.grid['teff'])
        self.grid['logg'] = np.asarray(self.grid['logg'])
        self.grid['feh'] = np.asarray(self.grid['feh'])

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


def _update_par(atmosphere_model='out.atm', line_list='linelist.moog', **kwargs):
    '''Update the parameter file (batch.par) with new linelists, atmosphere
    models, or others.

    Inputs
    -----
    atmosphere_model : str
      Path of the model atmosphere file for MOOG
    line_list : str
      Path of the line list

    Additional keyword arguments
    ----------------------------
    These additional keyword arguments allow the user to have full control
    over what is put into the MOOG input file. This has to be a dictionary
    The default values are:

    terminal        'x11'
    atmosphere      1
    molecules       1
    trudamp         1
    lines           1
    flux/int        0
    damping         2
    units           0
    iraf            0
    plot            0
    obspectrum      0       Unless obspectrum is provided to the function.
    opacit          0
    freeform        0
    strong          0       Unless a strong lines list is provided.
    plotpars        1       0.75 Gaussian smoothing by default. Show full
                            synthesized spectral range with y:[0, 1.2]
    histogram       0
    synlimits               Defaults to the wavelength range provided and
                            the given wavelength step size, and the delta
                            defaults to the wavelength step size.

    Outputs
    -------
    And updated configuration file
    '''

    # Path checks for input files
    if not os.path.exists(line_list):
        raise IOError('Line list file "%s" could not be found.' % (line_list))

    default_kwargs = {
        'driver': 'abfind',
        'atmosphere': 1,
        'molecules': 1,
        'trudamp': 1,  # Sure, why not? It's a black art anyway!
        'lines': 1,
        'terminal': 'x11',
        'flux/int': 0,
        'damping': 1,
        'units': 0,
        'iraf': 0,
        'plot': 0,
        'obspectrum': 0,
        'opacit': 0,
        'freeform': 0,
        'strong': 0,
        'summary': 'summary.out'}

    # Fill the keyword arguments with the defaults if they don't exist already
    for key, value in default_kwargs.iteritems():
        if key not in kwargs.keys():
            kwargs[key] = value
    # Generate a MOOG-compatible run file

    moog_contents = "%s\n"\
                    "terminal       %s\n"\
                    "model_in       '%s'\n"\
                    "summary_out    '%s'\n"\
                    "standard_out   '%s'\n"\
                    "lines_in       '%s'\n" % (kwargs['driver'], kwargs['terminal'], atmosphere_model,
                                               kwargs['summary'], 'result.out', line_list)

    settings = 'atmosphere,molecules,trudamp,lines,strong,flux/int,damping,'\
               'units,iraf,plot,opacit,freeform,obspectrum,histogram,'\
               'synlimits'.split(',')
    if 'plotpars' in kwargs:
        if kwargs['plotpars'] != 0:
            settings.append('plotpars')

    for setting in settings:
        if setting in kwargs:
            moog_contents += "%s %s\n" % (setting + ' ' * (14 - len(setting)), kwargs[setting])

    with open('batch.par', 'w') as moog:
        moog.writelines(moog_contents)


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

    # Generate a MOOG-compatible run file
    # out = '%s.spec' % line_list.rpartition('/')[2]
    out = 'smooth.out'
    moog_contents = "synth\n"\
                    "terminal          %s\n"\
                    "model_in          '%s'\n"\
                    "summary_out       '%s'\n"\
                    "smoothed_out      '%s'\n"\
                    "standard_out      'result.out'\n"\
                    "lines_in          'linelist.moog'\n"\
                    "plot              %s\n"\
                    "synlimits\n"\
                    "      %s      %s       %s      %s\n"\
                    "plotpars          %s\n"\
                    "damping        %s\n" % (default_kwargs['terminal'],
                                             default_kwargs['model_in'], default_kwargs['summary'], out,
                                             kwargs['options']['plotpars'], start_wave, end_wave,
                                             kwargs['options']['step_wave'], kwargs['options']['step_flux'],
                                             kwargs['options']['plotpars'], kwargs['options']['damping'])

    # Fill the keyword arguments with the defaults if they don't exist already
    for key, value in default_kwargs.iteritems():
        if key not in kwargs.keys():
            kwargs[key] = value

    settings = 'atmosphere,molecules,trudamp,lines,strong,flux/int,damping,'\
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


def _run_moog(par='batch.par', driver='abfind'):
    '''Run MOOGSILENT with the given parameter file

    Inputs
    ------
    par : str
      The input file for MOOG (default: batch.par)
    driver : str
      Which driver to use MOOG in. Choices are: 'abfind', 'synth'. (Default: 'abfind')

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


def fun_moog(x, atmtype, par='batch.par', results='summary.out', weights='null',
             driver='abfind', version=2014, r=None, fout=None, **options):
    '''Run MOOG and return slopes for abfind mode.

    Inputs
    ------
    x : tuple
      tuple/list with values (teff, logg, [Fe/H], vt) in that order
    atmtype : str
      The atmosphere type for the interpolation
    par : str
      The configuration file for MOOG (default: batch.par)
    results : str
      The summary file of MOOG
    weights : str
      The weights to be used in the slope calculation
    driver : str
      Which driver to use when running MOOG. Choices are: 'abfind' (default), 'synth'
    version : int
      The version of MOOG (default:2014)

    Output
    ------
    driver == 'synth' :
      Run MOOG once on silent mode
    driver == 'abfind' :
      The slopes and abundances for the different elements after a run with MOOG
    '''

    from interpolation import interpolator
    # Create an atmosphere model from input parameters
    teff, logg, feh, _ = x
    _, x = interpolator(x, atmtype=atmtype, result=True)

    # Run MOOG and get the slopes and abundaces
    _run_moog(par=par, driver=driver)
    if driver == 'abfind':
        m = Readmoog(params=x, fname=results, version=version)
        _, _, _, _, _, _, data, _ = m.fe_statistics()
        if version > 2013:
            EPs, _ = slope((data[:, 2], data[:, 6]), weights=weights)
            RWs, _ = slope((data[:, 5], data[:, 6]), weights=weights)
        else:
            EPs, _ = slope((data[:, 1], data[:, 5]), weights=weights)
            RWs, _ = slope((data[:, 4], data[:, 5]), weights=weights)
        m = Readmoog(params=x, fname=results, version=version)
        fe1, _, fe2, _, _, _, _, _ = m.fe_statistics()
        abundances = [fe1+7.47, fe2+7.47]
        res = EPs**2 + RWs**2 + np.diff(abundances)[0]**2
        return res, EPs, RWs, abundances, x
    elif driver == 'synth':
        # Create synthetic spectra
        spec = []
        # Run moog for each linelist file
        # N is number of intervals
        N = len(r)
        # Run moog for each linelist file
        for i in range(N):
            _update_par_synth('linelist/%s' % fout[i], r[i][0], r[i][1], options=options)
            _run_moog(driver='synth')
            x, y = _read_raw_moog('summary.out')
            # add broadening
            # This is done here so the x-axis is equidistant.
            # Currently, the wavelength array as to be regularly spaced.
            spec.append(broadening(x, y, resolution=options['resolution'], vsini=options['vsini'], epsilon=options['limb'], vmac=options['vmac']))
        # Gather all individual spectra to one
        w = np.column_stack(spec)[0]
        f = np.column_stack(spec)[1]
        return w, f


def fun_moog_synth(x, atmtype, par='batch.par', ranges=None, results='summary.out',
                   driver='synth', version=2014, **options):
    '''Run MOOG and create synthetic spectrum for the synth driver.

    :x: A tuple/list with values (teff, logg, [Fe/H], vt, vmic, vmac)
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
    logg  =  round(logg,4)
    feh   = round(feh,4)
    vsini = round(vsini, 3)
    vmac  = round(vmac, 3)

    interpolator(x[0:4], atmtype=atmtype)

    # Create synthetic spectrum
    start = ranges[0][0]
    end = ranges[-1][-1]
    _update_par_synth(start, end, options=options)
    _run_moog(driver='synth')
    x, y = _read_raw_moog('summary.out')

    spec = []
    for i, ri in enumerate(ranges):
        x_synth = x[(x > ri[0]) & (x < ri[1])]
        y_synth = y[(x > ri[0]) & (x < ri[1])]

        # Check for enough points for vmac
        # Define central wavelength
        lambda0 = (x_synth[0] + x_synth[-1]) / 2.0
        vmacro = vmac/(299792458.*1e-3)*lambda0
        n_wave = len(x_synth)
        dwave = x_synth[1]-x_synth[0]
        n_kernel = int(5*vmacro/dwave)
        if n_kernel % 2 == 0:
            n_kernel += 1
        # The kernel might be of too low resolution, or the the wavelength range
        # might be too narrow. In both cases, raise an appropriate error
        if n_kernel > n_wave:
            # Add extra points on both sides
            ex_points = n_kernel-n_wave
            print("Spectrum range too narrow for macroturbulent broadening. Adding %s points." % ex_points)
            if ex_points % 2 == 0:
                w_s = x_synth[0] - (dwave*((ex_points+2)/2))
                w_e = x_synth[-1] + (dwave*((ex_points+2)/2))
                _update_par_synth(w_s, w_e, options=options)
                _run_moog(driver='synth')
                x_synth, y_synth = _read_raw_moog('summary.out')

            else:
                ex_points += 1
                w_s = x_synth[0] - (dwave*((ex_points+2)/2))
                w_e = x_synth[-1] + (dwave*((ex_points+2)/2))
                _update_par_synth(w_s, w_e, options=options)
                _run_moog(driver='synth')
                x_synth, y_synth = _read_raw_moog('summary.out')
        spec.append(broadening(x_synth, y_synth, vsini, vmac, resolution=options['resolution'], epsilon=options['limb']))

    # Gather all individual spectra to one
    w = np.column_stack(spec)[0]
    f = np.column_stack(spec)[1]

    return w, f


class Readmoog:
    '''Read the output file from MOOG and return some useful informations

    Inputs
    ------
    params : list/tuple
      A list of the atmospheric parameters (Teff, logg, [Fe/H], vt). If not
      provided it is read from the output file.
    fname : str
      Path of the output file (default: 'summary.out')
    version : int
      Version of MOOG to be used (default: 2014)
    '''

    def __init__(self, params=None, fname='summary.out', version=2014):
        self.fname = fname
        self.nelements = 1
        self.idx = 1 if version > 2013 else 0
        self.version = version
        with open(self.fname, 'r') as f:
            self.lines = f.readlines()
        if params:
            self.teff = params[0]
            self.logg = params[1]
            self.feh = params[2]
            self.vt = params[3]
        else:
            self.parameters()

    def parameters(self):
        '''Get the atmospheric parameters

        Outputs
        -------
        params : tuple
          The atmospheric parameters (Teff, logg, [Fe/H], vt) in that order
        '''
        for line in self.lines:
            if 'Teff' in line:
                break
        line = line.split()
        self.teff = int(line[1])
        self.logg = float(line[4])
        self.vt = float(line[6])
        self.feh = float(line[-1].split('=')[-1])
        self.params = self.teff, self.logg, self.feh, self.vt
        return self.params

    def fe_statistics(self):
        '''Get statistics on Fe lines

        Outputs
        -------
        fe1 : float
          The abundance of FeI
        sigfe1 : float
          The deviation of FeI
        fe2 : float
          The abundance of FeII
        sigfe2 : float
          The deviation of FeII
        slopeEP : float
          The correlation between abundance and excitation potential (EP)
        slopeRW : float
          The correlation between abundance and reduced equivalent width (RW)
        linesFe1 : ndarray
          The structure of FeI from the output including
          (wavelength, ID [version>=2014], EP, loggf, EW, RW, abundance, deviation on abundace)
        linesFe2 : ndarray
          Same as for linesFe1 but for FeII
        '''
        self.readdata = False
        self.slopeEP = False
        self.slopeRW = False
        self.Fe1Lines = []
        self.Fe2Lines = []
        for line in self.lines:
            if '#lines' in line and self.nelements == 1:  # Statistics on FeI
                line = line.split()
                self.readdata = False
                self.nfe1 = int(line[-1])
                self.fe1 = float(line[3])
                self.sigfe1 = float(line[7])
            elif '#lines' in line and self.nelements == 2:  # Statistics on FeII
                line = line.split()
                self.readdata = False
                self.nfe2 = int(line[-1])
                self.fe2 = float(line[3])
                self.sigfe2 = float(line[7])
            elif 'E.P.' in line and self.nelements == 1:  # We only want information from FeI
                line = line.split()
                try:
                    self.slopeEP = float(line[4])
                except ValueError:
                    self.slopeEP = False
            elif 'R.W.' in line and self.nelements == 1:  # We only want information from FeI
                line = line.split()
                self.nelements += 1  # Done with this element, move to next one
                try:
                    self.slopeRW = float(line[4])
                except ValueError:
                    self.slopeRW = False
            else:
                if line.startswith('wavelength'):
                    self.readdata = True
                    continue
            if self.readdata:
                content = map(float, filter(None, line.split(' ')))
                if self.nelements == 1:
                    self.Fe1Lines.append(content)
                else:
                    self.Fe2Lines.append(content)

        # Store the line information in numpy arrays because lists are not for science!
        self.linesFe1 = np.zeros((len(self.Fe1Lines), 7+self.idx))
        self.linesFe2 = np.zeros((len(self.Fe2Lines), 7+self.idx))
        for i, f1 in enumerate(self.Fe1Lines):
            self.linesFe1[i, 0] = f1[0]
            self.linesFe1[i, 1] = f1[1]
            self.linesFe1[i, 2] = f1[2]
            self.linesFe1[i, 3] = f1[3]
            self.linesFe1[i, 4] = f1[4]
            self.linesFe1[i, 5] = f1[5]
            self.linesFe1[i, 6] = f1[6]
            if self.version > 2013:
                self.linesFe1[i, 7] = f1[7]
        for i, f2 in enumerate(self.Fe2Lines):
            self.linesFe2[i, 0] = f2[0]
            self.linesFe2[i, 1] = f2[1]
            self.linesFe2[i, 2] = f2[2]
            self.linesFe2[i, 3] = f2[3]
            self.linesFe2[i, 4] = f2[4]
            self.linesFe2[i, 5] = f2[5]
            self.linesFe2[i, 6] = f2[6]
            if self.version > 2013:
                self.linesFe2[i, 7] = f2[7]

        # If We don't have any RW slope, calculate it manually
        if not self.slopeRW:
            self.slopeRW, _ = np.polyfit(self.linesFe1[:, 4+self.idx], self.linesFe1[:, 5+self.idx], 1)
        if not self.slopeEP:
            self.slopeEP, _ = np.polyfit(self.linesFe1[:, 1+self.idx], self.linesFe1[:, 5+self.idx], 1)
        self.sigfe1 = self.sigfe1 / np.sqrt(self.nfe1)
        try:
            self.sigfe2 = self.sigfe2 / np.sqrt(self.nfe2)
        except AttributeError:
            raise ValueError('No FeII lines were measured.')
        return self.fe1-7.47, self.sigfe1, self.fe2-7.47, self.sigfe2, self.slopeEP, self.slopeRW, self.linesFe1, self.linesFe2

    def elements(self):
        '''Get the elements and abundances from the output file

        Outputs
        -------
        element : list
          The elements, e.g. FeI, TiII
        abundances : list
          The corresponding abundances to the elements
        '''
        abundances = []
        element = []
        for line in self.lines:
            # Get the average abundance
            if line.startswith('average abundance'):
                line = filter(None, line.split('abundance =')[1].split(' '))
                abundances.append(float(line[0]))
            # Get element
            elif line.startswith('Abundance'):
                line = filter(None, line.split(' '))
                element.append(str(line[4])+str(line[5]))
        return element, abundances


def _slopeSigma(x, y, weights):
    '''Sigma on a slope after fitting a straight line

    Inputs
    ------
    x : array_like
      Independent values
    y : array_like
      Dependent values
    weights : array_like
      Weights to be appplied for linear fit

    Output
    ------
    sigma : float
      The deviation on the slope calculated
    '''
    N = len(x)
    var = np.var(x) * N
    a, b = np.polyfit(x, y, 1, w=weights)
    chi2 = np.sum((y - a*x-b)**2)
    return np.sqrt(chi2/((N-2)*var))


def error(linelist, converged, params, atmtype, version=2014, weights='null'):
    '''Error estimation on a given line list

    Inputs
    ------
    linelist : str
      Line list (without the .moog/.ares) to find the result file
    converged : bool
      True if the linelist converged, False otherwise
    params : list/tuple
      A list of the atmospheric parameters (Teff, logg, [Fe/H], vt). If not
      provided it is read from the output file.
    atmtype : str
      The atmosphere type to be used for the error calculation
    version : int
      The version of MOOG (default: 2014)
    weights : str
      The weights to be applied for slope calculation (default: 'null')

    Outputs
    -------
    teff : int
      Teff
    errorteff : int
      Error on Teff
    logg : float
      logg
    errorlogg : float
      Error on logg
    feh : float
      [Fe/H]
    errorfeh : float
      Error on [Fe/H]
    vt : float
      microturbulence
    errormicro : float
      Error on microturbulence
    '''
    # Find the output file and read the current state of it
    idx = 1 if version > 2013 else 0
    if converged:
        m = Readmoog(params=params, fname='results/%s.out' % linelist, version=version)
        summary = m.fe_statistics()
    else:
        m = Readmoog(params=params, fname='results/%s.NC.out' % linelist, version=version)
        summary = m.fe_statistics()
    # Read the correct output file (error_summary.out).
    _update_par(line_list='linelist/%s' % linelist, summary='error_summary.out')
    data = summary[6]
    _, weights = slope((data[:, 1+idx], data[:, 5+idx]), weights=weights)

    # Prepare the different things we need
    teff, logg, feh, vt = m.parameters()
    Fe1 = summary[-2]
    sigmafe1 = summary[1]
    sigmafe2 = summary[3]

    siga1 = _slopeSigma(Fe1[:, 4+idx], Fe1[:, 5+idx], weights=weights)
    siga2 = _slopeSigma(Fe1[:, 1+idx], Fe1[:, 5+idx], weights=weights)

    # Error om microturbulence
    fun_moog((teff, logg, feh, vt+0.1), atmtype, results='error_summary.out', version=version)
    sumvt = Readmoog(params=(teff, logg, feh, vt+0.1), fname='error_summary.out', version=version).fe_statistics()
    slopeEP, slopeRW = sumvt[4], sumvt[5]
    if slopeRW == 0:
        errormicro = abs(siga1/0.001) * 0.10
    else:
        errormicro = abs(siga1/slopeRW) * 0.10

    # Contribution to [Fe/H]
    deltafe1micro = abs((errormicro/0.10) * (sumvt[0]-feh))

    # Error on Teff
    slopes = errormicro/0.10 * slopeEP
    errorslopeEP = np.hypot(slopes, siga2)
    fun_moog((teff+100, logg, feh, vt), atmtype, results='error_summary.out', version=version)
    sumteff = Readmoog(params=(teff+100, logg, feh, vt), fname='error_summary.out', version=version).fe_statistics()

    errorteff = abs(errorslopeEP/sumteff[4]) * 100
    # Contribution to [Fe/H]
    deltafe1teff = abs((errorteff/100) * (sumteff[0]-feh))
    # Error on logg
    fe2error = abs(errorteff/100 * (sumteff[2]-feh))
    sigmafe2total = np.hypot(sigmafe2, fe2error)
    fun_moog((teff, logg-0.20, feh, vt), atmtype, results='error_summary.out', version=version)
    sumlogg = Readmoog(params=(teff, logg-0.20, feh, vt), fname='error_summary.out', version=version).fe_statistics()
    errorlogg = abs(sigmafe2total/(sumlogg[2]-feh)*0.20)

    # Error on [Fe/H]
    errorfeh = np.sqrt(sigmafe1**2 + deltafe1teff**2 + deltafe1micro**2)

    errorteff = int(errorteff)
    errorlogg = round(errorlogg, 2)
    errorfeh = round(errorfeh, 2)
    errormicro = round(errormicro, 2)

    os.remove('error_summary.out')
    return teff, errorteff, logg, errorlogg, feh, errorfeh, vt, errormicro


def slope(data, weights='null'):
    '''Calculate the slope of a data set with weights.

    Inputs
    ------
    data : list
      x values in first element and y values in the second element
    weights : str
      The weight to be applied for the slope. Choices:
      'null' : All points have the same weight
      'sigma' : Weight based on the standard deviation
      'mad' : Weight based on the mean absolute deviation (mad)
      For 'sigma' and 'mad' all points more than 3sig/mad away form the mean
      will have weight: 0.01. Between 2 and 3: 0.10. Between 1 and 2: 0.25.
      Below 1: 1.00

    Outputs
    -------
    params : list
      The slope calculated
    w : ndarray
      The weights used
    '''
    import statsmodels.formula.api as sm
    weights = weights.lower()
    options = ['null', 'sigma', 'mad']
    if weights not in options:
        weights = None

    data = {'x': data[0], 'y': data[1]}
    fit = np.polyfit(data['x'], data['y'], 1)
    Y = np.poly1d(fit)(data['x'])
    dif = data['y'] - Y

    if not weights:
        w = 1/abs(dif)
        idx = np.isinf(w)
        w[~idx] /= w[~idx].max()
        w[idx] = 1
    if weights == 'null':
        w = np.ones(len(data['x']))
    elif weights == 'sigma':
        sig = np.std(dif)
        w = np.zeros(len(data['y'])) + 0.01
        mask3 = abs(data['y']-Y) < 3*sig
        w[mask3] = 0.10
        mask2 = abs(data['y'] - Y) < 2*sig
        w[mask2] = 0.25
        mask1 = abs(data['y'] - Y) < sig
        w[mask1] = 1.0
    elif weights == 'mad':
        mad = np.mean(np.absolute(dif - np.mean(dif, None)), None)
        w = np.zeros(len(data['y'])) + 0.01
        mask3 = abs(data['y']-Y) < 3*mad
        w[mask3] = 0.10
        mask2 = abs(data['y'] - Y) < 2*mad
        w[mask2] = 0.25
        mask1 = abs(data['y'] - Y) < mad
        w[mask1] = 1.0

    wls = sm.wls('y ~ x', data=data, weights=w).fit()
    return wls.params[1], w
