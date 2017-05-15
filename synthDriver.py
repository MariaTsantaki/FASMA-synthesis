#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os
import yaml
from utils import fun_moog_synth as func
from observations import read_obs_intervals, plot, snr
from minimization import minimize_synth, getMic,  getMac
from synthetic import read_linelist, save_synth_spec
import numpy as np


def _getSpt(spt):
    """Get the spectral type from a string like 'F5V'."""
    if len(spt) > 4:
        raise ValueError('Spectral type most be of the form: F8V')
    if '.' in spt:
        raise ValueError('Do not use half spectral types as %s' % spt)
    with open('SpectralTypes.yml', 'r') as f:
        d = yaml.safe_load(f)
    temp = spt[0:2]
    lum = spt[2:]
    try:
        line = d[lum][temp]
    except KeyError:
        print('Was not able to find the spectral type: %s' % spt)
        print('Teff=5777 and logg=4.44')
        return 5777, 4.44
    try:
        line = line.split()
        teff = int(line[0])
        logg = float(line[1])
    except AttributeError:
        teff = line
        logg = 4.44
    return teff, logg


def _options(options=None):
    '''Reads the options inside the config file'''
    defaults = {'spt':          False,
                'model':        'kurucz95',
                'MOOGv':        2014,
                'plotpars':     0,
                'save':         False,
                'fix_teff':     False,
                'fix_logg':     False,
                'fix_feh':      False,
                'fix_vt':       False,
                'fix_vmac':     False,
                'fix_vsini':    False,
                'flag_vt':      False,
                'flag_vmac':    False,
                'plot':         False,  # This is irrelevant with the batch.par value
                'plot_res':     False,
                'damping':      1,
                'step_wave':    0.01,
                'step_flux':    3.0,
                'minimize':     False,
                'refine':       False,
                'errors':       False,
                'observations': False,
                'inter_file':   'intervals.lst',
                'snr':          None,
                'resolution':   None,
                'limb':         0.6
                }
    if not options:
        return defaults
    else:
        for option in options.split(','):
            if ':' in option:
                option = option.split(':')
                defaults[option[0]] = option[1]
            else:
                # Clever way to change the boolean
                if option in ['teff', 'logg', 'feh', 'vt', 'vmac', 'vsini']:
                    option = 'fix_%s' % option
                defaults[option] = False if defaults[option] else True
        defaults['model']     = defaults['model'].lower()
        defaults['step_wave'] = float(defaults['step_wave'])
        defaults['step_flux'] = float(defaults['step_flux'])
        defaults['limb']      = float(defaults['limb'])
        defaults['MOOGv']     = int(defaults['MOOGv'])
        if defaults['observations'] and (defaults['snr'] is None):
            if os.path.isfile('spectra/%s' % defaults['observations']):
                defaults['snr'] = snr('spectra/%s' % defaults['observations'])
            elif os.path.isfile(str(defaults['observations'])):
                defaults['snr'] = snr(defaults['observations'])
        if defaults['inter_file']:
            defaults['inter_file'] = str(defaults['inter_file'])
        else:
            defaults['inter_file'] = 'intervals.lst'
        if defaults['resolution'] is not None:
            defaults['resolution'] = int(defaults['resolution'])
        return defaults


def _output(overwrite=None, header=None, parameters=None):
    """Create the output file 'synthresults.dat'

    Input
    -----
    overwrite - Overwrite the file
    header    - Only use True if this is for the file to be created
    """
    if header:
        hdr = ['linelist', 'observations', 'teff', 'tefferr', 'logg', 'loggerr', 'feh', 'feherr',
               'vt', 'vterr', 'vmac', 'ervmac', 'vsini', 'ervsini', 'chi2', 'status', 'errteff', 'errlogg',
               'errfeh', 'errvsini', 'time', 'model', 'resolution', 'snr']
        if overwrite:
            with open('synthresults.dat', 'w') as output:
                output.write('\t'.join(hdr)+'\n')
        else:
            if not os.path.isfile('synthresults.dat'):
                with open('synthresults.dat', 'w') as output:
                    output.write('\t'.join(hdr)+'\n')
    else:
        with open('synthresults.dat', 'a') as output:
            output.write('\t'.join(map(str, parameters))+'\n')


def synthdriver(starLines='StarMe_synth.cfg', overwrite=False):
    """The function that glues everything together

    Input
    -----
    starLines   -   Configuration file (default: StarMe.cfg)

    Output
    -----
    synthresults.dat             -   Easy readable table with results from many linelists
    """

    try:  # Cleaning from previous runs
        os.remove('captain.log')
    except OSError:
        pass
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler('captain.log')
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # Check if models exist
    if not os.path.isdir('models'):
        logger.error('Model folder is missing.')

    # Create results directory
    if not os.path.isdir('results'):
        os.mkdir('results')
        logger.info('results directory was created')

    # Creating the output file
    _output(overwrite=overwrite, header=True)

    with open(starLines, 'r') as lines:
        for line in lines:
            if not line[0].isalpha():
                logger.debug('Skipping header: %s' % line.strip())
                continue
            logger.info('Line processing: %s' % line.strip())
            line = line.strip()
            line = line.split(' ')

            # Check if configuration parameters are correct
            if len(line) not in [1, 2, 7, 8]:
                logger.error('Could not process this information: %s' % line)
                continue

            # Check if the linelist is inside the directory, if not, log it
            if not os.path.isfile('rawLinelist/%s' % line[0]):
                logger.error('Error: rawLinelist/%s not found.' % line[0])
                raise IOError('Linelist file does not exist in the folder!')

            logger.info('Line list: %s' % line[0])
            # Create synthetic spectrum with solar values and the default options
            if len(line) == 1:
                initial = [5777, 4.44, 0.00, 1.00, 3.21, 1.90]
                options = _options()
                ranges, atomic_data = read_linelist(line[0], intname=options['inter_file'])
                x_obs, y_obs = (None, None)
                # Create initial synthetic model
                logger.info('Getting initial model grid')
                x_initial, y_initial = func(initial, atmtype=options['model'],
                ranges=ranges, driver='synth', version=options['MOOGv'], **options)
                logger.info('Interpolation successful.')
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))
                if options['save']:
                    save_synth_spec(x_initial, y_initial, y_obs=None, initial=initial, final=None, fname='initial.spec', **options)
                    logger.info('Save initial synthetic spectrum')
                print('Synthetic spectrum contains %s points' % len(x_initial))

            # Create synthetic spectrum with solar values and the options defined by the user
            elif len(line) == 2:
                options = _options(line[1])
                ranges, atomic_data = read_linelist(line[0], intname=options['inter_file'])
                if options['spt']:
                    logger.info('Spectral type given: %s' % line[1])
                    Teff, logg = _getSpt(options['spt'])
                    mic = getMic(Teff, logg)
                    mac = getMac(Teff, logg)
                    initial = (Teff, logg, 0.00, mic, mac, 1.90)
                else:
                    initial = [5777, 4.44, 0.00, 1.00, 3.21, 1.90]

                # Create initial synthetic model
                logger.info('Getting initial model grid')
                x_initial, y_initial = func(initial, atmtype=options['model'],
                ranges=ranges, driver='synth', version=options['MOOGv'], **options)
                logger.info('Interpolation successful.')
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))
                if options['save']:
                    save_synth_spec(x_initial, y_initial, y_obs=None, initial=initial, final=None, fname='initial.spec', **options)
                    logger.info('Save initial synthetic spectrum')
                print('Synthetic spectrum contains %s points' % len(x_initial))

                if options['observations']:
                    # Check if observations exit, if not pass another line
                    if os.path.isfile('spectra/%s' % options['observations']):
                        x_obs, y_obs, delta_l = read_obs_intervals('spectra/%s' % options['observations'], ranges, snr=options['snr'])
                        print('This is your observed spectrum: %s' % options['observations'])
                        print('Observed spectrum contains %s points' % len(x_obs))
                    elif os.path.isfile(options['observations']):
                        x_obs, y_obs, delta_l = read_obs_intervals(options['observations'], ranges, snr=options['snr'])
                        print('This is your observed spectrum: %s' % options['observations'])
                        print('Observed spectrum contains %s points' % len(x_obs))
                    else:
                        logger.error('Error: %s not found.' % options['observations'])
                        continue
                    # If observations exists, then why not, find best fit?
                    if options['minimize']:
                        print('Starting minimization...')
                        logger.info('Starting the minimization procedure...')
                        params, x_final, y_final = minimize_synth(initial, x_obs, y_obs, x_initial, y_initial, delta_l, ranges=ranges, **options)
                        logger.info('Minimization done.')
                        tmp = [line[0]] + [options['observations']] + params + [options['model'], options['resolution'], options['snr']]
                        _output(parameters=tmp)
                        logger.info('Saved results to: synthresults.dat')
                        if options['save']:
                            save_synth_spec(x_final, y_final, y_obs=y_obs, initial=initial, final=(params[0],params[2],params[4],params[6],params[8],params[10]), fname='final.spec', **options)
                            logger.info('Save final synthetic spectrum')

                else:
                    x_obs, y_obs = (None, None)
                    x_final, y_final = (None, None)

                if options['plot']:  # if there in no observed only the synthetic will be plotted
                    plot(x_obs, y_obs, x_initial, y_initial, res=options['plot_res'])
                    if options['minimize']:
                        # Plot also final spectra
                        plot(x_obs, y_obs, x_final, y_final, res=options['plot_res'])

            # Create spectra with parameters defined by the user but the rest are set to default
            elif len(line) == 7:
                logger.info('Initial parameters given by the user.')
                initial = map(float, line[1::])
                initial[0] = int(initial[0])
                options = _options()
                ranges, atomic_data = read_linelist(line[0], intname=options['inter_file'])
                x_obs, y_obs = (None, None)  # No observed spectrum, no plots

                # Create initial synthetic model
                logger.info('Getting initial model grid')
                x_initial, y_initial = func(initial, atmtype=options['model'],
                ranges=ranges, driver='synth', version=options['MOOGv'], **options)
                if options['save']:
                    logger.info('Save initial synthetic spectrum')
                    save_synth_spec(x_initial, y_initial, y_obs=None, initial=initial, final=None, fname='initial.spec', **options)
                logger.info('Interpolation successful.')
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))
                print('Synthetic spectrum contains %s points' % len(x_initial))

            # Create synthetic spectra with values set by the user and options altered.
            elif len(line) == 8:
                logger.info('Initial parameters given by user.')
                initial = map(float, line[1:-1])
                initial[0] = int(initial[0])
                logger.info('Initial parameters: {0}, {1}, {2}, {3}'.format(*initial))
                options = _options(line[-1])
                ranges, atomic_data = read_linelist(line[0], intname=options['inter_file'])
                logger.info('Getting initial model grid')
                logger.info('Setting solar values {0}, {1}, {2}, {3}'.format(*initial))

                if options['save']:
                    # Create initial synthetic model
                    x_initial, y_initial = func(initial, atmtype=options['model'],ranges=ranges, driver='synth', version=options['MOOGv'], **options)
                    print('Synthetic spectrum contains %s points' % len(x_initial))
                    logger.info('Save initial synthetic spectrum')
                    save_synth_spec(x_initial, y_initial, y_obs=None, initial=initial, final=None, fname='initial.spec', **options)
                if options['observations']:
                    print('This is your observed spectrum: %s' % options['observations'])
                    # Check if observations exit, if not pass another line
                    if os.path.isfile('spectra/%s' % options['observations']):
                        x_obs, y_obs, delta_l = read_obs_intervals('spectra/%s' % options['observations'], ranges, snr=options['snr'])
                        print('Observed spectrum contains %s points' % len(x_obs))
                        logger.info('Observed spectrum read.')
                    elif  os.path.isfile(options['observations']):
                        x_obs, y_obs, delta_l = read_obs_intervals(options['observations'], ranges, snr=options['snr'])
                        print('Observed spectrum contains %s points' % len(x_obs))
                        logger.info('Observed spectrum read.')
                    else:
                        logger.error('Error: %s not found.' % options['observations'])
                        continue
                    # If observations exists, then why not, find best fit?
                    if options['minimize']:
                        print('Starting minimization...')
                        logger.info('Starting the minimization procedure...')
                        params, x_obs_final, y_obs_final = minimize_synth(initial, x_obs, y_obs, delta_l, ranges=ranges, **options)
                        logger.info('Minimization done.')
                        if options['save']:
                            parameters = [params[0], params[2], params[4], params[6], params[8], params[10]]
                            x_final_synth, y_final_synth = func(parameter, atmtype=options['model'],ranges=ranges, driver='synth', version=options['MOOGv'], **options)
                            save_synth_spec(x_final_synth, y_final_synth, y_obs=None, initial=initial, final=(params[0],params[2],params[4],params[6],params[8],params[10]), fname='final.spec', **options)
                            logger.info('Save final synthetic spectrum')
                        tmp = [line[0]] + [options['observations']] + params + [options['model'], options['resolution'], options['snr']]
                        _output(parameters=tmp)
                        logger.info('Saved results to synthresults.dat')
                else:
                    x_obs, y_obs     = (None, None)
                    x_final, y_final = (None, None)

                if options['plot']:  # if there in no observed only the synthetic will be plotted
                    x_initial, y_initial = func(initial, atmtype=options['model'],ranges=ranges, driver='synth', version=options['MOOGv'], **options)
                    plot(x_obs, y_obs, x_initial, y_initial, res=options['plot_res'])
                    if options['minimize']:
                        # Plot also final spectra
                        parameters = [params[0], params[2], params[4], params[6], params[8], params[10]]
                        x_final_synth, y_final_synth = func(parameters, atmtype=options['model'],ranges=ranges, driver='synth', version=options['MOOGv'], **options)
                        plot(x_obs_final, y_obs_final, x_final_synth, y_final_synth, res=options['plot_res'])

            else:
                logger.error('Could not process information for this line list: %s' % line)
                continue

            if options['model'] != 'kurucz95' and options['model'] != 'apogee_kurucz' and options['model'] != 'marcs':
                logger.error('Your request for type: %s is not available' % options['model'])
                continue

            # Options not in use will be removed
            if __name__ == '__main__':
                options['GUI'] = False  # Running batch mode
            else:
                options['GUI'] = True  # Running GUI mode
            options.pop('spt')

    return

if __name__ == '__main__':
    synthdriver()
