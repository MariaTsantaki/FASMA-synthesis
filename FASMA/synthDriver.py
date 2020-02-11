#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os
from .utils import fun_moog_synth as func
from .observations import read_obs_intervals, plot, snr
from .minimization import MinimizeSynth, getMic, getMac
from .synthetic import read_linelist, read_linelist_elem, save_synth_spec
import time
import yaml

class FASMA:
    def __init__(self, cfgfile='config.yml', overwrite=None, **kwargs):
        '''The function that glues everything together. A log file is created
        with the list of processes (fasma.log).

        Input
        -----
        cfgfile : str
          Configuration file (default: config.yml)
        overwrite : bool
          Overwrite the synthresults.dat file (default: False)

        Output
        ------
        synthresults.dat : file
          Easy readable table with results
        '''

        self.cfgfile = cfgfile
        self.overwrite = overwrite
        self.status = None
        self.kwargs = kwargs
        # Setup of logger
        if os.path.isfile('fasma.log'):  # Cleaning from previous runs
            os.remove('fasma.log')
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        handler = logging.FileHandler('fasma.log')
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)
        if kwargs:
            self.configure(**kwargs)
        self.synthdriver()
        result = self.result()

    def configure(cls, cfgfile='config.yml', **kwargs):
        '''Create configuration file from kwargs.
        Otherwise set to default.
        '''
        defaults = {
            'linelist': 'linelist.lst',
            'teff': 5777,
            'logg': 4.44,
            'feh': 0.0,
            'vt': 1.0,
            'vmac': 3.21,
            'vsini': 1.90,
            'model': 'apogee_kurucz',
            'MOOGv': 2014,
            'save': False,
            'element': False,
            'fix_teff': False,
            'fix_logg': False,
            'fix_feh': False,
            'fix_vt': False,
            'fix_vmac': False,
            'fix_vsini': False,
            'plot': False,
            'plot_res': False,
            'damping': 1,
            'step_wave': 0.01,
            'step_flux': 3.0,
            'minimize': False,
            'refine': False,
            'observations': False,
            'intervals_file': 'intervals.lst',
            'snr': None,
            'resolution': None,
            'limb': 0.6
        }

        defaults.update(kwargs)
        dic = {}
        dic['star'] = defaults
        with open(cfgfile, 'w') as f:
            yaml.dump(dic, f)

    def _setup(self, line):
        '''Do the setup with initial parameters and options.

        Input
        -----
        line : list
          Each line from the configuration file after being split at spaces
          The format is spaced separated: linelist_file (params) (options)
        '''
        self.linelist = line['linelist']
        self._options(line)
        self.initial = [self.options['teff'], self.options['logg'], self.options['feh'], self.options['vt'], self.options['vmac'], self.options['vsini']]

    def _genStar(self):
        '''A generator for the configuration file.
        '''
        with open(self.cfgfile, 'r') as stream:
            try:
                data = yaml.safe_load(stream)
                self.logger.info('Reading input parameters.')
            except yaml.scanner.ScannerError as exc:
                self.logger.error('Could not process this information, check input parameters.')
                raise

        if data is None:
            print('Could not process this information, check input parameters.')
            self.logger.error('Could not process this information, check input parameters.')
            return None

        for item, line in data.items():
            self._setup(line)
            yield self.initial, self.options

    def _prepare(self):
        '''Check if linelist exists and create the first synthetic spectrum with
        initial parameters.
        '''

        if not os.path.isfile(self.linelist):
            message = 'The line list file {0} does not exists!\n'.format(self.options['linelist'])
            print(message)
            self.logger.error(message)
            raise StopIteration

        if not os.path.isfile(self.options['intervals_file']):
            message = 'The intervals list file {0} does not exists!\n'.format(self.options['intervals_file'])
            print(message)
            self.logger.error(message)
            raise StopIteration

        if self.options['element']:
            el = ['Na', 'Mg', 'Al', 'Si', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Ni']
            if self.options['element'] not in el:
                message = 'This element is not in the line list:', self.options['element']
                print(message)
                self.logger.error(message)
                raise StopIteration

        if self.options is None:
            self.logger.error('The inputs were not inserted correctly.\n')
            raise StopIteration

        if self.options['minimize'] and self.options['element']:
            message = 'I am confused with the inputs! First derive parameters and then abundances'
            self.logger.error('Options minimize and element at the same time!\n')
            print(message)
            raise StopIteration

        if self.options['element']:
            self.ranges, atomic = read_linelist_elem(
                self.linelist,
                element=self.options['element'],
                intname=self.options['intervals_file'],
            )
            self.logger.info('Getting initial model grid')
            self.xspec, self.yspec = func(
                self.initial,
                atmtype=self.options['model'],
                abund=self.initial[2],
                elem=self.options['element'],
                ranges=self.ranges,
                driver='synth',
                version=self.options['MOOGv'],
                **self.options
            )
        else:
            self.ranges, atomic = read_linelist(
                self.linelist, intname=self.options['intervals_file']
            )
            self.logger.info('Getting initial model grid')
            self.xspec, self.yspec = func(
                self.initial,
                atmtype=self.options['model'],
                ranges=self.ranges,
                driver='synth',
                version=self.options['MOOGv'],
                **self.options
            )

        if __name__ in ('__main__', 'synthDriver'):
            self.options['GUI'] = False  # Running batch mode
        else:
            self.options['GUI'] = True

    def _options(self, options=None):
        '''Reads the options inside the config file otherwise set to defaults.
        '''

        defaults = {
            'linelist': 'linelist.lst',
            'teff': 5777,
            'logg': 4.44,
            'feh': 0.0,
            'vt': 1.0,
            'vmac': 3.21,
            'vsini': 1.90,
            'model': 'apogee_kurucz',
            'MOOGv': 2014,
            'save': False,
            'element': False,
            'fix_teff': False,
            'fix_logg': False,
            'fix_feh': False,
            'fix_vt': False,
            'fix_vmac': False,
            'fix_vsini': False,
            'plot': False,
            'plot_res': False,
            'damping': 1,
            'step_wave': 0.01,
            'step_flux': 3.0,
            'minimize': False,
            'refine': False,
            'observations': False,
            'intervals_file': 'intervals.lst',
            'snr': None,
            'resolution': None,
            'limb': 0.6,
        }

        if not options:
            self.options = defaults
        else:
            defaults.update(options)
            defaults['model'] = defaults['model'].lower()
            defaults['step_wave'] = float(defaults['step_wave'])
            defaults['step_flux'] = float(defaults['step_flux'])
            defaults['limb'] = float(defaults['limb'])
            defaults['MOOGv'] = int(defaults['MOOGv'])
            if defaults['observations'] and (defaults['snr'] is None):
                if os.path.isfile(str(defaults['observations'])):
                    defaults['snr'] = snr(defaults['observations'])
                elif os.path.isfile('spectra/' + str(defaults['observations'])):
                    defaults['snr'] = snr('spectra/' + str(defaults['observations']))
                else:
                    print('The SNR was not measured.')
                    self.logger.error('Error: %s not found.' % defaults['observations'])
            if defaults['intervals_file']:
                defaults['intervals_file'] = str(defaults['intervals_file'])
            if defaults['resolution'] is not None:
                defaults['resolution'] = int(float(defaults['resolution']))
            self.options = defaults

    def _output(
        self, overwrite=None, header=None, stellarparams=False, abundance=False
    ):
        '''Create the output file 'synthresults.dat'

        Input
        -----
        overwrite - Overwrite the file
        header    - Only use True if this is for the file to be created

        Output
        -----
        The 'synthresults.dat' file with the final parameters.
        '''
        if abundance:
            if header:
                hdr2 = [
                    'linelist',
                    'observations',
                    'Teff',
                    'erTeff',
                    'logg',
                    'erlogg',
                    '[M/H]',
                    'er[M/H]',
                    'vmic',
                    'ervmic',
                    'vmac',
                    'ervmac',
                    'vsini',
                    'ervsini',
                    'Na',
                    'erNa',
                    'Mg',
                    'erMg',
                    'Al',
                    'erAl',
                    'Si',
                    'erSi',
                    'Ca',
                    'erCa',
                    'Sc',
                    'erSc',
                    'Ti',
                    'erTi',
                    'V',
                    'erV',
                    'Cr',
                    'erCr',
                    'Mn',
                    'erMn',
                    'Ni',
                    'erNi',
                    'chi^2',
                    'time',
                    'model',
                    'resolution',
                    'snr',
                ]
                if overwrite:
                    with open('synthresults_elements.dat', 'w') as output:
                        output.write('\t'.join(hdr2) + '\n')
                else:
                    if not os.path.isfile('synthresults_elements.dat'):
                        with open('synthresults_elements.dat', 'w') as output:
                            output.write('\t'.join(hdr2) + '\n')
            else:
                with open('synthresults_elements.dat', 'a') as output:
                    output.write('\t'.join(map(str, self.parameters)) + '\n')

        if stellarparams:
            if header:
                hdr1 = [
                    'linelist',
                    'observations',
                    'Teff',
                    'erTeff',
                    'logg',
                    'erlogg',
                    '[M/H]',
                    'er[M/H]',
                    'vmic',
                    'ervmic',
                    'vmac',
                    'ervmac',
                    'vsini',
                    'ervsini',
                    'chi^2',
                    'time',
                    'model',
                    'resolution',
                    'snr',
                ]

                if overwrite:
                    with open('synthresults.dat', 'w') as output:
                        output.write('\t'.join(hdr1) + '\n')
                else:
                    if not os.path.isfile('synthresults.dat'):
                        with open('synthresults.dat', 'w') as output:
                            output.write('\t'.join(hdr1) + '\n')
            else:
                with open('synthresults.dat', 'a') as output:
                    output.write('\t'.join(map(str, self.parameters)) + '\n')

    def minizationRunner(self, p=None):
        '''A function to run the minimization routine

        Output
        ------
        params : output parameters
        '''

        print('Starting minimization...')
        self.logger.info('Starting the minimization procedure...')
        start_time = time.time()
        # Run the minimization routine first time
        if p is not None:
            function = MinimizeSynth(
                p, self.xobs, self.yobs, self.ranges, **self.options
            )
            self.params, self.xo, self.yo = function.minimize()
        else:
            function = MinimizeSynth(
                self.initial, self.xobs, self.yobs, self.ranges, **self.options
            )
            self.params, self.xo, self.yo = function.minimize()

        if self.options['refine']:
            print('Refining the parameters...')
            print('Patience is the key...')
            params2 = [
                self.params[0],
                self.params[2],
                self.params[4],
                self.params[6],
                self.params[8],
                self.params[10],
            ]
            params2[3] = getMic(self.params[0], self.params[2], self.params[4])
            params2[4] = getMac(self.params[0], self.params[2])
            if round(params2[1], 2) == 5.0:
                params2[1] = 4.9
                self.options['fix_logg'] = True
            self.options['fix_vt'] = True
            self.options['fix_vmac'] = True
            function = MinimizeSynth(
                params2, self.xo, self.yo, self.ranges, **self.options
            )
            self.params, xxo, yyo = function.minimize()

        self.end_time = int(time.time() - start_time)
        print('Minimization finished in %s sec' % int(self.end_time))
        self.logger.info('Minimization done.')
        self.status = 1
        return self.status

    def minizationElementRunner(self, p=None):
        '''A function to run the minimization routine for element abundances.

        Output
        ------
        params : output parameters
        '''
        import numpy as np

        print('Starting minimization for element abundance...')
        self.logger.info('Starting minimization for element abundance...')
        start_time = time.time()
        # Run the minimization routine first time
        if p is not None:
            function = MinimizeSynth(
                p, self.xobs, self.yobs, self.ranges, **self.options
            )
            self.elemabund, self.xo, self.yo = function.minimizeElement()
        else:
            function = MinimizeSynth(
                self.initial, self.xobs, self.yobs, self.ranges, **self.options
            )
            self.elemabund, self.xo, self.yo = function.minimizeElement()

        self.end_time = int(time.time() - start_time)
        print('Minimization finished in %s sec' % int(self.end_time))
        self.logger.info('Minimization done.')

        # Make the output array to save results.
        el = ['Na', 'Mg', 'Al', 'Si', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Ni']
        a = np.empty(2 * len(el) + 1)
        a[:] = np.nan
        ind1 = el.index(self.options['element'])
        ind2 = ind1 + 1
        a[ind1] = self.elemabund[0]
        a[ind2] = self.elemabund[1]
        a[-1] = self.elemabund[2]
        self.abund = a
        self.status = 1
        return self.status

    def plotRunner(self, x=None, y=None, xs=None, ys=None, xf=None, yf=None, res=False):
        '''A function to plot spectra with or without residuals.

        Output
        ------
        plots
        '''

        if self.options['minimize']:
            # Create final spectrum with final parameters.
            p1 = [
                self.params[0],
                self.params[2],
                self.params[4],
                self.params[6],
                self.params[8],
                self.params[10],
            ]
            xf, yf = func(
                p1,
                atmtype=self.options['model'],
                ranges=self.ranges,
                driver='synth',
                version=self.options['MOOGv'],
                **self.options
            )
        elif self.options['element']:
            # Create final spectrum with final parameters.
            xf, yf = func(
                self.initial,
                atmtype=self.options['model'],
                abund=self.elemabund[0],
                elem=self.options['element'],
                ranges=self.ranges,
                driver='synth',
                version=self.options['MOOGv'],
                **self.options
            )
        else:
            xf, yf = (None, None)

        if self.options['plot_res']:
            plot(x, y, xs, ys, xf, yf, res=True)
        else:
            plot(x, y, xs, ys, xf, yf)

    def saveRunner(self):
        '''A function to save spectra in a fits like format in the results/ folder.
        If initial spectrum exists, save it.

        Output
        ------
        saved spectrum
        '''

        if not os.path.isdir('results'):
            os.mkdir('results')
            self.logger.info('results directory was created')

        if self.xspec is not None:
            save_synth_spec(
                self.xspec, self.yspec, initial=self.initial, **self.options
            )
        else:
            print('No spectrum was created.')

    def synthdriver(self):
        '''A function to connect all. This function applies the options set by
        the user in the config.yml file.
        '''

        # Creating the output file
        self._output(header=True, stellarparams=True, abundance=True)

        # Define options
        for (self.initial, self.options) in self._genStar():
            self.logger.info(
                'Initial parameters: {:.0f}, {:.2f}, {:.2f}, {:.2f}'.format(
                    *self.initial
                )
            )

            try:
                self._prepare()
            except StopIteration:
                continue

            if self.options['observations']:
                if os.path.isfile(self.options['observations']):
                    self.xobs, self.yobs, delta_l = read_obs_intervals(self.options['observations'],
                    self.ranges, snr=self.options['snr'])
                    print('Observed spectrum contains %s points' % len(self.xobs))
                    self.logger.info('Observed spectrum read.')
                elif os.path.isfile('spectra/' + self.options['observations']):
                    print(
                        'This is your observed spectrum: %s'
                        % self.options['observations']
                    )
                    self.xobs, self.yobs, delta_l = read_obs_intervals(
                        'spectra/' + self.options['observations'],
                        self.ranges,
                        snr=self.options['snr']
                    )
                    print('Observed spectrum contains %s points' % len(self.xobs))
                    self.logger.info('Observed spectrum read.')
                else:
                    print(
                        'The observed spectrum %s was not found.'
                        % self.options['observations']
                    )
                    self.logger.error(
                        'Error: %s not found.' % self.options['observations']
                    )
                    continue
            else:
                self.xobs, self.yobs = (None, None)

            if self.options['minimize']:
                if self.xobs is None:
                    continue
                self.status = self.minizationRunner()
                if self.status is None:
                    self.logger.error(
                        'The minimization routine did not finish succesfully.'
                    )
                    continue  # Problem with the minimization routine
                else:
                    self.logger.info('The minimization routine finished succesfully.')
                self.parameters = (
                    [self.linelist]
                    + [self.options['observations']]
                    + self.params
                    + [self.end_time]
                    + [
                        self.options['model'],
                        self.options['resolution'],
                        self.options['snr']
                    ]
                )
                # Save the results in the output file.
                self._output(stellarparams=True)
                self.xobs, self.yobs = self.xo, self.yo

            if self.options['element']:
                if self.xobs is None:
                    continue
                self.status = self.minizationElementRunner()
                if self.status is None:
                    self.logger.error(
                        'The minimization routine did not finish succesfully.'
                    )
                    continue  # Problem with the minimization routine
                else:
                    self.logger.info('The minimization routine finished succesfully.')
                self.parameters = (
                    [self.linelist]
                    + [self.options['observations']]
                    + self.initial
                    + self.abund.tolist()
                    + [self.end_time]
                    + [
                        self.options['model'],
                        self.options['resolution'],
                        self.options['snr'],
                    ]
                )
                self._output(abundance=True)
            else:
                self.abund = (None, None)

            if self.options['save']:
                self.logger.info('Save synthetic spectrum.')
                self.saveRunner()

            if self.options['plot']:
                self.logger.info('Plotting results.')
                self.plotRunner(x=self.xobs, y=self.yobs, xs=self.xspec, ys=self.yspec)

    def result(self):
        '''If any, get the output parameters.

        Input
        ------
        status : if None, no minimization happened

        Output
        ------
        params : output parameters in dictionary
        '''
        if self.status is None:
            self.logger.info('No parameters are returned.\n')
            result = None
        else:
            if self.options['element']:
                result = {
                "teff": self.initial[0],
                "logg": self.initial[1],
                "feh": self.initial[2],
                "vt": self.initial[3],
                "vmac": self.initial[4],
                "vsini": self.initial[5],
                "element": self.options['element'],
                "abund": self.abund[0],
                "erabund": self.abund[1],
                "spectrum": {
                "wave" : self.xobs,
                "flux" : self.yobs},
                }

            if self.options['minimize']:
                result = {
                "teff": self.params[0],
                "erteff": self.params[1],
                "logg": self.params[2],
                "erlogg": self.params[3],
                "feh": self.params[4],
                "erfeh": self.params[5],
                "vt": self.params[6],
                "ervt": self.params[7],
                "vmac": self.params[8],
                "ervmac": self.params[9],
                "vsini": self.params[10],
                "ervsini": self.params[11],
                "chi2": self.params[12],
                "spectrum": {
                "wave" : self.xobs,
                "flux" : self.yobs},
                }
        return result


if __name__ == '__main__':
    import sys

    if len(sys.argv) > 1:
        cfgfile = sys.argv[1]
    else:
        cfgfile = 'config.yml'

    driver = fasma(cfgfile=cfgfile, overwrite=None)
