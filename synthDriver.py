#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os
import yaml
from utils import fun_moog_synth as func
from observations import read_obs_intervals, plot, snr
from minimization import MinimizeSynth, getMic, getMac
from synthetic import read_linelist, save_synth_spec
import time

class synthMethod:

    def __init__(self, cfgfile='StarMe_synth.cfg', overwrite=None):
        """The function that glues everything together.

        Input
        -----
        cfgfile : str
          Configuration file (default: StarMe_synth.cfg)
        overwrite : bool
          Overwrite the synthresults.dat file (default: False)

        Output
        ------
        synthresults.dat : file
          Easy readable table with results
        """

        self.cfgfile = cfgfile
        self.overwrite = overwrite
        # Setup of logger
        if os.path.isfile('captain.log'):  # Cleaning from previous runs
            os.remove('captain.log')
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.DEBUG)
        handler = logging.FileHandler('captain.log')
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

    def _setup(self, line):
        """Do the setup with initial parameters and options.

        Input
        -----
        line : list
          Each line from the configuration file after being split at spaces
          The format is spaced separated: linelist file (params) (options)
        """

        self.linelist = line[0]
        if len(line) == 1:
            self.initial = [5777, 4.44, 0.00, 1.00, 3.21, 1.90]
            self._options()
        elif len(line) == 2:
            self.initial = [5777, 4.44, 0.00, 1.00, 3.21, 1.90]
            self._options(line[-1])
        elif len(line) == 7:
            self.initial = list(map(float, line[1::]))
            self.initial[0] = int(self.initial[0])
            self._options()
        elif len(line) == 8:
            self.initial = list(map(float, line[1:-1]))
            self.initial[0] = int(self.initial[0])
            self._options(line[-1])

    def _genStar(self):
        """A generator for the configuration file."""

        lines = open(self.cfgfile, 'r')
        for line in lines:
            if not line[0].isalpha():
                self.logger.debug('Skipping header: %s' % line.strip())
                continue
            self.logger.info('Line processing: %s' % line.strip())
            line = line.strip()
            line = line.split(' ')
            # Check if configuration parameters are correct
            if len(line) not in [1, 2, 7, 8]:
                self.logger.error('Could not process this information: %s' % line)
                continue

            self._setup(line)
            yield self.initial, self.options, line

    def _prepare(self):
        """Prepare the run with setup and first interpolation for model atmospheres."""

        if not os.path.isfile('rawLinelist/%s' % self.linelist):
            print('The line list does not exists!\n')
            self.logger.error('The line list does not exists!\n')
            return None

        self.ranges, atomic = read_linelist(self.linelist, intname=self.options['inter_file'])
        self.logger.info('Getting initial model grid')
        self.xspec, self.yspec = func(self.initial, atmtype=self.options['model'], ranges=self.ranges, driver='synth', version=self.options['MOOGv'], **self.options)

        if __name__ in ('__main__', 'synthDriver'):
            self.options['GUI'] = False  # Running batch mode
        else:
            self.options['GUI'] = True  # Running GUI mode

    def _options(self, options=None):
        '''Reads the options inside the config file'''

        defaults = {'model':        'kurucz95',
                    'MOOGv':        2014,
                    'save':         False,
                    'element':      False,
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
                    'inter_file':   'intervals_hr10_15n.lst',
                    'snr':          None,
                    'resolution':   None,
                    'limb':         0.6
                    }
        if not options:
            self.options = defaults
        else:
            for option in options.split(','):
                if ':' in option:
                    option = option.split(':')
                    defaults[option[0]] = option[1]
                else:
                    # Clever way to change the boolean
                    if option in ['teff', 'logg', 'mh', 'vt', 'vmac', 'vsini']:
                        option = 'fix_%s' % option
                    defaults[option] = False if defaults[option] else True
            defaults['model']     = defaults['model'].lower()
            defaults['step_wave'] = float(defaults['step_wave'])
            defaults['step_flux'] = float(defaults['step_flux'])
            defaults['limb']      = float(defaults['limb'])
            defaults['MOOGv']     = int(defaults['MOOGv'])
            if defaults['observations'] and (defaults['snr'] is None):
                if os.path.isfile(str(defaults['observations'])):
                    defaults['snr'] = snr(defaults['observations'])
                elif os.path.isfile('spectra/' + str(defaults['observations'])):
                    defaults['snr'] = snr('spectra/' + str(defaults['observations']))
                else:
                    self.logger.error('Error: %s not found.' % defaults['observations'])
            if defaults['inter_file']:
                defaults['inter_file'] = str(defaults['inter_file'])
            if defaults['resolution'] is not None:
                defaults['resolution'] = int(float(defaults['resolution']))
            self.options = defaults

    def _output(self, overwrite=None, header=None):
        """Create the output file 'synthresults.dat'

        Input
        -----
        overwrite - Overwrite the file
        header    - Only use True if this is for the file to be created
        """

        if header:
            hdr = ['linelist', 'observations', 'Teff', 'errTeff', 'logg', 'errlogg', '[M/H]', 'err[M/H]',
                'vmic', 'errvmic', 'vmac', 'errvmac', 'vsini', 'errvsini', 'chi^2', 'time', 'model', 'resolution', 'snr']
            if overwrite:
                with open('synthresults.dat', 'w') as output:
                    output.write('\t'.join(hdr)+'\n')
            else:
                if not os.path.isfile('synthresults.dat'):
                    with open('synthresults.dat', 'w') as output:
                        output.write('\t'.join(hdr)+'\n')
        else:
            with open('synthresults.dat', 'a') as output:
                output.write('\t'.join(map(str, self.parameters))+'\n')

    def minizationRunner(self, p=None):
        """A function to run the minimization routine

        Output
        ------
        params : output parameters
        """

        print('Starting minimization...')
        self.logger.info('Starting the minimization procedure...')
        start_time = time.time()
        # Run the minimization routine first time
        if p is not None:
            function = MinimizeSynth(p, self.xobs, self.yobs, self.ranges, **self.options)
            self.params, self.xo, self.yo = function.minimize()
        else:
            function = MinimizeSynth(self.initial, self.xobs, self.yobs, self.ranges, **self.options)
            self.params, self.xo, self.yo = function.minimize()

        if self.options['refine']:
            print('Refining the parameters...')
            print('Patience is the key...')
            params2 = [self.params[0], self.params[2], self.params[4], self.params[6], self.params[8], self.params[10]]
            params2[3] = getMic(self.params[0], self.params[2], self.params[4])
            params2[4] = getMac(self.params[0], self.params[2])
            if round(params2[1], 2) > 4.99:
                params2[1] = 4.9
                self.options['fix_logg'] = True
            self.options['fix_vt'] = True
            self.options['fix_vmac'] = True
            function = MinimizeSynth(params2, self.xo, self.yo, self.ranges, **self.options)
            self.params, xxo, yyo = function.minimize()

        self.logger.info('Minimization done.')
        self.end_time = int(time.time() - start_time)
        print('Minimization finished in %s sec' % int(self.end_time))
        status = 1
        return status

    def minizationElementRunner(self, p=None):
        """A function to run the minimization routine

        Output
        ------
        params : output parameters
        """

        print('Starting Li minimization...')
        self.logger.info('Starting the Li minimization procedure...')
        start_time = time.time()
        # Run the minimization routine first time
        if p is not None:
            function = MinimizeSynth(p, self.xobs, self.yobs, self.ranges, **self.options)
            self.elemabund, self.xo, self.yo = function.minimizeElement()
        else:
            function = MinimizeSynth(self.initial, self.xobs, self.yobs, self.ranges, **self.options)
            self.elemabund, self.xo, self.yo = function.minimizeElement()

        self.logger.info('Minimization done.')
        self.end_time = int(time.time() - start_time)
        print('Minimization finished in %s sec' % int(self.end_time))
        status = 1
        return status

    def plotRunner(self, x=None, y=None, xs=None, ys=None, xf=None, yf=None, res=False):
        """A function to plot spectra with or without residuals.
        """

        if self.options['minimize']:
            p1 = [self.params[0], self.params[2], self.params[4], self.params[6], self.params[8], self.params[10]]
            xf, yf = func(p1, atmtype=self.options['model'], ranges=self.ranges, driver='synth', version=self.options['MOOGv'], **self.options)
        else:
            xf, yf = (None, None)
        if self.options['plot_res']:
            plot(x, y, xs, ys, xf, yf, res=True)
        else:
            plot(x, y, xs, ys, xf, yf, res=False)

    def saveRunner(self):
        """A function to save spectra.
        """
        if not os.path.isdir('results'):
            os.mkdir('results')
            self.logger.info('results directory was created')

        if (self.xspec is not None) and (self.yspec is not None):
            save_synth_spec(self.xspec, self.yspec, initial=self.initial, **self.options)
        else:
            print('No spectrum was created.')

    def synthdriver(self):
        """A function to synthesize spectra.
        """

        # Creating the output file
        self._output(header=True)

        for (self.initial, self.options, line) in self._genStar():
            self.logger.info('Initial parameters: {:.0f}, {:.2f}, {:.2f}, {:.2f}'.format(*self.initial))
            self._prepare()
            if self.options is None:
                self.logger.error('The line list does not exists!\n')
                continue
            if self.options['minimize'] and self.options['element']:
                print('I am confused! First derive parameters and then abundances')
                self.logger.error('Minimization of parameters and abundances at the same time!\n')
                continue

            if self.options['save']:
                self.logger.info('Save synthetic spectrum.')
                self.saveRunner()

            if self.options['observations']:
                if os.path.isfile(self.options['observations']):
                    self.xobs, self.yobs, self.delta_l = read_obs_intervals(self.options['observations'], self.ranges, snr=self.options['snr'])
                    print('Observed spectrum contains %s points' % len(self.xobs))
                    self.logger.info('Observed spectrum read.')
                elif os.path.isfile('spectra/' + self.options['observations']):
                    print('This is your observed spectrum: %s' % self.options['observations'])
                    self.xobs, self.yobs, self.delta_l = read_obs_intervals('spectra/' + self.options['observations'], self.ranges, snr=self.options['snr'])
                    print('Observed spectrum contains %s points' % len(self.xobs))
                    self.logger.info('Observed spectrum read.')
                else:
                    print('The observed spectrum %s was not found.' % self.options['observations'])
                    self.logger.error('Error: %s not found.' % self.options['observations'])
                    continue
            else:
                self.xobs, self.yobs = (None, None)

            if self.options['minimize']:
                self.logger.info('Starting the minimization routine...')
                status = self.minizationRunner()
                if status is None:
                    self.logger.error('The minimization routine did not finish succesfully.')
                    continue  # Problem with the minimization routine
                else:
                    self.logger.info('The minimization routine finished succesfully.')
                self.parameters = [self.linelist] + [self.options['observations']] + self.params + [self.end_time] + [self.options['model'], self.options['resolution'], self.options['snr']]
                self._output()
                self.xobs, self.yobs = self.xo, self.yo

            if self.options['element']:
                self.logger.info('Starting the minimization routine...')
                status = self.minizationElementRunner()
                if status is None:
                    self.logger.error('The minimization routine did not finish succesfully.')
                    continue  # Problem with the minimization routine
                else:
                    self.logger.info('The minimization routine finished succesfully.')
                self.parameters = [self.linelist] + [self.options['observations']] + self.initial + self.elemabund + [self.end_time] + [self.options['model'], self.options['resolution'], self.options['snr']]
                self._output()

            if self.options['plot']:
                self.logger.info('Plotting results.')
                self.plotRunner(x=self.xobs, y=self.yobs, xs=self.xspec, ys=self.yspec)


if __name__ == '__main__':

    import sys
    if len(sys.argv) > 1:
        cfgfile = sys.argv[1]
    else:
        cfgfile = 'StarMe_synth.cfg'
    driver = synthMethod(cfgfile=cfgfile, overwrite=None)
    driver.synthdriver()
