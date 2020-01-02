#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from copy import copy
from .mpfit import mpfit


class MinimizeSynth:
    '''Minimize the chi square function between a synthetic spectrum to an observed.

    '''

    def __init__(self, p0, xobs, yobs, ranges, **kwargs):
        '''
        Input
        -----
        p0 : list
        Initial parameters (teff, logg, feh, vt, vmac, vsini)
        xobs : ndarray
        Observed wavelength
        yobs : ndarray
        Observed flux
        ranges : ndarray
        ranges of the intervals
        '''

        self.p0 = p0
        self.xobs = xobs
        self.yobs = yobs
        self.ranges = ranges
        self.model = kwargs['model']
        self.elem = kwargs['element']
        self.kwargs = kwargs
        self.y_obserr = 0.1  # arbitary value, the error of the flux
        # Set which values are fixed
        self.fix_teff = 1 if kwargs['fix_teff'] else 0
        self.fix_logg = 1 if kwargs['fix_logg'] else 0
        self.fix_feh = 1 if kwargs['fix_feh'] else 0
        self.fix_vt = 1 if kwargs['fix_vt'] else 0
        self.fix_vmac = 1 if kwargs['fix_vmac'] else 0
        self.fix_vsini = 1 if kwargs['fix_vsini'] else 0

    def parinfo_limit(self):
        '''Smart way to calculate the bounds of each of parameters depending on
        the grid of model atmospheres.
        '''

        if self.model.lower() == 'apogee_kurucz':
            bounds = [3500, 30000, 0.0, 5.0, -5.0, 1.5]
        if self.model.lower() == 'marcs':
            bounds = [3500, 7500, 0.0, 5.5, -4.0, 1.0]
        return bounds

    def bounds(self, i, p):
        '''
        Function to check if parameters during minimization are within bounds.
        Input
        -----
        p: parameter we want to check;
        i: the index of the bounds we want to check.
        '''

        if self.model.lower() == 'apogee_kurucz':
            bounds = [3500, 30000, 0.0, 5.0, -5.0, 1.5, 0.0, 9.99, 0.0, 20.0, 0.0, 99.9]
        if self.model.lower() == 'marcs':
            bounds = [3500, 7500, 0.0, 5.5, -4.0, 1.0, 0.0, 9.99, 0.0, 20.0, 0.0, 99.9]

        if p[int((i - 1) / 2)] < bounds[i - 1]:
            p[int((i - 1) / 2)] = bounds[i - 1]
        elif p[int((i - 1) / 2)] > bounds[i]:
            p[int((i - 1) / 2)] = bounds[i]
        return p

    def convergence_info(self, res):
        ''' Information on convergence from mpfit function. All values greater
        than zero can represent success (however status == 5 may indicate failure
        to converge).
        If the fit is unweighted (i.e. no errors were given, or the weights
        were uniformly set to unity), then .perror will probably not represent
        the true parameter uncertainties.
        *If* you can assume that the true reduced chi-squared value is unity --
        meaning that the fit is implicitly assumed to be of good quality --
        then the estimated parameter uncertainties can be computed by scaling
        .perror by the measured chi-squared value.

        Input
        -----
        res : ndarray
        the result of the mpfit function

        Output
        ------
        parameters : ndarray
        stellar parameters or abundances with their errors
        '''

        if res.status == -16:
            print(
                'status = %s : A parameter or function value has become infinite or an undefined number.'
                % res.status
            )
        if -15 <= res.status <= -1:
            print(
                'status = %s : MYFUNCT or iterfunct functions return to terminate the fitting process. '
                % res.status
            )
        if res.status == 0:
            print('status = %s : Improper input parameters.' % res.status)
        if res.status == 1:
            print(
                'status = %s : Both actual and predicted relative reductions in the sum of squares are at most ftol.'
                % res.status
            )
        if res.status == 2:
            print(
                'status = %s : Relative error between two consecutive iterates is at most xtol.'
                % res.status
            )
        if res.status == 3:
            print(
                'status = %s : Conditions for status = 1 and status = 2 both hold.'
                % res.status
            )
        if res.status == 4:
            print(
                'status = %s : The cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value.'
                % res.status
            )
        if res.status == 5:
            print(
                'status = %s : The maximum number of iterations has been reached.'
                % res.status
            )
        if res.status == 6:
            print('status = %s : ftol is too small.' % res.status)
        if res.status == 7:
            print('status = %s : xtol is too small.' % res.status)
        if res.status == 8:
            print('status = %s : gtol is too small.' % res.status)

        xreduced = round((res.fnorm / self.dof), 4)
        print('Iterations: %s' % res.niter)
        print('Value of the summed squared residuals: %s' % res.fnorm)
        print('Reduced chi squared: %s' % xreduced)
        print('Fitted parameters with uncertainties:')
        # scaled uncertainties
        pcerror = res.perror * np.sqrt(res.fnorm / self.dof)
        # This should separate elements from parameters
        if len(res.params) > 1:
            teff = round(float(res.params[0]), 0)
            logg = round(float(res.params[1]), 3)
            feh = round(float(res.params[2]), 3)
            vt = round(float(res.params[3]), 2)
            vmac = round(float(res.params[4]), 2)
            vsini = round(float(res.params[5]), 2)
            # scaled error
            erteff = round(float(pcerror[0]), 0)
            erlogg = round(float(pcerror[1]), 3)
            erfeh = round(float(pcerror[2]), 3)
            ervt = round(float(pcerror[3]), 2)
            ervmac = round(float(pcerror[4]), 2)
            ervsini = round(float(pcerror[5]), 2)
            # Save only the scaled error
            parameters = [
                teff,
                erteff,
                logg,
                erlogg,
                feh,
                erfeh,
                vt,
                ervt,
                vmac,
                ervmac,
                vsini,
                ervsini,
                xreduced,
            ]
            for i, x in enumerate(res.params):
                print(
                    "\t%s: %s +- %s (scaled error +- %s)"
                    % (
                        self.parinfo[i]['parname'],
                        round(x, 3),
                        round(res.perror[i], 3),
                        round(pcerror[i], 3),
                    )
                )
        else:
            abund = round(float(res.params), 3)
            erabund = round(float(pcerror), 3)
            # Save only the scaled error
            parameters = [abund, erabund, xreduced]
            for i, x in enumerate(res.params):
                print(
                    "\t%s: %s +- %s (scaled error +- %s)"
                    % (
                        self.parinfo[i]['parname'],
                        round(x, 3),
                        round(res.perror[i], 3),
                        round(pcerror[i], 3),
                    )
                )
        return parameters

    def myfunct(self, p, **kwargs):
        '''Function that returns the weighted deviates (to be minimized).

        Input
        ----
        p : list
          Parameters to be minimized

        Output
        -----
        (y-ymodel)/err : ndarray
          Model deviation from observation
        '''
        from .utils import fun_moog_synth as func
        from scipy.interpolate import InterpolatedUnivariateSpline

        options = kwargs['options']

        # minimize for element abundance
        if options['element']:
            xs, ys = func(
                self.p0,
                atmtype=self.model,
                abund=p,
                elem=self.elem,
                driver='synth',
                ranges=self.ranges,
                **options
            )
            print('    [' + self.elem + '/H]: {:1.2f}'.format(*p))

        # minimize for stellar parameters
        else:
            # Check if parameters are within bounds
            for i in range(1, 12, 2):
                p = self.bounds(i, p)
            print(
                '    Teff:{:8.1f}   logg: {:1.2f}   [Fe/H]: {:1.2f}   vt: {:1.2f}   vmac: {:1.2f}   vsini: {:1.2f}'.format(
                    *p
                )
            )
            xs, ys = func(
                p, atmtype=self.model, driver='synth', ranges=self.ranges, **options
            )

        sl = InterpolatedUnivariateSpline(xs, ys, k=1)
        self.ymodel = sl(self.xobs)
        # Check if interpolation is done correctly
        if np.isnan(self.ymodel).any():
            print('Warning: Check overlapping intervals.')
        # Error on the flux
        err = np.zeros(len(self.yobs)) + self.y_obserr
        status = 0
        # Print parameters at each function call
        return [status, (self.yobs - self.ymodel) / err]

    def exclude_bad_points(self):
        '''Function to exclude points from the spectrum, e.g. lines which are not
        in line list by comparing observations with the best synthetic model.
        '''

        delta_y = abs(np.subtract(self.ymodel, self.yobs) / self.ymodel)
        self.yobs_lpts = self.yobs[np.where((delta_y < 0.03) | (self.ymodel < 0.99))]
        self.xobs_lpts = self.xobs[np.where((delta_y < 0.03) | (self.ymodel < 0.99))]
        return self.xobs_lpts, self.yobs_lpts

    def minimize(self):
        '''Function to glue the minimization together. Set the minimization options here.
        '''

        # Set PARINFO structure for all free parameters for mpfit
        # The limits are also cheched by the bounds function
        teff_info = {
            'parname': 'Teff',
            'limited': [1, 1],
            'limits': [self.parinfo_limit()[0], self.parinfo_limit()[1]],
            'step': 100,
            'mpside': 2,
            'fixed': self.fix_teff,
        }
        logg_info = {
            'parname': 'logg',
            'limited': [1, 1],
            'limits': [self.parinfo_limit()[2], self.parinfo_limit()[3]],
            'step': 0.15,
            'mpside': 2,
            'fixed': self.fix_logg,
        }
        feh_info = {
            'parname': '[Fe/H]',
            'limited': [1, 1],
            'limits': [self.parinfo_limit()[4], self.parinfo_limit()[5]],
            'step': 0.10,
            'mpside': 2,
            'fixed': self.fix_feh,
        }
        vt_info = {
            'parname': 'vt',
            'limited': [1, 1],
            'limits': [0.0, 9.99],
            'step': 0.2,
            'mpside': 2,
            'fixed': self.fix_vt,
        }
        vmac_info = {
            'parname': 'vmac',
            'limited': [1, 1],
            'limits': [0.0, 20.0],
            'step': 0.5,
            'mpside': 2,
            'fixed': self.fix_vmac,
        }
        vsini_info = {
            'parname': 'vsini',
            'limited': [1, 1],
            'limits': [0.0, 99.0],
            'step': 1.5,
            'mpside': 2,
            'fixed': self.fix_vsini,
        }
        self.parinfo = [teff_info, logg_info, feh_info, vt_info, vmac_info, vsini_info]

        # A dictionary which contains the parameters to be passed to the user-supplied function specified by myfunct.
        # This is the way you can pass additional data to your user-supplied function without using global variables.
        self.fa = {
            'xobs': self.xobs,
            'ranges': self.ranges,
            'model': self.model,
            'yobs': self.yobs,
            'y_obserr': self.y_obserr,
            'options': self.kwargs,
        }
        m = mpfit(
            self.myfunct,
            xall=self.p0,
            parinfo=self.parinfo,
            ftol=1e-4,
            xtol=1e-4,
            gtol=1e-4,
            functkw=self.fa,
            maxiter=20,
        )
        self.dof = len(self.yobs) - len(m.params)
        self.parameters = self.convergence_info(m)

        # Prepare observations for the next iteration because of the refine option.
        if self.kwargs['refine']:
            self.xobs_lpts, self.yobs_lpts = self.exclude_bad_points()
        else:
            self.xobs_lpts, self.yobs_lpts = self.xobs, self.yobs
        return self.parameters, self.xobs_lpts, self.yobs_lpts

    def minimizeElement(self):

        # Initial value is the initial metallicity.
        pelem = self.p0[2]
        # Set PARINFO structure for all free parameters for mpfit, Limits to be added..
        elem_info = {'parname': self.elem, 'step': 0.10, 'mpside': 2}
        self.parinfo = [elem_info]

        # A dictionary which contains the parameters to be passed to the user-supplied function specified by myfunct.
        # This is the way you can pass additional data to your user-supplied function without using global variables.
        self.fa = {
            'xobs': self.xobs,
            'ranges': self.ranges,
            'model': self.model,
            'yobs': self.yobs,
            'y_obserr': self.y_obserr,
            'params': self.p0,
            'element': self.elem,
            'options': self.kwargs,
        }
        m = mpfit(
            self.myfunct,
            xall=[pelem],
            parinfo=self.parinfo,
            ftol=1e-4,
            xtol=1e-4,
            gtol=1e-4,
            functkw=self.fa,
            maxiter=20,
        )
        self.dof = len(self.yobs) - len(m.params)
        self.parameters = self.convergence_info(m)
        return self.parameters, self.xobs, self.yobs


def getMic(teff, logg, feh):
    '''Calculate microturbulence.
    '''
    if logg >= 3.80:  # Dwarfs (Tsantaki et al. 2013)
        mic = (6.932 * teff * (10 ** (-4))) - (0.348 * logg) - 1.437
        # mic = 1.163 + (7.808 * (10**(-4)) * (teff - 5800.0)) - (0.494*(logg - 4.30)) - (0.050*feh)
    elif logg < 3.80:  # Giants (Adibekyan et al. 2015)
        mic = 2.72 - (0.457 * logg) + (0.072 * feh)

    # Take care of negative values
    if mic < 0.1:
        mic = 0.1
    return round(mic, 2)


def getMac(teff, logg):
    '''Calculate macroturbulence.
    For hotter dwarfs: Doyle et al. 2014
    5200 < teff < 6400 K
    4.0 < logg < 4.6 dex
    For cooler dwarfs: Valenti et al. 2005
    For subgiants and giants (logg < 3.80 dex): Hekker & Melendez 2007
    '''

    if logg > 3.80:
        if teff > 5200:
            mac = (
                3.21
                + (2.33 * (teff - 5777.0) * (10 ** (-3)))
                + (2.00 * ((teff - 5777.0) ** 2) * (10 ** (-6)))
                - (2.00 * (logg - 4.44))
            )
        else:
            mac = 3.98 + ((teff - 5770.0) / 650.0)
    elif 2.90 <= logg <= 3.80:  # subgiants
        mac = -8.426 + (0.00241 * teff)
    elif 1.50 <= logg < 2.90:  # giants
        mac = -3.953 + (0.00195 * teff)
    elif logg < 1.50:  # bright giants
        mac = -0.214 + (0.00158 * teff)

    # Take care of negative values
    if mac < 0.30:
        mac = 0.30
    return round(mac, 2)
