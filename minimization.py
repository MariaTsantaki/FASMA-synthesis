#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from copy import copy
import time


class Minimize:
    '''Minimize for best parameters given a line list'''

    def __init__(self, x0, func, model, weights='null',
                 fix_teff=False, fix_logg=False, fix_feh=False, fix_vt=False,
                 iterations=160, EPcrit=0.001, RWcrit=0.003, ABdiffcrit=0.01,
                 MOOGv=2014, GUI=True, **kwargs):
        self.x0 = x0
        self.func = func
        self.model = model
        self.weights = weights
        self.fix_teff = fix_teff
        self.fix_logg = fix_logg
        self.fix_feh = fix_feh
        self.fix_vt = fix_vt
        self.maxiterations = iterations
        self.iteration = 0
        self.EPcrit = EPcrit
        self.RWcrit = RWcrit
        self.ABdiffcrit = ABdiffcrit
        self.MOOGv = MOOGv
        self.GUI = GUI
        if self.model.lower() == 'kurucz95':
            self.bounds = [3750, 39000, 0.0, 5.0, -3, 1, 0, 9.99]
        if self.model.lower() == 'apogee_kurucz':
            self.bounds = [3500, 30000, 0.0, 5.0, -5, 1.5, 0, 9.99]
        if self.model.lower() == 'marcs':
            self.bounds = [2500, 8000, 0.0, 5.0, -5, 1.0, 0, 9.99]

    def _getMic(self):
        '''Get the microturbulence if this is fixed'''
        if self.x0[1] >= 3.95:
            self.x0[3] = 6.932*self.x0[0]/10000 - 0.348*self.x0[1] - 1.437
        else:
            self.x0[3] = 2.72 - 0.457*self.x0[1] + 0.072*self.x0[2]

    def print_format(self):
        '''Print the stellar atmospheric parameters in a nice format'''
        rest = self.x0 + list((self.slopeEP, self.slopeRW, self.Abdiff))
        if self.iteration == 0:
            if self.GUI:
                print(' i     Teff       logg     [Fe/H]    vt    EPslope    RWslope    |FeI-FeII|')
                print('-' * 99)
            else:
                print(' i    Teff    logg    [Fe/H]    vt    EPslope    RWslope    |FeI-FeII|')
                print('-' * 70)
        else:
            print '{:4d}{:>6d}{:>8.2f}{:>+9.2f}{:>8.2f}{:>+9.3f}{:>+11.3f}{:>11.2f}'.format(self.iteration, *rest)

    def check_bounds(self, i):
        '''
        Function which checks if parameters are within bounds.
        Input - parameter: what we want to check; bounds: ze bounds;
        i: the index of the bounds we want to check'''
        if self.x0[int((i-1)/2)] < self.bounds[i-1]:
            self.x0[int((i-1)/2)] = self.bounds[i-1]
        elif self.x0[int((i-1)/2)] > self.bounds[i]:
            self.x0[int((i-1)/2)] = self.bounds[i]

    def check_convergence(self, fe_input):
        '''Check the convergence criteria'''
        self.slopeEP = 0.00 if self.fix_teff else self.slopeEP
        self.slopeRW = 0.00 if self.fix_vt else self.slopeRW
        self.Abdiff = 0.00 if self.fix_logg else self.Abdiff
        fe_input = self.x0[2]+7.47 if self.fix_feh else fe_input

        cond1 = abs(self.slopeRW) <= self.RWcrit
        cond2 = abs(self.Abdiff) <= self.ABdiffcrit
        cond3 = abs(self.slopeEP) <= self.EPcrit
        cond4 = round(fe_input, 2) == round(self.x0[2]+7.47, 2)
        return cond1 and cond2 and cond3 and cond4

    def _bump(self, alpha):
        '''Bump to the values in the list, x'''
        for i, X in enumerate(zip(alpha, self.x0)):
            ai, xi = X
            sig = 0.01 if ai*xi == 0 else ai*xi
            if ai:
                self.x0[i] = np.random.normal(xi, abs(sig))

    def _format_x0(self):
        '''Format the values in x0, so first value is an integer'''
        self.x0[0] = int(self.x0[0])
        self.x0[1] = round(self.x0[1], 2)
        self.x0[2] = round(self.x0[2], 2)
        self.x0[3] = round(self.x0[3], 2)

    def minimize(self):
        self._format_x0()
        res, self.slopeEP, self.slopeRW, abundances, self.x0 = self.func(self.x0, self.model, version=self.MOOGv)
        self.Abdiff = np.diff(abundances)[0]
        self.x0 = list(self.x0)

        if self.check_convergence(abundances[0]):
            return self.x0, True

        parameters = [copy(self.x0)]
        best = {}
        # Print the header before starting
        self.print_format()

        while self.iteration < self.maxiterations:
            # Step for Teff
            if (abs(self.slopeEP) >= self.EPcrit) and not self.fix_teff:
                self.x0[0] += 2000*self.slopeEP
                self.check_bounds(1)

            # Step for VT
            if (abs(self.slopeRW) >= self.RWcrit) and not self.fix_vt:
                self.x0[3] += 1.5*self.slopeRW
                self.check_bounds(7)

            # Step for logg
            if (abs(self.Abdiff) >= self.ABdiffcrit) and not self.fix_logg:
                self.x0[1] -= self.Abdiff
                self.check_bounds(3)

            # Step for [Fe/H]
            if not self.fix_feh:
                self.x0[2] = abundances[0]-7.47
                self.check_bounds(5)

            if self.fix_vt:
                self._getMic()  # Reset the microturbulence
                self.check_bounds(7)

            if self.x0 in parameters:
                alpha = [0] * 4
                alpha[0] = abs(self.slopeEP) if not self.fix_teff else 0
                alpha[1] = abs(self.Abdiff) if not self.fix_logg else 0
                alpha[2] = 0.01 if not self.fix_feh else 0
                alpha[3] = abs(self.slopeRW) if not self.fix_vt else 0
                self._bump(alpha)
                self.check_bounds(1)
                self.check_bounds(3)
                self.check_bounds(5)
                self.check_bounds(7)
            parameters.append(copy(self.x0))

            self._format_x0()
            res, self.slopeEP, self.slopeRW, abundances, self.x0 = self.func(self.x0, self.model, weights=self.weights, version=self.MOOGv)
            self.Abdiff = np.diff(abundances)[0]
            self.iteration += 1
            self.print_format()
            best[res] = parameters[-1]
            if self.check_convergence(abundances[0]):
                print '\nStopped in %i iterations' % self.iteration
                return self.x0, True

        print '\nStopped in %i iterations' % self.iteration
        if self.check_convergence(abundances[0]):
            return self.x0, True
        else:
            # Return the best solution rather than the last iteration
            _ = self.func(best[min(best.keys())], self.model, weights=self.weights, version=self.MOOGv)
            return best[min(best.keys())], False


class Minimize_synth:
    '''Minimize a synthetic spectrum to an observed

    Input
    -----
    p0 : list
      Initial parameters (teff, logg, feh, vt)
    x_obs : ndarray
      Observed wavelength
    y_obs : ndarray
      Observed flux
    r : ndarray
      ranges of the intervals
    fout : str
      Input line list files

    Output
    -----
    params : list
      Final parameters
    x_final : ndarray
      Final wavelength
    y_final : ndarray
      Final synthetic flux
    '''

    from utils import fun_moog_synth as func
    from mpfit import mpfit
    from scipy.interpolate import InterpolatedUnivariateSpline
    from synthetic import save_synth_spec

    def __init__(self, p0, x_obs, y_obs, r, fout, model='kurucz95',
                 fix_teff=None, fix_logg=None, fix_feh=None, fix_vt=None,
                 fix_vmac=None, fix_vsini=None, **kwargs):

        self.p0 = p0
        self.x_obs = x_obs
        self.y_obs = y_obs
        self.r = r
        self.fout = fout
        self.model = model
        self.fix_teff = 1 if fix_teff else 0
        self.fix_logg = 1 if fix_logg else 0
        self.fix_feh = 1 if fix_feh else 0
        self.fix_vt = 1 if fix_vt else 0
        self.fix_vmac = 1 if fix_vmac else 0
        self.fix_vsini = 1 if fix_vsini else 0

        # Setting up the bounds
        if self.model.lower() == 'kurucz95':
            self.bounds = [3750, 39000, 0.0, 5.0, -3, 1, 0, 9.99, 0, 50, 0, 100]
        if self.model.lower() == 'apogee_kurucz':
            self.bounds = [3500, 30000, 0.0, 5.0, -5, 1.5, 0, 9.99, 0, 50, 0, 100]
        if self.model.lower() == 'marcs':
            self.bounds = [2500, 8000, 0.0, 5.0, -5, 1.0, 0, 9.99, 0, 50, 0, 100]

        # Setting up PARINFO for mpfit
        teff_info  = {'limited': [1, 1], 'limits': self.bounds[0:2],   'step': 30,   'mpside': 2, 'fixed': fix_teff}
        feh_info   = {'limited': [1, 1], 'limits': self.bounds[4:6],   'step': 0.05, 'mpside': 2, 'fixed': fix_feh}
        logg_info  = {'limited': [1, 1], 'limits': self.bounds[2:4],   'step': 0.2,  'mpside': 2, 'fixed': fix_logg}
        vt_info    = {'limited': [1, 1], 'limits': self.bounds[6:8],   'step': 0.3,  'mpside': 2, 'fixed': fix_vt}
        vmac_info  = {'limited': [1, 1], 'limits': self.bounds[8:10],  'step': 0.5,  'mpside': 2, 'fixed': fix_vmac}
        vsini_info = {'limited': [1, 1], 'limits': self.bounds[10:12], 'step': 0.5,  'mpside': 2, 'fixed': fix_vsini}
        self.parinfo = [teff_info, logg_info, feh_info, vt_info, vmac_info, vsini_info]

        # Setting up keyword arguments for myfunct
        self.fa = {'x_obs': x_obs, 'r': r, 'fout': fout, 'model': model,
                   'y': y_obs, 'options': kwargs}

    def bounds(self, i, p):
        if p[int((i-1)/2)] < self.bounds[i-1]:
            p[int((i-1)/2)] = self.bounds[i-1]
        elif p[int((i-1)/2)] > self.bounds[i]:
            p[int((i-1)/2)] = self.bounds[i]
        return p

    def myfunct(self, p, y=None, **kwargs):
        '''Function that return the weighted deviates (to be minimized).

        Input
        ----
        p : list
          Parameters for the model atmosphere
        x_obs : ndarray
          Wavelength
        r : ndarray
          ranges of the intervals
        fout : str
          Line list file
        model : str
          Model atmosphere type
        y : ndarray
          Observed flux


        Output
        -----
        (y-ymodel)/err : ndarray
          Model deviation from observation
        '''

        options = kwargs['options']
        for i in range(1, 12, 2):
            p = self.bounds(i, p)

        x_s, y_s = func(p, atmtype=model, driver='synth', r=self.r, fout=self.fout, **options)
        sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
        flux_s = sl(self.x_obs)
        ymodel = flux_s
        # Error on the flux #needs corrections
        err = np.zeros(len(y)) + 0.01
        status = 0
        return([status, (y-ymodel)/err])

    def minimize(self):
        start_time = time.time()
        m = mpfit(self.myfunct, xall=self.p0, parinfo=self.parinfo,
                  ftol=1e-5, xtol=1e-5, gtol=1e-10, functkw=self.fa)
        end_time = time.time()-start_time
        print('status = %s' % m.status)
        print('Iterations: %s' % m.niter)
        print('Fitted pars:%s' % m.params)
        print('Uncertainties: %s' % m.perror)  # TODO: We can use them we define a realistic error on the flux
        print('Value of the summed squared residuals: %s' % m.fnorm)
        print('Number of calls to the function: %s' % m.nfev)
        print('Calculations finished in %s sec' % int(end_time))
        x_s, y_s = func(m.params, atmtype=model, driver='synth', r=r, fout=fout, **kwargs)
        sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
        flux_final = sl(x_obs)
        save_synth_spec(x_obs, flux_final, fname='final.spec')
        chi = ((x_obs - flux_final)**2)
        chi2 = np.sum(chi)
        print('This is your chi2 value: '), chi2

        return m.params, x_s, flux_final


def minimize_synth(p0, x_obs, y_obs, x_s, y_s, ranges, **kwargs):
    '''Minimize a synthetic spectrum to an observed

     Input
     -----
     p0 : list
       Initial parameters (teff, logg, feh, vt)
     x_obs : ndarray
       Observed wavelength
     y_obs : ndarray
       Observed flux
     x_s : ndarray
       Synthetic wavelength
     y_s : ndarray
       Synthetic flux
     ranges : ndarray
       ranges of the intervals
     atomic_data : ndarray
       Atomic data

     Output
     -----
     params : list
       Final parameters
     x_final : ndarray
       Final wavelength
     y_final : ndarray
       Final synthetic flux
    '''

    from utils import fun_moog_synth as func
    from mpfit import mpfit
    from scipy.interpolate import InterpolatedUnivariateSpline
    from synthetic import save_synth_spec


    def wave_step(delta_l, step_wave=0.01):
        '''Find the step of synthesis in wavelength depending the observations
        '''

        if delta_l < step_wave:
            step_wave = delta_l
        elif delta_l > step_wave:
            step_wave = delta_l
        else:
            step_wave
        return round(step_wave,3)


    def exclude_bad_points(x_obs, y_obs, x_s, y_s):
        '''Exclude points from the spectrum as continuum or bad points
        '''
        # Exclude some continuum points
        y_obs_lpts = y_obs[np.where(y_obs < 1.0)]
        x_obs_lpts = x_obs[np.where(y_obs < 1.0)]

        # Exclude some bad points
        sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
        ymodel = sl(x_obs_lpts)
        # Check if interpolation is done correctly
        if np.isnan(ymodel).any():
            print('Warning: Check overlapping intervals.')

        delta_y    = (np.subtract(ymodel,y_obs_lpts)/ymodel)
        y_obs_lpts = y_obs_lpts[np.where((delta_y < 0.03) | (ymodel<0.95))]
        x_obs_lpts = x_obs_lpts[np.where((delta_y < 0.03) | (ymodel<0.95))]
        return x_obs_lpts, y_obs_lpts


    def bounds(i, p, model):
        '''Smart way to calculate the bounds of each of parameters'''
        if model.lower() == 'kurucz95':
            bounds = [3750, 39000, 0.0, 5.0, -3, 1, 0, 9.99, 0, 50, 0, 100]
        if model.lower() == 'apogee_kurucz':
            bounds = [3500, 30000, 0.0, 5.0, -5, 1.5, 0, 9.99, 0, 50, 0, 100]
        if model.lower() == 'marcs':
            bounds = [2500, 8000, 0.0, 5.0, -5, 1.0, 0, 9.99, 0, 50, 0, 100]

        if p[int((i-1)/2)] < bounds[i-1]:
            p[int((i-1)/2)] = bounds[i-1]
        elif p[int((i-1)/2)] > bounds[i]:
            p[int((i-1)/2)] = bounds[i]
        return p


    def parinfo_limit(model):
        '''Smart way to calculate the bounds of each of parameters'''
        if model.lower() == 'kurucz95':
            bounds = [3750, 39000, 0.0, 5.0, -3, 1]
        if model.lower() == 'apogee_kurucz':
            bounds = [3500, 30000, 0.0, 5.0, -5, 1.5]
        if model.lower() == 'marcs':
            bounds = [2500, 8000, 0.0, 5.0, -5, 1.0]
        return bounds


    def _getMic(teff, logg, feh):
        """Calculate micro turbulence."""
        if logg >= 3.95:  # Dwarfs Tsantaki 2013
            mic = 6.932 * teff * (10**(-4)) - 0.348 * logg - 1.437
            # Take care of negative values
            if mic < 0:
                return 0.3
            return round(mic, 2)
        else:  # Giants Adibekyan 2015
            mic = 2.72 - (0.457 * logg) + (0.072 * feh)
            # Take care of negative values
            if mic < 0:
                return 0.3
            return round(mic, 2)


    def _getMac(teff, logg):
        """Calculate macro turbulence."""
        # For Dwarfs: Doyle et al. 2014
        # 5200 < teff < 6400
        # 4.0 < logg < 4.6
        if logg > 3.95:
            mac = 3.21 + (2.33 * (teff - 5777.) * (10**(-3)))
            + (2.00 * ((teff - 5777.)**2) * (10**(-6))) - (2.00 * (logg - 4.44))
        # For subgiants and giants: Hekker & Melendez 2007
        elif 2.0 <= logg <= 3.95:
            mac = -8.426 + (0.00241*teff)
        elif 1.0 <= logg < 2.0:
            mac = -3.953 + (0.00195*teff)
        elif logg < 1.0:
            mac = -0.214 + (0.00158*teff)
        
        # For negative values, keep a minimum of 0.3 km/s
        if mac < 0:
            mac = 0.30
        return round(mac, 2)


    def convergence_info(res, parinfo, dof):
        """
        Information on convergence. All values greater than zero can
        represent success (however status == 5 may indicate failure to
        converge).
        If the fit is unweighted (i.e. no errors were given, or the weights
        were uniformly set to unity), then .perror will probably not represent
        the true parameter uncertainties.
        *If* you can assume that the true reduced chi-squared value is unity --
        meaning that the fit is implicitly assumed to be of good quality --
        then the estimated parameter uncertainties can be computed by scaling
        .perror by the measured chi-squared value.
        """

        if res.status == -16:
            print('status = %s : A parameter or function value has become infinite or an undefined number.' % res.status)
        if -15 <= res.status <= -1:
            print('status = %s : MYFUNCT or iterfunct functions return to terminate the fitting process. ' % res.status)
        if res.status == 0:
            print('status = %s : Improper input parameters.' % res.status)
        if res.status == 1:
            print('status = %s : Both actual and predicted relative reductions in the sum of squares are at most ftol.' % res.status)
        if res.status == 2:
            print('status = %s : Relative error between two consecutive iterates is at most xtol.' % res.status)
        if res.status == 3:
            print('status = %s : Conditions for status = 1 and status = 2 both hold.' % res.status)
        if res.status == 4:
            print('status = %s : The cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value.' % res.status)
        if res.status == 5:
            print('status = %s : The maximum number of iterations has been reached.' % res.status)
        if res.status == 6:
            print('status = %s : ftol is too small.' % res.status)
        if res.status == 7:
            print('status = %s : xtol is too small.' % res.status)
        if res.status == 8:
            print('status = %s : gtol is too small.' % res.status)

        print('Iterations: %s' % res.niter)
        print('Fitted parameters with uncertainties:')
        # scaled uncertainties
        pcerror = res.perror * np.sqrt(res.fnorm / dof)
        teff = round(float(res.params[0]),0)
        logg = round(float(res.params[1]),3)
        feh = round(float(res.params[2]),3)
        vt = round(float(res.params[3]),2)
        vmac = round(float(res.params[4]),2)
        vsini = round(float(res.params[5]),1)

        erteff = round(float(pcerror[0]),0)
        erlogg = round(float(pcerror[1]),3)
        erfeh = round(float(pcerror[2]),3)
        ervt = round(float(pcerror[3]),2)
        ervmac = round(float(pcerror[4]),2)
        ervsini = round(float(pcerror[5]),1)
        parameters = [teff, erteff, logg, erlogg, feh, erfeh, vt, ervt, vmac, ervmac, vsini, ervsini, res.status]
        for i, x in enumerate(res.params):
                    print( "\t%s: %s +- %s (scaled error +- %s)" % (parinfo[i]['parname'], round(x, 3), round(res.perror[i], 3), round(pcerror[i], 3)))
        print('Value of the summed squared residuals: %s' % res.fnorm)
        return parameters


    def myfunct(p, x_obs=None, ranges=None, model=None, y=None, y_obserr=0.01, **kwargs):
        '''Function that return the weighted deviates (to be minimized).

        Input
        ----
        p : list
          Parameters for the model atmosphere
        x_obs : ndarray
          Wavelength
        ranges : ndarray
          ranges of the intervals
        atomic_data : ndarray
          Atomic data
        model : str
          Model atmosphere type
        y : ndarray
          Observed flux


        Output
        -----
        (y-ymodel)/err : ndarray
          Model deviation from observation
        '''

        # Definition of the Model spectrum to be iterated
        options = kwargs['options']
        # Check for bounds
        p = bounds(1, p, model)
        p = bounds(3, p, model)
        p = bounds(5, p, model)
        p = bounds(7, p, model)
        p = bounds(9, p, model)
        p = bounds(11, p, model)

        if options['fix_vt'] and options['flag_vt']:
            p[3] = _getMic(p[0], p[1], p[2])
        if options['fix_vmac'] and options['flag_vmac']:
            p[4] = _getMac(p[0], p[1])
        x_s, y_s = func(p, atmtype=model, driver='synth', ranges=ranges, **options)
        sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
        ymodel = sl(x_obs)
        # Check if interpolation is done correctly
        if np.isnan(ymodel).any():
            print('Warning: Check overlapping intervals.')
        # Error on the flux #needs corrections
        err = np.zeros(len(y)) + y_obserr
        status = 0
        #Print parameters at each function call
        print('   Teff:{:8.1f}   logg: {:1.2f}   [Fe/H]: {:1.2f}   vt: {:1.2f}   vmac: {:1.2f}   vsini: {:1.2f}'.format(*p))
        return([status, (y-ymodel)/err])


    #Define step for synthesis according to observations
    delta_l = x_obs[1] - x_obs[0]
    kwargs['step_wave'] = wave_step(delta_l)
    # Define the observation points
    x_o, y_o = exclude_bad_points(x_obs, y_obs, x_s, y_s)

    model = kwargs['model']
    y_obserr = 1.0/(kwargs['snr']) #Gaussian noise
    fix_teff = 1 if kwargs['fix_teff'] else 0
    fix_logg = 1 if kwargs['fix_logg'] else 0
    fix_feh = 1 if kwargs['fix_feh'] else 0
    fix_vt = 1 if kwargs['fix_vt'] else 0
    fix_vmac = 1 if kwargs['fix_vmac'] else 0
    fix_vsini = 1 if kwargs['fix_vsini'] else 0

    # Set PARINFO structure for all 6 free parameters for mpfit
    # Teff, logg, feh, vt, vmac, vsini
    # The limits are also cheched by the bounds function
    teff_info  = {'parname':'Teff',   'limited': [1, 1], 'limits': [parinfo_limit(model)[0], parinfo_limit(model)[1]], 'step': 100,  'mpside': 2, 'fixed': fix_teff}
    logg_info  = {'parname':'logg',   'limited': [1, 1], 'limits': [parinfo_limit(model)[2], parinfo_limit(model)[3]], 'step': 0.1,  'mpside': 2, 'fixed': fix_logg}
    feh_info   = {'parname':'[Fe/H]', 'limited': [1, 1], 'limits': [parinfo_limit(model)[4], parinfo_limit(model)[5]], 'step': 0.05, 'mpside': 2, 'fixed': fix_feh}
    vt_info    = {'parname':'vt',     'limited': [1, 1], 'limits': [0.0, 9.99],  'step': 0.5,  'mpside': 2, 'fixed': fix_vt}
    vmac_info  = {'parname':'vmac',   'limited': [1, 1], 'limits': [0.0, 50.0],  'step': 2.0,  'mpside': 2, 'fixed': fix_vmac}
    vsini_info = {'parname':'vsini',  'limited': [1, 1], 'limits': [0.3, 100.0], 'step': 2.0,  'mpside': 2, 'fixed': fix_vsini}

    parinfo = [teff_info, logg_info, feh_info, vt_info, vmac_info, vsini_info]

    # A dictionary which contains the parameters to be passed to the
    # user-supplied function specified by myfunct via the standard Python
    # keyword dictionary mechanism. This is the way you can pass additional
    # data to your user-supplied function without using global variables.
    fa = {'x_obs': x_o, 'ranges': ranges, 'model': model, 'y': y_o, 'y_obserr': y_obserr, 'options': kwargs}

    # Minimization starts here
    # Measure time
    start_time = time.time()
    m = mpfit(myfunct, xall=p0, parinfo=parinfo, ftol=1e-5, xtol=1e-5, gtol=1e-4, functkw=fa)
    #Print results
    dof = len(y_o)-len(m.params)
    if kwargs['refine']:
        print('Refining the parameters...')
        kwargs['flag_vt'] = True
        kwargs['flag_vmac'] = True
        f = mpfit(myfunct, xall=m.params, parinfo=parinfo, ftol=1e-5, xtol=1e-5, gtol=1e-4, functkw=fa)
        parameters = convergence_info(f, parinfo, dof)
        end_time = time.time()-start_time
        print('Minimization finished in %s sec' % int(end_time))
        #Final synthetic spectrum
        x_s, y_s = func(f.params, atmtype=model, driver='synth', ranges=ranges, **kwargs)
        sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
        flux_final = sl(x_o)
    else:
        parameters = convergence_info(m, parinfo, dof)
        end_time = time.time()-start_time
        print('Minimization finished in %s sec' % int(end_time))
        #Final synthetic spectrum
        x_s, y_s = func(m.params, atmtype=model, driver='synth', ranges=ranges, **kwargs)
        sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
        flux_final = sl(x_o)

    x_init, y_init = func(p0, atmtype=model, driver='synth', ranges=ranges, **kwargs)
    sl = InterpolatedUnivariateSpline(x_init, y_init, k=1)
    flux_initial = sl(x_o)

    err = np.zeros(len(y_o)) + y_obserr
    chi = ((y_o - flux_final)**2/(err**2))
    chi2 = np.sum(chi)/dof
    print('This is your reduced chi2 value: '), round(chi2,2)
    for i, r in enumerate(ranges):
        wave   = x_o[np.where((x_o >= float(r[0])) & (x_o <= float(r[1])))]
        fm     = flux_final[np.where((x_o >= float(r[0])) & (x_o <= float(r[1])))]
        fminit = flux_initial[np.where((x_o >= float(r[0])) & (x_o <= float(r[1])))]
        fobs   = y_o[np.where((x_o >= float(r[0])) & (x_o <= float(r[1])))]

        err = np.zeros(len(fobs)) + y_obserr
        chi = ((fobs - fm)**2/(err**2))
        dof = len(fobs)-len(m.params)
        chi2final = np.sum(chi)/dof

        chiinit = ((fobs - fminit)**2/(err**2))
        chi2init = np.sum(chiinit)/dof
        print('%s This is your reduced chi2 value: initial: %s final: %s') % (i, round(chi2init,2), round(chi2final,2))

    parameters = parameters + [round(chi2,2)] + [int(end_time)]
    return parameters, x_o, flux_final


def mcmc_synth(x0, observed, limits):
    '''This could be cool if it worked'''

    import emcee
    from utils import interpol_synthetic, fun_moog as func

    def lnlike(theta, x, y, yerr):
        teff, logg = theta
        x0 = (teff, logg, 0.0, 1.0)
        func(x0, driver='synth')
        _, _, model = interpol_synthetic(x, y, limits[0], limits[1])
        inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2))
        return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

    def lnprior(theta):
        teff, logg = theta
        if 5500 < teff < 6000 and 4.0 < logg < 4.9:
            return 0.0
        return -np.inf

    def lnprob(theta, x, y, yerr):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(theta, x, y, yerr)

    x, y = np.loadtxt(observed, unpack=True, usecols=(0, 1))
    idx = (x >= limits[0]) & (x <= limits[1])
    x, y = x[idx], y[idx]
    y /= np.median(y)
    # Normalization (use first 50 points below 1.2 as constant continuum)
    maxes = y[(y < 1.2)].argsort()[-50:][::-1]
    y /= np.median(y[maxes])

    x0 = np.array(x0)

    ndim, nwalkers = 2, 8
    Teff_step, logg_step = 50, 0.1
    pos = [x0 + np.random.randn()*np.array([Teff_step, logg_step]) for i in range(nwalkers)]

    print('starting MCMC')
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, np.array([0.01]*len(x))))
    print('still doing MCMC I guess')
    sampler.run_mcmc(pos, 500)

    print('Are we done soon???')
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    print('Done!')
    # print(samples)
    print(samples.shape)
    import corner
    import matplotlib.pyplot as plt
    fig = corner.corner(samples, labels=["$Teff$", "$logg$"], truths=[5777, 4.44])
    plt.show()
