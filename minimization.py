#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import time


def getMic(teff, logg, feh):
    """Calculate micro turbulence.

    Inputs
    ------
    teff : float
      Effective temperature
    logg : float
      Surface gravity
    feh : float
      Metallicity, [Fe/H]

    Output
    ------
    mic : float
      Microturbulence
    """
    if logg >= 3.50:  # Dwarfs Tsantaki 2013
        mic = (6.932 * teff * (10**(-4))) - (0.348 * logg) - 1.437
    elif logg < 3.50:  # Giants Adibekyan 2015
        mic = 2.72 - (0.457 * logg) + (0.072 * feh)

    # Take care of negative values
    if mic < 0.1:
        mic = 0.1
    return round(mic, 2)


def getMac(teff, logg):
    """Calculate macro turbulence.

    Inputs
    ------
    teff : float
      Effective temperature
    logg : float
      Surface gravity

    Output
    ------
    mac : float
      Macroturbulence"""
    # For Dwarfs: Doyle et al. 2014
    # 5200 < teff < 6400, 4.0 < logg < 4.6
    if logg > 3.9:
        mac = 3.21 + (2.33 * (teff - 5777.) * (10**(-3))) + (2.00 * ((teff - 5777.)**2) * (10**(-6))) - (2.00 * (logg - 4.44))
    # For subgiants and giants: Hekker & Melendez 2007
    elif 2.9 <= logg <= 3.9: # subgiants
        mac = -8.426 + (0.00241*teff)
    elif 1.5 <= logg < 2.9: # giants
        mac = -3.953 + (0.00195*teff)
    elif logg < 1.5: # bright giants
        mac = -0.214 + (0.00158*teff)

    # For negative values, keep a minimum of 0.3 km/s
    if mac < 0.10:
        mac = 0.10
    return round(mac, 2)


def minimize_synth(p0, x_obs, y_obs, delta_l, ranges, **kwargs):
    '''Minimize a synthetic spectrum to an observed

     Input
     -----
     p0 : list
       Initial parameters (teff, logg, feh, vt)
     x_obs : np.ndarray
       Observed wavelength
     y_obs : np.ndarray
       Observed flux
     ranges : np.ndarray
       ranges of the intervals

     Output
     -----
     params : list
       Final parameters
     x_final : np.ndarray
       Final wavelength
     y_final : np.ndarray
       Final synthetic flux
    '''

    from utils import fun_moog_synth as func
    from mpfit import mpfit
    from scipy.interpolate import InterpolatedUnivariateSpline
    from synthetic import save_synth_spec


    def exclude_bad_points(x_obs, y_obs, x_s, y_s):
        '''Exclude points from the spectrum as continuum or bad points
        '''

        # Exclude some bad points
        sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
        ymodel = sl(x_obs)
        # Check if interpolation is done correctly
        if np.isnan(ymodel).any():
            print('Warning: Check overlapping intervals.')
        delta_y    = abs(np.subtract(ymodel,y_obs)/ymodel)
        y_obs_lpts = y_obs[np.where((delta_y < 0.03) | (ymodel < 0.98))]
        x_obs_lpts = x_obs[np.where((delta_y < 0.03) | (ymodel < 0.98))]
        return x_obs_lpts, y_obs_lpts


    def bounds(i, p, model):
        '''Smart way to calculate the bounds of each of parameters'''
        if model.lower() == 'apogee_kurucz':
            bounds = [3500, 30000, 0.0, 5.0, -5.0, 1.5, 0.0, 9.99, 0.0, 20.0, 0.0, 99.9]
        if model.lower() == 'marcs':
            bounds = [2500, 8000, 0.0, 5.0, -5.0, 1.0, 0.0, 9.99, 0.0, 20.0, 0.0, 99.9]

        if p[int((i-1)/2)] < bounds[i-1]:
            p[int((i-1)/2)] = bounds[i-1]
        elif p[int((i-1)/2)] > bounds[i]:
            p[int((i-1)/2)] = bounds[i]
        return p


    def parinfo_limit(model):
        '''Smart way to calculate the bounds of each of parameters'''
        if model.lower() == 'apogee_kurucz':
            bounds = [3500, 30000, 0.0, 5.0, -5.0, 1.5]
        if model.lower() == 'marcs':
            bounds = [2500, 8000, 0.0, 5.0, -5.0, 1.0]
        return bounds


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

        x_red = round((res.fnorm / dof),4)
        print('Iterations: %s' % res.niter)
        print('Value of the summed squared residuals: %s' % res.fnorm)
        print('Reduced chi squared: %s' % x_red)
        print('Fitted parameters with uncertainties:')
        # scaled uncertainties
        pcerror = res.perror * np.sqrt(res.fnorm / dof)
        teff = round(float(res.params[0]),0)
        logg = round(float(res.params[1]),3)
        feh = round(float(res.params[2]),3)
        vt = round(float(res.params[3]),2)
        vmac = round(float(res.params[4]),2)
        vsini = round(float(res.params[5]),1)
        #scaled error
        erteff = round(float(pcerror[0]),0)
        erlogg = round(float(pcerror[1]),3)
        erfeh = round(float(pcerror[2]),3)
        ervt = round(float(pcerror[3]),2)
        ervmac = round(float(pcerror[4]),2)
        ervsini = round(float(pcerror[5]),1)
        # Save only the scaled error
        parameters = [teff, erteff, logg, erlogg, feh, erfeh, vt, ervt, vmac, ervmac, vsini, ervsini, x_red, res.status]
        for i, x in enumerate(res.params):
                    print( "\t%s: %s +- %s (scaled error)" % (parinfo[i]['parname'], round(x, 3), round(pcerror[i], 3)))
                    #print( "\t%s: %s +- %s (scaled error +- %s)" % (parinfo[i]['parname'], round(x, 3), round(res.perror[i], 3), round(pcerror[i], 3)))
        return parameters


    def myfunct(p, x_obs=None, ranges=None, model=None, y=None, y_obserr=0.1, **kwargs):
        '''Function that return the weighted deviates (to be minimized).

        Input
        ----
        p : list
          Parameters for the model atmosphere
        x_obs : ndarray
          Wavelength
        ranges : ndarray
          ranges of the intervals
        model : str
          Model atmosphere type
        y : ndarray
          Observed flux

        Output
        -----
        (y-ymodel)/err : ndarray
          Model deviation from observations
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
            p[3] = getMic(p[0], p[1], p[2])
        if options['fix_vmac'] and options['flag_vmac']:
            p[4] = getMac(p[0], p[1])

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
        print('    Teff:{:8.1f}   logg: {:1.2f}   [Fe/H]: {:1.2f}   vt: {:1.2f}   vmac: {:1.2f}   vsini: {:1.2f}'.format(*p))
        return([status, (y-ymodel)/err])


    def error_synth(p, **kwargs):
        '''Function to calculate the errors.'''

        teff_up    = p[0] + 100.0
        teff_down  = p[0] - 100.0
        logg_up    = p[1] + 0.1
        logg_down  = p[1] - 0.1
        feh_up     = p[2] + 0.05
        feh_down   = p[2] - 0.05
        vsini_up   = p[5] + 0.5
        vsini_down = p[5] - 0.5
        #teff
        vt_info['fixed']    = 1
        vmac_info['fixed']  = 1
        teff_info['fixed']  = 1
        logg_info['fixed']  = 0
        feh_info['fixed']   = 0
        vsini_info['fixed'] = 0
        print('Calculation of the errors')
        print('    Fix Teff +100K')
        teff_up_fit = mpfit(myfunct, xall=(teff_up, p[1], p[2], p[3], p[4], p[5]), parinfo=parinfo, ftol=1e-2, xtol=1e-2, gtol=1e-2, functkw=fa, maxiter=2, nocovar=1, quiet=1, iterfunct=None)
        print('    Fix Teff -100K')
        teff_down_fit = mpfit(myfunct, xall=(teff_down, p[1], p[2], p[3], p[4], p[5]), parinfo=parinfo, ftol=1e-2, xtol=1e-2, gtol=1e-2, functkw=fa, maxiter=2, nocovar=1, quiet=1, iterfunct=None)
        #logg
        teff_info['fixed']  = 0
        logg_info['fixed']  = 1
        print('    Fix logg +0.1dex')
        logg_up_fit = mpfit(myfunct, xall=(p[0], logg_up, p[2], p[3], p[4], p[5]), parinfo=parinfo, ftol=1e-2, xtol=1e-2, gtol=1e-2, functkw=fa, maxiter=2, nocovar=1, quiet=1, iterfunct=None)
        print('    Fix logg -0.1dex')
        logg_down_fit = mpfit(myfunct, xall=(p[0], logg_down, p[2], p[3], p[4], p[5]), parinfo=parinfo, ftol=1e-2, xtol=1e-2, gtol=1e-2, functkw=fa, maxiter=2, nocovar=1, quiet=1, iterfunct=None)
        #feh
        logg_info['fixed']  = 0
        feh_info['fixed']   = 1
        print('    Fix [Fe/H] +0.05dex')
        feh_up_fit = mpfit(myfunct, xall=(p[0], p[1], feh_up, p[3], p[4], p[5]), parinfo=parinfo, ftol=1e-2, xtol=1e-2, gtol=1e-2, functkw=fa, maxiter=2, nocovar=1, quiet=1, iterfunct=None)
        print('    Fix [Fe/H] -0.05dex')
        feh_down_fit = mpfit(myfunct, xall=(p[0], p[1], feh_down, p[3], p[4], p[5]), parinfo=parinfo, ftol=1e-2, xtol=1e-2, gtol=1e-2, functkw=fa, maxiter=2, nocovar=1, quiet=1, iterfunct=None)
        #vsini
        feh_info['fixed']   = 0
        vsini_info['fixed'] = 1
        print('    Fix vsini +0.5km/s')
        vsini_up_fit = mpfit(myfunct, xall=(p[0], p[1], p[2], p[3], p[4], vsini_up), parinfo=parinfo, ftol=1e-2, xtol=1e-2, gtol=1e-2, functkw=fa, maxiter=2, nocovar=1, quiet=1, iterfunct=None)
        print('    Fix vsini -0.5km/s')
        vsini_down_fit = mpfit(myfunct, xall=(p[0], p[1], p[2], p[3], p[4], vsini_down), parinfo=parinfo, ftol=1e-2, xtol=1e-2, gtol=1e-2, functkw=fa, maxiter=2, nocovar=1, quiet=1, iterfunct=None)

        teff = [logg_up_fit.params[0], logg_down_fit.params[0], feh_up_fit.params[0], feh_down_fit.params[0], vsini_up_fit.params[0], vsini_down_fit.params[0]]
        logg = [teff_up_fit.params[1], teff_down_fit.params[1], feh_up_fit.params[1], feh_down_fit.params[1], vsini_up_fit.params[1], vsini_down_fit.params[1]]
        feh = [teff_up_fit.params[2], teff_down_fit.params[2], logg_up_fit.params[2], logg_down_fit.params[2], vsini_up_fit.params[2], vsini_down_fit.params[2]]
        vsini = [teff_up_fit.params[5], teff_down_fit.params[5], logg_up_fit.params[5], logg_down_fit.params[5], feh_up_fit.params[5], feh_down_fit.params[5]]

        # errors
        error_teff = (teff - p[0])**2
        error_teff_std = np.sqrt(np.sum(error_teff)/len(teff))

        error_logg = (logg - p[1])**2
        error_logg_std = np.sqrt(np.sum(error_logg)/len(logg))

        error_feh = (feh - p[2])**2
        error_feh_std = np.sqrt(np.sum(error_feh)/len(feh))

        error_vsini = (vsini - p[5])**2
        error_vsini_std = np.sqrt(np.sum(error_vsini)/len(vsini))
        return error_teff_std, error_logg_std, error_feh_std, error_vsini_std

    #Define step for synthesis according to observations
    kwargs['step_wave'] = round(float(delta_l),4)
    model = kwargs['model']
    y_obserr = 0.01 #arbitary value

    fix_teff  = 1 if kwargs['fix_teff']  else 0
    fix_logg  = 1 if kwargs['fix_logg']  else 0
    fix_feh   = 1 if kwargs['fix_feh']   else 0
    fix_vt    = 1 if kwargs['fix_vt']    else 0
    fix_vmac  = 1 if kwargs['fix_vmac']  else 0
    fix_vsini = 1 if kwargs['fix_vsini'] else 0

    # Set PARINFO structure for all 6 free parameters for mpfit
    # Teff, logg, feh, vt, vmac, vsini
    # The limits are also cheched by the bounds function
    teff_info  = {'parname':'Teff',   'limited': [1, 1], 'limits': [parinfo_limit(model)[0], parinfo_limit(model)[1]], 'step': 100,  'mpside': 2, 'fixed': fix_teff}
    logg_info  = {'parname':'logg',   'limited': [1, 1], 'limits': [parinfo_limit(model)[2], parinfo_limit(model)[3]], 'step': 0.1,  'mpside': 2, 'fixed': fix_logg}
    feh_info   = {'parname':'[Fe/H]', 'limited': [1, 1], 'limits': [parinfo_limit(model)[4], parinfo_limit(model)[5]], 'step': 0.05, 'mpside': 2, 'fixed': fix_feh}
    vt_info    = {'parname':'vt',     'limited': [1, 1], 'limits': [0.0, 9.99], 'step': 0.5,  'mpside': 2, 'fixed': fix_vt}
    vmac_info  = {'parname':'vmac',   'limited': [1, 1], 'limits': [0.0, 20.0], 'step': 0.5,  'mpside': 2, 'fixed': fix_vmac}
    vsini_info = {'parname':'vsini',  'limited': [1, 1], 'limits': [0.0, 99.0], 'step': 1.0,  'mpside': 2, 'fixed': fix_vsini}
    parinfo = [teff_info, logg_info, feh_info, vt_info, vmac_info, vsini_info]

    # A dictionary which contains the parameters to be passed to the
    # user-supplied function specified by myfunct via the standard Python
    # keyword dictionary mechanism. This is the way you can pass additional
    # data to your user-supplied function without using global variables.
    fa = {'x_obs': x_obs, 'ranges': ranges, 'model': model, 'y': y_obs, 'y_obserr': y_obserr, 'options': kwargs}

    # Minimization starts here.
    # Measure time
    start_time = time.time()
    m = mpfit(myfunct, xall=p0, parinfo=parinfo, ftol=1e-4, xtol=1e-4, gtol=1e-4, functkw=fa, maxiter=20)
    # Output
    dof = len(y_obs) - len(m.params)
    parameters_1 = convergence_info(m, parinfo, dof)

    if kwargs['refine']:
        print('Refining the parameters...')
        print('Patience is the key...')
        #kwargs['flag_vt']   = True
        #kwargs['flag_vmac'] = True
        m.params[3] = getMic(m.params[0], m.params[1], m.params[2])
        m.params[4] = getMac(m.params[0], m.params[1])
        x_s, y_s = func(m.params, atmtype=model, driver='synth', ranges=ranges, **kwargs)
        x_o, y_o = exclude_bad_points(x_obs, y_obs, x_s, y_s)

        fa = {'x_obs': x_o, 'ranges': ranges, 'model': model, 'y': y_o, 'y_obserr': y_obserr, 'options': kwargs}
        f = mpfit(myfunct, xall=m.params, parinfo=parinfo, ftol=1e-3, xtol=1e-3, gtol=1e-3, functkw=fa, maxiter=20)
        # Output
        dof = len(y_o) - len(f.params)
        parameters_2 = convergence_info(f, parinfo, dof)

        #error estimation
        if kwargs['errors']:
            kwargs['flag_vt']   = False
            kwargs['flag_vmac'] = False
            kwargs['refine']    = False
            teff_error, logg_error, feh_error, vsini_error = error_synth(f.params, **kwargs)
            end_time = time.time() - start_time
            print('Minimization finished in %s sec' % int(end_time))
            parameters = parameters_2 + parameters_1 + [int(teff_error)] + [round(logg_error,2)] + [round(feh_error,3)] + [round(vsini_error,2)] + [int(end_time)]
        else:
            end_time = time.time() - start_time
            print('Minimization finished in %s sec' % int(end_time))
            parameters = parameters_2 + parameters_1 + [0] + [0] + [0] + [0] + [int(end_time)]
    else:
        #error estimation
        if kwargs['errors']:
            kwargs['flag_vt']   = False
            kwargs['flag_vmac'] = False
            kwargs['refine']    = False
            teff_error, logg_error, feh_error, vsini_error = error_synth(m.params, **kwargs)
            end_time = time.time() - start_time
            print('Minimization finished in %s sec' % int(end_time))
            parameters = parameters_1 + parameters_1 + [round(teff_error,1)] + [round(logg_error,2)] + [round(feh_error,3)] + [round(vsini_error,2)] + [int(end_time)]
        else:
            end_time = time.time() - start_time
            print('Minimization finished in %s sec' % int(end_time))
            parameters = parameters_1 + parameters_1 + [0] + [0] + [0] + [0] + [int(end_time)]
        x_o, y_o = x_obs, y_obs

    return parameters, x_o, y_o
