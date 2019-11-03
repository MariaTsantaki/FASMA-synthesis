# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import os
from astropy.io import fits
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def mad(data, axis=None):
    return np.median(np.absolute(data - np.median(data, axis)), axis)

def local_norm(obs_fname, r, snr, method='linear', lol=1.0, plot=False):
    '''Local Normalisation function. Make a linear fit from the maximum points
    of each segment.
    Input
    -----
    obs_fname : observations file
    r : range of the interval

    Output
    ------
    new_flux : normalized flux
    '''

    # Define the area of Normalization
    start_norm = r[0] - 1.0
    end_norm   = r[1] + 1.0
    #Transform SNR to noise
    if snr is None:
        noise = 0.0
    else:
        snr = float(snr)
        noise = 1.0/snr
    #Read observations
    wave_obs, flux_obs, delta_l = read_observations(obs_fname, start_norm, end_norm)

    flux_obs = flux_obs/np.median(flux_obs)
    # Clean for cosmic rays
    med = np.median(flux_obs)
    sig = mad(flux_obs)
    flux_clean = np.where(flux_obs < (med + (sig*3.0)), flux_obs, med)
    flux_obs = flux_clean

    pol_fit = np.polyfit(wave_obs, flux_obs, 1)
    fit_line = np.poly1d(pol_fit)
    for i in range(5):
        condition = flux_obs - fit_line(wave_obs) + noise > 0
        cont_points_wl = wave_obs[condition]
        cont_points_fl = flux_obs[condition]
        pol_fit_new = np.polyfit(cont_points_wl, cont_points_fl, 1)
        fit_line = np.poly1d(pol_fit_new)

    new_flux = flux_obs/fit_line(wave_obs)
    wave     = wave_obs[np.where((wave_obs >= float(r[0])) & (wave_obs <= float(r[1])))]
    new_flux = new_flux[np.where((wave_obs >= float(r[0])) & (wave_obs <= float(r[1])))]

    if plot:
        plt.plot(wave_obs, flux_obs, label='raw spectrum')
        plt.plot(cont_points_wl, cont_points_fl, 'o')
        plt.xlabel(r'Wavelength $\AA{}$')
        plt.ylabel('Normalized flux')
        plt.legend(loc='best', frameon=False)
        plt.grid(True)
        plt.show()

        x = [start_norm, end_norm]
        y = [1.0, 1.0]
        plt.plot(x, y)
        plt.plot(wave, new_flux, label='normalized')
        plt.xlabel(r'Wavelength $\AA{}$')
        plt.ylabel('Normalized flux')
        plt.legend(loc='best', frameon=False)
        plt.grid(True)
        plt.show()
    return wave, new_flux, delta_l

def read_observations(fname, start_synth, end_synth):
    """Read observed spectrum of different types and return wavelength and flux.
    Input
    -----
    fname : filename of the spectrum. These are the approved formats: '.dat', '.txt',
    '.spec', '.fits'.
    start_synth : starting wavelength where the observed spectrum is cut
    end_synth : ending wavelength where the observed spectrum is cut

    Output
    -----
    wave_obs : raw observed wavelength
    flux_obs : raw observed flux
    """

    extension = ('.dat', '.txt', '.spec', '.fits')
    if fname.endswith(extension):
        if (fname[-4:] == '.dat') or (fname[-4:] == '.txt'):
            with open(fname, 'r') as f:
                lines = (line for line in f if not line[0].isalpha())  # skip header
                wave, flux = np.loadtxt(lines, unpack=True, usecols=(0, 1))
        elif fname[-5:] == '.fits':
            hdulist = fits.open(fname)
            header = hdulist[0].header
            # Only 1-D spectrum accepted.
            flux = hdulist[0].data  # flux data in the primary
            flux = np.array(flux, dtype=np.float64)
            start_wave = header['CRVAL1']  # initial wavelenght
            # step = header['CD1_1'] # step in wavelenght
            step = header['CDELT1']  # increment per pixel
            n = len(flux)
            w = start_wave + step * n
            wave = np.linspace(start_wave, w, n, endpoint=False)
        # These types are produced by FASMA (fits format).
        elif fname[-5:] == '.spec':
            hdulist = fits.open(fname)
            x = hdulist[1].data
            flux = x['flux']
            wave = x['wavelength']
        # Cut observations to the intervals of the synthesis
        delta_l = wave[1] - wave[0]
        wave_obs = wave[np.where((wave >= float(start_synth)) & (wave <= float(end_synth)))]
        flux_obs = flux[np.where((wave >= float(start_synth)) & (wave <= float(end_synth)))]
    else:
        print('Spectrum is not in acceptable format. Convert to ascii or fits.')
        wave_obs, flux_obs, delta_l = (None, None, None)
    return wave_obs, flux_obs, delta_l

def read_obs_intervals(obs_fname, r, snr=100, method='linear'):
    """Read only the spectral chunks from the observed spectrum
    This function does the same as read_observations but for the whole linelist.
    Input
    -----
    fname : filename of the spectrum.
    r : ranges of wavelength intervals where the observed spectrum is cut
    (starting and ending wavelength)

    Output
    -----
    xobs : observed normalized wavelength
    yobs : observed normalized flux
    """

    lol = 1.0 # This is here for no reason.
    # Obtain the normalized spectrum
    spec = [local_norm(obs_fname, ri, snr, method, lol) for ri in r]
    xobs = np.hstack(np.vstack(spec).T[0])
    yobs = np.hstack(np.vstack(spec).T[1])
    delta_l = spec[0][2]
    if any(i == 0 for i in yobs):
        print('Warning: Flux contains 0 values.')

    print('SNR: %s' % snr)
    return xobs, yobs, delta_l

def plot(xobs, yobs, x, y, res=False):
    """Function to plot synthetic spectrum.
    Input
    -----
    xobs : observed wavelength
    yobs : observed flux
    x : synthetic wavelength
    y : synthetic flux

    Output
    ------
    plots
    """

    # if nothing exists, pass
    if (xobs is None) and (x is None):
        pass
    # if there is not observed spectrum, plot only synthetic (case 1, 3)
    if xobs is None:
        plt.plot(x, y, label='synthetic')
    # if both exist
    else:
        plt.plot(x, y, label='synthetic')
        plt.plot(xobs, yobs, label='observed')
        if res:
            sl = InterpolatedUnivariateSpline(x, y, k=1)
            ymodel = sl(xobs)
            plt.plot(xobs, (yobs-ymodel)*10, label='residuals')
    plt.xlabel(r'Wavelength $\AA{}$')
    plt.ylabel('Normalized flux')
    plt.legend(loc='best', frameon=False)
    plt.grid(True)
    plt.show()
    return

def snr(fname, plot=False):
    """Calculate SNR using for various intervals.
    Input
    ----
    fname : spectrum
    plot : plot snr fit
    Output
    -----
    snr : snr value averaged from the continuum intervals
    """

    from PyAstronomy import pyasl

    def sub_snr(snr_region):
        '''Measure the SNR on a small interval
        Input
        -----
        interval : list
          Upper and lower limit on wavelength

        Output
        ------
        SNR : float
          The SNR in the small interval
        '''

        w1, w2 = snr_region
        wave_cut, flux_cut, l = read_observations(fname, w1, w2)

        # Clean for cosmic rays
        med = np.median(flux_cut)
        sig = mad(flux_cut)
        flux_obs = np.where(flux_cut < (med + (sig*3.0)), flux_cut, med)
        pol_fit = np.polyfit(wave_cut, flux_obs, 1)
        fit_line = np.poly1d(pol_fit)
        for i in range(10):
            condition = abs(flux_obs - fit_line(wave_cut)) < 3*sig
            cont_points_wl = wave_cut[condition]
            cont_points_fl = flux_obs[condition]

        num_points = int(len(cont_points_fl)/2.0)
        if num_points != 0:
            snrEstimate = pyasl.estimateSNR(cont_points_wl, cont_points_fl, num_points, deg=2, controlPlot=plot)
            return snrEstimate["SNR-Estimate"]
        else:
            return 0

    snr_regions = [[5744, 5746], [6048, 6052], [6068, 6076], [6682, 6686], [6649, 6652],
                [6614, 6616], [5438.5, 5440], [5449.5, 5051], [5458, 5459.25],
                [5498.3, 5500], [5541.5, 5542.5]]

    snr = [sub_snr(snr_region) for snr_region in snr_regions]

    if len(snr):
        snr = [value for value in snr if value != 0]
        snr_clean = [value for value in snr if not np.isnan(value)]
        snr_total = np.average(snr_clean)
        snr = int(snr_total)
        return snr
    else:
        return None
