# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import os
from astropy.io import fits
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def local_norm(obs_fname, r, snr=None, method='linear', plot=False):
    '''Local Normalisation function. Make a linear fit from the maximum points
    of each segment.

    Input
    -----
    obs_fname : str
      Observations file
    r : list
      range of the interval
    method : str
      'scalar' : Normalize with median value
      'linear' : Normalize with a linear fit
    plot : bool
      Show plot

    Output
    ------
    wave : ndarray
      Wavelength
    new_flux : ndarray
      Normalized flux
    delta_l : float
      Wavelength step
    '''

    # Define the area of Normalization
    start_norm = r[0]
    end_norm = r[1]
    # Transform SNR to noise
    if snr is None:
        noise = 0.0
    else:
        snr = float(snr)
        noise = 1.0/(snr)

    # Read observations
    wave_obs, flux_obs, delta_l = read_observations(obs_fname, start_norm, end_norm)

    # Divide in 2 and find the maximum points
    y = np.array_split(flux_obs, 2)
    x = np.array_split(wave_obs, 2)
    index_max1 = np.sort(np.argsort(y[0])[-8:])
    index_max2 = np.sort(np.argsort(y[1])[-8:])
    f_max1 = y[0][index_max1]
    f_max2 = y[1][index_max2]

    w_max1 = x[0][index_max1]
    w_max2 = x[1][index_max2]

    f_max = np.concatenate((f_max1, f_max2))
    w_max = np.concatenate((w_max1, w_max2))

    if method == 'scalar':
        # Divide with the median of maximum values.
        new_flux = flux_obs/np.median(f_max)
        if snr<20:
            new_flux =  new_flux + (2.0*noise)
        elif 20<=snr<200:
            new_flux =  new_flux + (1.5*noise)
        elif 200<=snr<350:
            new_flux =  new_flux + (1.0*noise)
        elif 350<=snr:
            new_flux =  new_flux + (0.0*noise)
    if method == 'linear':
        z = np.polyfit(w_max, f_max, 1)
        p = np.poly1d(z)
        new_flux = flux_obs/p(wave_obs)
        if snr<20:
            new_flux =  new_flux + (2.0*noise)
        elif 20<=snr<200:
            new_flux =  new_flux + (1.5*noise)
        elif 200<=snr<350:
            new_flux =  new_flux + (1.0*noise)
        elif 350<=snr:
            new_flux =  new_flux + (0.0*noise)

    # Exclude some continuum points which differ 0.5% from continuum level
    wave = wave_obs[np.where((wave_obs >= float(r[0])) & (wave_obs <= float(r[1])))]
    new_flux = new_flux[np.where((wave_obs >= float(r[0])) & (wave_obs <= float(r[1])))]

    if plot:
        plt.plot(wave_obs, flux_obs)
        x = [start_norm, end_norm]
        y = [np.median(f_max), np.median(f_max)]
        plt.plot(x, y)
        plt.plot(w_max, f_max, 'o')
        plt.show()  # Are you sure you want this here?

        y = [1.0, 1.0]
        plt.plot(x, y)
        plt.plot(wave, new_flux)
        plt.show()
    return wave, new_flux, delta_l


def read_observations(fname, start_synth, end_synth):
    '''Read observed spectrum of different types and return wavelength and flux.
    Input
    -----
    fname : str
      filename of the spectrum. Currently only fits and text files accepted.
    start_synth : float
      starting wavelength where the observed spectrum is cut
    end_synth : float
      ending wavelength where the observed spectrum is cut

    Output
    -----
    wavelength_obs : ndarray
      observed wavelength
    flux_obs : ndarray
      observed flux
    '''
    # These are the approved formats
    extension = ('.dat', '.txt', '.spec', '.fits')
    if fname.endswith(extension):
        if (fname[-4:] == '.dat') or (fname[-4:] == '.txt'):
            with open(fname, 'r') as f:
                lines = (line for line in f if not line[0].isalpha())  # skip header
                wave, flux = np.loadtxt(lines, unpack=True, usecols=(0, 1))  # Not sure this is right. Maybe np.array(lines)

        elif fname[-5:] == '.fits':
            hdulist = fits.open(fname)
            header = hdulist[0].header
            # Only 1-D spectrum accepted.
            flux = hdulist[0].data  # flux data in the primary
            flux = np.array(flux, dtype=np.float64)
            start_wave = header['CRVAL1']  # initial wavelenght
            step = header['CDELT1']  # increment per pixel
            w0, dw, n = start_wave, step, len(flux)
            w = start_wave + step * n
            wave = np.linspace(w0, w, n, endpoint=False)
        # These types are produced by FASMA (fits format).
        elif fname[-5:] == '.spec':
            hdulist = fits.open(fname)
            x = hdulist[1].data
            flux = x['flux']  # flux data in the primary
            wave = x['wavelength']
        # Cut observations to the intervals of the synthesis
        delta_l = wave[1] - wave[0]
        wavelength_obs = wave[np.where((wave >= float(start_synth)) & (wave <= float(end_synth)))]
        flux_obs = flux[np.where((wave >= float(start_synth)) & (wave <= float(end_synth)))]
    else:
        print('Spectrum is not in acceptable format. Convert to ascii of fits.')
        wavelength_obs, flux_obs, delta_l = (None, None, None)
    return wavelength_obs, flux_obs, delta_l


def read_obs_intervals(obs_fname, r, snr=100, method='linear'):
    '''Read only the spectral chunks from the observed spectrum
    This function does the same as read_observations but for the whole linelist.
    Input
    -----
    fname : filename of the spectrum. Currently only fits and text files accepted.
    r : ranges of wavelength intervals where the observed spectrum is cut
    (starting and ending wavelength)

    Output
    -----
    wavelength_obs : observed wavelength
    flux_obs : observed flux
    '''

    # N : number of intervals
    N = len(r)
    x_obs = []
    y_obs = []
    for i in range(N):
        # Obtain the normalized spectrum
        x_obs.append(local_norm(obs_fname, r[i], snr, method)[0])
        y_obs.append(local_norm(obs_fname, r[i], snr, method)[1])

    x_obs = np.hstack(x_obs)
    y_obs = np.hstack(y_obs)
    if any(i == 0 for i in y_obs):
        print('Warning: Flux contains 0 values.')

    delta_l = local_norm(obs_fname, r[0], snr, method)[2]
    print('SNR: %s' % snr)
    return x_obs, y_obs, delta_l


def plot(x_obs=None, y_obs=None, x=None, y=None, res=False):
    '''Function to plot synthetic spectrum.
    Input
    -----
    x_obs : ndarray
      observed wavelength
    y_obs : ndarray
      observed flux
    x : ndarray
      synthetic wavelength
    y : ndarray
      synthetic flux
    res : bool
      If True, then calculate the difference in flux

    Output
    ------
    plots
    '''

    # if nothing exists, exit
    if (x_obs is None) and (x is None):
        return
    # if there is not observed spectrum, plot only synthetic (case 1, 3)
    if x_obs is None:
        plt.plot(x, y, label='synthetic')
        if res:
            sl = InterpolatedUnivariateSpline(x, y, k=1)
            ymodel = sl(x_obs)
            plt.plot(x_obs, (y_obs-ymodel)*10, label='residuals')
        plt.legend(loc='best', frameon=False)
        plt.show()
    # if both exist
    else:
        plt.plot(x, y, label='synthetic')
        plt.plot(x_obs, y_obs, label='observed')
        if res:
            sl = InterpolatedUnivariateSpline(x, y, k=1)
            ymodel = sl(x_obs)
            plt.plot(x_obs, (y_obs-ymodel)*10, label='residuals')

        plt.legend(loc='best', frameon=False)
        plt.show()
    return


def snr(fname, plot=False):
    '''Calculate SNR using intervals depending on giraffe mode.
    Input
    ----
    fname : str
      Name of the spectrum
    plot : bool
      Plot snr fit

    Output
    -----
    snr : float
      snr value averaged from the continuum intervals
    '''
    from PyAstronomy import pyasl
    def sub_snr(interval):
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
        w1, w2 = interval
        wave_cut, flux_cut, l = read_observations(fname, w1, w2)
        num_points = int(len(flux_cut)/4)
        if num_points != 0:
            snrEsti1 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=2, controlPlot=plot)
            return snrEsti1["SNR-Estimate"]
        else:
            return 0


    intervals = [[5744, 5746],   [6048, 6052],    [6068, 6076],     [6682, 6686],
                 [6649, 6652],   [6614, 6616],    [5374.5, 5376.5], [5438.5, 5440],
                 [5449.5, 5051], [5458, 5459.25], [5498.3, 5500],   [5541.5, 5542.5]]

    snr = []
    for interval in intervals:
        snr.append(sub_snr(interval))

    if not len(snr):
        snr = None
    else:
        snr = [value for value in snr if value != 0]
        snr_clean = [value for value in snr if not np.isnan(value)]
        snr_total = np.average(snr_clean)
        snr = round(snr_total, 1)
    return snr
