# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
import os
from astropy.io import fits
import matplotlib.pyplot as plt


def local_norm(obs_fname, r, snr, method='linear', plot=False):
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

    # Define the area of Normalization, 10A around the center of each interval
    center = (r[0] + r[1])/2.0
    start_norm = center - 10.0
    end_norm = center + 10.0
    start_norm = r[0]
    end_norm = r[1]
    #Transform SNR to noise
    noise = 1.0/(2.0*snr)
    #Read observations
    wave_obs, flux_obs = read_observations(obs_fname, start_norm, end_norm)

    # Divide in 2 and find the maximum points
    y = np.array_split(flux_obs, 2)
    x = np.array_split(wave_obs, 2)
    index_max1 = np.sort(np.argsort(y[0])[-5:])  # this can be done better
    index_max2 = np.sort(np.argsort(y[1])[-5:])  # this can be done better
    # index_max = np.sort(np.argsort(flux_obs)[-10:]) # this can be done better
    f_max1 = y[0][index_max1]
    f_max2 = y[1][index_max2]

    w_max1 = x[0][index_max1]
    w_max2 = x[1][index_max2]

    # f_max = flux_obs[index_max]
    # w_max = wave_obs[index_max]
    f_max = np.concatenate((f_max1, f_max2))
    f_max = f_max - (f_max*noise)
    #std = np.std(f_max, ddof=1)
    w_max = np.concatenate((w_max1, w_max2))

    if method == 'scalar':
        # Divide with the median of maximum values.
        new_flux = flux_obs/np.median(f_max)

    if method == 'linear':
        z = np.polyfit(w_max, f_max, 1)
        p = np.poly1d(z)
        #f = p(wave_obs)
        new_flux = flux_obs/p(wave_obs)

    wave = wave_obs[np.where((wave_obs >= float(r[0])) & (wave_obs <= float(r[1])))]
    new_flux = new_flux[np.where((wave_obs >= float(r[0])) & (wave_obs <= float(r[1])))]
    if plot:
        plt.plot(wave_obs, flux_obs)
        x = [center, center]
        y = [np.median(f_max), np.median(f_max)]
        plt.plot(x, y)
        plt.plot(w_max, f_max, 'o')
        plt.show()

        plt.plot(wave, new_flux)
        plt.show()
    return wave, new_flux


def read_observations(fname, start_synth, end_synth):
    """Read observed spectrum of different types and return wavelength and flux.
    Input
    -----
    fname : filename of the spectrum. Currently only fits and text files accepted.
    start_synth : starting wavelength where the observed spectrum is cut
    end_synth : ending wavelength where the observed spectrum is cut

    Output
    -----
    wavelength_obs : observed wavelength
    flux_obs : observed flux
    """
    # These are the approved formats
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
            # step = header['CD1_1'] #step in wavelenght
            step = header['CDELT1']  # increment per pixel
            w0, dw, n = start_wave, step, len(flux)
            w = start_wave + step * n
            wave = np.linspace(w0, w, n, endpoint=False)
        # T hese types are produced by MOOGme (fits format).
        elif fname[-5:] == '.spec':
            hdulist = fits.open(fname)
            x = hdulist[1].data
            flux = x['flux']  # flux data in the primary
            wave = x['wavelength']
        # Cut observations to the intervals of the synthesis
        wavelength_obs = wave[np.where((wave >= float(start_synth)) & (wave <= float(end_synth)))]
        flux_obs = flux[np.where((wave >= float(start_synth)) & (wave <= float(end_synth)))]

    else:
        print('Spectrum is not in acceptable format. Convert to ascii of fits.')
        wavelength_obs, flux_obs = (None, None)
    return wavelength_obs, flux_obs


def read_obs_intervals(obs_fname, r, snr=100, method='scalar'):
    """Read only the spectral chunks from the observed spectrum
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
    """

    # N : number of intervals
    N = len(r)
    spec = []
    for i in range(N):
        # Obtain the normalized spectrum
        spec.append(local_norm(obs_fname, r[i], snr, method))

    x_obs = np.column_stack(spec)[0]
    y_obs = np.column_stack(spec)[1]

    if any(i == 0 for i in y_obs):
        print('Warning: Flux contains 0 values.')
    return x_obs, y_obs


def plot(x_obs, y_obs, x, y, res=False):
    """Function to plot synthetic spectrum.
    Input
    -----
    x_obs : observed wavelength
    y_obs : observed flux
    x : synthetic wavelength
    y : synthetic flux

    Output
    ------
    plots
    """

    # if nothing exists, pass
    if (x_obs is None) and (x is None):
        pass
    # if there is not observed spectrum, plot only synthetic (case 1, 3)
    if x_obs is None:
        plt.plot(x, y, label='synthetic')
        if res:
            sl = InterpolatedUnivariateSpline(x, y, k=1)
            ymodel = sl(x_obs)
            plt.plot(x_obs, (y_obs-ymodel)*10, label='residuals')
        plt.legend()
        plt.show()
    # if both exist
    else:
        plt.plot(x, y, label='synthetic')
        plt.plot(x_obs, y_obs, label='observed')
        if res:
            sl = InterpolatedUnivariateSpline(x, y, k=1)
            ymodel = sl(x_obs)
            plt.plot(x_obs, (y_obs-ymodel)*10, label='residuals')

        plt.legend()
        plt.show()
    return


def snr(fname, plot=False):
    """Calculate SNR using intervals depending on giraffe mode.
    Input
    ----
    fname : spectrum
    plot : plot snr fit
    Output
    -----
    snr : snr value averaged from the continuum intervals
    """
    # I know I know, there is a better way to do this!
    from PyAstronomy import pyasl

    snr = []
    wave_cut, flux_cut = read_observations(fname, 6681, 6690)
    num_points = int(len(flux_cut)/2)
    if num_points != 0:
        snrEsti1 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=1, controlPlot=plot)
        snr.append(snrEsti1["SNR-Estimate"])
    else:
        pass

    wave_cut, flux_cut = read_observations(fname, 6649, 6652)
    num_points = int(len(flux_cut)/2)
    if num_points != 0:
        snrEsti2 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=1, controlPlot=plot)
        snr.append(snrEsti2["SNR-Estimate"])
    else:
        pass

    wave_cut, flux_cut = read_observations(fname, 6614, 6623)
    num_points = int(len(flux_cut)/2)
    if num_points != 0:
        snrEsti3 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=1, controlPlot=plot)
        snr.append(snrEsti3["SNR-Estimate"])
    else:
        pass

    wave_cut, flux_cut = read_observations(fname, 5374.5, 5376.5)
    num_points = int(len(flux_cut)/2)
    if num_points != 0:
        snrEsti1 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=1, controlPlot=plot)
        snr.append(snrEsti1["SNR-Estimate"])
    else:
        pass

    wave_cut, flux_cut = read_observations(fname, 5438, 5440)
    num_points = int(len(flux_cut)/2)
    if num_points != 0:
        snrEsti1 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=1, controlPlot=plot)
        snr.append(snrEsti1["SNR-Estimate"])
    else:
        pass

    wave_cut, flux_cut = read_observations(fname, 5449.5, 5051)
    num_points = int(len(flux_cut)/2)
    if num_points != 0:
        snrEsti1 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=1, controlPlot=plot)
        snr.append(snrEsti1["SNR-Estimate"])
    else:
        pass

    wave_cut, flux_cut = read_observations(fname, 5458, 5459.5)
    num_points = int(len(flux_cut)/2)
    if num_points != 0:
        snrEsti1 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=1, controlPlot=plot)
        snr.append(snrEsti1["SNR-Estimate"])
    else:
        pass

    wave_cut, flux_cut = read_observations(fname, 5498.3, 5500)
    num_points = int(len(flux_cut)/2)
    if num_points != 0:
        snrEsti1 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=1, controlPlot=plot)
        snr.append(snrEsti1["SNR-Estimate"])
    else:
        pass

    wave_cut, flux_cut = read_observations(fname, 5541.5, 5542.5)
    num_points = int(len(flux_cut)/2)
    if num_points != 0:
        snrEsti1 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=1, controlPlot=plot)
        snr.append(snrEsti1["SNR-Estimate"])
    else:
        pass

    wave_cut, flux_cut = read_observations(fname, 8643.8, 8645.8)
    num_points = int(len(flux_cut)/2)
    if num_points != 0:
        snrEsti1 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=1, controlPlot=plot)
        snr.append(snrEsti1["SNR-Estimate"])
    else:
        pass

    wave_cut, flux_cut = read_observations(fname, 8760.0, 8762.2)
    num_points = int(len(flux_cut)/2)
    if num_points != 0:
        snrEsti1 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=1, controlPlot=plot)
        snr.append(snrEsti1["SNR-Estimate"])
    else:
        pass
    wave_cut, flux_cut = read_observations(fname, 8870.0, 8872.0)
    num_points = int(len(flux_cut)/2)
    if num_points != 0:
        snrEsti1 = pyasl.estimateSNR(wave_cut, flux_cut, num_points, deg=1, controlPlot=plot)
        snr.append(snrEsti1["SNR-Estimate"])
    else:
        pass

    snr_total = sum(snr)/len(snr)
    return round(snr_total,1)


# The rest are useless functions for now....
def plot_synth(fname):
    """Function to plot synthetic spectrum
    """
    from synthetic import _read_raw_moog, _read_moog

    if fname == 'smooth.out':
        x, y = _read_moog(fname='smooth.out')
        plt.plot(x, y)
        plt.show()
        plt.close()
    elif fname == 'summary.out':
        x, y = _read_raw_moog('summary.out')
        # z = pyasl.instrBroadGaussFast(x, y, 50000, edgeHandling="firstlast")
        plt.plot(x, y)
        # plt.plot(x, z)
        plt.show()
        plt.close()
    else:
        if os.path.isfile(fname):
            x, y = _read_moog(fname)
            plt.plot(x, y)
            plt.show()
            plt.close()
        else:
            print('Synthetic spectrum does not exist.')
    return


def interpol_synthetic(wave_obs, wave_synth, flux_synth):
    """Interpolation of the synthetic flux to the observed wavelength"""
    sl = InterpolatedUnivariateSpline(wave_synth, flux_synth, k=1)
    int_flux = sl(wave_obs)
    return int_flux


def normalize(x, y, num=20, k=1):
    index_max = np.sort(np.argsort(y)[-num:])
    f_obs_max = y[index_max]
    w_obs_max = x[index_max]
    sl = InterpolatedUnivariateSpline(w_obs_max, f_obs_max, k=1)
    continuum = sl(x)
    norm = y/continuum
    return x, norm


def chi2(wave_obs, flux_obs, wave_synth, flux_synth):
    """Some x2 statistics"""

    # Interpolation of the synthetic flux to the observed wavelength
    sl = InterpolatedUnivariateSpline(wave_synth, flux_synth, k=1)
    int_flux = sl(wave_obs)
    f = interp1d(wave_synth, flux_synth, kind='cubic')
    error = 1./200
    chi2 = np.sum(((flux_obs - int_flux)/error)**2)
    # which or the 2 is the correct format?
    chi = ((flux_obs - int_flux)**2)
    chi2 = np.sum(chi)
    print('This is your chi2 value: '), chi2
    return
