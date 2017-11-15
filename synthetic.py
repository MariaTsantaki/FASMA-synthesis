# -*- coding: utf8 -*-

# My imports
from __future__ import division
import os
import numpy as np
import pandas as pd
from astropy.io import fits
from PyAstronomy import pyasl
from scipy.integrate import quad
from scipy.signal import fftconvolve
from scipy.interpolate import InterpolatedUnivariateSpline

def save_synth_spec(x, y, y_obs=None, initial=None, final=None, fname='initial.spec', **options):
    '''Save synthetic spectrum of all intervals

    Input
    ----
    x : np.ndarray
      Wavelength
    y : np.ndarray
      Flux
    y_obs : np.ndarray
      Observed flux
    initial : list
      Initial parameters
    final : list
      Final parameters
    fname : str
      Filename of fits file
    options : dict
      Option dictionary from 'synthDriver'

    Output
    -----
    fname fits file
    '''
    # Create header
    header = fits.Header()
    header['CRVAL1']   = x[0]
    header['CDELT1']   = x[1] - x[0]
    header['Teff_in']  = initial[0]
    header['logg_in']  = initial[1]
    header['FeH_in']   = initial[2]
    header['vt_in']    = initial[3]
    header['vmac_in']  = initial[4]
    header['vsini_in'] = initial[5]
    header['Mod_atmo'] = options['model']
    header['Damping']  = options['damping']
    header['Interval'] = options['inter_file']
    header['Resol']    = options['resolution']

    if final:
        header['Teff_f']  = final[0]
        header['logg_f']  = final[1]
        header['FeH_f']   = final[2]
        header['vt_f']    = final[3]
        header['vmac_f']  = final[4]
        header['vsini_f'] = final[5]
        header['Obs']     = options['observations']
        header['SNR']     = options['snr']

    if options['observations'] and (final is None):
        fname = options['observations'].split('/')[-1]
        fname = fname.split('.')[0] + '_input.spec'
    elif final:
        fname = options['observations'].split('/')[-1]
        fname = fname.split('.')[0] + '_output.spec'
    else:
        fname = str(initial[0]) + '_' + str(initial[1]) + '_' + str(initial[2]) + '_' + str(initial[3]) + '_' + str(initial[4]) + '_' + str(initial[5]) +  '_' + str(options['resolution']) + '.spec'

    tbhdu = fits.BinTableHDU.from_columns([fits.Column(name='wavelength', format='D', array=x),
                                           fits.Column(name='flux', format='D', array=y),
                                           fits.Column(name='y_obs', format='D', array=y_obs)], header=header)
    tbhdu.writeto('results/%s' % fname, clobber=True)
    print('Synthetic spectrum saved: results/%s' % fname)


def broadening(x, y, vsini, vmac, resolution=None, epsilon=0.60):
    '''This function broadens the given data using velocity kernels,
    e.g. instrumental profile, vsini and vmac.
    Based on http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/broadening.html
    Input
    ----
    x : np.ndarray
      wavelength
    y : np.ndarray
      flux
    vsini : float
      vsini in km/s
    vmac : float
      vmac in km/s
    resolution : float
      Instrumental resolution (lambda /delta lambda)
    epsilon : float
      Linear limb-darkening coefficient (0-1).  Linear limb-darkening coefficient (0-1).

    Output
    -----
    x : np.ndarray
      Same wavelength
    y_broad : np.ndarray
      Broadened flux
    '''

    def instrumental_profile(x, y, resolution):
        '''
        Inputs
        -----
        x, y : The abscissa and ordinate of the data.
        sigma : The width (i.e., standard deviation) of the Gaussian profile
        used in the convolution.
        edgeHandling : None, "firstlast". Determines the way edges will be
        handled. If None, nothing will be done about it. If set to "firstlast",
        the spectrum will be extended by using the first and last value at the
        start or end. Note that this is not necessarily appropriate.
        The default is None.
        maxsig : The extent of the broadening kernel in terms of standrad
        deviations. By default, the Gaussian broadening kernel will be extended
        over the entire given spectrum, which can cause slow evaluation in the
        case of large spectra. A reasonable choice could, e.g., be five.

        Output
        -----
        y_inst : convolved flux
        '''

        # Deal with zero or None values seperately
        if (resolution is None) or (resolution == 0):
            return y
        else:
            y_inst = pyasl.instrBroadGaussFast(x, y, resolution, edgeHandling="firstlast", fullout=False, maxsig=None)
            return y_inst

    def vsini_broadening(x, y, epsilon, vsini):
        '''
        Apply rotational broadening to a spectrum assuming a linear limb darkening
        law. The adopted limb darkening law is the linear one, parameterize by the
        linear limb darkening parameter: epsilon = 0.6.
        The effect of rotational broadening on the spectrum is
        wavelength dependent, because the Doppler shift depends
        on wavelength. This function neglects this dependence, which
        is weak if the wavelength range is not too large.
        Code from: http://www.phoebe-project.org/2.0/
        .. note:: numpy.convolve is used to carry out the convolution
              and "mode = same" is used. Therefore, the output
              will be of the same size as the input, but it
              will show edge effects.
        Input
        -----
        wvl : The wavelength
        flux : The flux
        epsilon : Linear limb-darkening coefficient (0-1).
        vsini : Projected rotational velocity in km/s.
        effWvl : The wavelength at which the broadening kernel is evaluated.
        If not specified, the mean wavelength of the input will be used.

        Output
        ------
        y_rot : convolved flux
        '''

        if vsini == 0:
            return y
        else:
            y_rot = pyasl.rotBroad(x, y, epsilon, vsini, edgeHandling='firstlast')
            return y_rot

    def vmacro_kernel(dlam, Ar, At, Zr, Zt):
        '''
        Macroturbulent velocity kernel.
        Zr == Zt = vmac
        Ar == At == 1.0
        '''
        dlam[dlam == 0] = 1e-8
        return np.array([(2*Ar*idlam/(np.sqrt(np.pi)*Zr**2) + 2*At*idlam/(np.sqrt(np.pi)*Zt**2)) *
                         quad(lambda u: np.exp(-1/u**2), 0, Zr/idlam)[0]
                         for idlam in dlam])

    def vmac_broadening(wave, flux, vmacro_rad):
        '''
        Apply macroturbulent broadening.
        The macroturbulent kernel is defined in Gray 2005.
        Same functions are used in iSpec (Blanco-Cuaresma et al. 2014)

        Input
        -----
        :parameter wave: Wavelength of the spectrum
        :parameter flux: Flux of the spectrum
        :parameter vmacro_rad: macroturbulent broadening, radial component

        Output
        ------
        y_mac : broadened flux
        '''

        # radial component is equal to the tangential component
        vmacro_tan = vmacro_rad

        if vmacro_rad == vmacro_tan == 0:
            return flux

        # Define central wavelength
        lambda0  = (wave[0] + wave[-1]) / 2.0
        vmac_rad = vmacro_rad/(299792458.*1e-3)*lambda0
        vmac_tan = vmac_rad

        # Make sure the wavelength range is equidistant before applying the convolution
        delta_wave = np.diff(wave).min()
        range_wave = wave.ptp()
        n_wave = int(range_wave/delta_wave)+1
        wave_ = np.linspace(wave[0], wave[-1], n_wave)
        flux_ = np.interp(wave_, wave, flux)
        dwave = wave_[1]-wave_[0]
        n_kernel = int(5*max(vmac_rad, vmac_tan)/dwave)
        if n_kernel % 2 == 0:
            n_kernel += 1
        # The kernel might be of too low resolution, or the the wavelength range
        # might be too narrow. In both cases, raise an appropriate error
        if n_kernel == 0:
            raise ValueError("Spectrum resolution too low for macroturbulent broadening")
        elif n_kernel > n_wave:
            raise ValueError("Spectrum range too narrow for macroturbulent broadening")
        # Construct the broadening kernel
        wave_k = np.arange(n_kernel)*dwave
        wave_k -= wave_k[-1]/2.
        kernel = vmacro_kernel(wave_k, 1.0, 1.0, vmac_rad, vmac_tan)
        kernel /= sum(kernel)

        flux_conv = fftconvolve(1.0-flux_, kernel, mode='same')
        # And interpolate the results back on to the original wavelength array,
        # taking care of even vs. odd-length kernels
        if n_kernel % 2 == 1:
            offset = 0.0
        else:
            offset = dwave / 2.0
        flux = np.interp(wave+offset, wave_, 1.0-flux_conv)
        return flux

    # Instrumental broadening
    y_inst = instrumental_profile(x, y, resolution)
    # vsini broadening
    y_rot = vsini_broadening(x, y_inst, epsilon, vsini)
    # vmac broadening
    y_broad = vmac_broadening(x, y_rot, vmac)
    return x, y_broad


def _read_raw_moog(fname='summary.out'):
    '''Read the summary.out and return them

    Inputs
    ------
    fname : str (default: summary.out)
      Filename of the output file from MOOG from summary_out

    Output
    ------
    wavelenth : ndarray
      The wavelenth vector
    flux : ndarray
      The flux vector
    '''

    with open(fname, 'r') as f:
        f.readline()
        f.readline()
        start_wave, end_wave, step, flux_step = map(float, f.readline().split())
        lines = f.readlines()

    # Remove trailing '\n' from every line in lines
    data = map(lambda s: s.strip().replace('-', ' -'), lines)
    # Convert every element to a float
    flux = map(float, ' '.join(data).split(' '))
    flux = 1.0 - np.array(flux)

    w0, dw, n = float(start_wave), float(step), len(flux)
    w = w0 + dw * n
    wavelength = np.linspace(w0, w, n, endpoint=False)
    return wavelength, flux


def _read_moog(fname='smooth.out'):
    '''Read the output of moog - synthetic spectra.

    Input
    -----
    fname : str (default: smooth.out)
      The smoothed spectrum from MOOG

    Output
    ------
    wavelength : np.ndarray
      The wavelenth vector
    flux : np.ndarray
      The flux vector
    '''

    wavelength, flux = np.loadtxt(fname, skiprows=2, usecols=(0, 1), unpack=True)
    return wavelength, flux


def read_linelist(fname, intname='intervals_hr10_15n.lst'):
    '''Read the line list return atomic data and ranges

    Input
    -----
    fname : str
      File that contains the linelist
    intname : str
      File that contains the intervals

    Output
    ------
    ranges : wavelength ranges of the linelist
    atomic : atomic data
    '''

    if not os.path.isfile('rawLinelist/%s' % intname):
        raise IOError('The interval list is not in the correct place.')
    if not os.path.isfile('rawLinelist/%s' % fname):
        raise IOError('The line list is not in the correct place.')
    lines = pd.read_csv('rawLinelist/%s' % fname,
                skiprows=1, comment='#',
                delimiter='\t', usecols=range(6),
                names=['wl', 'elem', 'excit', 'loggf', 'vdwaals', 'Do'],
                converters={'Do': lambda x : x.replace("nan"," "), 'vdwaals': lambda x : float(x)})
    lines.sort_values(by='wl', inplace=True)

    intervals = pd.read_csv('rawLinelist/%s' % intname, comment='#',
                    names=['start', 'end'], delimiter='\t')

    ranges = intervals.values
    atomic = []
    N = []
    for ri in intervals.values:
        idx = (lines.wl > ri[0]) & (lines.wl < ri[1])
        a = lines[idx]
        N.append(len(a))
        a = a.as_matrix()
        atomic.append(a)
    N = sum(N)
    atomic = np.vstack(atomic)
    print('Linelist contains %s lines in %s intervals' % (N, len(ranges)))

    # Create line list for MOOG
    fmt = ["%8s", "%9s", "%10s", "%10s", "%6s", "%6s"]
    header = 'Wavelength     ele       EP      loggf   vdwaals   Do'
    np.savetxt('linelist.moog', atomic, fmt=fmt, header=header)
    return ranges, atomic


def interpol_synthetic(wave_obs, wave_synth, flux_synth):
    '''Interpolation of the synthetic flux to the observed wavelength.
    Input
    -----
    wave_obs : np.ndarray
      Observed wavelength
    wave_synth : np.ndarray
      Synthetic wavelength
    flux_synth : np.ndarray
      Synthetic flux

    Output
    ------
    int_flux : np.ndarray
      Interpolated synthetic flux
    '''
    sl = InterpolatedUnivariateSpline(wave_synth, flux_synth, k=1)
    int_flux = sl(wave_obs)
    return int_flux
