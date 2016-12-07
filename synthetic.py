# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import os
import pandas as pd
from astropy.io import fits


def save_synth_spec(x, y, y_obs=None, initial=None, final=None, fname='initial.spec', **options):
    '''Save synthetic spectrum of all intervals

    Input
    ----
    x : ndarray
      Wavelength
    y : ndarray
      Flux
    fname : str
      Filename of fits file

    Output
    -----
    fname fits file
    '''
    #Create header
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
    return


def broadening(x, y, vsini, vmac, resolution=None, epsilon=0.60):
    '''This function broadens the given data using velocity kernels,
    e.g. instrumental profile, vsini and vmac.
    Based on http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/broadening.html
    Input
    ----
    x : ndarray
      wavelength
    y : ndarray
      flux
    resolution : float
      Instrumental resolution (lambda /delta lambda)
    vsini : float
      vsini in km/s
    vmac : float
      vmac in km/s

    Output
    -----
    y_broad : ndarray
      Broadened flux
    x : ndarray
      Same wavelength
    '''

    from PyAstronomy import pyasl
    from scipy.signal import fftconvolve
    from scipy.integrate import quad

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
            y_inst = y
        else:
            y_inst = pyasl.instrBroadGaussFast(x, y, resolution, edgeHandling="firstlast", fullout=False, maxsig=None)
        return y_inst

    def vsini_broadening(x, y, epsilon, vsini):
        """
        From PyAstronomy:
        -----------------
        Apply rotational broadening to a spectrum.

        This function applies rotational broadening to a given
        spectrum using the formulae given in Gray's "The Observation
        and Analysis of Stellar Photospheres". It allows for
        limb darkening parameterized by the linear limb-darkening law.

        The `edgeHandling` parameter determines how the effects at
        the edges of the input spectrum are handled. If the default
        option, "firstlast", is used, the input spectrum is internally
        extended on both sides; on the blue edge of the spectrum, the
        first flux value is used and on the red edge, the last value
        is used to extend the flux array. The extension is neglected
        in the return array. If "None" is specified, no special care
        will be taken to handle edge effects.

        .. note:: Currently, the wavelength array as to be regularly spaced.

        Parameters
        ----------
        wvl : array
        The wavelength array [A]. Note that a
        regularly spaced array is required.
        flux : array
        The flux array.
        vsini : float
        Projected rotational velocity [km/s].
        epsilon : float
        Linear limb-darkening coefficient (0-1).
        edgeHandling : string, {"firstlast", "None"}
        The method used to handle edge effects.

        Returns
        -------
        Broadened spectrum : array
        An array of the same size as the input flux array,
        which contains the broadened spectrum.
        """
        if vsini == 0:
            y_rot = y
        else:
            y_rot = pyasl.rotBroad(x, y, epsilon, vsini, edgeHandling='firstlast')
        return y_rot

    def _vmacro_kernel(dlam, Ar, At, Zr, Zt):
        '''
        Macroturbulent velocity kernel.
        '''
        dlam[dlam == 0] = 1e-8
        if Zr != Zt:
            return np.array([(2*Ar*idlam/(np.sqrt(np.pi)*Zr**2) * quad(lambda u: np.exp(-1/u**2), 0, Zr/idlam)[0] +
                              2*At*idlam/(np.sqrt(np.pi)*Zt**2) * quad(lambda u: np.exp(-1/u**2), 0, Zt/idlam)[0])
                             for idlam in dlam])
        else:
            return np.array([(2*Ar*idlam/(np.sqrt(np.pi)*Zr**2) + 2*At*idlam/(np.sqrt(np.pi)*Zt**2)) *
                             quad(lambda u: np.exp(-1/u**2), 0, Zr/idlam)[0]
                             for idlam in dlam])

    def _broadening_macroturbulent(wave, flux, vmacro_rad, vmacro_tan=None,
                                   return_kernel=False):
        '''
        Apply macroturbulent broadening.
        The macroturbulent kernel is defined as in [Gray2005]:

        .. math::
            K_\mathrm{macro}(\Delta\lambda) = \frac{2A_R\Delta\lambda}{\sqrt{\pi}\zeta_R^2}\int_0^{\zeta_R/\Delta\lambda}e^{-1/u^2}du

             & + \frac{2A_T\Delta\lambda}{\sqrt{\pi}\zeta_T^2}\int_0^{\zeta_T/\Delta\lambda}e^{-1/u^2}du

        If :envvar:`vmacro_tan` is :envvar:`None`, then the value will be put equal
        to the radial component :envvar:`vmacro_rad`.

        Input
        -----
        :parameter wave: Wavelength of the spectrum
        :type wave: array
        :parameter flux: Flux of the spectrum
        :type flux: array
        :parameter vmacro_rad: macroturbulent broadening, radial component
        :type vmacro_rad: float
        :parameter vmacro_tan: macroturbulent broadening, tangential component
        :type vmacro_tan: float
        :parameter return_kernel: return kernel
        :type return_kernel: bool

        Output
        ------
        y_mac : broadened flux [, (wavelength, kernel)]
        '''

        if vmacro_tan is None:
            vmacro_tan = vmacro_rad

        if vmacro_rad == vmacro_tan == 0:
            return flux

        # Define central wavelength
        lambda0 = (wave[0] + wave[-1]) / 2.0
        vmac_rad = vmacro_rad/(299792458.*1e-3)*lambda0
        vmac_tan = vmacro_tan/(299792458.*1e-3)*lambda0

        # Make sure the wavelength range is equidistant before applying the
        # convolution
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
            raise ValueError(("Spectrum resolution too low for macroturbulent broadening"))
        elif n_kernel > n_wave:
            raise ValueError(("Spectrum range too narrow for macroturbulent broadening"))
        # Construct the broadening kernel
        wave_k = np.arange(n_kernel)*dwave
        wave_k -= wave_k[-1]/2.
        kernel = _vmacro_kernel(wave_k, 1.0, 1.0, vmac_rad, vmac_tan)
        kernel /= sum(kernel)

        flux_conv = fftconvolve(1-flux_, kernel, mode='same')
        # And interpolate the results back on to the original wavelength array,
        # taking care of even vs. odd-length kernels
        if n_kernel % 2 == 1:
            offset = 0.0
        else:
            offset = dwave / 2.0
        flux = np.interp(wave+offset, wave_, 1-flux_conv)

        # Return the results.
        if return_kernel:
            return flux, (wave_k, kernel)
        else:
            return flux

    # Instrumental broadening
    y_inst = instrumental_profile(x, y, resolution)
    # vsini broadening
    y_rot = vsini_broadening(x, y_inst, epsilon, vsini)
    # vmac broadening
    y_broad = _broadening_macroturbulent(x, y_rot, vmacro_rad=vmac,
                                         vmacro_tan=None, return_kernel=False)
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
    import itertools

    with open('summary.out', 'r') as f:
        f.readline()
        f.readline()
        start_wave, end_wave, step, flux_step = map(float, f.readline().split())
        lines = f.readlines()

    data = []
    for line in lines:
        line = line.replace('-',' ')
        line = line.replace('\n','').split(' ')
        line = filter(None, line)
        data.append(line)

    flux = list(itertools.chain(*data))
    flux = np.array(flux)
    flux = flux.astype(np.float)
    flux = 1.0 - flux

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
    wavelength : ndarray
      The wavelenth vector
    flux : ndarray
      The flux vector
    '''

    wavelength, flux = np.loadtxt(fname, skiprows=2, usecols=(0, 1), unpack=True)
    return wavelength, flux


def read_linelist(fname, intname='intervals.lst'):
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
        raise IOError('The interval list is not in the correct place!')

    lines = pd.read_csv('rawLinelist/%s' % fname, skiprows=1, comment='#', delimiter='\t', usecols=range(6),
    names=['wl', 'elem', 'excit', 'loggf', 'vdwaals', 'Do'],
    converters={'Do': lambda x : x.replace("nan"," "), 'vdwaals': lambda x : float(x)})
    lines.sort_values(by='wl', inplace=True)

    intervals = pd.read_csv('rawLinelist/%s' % intname, comment='#', names=['start', 'end'], delimiter='\t')
    ranges = intervals.values
    atomic = []
    N = []
    for i, ri in enumerate(intervals.values):
        a = lines[(lines.wl>ri[0]) & (lines.wl<ri[1])]
        a = a.as_matrix()
        atomic.append(a)
        N.append(len(lines[(lines.wl>ri[0]) & (lines.wl<ri[1])]))
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
    wave_obs : ndarray
      Observed wavelength
    wave_synth : ndarray
      Synthetic wavelength
    flux_synth : ndarray
      Synthetic flux

    Output
    ------
    int_flux : ndarray
      Interpolated synthetic flux
    '''

    from scipy.interpolate import InterpolatedUnivariateSpline

    sl = InterpolatedUnivariateSpline(wave_synth, flux_synth, k=1)
    int_flux = sl(wave_obs)
    return int_flux
