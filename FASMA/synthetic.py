# -*- coding: utf8 -*-
# My imports
from __future__ import division
import numpy as np
import os
import pandas as pd
from astropy.io import fits

def save_synth_spec(xobs=None, yobs=None, xs=None, ys=None, xf=None, yf=None, params=None, star=None, **options):
    '''Save synthetic spectrum of all intervals

    Input
    ----
    x : ndarray
      Wavelength
    y : ndarray
      Flux
    initial : list
      Set of parameters to name the new file, else it is named 'synthetic.spec'.

    Output
    -----
    fname fits file
    '''
    # Create header
    header = fits.Header()
    header['CRVAL1'] = xs[0]
    header['CDELT1'] = xs[1] - xs[0]

    fname = (
            str(star) 
            + '_'
            + str(params[0])
            + '_'
            + str(params[1])
            + '_'
            + str(params[2])
            + '_'
            + str(params[3])
            + '_'
            + str(params[4])
            + '_'
            + str(params[5])
            + '_'
            + str(options['resolution'])
            + '.spec'
        )        
    if xf is None: 
        xf = np.empty(len(xs))
        xf[:] = np.nan
        yf = np.empty(len(xs))
        yf[:] = np.nan
    if xobs is None: 
        xobs = np.empty(len(xs))
        xobs[:] = np.nan
        yobs = np.empty(len(xs))
        yobs[:] = np.nan
        
    tbhdu = fits.BinTableHDU.from_columns(
        [
            fits.Column(name='wavelength_synthetic_initial', format='D', array=xs),
            fits.Column(name='flux_synthetic_initial', format='D', array=ys),
            fits.Column(name='wavelength_synthetic_final', format='D', array=xf),
            fits.Column(name='flux_synthetic_final', format='D', array=yf),
            fits.Column(name='wavelength_observed', format='D', array=xobs),
            fits.Column(name='flux_observed', format='D', array=yobs),
        ],
        header=header,
    )
    tbhdu.writeto('results/%s' % fname, overwrite=True)
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
    epsilon : limb-darkening parameter

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
            y_inst = pyasl.instrBroadGaussFast(
                x, y, resolution, edgeHandling="firstlast", fullout=False, maxsig=None
            )
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
            y_rot = y
        else:
            y_rot = pyasl.rotBroad(x, y, epsilon, vsini, edgeHandling='firstlast')
        return y_rot

    def vmacro_kernel(dlam, Ar, At, Zr, Zt):
        '''
        Macroturbulent velocity kernel.
        '''
        dlam[dlam == 0] = 1e-8
        if Zr != Zt:
            return np.array(
                [
                    (
                        2
                        * Ar
                        * idlam
                        / (np.sqrt(np.pi) * Zr ** 2)
                        * quad(lambda u: np.exp(-1 / u ** 2), 0, Zr / idlam)[0]
                        + 2
                        * At
                        * idlam
                        / (np.sqrt(np.pi) * Zt ** 2)
                        * quad(lambda u: np.exp(-1 / u ** 2), 0, Zt / idlam)[0]
                    )
                    for idlam in dlam
                ]
            )
        else:
            return np.array(
                [
                    (
                        2 * Ar * idlam / (np.sqrt(np.pi) * Zr ** 2)
                        + 2 * At * idlam / (np.sqrt(np.pi) * Zt ** 2)
                    )
                    * quad(lambda u: np.exp(-1 / u ** 2), 0, Zr / idlam)[0]
                    for idlam in dlam
                ]
            )

    def vmac_broadening(wave, flux, vmacro_rad):
        '''
        Apply macroturbulent broadening.
        The macroturbulent kernel is defined as in [Gray2005].
        These functions are taken from iSpec (Blanco-Cuaresma et al. 2014)

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
        lambda0 = (wave[0] + wave[-1]) / 2.0
        vmac_rad = vmacro_rad / (299792458.0 * 1e-3) * lambda0
        vmac_tan = vmac_rad

        # Make sure the wavelength range is equidistant before applying the
        # convolution
        delta_wave = np.diff(wave).min()
        range_wave = wave.ptp()
        n_wave = int(range_wave / delta_wave) + 1
        wave_ = np.linspace(wave[0], wave[-1], n_wave)
        flux_ = np.interp(wave_, wave, flux)
        dwave = wave_[1] - wave_[0]
        n_kernel = int(5 * max(vmac_rad, vmac_tan) / dwave)
        if n_kernel % 2 == 0:
            n_kernel += 1
        # The kernel might be of too low resolution, or the the wavelength range
        # might be too narrow. In both cases, raise an appropriate error
        if n_kernel == 0:
            raise ValueError(
                ("Spectrum resolution too low for macroturbulent broadening")
            )
        elif n_kernel > n_wave:
            raise ValueError(
                ("Spectrum range too narrow for macroturbulent broadening")
            )
        # Construct the broadening kernel
        wave_k = np.arange(n_kernel) * dwave
        wave_k -= wave_k[-1] / 2.0
        kernel = vmacro_kernel(wave_k, 1.0, 1.0, vmac_rad, vmac_tan)
        kernel /= sum(kernel)

        flux_conv = fftconvolve(1 - flux_, kernel, mode='same')
        # And interpolate the results back on to the original wavelength array,
        # taking care of even vs. odd-length kernels
        if n_kernel % 2 == 1:
            offset = 0.0
        else:
            offset = dwave / 2.0
        flux = np.interp(wave + offset, wave_, 1 - flux_conv)
        return flux

    # vmac broadening
    y_mac = vmac_broadening(x, y, vmacro_rad=vmac)
    # vsini broadening
    y_rot = vsini_broadening(x, y_mac, epsilon, vsini)
    # Instrumental broadening
    y_inst = instrumental_profile(x, y_rot, resolution)
    return x, y_inst

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
        start_wave, end_wave, step, flux_step = list(map(float, f.readline().split()))
        lines = f.readlines()

    data = []
    for line in lines:
        line = line.replace('-', ' ')
        line = line.replace('\n', '').split(' ')
        line = filter(None, line)
        data.append(line)

    flux = list(itertools.chain(*data))
    flux = np.array(flux)
    flux = flux.astype(float)
    flux = 1.0 - flux

    w0, dw, n = float(start_wave), float(step), len(flux)
    w = w0 + dw * n
    wavelength = np.linspace(w0, w, n, endpoint=False)
    return wavelength, flux

def read_linelist(fname, intname='intervals.lst'):
    '''Read the line list (atomic data) and the file which includes the ranges
    where the synthesis will happen.

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

    lines = pd.read_csv(
        fname,
        skiprows=1,
        comment='#',
        delimiter='\t',
        usecols=range(6),
        names=['wl', 'elem', 'excit', 'loggf', 'vdwaals', 'Do'],
        converters={
            'Do': lambda x: x.replace("nan", " "),
            'vdwaals': lambda x: float(x),
        },
    )
    lines.sort_values(by='wl', inplace=True)

    intervals = pd.read_csv(
        intname, comment='#', names=['start', 'end'], delimiter='\t'
    )
    ranges = intervals.values
    atomic = []
    N = []
    for i, ri in enumerate(intervals.values):
        a = lines[(lines.wl > ri[0]) & (lines.wl < ri[1])]
        atomic.append(a.values)
        N.append(len(a))
    N = sum(N)
    atomic = np.vstack(atomic)
    print('Linelist contains %s lines in %s intervals' % (N, len(ranges)))

    # Create line list for MOOG
    fmt = ['%9.3f', '%10.1f', '%9.2f', '%9.3f', '%9.3f', '%7.4s']
    header = 'Wavelength     ele       EP      loggf   vdwaals   Do'
    np.savetxt('linelist.moog', atomic, fmt=fmt, header=header)
    return ranges, atomic

def read_linelist_elem(fname, element=None, intname='intervals_elements.lst'):
    '''Read the line list (atomic data) and the file which includes the ranges
    where the synthesis will happen for the element abundances.

    Input
    -----
    fname : str
      File that contains the linelist
    element : str
      The element to be searched in the line list
    intname : str
      File that contains the central line where -+2.0 \AA{} are added to create
      the interval around each line.

    Output
    ------
    ranges : wavelength ranges of the linelist
    atomic : atomic data
    '''

    if not os.path.isfile(intname):
        raise IOError('The interval list is not in the correct place!')
    print('Line list:', fname)
    print('Intervals list:', intname)
    lines = pd.read_csv(
        fname,
        skiprows=1,
        comment='#',
        delimiter='\t',
        usecols=range(6),
        names=['wl', 'elem', 'excit', 'loggf', 'vdwaals', 'Do'],
    )
    lines.sort_values(by='wl', inplace=True)
    intervals = pd.read_csv(
        intname,
        comment='#',
        usecols=(0, 1),
        names=['El', 'wave'],
        delimiter='\t',
    )

    try:
        intervals['El'] = intervals['El'].map(lambda x: x.strip().strip("I"))
    except AttributeError:
        print('The format of the line list is not correct.')
    intervals = intervals[intervals['El'] == element]
    intervals['wave'] = intervals['wave'].apply(pd.to_numeric)
    intervals.sort_values(by='wave', inplace=True)

    ranges = []
    merged = [[intervals.wave.iloc[0] - 2.0, intervals.wave.iloc[0] + 2.0]]
    for i, ri in enumerate(intervals.wave):
        r1 = float(ri) - 2.0
        r2 = float(ri) + 2.0
        previous = merged[-1]
        if r1 <= previous[1]:
            previous[1] = max(previous[1], r2)
        else:
            merged.append([r1, r2])
        ranges.append([r1, r2])

    atomic = []
    N = []
    for ri in merged[:]:
        a = lines[(lines.wl > ri[0]) & (lines.wl < ri[1])]
        atomic.append(a)
        N.append(len(a))

    atomic = pd.concat(atomic)
    atomic = atomic.fillna(' ')
    atomic = atomic.drop_duplicates()
    atomic.sort_values(by='wl', inplace=True)
    N = sum(N)
    merged = [[x[0]-3., x[1]+3.] for x in merged]
    # Create line list for MOOG
    print('Linelist contains %s lines in %s interval(s).' % (N, len(merged)))
    fmt = ['%9.3f', '%10.1f', '%9.2f', '%9.3f', '%9.3f', '%7.4s']
    header = 'Wavelength     ele       EP      loggf   vdwaals   Do'
    np.savetxt('linelist.moog', atomic.values, fmt=fmt, header=header)
    return merged, atomic
