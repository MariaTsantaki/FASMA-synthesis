#!/usr/bin/python
from __future__ import division
from periodictable import elements
from scipy.interpolate import griddata
import numpy as np
import gzip
from .solar_abundance import solar
from .utils import GetModels

def read_model(fname):
    '''Read the KURUCZ model atmosphere.

    Input
    -----
    fname : str
      The gz file of the model atmosphere.

    Output
    ------
    model : ndarray
      The correct atmosphere, the columns and tauross in a tuple
    '''

    f = gzip.open(fname, compresslevel=1)
    data = f.readlines()
    model = np.loadtxt(data[23:-2])
    return model

def solar_abundance(element):
    '''For a given atomic number, return solar abundance from Asplund et al. 2009.

    Input
    -----
    element : str
      The name of an element

    Output
    ------
    index : int
      The atomic number of the element
    abundance : float
      The solar abundance of the atom in dex
    '''
    abundance = solar.get(element, None)
    if abundance is None:
        return None, None
    index = elements.symbol(element).number
    return index, abundance

def interpolator_kurucz(params, atmtype='apogee_kurucz'):
    '''Interpolation for Kurucz model atmospheres.
    Input
    -----
    params : ndarray
      Stellar parameters for the Interpolation.

    Output
    ------
    newatm : ndarray
      The interpolated atmosphere, the columns and tauross in a tuple
    '''

    m = GetModels(params[0], params[1], params[2], atmtype=atmtype)
    mdict = m.getmodels()
    if mdict is False:
        return False

    params[0] = mdict['teff'][0]
    params[1] = mdict['logg'][0]
    params[2] = mdict['feh'][0]
    mnames = mdict['models']
    teff = mdict['teff']
    logg = mdict['logg']
    feh = mdict['feh']
    # Making the surrounding grid points
    gridpoints = []
    for temp in teff[1]:
        for grav in logg[1]:
            for metal in feh[1]:
                gridpoints.append((temp, grav, metal))
    gridpoints = np.asarray(gridpoints)
    # Define the points to obtain at the end
    teff = teff[0]
    logg = logg[0]
    feh = feh[0]
    options = {'method': 'linear', 'rescale': True}

    # Reading the models
    models = []
    for mname in mnames:
        tatm = read_model(mname)
        models.append(tatm)

    layers = range(min([model.shape[0] for model in models]))
    columns = range(6)
    newatm = np.zeros((len(layers), len(columns)))
    for layer in layers:
        for column in columns:
            tlayer = [model[layer, column] for model in models]
            newatm[layer, column] = griddata(
                gridpoints, tlayer, (teff, logg, feh), **options
            )
    vt_array = np.zeros(len(layers)) + params[-1] * 1e5
    newatm = np.hstack((newatm, vt_array[:, np.newaxis]))
    return newatm

def interpolator_marcs(params, atmtype='marcs'):
    '''Interpolation for Kurucz model atmospheres.
    Input
    -----
    params : ndarray
      Stellar parameters for the Interpolation.

    Output
    ------
    newatm : ndarray
      The interpolated atmosphere, the columns and tauross in a tuple
    '''

    m = GetModels(params[0], params[1], params[2], atmtype=atmtype)
    mdict = m.getmodels()
    if mdict is False:
        return False

    params[0] = mdict['teff'][0]
    params[1] = mdict['logg'][0]
    params[2] = mdict['feh'][0]
    mnames = mdict['models']
    teff = mdict['teff']
    logg = mdict['logg']
    feh = mdict['feh']
    # Making the surrounding grid points
    gridpoints = []
    for temp in teff[1]:
        for grav in logg[1]:
            for metal in feh[1]:
                gridpoints.append((temp, grav, metal))
    gridpoints = np.asarray(gridpoints)
    # Define the points to obtain at the end
    teff = teff[0]
    logg = logg[0]
    feh = feh[0]
    options = {'method': 'linear', 'rescale': True}

    # Reading the models
    models = []
    for mname in mnames:
        tatm = read_model(mname)
        models.append(tatm)

    layers = range(min([model.shape[0] for model in models]))
    columns = range(6)
    newatm = np.zeros((len(layers), len(columns)))
    for layer in layers:
        for column in columns:
            tlayer = [model[layer, column] for model in models]
            newatm[layer, column] = griddata(
                gridpoints, tlayer, (teff, logg, feh), **options
            )
    vt_array = np.zeros(len(layers)) + params[-1] * 1e5
    newatm = np.hstack((newatm, vt_array[:, np.newaxis]))
    return newatm

def interpolator(params, abund=0.0, elem=False, save=True, atmtype='apogee_kurucz', result=None):
    '''Function to connect all. For a given set of params, return a model atmosphere.

    Input
    -----
    params : list of length 3
      Teff, logg, [Fe/H].
    abund : float
      abundance of a given element to added in the atmosphere.
    element : float
      any element to change abundance in the atmosphere.
    save : bool
      Whether the new atmosphere should be saved. Default is True.
    atmtype : str
      The atmosphere models being used. Default is apogee_kurucz.
    result : bool
      return the new atmosphere. Default is False.

    Output
    ------
    newatm : ndarray
      New interpolated atmosphere.
    '''

    params = list(params)
    if atmtype == 'marcs':
        # newatm = interpolator_marcs(params, fesun=7.41, microlim=3.0)
        newatm = interpolator_marcs(params, atmtype=atmtype)
    elif atmtype == 'apogee_kurucz':
        newatm = interpolator_kurucz(params, atmtype=atmtype)
    else:
        raise NameError('Could not find model atmosphere: ', atmtype)

    if newatm is False:
        raise NameError('Could not find model atmosphere: ', atmtype)
    if save:
        save_model(newatm, params, abund=abund, elem=elem, type=atmtype)
    if result:
        return newatm, params

def save_model(model, params, abund=0.0, elem=False, type='apogee_kurucz', fout='out.atm'):
    '''Save the model atmosphere in the right format.

    Input
    -----
    model : ndarray
      The interpolated model atmosphere.
    params : list
      Teff, logg, [M/H], vt of the interpolated atmosphere.
    abund : float
      abundance of a given element to added in the atmosphere.
    element : float
      any element to change abundance in the atmosphere.
    type : str
      Type of atmospheric parameters. Default is apogee_kurucz
    fout : str
      Name of the saved atmosphere. Default is out.atm

    Output
    ------
    Saved atmospheric model in file.
    '''

    teff, logg, feh, vt = params
    if type in ['apogee_kurucz', 'marcs']:
        # The name in the header shows the format of the model not the type.
        header = (
            'KURUCZ\n'
            'Teff= %i   log g= %.2f\n'
            'NTAU        %i' % (teff, logg, model.shape[0])
        )

    if elem:
        # Get atomic number and solar abundance, only one element per time.
        num, solabund = solar_abundance(elem)
        num_fe, solabund_fe = solar_abundance('Fe')
        footer = (
            '    %.3e\n'
            'NATOMS     2  %.2f\n'
            '      26.0    %.2f\n'
            '      %.1f    %.2f\n'
            'NMOL      19\n'
            '      606.0    106.0    607.0    608.0    107.0    108.0    112.0  707.0\n'
            '       708.0    808.0     12.1  60808.0  10108.0    101.0     6.1    7.1\n'
            '         8.1    822.0     22.1'
            % (vt * 1e5, feh, solabund_fe + feh, num, solabund + abund)
        )

    else:
        num_fe, solabund_fe = solar_abundance('Fe')
        footer = (
            '    %.3e\n'
            'NATOMS     1  %.2f\n'
            '      26.0    %.2f\n'
            'NMOL      19\n'
            '      606.0    106.0    607.0    608.0    107.0    108.0    112.0  707.0\n'
            '       708.0    808.0     12.1  60808.0  10108.0    101.0     6.1    7.1\n'
            '         8.1    822.0     22.1' % (vt * 1e5, feh, solabund_fe + feh)
        )

    _fmt = ('%15.8E', '%8.1f', '%.3E', '%.3E', '%.3E', '%.3E', '%.3E')
    while model.shape[1] < len(_fmt):
        model = np.column_stack((model, np.zeros_like(model[:, 0])))
    np.savetxt(
        fout, model, header=header, footer=footer, comments='', delimiter=' ', fmt=_fmt
    )
