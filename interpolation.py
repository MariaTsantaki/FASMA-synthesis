#!/usr/bin/python
from __future__ import division
import numpy as np
import gzip
from scipy.interpolate import griddata
from utils import GetModels

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

def solar_abundance(elem):
    '''For a given atomic number, return solar abundance from Asplund et al. 2009.

    Input
    -----
    atom : int
      The atomic number

    Output
    ------
    abundance : float
      The solar abundance of the atom in dex
    '''

    element = ['H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne' ,
    'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', 'Sc', 'Ti',
    'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se',
    'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd',
    'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce',
    'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb',
    'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U',  'Np', 'Pu',
    'Am', 'Cm', 'Bk', 'Cf', 'Es'] #periodic table

    solar = [12.00, 10.93, 0.96, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93,
             6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40, 5.03, 6.34,
             3.15, 4.95, 3.93, 5.64, 5.43, 7.47, 4.99, 6.22, 4.19, 4.56,
             3.04, 3.65, 2.30, 3.34, 2.54, 3.25, 2.52, 2.87, 2.21, 2.58,
             1.46, 1.88, -5.00, 1.75, 0.91, 1.57, 0.94, 1.71, 0.80, 2.04,
             1.01, 2.18, 1.55, 2.24, 1.08, 2.18, 1.10, 1.58, 0.72, 1.42,
             -5.00, 0.96, 0.52, 1.07, 0.30, 1.10, 0.48, 0.92, 0.10, 0.84,
             0.10, 0.85, -0.12, 0.85, 0.26, 1.40, 1.38, 1.62, 0.92, 1.17,
             0.90, 1.75, 0.65, -5.00, -5.00, -5.00, -5.00, -5.00, -5.00,
             0.02, -5.00, -0.54, -5.00, -5.00, -5.00] #solar abundances

    if elem in element:
        ind = element.index(str(elem))
    else:
        print('Element does not exist in the periodic table.')
        ind, solar[ind] = None, None
    return ind+1, solar[ind]

def interpolator_kurucz(params, atmtype='kurucz95'):
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
    feh  = feh[0]
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
            newatm[layer, column] = griddata(gridpoints, tlayer, (teff, logg, feh), **options)
    vt_array = np.zeros(len(layers))+params[-1]*1e5
    newatm = np.hstack((newatm, vt_array[:, np.newaxis]))
    return newatm

def interpolator_marcs(params, fesun=7.47, microlim=3.0):
    '''Interpolation for marcs models. The function is taken from STEPAR
    (Tabernero et al. 2019) to deal the gaps in the grid.
    Input
    -----
    params : ndarray
      Stellar parameters for the Interpolation.

    Output
    ------
    newatm : ndarray
      The interpolated atmosphere, the columns and tauross in a tuple
    '''

    import _pickle as pic
    import os

    fname = "models/marcs/MARCS1M.bin"
    if not os.path.isfile(fname):
        return False

    gridMODS = open(fname,"rb")
    tmod     = pic.load(gridMODS)
    gmod     = pic.load(gridMODS)
    mmod     = pic.load(gridMODS)
    ltaumod  = pic.load(gridMODS)
    Temod    = pic.load(gridMODS)
    lpgmod   = pic.load(gridMODS)
    lpemod   = pic.load(gridMODS)
    rhoxmod  = pic.load(gridMODS)
    kmod     = pic.load(gridMODS)
    gridMODS.close()

    x = list(params)
    Teff  = np.round(x[0],0)
    logg  = np.round(x[1],2)
    metal = np.round(x[2],2)
    micro = np.round(x[3],2)
    sel = np.where((np.abs(tmod-Teff) <= 251.) & (np.abs(mmod-metal) <= 0.251) & (np.abs(gmod-logg) <= 0.5))
    lT = np.shape(tmod[sel])[0]
    lm = np.shape(mmod[sel])[0]
    lg = np.shape(gmod[sel])[0]
    if micro > microlim:
        micro = microlim
    if lT > 1 and lm > 1 and lg > 1 and micro <= microlim and micro >= 0.:
        dT = max(tmod[sel])-min(tmod[sel])
        dm = max(mmod[sel])-min(mmod[sel])
        dg = max(gmod[sel])-min(gmod[sel])
        testt = (min(tmod[sel]) <= Teff) and (max(tmod[sel]) >= Teff) and dT > 0.
        testm = (min(mmod[sel]) <= metal) and (max(mmod[sel]) >= metal) and dm > 0.
        testg = (min(gmod[sel]) <= logg) and (max(gmod[sel]) >= logg) and dg > 0.
        if testt and testm and testg:
            testint = True
        else:
            testint = False
    else:
        testint = False
    if testint:
        teint = griddata((tmod[sel], gmod[sel], mmod[sel]), Temod[sel], (Teff, logg, metal), method='linear', fill_value=np.nan, rescale=True)
        if np.isfinite(teint[0]):
            lpgint  = griddata((tmod[sel], gmod[sel], mmod[sel]), lpgmod[sel], (Teff, logg, metal),  method='linear', fill_value=np.nan, rescale=True)
            lpeint  = griddata((tmod[sel], gmod[sel], mmod[sel]), lpemod[sel], (Teff, logg, metal),  method='linear', fill_value=np.nan, rescale=True)
            rhoxint = griddata((tmod[sel], gmod[sel], mmod[sel]), rhoxmod[sel], (Teff, logg, metal), method='linear', fill_value=np.nan, rescale=True)
            kint    = griddata((tmod[sel], gmod[sel], mmod[sel]), kmod[sel], (Teff, logg, metal),    method='linear', fill_value=np.nan, rescale=True)
            nlayers = len(kint)
            newatm = []
            for i in range(nlayers):
                newatm.append([rhoxint[i], teint[i], 10.**lpgint[i], 10.**lpeint[i], kint[i]])
            newatm = np.array(newatm)
            vt_array = np.zeros(nlayers)+params[-1]*1e5
            newatm = np.hstack((newatm, vt_array[:, np.newaxis]))
            return newatm
        else:
            return False
    else:
        return False

def interpolator(params, abund=0.0, elem=False, save=True, atmtype='kurucz95', result=None):
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
      The atmosphere models being used. Default is Kurucz95.
    result : bool
      return the new atmosphere. Default is False.

    Output
    ------
    newatm : ndarray
      New interpolated atmosphere.
    '''

    params = list(params)
    if atmtype == 'marcs':
        newatm = interpolator_marcs(params, fesun=7.47, microlim=3.0)
        if newatm is False:
            raise NameError('Could not find %s models' % atmtype)
    elif atmtype == 'kurucz95':
        newatm = interpolator_kurucz(params, atmtype=atmtype)
        if newatm is False:
            raise NameError('Could not find %s models' % atmtype)
    elif atmtype == 'apogee_kurucz':
        newatm = interpolator_kurucz(params, atmtype=atmtype)
        if newatm is False:
            raise NameError('Could not find %s models' % atmtype)
    else:
        raise NameError('Could not find %s models' % atmtype)

    if save:
            save_model(newatm, params, abund=abund, elem=elem, type=atmtype)
    if result:
        return newatm, params

def save_model(model, params, abund=0.0, elem=False, type='kurucz95', fout='out.atm'):
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
      Type of atmospheric parameters. Default is Kurucz95
    fout : str
      Name of the saved atmosphere. Default is out.atm

    Output
    ------
    Saved atmospheric model in file.
    '''

    teff, logg, feh, vt = params
    if type in ['kurucz95', 'apogee_kurucz', 'marcs']:
        # The name in the header shows the format of the model not the type.
        header = 'KURUCZ\n'\
                 'Teff= %i   log g= %.2f\n'\
                 'NTAU        %i' % (teff, logg, model.shape[0])
    else:
        raise NameError('Could not find %s models' % type)

    if elem:
        # Get atomic number and solar abundance, only one element per time.
        num, solabund = solar_abundance(elem)
        footer = '    %.3e\n'\
        'NATOMS     2  %.2f\n'\
        '      26.0    %.2f\n'\
        '      %.1f    %.2f\n'\
        'NMOL      19\n'\
        '      606.0    106.0    607.0    608.0    107.0    108.0    112.0  707.0\n'\
        '       708.0    808.0     12.1  60808.0  10108.0    101.0     6.1    7.1\n'\
        '         8.1    822.0     22.1' % (vt*1e5, feh, 7.47+feh, num, solabund+abund)

    else:
        footer = '    %.3e\n'\
        'NATOMS     1  %.2f\n'\
        '      26.0    %.2f\n'\
        'NMOL      19\n'\
        '      606.0    106.0    607.0    608.0    107.0    108.0    112.0  707.0\n'\
        '       708.0    808.0     12.1  60808.0  10108.0    101.0     6.1    7.1\n'\
        '         8.1    822.0     22.1' % (vt*1e5, feh, 7.47+feh)

    _fmt = ('%15.8E', '%8.1f', '%.3E', '%.3E', '%.3E', '%.3E', '%.3E')
    while model.shape[1] < len(_fmt):
        model = np.column_stack((model, np.zeros_like(model[:, 0])))
    np.savetxt(fout, model, header=header, footer=footer, comments='', delimiter=' ', fmt=_fmt)

if __name__ == '__main__':
    import argparse
    args = argparse.ArgumentParser(description='Get a model atmosphere.')
    args.add_argument('teff',   type=int,       help='Effective temperature')
    args.add_argument('logg',   type=float,     help='Surface gravity')
    args.add_argument('feh',    type=float,     help='Metallicity, [Fe/H]')
    args.add_argument('vt',     type=float,     help='Microturbulence')
    args.add_argument('-elem',  type=str,       help='Element', default=False)
    args.add_argument('-abund', type=float,     help='Abundance', default=0.0)
    args.add_argument('-o',     '--out',        help='Output atmosphere', default='out.atm')
    args.add_argument('-a',     '--atmosphere', help='Model atmosphere', choices=['kurucz95', 'apogee_kurucz', 'marcs'], default='kurucz95')
    args = args.parse_args()

    params = [args.teff, args.logg, args.feh, args.vt]
    atmosphere, p = interpolator(params, elem=args.elem, abund=args.abund, save=False, atmtype=args.atmosphere, result=True)
    save_model(atmosphere, params, elem=args.elem, abund=args.abund, type=args.atmosphere, fout=args.out)
    print('Atmosphere model sucessfully saved in: %s' % args.out)
