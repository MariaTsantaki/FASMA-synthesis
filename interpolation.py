#!/usr/bin/python
from __future__ import division
import numpy as np
import gzip
from scipy.interpolate import griddata
from utils import GetModels


def read_model(fname):
    '''Read the model atmosphere

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


def interpolator(params, save=True, atmtype='marcs', result=None):
    '''This is a new approach based on a scipy interpolator.
    Resembles the original interpolator we used but with a change

    Input
    -----
    params : list of length 3
      Teff, logg, [Fe/H] desired.
    save : bool
      Wether the new atmosphere should be saved. Default is True.
    atmtype : str
      The atmosphere models being used. Default is MARCS.

    Output
    ------
    newatm : ndarray
      New interpolated atmosphere.
    '''

    params = list(params)
    m = GetModels(params[0], params[1], params[2], atmtype=atmtype)
    mdict = m.getmodels()
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
            newatm[layer, column] = griddata(gridpoints, tlayer, (teff, logg, feh), **options)
    vt_array = np.zeros(len(layers))+params[-1]*1e5
    newatm = np.hstack((newatm, vt_array[:, np.newaxis]))
    if save:
        save_model(newatm, params, type=atmtype)
    if result:
        return newatm, params


def save_model(model, params, type='marcs', fout='out.atm'):
    '''Save the model atmosphere in the right format

    Input
    -----
    model : ndarray
      The interpolated model atmosphere.
    params : list
      Teff, logg, [Fe/H], vt of the interpolated atmosphere.
    type : str
      Type of atmospheric parameters. Default is MARCS
    fout : str
      Name of the saved atmosphere. Default is out.atm

    Output
    ------
    Atmospheric model.
    '''
    teff, logg, feh, vt = params
    if type in ['apogee_kurucz', 'marcs']:
        header = 'KURUCZ\n'\
                 'Teff= %i   log g= %.2f\n'\
                 'NTAU        %i' % (teff, logg, model.shape[0])
    else:
        raise NameError('Could not find %s models' % type)

    footer = '    %.3e\n'\
             'NATOMS     1  %.2f\n'\
             '      26.0   %.2f\n'\
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
    args.add_argument('teff', type=int, help='Effective temperature')
    args.add_argument('logg', type=float, help='Surface gravity')
    args.add_argument('feh', type=float, help='Metallicity, [Fe/H]')
    args.add_argument('vt', type=float, help='Microturbulence')
    args.add_argument('-o', '--out', help='Output atmosphere', default='out.atm')
    args.add_argument('-a', '--atmosphere', help='Model atmosphere', choices=['apogee_kurucz', 'marcs'], default='marcs')

    args = args.parse_args()

    params = [args.teff, args.logg, args.feh, args.vt]

    atmosphere, p = interpolator(params, save=False, atmtype=args.atmosphere, result=True)
    save_model(atmosphere, params, type=args.atmosphere, fout=args.out)

    print 'Atmosphere model sucessfully saved in: %s' % args.out
