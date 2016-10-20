#!/usr/bin/python

from __future__ import division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import argparse

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
try:
    import seaborn as sns
    sns.set_style('dark')
    sns.set_context('talk', font_scale=1.2)
    colorSB = sns.color_palette()
except ImportError:
    print 'Please install seaborn: pip install seaborn'


def massTorres(teff, erteff, logg, erlogg, feh, erfeh):
    """Calculate a mass using the Torres calibration"""
    ntrials = 100
    randomteff = teff + erteff * np.random.randn(ntrials)
    randomlogg = logg + erlogg * np.random.randn(ntrials)
    randomfeh = feh + erfeh * np.random.randn(ntrials)

    # Parameters for the Torres calibration:
    a1, a2, a3 = 1.5689, 1.3787, 0.4243
    a4, a5, a6 = 1.139, -0.1425, 0.01969
    a7 = 0.1010

    logM = np.zeros(ntrials)
    for i in xrange(ntrials):
        X = np.log10(randomteff[i]) - 4.1
        logM[i] = a1 + a2*X + a3*X**2 + a4*X**3 + a5*randomlogg[i]**2 + a6*randomlogg[i]**3 + a7*randomfeh[i]

    meanMasslog = np.mean(logM)
    sigMasslog = np.sqrt(np.sum(logM-meanMasslog)**2)/(ntrials-1)
    sigMasslogTot = np.sqrt(0.027**2 + sigMasslog**2)

    meanMass = 10**meanMasslog
    sigMass = 10**(meanMasslog + sigMasslogTot) - meanMass
    return meanMass, sigMass


def radTorres(teff, erteff, logg, erlogg, feh, erfeh):
    ntrials = 100
    randomteff = teff + erteff*np.random.randn(ntrials)
    randomlogg = logg + erlogg*np.random.randn(ntrials)
    randomfeh = feh + erfeh*np.random.randn(ntrials)

    # Parameters for the Torres calibration:
    b1, b2, b3 = 2.4427, 0.6679, 0.1771
    b4, b5, b6 = 0.705, -0.21415, 0.02306
    b7 = 0.04173

    logR = np.zeros(ntrials)
    for i in xrange(ntrials):
        X = np.log10(randomteff[i]) - 4.1
        logR[i] = b1 + b2*X + b3*X**2 + b4*X**3 + b5*randomlogg[i]**2 + b6*randomlogg[i]**3 + b7*randomfeh[i]

    meanRadlog = np.mean(logR)
    sigRadlog = np.sqrt(np.sum((logR-meanRadlog)**2))/(ntrials-1)
    sigRadlogTot = np.sqrt(0.014**2 + sigRadlog**2)

    meanRad = 10**meanRadlog
    sigRad = 10**(meanRadlog + sigRadlogTot) - meanRad
    return meanRad, sigRad


def _parser():
    parser = argparse.ArgumentParser(description='Preprocess the results')
    p = ['teff', 'tefferr', 'logg', 'loggerr', 'feh', 'feherr', 'vt', 'vterr']
    p += ['lum', 'mass', 'masserr', 'radius', 'radiuserr', 'age']
    parser.add_argument('x', choices=p)
    parser.add_argument('y', choices=p)
    parser.add_argument('-z', help='Color scale', choices=p, default=None)
    parser.add_argument('-i', '--input', help='File name of result file', default='results.csv')
    parser.add_argument('-c', '--convergence', help='Only plot converged results', default=True, action='store_false')
    parser.add_argument('-ix', help='Inverse x axis', default=False, action='store_true')
    parser.add_argument('-iy', help='Inverse y axis', default=False, action='store_true')
    parser.add_argument('-iz', help='Inverse z axis', default=False, action='store_true')
    parser.add_argument('-lx', help='Logarithmic x axis', default=False, action='store_true')
    parser.add_argument('-ly', help='Logarithmic y axis', default=False, action='store_true')
    parser.add_argument('-s', help='Place Solar values in the plot', default=False, action='store_true')
    parser.add_argument('-l', help='Fit a linear regression', default=False, action='store_true')
    args = parser.parse_args()
    return args


if __name__ == '__main__':

    args = _parser()

    df = pd.read_csv(args.input, delimiter=r'\s+', comment='#')
    df = df[(df.convergence) | (~df.convergence)]  # Remove blank lines and comments
    df.teff = pd.to_numeric(df.teff, errors='coarse')
    df.tefferr = pd.to_numeric(df.tefferr, errors='coarse')
    df.logg = pd.to_numeric(df.logg, errors='coarse')
    df.loggerr = pd.to_numeric(df.loggerr, errors='coarse')
    df.feh = pd.to_numeric(df.feh, errors='coarse')
    df.feherr = pd.to_numeric(df.feherr, errors='coarse')
    df.vt = pd.to_numeric(df.vt, errors='coarse')
    df.vterr = pd.to_numeric(df.vterr, errors='coarse')

    m_ = ['mass', 'masserr', 'lum', 'radius', 'radiuserr', 'age']
    if (args.x in m_) or (args.y in m_) or (args.z in m_):
        params = zip(df.teff, df.tefferr, df.logg, df.loggerr, df.feh, df.feherr)
        m = [massTorres(t, et, l, el, f, ef) for t, et, l, el, f, ef in params]
        r = [radTorres(t, et, l, el, f, ef) for t, et, l, el, f, ef in params]
        df['mass'] = pd.Series(np.asarray(m)[:, 0])
        df['masserr'] = pd.Series(np.asarray(m)[:, 1])
        df['radius'] = pd.Series(np.asarray(r)[:, 0])
        df['radiuserr'] = pd.Series(np.asarray(r)[:, 1])
        df['lum'] = (df.teff/5777)**4 * df.radius**2

    if (args.x == 'age') or (args.y == 'age') or (args.z == 'age'):
        from isochrones.dartmouth import Dartmouth_Isochrone
        dar = Dartmouth_Isochrone()
        age = np.zeros(df.shape[0])
        for i, (mass, feh) in enumerate(df[['mass', 'feh']].values):
            age[i] = np.mean(dar.agerange(mass, feh))
        age = (10**age)/1e9
        df['age'] = pd.Series(age)

    df1 = df[df.convergence]
    df2 = df[~df.convergence]

    # Plot the results
    plt.figure()
    if args.z:
        if args.iz:
            z = 1/df1[args.z].values
        else:
            z = df1[args.z].values
        color = df1[args.z].values
        u = z[~np.isnan(z)]
        size = (z-u.min())/(u.max()-u.min())*100
        size[np.argmin(size)] = 10  # Be sure to show the "smallest" point
        plt.scatter(df1[args.x], df1[args.y], c=color, s=size, cmap=cm.seismic, label='Converged')
    else:
        plt.scatter(df1[args.x], df1[args.y], c=colorSB[0], s=40, label='Converged')
    if not args.convergence:
        if args.z:
            plt.scatter(df2[args.x], df2[args.y], c=df2[args.z].values, cmap=cm.seismic, s=55, marker='x', label='Not converged')
        else:
            plt.scatter(df2[args.x], df2[args.y], c=color[2], s=9, marker='d', label='Not converged')
        plt.legend(loc='best', frameon=False)

    if args.l:
        p = np.polyfit(df1[args.x], df1[args.y], deg=1)
        print '  y=%.3f*x+%.3f' % (p[0], p[1])
        yfit = np.poly1d(p)(df1[args.x])
        plt.plot(df1[args.x], yfit, '-k')

    labels = {'teff': r'$T_\mathrm{eff}$ [K]',
              'tefferr': r'$\sigma T_\mathrm{eff}$ [K]',
              'logg': r'$\log(g)$ [cgs]',
              'loggerr': r'$\sigma \log(g)$ [cgs]',
              'feh': '[Fe/H]',
              'feherr': r'$\sigma$ [Fe/H]',
              'vt': r'$\xi_\mathrm{micro}$ [km/s]',
              'vterr': r'$\sigma\xi_\mathrm{micro}$ [km/s]',
              'lum': r'$L_\odot$',
              'mass': r'$M_\odot$',
              'masserr': r'$\sigma M_\odot$',
              'radius': r'$R_\odot$',
              'radiuserr': r'$\sigma R_\odot$',
              'age': r'Age $[Gyr]$'}

    plt.xlabel(labels[args.x])
    plt.ylabel(labels[args.y])
    if args.z:
        cbar = plt.colorbar()
        cbar.set_label(labels[args.z])
    if args.s:
        sun = {'teff': 5777,
               'tefferr': 1,
               'logg': 4.44,
               'loggerr': 0.01,
               'feh': 0.00,
               'feherr': 0.01,
               'vt': 1.00,
               'vterr': 0.01,
               'lum': 1,
               'mass': 1,
               'masserr': 0.01,
               'radius': 1,
               'radiuserr': 0.01,
               'age': 4.567}
        plt.plot(sun[args.x], sun[args.y], color=colorSB[4], marker='*', ms=20, alpha=0.8)
    if args.ix:
        plt.xlim(plt.xlim()[::-1])
    if args.iy:
        plt.ylim(plt.ylim()[::-1])

    if args.lx:
        plt.xscale('log')
    if args.ly:
        plt.yscale('log')

    plt.grid(True)
    plt.tight_layout()
    plt.show()
