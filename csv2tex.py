#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import pandas as pd
import argparse


def paramError(s1, s2):
    """Add together s1 and s2 with the right delimiter \pm and surroundings $$

    Inputs
    ------
    s1 : pandas Series
      A pandas Series with main parameter, e.g. Teff
    s2 : pandas Series
      A pandas Series with error, e.g. Tefferr

    Output
    ------
    sout : pandas Series
      A new series with format like this: $5777 \\pm 42$
    """
    sout = '$' + s1.astype('str') + ' \pm ' + s2.astype('str') + '$'
    return sout


def _parser():
    parser = argparse.ArgumentParser(description='Convert the result file to a TeX table')
    parser.add_argument('-i', help='Input file', default='results.csv')
    parser.add_argument('-o', help='Output TeX table', default=None)
    return parser.parse_args()


if __name__ == '__main__':

    args = _parser()
    df = pd.read_csv(args.i, delimiter=r'\s+')
    if args.o is None:
        output = args.i.rpartition('.')[0] + '.tex'
    else:
        output = args.o

    junk = ['fixteff', 'fixlogg', 'fixfeh', 'fixvt', 'outlier', 'weights',
            'model', 'refine', 'EPcrit', 'RWcrit', 'ABdiffcrit']
    df.drop(junk, inplace=True, axis=1)

    df['teff'] = paramError(df['teff'], df['tefferr'])
    df['logg'] = paramError(df['logg'], df['loggerr'])
    df['feh'] = paramError(df['feh'], df['feherr'])
    df['vt'] = paramError(df['vt'], df['vterr'])
    df['loggastero'] = paramError(df['loggastero'], df['dloggastero'])
    df['loggLC'] = paramError(df['loggLC'], df['dloggLC'])
    junk = ['tefferr', 'loggerr', 'feherr', 'vterr', 'dloggastero', 'dloggLC', 'convergence']
    df.drop(junk, inplace=True, axis=1)

    df.to_latex(buf=output, index=False, escape=False)
