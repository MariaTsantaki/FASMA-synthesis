#!/usr/bin/env python
# -*- coding: utf8 -*-
from __future__ import print_function
from synthDriver import synthMethod
from gooey import Gooey, GooeyParser

def synth(args):
    """Driver for the synthesis method"""
    fout = ''
    linelist = args.linelist.rpartition('/')[-1]
    inter_file = args.inter_file.rpartition('/')[-1]
    if not args.temperature:
        print('Temperature was set to 5777 K')
        args.temperature = 5777
    if not args.surfacegravity:
        print('Surface gravity was set to 4.44 dex')
        args.surfacegravity = 4.44
    if not args.FeH:
        args.FeH = 0.0
    if not args.microturbulence:
        print('Microturbulence was set to 1.00 km/s')
        args.microturbulence = 1.00
    if not args.macroturbulence:
        print('Macroturbulence was set to 3.21 km/s')
        args.macroturbulence = 3.21
    if not args.vsini:
        print('Rotation was set to 1.90 km/s')
        args.microturbulence = 1.90
    fout += '%s %s %s %s %s %s %s ' % (linelist, args.temperature, args.surfacegravity, args.FeH, args.microturbulence, args.macroturbulence, args.vsini)

    fout += 'model:%s,step_wave:%s,step_flux:%s,inter_file:%s,limb:%s,damping:%s' % (args.model, args.step_wave, args.step_flux, inter_file, args.limb, args.damping)
    if args.observations:
        fout += ',observations:%s' % args.observations
    if args.resolution:
        fout += ',resolution:%s' % args.resolution
    if args.snr:
        fout += ',snr:%s' % args.snr
    if args.refine:
        fout += ',refine'
    if args.plot:
        fout += ',plot'
    if args.plot_res:
        fout += ',plot_res'
    if args.save:
        fout += ',save'
    if args.minimize:
        fout += ',minimize'
    if args.Fixteff:
        fout += ',teff'
    if args.Fixgravity:
        fout += ',logg'
    if args.FixFeH:
        fout += ',feh'
    if args.Fixmicroturbulence:
        fout += ',vt'
    if args.Fixmacroturbulence:
        fout += ',vmac'
    if args.Fixvsini:
        fout += ',vsini'
    with open('StarMe_synth.cfg', 'w') as f:
        f.writelines(fout)
    driver = synthMethod(cfgfile='StarMe_synth.cfg', overwrite=None)
    driver.synthdriver()


@Gooey(program_name='FASMA - Spectral Synthesis with MOOG',
       default_size=(900, 1000),
       image_dir='./img')

def main():
    '''Take care of all the argparse stuff.
    :returns: the args
    '''
    parser = GooeyParser(description='Set parameters')

    subparsers = parser.add_subparsers()

    parent_parser = GooeyParser(add_help=False)
    parent_parser.add_argument('--temperature',     help='(in K)',    default=5777,   type=int,   metavar='Effective temperature')
    parent_parser.add_argument('--surfacegravity',  help='(in dex)',  default=4.44,   type=float, metavar='Surface gravity')
    parent_parser.add_argument('--FeH',             help='(in dex)',  default='0.00', type=float, metavar='Metallicity')
    parent_parser.add_argument('--microturbulence', help='(in km/s)', default=1.00,   type=float, metavar='Microturbulence')
    parent_parser.add_argument('--macroturbulence', help='(in km/s)', default=3.21,   type=float, metavar='Macroturbulence')
    parent_parser.add_argument('--vsini',           help='(in km/s)', default=1.90,   type=float, metavar='Rotational Velocity')
    parent_parser.add_argument('--model',           help='Grid of models', default='kurucz95', choices=['kurucz95', 'apogee_kurucz', 'marcs'], metavar='Model atmosphere')

    synth_parser = subparsers.add_parser('synth', parents=[parent_parser], help='Spectal synthesis')
    synth_parser.add_argument('linelist',             help='Input line list file', widget='FileChooser', metavar='Line list')
    synth_parser.add_argument('--inter_file',         help='Input wavelength intervals', default='intervals.lst', widget='FileChooser', metavar='Intervals file')
    synth_parser.add_argument('--observations',       help='Input spectrum in .fits, .dat, .spec', widget='FileChooser', metavar='Observed spectrum')
    synth_parser.add_argument('--damping',            help='Parameter for MOOG', default='1', choices=['0', '1', '2'], type=str, metavar='Damping')
    synth_parser.add_argument('--step_wave',          help='Step in wavelength for synthesis', default=0.01, type=float, metavar='Step Wave')
    synth_parser.add_argument('--step_flux',          help='Step in flux for synthesis', default=3.0, type=float, metavar='Step Flux')
    synth_parser.add_argument('--minimize',           help='Start minimization', action='store_true', metavar='Minimization')
    synth_parser.add_argument('--refine',             help='Refine parameters',   action='store_true', metavar='Refine parameters')
    synth_parser.add_argument('--Iterations',         help='Maximum number of iterations', default=100, type=int)
    synth_parser.add_argument('--Fixteff',            help='Fix Teff',   action='store_true', metavar='Fix temperature')
    synth_parser.add_argument('--Fixgravity',         help='Fix logg',   action='store_true', metavar='Fix gravity')
    synth_parser.add_argument('--FixFeH',             help='Fix [Fe/H]', action='store_true', metavar='Fix metallicity')
    synth_parser.add_argument('--Fixmicroturbulence', help='Fix vt',     action='store_true', metavar='Fix microturbulence')
    synth_parser.add_argument('--Fixmacroturbulence', help='Fix vmac',   action='store_true', metavar='Fix macroturbulence')
    synth_parser.add_argument('--Fixvsini',           help='Fix vsini',  action='store_true', metavar='Fix rotational velocity')
    synth_parser.add_argument('--resolution',         help='Instrumental resolution', default=None, type=int, metavar='Resolution')
    synth_parser.add_argument('--snr',                help='Signal-to-noise ratio of spectrum', default=None, type=float, metavar='SNR')
    synth_parser.add_argument('--plot',               help='Plot output', action='store_true', metavar='Plot')
    synth_parser.add_argument('--plot_res',           help='Add residuals in the plot', action='store_true', metavar='Plot residuals')
    synth_parser.add_argument('--limb',               help='Limb darkening coefficient', default=0.6, type=float, metavar='Limb Darkening')
    synth_parser.add_argument('--save',               help='Save synthetic spectrum', action='store_true', metavar='Save output')

    synth_parser.set_defaults(driver=synth)
    args = parser.parse_args()
    return args.driver(args)


if __name__ == '__main__':
    main()
