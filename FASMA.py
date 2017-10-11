#!/usr/bin/env python
# -*- coding: utf8 -*-

from __future__ import print_function
from synthDriver import synthdriver
from gooey import Gooey, GooeyParser


def synth(args):
    """Driver for the synthesis method"""
    fout = ''
    linelist = args.linelist.rpartition('/')[-1]
    inter_file = args.inter_file.rpartition('/')[-1]
    if args.spectralType:
        if not args.temperature and not args.surfacegravity:
            fout += '%s spt:%s,' % (linelist, args.spectralType)
        else:
            print('Temperature and/or surface gravity set. Will not use spectral type.')
    else:
        if not args.temperature:
            print('Warning: Temperature was set to 5777 K')
            args.temperature = 5777
        if not args.surfacegravity:
            print('Warning: Surface gravity was set to 4.44 dex')
            args.surfacegravity = 4.44
        if not args.FeH:
            args.FeH = 0.0
        if not args.microturbulence:
            print('Warning: Microturbulence was set to 1.00 km/s')
            args.microturbulence = 1.00
        if not args.macroturbulence:
            print('Warning: Macroturbulence was set to 3.21 km/s')
            args.macroturbulence = 3.21
        if not args.rotation:
            print('Warning: Rotation was set to 1.90 km/s')
            args.microturbulence = 1.90
        fout += '%s %s %s %s %s %s %s ' % (linelist, args.temperature, args.surfacegravity, args.FeH, args.microturbulence, args.macroturbulence, args.rotation)

    fout += 'model:%s,MOOGv:%s,step_wave:%s,step_flux:%s,inter_file:%s,limb:%s,damping:%s' % (args.model, args.MOOGv, args.step_wave, args.step_flux, inter_file, args.limb, args.damping)
    if args.observations:
        fout += ',observations:%s' % args.observations
    if args.resolution:
        fout += ',resolution:%s' % args.resolution
    if args.snr:
        fout += ',snr:%s' % args.snr
    if args.refine:
        fout += ',refine'
    if args.errors:
        fout += ',errors'
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
    if args.Fixrotation:
        fout += ',vsini'
    if args.flag_vt:
        fout += ',flag_vt'
    if args.flag_vmac:
        fout += ',flag_vmac'
    with open('StarMe_synth.cfg', 'w') as f:
        f.writelines(fout)
    synthdriver(overwrite=args.overwrite)


@Gooey(program_name='FASMA - spectral analysis',
       default_size=(900, 1000),
       image_dir='./img')
def main():
    '''Take care of all the argparse stuff.
    :returns: the args
    '''
    parser = GooeyParser(description='Set parameters')

    subparsers = parser.add_subparsers()

    # Common to all
    parent_parser = GooeyParser(add_help=False)
    parent_parser.add_argument('--temperature',     help='(in K)',    default=5777,  type=int,   metavar='Temperature')
    parent_parser.add_argument('--surfacegravity',  help='(in dex)',  default=4.44,  type=float, metavar='logg')
    parent_parser.add_argument('--FeH',             help='(in dex)',  default='0.00',type=float, metavar='[Fe/H]')
    parent_parser.add_argument('--microturbulence', help='(in km/s)', default=1.0,   type=float, metavar='Microturbulence')
    parent_parser.add_argument('--MOOGv',           default='2014', choices=['2013', '2014', '2016'], type=str, metavar='MOOG version')
    parent_parser.add_argument('--model',           help='Grid of models', default='marcs', choices=['apogee_kurucz', 'kurucz08', 'marcs'], metavar='Model atmosphere')

    # For the synhtesis method
    synth_parser = subparsers.add_parser('synth', parents=[parent_parser], help='Synthesis method')
    synth_parser.add_argument('linelist',             help='Input linelist file', widget='FileChooser')
    synth_parser.add_argument('--macroturbulence',    help='(in km/s)',  default=3.21, type=float, metavar='Macroturbulence')
    synth_parser.add_argument('--rotation',           help='(in km/s)',  default=1.90, type=float, metavar='Rotational velocity')
    synth_parser.add_argument('--spectralType',       help='(optional)', default=False, metavar='Spectral type')
    synth_parser.add_argument('--minimize',           help='Start minimization', action='store_true', metavar='Minimization procedure')
    synth_parser.add_argument('--Fixteff',            help='Fix Teff',   action='store_true', metavar='Fix temperature')
    synth_parser.add_argument('--Fixgravity',         help='Fix logg',   action='store_true', metavar='Fix gravity')
    synth_parser.add_argument('--FixFeH',             help='Fix [Fe/H]', action='store_true', metavar='Fix metallicity')
    synth_parser.add_argument('--Fixmicroturbulence', help='Fix vt',     action='store_true', metavar='Fix microturbulence')
    synth_parser.add_argument('--Fixmacroturbulence', help='Fix vmac',   action='store_true', metavar='Fix macroturbulence')
    synth_parser.add_argument('--Fixrotation',        help='Fix vsini',  action='store_true', metavar='Fix rotation')
    synth_parser.add_argument('--flag_vt',            help='For each iteration', action='store_true', metavar='Change vt in minimization')
    synth_parser.add_argument('--flag_vmac',          help='For each iteration', action='store_true', metavar='Change vmac in minimization')
    synth_parser.add_argument('--step_wave',          help='(in Angstroms)', default=0.01, type=float, metavar='Wavelength step for synthesis')
    synth_parser.add_argument('--step_flux',          help='(in Angstroms)', default=10.0, type=float, metavar='Flux step for synthesis')
    synth_parser.add_argument('--inter_file',         help='File with the intervals', metavar='Intervals', widget='FileChooser')
    synth_parser.add_argument('--limb',               help='Coefficient for vsini broadening', default=0.6, type=float, metavar='Limb darkening')
    synth_parser.add_argument('--damping',            help='Van der waals damping', default='1', metavar='Damping option',  choices=['0', '1', '2'])
    synth_parser.add_argument('--observations',       help='File with the observations', widget='FileChooser', metavar='Observations')
    synth_parser.add_argument('--resolution',         help='Instrumental resolution', default=None, metavar='Resolution')
    synth_parser.add_argument('--snr',                help='Signal-to-noise ratio',   default=None, metavar='SNR')
    synth_parser.add_argument('--save',               help='Save spectrum', action='store_true', metavar='Save output')
    synth_parser.add_argument('--refine',             help='Refine results', action='store_true', metavar='Refine results')
    synth_parser.add_argument('--errors',             help='Calculate errors', action='store_true', metavar='Errors')
    synth_parser.add_argument('--plot',               help='Plot spectrum', action='store_true', metavar='Plot spectrum')
    synth_parser.add_argument('--plot_res',           help='Plot residuals', action='store_true', metavar='Plot residuals')
    synth_parser.add_argument('--overwrite',          help='Overwrite results.csv', action='store_true', default=False)
    synth_parser.set_defaults(driver=synth)

    args = parser.parse_args()
    return args.driver(args)


if __name__ == '__main__':
    main()
