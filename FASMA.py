#!/usr/bin/env python
# -*- coding: utf8 -*-
from __future__ import print_function
from synthDriver import synthMethod
from gooey import Gooey, GooeyParser

def synth(args):
    '''Driver for the synthesis method for stellar parameters.
    '''

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

def abund(args):
    '''Driver for the synthesis method for element abundances.
    '''

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
    if args.element:
        fout += ',element:%s' % args.element
    if args.observations:
        fout += ',observations:%s' % args.observations
    if args.resolution:
        fout += ',resolution:%s' % args.resolution
    if args.snr:
        fout += ',snr:%s' % args.snr
    if args.plot:
        fout += ',plot'
    if args.plot_res:
        fout += ',plot_res'
    if args.save:
        fout += ',save'
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

    synth_parser = subparsers.add_parser('Parameters', parents=[parent_parser], help='Spectal Synthesis')
    synth_parser.add_argument('linelist',             help='List with atomic data', widget='FileChooser', metavar='Line List')
    synth_parser.add_argument('--inter_file',         help='List of wavelength intervals', default='intervals.lst', widget='FileChooser', metavar='Intervals File')
    synth_parser.add_argument('--observations',       help='Input formats: fits, dat, spec', widget='FileChooser', metavar='Observed Spectrum')
    synth_parser.add_argument('--damping',            help='Parameter for MOOG', default='1', choices=['0', '1', '2'], type=str, metavar='Damping')
    synth_parser.add_argument('--step_wave',          help='for synthetic spectrum', default=0.01, type=float, metavar='Step Wavelength')
    synth_parser.add_argument('--step_flux',          help='for synthetic spectrum', default=3.0, type=float, metavar='Step Flux')
    synth_parser.add_argument('--minimize',           help='Start minimization process', action='store_true', metavar='Minimization')
    synth_parser.add_argument('--refine',             help='Refine parameters',   action='store_true', metavar='Refine Parameters')
    synth_parser.add_argument('--Iterations',         help='Maximum number', default=100, type=int)
    synth_parser.add_argument('--Fixteff',            help='Fix Teff',   action='store_true', metavar='Fix Temperature')
    synth_parser.add_argument('--Fixgravity',         help='Fix logg',   action='store_true', metavar='Fix Gravity')
    synth_parser.add_argument('--FixFeH',             help='Fix [Fe/H]', action='store_true', metavar='Fix Metallicity')
    synth_parser.add_argument('--Fixmicroturbulence', help='Fix vt',     action='store_true', metavar='Fix Microturbulence')
    synth_parser.add_argument('--Fixmacroturbulence', help='Fix vmac',   action='store_true', metavar='Fix Macroturbulence')
    synth_parser.add_argument('--Fixvsini',           help='Fix vsini',  action='store_true', metavar='Fix Rotational Velocity')
    synth_parser.add_argument('--resolution',         help='of the spectrograph', default=None, type=int, metavar='Resolution')
    synth_parser.add_argument('--snr',                help='ratio of spectrum', default=None, type=float, metavar='Signal-to-noise')
    synth_parser.add_argument('--plot',               help='Plot output', action='store_true', metavar='Plot')
    synth_parser.add_argument('--plot_res',           help='Add residuals to the plot', action='store_true', metavar='Plot Residuals')
    synth_parser.add_argument('--limb',               help='coefficient for rotation', default=0.6, type=float, metavar='Limb Darkening')
    synth_parser.add_argument('--save',               help='Save synthetic spectrum', action='store_true', metavar='Save Output')
    synth_parser.set_defaults(driver=synth)

    abund_parser = subparsers.add_parser('Abundances', parents=[parent_parser], help='Spectal synthesis for chemical abundance')
    abund_parser.add_argument('linelist',             help='List with atomic data', widget='FileChooser', metavar='Line List')
    abund_parser.add_argument('--inter_file',         help='List of wavelength intervals', default='intervals_elements.lst', widget='FileChooser', metavar='Intervals File')
    abund_parser.add_argument('--element',            help='for abundance determination', choices=['Na', 'Mg', 'Al', 'Si', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Ni'], default='Na', metavar='Element')
    abund_parser.add_argument('--observations',       help='Input formats: fits, dat, spec', widget='FileChooser', metavar='Observed Spectrum')
    abund_parser.add_argument('--damping',            help='Parameter for MOOG', default='1', choices=['0', '1', '2'], type=str, metavar='Damping')
    abund_parser.add_argument('--step_wave',          help='for synthetic spectrum', default=0.01, type=float, metavar='Step Wavelength')
    abund_parser.add_argument('--step_flux',          help='for synthetic spectrum', default=3.0, type=float, metavar='Step Flux')
    abund_parser.add_argument('--Iterations',         help='Maximum number', default=100, type=int)
    abund_parser.add_argument('--resolution',         help='of the spectrograph', default=None, type=int, metavar='Resolution')
    abund_parser.add_argument('--snr',                help='ratio of spectrum', default=None, type=float, metavar='Signal-to-noise')
    abund_parser.add_argument('--plot',               help='Plot output', action='store_true', metavar='Plot')
    abund_parser.add_argument('--plot_res',           help='Add residuals to the plot', action='store_true', metavar='Plot Residuals')
    abund_parser.add_argument('--limb',               help='coefficient for rotation', default=0.6, type=float, metavar='Limb Darkening')
    abund_parser.add_argument('--save',               help='Save synthetic spectrum', action='store_true', metavar='Save Output')
    abund_parser.set_defaults(driver=abund)

    args = parser.parse_args()
    return args.driver(args)


if __name__ == '__main__':
    main()
