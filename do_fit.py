#!/usr/bin/env python

'''
Fit BEC C2 and R2 functions at 13 TeV
'''


from argparse import ArgumentParser
from os import system
import json
import sys
from datetime import datetime
from utils import remove_output_file, erase_line_in_table_file
from utils import add_histogram
from utils import add_latex_table_line, add_csv_table_line, add_csv_table_line_averages
import yaml

def main():
    """ Main fitting script. """
    parser = ArgumentParser(description='Yay, fit multiplicity or pair '
                                        'transverse momentum! (use reference sample option)')

    parser.add_argument('-u', '--uls', action='store_true',
                        default=False,
                        help='use ULS reference sample')
    parser.add_argument('-o', '--ohp', action='store_true',
                        default=False,
                        help='use OHP reference sample')
    parser.add_argument('-m', '--mix', action='store_true',
                        default=False,
                        help='use MIX reference sample')
    parser.add_argument('-n', '--mult', action='store_true',
                        default=False,
                        help='fit all multiplicities')
    parser.add_argument('-k', '--kt', action='store_true',
                        default=False,
                        help='fit all kts')
    parser.add_argument('--tt', action='store_true',
                        default=False,
                        help='fit all tts')
    parser.add_argument('--mult-index',
                        type=int, default=0,
                        help='Multiplicity index')
    parser.add_argument('--kt-index',
                        type=int, default=0,
                        help='kT index')
    parser.add_argument('--tt-index',
                        type=int, default=0,
                        help='tt index')
    parser.add_argument('-t', '--test', action='store_true',
                        default=False,
                        help='test fit')
    parser.add_argument('--c2-index',
                        type=int, default=2,
                        help='C2 index')
    parser.add_argument('--dont-use-eps', action='store_true',
                        default=False,
                        help='Don\'t use epsilon in fit')
    parser.add_argument('--dont-show-c0', action='store_true',
                        default=False,
                        help='Don\'t show C0 in plots')
    parser.add_argument('--fix-c0', action='store_true',
                        default=False,
                        help='Fix C0 in fits')
    parser.add_argument('--unf', action='store_true',
                        default=False,
                        help='Q spectra are unfolded')
    parser.add_argument('--style-8tev', action='store_true',
                        default=False,
                        help='8 TeV histogram style')
    parser.add_argument('--style-13tev', action='store_true',
                        default=False,
                        help='13 TeV histogram style')
    parser.add_argument('--style-dubna', action='store_true',
                        default=False,
                        help='Dubna histogram style')
    parser.add_argument('--enlarge-uncert', action='store_true',
                        default=False,
                        help='Use uncertainty enlargement')
    parser.add_argument('--unc-from',
                        type=float, default=1.,
                        help='uncertainty enlarged from Q')
    parser.add_argument('--unc-to',
                        type=float, default=-1.,
                        help='uncertainty enlarged to Q')
    parser.add_argument('--unc-factor',
                        type=int, default=-1.,
                        help='uncertainty in Q enlarged by')
    parser.add_argument('--unc2-from',
                        type=float, default=1.,
                        help='2nd uncertainty enlarged from Q')
    parser.add_argument('--unc2-to',
                        type=float, default=-1.,
                        help='2nd uncertainty enlarged to Q')
    parser.add_argument('--unc2-factor',
                        type=int, default=-1.,
                        help='2nd uncertainty in Q enlarged by')
    parser.add_argument('--q-min',
                        type=float, default=-1,
                        help='Q min')
    parser.add_argument('--q-max',
                        type=float, default=-1,
                        help='Q max')
    parser.add_argument('--rej-fin', action='store_true',
                        default=False,
                        help='Use final rejection')
    parser.add_argument('--rej-from',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej-to',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej2-from',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej2-to',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej3-from',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej3-to',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej-from-out',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej-to-out',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej-from-side',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej-to-side',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej-from-long',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej-to-long',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej2-from-out',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej2-to-out',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej2-from-side',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej2-to-side',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej2-from-long',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej2-to-long',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej3-from-out',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej3-to-out',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej3-from-side',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej3-to-side',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej3-from-long',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej3-to-long',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej4-from-out',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej4-to-out',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej4-from-side',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej4-to-side',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej4-from-long',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej4-to-long',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej5-from-out',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej5-to-out',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej5-from-side',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej5-to-side',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--rej5-from-long',
                        type=float, default=1.,
                        help='Excluded from Q')
    parser.add_argument('--rej5-to-long',
                        type=float, default=-1.,
                        help='Excluded from Q')
    parser.add_argument('--alpha-out',
                        type=float, default=-1.,
                        help='alpha out')
    parser.add_argument('--alpha-side',
                        type=float, default=-1.,
                        help='alpha side')
    parser.add_argument('--alpha-long',
                        type=float, default=-1.,
                        help='alpha long')
    parser.add_argument('--non-closure', action='store_true',
                        default=False,
                        help='Do non-closure correction')
    parser.add_argument('-i', '--input-yaml',
                    type=str, default='input/input.yaml',
                    help='Input yaml file location')
    # parser.add_argument('-i', '--input-json',
                        # type=str, default='input/file_list.json',
                        # help='Input JSON file location')

    args = parser.parse_args()

    rej_fin = args.rej_fin
    uls = args.uls
    ohp = args.ohp
    mix = args.mix
    input_yaml_path = args.input_yaml

    if args.test:
        uls = True
        args.mult = True

    if not uls and not ohp and not mix:
        parser.print_help()
        sys.exit(1)

    if args.mult and args.kt:
        parser.print_help()
        sys.exit(1)

    args.mult_index = args.mult_index - 1
    args.kt_index = args.kt_index - 1
    args.tt_index = args.tt_index - 1
    ttsplit = ''
    ttname = ''
    if args.mult_index >= 0 and args.kt_index >= 0:
        parser.print_help()
        sys.exit(1)
    if args.tt_index >= 0 and (args.kt_index >= 0 or args.mult_index >= 0):
        parser.print_help()
        sys.exit(1)

    if args.mult and args.mult_index >= 0:
        parser.print_help()
        sys.exit(1)

    if args.kt and args.kt_index >= 0:
        parser.print_help()
        sys.exit(1)

    if not args.style_8tev and not args.style_13tev and not args.style_dubna:
        args.style_13tev = True

    # system('make')

    try:
        with open(input_yaml_path, 'r') as stream:
            input_files = yaml.full_load(stream)
    except IOError:
        print('ERROR: Input file list not found!')
        exit(1)
    # try:
        # with open(file_list_path, 'r') as stream:
            # input_files = json.load(stream)
    # except IOError:
        # print('ERROR: Input file list not found!')
        # sys.exit(1)
    # except json.decoder.JSONDecodeError:
        # print('ERROR: Syntax error in input file list!')
        # sys.exit(1)
    file_data = ''
    file2_data = ''
    file_mc = ''
    file2_mc = ''
    # Settings
    if uls:
        try:
            file_data = input_files['uls']['file_data']
        except KeyError:
            pass
        try:
            file2_data = input_files['uls']['file2_data']
        except KeyError:
            pass
        try:
            file_mc = input_files['uls']['file_mc']
        except KeyError:
            pass
        try:
            file2_mc = input_files['uls']['file2_mc']
        except KeyError:
            pass
        plotname_prefix = 'bec_fit_uls'
        if args.style_8tev:
            hist_prefix = 'ppmm_q'
            hist2_prefix = 'pm_q'
            hist_suffix = ''
            hist2_suffix = ''
            hist_truth_suffix = ''
            hist2_truth_suffix = ''
            plotname_prefix += '_8tev'
        if args.style_13tev:
            hist_prefix = ''
            hist2_prefix = ''
            hist_suffix = 'Qosl_ppmm'
            hist2_suffix = 'Qosl_pm'
            hist_truth_suffix = 'truth_Q_ppmm'
            hist2_truth_suffix = 'truth_Q_pm'
            plotname_prefix += '_13tev'
        if args.style_dubna:
            hist_prefix = 'UCP_phN_ppmm'
            hist2_prefix = 'UCP_phN_pm'
            hist_suffix = '_W'
            hist2_suffix = '_W'
            hist_truth_suffix = '_truth'
            hist2_truth_suffix = '_truth'
            plotname_prefix += '_dubna'
        comment_prefix = 'ULS, variable binning;'
    if ohp:
        hist_prefix = ''
        hist2_prefix = ''
        try:
            file_data = input_files['ohp']['file_data']
        except KeyError:
            pass
        try:
            file2_data = input_files['ohp']['file2_data']
        except KeyError:
            pass
        try:
            file_mc = input_files['ohp']['file_mc']
        except KeyError:
            pass
        try:
            file2_mc = input_files['ohp']['file2_mc']
        except KeyError:
            pass
        # hist_suffix = 'Q_ppmm'
        # hist2_suffix = 'Q_ppmm'
        hist_suffix = 'Qosl_ppmm'
        hist2_suffix = 'Qosl_ppmm'
        hist_truth_suffix = 'truth_Q_ppmm'
        hist2_truth_suffix = 'truth_Q_ppmm'
        plotname_prefix = 'bec_fit_ohp'
        if args.style_13tev:
            plotname_prefix += '_13tev'
        comment_prefix = 'OHP, variable binning;'
    if mix:
        hist_prefix = ''
        hist2_prefix = ''
        try:
            file_data = input_files['mix']['file_data']
        except KeyError:
            pass
        try:
            file2_data = input_files['mix']['file2_data']
        except KeyError:
            pass
        try:
            file_mc = input_files['mix']['file_mc']
        except KeyError:
            pass
        try:
            file2_mc = input_files['mix']['file2_mc']
        except KeyError:
            pass
        hist_suffix = 'Qosl_ppmm'
        hist2_suffix = 'Qosl_ppmm'
        hist_truth_suffix = 'truth_Q_ppmm'
        hist2_truth_suffix = 'truth_Q_ppmm'
        plotname_prefix = 'bec_fit_mix'
        if args.style_13tev:
            plotname_prefix += '_13tev'
        comment_prefix = 'MIX, variable binning;'
    title = ''
    # Extract mixings from file name
    mixings = ''
    if mix:
        try:
            if 'MIX_' in file2_data:
                mixings = file2_data.split('MIX_')[1].split('_coulomb')[0]
        except NameError:
            pass
        try:
            if 'MIX_' in file2_mc:
                mixings = file2_mc.split('MIX_')[1].split('_coulomb')[0]
        except NameError:
            pass
        # Check if file2_data or file2_mc is not empty
    if not mix:
        try:
            if file2_data:
                if 'mix_thrust' in file2_data:
                    mixings = 'mix_thrust'
        except NameError:
            pass
        try:
            if file_data:
                if 'mix_thrust' in file_data:
                    mixings = 'mix_thrust'
        except NameError:
            pass
        try:
            if file2_mc:
                if 'mix_thrust' in file2_mc:
                    mixings = 'mix_thrust'
        except NameError:
            pass 
        try:
            if file_mc:
                if 'mix_thrust' in file_mc:
                    mixings = 'mix_thrust'
        except NameError:
            pass
    mix_thrust = False
    if mixings.find('mix_thrust') >= 0:
        mix_thrust = True
    plotname_prefix += mixings
    print('Mixings:', mixings)
    rej_from = args.rej_from
    rej_to = args.rej_to
    rej_from_Osl = [args.rej_from_out, args.rej_from_side, args.rej_from_long]
    rej2_from_Osl = [args.rej2_from_out, args.rej2_from_side, args.rej2_from_long]
    rej3_from_Osl = [args.rej3_from_out, args.rej3_from_side, args.rej3_from_long]
    rej4_from_Osl = [args.rej4_from_out, args.rej4_from_side, args.rej4_from_long]
    rej5_from_Osl = [args.rej5_from_out, args.rej5_from_side, args.rej5_from_long]
    rej_to_Osl = [args.rej_to_out, args.rej_to_side, args.rej_to_long]
    rej2_to_Osl = [args.rej2_to_out, args.rej2_to_side, args.rej2_to_long]
    rej3_to_Osl = [args.rej3_to_out, args.rej3_to_side, args.rej3_to_long]
    rej4_to_Osl = [args.rej4_to_out, args.rej4_to_side, args.rej4_to_long]
    rej5_to_Osl = [args.rej5_to_out, args.rej5_to_side, args.rej5_to_long]
    rej2_from = args.rej2_from
    rej2_to = args.rej2_to
    rej3_from = args.rej3_from
    rej3_to = args.rej3_to
    alpha = [args.alpha_out, args.alpha_side, args.alpha_long]
    unc_from = args.unc_from
    unc_to = args.unc_to
    unc_factor = args.unc_factor
    unc2_from = args.unc2_from
    unc2_to = args.unc2_to
    unc2_factor = args.unc2_factor
    # q_min = 0.02 #skip first bin, -1 for min
    # q_max = 1 # -1 for max (now 2.02)
    q_min = args.q_min
    q_max = args.q_max
    proj_min = 1
    proj_max = 30 #14 #10 for 200 MeV (14 for 280MeV, 22 for 600MeV)
    plotname_suffix = '_proj_%i_%i' % (proj_min, proj_max)
    pt = ""
    if file_data:
        pt = file_data.split('trk_pt_')[1].split('_')[0]
    else:
        pt = file_mc.split('trk_pt_')[1].split('_')[0]
    if file2_data and file2_mc:
        c2_type = "r2"
    elif file2_data:
        c2_type = "c2_data"
    elif file2_mc:
        c2_type = "c2_mc"
    else:
        print("ERROR: No file2 data or mc found!")
    plotname_prefix += '_pt_%s' % pt
    c2_index = args.c2_index
    plotname_prefix = plotname_prefix + '_%s_%i' % (c2_type, c2_index)
    if q_max > 0:
        plotname_suffix += '_qMax_%i' % (q_max * 1000.)
    rej_names = ['out', 'side', 'long']
    bin_start = 21
    bin_stop = 151


    if rej_to > rej_from:
        plotname_suffix += '_rej_%i_%i' % (rej_from * 1000., rej_to * 1000.)
    if rej2_to > rej2_from:
        plotname_suffix += '_rej2_%i_%i' % (rej2_from * 1000., rej2_to * 1000.)
    if rej3_to > rej3_from:
        plotname_suffix += '_rej3_%i_%i' % (rej3_from * 1000., rej3_to * 1000.)
    rej_Osl_lists = [
        (rej_from_Osl, rej_to_Osl, 'rej'),
        (rej2_from_Osl, rej2_to_Osl, 'rej2'),
        (rej3_from_Osl, rej3_to_Osl, 'rej3'),
        (rej4_from_Osl, rej4_to_Osl, 'rej4'),
        (rej5_from_Osl, rej5_to_Osl, 'rej5'),

    ]
    # for rej_from_temp, rej_to_temp, prefix in rej_Osl_lists:
    #     if all(a > b for a, b in zip(rej_to_temp, rej_from_temp)):
    #         plotname_suffix += ''.join(['_%s_%s_%i_%i' % (prefix, rej_names[i],
    #                             rej_from_temp[i] * 1000., rej_to_temp[i] * 1000.) 
    #                             for i in range(len(rej_names))])
    # Initialize a dictionary to store concatenated rejection regions
    rej_regions = {"out": [], "side": [], "long": []}

    # Concatenate rejection regions by type
    for (rej_from_osl, rej_to_osl, _) in rej_Osl_lists:
        if all(a > b for a, b in zip(rej_to_osl, rej_from_osl)):
            for key, idx in zip(rej_regions.keys(), range(3)):
                rej_regions[key].extend([int(rej_from_osl[idx] * 1000), int(rej_to_osl[idx] * 1000)])

    # Update plotname_suffix with concatenated rejection regions
    plotname_suffix += "".join(
        f"_rej_{key}_" + "_".join(f"{vals[i]}-{vals[i+1]}" for i in range(0, len(vals), 2))
        for key, vals in rej_regions.items()
    )

    if (unc_to > unc_from) and unc_factor >= 0.:
        plotname_suffix += '_unc_%i_%i_by_%i' % (unc_from * 1000.,
                                                 unc_to * 1000.,
                                                 unc_factor)
        comment_prefix += '%0.2f-%0.2f GeV uncert. x%i;' % (unc_from,
                                                            unc_to,
                                                            unc_factor)
    if (unc2_to > unc2_from) and unc2_factor >= 0.:
        plotname_suffix += '_unc2_%i_%i_by_%i' % (unc2_from * 1000.,
                                                  unc2_to * 1000.,
                                                  unc2_factor)
        comment_prefix += '%0.2f-%0.2f GeV uncert. x%i;' % (unc2_from,
                                                            unc2_to,
                                                            unc2_factor)

    if rej_fin:
        plotname_suffix += '_REJ'
        comment_prefix += 'REJ;'
    mults = (2, 10), (11, 20), (21, 30), (31, 40), \
            (41, 50), (51, 60), (61, 70), (71, 80), \
            (81, 90), (91, 100), (101, 125), (126, 150), \
            (151, 200), (201, 250), (21, 250)
    # mults = (2, 9), (10, 19), (20, 29), (30, 39), \
    #         (40, 49), (50, 59), (60, 69), (70, 79), \
    #         (80, 89), (90, 99), (100, 124), \
    #         (125, 149), (150, 199), (200, 249)

    kts = (100, 200), (200, 300), (300, 400), (400, 500), (500, 600), \
          (600, 700), (700, 1000), (1000, 1500), (1500, 2000)

    tts = (0.63, 0.71), (0.71, 0.78), (0.78, 0.85), (0.85, 0.93), (0.93, 1.0)
    
    if args.test:
        mults = ((51, 60), )
        # mults = ((50, 60), )

    if args.mult_index > len(mults):
        print("ERROR: Multiplicity index out of range!")
        sys.exit(1)
    if args.kt_index > len(kts):
        print("ERROR: kT index out of range!")
        sys.exit(1)
    if args.tt_index > len(tts):
        print("ERROR: tt index out of range!")
        sys.exit(1)
    
    args_prefix = ''
    try:
        args_prefix += ' --file-data "' + file_data +'"'
    except NameError:
        pass
    try:
        args_prefix += ' --file-mc "' + file_mc +'"'
    except NameError:
        pass
    try:
        args_prefix += ' --file2-data "' + file2_data +'"'
    except NameError:
        pass
    try:
        args_prefix += ' --file2-mc "' + file2_mc +'"'
    except NameError:
        pass

    args_prefix += ' --rej-from %f' % rej_from
    args_prefix += ' --rej-to %f' % rej_to
    for i, (rej_from_temp, rej_to_temp, prefix) in enumerate(rej_Osl_lists):
        args_prefix += ' --%s-from-out %f' % (prefix, rej_from_temp[0])
        args_prefix += ' --%s-to-out %f' % (prefix, rej_to_temp[0])
        args_prefix += ' --%s-from-side %f' % (prefix, rej_from_temp[1])
        args_prefix += ' --%s-to-side %f' % (prefix, rej_to_temp[1])
        args_prefix += ' --%s-from-long %f' % (prefix, rej_from_temp[2])
        args_prefix += ' --%s-to-long %f' % (prefix, rej_to_temp[2])
    # args_prefix += ' --rej-from-out %f' % rej_from_Osl[0]
    # args_prefix += ' --rej-to-out %f' % rej_to_Osl[0]
    # args_prefix += ' --rej-from-side %f' % rej_from_Osl[1]
    # args_prefix += ' --rej-to-side %f' % rej_to_Osl[1]
    # args_prefix += ' --rej-from-long %f' % rej_from_Osl[2]
    # args_prefix += ' --rej-to-long %f' % rej_to_Osl[2]
    # args_prefix += ' --rej2-from-out %f' % rej2_from_Osl[0]
    # args_prefix += ' --rej2-to-out %f' % rej2_to_Osl[0]
    # args_prefix += ' --rej2-from-side %f' % rej2_from_Osl[1]
    # args_prefix += ' --rej2-to-side %f' % rej2_to_Osl[1]
    # args_prefix += ' --rej2-from-long %f' % rej2_from_Osl[2]
    # args_prefix += ' --rej2-to-long %f' % rej2_to_Osl[2]
    # args_prefix += ' --rej3-from-out %f' % rej3_from_Osl[0]
    # args_prefix += ' --rej3-to-out %f' % rej3_to_Osl[0]
    # args_prefix += ' --rej3-from-side %f' % rej3_from_Osl[1]
    # args_prefix += ' --rej3-to-side %f' % rej3_to_Osl[1]
    # args_prefix += ' --rej3-from-long %f' % rej3_from_Osl[2]
    # args_prefix += ' --rej3-to-long %f' % rej3_to_Osl[2]
    args_prefix += ' --rej2-from %f' % rej2_from
    args_prefix += ' --rej2-to %f' % rej2_to
    args_prefix += ' --rej3-from %f' % rej3_from
    args_prefix += ' --rej3-to %f' % rej3_to
    args_prefix += ' --alpha-out %f' % alpha[0]
    args_prefix += ' --alpha-side %f' % alpha[1]
    args_prefix += ' --alpha-long %f' % alpha[2]
    args_prefix += ' --q-min %f' % q_min
    args_prefix += ' --q-max %f' % q_max
    args_prefix += ' --proj-min %d' % proj_min
    args_prefix += ' --proj-max %d' % proj_max
    args_prefix += ' --title "%s"' % title
    args_prefix += ' --c2-index %d' % c2_index
    if rej_fin:
        args_prefix += ' --rej-fin'
    if not args.dont_use_eps:
        args_prefix += " --use-eps"
    if args.fix_c0:
        args_prefix += " --fix-c0"
    if not args.dont_show_c0:
        args_prefix += " --show-c0"
    if args.enlarge_uncert:
        args_prefix += " --enlarge-uncert"
        args_prefix += ' --unc-from %f' % unc_from
        args_prefix += ' --unc-to %f' % unc_to
        args_prefix += ' --unc-factor %d' % unc_factor
        args_prefix += ' --unc2-from %f' % unc2_from
        args_prefix += ' --unc2-to %f' % unc2_to
        args_prefix += ' --unc2-factor %d' % unc2_factor
    c2_index = c2_type + "_" + str(c2_index)
    #
    # All multiplicity, all kT, all tt
    #
    if not args.mult and not args.kt and not args.tt:
        if args.mult_index >= 0:
            mult1, mult2 = mults[args.mult_index]
        else:
            mult1 = mults[0][0]
            mult2 = mults[-1][-1]
        if args.kt_index >= 0:
            kt1, kt2 = kts[args.kt_index]
        else:
            kt1 = kts[0][0]
            kt2 = kts[-1][-1]
        if mix_thrust:
            if args.tt_index >= 0:
                tt1, tt2 = tts[args.tt_index]
                ttsplit = 'tt_%0.2f_%0.2f/' % (tt1, tt2)
                if tt2 == 1.0:
                    ttsplit = 'tt_%0.2f_%0.1f/' % (tt1, tt2)
                ttname = '_tt_%0.2f_%0.2f' % (tt1, tt2)
            else:
                ttsplit = 'no_tt_split/'
                tt1 = tts[0][0]
                tt2 = tts[-1][-1]
                ttname = '_tt_%0.2f_%0.2f' % (tt1, tt2)
        if ttname != '':
            ttname = ttname.replace('.', '_')
        if args.kt_index >= 0:
            tablename_mult = plotname_prefix + '_mult_kt_%i_%i' % (kt1, kt2) + \
                             plotname_suffix
        else:
            tablename_mult = plotname_prefix + '_mult' + plotname_suffix
        if args.mult_index >= 0:
            tablename_kt = plotname_prefix + \
                           '_kt_mult_%i_%i' % (mult1, mult2) + \
                           plotname_suffix
        else:
            tablename_kt = plotname_prefix + '_kt' + plotname_suffix
        if mix_thrust and args.tt_index >= 0:
            tablename_tt = plotname_prefix + '_mult' + ttname + plotname_suffix
        else:
            tablename_tt = plotname_prefix + '_tt' + plotname_suffix
        erase_line_in_table_file(tablename_mult, 'tex',
                                 '%i--%i ' % (mult1, mult2), c2_index)
        erase_line_in_table_file(tablename_mult, 'csv',
                                 '%i,%i,' % (mult1, mult2), c2_index)
        erase_line_in_table_file(tablename_kt, 'tex', '%i--%i ' % (kt1, kt2), c2_index)
        erase_line_in_table_file(tablename_kt, 'csv', '%i,%i,' % (kt1, kt2), c2_index)
        if mix_thrust:
            erase_line_in_table_file(tablename_tt, 'tex', '%.2f--%.2f ' % (tt1, tt2), c2_index)
            erase_line_in_table_file(tablename_tt, 'csv', '%.2f,%.2f,' % (tt1, tt2), c2_index)
        if args.style_8tev:
            hist_data = '_g_n_%i_%i' % (mult1, mult2)
            hist_mc = '_n_%i_%i' % (mult1, mult2)
        if args.style_13tev:
            hist_data = 'no_multiplicity_split/no_kt_split/%s' % ttsplit
            hist_mc = 'no_multiplicity_split/no_kt_split/%s' % ttsplit
            if args.kt_index >= 0:
                hist_data = 'no_multiplicity_split/kt_%i_%i/%s' % (kt1, kt2, ttsplit)
                hist_mc = 'no_multiplicity_split/kt_%i_%i/%s' % (kt1, kt2, ttsplit)
            if args.mult_index >= 0:
                hist_data = 'multiplicity_%i_%i/no_kt_split/%s' % (mult1, mult2, ttsplit)
                hist_mc = 'multiplicity_%i_%i/no_kt_split/%s' % (mult1, mult2, ttsplit)
            #for now we need hist2data to be without tt split
            hist2_data = hist_data
            hist2_mc = hist_mc
            # hist_data = 'no_multiplicity_split/no_kt_split/'
            # hist_mc = 'no_multiplicity_split/no_kt_split/'
        if args.style_dubna:
            hist_data = '_kT%i-%i_Q20_N%i-%i_G' % (kt1, kt2, mult1, mult2)
            hist_mc = '_kT%i-%i_Q20_N%i-%i' % (kt1, kt2, mult1, mult2)
        plotname = plotname_prefix + \
                   '_mult_%i_%i_kt_%i_%i%s' % (mult1, mult2, kt1, kt2, ttname) + \
                   plotname_suffix
        comment = comment_prefix + '%i #leq n_{ch} #leq %i' % (mult1, mult2)

        # fit_command = 'gdbserver :1234 ./fit' + args_prefix + \
        fit_command = './fit' + args_prefix + \
            ' --hist-data ' + hist_prefix + hist_data + hist_suffix + \
            ' --hist2-data ' + hist2_prefix + hist2_data + hist2_suffix + \
            ' --hist-mc ' + hist_prefix + hist_mc + hist_suffix + \
            ' --hist2-mc ' + hist2_prefix + hist2_mc + hist2_suffix
        if args.non_closure:
            fit_command += ' --hist-truth ' + \
                           hist_prefix + hist_mc + hist_truth_suffix
            fit_command += ' --hist2-truth ' + \
                           hist2_prefix + hist_mc + hist2_truth_suffix
        fit_command += ' --plot-name ' + plotname + \
            ' --comment "' + comment + '"'

        print(fit_command)
        return_val = system(fit_command)
        if return_val != 0:
            sys.exit(1)
        if args.mult_index < 0:
            #make_filepath(filename, extension):
            add_latex_table_line(tablename_mult, plotname,
                                 '%i--%i ' % (mult1, mult2),
                                 c2_index, args.kt_index >= 0)
            add_csv_table_line(tablename_mult, plotname,
                               '%i,%i,' % (mult1, mult2),
                               c2_index, args.kt_index >= 0)
        if args.kt_index < 0:
            add_latex_table_line(tablename_kt, plotname, '%i--%i ' % (kt1, kt2), c2_index, args.tt_index >= 0)
            add_csv_table_line(tablename_kt, plotname, '%i,%i,' % (kt1, kt2), c2_index, args.tt_index >= 0)

        if mix_thrust and args.tt_index < 0:
            add_latex_table_line(tablename_tt, plotname, '%.2f--%.2f ' % (tt1, tt2), c2_index)
            add_csv_table_line(tablename_tt, plotname, '%.2f,%.2f,' % (tt1, tt2), c2_index)

        add_csv_table_line_averages(tablename_mult, bin_start, bin_stop, c2_index)


    if args.unf and args.style_13tev:
        hist_suffix += '_unf'
        hist2_suffix += '_unf'


    #
    # Fit multiplicity bins
    #
    if args.mult and not args.kt and not args.tt:
        if args.kt_index >= 0:
            kt1, kt2 = kts[args.kt_index]
        # only do, if mixings contains "mix_thrust"
        if mix_thrust:
            if args.tt_index >= 0:
                tt1, tt2 = tts[args.tt_index]
                ttsplit = 'tt_%0.2f_%0.2f/' % (tt1, tt2)
                if tt2 == 1.0:
                    ttsplit = 'tt_%0.2f_%0.1f/' % (tt1, tt2)
                ttname = '_tt_%0.2f_%0.2f' % (tt1, tt2)
            else:
                ttsplit = 'no_tt_split/'

        colname = plotname_prefix + '_mult' + plotname_suffix
        if args.kt_index >= 0:
            colname = plotname_prefix + '_mult_kt_%i_%i%s' % (kt1, kt2, ttname) + \
                      plotname_suffix
        if args.tt_index >= 0:
            colname = plotname_prefix + '_mult%s' % ttname + plotname_suffix
        remove_output_file(colname, 'tex', c2_index)
        remove_output_file(colname, 'csv', c2_index)
        remove_output_file(colname, 'root', c2_index)

        for mult1, mult2 in mults:
            if args.style_8tev:
                hist_data = '_g_n_%i_%i' % (mult1, mult2)
                hist_mc = '_n_%i_%i' % (mult1, mult2)
            if args.style_13tev:
                hist_data = 'multiplicity_%i_%i/no_kt_split/%s' % (mult1, mult2, ttsplit)
                hist_mc = 'multiplicity_%i_%i/no_kt_split/%s' % (mult1, mult2, ttsplit)
                if args.kt_index >= 0:
                    hist_data = 'multiplicity_%i_%i/kt_%i_%i/%s' % (mult1, mult2,
                                                                    kt1, kt2, ttsplit)
                    hist_mc = 'multiplicity_%i_%i/kt_%i_%i/%s' % (mult1, mult2,
                                                                  kt1, kt2, ttsplit)
                hist2_data = hist_data
                hist2_mc = hist_mc
                # hist_data = 'multiplicity_%i_%i/no_kt_split/' % (mult1, mult2)
                # hist_mc = 'multiplicity_%i_%i/no_kt_split/' % (mult1, mult2)
            if args.style_dubna:
                hist_data = '_pT100_Q20_Nch%i-%i_G' % (mult1, mult2)
                hist_mc = '_pT100_Q20_Nch%i-%i' % (mult1, mult2)
               # hist_data = '_kT100-2000_Q20_N%i-%i_G' % (mult1, mult2)
               # hist_mc = '_kT100-2000_Q20_N%i-%i' % (mult1, mult2)
                if args.kt_index >= 0:
                    hist_data = '_kT%i-%i_Q20_N%i-%i_G' % (kt1, kt2,
                                                           mult1, mult2)
                    hist_mc = '_kT%i-%i_Q20_N%i-%i' % (kt1, kt2,
                                                       mult1, mult2)

            plotname = plotname_prefix + '_mult_%i_%i' % (mult1, mult2) + \
                       plotname_suffix
            if args.kt_index >= 0:
                plotname = plotname_prefix + \
                           '_mult_%i_%i_kt_%i_%i%s' % (mult1, mult2, kt1, kt2, ttname) + \
                           plotname_suffix
            if args.tt_index >= 0:
                plotname = plotname_prefix + \
                           '_mult_%i_%i%s' % (mult1, mult2, ttname) + \
                           plotname_suffix
            comment = comment_prefix + '%i #leq n_{ch} #leq %i' % (mult1, mult2)

            # fit_command = 'gdbserver :1234 ./fit' + args_prefix + \
            fit_command = './fit' + args_prefix + \
                ' --hist-data ' + hist_prefix + hist_data + hist_suffix + \
                ' --hist2-data ' + hist2_prefix + hist2_data + hist2_suffix + \
                ' --hist-mc ' + hist_prefix + hist_mc + hist_suffix + \
                ' --hist2-mc ' + hist2_prefix + hist2_mc + hist2_suffix
            if args.non_closure:
                fit_command += ' --hist-truth ' + \
                               hist_prefix + hist_mc + hist_truth_suffix
                fit_command += ' --hist2-truth ' + \
                               hist2_prefix + hist_mc + hist2_truth_suffix
            fit_command += ' --plot-name ' + plotname + \
                ' --comment "' + comment + '"'

            print(fit_command)
            return_val = system(fit_command)
            if return_val != 0:
                sys.exit(1)

            add_latex_table_line(colname, plotname, '%i--%i ' % (mult1, mult2), c2_index)
            print(colname, plotname, '%i--%i ' % (mult1, mult2), c2_index)
            add_csv_table_line(colname, plotname, '%i,%i,' % (mult1, mult2), c2_index)
            add_histogram(colname, plotname, c2_index)


    #
    # Fit kT bins
    #
    if args.kt and not args.mult and not args.tt:
        if args.mult_index >= 0:
            mult1, mult2 = mults[args.mult_index]
        if mix_thrust:
            if args.tt_index >= 0:
                tt1, tt2 = tts[args.tt_index]
                ttsplit = 'tt_%0.2f_%0.2f/' % (tt1, tt2)
                if tt2 == 1.0:
                    ttsplit = 'tt_%0.2f_%0.1f/' % (tt1, tt2)
                ttname = '_tt_%0.2f_%0.2f' % (tt1, tt2)
            else:
                ttsplit = 'no_tt_split/'

        colname = plotname_prefix + '_kt' + plotname_suffix
        if args.mult_index >= 0:
            colname = plotname_prefix + '_kt_mult_%i_%i%s' % (mult1, mult2, ttname) + \
                      plotname_suffix
        elif args.tt_index >= 0:
            colname = plotname_prefix + '_kt%s' % ttname + plotname_suffix

        remove_output_file(colname, 'tex', c2_index)
        remove_output_file(colname, 'csv', c2_index)
        remove_output_file(colname, 'root', c2_index)

        for kt1, kt2 in kts:
            if args.style_8tev:
                hist_data = '_g_kt_%i_%i' % (kt1, kt2)
                hist_mc = '_kt_%i_%i' % (kt1, kt2)
            if args.style_13tev:
                hist_data = 'no_multiplicity_split/kt_%i_%i/%s' % (kt1, kt2, ttsplit)
                hist_mc = 'no_multiplicity_split/kt_%i_%i/%s' % (kt1, kt2, ttsplit)
                if args.mult_index >= 0:
                    hist_data = 'multiplicity_%i_%i/kt_%i_%i/%s' % (mult1, mult2,
                                                                    kt1, kt2, ttsplit)
                    hist_mc = 'multiplicity_%i_%i/kt_%i_%i/%s' % (mult1, mult2,
                                                                  kt1, kt2, ttsplit)
                hist2_data = hist_data
                hist2_mc = hist_mc
                # hist_data = 'no_multiplicity_split/kt_%i_%i/' % (kt1, kt2)
                # hist_mc = 'no_multiplicity_split/kt_%i_%i/' % (kt1, kt2)
            if args.style_dubna:
                hist_data = '_kT%i-%i_Q20_N2-250_G' % (kt1, kt2)
                hist_mc = '_kT%i-%i_Q20_N2-250' % (kt1, kt2)
                if args.mult_index >= 0:
                    hist_data = '_kT%i-%i_Q20_N%i-%i_G' % (kt1, kt2,
                                                           mult1, mult2)
                    hist_mc = '_kT%i-%i_Q20_N%i-%i' % (kt1, kt2,
                                                       mult1, mult2)
            plotname = plotname_prefix + '_kt_%i_%i' % (kt1, kt2) + \
                       plotname_suffix
            if args.mult_index >= 0:
                plotname = plotname_prefix + \
                           '_kt_%i_%i_mult_%i_%i' % (kt1, kt2, mult1, mult2) + \
                           plotname_suffix
            comment = comment_prefix + '%i #leq k_{T} < %i [MeV]' % (kt1, kt2)

            # fit_command = 'gdbserver :1234 ./fit' + args_prefix + \
            fit_command = './fit' + args_prefix + \
                ' --hist-data ' + hist_prefix + hist_data + hist_suffix + \
                ' --hist2-data ' + hist2_prefix + hist2_data + hist2_suffix + \
                ' --hist-mc ' + hist_prefix + hist_mc + hist_suffix + \
                ' --hist2-mc ' + hist2_prefix + hist2_mc + hist2_suffix
            if args.non_closure:
                fit_command += ' --hist-truth ' + \
                               hist_prefix + hist_mc + hist_truth_suffix
                fit_command += ' --hist2-truth ' + \
                               hist2_prefix + hist_mc + hist2_truth_suffix
            fit_command += ' --plot-name ' + plotname + \
                ' --comment "' + comment + '"'
            print(c2_index)
            print(fit_command)
            return_val = system(fit_command)
            if return_val != 0:
                sys.exit(1)

            add_latex_table_line(colname, plotname, '%i--%i ' % (kt1, kt2), c2_index)
            add_csv_table_line(colname, plotname, '%i,%i, ' % (kt1, kt2), c2_index)
            add_histogram(colname, plotname, c2_index)

    #
    # Fit tt bins
    #
    if mix_thrust and args.tt and not args.mult and not args.kt:
        if args.mult_index >= 0:
            mult1, mult2 = mults[args.mult_index]
        if args.kt_index >= 0:
            kt1, kt2 = kts[args.kt_index]

        colname = plotname_prefix + '_tt' + plotname_suffix
        if args.mult_index >= 0 and args.kt_index >= 0:
            colname = plotname_prefix + '_tt_mult_%i_%i_kt_%i_%i' % (mult1, mult2, kt1, kt2) + \
                      plotname_suffix
        elif args.mult_index >= 0:
            colname = plotname_prefix + '_tt_mult_%i_%i' % (mult1, mult2) + plotname_suffix
        elif args.kt_index >= 0:
            colname = plotname_prefix + '_tt_kt_%i_%i' % (kt1, kt2) + plotname_suffix

        colname = colname.replace('.', '_')
        remove_output_file(colname, 'tex', c2_index)
        remove_output_file(colname, 'csv', c2_index)
        remove_output_file(colname, 'root', c2_index)

        for tt1, tt2 in tts:
            ttsplit = 'tt_%0.2f_%0.2f/' % (tt1, tt2)
            if tt2 == 1.0:
                ttsplit = 'tt_%0.2f_%0.1f/' % (tt1, tt2)
            ttname = '_tt_%0.2f_%0.2f' % (tt1, tt2)
            ttname = ttname.replace('.', '_')
            if args.style_8tev:
                hist_data = '_g_tt_%0.2f_%0.2f' % (tt1, tt2)
                hist_mc = '_tt_%0.2f_%0.2f' % (tt1, tt2)
            if args.style_13tev:
                hist_data = 'no_multiplicity_split/no_kt_split/%s' % ttsplit
                hist_mc = 'no_multiplicity_split/no_kt_split/%s' % ttsplit
                if args.mult_index >= 0 and args.kt_index >= 0:
                    hist_data = 'multiplicity_%i_%i/kt_%i_%i/%s' % (mult1, mult2,
                                                                    kt1, kt2, ttsplit)
                    hist_mc = 'multiplicity_%i_%i/kt_%i_%i/%s' % (mult1, mult2,
                                                                  kt1, kt2, ttsplit)
                elif args.mult_index >= 0:
                    hist_data = 'multiplicity_%i_%i/no_kt_split/%s' % (mult1, mult2, ttsplit)
                    hist_mc = 'multiplicity_%i_%i/no_kt_split/%s' % (mult1, mult2, ttsplit)
                elif args.kt_index >= 0:
                    hist_data = 'no_multiplicity_split/kt_%i_%i/%s' % (kt1, kt2, ttsplit)
                    hist_mc = 'no_multiplicity_split/kt_%i_%i/%s' % (kt1, kt2, ttsplit)
                hist2_data = hist_data
                hist2_mc = hist_mc
                # hist_data = 'no_multiplicity_split/no_kt_split/'
                # hist_mc = 'no_multiplicity_split/no_kt_split/'

            # if args.style_dubna:
            plotname = plotname_prefix + ttname + plotname_suffix
            if args.mult_index >= 0 and args.kt_index >= 0:
                plotname = plotname_prefix + ttname + \
                           '_mult_%i_%i_kt_%i_%i' % (mult1, mult2, kt1, kt2) + \
                           plotname_suffix
            elif args.mult_index >= 0:
                plotname = plotname_prefix + ttname + '_mult_%i_%i' % (mult1, mult2) + \
                           plotname_suffix
            elif args.kt_index >= 0:
                plotname = plotname_prefix + ttname + '_kt_%i_%i' % (kt1, kt2) + \
                           plotname_suffix
            comment = comment_prefix + '%.2f #leq tt < %.2f' % (tt1, tt2)

            fit_command = './fit' + args_prefix + \
                ' --hist-data ' + hist_prefix + hist_data + hist_suffix + \
                ' --hist2-data ' + hist2_prefix + hist2_data + hist2_suffix + \
                ' --hist-mc ' + hist_prefix + hist_mc + hist_suffix + \
                ' --hist2-mc ' + hist2_prefix + hist2_mc + hist2_suffix
            if args.non_closure:
                fit_command += ' --hist-truth ' + hist_prefix + hist_mc + hist_truth_suffix
                fit_command += ' --hist2-truth ' + hist2_prefix + hist_mc + hist2_truth_suffix
            fit_command += ' --plot-name ' + plotname + ' --comment "' + comment + '"'

            print(fit_command)
            return_val = system(fit_command)
            if return_val != 0:
                sys.exit(1)

            add_latex_table_line(colname, plotname, '%.2f--%.2f ' % (tt1, tt2), c2_index)
            add_csv_table_line(colname, plotname, '%.2f,%.2f, ' % (tt1, tt2), c2_index)
            add_histogram(colname, plotname, c2_index)


if __name__ == '__main__':
    main()
