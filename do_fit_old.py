#!/usr/bin/env python

from argparse import ArgumentParser
from os import system
from os import path
from sys import exit
from datetime import datetime
import yaml


def make_filepath(filename, extension):
    filepath = 'output/' + extension
    filepath += '/%s_%s_%s' % (datetime.now().year,
                               datetime.now().strftime('%h'),
                               datetime.now().day)
    system('mkdir -p ' + filepath)
    filepath = filepath + '/' + filename + '.' + extension

    return filepath


def add_table_line(filename, line_prefix=''):
    filepath = make_filepath(filename, 'tex')

    if line_prefix:
        in_file = open('fit.out', 'r')
        line = in_file.read()
        in_file.close()
        system('rm -f fit.out')

        out_file = open(filepath, 'a+')
        out_file.write(line_prefix + line)
        out_file.close()
    else:
        system('rm -f ' + filepath)


def add_csv_line(filename, line_prefix=''):
    filepath = make_filepath(filename, 'csv')

    if line_prefix:
        in_file = open('fit.csv', 'r')
        lines = in_file.readlines()
        in_file.close()
        system('rm -f fit.csv')

        if not path.isfile(filepath):
            out_file = open(filepath, 'w')
            out_file.write('bin_start,bin_stop,' + lines[0])
            out_file.close()

        out_file = open(filepath, 'a+')
        out_file.write(line_prefix + lines[1])
        out_file.close()
    else:
        system('rm -f ' + filepath)


parser = ArgumentParser(description='Yay, fit multiplicity or pair '
                                    'transverse momentum!')

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
parser.add_argument('-t', '--test', action='store_true',
                    default=False,
                    help='test fit')
parser.add_argument('--c2-index',
                    type=int, default=1,
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
parser.add_argument('--enlarge-uncert', action='store_true',
                    default=False,
                    help='Use uncertainity enlargement')
parser.add_argument('--rej-from',
                    type=int, default=1,
                    help='Excluded from Q')
parser.add_argument('--rej-to',
                    type=int, default=1,
                    help='Excluded from Q')
parser.add_argument('-i', '--input-yaml',
                    type=str, default='input/input.yaml',
                    help='Input yaml file location')
parser.add_argument('--hist-type',
                    type=str, default='_qosl_g',
                    help='Histogram type')


args = parser.parse_args()

uls = args.uls
ohp = args.ohp
mix = args.mix
mult = args.mult
kt = args.kt
test = args.test
input_yaml_path = args.input_yaml
hist_type = args.hist_type

if test:
    uls = True
    mult = True

if not uls and not ohp and not mix:
    parser.print_help()
    exit(1)

if not mult and not kt:
    parser.print_help()
    exit(1)

system('make')

try:
    with open(input_yaml_path, 'r') as stream:
        input_files = yaml.load(stream)
except IOError:
    print('ERROR: Input file list not found!')
    exit(1)

# Settings
if uls:
    file_data = input_files['uls']['file_data']
    file2_data = input_files['uls']['file2_data']
    file_mc = input_files['uls']['file_mc']
    file2_mc = input_files['uls']['file2_mc']
    hist_prefix = 'ppmm'
    hist2_prefix = 'pm'
    plotname_prefix = 'bec_3d_fit_uls'
    comment_prefix = 'ULS, variable binning;'
if ohp:
    file_data = input_files['ohp']['file_data']
    file2_data = input_files['ohp']['file2_data']
    file_mc = input_files['ohp']['file_mc']
    file2_mc = input_files['ohp']['file2_mc']
    hist_prefix = 'ppmm'
    hist2_prefix = 'ppmm'
    plotname_prefix = 'bec_3d_fit_ohp'
    comment_prefix = 'OHP, variable binning;'
if mix:
    file_data = input_files['mix']['file_data']
    file2_data = input_files['mix']['file2_data']
    file_mc = input_files['mix']['file_mc']
    file2_mc = input_files['mix']['file2_mc']
    hist_prefix = 'ppmm'
    hist2_prefix = 'ppmm'
    plotname_prefix = 'bec_3d_fit_mix'
    comment_prefix = 'MIX, variable binning;'
title = ''
rej_from = args.rej_from 
rej_to = args.rej_to
q_min = 20
q_max = -1
proj_min = 1
proj_max = 10
plotname_suffix = '_proj_%i_%i_rej_%i_%i' % (proj_min, proj_max, rej_from, rej_to)
c2_index = args.c2_index
plotname_prefix = plotname_prefix + '_c2_%i' % c2_index


ns = (2, 10), (11, 20), (21, 30), (31, 40), \
     (41, 50), (51, 60), (61, 70), (71, 80), \
     (81, 100), (101, 0)
# before was (0, 10)
# in 1D BEC mults = (2, 10), (11, 20), (21, 30), (31, 40), \
#                   (41, 50), (51, 60), (61, 70), (71, 80), \
#                   (81, 90), (91, 100), (101, 125), \
#                   (126, 150), (151, 200), (201, 250)
kts = (100, 200), (200, 300), (300, 400), (400, 500), \
      (500, 600), (600, 700), (700, 1000), (1000, 1500)

if test:
    ns = (51, 60),


args_prefix = ' --file-data ' + file_data
if file_mc is not None:
    args_prefix += ' --file-mc ' + file_mc
if file2_data is not None:
    args_prefix += ' --file2-data ' + file2_data
if file2_mc is not None:
    args_prefix += ' --file2-mc ' + file2_mc
args_prefix += ' --rej-from %d' % rej_from
args_prefix += ' --rej-to %d' % rej_to
args_prefix += ' --q-min %d' % q_min
args_prefix += ' --q-max %d' % q_max
args_prefix += ' --proj-min %d' % proj_min
args_prefix += ' --proj-max %d' % proj_max
args_prefix += ' --title "%s"' % title
args_prefix += ' --c2-index %d' % c2_index
if not args.dont_use_eps:
    args_prefix += " --use-eps"
if args.fix_c0:
    args_prefix += " --fix-c0"
if not args.dont_show_c0:
    args_prefix += " --show-c0"
if args.enlarge_uncert:
    args_prefix += " --enlarge-uncert"

if mult and not kt:
    plotname_prefix += '_mult'
    add_table_line(plotname_prefix)
    add_csv_line(plotname_prefix)
    for n1, n2 in ns:
        hist_suffix = hist_type + '_n_%i_%i' % (n1, n2)
        plotname = plotname_prefix + '_%i_%i' % (n1, n2) + plotname_suffix
        comment = comment_prefix + '%i < n_{ch} #leq %i' % (n1, n2)

        args = args_prefix
        args += ' --hist-data ' + hist_prefix + hist_suffix
        args += ' --hist2-data ' + hist2_prefix + hist_suffix
        args += ' --hist-mc ' + hist_prefix + hist_suffix
        args += ' --hist2-mc ' + hist2_prefix + hist_suffix
        args += ' --plot-name ' + plotname
        args += ' --comment "' + comment + '"'

        print('./fit' + args)
        return_val = system('./fit' + args)
        if return_val > 0:
            exit(1)
        if n2 == 0:
            n2 = 150
        if n1 == 0:
            n1 = 2
        add_table_line(plotname_prefix, '%i--%i ' % (n1, n2))
        add_csv_line(plotname_prefix, '%i,%i,' % (n1, n2))

if kt and not mult:
    plotname_prefix += '_kt'
    add_table_line(plotname_prefix)
    add_csv_line(plotname_prefix)
    for kt1, kt2 in kts:
        hist_suffix = hist_type + '_kt_%i_%i' % (kt1, kt2)
        plotname = plotname_prefix + '_%i_%i' % (kt1, kt2) + plotname_suffix
        comment = '%i #leq k_{T} < %i [MeV]' % (kt1, kt2)

        args = args_prefix
        args += ' --hist-data ' + hist_prefix + hist_suffix
        args += ' --hist2-data ' + hist2_prefix + hist_suffix
        args += ' --hist-mc ' + hist_prefix + hist_suffix
        args += ' --hist2-mc ' + hist2_prefix + hist_suffix
        args += ' --plot-name ' + plotname
        args += ' --comment "' + comment + '"'

        print('./fit' + args)
        system('./fit' + args)
        if return_val > 0:
            exit(1)
        add_table_line(plotname_prefix, '%i--%i ' % (kt1, kt2))
        add_csv_line(plotname_prefix, '%i,%i, ' % (kt1, kt2))

if mult and kt:
    for n1, n2 in ns:
        for kt1, kt2 in kts:
            hist_suffix = hist_type + '_kt_%i_%i_n_%i_%i' % (kt1, kt2, n1, n2)
            plotname = plotname_prefix
            plotname += '_kt_%i_%i_mult_%i_%i' % (kt1, kt2, n1, n2)
            plotname += plotname_suffix
            comment = '%i #leq n_{obs} < %i;' % (n1, n2)
            comment += '%i #leq k_{T} < %i [MeV]' % (kt1, kt2)

            args = args_prefix
            args += ' --hist-data ' + hist_prefix + hist_suffix
            args += ' --hist2-data ' + hist2_prefix + hist_suffix
            args += ' --hist-mc ' + hist_prefix + hist_suffix
            args += ' --hist2-mc ' + hist2_prefix + hist_suffix
            args += ' --plot-name ' + plotname
            args += ' --comment "' + comment + '"'

            print('./fit' + args)
            system('./fit' + args)
            if return_val > 0:
                exit(1)
            add_table_line(plotname_prefix + '_mult_%i_%i_kt' % (n1, n2), '%i--%i ' % (kt1, kt2))
            add_csv_line(plotname_prefix + '_mult_%i_%i_kt' % (n1, n2), '%i,%i, ' % (kt1, kt2))
