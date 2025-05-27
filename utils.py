'''
Utils for do_fit
'''

from os import system
from os import path
from datetime import datetime
import pandas as pd
import numpy as np


def get_filepath(plotname, extension, c2_index):
    """ Determine full path to the output file. """
    day = datetime.now().day
    day_str = str(day)
    if day < 10:
        day_str = '0' + day_str
    filepath = 'output/%s_%s_%s' % (datetime.now().year,
                                    datetime.now().strftime('%h'),
                                    day_str)
    # system('mkdir -p ' + filepath + '/c2_' + str(c2_index) + '/' + extension)
    # filepath = filepath + '/c2_' + str(c2_index) + '/' + extension + '/' + plotname + '.' + extension
    system('mkdir -p ' + filepath + '/' + c2_index + '/' + extension)
    filepath = filepath + '/' + c2_index + '/' + extension + '/' + plotname + '.' + extension

    return filepath


def remove_output_file(filename, extension, c2_index):
    """ Remove output file. """
    filepath = get_filepath(filename, extension, c2_index)

    system('rm -f ' + filepath)


def erase_line_in_table_file(tablename, extension, string, c2_index):
    """ Erase line in table file. """
    filepath = get_filepath(tablename, extension, c2_index)

    try:
        with open(filepath, 'r+') as outfile:
            buff = outfile.readlines()
            outfile.seek(0)
            for line in buff:
                if string not in line:
                    outfile.write(line)
            outfile.truncate()
    except FileNotFoundError:
        pass


def add_latex_table_line(tablename, plotname, line_prefix, c2_index, erase=True):
    """ Add line to the tex table. """
    tablefile = get_filepath(tablename, 'tex', c2_index)
    linefile = get_filepath(plotname, 'tex', c2_index)

    in_file = open(linefile, 'r')
    line = in_file.read()
    in_file.close()
    if erase:
        system('rm -f ' + linefile)

    out_file = open(tablefile, 'a+')
    out_file.write(line_prefix + line)
    out_file.close()


def add_csv_table_line(tablename, plotname, line_prefix, c2_index, erase=True):
    """ Add line to the csv table. """
    tablefile = get_filepath(tablename, 'csv', c2_index)
    linefile = get_filepath(plotname, 'csv', c2_index)

    in_file = open(linefile, 'r')
    lines = in_file.readlines()
    in_file.close()
    if erase:
        system('rm -f ' + linefile)

    if not path.isfile(tablefile):
        out_file = open(tablefile, 'w')
        out_file.write('bin_start,bin_stop,' + lines[0])
        out_file.close()

    out_file = open(tablefile, 'a+')
    out_file.write(line_prefix + lines[1])
    out_file.close()

def add_csv_table_line_averages(tablename, bin_start, bin_stop , c2_index):
    """ Add line to the csv table. """
    tablefile = get_filepath(tablename, 'csv', c2_index)

    data = pd.read_csv(tablefile, sep=',')
    # Filter the DataFrame
    filtered_data = data[(data['bin_start'] >= bin_start) & (data['bin_stop'] <= bin_stop)]
    # Calculate the average for columns that don't start with '#sigma'
    average_row = filtered_data[[col for col in filtered_data.columns if not col.startswith('#sigma')]].mean()

    # Calculate the mean of the square of the values for columns that start with '#sigma'
    for col in filtered_data.columns:
        if col.startswith('#sigma'):
            average_row[col] = np.sqrt(np.sum(filtered_data[col]**2)/len(filtered_data))

    average_row['bin_start'] = bin_start
    average_row['bin_stop'] = bin_stop

    empty_row = pd.Series(['' for _ in range(len(data.columns))], index=data.columns)
    # Append the average row to the DataFrame
    data = pd.concat([data, pd.DataFrame([empty_row])], ignore_index=True)
    data = pd.concat([data, pd.DataFrame([average_row])], ignore_index=True)

    # Write the DataFrame back to the file
    data.to_csv(tablefile, sep=',', index=False)


def add_histogram(colname, plotname, c2_index, erase=True):
    """ Add fitted histogram + fit function to common root file. """
    colfile = get_filepath(colname, 'root', c2_index)
    histfile = get_filepath(plotname, 'root', c2_index)

    if not path.isfile(colfile):
        system('mv ' + histfile + ' ' + colfile)
    else:
        system('hadd -v 0 ' + colfile + '_new ' + colfile + ' ' + histfile)
        system('mv ' + colfile + '_new ' + colfile)

    if erase:
        system('rm -f ' + histfile)
