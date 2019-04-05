import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
import os
import time
import numpy as np

from typing import Union, Dict, List
from typing import TyperVar
T = TyperVar('T')
from typing import AnyStr
Array = Union[AnyStr[T], np.ndarray[T]]

from . import utilities as utili
from . import thundobj
from .background import background_removal as bg_remove
from . import peak_finding
from . import peak_fitting


def parse_user_args():
    import argparse

    parser = argparse.ArgumentParser(
        description='fit peaks and background to the given data given a set of parameter'
    )
    parser.add_argument('--param_file_path', type=str, default='./params.txt',
                        help='input filepath to param file, if you want to use it')
    parser.add_argument('--datapath', type=str, default='./data.txt',
                        help='relative path to the datafile from where python script is called')
    parser.add_argument('--x_ind', type=int, default=0,
                        help='the column in data which is the independent data')
    parser.add_argument('--y_ind', type=int, default=1,
                        help='the column in data which is the dependent data')
    parser.add_argument('--e_ind', type=Union[int, None], default=None,
                        help='NOT IMPLEMENTED YET. the column in data which is the independent data uncertainties')
    parser.add_argument('--no_peaks', type=int, default=None,
                        help='the number of peaks you would like fitted. If you have specified bounds or peak infomation'
                             '\n e.g. centres then please make sure this is the same length as that list')
    parser.add_argument('--background', type=Array, default="SCARF",
                        help="The stype of background you'd like to fit. 'SCARF' is a rolling ball solgay_filter "
                             "background subtraction. \n 'OLD' uses a soon-to-be not implemented numerical method"
                             "which doesn't work too well. \n 'no' specifies that you would like no background fitted."
                             " An np array of background can also be passed, but please note that it must be the same"
                             "length as the datafile (once rows containing nan values have been removed).")
    parser.add_argument('--scarf_params', type=Union[None, Dict], default=None,
                        help='a dictionary (or None) of parameters for scarf algorithm. If none an interactive routine'
                             'will be used. if the dictionary is specified it should be of the form: \n'
                             '{"rad":70, "b":90, "window_length":51, "poly_order":3}'
                             '\n where window length must be odd and greater than poly_order, and all must be integers')
    parser.add_argument('--peak_types', type=Union[None, List], default=None,
                        help='a list (or none) or the types of peak to be fitted. '
                             '\n e.g. ["LorentzianModel", "GaussianModel"] as strings! they must be implemented in lmfit')
    parser.add_argument('--peak_centres', type=Union[None, List], default=None,
                        help='a list (or none) or the centres values (x values) of the peaks to be fitted. '
                             '\n e.g. [488, 365] as integers!')
    parser.add_argument('--peak_widths', type=Union[None, List], default=None,
                        help='a list (or none) or the width values (x values) of the peaks to be fitted. '
                             '\n e.g. [10, 1] as integers!')
    parser.add_argument('--peak_amps', type=Union[None, List], default=None,
                        help='a list (or none) or the amplitude values (x values) of the peaks to be fitted. '
                             '\n e.g. [10, 5] as integers!')
    parser.add_argument('--tightness', type=str, default="med",
                        help='a string indicating how tight the auto-generated bounds should be. not used if bounds'
                             'are supplied. valid values are "low", "med" and "high", any other value will cause '
                             'default to be used.')
    parser.add_argument('--bounds', type=Union[None, Dict], default=None,
                        help='a dictionary of the bounds for peak values. if none passed then will be auto-generated.'
                             '\n of the form: {"centers":[(365, 390), (283,285)],"widths":[(2, 3), (1, 4)],'
                             '"amps":[(2, 3), (1, 4)]}'
                             '\n the list elements correspond to the peaks supplied earlier, the tuple elements'
                             'correspond to the low and high bounds on that specific value')
    parser.add_argument('--normalise', type=bool, default=False,
                        help='bool, True or False for should I normalise data or not')

    args = parser.parse_args()  # this allows us to now use them all

    return args

def normalise_all(y_bg_rem, bg, y_raw):
    y_data_bg_rm, (mean_y_data, std_dev) = normalisation.svn(y_bg_rem) # normalise the data
    background, _ = normalisation.svn(bg, mean_y_data, std_dev) #normalise with data from bg subtracted data
    y_data_norm, _ = normalisation.svn(y_raw, mean_y_data, std_dev) #normalise with data from bg subtracted data

    return y_data_bg_rm, background, y_data_norm

def main():
    args = parse_user_args()

    if args.param_file_path:  # if there is a params file then use it
        LOGGER.info('Using params file')
        arguments = utili.parse_param_file(args.param_file_path)  # parse it
    else:
        print('not using params file')
        arguments = utili.parse_args(args)  # else use argparse but put in dictionary form

    curr_time = time.localtime(time.time())  # reformat this so its nice
    dirname = utili.make_dir(f'analysed_{curr_time}')  # make a dict for the processed data to be saved in

    thunder = thundobj.main(arguments) # create a Thunder object

    thunder.background, thunder.y_data_bg_rm = bg_remove.background_finder(thunder.x_data, thunder.y_data,
                                                                           thunder.background, thunder.scarf_params)
                                                                           # determine the background
    if args.normalise:
        thunder.y_data_bg_rm, thunder.background, thunder.y_data_norm = \
                                                 normalise_all(thunder.y_data_bg_rm, thunder.background, thunder.y_data)

    thunder.peak_no, thunder.peak_centres, thunder.peak_amps, thunder.peak_widths, thunder.peak_types = \
                   peak_finding.peaks_unspecified(thunder.x_data, thunder.y_data_bg_rm, thunder.peak_no,
                                                  thunder.peak_centres, thunder.peak_amps, thunder.peak_widths,
                                                  thunder.peak_types) # find peaks/use them if supplied

    thunder.bounds = peak_fitting.make_bounds(thunder.tightness, thunder.no_peaks, thunder.bounds, thunder.peak_widths,
                                              thunder.peak_centres, thunder.peak_amps) # make bounds

    thunder.specs, thunder.model, thunder.peak_params, thunder.peaks = \
        peak_fitting.fit_peaks(thunder.x_data, thunder.y_data_bg_rm, thunder.peak_types, thunder.peak_centres,
                               thunder.peak_amps, thunder.peak_widths, thunder.bounds) # fit peaks

    thunder.chi_sq = thunder.peaks.chisqr # set the stats from the fits
    reduced_chi_sq = thunder.peaks.redchi
    thunder.free_params = round(thunder.chi_sq / reduced_chi_sq)

    thunder.plot_all() # plot the data in full and save as an object
    thunder.fit_report() # generate a fit report

    # save a plot of the figure and the thunder object
    dataname = os.path.basename(arguments['datapath'])
    utili.save_plot(thunder.plot, path=dirname, figname=f"{dataname}.svg")
    utili.save_thunder(thunder, path=dirname, filename=f"{dataname}.p")
    utili.save_fit_report(thunder.fit_report, path=dirname, filename=f"{dataname}_report.json")