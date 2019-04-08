import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
import os
import time

from typing import Union, Dict, List

from . import utilities as utili
from . import multi_obj
from .background import background_removal as bg_remove
from . import peak_finding
from . import peak_fitting

def str_or_none(value):
    try:
        return str(value)
    except:
        return None

def parse_user_args():
    import argparse

    parser = argparse.ArgumentParser(
        description='fit peaks and background to the given data given a set of parameter'
    )
    parser.add_argument('--param_file_path', type=str_or_none, default=None,
                        help='input filepath to param file, if you want to use it')
    parser.add_argument('--datapath', type=str_or_none, default=None,
                        help='relative path to the datafile from where python script is called')
    parser.add_argument('--x_ind', type=str_or_none, default=None,
                        help='the column in data which is the independent data')
    parser.add_argument('--y_ind', type=str_or_none, default=None,
                        help='the column in data which is the dependent data')
    parser.add_argument('--e_ind', type=Union[int, None], default=None,
                        help='NOT IMPLEMENTED YET. the column in data which is the independent data uncertainties')
    parser.add_argument('--no_peaks', type=int, default=None,
                        help='the number of peaks you would like fitted. If you have specified bounds or peak infomation'
                             '\n e.g. centres then please make sure this is the same length as that list')
    parser.add_argument('--background', type=str, default="SCARF",
                        help="The stype of background you'd like to fit. 'SCARF' is a rolling ball solgay_filter "
                             "background subtraction. \n 'OLD' uses a soon-to-be not implemented numerical method"
                             "which doesn't work too well. \n 'no' specifies that you would like no background fitted."
                             "NOT IMPLEMENTED YET: An np array of background can also be passed by passing the path to the file, but please note that it must be the same"
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
        arguments['datapath'] = args.datapath
    else:
        print('not using params file')
        arguments = utili.parse_args(args)  # else use argparse but put in dictionary form

    curr_time = time.strftime('%d_%m_%Y_%l:%M%p') #name directory with the current time
    dirname = utili.make_dir(f'analysed_{curr_time}')  # make a dict for the processed data to be saved in

    bag = multi_obj.main(arguments) # create a Thunder object


    bag.bag_iterator(bag.thunder_bag, bg_remove.background_finder, ('x_data', 'y_data',
                                                                   'background', 'scarf_params'), ('background', 'y_data_bg_rm', 'params')) # determine the background

    if args.normalise:
        bag.bag_iterator(bag.thunder_bag, normalise_all, ('y_data_bg_rm', 'background', 'y_data'), ('y_data_bg_rm', 'background', 'y_data_norm'))

    import ipdb
    ipdb.set_trace()
    bag.bag_iterator(bag.thunder_bag, peak_finding.peaks_unspecified, ('x_data', 'y_data_bg_rm', 'no_peaks',
                                                  'peak_centres', 'peak_amps', 'peak_widths',
                                                  'peak_types'), ('no_peaks', 'peak_centres', 'peak_amps', 'peak_widths', 'peak_types', 'prominence')) # find peaks/use them if supplied

    import ipdb
    ipdb.set_trace()
    bag.bag_iterator(bag.thunder_bag, peak_fitting.make_bounds, 'tightness', 'no_peaks', 'bounds', 'peak_widths',
                                              'peak_centres', 'peak_amps') # make bounds

    import ipdb
    ipdb.set_trace()
    bag.bag_iterator(bag.thunder_bag, peak_fitting.fit_peaks, 'x_data', 'y_data_bg_rm', 'peak_types', 'peak_centres',
                               'peak_amps', 'peak_widths', 'bounds') # fit peaks

    import ipdb
    ipdb.set_trace()
    bag.chi_sq = bag.peaks.chisqr # set the stats from the fits
    reduced_chi_sq = bag.peaks.redchi
    bag.free_params = round(bag.chi_sq / reduced_chi_sq)

    bag.plot_all() # plot the data in full and save as an object
    bag.gen_fit_report() # generate a fit report

    # save a plot of the figure and the thunder object
    dataname = os.path.basename(arguments['datapath'])
    utili.save_plot(bag.plot, path=dirname, figname=f"{dataname}.svg")
    utili.save_thunder(bag, path=dirname, filename=f"{dataname}.p")
    utili.save_fit_report(bag.fit_report, path=dirname, filename=f"{dataname}_report.json")