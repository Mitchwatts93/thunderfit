import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
import os
import time
import numpy as np
from typing import Union, Dict, List
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
import ast
from scipy import sparse

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
    file_name = os.path.basename(ast.literal_eval(arguments['datapath'])[0]) # the name
    file_name = file_name.split('.')[0] # the name of the file
    dirname = utili.make_dir(f'{file_name}_analysed_{curr_time}')  # make a dict for the processed data to be saved in

    bag = multi_obj.main(arguments) # create a Thunder object

    ###### clip the data if weird edges
    if arguments.get('clip_data', False):
        first_thunder = bag.thunder_bag[sorted(bag.thunder_bag.keys())[0]]
        clip_left, clip_right = utili.clip_data(first_thunder.x_data, first_thunder.y_data)
        for thund in bag.thunder_bag.values():
            setattr(thund, 'x_data', thund.x_data[clip_left:clip_right])
            setattr(thund, 'y_data', thund.y_data[clip_left:clip_right])

    ###### cosmic ray removal goes here
    ####################################################################################################################

    ###### remove background
    if arguments.get('bg_first_only', False):
        # add step to find bg parameters for first one and use for the rest.
        first_thunder = bag.thunder_bag[sorted(bag.thunder_bag.keys())[0]]
        _, _, params = bg_remove.background_finder(first_thunder.x_data, first_thunder.y_data, first_thunder.background, first_thunder.scarf_params)
        [param.pop('b', None) for param in params]
        for thund in bag.thunder_bag.values():
            setattr(thund, 'scarf_params', params) # set all the values to this
    bag.bag_iterator(bag.thunder_bag, bg_remove.background_finder, ('x_data', 'y_data',
                                                                   'background', 'scarf_params'), ('background', 'y_data_bg_rm', 'params')) # determine the background

    ###### normalisation
    if args.normalise:
        pool = bag.bag_iterator(bag.thunder_bag, normalise_all, ('y_data_bg_rm', 'background', 'y_data'), ('y_data_bg_rm', 'background', 'y_data_norm'), pool)
        pool.restart()

    ###### find peaks
    if arguments.get('peakf_first_only', False):
        # add step to find bg parameters for first one and use for the rest.
        first_thunder = bag.thunder_bag[sorted(bag.thunder_bag.keys())[0]]
        no_peaks, peak_centres, peak_amps, peak_widths, peak_types, prominence = \
            peak_finding.peaks_unspecified(first_thunder.x_data, first_thunder.y_data_bg_rm, first_thunder.no_peaks,
                                           first_thunder.peak_centres, first_thunder.peak_amps, first_thunder.peak_widths,
                                           first_thunder.peak_types)
        for thund in bag.thunder_bag.values(): # set these first values for all of them
            setattr(thund, 'no_peaks', no_peaks)  # set values
            setattr(thund, 'peak_centres', peak_centres)  # set values
            setattr(thund, 'peak_amps', peak_amps)  # set values
            setattr(thund, 'peak_widths', peak_widths)  # set values
            setattr(thund, 'peak_types', peak_types)  # set values
            setattr(thund, 'prominence', prominence)  # set values
    else:
        bag.bag_iterator(bag.thunder_bag, peak_finding.peaks_unspecified, ('x_data', 'y_data_bg_rm', 'no_peaks',
                                                  'peak_centres', 'peak_amps', 'peak_widths',
                                                  'peak_types'), ('no_peaks', 'peak_centres', 'peak_amps', 'peak_widths', 'peak_types', 'prominence')) # find peaks/use them if supplied

    ###### find bounds
    if arguments.get('peakf_first_only', False):
        # add step to find bg parameters for first one and use for the rest.
        first_thunder = bag.thunder_bag[sorted(bag.thunder_bag.keys())[0]]
        bounds = peak_fitting.make_bounds(first_thunder.tightness, first_thunder.no_peaks, first_thunder.bounds,
                                          first_thunder.peak_widths, first_thunder.peak_centres, first_thunder.peak_amps)
        for thund in bag.thunder_bag.values():  # set these first values for all of them
            setattr(thund, 'bounds', bounds)  # set values
    else:
        bag.bag_iterator(bag.thunder_bag, peak_fitting.make_bounds, ('tightness', 'no_peaks', 'bounds', 'peak_widths',
                                              'peak_centres', 'peak_amps'), ('bounds')) # make bounds

    ###### fit peaks
    bag.bag_iterator(bag.thunder_bag, peak_fitting.fit_peaks, ('x_data', 'y_data_bg_rm', 'peak_types', 'peak_centres',
                               'peak_amps', 'peak_widths', 'bounds'), ('specs', 'model', 'peak_params', 'peaks')) # fit peaks

    ###### store important values of fits in dicts
    # store all the peak parameters in a dictionary, so the keys are e.g. sigma, center, amplitude, and the values are
    # dictionaries with keys as the run number with values as lists of values for all the peaks for that run
    # this will for now assume the same types of peak for all fits!
    fit_params = {}
    first_thunder = bag.thunder_bag[sorted(bag.thunder_bag.keys())[0]]
    params = list(first_thunder.peak_params.keys())[: len(first_thunder.peak_params) // first_thunder.no_peaks] # what are the peak params for the first peak
    params = [param.split('_')[1] for param in params] # keep only the type of param
    for param in params:
        fit_params[param] = {} # e.g. 'center'
        for key in bag.thunder_bag.keys():
            fit_details = bag.thunder_bag[key].peak_params
            fit_details = [fit_details[key_] for key_ in fit_details.keys() if param in key_]
            fit_params[param][key] = fit_details

    ###### fetch stats etc
    stats = {'chisq':{}, 'reduced_chi_sq':{}, 'free_params':{}}
    for key in bag.thunder_bag.keys():
        chisq = bag.thunder_bag[key].peaks.chisqr
        reduced_chi_sq = bag.thunder_bag[key].peaks.redchi
        free_params = round(chisq / reduced_chi_sq)
        stats['chisq'][key] = chisq
        stats['reduced_chi_sq'][key] = reduced_chi_sq
        stats['free_params'][key] = free_params

    ###### plot map scan
    def shift_map_matr(coordinates_array):
        coordinates_array[:,0] = coordinates_array[:,0] - min(coordinates_array[:,0])
        coordinates_array[:, 1] = coordinates_array[:, 1] - min(coordinates_array[:, 1])
        return coordinates_array
    def map_scan_plot(coordinates, values):
        no_fits = len(list(values.values())[0])
        figs = []
        data = []
        for i in range(no_fits):
            X = []
            Y = []
            Z = []
            for key in values.keys():
                x, y = coordinates[key]
                z = values[key][i]
                X.append(x)
                Y.append(y)
                Z.append(z)
            x_step = np.unique(X)[1] - np.unique(X)[0]
            xx = (X/x_step).astype(int)
            y_step = np.unique(Y)[1] - np.unique(Y)[0]
            yy = (Y / y_step).astype(int)
            data_ = sparse.coo_matrix((Z, (xx, yy))).toarray()
            data.append(data_)
            f = plt.figure()
            ax = plt.gca()
            im = plt.imshow(data_, cmap='magma', extent=[min(X), max(X), max(Y), min(Y)], vmin=data_.min(), vmax=data_.max())
            plt.xlabel('x coordinates')
            plt.ylabel('y coordinates')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            plt.colorbar(im, cax=cax)
            f.tight_layout()
            figs.append(f)
        return figs

    coordinates_array = np.array(list(bag.coordinates.values())) # convert coordinates for each point into an array
    coordinates_array = shift_map_matr(coordinates_array) # shift so that each coordinate starts at 0, not normalised
    for i, key in enumerate(bag.coordinates):
        bag.coordinates[key] = coordinates_array.tolist()[i] # reassign in the correct format
    p = ''
    plot = None
    while True:
        ans = input("making map scans, please input which property you would like to scan. options are:"
              f"\n {[p_ for p_ in fit_params.keys()]}. or type e to exit")
        if ans == 'y':
            break
        else:
            try:
                p = ans
                plot = map_scan_plot(bag.coordinates, fit_params[p])
                for pt in plot:
                    pt.suptitle(f'{p}_heatmap')
            except KeyError:
                p = ''
                print('wrong answer entered, trying again!')
            try:
                for i, pt in enumerate(plot):
                    utili.save_plot(pt, path=dirname, figname=f"{p}_{i}.svg")
            except AttributeError:
                print("Tried to save plot but there is no plot yet! Somethuing wen't wrong in making the plot")

    ###### put here some code for cluster analysis and pca
    ####################################################################################################################


    ############## figure out what to save
    # save a plot of the figure and the bag object
    utili.save_fit_report(stats, path=dirname, filename=f"{file_name}_report.json")
    utili.save_fit_report(fit_params, path=dirname, filename=f"{file_name}_peak_info.json")
    utili.save_thunder(bag, path=dirname, filename=f"{file_name}.p")