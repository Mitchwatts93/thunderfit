import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
import os
import json
import math

from scipy.signal import find_peaks
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import leastsq, least_squares
import numpy as np
import pandas as pd
import tables
from typing import Dict, Union
import peakutils
import copy


import funcs


class BagmanThunder():
    """
    bagman object with all the methods we love inside it. Name generated using WuTang Clan name generator.
    """
    def __init__(self, input):
        self.input: Union[BagmanThunder, Dict] = input
        self.data: pd.DataFrame = pd.DataFrame([])  # this is what we will create later
        self.x_ind: int = 0
        self.y_ind: int = 1
        self.e_ind: Union[int, None] = None
        self.x_label: str = 'x_axis'
        self.y_label: str = 'y_axis'
        self.e_label: Union[str, None] = None
        self.datapath: str = 'data.txt'
        self.fit_data: Dict = {'yfit': None, 'background': None, 'peak_types': [], 'peak_centres': [], 'peak_widths':[],
                         'peak_amps': [], 'chisq': None, 'free_params': None, 'p_value':None}

        if isinstance(input, BagmanThunder):  # if only pass one but its already a spec1d object then just use that
            self.overwrite_bagman(input)  # add all the details in depending on args
        elif isinstance(input, dict):
            self.create_bagman(input)  # add all the details in depending on args
        else:
            raise TypeError('Cannot convert input to BagmanThunder object')

        self.data = self.load_data(self.datapath, self.x_ind, self.y_ind, self.x_label, self.y_label, self.e_ind,
                                   self.e_label) # load the data


    def overwrite_bagman(self, inp):
        bag = inp
        self.x_ind = bag.x_ind
        self.y_ind = bag.y_ind
        self.e_ind = bag.e_ind
        self.x_label = bag.x_label
        self.y_label = bag.y_label
        self.datapath = bag.datapath
        self.fit_data = bag.fit_data

    def create_bagman(self, inp: Dict):
        """
        Used to create a bagman object given different input types
        :param args: a,b,c depending on type of input and
        :return: None, we modify the object unless a spec1d object is passed, in which case we return that
        """
        self.x_label = inp.get('x_label', self.x_label)  # if the key exists then set as that, otherwise make it
        self.y_label = inp.get('y_label', self.y_label)
        self.e_label = inp.get('e_label', self.e_label)
        try:
            self.x_ind = inp['x_ind']
            self.y_ind = inp['y_ind']
            self.e_ind = inp['e_ind']
            self.datapath = inp['datapath']
        except KeyError as e:
            LOGGER.info(f"KeyError: Missing field in the data dictionary: {e}")
        self.fit_data = inp.get('fit_data', self.fit_data)

    @staticmethod
    def load_data(datapath, x_ind, y_ind, x_label, y_label, e_ind=None, e_label=None):
        """
        load in data as a pandas df - save by modifying self.data, use object params to load
        :return: None
        """
        if '.h5' in datapath: # if the data is already stored as a pandas df
            store = pd.HDFStore(datapath)
            keys = store.keys()
            if len(keys) > 1:
                LOGGER.warning("Too many keys in the hdfstore, will assume all should be concated")
                LOGGER.warning("not sure this concat works yet")
                data = store.concat([store[key] for key in keys]) # not sure this will work! concat all keys dfs together
            else:
                data = store[keys[0]] # if only one key then we use it as the datafile
        else: # its a txt or csv file
            data = pd.read_csv(datapath, header=None, sep='\t') # load in, works for .txt and .csv
            # this needs to be made more flexible/user defined

        col_ind = [x_ind, y_ind]
        col_lab = [x_label, y_label]
        if e_ind: # if we have specified this column then we use it, otherwise just x and y
            assert (len(data.columns) >= 2), "You have specified an e_ind but there are less than 3 columns in the data"
            col_ind.append(e_ind)
            col_lab.append(e_label)
        data = data[col_ind]  # keep only these columns, don't want to waste memory
        data.columns = col_lab   # rename the columns
        dropped = data.dropna() # drop any rows with NaN etc in them
        return data

    @staticmethod
    def peak_finder(data, prominence):
        peaks, _ = find_peaks(data, prominence) # find the peak positions in the data
        peaks = list(peaks) # convert to a list
        return peaks

    def remove_background(self, data):
        params = np.array([0.01, 10 ** 5])
        bounds = [np.array([0.001, 10 ** 5]), np.array([0.1, 10 ** 9])]
        baseline_values = least_squares(self.residual_baseline, params[:], args=(data.values,),
                                  bounds=bounds)
        #baseline_values = peakutils.baseline(data, deg=20) # use peakutils to find the baseline

        p, lam = baseline_values['x']
        baseline_values = self.baseline_als(data.values, lam, p, niter=10)
        return baseline_values

    def residual_baseline(self, params, y):
        p, lam = params
        niter = 10
        baseline = self.baseline_als(y, lam, p, niter)
        residual = y - baseline
        return residual

    @staticmethod
    def baseline_als(y, lam, p, niter=10):
        L = len(y)
        D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
        w = np.ones(L)
        if niter < 1:
            raise ValueError("n iter is too small!")
        for i in range(niter):
            W = sparse.spdiags(w, 0, L, L)
            Z = W + lam * D.dot(D.transpose())
            z = spsolve(Z, w * y)
            w = p * (y > z) + (1 - p) * (y < z)
        return z

    def fit_peaks(self, fit_data, data, x_label, y_label):
        num_peaks, params, bounds, peak_funcs = self.make_params(data, fit_data, x_label, y_label)
        optimised = least_squares(self.peak_resid, params[:], args=(data[x_label].values, data[y_label].values,
                                                                    num_peaks, peak_funcs, fit_data['background']),
                                                                    bounds=bounds)
        peak_sum = self.peak_func(data[x_label].values, num_peaks, peak_funcs, *optimised['x']) # the fitted curve
        peak_sum += fit_data['background']
        return optimised, peak_sum

    @staticmethod
    def make_params(data, fit_data, x_label, y_label):
        peak_funcs = []
        peak_cents = []
        peak_widths = []
        peak_amps = []
        for i, type in enumerate(fit_data['peak_types']):
            peak_funcs.append(eval(f'funcs.{type}'))  # append the functions included in funcs
            peak_cents.append(fit_data['peak_centres'][i])
            peak_widths.append(fit_data['peak_widths'][i])
            peak_amps.append(fit_data['peak_amps'][i])
        num_peaks = len(peak_funcs)
        params = np.array(peak_cents + peak_widths + peak_amps)

        cent_lb = data[x_label].min()
        cent_ub = data[x_label].max()
        l_cent_bounds = [cent_lb for _ in peak_cents]
        u_cent_bounds = [cent_ub for _ in peak_cents]
        width_ub = (data[x_label].max() - data[x_label].min()) / 10 # width of data
        l_width_bounds = [0 for _ in peak_widths]
        u_width_bounds = [width_ub for _ in peak_widths]
        amps_lb = data[y_label].min()
        amps_ub = data[y_label].max()
        l_amp_bounds = [amps_lb for _ in peak_amps]
        u_amp_bounds = [amps_ub for _ in peak_amps]

        l_bounds = np.array(l_cent_bounds + l_width_bounds + l_amp_bounds)
        u_bounds = np.array(u_cent_bounds + u_width_bounds + u_amp_bounds)
        bounds = [l_bounds, u_bounds]

        return num_peaks, params, bounds, peak_funcs

    def peak_resid(self, params, x_data, y_data, num_peaks, peak_funcs, background):
        predicted_y = self.peak_func(x_data, num_peaks, peak_funcs, *params) # params is in format:p peak funcs,
                                                                            # peak cents, peak widths, peak amps
        residual = y_data - predicted_y + background
        return residual

    @staticmethod
    def peak_func(x_data, num_peaks, peak_funcs, *params): # args is in format: peak cents, peak widths, peak amps
        peak_centres = params[0:num_peaks]
        peak_widths = params[num_peaks:2 * num_peaks]
        peak_amps = params[2 * num_peaks:3 * num_peaks]

        peak_sum = x_data.copy() # copy so we don't change it
        peak_sum[:] = 0 # then set all the values to zero
        for i, fun in enumerate(peak_funcs):
            cent = peak_centres[i]
            width = peak_widths[i]
            amp = peak_amps[i]
            peak_sum += fun(x_data, cent, width, amp) # add together the pd dfs of y values at all x for all of the functions

        return peak_sum


def main(arguments):
    bagman = BagmanThunder(copy.deepcopy(arguments)) # load object

    if not bagman.fit_data['peak_types']:  # i.e. if no peaks were specified, then we detect the centres
        y_data = bagman.data[bagman.y_label]
        largest_y = (y_data.max() - y_data.min())  # what is the largest y value
        prominence = 10 ** (math.floor(math.log10(largest_y)))  # set the prominence of a peak to be the order of
        # magnitude of the largest peak
        # TODO: fix the way to guess prominence, order of mag is probably too big
        peak_centres = bagman.peak_finder(bagman.data[bagman.y_label], prominence)  # run peak finder here
        bagman.fit_data['peak_centres'] = peak_centres
        bagman.fit_data['peak_types'] = ['gaussian' for _ in peak_centres]  # we assume all the types are gaussian
        y_peaks = bagman.data[bagman.data[bagman.x_label] == bagman.fit_data['peak_centres']][bagman.y_label]
        import ipdb
        ipdb.set_trace()
        bagman.fit_data['peak_amps'] = list(y_peaks) # all peak amps are the order of mag of largest y
        x_data = bagman.data[bagman.x_label]
        def_width = (x_data.max() - x_data.min()) / 1000
        bagman.fit_data['peak_widths'] = [def_width for _ in peak_centres]

    if bagman.fit_data['background'] is None: # then determine the background
        bagman.fit_data['background'] = bagman.remove_background(bagman.data[bagman.y_label])
    elif bagman.fit_data['background'] == 'no': # then user doesn't want to make a background
        LOGGER.warning("not using a background, this may prevent algorithm from converging if a background is present")
        bagman.fit_data['background'] = np.array([0 for _ in bagman.data[bagman.y_label]])  # set the background as 0 everywhere
    elif not isinstance(bagman.fit_data['background'], np.ndarray):
        raise TypeError('the background passed is in the incorrect format, please pass as type np array')
    else:  # then it is the correct type and has been passed, so check the length
        assert len(bagman.fit_data['background']) == len(bagman.data[bagman.y_label]), "the background generated or passed" \
                                                                                 "is of incorrect length"

    # now fit peaks
    peak_info, peak_sum = bagman.fit_peaks(bagman.fit_data, bagman.data, bagman.x_label, bagman.y_label)
    fit_info = {'bagman':bagman, 'peak_info':peak_info, 'peak_sum':peak_sum}
    return fit_info



# for saving and parsing
def parse_param_file(filepath='./params.txt'):
    """
    parse a params file which we assume is a dictionary
    :param filepath: str: path to params file
    :return: dictionary of paramters
    """
    # maybe use json loads if you end up writing parameter files non-manually

    with open(filepath, 'r') as f:
        arguments = json.load(f)
        f.close()

    # TODO: add some checks to user passed data
    return arguments


def parse_args(arg):
    """
    convert argparse arguments into a dictionary for consistency later
    :param arg: argparse parsed args
    :return: dictionary of parameters
    """
    arguments = {}
    arguments['x_label'] = arg.x_label
    arguments['y_label'] = arg.y_label
    arguments['e_label'] = arg.y_label
    arguments['x_ind'] = arg.x_ind
    arguments['y_ind'] = arg.y_ind
    arguments['e_ind'] = arg.e_ind
    arguments['datapath'] = arg.datapath
    arguments['fit_data'] = arg.fit_data

    # TODO: add some checks to user passed data

    return arguments


def make_dir(dirname, i=1):
    """
    function to make a directory, recursively adding _new if that name already exists
    :param dirname: str: name of directory to create
    :param i: the run number we are on
    :return: str: the directory name which was available, and all subsequent data should be saved in
    """
    try:
        os.mkdir(f'{dirname}')
    except FileExistsError as e:
        dirname = make_dir(f'{dirname}_new', i+1)
        if i ==1:
            print(e, f' so I named the file: {dirname}')
        return dirname
    return dirname


if __name__ == '__main__':
    # i.e. called from bash
    import argparse

    parser = argparse.ArgumentParser(
        description='Compute the cumulative counts at various radii given the beam data and background files'
    )
    parser.add_argument('--param_file_path', type=str, default='./params.txt',
                        help='input filepath to param file, if you want to use it')
    parser.add_argument('--x_label', type=str, default='x_axis',
                        help='the label for independent variables')
    parser.add_argument('--y_label', type=str, default='y_axis',
                        help='the label for dependent variables')
    parser.add_argument('--e_label', type=str, default='y_error',
                        help='the label for uncertainties in y')
    parser.add_argument('--x_ind', type=int, default=0,
                        help='the column in data which is the independent data')
    parser.add_argument('--y_ind', type=int, default='',
                        help='the column in data which is the dependent data')
    parser.add_argument('--e_ind', type=Union[int, None], default=None,
                        help='the column in data which is the independent data uncertainties')
    parser.add_argument('--datapath', type=str, default='./data.txt',
                        help='relative path to the datafile from where python script is called')
    parser.add_argument('--fit_data', type=Dict, default={'yfit': None, 'background': None, 'peak_types': [],
                            'peak_centres': [], 'peak_widths':[], 'peak_amps': [], 'chisq': None, 'free_params': None,
                                                          'p_value':None},
                        help='the fit data as specified in the bagmanthunder __init__')
    args = parser.parse_args()  # this allows us to now use them all

    if args.param_file_path: # if there is a params file then use it
        LOGGER.info('Using params file')
        arguments = parse_param_file(args.param_file_path) # parse it
    else:
        print('not using params file')
        arguments = parse_args(args) # else use argparse but put in dictionary form

    dirname = make_dir('analysed')  # make a dict for the processed data to be saved in

    main(arguments)