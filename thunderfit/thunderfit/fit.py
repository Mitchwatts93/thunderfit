import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
import os
import json
import dill
import operator
import difflib
import re
from typing import Dict, Union
import copy

from scipy.signal import find_peaks as peak_find
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import least_squares
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import lmfit
from lmfit import models

import scarf
import utilities as utili
import plotting
import background_removal as bg_remove


# TODO: need to fail if peak fitting doesn't work!
# fix peak finder
# add argument for number of peaks to use pick these based on some paramter e.g. prominence
# sort peaks in this order
# add try for matplotlib and set interactive background to not run if it fails

class Thunder():
    """
    thunder object with all the methods we love inside it. Name generated using WuTang Clan name generator.
    """
    def __init__(self, input):
        self.input: Union[Thunder, Dict] = input
        self.data: pd.DataFrame = pd.DataFrame([])  # this is what we will create later
        self.data_bg_rm: pd.DataFrame = pd.DataFrame([]) # later we will fill this with background remove data
        self.x_ind: int = 0
        self.y_ind: int = 1
        self.e_ind: Union[int, None] = None
        self.x_label: str = 'x_axis'
        self.y_label: str = 'y_axis'
        self.e_label: Union[str, None] = None
        self.datapath: str = 'data.txt'
        self.peaks: lmfit.model.ModelResult
        self.plot: plt = None
        self.fit_data: {} = {}

        self.user_params: Dict = {'yfit': None, 'background': None, 'peak_types': [], 'peak_centres': [], 'peak_widths':[],
                                'peak_amps': [], 'chisq': None, 'free_params': None, 'p_value':None, 'tightness':None,
                                'bounds' : {'centers':None, 'widths':None, 'amps':None}}

        if isinstance(input, Thunder):  # if only pass one but its already a thunder object then just use that
            self.overwrite_thunder(input)  # add all the details in depending on args
        elif isinstance(input, dict):
            self.create_thunder(input)  # add all the details in depending on args
        else:
            raise TypeError('Cannot convert input to Thunder object')

        self.data = utili.load_data(self.datapath, self.x_ind, self.y_ind, self.x_label, self.y_label, self.e_ind,
                                   self.e_label) # load the data

        self.tightness = utili.tightness_setter(self.user_params['tightness'])

    #### loading thunder object
    def overwrite_thunder(self, inp):
        thun = inp
        self.x_ind = thun.x_ind
        self.y_ind = thun.y_ind
        self.e_ind = thun.e_ind
        self.x_label = thun.x_label
        self.y_label = thun.y_label
        self.datapath = thun.datapath
        self.user_params = thun.user_params

    def create_thunder(self, inp: Dict):
        """
        Used to create a thunder object given different input types
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
        self.user_params = inp.get('fit_params', self.user_params)
    #### end loading

    #### background
    def background_finder(self):
        y_label = self.y_label
        x_label = self.x_label
        y_data = self.data[y_label]
        x_data = self.data[x_label]
        bg = self.user_params['background']
        data_bg_rm = self.data_bg_rm

        if bg == 'no':  # then user doesn't want to make a background
            LOGGER.warning(
                "Warning: no background specified, so not using a background,"
                " this may prevent algorithm from converging")
            bg = np.array([0 for _ in y_data])  # set the background as 0 everywhere
            data_bg_rm[y_label] = y_data # no background subtracted
            data_bg_rm[x_label] = x_data

        elif bg == 'SCARF':
            bg = np.array([0 for _ in y_data], dtype=np.float64)
            rad = 20
            b = 0
            window_length, poly_order = 51, 3
            L_sg = 0
            data_bg_rm[y_label] = y_data

            while True:
                while True:
                    D = scarf.rcf(data_bg_rm[y_label], rad)
                    fig, ax = plt.subplots()
                    ax.plot(x_data, D)
                    ax.plot(x_data, data_bg_rm[y_label])
                    print(f"SCARF background removal requires user input. Please look at the following bg with rad={rad}")
                    plt.show(block=True)
                    ans = input("If you are happy with the plot, type y. if not then please type a new rad")
                    if ans == 'y':
                        break
                    else:
                        try:
                            rad = int(ans)
                        except ValueError:
                            print("You entered an incorrect answer! Trying again...")

                L = D + b
                while True: # now estimate a baseline to add to D to get L
                    fig, ax = plt.subplots()
                    ax.plot(x_data, L)
                    ax.plot(x_data, data_bg_rm[y_label])
                    print(f"Please look at the following bg with a shift={b}")
                    plt.show(block=True)
                    ans = input("If you are happy with the plot, type y. if not then please type a new background value. \n"
                                "Please note that the background should NOT intercept the data. Ideally it would pass through"
                                "the mean of the noise for the correct bg already fit")
                    if ans == 'y':
                        L = D + b
                        break
                    else:
                        try:
                            b = int(ans)
                            L = D + b
                        except ValueError:
                            print("You entered an incorrect answer! Trying again...")

                # then apply SG filter to L
                while True:
                    try:
                        L_sg = scarf.smooth(L, window_length, poly_order)
                        fig, ax = plt.subplots()
                        ax.plot(x_data, L_sg)
                        ax.plot(x_data, data_bg_rm[y_label])
                        print(f"Please look at the following bg with Sg filter parameters (window length, polynomial order): "
                              f"{window_length}, {poly_order}")
                        plt.show(block=True)
                    except ValueError as e:
                        print(
                            "Incorrect values for window_length and poly_order have been entered. Poly order must be less than window length and window length must be odd")
                    ans = input("please enter y if you are happy with these values, or enter two integers with a space "
                                    "for window_length and poly_order")
                    if ans == 'y':
                        L = L_sg
                        break
                    else:
                        try:
                            ans = ans.split(' ')
                            if len(ans) != 2:
                                raise ValueError("The tuple was more than two elements long")
                            window_length = int(ans[0])
                            poly_order = int(ans[1])
                        except ValueError:
                            print("You entered an incorrect answer! Trying again...")

                # final question before exiting
                fig, ax = plt.subplots()
                ax.plot(x_data, L)
                ax.plot(x_data, data_bg_rm[y_label])
                print(f"Please look at the following bg with selected parameters")
                plt.show(block=True)
                ans = input("Are you happy with this bg? If yes, type y, else type n. n will restart the fitting. \n"
                            "typing repeat will add an additional bg subtraction to this one")
                if ans == 'y':
                    bg += L
                    break
                elif ans == 'n':
                    pass
                elif ans =='repeat':
                    bg += L
                    print("apply two bg removal steps, this will mean the background just specified will be removed "
                          "from the data")
                    data_bg_rm[y_label] -= L # remove the bg found here from the original data and go again
                else:
                    print("You entered an incorrect answer! Trying whole fitting routine again...")

            data_bg_rm[y_label] -= L  # subtract background from the data
            data_bg_rm[x_label] = x_data

        elif isinstance(bg, np.ndarray):
            assert len(self.user_params['background']) == len(y_data), \
                    "the background generated or passed is of incorrect length"
            data_bg_rm[y_label] = y_data - bg # subtract user supplied background from the data
            data_bg_rm[x_label] = x_data

        elif bg == 'OLD':
            bg = bg_remove.find_background(y_data, bg_remove.residual_baseline, bg_remove.baseline_als) # find a background the old way
            data_bg_rm[y_label] = y_data - bg  # subtract background from the data
            data_bg_rm[x_label] = x_data

        else:  # then it is the incorrect type
            raise TypeError('the background passed is in the incorrect format, please pass as type np array')

        #y_min = data_bg_rm[y_label].min()
        #if y_min < 0:
        #    data_bg_rm[y_label] += abs(y_min)  # then shift all the data up so no points are below zero
        #    bg -= abs(y_min)  # and lower the bg we have calculated by that shift too

        self.user_params['background'] = bg
        self.data_bg_rm = data_bg_rm
    ##### background end

    #### normalise
    @staticmethod
    def normalisation(y_data):
        """normalise using std variance normalisation"""
        mean_y_data = np.mean(y_data)
        shifted_y_data = y_data - mean_y_data
        std_dev = np.std(y_data)
        normalised_y = shifted_y_data / std_dev
        return normalised_y

    #### normalise end

    ##### peak finding
    def peaks_unspecified(self, specified_dict):
        x_data = self.data_bg_rm[self.x_label]
        # fix this up since only cents_specified works now

        if not specified_dict['cents_specified']:

            #width_ranges = [50, len(x_data) / 2]  # these are index widths TODO make this a variable...
            prominence = 1.6
            # do question asking routine for picking prominence to find peaks
            peak_info = self.peak_finder(self.data_bg_rm[self.y_label],
                                                    prominence)  # find the peak centers

            # use peak info dict and store the heights and widths of peaks
            self.user_params['peak_centres'] = x_data[peak_info['center_indices']].values  # these are the indices of the centres
            self.user_params['peak_widths'] = x_data[peak_info['right_edges']].values - x_data[peak_info['left_edges']].values
            self.user_params['peak_amps'] = peak_info['amps']
            import ipdb
            ipdb.set_trace()

            # set bounds from these too


        if not specified_dict['amps_specified']: # find a faster way to do this
            xcents = self.user_params['peak_centres'] # this is x data
            peak_centres_indices = [self.data_bg_rm[self.x_label].iloc[(self.data_bg_rm[self.x_label] - xval)
                                 .abs().argsort()[:1]].index for xval in xcents] #find the indices for these xvalues
            peak_centres_indices = [ind.tolist()[0] for ind in peak_centres_indices] # stupid pandas index type

            y_peaks = self.data_bg_rm[self.y_label][peak_centres_indices]  # get the y values from the indices
            self.user_params['peak_amps'] = list(y_peaks)  # all peak amps are the order of mag of largest y

        if not specified_dict['widths_specified']:
            width = x_data.max() - x_data.min()
            self.user_params['peak_widths'] = [(width / self.tightness['width']) * np.random.random() for _ in self.user_params['peak_centres']]

        if not specified_dict['types_specified']:
            self.user_params['peak_types'] = ['LorentzianModel' for _ in
                                              self.user_params['peak_centres']]  # we assume all the types are gaussian


        len_ord_specified = sorted(specified_dict.items(), key=operator.itemgetter(1))  # get the shortest
        len_ord_specified = filter(lambda tup: tup[1] > 0, len_ord_specified)
        try:
            shortest_specified = next(len_ord_specified)[0]  # this is the dict key with the shortest specified data

            for param in ['peak_amps', 'peak_centres', 'peak_widths', 'peak_types']:
                if len(self.user_params[param]) > specified_dict[shortest_specified]: # then we need to trim it
                    LOGGER.warning("Some of the specified peak parameters differ in length. Choosing peak paramters"
                                   "as the first n parameters where n is the length of the shortest set of parameters")
                    self.user_params[param] = self.user_params[param][:specified_dict[shortest_specified]]
        except StopIteration:
            pass

    @staticmethod
    def peak_finder(data, prominence):
        # do a routine looping through until the right number of peaks is found

        peaks, properties = peak_find(data, prominence=prominence) # find the peak positions in the data

        peaks = list(peaks) # convert to a list
        amps = list(properties['prominences']) # store the heights

        peak_info = {'center_indices':peaks, 'right_edges':list(properties['right_bases']),
                     'left_edges':list(properties['left_bases']), 'amps':amps}
        return peak_info
    ##### peak finding end

    ##### peak fitting
    def fit_peaks(self):
        self.user_params = self.make_bounds(self.data_bg_rm, self.user_params, self.y_label)
        self.specs = self.build_specs(self.data_bg_rm[self.x_label].values, self.data_bg_rm[self.y_label].values, self.user_params)

        self.model, self.peak_params = self.generate_model(self.specs)
        self.peaks = self.model.fit(self.specs['y_bg_rm'], self.peak_params, x=self.specs['x_bg_rm'])
        if not self.peaks.success:
            logging.warning('The fitting routine failed! exiting programme. Try lowering tightness settings or manually '
                         'inputting a background, peak bounds and peak info.')
        self.peak_params = self.peaks.best_values

    def make_bounds(self, data_bg_rm, user_params, y_label):
        if user_params['bounds']['centers'] is None:
            l_cent_bounds = [cent - self.tightness['centre_bounds'] *
                             user_params['peak_widths'][i] for i, cent in enumerate(user_params['peak_centres'])]
            u_cent_bounds = [cent + self.tightness['centre_bounds'] *
                             user_params['peak_widths'][i] for i, cent in enumerate(user_params['peak_centres'])]
            cent_bounds = list(zip(l_cent_bounds, u_cent_bounds))
            user_params['bounds']['centers'] = cent_bounds

        if user_params['bounds']['widths'] is None:
            peak_widths = user_params['peak_widths']
            l_width_bounds = [width / self.tightness['width_bounds'][0] for width in peak_widths]
            u_width_bounds = [width * self.tightness['width_bounds'][1] for width in peak_widths]
            width_bounds = list(zip(l_width_bounds, u_width_bounds))
            user_params['bounds']['widths'] = width_bounds

        if user_params['bounds']['amps'] is None:
            peak_amps = user_params['peak_amps']
            amps_lb = data_bg_rm[y_label].mean() # maybe change this to min
            amps_ub = data_bg_rm[y_label].max()
            l_amp_bounds = [amps_lb for _ in peak_amps]
            u_amp_bounds = [amps_ub for _ in peak_amps]
            amp_bounds = list(zip(l_amp_bounds, u_amp_bounds))
            user_params['bounds']['amps'] = amp_bounds

        # todo currently our bounds are set by the data ranges. It may make sense to define
        # narrower ranges around the peaks themselves
        return user_params

    @staticmethod
    def build_specs(x_bg_rm, y_bg_rm, user_params):
        specs = {'x_bg_rm':x_bg_rm, 'y_bg_rm':y_bg_rm,
                'model': [
                    {'type': user_params['peak_types'][i],
                    'params': {'center': user_params['peak_centres'][i], 'amp': user_params['peak_amps'][i],
                               'sigma': user_params['peak_widths'][i], 'gamma':user_params['peak_widths'][i]},
                     'bounds': {'centers': user_params['bounds']['centers'][i], 'amps': user_params['bounds']['amps'][i],
                                'widths': user_params['bounds']['widths'][i]}
                    }
                for i, _ in enumerate(user_params['peak_centres'])]
                }
        return specs

    @staticmethod
    def generate_model(spec):
        """
        https://chrisostrouchov.com/post/peak_fit_xrd_python/
        :param spec:
        :return:
        """
        composite_model = None
        params = None
        for i, basis_func in enumerate(spec['model']):
            prefix = f'm{i}_'
            model = getattr(models, basis_func['type'])(prefix=prefix)
            if basis_func['type'] in ['GaussianModel', 'LorentzianModel','VoigtModel']:
                # for now VoigtModel has gamma constrained to sigma
                w_min = basis_func['bounds']['widths'][0]
                w_max = basis_func['bounds']['widths'][1]
                x_min = basis_func['bounds']['centers'][0]
                x_max = basis_func['bounds']['centers'][1]
                y_min = basis_func['bounds']['amps'][0]
                y_max = basis_func['bounds']['amps'][1]

                model.set_param_hint('sigma', min=w_min, max=w_max)
                model.set_param_hint('center', min=x_min, max=x_max)
                model.set_param_hint('height', min=y_min, max=1.1 * y_max)
                model.set_param_hint('amplitude', min=1e-6)

                # default guess is horrible!! do not use guess()
                default_params = {
                    prefix + 'center': basis_func['params']['center'],
                    prefix + 'height': basis_func['params']['amp'],
                    prefix + 'sigma': basis_func['params']['sigma']
                }
            else:
                raise NotImplemented(f'model {basis_func["type"]} not implemented yet')

            model_params = model.make_params(**default_params, **basis_func.get('params', {}))

            if params is None: # first loop
                params = model_params
                composite_model = model
            else: # subsequent loops
                params.update(model_params)
                composite_model = composite_model + model

        return composite_model, params
    ##### peak fitting end

    def plot_all(self):
        ax = plotting.plot_fits(self.data[self.x_label], self.peaks.eval_components()) # plot each component of the model
        ax = plotting.plot_background(self.data[self.x_label], self.user_params['background'], ax) #plot the background supplied by user
        ax = plotting.plot_fit_sum(self.data[self.x_label], self.peaks.best_fit, self.user_params['background'], ax) # plot the fitted data
        try:
            ax = plotting.plot_uncertainty_curve(self.data[self.x_label], self.peaks.eval_uncertainty(sigma=3),
                                         self.peaks.best_fit, ax) #plot a band of uncertainty
        except TypeError:
            logging.warning('There are not uncertainties available for some reason - '
                         'try lowering the tightness of automatic bounds')
        ax = plotting.plot_data(self.data[self.x_label], self.data[self.y_label], ax)  # plot the raw data

        ax.minorticks_on()
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)

        self.plot = plt

    #this needs some fixing
    def fit_report(self):
        self.fit_data = {mod_no:{} for mod_no in range(len(self.user_params['peak_types']))}

        ## total fit data
        chi_sq = self.peaks.chisqr
        reduced_chi_sq = self.peaks.redchi
        free_params = round(chi_sq / reduced_chi_sq)

        ## individual parameter data
        param_info = {"center":"centers", "amplitude":"amps", "sigma":"widths", "fwhm":False, "height":False}
        for parameter, param_obj in self.peaks.params.items():
            model_no = int(re.findall(r'\d+', parameter)[0])
            param_type = param_info[difflib.get_close_matches(parameter, param_info.keys())[0]]

            if param_type:
                value =param_obj.value
                err = param_obj.stderr
                type = self.user_params['peak_types'][model_no]
                bounds = self.user_params['bounds'][param_type][model_no]

                fit_info = {"value":value,
                            "stderr":err,
                            "peak_type":type,
                            "bounds":bounds}

                self.fit_data[model_no][param_type] = fit_info


# TODO move bounds making into a new function in main
def main(arguments):
    thunder = Thunder(copy.deepcopy(arguments)) # load object

    thunder.background_finder() # then determine the background
    thunder.data_bg_rm[thunder.y_label] = thunder.normalisation(thunder.data_bg_rm[thunder.y_label]) # normalise the data

    specified_dict = utili.peak_details(thunder.user_params)
    thunder.peaks_unspecified(specified_dict)

    # now fit peaks
    thunder.fit_peaks()
    thunder.plot_all()
    thunder.fit_report()

    return thunder

if __name__ == '__main__':
    ##### for saving and parsing
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
        arguments['user_params'] = arg.user_params

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
            dirname = make_dir(f'{dirname}_new', i + 1)
            if i == 1:
                print(e, f'. So I named the file: {dirname}')
            return dirname
        return dirname
    #####

    # i.e. called from bash
    import argparse

    parser = argparse.ArgumentParser(
        description='fit peaks and background to the given data given a set of parameter'
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
    parser.add_argument('--y_ind', type=int, default=1,
                        help='the column in data which is the dependent data')
    parser.add_argument('--e_ind', type=Union[int, None], default=None,
                        help='the column in data which is the independent data uncertainties')
    parser.add_argument('--datapath', type=str, default='./data.txt',
                        help='relative path to the datafile from where python script is called')
    parser.add_argument('--user_params', type=Dict, default={'yfit': None, 'background': None, 'peak_types': [],
                            'peak_centres': [], 'peak_widths':[], 'peak_amps': [], 'chisq': None, 'free_params': None,
                                                          'p_value':None, 'tightness':None},
                        help='the fit data as specified in the Thunder __init__')
    args = parser.parse_args()  # this allows us to now use them all

    if args.param_file_path: # if there is a params file then use it
        LOGGER.info('Using params file')
        arguments = utili.parse_param_file(args.param_file_path) # parse it
    else:
        print('not using params file')
        arguments = parse_args(args) # else use argparse but put in dictionary form

    dirname = make_dir('analysed')  # make a dict for the processed data to be saved in

    thunder = main(arguments)

    # save a plot of the figure and the thunder object
    dataname = os.path.basename(arguments['datapath'])
    utili.save_plot(thunder.plot, path=dirname, figname=f"{dataname}.svg")
    utili.save_thunder(thunder, path=dirname, filename=f"{dataname}.p")
    utili.save_fit_report(thunder.fit_data, path=dirname, filename=f"{dataname}_report.json")