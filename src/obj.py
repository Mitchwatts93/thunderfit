import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
import os
import json

from scipy.signal import find_peaks_cwt as peak_find
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import least_squares
import numpy as np
import pandas as pd
from matplotlib import axes
import matplotlib.pyplot as plt

from typing import Dict, Union
import copy

import lmfit
from lmfit import models







class BagmanThunder():
    """
    bagman object with all the methods we love inside it. Name generated using WuTang Clan name generator.
    """
    def __init__(self, input):
        self.input: Union[BagmanThunder, Dict] = input
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

        self.user_params: Dict = {'yfit': None, 'background': None, 'peak_types': [], 'peak_centres': [], 'peak_widths':[],
                         'peak_amps': [], 'chisq': None, 'free_params': None, 'p_value':None, 'bounds' : {'centers':None,
                    'widths':None,
                    'amps':None}}

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
        self.user_params = bag.user_params

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
        self.user_params = inp.get('fit_params', self.user_params)

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

    ##### peak finding
    def peaks_unspecified(self):
        x_data = self.data_bg_rm[self.x_label]

        width_ranges = [50, len(x_data) / 2]  # these are index widths TODO make this a variable...
        peak_centres_indices = self.peak_finder(self.data_bg_rm[self.y_label],
                                                  width_ranges)  # run peak finder here
        self.user_params['peak_centres'] = x_data[peak_centres_indices].values  # these are the indices of the centres

        if not self.user_params['peak_types']:
            self.user_params['peak_types'] = ['GaussianModel' for _ in
                                             peak_centres_indices]  # we assume all the types are gaussian
        elif len(self.user_params['peak_types']) != len(peak_centres_indices):
            self.user_params['peak_types'] = [self.user_params['peak_types'][0] for _ in peak_centres_indices]

        y_peaks = self.data_bg_rm[self.y_label][peak_centres_indices]  # get the y values from the indices
        self.user_params['peak_amps'] = list(y_peaks)  # all peak amps are the order of mag of largest y

        width = x_data.max() - x_data.min()
        self.user_params['peak_widths'] = [(width / 10) * np.random.random() for _ in self.user_params['peak_centres']]

    @staticmethod
    def peak_finder(data, width_range):
        peaks = peak_find(data, widths=width_range) # find the peak positions in the data
        peaks = list(peaks) # convert to a list
        return peaks
    ##### peak finding end

    ##### background
    def background_finder(self):
        if self.user_params['background'] is None:
            self.user_params['background'] = self.find_background(self.data[self.y_label])

            self.data_bg_rm[self.y_label] = self.data[self.y_label] - self.user_params[
                'background']  # subtract background from the data
            self.data_bg_rm[self.x_label] = self.data[self.x_label]
        elif self.user_params['background'] == 'no':  # then user doesn't want to make a background
            LOGGER.warning(
                "not using a background, this may prevent algorithm from converging if a background is present")
            self.user_params['background'] = np.array(
                [0 for _ in self.data[self.y_label]])  # set the background as 0 everywhere

            self.data_bg_rm[self.y_label] = self.data[self.y_label] - 0 # subtract background from the data
            self.data_bg_rm[self.x_label] = self.data[self.x_label]
        elif not isinstance(self.user_params['background'], np.ndarray):
            raise TypeError('the background passed is in the incorrect format, please pass as type np array')
        else:  # then it is the correct type and has been passed, so check the length
            assert len(self.user_params['background']) == len(
                self.data[self.y_label]), "the background generated or passed" \
                                              "is of incorrect length"
            self.data_bg_rm[self.y_label] = self.data[self.y_label] - self.user_params[
                'background']  # subtract background from the data
            self.data_bg_rm[self.x_label] = self.data[self.x_label]

    def find_background(self, data):
        params = np.array([0.01, 10 ** 5])
        bounds = [np.array([0.001, 10 ** 5]), np.array([0.1, 10 ** 9])]
        baseline_values = least_squares(self.residual_baseline, params[:], args=(data.values,),
                                  bounds=bounds)

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
    ##### background end

    ##### peak fitting
    def fit_peaks(self):
        self.user_params = self.make_bounds(self.data_bg_rm, self.user_params, self.x_label, self.y_label)
        self.specs = self.build_specs(self.data_bg_rm[self.x_label].values, self.data_bg_rm[self.y_label].values, self.user_params)

        self.model, self.peak_params = self.generate_model(self.specs)

        self.peaks = self.model.fit(self.specs['y_bg_rm'], self.peak_params, x=self.specs['x_bg_rm'])
        return self.peaks

    def make_bounds(self, data_bg_rm, user_params, x_label, y_label):
        if user_params['bounds']['centers'] is None:
            peak_cents = user_params['peak_centres']
            cent_lb = data_bg_rm[x_label].min()
            cent_ub = data_bg_rm[x_label].max()
            l_cent_bounds = [cent_lb for _ in peak_cents]
            u_cent_bounds = [cent_ub for _ in peak_cents]
            cent_bounds = list(zip(l_cent_bounds, u_cent_bounds))
            user_params['bounds']['centers'] = cent_bounds

        if user_params['bounds']['widths'] is None:
            peak_widths = user_params['peak_widths']
            width_ub = (data_bg_rm[x_label].max() - data_bg_rm[x_label].min()) # width of data
            l_width_bounds = [1e-6 for _ in peak_widths]
            u_width_bounds = [width_ub for _ in peak_widths] # TODO maybe use upper bound as e.g. 2x height
            width_bounds = list(zip(l_width_bounds, u_width_bounds))
            user_params['bounds']['widths'] = width_bounds

        if user_params['bounds']['amps'] is None:
            peak_amps = user_params['peak_amps']
            amps_lb = data_bg_rm[y_label].min()
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

    ##### plotting
    #todo fix the assertions in these
    # for all of these take a figure as optional input so can be plotted on same axis
    @staticmethod
    def plot_data(x, y, ax=False, params=('r-',)):
        if ax:
            #assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
                                                                " an axes object"
        else:
            fig, ax = plt.subplots()

        if params:
            assert isinstance(params, tuple), "invalid plot params passed"

        ax.plot(x, y, *params)
        return ax

    @staticmethod
    def plot_fits(x, peaks, ax=False):
        if ax:
            #assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
                                                                " an axes object"
        else:
            fig, ax = plt.subplots()

        for i, peak in enumerate(peaks):
            ax.plot(x, peaks[peak])
        return ax

    @staticmethod
    def plot_background(x, background_data, ax=False, params=('b-',)):
        if ax:
            #assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
                                                                " an axes object"
        else:
            fig, ax = plt.subplots()

        if params:
            assert isinstance(params, tuple), "invalid plot params passed"

        ax.plot(x, background_data, *params)
        return ax

    @staticmethod
    def plot_fit_sum(x, peak_sum, ax=False, params=('k-',)): # option of including background
        if ax:
            #assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
                                                                " an axes object"
        else:
            fig, ax = plt.subplots()

        if params:
            assert isinstance(params, tuple), "invalid plot params passed"

        ax.plot(x, peak_sum, *params)
        return ax

    def plot_all(self):
        ax = self.plot_data(self.data[self.x_label], self.data[self.y_label])
        ax = self.plot_fits(self.data[self.x_label], self.peaks.eval_components(), ax)
        ax = self.plot_background(self.data[self.x_label], self.user_params['background'], ax)
        ax = self.plot_fit_sum(self.data[self.x_label], self.peaks.eval(), ax)

        plt.show()
    ##### plotting end


# TODO move bounds making into a new function in main
def main(arguments):
    bagman = BagmanThunder(copy.deepcopy(arguments)) # load object

    bagman.background_finder() # then determine the background

    if not bagman.user_params['peak_centres']:  # i.e. if no peak centres were specified, then we detect the centres
        bagman.peaks_unspecified()

    # now fit peaks
    bagman.fit_peaks()

    return bagman


if __name__ == '__main__':

    ##### for saving and parsing
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
                print(e, f' so I named the file: {dirname}')
            return dirname
        return dirname
    #####


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
    parser.add_argument('--user_params', type=Dict, default={'yfit': None, 'background': None, 'peak_types': [],
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