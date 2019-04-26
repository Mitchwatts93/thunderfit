import logging
from copy import deepcopy
from difflib import get_close_matches
from re import findall
from typing import Dict, Union

import matplotlib.pyplot as plt
from lmfit.model import ModelResult
from numpy import ndarray

from . import plotting
from . import utilities as utili
from .background import background_removal as bg_remove
from . import peak_finding
from . import peak_fitting

# TODO: need to fail if peak fitting doesn't work!

class Thunder():
    """
    thunder object with all the methods we love inside it. Name generated using WuTang Clan name generator.
    """

    def __init__(self, input, x_data=None, y_data=None, e_data=None):
        self.input: Union[Thunder, Dict] = input

        self.x_ind: int = 0
        self.y_ind: int = 1
        self.e_ind: Union[int, None] = None

        self.x_data: Union[None, ndarray] = x_data
        self.y_data: Union[None, ndarray] = y_data
        self.e_data: Union[None, ndarray] = e_data

        self.y_data_bg_rm = None
        self.y_data_norm = None

        self.datapath: str = './data.txt'

        self.no_peaks: int = 0
        self.background: str = "SCARF"  # fix this
        self.scarf_params: Union[None, Dict] = None

        self.peak_info_dict: Dict = {}
        self.peak_finder_type: str = 'auto'

        self.bounds: Union[None, Dict] = None

        self.peaks: ModelResult
        self.plot: plt = None
        self.fit_report: {} = {}
        self.peak_params: {} = {}
        self.model = None  # give type!

        self.free_params: int = 0
        self.p_value: int = 0
        self.chisq: int = 0

        self.method: str = 'leastsq'
        self.tol: float = 0.0000001

        if isinstance(input, Thunder):  # if only pass one but its already a thunder object then just use that
            self.overwrite_thunder(input)  # add all the details in depending on args
        elif isinstance(input, dict):
            self.create_thunder(input)  # add all the details in depending on args
        else:
            raise TypeError('Cannot convert input to Thunder object')

        if isinstance(self.x_data, ndarray) and isinstance(self.y_data, ndarray):
            pass  # they're already loaded as they've been passed
        else:
            self.x_data, self.y_data, self.e_data = utili.load_data(self.datapath, self.x_ind,
                                                                    self.y_ind)  # load the data

    def overwrite_thunder(self, inp):
        logging.debug('overwriting thund obj')
        thun = inp

        if thun.x_data and thun.y_data:
            self.x_data = thun.x_data
            self.y_data = thun.y_data
        else:
            self.x_ind = thun.x_ind
            self.y_ind = thun.y_ind
            self.e_ind = thun.e_ind
            self.datapath = thun.datapath

        if thun.y_data_bg_rm:
            self.y_data_bg_rm = thun.y_data_bg_rm
        if thun.y_data_norm:
            self.y_data_norm = thun.y_data_norm

        self.no_peaks = thun.no_peaks
        self.background = thun.background
        self.scarf_params = thun.scarf_params

        self.peak_info_dict = thun.peak_info_dict
        self.peak_finder_type = thun.peak_finder_type

        self.bounds = thun.bounds

        self.method = thun.method
        self.tol = thun.tol

    def create_thunder(self, inp: Dict):
        """
        Used to create a thunder object given different input types
        :param args: a,b,c depending on type of input and
        :return: None, we modify the object unless a spec1d object is passed, in which case we return that
        """
        logging.debug('creating thund obj')
        try:  # only continue if its e_ind missing
            self.e_ind = inp['e_ind']
        except KeyError as e:
            logging.info(f"KeyError: Missing field in the data dictionary: {e}")

        try:  # If its the others you need to fail here
            self.datapath = inp['datapath']
            self.x_ind = inp['x_ind']
            self.y_ind = inp['y_ind']
        except KeyError as e:
            raise KeyError(f"Missing vital information to load object: {e}")

        self.no_peaks = inp.get('no_peaks', self.no_peaks)
        self.background = inp.get('background', self.background)
        # do some check on background here to set it to an np array
        self.scarf_params = inp.get('scarf_params', self.scarf_params)

        self.peak_info_dict = inp.get('peak_info_dict', self.peak_info_dict)
        self.peak_finder_type = inp.get('peak_finder_type', self.peak_finder_type)

        self.bounds = inp.get('bounds', self.bounds)
        self.method = inp.get('method', self.method)
        self.tol = inp.get('tol', self.tol)

    def clip_data(self):
        clip_left, clip_right = utili.clip_data(getattr(self, 'x_data'), getattr(self, 'y_data'))
        setattr(self, 'x_data', getattr(self, 'x_data')[clip_left:clip_right])
        setattr(self, 'y_data', getattr(self, 'y_data')[clip_left:clip_right])

    def remove_bg(self):
        background, y_data_bg_rm, params = bg_remove.background_finder(getattr(self, 'x_data'),
                                                                       getattr(self, 'y_data'),
                                                                       getattr(self, 'background'),
                                                                       getattr(self, 'scarf_params'))
        setattr(self, 'background', background)
        setattr(self, 'y_data_bg_rm', y_data_bg_rm)
        setattr(self, 'params', params)

    def normalise(self):
        y_data_bg_rm, background, y_data_norm = utili.normalise_all('y_data_bg_rm', 'background', 'y_data')
        setattr(self, 'background', background)
        setattr(self, 'y_data_bg_rm', y_data_bg_rm)
        setattr(self, 'y_data_norm', y_data_norm)

    def find_peaks(self):
        no_peaks, peak_info_dict, prominence = \
            peak_finding.find_peak_details(getattr(self, 'x_data'), getattr(self, 'y_data_bg_rm'),
                                           getattr(self, 'no_peaks'),
                                           getattr(self, 'peak_finder_type', 'auto'))
        setattr(self, 'no_peaks', no_peaks)  # set values
        center_indices = utili.find_closest_indices(self.x_data, peak_info_dict['center'])  # get the indices from the x centres
        center_indices = peak_finding.match_peak_centres(center_indices, self.y_data)  # match to the peakfinding cents
        peak_centres = self.x_data[center_indices]  # convert back to x values
        peak_info_dict['center'] = peak_centres  # set values
        setattr(self, peak_info_dict, peak_info_dict)

    def bound_setter(self, bounds=None):
        if not bounds:
            bounds = peak_finding.make_bounds(getattr(self, 'x_data'), getattr(self, 'y_data'), getattr(self, 'no_peaks'),
                                              self.peak_info_dict)
        setattr(self, 'bounds', bounds)  # set values

    def fit_peaks(self):
        specs, model, peak_params, peaks = peak_fitting.fit_peaks(getattr(self, 'x_data'),
                                                                          getattr(self, 'y_data_bg_rm'),
                                                                          getattr(self, 'peak_info_dict'),
                                                                          getattr(self, 'bounds'),
                                                                          getattr(self, 'method'),
                                                                          getattr(self, 'tol'))
        setattr(self, 'specs', specs)
        setattr(self, 'model', model)
        setattr(self, 'peak_params', peak_params)
        setattr(self, 'peaks', peaks)

    ## plot_all and fit_report need imporovements e.g. to check which attributes exists in the object
    def plot_all(self):
        logging.debug('plotting all for thund obj')
        ax, plt = plotting.plot_fits(self.x_data, self.peaks.eval_components())  # plot each component of the model
        ax, plt = plotting.plot_background(self.x_data, self.background, ax)  # plot the background supplied by user
        ax, plt = plotting.plot_fit_sum(self.x_data, self.peaks.best_fit, self.background, ax)  # plot the fitted data
        try:
            ax, plt = plotting.plot_uncertainty_curve(self.x_data, self.peaks.eval_uncertainty(sigma=3),
                                                      self.peaks.best_fit, ax)  # plot a band of uncertainty
        except TypeError:
            logging.warning('There are not uncertainties available for some reason - '
                            'try lowering the tightness of automatic bounds')
        ax, plt = plotting.plot_data(self.x_data, self.y_data, ax)  # plot the raw data

        ax.minorticks_on()
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)

        self.plot = plt

    def gen_fit_report(self):
        logging.debug('genertaing fit report for thund obj')
        self.fit_report = {mod_no: {} for mod_no in range(len(self.peak_info_dict['peak_types']))}

        ## total fit data
        self.fit_report['chi_sq'] = self.chi_sq
        self.fit_report['free_params'] = self.free_params
        self.fit_report['p_value'] = 'not implemented'

        ## individual parameter data
        param_info = {"center": "center", "amplitude": "amplitude", "sigma": "sigma", "fwhm": False, "height": False}
        for parameter, param_obj in self.peaks.params.items():
            model_no = int(findall(r'\d+', parameter)[0])
            param_type = param_info[get_close_matches(parameter, param_info.keys())[0]]

            if param_type:
                value = param_obj.value
                err = param_obj.stderr
                type = self.peak_info_dict['peak_types'][model_no]
                bounds = self.bounds[param_type][model_no]

                fit_info = {"value": value,
                            "stderr": err,
                            "peak_type": type,
                            "bounds": bounds}

                self.fit_report[model_no][param_type] = fit_info


def main(arguments):
    thunder = Thunder(deepcopy(arguments))  # load object
    return thunder
