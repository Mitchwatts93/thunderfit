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
        self.peak_types: Union[None, list] = []
        self.peak_centres: Union[None, list] = []
        self.peak_widths: Union[None, list] = []
        self.peak_amps: Union[None, list] = []
        self.tightness: str = "med"
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
        self.amp_bounds: bool = False

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

        self.tightness = utili.tightness_setter(self.tightness)

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
        self.peak_types = thun.peak_types
        self.peak_centres = thun.peak_centres
        self.peak_widths = thun.peak_widths
        self.peak_amps = thun.peak_amps
        self.tightness = thun.tightness
        self.bounds = thun.bounds

        self.method = thun.method
        self.tol = thun.tol
        self.amp_bounds = thun.amp_bounds

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
        self.peak_types = inp.get('peak_types', self.peak_types)
        self.peak_centres = inp.get('peak_centres', self.peak_centres)
        self.peak_widths = inp.get('peak_widths', self.peak_widths)
        self.peak_amps = inp.get('peak_amps', self.peak_amps)
        self.tightness = inp.get('tightness', self.tightness)
        self.bounds = inp.get('bounds', self.bounds)
        self.method = inp.get('method', self.method)
        self.tol = inp.get('tol', self.tol)
        self.amp_bounds = inp.get('amp_bound', self.amp_bounds)

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
        self.fit_report = {mod_no: {} for mod_no in range(len(self.peak_types))}

        ## total fit data
        self.fit_report['chi_sq'] = self.chi_sq
        self.fit_report['free_params'] = self.free_params
        self.fit_report['p_value'] = 'not implemented'

        ## individual parameter data
        param_info = {"center": "centers", "amplitude": "amps", "sigma": "widths", "fwhm": False, "height": False}
        for parameter, param_obj in self.peaks.params.items():
            model_no = int(findall(r'\d+', parameter)[0])
            param_type = param_info[get_close_matches(parameter, param_info.keys())[0]]

            if param_type:
                value = param_obj.value
                err = param_obj.stderr
                type = self.peak_types[model_no]
                bounds = self.bounds[param_type][model_no]

                fit_info = {"value": value,
                            "stderr": err,
                            "peak_type": type,
                            "bounds": bounds}

                self.fit_report[model_no][param_type] = fit_info


def main(arguments):
    thunder = Thunder(deepcopy(arguments))  # load object
    return thunder
