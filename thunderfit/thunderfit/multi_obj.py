from ast import literal_eval
from copy import deepcopy
from glob import glob

import matplotlib
import matplotlib.pyplot as plt
from numpy import round
from pandas.errors import ParserError
from tqdm import tqdm

matplotlib.use('TkAgg')
from typing import Union
import logging

from .thundobj import Thunder
from . import utilities as utili
from .background import background_removal as bg_remove
from . import peak_finding
from . import peak_fitting


############## NOT IMPLEMENTED!
# TODO
# make option of passing in many params files - one for each data file

class ThunderBag():

    def __init__(self, input):
        # initialise everything first
        self.thunder_bag: {} = {}
        self.coordinates: {} = {}
        self.first: str = ''
        self.stats: {} = {}
        self.fit_params: {} = {}
        self.x_ind: Union[None, int] = None
        self.y_ind: Union[None, int] = None
        self.e_ind: Union[None, int] = None
        self.img_path: Union[None, str] = None
        self.map: Union[None, str] = None
        self.datapath: Union[None, str] = None
        self.x_coord_ind: Union[None, int] = None
        self.y_coord_ind: Union[None, int] = None

        if isinstance(input, Thunder):  # if only pass one but its already a thunder object then just use that
            self.thunder_bag[0] = Thunder(input)  # add all the details in depending on args
        elif isinstance(input, dict):
            self.create_bag(input)  # add all the details in depending on args
        else:
            raise TypeError('Cannot convert input to ThunderBag object')

    def create_bag(self, inp):
        logging.debug('creating bag object')
        self.x_ind = inp.get('x_ind', 2)
        self.y_ind = inp.get('y_ind', 3)
        self.e_ind = inp.get('e_ind', None)
        self.img_path = inp.get('imgpath', None)
        self.coordinates = inp.get('coords', {})

        self.map = inp.get('map', None)  # if user passes map as True then the file will be treated as a map file

        data_paths = inp.get('datapath', None)
        self.datapath = literal_eval(data_paths)  # this is a bit dangerous!!!!!

        for i, data in tqdm(enumerate(self.datapath)):
            if len(self.datapath):
                prefix = f'{i}_'  # if more than one datapath then we name them with i_j
            if isinstance(data, Thunder):
                self.thunder_bag[i] = data
            elif isinstance(data, str):
                # then read the data file
                if self.map == True:
                    prefix = ''
                    self.x_coord_ind, self.y_coord_ind = inp.get('x_coord_ind', 0), inp.get('y_coord_ind', 1)
                    map_path = glob(data)[0]  # save the filepath to the mapscan as self.map for later
                    x_data, y_data, x_coords, y_coords = self.read_map(map_path, self.x_ind, self.y_ind,
                                                                       self.x_coord_ind, self.y_coord_ind)

                    for j in range(len(x_data)):  # go through the list of x_data
                        x_data_, y_data_ = x_data[j], y_data[j]  # the x and y data for each coordinate set
                        self.thunder_bag[f'{prefix}{j}'] = Thunder(inp, x_data=x_data_,
                                                                   y_data=y_data_)  # make a thunder obj
                        # with this data
                        x_coords_, y_coords_ = x_coords[j], y_coords[j]
                        self.coordinates[f'{prefix}{j}'] = (
                            x_coords_, y_coords_)  # for each i we will have a list of tuples of x and y coords
                elif '*' in data:
                    filematches = glob(data)
                    for j, file in enumerate(filematches):
                        try:
                            self.thunder_bag[f'{prefix}{j}'] = self.create_thunder(file,
                                                                                   inp)  # make a thunder object for
                            # each file
                        except ParserError as e:
                            logging.warning(f"A Thunder object could not be created for the datafile: {file}, skipping")
                else:
                    try:
                        self.thunder_bag[str(i)] = self.create_thunder(data, inp)
                    except ParserError as e:
                        logging.warning(f"A Thunder object could not be created for the datafile: {file}, skipping")
            else:
                logging.warning(f"wrong format in data list detected for {i}th element: {data}. Skipping element")
                pass

    @staticmethod
    def create_thunder(file, inp):
        logging.debug('creating thunder object')
        arguments = deepcopy(inp)
        arguments['datapath'] = file
        thund_obj = Thunder(arguments)
        return thund_obj

    @staticmethod
    def read_map(file_address, x_ind, y_ind, x_coord_ind, y_coord_ind):
        logging.debug('reading in mapscan file')
        x_data, y_data, _ = utili.load_data(file_address, x_ind, y_ind)  # load the data. note these drop nan rows but
        # does that for the whole filepath so will be consistent for data and coordinates
        x_coords, y_coords, _ = utili.load_data(file_address, x_coord_ind, y_coord_ind)  # load the coordinates
        x_data, y_data, x_coords, y_coords = utili.map_unique_coords(x_data, y_data, x_coords, y_coords)  #

        return x_data, y_data, x_coords, y_coords

    @staticmethod
    def bag_iterator(bag, func, input_args, sett_args):
        bagkeys = tqdm(bag.keys())
        bagkeys.set_description(f"Operating with: {func.__name__}, to find: {sett_args}")
        for key in bagkeys:
            thund = bag.get(key)  # bag[key]
            kwargs_ = [getattr(thund, arg) for arg in input_args]
            _, val = utili.apply_func((key, kwargs_), func)
            for i, arg in enumerate(sett_args):
                try:
                    setattr(thund, arg, val[i])
                except KeyError as e:
                    if isinstance(val, dict):
                        setattr(thund, arg, val)
                    else:
                        print(f'Weird KeyError encountered: {e}')

    def choose_spectrum(self):
        logging.debug('choosing which thunder object will be the user specified data for bg etc')
        # then we have to choose which spectrum we want
        first = next(iter(self.thunder_bag.keys()))  # changed from list to iter
        while True:
            try:
                first_thunder = self.thunder_bag[first]
                fig, ax = plt.subplots()
                ax.plot(getattr(first_thunder, 'x_data'), getattr(first_thunder, 'y_data'))
                print(f"Need a decision on which plot is representitive of data, the following is for index {first}")
                plt.show(block=True)
                ans = input("If you are happy with using this data file, type y, otherwise enter a new index")
                if ans == 'y':
                    break
                else:
                    try:
                        first = str(ans)
                    except ValueError:
                        print("You entered an incorrect answer! Trying again...")
            except KeyError:
                print('incorrect key, please enter a lower index value')
                first = next(iter(self.thunder_bag.keys()))
        self.first = first

    def clip_data(self):
        logging.debug('clipping data based on user specified plot')
        first_thunder = self.thunder_bag[self.first]
        clip_left, clip_right = utili.clip_data(getattr(first_thunder, 'x_data'), getattr(first_thunder, 'y_data'))
        for thund in self.thunder_bag.values():
            setattr(thund, 'x_data', getattr(thund, 'x_data')[clip_left:clip_right])
            setattr(thund, 'y_data', getattr(thund, 'y_data')[clip_left:clip_right])

    def bg_param_setter(self):
        logging.debug('setting backgrounds for all based on background of user specified plot')
        # add step to find bg parameters for first one and use for the rest.
        first_thunder = self.thunder_bag[self.first]
        if isinstance(getattr(first_thunder, 'background'), str) and getattr(first_thunder, 'background') == "SCARF":
            _, _, params = bg_remove.background_finder(getattr(first_thunder, 'x_data'),
                                                       getattr(first_thunder, 'y_data'),
                                                       getattr(first_thunder, 'background'),
                                                       getattr(first_thunder, 'scarf_params'))
            [param.pop('b', None) for param in params]  # we want to find b each time so don't set it for all others
            for thund in self.thunder_bag.values():
                setattr(thund, 'scarf_params', params)  # set all the values to this

    def remove_backgrounds(self):
        # do some checks
        self.bag_iterator(getattr(self, 'thunder_bag'), bg_remove.background_finder,
                         ('x_data', 'y_data', 'background', 'scarf_params'),
                         ('background', 'y_data_bg_rm', 'params'))

    def normalise_data(self):
        self.bag_iterator(getattr(self, 'thunder_bag'), utili.normalise_all, ('y_data_bg_rm', 'background', 'y_data'),
                         ('y_data_bg_rm', 'background', 'y_data_norm'))

    def find_peaks(self):
        logging.debug('setting peak no, centres and types based on user specified plot details')
        # add step to find bg parameters for first one and use for the rest.
        first_thunder = self.thunder_bag[self.first]
        no_peaks, peak_info_dict, prominence = \
            peak_finding.find_peak_details(getattr(first_thunder, 'x_data'), getattr(first_thunder, 'y_data_bg_rm'),
                                           getattr(first_thunder, 'no_peaks'), getattr(first_thunder, 'peak_finder_type', 'auto'))
        for thund in self.thunder_bag.values():  # set these first values for all of them
            setattr(thund, 'no_peaks', no_peaks)  # set values
            center_indices = utili.find_closest_indices(thund.x_data, peak_info_dict['center']) # get the indices from the x centres
            center_indices = peak_finding.match_peak_centres(center_indices, thund.y_data) # match to the peakfinding cents
            peak_centres = thund.x_data[center_indices] # convert back to x values
            peak_info_dict['center'] = peak_centres # set values
            thund.peak_info_dict = peak_info_dict

    def peaks_adj_params(self):
        for thund in self.thunder_bag.values():  # set these first values for all of them
            center_indices = utili.find_closest_indices(thund.x_data, thund.peak_info_dict['center']) # get the indices from the x centres
            center_indices = peak_finding.match_peak_centres(center_indices, thund.y_data) # match to the peakfinding cents
            peak_centres = thund.x_data[center_indices] # convert back to x values
            thund.peak_info_dict['center'] = peak_centres # set values

    def bound_setter(self, bounds=None):
        logging.debug('setting bounds based on user provided bounds or found for user specified plot')
        if not bounds:
            first_thunder = self.thunder_bag[self.first]
            bounds = peak_finding.make_bounds(getattr(first_thunder, 'x_data'), getattr(first_thunder, 'y_data'), getattr(first_thunder, 'no_peaks'),
                                              first_thunder.peak_info_dict)
        for thund in self.thunder_bag.values():  # set these first values for all of them
            setattr(thund, 'bounds', bounds)  # set values

    def fit_peaks(self):
        self.bag_iterator(getattr(self, 'thunder_bag'), peak_fitting.fit_peaks,
                         ('x_data', 'y_data_bg_rm', 'peak_info_dict', 'bounds', 'method', 'tol'),
                         ('specs', 'model', 'peak_params', 'peaks'))  # fit peaks

    def make_fit_params(self):
        logging.debug('generating fit params')
        fit_params = {}
        first_thunder = self.thunder_bag.get(self.first)
        params = list(getattr(first_thunder, 'peak_params').keys())[
                 : len(
                     getattr(first_thunder, 'peak_params')) // getattr(first_thunder, 'no_peaks')]  # what are the peak params for the first peak
        params = [param.split('_')[1] for param in params]  # keep only the type of param
        for param in params:
            fit_params[param] = {}  # e.g. 'center'
            for key in self.thunder_bag.keys():
                fit_details = getattr(self.thunder_bag.get(key), 'peak_params')
                fit_details = [fit_details.get(key_) for key_ in fit_details.keys() if param in key_]
                fit_params[param][key] = fit_details
        self.fit_params = fit_params

    def get_fit_stats(self):
        logging.debug('generating fit stats')
        stats = {'chisq': {}, 'reduced_chi_sq': {}, 'free_params': {}}
        for key, thund in self.thunder_bag.items():
            chisq = getattr(getattr(thund, 'peaks'), 'chisqr')
            reduced_chi_sq = getattr(getattr(thund, 'peaks'), 'redchi')
            free_params = round(chisq / reduced_chi_sq)
            stats['chisq'][key] = chisq
            stats['reduced_chi_sq'][key] = reduced_chi_sq
            stats['free_params'][key] = free_params
        self.stats = stats

    def save_failed_plots(self, dirname):
        logging.debug('saving failed plots')
        for key, thund in self.thunder_bag.items():
            if not getattr(getattr(thund, 'peaks'), 'success'):
                thund.plot_all()
                utili.save_plot(thund.plot, path=dirname,
                                figname=f"failed_plot_{key}_at_position_{self.coordinates.get(key)}.svg")
                thund.plot.close()  # close so memory is conserved.


def main(arguments):
    bag = ThunderBag(deepcopy(arguments))  # load object
    return bag
