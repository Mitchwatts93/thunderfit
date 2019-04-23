import logging

from os.path import join
from os import mkdir
from json import dump as j_dump
from json import load as j_load
from dill import dump as d_dump
from dill import load as d_load
import pandas as pd
from numpy import vstack
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from . import normalisation


#### tools
def save_thunder(obj, path, filename='thunder.d'):
    d_dump(obj, open(join(path, filename), 'wb'))

def load_thunder(path):
    obj = d_load(open(path, 'rb'))
    return obj

def save_plot(plot, path='.', figname='figure.png'):
    plot.savefig(join(path, figname), transparent=True, format='svg')

def save_fit_report(obj, path, filename="report.json"):
    j_dump(obj, open(join(path, filename), 'w'))

def find_closest_indices(list1, list2):
    try:
        list_of_matching_indices = [min(range(len(list1)), key=lambda i: abs(list1[i] - cent))
                                for cent in list2]
    except ValueError:
        print('this dataset has no peaks!')
        return
    return list_of_matching_indices

def normalise_all(y_bg_rem, bg, y_raw):
    y_data_bg_rm, (mean_y_data, std_dev) = normalisation.svn(y_bg_rem) # normalise the data
    background, _ = normalisation.svn(bg, mean_y_data, std_dev) #normalise with data from bg subtracted data
    y_data_norm, _ = normalisation.svn(y_raw, mean_y_data, std_dev) #normalise with data from bg subtracted data

    return y_data_bg_rm, background, y_data_norm
#### tools

### user inputs and loading etc
def tightness_setter(tightness):
    tight_dict = {}
    tight_dict['centre_bounds'] = 1
    tight_dict['width_bounds'] = (5, 2)
    tight_dict['amps_bounds'] = (2, 2)

    if tightness == 'low':
        tight_dict['centre_bounds'] = 10
        tight_dict['width_bounds'] = (20, 3)
        tight_dict['amps_bounds'] = (5, 3)
    elif tightness == "med":
        pass
    elif tightness == 'high':
        tight_dict['centre_bounds'] = 0.5
        tight_dict['width_bounds'] = (2, 1)
        tight_dict['amps_bounds'] = (1.2, 1.2)

    else:
        logging.warning(
            'The tightness defined was incorrect format, use low, med or high. Using default med settings')

    return tight_dict

def load_data(datapath, x_ind, y_ind, e_ind=None, LOGGER=logging.getLogger()):
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

    if e_ind: # if we have specified this column then we use it, otherwise just x and y
        assert (len(data.columns) >= 2), "You have specified an e_ind but there are less than 3 columns in the data"
        e_data = data[e_ind].values
    else:
        e_data = None

    data.dropna() # drop any rows with NaN etc in them

    x_data = data[x_ind].values
    y_data = data[y_ind].values

    return x_data, y_data, e_data

def map_unique_coords(x_data, y_data, x_coords, y_coords):
    data = vstack((x_coords, y_coords, x_data, y_data)).transpose() # now have columns as the data
    df = pd.DataFrame(data=data, columns=['x_coords', 'y_coords', 'x_data', 'y_data'])
    unique_dict = dict(tuple(df.groupby(['x_coords', 'y_coords']))) # get a dictionary of the unique values for
                                        # coordinates (as tuples of (x,y)) and then the whole df rows for these values

    x_data, y_data, x_coords, y_coords = [], [], [], []
    for key in unique_dict.keys():
        x_data_ = unique_dict[key]['x_data'].values # get the x_data
        x_data.append(x_data_)
        y_data_ = unique_dict[key]['y_data'].values
        y_data.append(y_data_)
        x_coords.append(key[0])
        y_coords.append(key[1])

    return x_data, y_data, x_coords, y_coords

def parse_param_file(filepath='./params.txt'):
    """
    parse a params file which we assume is a dictionary
    :param filepath: str: path to params file
    :return: dictionary of paramters
    """
    # maybe use json loads if you end up writing parameter files non-manually

    with open(filepath, 'r') as f:
        arguments = j_load(f)
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
    arguments['x_ind'] = arg.x_ind
    arguments['y_ind'] = arg.y_ind
    arguments['e_ind'] = arg.e_ind
    arguments['datapath'] = arg.datapath
    arguments['no_peaks'] = arg.no_peaks
    arguments['background'] = arg.background
    arguments['scarf_params'] = arg.scarf_params
    arguments['peak_types'] = arg.peak_types
    arguments['peak_centres'] = arg.peak_centres
    arguments['peak_widths'] = arg.peak_widths
    arguments['peak_amps'] = arg.peak_amps
    arguments['tightness'] = arg.tightness
    arguments['bounds'] = arg.bounds

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
        mkdir(f'{dirname}')
    except FileExistsError as e:
        dirname = make_dir(f'{dirname}_new', i + 1)
        if i == 1:
            print(e, f'. So I named the file: {dirname}')
        return dirname
    return dirname

def clip_data(x_data, y_data):

    clip_left, clip_right = 0, len(x_data) - 1
    while True:
        fig, ax = plt.subplots()
        ax.plot(x_data[clip_left:clip_right], y_data[clip_left:clip_right])
        print(f"Removing background, please type two x values seperated by a space for the clips. \n"
              f"Current values are: {x_data[clip_left]}, {x_data[clip_right]}. \n"
              f"PLEASE MAKE SURE YOU ENTER IN THE SAME ORDER AS HERE. i.e. if first value is larger than right then the "
              f"first value will be the large x_clip second small")
        plt.show(block=True)
        ans = input("If you are happy with the clips type y. If not then please type a new pair of values ")
        if ans == 'y':
            break
        else:
            try:
                ans = ans.split(' ')
                if len(ans) != 2:
                    raise ValueError("The tuple was more than two elements long")
                clip_left = float(ans[0])
                clip_left = find_closest_indices(list(x_data), [clip_left])[0]
                clip_right = float(ans[1])
                clip_right = find_closest_indices(list(x_data), [clip_right])[0]
            except ValueError:
                print("You entered an incorrect answer! Trying again...")

    plt.close()
    return clip_left, clip_right

def apply_func(key_kwargs_, func):
    key = key_kwargs_[0]
    kwargs_ = key_kwargs_[1]
    val = func(*kwargs_)
    return key, val
####


