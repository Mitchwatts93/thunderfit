import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
import os
import json
import dill
import pandas as pd

#### tools
def save_thunder(obj, path, filename='thunder.p'):
    dill.dump(obj, open(os.path.join(path, filename), 'wb'))

def load_thunder(path):
    obj = dill.load(open(path, 'rb'))
    return obj

def save_plot(plot, path='.', figname='figure.png'):
    plot.savefig(os.path.join(path, figname), transparent=True, format='svg')

def save_fit_report(obj, path, filename="report.json"):
    json.dump(obj, open(os.path.join(path, filename), 'w'))

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
#### tools

#### parsing user params

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
####
#### loading data
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
    data.dropna() # drop any rows with NaN etc in them
    return data

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
####