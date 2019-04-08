from obj import *
import glob
import copy
import pandas
from tqdm import tqdm

############## NOT IMPLEMENTED!
# TODO
# make option of passing in many params files - one for each data file

class ThunderBag():

    def __init__(self, input):
        # initialise everything first
        self.thunder_bag: {} = {}

        if isinstance(input, Thunder):  # if only pass one but its already a thunder object then just use that
            self.overwrite_thunder(input)  # add all the details in depending on args
        elif isinstance(input, dict):
            self.create_bag(input)  # add all the details in depending on args
        else:
            raise TypeError('Cannot convert input to Thunder object')

    def overwite_thunder(self, inp):
        thun = inp
        self.x_ind = thun.x_ind
        self.y_ind = thun.y_ind
        self.e_ind = thun.e_ind
        self.x_label = thun.x_label
        self.y_label = thun.y_label
        self.datapath = thun.datapath

    def create_bag(self, inp):
        self.x_label = inp['x_label']
        self.y_label = inp['y_label']
        self.e_label = inp['e_label']
        self.x_ind =  inp['x_ind']
        self. y_ind = inp['y_ind']
        self.e_ind = inp['e_ind']
        self.img_path = inp['imgpath']

        assert isinstance(inp['datapath'], list), "Wrong format for datapath, should be a list"
        self.datapath = inp['datapath']

        for i, data in tqdm(enumerate(self.datapath)):
            if isinstance(data, Thunder):
                self.thunder_bag[i] = data
            elif isinstance(data, str):
                # then read the data file
                if '*' in data:
                    filematches = glob.glob(data)
                    for j, file in enumerate(filematches):
                        try:
                            self.thunder_bag[f'{i}_{j}'] = self.create_thunder(file, inp) # make a thunder object for each file
                        except pandas.errors.ParserError as e:
                            logging.warn(f"A Thunder object could not be created for the datafile: {file}, skipping")
                else:
                    try:
                        self.thunder_bag[str(i)] = self.create_thunder(data, inp)
                    except pandas.errors.ParserError as e:
                        logging.warn(f"A Thunder object could not be created for the datafile: {file}, skipping")
            else:
                logging.warn(f"wrong format in data list detected for {i}th element: {data}. Skipping element")
                pass

    @staticmethod
    def create_thunder(file, inp):
        arguments = copy.deepcopy(inp)
        arguments['datapath'] = file
        thund_obj = Thunder(arguments)
        return thund_obj

    @staticmethod
    def fit_bag(bag_dict):
        for baglabel, thund in tqdm(bag_dict.items()):
            thund.background_finder()  # then determine the background
            specified_dict = peak_details(thund.user_params)
            thund.peaks_unspecified(specified_dict)

            # now fit peaks
            thund.fit_peaks()
            #thund.plot_all()
            thund.fit_report()

        return bag_dict

def main(arguments):

    bag = ThunderBag(copy.deepcopy(arguments)) # load object

    bag.fit_bag(bag.thunder_bag)
    import ipdb
    ipdb.set_trace()


    return bag


def parse_param_file(filepath='./bag_params.txt'):
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
