from obj import *

# here we want to fit multiple datasets.

class ThunderBag():

    def __init__(self, input):
        # initialise everything first

        if isinstance(input, Thunder):  # if only pass one but its already a thunder object then just use that
            self.overwrite_thunder(input)  # add all the details in depending on args
        elif isinstance(input, dict):
            self.create_thunder(input)  # add all the details in depending on args
        else:
            raise TypeError('Cannot convert input to Thunder object')

        self.data = self.load_data(self.datapath, self.x_ind, self.y_ind, self.x_label, self.y_label, self.e_ind,
                                   self.e_label) # load the data




def main(arguments):
    bag = ThunderBag(copy.deepcopy(arguments)) # load object



    return bagman


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
        arguments['imgpath'] = arg.imgpath
        arguments['info_path'] = arg.info_path

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
        description='fit peaks and background to the given data given a set of parameter'
    )
    parser.add_argument('--param_file_path', type=str, default='./bag_params.txt',
                        help='input filepath to param file, if you want to use it, if you use this no other inputs will '
                             'be parsed')
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
    parser.add_argument('--datapath', type=str, default=['./data.txt'],
                        help='relative paths to the datafiles from where python script is called, a list of paths')
    parser.add_argument('--imgpath', type=str, default=['./data.txt'],
                        help='relative path to the image file for the scans - if it exists')
    parser.add_argument('--info_path', type=str, default=['./data.txt'],
                        help='relative path to the scan infomation - if it exists')
    parser.add_argument('--user_params', type=Dict, default={'yfit': None, 'background': None, 'peak_types': [],
                                                             'peak_centres': [], 'peak_widths': [], 'peak_amps': [],
                                                             'chisq': None, 'free_params': None,
                                                             'p_value': None, 'tightness': None},
                        help='the fit data as specified in the Thunder __init__')
    args = parser.parse_args()  # this allows us to now use them all

    if args.param_file_path: # if there is a params file then use it
        LOGGER.info('Using params file')
        arguments = parse_param_file(args.param_file_path) # parse it
    else:
        print('not using params file')
        arguments = parse_args(args) # else use argparse but put in dictionary form



    # now save things!