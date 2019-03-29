from obj import *

# here we want to fit multiple datasets.

class ThunderBag():




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
    args = parser.parse_args()  # this allows us to now use them all

    if args.param_file_path: # if there is a params file then use it
        LOGGER.info('Using params file')
        arguments = parse_param_file(args.param_file_path) # parse it
    else:
        print('not using params file')
        arguments = parse_args(args) # else use argparse but put in dictionary form



    # now save things!