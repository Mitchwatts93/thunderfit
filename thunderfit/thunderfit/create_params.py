import matplotlib.pyplot as plt
from os.path import basename
from ast import literal_eval
import numpy as np

from . import utilities
from . import parsing
from . import multi_obj
from . import thundobj
from . import peak_finding

def main():
    args = parsing.parse_user_args()

    arguments = parsing.using_user_args(args)

    if arguments['map']:
        try:  # a user can pass in a list of filenames or just one
            file_name = basename(literal_eval(arguments['datapath'])[0])
        except SyntaxError:  # assume its just a string and not a list passed
            file_name = None
            log_name = arguments['datapath']
            arguments['datapath'] = f"['{arguments['datapath']}',]"  # as this is what multiobj needs
        bag = multi_obj.main(arguments)  # create a Thunder object

        bag.choose_spectrum() # get the user to choose which spectrum to create params from
        thund = bag.thunder_bag[bag.first]
    else:
        thund = thundobj.main(arguments)  # create a Thunder object


    thund.y_data = utilities.sharpening_routine(thund.x_data, thund.y_data)

    peak_info, _ = peak_finding.interactive_peakfinder(thund.x_data, thund.y_data, type='user')
    utilities.save_fit_report(peak_info, './', 'generated_params.txt')

