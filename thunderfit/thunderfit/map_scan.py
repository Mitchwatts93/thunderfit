import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
from os.path import basename
from ast import literal_eval

from . import utilities as utili
from . import multi_obj
from .background import background_removal as bg_remove
from . import peak_finding
from . import peak_fitting
from . import parsing
from . import map_scan_tools

def main():
    args = parsing.parse_user_args()

    arguments, LOGGER = parsing.using_user_args(args, LOGGER)

    # fix this. can be moved to then end near dir creation if you sort out the list requirement in multiobj.
    try:  # a user can pass in a list of filenames or just one
        file_name = basename(literal_eval(arguments['datapath'])[0])
    except SyntaxError:  # assume its just a string and not a list passed
        file_name = None
        arguments['datapath'] = f"['{arguments['datapath']}',]"  # as this is what multiobj needs

    bag = multi_obj.main(arguments) # create a Thunder object

    first = map_scan_tools.choose_spectrum(bag, arguments) # choose which spectrum to base everything off of

    ###### clip the data if weird edges
    if arguments.get('clip_data', False):
        bag = map_scan_tools.clip_data(bag, first)

    ###### cosmic ray removal goes here
    ####################################################################################################################

    ###### remove background
    if arguments.get('bg_first_only', False):
        bag = map_scan_tools.bg_param_setter(bag, first)
    bag.bag_iterator(bag.thunder_bag, bg_remove.background_finder, ('x_data', 'y_data',
                                                                   'background', 'scarf_params'), ('background', 'y_data_bg_rm', 'params')) # determine the background

    ###### normalisation
    if args.normalise:
        bag.bag_iterator(bag.thunder_bag, utili.normalise_all, ('y_data_bg_rm', 'background', 'y_data'), ('y_data_bg_rm', 'background', 'y_data_norm'))

    ###### find peaks
    if arguments.get('peakf_first_only', False):
        bag = map_scan_tools.peak_info_setter(bag, first)
    bag.bag_iterator(bag.thunder_bag, peak_finding.find_peak_details, ('x_data', 'y_data_bg_rm', 'no_peaks',
                                                  'peak_centres', 'peak_amps', 'peak_widths',
                                                  'peak_types'), ('no_peaks', 'peak_centres', 'peak_amps', 'peak_widths', 'peak_types', 'prominence')) # find peaks/use them if supplied

    ###### find bounds
    if arguments.get('bounds_first_only', False):
        bag = map_scan_tools.bound_setter(bag, first)
    else:
        bounds = {'amps': False, 'centers': False, 'widths': False}
        bag = map_scan_tools.bound_setter(bag, first, bounds)
    bag.bag_iterator(bag.thunder_bag, peak_finding.make_bounds, ('tightness', 'no_peaks', 'bounds', 'peak_widths',
                                              'peak_centres', 'peak_amps'), ('bounds',)) # make bounds

    ###### fit peaks
    bag.bag_iterator(bag.thunder_bag, peak_fitting.fit_peaks, ('x_data', 'y_data_bg_rm', 'peak_types', 'peak_centres',
                               'peak_amps', 'peak_widths', 'bounds'), ('specs', 'model', 'peak_params', 'peaks')) # fit peaks

    ##### fit params dictionary
    # store all the peak parameters in a dictionary, so the keys are e.g. sigma, center, amplitude, and the values are
    # dictionaries with keys as the run number with values as lists of values for all the peaks for that run
    # this will for now assume the same types of peak for all fits!
    fit_params = map_scan_tools.make_fit_params(first, bag)

    ###### fetch stats etc
    stats = map_scan_tools.get_fit_stats(bag)

    ###### make directory to save everything in
    file_name, dirname = parsing.make_user_files(arguments, file_name)

    ###### plot map scan
    map_scan_tools.plot_map_scan(bag, fit_params, dirname)

    # save individual plots for each of the failed fits
    map_scan_tools.save_failed_plots(bag, dirname)

    ###### put here some code for cluster analysis and pca
    ####################################################################################################################

    # save the bag object and it reports
    utili.save_fit_report(stats, path=dirname, filename=f"{file_name}_report.json")
    utili.save_fit_report(fit_params, path=dirname, filename=f"{file_name}_peak_info.json")
    utili.save_thunder(bag, path=dirname, filename=f"{file_name}.d")
