import logging
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

    arguments = parsing.using_user_args(args)

    LOGGER = logging.getLogger()

    # fix this. can be moved to then end near dir creation if you sort out the list requirement in multiobj.
    try:  # a user can pass in a list of filenames or just one
        file_name = basename(literal_eval(arguments['datapath'])[0])
    except SyntaxError:  # assume its just a string and not a list passed
        file_name = None
        arguments['datapath'] = f"['{arguments['datapath']}',]"  # as this is what multiobj needs

    bag = multi_obj.main(arguments) # create a Thunder object

    bag.first = list(bag.thunder_bag.keys())[0]
    if arguments.get('clip_data', False) or arguments.get('bg_first_only', False) or arguments.get(
            'peakf_first_only', False) or arguments.get('bounds_first_only', False):
        bag.choose_spectrum() # choose which spectrum to base everything off of if user wants to use one spectra to choose parameters

    ###### clip the data if weird edges
    if arguments.get('clip_data', False):
        LOGGER.info('clipping data')
        bag.clip_data()

    ###### cosmic ray removal goes here
    ####################################################################################################################

    ###### remove background
    if arguments.get('bg_first_only', False):
        LOGGER.info('determining background conditions for all based on user guided for first')
        bag.bg_param_setter()
    LOGGER.info('removing background from data for all thunder objects')
    bag.bag_iterator(getattr(bag, 'thunder_bag'), bg_remove.background_finder, ('x_data', 'y_data',
                                                                   'background', 'scarf_params'), ('background', 'y_data_bg_rm', 'params')) # determine the background

    ###### normalisation
    if args.normalise:
        LOGGER.info('normalising data using svn normalisation')
        bag.bag_iterator(getattr(bag, 'thunder_bag'), utili.normalise_all, ('y_data_bg_rm', 'background', 'y_data'), ('y_data_bg_rm', 'background', 'y_data_norm'))

    ###### find peaks
    if arguments.get('peakf_first_only', False):
        LOGGER.info('running user guided routine to determine peak information')
        bag.peak_info_setter()
    LOGGER.info('setting peak information for all thunder objects')
    bag.bag_iterator(getattr(bag, 'thunder_bag'), peak_finding.find_peak_details, ('x_data', 'y_data_bg_rm', 'no_peaks',
                                                  'peak_centres', 'peak_amps', 'peak_widths',
                                                  'peak_types'), ('no_peaks', 'peak_centres', 'peak_amps', 'peak_widths', 'peak_types', 'prominence')) # find peaks/use them if supplied

    ###### find bounds
    if arguments.get('bounds_first_only', False):
        LOGGER.info('setting bounds based on first')
        bag.bound_setter()
    else:
        LOGGER.info('setting all bounds to preset')
        bounds = {'amps': False, 'centers': False, 'widths': False} # should really do this in the thunderobj
        bag.bound_setter(bounds)
    LOGGER.info('finding bounds for all data sets')
    bag.bag_iterator(getattr(bag, 'thunder_bag'), peak_finding.make_bounds, ('tightness', 'no_peaks', 'bounds', 'peak_widths',
                                              'peak_centres', 'peak_amps'), ('bounds',)) # make bounds

    ###### fit peaks
    LOGGER.info('fitting peaks for all')
    bag.bag_iterator(getattr(bag, 'thunder_bag'), peak_fitting.fit_peaks, ('x_data', 'y_data_bg_rm', 'peak_types', 'peak_centres',
                               'peak_amps', 'peak_widths', 'bounds'), ('specs', 'model', 'peak_params', 'peaks')) # fit peaks

    ##### fit params dictionary
    # store all the peak parameters in a dictionary, so the keys are e.g. sigma, center, amplitude, and the values are
    # dictionaries with keys as the run number with values as lists of values for all the peaks for that run
    # this will for now assume the same types of peak for all fits!
    LOGGER.info('making fit parameters dictionary')
    bag.make_fit_params()

    ###### fetch stats etc
    LOGGER.info('making stats dictionary')
    bag.get_fit_stats()

    ###### make directory to save everything in
    LOGGER.info(f'creating directory {file_name}')
    file_name, dirname = parsing.make_user_files(arguments, file_name)

    ###### plot map scan
    LOGGER.info('plotting map scans')
    map_scan_tools.plot_map_scan(bag, getattr(bag, 'fit_params'), dirname)

    # save individual plots for each of the failed fits
    LOGGER.info('saving failed fit plots')
    bag.save_failed_plots(dirname)

    ###### put here some code for cluster analysis and pca
    LOGGER.info('not currently doing cluster analysis or pca')
    ####################################################################################################################

    # save the bag object and it reports
    LOGGER.info('saving fit reports on stats and fitting parameters')
    utili.save_fit_report(getattr(bag, 'stats'), path=dirname, filename=f"{file_name}_report.json")
    utili.save_fit_report(getattr(bag, 'fit_params'), path=dirname, filename=f"{file_name}_peak_info.json")
    LOGGER.info('saving thunderbag object')
    utili.save_thunder(bag, path=dirname, filename=f"{file_name}.d")
