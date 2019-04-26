import logging
from ast import literal_eval
from os.path import basename
from os.path import join
from os import rename
from time import strftime

from . import map_scan_tools
from . import multi_obj
from . import parsing
from . import utilities as utili

def setup_logger(log_name):
    curr_time = strftime('%d_%m_%Y_%l:%M%p')
    log_filename = f"{log_name}_{curr_time}.log"
    logging.getLogger().setLevel(logging.DEBUG)
    logger = logging.getLogger('')
    logger.handlers = []
    logging.basicConfig(filename=log_filename, level=logging.DEBUG)
    logging.info('have read in user arguments')
    return log_filename

def main():
    args = parsing.parse_user_args()

    arguments = parsing.using_user_args(args)

    try:  # a user can pass in a list of filenames or just one
        file_name = basename(literal_eval(arguments['datapath'])[0])
        log_name = file_name
    except SyntaxError:  # assume its just a string and not a list passed
        file_name = None
        log_name = arguments['datapath']
        arguments['datapath'] = f"['{arguments['datapath']}',]"  # as this is what multiobj needs

    # setup logger
    log_filename = utili.setup_logger(log_name)

    logging.info('creating multi_obj object')
    bag = multi_obj.main(arguments)  # create a Thunder object

    bag.first = next(iter((bag.thunder_bag.keys())))

    if arguments.get('clip_data', False) or arguments.get('bg_first_only', False) or arguments.get(
            'peakf_first_only', False) or arguments.get('bounds_first_only', False):
        logging.info('choosing spectrum for data')
        bag.choose_spectrum()  # choose which spectrum to base everything off of if user wants to use one spectra to
        # choose parameters
        logging.info(f'using spectra {bag.first} as spectra to set any user variables from')

    ###### clip the data if weird edges
    if arguments.get('clip_data', False):
        logging.info('clipping data')
        bag.clip_data()

    ###### cosmic ray removal goes here
    ####################################################################################################################

    ###### remove background
    if arguments.get('bg_first_only', False):
        logging.info('determining background conditions for all based on user guided for first')
        bag.bg_param_setter()
    logging.info('removing background from data for all thunder objects')
    bag.remove_backgrounds()

    ###### normalisation
    if arguments.get('normalise', False):
        logging.info('normalising data using svn normalisation')
        bag.normalise_data()

    ###### find peaks
    if arguments.get('find_peaks', False):
        print("Warning: Finding peaks automatically will overwirte and peak_info_dict supplied")
        logging.info('running user guided routine to determine peak information')
        bag.find_peaks()
    elif arguments.get('adj_params', False):
        bag.peaks_adj_params()
    logging.info('setting peak information for all thunder objects')

    ###### find bounds
    if arguments.get('find_bounds', False):
        logging.info('finding bounds via user guided routine')
        bag.bound_setter()
    else:
        logging.info("setting all bounds to either user supplied or preset: {'amplitude': False, 'center': False, 'sigma': False}")
        bounds = arguments.get('bounds', {'amplitude': False, 'center': False, 'sigma': False})  # should really do this in the thunderobj
        bag.bound_setter(bounds)

    ###### fit peaks
    logging.info('fitting peaks for all')
    bag.fit_peaks()

    ##### fit params dictionary
    # store all the peak parameters in a dictionary, so the keys are e.g. sigma, center, amplitude, and the values are
    # dictionaries with keys as the run number with values as lists of values for all the peaks for that run
    # this will for now assume the same types of peak for all fits!
    logging.info('making fit parameters dictionary')
    bag.make_fit_params()

    ###### fetch stats etc
    logging.info('making stats dictionary')
    bag.get_fit_stats()

    ###### make directory to save everything in
    file_name, dirname = parsing.make_user_files(arguments, file_name)

    ###### plot map scan
    logging.info('plotting map scans')
    map_scan_tools.plot_map_scan(bag, getattr(bag, 'fit_params'), dirname)

    # save individual plots for each of the failed fits
    logging.info('saving failed fit plots')
    bag.save_failed_plots(dirname)

    ###### put here some code for cluster analysis and pca
    logging.info('not currently doing cluster analysis or pca')
    ####################################################################################################################

    # save the bag object and it reports
    logging.info('saving fit reports on stats and fitting parameters')
    utili.save_fit_report(getattr(bag, 'stats'), path=dirname, filename=f"{file_name}_report.json")
    utili.save_fit_report(getattr(bag, 'fit_params'), path=dirname, filename=f"{file_name}_peak_info.json")
    logging.info('saving thunderbag object')
    utili.save_thunder(bag, path=dirname, filename=f"{file_name}.d")

    # move the log file in with all the rest of it
    log_filename_ = str(join(dirname, f'{file_name}.log'))
    rename(log_filename, log_filename_) # use os.rename to move the log file to the final destination
    utili.save_fit_report(arguments, path=dirname, filename=f"{file_name}_inpargs.json")


