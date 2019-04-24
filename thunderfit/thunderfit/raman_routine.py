import logging
from os import rename
from os.path import join
from time import strftime

from . import parsing
from . import peak_finding
from . import peak_fitting
from . import thundobj
from . import utilities as utili
from .background import background_removal as bg_remove


def main():

    args = parsing.parse_user_args()

    arguments = parsing.using_user_args(args)

    # save a plot of the figure and the thunder object

    curr_time = strftime('%d_%m_%Y_%l:%M%p')
    file_name = arguments['datapath']
    log_filename = f'{file_name}_{curr_time}.log'
    logging.getLogger().setLevel(logging.DEBUG)
    logger = logging.getLogger('')
    logger.handlers = []
    logging.FileHandler(log_filename, mode="w", encoding=None, delay=False)
    logging.basicConfig(filename=log_filename, level=logging.DEBUG)

    logging.info('creating thunder obj')
    thunder = thundobj.main(arguments) # create a Thunder object

    logging.info('setting and subtracting bg')
    thunder.background, thunder.y_data_bg_rm, _ = bg_remove.background_finder(thunder.x_data, thunder.y_data,
                                                                           thunder.background, thunder.scarf_params)
                                                                           # determine the background
    if args.normalise:
        logging.info('normalising data')
        thunder.y_data_bg_rm, thunder.background, thunder.y_data_norm = \
                                                 utili.normalise_all(thunder.y_data_bg_rm, thunder.background, thunder.y_data)

    logging.info('setting peak info')
    thunder.no_peaks, thunder.peak_centres, thunder.peak_amps, thunder.peak_widths, thunder.peak_types, _ = \
                   peak_finding.find_peak_details(thunder.x_data, thunder.y_data_bg_rm, thunder.no_peaks,
                                                  thunder.peak_centres, thunder.peak_amps, thunder.peak_widths,
                                                  thunder.peak_types) # find peaks/use them if supplied

    logging.info('setting bounds')
    thunder.bounds = peak_finding.make_bounds(thunder.tightness, thunder.no_peaks, thunder.bounds, thunder.peak_widths,
                                              thunder.peak_centres, thunder.peak_amps) # make bounds

    logging.info('fitting peaks')
    thunder.specs, thunder.model, thunder.peak_params, thunder.peaks = \
        peak_fitting.fit_peaks(thunder.x_data, thunder.y_data_bg_rm, thunder.peak_types, thunder.peak_centres,
                               thunder.peak_amps, thunder.peak_widths, thunder.bounds) # fit peaks

    logging.info('setting stats etc')
    thunder.chi_sq = thunder.peaks.chisqr # set the stats from the fits
    reduced_chi_sq = thunder.peaks.redchi
    thunder.free_params = round(thunder.chi_sq / reduced_chi_sq)

    logging.info('plotting and genertaing fit reports')
    thunder.plot_all() # plot the data in full and save as an object
    thunder.gen_fit_report() # generate a fit report

    # create a directory to save everything in
    file_name, dirname = parsing.make_user_files(arguments, file_name=None)

    logging.info('saving plots, reports and thund obj')
    utili.save_plot(thunder.plot, path=dirname, figname=f"{file_name}.svg")
    utili.save_thunder(thunder, path=dirname, filename=f"{file_name}.d")
    utili.save_fit_report(thunder.fit_report, path=dirname, filename=f"{file_name}_report.json")

    # move the log file in with all the rest of it
    log_filename_ = str(join(dirname, f'{file_name}.log'))
    rename(log_filename, log_filename_)  # use os.rename to move the log file to the final destination