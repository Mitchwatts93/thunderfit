import logging

from . import thundobj
from . import utilities as utili
from .background import background_removal as bg_remove
from . import peak_finding
from . import peak_fitting
from . import parsing


def main():
    LOGGER = logging.getLogger(__name__)
    LOGGER.setLevel(logging.INFO)

    args = parsing.parse_user_args()

    arguments, LOGGER = parsing.using_user_args(args, LOGGER)

    thunder = thundobj.main(arguments) # create a Thunder object

    thunder.background, thunder.y_data_bg_rm, _ = bg_remove.background_finder(thunder.x_data, thunder.y_data,
                                                                           thunder.background, thunder.scarf_params)
                                                                           # determine the background
    if args.normalise:
        thunder.y_data_bg_rm, thunder.background, thunder.y_data_norm = \
                                                 utili.normalise_all(thunder.y_data_bg_rm, thunder.background, thunder.y_data)

    thunder.no_peaks, thunder.peak_centres, thunder.peak_amps, thunder.peak_widths, thunder.peak_types, _ = \
                   peak_finding.find_peak_details(thunder.x_data, thunder.y_data_bg_rm, thunder.no_peaks,
                                                  thunder.peak_centres, thunder.peak_amps, thunder.peak_widths,
                                                  thunder.peak_types) # find peaks/use them if supplied

    thunder.bounds = peak_finding.make_bounds(thunder.tightness, thunder.no_peaks, thunder.bounds, thunder.peak_widths,
                                              thunder.peak_centres, thunder.peak_amps) # make bounds

    thunder.specs, thunder.model, thunder.peak_params, thunder.peaks = \
        peak_fitting.fit_peaks(thunder.x_data, thunder.y_data_bg_rm, thunder.peak_types, thunder.peak_centres,
                               thunder.peak_amps, thunder.peak_widths, thunder.bounds) # fit peaks

    thunder.chi_sq = thunder.peaks.chisqr # set the stats from the fits
    reduced_chi_sq = thunder.peaks.redchi
    thunder.free_params = round(thunder.chi_sq / reduced_chi_sq)

    thunder.plot_all() # plot the data in full and save as an object
    thunder.gen_fit_report() # generate a fit report

    # save a plot of the figure and the thunder object
    file_name, dirname = parsing.make_user_files(arguments, file_name=None)

    utili.save_plot(thunder.plot, path=dirname, figname=f"{file_name}.svg")
    utili.save_thunder(thunder, path=dirname, filename=f"{file_name}.d")
    utili.save_fit_report(thunder.fit_report, path=dirname, filename=f"{file_name}_report.json")