from numpy import unique, round, array
from scipy.sparse import coo_matrix
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
from mpl_toolkits.axes_grid1 import make_axes_locatable

from . import utilities as utili
from .background import background_removal as bg_remove
from . import peak_finding



##### funcs for analysing each spectra
def shift_map_matr(coordinates_array):
    coordinates_array[:, 0] = coordinates_array[:, 0] - min(coordinates_array[:, 0])
    coordinates_array[:, 1] = coordinates_array[:, 1] - min(coordinates_array[:, 1])
    return coordinates_array

def map_scan_plot(coordinates, values):
    no_fits = len(list(values.values())[0])
    figs = []
    data = []
    for i in range(no_fits):
        X = []
        Y = []
        Z = []
        for key in values.keys():
            x, y = coordinates[key]
            z = values[key][i]
            X.append(x)
            Y.append(y)
            Z.append(z)
        x_step = unique(X)[1] - unique(X)[0]
        xx = (round(X / x_step)).astype(int)
        y_step = unique(Y)[1] - unique(Y)[0]
        yy = round((Y / y_step)).astype(int)
        data_ = coo_matrix((Z, (xx, yy))).toarray()
        data.append(data_)
        f = plt.figure()
        ax = plt.gca()
        im = plt.imshow(data_, cmap='magma', extent=[min(X), max(X), max(Y), min(Y)], vmin=data_.min(),
                        vmax=data_.max())
        plt.xlabel('x coordinates')
        plt.ylabel('y coordinates')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        # f.tight_layout()
        figs.append(f)
    return figs

def choose_spectrum(bag, arguments):
    first = list(bag.thunder_bag.keys())[0]
    if  arguments.get('clip_data', False) or arguments.get('bg_first_only', False) or arguments.get('peakf_first_only', False) or arguments.get('bounds_first_only', False):
        # then we have to choose which spectrum we want
        first = list(bag.thunder_bag.keys())[0]
        while True:
            try:
                first_thunder = bag.thunder_bag[first]
            except KeyError:
                print('incorrect key, please enter a lower index value')
                first = list(bag.thunder_bag.keys())[0]
                first_thunder = bag.thunder_bag[first]
            fig, ax = plt.subplots()
            ax.plot(first_thunder.x_data, first_thunder.y_data)
            print(f"Need a decision on which plot is representitive of data, the following is for index {first}")
            plt.show(block=True)
            ans = input("If you are happy with using this data file, type y, otherwise enter a new index")
            if ans == 'y':
                break
            else:
                try:
                    first = str(ans)
                except ValueError:
                    print("You entered an incorrect answer! Trying again...")

    return first

def clip_data(bag, first):
    first_thunder = bag.thunder_bag[first]
    clip_left, clip_right = utili.clip_data(first_thunder.x_data, first_thunder.y_data)
    for thund in bag.thunder_bag.values():
        setattr(thund, 'x_data', thund.x_data[clip_left:clip_right])
        setattr(thund, 'y_data', thund.y_data[clip_left:clip_right])
    return bag

def bg_param_setter(bag, first):
    # add step to find bg parameters for first one and use for the rest.
    first_thunder = bag.thunder_bag[first]
    _, _, params = bg_remove.background_finder(first_thunder.x_data, first_thunder.y_data, first_thunder.background,
                                               first_thunder.scarf_params)
    [param.pop('b', None) for param in params]  # we want to find b each time so don't set it for all others
    for thund in bag.thunder_bag.values():
        setattr(thund, 'scarf_params', params)  # set all the values to this
    return bag

def peak_info_setter(bag, first):
    # add step to find bg parameters for first one and use for the rest.
    first_thunder = bag.thunder_bag[first]
    no_peaks, peak_centres, peak_amps, peak_widths, peak_types, prominence = \
        peak_finding.find_peak_details(first_thunder.x_data, first_thunder.y_data_bg_rm, first_thunder.no_peaks,
                                       first_thunder.peak_centres, first_thunder.peak_amps, first_thunder.peak_widths,
                                       first_thunder.peak_types)
    for thund in bag.thunder_bag.values():  # set these first values for all of them
        setattr(thund, 'no_peaks', no_peaks)  # set values
        setattr(thund, 'peak_centres', peak_centres)  # set values
        setattr(thund, 'peak_types', peak_types)  # set values

    return bag

def bound_setter(bag, first, bounds=None):
    if not bounds:
        first_thunder = bag.thunder_bag[first]
        bounds = peak_finding.make_bounds(first_thunder.tightness, first_thunder.no_peaks, first_thunder.bounds,
                                      first_thunder.peak_widths, first_thunder.peak_centres, first_thunder.peak_amps)
    for thund in bag.thunder_bag.values():  # set these first values for all of them
        setattr(thund, 'bounds', bounds)  # set values
    return bag
#####

##### funcs for making info dicts
def make_fit_params(first, bag):
    fit_params = {}
    first_thunder = bag.thunder_bag[first]
    params = list(first_thunder.peak_params.keys())[
             : len(first_thunder.peak_params) // first_thunder.no_peaks]  # what are the peak params for the first peak
    params = [param.split('_')[1] for param in params]  # keep only the type of param
    for param in params:
        fit_params[param] = {}  # e.g. 'center'
        for key in bag.thunder_bag.keys():
            fit_details = bag.thunder_bag[key].peak_params
            fit_details = [fit_details[key_] for key_ in fit_details.keys() if param in key_]
            fit_params[param][key] = fit_details
    return fit_params

def get_fit_stats(bag):
    stats = {'chisq': {}, 'reduced_chi_sq': {}, 'free_params': {}}
    for key in bag.thunder_bag.keys():
        chisq = bag.thunder_bag[key].peaks.chisqr
        reduced_chi_sq = bag.thunder_bag[key].peaks.redchi
        free_params = round(chisq / reduced_chi_sq)
        stats['chisq'][key] = chisq
        stats['reduced_chi_sq'][key] = reduced_chi_sq
        stats['free_params'][key] = free_params
    return stats
#####

##### funcs for plotting
def plot_map_scan(bag, fit_params, dirname):
    coordinates_array = array(list(bag.coordinates.values()))  # convert coordinates for each point into an array
    coordinates_array = shift_map_matr(
        coordinates_array)  # shift so that each coordinate starts at 0, not normalised
    for i, key in enumerate(bag.coordinates):
        bag.coordinates[key] = coordinates_array.tolist()[i]  # reassign in the correct format
    while True:
        ans = input("making map scans, please input which property you would like to scan. options are:"
                    f"\n {[p_ for p_ in fit_params.keys()]}. or type e to exit")
        if ans == 'e':
            break
        else:
            try:
                p = ans
                plot = map_scan_plot(bag.coordinates, fit_params[p])
                for i, pt in enumerate(plot):
                    try:
                        cents = next(iter(fit_params['center'].values()))
                        pt.suptitle(f'{p}_{i}_heatmap. peak {i} is centered at: {cents[i]}')
                    except KeyError:
                        print("tried to add label for peak centers onto graph, but couldn't fetch the right variable")
                        pt.suptitle(f'{p}_{i}_heatmap')
            except KeyError:
                p = ''
                plot = None
                print('wrong answer entered, trying again!')
            try:
                for i, pt in enumerate(plot):
                    utili.save_plot(pt, path=dirname, figname=f"{p}_{i}.svg")
                    pt.close()
            except AttributeError:
                print("Tried to save plot but there is no plot yet! Somethuing wen't wrong in making the plot")

def save_failed_plots(bag, dirname):
    for key, thund in bag.thunder_bag.items():
        if not thund.peaks.success:
            thund.plot_all()
            utili.save_plot(thund.plot, path=dirname, figname=f"failed_plot_{key}.svg")
            thund.plot.close() # close so memory is conserved.
#####