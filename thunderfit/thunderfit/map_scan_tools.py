import matplotlib
import matplotlib.pyplot as plt
from numpy import unique, round, array, nanmin, nanmax, nan
from scipy.sparse import coo_matrix

matplotlib.use('TkAgg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
import logging

from . import utilities as utili


# funcs for plotting
def shift_map_matr(coordinates_array):
    logging.debug('shifting coordinates array')
    coordinates_array[:, 0] = coordinates_array[:, 0] - min(coordinates_array[:, 0])
    coordinates_array[:, 1] = coordinates_array[:, 1] - min(coordinates_array[:, 1])
    return coordinates_array


def map_scan_plot(coordinates, values):
    logging.debug('plotting mapscans')
    no_fits = len(list(values.values())[0])
    figs = []
    axs = []
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
        data_ = coo_matrix((Z, (xx, yy))).toarray() # out=nanmat should mean that any missing values are nan
        data_[data_ == 0] = nan
        f = plt.figure()
        ax = plt.gca()
        magma_cmap = matplotlib.cm.get_cmap('magma')
        magma_cmap.set_bad(color='green')
        im = plt.imshow(data_, cmap=magma_cmap, extent=[min(X), max(X), max(Y), min(Y)], vmin=nanmin(data_),
                        vmax=nanmax(data_))
        plt.xlabel('x coordinates')
        plt.ylabel('y coordinates')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        figs.append(f)
        axs.append(ax)
    return figs, axs


def plot_map_scan(bag, fit_params, dirname):
    logging.debug('runnning user input routine to generate/save user chosen variables in maps')
    coordinates_array = array(list(getattr(bag, 'coordinates').values()))  # convert coordinates for each point into an
    # array
    coordinates_array = shift_map_matr(
        coordinates_array)  # shift so that each coordinate starts at 0, not normalised
    for i, key in enumerate(getattr(bag, 'coordinates')):
        getattr(bag, 'coordinates')[key] = coordinates_array.tolist()[i]  # reassign in the correct format
    while True:
        plt.close()
        plt.clf()
        ans = input("making map scans, please input which property you would like to scan. options are:"
                    f"\n {[p_ for p_ in fit_params.keys()]}, or type 'all' to plot all, or type e to exit")
        if ans == 'e':
            break
        elif ans == 'all':
            for p in fit_params.keys():
                plot, ax = map_scan_plot(getattr(bag, 'coordinates'), fit_params.get(p))
                for i, pt in enumerate(plot):
                    try:
                        cents = next(iter(fit_params.get('center').values()))
                        pt.suptitle(f'{p}_{i}_heatmap. peak {i} is centered at: {cents[i]}')
                    except KeyError:
                        print("tried to add label for peak centers onto graph, but couldn't fetch the right variable")
                        pt.suptitle(f'{p}_{i}_heatmap')
                for i, pt in enumerate(plot):
                    utili.save_plot(pt, path=dirname, figname=f"{p}_{i}.svg")
            break
        else:
            try:
                p = ans
                plot, ax = map_scan_plot(getattr(bag, 'coordinates'), fit_params.get(p))
                for i, pt in enumerate(plot):
                    try:
                        cents = next(iter(fit_params.get('center').values()))
                        pt.suptitle(f'{p}_{i}_heatmap. peak {i} is centered at: {cents[i]}')
                    except KeyError:
                        print("tried to add label for peak centers onto graph, but couldn't fetch the right variable")
                        pt.suptitle(f'{p}_{i}_heatmap')
            except KeyError:
                p = ''
                ax = None
                plot = None
                print('wrong answer entered, trying again!')
            try:
                for i, pt in enumerate(plot):
                    utili.save_plot(pt, path=dirname, figname=f"{p}_{i}.svg")
            except AttributeError:
                print("Tried to save plot but there is no plot yet! Something wen't wrong in making the plot")
#
