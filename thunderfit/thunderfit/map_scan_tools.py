import matplotlib
import matplotlib.pyplot as plt
from numpy import unique, round, array, nanmin, nanmax, nan, nanpercentile
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

def map_scan_matrices(coordinates, values):
    logging.debug('generating map matrices')
    no_maps = len(list(values.values())[0])
    data = {}
    X_ = {}
    Y_ = {}
    for i in range(no_maps):
        X = []
        Y = []
        Z = []
        for key in values.keys():
            x, y = coordinates[key]
            z = utili.safe_list_get(values[key], i, nan)
            if type(z) == str:
                return [], []
            X.append(x)
            Y.append(y)
            Z.append(z)
        x_step = unique(X)[1] - unique(X)[0]
        xx = (round(X / x_step)).astype(int)
        y_step = unique(Y)[1] - unique(Y)[0]
        yy = round((Y / y_step)).astype(int)
        data_ = coo_matrix((Z, (yy, xx))).toarray()  # out=nanmat should mean that any missing values are nan
        data_[data_ == 0] = nan
        data[i] = data_
        X_[i] = X
        Y_[i] = Y
    return data, X_, Y_

def map_scan_plot(data_mat, X_coords, Y_coords):
    logging.debug('plotting mapscans')
    no_fits = len(data_mat.keys())
    figs = []
    axs = []
    for i in range(no_fits):
        data = data_mat[i]
        X = X_coords[i]
        Y = Y_coords[i]
        f = plt.figure()
        ax = plt.gca()
        magma_cmap = matplotlib.cm.get_cmap('magma')
        magma_cmap.set_bad(color='green')
        im = plt.imshow(data, cmap=magma_cmap, vmin=nanpercentile(data,1), extent=[min(X), max(X), max(Y), min(Y)],
                        vmax=nanpercentile(data,99)) # we plot the 99th percentile and 1st percentil as the max and min
        # colours
        plt.xlabel('x coordinates')
        plt.ylabel('y coordinates')
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)
        figs.append(f)
        axs.append(ax)
    return figs, axs

def plot_map_scan(bag, fit_params, map_matrices, X_coords, Y_coords, dirname):
    logging.debug('runnning user input routine to generate/save user chosen variables in maps')
    while True:
        plt.close()
        plt.clf()
        ans = input("making map scans, please input which property you would like to scan. options are:"
                    f"\n {[p_ for p_ in fit_params.keys()]}, or type 'all' to plot all, or type e to exit")
        if ans == 'e':
            break
        elif ans == 'all':
            for p in fit_params.keys():
                data_mat = map_matrices[p]
                plot, ax = map_scan_plot(data_mat, X_coords[p], Y_coords[p])
                for i, pt in enumerate(plot):
                    try:
                        cents = next(iter(fit_params.get('center').values()))
                        pt.suptitle(f"{p}_{i}_heatmap. peak {i} is centered at: {utili.safe_list_get(cents, i, 'na')}")
                    except (KeyError, AttributeError):
                        print("tried to add label for peak centers onto graph, but couldn't fetch the right variable")
                        pt.suptitle(f'{p}_{i}_heatmap')
                for i, pt in enumerate(plot):
                    utili.save_plot(pt, path=dirname, figname=f"{p}_{i}.svg")
            break
        else:
            try:
                p = ans
                data_mat = map_matrices[p]
                plot, ax = map_scan_plot(data_mat, X_coords[p], Y_coords[p])
                for i, pt in enumerate(plot):
                    try:
                        cents = next(iter(fit_params.get('center').values()))
                        pt.suptitle(f'{p}_{i}_heatmap. peak {i} is centered at: {cents[i]}')
                    except (KeyError, AttributeError):
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
    return map_matrices