import matplotlib
import matplotlib.pyplot as plt
from numpy import unique, round, array, nanmin, nanmax, nan, nanpercentile, nanmean
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

def generate_map_matrix(coordinates, peak_label_vals):
    X = []
    Y = []
    Z = []

    for key in peak_label_vals.keys():
        x, y = coordinates[key]
        z = peak_label_vals[key]
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

    return X, Y, data_

def map_scan_matrices_from_dicts(coordinates, values):
    logging.debug('generating map matrices')
    peak_labels = list(values.values())[0].keys()
    data = {}
    X_ = {}
    Y_ = {}
    for peak_label in peak_labels:
        peak_label_vals = {key: values[key].get(peak_label, nan) for key in values.keys()}
        X, Y, data_ = generate_map_matrix(coordinates, peak_label_vals)
        data_[data_ == 0] = nan
        data[peak_label] = data_
        X_[peak_label] = X
        Y_[peak_label] = Y
    return data, X_, Y_

def map_plotter(data, X, Y):
    f = plt.figure()
    ax = plt.gca()
    magma_cmap = matplotlib.cm.get_cmap('magma')
    magma_cmap.set_bad(color='green')
    im = plt.imshow(data, cmap=magma_cmap, vmin=nanpercentile(data, 3), extent=[min(X), max(X), max(Y), min(Y)],
                    vmax=nanpercentile(data, 97))  # we plot the 99th percentile and 1st percentil as the max and min
    # colours
    plt.xlabel('x coordinates')
    plt.ylabel('y coordinates')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

    return f, ax

def map_scan_plot_dicts(data_mat, X_coords, Y_coords):
    logging.debug('plotting mapscans')
    peak_labels = data_mat.keys()
    figs = {}
    axs = {}
    for peak_label in peak_labels:
        data = data_mat[peak_label]
        X = X_coords[peak_label]
        Y = Y_coords[peak_label]
        f, ax = map_plotter(data, X, Y)
        figs[peak_label] = f
        axs[peak_label] = ax
    return figs, axs

def get_cent_mean(peak_label, fit_params):
    cents = [cent.get(peak_label, 0) for cent in fit_params.get('center', nan).values()]
    cent_mean = nanmean(cents)
    return cent_mean

def save_mapscan(peak_label, fit_params, plot, dirname, p):
    cent_mean = get_cent_mean(peak_label, fit_params)
    plot[peak_label].suptitle(f"{p}_{peak_label}_heatmap. peak {peak_label} is centered at:  {cent_mean}")
    utili.save_plot(plot[peak_label], path=dirname, figname=f"{p}_{peak_label}.svg")

def mapscans_for_parameter(map_matrices, X_coords, Y_coords, p, fit_params, dirname):
    data_mat = map_matrices[p]
    plot, ax = map_scan_plot_dicts(data_mat, X_coords[p], Y_coords[p])
    for peak_label in plot.keys():
        save_mapscan(peak_label, fit_params, plot, dirname, p)
    plt.close()


def plot_map_scan(fit_params, map_matrices, X_coords, Y_coords, dirname):
    logging.debug('runnning user input routine to generate/save user chosen variables in maps')
    while True:
        ans = input("making map scans, please input which property you would like to scan. options are:"
                    f"\n {[p_ for p_ in fit_params.keys()]}, or type 'all' to plot all, or type e to exit")
        if ans == 'e':
            break
        elif ans == 'all':
            for p in fit_params.keys():
                mapscans_for_parameter(map_matrices, X_coords, Y_coords, p, fit_params, dirname)
            break #break the while loop
        else:
            try:
                p = ans
                mapscans_for_parameter(map_matrices, X_coords, Y_coords, p, fit_params, dirname)
            except KeyError:
                print('wrong answer entered, trying again!')
#
    return map_matrices