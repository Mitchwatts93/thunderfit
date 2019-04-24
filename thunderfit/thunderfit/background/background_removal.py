import logging

# from scipy.sparse import diags, spdiags
# from scipy.sparse.linalg import spsolve
# from scipy.optimize import least_squares
from numpy import array, ndarray

from . import scarf

"""
#### old method
def find_background(data, residual_baseline_func, baseline_asl_func):
    params = (array([0.01, 10 ** 5]))
    bounds = [array([0.001, 10 ** 5]), array([0.1, 10 ** 9])]
    baseline_values = least_squares(residual_baseline_func, params[:], args=(data.values,), bounds=bounds)
    p, lam = baseline_values['x']
    baseline_values = baseline_asl_func(data.values, lam, p, niter=10)
    return baseline_values

def residual_baseline(params, y):
    p, lam = params
    niter = 10
    baseline = baseline_als(y, lam, p, niter)
    residual = y - baseline
    return residual

def baseline_als(y, lam, p, niter=10):
    L = len(y)
    D = diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    w = ones(L)
    if niter < 1:
        raise ValueError("n iter is too small!")
    for i in range(niter):
        W = spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z
####
"""


def correct_negative_bg(y_bg_rm, bg):
    y_min = y_bg_rm.min()
    if y_min < 0:
        y_bg_rm += abs(y_min)  # then shift all the data up so no points are below zero
        bg -= abs(y_min)  # and lower the bg we have calculated by that shift too
    return y_bg_rm, bg


def background_finder(x_data, y_data, bg, scarf_params):
    if bg == 'no':  # then user doesn't want to make a background
        logging.warning(
            "Warning: no background specified, so not using a background,"
            " this may prevent algorithm from converging")
        bg = array([0 for _ in y_data])  # set the background as 0 everywhere
        data_bg_rm_y = y_data  # no background subtracted
        params = 'no'

    elif bg == 'SCARF':
        logging.debug('using SCARF method for background subtraction')
        data_bg_rm_y, bg, params = scarf.perform_scarf(x_data, y_data, scarf_params)

    elif isinstance(bg, ndarray):
        assert len(bg) == len(y_data), \
            "the background generated or passed is of incorrect length"
        logging.debug('using numpy array as supplied by user')
        data_bg_rm_y = y_data - bg  # subtract user supplied background from the data
        params = 'user_specified'

    elif bg == 'OLD':
        logging.warning('user specified old bg subtraction method which is no longer supported')
        # bg = find_background(y_data, residual_baseline, baseline_als) # find a background the old way
        # data_bg_rm_y = y_data - bg  # subtract background from the data
        # params = 'old_method'
        raise TypeError('old method is no longer supported')

    else:  # then it is the incorrect type
        logging.warning('user specified unknown bg subtraction method')
        raise TypeError('the background passed is in the incorrect format, please pass as type np array')

    return bg, data_bg_rm_y, params
