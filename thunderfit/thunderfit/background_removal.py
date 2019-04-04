import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import least_squares
import numpy as np

#### old method
def find_background(data, residual_baseline_func, baseline_asl_func):
    params = (np.array([0.01, 10 ** 5]))
    bounds = [np.array([0.001, 10 ** 5]), np.array([0.1, 10 ** 9])]
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
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    w = np.ones(L)
    if niter < 1:
        raise ValueError("n iter is too small!")
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z
####

def correct_negative_bg(y_bg_rm, bg):
    y_min = y_bg_rm.min()
    if y_min < 0:
        y_bg_rm += abs(y_min)  # then shift all the data up so no points are below zero
        bg -= abs(y_min)  # and lower the bg we have calculated by that shift too
    return y_bg_rm, bg