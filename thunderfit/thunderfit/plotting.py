import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import logging


# todo fix the assertions in these
def plot_data(x, y, ax=False, line='r-', linethickness=0.5):
    logging.debug('plotting data')
    if ax:
        # assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
        " an axes object"
    else:
        fig, ax = plt.subplots()

    ax.plot(x, y, line, linewidth=linethickness, alpha=0.5)
    return ax, plt


def plot_fits(x, peaks, ax=False, linethickness=0.5):
    logging.debug('plotting fits')
    if ax:
        # assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
        " an axes object"
    else:
        fig, ax = plt.subplots()

    for i, peak in enumerate(peaks):
        ax.plot(x, peaks[peak], linewidth=linethickness)
    return ax, plt


def plot_background(x, background_data, ax=False, line='b--', linethickness=0.5):
    logging.debug('plotting bg')
    if ax:
        # assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
        " an axes object"
    else:
        fig, ax = plt.subplots()

    ax.plot(x, background_data, line, linewidth=linethickness)
    return ax, plt


def plot_fit_sum(x, peak_sum, background, ax=False, line='k-', linethickness=0.5):  # option of including background
    logging.debug('plotting fit sum')
    if ax:
        # assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
        " an axes object"
    else:
        fig, ax = plt.subplots()

    sum = peak_sum + background

    ax.plot(x, sum, line, linewidth=linethickness)
    return ax, plt


def plot_uncertainty_curve(x, eval_unc, peak_sum, ax=False, color="#ABABAB"):
    logging.debug('plotting fit uncertainty curve')
    if ax:
        # assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
        " an axes object"
    else:
        fig, ax = plt.subplots()

    ax.fill_between(x, peak_sum - eval_unc, peak_sum + eval_unc, color=color)  # plot a grey band of uncertainty

    return ax, plt
