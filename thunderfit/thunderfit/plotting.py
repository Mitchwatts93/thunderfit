import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
import matplotlib.pyplot as plt

#todo fix the assertions in these
def plot_data(x, y, ax=False, line='r-', linethickness=0.5):
        if ax:
            #assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
                                                                 " an axes object"
        else:
            fig, ax = plt.subplots()

        ax.plot(x, y, line, linewidth=linethickness, alpha=0.5)
        return ax

def plot_fits(x, peaks, ax=False, linethickness=0.5):
        if ax:
            #assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
                                                                " an axes object"
        else:
            fig, ax = plt.subplots()

        for i, peak in enumerate(peaks):
            ax.plot(x, peaks[peak], linewidth=linethickness)
        return ax

def plot_background(x, background_data, ax=False, line='b--', linethickness=0.5):
        if ax:
            #assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
                                                                 " an axes object"
        else:
            fig, ax = plt.subplots()

        ax.plot(x, background_data, line, linewidth=linethickness)
        return ax

def plot_fit_sum(x, peak_sum, background, ax=False, line='k-', linethickness=0.5): # option of including background
        if ax:
            #assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
                                                                 " an axes object"
        else:
            fig, ax = plt.subplots()

        sum = peak_sum + background

        ax.plot(x, sum, line, linewidth=linethickness)
        return ax

def plot_uncertainty_curve(x, eval_unc, peak_sum, ax=False, color="#ABABAB"):
        if ax:
            #assert isinstance(ax, axes._subplots.AxesSubplot), "the figure passed isn't the correct format, please pass" \
                                                                 " an axes object"
        else:
            fig, ax = plt.subplots()

        ax.fill_between(x, peak_sum - eval_unc, peak_sum + eval_unc, color=color) #plot a grey band of uncertainty

        return ax