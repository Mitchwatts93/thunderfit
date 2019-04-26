import logging

from numpy import mean, std


def svn(y_data, mean_y_data=False, std_dev=False):
    logging.debug('normalising using svn normalisation')
    """normalise using std variance normalisation"""
    if not mean_y_data and not std_dev:
        mean_y_data = mean(y_data)
        std_dev = std(y_data)

    shifted_y_data = y_data - mean_y_data
    normalised_y = shifted_y_data / std_dev
    return normalised_y, (mean_y_data, std_dev)

