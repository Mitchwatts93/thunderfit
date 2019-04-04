import numpy as np

def svn(y_data, mean_y_data=False, std_dev=False):
    """normalise using std variance normalisation"""
    if not mean_y_data and not std_dev:
        mean_y_data = np.mean(y_data)
        std_dev = np.std(y_data)

    shifted_y_data = y_data - mean_y_data
    normalised_y = shifted_y_data / std_dev
    return normalised_y, (mean_y_data, std_dev)
