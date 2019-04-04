import numpy as np

def svn(y_data):
    """normalise using std variance normalisation"""
    mean_y_data = np.mean(y_data)
    shifted_y_data = y_data - mean_y_data
    std_dev = np.std(y_data)
    normalised_y = shifted_y_data / std_dev
    return normalised_y