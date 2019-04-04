from scipy.signal import find_peaks as peak_find
from itertools import product


def peak_finder(data, prominence):
    # do a routine looping through until the right number of peaks is found

    peaks, properties = peak_find(data, prominence=prominence)  # find the peak positions in the data

    peaks = list(peaks)  # convert to a list
    amps = list(properties['prominences'])  # store the heights

    peak_info = {'center_indices': peaks, 'right_edges': list(properties['right_bases']),
                 'left_edges': list(properties['left_bases']), 'amps': amps}
    return peak_info

def find_cents(prominence, y_data, find_all=False):
    peak_info = peak_finder(y_data, prominence)  # find the peak centers
    if find_all:
        return peak_info
    center_indices = peak_info['center_indices']
    return center_indices

def find_peak_properties(prominence, center_list, y_data, peak_info_key):
    peak_info = peak_finder(y_data, prominence)
    center_indices = peak_info['center_indices']

    matching_indices = find_closest_indices(center_indices, center_list)

    if peak_info_key=='widths':
        peak_properties = ([peak_info['left_edges'][i] for i in matching_indices],
                           [peak_info['right_edges'][i] for i in matching_indices])
    else:
        peak_properties = [peak_info[peak_info_key][i] for i in matching_indices]
    return peak_properties

def find_closest_indices(list1, list2):
    list_of_matching_indices = [min(range(len(list1)), key=lambda i: abs(list1[i] - cent))
                                for cent in list2]
    return list_of_matching_indices
