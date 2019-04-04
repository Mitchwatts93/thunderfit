from scipy.signal import find_peaks as peak_find
import numpy as np

def peak_finder(data, prominence, height=0, width=0):
    # do a routine looping through until the right number of peaks is found

    peaks, properties = peak_find(data, prominence=prominence, height=height, width=width)  # find the peak positions in the data

    peaks = list(peaks)  # convert to a list
    amps = list(properties['peak_heights'])  # store the heights
    sorted_indices = np.argsort(amps)[::-1] # we will sort below in order of amplitudes

    peak_info = {'center_indices': sort_lists(sorted_indices, peaks), 'right_edges': sort_lists(sorted_indices, list(properties['right_bases'])),
                 'left_edges': sort_lists(sorted_indices, list(properties['left_bases'])), 'amps': sort_lists(sorted_indices, amps)}
    return peak_info

def sort_lists(sorted_indices, list_to_sort):
    return [list_to_sort[i] for i in sorted_indices]

def find_cents(prominence, y_data, find_all=False):
    peak_info = peak_finder(y_data, prominence, height=0, width=0)  # find the peak centers
    if find_all:
        return peak_info
    center_indices = peak_info['center_indices']
    return center_indices

def find_peak_properties(prominence, center_list, y_data, peak_info_key):
    peak_info = peak_finder(y_data, prominence, height=0, width=0)
    center_indices = peak_info['center_indices']

    matching_indices = find_closest_indices(center_indices, center_list)

    if peak_info_key=='widths':
        peak_properties = ([peak_info['right_edges'][i] for i in matching_indices],
                           [peak_info['left_edges'][i] for i in matching_indices])
    else:
        peak_properties = [peak_info[peak_info_key][i] for i in matching_indices]
    return peak_properties

def find_closest_indices(list1, list2):
    list_of_matching_indices = [min(range(len(list1)), key=lambda i: abs(list1[i] - cent))
                                for cent in list2]
    return list_of_matching_indices
