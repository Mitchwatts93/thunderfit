import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
from scipy.signal import find_peaks as peak_find
import numpy as np

from . import utilities as utili

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

    matching_indices = utili.find_closest_indices(center_indices, center_list)

    if peak_info_key=='widths':
        peak_properties = ([peak_info['right_edges'][i] for i in matching_indices],
                           [peak_info['left_edges'][i] for i in matching_indices])
    else:
        peak_properties = [peak_info[peak_info_key][i] for i in matching_indices]
    return peak_properties

#THIS NEEDS A BIT OF TIDYING UP?
def peaks_unspecified(x_data, y_data, peak_no, peak_centres, peak_amps, peak_widths, peak_types):
    #todo: pass in prominence x2 values
    if len(peak_centres) == 0 or len(peak_centres) < peak_no:
        if peak_no and len(peak_centres) < peak_no and len(peak_centres):
            logging.warning("you specified less peak centers than peak_numbers."
                 " Currently only finding all peaks based on tightness criteria or using all supplied is possible")
        prominence = 1.6
        if not peak_no: # then they don't know so we can find everything in one go and save some time
            peak_info = find_cents(prominence, y_data, find_all=True)
            center_indices = peak_info['center_indices']
            peak_amps = peak_info['amps']
            peak_left_edges, peak_right_edges = peak_info['left_edges'], peak_info['left_edges']
            peak_widths = x_data[peak_right_edges] - x_data[peak_left_edges] # the xvalues can be
                                                                                            # indexed from the data
            peak_centres = x_data[center_indices]
            peak_no = len(center_indices)
        else: # just find the centers
            center_indices = find_cents(prominence, y_data)
            center_indices = center_indices[:peak_no] # take the first n as user has specified how many peaks
            peak_centres = x_data[center_indices]
    elif len(peak_centres) > peak_no:
        logging.warning("specified more peak centers than no_peaks. cutting the peaks supplied as [:no_peaks]")
        peak_centres = peak_centres[:peak_no]

    if len(peak_amps) == 0 or len(peak_amps) < peak_no:
        if peak_no and len(peak_amps) < peak_no and len(peak_amps):
            logging.warning("you specified less peak amps than peak_numbers."
                " Currently only finding all peaks based on tightness criteria or using all supplied is possible")
        peak_amps = find_peak_properties(1, peak_centres, y_data, 'amps')
    elif len(peak_amps) > peak_no:
        logging.warning("specified more peak amps than no_peaks. cutting the peaks supplied as [:no_peaks]")
        peak_amps = peak_amps[:peak_no]

    if len(peak_widths) == 0 or len(peak_widths) < peak_no:
        if peak_no and len(peak_widths) < peak_no and len(peak_widths):
            logging.warning("you specified less peak widths than peak_numbers."
                " Currently only finding all peaks based on tightness criteria or using all supplied is possible")
        peak_left_edges, peak_right_edges = find_peak_properties(1, peak_centres, y_data, 'widths') # get edge indices
        peak_widths = x_data[peak_right_edges] - x_data[peak_left_edges]
    elif len(peak_widths) > peak_no:
        logging.warning("specified more peak widths than no_peaks. cutting the peaks supplied as [:no_peaks]")
        peak_widths = peak_widths[:peak_no]

    if len(peak_types) == 0 or len(peak_types) < peak_no:
        if peak_no and len(peak_types) < peak_no and len(peak_types):
            logging.warning("you specified less peak types than peak_numbers."
                " Currently only finding all peaks based on tightness criteria or using all supplied is possible")
        peak_types = ['LorentzianModel' for _ in peak_centres]  # we assume all the types are Lorentzian for now
    elif len(peak_types) > peak_no:
        logging.warning("specified more peak types than no_peaks. cutting the peaks supplied as [:no_peaks]")
        peak_types = peak_widths[:peak_no]

    return peak_no, peak_centres, peak_amps, peak_widths, peak_types