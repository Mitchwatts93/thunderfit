import matplotlib
import matplotlib.pyplot as plt
from numpy import argsort
from scipy.signal import find_peaks as peak_find

matplotlib.use('TkAgg')
import logging

from . import utilities as utili

def peak_finder(data, prominence, height=0, width=0):
    logging.debug(f'finding peaks based on prominence of {prominence}')
    # do a routine looping through until the right number of peaks is found

    peaks, properties = peak_find(data, prominence=prominence, height=height, width=width)  # find the peak positions in the data

    peaks = list(peaks)  # convert to a list
    amps = list(properties['peak_heights'])  # store the heights
    sorted_indices = argsort(amps)[::-1] # we will sort below in order of amplitudes

    peak_info = {'center_indices': sort_lists(sorted_indices, peaks), 'right_edges': sort_lists(sorted_indices, list(properties['right_bases'])),
                 'left_edges': sort_lists(sorted_indices, list(properties['left_bases'])), 'amplitude': sort_lists(sorted_indices, amps)}
    return peak_info

def sort_lists(sorted_indices, list_to_sort):
    return [list_to_sort[i] for i in sorted_indices]

def find_cents(prominence, y_data, find_all=False):
    logging.debug(f'finding centre values based on prominence of {prominence}')
    peak_info = peak_finder(y_data, prominence, height=0, width=0)  # find the peak centers
    if find_all:
        return peak_info
    #center_indices = peak_info['center_indices']
    return peak_info

def auto_peak_finder(prominence, x_data, y_data):
    while True:
        peak_info = find_cents(prominence, y_data, find_all=True)
        plt.plot(x_data, y_data)
        peak_coordinates = [x_data[ind] for ind in peak_info['center_indices']]
        for xc in peak_coordinates:
            plt.axvline(x=xc)
        print(f"Peak finder requires user input, please look at the following plot with prominence={prominence}")
        plt.show()
        ans = input("If you are happy with the plot, type y. If not then please type a new prominence ")
        if ans == 'y':
            break
        else:
            try:
                prominence = float(ans)
            except ValueError:
                print("You entered an incorrect answer! Trying again...")
    plt.close()
    return peak_info, prominence

def user_peak_finder(x_data, y_data):
    peak_info = {'center_values':(), 'amplitude':(), 'sigma':()}
    for key, peak_value in peak_info.items():
        while True:
            to_plot = peak_value
            if key == 'sigma':
                cents = peak_info['center_values']
                to_plot = [(cents[i] - width/2, cents[i] + width/2) for i, width in enumerate(to_plot)]
            print(f"Peak finder requires user input, please look at the following plot of values for {key}: {peak_value}")
            plot_values(key, to_plot, x_data, y_data)
            ans = input("If you are happy with the plot, type y. If not then please type a list of values seperated by commas ")
            if ans == 'y':
                break
            else:
                try:
                    peak_value_ = ans.split(',')
                    peak_value = [float(i) for i in peak_value_]
                except ValueError:
                    print("You entered an incorrect answer! Trying again...")
        plt.close()
        peak_info[key] = peak_value
    return peak_info

def plot_values(key, list_of_values, x_data, y_data):
    plt.plot(x_data, y_data)
    if key == 'center_values':
        for xc in list_of_values:
            plt.axvline(x=xc)
    elif key == 'amplitude':
        for yc in list_of_values:
            plt.axhline(y=yc)
    elif key =='sigma': # these have to be passed in as band edges
        for band in list_of_values:
            plt.axvspan(band[0], band[1], alpha=0.1, color='red')
    plt.show()

def interactive_peakfinder(x_data, y_data, type='auto', prominence=1):
    logging.debug('finding centres of peaks with user guided routine')
    if type == 'user':
        peak_info = user_peak_finder(x_data, y_data)
        prominence = None
    else:
        peak_info, prominence = auto_peak_finder(prominence, x_data, y_data)
        center_indices = peak_info['center_indices']
        peak_info['center_values'] = x_data[center_indices]
    return peak_info, prominence

def match_peak_centres(center_indices, y_data, prominence=1):
    while True:
        peak_info_ = find_cents(prominence, y_data, find_all=True)
        center_indices_ = utili.find_closest_indices(peak_info_['center_indices'],
                                                     center_indices)  # the indices might be slightly off so fix that
        if len(center_indices_) == len(center_indices):
            break
        elif len(center_indices_) > len(center_indices):
            center_indices_ = center_indices_[:len(center_indices)]  # as they're in order of prominence
            break
        else:
            prominence *= 10  # increase prominence until we get more than or equal to number.
        # do something smarter!
    center_indices = [peak_info_['center_indices'][i] for i in center_indices_]
    return center_indices

def find_peak_details(x_data, y_data, no_peaks, type='auto', prominence=1):
    logging.debug(f'finding peak details based on prominence of {prominence}, and user provided details:'
                  f'no_peaks:{no_peaks}')

    peak_info, prominence = interactive_peakfinder(x_data, y_data, type, prominence)
    no_peaks = max(len(value) for value in peak_info.values())
    center = peak_info['center_values']
    amplitude = peak_info['amplitude']
    sigma = peak_info.get('sigma', ())
    for i, sig in enumerate(sigma):
        x1 = x_data[0] + sig
        sigma[i] = utili.find_closest_indices(x_data, x1)
        import ipdb
        ipdb.set_trace()
    type = ["LorentzianModel" for i in range(no_peaks)]
    peak_info_dict = {}

    peak_info_dict['center'], peak_info_dict['amplitude'], peak_info_dict['sigma'], \
    peak_info_dict['type'] = center, amplitude, sigma, type

    return no_peaks, peak_info_dict, prominence

def make_bounds(x_data, y_data, no_peaks, peak_info_dict):
    logging.debug(f'making bounds based on: no_peaks:{no_peaks}, '
                  f'peak_info_dict: {peak_info_dict}')
    bounds = {}
    for key in peak_info_dict:
        while True:
            ans = input(f"Do you want to create bounds for {key}? \n"
                        f"Current values are: {peak_info_dict[key]} \n"
                        f"type y or n")
            if ans == 'y':
                bounds[key] = interactive_bounds(x_data, y_data, len(peak_info_dict[key]))
                break
            elif ans == 'n':
                break
            else:
                print("incorrect answer, please type y or n")
    return bounds

def interactive_bounds(x_data, y_data, no_bounds:int):
    bounds = [(0, 0) for i in range(no_bounds)]
    while True:
        plt.plot(x_data, y_data)
        for bound in bounds:
            plt.axvspan(bound[0], bound[1], alpha=0.1, color='red')
        print(f"current bounds are: {bounds}")
        plt.show()
        ans = input(f"Please enter a list of values for the bounds, comma delimited. You must enter {2*no_bounds}"
                    f" numbers. type y if you are happy with these bounds")
        if ans == 'y':
            break
        else:
            try:
                bounds = ans.split(',')
                bounds = [float(i) for i in bounds]
                if len(bounds) != 2 * no_bounds:
                    raise ValueError("incorrect no of bounds passed")
                else:
                    bounds = [(bounds[2*i], bounds[2*i+1]) for i in range(no_bounds)]
            except Exception as e:
                print("You entered an incorrect answer! Trying again...")
    plt.close()
    return bounds