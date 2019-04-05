import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
import logging
LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

def rcf(A, rad):
    """
    Popular rolling circle filter (RCF) routine as described by James et al.: https://doi.org/10.1366%2F12-06766
    :param A: np matrix of data
    :param r: radius of circle to filter with
    :return: L, a numpy array of length len(A) which is the locus of the RCF
    """
    if isinstance(A, pd.core.frame.DataFrame) or isinstance(A, pd.core.series.Series):
        A = A.values # convert to numpy array
    elif not isinstance(A, np.ndarray):
        raise TypeError("Incorrect format passed for A")

    assert isinstance(rad, int) and abs(rad) == rad, "rad needs to be a positive integer"
    assert len(A.shape) == 1, "A needs to be a single column"

    RC = semi_circle_array(rad)
    n = len(A)

    if len(RC) >= n:
        print('something has gone wrong')
        raise RuntimeError("The data has less points than the selected radius!")

    A_sub, RC_sub = generate_sub_matrices(A, rad, RC, n)

    D = generate_diff_matrix(A_sub, RC_sub)

    return D

def semi_circle_array(rad):
    """Generate an array of length 2r+1 with values corresponding to points on a semi circle"""
    RC = np.array([np.sqrt(rad**2 - (i - rad)**2) for i in range(2*rad + 1)])
    return RC

def generate_sub_matrices(A, rad, RC, n):
    """
    Generate the sub matrices used in RCF routine
    :param A: Data matrix, length n
    :param rad: radius of RCF filter
    :param RC: np array of length 2r+1 with
    :param n:
    :return: A_sub, RC_sub, lists of same dimensions, sampled according to paper: https://doi.org/10.1366%2F12-06766
    """
    ni = n
    A_sub = []
    RC_sub = []

    for i in range(len(A)): # vectorise this when you can be bothered
        if i < (rad + 1):
            a_sub = A[:(i + rad + 1)]
            rc_sub = RC[(rad - i):(2*rad + 1)]
            A_sub.append(a_sub)
            RC_sub.append(rc_sub)

        elif (rad + 1) <= i <= (n - (rad + 1)):
            a_sub = A[(i - rad):(i + rad + 1)]
            rc_sub = RC
            A_sub.append(a_sub)
            RC_sub.append(rc_sub)

        elif i > (n - (rad + 1)):
            a_sub = A[(i - rad):(n + 1)]
            rc_sub = RC[:(rad + (n - i))]
            A_sub.append(a_sub)
            RC_sub.append(rc_sub)

        else:
            print('something has gone wrong')
            import ipdb
            ipdb.set_trace()
            raise RuntimeError("Something has gone wrong and I don't know why")

    return A_sub, RC_sub

def generate_diff_matrix(A_sub, RC_sub):
    """
    Generate a difference matrix, which will be len(A_sub)=len(RC_sub). Each element will be the minimum of matrix
    formed by the difference of AC_sub_i and RC_sub_i at each i (element).
    :param A_sub: List of np arrays
    :param RC_sub: List of np arrays
    :return: D, Difference matrix of shape (len(A_sub),), elements as described above
    """
    D = []
    for i in range(len(A_sub)):
        diff = A_sub[i] - RC_sub[i]
        D.append(np.min(diff))

    D = np.array(D)
    return D

def smooth(L, window_length, polyorder):
    """
    Savgol filter applied to data
    :param L: the data to be smoothed
    :param window_length: positive odd integer
    :param polyorder: must be less than window length
    :return:
    """
    return savgol_filter(L, window_length, polyorder, mode='mirror')

def perform_scarf(x_data, y_data, scarf_params=False):
    bg = np.array([0 for _ in y_data], dtype=np.float64)
    rad = 70
    b = 90
    window_length, poly_order = 51, 3
    L_sg = 0
    data_bg_rm_y = y_data #.copy????

    if not scarf_params:
        while True:
            while True:
                D = rcf(data_bg_rm_y, rad)
                fig, ax = plt.subplots()
                ax.plot(x_data, D)
                ax.plot(x_data, data_bg_rm_y)
                print(f"SCARF background removal requires user input. Please look at the following bg with rad={rad}")
                plt.show(block=True)
                ans = input("If you are happy with the plot, type y. if not then please type a new rad")
                if ans == 'y':
                    break
                else:
                    try:
                        rad = int(ans)
                    except ValueError:
                        print("You entered an incorrect answer! Trying again...")

            L = D + b
            while True:  # now estimate a baseline to add to D to get L
                fig, ax = plt.subplots()
                ax.plot(x_data, L)
                ax.plot(x_data, data_bg_rm_y)
                print(f"Please look at the following bg with a shift={b}")
                plt.show(block=True)
                ans = input("If you are happy with the plot, type y. if not then please type a new background value. \n"
                            "Please note that the background should NOT intercept the data. Ideally it would pass through"
                            "the mean of the noise for the correct bg already fit")
                if ans == 'y':
                    L = D + b
                    break
                else:
                    try:
                        b = int(ans)
                        L = D + b
                    except ValueError:
                        print("You entered an incorrect answer! Trying again...")

            # then apply SG filter to L
            while True:
                try:
                    L_sg = smooth(L, window_length, poly_order)
                    fig, ax = plt.subplots()
                    ax.plot(x_data, L_sg)
                    ax.plot(x_data, data_bg_rm_y)
                    print(f"Please look at the following bg with Sg filter parameters (window length, polynomial order): "
                          f"{window_length}, {poly_order}")
                    plt.show(block=True)
                except ValueError as e:
                    print(
                        "Incorrect values for window_length and poly_order have been entered. Poly order must be less than window length and window length must be odd")
                ans = input("please enter y if you are happy with these values, or enter two integers with a space "
                            "for window_length and poly_order")
                if ans == 'y':
                    L = L_sg
                    break
                else:
                    try:
                        ans = ans.split(' ')
                        if len(ans) != 2:
                            raise ValueError("The tuple was more than two elements long")
                        window_length = int(ans[0])
                        poly_order = int(ans[1])
                    except ValueError:
                        print("You entered an incorrect answer! Trying again...")

            # final question before exiting
            fig, ax = plt.subplots()
            ax.plot(x_data, L)
            ax.plot(x_data, data_bg_rm_y)
            print(f"Please look at the following bg with selected parameters")
            plt.show(block=True)
            ans = input("Are you happy with this bg? If yes, type y, else type n. n will restart the fitting. \n"
                        "typing repeat will add an additional bg subtraction to this one")
            if ans == 'y':
                bg += L
                data_bg_rm_y -= L
                break
            elif ans == 'n':
                pass
            elif ans == 'repeat':
                bg += L
                print("apply two bg removal steps, this will mean the background just specified will be removed "
                      "from the data")
                data_bg_rm_y -= L  # remove the bg found here from the original data and go again
            else:
                print("You entered an incorrect answer! Trying whole fitting routine again...")

    else:
        rad, b, window_length, poly_order = \
            scarf_params['rad'], scarf_params['b'], scarf_params['window_length'], scarf_params['poly_order']
        D = rcf(data_bg_rm_y, rad)
        L = D + b
        L = smooth(L, window_length, poly_order)
        bg += L
        data_bg_rm_y -= L

    return data_bg_rm_y, bg