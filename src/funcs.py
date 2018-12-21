import numpy as np

def gaussian(x_data, cent, width, amp):
    return amp * np.exp(-(x_data - cent) ** 2 / (2 * width **2))