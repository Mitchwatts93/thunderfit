import matplotlib.pyplot as plt

from . import utilities
from . import parsing
from . import multi_obj

def main():
    args = parsing.parse_user_args()

    arguments = parsing.using_user_args(args)

    bag = multi_obj.main(arguments)  # create a Thunder object

    bag.choose_spectrum() # get the user to choose which spectrum to create params from
    thund = bag.thunder_bag[bag.first]

    cents = []
    while True:
        plt.plot(thund.x_data, thund.y_data)
        for xc in cents:
            plt.axvline(x=xc)
        print(f"First looking at peak center values, look at the following plot with centers at values:"
              f"{cents}")
        plt.show()
        ans = input("Please enter a list of values e.g. [1,2,3] etc. or type y if happy with these values")
        if ans == 'y':
            break
        else:
            try:
                cents = [float(i) for i in ans]
            except:
                print("You entered an incorrect answer! Trying again...")
    plt.close()

    amps = []
    while True:
        plt.plot(thund.x_data, thund.y_data)
        for yc in amps:
            plt.axhline(y=yc)
        print(f"First looking at peak center values, look at the following plot with amps at values:"
              f"{amps}")
        plt.show()
        ans = input("Please enter a list of values e.g. [1,2,3] etc. or type y if happy with these values")
        if ans == 'y':
            break
        else:
            try:
                amps = [float(i) for i in ans]
            except:
                print("You entered an incorrect answer! Trying again...")
    plt.close()

    widths = []
    while True:
        plt.plot(thund.x_data, thund.y_data)
        for width in widths:
            plt.axvline(x=width)
        print(f"First looking at peak center values, look at the following plot with widths at values:"
              f"{widths}")
        plt.show()
        ans = input("Please enter a list of values e.g. [1,2,3] etc. or type y if happy with these values")
        if ans == 'y':
            break
        else:
            try:
                widths = [float(i) for i in ans]
            except:
                print("You entered an incorrect answer! Trying again...")
    plt.close()

    params = {"peak_centres":cents, "peak_amps":amps, "peak_widths":widths}
    utilities.save_fit_report(params, './', 'generated_params.txt')

