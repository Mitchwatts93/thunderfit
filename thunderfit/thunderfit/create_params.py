import matplotlib.pyplot as plt
from os.path import basename
from ast import literal_eval

from . import utilities
from . import parsing
from . import multi_obj
from . import thundobj

def main():
    args = parsing.parse_user_args()

    arguments = parsing.using_user_args(args)

    if arguments['map']:
        try:  # a user can pass in a list of filenames or just one
            file_name = basename(literal_eval(arguments['datapath'])[0])
        except SyntaxError:  # assume its just a string and not a list passed
            file_name = None
            log_name = arguments['datapath']
            arguments['datapath'] = f"['{arguments['datapath']}',]"  # as this is what multiobj needs
        bag = multi_obj.main(arguments)  # create a Thunder object

        bag.choose_spectrum() # get the user to choose which spectrum to create params from
        thund = bag.thunder_bag[bag.first]
    else:
        thund = thundobj.main(arguments)  # create a Thunder object

    cents = []
    while True:
        plt.plot(thund.x_data, thund.y_data)
        for xc in cents:
            plt.axvline(x=xc)
        print(f"First looking at peak center values, look at the following plot with centers at values:"
              f"{cents}")
        plt.show()
        ans = input("Please enter a list of values e.g. 1,2,3 etc. or type y if happy with these values")
        if ans == 'y':
            break
        else:
            try:
                cents = ans.split(',')
                cents = [float(i) for i in cents]
            except:
                print("You entered an incorrect answer! Trying again...")
    plt.close()

    amps = []
    while True:
        plt.plot(thund.x_data, thund.y_data)
        if len(amps) == len(cents):
            plt.plot(cents, amps, 'r*')
        else:
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
                amps = ans.split(',')
                amps = [float(i) for i in amps]
            except:
                print("You entered an incorrect answer! Trying again...")
    plt.close()

    widths = []
    while True:
        plt.plot(thund.x_data, thund.y_data)
        if len(widths)%2:
            for width in widths:
                plt.axvline(x=width)
        else:
            for band in range(len(widths)//2):
                plt.axvspan(widths[2*band], widths[2*band + 1], alpha=0.1, color='red')
        print(f"First looking at peak center values, look at the following plot with widths at values:"
              f"{widths}")
        plt.show()
        ans = input("Please enter a list of values e.g. [1,2,3] etc. or type y if happy with these values")
        if ans == 'y':
            break
        else:
            try:
                widths = ans.split(',')
                widths = [float(i) for i in widths]
            except:
                print("You entered an incorrect answer! Trying again...")
    plt.close()
    widths = [abs(widths[2*i+1]-widths[2*i]) for i in range(len(widths)//2)]

    params = {"peak_centres":cents, "peak_amps":amps, "peak_widths":widths}
    utilities.save_fit_report(params, './', 'generated_params.txt')

