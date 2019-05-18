Quick install
=============
Python version >=3.6 needed for this package

Mac/Linux:
----------
 
Using pip:
^^^^^^^^^^
First have pip installed: https://www.makeuseof.com/tag/install-pip-for-python/
Be sure to install python >=3.6

(optional - recommended): 
    Create a new environment in python so that packages aren't corrupted. Maintenance of this package won't be great so dependencies are set to specific releases.

    1. Choose a directory you will store your python environment in. recommended to be somewhere convenient to access
    2. `python3 -m pip install --user virtualenv`
    3. `python3 -m virtualenv thunder`
    4. When it comes time to use your environment (when installing the package or when using it):

        i. `source path_to_env/thunder/bin/activate`
        ii. (Do this step whenever you're finished using thunderfit) to deactivate just type `deactivate`

Now with you environment active (if using one) type::

    pip install thunderfit

You can check the correct script for ramananalyse (or any other script in future releases) is present by typing::

    command -v ramananalyse

To use ramananalyse::
    ramananalyse --param_file_path path_to_param_file --datapath path_to_data_file

Using anaconda:
^^^^^^^^^^^^^^^

ANACONDA NOT CURRENTLY SUPPORTED.

Windows:
--------

Currenlty untested, coming soon. Below may or may not work.

First have pip installed: https://www.makeuseof.com/tag/install-pip-for-python/
Be sure to install python >=3.6

(optional - recommended): 
    Create a new environment in python so that packages aren't corrupted. Maintenance of this package won't be great so dependencies are set to specific releases.

    1. Choose a directory you will store your python environment in. recommended to be somewhere convenient to access
    2. `py -m pip install --user virtualenv`
    3. `py -m virtualenv thunder`
    4. When it comes time to use your environment (when installing the package or when using it):

        i. `.\path_to_env\thunder\Scripts\activate`
        ii. (Do this step whenever you're finished using thunderfit) to deactivate just type `deactivate`

Now with you environment active (if using one) type::

    pip install thunderfit

scripts in windows install as .exe so check inside you env inside the thunderfit directory and see if the .exe exists


Using windows subsystem for linux (WSL):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Follow instructions for Mac/Linux


Using Thunderfit
================

To create a thunderfit object, call the Thunder class in the thundobj module with correct inputs. See code for this.

Alternatively a thunderfit object can be created by passing a thunderfit object to the Thunder class, and all attributes will be copied into a new object.

The param file
--------------

The param file is in json format and an example is below:::

    {"x_ind": 2, "y_ind": 3,"x_coord_ind":1, "y_coord_ind":0, "map": true, "background": "no",
    "clip_data": true, "clips":[3100, 1100], "method": "leastsq", "tol": 0.001, "no_peaks": 4,
    "peak_info_dict":
    {
    "type": ["LorentzianModel", "LorentzianModel", "LorentzianModel", "PowerLawModel", "LinearModel"],
    "center": [1350, 1590, 2700],
    "sigma": [15, 10, 30],
    "height": [100, 1000, 500],
    "amplitude":[null,null,null,0],
    "exponent":[null,null,null,1],
    "slope":[null,null,null,null,0],
    "intercept":[null,null,null,null,0]
    },
    "bounds":
    {
    "center": [[1330, 1370], [1570, 1610], [2680, 2730]],
    "sigma": [[8, 40], [10, 30], [30, 60]],
    "amplitude": [[0.0001,null], [0.0001,null], [0.0001,null]]
    }
    }

Possible Arguments are:

1. x_ind - the data should be in a csv format only currently. x_ind speicifies which column of the csv data is the x data
2. y_ind - the data should be in a csv format only currently. y_ind speicifies which column of the csv data is the y data
3. e_ind - (optional) the data should be in a csv format only currently. e_ind speicifies which column of the csv data is the e data. If not specified then only x and y data will be loaded

4. x_coord_ind - which column of the map has the x coordinates
5. y_coord_ind - which column of the map has the x coordinates
6. map - is this a mapscan? defaults to no
7. background - either "SCARF" or "no" to subtract either a scarf generated background or no background before fitting (note using e.g. linear models and powerlaw models is a good way to do a background simultaneously with the peaks)
8. clip_data - true or false. should the data be clipped? defaults to false
9. clips - if the data is being clipped this will be read. should be a list, e.g. `[10,20]` where the two elements are the left clip and right clip of the data. Note the order is important and if the data file has x read in backwards then the first number should be the right clip
10. method - what type of fitting method to use. uses same names as lmfit methods
11. tol - what tolerance to use. currently defaults to same as lmfit and tol is set for xtol and ftol
12. no_peaks - how many peaks to fit (will be depreciated soon)
13. peak_info_dict - this is a dictionary of information about the models to fit. the very minimum is to include type as a key. pass in the format {"key": value}. the value for all should be a list `[]` which is comma seperated. note that the element number will correspond to the model. if its not appropriate for that model type null. default is to not set parameters for models unless specified

    a. type - a key to specify models. the value. currently most of lmfits models are supported. expression model and split lorentzian currently aren't
    b. model parameters - see lmfit built in models to see which parameters can be passed

14. bounds - this has the same format as peak_info_dict except the values should be a list of list, with each sublist being two elements for a lower and upper bound on that parameter

    a.model parameters - [[low,upp],[low,upp]] replace low and upp with numerical bound values

15. datapath - the relative path to the data. Data should be in csv format. note and nan rows will be removed. - if passed into command line then that always takes precedence.
16. scarf_params - a dictionary containing parameters for the "SCARF" background method. if null then it will launch an interactive procedure for choosing the parameters which could be passed in here.

    a. rad - a number which corresponds to the radius of the rolling ball
    b. b - a number which corresponds to the shift in the background generated by rolling ball method
    c. window_length - a parameter for Savgol filter (current implementation uses scipy savgol_filter from signal)
    d. poly_order - a parameter for Savgol filter (current implementation uses scipy savgol_filter from signal)

17. normalise - bool - should the data be normalised
18. bg_first_only - bool- if finding a background with a user guided routine, should the routine only be for the first spectra
19. bounds_first_only - bool - if finding bounds interactively should do this for only first
20. peakf_first_only - bool - find peaks interactively only for first
21. find_peaks - find peaks interactively
22. adj_params - should the parameters be adjusted for each spectra e.g. by peak finding to slightly move the guess and improve covergence time?
23. find_bounds - interactively find the bounds
24. make_gif - make a gif of all the fits
25. peak_finder_type - what type of peak finding should be performed?


Scripts
=======

The below scripts will install with Thunderfit by default. They are useful for either analysing a single Raman spectra, a mapscan or generating a parameters file with user guided routine.

The ramananalyse script
-----------------------

needs user input for the param file location at a minimum

Currently this script processes user inputs and parses everything, it then creates a new directory in the current directory named analysed_{time}. This will contain all the analysis data . Then it creates a Thunder object based on input and params file. The background and the data with the background removed are then saved as variables in the object. Then peaks are fitted to the data using the peak information and the bounds information (and of course the y data with the bg removed). Then the original data, fitted peaks, background, the fit sum and the uncertainties on the fitted peaks are all plotted using matplot lib and the plot object returned. A fit report is then generated. The plots are then saved in the generated directory from earlier, as is the fit report and the Thunder object (using dill).

The map_scan script
-------------------

same to run as ramananalyse

Further details coming soon. Run like:

mapscan --param_file_path ../bag_params.txt --datapath './map.txt'
