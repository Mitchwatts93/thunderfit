Quick install
=============
Python version >=3.6 needed for this package

Mac/Linux:
----------
 
Using pip:
~~~~~~~~~~
First have pip installed: https://www.makeuseof.com/tag/install-pip-for-python/
Be sure to install python >=3.6

(optional - recommended): 
	Create a new environment in python so that packages aren't corrupted. Maintenance of this package won't be great so dependencies are set to specific releases.
	1. Choose a directory you will store your python environment in. recommended to be somewhere convenient to access
	2. :: python3 -m pip install --user virtualenv
	3. :: python3 -m virtualenv thunder
	4. When it comes time to use your environment (when installing the package or when using it):
		i. :: source path_to_env/thunder/bin/activate
		ii. (Do this step whenever you're finished using thunderfit) to deactivate just type ::deactivate

Now with you environment active (if using one) type::
	pip install thunderfit

You can check the correct script for ramananalyse (or any other script in future releases) is present by typing::
	command -v ramananalyse

To use ramananalyse::
	ramananalyse --param_file_path path_to_param_file --datapath path_to_data_file

Using anaconda:
~~~~~~~~~~~~~~~
First have conda installed: https://www.anaconda.com/distribution/
Be sure to install python >=3.6

(optional - recommended): 
	Create a new environment in python so that packages aren't corrupted. Maintenance of this package won't be great so dependencies are set to specific releases.
	1. Choose a directory you will store your python environment in. recommended to be somewhere convenient to access
	2. :: conda create --n thunder python=3.6
	3. When it comes time to use your environment (when installing the package or when using it):
		i. :: conda activate thunder
		ii. (Do this step whenever you're finished using thunderfit) to deactivate just type ::conda deactivate 

Now with you environment active (if using one) type::
	conda skeleton pypi thunderfit
	conda build thunderfit

You can check the correct script for ramananalyse (or any other script in future releases) is present by typing::
	command -v ramananalyse

To use ramananalyse::
	ramananalyse --param_file_path path_to_param_file --datapath path_to_data_file

Windows:
--------

Currenlty untested, coming soon. Below may or may not work.

First have pip installed: https://www.makeuseof.com/tag/install-pip-for-python/
Be sure to install python >=3.6

(optional - recommended): 
	Create a new environment in python so that packages aren't corrupted. Maintenance of this package won't be great so dependencies are set to specific releases.
	1. Choose a directory you will store your python environment in. recommended to be somewhere convenient to access
	2. :: py -m pip install --user virtualenv
	3. :: py -m virtualenv thunder
	4. When it comes time to use your environment (when installing the package or when using it):
		i. :: .\path_to_env\thunder\Scripts\activate
		ii. (Do this step whenever you're finished using thunderfit) to deactivate just type ::deactivate. 

Now with you environment active (if using one) type::
	pip install thunderfit

scripts in windows install as .exe so check inside you env inside the thunderfit directory and see if the .exe exists


Using windows subsystem for linux (WSL):
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Follow instructions for Mac/Linux


Creating a thunderfit object
============================

To create a thunderfit object, call the Thunder class in the thundobj module. This class has many inputs currently, but is going to be cleaned up in the near future. Currently the primary usage of this code is via the ramananalyse command which is created. This takes in two main arguments: --param_file_path and --datapath (it also takes in many more object but this will also be cleaned up soon). The datapath argument is self explanatory, and is just a relative path to the data to be analysed. The --param_file_path argu,emt specifies a relative path to a parameters file, which contains all the Thunder class variables needed. NOTE: not all inputs are passified yet so please be careful with input and follow guidelines below.

Alternatively a thunderfit object can be created by passing a thunderfit object to the Thunder class, and all attributes will be copied into a new object.

The param file
~~~~~~~~~~~~~~

The param file is in json format and an example is below:::

	{"x_ind" : 0, "y_ind" : 1, "e_ind" : null, "datapath": "data.txt", "no_peaks":3, "background": "SCARF", "scarf_params":{"rad":70, "b":90, "window_length":51, "poly_order":3}, "peak_types": [], "peak_centres": [], "peak_widths":[], "peak_amps":[], "tightness":"med", "bounds" : {"centers":null,"widths":null,"amps":null}}

Arguments are:
1. x_ind - the data should be in a csv format only currently. x_ind speicifies which column of the csv data is the x data
2. y_ind - the data should be in a csv format only currently. y_ind speicifies which column of the csv data is the y data
3. e_ind - (optional) the data should be in a csv format only currently. e_ind speicifies which column of the csv data is the e data. If not specified then only x and y data will be loaded
4. datapath - the relative path to the data. Data should be in csv format. note and nan rows will be removed.
5. no_peaks - the number of peaks to be fitted. Can either be null or can be an integer. Note that increasing the number of peaks without specifying bounds etc may result in a bad fit
6. yfit - null (I think no longer used - will be cleaned up in future version)
7. background - either a numpy array (of same length as the data! - data with nan rows removed!) containing numerical values or a string: "SCARF" for a rolling ball-GS smoothed background to be fitted (see parameters below or will be interactive) reference in relevant functions or "OLD" which uses a numerical method to fit the background, which usually struggles somewhat - need to put reference in for this method soon. Can also be "no" if no backgorund subtraction is wanted.
8. scarf_params - a dictionary containing parameters for the "SCARF" background method. if null then it will launch an interactive procedure for choosing the parameters which could be passed in here.
	a. rad - a number which corresponds to the radius of the rolling ball
	b. b - a number which corresponds to the shift in the background generated by rolling ball method
	c. window_length - a parameter for Savgol filter (current implementation uses scipy savgol_filter from signal)
	d. poly_order - a parameter for Savgol filter (current implementation uses scipy savgol_filter from signal)
9. peak_types - a list of peak types, these will be models used by lmfit, see documentation for lmfit for supported models, currently using "LorentzianModel", "GaussianModel" or "VoigtModel" only others not implemented yet. VoigtModel will have gamma set as sigma for now. if less are specified than no_peaks then these will be ignored. if more are specified then will clip the list to [:no_peaks]
10. peak_centres - a list of peak centre values, these must be the x values, not the indices of the peaks from y_data. if less are specified than no_peaks then these will be ignored. if more are specified then will clip the list to [:no_peaks]
11. peak_widths - a list of widths. if less are specified than no_peaks then these will be ignored. if more are specified then will clip the list to [:no_peaks]
12. peak_amps - a list of amplitudes. if less are specified than no_peaks then these will be ignored. if more are specified then will clip the list to [:no_peaks]
13. tightness - either "low" "med" or "high" which controls how relaxed the bounds will be around the peaks found or specified if bounds are not given
14. bounds - a dictionary of the following:
	a. centers - peak center bounds. a list of tuples, which each tuple is length=2. Each tuple contains the low and high bounds for each peak - if an empty list, list not the same length as no_peaks or null then bounds will automatically be generated.
	b. widths - peak width bounds. a list of tuples, which each tuple is length=2. Each tuple contains the low and high bounds for each peak - if an empty list, list not the same length as no_peaks or null then bounds will automatically be generated.
	c. amps - peak amplitude bounds. a list of tuples, which each tuple is length=2. Each tuple contains the low and high bounds for each peak - if an empty list, list not the same length as no_peaks or null then bounds will automatically be generated.

The very minimum which can be supplied is the datapath, x_ind and y_ind so that the data can be loaded. If this isn't specified in param file (or on command line by --datapath --x_ind and --y_ind) then it will fail.

The ramananalyse script
~~~~~~~~~~~~~~~~~~~~~~~

Currently this script processes user inputs and parses everything, it then creates a new directory in the current directory named analysed_{time}. This will contain all the analysis data (and in a future version also a log file - currently logs are output directly to user). Then it creates a Thunder object based on input and params file. The background and the data with the background removed are then saved as variables in the object. Currenly it doesn't normalise but in the future there will be an option to perform normalisation on the background subtracted data, and then on the generated background and original data in order to make a nice plot at the end (currently only svn normalisation is implemented and bg and original data use the mean and stddev from the background subtracted data).  Then it determines if peak information has been passed by the user, and finds the peak information automatically or just uses the information if correct. Then bounds are either used or generated for these peaks. Then peaks are fitted to the data using the peak information and the bounds information (and of course the y data with the bg removed). Then the original data, fitted peaks, background, the fit sum and the uncertainties on the fitted peaks (if available - will be improved in future release) are all plotted using matplot lib and the plot object returned. A fit report is then generated. The plots are then saved in the generated directory from earlier, as is the fit report and the Thunder object (using dill).

