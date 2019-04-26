import logging

from lmfit import models
from numpy import nan

from . import utilities as utils

def get_tols(method, tol):
    if method in ['leastsq', 'least_squares']:
        tols = {'ftol': tol, 'xtol': tol}  # warning: need to experiment with these/allow control maybe?
    else:
        tols = {'ftol': tol}
    return tols

def prep_algo(method, tols):
    if method == 'basinhopping' or method == 'ampgo':
        print('Warning: This is a very slow but thorough algorithm')
    if method == 'differential_evolution':
        all_bounds = True
    if method in ['leastsq', 'least_squares', 'nelder', 'cobyla']:
        tols  = get_tols(method, tols)
    return tols

def fit_peaks(x_data, y_data, peak_info_dict, bounds, method='leastsq', tol=0.0000001):
    impl_methods = ['leastsq', 'least_squares', 'nelder', 'lbfgsb', 'powell', 'cg', 'cobyla', 'bfgsb',
                    'differential_evolution', 'basinhopping', 'ampgo']
    if not method in impl_methods:
        raise ValueError(f"The method supplied is not supported. Available methods: {impl_methods}")
    logging.debug(f'fitting peaks:  {peak_info_dict}')

    tols = prep_algo(method, tol)
    model_specs = build_specs(peak_info_dict, bounds)
    model, peak_params = generate_model(model_specs)

    if method in ['leastsq', 'least_squares', 'nelder', 'cobyla']:
        peaks = model.fit(y_data, peak_params, x=x_data, method=method, fit_kws=tols)
        if not peaks.success:
            print('peaks failed to fit, raising tolerance by one order magnitude and trying again')
            tols = {key:value*10 for key, value in tols.items()} # try raising the tolerance if it fails by one order
            peaks = model.fit(y_data, peak_params, x=x_data, method=method, fit_kws=tols)
    else:
        peaks = model.fit(y_data, peak_params, x=x_data, method=method)

    peak_params = peaks.best_values
    if not peaks.success:
        print('peaks failed to fit')
        peak_params = {key: nan for key in peak_params}

    return model_specs, model, peak_params, peaks

def safe_list_get (l, idx, default):
    """fetch items safely from a list, if it isn't long enough return a default value"""
    try:
        return l[idx]
    except IndexError:
        return default

##### Call this with correct arguments, should always call with all the arguments!
def build_specs(peak_info_dict: dict={}, bounds: dict={}):
    """Build a specs list, which has the peak specification details for each element corresponding to a peak to fit.
    Note that each element is a dictionary containing either None, or a value for that peak i.e. the parameter guess,
    or the bounds for that value. Only some will be valid for each peak type, which will be handled elsewhere"""
    logging.debug('building specs')
    type = peak_info_dict.get('type', ())
    if 'PolynomialModel' in type:
        raise NotImplementedError("No support for polynomial models yet")
    center = peak_info_dict.get('center', ())
    amplitude = peak_info_dict.get('amplitude', ())
    sigma = peak_info_dict.get('sigma', ())
    sigma_r = peak_info_dict.get('sigma_r', ())
    gamma = peak_info_dict.get('gamma', ())
    fraction = peak_info_dict.get('fraction', ())
    beta = peak_info_dict.get('beta', ())
    exponent = peak_info_dict.get('exponent', ())
    q = peak_info_dict.get('q', ())
    c = peak_info_dict.get('c', ())
    intercept = peak_info_dict.get('intercept', ())
    slope = peak_info_dict.get('slope', ())
    a = peak_info_dict.get('a', ())
    b = peak_info_dict.get('b', ())
    degree = peak_info_dict.get('degree', ())
    #poly_consts = peak_info_dict.get('poly_consts', ())
    form = peak_info_dict.get('form', ())
    center1 = peak_info_dict.get('center1', ())
    center2 = peak_info_dict.get('center2', ())
    sigma1 = peak_info_dict.get('sigma1', ())
    sigma2 = peak_info_dict.get('sigma2', ())
    decay = peak_info_dict.get('decay', ())

    specs = [
        {
            'type': type[i],
            'params': {'center': utils.safe_list_get(center, i, None),
                       'amplitude': utils.safe_list_get(amplitude, i, None),
                       'sigma': utils.safe_list_get(sigma, i, None),
                       'sigma_r': utils.safe_list_get(sigma_r, i, None),
                       'gamma': utils.safe_list_get(gamma, i, None),
                       'fraction': utils.safe_list_get(fraction, i, None),
                       'beta': utils.safe_list_get(beta, i, None),
                       'exponent': utils.safe_list_get(exponent, i, None),
                       'q': utils.safe_list_get(q, i, None),
                       'c': utils.safe_list_get(c, i, None),
                       'intercept': utils.safe_list_get(intercept, i, None),
                       'slope': utils.safe_list_get(slope, i, None),
                       'a': utils.safe_list_get(a, i, None),
                       'b': utils.safe_list_get(b, i, None),
                       'degree': utils.safe_list_get(degree, i, None),

                       'form': utils.safe_list_get(form, i, None),
                       'center1': utils.safe_list_get(center1, i, None),
                       'center2': utils.safe_list_get(center2, i, None),
                       'sigma1': utils.safe_list_get(sigma1, i, None),
                       'sigma2': utils.safe_list_get(sigma2, i, None),
                       'decay': utils.safe_list_get(decay, i, None)
                       },
            'bounds': {'center': utils.safe_list_get(bounds.get('center', []), i,(None, None)),
                       'amplitude': utils.safe_list_get(bounds.get('amplitude', []), i, (None, None)),
                       'sigma': utils.safe_list_get(bounds.get('sigma', []), i, (None, None)),
                       'sigma_r': utils.safe_list_get(bounds.get('sigma_r', []), i, (None, None)),
                       'gamma': utils.safe_list_get(bounds.get('gamma', []), i, (None, None)),
                       'fraction': utils.safe_list_get(bounds.get('fraction', []), i, (None, None)),
                       'beta': utils.safe_list_get(bounds.get('beta', []), i, (None, None)),
                       'exponent': utils.safe_list_get(bounds.get('exponent', []), i, (None, None)),
                       'q': utils.safe_list_get(bounds.get('q', []), i, (None, None)),
                       'c': utils.safe_list_get(bounds.get('c', []), i, (None, None)),
                       'intercept': utils.safe_list_get(bounds.get('intercept', []), i, (None, None)),
                       'slope': utils.safe_list_get(bounds.get('slope', []), i, (None, None)),
                       'a': utils.safe_list_get(bounds.get('a', []), i, (None, None)),
                       'b': utils.safe_list_get(bounds.get('b', []), i, (None, None)),
                       'degree': utils.safe_list_get(bounds.get('degree', []), i, (None, None)),

                       'form': utils.safe_list_get(bounds.get('form', []), i, (None, None)),
                       'center1': utils.safe_list_get(bounds.get('center1', []), i, (None, None)),
                       'center2': utils.safe_list_get(bounds.get('center2', []), i, (None, None)),
                       'sigma1': utils.safe_list_get(bounds.get('sigma1', []), i, (None, None)),
                       'sigma2': utils.safe_list_get(bounds.get('sigma2', []), i, (None, None)),
                       'decay': utils.safe_list_get(bounds.get('decay', []), i, (None, None))
                       }
        }
        for i in range(len(type))
    ]

    #'poly_consts': utils.safe_list_get(poly_consts, i, None), and
    # 'poly_consts': utils.safe_list_get(bounds.get('poly_consts', []), i, [None, None]),
    # not yet implemented!

    return specs

def generate_model(model_specs):
    """
    https://chrisostrouchov.com/post/peak_fit_xrd_python/
    :param model_specs:
    :return:
    """
    logging.debug('generating model specs')
    composite_model = None
    params = None
    for i, spec in enumerate(model_specs):
        prefix = f'm{i}_'
        model = getattr(models, spec['type'])(prefix=prefix) # generate the lmfit model based on the type specified
        model = decide_model_actions(spec, model) # call another function to decide what to do
        model_params = model.make_params() # make the params object
        if params is None: # first loop
            params = model_params
            composite_model = model
        else: # subsequent loops
            params.update(model_params)
            composite_model = composite_model + model
    return composite_model, params

def decide_model_actions(spec, model):
    for param_key, param_value in spec['params'].items():
        if param_value: # then set this value
            model.set_param_hint(param_key, value=param_value)
    for bound_key, bound_value in spec['bounds'].items():
        if bound_value[0]: # then set lower bound
            model.set_param_hint(bound_key, min=bound_value[0])
        if bound_value[1]: # then set upper bound
            model.set_param_hint(bound_key, max=bound_value[1])
    return model

