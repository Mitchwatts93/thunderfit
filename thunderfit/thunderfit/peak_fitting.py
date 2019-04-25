import logging

from lmfit import models
from numpy import nan

def get_tols(method, tol):
    if method in ['least_sq', 'least_squares']:
        tols = {'ftol': tol, 'xtol': tol}  # warning: need to experiment with these/allow control maybe?
    else:
        tols = {'ftol': tol}
    return tols

def fit_peaks(x_data, y_data, peak_types, peak_centres, peak_amps, peak_widths, bounds, method='leastsq', tol=0.0000001,
              amp_bounds=False):
    impl_methods = ['leastsq', 'least_squares', 'nelder', 'lbfgsb', 'powell', 'cg', 'cobyla', 'bfgsb',
                    'differential_evolution', 'basinhopping', 'ampgo']
    if not method in impl_methods:
        raise ValueError(f"The method supplied is not supported. Available methods: {impl_methods}")

    if method == 'basinhopping' or method == 'ampgo':
        print('Warning: This is a very slow but thorough algorithm')

    logging.debug(f'fitting peaks:  peak_centres:{peak_centres}, peak_amps:{peak_amps}, peak_widths:{peak_widths}, '
                  f'peak_types:{peak_types}, bounds:{bounds}')

    model_specs = build_specs(peak_types, peak_centres, peak_amps, peak_widths, bounds)
    if method == 'differential_evolution':
        amp_bounds = True
    model, peak_params = generate_model(model_specs, amp_bounds)


    if method in ['least_sq', 'least_squares', 'nelder', 'cobyla']:
        tols  = get_tols(method, tol)
        peaks = model.fit(y_data, peak_params, x=x_data, method=method, fit_kws=tols)
    else:
        peaks = model.fit(y_data, peak_params, x=x_data, method=method)

    peak_params = peaks.best_values
    if not peaks.success:
        print('peaks failed to fit')
        peak_params = {key: nan for key in peak_params}
    return model_specs, model, peak_params, peaks

def build_specs(peak_types, peak_centres, peak_amps, peak_widths, bounds):
    logging.debug('building specs')
    specs = [
        {
            'type': peak_types[i],
            'params': {'center': peak_centres[i], 'amp': peak_amps[i] ,'sigma': peak_widths[i]
                       },
            'bounds': {'centers': bounds['centers'][i], 'amps': bounds['amps'][i],
                       'widths': bounds['widths'][i]
                       }
        }
        for i, _ in enumerate(peak_centres)
    ]

    return specs

def generate_model(model_specs, amp_bounds=False):
    """
    https://chrisostrouchov.com/post/peak_fit_xrd_python/
    :param model_specs:
    :return:
    """
    logging.debug('generating model specs')
    composite_model = None
    params = None
    for i, basis_func in enumerate(model_specs):
        prefix = f'm{i}_'
        model = getattr(models, basis_func['type'])(prefix=prefix)
        if basis_func['type'] in ['GaussianModel', 'LorentzianModel' ,'VoigtModel']:
            # for now VoigtModel has gamma constrained to sigma
            w_min = basis_func['bounds']['widths'][0]
            w_max = basis_func['bounds']['widths'][1]
            x_min = basis_func['bounds']['centers'][0]
            x_max = basis_func['bounds']['centers'][1]

            model.set_param_hint('sigma', value=basis_func['params']['sigma'], min=w_min, max=w_max)
            model.set_param_hint('center', value=basis_func['params']['center'], min=x_min, max=x_max)
            if amp_bounds: # for differential_evolution algorithm this is needed
                y_min = basis_func['bounds']['amps'][0]
                y_max = basis_func['bounds']['amps'][1]
                model.set_param_hint('amplitude', value=basis_func['params']['amp'], min=y_min, max=y_max)
            else: # otherwise we don't put any bounds on
                model.set_param_hint('height', value=basis_func['params']['amp'])
        else:
            raise NotImplemented(f'model {basis_func["type"]} not implemented yet')

        model_params = model.make_params()

        if params is None: # first loop
            params = model_params
            composite_model = model
        else: # subsequent loops
            params.update(model_params)
            composite_model = composite_model + model

    return composite_model, params