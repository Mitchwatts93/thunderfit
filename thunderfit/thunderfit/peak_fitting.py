import logging

from lmfit import models
from numpy import nan


def fit_peaks(x_data, y_data, peak_types, peak_centres, peak_amps, peak_widths, bounds):
    logging.debug(f'fitting peaks:  peak_centres:{peak_centres}, peak_amps:{peak_amps}, peak_widths:{peak_widths}, '
                  f'peak_types:{peak_types}, bounds:{bounds}')
    model_specs = build_specs(peak_types, peak_centres, peak_amps, peak_widths, bounds)
    model, peak_params = generate_model(model_specs)
    peaks = model.fit(y_data, peak_params, x=x_data)
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
            'params': {'center': peak_centres[i], 'amp': peak_amps[i] ,'sigma': peak_widths[i],
                       'gamma' :peak_widths[i]
                       },
            'bounds': {'centers': bounds['centers'][i], 'amps': bounds['amps'][i],
                       'widths': bounds['widths'][i]
                       }
        }
        for i, _ in enumerate(peak_centres)
    ]

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
    for i, basis_func in enumerate(model_specs):
        prefix = f'm{i}_'
        model = getattr(models, basis_func['type'])(prefix=prefix)
        if basis_func['type'] in ['GaussianModel', 'LorentzianModel' ,'VoigtModel']:
            # for now VoigtModel has gamma constrained to sigma
            w_min = basis_func['bounds']['widths'][0]
            w_max = basis_func['bounds']['widths'][1]
            x_min = basis_func['bounds']['centers'][0]
            x_max = basis_func['bounds']['centers'][1]
            y_min = basis_func['bounds']['amps'][0]
            y_max = basis_func['bounds']['amps'][1]

            model.set_param_hint('sigma', min=w_min, max=w_max)
            model.set_param_hint('center', min=x_min, max=x_max)
            model.set_param_hint('height', min=y_min, max=y_max)
            model.set_param_hint('amplitude', min=1e-6)

            # default guess is horrible!! do not use guess()
            default_params = {
                prefix + 'center': basis_func['params']['center'],
                prefix + 'height': basis_func['params']['amp'],
                prefix + 'sigma': basis_func['params']['sigma']
            }
        else:
            raise NotImplemented(f'model {basis_func["type"]} not implemented yet')

        model_params = model.make_params(**default_params, **basis_func.get('params', {}))

        if params is None: # first loop
            params = model_params
            composite_model = model
        else: # subsequent loops
            params.update(model_params)
            composite_model = composite_model + model

    return composite_model, params