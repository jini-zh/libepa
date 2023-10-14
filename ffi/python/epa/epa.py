import enum

from epa._epa_cffi import lib, ffi
import epa._epa_functions

lib.epa_init()

_python_error_layer = 0x6e6f6874 # "Pyth"

from epa._epa_vars import *
proton_mass                       = lib.epa_proton_mass()
proton_magnetic_moment            = lib.epa_proton_magnetic_moment()
proton_dipole_form_factor_lambda2 = lib.epa_proton_dipole_form_factor_lambda2()

pp_elastic_slope                  = lib.epa_pp_elastic_slope

class Error(Exception):
    def __init__(self, message):
        super(Exception, self).__init__(message)

def _clear_error():
    global _error_handle
    lib.epa_clear_error()
    _error_handle = None

def _set_error(error):
    global _error_handle
    _error_handle = ffi.new_handle(error)
    lib.epa_set_error(_error_handle, _python_error_layer)

def _fail(error = None):
    if not error:
        error = lib.epa_get_error()
    layer = lib.epa_get_error_layer()
    if layer == _python_error_layer:
        exception = ffi.from_handle(error)
    elif layer == lib.epa_cpp_error_layer:
        exception = Error(ffi.string(lib.epa_cpp_error_message(error)))
    else:
        raise Exception('Unexpected EPA error layer: ', layer)
    _clear_error()
    raise exception

def _check_error():
    error = lib.epa_get_error()
    if error:
        _fail(error)

class IntegrationMethod(enum.IntEnum):
    gauss15 = lib.GSL_INTEG_GAUSS15,
    gauss21 = lib.GSL_INTEG_GAUSS21,
    gauss31 = lib.GSL_INTEG_GAUSS31,
    gauss41 = lib.GSL_INTEG_GAUSS41,
    gauss51 = lib.GSL_INTEG_GAUSS51,
    gauss61 = lib.GSL_INTEG_GAUSS61

def get_default_integration_method():
    return IntegrationMethod(lib.epa_get_default_integration_method())

def set_default_integration_method(method):
    lib.epa_set_default_integration_method(method.value)

def QAGWorkspace(limit = 1000):
    return ffi.gc(
            lib.epa_make_qag_integration_workspace(limit),
            lib.epa_destroy_qag_integration_workspace
    )

def CQuadWorkspace(limit = 100):
    return ffi.gc(
            lib.epa_make_cquad_integration_workspace(limit),
            lib.epa_destroy_cquad_integration_workspace
    )

class Function:
    functions = epa._epa_functions.functions
    callbacks = {}
    types     = {}

    def __init__(self, function, type = None, handles = None):
        if function == ffi.NULL:
            _fail()

        if not type:
            if isinstance(function, ffi.CData):
                type = ffi.typeof(function)
            else:
                raise Exception('epa.Function: no type specified for ' + function)
        elif isinstance(type, str):
            type = Function._typeof(type)

        if type.kind == 'pointer' and Function._is_epa_function_type(type):
            ftype = type.item.fields[0][1].type
        elif type.kind == 'function':
            ftype = type
        else:
            raise Exception('Invalid function type: ' + type.cname)

        if not isinstance(function, ffi.CData):
            handle = ffi.new_handle(function)
            callback = ffi.cast('void (*)()', Function._callback(ftype))
            function = lib.epa_make_function(callback, handle, ffi.NULL)
            if not handles:
                handles = []
            handles.append(handle)
            if type.kind == 'pointer':
                handles.append(function)
                function = ffi.cast(type, function)

        function = ffi.gc(function, Function._destroy_function)

        lift_result = Function._is_epa_function_type(ftype.result)
        wrap = lift_result
        if not wrap:
            for arg in ftype.args:
                if Function._is_epa_function_type(arg):
                    wrap = True
                    break

        # no function arguments must be lowered, nor the result must be lifted
        def call_simple(*args):
            return function.function(*args, function.data)

        # some of the arguments must be lowered, or the result must be lifted
        def call_complicated(*args):
            handles = []
            args = (
                    _lower(type, arg, handles) \
                    for type, arg in zip(ftype.args, args)
            )
            result = call_simple(*args)
            if lift_result:
                return Function(result, ftype.result, handles)
            return result

        if wrap:
            call = call_complicated
        else:
            call = call_simple

        self.epa_function = function
        self.function     = call
        self.handles      = handles

    def __call__(self, *args):
        result = self.function(*args)
        _check_error()
        return result

    def _typeof(type):
        if isinstance(type, str):
            try:
                return Function.types[type]
            except KeyError:
                t = ffi.typeof(type)
                if t.kind == 'function':
                    t = ffi.typeof(type[:-1] + ', void*)')
                elif not Function._is_epa_function_type(t):
                    raise Exception('Invalid function type: ' + type)
                Function.types[type] = t
                return t
        return type

    def _callback(type):
        def callback(*args):
            try:
                return ffi.from_handle(args[-1])(*args[:-1])
            except Exception as e:
                _set_error(e)
                return ffi.cast(type.result, 0)
        try:
            return Function.callbacks[type.cname]
        except KeyError:
            cb = ffi.callback(type, callback)
            Function.callbacks[type.cname] = cb
            return cb

    def _is_epa_function_type(type):
        return      type.kind == 'pointer' \
                and type.item.kind == 'struct' \
                and type.item.cname[7:] in Function.functions

    def _destroy_function(epa_function):
        lib.epa_destroy_function(ffi.cast('epa_function*', epa_function))

def _lower(type, arg, handles):
    if isinstance(arg, Function):
        handles.append(arg)
        return arg.epa_function
    if callable(arg):
        f = Function(arg, type)
        handles.append(f)
        return f.epa_function
    return arg

def _lower_integrator(integrator, handles):
    if integrator:
        return _lower('epa_integrator*', integrator, handles)
    return ffi.NULL

def _lower_integrator_generator(generator, handles):
    if generator:
        def generate(level):
            return _lower('epa_integrator*', generator(level), handles)
        return Function(generate, 'epa_integrator_generator*')
    return ffi.NULL

def _lower_spectra(type, spectrum1, spectrum2, handles):
    spectrum1_f = _lower(type, spectrum1, handles)
    if spectrum2 and spectrum2 is not spectrum1:
        spectrum2_f = _lower(type, spectrum2, handles)
    else:
        spectrum2_f = spectrum1_f
    return spectrum1_f, spectrum2_f

def _integrator_generator(generator, level):
    return generator(
            relative_error = get_default_relative_error() \
                           * get_default_error_step() ** level
    )

def qag_integrator(
        absolute_error = None,
        relative_error = None,
        method         = None,
        workspace      = None
):
    if absolute_error is None:
        absolute_error = get_default_absolute_error()
    if relative_error is None:
        relative_error = get_default_relative_error()
    if method is None:
        method = get_default_integration_method()
    if not workspace:
        workspace = ffi.NULL
    return Function(
            lib.epa_qag_integrator(
                absolute_error, relative_error, method, workspace
            )
    )

def qag_integrator_generator(level):
    return _integrator_generator(qag_integrator, level)

def cquad_integrator(
        absolute_error = None,
        relative_error = None,
        workspace      = None
):
    if absolute_error is None:
        absolute_error = get_default_absolute_error()
    if relative_error is None:
        relative_error = get_default_relative_error()
    if not workspace:
        workspace = ffi.NULL
    return Function(
            lib.epa_cquad_integrator(absolute_error, relative_error, workspace)
    )

def cquad_integrator_generator(level):
    return _integrator_generator(cquad_integrator, level)

default_integrator = qag_integrator_generator

def form_factor_monopole(lambda2):
    return Function(lib.epa_form_factor_monopole(lambda2))

def form_factor_dipole(lambda2):
    return Function(lib.epa_form_factor_dipole(lambda2))

def _spectrum(epa_spectrum, Z, gamma, form_factor, integrator):
    handles = []
    form_factor = _lower('epa_function1d*', form_factor, handles)
    integrator = _lower_integrator(integrator, handles)
    return Function(
            epa_spectrum(Z, gamma, form_factor, integrator),
            handles = handles
    )

def spectrum(Z, gamma, form_factor, integrator = None):
    return _spectrum(lib.epa_spectrum, Z, gamma, form_factor, integrator)

def spectrum_monopole(Z, gamma, lambda2):
    return Function(lib.epa_spectrum_monopole(Z, gamma, lambda2))

def spectrum_dipole(Z, gamma, lambda2):
    return Function(lib.epa_spectrum_dipole(Z, gamma, lambda2))

def spectrum_b(Z, gamma, form_factor, integrator = None):
    return _spectrum(lib.epa_spectrum_b, Z, gamma, form_factor, integrator)

def spectrum_b_point(Z, gamma):
    return Function(lib.epa_spectrum_b_point(Z, gamma))

def spectrum_b_monopole(Z, gamma, lambda2):
    return Function(lib.epa_spectrum_b_monopole(Z, gamma, lambda2))

def spectrum_b_dipole(Z, gamma, lambda2):
    return Function(lib.epa_spectrum_b_dipole(Z, gamma, lambda2))

def _luminosity(epa_luminosity, spectrum1, spectrum2, integrator):
    handles = []
    spectrum1, spectrum2 = _lower_spectra(
            'epa_function1d*', spectrum1, spectrum2, handles
    )
    integrator = _lower_integrator(integrator, handles)
    return Function(
            epa_luminosity(spectrum1, spectrum2, integrator),
            handles = handles
    )

def luminosity(spectrum1, spectrum2 = None, integrator = None):
    return _luminosity(lib.epa_luminosity, spectrum1, spectrum2, integrator)

def luminosity_y(spectrum1, spectrum2 = None):
    handles = []
    spectrum1, spectrum2 = _lower_spectra(
            'epa_function1d*', spectrum1, spectrum2, handles
    )
    return Function(lib.epa_luminosity_y(spectrum1, spectrum2), handles = handles)

def luminosity_fid(spectrum1, spectrum2 = None, integrator = None):
    return _luminosity(lib.epa_luminosity_fid, spectrum1, spectrum2, integrator)

def _luminosity_b(
        epa_luminosity,
        upc_probability,
        spectrum1,
        spectrum2,
        integrator_generator,
        integration_level
):
    handles = []
    upc_probability = _lower('epa_function1d*', upc_probability, handles)
    spectrum1, spectrum2 = _lower_spectra(
            'epa_function2d*', spectrum1, spectrum2, handles
    )
    generator = _lower_integrator_generator(integrator_generator, handles)
    return Function(
           epa_luminosity(
                spectrum1,
                spectrum2,
                upc_probability,
                generator.epa_function,
                integration_level
            ),
            handles = handles
    )

def luminosity_b(
        upc_probability,
        spectrum1,
        spectrum2 = None,
        integrator_generator = default_integrator,
        integration_level = 0
):
    return _luminosity_b(
            lib.epa_luminosity_b,
            upc_probability,
            spectrum1,
            spectrum2,
            integrator_generator,
            integration_level
    )

def luminosity_y_b(
        upc_probability,
        spectrum1,
        spectrum2 = None,
        integrator_generator = default_integrator,
        integration_level = 0,
        integrators = None
):
    return _luminosity_b(
            lib.epa_luminosity_y_b,
            upc_probability,
            spectrum1,
            spectrum2,
            integrator_generator,
            integration_level,
            integrators
    )

def luminosity_fid_b(
        upc_probability,
        spectrum1,
        spectrum2 = None,
        integrator_generator = default_integrator,
        integration_level = 0,
        integrators = None
):
    return _luminosity_b(
            lib.epa_luminosity_fid_b,
            upc_probability,
            spectrum1,
            spectrum2,
            integrator_generator,
            integration_level,
            integrators
    )

def xsection(photons_xsection, luminosity):
    handles = []
    photons_xsection = _lower('epa_function1d*', photons_xsection, handles)
    luminosity       = _lower('epa_function1d*', luminosity,       handles)
    return Function(
            lib.epa_xsection(photons_xsection, luminosity),
            handles = handles
    )

def xsection_b(photons_xsection, luminosity):
    handles = []
    photons_xsection = _lower('epa_xsection_b_f*',   photons_xsection, handles)
    luminosity       = _lower('epa_luminosity_b_f*', luminosity,       handles)
    return Function(
            lib.epa_xsection_b(photons_xsection, luminosity),
            handles = handles
    )

def xsection_fid(
        photons_xsection_pT,
        luminosity_fid,
        mass    = 0,
        pT_min  = 0,
        eta_max = infinity,
        w1_min  = 0,
        w1_max  = infinity,
        w2_min  = 0,
        w2_max  = infinity,
        integrator = None
):
    handles = []
    photons_xsection_pT = _lower(
            'epa_function2d*', photons_xsection_pT, handles
    )
    luminosity_fid = _lower('epa_function3d*', luminosity_fid, handles)
    integrator = _lower_integrator(integrator, handles)
    return Function(
            lib.epa_xsection_fid(
                photons_xsection_pT,
                luminosity_fid,
                mass,
                pT_min,
                eta_max,
                w1_min,
                w1_max,
                w2_min,
                w2_max,
                integrator
            ),
            handles = handles
    )

def xsection_fid_b(
        photons_xsection_pT,
        luminosity_fid,
        mass    = 0,
        pT_min  = 0,
        eta_max = infinity,
        w1_min  = 0,
        w1_max  = infinity,
        w2_min  = 0,
        w2_max  = infinity,
        integrator = None
):
    handles = []
    photons_xsection_pT = _lower(
            'epa_xsection_pT_b*', photons_xsection_pT, handles
    )
    luminosity_fid = _lower('epa_luminosity_fid_b_f*', luminosity_fid, handles)
    integrator = _lower_integrator(integrator, handles)
    return Function(
            lib.epa_xsection_fid_b(
                photons_xsection_pT,
                luminosity_fid,
                mass,
                pT_min,
                eta_max,
                w1_min,
                w1_max,
                w2_min,
                w2_max,
                integrator
            ),
            handles = handles
    )

def photons_to_fermions(mass, charge = 1):
    return Function(lib.epa_photons_to_fermions(mass, charge))

def photons_to_fermions_pT(mass, charge = 1):
    return Function(lib.epa_photons_to_fermions_pT(mass, charge))

def photons_to_fermions_b(mass, charge = 1):
    return Function(lib.epa_photons_to_fermions_b(mass, charge))

def photons_to_fermions_pT_b(mass, charge = 1):
    return Function(lib.epa_photons_to_fermions_pT_b(mass, charge))


def proton_dipole_form_factor(lambda2 = proton_dipole_form_factor_lambda2):
    return Function(lib.epa_proton_dipole_form_factor(lambda2))

def proton_dipole_spectrum(energy, lambda2 = proton_dipole_form_factor_lambda2):
    return Function(lib.epa_proton_dipole_spectrum(energy, lambda2))

def proton_dipole_spectrum_Dirac(
        energy, lambda2 = proton_dipole_form_factor_lambda2
):
    return Function(lib.epa_proton_dipole_spectrum_Dirac(energy, lambda2))

def proton_dipole_spectrum_b_Dirac(
        energy, lambda2 = proton_dipole_form_factor_lambda2
):
    return Function(lib.epa_proton_dipole_spectrum_b_Dirac(energy, lambda2))

def pp_upc_probability(collision_energy):
    return Function(lib.epa_pp_upc_probability(collision_energy))

def pp_luminosity(collision_energy, integrator = None):
    if integrator:
        handles = []
        integrator = _lower('epa_integrator*', integrator, handles)
    else:
        handles = None
        integrator = ffi.NULL
    return Function(
            lib.epa_pp_luminosity(collision_energy, integrator),
            handles = handles
    )

def pp_luminosity_y(collision_energy):
    return Function(lib.epa_pp_luminosity_y(collision_energy))

def pp_luminosity_fid(collision_energy, integrator = None):
    if integrator:
        handles = []
        integrator = _lower('epa_integrator*', integrator, handles)
    else:
        handles = None
        integrator = ffi.NULL
    return Function(
            lib.epa_pp_luminosity_fid(collision_energy, integrator),
            handles = handles
    )

def _ppx_luminosity_b(
        epa_luminosity,
        spectrum,
        spectrum_b,
        B,
        integrator_generator,
        integration_level
):
    handles = []
    if spectrum:
        spectrum = _lower('epa_function1d*', spectrum, handles)
    else:
        spectrum = ffi.NULL
    spectrum_b = _lower('epa_function1d*', spectrum_b, handles)
    generator = _lower_integrator_generator(integrator_generator, handles)
    return Function(
            epa_luminosity(
                spectrum,
                spectrum_b,
                B,
                generator.epa_function,
                integration_level
            ),
            handles = handles
    )

def ppx_luminosity_b(
        spectrum,
        spectrum_b,
        B,
        integrator_generator = default_integrator,
        integration_level = 0
):
    return _ppx_luminosity_b(
            lib.epa_ppx_luminosity_b,
            spectrum,
            spectrum_b,
            B,
            integrator_generator,
            integration_level
    )

def ppx_luminosity_y_b(
        spectrum_b,
        B,
        integrator_generator = default_integrator,
        integration_level = 0
):
    handles = []
    spectrum_b = _lower('epa_function2d*', spectrum_b, handles)
    generator = _lower_integrator_generator(integrator_generator, handles)
    return Function(
            lib.epa_ppx_luminosity_y_b(
                spectrum_b,
                B,
                generator.epa_function,
                integation_level
            ),
            handles = handles
    )

def ppx_luminosity_fid_b(
        spectrum,
        spectrum_b,
        B,
        integrator_generator = default_integrator,
        integration_level = 0
):
    return _ppx_luminosity_b(
            lib.epa_ppx_luminosity_fid_b,
            spectrum,
            spectrum_b,
            B,
            integrator_generator,
            integration_level
    )

def _pp_luminosity_b(
        epa_luminosity,
        collision_energy,
        integrator_generator,
        integration_level
):
    handles = []
    generator = _lower_integrator_generator(integrator_generator, handles)
    return Function(
            epa_luminosity(
                collision_energy, generator.epa_function, integration_level
            ),
            handles = handles
    )

def pp_luminosity_b(
        collision_energy,
        integrator_generator = default_integrator,
        integration_level = 0
):
    return _pp_luminosity_b(
            lib.epa_pp_luminosity_b,
            collision_energy,
            integrator_generator,
            integration_level
    )

def pp_luminosity_y_b(
        collision_energy,
        integrator_generator = default_integrator,
        integration_level = 0
):
    return _pp_luminosity_b(
            lib.epa_pp_luminosity_y_b,
            collision_energy,
            integrator_generator,
            integration_level
    )

def pp_luminosity_fid_b(
        collision_energy,
        integrator_generator = default_integrator,
        integration_level = 0
):
    return _pp_luminosity_b(
            lib.epa_pp_luminosity_fid_b,
            collision_energy,
            integrator_generator,
            integration_level
    )

def _pp_to_ppll(
        epa_pp_to_ppll,
        collision_energy,
        mass,
        charge,
        pT_min,
        eta_max,
        integrator_generator,
        integration_level
):
    handles = []
    generator = _lower_integrator_generator(integrator_generator, handles)
    return Function(
            epa_pp_to_ppll(
                collision_energy,
                mass,
                charge,
                pT_min,
                eta_max,
                generator.epa_function,
                integration_level
            ),
            handles = handles
    )

def pp_to_ppll(
        collision_energy,
        mass                 = 0,
        charge               = 1,
        pT_min               = 0,
        eta_max              = infinity,
        integrator_generator = default_integrator,
        integration_level    = 0
):
    return _pp_to_ppll(
            lib.epa_pp_to_ppll,
            collision_energy,
            mass,
            charge,
            pT_min,
            eta_max,
            integrator_generator,
            integration_level
    )

def pp_to_ppll_b_integrator(integration_level = 0):
    return Function(lib.epa_pp_to_ppll_b_integrator(integration_level))

def pp_to_ppll_b(
        collision_energy,
        mass                 = 0,
        charge               = 1,
        pT_min               = 0,
        eta_max              = infinity,
        integrator_generator = None,
        integration_level    = 0
):
    if not integrator_generator:
        if pT_min == 0 and eta_max == infinity:
            integrator_generator = default_integrator
        else:
            integrator_generator = pp_to_ppll_b_integrator(integration_level)
    return _pp_to_ppll(
            lib.epa_pp_to_ppll_b,
            collision_energy,
            mass,
            charge,
            pT_min,
            eta_max,
            integrator_generator,
            integration_level
    )
