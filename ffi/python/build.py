#!/usr/bin/python

import cffi
import re
import subprocess

ffi = cffi.FFI()

ffi.cdef(
        subprocess.run(
            ['cpp', '-I', '../c', '-'],
            input = '\n'.join('#include "' + module + '.h"' for module in [ 'epa', 'proton' ]),
            capture_output = True,
            text = True
        ).stdout
)

functions = { 'epa_function' }
with open('../c/epa.h') as file:
    for line in file:
        m = re.match(r'^\s*defun\s*\(\s*(\w+)', line)
        if m:
            functions.add(m[1])

ffi.cdef(r'''
enum {
    GSL_INTEG_GAUSS15,
    GSL_INTEG_GAUSS21,
    GSL_INTEG_GAUSS31,
    GSL_INTEG_GAUSS41,
    GSL_INTEG_GAUSS51,
    GSL_INTEG_GAUSS61,
    ...
};
''')

ffi.set_source(
        'epa._epa_cffi',
        r'''
#include <gsl/gsl_integration.h>
#include "../c/epa.h"
#include "../c/proton.h"
''',
        include_dirs = [ '../c' ],
        libraries = [ 'epa', 'gsl' ],
        library_dirs = [ '../..' ]
)

with open('epa/_epa_functions.py', mode = 'w') as file:
    print('functions =', functions, file = file)

with open('epa/_epa_vars.py', mode = 'w') as file:
    print('from epa._epa_cffi import lib', file = file)
    def cpp(const, var):
        with subprocess.Popen(
                [ 'cpp', '-P', '-D', const, '-D', var, '../c/epa_vars.h' ],
                text = True,
                stdout = subprocess.PIPE
        ) as cpp:
            for line in cpp.stdout:
                print(line.replace(';', ''), end = '', file = file)
    cpp(
            'defconst(type, name)=name = lib.epa_ ## name()',
            'defvar(type, name)=get_ ## name = lib.epa_get_ ## name'
    )
    cpp(
            'defconst(type, name)=',
            'defvar(type, name)=set_ ## name = lib.epa_set_ ## name'
    )

ffi.compile(verbose = True)
