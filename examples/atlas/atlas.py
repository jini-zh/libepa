#!/usr/bin/python

import argparse
import os
import sys
import threading

sys.path.append(os.getcwd() + '/../../ffi/python')

import epa

def usage():
    print('''This program calculates the fiducial differential cross section for muon pair production in ultraperipheral proton-proton collisions. The fiducial region is defined by the cuts on transverse momentum pT and pseudorapidity eta of each muon:
pT >  6 GeV for 12 GeV < sqrt(s) < 30 GeV,
pT > 10 GeV for 30 GeV < sqrt(s) < 70 GeV;
abs(eta) < 2.4.
Parameters of this calculation are the same as were set in the ATLAS experiment for the paper 1708.04053.
Usage: {:s} options...
Allowed options:
  -h or --help:     print this message
  -n or --npoints:  number of points in the invariant mass range
  -s or --step:     points step in the invariant mass range (GeV)
  -v or --verbose:  be verbose while calculating
  -S or --survival: account for non-electromagnetic interactions
'''.format(sys.argv[0])
)

class Function1d:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __call__(self, x):
        i = bisect(self.x, x)
        return self.y[i] \
             + (self.y[i+1] - self.y[i]) \
               * (x - self.x[i]) / (self.x[i+1] - self.x[i])

    def dump(self, output = sys.stdout):
        for x, y in zip(self.x, self.y):
            print('{:19.12e} {:19.12e}'.format(x, y))

def bisect(array, value):
    if value < array[0] or value > array[-1]:
        return None
    left  = 0
    right = len(array)
    while right - left > 1:
        middle = (left + right) // 2
        if value < array[middle]:
            right = middle
        else:
            left = middle
    return left

def grid(start, end, n, log = False):
    n -= 1
    x = start
    if log:
        step = (end / start) ** (1 / n)
        for i in range(n):
            yield x
            x *= step
    else:
        step = (end - start) / n
        for i in range(n):
            yield x
            x += step
    yield end

def make_function1d(f, x):
    return Function1d(x, list(map(f, x)))

def make_function1d_async(generator, x, nthreads = None, verbose = None):
    if not nthreads:
        nthreads = len(os.sched_getaffinity(0))
    if verbose:
        print(f'Using {nthreads} threads', file = sys.stderr)

    y = [0] * len(x)
    index = 0

    m_queue = threading.Lock()
    m_output = threading.Lock()
    def calculate():
        f = generator()
        nonlocal index
        while index < len(x):
            with m_queue:
                i = index
                index += 1
            y[i] = f(x[i])
            if verbose:
                with m_output:
                    print(verbose, x[i], '=>', y[i], file = sys.stderr)

    threads = []
    for i in range(nthreads):
        thread = threading.Thread(target = calculate, daemon = True)
        thread.start()
        threads.append(thread)
    for thread in threads:
        thread.join()

    return Function1d(x, y)

muon_mass = 105.6583745e-3
def atlas_xsection(epa_xsection):
    x0 = epa_xsection(13e3, mass = muon_mass, pT_min =  6, eta_max = 2.4)
    x1 = epa_xsection(13e3, mass = muon_mass, pT_min = 10, eta_max = 2.4)
    return lambda rs: x0(rs) if rs < 30 else x1(rs)

argparser = argparse.ArgumentParser(add_help = False)
argparser.add_argument('-h', '--help',     action = 'store_true')
argparser.add_argument('-n', '--npoints',  type = int)
argparser.add_argument('-s', '--step',     type = float)
argparser.add_argument('-v', '--verbose',  action = 'store_true')
argparser.add_argument('-S', '--survival', action = 'store_true')

args = argparser.parse_args(sys.argv[1:])

if args.help:
    usage()
    sys.exit(0)

if args.npoints is not None and args.step is not None:
    print('--npoints and --step parameters are incompatible', file = sys.stderr)
    sys.exit(1)

if args.npoints is None:
    if args.step is None:
        args.npoints = 100
    else:
        args.npoints = (70 - 30) / step

x = list(grid(12, 70, args.npoints))
i = bisect(x, 30)
if x[i] != 30:
    x.insert(i + 1, 30)

epa.set_default_error_step(1 / 3)

if args.survival:
    xsection = make_function1d_async(
            lambda: atlas_xsection(epa.pp_to_ppll_b),
            x,
            verbose = 'xsection_b' if args.verbose else None
    )
else:
    xsection = make_function1d(atlas_xsection(epa.pp_to_ppll), x)

xsection.dump()

if args.verbose:
    integrate = epa.qag_integrator(relative_error = 1e-2)
    print(
            'Integrated cross section:',
            integrate(xsection, x[0], 30) + integrate(xsection, 30, x[-1])
    )
