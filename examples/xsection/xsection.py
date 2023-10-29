#!/usr/bin/python

import argparse
import os
import sys
import threading

sys.path.append(os.getcwd() + '/../../ffi/python')

import epa

def usage():
    print(f'''This program plots the cross section for the production of a pair of fermions in ultraperipheral proton-proton collisions with respect to the fermion mass
Usage: {sys.argv[0]} options...
Allowed options:
  -h or --help: print this message and exit (no argument)
  -E or --collision-energy: specify proton-proton collision energy (in GeV). Default is 13e3
  -f or --from:             beginning of the mass range (GeV)
  -t or --to:               ending of the mass range (GeV)
  -n or --npoints:          points step in the mass range (GeV)
  -l or --log:              assume logarithmic scale in the mass range (no argument)
  -v or --verbose:          be verbose while calculating
  -S or --survival:         account for non-electromagnetic interactions
  --rs-max:                 cut off integration over invariant mass (sqrt(s)) at this value (GeV)''')

def die(message):
    print(message, file = sys.stderr)
    sys.exit(1)

def get_nprocs():
    return len(os.sched_getaffinity(0))

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
    if n < 2:
        if n == 1:
            yield start
        return

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
        nthreads = get_nprocs()
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

argparser = argparse.ArgumentParser(add_help = False)
argparser.add_argument('-h', '--help',    action = 'store_true')
argparser.add_argument('-E', '--collision-energy', type = float, default = 13e3)
argparser.add_argument('-f', '--from',             type = float)
argparser.add_argument('-t', '--to',               type = float)
argparser.add_argument('-n', '--npoints',          type = int)
argparser.add_argument('-s', '--step',             type = float)
argparser.add_argument('-l', '--log',      action = 'store_true')
argparser.add_argument('-v', '--verbose',  action = 'store_true')
argparser.add_argument('-S', '--survival', action = 'store_true')
argparser.add_argument('--rs-max',                 type = float)
args = argparser.parse_args(sys.argv[1:])

if args.help:
    usage()
    sys.exit(0)

args.start = getattr(args, 'from')
args.end   = args.to

i = 0
for value in [ 'start', 'end', 'step', 'npoints' ]:
    if getattr(args, value) is not None:
        i += 1
if i != 3:
    die('Exactly 3 of --from, --to, --step, --npoints must be specified')

if args.start is None:
    args.start = args.end - args.step * args.npoints
elif args.npoints is None:
    args.npoints = (args.start - args.end) / args.step + 1

if args.start < 0:
    die('Beginning of the invariant mass range cannot be negative')
if args.end <= args.start:
    die('End of the invariant mass range must be greater than the beginning')
if args.npoints < 1:
    sys.exit(0)

if args.rs_max is None:
    args.rs_max = args.collision_energy / 2

epa.set_default_relative_error(1e-2)

x = list(grid(args.start, args.end, args.npoints, args.log))

if args.survival:
    def xsection_():
        integrate = epa.default_integrator(0)
        luminosity = epa.pp_luminosity_b(
                args.collision_energy, integration_level = 1
        )
        return lambda mass: integrate(
                epa.xsection_b(epa.photons_to_fermions_b(mass), luminosity),
                2 * mass,
                args.rs_max
        )
    xsection = make_function1d_async(
            xsection_, x, verbose = args.verbose and 'cross section'
    )
else:
    integrate = epa.default_integrator(0)
    luminosity = epa.pp_luminosity(
            args.collision_energy, integrator = epa.default_integrator(1)
    )
    xsection = make_function1d(
            lambda mass: integrate(
                epa.xsection(epa.photons_to_fermions(mass), luminosity),
                2 * mass,
                args.rs_max
            ),
            x
    )

xsection.dump()
