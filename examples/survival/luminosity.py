#!/usr/bin/env python3

import argparse
import os
import sys
import threading

sys.path.append(os.getcwd() + '/../../ffi/python')

import epa

def usage():
    print(f'''This program calculates photon-photon luminosity in ultraperipheral proton-proton collisions
Usage: {sys.argv[0]} options... > output
Allowed options:
  -h or --help:             print this message and exit (no argument)
  -E or --collision-energy: specify proton-proton collision energy (in GeV). Default is 13e3
  -f or --from:             beginning of the invariant mass range (GeV)
  -t or --to:               ending of the invariant mass range (GeV)
  -n or --npoints:          number of points in the invariant mass range
  -s or --step:             points step in the invariant mass range (GeV)
  -l or --log:              assume logarithmic scale in the invariant mass range (no argument)
  -T or --nthreads:         use this many threads for calculation. Default is {get_nprocs()}
  -v or --verbose:          be verbose while calculating.
  -S or --survival:         account for non-electromagnetic interactions
  -p or --polarization parallel|perpendicular: assume polarized photons (no argument). Implies -S.''')

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
argparser.add_argument('-h', '--help',     action = 'store_true')
argparser.add_argument('-E', '--collision-energy', type = float, default = 13e3)
argparser.add_argument('-f', '--from',             type = float)
argparser.add_argument('-t', '--to',               type = float)
argparser.add_argument('-n', '--npoints',          type = int)
argparser.add_argument('-s', '--step',             type = float)
argparser.add_argument('-l', '--log',      action = 'store_true')
argparser.add_argument(
        '-T', '--nthreads', type = int, default = get_nprocs()
)
argparser.add_argument('-v', '--verbose',  action = 'store_true')
argparser.add_argument('-S', '--survival', action = 'store_true')
argparser.add_argument('-p', '--polarization')
args = argparser.parse_args(sys.argv[1:])

if args.help:
    usage()
    sys.exit(0)

if args.polarization:
    if args.polarization == 'parallel':
        polarization = (1, 0)
    elif args.polarization == 'perpendicular':
        polarization = (0, 1)
    else:
        die(f"Invalid polarization `{args.polarization}', expected `parallel' or `perpendicular'")
    survival = True
else:
    polarization = (1, 1)
    survival = args.survival

collision_energy = args.collision_energy

start   = getattr(args, 'from')
end     = args.to
step    = args.step
npoints = args.npoints

i = 0
if start   is not None: i += 1
if end     is not None: i += 1
if step    is not None: i += 1
if npoints is not None: i += 1
if i != 3:
    die('Exactly 3 of --from, --to, --step, --npoints must be specified')

if start is None:
    start = end - step * npoints
elif npoints is None:
    npoints = (end - start) / step + 1

if start < 0:
    die('Beginning of the invariant mass range cannot be negative')
if end <= start:
    die('End of the invariant mass range must be greater than the beginning')
if npoints < 1:
    sys.exit(0)

nthreads = args.nthreads
log = args.log
verbose = args.verbose

x = list(grid(start, end, npoints, log))

def luminosity_b():
    luminosity = epa.pp_luminosity_b(collision_energy)
    return lambda rs: luminosity(rs, polarization)

if survival:
    luminosity = make_function1d_async(
            luminosity_b,
            x,
            nthreads,
            verbose and 'luminosity_b'
    )
else:
    luminosity = make_function1d(epa.pp_luminosity(collision_energy), x)

luminosity.dump()
