import setuptools

setuptools.setup(
        name = 'libepa',
        version = '1.0.0',
        description = 'Library for calculations of cross sections of ultraperipheral collisions of high energy particles under the equivalent photons approximation',
        url = 'https://github.com/jini-zh/libepa',
        author = 'Evgenii Zhemchugov',
        author_email = 'jini.zh@gmail.com',
        license = 'GPLv3',
        packages = [ 'epa' ],
        cffi_modules = ["./build.py:ffi"]
)
