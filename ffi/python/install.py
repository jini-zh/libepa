import os
import os.path
import shutil
import sys
import sysconfig

action      = sys.argv[1]
prefix      = sys.argv[2] if len(sys.argv) > 2 else '/usr'
destination = sysconfig.get_path('platlib', vars = { 'platbase': prefix })
destination = os.path.join(destination, 'epa')

if action == 'u':
    try:
        shutil.rmtree(destination)
        print('removed', destination)
    except FileNotFoundError:
        pass
elif action == 'i':
    source      = os.path.dirname(sys.argv[0])
    source      = os.path.join(source, 'epa')

    try:
        os.makedirs(destination)
    except FileExistsError:
        pass

    files = [
            '__init__.py',
            '_epa_cffi.so',
            '_epa_functions.py',
            '_epa_vars.py',
            'epa.py'
    ]

    for file in files:
        src = os.path.join(source,      file)
        dst = os.path.join(destination, file)
        print(src, '->', dst)
        shutil.copyfile(src, dst)
else:
    raise Exception(sys.argv[0] + ': invalid action: ' + action)
