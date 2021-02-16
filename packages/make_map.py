
import pathlib
from pathlib import PurePath

"""
Print the map of files being called throughout the process.
"""

not_packs = (
    'os', 'sys', 'gc', 'time', 'pickle', 'argparse', 'traceback', 'astropy',
    'numpy', 'scipy', 'matplotlib', 'emcee', 'warnings', 'itertools',
    'urllib', 'collections', 'shutil', 'multiprocessing', '_version',
    'prep_plots', 'subprocess', 'functools',)

mypath = pathlib.Path().absolute()


def printImport(N, file):
    """
    """
    line = str(file).split('/')
    idx = line.index('packages') + 1
    print("    " * N, "{}-->".format(N), '/'.join(line[idx:])[:-3])


def fileCheck(line):
    """
    """
    try:
        f0 = PurePath(mypath, '/'.join(line) + '.py')
        with open(f0):
            pass
        return f0
    except:
        pass

    try:
        f1 = PurePath(mypath, '/'.join(line[1:]) + '.py')
        with open(f1):
            pass
        return f1
    except:
        pass

    try:
        f2 = PurePath(mypath, line[0] + '.py')
        with open(f2):
            pass
        return f2
    except:
        pass


def allImports(file):
    file_imports = []
    with open(file, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if any(x in line for x in not_packs):
                continue
            if line.startswith('from'):
                line = line.replace('packages', '').replace(
                    'from', '').replace('import', '').replace(
                    '.', ' ').replace(',', ' ').split()

                line0 = str(file).split('/')
                idx = line0.index('packages') + 1
                line = line0[idx:-1] + line

                fline = fileCheck(line)
                if fline is not None:
                    file_imports.append(fline)

    return file_imports


print("/ asteca_run")
file_imports = allImports(PurePath(mypath, 'asteca_run.py'))

for file in file_imports:
    print("")
    printImport(1, file)
    for file in allImports(file):
        printImport(2, file)
        for file in allImports(file):
            printImport(3, file)
            for file in allImports(file):
                printImport(4, file)
                for file in allImports(file):
                    printImport(5, file)
