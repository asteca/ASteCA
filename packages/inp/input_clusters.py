
from os.path import join
from os import walk
import re


def main(mypath, file_end):
    """
    Store the paths and names of all the input clusters stored in the
    'input/' folder.
    """

    # Path where cluster data files are stored.
    input_dir = 'input' + file_end
    full_in_dir = join(mypath, input_dir)

    # Store subdir names and file names in cl_files.
    cl_files = []
    for root, dirs, files in walk(full_in_dir):

        # Don't read this sub-folder so it can be used as a container.
        if 'dont_read' in dirs:
            dirs.remove('dont_read')

        for f in files:
            # Remove input_dir from sub-dir path.
            subdir0 = re.sub(r'^' + re.escape(full_in_dir), '', root)
            if subdir0.startswith('/'):
                # Remove possible extra '/' char at beginning.
                subdir = subdir0[1:]
            else:
                subdir = subdir0
            # Don't attempt to read certain files
            if not f.endswith(('.md', '.pickle')) and not\
                    f.startswith('.'):
                cl_files.append([mypath, input_dir, subdir, f])

    # Return sorted list by cluster file name.
    cl_files.sort(key=lambda x: x[-1].lower())

    return cl_files
