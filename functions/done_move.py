
from os.path import join, exists, isfile
from os import rmdir, makedirs
import shutil
import get_in_params as g


def done_move(mypath, cl_file, data_file, memb_file):
    '''
    Move cluster data file to 'done' dir if flag is set.
    '''

    input_dir = join(mypath, 'input')

    if g.flag_move_file:

        dst_dir = join(g.done_dir, cl_file[0])
        # If the sub-dir doesn't exist, create it before moving the file.
        if not exists(dst_dir):
            makedirs(dst_dir)
        try:
            shutil.move(data_file, dst_dir)
            print 'Photometric data file moved.'
        except:
            print "Data file already exists in 'done' destination folder."

        # Also move *memb_data.dat file if it exists.
        if isfile(memb_file):
            try:
                shutil.move(memb_file, dst_dir)
                print 'Membership data file moved.'
            except:
                print ("Membership file already exists in 'done' "
                    "destination folder.")

        # If sub-dir left behind is empty, remove it.
        if cl_file[0] != '':
            try:
                rmdir(join(input_dir, cl_file[0]))
            except OSError:
                # Sub-dir not empty, skip.
                pass
    else:
        pass