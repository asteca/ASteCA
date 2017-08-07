
from os.path import exists, isfile
from os import makedirs
import shutil


def main(pd, dst_dir, data_file, memb_file, **kwargs):
    '''
    Move cluster data file to 'done' dir if flag is set.
    '''

    if pd['flag_move_file']:

        # If the sub-dir doesn't exist, create it before moving the file.
        if not exists(dst_dir):
            makedirs(dst_dir)
        try:
            shutil.move(data_file, dst_dir)
            print('Photometric data file moved.')
        except:
            print("Data file already exists in 'done' destination folder.")

        # Also move *memb_data.dat file if it exists.
        if isfile(memb_file):
            try:
                shutil.move(memb_file, dst_dir)
                print('Membership data file moved.')
            except:
                print("Membership file already exists in 'done' "
                      "destination folder.")
