
import os
from shutil import copyfile
import traceback


def check(mypath):
    """
    Check if this is the first run of the code or not. If it is, attempt to
    copy input .dat files from packages/ folder into root folder.
    """

    try:
        fr_file = mypath + 'packages/.first_run'
        if os.path.isfile(fr_file):

            success = 0
            with open(fr_file) as f:
                if f.read() == 'y':
                    print("* First run of the code detected *\n")

                    print("Copy input .dat files to the root folder.\n")
                    # Check if the input .dat files already exist.
                    for i_f in ['params', 'semi']:
                        i_f_path = mypath + i_f + '_input.dat'
                        o_f_path = mypath + i_f + '_input_OLD.dat'
                        n_f_path = mypath + 'packages/.' + i_f + '_input.dat'
                        # First, check if the file exists in packages/
                        if os.path.isfile(n_f_path):
                            # Now check if it exists in the root folder.
                            if os.path.isfile(i_f_path):
                                # File exist in root folder. Rename it.
                                print("File {}\nalready present in"
                                      " folder.".format(i_f_path))
                                print("Rename to {}.".format(
                                    i_f + '_input_OLD.dat'))
                                os.rename(i_f_path, o_f_path)
                            # Copy into root folder.
                            print("Copy {} to root folder.\n".format(
                                i_f + '_input.dat'))
                            copyfile(n_f_path, i_f_path)
                            success += 1
                        else:
                            print("WARNING: file {}\nis not present.\n".format(
                                n_f_path))

            # If both files could be copied to the root folder,
            # change .first_run file so it won't be done again.
            if success == 2:
                with open(fr_file, 'w') as fw:
                    fw.write("n")
        else:
            print("WARNING: No 'packages/.first_run' file present.\n")
    except:
        print("WARNING: Could not check/copy input .dat files\n"
              "to root folder.\n")
        print(traceback.format_exc())
