
import os
from shutil import copyfile, copytree
import traceback


def main(mypath):
    """
    Check if this is the first run of the code or not. If it is, attempt to
    copy input .dat files and input/, isochrones/ folders from packages/ folder
    into root folder. If these files/folders already exist, halt.
    """

    # Generate output dir (if it doesn't exist).
    if not os.path.exists('output/'):
        os.mkdir('output/')

    first_run_flag = False
    try:
        fr_file = mypath + 'packages/.first_run'
        if os.path.isfile(fr_file):
            # Read .first_run file.
            with open(fr_file) as f:
                N = str(f.read()).rstrip()
            # If this is the first run.
            if N == '0':
                print("* First run of the code detected *\n")
                files_folders_check(mypath)
                # Change value to indicate that the first run was processed.
                with open(fr_file, 'w') as fw:
                    fw.write("1")
                first_run_flag = True
        else:
            print("ERROR: No 'packages/.first_run' file present.\n")
    except Exception:
        print(traceback.format_exc())
        print("ERROR: Could not check/copy input data files and folders.\n")

    if first_run_flag:
        print("First run processes completed.\n")

    return first_run_flag


def files_folders_check(mypath):
    """
    Check if input files and input folders already exist in root folder.
    """
    files_lst = ['params_input.dat', 'CMD_systs.dat']
    folds_lst = ['isochrones', 'input']

    success = 0
    # Check files
    for i_f in files_lst:
        i_f_path = mypath + i_f
        n_f_path = mypath + 'packages/defvals/' + i_f
        # First, check if the file exists in packages/ (it should)
        if os.path.isfile(n_f_path):
            # Now check if it exists in the root folder.
            if os.path.isfile(i_f_path):
                # File exist in root folder. Rename it.
                raise ValueError(
                    "File '{}' is already present in root folder.\nPlease"
                    " rename/delete it and try again.\n".format(i_f))
            else:
                success += 1
        else:
            print("WARNING: file {}\nis not present.\n".format(n_f_path))
    # Check folders.
    for i_fo in folds_lst:
        i_fo_path = mypath + i_fo + '/'
        n_fo_path = mypath + 'packages/defvals/' + i_fo + '/'
        # First, check if the folder exists in packages/ (it should)
        if os.path.isdir(n_fo_path):
            # Now, check if the folder exists in the root.
            if os.path.isdir(i_fo_path):
                raise ValueError(
                    "Folder '/{}' is already present in root folder.\n"
                    "Please rename/delete it and try again\n.".format(i_fo))
            else:
                success += 1
        else:
            print("WARNING: folder {}\nis not present.\n".format(n_fo_path))

    if success == (len(files_lst) + len(folds_lst)):
        # Copy files into root folder.
        for i_f in files_lst:
            i_f_path = mypath + i_f
            n_f_path = mypath + 'packages/defvals/' + i_f
            copyfile(n_f_path, i_f_path)
            print("File {} copied into root folder".format(i_f))
        # Copy folders into root folder.
        for i_fo in folds_lst:
            i_fo_path = mypath + i_fo + '/'
            n_fo_path = mypath + 'packages/defvals/' + i_fo + '/'
            # Copy folders into root folder.
            copytree(n_fo_path, i_fo_path)
            print("Folder {}/ copied into root folder".format(i_fo))
        print("")
