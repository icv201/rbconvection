'''
Instructions for Use:
Change the base direc to be where you are working from in the terminal - it should equal thisdir as a sanity check.
results_direc is the folder name that the raw_data2 folder is saved into. 
'''

from dedalus import public as de
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import h5py
import numpy as np
import sys
import pathlib
import os
from dedalus.tools import post

# Getting the current work directory (cwd)
thisdir = os.getcwd()

#base_direc = "Users/tomaddison/Desktop/astro_project/MattExplainedCode/"
results_direc = "Results3j"

# This is the directory ending (including) raw_data2. Change raw_data2 if that is not the folder where the results are. 
direc = os.path.join(thisdir, results_direc, "raw_data2")


# Change this to analysis, run_parameters and snapshots depending on where we are at. 
section = "run_parameters"

specificdir =os.path.join(direc, section)
# r=root, d=directories, f = files
list_of_run_parameter_files = []
for r, d, f in os.walk(specificdir): #you need thisdir not specificdir
    for file in f:
        if file.endswith(".h5"):
            list_of_run_parameter_files.append(os.path.join(r,file))
            print(os.path.join(r, file))

# Useful for error hunting if the loop above doesn't appear to work as it should
#print(list_of_run_parameter_files)

# If there are 0 or 1 files nothing needs to be merged. But useful to know how many files there are. 
if len(list_of_run_parameter_files) == 1:
    final_file = list_of_run_parameter_files[0]
    print(final_file)
    # are we just going to pass?
elif len(list_of_run_parameter_files) == 0:
    print("No Run Parameter files found")
else:
    # Do the merge for run_parameters here.
    final_file = post.merge_sets(results_direc + "/raw_data2/" + section + ".h5", list_of_run_parameter_files, cleanup=True) 


section = "analysis"
specificdir =specificdir =os.path.join(direc, section)
# r=root, d=directories, f = files
list_of_analysis_files = []
for r, d, f in os.walk(specificdir): #you need thisdir not specificdir
    for file in f:
        if file.endswith(".h5"):
            list_of_analysis_files.append(os.path.join(r,file))
            print(os.path.join(r, file))

# Useful for error hunting if the loop above doesn't appear to work as it should
#print(list_of_analysis_files)

# If there are 0 or 1 files nothing needs to be merged. But useful to know how many files there are. 
if len(list_of_analysis_files) == 1:
    final_analysis_file = list_of_analysis_files[0]
    print(final_analysis_file)
elif len(list_of_analysis_files) == 0:
    print("No analysis files found")
else:
    # Do the merge for run_parameters here.
    final_analysis_file = post.merge_sets(results_direc + "/raw_data2/" + section + ".h5", list_of_analysis_files, cleanup=True)



section = "snapshots"
specificdir =specificdir =os.path.join(direc, section)
# r=root, d=directories, f = files
list_of_snapshot_files = []
for r, d, f in os.walk(specificdir): #you need thisdir not specificdir
    for file in f:
        if file.endswith(".h5"):
            list_of_snapshot_files.append(os.path.join(r,file))
            print(os.path.join(r, file))

# Useful for error hunting if the loop above doesn't appear to work as it should
# print(list_of_snapshot_files)

# If there are 0 or 1 files nothing needs to be merged. But useful to know how many files there are. 
if len(list_of_snapshot_files) == 1:
    final_snapshot_file = list_of_snapshot_files[0]
    print(final_snapshot_file)
elif len(list_of_snapshot_files) == 0:
    print("No snapshot files found")
else:
    # Do the merge for run_parameters here.
    final_snapshot_file = post.merge_sets(results_direc + "/raw_data2/" + section + ".h5", list_of_snapshot_files, cleanup=True) 
