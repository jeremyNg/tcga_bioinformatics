#!/usr/local/bin/python3

# load the modules first
import numpy as np
import os
import os.path
import sys

if len(sys.argv)<2:
    print ("Not enough arguments\nUsage: pythonscript.py [file]") # print error msg if no path is specified
else:
    filepath=sys.argv[1] # assigns variable as file path
    try: # to check if the file that is selected does exist
        with open(filepath) as file:
            pass
    except IOError as e:
        print ("File %s does not exist" %sys.argv[1])
        exit(0) # exits the script if the file does not exist
#print ("Reading in data from %s" %filepath)
rnaseq=np.genfromtxt(filepath,delimiter=",") # reads in the array using numpy
#print ("Dimensions of array extent: %d  genes and %d samples" %rnaseq.shape)

exit(0) # exits the script
