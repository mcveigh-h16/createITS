# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:21:36 2022

@author: mcveigh
"""

import sys

inputfile = sys.argv[1]
outputfile = sys.argv[2]

acclist = inputfile
print ("The original list is : " +  str(acclist))

# using set()
# to remove duplicated 
# from list 
acclist = list(set(acclist))

print ("The list after removing duplicates : " + str(acclist))
