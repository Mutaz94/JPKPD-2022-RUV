#--!encoding-utf8
"""
Project: Investigating the contribution of residual unexplaind variability components in NLME approach
Program: Create NONMEM files and attach them to corresponding dataset
Author: Mutaz M. Jaber <jaber038@umn.edu>
Date created: 9/5/21
Date modified: 9/8/21
"""

import numpy
import os
import sys
import string



def create_control():
    """
    Function to create NONMEM control streams based 
    on template file. 
    
    Example
    ------------------------
    > createpy dat1 tmp1 
    # This will produce something!! 
    """
    with open(f'{sys.argv[1]}') as nmfi:
        nmfi = nmfi.readlines()
        for line in nmfi:
            print(line)




