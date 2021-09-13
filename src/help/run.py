
"""
Project: Investigating the contribution of residual unexplaind variability components in NLME approach
Program: invoke nonmem models 
Author: Mutaz M. Jaber <jaber038@umn.edu>
Date created: 9/5/21
Date modified: 9/12/21
"""

import sys
import os 
import subprocess
import glob 



def run():
    CALL_NM = 'nmfe75'
    infile=[]
    outfile=[]
    ARGS='-prdefault'
    for file in 


def clean():
    RM_FILE=['FDATA.csv', 'FCON', 'LINK', 'gfortran.txt']
    pass

