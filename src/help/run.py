
"""
Project: Investigating the contribution of residual unexplaind variability components in NLME approach
Program: invoke nonmem models 
Author: Mutaz M. Jaber <jaber038@umn.edu>
Date created: 9/5/21
Date modified: 9/15/21
"""

import sys
import os 
import subprocess
import glob 


def run(n):
    CALL_NM = 'nmfe75'
    NSIMS = n 
    infile=[]
    outfile=[]
    for i in range(1, NSIMS+1):
        infile.append(f'm{i}.ctl')
        outfile.append(f'm{i}.res')

    ARGS='-prdefault'
    with open('m_run.out', 'w') as tmpfile:
        for i in range(NSIMS):
            subprocess.call([CALL_NM, infile[i], outfile[i], ARGS, '-background'], 
                    stdout=tmpfile, stderr=subprocess.STDOUT)




