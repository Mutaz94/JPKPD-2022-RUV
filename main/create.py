
"""
Project: Investigating the contribution of residual unexplaind variability components in NLME approach
Program: Main program to create, estimate, and produce results  
Author: Mutaz M. Jaber <jaber038@umn.edu>
Date created: 9/5/21
Date modified: 9/26/21
"""
import os 
import sys
import subprocess
import glob

sys.path.append('src/help')
import parse
import run
import read 

CMT='2' 
TYPE=['SD1', 'SD2', 'SD3', 'SD4']
PER=['B','M', 'A1', 'A2', 'A3', 'S1', 'SL1', 'SL2', 'SL3', 'S2', 'TD1', 'TD2', 'D', 'All']
print('Creating control streams....')
os.chdir('src/tmp')
for type in TYPE:
    for P in PER:
        if P == 'M':
                N = len(glob.glob(f'../../data/{type}/B/*.csv'))
        else:
                N = len(glob.glob(f'../../data/{type}/{P}/*.csv'))
        for s in range(1, N+1, 1):
            if P == 'M' or P == 'All':
                parse.create_control(comp='1', design=type, per=P, n=s, dir=None)
            else:
                parse.create_control(comp=CMT, design=type,per=P,  n=s, dir=None)

