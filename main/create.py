
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

sys.path.append('src/help')
import parse
import run
import read 

CMT='2' 
NSUBS=100
TYPE=['int', 'spa', 'spa1', 'spa2']
PER=['B','M', 'A1', 'A2', 'A3', 'S1', 'SL1', 'SL2', 'SL3', 'S2', 'TD1', 'TD2', 'D', 'All']
print('Creating control streams....')
os.chdir('src/tmp')
for type in TYPE:
    for P in PER:
        for s in range(1, NSUBS+1, 1):
            if P == 'M' or P == 'All':
                parse.create_control(comp='1', design=type, per=P, nsub=s, dir=None)
            else:
                parse.create_control(comp=CMT, design=type,per=P,  nsub=s, dir=None)

