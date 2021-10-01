
"""
Project: Investigating the contribution of residual unexplaind variability components in NLME approach
Program: Main program to create, estimate, and produce results  
Author: Mutaz M. Jaber <jaber038@umn.edu>
Date created: 9/5/21
Date modified: 9/28/21
"""
import os 
import sys
import subprocess

sys.path.append('src/help')
import parse
import run
import read 
print('------------------------------------------------------')
print('Investigating the contribution of residual unexplained')
print('variability components in NLME approach') 
print('------------------------------------------------------') 

if not os.path.exists('results'):
        os.mkdir('results')

TYPE=['SD1', 'SD2', 'SD3', 'SD4']
PER=['B','M', 'A1', 'A2', 'A3', 'S1', 'SL1', 'SL2', 'SL3', 'S2', 'TD1', 'TD2', 'D', 'All']
for type in TYPE:
    for P in PER:
        if P == 'M' or P == 'All':
            EXT = read.GetExt(type,P,1)
        else:
            EXT = read.GetExt(type,P, 2)

        read.GetValues(EXT, type, P) 
        read.GetShrinkage(type, P) 
