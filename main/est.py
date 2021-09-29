
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
print('------------------------------------------------------')
print('Investigating the contribution of residual unexplained')
print('variability components in NLME approach') 
print('------------------------------------------------------') 

# os.chdir('../../')
TYPE=['int', 'spa', 'spa1', 'spa2']
PER=['B','M', 'A1', 'A2', 'A3', 'S1', 'SL1', 'SL2', 'SL3', 'S2', 'TD1', 'TD2', 'D', 'All']
print('Moving to the estimation step. This step going to take a lot of time. Relax!')
NSIMS=100
for type in TYPE:
    for P in PER:
        os.chdir(f'src/est/{type}/{P}')
        run.run(n=NSIMS)
        print('Back to Home')
        os.chdir(f'../../../../')
