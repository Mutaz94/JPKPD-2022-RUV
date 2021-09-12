#!/usr/bin/python 

import os 
import sys
import subprocess
sys.path.append('src/help')
import parse

print('------------------------------------------------')
print('Investigating the contribution of residual unexplained')
print('variability components in NLME approach') 
print('------------------------------------------------') 

DOSE='120'
TDOSE='0'
NSIMS='100'
BASE='base.cpp'
TYPE=['int', 'spa']
PER=['B','A1', 'A2', 'A3', 'S1', 'SL1', 'SL2', 'SL3', 'S2', 'TD1', 'TD2', 'D', 'All']
NSUBS='100'
CMT=2 
print('Creating datasets...')
os.chdir('src/sim')

for type in TYPE:
    for P in PER:
       args = [
            'Rscript', 
            '--vanilla', 
            '-e', 
            'source("generate.R")', 
             NSUBS, NSIMS, TDOSE, DOSE, type, BASE, P 
            ]
       subprocess.call(args, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL) 

print('Data has been generated') 
print('Back to main dir')
os.chdir('../../')






NSUBS=100

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

