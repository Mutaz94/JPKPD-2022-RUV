"""
Project: Investigating the contribution of residual unexplaind variability components in NLME approach
Program: Main program to create, estimate, and produce results  
Author: Mutaz M. Jaber <jaber038@umn.edu>
Date created: 9/5/21
Date modified: 9/20/21
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

# Adding M compartment model misspecification 
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


os.chdir('../../')
print('Moving to the estimation step. This step going to take a lot of time. Relax!')
NSIMS=100
for type in TYPE:
    for P in PER:
        os.chdir(f'src/est/{type}/{P}')
        run.run(n=NSIMS)
        print('Back to Home')
        os.chdir(f'../../../../')


# test=['B','M', 'A1', 'A2', 'A3', 'S1', 'S2', 'TD1', 'TD2', 'D', 'All']
print('Done with execution, moving to generate results')
os.mkdir('results') 

for type in TYPE:
    for P in PER:
        if P == 'M' or P == 'All':
            EXT = read.GetExt(type,P,1)
        else:
            EXT = read.GetExt(type,P, 2)

        read.GetValues(EXT, type, P) 

subprocess.call(['bash', 'clean.sh']) 
