
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
TYPE=['SD1', 'SD2', 'SD3', 'SD4']
PER=['B','A1', 'A2', 'A3', 'S1', 'SL1', 'SL2', 'SL3','TD1', 'TD2', 'D', 'D2','All']
NSUBS='100'
CMT=2 
print('Creating datasets...')
os.chdir('src/sim')
with open('tmp.Rout', 'w') as tmpfile:
    for type in TYPE:
        for P in PER:
           args = [
                'Rscript', 
                '--vanilla', 
                '-e', 
                'source("generate.R")', 
                 NSUBS, NSIMS, TDOSE, DOSE, type, BASE, P
                ]
           subprocess.call(args, stderr=subprocess.STDOUT, stdout=tmpfile) 

print('Data has been generated') 
