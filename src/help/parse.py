#--!encoding-utf8
"""
Project: Investigating the contribution of residual unexplaind variability components in NLME approach
Program: Create NONMEM files and attach them to corresponding dataset
Author: Mutaz M. Jaber <jaber038@umn.edu>
Date created: 9/5/21
Date modified: 9/11/21
"""

import numpy
import os
import sys
import string



def create_control(comp, design, per, n, dir=None):
    """
    Function to create NONMEM control streams based 
    on template file. 
    
    Example
    ------------------------
    > createpy dat1 tmp1 
    # This will produce something!! 
    """
    if comp == '1':
        val = 'tmp2.ctl'
    else:
        val = 'tmp1.ctl' 
    
    with open(f'{val}') as nmfi:
        nmfi = nmfi.readlines()
    
        # design=sys.argv[2]
        # per=sys.argv[3]
        # nsub=sys.argv[4]
        newfile=""
        for line in nmfi:
            sl = line.strip()
            if per == 'M':
                newline=sl.replace('../../data/SD1/B/dat1.csv', f'../../../../data/{design}/B/dat{n}.csv')
            else:
                newline=sl.replace('../../data/SD1/B/dat1.csv', f'../../../../data/{design}/{per}/dat{n}.csv')

            newfile += newline + '\n'
        
        # Create directory
        dir = f'../est/{design}/{per}/'
        os.makedirs(dir, exist_ok=True) 
        control = open(f'{dir}/m{n}.ctl', 'w')
        control.write(newfile)
        control.close() 
        print(f'Model {design}-{per}-{n} has been created')
    
# if __name__ == '__main__':
#     create_control(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]) 


