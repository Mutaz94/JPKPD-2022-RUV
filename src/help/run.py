
"""
Project: Investigating the contribution of residual unexplaind variability components in NLME approach
Program: invoke nonmem models 
Author: Mutaz M. Jaber <jaber038@umn.edu>
Date created: 9/5/21
Date modified: 9/11/21
"""

import sys
import os 
import subprocess

CALL_NM = subprocess.call('nmfe75')
RM_FILES=['FDATA', 'FSUB', 'LINK']
RM_FILES
