
"""
Project: Investigating the contribution of residual unexplaind variability components in NLME approach
Program: Read ext file and generate result table 
Author: Mutaz M. Jaber <jaber038@umn.edu>
Date created: 9/18/21
Date modified: 9/19/21
"""
import pandas
import os 
import sys
import glob
import statistics
import math 

from scipy import stats


def GetExt(design, model, comp):
    Final = []
    if comp == 2:
        names = ['CL', 'V', 'Q', 'Vp', 'Ka', 'BCL', 'BV', 'BQ', 'BVp', 'BKa','RUV']
    else:
        names = ['CL', 'V', 'Ka','BCL', 'BV', 'BKa','RUV']
        
    for file in glob.glob(f'src/est/{design}/{model}/*.ext'):
        dat = pandas.read_table(file, skiprows=1, delim_whitespace=True)
        dat = dat[dat['ITERATION'] == -1E9]
        thetaNames = [col for col in dat.columns if 'THETA' in col]
        EST = dat[thetaNames]
        EST.columns = names
        Final.append(EST)
    ESTi = pandas.concat(Final)
    return ESTi

def GetValues(value, design, model):
    Med= []
    Max = []
    Min = []
    up95, lo95 = [], []
    for nam in value.columns:
        Med.append(statistics.median(value[nam]))
        Max.append(max(value[nam]))
        Min.append(min(value[nam]))
        up95.append(statistics.mean(value[nam] + stats.norm.ppf(0.975) * statistics.stdev(value[nam])/math.sqrt(len(value[nam]))))
        lo95.append(statistics.mean(value[nam] - stats.norm.ppf(0.975) * statistics.stdev(value[nam])/math.sqrt(len(value[nam]))))
    Data = pandas.DataFrame({'Median': Med, 'Min': Min, 'Max': Max, '95CIlo': lo95, '95CIup': up95},value.columns)
    Data.to_csv(f'results/{design}-{model}.csv')
    return Data


# def CSV(Data, design, model):
#     Data.to_csv(f'{design}-{model}.csv', index=False)

#    # Create results directory
#    # Dump files in that directory 
