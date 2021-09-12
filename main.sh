#!/usr/bin/sh

echo ------------------------------------------------------- \n
echo Investigating the contribution of residual unexplained \n
echo variability components in NLME approach \n
echo ------------------------------------------------------- \n

# declare -a CMT=(1 2) 
declare -a TYPE=("int" "spa")
declare -a PER=("B" "A1" "A2" "A3" "S1" "SL1" "SL2" "SL3" "S2" "TD1" "TD2" "D" "All") 
NSUBS=100
CMT=2
echo Creating control streams...
cd src/tmp 

for c in ${CMT[@]}; do
    for i in ${TYPE[@]}; do 
        for p in ${PER[@]}; do 
            for s in $(seq 1 $NSUBS); do 
             python ../help/parse.py $c $i $p $s
            done
        done 
    done
done




