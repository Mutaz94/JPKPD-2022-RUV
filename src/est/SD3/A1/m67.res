Sat Oct 23 21:52:56 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/SD3/A1/dat67.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2

$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr
$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      600
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      500
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m67.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1546.02694478340        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8467E+02 -9.5393E+01  1.8153E+02 -1.3203E+02  1.2115E+02  3.1585E+01 -2.1990E+01 -4.6257E+02 -3.5431E+01 -4.1447E+01
            -5.2522E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1840.68618327976        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.5129E-01  1.2374E+00  7.0049E-01  1.0449E+00  8.3956E-01  9.7013E-01  1.0283E+00  2.0343E+00  9.6822E-01  8.9868E-01
             1.2824E+00
 PARAMETER:  5.0060E-02  3.1303E-01 -2.5598E-01  1.4396E-01 -7.4876E-02  6.9679E-02  1.2790E-01  8.1013E-01  6.7708E-02 -6.8252E-03
             3.4873E-01
 GRADIENT:   1.0920E+02  1.3363E+02  4.2158E+01  1.4735E+02 -1.1442E+01  5.1489E+00 -4.3579E+00 -1.9844E+01  1.5811E+00  1.6006E+01
            -8.5185E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1846.51203074622        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      200
 NPARAMETR:  9.8434E-01  1.2048E+00  4.2541E-01  1.0404E+00  6.6119E-01  8.8001E-01  1.0823E+00  1.7858E+00  9.3109E-01  5.5685E-01
             1.2764E+00
 PARAMETER:  8.4220E-02  2.8631E-01 -7.5469E-01  1.3960E-01 -3.1371E-01 -2.7827E-02  1.7910E-01  6.7985E-01  2.8597E-02 -4.8547E-01
             3.4401E-01
 GRADIENT:  -4.0620E+01  8.3190E+01 -3.1079E+00  1.3397E+02 -3.9400E+00 -6.1095E+01  4.1194E+00 -4.8761E+01 -1.4686E+01 -2.2238E+00
            -1.0229E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1871.47737534925        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  9.9846E-01  1.3070E+00  4.5732E-01  9.0530E-01  7.5099E-01  1.0097E+00  9.5795E-01  2.2772E+00  9.6435E-01  4.2202E-01
             1.5430E+00
 PARAMETER:  9.8456E-02  3.6771E-01 -6.8238E-01  5.0667E-04 -1.8636E-01  1.0968E-01  5.7035E-02  9.2295E-01  6.3700E-02 -7.6270E-01
             5.3375E-01
 GRADIENT:  -2.8820E+00 -1.9408E+01  3.4607E+00 -1.5665E+00 -3.8235E+01  3.5269E+00  2.0032E-01 -1.8554E+01 -1.9491E+00  1.4965E+00
             2.5047E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1882.92987784646        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  9.9601E-01  1.6518E+00  5.0783E-01  7.3638E-01  9.8412E-01  9.9399E-01  9.0685E-01  2.8768E+00  1.0461E+00  6.1927E-01
             1.6230E+00
 PARAMETER:  9.5999E-02  6.0185E-01 -5.7761E-01 -2.0602E-01  8.3997E-02  9.3967E-02  2.2220E-03  1.1567E+00  1.4504E-01 -3.7922E-01
             5.8427E-01
 GRADIENT:   1.2624E+02  1.8126E+02  4.9434E+00  5.6641E+01 -2.6029E-02  1.0258E+01  8.8610E+00  4.0645E+01 -3.7265E-01  1.7925E+00
             5.5057E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1883.29862066654        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      725
 NPARAMETR:  1.0010E+00  1.6517E+00  5.0940E-01  7.3636E-01  9.8413E-01  1.0036E+00  8.5184E-01  2.8774E+00  1.1143E+00  5.7695E-01
             1.6226E+00
 PARAMETER:  1.0105E-01  6.0178E-01 -5.7452E-01 -2.0604E-01  8.4005E-02  1.0364E-01 -6.0354E-02  1.1569E+00  2.0824E-01 -4.5001E-01
             5.8402E-01
 GRADIENT:   3.5114E-01  9.4677E+00  1.1925E-01  3.4825E+01  4.4434E-01  8.0268E-01 -5.7470E-02 -6.4092E+00  9.4244E-02 -8.8803E-02
             4.6885E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1886.65411192913        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      903
 NPARAMETR:  9.9590E-01  1.6863E+00  4.9081E-01  6.8456E-01  1.0081E+00  1.0010E+00  8.5187E-01  2.8965E+00  1.1473E+00  6.8719E-01
             1.4941E+00
 PARAMETER:  9.5887E-02  6.2254E-01 -6.1169E-01 -2.7898E-01  1.0810E-01  1.0099E-01 -6.0320E-02  1.1635E+00  2.3738E-01 -2.7515E-01
             5.0152E-01
 GRADIENT:  -7.7854E+00 -2.4157E+01  1.0238E+00  4.3051E+00 -3.4548E+00 -7.1818E-01  1.8900E-01 -3.4069E+00 -7.9600E-01 -8.9058E-02
            -7.4990E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1886.90844824774        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1065
 NPARAMETR:  9.9549E-01  1.7101E+00  4.8860E-01  6.8052E-01  1.0108E+00  1.0015E+00  8.4771E-01  2.9047E+00  1.1504E+00  6.9174E-01
             1.4918E+00
 PARAMETER:  9.5482E-02  6.3658E-01 -6.1621E-01 -2.8490E-01  1.1079E-01  1.0147E-01 -6.5217E-02  1.1663E+00  2.4013E-01 -2.6854E-01
             5.0000E-01
 GRADIENT:  -8.7426E+00 -2.8735E+00  2.3521E+00  1.2220E+01 -1.0135E+01 -7.2378E-01 -1.2827E-01 -3.0830E+00 -1.5527E+00 -5.9352E-02
            -9.6908E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1887.31073605229        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1243
 NPARAMETR:  1.0030E+00  1.7067E+00  4.9533E-01  6.6945E-01  1.0304E+00  1.0052E+00  8.4449E-01  3.0208E+00  1.1683E+00  6.8976E-01
             1.5121E+00
 PARAMETER:  1.0298E-01  6.3459E-01 -6.0253E-01 -3.0129E-01  1.2992E-01  1.0521E-01 -6.9020E-02  1.2055E+00  2.5552E-01 -2.7141E-01
             5.1352E-01
 GRADIENT:   7.7215E+00 -3.2354E+01 -9.5412E-01  2.9097E+00  5.8019E-01  9.2970E-01 -3.5268E-02  1.3242E+00  9.0858E-02 -5.5934E-01
            -8.9819E-01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1887.35083649943        NO. OF FUNC. EVALS.: 153
 CUMULATIVE NO. OF FUNC. EVALS.:     1396
 NPARAMETR:  9.9963E-01  1.7039E+00  5.0521E-01  6.6712E-01  1.0374E+00  1.0030E+00  8.4619E-01  3.0653E+00  1.1673E+00  7.0811E-01
             1.5113E+00
 PARAMETER:  9.9667E-02  6.3470E-01 -5.8273E-01 -3.0460E-01  1.3691E-01  1.0328E-01 -6.7016E-02  1.2081E+00  2.5515E-01 -2.4502E-01
             5.1271E-01
 GRADIENT:   1.5246E-01  2.1717E+03  9.9681E-03  2.1890E-01  3.0305E-01  1.9101E-01 -3.0150E-03 -2.0390E+00  8.1052E-02  1.1126E-02
            -3.6308E-01
 NUMSIGDIG:         3.2         2.3         3.8         3.0         2.7         2.3         3.8         1.8         2.5         3.0
                    3.0

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1396
 NO. OF SIG. DIGITS IN FINAL EST.:  1.8

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.5498E-04 -7.7348E-03 -3.4451E-02  2.3572E-02 -4.8063E-02
 SE:             2.9717E-02  2.3795E-02  2.1045E-02  2.1320E-02  1.5581E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7705E-01  7.4514E-01  1.0163E-01  2.6888E-01  2.0373E-03

 ETASHRINKSD(%)  4.4374E-01  2.0283E+01  2.9496E+01  2.8574E+01  4.7802E+01
 ETASHRINKVR(%)  8.8552E-01  3.6453E+01  5.0292E+01  4.8984E+01  7.2754E+01
 EBVSHRINKSD(%)  7.5282E-01  1.9899E+01  3.0575E+01  3.1273E+01  4.8820E+01
 EBVSHRINKVR(%)  1.5000E+00  3.5839E+01  5.1802E+01  5.2766E+01  7.3806E+01
 RELATIVEINF(%)  9.8466E+01  1.1294E+01  1.9975E+01  8.4043E+00  8.7583E+00
 EPSSHRINKSD(%)  3.4799E+01
 EPSSHRINKVR(%)  5.7489E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1887.3508364994332     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -968.41230329476048     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1887.351       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.71E+00  5.05E-01  6.67E-01  1.04E+00  1.00E+00  8.46E-01  3.03E+00  1.17E+00  7.08E-01  1.51E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,      116.707
Stop Time:
Sat Oct 23 21:53:14 CDT 2021
