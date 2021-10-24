Sat Oct 23 13:20:29 CDT 2021
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
$DATA ../../../../data/SD1/B/dat38.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3473.71692737842        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6179E+02  7.6035E+01  1.2818E+02  4.6581E+01  5.7925E+01  1.5623E+01 -4.4018E+01 -4.4156E+02 -1.2320E+02 -5.6294E+00
            -8.9409E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3672.05815609147        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      171
 NPARAMETR:  9.4162E-01  1.0863E+00  9.2562E-01  9.7896E-01  1.1340E+00  1.0548E+00  1.2057E+00  3.3909E+00  1.1842E+00  9.3157E-01
             1.0260E+00
 PARAMETER:  3.9846E-02  1.8278E-01  2.2713E-02  7.8733E-02  2.2576E-01  1.5339E-01  2.8707E-01  1.3211E+00  2.6910E-01  2.9114E-02
             1.2567E-01
 GRADIENT:  -9.9532E+01 -2.3574E+01 -2.8166E+01  7.3936E+00 -2.9604E+01 -1.5576E+01 -7.2477E+00 -1.7015E+01  3.6462E+01 -2.5534E+00
             2.5541E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3677.75310712299        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      357
 NPARAMETR:  9.5472E-01  1.0972E+00  9.8022E-01  9.6022E-01  1.1665E+00  1.0219E+00  1.2255E+00  3.4935E+00  1.1307E+00  9.5185E-01
             1.0143E+00
 PARAMETER:  5.3659E-02  1.9275E-01  8.0017E-02  5.9411E-02  2.5403E-01  1.2166E-01  3.0335E-01  1.3509E+00  2.2283E-01  5.0655E-02
             1.1418E-01
 GRADIENT:  -7.6926E+01 -3.0939E+01 -2.1527E+01 -8.8145E+00 -1.4802E+01 -2.7186E+01 -6.2029E+00 -1.3042E+01  2.9672E+01 -2.0779E+00
             1.1085E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3678.70060086045        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:      505
 NPARAMETR:  9.5534E-01  1.0984E+00  9.8409E-01  9.5916E-01  1.1669E+00  1.0861E+00  1.2257E+00  3.4670E+00  1.1306E+00  9.5393E-01
             1.0149E+00
 PARAMETER:  5.4313E-02  1.9382E-01  8.3962E-02  5.8300E-02  2.5432E-01  1.8255E-01  3.0349E-01  1.3433E+00  2.2279E-01  5.2837E-02
             1.1481E-01
 GRADIENT:   3.7932E+02  1.6276E+02 -1.5280E+01  7.6807E+01  9.6736E+01  1.1702E+02  4.1985E+01  4.4980E+01  4.5735E+01  1.9711E+00
             4.3721E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3679.59276717639        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      651
 NPARAMETR:  9.5726E-01  1.1003E+00  9.9765E-01  9.5953E-01  1.1677E+00  1.0850E+00  1.2280E+00  3.4648E+00  1.1141E+00  9.5505E-01
             1.0148E+00
 PARAMETER:  5.6321E-02  1.9563E-01  9.7643E-02  5.8686E-02  2.5500E-01  1.8158E-01  3.0535E-01  1.3426E+00  2.0800E-01  5.4009E-02
             1.1469E-01
 GRADIENT:  -6.2974E+01 -2.8408E+01 -1.9045E+01 -9.4626E+00 -1.6761E+01 -1.0288E+00 -5.9897E+00 -1.5672E+01  2.6694E+01 -1.7501E+00
             1.8013E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3681.81556263137        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      833
 NPARAMETR:  9.5750E-01  1.0998E+00  9.9741E-01  9.5930E-01  1.1670E+00  1.0892E+00  1.3147E+00  3.4547E+00  9.4926E-01  9.7781E-01
             1.0151E+00
 PARAMETER:  5.6567E-02  1.9517E-01  9.7408E-02  5.8450E-02  2.5441E-01  1.8547E-01  3.7363E-01  1.3397E+00  4.7932E-02  7.7563E-02
             1.1496E-01
 GRADIENT:  -6.2111E+01 -2.1624E+01 -1.9891E+01 -1.7561E+01 -2.5903E+01  6.5489E-01 -2.4175E-01 -2.0218E+01  4.3335E-01  3.4824E-01
             4.2790E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3681.88062456737        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1023
 NPARAMETR:  9.5754E-01  1.0999E+00  9.9948E-01  9.5931E-01  1.1670E+00  1.0893E+00  1.3210E+00  3.4590E+00  9.4274E-01  9.8344E-01
             1.0151E+00
 PARAMETER:  5.6616E-02  1.9520E-01  9.9484E-02  5.8454E-02  2.5444E-01  1.8550E-01  3.7839E-01  1.3410E+00  4.1035E-02  8.3297E-02
             1.1497E-01
 GRADIENT:  -6.1988E+01 -2.1524E+01 -1.9549E+01 -1.7708E+01 -2.7034E+01  6.9602E-01  3.1023E-01 -2.0015E+01 -1.9591E-01  9.5602E-01
             4.5582E+00

0ITERATION NO.:   32    OBJECTIVE VALUE:  -3681.91230148143        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1093
 NPARAMETR:  9.5753E-01  1.0999E+00  1.0011E+00  9.5930E-01  1.1669E+00  1.0892E+00  1.3208E+00  3.4573E+00  9.4270E-01  9.8245E-01
             1.0151E+00
 PARAMETER:  5.6616E-02  1.9520E-01  1.0111E-01  5.8454E-02  2.5444E-01  1.8550E-01  3.7840E-01  1.3410E+00  4.1034E-02  8.3297E-02
             1.1497E-01
 GRADIENT:   1.9669E+04  2.8361E+03 -5.6953E+03  5.5764E+03  1.3148E+04  7.5357E+03  6.6902E+03  2.4474E+03  2.5326E+04  9.3345E-01
            -2.9195E+04
 NUMSIGDIG:         2.8         3.1         3.2         3.1         2.3         2.3         2.3         2.3         2.3         0.9
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1093
 NO. OF SIG. DIGITS IN FINAL EST.:  0.9
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1544E-02 -1.8405E-02  5.5869E-03  3.2830E-02 -4.9122E-02
 SE:             2.9758E-02  2.4874E-02  2.7481E-02  2.5416E-02  2.1745E-02
 N:                     100         100         100         100         100

 P VAL.:         2.8914E-01  4.5934E-01  8.3890E-01  1.9645E-01  2.3883E-02

 ETASHRINKSD(%)  3.0616E-01  1.6669E+01  7.9348E+00  1.4855E+01  2.7152E+01
 ETASHRINKVR(%)  6.1138E-01  3.0559E+01  1.5240E+01  2.7503E+01  4.6932E+01
 EBVSHRINKSD(%)  2.3534E-01  1.6643E+01  1.1044E+01  1.7563E+01  2.7505E+01
 EBVSHRINKVR(%)  4.7013E-01  3.0516E+01  2.0869E+01  3.2041E+01  4.7445E+01
 RELATIVEINF(%)  9.9529E+01  4.0762E+01  7.3132E+01  4.2862E+01  3.2310E+01
 EPSSHRINKSD(%)  2.3459E+01
 EPSSHRINKVR(%)  4.1415E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3681.9123014814322     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2027.8229417130215     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.33
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3681.912       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.58E-01  1.10E+00  1.00E+00  9.59E-01  1.17E+00  1.09E+00  1.32E+00  3.46E+00  9.43E-01  9.83E-01  1.02E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       78.989
Stop Time:
Sat Oct 23 13:20:42 CDT 2021
