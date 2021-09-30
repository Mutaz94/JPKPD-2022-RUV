Wed Sep 29 14:09:41 CDT 2021
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
$DATA ../../../../data/spa/S1/dat28.csv ignore=@
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
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       29 SEP 2021
Days until program expires : 200
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
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1648.12027112526        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1345E+02 -8.3157E+01 -5.0989E+01 -6.1940E+01  9.2803E+01  4.0456E+01  2.7590E+00  4.9829E+00 -1.1605E+01  2.1682E+01
            -8.7425E-01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1660.30713111149        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.8952E-01  1.0331E+00  1.1099E+00  1.0735E+00  9.8039E-01  9.4258E-01  9.4453E-01  9.8359E-01  1.0772E+00  8.1211E-01
             1.0168E+00
 PARAMETER:  8.9463E-02  1.3254E-01  2.0430E-01  1.7090E-01  8.0196E-02  4.0869E-02  4.2936E-02  8.3451E-02  1.7438E-01 -1.0812E-01
             1.1665E-01
 GRADIENT:   1.7800E+01  2.5396E+00  2.4545E+00  5.0253E+00  2.3951E+01 -1.2247E+01  5.5967E+00 -6.5688E+00  3.6520E+00 -6.6204E+00
            -9.5345E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1661.39202733518        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.8978E-01  9.3673E-01  1.1659E+00  1.1469E+00  9.4891E-01  9.6474E-01  8.0563E-01  1.0628E+00  1.0526E+00  7.9093E-01
             1.0152E+00
 PARAMETER:  8.9732E-02  3.4641E-02  2.5347E-01  2.3707E-01  4.7554E-02  6.4103E-02 -1.1613E-01  1.6094E-01  1.5130E-01 -1.3454E-01
             1.1506E-01
 GRADIENT:   1.9244E+01  1.6758E+01  7.9100E+00  2.2501E+01 -4.9942E-01 -2.4710E+00  3.9224E-01 -5.1390E+00  2.4848E+00 -8.0086E+00
            -1.8761E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1662.28026529479        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  9.7820E-01  8.2341E-01  1.3015E+00  1.2113E+00  9.5698E-01  9.7397E-01  6.6111E-01  1.2232E+00  1.0090E+00  8.6783E-01
             1.0116E+00
 PARAMETER:  7.7955E-02 -9.4299E-02  3.6350E-01  2.9170E-01  5.6031E-02  7.3630E-02 -3.1384E-01  3.0146E-01  1.0894E-01 -4.1754E-02
             1.1150E-01
 GRADIENT:  -5.1755E+00  9.5474E+00  3.9556E+00  8.7087E+00 -1.2460E+01  2.0054E+00  2.1547E-01  1.0210E+00 -2.0599E+00  2.4811E+00
            -4.4314E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1662.76607831365        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  9.7883E-01  5.5526E-01  1.3902E+00  1.3806E+00  9.0948E-01  9.6456E-01  4.4080E-01  1.2031E+00  9.1253E-01  8.3840E-01
             1.0136E+00
 PARAMETER:  7.8607E-02 -4.8831E-01  4.2946E-01  4.2252E-01  5.1216E-03  6.3920E-02 -7.1916E-01  2.8491E-01  8.4669E-03 -7.6256E-02
             1.1346E-01
 GRADIENT:   2.7320E+00  3.6808E+00  1.2579E+00  7.4593E+00 -1.5943E-01 -9.7802E-01 -9.8807E-02 -6.9794E-01 -6.9624E-01 -7.0795E-01
            -3.9628E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1662.83260893152        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      877
 NPARAMETR:  9.7617E-01  4.2082E-01  1.3984E+00  1.4659E+00  8.7297E-01  9.6545E-01  2.9793E-01  1.1987E+00  8.6058E-01  8.2017E-01
             1.0138E+00
 PARAMETER:  7.5883E-02 -7.6554E-01  4.3529E-01  4.8249E-01 -3.5856E-02  6.4840E-02 -1.1109E+00  2.8123E-01 -5.0154E-02 -9.8246E-02
             1.1375E-01
 GRADIENT:   1.7711E-01  3.6599E+00  1.7789E+00  1.0300E+01 -2.5974E+00 -2.3279E-01 -3.6709E-02 -5.4410E-01 -9.7087E-01 -1.2719E-01
            -2.8847E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1662.86380280519        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1052
 NPARAMETR:  9.7417E-01  3.1247E-01  1.3899E+00  1.5312E+00  8.4079E-01  9.6575E-01  1.7617E-01  1.1993E+00  8.2063E-01  7.9980E-01
             1.0142E+00
 PARAMETER:  7.3831E-02 -1.0633E+00  4.2926E-01  5.2604E-01 -7.3418E-02  6.5152E-02 -1.6363E+00  2.8170E-01 -9.7682E-02 -1.2340E-01
             1.1409E-01
 GRADIENT:  -1.0393E+00  2.5699E+00  1.4362E+00  7.8277E+00 -3.3901E+00  1.9925E-01 -5.0970E-03 -1.7586E-01 -8.6139E-01  2.7777E-01
            -1.1796E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1662.90776620119        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1238             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7457E-01  2.8904E-01  1.3895E+00  1.5384E+00  8.3550E-01  9.6530E-01  1.6251E-01  1.2022E+00  8.1390E-01  7.9518E-01
             1.0142E+00
 PARAMETER:  7.4239E-02 -1.1412E+00  4.2894E-01  5.3073E-01 -7.9727E-02  6.4679E-02 -1.7170E+00  2.8412E-01 -1.0592E-01 -1.2918E-01
             1.1415E-01
 GRADIENT:   3.5846E+02  2.9891E+01  9.5263E+00  6.2789E+02  6.7212E+00  2.9765E+01  1.0439E-01  1.0029E+00  8.8514E+00  8.6353E-01
             9.2270E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1662.90836865617        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:     1370
 NPARAMETR:  9.7463E-01  2.8866E-01  1.3887E+00  1.5381E+00  8.3563E-01  9.6523E-01  1.6532E-01  1.2026E+00  8.1391E-01  7.9427E-01
             1.0142E+00
 PARAMETER:  7.4305E-02 -1.1425E+00  4.2834E-01  5.3056E-01 -7.9567E-02  6.4610E-02 -1.6999E+00  2.8445E-01 -1.0591E-01 -1.3033E-01
             1.1414E-01
 GRADIENT:   1.1086E+00  1.0758E-01  8.7271E-02 -8.3261E+00  6.1933E-01  8.4543E-02 -1.3512E-03 -1.0813E-02  7.4003E-02  2.4951E-02
            -5.9522E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1370
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.8016E-04 -1.6623E-03 -2.8059E-02 -4.1846E-03 -3.2995E-02
 SE:             2.9843E-02  9.0814E-04  1.8691E-02  2.9277E-02  1.9958E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8984E-01  6.7185E-02  1.3330E-01  8.8634E-01  9.8280E-02

 ETASHRINKSD(%)  2.1367E-02  9.6958E+01  3.7384E+01  1.9178E+00  3.3139E+01
 ETASHRINKVR(%)  4.2729E-02  9.9907E+01  6.0793E+01  3.7987E+00  5.5296E+01
 EBVSHRINKSD(%)  4.3842E-01  9.7258E+01  3.9454E+01  2.4085E+00  3.1485E+01
 EBVSHRINKVR(%)  8.7491E-01  9.9925E+01  6.3342E+01  4.7590E+00  5.3057E+01
 RELATIVEINF(%)  9.7724E+01  4.7059E-03  7.0935E+00  8.0469E+00  5.4137E+00
 EPSSHRINKSD(%)  4.4741E+01
 EPSSHRINKVR(%)  6.9465E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1662.9083686561653     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -927.75754209242712     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.03
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.47
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1662.908       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  2.89E-01  1.39E+00  1.54E+00  8.36E-01  9.65E-01  1.65E-01  1.20E+00  8.14E-01  7.94E-01  1.01E+00
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.24E+03
 
 TH 2
+       -2.67E+01  3.86E+02
 
 TH 3
+        1.47E+00  7.88E+01  1.52E+02
 
 TH 4
+       -1.01E+01  4.30E+02 -1.57E+01  6.75E+02
 
 TH 5
+        2.29E+00 -3.43E+02 -3.58E+02 -7.43E+01  1.21E+03
 
 TH 6
+        2.36E-01 -2.23E+00  5.79E-01 -2.16E+00 -9.44E-01  2.10E+02
 
 TH 7
+        3.06E-02 -8.98E-01 -1.05E-02 -4.53E-01  3.15E-01 -1.83E-02 -9.97E-02
 
 TH 8
+        2.22E-01  1.72E+00 -3.17E+01 -2.79E+00 -1.16E+01  1.69E-01  4.44E-02  3.70E+01
 
 TH 9
+        2.27E+00 -1.00E+02  4.95E+00 -1.74E-02  2.76E+00 -1.00E+00  1.07E+00 -1.54E-01  2.74E+02
 
 TH10
+        3.83E-01  9.35E+00 -3.89E-01 -5.34E-01 -1.04E+02  2.74E-02  2.02E-01  1.84E+01  1.62E+00  9.43E+01
 
 TH11
+       -7.78E+00 -7.07E+00 -4.42E+00 -7.14E+00 -1.13E+01  2.42E+00  8.98E-02  6.94E+00  1.07E+01  1.57E+01  2.05E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       22.528
Stop Time:
Wed Sep 29 14:10:05 CDT 2021
