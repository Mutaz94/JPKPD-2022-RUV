Sat Oct 23 21:51:45 CDT 2021
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
$DATA ../../../../data/SD3/A1/dat61.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m61.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1733.94255740101        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8999E+02  1.7789E+01  3.9424E+01  3.4507E+01  8.5527E+01  6.2483E+01 -4.6012E+01 -2.4284E+01  9.4535E+00 -3.5073E+01
            -6.8064E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1898.03340002543        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0330E+00  1.0071E+00  9.7109E-01  1.0473E+00  9.3075E-01  1.0340E+00  1.2121E+00  7.7362E-01  9.0258E-01  9.1716E-01
             1.8270E+00
 PARAMETER:  1.3249E-01  1.0705E-01  7.0662E-02  1.4618E-01  2.8235E-02  1.3340E-01  2.9235E-01 -1.5667E-01 -2.4937E-03  1.3528E-02
             7.0265E-01
 GRADIENT:   1.8504E+02  3.1693E+01  1.5524E+01  3.4542E+01 -1.0758E+01  4.3903E+01  2.6876E+00  5.2499E+00  2.5632E-01 -5.1554E+00
             6.1831E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1907.74033170987        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0238E+00  7.8628E-01  5.0185E-01  1.1481E+00  5.8979E-01  9.5726E-01  1.5899E+00  2.2850E-01  8.1119E-01  7.8557E-01
             1.7345E+00
 PARAMETER:  1.2357E-01 -1.4044E-01 -5.8946E-01  2.3811E-01 -4.2799E-01  5.6319E-02  5.6364E-01 -1.3762E+00 -1.0925E-01 -1.4134E-01
             6.5070E-01
 GRADIENT:   1.4951E+02  3.3932E+01 -4.9965E+01  1.8704E+02  9.3103E+01  8.5624E+00  2.0660E+01  1.7066E+00 -4.3808E+00  2.4611E+01
             6.4146E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1911.56868238358        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      341
 NPARAMETR:  1.0253E+00  6.8295E-01  6.2379E-01  1.2377E+00  6.0785E-01  9.7576E-01  1.7566E+00  3.5872E-01  8.3117E-01  7.5240E-01
             1.7052E+00
 PARAMETER:  1.2503E-01 -2.8134E-01 -3.7195E-01  3.1325E-01 -3.9783E-01  7.5463E-02  6.6336E-01 -9.2522E-01 -8.4923E-02 -1.8449E-01
             6.3371E-01
 GRADIENT:  -4.5873E+00  4.3853E+01 -5.9413E+00  9.6386E+01 -4.9452E+00  2.2211E+00  7.1311E+00  3.0322E+00 -3.8311E-01  9.1700E+00
             3.5236E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1921.09267230696        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      519
 NPARAMETR:  1.0289E+00  5.7391E-01  4.1266E-01  1.1641E+00  4.5120E-01  9.7427E-01  1.6402E+00  1.4198E-01  8.4552E-01  5.3429E-01
             1.6112E+00
 PARAMETER:  1.2851E-01 -4.5529E-01 -7.8514E-01  2.5197E-01 -6.9586E-01  7.3929E-02  5.9484E-01 -1.8521E+00 -6.7799E-02 -5.2682E-01
             5.7696E-01
 GRADIENT:   1.3915E+00  1.6305E+00  1.1410E+00 -1.4918E+00 -2.1744E+00  2.4932E-01  8.0216E-02  2.9355E-01 -3.0530E-02  2.5580E-01
             1.6197E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1921.20785146503        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      696
 NPARAMETR:  1.0299E+00  5.6434E-01  4.1729E-01  1.1703E+00  4.5263E-01  9.7198E-01  1.6629E+00  6.2223E-02  8.4352E-01  5.4320E-01
             1.6101E+00
 PARAMETER:  1.2947E-01 -4.7209E-01 -7.7398E-01  2.5723E-01 -6.9269E-01  7.1583E-02  6.0856E-01 -2.6770E+00 -7.0177E-02 -5.1028E-01
             5.7628E-01
 GRADIENT:   4.2232E+00 -5.5775E-01 -1.1329E+00  3.5827E-01  1.7699E+00 -5.2741E-01 -1.0631E-01  5.4546E-02  5.1377E-02  3.1857E-01
            -1.0435E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1921.23859525443        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      858
 NPARAMETR:  1.0282E+00  5.6427E-01  4.1676E-01  1.1700E+00  4.5209E-01  9.7303E-01  1.6648E+00  1.0000E-02  8.4327E-01  5.4180E-01
             1.6125E+00
 PARAMETER:  1.2783E-01 -4.7222E-01 -7.7525E-01  2.5700E-01 -6.9388E-01  7.2657E-02  6.0972E-01 -4.6089E+00 -7.0470E-02 -5.1286E-01
             5.7779E-01
 GRADIENT:   3.6041E-01 -9.7765E-02 -1.9130E-01 -7.1458E-02  2.5574E-01 -1.0894E-01  2.6868E-02  0.0000E+00  3.6025E-02  2.9955E-02
             9.8679E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      858
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1911E-03  2.4186E-02 -4.4922E-04 -1.4881E-02  1.3695E-02
 SE:             2.9719E-02  2.1438E-02  2.2409E-04  2.6966E-02  1.8096E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6803E-01  2.5925E-01  4.4995E-02  5.8107E-01  4.4916E-01

 ETASHRINKSD(%)  4.3705E-01  2.8179E+01  9.9249E+01  9.6597E+00  3.9378E+01
 ETASHRINKVR(%)  8.7218E-01  4.8418E+01  9.9994E+01  1.8386E+01  6.3249E+01
 EBVSHRINKSD(%)  8.7010E-01  2.8328E+01  9.9212E+01  9.3396E+00  3.8649E+01
 EBVSHRINKVR(%)  1.7326E+00  4.8632E+01  9.9994E+01  1.7807E+01  6.2361E+01
 RELATIVEINF(%)  9.8084E+01  1.0824E+01  6.3327E-04  3.7413E+01  2.6139E+00
 EPSSHRINKSD(%)  3.1211E+01
 EPSSHRINKVR(%)  5.2680E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1921.2385952544325     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1002.3000620497598     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.57
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1921.239       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  5.64E-01  4.17E-01  1.17E+00  4.52E-01  9.73E-01  1.66E+00  1.00E-02  8.43E-01  5.42E-01  1.61E+00
 


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
 #CPUT: Total CPU Time in Seconds,       65.567
Stop Time:
Sat Oct 23 21:51:56 CDT 2021
