Sat Sep 18 08:44:50 CDT 2021
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
$DATA ../../../../data/spa/B/dat80.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1716.86235666759        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.2645E+01 -2.6931E+01  9.6198E+00 -4.8295E+01 -3.7800E+01  3.7088E+01  1.2508E+01  6.7117E+00  3.4051E+01  1.9666E+01
            -3.2703E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1725.12667830510        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0387E+00  1.0624E+00  1.0361E+00  9.9799E-01  1.0598E+00  8.5014E-01  9.1128E-01  9.4975E-01  7.9464E-01  8.7993E-01
             1.1429E+00
 PARAMETER:  1.3797E-01  1.6053E-01  1.3550E-01  9.7984E-02  1.5808E-01 -6.2349E-02  7.0982E-03  4.8444E-02 -1.2987E-01 -2.7917E-02
             2.3357E-01
 GRADIENT:   7.9653E+01  5.7552E+00  1.1661E+01 -8.6636E+00 -6.7909E-01 -2.4505E+01 -2.6280E+00 -1.7573E+00 -1.1022E+01 -1.7692E+00
             1.2413E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1726.10934936142        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0312E+00  1.0490E+00  7.8431E-01  9.8794E-01  9.0462E-01  8.7853E-01  9.9328E-01  7.8929E-01  7.4378E-01  6.6136E-01
             1.1154E+00
 PARAMETER:  1.3077E-01  1.4784E-01 -1.4295E-01  8.7862E-02 -2.3888E-04 -2.9504E-02  9.3254E-02 -1.3662E-01 -1.9601E-01 -3.1345E-01
             2.0925E-01
 GRADIENT:   5.0324E+01  3.8842E+00  7.8752E+00 -6.8769E+00 -4.0506E+00 -1.1118E+01 -3.9463E+00  5.7571E-01 -1.5171E+01 -6.2545E+00
             4.6726E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1727.91539111154        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.0178E+00  1.1453E+00  6.1176E-01  9.1591E-01  8.4841E-01  9.0078E-01  9.3345E-01  4.3364E-01  8.3291E-01  6.8283E-01
             1.0899E+00
 PARAMETER:  1.1761E-01  2.3570E-01 -3.9142E-01  1.2162E-02 -6.4397E-02 -4.4932E-03  3.1128E-02 -7.3555E-01 -8.2828E-02 -2.8151E-01
             1.8605E-01
 GRADIENT:   7.2640E+00  8.0696E+00 -1.2133E+00  5.3334E+00 -4.8726E+00 -2.1781E+00 -5.6581E-02  1.2900E+00 -1.0754E+00  1.4520E+00
             2.1860E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1727.97986779945        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0142E+00  1.1022E+00  6.4396E-01  9.4136E-01  8.5319E-01  9.0592E-01  9.6328E-01  3.6299E-01  8.3013E-01  7.0781E-01
             1.0858E+00
 PARAMETER:  1.1410E-01  1.9728E-01 -3.4012E-01  3.9571E-02 -5.8778E-02  1.1929E-03  6.2591E-02 -9.1337E-01 -8.6175E-02 -2.4558E-01
             1.8235E-01
 GRADIENT:  -1.6508E+00 -1.3470E+00 -2.5665E+00  7.2765E-01  2.4137E+00  1.2411E-01  4.6232E-01  6.6847E-01  5.1146E-01  6.9845E-01
             7.2128E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1727.99416672957        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.0140E+00  1.0820E+00  6.3442E-01  9.5110E-01  8.3719E-01  9.0660E-01  9.8239E-01  2.7288E-01  8.2149E-01  7.0349E-01
             1.0843E+00
 PARAMETER:  1.1392E-01  1.7886E-01 -3.5504E-01  4.9863E-02 -7.7699E-02  1.9432E-03  8.2230E-02 -1.1987E+00 -9.6641E-02 -2.5170E-01
             1.8094E-01
 GRADIENT:  -2.2337E+00 -2.9434E+00 -2.5399E+00  2.0005E-01  2.9202E+00  3.9908E-01  7.2800E-01  3.7147E-01  7.3045E-01  7.9386E-01
             4.4320E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1728.08410811551        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  1.0151E+00  1.0945E+00  5.9117E-01  9.3740E-01  8.1202E-01  9.0613E-01  9.7528E-01  7.9073E-02  8.1772E-01  6.7060E-01
             1.0832E+00
 PARAMETER:  1.1502E-01  1.9029E-01 -4.2564E-01  3.5350E-02 -1.0823E-01  1.4306E-03  7.4970E-02 -2.4374E+00 -1.0124E-01 -2.9958E-01
             1.7993E-01
 GRADIENT:   1.0469E-01 -8.4431E-01  3.7011E-01 -9.6524E-01 -1.1725E+00  3.6753E-02 -2.3269E-01  3.4457E-02 -3.1814E-02  1.1950E-01
            -5.1483E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1728.71577322632        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      588
 NPARAMETR:  1.0298E+00  1.1864E+00  5.9574E-01  8.9039E-01  8.6395E-01  9.1389E-01  9.0444E-01  1.2834E-02  8.5940E-01  7.0749E-01
             1.0875E+00
 PARAMETER:  1.2933E-01  2.7092E-01 -4.1795E-01 -1.6095E-02 -4.6241E-02  9.9567E-03 -4.3621E-04 -4.2556E+00 -5.1516E-02 -2.4604E-01
             1.8387E-01
 GRADIENT:  -3.7911E+00 -2.7896E-01  1.3397E+00  5.8471E-01 -2.1456E+00  1.2660E-01 -1.5478E+00  4.2366E-04 -9.5271E-01 -6.7312E-01
            -1.2347E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1728.85438514455        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      766
 NPARAMETR:  1.0306E+00  1.3266E+00  5.5319E-01  8.0373E-01  9.1596E-01  9.1430E-01  8.3464E-01  1.0000E-02  9.3206E-01  7.2233E-01
             1.0890E+00
 PARAMETER:  1.3014E-01  3.8259E-01 -4.9206E-01 -1.1849E-01  1.2219E-02  1.0401E-02 -8.0750E-02 -5.5651E+00  2.9643E-02 -2.2527E-01
             1.8526E-01
 GRADIENT:  -2.2073E+00  1.1450E+00  9.6118E-02  1.1185E+00 -2.6367E-01  7.6025E-02  1.3563E-01  0.0000E+00  1.3371E-01 -3.4777E-02
            -2.6489E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1728.85781403744        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      928
 NPARAMETR:  1.0315E+00  1.3462E+00  5.4450E-01  7.9045E-01  9.2203E-01  9.1417E-01  8.2407E-01  1.0000E-02  9.4202E-01  7.2269E-01
             1.0893E+00
 PARAMETER:  1.3101E-01  3.9730E-01 -5.0789E-01 -1.3516E-01  1.8822E-02  1.0262E-02 -9.3499E-02 -5.7374E+00  4.0267E-02 -2.2478E-01
             1.8554E-01
 GRADIENT:   1.7864E-02 -4.7032E-03  8.2938E-03 -1.5072E-02 -1.4108E-02 -1.3348E-03 -1.6553E-03  0.0000E+00  2.1285E-03 -1.4442E-03
            -1.1683E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      928
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.4361E-05 -1.3531E-02 -3.4154E-04  8.2439E-03 -2.1422E-02
 SE:             2.9809E-02  2.3252E-02  1.5666E-04  2.3828E-02  2.0782E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9908E-01  5.6063E-01  2.9242E-02  7.2937E-01  3.0264E-01

 ETASHRINKSD(%)  1.3754E-01  2.2101E+01  9.9475E+01  2.0172E+01  3.0377E+01
 ETASHRINKVR(%)  2.7490E-01  3.9318E+01  9.9997E+01  3.6275E+01  5.1526E+01
 EBVSHRINKSD(%)  5.7363E-01  2.1988E+01  9.9511E+01  2.0616E+01  3.0378E+01
 EBVSHRINKVR(%)  1.1440E+00  3.9141E+01  9.9998E+01  3.6981E+01  5.1528E+01
 RELATIVEINF(%)  9.8787E+01  2.6835E+00  1.5651E-04  3.0310E+00  4.1730E+00
 EPSSHRINKSD(%)  4.2729E+01
 EPSSHRINKVR(%)  6.7200E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1728.8578140374389     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -993.70698747370068     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.71
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1728.858       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.35E+00  5.44E-01  7.90E-01  9.22E-01  9.14E-01  8.24E-01  1.00E-02  9.42E-01  7.23E-01  1.09E+00
 


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
+       -1.04E+01  5.65E+02
 
 TH 3
+        1.32E+01  3.41E+02  9.52E+02
 
 TH 4
+       -2.51E+01  3.95E+02 -5.60E+02  1.28E+03
 
 TH 5
+       -6.24E+00 -4.56E+02 -9.38E+02  4.81E+02  1.21E+03
 
 TH 6
+       -5.33E-01 -1.69E+00  3.28E+00 -3.87E+00 -2.20E+00  2.36E+02
 
 TH 7
+       -5.62E-01  2.51E+01 -2.87E+01 -1.17E+01 -7.28E-01 -6.91E-02  1.20E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.97E+00 -2.44E+01 -4.50E+01  5.72E+01 -6.92E+00 -1.33E-01  2.06E+01  0.00E+00  1.02E+02
 
 TH10
+       -8.23E-01 -1.25E+01 -5.84E+01 -2.40E+01 -6.70E+01 -5.79E-01  3.09E+01  0.00E+00  1.41E+01  9.59E+01
 
 TH11
+       -8.29E+00 -1.95E+01 -3.39E+01 -4.32E+00 -7.31E+00  2.75E+00  1.00E+01  0.00E+00  1.36E+01  2.52E+01  1.83E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.170
Stop Time:
Sat Sep 18 08:45:06 CDT 2021
