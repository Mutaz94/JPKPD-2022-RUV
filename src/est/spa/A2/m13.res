Sat Sep 18 09:39:44 CDT 2021
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
$DATA ../../../../data/spa/A2/dat13.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m13.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1032.91641400534        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9117E+01 -4.6022E+01  6.2656E+01 -1.2237E+02  3.3821E+00  1.4109E+00 -6.9465E+00 -4.1122E+01 -3.4762E-01 -2.1972E+00
            -1.1742E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1424.05119884192        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0387E+00  1.0888E+00  9.8278E-01  1.0508E+00  1.0701E+00  8.8076E-01  8.7999E-01  1.0835E+00  7.0104E-01  7.4553E-01
             3.1209E+00
 PARAMETER:  1.3797E-01  1.8510E-01  8.2626E-02  1.4960E-01  1.6772E-01 -2.6975E-02 -2.7845E-02  1.8018E-01 -2.5519E-01 -1.9365E-01
             1.2381E+00
 GRADIENT:   6.7678E+01 -1.8209E+01 -1.3529E+01 -1.0646E+01  7.0971E+00 -2.3570E+01  2.4934E+00  2.7885E+00  4.7119E+00  1.1911E+01
             9.6343E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1431.96748118743        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0411E+00  1.2274E+00  8.3322E-01  9.6147E-01  1.0029E+00  9.6365E-01  7.7221E-01  1.4921E+00  8.9698E-01  3.4102E-01
             2.9337E+00
 PARAMETER:  1.4031E-01  3.0492E-01 -8.2457E-02  6.0707E-02  1.0288E-01  6.2975E-02 -1.5850E-01  5.0022E-01 -8.7248E-03 -9.7582E-01
             1.1763E+00
 GRADIENT:   6.2767E+01 -6.7389E-01 -2.0840E+00  4.7718E+00 -1.4165E+01  2.9244E+00  3.3051E+00  6.3329E+00  1.1974E+01  2.5625E+00
             8.4307E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1443.04805929321        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.0053E+00  1.0795E+00  3.4183E-01  9.7184E-01  5.8414E-01  9.7940E-01  9.9975E-01  6.2108E-01  7.1205E-01  1.7906E-01
             2.3213E+00
 PARAMETER:  1.0531E-01  1.7650E-01 -9.7344E-01  7.1437E-02 -4.3762E-01  7.9185E-02  9.9750E-02 -3.7630E-01 -2.3961E-01 -1.6201E+00
             9.4213E-01
 GRADIENT:  -3.6979E+00  4.4539E+01  1.9031E+00  6.2311E+01 -3.5236E+00  1.1928E+00  2.7166E+00 -2.7551E+00  6.4107E-02  1.0906E+00
            -2.1877E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1447.42812398326        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  9.8999E-01  1.0497E+00  2.0779E-01  8.6697E-01  4.8350E-01  9.7723E-01  8.9363E-01  4.1185E-01  7.5701E-01  6.5501E-02
             2.2937E+00
 PARAMETER:  8.9935E-02  1.4849E-01 -1.4712E+00 -4.2749E-02 -6.2670E-01  7.6964E-02 -1.2466E-02 -7.8710E-01 -1.7838E-01 -2.6257E+00
             9.3016E-01
 GRADIENT:  -6.4462E+00 -2.3269E+00  4.2361E+00  1.9956E-02 -1.2360E+00  8.9365E-01  8.3317E-01 -3.0876E+00  1.7688E+00  2.1119E-01
            -8.7296E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1470.09258750883        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      388
 NPARAMETR:  9.8402E-01  1.4157E+00  1.7304E-01  6.7937E-01  6.4671E-01  9.4027E-01  7.1566E-01  1.8473E+00  8.7750E-01  1.9600E-02
             2.3150E+00
 PARAMETER:  8.3894E-02  4.4765E-01 -1.6542E+00 -2.8660E-01 -3.3585E-01  3.8410E-02 -2.3456E-01  7.1374E-01 -3.0683E-02 -3.8322E+00
             9.3940E-01
 GRADIENT:  -1.2317E+01 -5.4725E+01 -1.8065E+00  3.1304E+01  4.9643E+01 -6.4049E+00 -5.7612E+00 -3.3123E+00  7.9974E+00  2.5238E-02
             7.5599E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1481.93802980400        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  9.7543E-01  1.6732E+00  1.2873E-01  5.3047E-01  7.5504E-01  9.5467E-01  7.2293E-01  2.1097E+00  9.5811E-01  1.0000E-02
             2.0014E+00
 PARAMETER:  7.5127E-02  6.1474E-01 -1.9500E+00 -5.3399E-01 -1.8098E-01  5.3612E-02 -2.2444E-01  8.4655E-01  5.7203E-02 -4.9323E+00
             7.9384E-01
 GRADIENT:  -1.9196E+01 -4.5098E+00 -7.4498E+00  1.8327E+01  2.7264E+01 -2.2716E+00 -4.7587E-01 -1.2345E+00 -8.3675E+00  0.0000E+00
             1.0845E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1482.03337008318        NO. OF FUNC. EVALS.: 111
 CUMULATIVE NO. OF FUNC. EVALS.:      572
 NPARAMETR:  9.7574E-01  1.6750E+00  1.2830E-01  5.2938E-01  7.5549E-01  9.5498E-01  7.2250E-01  2.1091E+00  9.6555E-01  1.0000E-02
             1.9982E+00
 PARAMETER:  7.5442E-02  6.1579E-01 -1.9534E+00 -5.3604E-01 -1.8040E-01  5.3935E-02 -2.2504E-01  8.4626E-01  6.4946E-02 -4.9461E+00
             7.9222E-01
 GRADIENT:  -2.7250E+01 -1.6547E+01 -8.4450E+00  1.5915E+01  2.5765E+01 -2.9382E+00 -6.4620E-01 -1.5824E+00 -8.0896E+00  0.0000E+00
             9.8632E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1482.08707547790        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      734
 NPARAMETR:  9.7593E-01  1.6747E+00  1.2836E-01  5.2915E-01  7.5533E-01  9.6247E-01  7.2257E-01  2.1095E+00  9.6804E-01  1.0000E-02
             1.9975E+00
 PARAMETER:  7.5632E-02  6.1563E-01 -1.9529E+00 -5.3648E-01 -1.8060E-01  6.1751E-02 -2.2494E-01  8.4644E-01  6.7516E-02 -4.9461E+00
             7.9191E-01
 GRADIENT:  -1.7268E+01 -4.6334E+00 -7.1219E+00  1.7227E+01  2.5675E+01  7.8907E-01 -3.8190E-01 -1.1189E+00 -7.8747E+00  0.0000E+00
             1.0672E+01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1482.08707547790        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      813
 NPARAMETR:  9.7592E-01  1.6747E+00  1.2836E-01  5.2914E-01  7.5534E-01  9.6247E-01  7.2273E-01  2.1094E+00  9.6804E-01  1.0000E-02
             1.9976E+00
 PARAMETER:  7.5632E-02  6.1563E-01 -1.9529E+00 -5.3648E-01 -1.8060E-01  6.1751E-02 -2.2494E-01  8.4644E-01  6.7516E-02 -4.9461E+00
             7.9191E-01
 GRADIENT:   1.5511E+05  2.5176E+04  7.9463E+03  1.4468E+04 -8.5875E+04  1.2194E-02 -5.9158E-01  1.8297E+04 -7.7573E+04  0.0000E+00
            -1.9584E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      813
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1695E-02 -7.1456E-03  4.1469E-03  2.8501E-03 -5.0357E-04
 SE:             2.9418E-02  2.6049E-02  1.7139E-02  2.0709E-02  2.5177E-04
 N:                     100         100         100         100         100

 P VAL.:         6.9096E-01  7.8384E-01  8.0881E-01  8.9054E-01  4.5487E-02

 ETASHRINKSD(%)  1.4472E+00  1.2733E+01  4.2583E+01  3.0621E+01  9.9157E+01
 ETASHRINKVR(%)  2.8735E+00  2.3845E+01  6.7033E+01  5.1866E+01  9.9993E+01
 EBVSHRINKSD(%)  1.5585E+00  1.2983E+01  4.2456E+01  3.3848E+01  9.9038E+01
 EBVSHRINKVR(%)  3.0927E+00  2.4281E+01  6.6886E+01  5.6240E+01  9.9991E+01
 RELATIVEINF(%)  9.5233E+01  9.3524E+00  1.4726E+01  6.2859E+00  1.1527E-03
 EPSSHRINKSD(%)  3.8173E+01
 EPSSHRINKVR(%)  6.1775E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1482.0870754779046     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -746.93624891416641     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.98
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1482.087       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  1.67E+00  1.28E-01  5.29E-01  7.55E-01  9.62E-01  7.23E-01  2.11E+00  9.68E-01  1.00E-02  2.00E+00
 


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
+        4.07E+08
 
 TH 2
+        3.65E+02  3.65E+06
 
 TH 3
+        1.52E+03  4.11E+04  6.19E+07
 
 TH 4
+       -9.71E+04 -8.91E+03 -3.96E+04  4.81E+07
 
 TH 5
+       -3.05E+03 -5.38E+02 -1.14E+08  7.02E+04  2.08E+08
 
 TH 6
+        3.76E+03  3.46E+02  1.46E+03 -1.42E+08 -2.68E+03  2.01E+02
 
 TH 7
+        1.65E+04  1.58E+03  6.36E+03 -8.40E+07 -1.18E+04 -9.93E-01  1.47E+08
 
 TH 8
+        2.30E+02 -6.80E+03  2.34E+04 -5.30E+03  3.64E+01  2.06E+02  8.99E+02  1.21E+06
 
 TH 9
+        6.92E+03  6.40E+02  2.74E+03 -1.41E+08 -4.96E+03 -3.77E+03 -1.66E+04  3.82E+02  4.14E+08
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.71E+02 -1.21E+03 -2.65E+04  5.99E+03 -9.96E+01 -2.28E+02 -9.99E+02 -6.84E+02 -4.09E+02  0.00E+00  1.55E+06
 
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
 #CPUT: Total CPU Time in Seconds,       15.577
Stop Time:
Sat Sep 18 09:40:01 CDT 2021
