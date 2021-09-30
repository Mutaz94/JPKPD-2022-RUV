Wed Sep 29 12:23:30 CDT 2021
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
$DATA ../../../../data/spa/A1/dat82.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m82.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1274.46768226971        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2740E+02 -6.8686E+00 -5.6361E+01  9.0891E+01  1.8033E+02  6.7085E+01 -3.3442E+01  8.4007E-01 -3.4306E+01 -5.5562E+01
            -6.4507E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1465.05680983259        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0252E+00  1.0200E+00  1.0777E+00  1.0007E+00  9.3533E-01  8.7177E-01  1.1733E+00  9.3733E-01  1.1278E+00  1.0499E+00
             1.9248E+00
 PARAMETER:  1.2486E-01  1.1978E-01  1.7483E-01  1.0068E-01  3.3143E-02 -3.7226E-02  2.5980E-01  3.5280E-02  2.2024E-01  1.4868E-01
             7.5483E-01
 GRADIENT:   2.0769E+02  7.1218E+00 -7.8586E+00  2.0536E+01  1.1033E+01 -2.5704E+00  8.6711E-01  4.4084E+00  1.3143E+01 -4.1012E+00
            -5.4055E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1467.34338617296        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.0057E+00  9.0979E-01  1.1658E+00  1.0751E+00  9.2672E-01  8.8405E-01  1.1461E+00  4.0898E-01  1.1200E+00  1.1105E+00
             1.9221E+00
 PARAMETER:  1.0570E-01  5.4615E-03  2.5344E-01  1.7244E-01  2.3898E-02 -2.3238E-02  2.3637E-01 -7.9408E-01  2.1337E-01  2.0483E-01
             7.5341E-01
 GRADIENT:   1.3425E+02  7.3933E+00 -1.4876E+00  4.1194E+01  8.0161E+00  5.1205E+00 -1.9719E+00  3.4678E-01  1.5175E+01 -1.5334E+00
            -7.0428E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1468.32676920460        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      268
 NPARAMETR:  1.0006E+00  8.9270E-01  1.1238E+00  1.0830E+00  9.0401E-01  8.9795E-01  1.3048E+00  4.1279E-01  1.0144E+00  1.1016E+00
             1.9509E+00
 PARAMETER:  1.0058E-01 -1.3506E-02  2.1676E-01  1.7977E-01 -9.1523E-04 -7.6436E-03  3.6609E-01 -7.8482E-01  1.1430E-01  1.9679E-01
             7.6827E-01
 GRADIENT:  -2.1783E+00  2.3812E+00 -1.2325E+00  1.0203E+00  3.6656E-01  2.1883E+00 -2.4520E+00  4.6397E-01 -1.0457E+00 -7.5047E-01
            -3.9545E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1469.56688843672        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  9.9579E-01  5.3105E-01  1.2033E+00  1.3187E+00  8.0959E-01  8.9343E-01  1.7826E+00  5.4594E-02  9.0729E-01  1.0906E+00
             1.9576E+00
 PARAMETER:  9.5778E-02 -5.3290E-01  2.8503E-01  3.7663E-01 -1.1123E-01 -1.2684E-02  6.7807E-01 -2.8078E+00  2.7061E-03  1.8672E-01
             7.7173E-01
 GRADIENT:  -7.2845E+00  7.0799E+00  2.3046E+00  1.8231E+01 -7.2789E+00  1.5384E+00 -7.7309E-01  4.4961E-03 -8.8213E-01 -4.7532E-01
            -4.6784E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1470.20547929470        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      621
 NPARAMETR:  9.9584E-01  2.9475E-01  1.3826E+00  1.4715E+00  8.0849E-01  8.8271E-01  2.3695E+00  1.0000E-02  8.6213E-01  1.1341E+00
             1.9925E+00
 PARAMETER:  9.5835E-02 -1.1216E+00  4.2398E-01  4.8625E-01 -1.1259E-01 -2.4756E-02  9.6270E-01 -5.9383E+00 -4.8352E-02  2.2584E-01
             7.8937E-01
 GRADIENT:   1.2954E+00  3.8478E+00  1.6524E+00  1.7507E+01 -4.4475E+00 -1.1644E+00  6.0503E-01  0.0000E+00  2.2395E-01  3.0486E-01
             1.6209E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1470.63480064547        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      800
 NPARAMETR:  9.9294E-01  1.0832E-01  1.3942E+00  1.5793E+00  7.6793E-01  8.8304E-01  3.5828E+00  1.0000E-02  8.2364E-01  1.1396E+00
             1.9805E+00
 PARAMETER:  9.2920E-02 -2.1227E+00  4.3233E-01  5.5697E-01 -1.6405E-01 -2.4385E-02  1.3762E+00 -1.2520E+01 -9.4020E-02  2.3068E-01
             7.8333E-01
 GRADIENT:   1.0782E+00  9.3756E-01  8.8626E-01  1.3599E+01 -4.3768E+00 -2.2002E-01 -4.8406E-01  0.0000E+00 -1.9746E+00  8.4066E-01
            -1.9132E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1470.79172074710        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      976
 NPARAMETR:  9.9120E-01  4.3634E-02  1.4507E+00  1.6183E+00  7.7069E-01  8.8235E-01  5.3822E+00  1.0000E-02  8.1545E-01  1.1438E+00
             1.9902E+00
 PARAMETER:  9.1162E-02 -3.0319E+00  4.7204E-01  5.8138E-01 -1.6047E-01 -2.5161E-02  1.7831E+00 -1.8781E+01 -1.0402E-01  2.3432E-01
             7.8825E-01
 GRADIENT:  -9.0921E-01  2.5729E-01  9.5027E-01  7.0362E+00 -2.6704E+00 -5.2589E-02 -2.0469E-01  0.0000E+00 -2.9430E-01  1.6270E-01
            -7.9776E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1470.85509348291        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1155
 NPARAMETR:  9.9085E-01  1.0000E-02  1.4830E+00  1.6407E+00  7.7330E-01  8.8225E-01  1.1939E+01  1.0000E-02  8.0813E-01  1.1470E+00
             1.9923E+00
 PARAMETER:  9.0806E-02 -4.6689E+00  4.9405E-01  5.9512E-01 -1.5709E-01 -2.5283E-02  2.5798E+00 -3.0267E+01 -1.1303E-01  2.3713E-01
             7.8930E-01
 GRADIENT:  -5.3472E-01  0.0000E+00  5.8374E-01  6.8716E+00 -1.4016E+00  9.2830E-02 -2.1273E-02  0.0000E+00 -4.1384E-01 -2.0669E-01
            -8.1892E-02

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1470.88304556827        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1212
 NPARAMETR:  9.9104E-01  1.0000E-02  1.4684E+00  1.6347E+00  7.7099E-01  8.8216E-01  1.2083E+01  1.0000E-02  8.0925E-01  1.1458E+00
             1.9904E+00
 PARAMETER:  9.1002E-02 -4.6784E+00  4.8417E-01  5.9144E-01 -1.6008E-01 -2.5380E-02  2.5918E+00 -3.0366E+01 -1.1165E-01  2.3612E-01
             7.8835E-01
 GRADIENT:   4.7087E-01  0.0000E+00 -4.1870E-01 -2.4709E+00  1.3479E+00  8.6094E-02 -5.4748E-03  0.0000E+00 -2.8552E-02 -2.3835E-01
            -8.1264E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1212
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -8.5064E-04  2.6493E-04  9.7578E-06 -1.2345E-02 -2.9615E-02
 SE:             2.9314E-02  1.9099E-03  1.0700E-04  2.7865E-02  2.2726E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7685E-01  8.8968E-01  9.2734E-01  6.5774E-01  1.9253E-01

 ETASHRINKSD(%)  1.7935E+00  9.3602E+01  9.9642E+01  6.6488E+00  2.3864E+01
 ETASHRINKVR(%)  3.5548E+00  9.9591E+01  9.9999E+01  1.2856E+01  4.2033E+01
 EBVSHRINKSD(%)  1.8671E+00  9.4086E+01  9.9631E+01  6.2736E+00  2.2437E+01
 EBVSHRINKVR(%)  3.6994E+00  9.9650E+01  9.9999E+01  1.2154E+01  3.9840E+01
 RELATIVEINF(%)  9.1246E+01  7.9095E-03  9.8381E-05  2.3956E+00  3.7102E+00
 EPSSHRINKSD(%)  3.5895E+01
 EPSSHRINKVR(%)  5.8906E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1470.8830455682744     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -735.73221900453620     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.20
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1470.883       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.91E-01  1.00E-02  1.47E+00  1.63E+00  7.71E-01  8.82E-01  1.21E+01  1.00E-02  8.09E-01  1.15E+00  1.99E+00
 


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
+        1.40E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        1.58E+00  0.00E+00  7.10E+01
 
 TH 4
+       -3.64E+01  0.00E+00 -2.36E+01  5.60E+02
 
 TH 5
+        1.85E+00  0.00E+00 -2.29E+02 -6.62E+01  8.57E+02
 
 TH 6
+       -8.59E-01  0.00E+00  1.43E+00 -7.27E+00 -2.31E+00  2.36E+02
 
 TH 7
+        2.64E-03  0.00E+00 -1.04E-02  7.68E+00  3.40E-02 -8.62E-04 -1.85E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.46E+00  0.00E+00  5.10E+00 -6.98E+00 -4.42E+00 -2.14E+00 -4.64E-02  0.00E+00  2.34E+02
 
 TH10
+       -3.49E+00  0.00E+00  4.32E+00 -2.25E+00 -5.65E+01  8.09E-01 -2.64E-02  0.00E+00  1.72E+00  5.98E+01
 
 TH11
+       -1.47E+01  0.00E+00 -4.96E+00 -1.10E+01  2.32E+00  4.42E+00 -7.78E-04  0.00E+00  1.14E+01  1.52E+01  6.70E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.655
Stop Time:
Wed Sep 29 12:23:55 CDT 2021
