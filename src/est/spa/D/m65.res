Wed Sep 29 20:14:15 CDT 2021
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
$DATA ../../../../data/spa/D/dat65.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m65.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   13640.6713057158        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.2380E+02  3.5318E+02 -1.3593E+01  8.3654E+01  4.0118E+02 -2.3438E+03 -8.7602E+02 -1.0319E+02 -1.6951E+03 -1.0037E+03
            -2.4257E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -609.070357341007        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3238E+00  9.9073E-01  9.5598E-01  1.4193E+00  1.2088E+00  1.7331E+00  1.1971E+00  9.8562E-01  1.3247E+00  1.2087E+00
             1.4854E+01
 PARAMETER:  3.8050E-01  9.0684E-02  5.4986E-02  4.5018E-01  2.8964E-01  6.4990E-01  2.7993E-01  8.5513E-02  3.8116E-01  2.8958E-01
             2.7983E+00
 GRADIENT:  -1.1923E+01  4.1043E+00 -4.1043E+00  2.0955E+00 -8.9406E+00  2.8412E+01  1.2964E+00  4.0567E+00  1.3277E+01  4.2863E+00
             1.2075E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -618.083093053773        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3361E+00  8.4901E-01  1.1065E+00  1.6636E+00  2.8538E+00  1.5331E+00  3.4072E+00  5.3421E-01  1.3092E+00  6.3568E+00
             1.3067E+01
 PARAMETER:  3.8976E-01 -6.3687E-02  2.0123E-01  6.0901E-01  1.1487E+00  5.2728E-01  1.3259E+00 -5.2698E-01  3.6945E-01  1.9495E+00
             2.6701E+00
 GRADIENT:   2.1384E+01  2.2619E+01 -3.0179E+00  3.7530E+01 -1.7465E+01 -4.3372E+01  8.1110E+00  3.6569E-01  2.0215E+01  1.0975E+01
             2.6632E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -645.547525260147        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.1325E+00  2.6009E-01  7.7265E-01  1.6296E+00  4.1335E+00  1.5662E+00  2.3756E+00  2.2546E-02  4.5093E-01  7.3827E+00
             1.2947E+01
 PARAMETER:  2.2447E-01 -1.2467E+00 -1.5793E-01  5.8832E-01  1.5191E+00  5.4868E-01  9.6525E-01 -3.6922E+00 -6.9644E-01  2.0991E+00
             2.6608E+00
 GRADIENT:  -3.5628E+01  8.4039E+00  1.2155E+01  2.9235E+01 -3.4890E+00 -2.2725E+00  8.0246E-01 -2.4641E-04  2.9786E+00  9.3146E+00
             2.7288E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -713.119041855283        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.4268E-01  1.0000E-02  1.2560E-01  1.0016E+00  4.5799E+00  1.6593E+00  3.7109E-01  1.0000E-02  2.3375E-02  1.6000E+00
             1.0949E+01
 PARAMETER:  4.0970E-02 -5.7651E+00 -1.9746E+00  1.0160E-01  1.6217E+00  6.0641E-01 -8.9131E-01 -1.2859E+01 -3.6561E+00  5.7003E-01
             2.4933E+00
 GRADIENT:   6.8553E+01  0.0000E+00 -4.3792E+01  1.6204E+02 -7.6594E+00 -1.1003E+02  3.4856E-04  0.0000E+00  1.3541E-03  3.5960E+00
            -8.1972E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -740.643677155775        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      399
 NPARAMETR:  6.7996E-01  1.0000E-02  4.4048E-02  4.6356E-01  5.8990E+00  1.8408E+00  2.9591E-02  1.0000E-02  1.0000E-02  1.0875E+00
             8.7703E+00
 PARAMETER: -2.8572E-01 -1.1162E+01 -3.0225E+00 -6.6882E-01  1.8748E+00  7.1019E-01 -3.4203E+00 -2.3284E+01 -7.9653E+00  1.8388E-01
             2.2714E+00
 GRADIENT:   9.4910E+01  0.0000E+00 -1.9569E+01  2.6936E+01  1.3890E+01 -6.9379E+01 -2.5425E-05  0.0000E+00  0.0000E+00 -1.0940E+01
            -2.3815E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -782.162359674812        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      577
 NPARAMETR:  4.1031E-01  1.0000E-02  1.5639E-02  2.0741E-01  7.2431E+00  1.9969E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.1273E+00
             1.0459E+01
 PARAMETER: -7.9085E-01 -1.6840E+01 -4.0580E+00 -1.4731E+00  2.0800E+00  7.9160E-01 -6.8939E+00 -3.9753E+01 -1.4455E+01  2.1980E-01
             2.4475E+00
 GRADIENT:   1.9351E+00  0.0000E+00 -9.3305E+00  3.3882E+00 -7.8688E+00  2.0235E+00  0.0000E+00  0.0000E+00  0.0000E+00  3.6840E+00
            -1.3705E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -782.370954269769        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      754
 NPARAMETR:  4.3977E-01  1.0000E-02  1.8770E-02  2.3958E-01  7.4702E+00  1.9706E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.1263E+00
             1.0736E+01
 PARAMETER: -7.2150E-01 -1.5996E+01 -3.8755E+00 -1.3289E+00  2.1109E+00  7.7832E-01 -6.5017E+00 -3.7955E+01 -1.3651E+01  2.1892E-01
             2.4736E+00
 GRADIENT:  -1.2534E+00  0.0000E+00 -8.4219E-01  3.1072E+00  1.1623E+00 -5.4867E+00  0.0000E+00  0.0000E+00  0.0000E+00 -4.5165E-01
             4.1583E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -783.270747657148        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      892
 NPARAMETR:  4.3320E-01  1.0000E-02  1.8530E-02  2.3683E-01  3.7891E+03  1.9904E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.1587E+00
             1.0642E+01
 PARAMETER: -7.3655E-01 -1.6050E+01 -3.8884E+00 -1.3404E+00  8.3399E+00  7.8836E-01 -6.5257E+00 -3.8145E+01 -1.3701E+01  2.4732E-01
             2.4648E+00
 GRADIENT:   6.9868E-01  0.0000E+00 -1.5888E+01  2.1388E+01  2.4928E-04 -3.9016E+00  0.0000E+00  0.0000E+00  0.0000E+00 -5.0104E-09
            -5.3115E+00

0ITERATION NO.:   42    OBJECTIVE VALUE:  -783.361913845728        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      949
 NPARAMETR:  4.3331E-01  1.0000E-02  1.8691E-02  2.3612E-01  3.4859E+03  2.0006E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.1586E+00
             1.0686E+01
 PARAMETER: -7.3631E-01 -1.6050E+01 -3.8797E+00 -1.3434E+00  8.2565E+00  7.9345E-01 -6.5257E+00 -3.8145E+01 -1.3701E+01  2.4720E-01
             2.4689E+00
 GRADIENT:  -9.2775E-01  0.0000E+00  2.7132E+00 -8.1503E-01  1.6345E-04 -1.3432E+00  0.0000E+00  0.0000E+00  0.0000E+00 -5.7486E-09
            -4.8832E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      949
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8236E-03  2.9807E-06  7.5253E-05 -1.9270E-04  6.5229E-08
 SE:             2.9451E-02  9.6797E-07  3.3483E-04  3.9811E-04  3.2741E-07
 N:                     100         100         100         100         100

 P VAL.:         9.2362E-01  2.0748E-03  8.2217E-01  6.2836E-01  8.4208E-01

 ETASHRINKSD(%)  1.3346E+00  9.9997E+01  9.8878E+01  9.8666E+01  9.9999E+01
 ETASHRINKVR(%)  2.6515E+00  1.0000E+02  9.9987E+01  9.9982E+01  1.0000E+02
 EBVSHRINKSD(%)  2.0289E+00  9.9995E+01  9.8932E+01  9.8722E+01  9.9999E+01
 EBVSHRINKVR(%)  4.0166E+00  1.0000E+02  9.9989E+01  9.9984E+01  1.0000E+02
 RELATIVEINF(%)  9.6797E-01  9.2376E-08  2.4023E-05  2.4307E-05  1.3449E-10
 EPSSHRINKSD(%)  7.3839E+00
 EPSSHRINKVR(%)  1.4223E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -783.36191384572817     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -48.211087281989990     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.60
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -783.362       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.33E-01  1.00E-02  1.87E-02  2.36E-01  3.49E+03  2.00E+00  1.00E-02  1.00E-02  1.00E-02  1.16E+00  1.07E+01
 


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
+        1.41E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -9.85E+03  0.00E+00  2.09E+06
 
 TH 4
+       -2.92E+02  0.00E+00 -1.95E+05  1.95E+04
 
 TH 5
+        2.24E-06  0.00E+00 -6.89E-05  5.14E-06 -6.78E-12
 
 TH 6
+        3.65E+00  0.00E+00  5.30E+02 -7.42E+01  3.97E-08  4.56E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        1.03E-04  0.00E+00  6.13E-04 -3.33E-04  1.74E-09 -1.89E-05  0.00E+00  0.00E+00  0.00E+00 -7.26E-04
 
 TH11
+       -1.38E+01  0.00E+00  2.84E+02 -1.96E+01 -3.09E-08  7.02E-01  0.00E+00  0.00E+00  0.00E+00  4.93E-06  3.54E+00
 
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
 #CPUT: Total CPU Time in Seconds,       18.787
Stop Time:
Wed Sep 29 20:14:35 CDT 2021
