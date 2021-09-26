Fri Sep 24 21:45:43 CDT 2021
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
$DATA ../../../../data/int/A2/dat79.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1859.35808601689        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8304E+01  1.4964E+02  1.9124E+02 -1.1431E+02  1.3032E+02  3.0682E+01 -1.1353E+02 -4.5081E+02 -8.2009E+01 -8.4596E+01
            -3.3247E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3070.74279583232        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0816E+00  1.0584E+00  8.2436E-01  1.0303E+00  9.5573E-01  9.0761E-01  9.8234E-01  4.3973E-01  8.6870E-01  9.9013E-01
             2.1055E+00
 PARAMETER:  1.7840E-01  1.5674E-01 -9.3153E-02  1.2984E-01  5.4716E-02  3.0575E-03  8.2187E-02 -7.2160E-01 -4.0756E-02  9.0085E-02
             8.4455E-01
 GRADIENT:   1.5895E+02  8.2791E+01 -7.5248E+01  5.2114E+01  1.9179E+01 -1.1319E+01 -9.6326E+00  3.2202E+00 -1.9051E+01  1.5125E+00
            -1.1452E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3077.38433288358        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0892E+00  7.8992E-01  6.9981E-01  1.1540E+00  7.3114E-01  9.4108E-01  9.8951E-01  3.3666E-01  9.3436E-01  9.5459E-01
             2.1101E+00
 PARAMETER:  1.8547E-01 -1.3582E-01 -2.5694E-01  2.4323E-01 -2.1315E-01  3.9275E-02  8.9456E-02 -9.8867E-01  3.2111E-02  5.3529E-02
             8.4672E-01
 GRADIENT:   1.6949E+02  5.3757E+01 -8.6633E+01  5.1750E+01  2.6323E+01  1.0457E+00  1.7445E+00  3.5620E+00  2.9226E+00  6.5290E+00
             3.6716E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3088.25775514655        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0175E+00  7.3487E-01  7.6509E-01  1.1519E+00  7.2469E-01  9.1812E-01  1.0240E+00  3.4760E-01  9.0626E-01  9.2320E-01
             2.0815E+00
 PARAMETER:  1.1733E-01 -2.0806E-01 -1.6777E-01  2.4137E-01 -2.2201E-01  1.4576E-02  1.2369E-01 -9.5671E-01  1.5695E-03  2.0088E-02
             8.3310E-01
 GRADIENT:  -3.0683E+00 -1.2807E+00  4.8942E-01  1.9636E+00 -1.3798E+00  1.1776E+00  6.6248E-01  3.1001E+00  6.5514E-01 -2.5991E-01
             1.7472E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3089.71014860223        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.0209E+00  7.1848E-01  7.4415E-01  1.1540E+00  7.1060E-01  9.1339E-01  1.0144E+00  6.7014E-02  9.1498E-01  9.4901E-01
             2.0594E+00
 PARAMETER:  1.2069E-01 -2.3061E-01 -1.9552E-01  2.4323E-01 -2.4164E-01  9.4118E-03  1.1431E-01 -2.6028E+00  1.1150E-02  4.7659E-02
             8.2241E-01
 GRADIENT:   7.2354E+00 -6.0985E+00 -2.9222E+00 -2.7602E+00  7.7345E+00 -7.0662E-01  2.2996E+00  1.3152E-01  2.1277E+00  1.5690E+00
            -4.5616E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3089.83346348522        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0182E+00  7.2415E-01  7.4675E-01  1.1547E+00  7.1081E-01  9.1513E-01  9.7027E-01  4.1871E-02  9.1157E-01  9.5722E-01
             2.0639E+00
 PARAMETER:  1.1808E-01 -2.2276E-01 -1.9203E-01  2.4384E-01 -2.4135E-01  1.1310E-02  6.9824E-02 -3.0732E+00  7.4138E-03  5.6278E-02
             8.2460E-01
 GRADIENT:  -2.0248E-01  7.1271E-01  2.9105E-01  9.0248E-01 -4.9277E-01  1.0801E-02 -2.8641E-02  5.1732E-02 -7.8238E-02  1.7480E-02
            -2.8835E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3089.83554732267        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  1.0183E+00  7.2241E-01  7.4513E-01  1.1550E+00  7.0928E-01  9.1508E-01  9.7117E-01  2.6461E-02  9.1164E-01  9.5638E-01
             2.0639E+00
 PARAMETER:  1.1812E-01 -2.2516E-01 -1.9420E-01  2.4407E-01 -2.4351E-01  1.1261E-02  7.0750E-02 -3.5321E+00  7.4864E-03  5.5400E-02
             8.2460E-01
 GRADIENT:  -7.2763E-02  2.3349E-01  1.1092E-01  2.8176E-01 -1.8864E-01  1.1046E-03 -6.2936E-03  2.0870E-02 -2.3206E-02 -3.0535E-03
            -9.2414E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3089.83808139354        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      520
 NPARAMETR:  1.0183E+00  7.2167E-01  7.4441E-01  1.1551E+00  7.0863E-01  9.1506E-01  9.7128E-01  1.0892E-02  9.1169E-01  9.5620E-01
             2.0639E+00
 PARAMETER:  1.1814E-01 -2.2619E-01 -1.9517E-01  2.4416E-01 -2.4442E-01  1.1240E-02  7.0856E-02 -4.4198E+00  7.5413E-03  5.5217E-02
             8.2461E-01
 GRADIENT:  -7.2289E-03  3.9233E-03 -6.6766E-03  3.4802E-05 -2.9071E-03 -2.2455E-03 -1.0610E-03  3.5762E-03 -2.2409E-03  4.5266E-03
             1.5040E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3089.99417746716        NO. OF FUNC. EVALS.: 154
 CUMULATIVE NO. OF FUNC. EVALS.:      674
 NPARAMETR:  1.0232E+00  7.5040E-01  7.7197E-01  1.1480E+00  7.3826E-01  9.1694E-01  9.6680E-01  1.0000E-02  9.1363E-01  9.6849E-01
             2.0681E+00
 PARAMETER:  1.2297E-01 -1.8715E-01 -1.5881E-01  2.3806E-01 -2.0345E-01  1.3284E-02  6.6236E-02 -4.5537E+00  9.6737E-03  6.7984E-02
             8.2663E-01
 GRADIENT:   1.0699E-02 -4.7210E-05  3.6146E-02  1.6947E-02 -2.5307E-02  1.1814E-02  6.6816E-03  0.0000E+00  2.8090E-03 -1.1075E-02
             2.0514E-02

0ITERATION NO.:   42    OBJECTIVE VALUE:  -3089.99418019573        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      731
 NPARAMETR:  1.0232E+00  7.5041E-01  7.7194E-01  1.1480E+00  7.3826E-01  9.1691E-01  9.6659E-01  1.0000E-02  9.1365E-01  9.6860E-01
             2.0681E+00
 PARAMETER:  1.2297E-01 -1.8714E-01 -1.5884E-01  2.3805E-01 -2.0346E-01  1.3251E-02  6.6016E-02 -4.5538E+00  9.6887E-03  6.8094E-02
             8.2662E-01
 GRADIENT:   3.6070E-03  6.1751E-04  8.5461E-03  5.2854E-03 -7.0711E-03 -4.3694E-04  1.3997E-03  0.0000E+00 -1.1584E-03 -2.8355E-03
            -1.1420E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      731
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0986E-03 -1.4837E-02 -1.2650E-04  1.7303E-03 -1.0892E-02
 SE:             2.9518E-02  1.9207E-02  1.7235E-04  2.7681E-02  2.5645E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7031E-01  4.3983E-01  4.6295E-01  9.5016E-01  6.7104E-01

 ETASHRINKSD(%)  1.1105E+00  3.5656E+01  9.9423E+01  7.2662E+00  1.4086E+01
 ETASHRINKVR(%)  2.2087E+00  5.8598E+01  9.9997E+01  1.4004E+01  2.6187E+01
 EBVSHRINKSD(%)  1.2456E+00  3.5362E+01  9.9337E+01  7.3939E+00  1.5004E+01
 EBVSHRINKVR(%)  2.4757E+00  5.8219E+01  9.9996E+01  1.4241E+01  2.7757E+01
 RELATIVEINF(%)  9.7488E+01  1.0220E+01  1.2349E-03  5.4207E+01  1.0638E+01
 EPSSHRINKSD(%)  1.7640E+01
 EPSSHRINKVR(%)  3.2169E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3089.9941801957339     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1435.9048204273231     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.33
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3089.994       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  7.50E-01  7.72E-01  1.15E+00  7.38E-01  9.17E-01  9.67E-01  1.00E-02  9.14E-01  9.69E-01  2.07E+00
 


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
+        1.23E+03
 
 TH 2
+       -9.04E+00  8.76E+02
 
 TH 3
+        1.83E+00  3.97E+01  8.01E+02
 
 TH 4
+       -1.53E+01  2.74E+02 -9.65E+01  8.64E+02
 
 TH 5
+       -3.32E+00 -7.46E+02 -7.66E+02  1.42E+02  1.55E+03
 
 TH 6
+        7.88E+00 -4.35E+00  5.14E+00 -6.13E+00 -2.93E+00  2.25E+02
 
 TH 7
+       -8.80E-01  5.31E+00 -6.99E+00  1.57E+00 -8.91E+00 -3.90E-01  3.39E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.15E+00 -2.52E+01  2.78E+01  1.12E+01 -4.75E+00  1.12E-02  1.46E+01  0.00E+00  1.76E+02
 
 TH10
+       -1.84E+00 -3.44E-01 -3.71E+01 -7.77E-01 -3.06E+00 -1.47E-01  3.59E+01  0.00E+00 -6.09E+00  1.08E+02
 
 TH11
+       -1.64E+01 -2.13E+01 -1.57E+01 -1.25E+01  4.83E+00  3.29E+00  2.26E+00  0.00E+00  7.36E+00  6.71E+00  2.72E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       25.985
Stop Time:
Fri Sep 24 21:46:11 CDT 2021
