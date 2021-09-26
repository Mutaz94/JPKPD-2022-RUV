Sat Sep 25 11:38:52 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat31.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m31.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1593.94002852630        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0075E+02 -5.2714E+01 -5.8155E+01  1.8539E+00  1.2245E+02  2.2359E+00  3.6311E+00  8.4846E+00  4.7726E-01 -1.4118E+01
            -9.7360E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1607.88141311838        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9839E-01  1.0000E+00  1.0810E+00  1.0060E+00  9.4556E-01  9.7430E-01  9.2174E-01  9.1817E-01  9.8490E-01  9.5191E-01
             1.3160E+00
 PARAMETER:  9.8393E-02  1.0005E-01  1.7789E-01  1.0597E-01  4.4019E-02  7.3967E-02  1.8508E-02  1.4623E-02  8.4787E-02  5.0717E-02
             3.7463E-01
 GRADIENT:   7.7935E+01 -1.1582E+01 -2.6001E-01 -1.8935E+01  7.1348E+00 -7.1204E+00  2.8912E+00  3.7628E+00 -2.4579E+00 -2.8563E+00
             3.0471E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1609.57114224690        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.9445E-01  9.2550E-01  1.0624E+00  1.0557E+00  9.1731E-01  9.7479E-01  7.6884E-01  6.5447E-01  1.0143E+00  1.0649E+00
             1.2890E+00
 PARAMETER:  9.4431E-02  2.2582E-02  1.6056E-01  1.5417E-01  1.3685E-02  7.4466E-02 -1.6287E-01 -3.2393E-01  1.1423E-01  1.6293E-01
             3.5383E-01
 GRADIENT:   7.0284E+01 -8.9272E+00 -5.5480E+00 -5.2267E+00  5.3216E+00 -6.7512E+00  1.0842E+00  1.7156E+00  1.7452E+00  8.8933E+00
             2.7312E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1611.44603604696        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.6678E-01  9.9536E-01  8.8602E-01  1.0052E+00  8.6278E-01  9.8823E-01  9.1055E-01  4.2774E-01  1.0073E+00  9.5532E-01
             1.2098E+00
 PARAMETER:  6.6217E-02  9.5353E-02 -2.1010E-02  1.0514E-01 -4.7595E-02  8.8163E-02  6.2987E-03 -7.4925E-01  1.0726E-01  5.4294E-02
             2.9048E-01
 GRADIENT:   6.9889E+00 -3.7297E+00 -2.2211E+00 -3.3936E+00  1.2910E+00 -1.1640E+00  5.4762E-01  8.2489E-01  7.7957E-01  1.6791E+00
             1.7115E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1611.49257783990        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.6243E-01  1.0063E+00  8.4762E-01  9.9932E-01  8.4792E-01  9.9300E-01  9.2409E-01  2.9162E-01  1.0017E+00  9.3456E-01
             1.2039E+00
 PARAMETER:  6.1709E-02  1.0632E-01 -6.5318E-02  9.9322E-02 -6.4972E-02  9.2974E-02  2.1053E-02 -1.1323E+00  1.0169E-01  3.2321E-02
             2.8554E-01
 GRADIENT:  -3.4063E+00  1.1954E+00 -5.6694E-01  1.7222E+00 -3.9483E-01  3.4431E-01 -1.3973E-02  3.0655E-01 -2.8059E-01  1.0784E-01
            -7.7518E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1611.64209764145        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.6464E-01  1.0092E+00  8.2384E-01  9.9414E-01  8.3883E-01  9.9180E-01  9.3207E-01  7.7242E-02  1.0017E+00  9.2871E-01
             1.2068E+00
 PARAMETER:  6.3998E-02  1.0917E-01 -9.3776E-02  9.4120E-02 -7.5743E-02  9.1765E-02  2.9653E-02 -2.4608E+00  1.0171E-01  2.6044E-02
             2.8800E-01
 GRADIENT:   9.7885E-01 -8.6741E-01 -4.6965E-01 -5.7713E-01  5.7131E-01 -1.1181E-01 -5.8487E-02  1.7268E-02  7.5159E-02  1.5968E-01
             2.4128E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1611.84590338725        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  9.7739E-01  1.0440E+00  8.2294E-01  9.7717E-01  8.5291E-01  9.9814E-01  9.0896E-01  2.2087E-02  1.0207E+00  9.3786E-01
             1.2100E+00
 PARAMETER:  7.7132E-02  1.4305E-01 -9.4876E-02  7.6903E-02 -5.9097E-02  9.8137E-02  4.5439E-03 -3.7128E+00  1.2053E-01  3.5850E-02
             2.9061E-01
 GRADIENT:   1.6939E+00  7.3517E-01  1.9557E-01  9.4509E-01 -6.0612E-01 -2.8639E-01  3.0333E-02  9.8513E-04  1.6188E-02  1.2153E-01
             5.2555E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1611.85084141739        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  9.7673E-01  1.0785E+00  8.1424E-01  9.5465E-01  8.6426E-01  9.9893E-01  8.8592E-01  1.2710E-02  1.0404E+00  9.4132E-01
             1.2088E+00
 PARAMETER:  7.6454E-02  1.7560E-01 -1.0550E-01  5.3589E-02 -4.5885E-02  9.8934E-02 -2.1124E-02 -4.2653E+00  1.3964E-01  3.9524E-02
             2.8961E-01
 GRADIENT:   2.6466E-02  2.9636E-03 -1.0300E-02  6.2904E-03  1.0325E-02 -1.1755E-03  2.9664E-03  3.4854E-04  4.3417E-03  4.6555E-03
             1.0397E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1611.85091181156        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  9.7671E-01  1.0773E+00  8.1453E-01  9.5539E-01  8.6386E-01  9.9894E-01  8.8666E-01  1.0000E-02  1.0397E+00  9.4117E-01
             1.2087E+00
 PARAMETER:  7.6439E-02  1.7450E-01 -1.0515E-01  5.4370E-02 -4.6341E-02  9.8936E-02 -2.0295E-02 -4.5554E+00  1.3898E-01  3.9364E-02
             2.8958E-01
 GRADIENT:   4.0068E-04 -6.4646E-04  4.6752E-06 -1.1857E-03 -1.1484E-04 -3.7423E-05 -2.6448E-05  0.0000E+00 -3.4581E-04  1.7026E-04
             6.5478E-05

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1611.85091181156        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      908
 NPARAMETR:  9.7671E-01  1.0773E+00  8.1453E-01  9.5539E-01  8.6386E-01  9.9894E-01  8.8666E-01  1.0000E-02  1.0397E+00  9.4117E-01
             1.2087E+00
 PARAMETER:  7.6439E-02  1.7450E-01 -1.0515E-01  5.4370E-02 -4.6341E-02  9.8936E-02 -2.0295E-02 -4.5554E+00  1.3898E-01  3.9364E-02
             2.8958E-01
 GRADIENT:   4.0068E-04 -6.4646E-04  4.6752E-06 -1.1857E-03 -1.1484E-04 -3.7423E-05 -2.6448E-05  0.0000E+00 -3.4581E-04  1.7026E-04
             6.5478E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      908
 NO. OF SIG. DIGITS IN FINAL EST.:  4.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.4372E-05 -1.6677E-02 -2.8947E-04  4.1302E-03 -2.2084E-02
 SE:             2.9751E-02  1.7776E-02  1.4995E-04  2.5260E-02  2.3866E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9854E-01  3.4816E-01  5.3556E-02  8.7012E-01  3.5480E-01

 ETASHRINKSD(%)  3.2974E-01  4.0447E+01  9.9498E+01  1.5376E+01  2.0046E+01
 ETASHRINKVR(%)  6.5839E-01  6.4534E+01  9.9997E+01  2.8387E+01  3.6074E+01
 EBVSHRINKSD(%)  6.2670E-01  4.0129E+01  9.9529E+01  1.5492E+01  1.8962E+01
 EBVSHRINKVR(%)  1.2495E+00  6.4155E+01  9.9998E+01  2.8585E+01  3.4329E+01
 RELATIVEINF(%)  9.8437E+01  1.1586E+00  2.4703E-04  3.1478E+00  5.9995E+00
 EPSSHRINKSD(%)  4.2170E+01
 EPSSHRINKVR(%)  6.6557E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1611.8509118115585     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -876.70008524782031     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.05
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1611.851       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.77E-01  1.08E+00  8.15E-01  9.55E-01  8.64E-01  9.99E-01  8.87E-01  1.00E-02  1.04E+00  9.41E-01  1.21E+00
 


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
+        1.15E+03
 
 TH 2
+       -9.04E+00  4.57E+02
 
 TH 3
+        1.41E+01  2.06E+02  3.96E+02
 
 TH 4
+       -1.47E+01  4.01E+02 -1.45E+02  8.05E+02
 
 TH 5
+       -1.87E+00 -3.97E+02 -5.58E+02  1.62E+02  1.10E+03
 
 TH 6
+        2.20E-01 -1.76E+00  8.90E-01 -5.01E+00 -3.44E+00  1.93E+02
 
 TH 7
+        1.99E+00  6.69E+00  5.74E+00 -7.43E+00 -1.77E+01  6.40E-01  3.37E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.30E+00 -2.67E+01 -1.72E+01  2.68E+01 -7.52E-01  1.08E+00  2.56E+01  0.00E+00  9.68E+01
 
 TH10
+       -9.51E-01 -1.76E+00 -4.41E+01 -1.28E+01 -5.28E+01  5.64E-02  1.46E+01  0.00E+00  3.52E+00  9.52E+01
 
 TH11
+       -9.72E+00 -1.79E+01 -3.03E+01 -4.50E+00  7.22E+00  3.44E+00  5.68E+00  0.00E+00  9.43E+00  2.12E+01  1.51E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.748
Stop Time:
Sat Sep 25 11:39:08 CDT 2021
