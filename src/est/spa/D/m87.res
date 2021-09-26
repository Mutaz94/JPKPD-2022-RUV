Sat Sep 25 14:44:34 CDT 2021
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
$DATA ../../../../data/spa/D/dat87.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   22331.2649786302        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2049E+02  3.9745E+02 -8.9695E+01  3.0633E+02  2.1125E+02 -2.7404E+03 -1.1121E+03 -3.5626E+01 -1.7832E+03 -4.7196E+02
            -4.1584E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -507.808806905796        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2562E+00  1.1877E+00  9.5355E-01  1.5535E+00  1.3602E+00  1.7160E+00  1.0498E+00  9.5756E-01  9.2531E-01  9.0509E-01
             1.5139E+01
 PARAMETER:  3.2812E-01  2.7205E-01  5.2438E-02  5.4052E-01  4.0765E-01  6.4002E-01  1.4858E-01  5.6628E-02  2.2370E-02  2.7512E-04
             2.8172E+00
 GRADIENT:  -3.7268E+01  3.3513E+01 -2.2094E+00  5.7559E+01 -7.3955E+00  1.6046E+01 -1.7481E+00  3.3393E+00  2.3958E-01  1.1966E+00
            -1.7459E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -517.612492037045        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.2847E+00  1.0646E+00  1.1949E+00  1.4257E+00  1.5500E+00  1.4812E+00  1.1658E+00  4.0188E-01  6.5463E-01  4.5393E-01
             1.5917E+01
 PARAMETER:  3.5056E-01  1.6264E-01  2.7807E-01  4.5467E-01  5.3823E-01  4.9286E-01  2.5338E-01 -8.1161E-01 -3.2368E-01 -6.8981E-01
             2.8674E+00
 GRADIENT:  -1.3897E+01  1.1298E+01  1.5063E+00  1.1102E+01 -2.4366E+00  4.6573E+00  1.1414E+00  3.2275E-01  1.4308E+00  1.6588E-01
             3.0778E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -523.021375866363        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.2316E+00  5.9094E-01  1.1184E+00  1.6275E+00  2.1336E+00  1.4738E+00  1.7400E+00  8.0748E-02  3.9189E-01  2.3014E+00
             1.5231E+01
 PARAMETER:  3.0829E-01 -4.2604E-01  2.1193E-01  5.8703E-01  8.5780E-01  4.8786E-01  6.5388E-01 -2.4164E+00 -8.3678E-01  9.3350E-01
             2.8233E+00
 GRADIENT:  -5.6787E+00  7.2229E+00  1.2250E+00  1.0987E+01 -3.1581E+00  2.5770E+00  6.6828E-01  1.1667E-02  2.6850E-01  1.3182E+00
             2.4828E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -567.036828354996        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  9.4144E-01  3.3279E-02  1.4067E-01  1.0939E+00  4.7290E+00  1.4584E+00  1.6882E-01  1.0000E-02  1.5184E-02  1.0000E-02
             1.7145E+01
 PARAMETER:  3.9651E-02 -3.3028E+00 -1.8614E+00  1.8974E-01  1.6537E+00  4.7735E-01 -1.6789E+00 -9.6779E+00 -4.0875E+00 -5.0588E+00
             2.9417E+00
 GRADIENT:   5.5729E+00 -1.5968E+00 -3.5431E+01  1.0317E+02  1.0923E+01 -5.0523E+01  1.0566E-04  0.0000E+00  9.3526E-03  0.0000E+00
             6.1869E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -594.015268642168        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  6.1526E-01  1.0000E-02  5.4573E-02  5.3851E-01  9.8064E+00  1.7154E+00  2.2471E+00  1.0000E-02  1.0000E-02  1.0000E-02
             1.2669E+01
 PARAMETER: -3.8571E-01 -4.5979E+00 -2.8082E+00 -5.1894E-01  2.3830E+00  6.3967E-01  9.0963E-01 -1.3027E+01 -4.9682E+00 -7.1537E+00
             2.6392E+00
 GRADIENT:   9.4833E+00  0.0000E+00 -2.6084E+00  1.6202E+01  2.1829E-01  1.3004E+01  1.8369E-03  0.0000E+00  0.0000E+00  0.0000E+00
            -5.1114E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -597.533831895519        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      491
 NPARAMETR:  5.2868E-01  1.0000E-02  3.5368E-02  3.9158E-01  1.7264E+01  1.5886E+00  6.5302E-01  1.0000E-02  1.0000E-02  1.0000E-02
             1.3460E+01
 PARAMETER: -5.3737E-01 -4.9326E+00 -3.2419E+00 -8.3756E-01  2.9486E+00  5.6284E-01 -3.2615E-01 -1.3604E+01 -5.8062E+00 -8.5605E+00
             2.6997E+00
 GRADIENT:  -6.3892E+00  0.0000E+00  1.0916E+01 -1.0746E+01 -1.1433E-02 -3.7459E+00  2.0330E-04  0.0000E+00  0.0000E+00  0.0000E+00
             1.9181E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -597.660584563884        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      672            RESET HESSIAN, TYPE II
 NPARAMETR:  5.2229E-01  1.0000E-02  3.3363E-02  3.7753E-01  9.7115E+00  1.6073E+00  5.5552E-01  1.0000E-02  1.0000E-02  1.0000E-02
             1.3474E+01
 PARAMETER: -5.4953E-01 -4.9616E+00 -3.3003E+00 -8.7411E-01  2.3733E+00  5.7456E-01 -4.8785E-01 -1.3727E+01 -5.8613E+00 -8.7298E+00
             2.7008E+00
 GRADIENT:   7.7221E-01  0.0000E+00  7.2075E+00 -2.8740E+00  4.3592E-02  7.7552E-01  2.1133E-04  0.0000E+00  0.0000E+00  0.0000E+00
             3.7009E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -597.682671725902        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      768
 NPARAMETR:  5.1982E-01  1.0000E-02  3.2914E-02  3.7642E-01  7.3186E+00  1.6043E+00  5.4772E-01  1.0000E-02  1.0000E-02  1.0000E-02
             1.3374E+01
 PARAMETER: -5.5428E-01 -4.9616E+00 -3.3138E+00 -8.7706E-01  2.0904E+00  5.7269E-01 -5.0199E-01 -1.3727E+01 -5.8613E+00 -8.7298E+00
             2.6933E+00
 GRADIENT:  -2.7311E+00  0.0000E+00 -4.1821E-01  1.8293E+00  4.3531E-02 -1.2919E+00  7.6919E-04  0.0000E+00  0.0000E+00  0.0000E+00
            -3.3950E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -597.703863334108        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      948
 NPARAMETR:  5.2394E-01  1.0000E-02  3.3129E-02  3.7811E-01  7.2341E+00  1.6091E+00  3.8031E-01  1.0000E-02  1.0000E-02  1.0000E-02
             1.3457E+01
 PARAMETER: -5.4637E-01 -4.9616E+00 -3.3073E+00 -8.7258E-01  2.0788E+00  5.7567E-01 -8.6678E-01 -1.3727E+01 -5.8613E+00 -8.7298E+00
             2.6995E+00
 GRADIENT:  -8.9697E-02  0.0000E+00  6.3070E-01 -9.7292E-01 -3.0487E-02 -8.1138E-02  3.9660E-04  0.0000E+00  0.0000E+00  0.0000E+00
             4.7979E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -597.704932540896        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1125
 NPARAMETR:  5.2640E-01  1.0000E-02  3.3522E-02  3.8171E-01  7.2008E+00  1.6105E+00  1.1228E-01  1.0000E-02  1.0000E-02  1.0000E-02
             1.3461E+01
 PARAMETER: -5.4170E-01 -4.9616E+00 -3.2956E+00 -8.6310E-01  2.0742E+00  5.7654E-01 -2.0868E+00 -1.3727E+01 -5.8613E+00 -8.7298E+00
             2.6998E+00
 GRADIENT:   1.0417E-02  0.0000E+00 -2.4013E-02  3.2874E-02  1.7793E-03 -1.3903E-03  3.3866E-05  0.0000E+00  0.0000E+00  0.0000E+00
            -1.3547E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -597.704953266360        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1303
 NPARAMETR:  5.2631E-01  1.0000E-02  3.3506E-02  3.8157E-01  7.2032E+00  1.6104E+00  3.9467E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.3462E+01
 PARAMETER: -5.4187E-01 -4.9616E+00 -3.2960E+00 -8.6346E-01  2.0745E+00  5.7649E-01 -3.1323E+00 -1.3727E+01 -5.8613E+00 -8.7298E+00
             2.6998E+00
 GRADIENT:   1.3497E-02  0.0000E+00 -2.3129E-02  2.3414E-02  2.5729E-03 -6.9257E-03  4.1741E-06  0.0000E+00  0.0000E+00  0.0000E+00
             7.1377E-03

0ITERATION NO.:   57    OBJECTIVE VALUE:  -597.704956349766        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1360
 NPARAMETR:  5.2626E-01  1.0000E-02  3.3501E-02  3.8152E-01  7.2006E+00  1.6105E+00  2.7865E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.3461E+01
 PARAMETER: -5.4196E-01 -4.9616E+00 -3.2962E+00 -8.6359E-01  2.0742E+00  5.7652E-01 -3.4804E+00 -1.3727E+01 -5.8613E+00 -8.7298E+00
             2.6998E+00
 GRADIENT:  -1.0063E-02  0.0000E+00  1.3623E-02 -1.7645E-02 -1.0859E-03  2.8552E-03  2.0973E-06  0.0000E+00  0.0000E+00  0.0000E+00
             1.1828E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1360
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8189E-03  7.0455E-06  1.1183E-04 -2.3216E-04  2.1551E-06
 SE:             2.8595E-02  3.8312E-06  2.2881E-04  2.9730E-04  1.4997E-05
 N:                     100         100         100         100         100

 P VAL.:         9.2147E-01  6.5918E-02  6.2501E-01  4.3485E-01  8.8574E-01

 ETASHRINKSD(%)  4.2042E+00  9.9987E+01  9.9233E+01  9.9004E+01  9.9950E+01
 ETASHRINKVR(%)  8.2316E+00  1.0000E+02  9.9994E+01  9.9990E+01  1.0000E+02
 EBVSHRINKSD(%)  4.3157E+00  9.9984E+01  9.9159E+01  9.8917E+01  9.9939E+01
 EBVSHRINKVR(%)  8.4452E+00  1.0000E+02  9.9993E+01  9.9988E+01  1.0000E+02
 RELATIVEINF(%)  5.2752E+00  2.7924E-07  5.4311E-05  9.3850E-05  5.1105E-06
 EPSSHRINKSD(%)  5.3483E+00
 EPSSHRINKVR(%)  1.0411E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -597.70495634976567     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       137.44587021397251     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.92
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.27
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -597.705       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.26E-01  1.00E-02  3.35E-02  3.82E-01  7.20E+00  1.61E+00  2.79E-02  1.00E-02  1.00E-02  1.00E-02  1.35E+01
 


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
+       -5.94E+03  0.00E+00  4.22E+05
 
 TH 4
+       -2.40E+02  0.00E+00 -4.60E+04  5.70E+03
 
 TH 5
+        1.92E+00  0.00E+00 -1.03E+02  1.12E+01  6.84E-02
 
 TH 6
+        2.13E+00  0.00E+00  4.82E+02 -7.25E+01  1.87E-02  6.28E+01
 
 TH 7
+       -7.71E-01  0.00E+00 -4.78E-01  1.90E-01  1.32E-02 -3.37E-03 -3.23E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.75E+01  0.00E+00  1.99E+02 -1.36E+01 -6.01E-02  9.78E-01 -1.93E-03  0.00E+00  0.00E+00  0.00E+00  2.02E+00
 
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
 #CPUT: Total CPU Time in Seconds,       22.255
Stop Time:
Sat Sep 25 14:45:02 CDT 2021
