Sat Sep 25 08:44:26 CDT 2021
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
$DATA ../../../../data/spa/A2/dat57.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1193.43494971493        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.8063E+01 -5.6650E+01 -9.2006E-01 -7.9144E+01  8.0075E+01  3.8545E+01 -6.9176E+00  1.5635E+00 -1.3524E+01 -1.9051E+01
            -9.6036E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1488.47062624166        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0307E+00  1.0227E+00  1.0755E+00  1.0334E+00  9.9416E-01  8.0848E-01  9.1100E-01  9.1896E-01  9.2456E-01  8.1393E-01
             1.8626E+00
 PARAMETER:  1.3023E-01  1.2243E-01  1.7280E-01  1.3283E-01  9.4142E-02 -1.1260E-01  6.7848E-03  1.5490E-02  2.1558E-02 -1.0588E-01
             7.2198E-01
 GRADIENT:   3.0737E+00 -2.2478E+01  2.8447E+00 -4.0357E+01  2.8066E+01 -3.5176E+01 -3.4366E-01  8.4899E-01 -1.1116E+01 -4.3054E+00
            -1.4945E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1501.92238901010        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0404E+00  1.3063E+00  5.9659E-01  8.7529E-01  8.0568E-01  8.8665E-01  8.7692E-01  6.2305E-01  1.0226E+00  4.7719E-01
             2.0465E+00
 PARAMETER:  1.3960E-01  3.6717E-01 -4.1653E-01 -3.3197E-02 -1.1606E-01 -2.0301E-02 -3.1335E-02 -3.7313E-01  1.2232E-01 -6.3984E-01
             8.1612E-01
 GRADIENT:   8.3618E+00  9.6084E+01  4.6933E+01  4.7086E+00 -9.4264E+01  5.6882E-02  2.8261E+00  2.1127E-01 -1.0002E+00 -2.6211E-01
            -8.1048E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1512.04659466765        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.0407E+00  1.0630E+00  5.6362E-01  9.8968E-01  7.1332E-01  8.7483E-01  9.8729E-01  4.3299E-01  8.5363E-01  3.2481E-01
             2.3294E+00
 PARAMETER:  1.3990E-01  1.6111E-01 -4.7337E-01  8.9627E-02 -2.3782E-01 -3.3722E-02  8.7209E-02 -7.3704E-01 -5.8262E-02 -1.0245E+00
             9.4561E-01
 GRADIENT:  -1.8371E+00  7.6103E+00  4.8543E+00 -4.9281E+00 -9.9449E+00 -2.8699E+00  1.9908E+00  1.1416E+00 -1.3734E+00  1.8939E+00
            -3.9916E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1512.82291251631        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      324
 NPARAMETR:  1.0422E+00  9.3090E-01  4.4502E-01  1.0327E+00  5.8408E-01  8.8830E-01  1.0681E+00  1.8547E-01  8.1326E-01  2.0121E-01
             2.2970E+00
 PARAMETER:  1.4134E-01  2.8398E-02 -7.0964E-01  1.3218E-01 -4.3771E-01 -1.8446E-02  1.6587E-01 -1.5849E+00 -1.0670E-01 -1.5034E+00
             9.3160E-01
 GRADIENT:  -9.4069E+00  1.5192E+00 -5.5432E-01  1.2062E+00 -2.2635E+00  1.0575E-01  7.2463E-01  1.6989E-01 -1.4578E-01  5.8323E-01
             4.9485E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1512.97449485065        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      500
 NPARAMETR:  1.0455E+00  9.3147E-01  4.9552E-01  1.0467E+00  6.2018E-01  8.8489E-01  1.0977E+00  1.3894E-01  8.1336E-01  1.7211E-01
             2.3356E+00
 PARAMETER:  1.4453E-01  2.9010E-02 -6.0215E-01  1.4568E-01 -3.7775E-01 -2.2293E-02  1.9322E-01 -1.8737E+00 -1.0658E-01 -1.6596E+00
             9.4826E-01
 GRADIENT:  -1.0202E+00 -1.6503E+00 -1.8017E+00  1.0162E+00  2.8919E+00 -9.3965E-02  9.9283E-01  3.5497E-02  4.1815E-01  1.3509E-01
             1.0501E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1513.20709117252        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      677
 NPARAMETR:  1.0447E+00  1.0653E+00  4.0788E-01  9.4914E-01  6.1681E-01  8.8718E-01  9.5624E-01  1.0000E-02  8.6110E-01  5.7284E-02
             2.3156E+00
 PARAMETER:  1.4370E-01  1.6322E-01 -7.9679E-01  4.7797E-02 -3.8320E-01 -1.9706E-02  5.5258E-02 -4.6370E+00 -4.9545E-02 -2.7597E+00
             9.3965E-01
 GRADIENT:  -6.8828E-01 -8.6951E-01 -8.2493E-01  8.4634E-03  1.5603E+00 -6.8169E-02  8.4777E-01  0.0000E+00  1.5486E-01  4.7820E-02
             6.2117E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1513.74904896241        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      856
 NPARAMETR:  1.0407E+00  1.1210E+00  2.3480E-01  8.3448E-01  5.1671E-01  8.9234E-01  8.2898E-01  1.0000E-02  8.8415E-01  1.0000E-02
             2.2403E+00
 PARAMETER:  1.3993E-01  2.1425E-01 -1.3490E+00 -8.0943E-02 -5.6026E-01 -1.3908E-02 -8.7554E-02 -1.3258E+01 -2.3127E-02 -6.4468E+00
             9.0660E-01
 GRADIENT:   1.5778E+01  5.0918E+00  4.7246E+00 -6.2087E+00 -1.0447E+01  2.9694E-01  2.5791E+00  0.0000E+00 -3.7242E+00  0.0000E+00
             8.1634E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1517.89899346499        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1034
 NPARAMETR:  1.0125E+00  1.5551E+00  1.3037E-01  5.7784E-01  6.4580E-01  8.8285E-01  6.6342E-01  1.0000E-02  1.1913E+00  1.0000E-02
             2.1375E+00
 PARAMETER:  1.1240E-01  5.4156E-01 -1.9374E+00 -4.4846E-01 -3.3726E-01 -2.4603E-02 -3.1034E-01 -4.3040E+01  2.7505E-01 -1.8914E+01
             8.5963E-01
 GRADIENT:  -2.4989E+01  2.2839E+01 -4.3231E-01  3.0053E+01 -3.2469E+00 -9.8128E-01 -1.2805E+00  0.0000E+00 -3.1672E-01  0.0000E+00
            -1.0932E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1528.07524749887        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1210
 NPARAMETR:  1.0271E+00  2.0596E+00  5.2155E-02  3.1059E-01  8.9713E-01  8.7842E-01  6.2835E-01  1.0000E-02  1.8683E+00  1.0000E-02
             2.1857E+00
 PARAMETER:  1.2674E-01  8.2249E-01 -2.8535E+00 -1.0693E+00 -8.5504E-03 -2.9625E-02 -3.6465E-01 -1.0105E+02  7.2502E-01 -4.3488E+01
             8.8193E-01
 GRADIENT:   1.1658E+01  3.5752E+01 -4.7837E+00  3.1772E+01  1.1886E+01 -6.4987E-01  6.7847E+00  0.0000E+00 -6.3454E+00  0.0000E+00
             5.8765E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1528.91206483282        NO. OF FUNC. EVALS.: 205
 CUMULATIVE NO. OF FUNC. EVALS.:     1415             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0265E+00  2.1032E+00  4.7636E-02  2.9172E-01  9.2305E-01  8.7822E-01  6.2170E-01  1.0000E-02  1.9791E+00  1.0000E-02
             2.1807E+00
 PARAMETER:  1.2616E-01  8.4348E-01 -2.9442E+00 -1.1320E+00  1.9931E-02 -2.9860E-02 -3.7530E-01 -1.0675E+02  7.8267E-01 -4.5872E+01
             8.7963E-01
 GRADIENT:   1.9040E+01  7.2209E+01 -2.9970E+00  3.2031E+01  1.0592E+01 -1.0747E-01  5.5687E+00  0.0000E+00 -6.8945E+00  0.0000E+00
             4.3531E+00

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1528.91206483282        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1488
 NPARAMETR:  1.0265E+00  2.1033E+00  4.7638E-02  2.9173E-01  9.2305E-01  8.7831E-01  6.2170E-01  1.0000E-02  1.9792E+00  1.0000E-02
             2.1806E+00
 PARAMETER:  1.2616E-01  8.4348E-01 -2.9442E+00 -1.1320E+00  1.9931E-02 -2.9860E-02 -3.7530E-01 -1.0675E+02  7.8267E-01 -4.5872E+01
             8.7963E-01
 GRADIENT:  -2.2644E+05 -3.3823E+04 -4.8433E+03 -2.5221E+04  1.4285E+05 -7.3189E-01 -7.6113E+04  0.0000E+00 -1.8285E+04  0.0000E+00
             1.6239E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1488
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.1704E-03 -2.8092E-02  1.5104E-04 -1.9302E-02 -3.2821E-04
 SE:             2.9326E-02  2.5700E-02  7.3259E-05  2.1037E-02  2.5202E-04
 N:                     100         100         100         100         100

 P VAL.:         9.1391E-01  2.7436E-01  3.9237E-02  3.5887E-01  1.9282E-01

 ETASHRINKSD(%)  1.7529E+00  1.3903E+01  9.9755E+01  2.9524E+01  9.9156E+01
 ETASHRINKVR(%)  3.4751E+00  2.5873E+01  9.9999E+01  5.0331E+01  9.9993E+01
 EBVSHRINKSD(%)  2.0810E+00  1.3327E+01  9.9783E+01  2.7126E+01  9.9072E+01
 EBVSHRINKVR(%)  4.1186E+00  2.4878E+01  1.0000E+02  4.6893E+01  9.9991E+01
 RELATIVEINF(%)  9.3180E+01  1.5169E+01  1.9925E-04  1.2346E+01  1.6242E-03
 EPSSHRINKSD(%)  3.1953E+01
 EPSSHRINKVR(%)  5.3696E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1528.9120648328239     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -793.76123826908577     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.92
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0T MATRIX UNOBTAINABLE
 Elapsed covariance  time in seconds:     6.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1528.912       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.10E+00  4.76E-02  2.92E-01  9.23E-01  8.78E-01  6.22E-01  1.00E-02  1.98E+00  1.00E-02  2.18E+00
 


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
+        4.26E+08
 
 TH 2
+       -9.28E+02  2.27E+06
 
 TH 3
+       -1.16E+04  2.86E+07  3.61E+08
 
 TH 4
+        1.83E+05  1.35E+04  1.66E+05  6.56E+07
 
 TH 5
+       -8.18E+03 -1.03E+03 -8.03E+03 -2.34E+08  8.38E+08
 
 TH 6
+       -1.20E+04 -8.84E+02 -1.11E+04 -4.73E+03  1.68E+04  2.39E+02
 
 TH 7
+       -8.55E+03 -6.56E+02 -7.97E+03 -3.34E+03  1.20E+04 -6.66E+03  1.31E+08
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.08E+05  7.86E+03  9.96E+04  1.40E+07 -5.00E+07 -1.00E+03 -7.10E+02  0.00E+00  2.99E+06
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        7.92E+02  8.86E+02  1.14E+04 -1.24E+04  5.28E+02  8.13E+02  5.97E+02  0.00E+00 -7.26E+03  0.00E+00  1.94E+06
 
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
 #CPUT: Total CPU Time in Seconds,       24.440
Stop Time:
Sat Sep 25 08:44:52 CDT 2021
