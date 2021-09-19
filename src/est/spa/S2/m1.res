Sat Sep 18 13:11:30 CDT 2021
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
$DATA ../../../../data/spa/S2/dat1.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1671.92268951276        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4592E+01 -5.7936E+00  3.0398E+01 -3.6125E+01 -3.8291E+01  8.7311E+00  2.5176E-02 -1.1031E+00  6.8767E+00 -2.6725E+01
             2.4436E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1676.26825331587        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0138E+00  1.0439E+00  9.5751E-01  9.9051E-01  1.0810E+00  9.6408E-01  9.9802E-01  9.5655E-01  9.3986E-01  1.3489E+00
             8.9865E-01
 PARAMETER:  1.1372E-01  1.4297E-01  5.6576E-02  9.0461E-02  1.7787E-01  6.3422E-02  9.8021E-02  5.5581E-02  3.7980E-02  3.9926E-01
            -6.8668E-03
 GRADIENT:   9.5018E+01 -1.1351E+01 -2.1401E+01  1.2710E-01  9.3279E+00 -5.7174E+00  1.9354E+00  5.4837E+00 -8.1943E+00  1.2625E+01
            -1.1699E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1678.57524363553        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.9946E-01  1.2351E+00  9.6928E-01  8.8983E-01  1.1769E+00  9.7849E-01  8.2109E-01  1.0165E+00  1.1188E+00  1.3810E+00
             8.9700E-01
 PARAMETER:  9.9457E-02  3.1117E-01  6.8798E-02 -1.6720E-02  2.6288E-01  7.8250E-02 -9.7128E-02  1.1641E-01  2.1229E-01  4.2278E-01
            -8.7038E-03
 GRADIENT:   5.3456E+01  3.3671E+01 -1.0207E+00  2.5315E+01  2.2129E-01  8.3627E-01  2.1384E+00 -5.1030E-01  2.3262E+00  4.9131E+00
            -1.3267E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1679.25502160264        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      316
 NPARAMETR:  1.0024E+00  1.3178E+00  9.9213E-01  8.2685E-01  1.2293E+00  9.9794E-01  7.2994E-01  1.1938E+00  1.1914E+00  1.3632E+00
             9.3708E-01
 PARAMETER:  1.0238E-01  3.7593E-01  9.2099E-02 -9.0132E-02  3.0648E-01  9.7937E-02 -2.1479E-01  2.7715E-01  2.7511E-01  4.0987E-01
             3.5014E-02
 GRADIENT:   4.9742E+00  4.3530E+00  3.5227E+00  3.7714E+00 -2.7387E+00  2.8968E+00 -7.5941E-01 -6.6966E-01 -9.0829E-01 -2.6050E+00
             4.8408E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1680.11265754016        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      492
 NPARAMETR:  9.9919E-01  1.5851E+00  7.5108E-01  6.4721E-01  1.2838E+00  9.8603E-01  6.8750E-01  1.0775E+00  1.4187E+00  1.4063E+00
             9.1228E-01
 PARAMETER:  9.9189E-02  5.6067E-01 -1.8624E-01 -3.3509E-01  3.4984E-01  8.5931E-02 -2.7469E-01  1.7460E-01  4.4971E-01  4.4099E-01
             8.1917E-03
 GRADIENT:  -5.0640E+00  4.4476E+00 -7.0908E-03  5.3796E+00 -2.5151E+00 -2.6330E+00 -8.7098E-01  9.4750E-02  2.0045E-01  2.3000E+00
            -5.6340E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1680.28544570169        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      678
 NPARAMETR:  9.9978E-01  1.6337E+00  7.0232E-01  6.1436E-01  1.2906E+00  9.8753E-01  6.8656E-01  1.0275E+00  1.4633E+00  1.3989E+00
             9.1452E-01
 PARAMETER:  9.9783E-02  5.9087E-01 -2.5336E-01 -3.8718E-01  3.5508E-01  8.7448E-02 -2.7606E-01  1.2714E-01  4.8071E-01  4.3570E-01
             1.0648E-02
 GRADIENT:  -4.2645E+00  5.2614E+00 -5.2256E-02  5.3704E+00 -2.3849E+00 -2.1434E+00 -7.1533E-01  1.4602E-01  1.3459E-01  1.9172E+00
            -4.6968E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1680.71337563538        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      858
 NPARAMETR:  1.0077E+00  1.8086E+00  5.3908E-01  4.9302E-01  1.3369E+00  1.0040E+00  6.7928E-01  7.8362E-01  1.6569E+00  1.3689E+00
             9.3853E-01
 PARAMETER:  1.0770E-01  6.9255E-01 -5.1790E-01 -6.0720E-01  3.9033E-01  1.0396E-01 -2.8672E-01 -1.4383E-01  6.0493E-01  4.1400E-01
             3.6559E-02
 GRADIENT:   1.1613E+01  1.5873E+00 -8.4122E-01  3.1091E+00  5.2427E+00  4.0428E+00 -1.3741E-01  1.7359E-01 -1.0751E+00 -2.1113E+00
             5.3474E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1680.97515528163        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1033
 NPARAMETR:  1.0001E+00  1.9149E+00  4.5455E-01  4.1968E-01  1.3684E+00  9.8926E-01  6.6815E-01  6.2892E-01  1.8318E+00  1.3942E+00
             9.2410E-01
 PARAMETER:  1.0015E-01  7.4966E-01 -6.8844E-01 -7.6827E-01  4.1362E-01  8.9205E-02 -3.0324E-01 -3.6375E-01  7.0528E-01  4.3235E-01
             2.1070E-02
 GRADIENT:  -5.0144E+00  4.2542E+00 -6.3460E-01  1.9773E+00  7.1433E-01 -1.8410E+00  1.6729E-01  1.7392E-01 -4.9426E-01  6.9137E-03
            -8.8186E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1681.00587987247        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1209
 NPARAMETR:  1.0025E+00  1.9467E+00  4.3130E-01  3.9506E-01  1.3818E+00  9.9419E-01  6.6290E-01  5.6322E-01  1.9074E+00  1.3999E+00
             9.2621E-01
 PARAMETER:  1.0247E-01  7.6613E-01 -7.4094E-01 -8.2872E-01  4.2337E-01  9.4171E-02 -3.1113E-01 -4.7409E-01  7.4574E-01  4.3643E-01
             2.3343E-02
 GRADIENT:   2.8155E-01 -2.1527E-01 -2.0888E-01 -3.5104E-02 -1.5131E-01  1.3768E-01  8.2956E-03  1.2820E-01  1.5021E-01  7.7059E-02
             1.8585E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1681.04176377946        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1387
 NPARAMETR:  1.0024E+00  1.9420E+00  4.2047E-01  3.9743E-01  1.3705E+00  9.9382E-01  6.6819E-01  2.3029E-01  1.8886E+00  1.3890E+00
             9.2736E-01
 PARAMETER:  1.0238E-01  7.6373E-01 -7.6639E-01 -8.2273E-01  4.1521E-01  9.3802E-02 -3.0318E-01 -1.3684E+00  7.3582E-01  4.2856E-01
             2.4585E-02
 GRADIENT:  -7.3691E-02 -4.2344E-01 -6.0662E-02 -6.7006E-02  3.4262E-01 -2.4378E-02  5.2347E-03  8.8286E-03 -1.2469E-01 -2.9229E-02
            -4.8674E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1681.04546388178        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1562
 NPARAMETR:  1.0024E+00  1.9445E+00  4.1621E-01  3.9575E-01  1.3690E+00  9.9391E-01  6.6867E-01  5.1460E-02  1.8930E+00  1.3872E+00
             9.2778E-01
 PARAMETER:  1.0242E-01  7.6499E-01 -7.7658E-01 -8.2698E-01  4.1411E-01  9.3888E-02 -3.0247E-01 -2.8669E+00  7.3815E-01  4.2727E-01
             2.5043E-02
 GRADIENT:  -3.7763E-03 -5.0713E-02 -3.9585E-02  3.4151E-02  3.3082E-02  2.7312E-03  2.2977E-02  3.7482E-04  4.5284E-03  5.8711E-03
             1.1636E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1681.04557816591        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1738
 NPARAMETR:  1.0024E+00  1.9452E+00  4.1579E-01  3.9524E-01  1.3694E+00  9.9390E-01  6.6852E-01  1.0000E-02  1.8946E+00  1.3874E+00
             9.2779E-01
 PARAMETER:  1.0242E-01  7.6535E-01 -7.7757E-01 -8.2826E-01  4.1436E-01  9.3879E-02 -3.0269E-01 -4.6477E+00  7.3899E-01  4.2741E-01
             2.5046E-02
 GRADIENT:  -3.4449E-04 -9.5409E-03 -4.1927E-03  1.9391E-03  2.9691E-03 -3.7235E-04  1.4325E-03  0.0000E+00  8.2655E-04  5.4887E-04
             1.7678E-03

0ITERATION NO.:   56    OBJECTIVE VALUE:  -1681.04557816591        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1760
 NPARAMETR:  1.0024E+00  1.9452E+00  4.1579E-01  3.9524E-01  1.3694E+00  9.9390E-01  6.6852E-01  1.0000E-02  1.8946E+00  1.3874E+00
             9.2779E-01
 PARAMETER:  1.0242E-01  7.6535E-01 -7.7757E-01 -8.2826E-01  4.1436E-01  9.3879E-02 -3.0269E-01 -4.6477E+00  7.3899E-01  4.2741E-01
             2.5046E-02
 GRADIENT:  -3.4449E-04 -9.5409E-03 -4.1927E-03  1.9391E-03  2.9691E-03 -3.7235E-04  1.4325E-03  0.0000E+00  8.2655E-04  5.4887E-04
             1.7678E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1760
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5120E-04 -5.0469E-02 -2.6692E-04  4.3182E-02 -5.1760E-02
 SE:             2.9891E-02  2.2890E-02  8.7955E-05  2.2858E-02  2.2984E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9329E-01  2.7463E-02  2.4075E-03  5.8873E-02  2.4322E-02

 ETASHRINKSD(%)  1.0000E-10  2.3316E+01  9.9705E+01  2.3423E+01  2.3001E+01
 ETASHRINKVR(%)  1.0000E-10  4.1196E+01  9.9999E+01  4.1360E+01  4.0711E+01
 EBVSHRINKSD(%)  3.8010E-01  2.1374E+01  9.9766E+01  2.6009E+01  1.9097E+01
 EBVSHRINKVR(%)  7.5876E-01  3.8179E+01  9.9999E+01  4.5253E+01  3.4546E+01
 RELATIVEINF(%)  9.9217E+01  6.4137E+00  1.0034E-04  5.9510E+00  2.5328E+01
 EPSSHRINKSD(%)  4.5350E+01
 EPSSHRINKVR(%)  7.0134E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1681.0455781659118     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -945.89475160217364     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.08
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1681.046       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.95E+00  4.16E-01  3.95E-01  1.37E+00  9.94E-01  6.69E-01  1.00E-02  1.89E+00  1.39E+00  9.28E-01
 


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
+        1.12E+03
 
 TH 2
+       -6.11E+00  4.52E+02
 
 TH 3
+        5.86E+00  1.89E+02  4.38E+02
 
 TH 4
+       -1.30E+01  3.59E+02 -3.16E+02  1.13E+03
 
 TH 5
+       -1.82E+00 -9.47E+01 -1.50E+02  1.35E+02  2.06E+02
 
 TH 6
+        1.57E+00 -2.52E-01  1.17E+00 -3.33E+00  1.61E+00  2.09E+02
 
 TH 7
+        1.43E+00 -1.29E+01 -1.73E+01 -1.21E+01 -1.48E+01 -1.05E+00  1.98E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.03E+00 -2.14E+01 -3.54E+01  6.79E+01 -2.84E-01 -2.69E-01  1.49E+01  0.00E+00  2.43E+01
 
 TH10
+       -2.27E-01 -1.50E+01 -2.91E+01  3.57E+00 -3.66E+01  1.03E+00  7.44E+00  0.00E+00  3.50E+00  4.97E+01
 
 TH11
+       -5.46E+00 -2.16E+01 -3.05E+01  2.50E+00  1.11E+00  3.16E+00  1.49E+01  0.00E+00  3.32E+00  1.05E+01  2.51E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.545
Stop Time:
Sat Sep 18 13:12:00 CDT 2021
