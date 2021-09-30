Thu Sep 30 10:03:52 CDT 2021
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
$DATA ../../../../data/spa2/D/dat91.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   42820.6957070603        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2253E+03  6.5769E+02 -1.9559E+01  7.6907E+02 -7.3198E+01 -3.1278E+03 -2.1311E+03 -2.8635E+01 -2.2333E+03 -4.7582E+02
            -8.1418E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -289.718326753186        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.1177E+00  1.2265E+00  1.0361E+00  1.2588E+00  9.7066E-01  1.7864E+00  1.2299E+00  9.8161E-01  1.1389E+00  9.5857E-01
             1.4647E+01
 PARAMETER:  2.1125E-01  3.0413E-01  1.3547E-01  3.3018E-01  7.0217E-02  6.8019E-01  3.0697E-01  8.1436E-02  2.3008E-01  5.7686E-02
             2.7843E+00
 GRADIENT:  -2.7728E+01  1.8481E+01  2.4545E+00  2.9886E+01 -1.5734E+01  3.8635E+00 -2.9933E+01  2.5394E+00 -1.9016E+01  7.5391E+00
            -6.6764E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -372.338803089146        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.1299E+00  1.0023E+00  5.0353E+00  1.6050E+00  1.9321E+00  2.2158E+00  5.2533E+00  3.6430E-01  1.8604E+00  3.7463E-02
             1.4480E+01
 PARAMETER:  2.2212E-01  1.0231E-01  1.7165E+00  5.7314E-01  7.5863E-01  8.9562E-01  1.7589E+00 -9.0978E-01  7.2078E-01 -3.1844E+00
             2.7728E+00
 GRADIENT:  -2.1193E+01 -5.1933E+00 -5.7924E+00 -1.1398E+01  1.1143E+01  1.5129E+01  3.0525E+01  4.9072E-02  2.5169E+01  1.0005E-02
             7.2118E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -381.473006638662        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      277
 NPARAMETR:  1.1662E+00  9.4241E-01  5.5943E+00  1.7148E+00  1.9421E+00  2.0995E+00  5.8664E+00  3.0516E-01  1.2398E+00  7.1988E-02
             1.4519E+01
 PARAMETER:  2.5372E-01  4.0683E-02  1.8217E+00  6.3930E-01  7.6378E-01  8.4170E-01  1.8692E+00 -1.0869E+00  3.1498E-01 -2.5313E+00
             2.7755E+00
 GRADIENT:  -1.3275E+01  3.2471E+00 -3.1571E+00  1.2447E+01  3.8690E+00 -1.0919E+01 -1.2289E-01  3.3589E-02 -1.0386E+00  3.2205E-02
             4.6565E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -382.874664970397        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      454
 NPARAMETR:  1.2050E+00  7.6776E-01  1.0356E+01  1.7812E+00  2.0286E+00  2.1802E+00  6.2001E+00  2.4364E-01  1.3164E+00  6.6875E-02
             1.4377E+01
 PARAMETER:  2.8649E-01 -1.6428E-01  2.4376E+00  6.7729E-01  8.0736E-01  8.7940E-01  1.9246E+00 -1.3121E+00  3.7493E-01 -2.6049E+00
             2.7656E+00
 GRADIENT:   2.3129E+00 -2.2283E-01 -3.7019E-01  4.2373E+00  2.5741E-01  2.3145E-01 -1.0268E+00  5.9511E-03 -2.9168E-01  2.8220E-02
            -4.3687E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -382.945154155085        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      613
 NPARAMETR:  1.1957E+00  8.1130E-01  1.3330E+01  1.7450E+00  2.0357E+00  2.1800E+00  6.2633E+00  1.1136E-01  1.2840E+00  2.5585E-02
             1.4339E+01
 PARAMETER:  2.7869E-01 -1.0912E-01  2.6900E+00  6.5677E-01  8.1086E-01  8.7933E-01  1.9347E+00 -2.0950E+00  3.4997E-01 -3.5658E+00
             2.7630E+00
 GRADIENT:   8.8888E+00  1.6463E+00  1.7556E-01  1.2066E+01 -1.3644E+00  1.5615E+01  4.1174E+01  6.8532E-04  3.5477E-01  4.4840E-03
             2.9027E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -382.984220252263        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      789
 NPARAMETR:  1.1977E+00  7.8710E-01  1.3210E+01  1.7512E+00  2.0700E+00  2.1814E+00  6.2510E+00  8.3781E-02  1.2874E+00  1.9118E-02
             1.4383E+01
 PARAMETER:  2.8040E-01 -1.3940E-01  2.6810E+00  6.6027E-01  8.2753E-01  8.7998E-01  1.9327E+00 -2.3796E+00  3.5260E-01 -3.8571E+00
             2.7661E+00
 GRADIENT:   4.3237E-02 -5.3155E-02 -2.0806E-02 -3.7919E-01  1.1267E-01  3.8231E-02  2.0271E+00  3.8555E-04  4.9757E-02  2.3092E-03
             2.3557E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -382.987124565901        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      972             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1980E+00  7.8658E-01  1.3370E+01  1.7529E+00  2.0712E+00  2.1874E+00  6.2666E+00  4.7919E-02  1.2886E+00  1.0000E-02
             1.4378E+01
 PARAMETER:  2.8066E-01 -1.4005E-01  2.6930E+00  6.6128E-01  8.2813E-01  8.8274E-01  1.9352E+00 -2.9382E+00  3.5354E-01 -4.5844E+00
             2.7657E+00
 GRADIENT:   9.1867E+00  5.4000E-01  1.5526E-02  1.0008E+01  6.2890E-01  1.6822E+01  4.0563E+01  1.2770E-04  7.3341E-01  0.0000E+00
             3.2980E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -382.987415607904        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:     1091
 NPARAMETR:  1.1975E+00  7.8282E-01  1.3487E+01  1.7539E+00  2.0709E+00  2.1879E+00  6.2675E+00  2.3631E-02  1.2879E+00  1.0000E-02
             1.4374E+01
 PARAMETER:  2.8023E-01 -1.4485E-01  2.7017E+00  6.6186E-01  8.2797E-01  8.8295E-01  1.9354E+00 -3.6452E+00  3.5302E-01 -4.5844E+00
             2.7654E+00
 GRADIENT:   9.0866E+00  4.5589E-01  2.8261E-02  1.0139E+01  5.5409E-01  1.6941E+01  4.0416E+01  3.0841E-05  6.3703E-01  0.0000E+00
             3.2775E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -382.987614290830        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1255             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1975E+00  7.8151E-01  1.3512E+01  1.7550E+00  2.0714E+00  2.1880E+00  6.2738E+00  1.7766E-02  1.2884E+00  1.0000E-02
             1.4374E+01
 PARAMETER:  2.8027E-01 -1.4653E-01  2.7036E+00  6.6249E-01  8.2824E-01  8.8299E-01  1.9364E+00 -3.9304E+00  3.5343E-01 -4.5844E+00
             2.7654E+00
 GRADIENT:   9.1007E+00  4.7351E-01  2.7397E-02  1.0179E+01  5.6259E-01  1.6949E+01  4.0539E+01  1.7460E-05  6.3692E-01  0.0000E+00
             3.2736E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -382.987719207463        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1441
 NPARAMETR:  1.1975E+00  7.8132E-01  1.3527E+01  1.7559E+00  2.0716E+00  2.1879E+00  6.2781E+00  1.0000E-02  1.2890E+00  1.0000E-02
             1.4374E+01
 PARAMETER:  2.8028E-01 -1.4677E-01  2.7047E+00  6.6301E-01  8.2830E-01  8.8293E-01  1.9371E+00 -5.3082E+00  3.5384E-01 -4.5844E+00
             2.7654E+00
 GRADIENT:   9.0459E-02  2.3834E-02  4.8743E-03  1.6896E-02 -7.0186E-02  8.0920E-01  2.2205E+00  0.0000E+00 -6.8424E-02  0.0000E+00
            -5.3700E-01

0ITERATION NO.:   51    OBJECTIVE VALUE:  -382.987719207463        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1463
 NPARAMETR:  1.1975E+00  7.8132E-01  1.3527E+01  1.7559E+00  2.0716E+00  2.1879E+00  6.2781E+00  1.0000E-02  1.2890E+00  1.0000E-02
             1.4374E+01
 PARAMETER:  2.8028E-01 -1.4677E-01  2.7047E+00  6.6301E-01  8.2830E-01  8.8293E-01  1.9371E+00 -5.3082E+00  3.5384E-01 -4.5844E+00
             2.7654E+00
 GRADIENT:   9.0459E-02  2.3834E-02  4.8743E-03  1.6896E-02 -7.0186E-02  8.0920E-01  2.2205E+00  0.0000E+00 -6.8424E-02  0.0000E+00
            -5.3700E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1463
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6985E-02  4.4505E-02  4.5085E-06 -7.6085E-02 -5.7297E-06
 SE:             2.7807E-02  2.2016E-02  2.5338E-06  1.4235E-02  3.6865E-05
 N:                     100         100         100         100         100

 P VAL.:         5.4132E-01  4.3228E-02  7.5184E-02  9.0631E-08  8.7649E-01

 ETASHRINKSD(%)  6.8417E+00  2.6243E+01  9.9992E+01  5.2311E+01  9.9876E+01
 ETASHRINKVR(%)  1.3215E+01  4.5600E+01  1.0000E+02  7.7258E+01  1.0000E+02
 EBVSHRINKSD(%)  7.6265E+00  2.4234E+01  9.9981E+01  5.1240E+01  9.9783E+01
 EBVSHRINKVR(%)  1.4671E+01  4.2596E+01  1.0000E+02  7.6224E+01  1.0000E+02
 RELATIVEINF(%)  8.3995E+01  2.6658E+01  4.3031E-07  1.0233E+01  5.3325E-05
 EPSSHRINKSD(%)  4.7496E+00
 EPSSHRINKVR(%)  9.2735E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -382.98771920746310     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       719.73852063814400     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.50
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    11.96
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -382.988       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.20E+00  7.81E-01  1.35E+01  1.76E+00  2.07E+00  2.19E+00  6.28E+00  1.00E-02  1.29E+00  1.00E-02  1.44E+01
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.18E+02
 
 TH 2
+       -8.76E+01  5.09E+01
 
 TH 3
+       -1.24E-01  2.12E-02  1.79E-04
 
 TH 4
+       -1.57E+02  5.84E+01  1.32E-01  1.41E+02
 
 TH 5
+        2.09E+01 -7.30E+00 -2.11E-02 -2.01E+01  2.99E+00
 
 TH 6
+       -3.62E+01 -1.08E+00  6.03E-02  2.85E+01 -5.70E+00  3.40E+01
 
 TH 7
+       -1.09E+01  1.58E+01 -2.51E-02 -3.85E+00  1.22E+00 -9.01E+00  1.14E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.57E+01 -1.08E+01 -5.77E-02 -4.97E+01  7.42E+00 -1.39E+01  7.65E+00  0.00E+00  2.10E+01
 
 TH10
+        4.30E-65 -6.66E-65  1.06E-67  8.09E-66 -4.61E-66  4.95E-65 -4.56E-65  0.00E+00 -2.77E-65  4.5E-111
 
 TH11
+        6.00E-01 -2.26E+00 -4.43E-03 -5.40E+00  7.73E-01  1.20E-01  1.66E-01  0.00E+00  1.92E+00  1.02E-47  6.08E-01
 
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
+        1.33E+02
 
 TH 2
+       -3.44E+00  2.78E+01
 
 TH 3
+        9.35E-04  2.24E-02  3.86E-03
 
 TH 4
+       -4.92E+00  2.14E+01  3.01E-02  6.63E+01
 
 TH 5
+       -5.95E-01 -4.16E+00 -1.53E-01 -7.87E+00  1.14E+01
 
 TH 6
+       -1.63E+00 -2.41E+00  1.02E-02  1.62E+00 -7.26E-01  2.79E+01
 
 TH 7
+        7.17E-02  3.85E+00 -4.81E-03 -3.87E+00  2.30E-01 -7.95E-01  2.86E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.83E-02 -2.43E+00 -3.37E-02 -2.25E+01  3.01E+00 -8.44E-01  1.60E+00  0.00E+00  1.96E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  5.49E+00
 
 TH11
+       -5.40E+00 -2.16E+00 -1.36E-03 -5.44E+00  6.23E-01  6.72E-01  1.72E-01  0.00E+00  1.63E+00  0.00E+00  2.75E+00
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.39E+02
 
 TH 2
+        3.97E+01  2.81E+01
 
 TH 3
+        5.45E-02  2.65E-02  6.02E-04
 
 TH 4
+        5.95E+01  2.28E+01  5.25E-02  6.94E+01
 
 TH 5
+       -1.07E+01 -2.38E+00 -4.11E-02 -9.56E+00  4.08E+00
 
 TH 6
+        1.59E+01  6.93E+00  7.12E-03 -8.52E+00 -1.42E+00  3.31E+01
 
 TH 7
+       -1.36E+00  4.34E+00 -3.88E-03 -4.51E+00  1.46E+00  1.56E+00  2.82E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -8.90E+00 -1.94E+00 -9.52E-03 -1.92E+01  2.82E+00  6.22E+00  1.29E+00  0.00E+00  1.42E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.72E+01 -3.72E+00 -2.03E-02 -7.20E+00  1.99E+00  4.22E-01  6.67E-01  0.00E+00  1.07E+00  0.00E+00  2.83E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       46.570
Stop Time:
Thu Sep 30 10:04:41 CDT 2021
