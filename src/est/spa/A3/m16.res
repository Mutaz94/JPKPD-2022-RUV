Wed Sep 29 13:20:08 CDT 2021
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
$DATA ../../../../data/spa/A3/dat16.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m16.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   95.3320357937159        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7384E+02 -2.0265E+00  5.5988E+01 -6.4481E+01  2.4603E+02  8.2144E+01 -5.4983E+01 -5.2911E+01 -1.4341E+02 -1.8379E+02
            -3.0970E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1025.24428403969        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0646E+00  9.5579E-01  8.0186E-01  1.1370E+00  7.7900E-01  8.3345E-01  1.1124E+00  1.0520E+00  1.2780E+00  1.1889E+00
             1.8278E+00
 PARAMETER:  1.6257E-01  5.4785E-02 -1.2083E-01  2.2835E-01 -1.4974E-01 -8.2181E-02  2.0648E-01  1.5070E-01  3.4530E-01  2.7301E-01
             7.0312E-01
 GRADIENT:   3.3301E+02  4.3240E+01  1.5514E+01  7.2009E+01  4.0708E+01 -7.0932E-01 -7.9272E+00 -2.7102E+00 -9.4915E+00 -1.7758E+01
            -7.7189E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1040.41313382287        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      200
 NPARAMETR:  1.0654E+00  7.4362E-01  4.1016E-01  1.1450E+00  3.8358E-01  7.4634E-01  1.2099E+00  1.0206E+00  1.2282E+00  5.5286E-01
             1.8314E+00
 PARAMETER:  1.6339E-01 -1.9622E-01 -7.9120E-01  2.3543E-01 -8.5822E-01 -1.9257E-01  2.9053E-01  1.2037E-01  3.0559E-01 -4.9266E-01
             7.0508E-01
 GRADIENT:   1.7976E+02  1.6701E+02  2.3715E+02 -1.1648E+01 -2.5875E+02 -5.2412E+01 -2.0352E+01 -5.9367E+01 -2.7105E+01 -2.3520E+01
            -6.8255E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1173.50903924914        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  1.0515E+00  7.4773E-01  3.5547E-01  1.0756E+00  3.6210E-01  8.0809E-01  5.9177E-01  2.1343E+00  1.3204E+00  4.6815E-01
             2.1059E+00
 PARAMETER:  1.5024E-01 -1.9071E-01 -9.3432E-01  1.7291E-01 -9.1583E-01 -1.1308E-01 -4.2464E-01  8.5816E-01  3.7792E-01 -6.5897E-01
             8.4473E-01
 GRADIENT:   1.3691E+02  1.0027E+02  1.4135E+02  7.3516E+00 -1.0870E+02 -1.1568E+01 -1.2879E+01 -4.8319E+01 -3.4486E+01 -2.4341E+00
            -3.1058E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1277.73831478856        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      554
 NPARAMETR:  1.0204E+00  5.3534E-01  1.7542E-01  1.0483E+00  2.6459E-01  8.1373E-01  3.2105E-01  1.6476E+00  1.4922E+00  1.8546E-01
             3.1136E+00
 PARAMETER:  1.2022E-01 -5.2486E-01 -1.6406E+00  1.4718E-01 -1.2296E+00 -1.0613E-01 -1.0361E+00  5.9934E-01  5.0022E-01 -1.5849E+00
             1.2358E+00
 GRADIENT:   4.5265E+00 -4.0461E+00  4.3541E+00 -8.7869E+00 -6.5946E+00 -5.6766E+00 -4.4105E+00  4.1264E-01 -1.1370E+00 -2.9365E+00
            -2.4569E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1284.07180109059        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  1.0194E+00  6.5639E-01  1.7140E-01  1.0139E+00  2.8842E-01  8.1261E-01  7.6519E-01  1.2575E+00  1.3700E+00  1.7702E-01
             3.3760E+00
 PARAMETER:  1.1918E-01 -3.2100E-01 -1.6637E+00  1.1379E-01 -1.1433E+00 -1.0751E-01 -1.6763E-01  3.2916E-01  4.1481E-01 -1.6315E+00
             1.3167E+00
 GRADIENT:  -4.6312E+00  2.1071E+01  1.8333E+01 -2.6067E+00 -3.1423E+01 -7.6318E+00 -6.4042E+00 -3.3557E+00 -6.7217E+00  9.1510E-02
             1.5771E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1285.37762037783        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      906
 NPARAMETR:  1.0165E+00  6.3336E-01  1.6322E-01  1.0209E+00  2.7922E-01  8.3430E-01  8.9278E-01  1.1263E+00  1.4348E+00  1.5053E-01
             3.2648E+00
 PARAMETER:  1.1637E-01 -3.5672E-01 -1.7127E+00  1.2071E-01 -1.1757E+00 -8.1165E-02 -1.3420E-02  2.1895E-01  4.6099E-01 -1.7936E+00
             1.2832E+00
 GRADIENT:  -4.6191E+00  1.0030E+01  1.0765E+01  4.8077E+00 -1.9260E+01 -2.3535E-01 -3.2408E+00 -5.5679E+00  1.7178E-01 -2.4405E-03
             4.5910E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1286.33700243220        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1083
 NPARAMETR:  1.0123E+00  6.4856E-01  1.5607E-01  9.9729E-01  2.8140E-01  8.3308E-01  9.1223E-01  1.3740E+00  1.4507E+00  1.2220E-01
             3.1439E+00
 PARAMETER:  1.1220E-01 -3.3300E-01 -1.7574E+00  9.7284E-02 -1.1680E+00 -8.2631E-02  8.1316E-03  4.1775E-01  4.7208E-01 -2.0021E+00
             1.2455E+00
 GRADIENT:   4.8488E-01 -2.2726E+00 -5.5222E-01 -4.0112E+00  2.9334E+00  8.5369E-03 -9.6374E-02 -1.0315E-01 -1.3989E+00  1.5472E-01
            -9.6298E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1286.38244817985        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1265
 NPARAMETR:  1.0124E+00  6.4742E-01  1.5612E-01  1.0030E+00  2.8069E-01  8.3314E-01  9.1099E-01  1.3747E+00  1.4615E+00  7.7501E-02
             3.1450E+00
 PARAMETER:  1.1232E-01 -3.3476E-01 -1.7571E+00  1.0295E-01 -1.1705E+00 -8.2554E-02  6.7755E-03  4.1825E-01  4.7948E-01 -2.4575E+00
             1.2458E+00
 GRADIENT:   2.5302E-01  9.0034E-01  4.0223E-01  1.7619E-01 -3.6560E+00 -1.4243E-01 -8.1457E-01 -3.6196E-01 -2.0905E-01  2.0140E-02
            -1.4838E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1286.39817204323        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1442
 NPARAMETR:  1.0117E+00  6.4785E-01  1.5572E-01  1.0018E+00  2.8085E-01  8.3327E-01  9.2640E-01  1.3808E+00  1.4600E+00  3.0292E-02
             3.1469E+00
 PARAMETER:  1.1164E-01 -3.3410E-01 -1.7597E+00  1.0184E-01 -1.1699E+00 -8.2399E-02  2.3549E-02  4.2268E-01  4.7841E-01 -3.3969E+00
             1.2464E+00
 GRADIENT:  -1.4520E+00 -3.9174E-01 -6.3620E-01 -2.1867E-01 -2.2080E+00 -7.0636E-02  5.1410E-02  2.0416E-04 -2.1526E-01  2.9152E-03
            -3.3479E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1286.40574585445        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1618
 NPARAMETR:  1.0119E+00  6.4906E-01  1.5621E-01  1.0023E+00  2.8154E-01  8.3304E-01  9.2563E-01  1.3804E+00  1.4577E+00  1.0000E-02
             3.1502E+00
 PARAMETER:  1.1184E-01 -3.3224E-01 -1.7566E+00  1.0233E-01 -1.1675E+00 -8.2669E-02  2.2722E-02  4.2234E-01  4.7683E-01 -4.6170E+00
             1.2475E+00
 GRADIENT:  -1.1932E+00 -8.4804E-01 -6.8456E-01  3.0928E-03 -8.4586E-01 -8.4314E-02  1.8174E-02  7.8050E-03 -2.5434E-01  0.0000E+00
            -1.4645E-01

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1286.40574585445        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1640
 NPARAMETR:  1.0119E+00  6.4906E-01  1.5621E-01  1.0023E+00  2.8154E-01  8.3304E-01  9.2563E-01  1.3804E+00  1.4577E+00  1.0000E-02
             3.1502E+00
 PARAMETER:  1.1184E-01 -3.3224E-01 -1.7566E+00  1.0233E-01 -1.1675E+00 -8.2669E-02  2.2722E-02  4.2234E-01  4.7683E-01 -4.6170E+00
             1.2475E+00
 GRADIENT:  -1.1932E+00 -8.4804E-01 -6.8456E-01  3.0928E-03 -8.4586E-01 -8.4314E-02  1.8174E-02  7.8050E-03 -2.5434E-01  0.0000E+00
            -1.4645E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1640
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0325E-03 -5.6074E-03  2.3287E-04 -9.1629E-03  3.1093E-04
 SE:             2.8329E-02  1.9188E-02  1.6993E-02  2.5719E-02  3.3090E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4280E-01  7.7011E-01  9.8907E-01  7.2164E-01  3.4740E-01

 ETASHRINKSD(%)  5.0942E+00  3.5718E+01  4.3072E+01  1.3839E+01  9.8891E+01
 ETASHRINKVR(%)  9.9289E+00  5.8678E+01  6.7592E+01  2.5762E+01  9.9988E+01
 EBVSHRINKSD(%)  5.0787E+00  3.5241E+01  4.2710E+01  1.2535E+01  9.8951E+01
 EBVSHRINKVR(%)  9.8995E+00  5.8063E+01  6.7178E+01  2.3499E+01  9.9989E+01
 RELATIVEINF(%)  8.3791E+01  1.5675E+00  1.0291E+01  4.7099E+01  3.3623E-04
 EPSSHRINKSD(%)  3.3673E+01
 EPSSHRINKVR(%)  5.6008E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1286.4057458544521     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -551.25491929071393     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.68
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.07
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1286.406       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  6.49E-01  1.56E-01  1.00E+00  2.82E-01  8.33E-01  9.26E-01  1.38E+00  1.46E+00  1.00E-02  3.15E+00
 


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
+        1.42E+03
 
 TH 2
+       -4.10E+01  1.29E+03
 
 TH 3
+       -4.66E+02  1.93E+03  7.42E+03
 
 TH 4
+       -3.92E+01  1.26E+02 -4.64E+02  3.99E+02
 
 TH 5
+        3.64E+02 -4.31E+03 -8.00E+03 -2.08E+01  1.76E+04
 
 TH 6
+        6.10E-01 -1.58E+01  2.84E+01 -1.31E+01  9.59E+01  2.35E+02
 
 TH 7
+       -4.98E+00 -9.15E+00 -6.11E+01 -1.09E+00  1.08E+02  1.23E+00  3.92E+01
 
 TH 8
+        2.49E+00 -7.68E+00 -1.23E+01 -2.51E+00  2.23E+01  2.06E+00  5.88E+00  1.12E+01
 
 TH 9
+        1.68E+01 -2.89E+01  9.09E+01 -6.07E+00  1.02E+02  5.00E+00  7.01E+00  5.96E-01  4.61E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.46E+01 -2.42E+00 -4.43E+01 -2.43E+00  2.01E+01  2.12E+00  9.72E+00  6.72E+00  6.30E+00  0.00E+00  2.60E+01
 
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
 #CPUT: Total CPU Time in Seconds,       29.791
Stop Time:
Wed Sep 29 13:20:39 CDT 2021
