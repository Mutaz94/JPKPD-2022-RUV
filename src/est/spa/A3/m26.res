Wed Sep 29 13:24:24 CDT 2021
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
$DATA ../../../../data/spa/A3/dat26.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -55.8037417547771        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2013E+02  1.8853E+01  1.4589E+02 -1.0626E+02  5.9081E+01  6.3805E+01 -2.6778E+01 -1.2879E+02 -1.2211E+02 -1.1735E+02
            -2.7772E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1256.47023780899        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0486E+00  1.0302E+00  9.4795E-01  1.1070E+00  1.0393E+00  7.6461E-01  9.4753E-01  1.0234E+00  9.5109E-01  9.8357E-01
             4.1714E+00
 PARAMETER:  1.4741E-01  1.2971E-01  4.6547E-02  2.0169E-01  1.3855E-01 -1.6838E-01  4.6102E-02  1.2313E-01  4.9858E-02  8.3432E-02
             1.5283E+00
 GRADIENT:   1.5360E+02 -2.2433E+01 -1.7448E+01 -1.3578E+01  4.6117E+00 -1.1921E+01  9.8077E+00  6.3727E+00  1.7014E+01  1.5781E+01
             5.6932E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1268.47672989883        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0025E+00  1.2452E+00  1.3767E+00  1.0026E+00  1.3602E+00  7.6516E-01  3.9238E-01  1.1751E+00  9.8373E-01  7.4211E-01
             4.2116E+00
 PARAMETER:  1.0251E-01  3.1931E-01  4.1971E-01  1.0264E-01  4.0761E-01 -1.6767E-01 -8.3553E-01  2.6139E-01  8.3600E-02 -1.9826E-01
             1.5378E+00
 GRADIENT:  -4.1202E-01  2.3957E+01 -2.4557E+00  2.7203E+01 -7.0853E+00 -4.9697E+00  2.2341E-01  2.8254E+00  4.6418E+00  3.9002E+00
             4.7838E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1272.10760519789        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.9722E-01  1.0166E+00  1.2300E+00  1.1073E+00  1.1664E+00  7.7823E-01  6.0934E-01  3.0988E-01  8.4202E-01  3.6091E-01
             3.9547E+00
 PARAMETER:  9.7217E-02  1.1645E-01  3.0705E-01  2.0196E-01  2.5395E-01 -1.5073E-01 -3.9538E-01 -1.0716E+00 -7.1955E-02 -9.1913E-01
             1.4749E+00
 GRADIENT:   1.4780E+01 -6.6855E-01  1.6004E-01 -5.6277E-01 -2.2913E+00 -3.3096E-01 -1.7073E-01  3.8288E-01 -9.0614E-01  7.2960E-01
            -3.7884E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1272.62124007675        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      322
 NPARAMETR:  9.9473E-01  1.0564E+00  1.2573E+00  1.0834E+00  1.2077E+00  7.8191E-01  6.1543E-01  5.5039E-02  8.6459E-01  9.8851E-02
             3.9712E+00
 PARAMETER:  9.4717E-02  1.5483E-01  3.2896E-01  1.8014E-01  2.8869E-01 -1.4601E-01 -3.8543E-01 -2.7997E+00 -4.5498E-02 -2.2141E+00
             1.4791E+00
 GRADIENT:  -1.6610E+01 -3.6440E+00 -6.6457E-02 -7.3888E+00 -6.1523E-01 -4.0469E-01 -7.5340E-02  1.0796E-02  3.6392E-01  3.8148E-02
            -1.2830E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1273.00469016479        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      499
 NPARAMETR:  9.9981E-01  8.4471E-01  1.4334E+00  1.2360E+00  1.1859E+00  7.7889E-01  7.6642E-01  1.4322E-02  7.3735E-01  3.6163E-02
             4.0401E+00
 PARAMETER:  9.9813E-02 -6.8757E-02  4.6005E-01  3.1187E-01  2.7047E-01 -1.4988E-01 -1.6603E-01 -4.1460E+00 -2.0469E-01 -3.2197E+00
             1.4963E+00
 GRADIENT:  -1.3196E+00  4.0998E+00  1.3714E+00  5.7590E+00 -2.7878E+00 -1.1651E+00 -1.2278E-01  6.6055E-04 -7.7574E-01  5.4064E-03
            -1.7006E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1273.14278620942        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      675
 NPARAMETR:  9.9877E-01  6.7607E-01  1.3565E+00  1.3345E+00  1.0884E+00  7.8151E-01  8.5152E-01  1.0000E-02  7.1231E-01  1.0000E-02
             4.0388E+00
 PARAMETER:  9.8766E-02 -2.9146E-01  4.0494E-01  3.8853E-01  1.8469E-01 -1.4653E-01 -6.0730E-02 -6.3276E+00 -2.3925E-01 -4.7756E+00
             1.4959E+00
 GRADIENT:   6.2510E-01  6.4923E-01  1.0327E-01  4.6813E-01 -4.6447E-01 -1.2805E-01 -7.0188E-02  0.0000E+00  4.2719E-03  0.0000E+00
             2.0119E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1273.25695298375        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      850
 NPARAMETR:  9.9553E-01  4.4863E-01  1.4036E+00  1.4747E+00  1.0266E+00  7.8087E-01  1.1718E+00  1.0000E-02  6.5429E-01  1.0000E-02
             4.0372E+00
 PARAMETER:  9.5523E-02 -7.0155E-01  4.3907E-01  4.8843E-01  1.2626E-01 -1.4734E-01  2.5852E-01 -1.1571E+01 -3.2421E-01 -8.7300E+00
             1.4956E+00
 GRADIENT:  -2.3687E-01  2.9389E-01  1.4769E-01  2.1182E-01  1.4554E-01 -5.0521E-02  1.0026E-01  0.0000E+00  2.8624E-01  0.0000E+00
             7.1202E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1273.29951464754        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1027
 NPARAMETR:  9.9343E-01  3.1642E-01  1.3390E+00  1.5503E+00  9.5091E-01  7.8126E-01  1.3918E+00  1.0000E-02  6.3643E-01  1.0000E-02
             4.0289E+00
 PARAMETER:  9.3404E-02 -1.0507E+00  3.9196E-01  5.3845E-01  4.9667E-02 -1.4685E-01  4.3063E-01 -1.6873E+01 -3.5189E-01 -1.2698E+01
             1.4935E+00
 GRADIENT:  -1.1912E+00  9.2779E-01  1.0272E+00  4.4298E+00 -2.4493E+00 -9.5843E-04 -8.3186E-02  0.0000E+00 -2.8059E-01  0.0000E+00
            -8.2582E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1273.33294158514        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1206
 NPARAMETR:  9.9175E-01  2.1387E-01  1.4048E+00  1.6156E+00  9.5197E-01  7.7991E-01  1.6813E+00  1.0000E-02  6.1641E-01  1.0000E-02
             4.0311E+00
 PARAMETER:  9.1719E-02 -1.4424E+00  4.3988E-01  5.7974E-01  5.0782E-02 -1.4858E-01  6.1958E-01 -2.3376E+01 -3.8384E-01 -1.7649E+01
             1.4940E+00
 GRADIENT:  -9.7686E-01  4.8234E-01  2.3697E-01  4.2190E+00 -7.4053E-01 -2.2460E-01 -4.4747E-02  0.0000E+00 -2.2136E-01  0.0000E+00
            -4.3000E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1273.34944440039        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1385
 NPARAMETR:  9.9048E-01  1.3667E-01  1.4274E+00  1.6627E+00  9.3952E-01  7.7947E-01  2.1592E+00  1.0000E-02  6.0372E-01  1.0000E-02
             4.0304E+00
 PARAMETER:  9.0431E-02 -1.8902E+00  4.5585E-01  6.0843E-01  3.7619E-02 -1.4915E-01  8.6972E-01 -3.1380E+01 -4.0464E-01 -2.3733E+01
             1.4939E+00
 GRADIENT:  -7.3800E-01  3.0325E-01  1.1209E-01  4.3614E+00 -4.7919E-01 -2.2864E-01 -1.9331E-02  0.0000E+00 -2.1091E-01  0.0000E+00
            -4.8826E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1273.36481565302        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1563
 NPARAMETR:  9.8933E-01  8.2402E-02  1.3844E+00  1.6894E+00  9.0532E-01  7.7990E-01  2.9636E+00  1.0000E-02  5.9859E-01  1.0000E-02
             4.0294E+00
 PARAMETER:  8.9270E-02 -2.3961E+00  4.2530E-01  6.2436E-01  5.2806E-04 -1.4859E-01  1.1864E+00 -4.0914E+01 -4.1319E-01 -3.0966E+01
             1.4936E+00
 GRADIENT:  -5.1875E-01  1.2493E-01  2.8088E-01  1.9628E+00 -8.4404E-01  3.2887E-02 -1.6366E-03  0.0000E+00 -2.3600E-01  0.0000E+00
            -1.3576E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1273.36763065506        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1740
 NPARAMETR:  9.8886E-01  5.9843E-02  1.3895E+00  1.7033E+00  9.0152E-01  7.7970E-01  3.5874E+00  1.0000E-02  5.9813E-01  1.0000E-02
             4.0289E+00
 PARAMETER:  8.8801E-02 -2.7160E+00  4.2892E-01  6.3254E-01 -3.6756E-03 -1.4885E-01  1.3774E+00 -4.7029E+01 -4.1395E-01 -3.5621E+01
             1.4935E+00
 GRADIENT:  -7.3410E-01  1.0247E-01  2.4774E-01  2.4638E+00 -7.0988E-01  2.7705E-03  1.6589E-02  0.0000E+00  2.7732E-01  0.0000E+00
             1.5629E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1273.37057892530        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1918
 NPARAMETR:  9.8854E-01  3.7537E-02  1.4170E+00  1.7184E+00  9.0747E-01  7.7928E-01  4.7728E+00  1.0000E-02  5.9257E-01  1.0000E-02
             4.0292E+00
 PARAMETER:  8.8471E-02 -3.1824E+00  4.4857E-01  6.4139E-01  2.9066E-03 -1.4939E-01  1.6629E+00 -5.6052E+01 -4.2328E-01 -4.2500E+01
             1.4936E+00
 GRADIENT:  -5.1341E-01  6.6973E-02  1.2237E-01  2.1952E+00 -3.7330E-01 -6.6558E-02  2.3796E-02  0.0000E+00  1.1406E-01  0.0000E+00
            -1.7710E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1273.37291316419        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2096
 NPARAMETR:  9.8829E-01  2.1315E-02  1.4292E+00  1.7283E+00  9.0844E-01  7.7912E-01  6.5771E+00  1.0000E-02  5.8938E-01  1.0000E-02
             4.0293E+00
 PARAMETER:  8.8222E-02 -3.7483E+00  4.5710E-01  6.4713E-01  3.9781E-03 -1.4959E-01  1.9836E+00 -6.7099E+01 -4.2868E-01 -5.0912E+01
             1.4936E+00
 GRADIENT:  -1.7559E-01  3.3351E-02  2.3748E-02  1.0615E+00 -7.5283E-02 -5.3598E-02  2.0708E-02  0.0000E+00  1.0770E-02  0.0000E+00
            -6.1823E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1273.37494914032        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2274
 NPARAMETR:  9.8803E-01  1.0131E-02  1.4234E+00  1.7344E+00  9.0284E-01  7.7916E-01  9.7616E+00  1.0000E-02  5.8852E-01  1.0000E-02
             4.0291E+00
 PARAMETER:  8.7962E-02 -4.4921E+00  4.5303E-01  6.5066E-01 -2.2106E-03 -1.4954E-01  2.3785E+00 -8.1683E+01 -4.3015E-01 -6.2007E+01
             1.4936E+00
 GRADIENT:  -3.0684E-01  4.1787E-02  8.6070E-02  1.5269E+00 -2.9557E-01 -3.6427E-02  1.2079E-02  0.0000E+00 -4.3119E-02  0.0000E+00
            -8.2447E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1273.37638414425        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2459
 NPARAMETR:  9.8806E-01  1.0000E-02  1.4216E+00  1.7331E+00  9.0314E-01  7.7915E-01  8.4662E+00  1.0000E-02  5.8897E-01  1.0000E-02
             4.0289E+00
 PARAMETER:  8.7989E-02 -4.5676E+00  4.5180E-01  6.4993E-01 -1.8824E-03 -1.4955E-01  2.2361E+00 -8.1683E+01 -4.2938E-01 -6.2007E+01
             1.4935E+00
 GRADIENT:   3.0613E-01  0.0000E+00 -2.2959E-01 -1.3993E+00  5.1991E-01  1.0980E-02  8.1214E-04  0.0000E+00  4.5636E-02  0.0000E+00
             4.1671E-02

0ITERATION NO.:   83    OBJECTIVE VALUE:  -1273.37662918553        NO. OF FUNC. EVALS.: 110
 CUMULATIVE NO. OF FUNC. EVALS.:     2569
 NPARAMETR:  9.8800E-01  1.0000E-02  1.4207E+00  1.7330E+00  9.0272E-01  7.7911E-01  8.4154E+00  1.0000E-02  5.8898E-01  1.0000E-02
             4.0286E+00
 PARAMETER:  8.7952E-02 -4.5219E+00  4.5248E-01  6.4999E-01 -2.7593E-03 -1.4956E-01  2.2317E+00 -8.1683E+01 -4.2954E-01 -6.2007E+01
             1.4934E+00
 GRADIENT:   3.6093E-02  5.0305E-04  1.1262E-01  2.1307E-01 -9.6065E-02  5.3623E-03  3.1508E-05  0.0000E+00 -1.0634E-02  0.0000E+00
             1.3990E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2569
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.7051E-04 -9.9049E-04  7.9095E-05 -1.5021E-02 -3.4484E-06
 SE:             2.7825E-02  1.2245E-03  7.3871E-05  2.2318E-02  1.3069E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7504E-01  4.1857E-01  2.8430E-01  5.0091E-01  9.7895E-01

 ETASHRINKSD(%)  6.7832E+00  9.5898E+01  9.9753E+01  2.5232E+01  9.9562E+01
 ETASHRINKVR(%)  1.3106E+01  9.9832E+01  9.9999E+01  4.4097E+01  9.9998E+01
 EBVSHRINKSD(%)  6.6494E+00  9.5888E+01  9.9711E+01  2.5323E+01  9.9540E+01
 EBVSHRINKVR(%)  1.2857E+01  9.9831E+01  9.9999E+01  4.4234E+01  9.9998E+01
 RELATIVEINF(%)  4.3980E+01  4.2615E-04  2.8044E-05  1.7306E-01  5.1290E-05
 EPSSHRINKSD(%)  1.7802E+01
 EPSSHRINKVR(%)  3.2434E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1273.3766291855341     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -538.22580262179588     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    32.58
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1273.377       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.88E-01  1.00E-02  1.42E+00  1.73E+00  9.02E-01  7.79E-01  8.43E+00  1.00E-02  5.89E-01  1.00E-02  4.03E+00
 


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
+        1.87E+03
 
 TH 2
+       -4.19E+01  4.47E+01
 
 TH 3
+       -2.95E+00  1.22E+01  4.40E+01
 
 TH 4
+       -1.63E+02  1.53E+02  1.97E+00  5.65E+02
 
 TH 5
+        3.31E+01 -6.43E+01 -1.11E+02 -1.29E+02  3.06E+02
 
 TH 6
+        2.65E+02  2.35E+00  8.17E+00 -5.66E+00 -1.72E+01  2.76E+02
 
 TH 7
+        3.06E-03 -2.69E-03 -2.76E-03 -7.31E-03  8.56E-03  2.94E-03  3.10E-07
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.78E+00 -1.57E+01 -7.04E+00 -5.11E+01  2.90E+01  6.93E+00  1.21E-03  0.00E+00  6.07E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.22E+01 -4.07E+00 -1.16E+00 -1.41E+01  6.25E+00  9.27E+00  4.05E-04  0.00E+00  1.96E+00  0.00E+00  1.10E+00
 
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
+        1.64E+03
 
 TH 2
+       -4.66E+01  2.28E+02
 
 TH 3
+       -2.06E+00  9.31E+00  4.16E+01
 
 TH 4
+       -1.79E+02  1.22E+02  4.23E+00  5.86E+02
 
 TH 5
+        3.26E+01 -4.92E+01 -1.01E+02 -1.36E+02  2.84E+02
 
 TH 6
+       -8.98E-01 -6.19E+00  4.84E+00 -3.15E+01 -2.11E+00  2.48E+02
 
 TH 7
+        3.40E-04  3.83E-02  7.25E-05 -1.09E-02  1.10E-02  2.10E-03  2.90E-04
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.48E+00 -9.42E+00  3.19E+00 -3.09E+01  2.49E+01  4.75E+00  3.36E-02  0.00E+00  1.70E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.38E+01 -2.85E+00  1.75E-01 -1.29E+01  3.50E+00  5.82E+00  3.78E-03  0.00E+00  1.85E+01  0.00E+00  2.85E+01
 
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
+        1.63E+03
 
 TH 2
+       -4.14E+01  2.80E+01
 
 TH 3
+        1.83E+00  7.93E+00  3.41E+01
 
 TH 4
+       -1.65E+02  1.26E+02  3.98E+00  5.98E+02
 
 TH 5
+        2.07E-02 -4.71E+01 -8.80E+01 -1.37E+02  2.62E+02
 
 TH 6
+       -2.77E+02 -7.17E+00  5.38E+00 -4.61E+01  5.16E+00  2.55E+02
 
 TH 7
+       -9.83E-03 -4.79E-03 -2.73E-03 -1.78E-02  1.50E-02  5.30E-03  6.31E-06
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -7.68E+00 -1.68E+01  4.50E-01 -6.79E+01  2.38E+01  9.59E+00  2.57E-02  0.00E+00  1.29E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.02E+02 -1.33E+01 -7.33E+00 -5.32E+01  2.59E+01  1.42E+01  3.50E-03  0.00E+00  2.00E+01  0.00E+00  1.02E+02
 
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
 #CPUT: Total CPU Time in Seconds,       38.732
Stop Time:
Wed Sep 29 13:25:04 CDT 2021
