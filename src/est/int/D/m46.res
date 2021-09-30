Wed Sep 29 08:57:18 CDT 2021
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
$DATA ../../../../data/int/D/dat46.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m46.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25357.5540947270        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.8008E+02  4.4372E+02 -9.9568E+00  3.0181E+02  1.8385E+02 -1.8825E+03 -9.9774E+02 -7.1540E+01 -1.4807E+03 -5.5705E+02
            -5.2742E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -838.674870682893        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  8.6411E-01  1.8053E+00  8.8833E-01  2.0052E+00  9.4597E-01  4.0894E+00  4.8647E+00  9.9105E-01  2.8715E+00  1.6728E+00
             1.2470E+01
 PARAMETER: -4.6057E-02  6.9074E-01 -1.8414E-02  7.9572E-01  4.4457E-02  1.5084E+00  1.6820E+00  9.1006E-02  1.1548E+00  6.1447E-01
             2.6233E+00
 GRADIENT:  -4.8424E+01  2.7133E+01 -4.2988E+01  9.6740E+01 -9.7370E+00  1.6901E+02  9.6631E+01  4.2710E+00  6.8430E+01  3.8690E+01
             5.6037E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -939.000623565038        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  8.6094E-01  1.0543E+00  4.8271E+01  2.9737E+00  2.6398E+00  3.6764E+00  7.5178E+00  8.1834E-01  3.2272E+00  2.2176E+00
             1.1804E+01
 PARAMETER: -4.9725E-02  1.5283E-01  3.9768E+00  1.1898E+00  1.0707E+00  1.4019E+00  2.1173E+00 -1.0048E-01  1.2716E+00  8.9642E-01
             2.5684E+00
 GRADIENT:  -4.9158E+01  2.1227E+01 -2.9654E+00  9.4169E+01  9.7026E+00  1.6405E+02  6.9580E+01  3.7349E-02  6.5109E+01  6.5647E+01
             5.3809E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1117.69618327246        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0456E+00  1.0888E+00  2.7544E+01  1.1574E+00  2.5666E+00  1.7354E+00  3.3235E+00  2.9923E+00  1.2607E+00  9.1542E-01
             9.8383E+00
 PARAMETER:  1.4462E-01  1.8510E-01  3.4158E+00  2.4621E-01  1.0426E+00  6.5123E-01  1.3010E+00  1.1960E+00  3.3169E-01  1.1625E-02
             2.3863E+00
 GRADIENT:  -3.4768E+01 -5.3992E+01 -1.3311E+00 -2.4935E+01  4.2242E+01 -3.1255E+01 -6.0545E+01  2.0665E-01  9.3609E+00  1.3389E+01
             3.9856E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1191.67740722123        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0528E+00  1.7906E+00  7.9772E+00  7.1198E-01  2.1580E+00  2.0359E+00  2.9356E+00  4.6406E+00  1.0725E+00  5.1604E-01
             7.8398E+00
 PARAMETER:  1.5142E-01  6.8252E-01  2.1766E+00 -2.3970E-01  8.6920E-01  8.1096E-01  1.1769E+00  1.6348E+00  1.7003E-01 -5.6157E-01
             2.1592E+00
 GRADIENT:   6.1627E+00  2.3305E+00  2.3422E+00 -3.9089E+00 -1.5067E+00  1.0577E+01 -1.1671E+01  1.4806E-01  4.4777E+00  5.5143E+00
             5.8574E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1194.52253471801        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0432E+00  1.6304E+00  2.3213E+00  7.3219E-01  1.8506E+00  1.9663E+00  3.1154E+00  2.1460E+00  8.2419E-01  2.1296E-01
             7.4464E+00
 PARAMETER:  1.4230E-01  5.8881E-01  9.4214E-01 -2.1171E-01  7.1551E-01  7.7615E-01  1.2364E+00  8.6358E-01 -9.3353E-02 -1.4467E+00
             2.1077E+00
 GRADIENT:   7.8367E+00 -2.9912E+00 -6.0375E+00  2.9837E+00  1.6316E+01 -5.0953E+00  9.5882E+00  3.5470E+00  2.9622E+00  1.0843E+00
            -4.9105E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1199.14711448769        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      470
 NPARAMETR:  1.0394E+00  1.6331E+00  2.1748E+00  7.5193E-01  1.7820E+00  2.0043E+00  3.1136E+00  1.0244E+00  5.5329E-01  1.2295E-01
             7.7438E+00
 PARAMETER:  1.3864E-01  5.9050E-01  8.7696E-01 -1.8511E-01  6.7773E-01  7.9530E-01  1.2358E+00  1.2415E-01 -4.9188E-01 -1.9959E+00
             2.1469E+00
 GRADIENT:  -1.0120E+01 -1.5485E+01 -9.5718E-01 -1.6829E+00 -6.3142E+00 -2.9231E+01 -5.1547E+01  4.3291E-01  1.0785E+00  2.9103E-01
            -1.3561E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1208.04782261823        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      647
 NPARAMETR:  1.0596E+00  1.4518E+00  3.7426E+00  8.9126E-01  1.9250E+00  2.1749E+00  4.1122E+00  5.5776E-01  3.9785E-01  8.4226E-02
             7.8067E+00
 PARAMETER:  1.5792E-01  4.7280E-01  1.4198E+00 -1.5122E-02  7.5491E-01  8.7698E-01  1.5140E+00 -4.8383E-01 -8.2168E-01 -2.3742E+00
             2.1550E+00
 GRADIENT:  -1.0713E-01  9.6771E-02 -2.7237E-02  7.1482E-01 -1.0391E-01 -2.7701E-01 -1.3408E-01  1.3064E-01 -2.4274E-01  1.1089E-01
            -7.4492E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1208.07234061468        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      824
 NPARAMETR:  1.0607E+00  1.4750E+00  3.6019E+00  8.7715E-01  1.9212E+00  2.1802E+00  4.0822E+00  3.7329E-01  3.4873E-01  6.3488E-02
             7.8084E+00
 PARAMETER:  1.5891E-01  4.8868E-01  1.3814E+00 -3.1072E-02  7.5297E-01  8.7941E-01  1.5066E+00 -8.8540E-01 -9.5346E-01 -2.6569E+00
             2.1552E+00
 GRADIENT:   3.8172E-01  9.9844E-02 -1.4524E-01  1.9923E-01 -2.6611E-02  5.1127E-01  3.0191E-01  6.8894E-02 -2.8452E-01  6.3095E-02
            -3.0655E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1208.07441520540        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1000
 NPARAMETR:  1.0594E+00  1.4815E+00  3.6474E+00  8.7448E-01  1.9287E+00  2.1757E+00  4.0650E+00  3.0990E-01  3.4295E-01  5.6653E-02
             7.8059E+00
 PARAMETER:  1.5770E-01  4.9306E-01  1.3940E+00 -3.4124E-02  7.5683E-01  8.7735E-01  1.5024E+00 -1.0715E+00 -9.7017E-01 -2.7708E+00
             2.1549E+00
 GRADIENT:  -1.0683E-01  1.7598E-01  1.1523E-01  3.1342E-01 -2.2576E-01 -2.5411E-01 -4.7072E-01  4.6026E-02 -2.8586E-01  4.9615E-02
            -9.4141E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1208.18168041786        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1182
 NPARAMETR:  1.0611E+00  1.4377E+00  3.6740E+00  8.9852E-01  1.9146E+00  2.1806E+00  4.1381E+00  1.3101E-02  4.2384E-01  1.0000E-02
             7.8128E+00
 PARAMETER:  1.5926E-01  4.6302E-01  1.4013E+00 -7.0057E-03  7.4950E-01  8.7960E-01  1.5202E+00 -4.2351E+00 -7.5841E-01 -4.5236E+00
             2.1558E+00
 GRADIENT:   4.1563E-01  5.1766E-03 -2.8791E-01  5.8946E-01  2.8097E-01  6.9783E-01  4.7492E-01  7.0490E-05 -3.6991E-02  4.9997E-04
             6.0561E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1208.18964892915        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1363
 NPARAMETR:  1.0599E+00  1.4355E+00  3.7005E+00  8.9934E-01  1.9161E+00  2.1903E+00  4.1545E+00  1.0000E-02  4.2425E-01  1.0000E-02
             7.8075E+00
 PARAMETER:  1.5821E-01  4.6151E-01  1.4085E+00 -6.0938E-03  7.5031E-01  8.8402E-01  1.5242E+00 -4.7000E+00 -7.5744E-01 -4.5706E+00
             2.1551E+00
 GRADIENT:   9.8679E-02  1.1082E-01 -1.9246E-01  1.4309E-01  2.1227E-02  2.1646E+00  1.2454E+00  0.0000E+00 -3.1205E-02  0.0000E+00
            -2.7158E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1208.19118041032        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1557             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0598E+00  1.4321E+00  3.7267E+00  9.0011E-01  1.9172E+00  2.1905E+00  4.1623E+00  1.0000E-02  4.2956E-01  1.0000E-02
             7.8077E+00
 PARAMETER:  1.5811E-01  4.5911E-01  1.4155E+00 -5.2345E-03  7.5088E-01  8.8411E-01  1.5261E+00 -4.7000E+00 -7.4499E-01 -4.5706E+00
             2.1551E+00
 GRADIENT:   1.2688E+01  1.0287E+01  3.5451E-01  1.2605E+00  4.0304E+00  3.9762E+01  7.0394E+01  0.0000E+00  2.4560E-01  0.0000E+00
             3.8910E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1208.19183389367        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1744
 NPARAMETR:  1.0598E+00  1.4307E+00  3.7377E+00  9.0099E-01  1.9180E+00  2.1905E+00  4.1647E+00  1.0000E-02  4.3206E-01  1.0000E-02
             7.8077E+00
 PARAMETER:  1.5812E-01  4.5817E-01  1.4185E+00 -4.2614E-03  7.5129E-01  8.8411E-01  1.5267E+00 -4.7000E+00 -7.3920E-01 -4.5706E+00
             2.1551E+00
 GRADIENT:   7.1808E-02 -1.5591E-02 -6.9952E-02 -6.0192E-01 -2.0194E-01  2.1895E+00  1.6194E+00  0.0000E+00  3.6573E-02  0.0000E+00
             8.0311E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1208.19270199671        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1938             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0598E+00  1.4292E+00  3.7511E+00  9.0258E-01  1.9190E+00  2.1905E+00  4.1673E+00  1.0000E-02  4.3116E-01  1.0000E-02
             7.8074E+00
 PARAMETER:  1.5812E-01  4.5712E-01  1.4221E+00 -2.4936E-03  7.5180E-01  8.8411E-01  1.5273E+00 -4.7000E+00 -7.4127E-01 -4.5706E+00
             2.1551E+00
 GRADIENT:   1.2678E+01  1.0321E+01  3.5072E-01  2.0024E+00  4.1670E+00  3.9752E+01  7.0267E+01  0.0000E+00  1.9308E-01  0.0000E+00
             3.8411E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1208.19307541527        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2125
 NPARAMETR:  1.0598E+00  1.4279E+00  3.7606E+00  9.0330E-01  1.9196E+00  2.1905E+00  4.1696E+00  1.0000E-02  4.3247E-01  1.0000E-02
             7.8074E+00
 PARAMETER:  1.5813E-01  4.5622E-01  1.4246E+00 -1.7048E-03  7.5212E-01  8.8411E-01  1.5278E+00 -4.7000E+00 -7.3825E-01 -4.5706E+00
             2.1551E+00
 GRADIENT:   6.7550E-02  8.8518E-02 -6.6512E-02  1.3649E-01 -1.0597E-01  2.1904E+00  1.3838E+00  0.0000E+00 -2.7285E-02  0.0000E+00
            -4.0985E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1208.19333007086        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2319             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0599E+00  1.4256E+00  3.7734E+00  9.0361E-01  1.9203E+00  2.1905E+00  4.1743E+00  1.0000E-02  4.3679E-01  1.0000E-02
             7.8083E+00
 PARAMETER:  1.5814E-01  4.5460E-01  1.4280E+00 -1.3616E-03  7.5247E-01  8.8411E-01  1.5290E+00 -4.7000E+00 -7.2830E-01 -4.5706E+00
             2.1552E+00
 GRADIENT:   1.2684E+01  1.0104E+01  3.9142E-01  1.3948E+00  4.1625E+00  3.9746E+01  7.0781E+01  0.0000E+00  2.5506E-01  0.0000E+00
             3.8917E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1208.19352872971        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2504
 NPARAMETR:  1.0599E+00  1.4248E+00  3.7789E+00  9.0411E-01  1.9206E+00  2.1905E+00  4.1758E+00  1.0000E-02  4.3773E-01  1.0000E-02
             7.8083E+00
 PARAMETER:  1.5814E-01  4.5404E-01  1.4294E+00 -8.0052E-04  7.5265E-01  8.8411E-01  1.5293E+00 -4.7000E+00 -7.2616E-01 -4.5706E+00
             2.1552E+00
 GRADIENT:   7.0082E-02 -4.7238E-02 -3.2249E-02 -4.6236E-01 -1.0997E-01  2.1895E+00  1.6764E+00  0.0000E+00  3.2421E-02  0.0000E+00
             9.9564E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1208.19368112377        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2698             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0599E+00  1.4245E+00  3.7857E+00  9.0538E-01  1.9211E+00  2.1905E+00  4.1757E+00  1.0000E-02  4.3592E-01  1.0000E-02
             7.8075E+00
 PARAMETER:  1.5814E-01  4.5382E-01  1.4312E+00  6.0104E-04  7.5289E-01  8.8411E-01  1.5293E+00 -4.7000E+00 -7.3029E-01 -4.5706E+00
             2.1551E+00
 GRADIENT:   1.2679E+01  1.0224E+01  3.8352E-01  2.2981E+00  4.2255E+00  3.9740E+01  7.0445E+01  0.0000E+00  1.8070E-01  0.0000E+00
             3.8234E+01

0ITERATION NO.:   94    OBJECTIVE VALUE:  -1208.19375689896        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     2838
 NPARAMETR:  1.0599E+00  1.4240E+00  3.7899E+00  9.0574E-01  1.9213E+00  2.1905E+00  4.1767E+00  1.0000E-02  4.3639E-01  1.0000E-02
             7.8074E+00
 PARAMETER:  1.5814E-01  4.5314E-01  1.4318E+00 -6.2917E-06  7.5294E-01  8.8411E-01  1.5298E+00 -4.7000E+00 -7.2298E-01 -4.5706E+00
             2.1552E+00
 GRADIENT:   4.5583E-04 -3.7571E-02 -2.3495E-02 -2.8180E-01 -6.5021E-02 -4.7925E-04  1.0060E-01  0.0000E+00  1.9542E-02  0.0000E+00
             2.6910E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2838
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.9216E-03  1.0654E-02 -2.6603E-05 -3.1996E-02  2.8652E-05
 SE:             2.8680E-02  2.7221E-02  3.7981E-05  7.3914E-03  1.5671E-04
 N:                     100         100         100         100         100

 P VAL.:         7.2939E-01  6.9552E-01  4.8366E-01  1.5000E-05  8.5493E-01

 ETASHRINKSD(%)  3.9185E+00  8.8051E+00  9.9873E+01  7.5238E+01  9.9475E+01
 ETASHRINKVR(%)  7.6835E+00  1.6835E+01  1.0000E+02  9.3868E+01  9.9997E+01
 EBVSHRINKSD(%)  3.2316E+00  5.1631E+00  9.9868E+01  8.1265E+01  9.9431E+01
 EBVSHRINKVR(%)  6.3588E+00  1.0060E+01  1.0000E+02  9.6490E+01  9.9997E+01
 RELATIVEINF(%)  9.3510E+01  4.5130E+01  3.4893E-05  1.5966E+00  6.8933E-04
 EPSSHRINKSD(%)  6.7246E+00
 EPSSHRINKVR(%)  1.2997E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1208.1937568989617     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       445.89560286944902     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    87.68
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    14.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1208.194       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.42E+00  3.79E+00  9.05E-01  1.92E+00  2.19E+00  4.18E+00  1.00E-02  4.39E-01  1.00E-02  7.81E+00
 


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
+        1.63E+02
 
 TH 2
+       -2.96E+00  7.16E+00
 
 TH 3
+        7.58E-01  7.88E-01  1.38E+00
 
 TH 4
+       -5.62E+01  3.89E+01 -6.86E+00  3.14E+02
 
 TH 5
+       -3.74E+00 -9.25E+00 -1.19E+01  4.32E+01  1.03E+02
 
 TH 6
+        1.29E+01 -1.14E+00 -4.90E-01 -5.45E+00  4.66E+00  1.31E+00
 
 TH 7
+        6.37E+00 -3.90E+00 -2.30E-01 -2.38E+01  3.40E+00  9.24E-01  2.28E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.54E+00 -6.64E+00 -2.87E-01 -4.09E+01  4.89E+00  1.35E+00  3.84E+00  0.00E+00  6.49E+00
 
 TH10
+        2.29E-62 -1.45E-62 -1.06E-62 -5.56E-63  9.48E-62  6.75E-63  7.18E-63  0.00E+00  1.11E-62  3.2E-106
 
 TH11
+       -2.40E-01 -2.33E+00  2.48E-02 -1.48E+01  6.01E-01  2.88E-01  1.27E+00  0.00E+00  2.20E+00  8.17E-45  1.26E+00
 
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
+        1.88E+02
 
 TH 2
+       -9.38E-01  2.83E+01
 
 TH 3
+        3.79E-01  8.81E-01  1.52E+00
 
 TH 4
+       -6.49E+00  4.27E+01 -5.94E+00  3.21E+02
 
 TH 5
+       -2.66E+00 -9.48E+00 -1.02E+01  3.38E+01  9.05E+01
 
 TH 6
+        1.78E+00 -1.11E-01  9.26E-02  1.18E+00 -1.63E+00  3.53E+01
 
 TH 7
+        4.13E-01  2.28E+00 -4.62E-01 -2.44E+01  2.45E+00 -3.51E-01  8.10E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -2.03E-02 -2.34E+00 -5.10E-01 -4.18E+01  4.37E+00 -1.11E-01  3.16E+00  0.00E+00  1.51E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  2.72E+01
 
 TH11
+       -6.71E+00 -3.26E+00 -9.77E-02 -1.83E+01  7.87E-01  1.72E+00  1.40E+00  0.00E+00  2.69E+00  0.00E+00  1.81E+01
 
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
+        1.93E+02
 
 TH 2
+        5.48E+01  2.84E+01
 
 TH 3
+        1.69E+00  9.72E-01  1.22E+00
 
 TH 4
+        8.12E+01  4.76E+01 -4.70E+00  3.53E+02
 
 TH 5
+       -1.66E+01 -5.46E+00 -8.96E+00  2.81E+01  7.53E+01
 
 TH 6
+        3.36E+01  7.42E+00 -1.14E-01 -1.58E+01  9.34E-01  4.03E+01
 
 TH 7
+        8.62E-01  3.66E+00 -5.26E-01 -2.95E+01  8.77E+00  6.45E+00  1.22E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -3.43E+00 -2.99E+00  4.87E-02 -4.52E+01  1.05E+00  1.68E+00  3.51E+00  0.00E+00  1.42E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.29E+02 -3.80E+01 -7.50E-01 -1.05E+02  1.20E+01  1.99E+01  1.05E+01  0.00E+00 -3.37E+00  0.00E+00  6.90E+02
 
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
 #CPUT: Total CPU Time in Seconds,      102.431
Stop Time:
Wed Sep 29 08:59:02 CDT 2021
