Thu Sep 30 08:44:27 CDT 2021
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
$DATA ../../../../data/spa2/D/dat22.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m22.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   1489.29975439640        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4340E+02 -4.7236E+00 -4.8897E+01 -1.7959E+01  3.8999E+02 -1.2696E+03 -3.6452E+02 -4.5636E+01 -5.8335E+02 -5.7953E+02
            -5.0491E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1390.81536969049        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:      106
 NPARAMETR:  1.0333E+00  8.9317E-01  9.3469E-01  1.2422E+00  7.8133E-01  1.9982E+00  1.6247E+00  1.0645E+00  2.6155E+00  2.3229E+00
             2.9745E+00
 PARAMETER:  1.3276E-01 -1.2980E-02  3.2457E-02  3.1691E-01 -1.4676E-01  7.9224E-01  5.8530E-01  1.6249E-01  1.0614E+00  9.4283E-01
             1.1901E+00
 GRADIENT:  -6.1801E+01 -7.6808E+00 -3.2534E+01  3.1529E+01  1.3857E+01 -2.0562E+02  8.6227E+00  9.5634E+00  4.5002E+01 -6.9242E+00
            -7.7482E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1409.80563458428        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      285
 NPARAMETR:  1.3383E+00  8.0856E-01  1.7834E+00  1.3545E+00  1.1568E+00  2.5924E+00  1.1424E+00  9.5377E-01  2.6350E+00  2.7451E+00
             2.6367E+00
 PARAMETER:  3.9141E-01 -1.1250E-01  6.7854E-01  4.0345E-01  2.4565E-01  1.0526E+00  2.3311E-01  5.2672E-02  1.0689E+00  1.1098E+00
             1.0695E+00
 GRADIENT:   5.4486E+01 -8.1419E+00 -5.7053E+01  3.2326E+01  3.7782E+01 -4.5365E+01  1.3644E+00  3.6273E+00  1.6640E+01  2.9722E+01
            -2.6897E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1436.66138294678        NO. OF FUNC. EVALS.: 109
 CUMULATIVE NO. OF FUNC. EVALS.:      394
 NPARAMETR:  1.2436E+00  8.2761E-01  3.0296E+00  1.2882E+00  1.1056E+00  2.6712E+00  9.7744E-01  7.4666E-01  2.6314E+00  2.7411E+00
             2.6418E+00
 PARAMETER:  3.1803E-01 -8.9214E-02  1.2084E+00  3.5327E-01  2.0038E-01  1.0825E+00  7.7181E-02 -1.9214E-01  1.0675E+00  1.1084E+00
             1.0715E+00
 GRADIENT:   2.2534E+02  2.2974E+01 -1.4639E+01  7.0926E+01 -2.4672E+01  3.4082E+02  1.2798E+00  1.4692E+00  8.6464E+01  8.3252E+01
            -2.4182E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1436.94496588677        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      563
 NPARAMETR:  1.2419E+00  8.2632E-01  3.0400E+00  1.2871E+00  1.1084E+00  2.6725E+00  9.7222E-01  7.4179E-01  2.6312E+00  2.7411E+00
             2.6421E+00
 PARAMETER:  3.1664E-01 -9.0773E-02  1.2119E+00  3.5243E-01  2.0294E-01  1.0830E+00  7.1832E-02 -1.9869E-01  1.0674E+00  1.1083E+00
             1.0716E+00
 GRADIENT:   2.7825E+01  1.7146E+01 -1.9624E+01  1.3988E+01 -3.0544E+01 -2.1362E+01  9.2081E-01  1.4036E+00  1.7216E-01  3.1859E+01
            -2.5713E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1441.65258934452        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      664
 NPARAMETR:  1.0529E+00  6.6314E-01  3.0313E+00  1.2857E+00  1.1078E+00  2.6528E+00  4.0004E-01  7.4219E-01  2.6203E+00  2.7457E+00
             2.6717E+00
 PARAMETER:  1.5153E-01 -3.1076E-01  1.2090E+00  3.5131E-01  2.0240E-01  1.0756E+00 -8.1619E-01 -1.9814E-01  1.0633E+00  1.1100E+00
             1.0827E+00
 GRADIENT:   6.6448E+01  6.4055E+00 -1.7708E+01  5.1882E+01 -1.3362E+01  3.6614E+02  1.2047E+00  1.6322E+00  1.0461E+02  8.7071E+01
            -2.1703E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1442.14120948897        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      735
 NPARAMETR:  9.6668E-01  6.1203E-01  3.0355E+00  1.2849E+00  1.1078E+00  2.5931E+00  1.4264E-01  7.4219E-01  2.6031E+00  2.7303E+00
             2.7347E+00
 PARAMETER:  6.6116E-02 -3.9097E-01  1.2104E+00  3.5071E-01  2.0242E-01  1.0529E+00 -1.8475E+00 -1.9815E-01  1.0567E+00  1.1044E+00
             1.1060E+00
 GRADIENT:   2.8188E+00 -4.8354E-01 -1.8227E+01  4.2089E+01 -1.0020E+01  3.5365E+02  2.2891E-01  1.7189E+00  1.0769E+02  8.4589E+01
            -1.7570E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1459.55644337864        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      896
 NPARAMETR:  1.1425E+00  7.3342E-01  3.0634E+00  1.2839E+00  1.1081E+00  2.4945E+00  1.0000E-02  7.4218E-01  2.5600E+00  2.6749E+00
             3.1067E+00
 PARAMETER:  2.3322E-01 -2.1004E-01  1.2195E+00  3.4989E-01  2.0264E-01  1.0141E+00 -8.2897E+00 -1.9817E-01  1.0400E+00  1.0839E+00
             1.2336E+00
 GRADIENT:  -2.3737E+00  7.9021E+00 -1.8493E+01 -2.1821E+00 -2.4359E+01 -5.1379E+01  0.0000E+00  1.6220E+00 -2.9163E+00  2.7881E+01
            -1.8521E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1460.65048375792        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1056
 NPARAMETR:  1.1732E+00  7.0964E-01  3.0786E+00  1.3019E+00  1.1074E+00  2.5048E+00  1.0000E-02  3.8093E-01  2.5548E+00  2.6837E+00
             3.1187E+00
 PARAMETER:  2.5972E-01 -2.4299E-01  1.2245E+00  3.6383E-01  2.0203E-01  1.0182E+00 -8.3086E+00 -8.6514E-01  1.0380E+00  1.0872E+00
             1.2374E+00
 GRADIENT:   1.2374E+02  1.1836E+01 -1.4340E+01  4.3882E+01 -1.8199E+01  1.9977E+02  0.0000E+00  4.8139E-01  6.4006E+01  6.5274E+01
             1.9130E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1463.11030905606        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1238
 NPARAMETR:  1.1731E+00  7.0918E-01  3.1582E+00  1.3540E+00  1.1084E+00  2.5691E+00  1.0000E-02  1.0000E-02  2.5499E+00  2.6315E+00
             3.1683E+00
 PARAMETER:  2.5964E-01 -2.4365E-01  1.2500E+00  4.0307E-01  2.0292E-01  1.0435E+00 -8.3086E+00 -6.6382E+00  1.0361E+00  1.0676E+00
             1.2532E+00
 GRADIENT:   6.4351E+00  1.0342E+01 -1.6638E+01  8.9591E+00 -2.5223E+01 -3.7055E+01  0.0000E+00  0.0000E+00  2.0989E+00  2.3624E+01
             5.7678E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1463.65098886018        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1429
 NPARAMETR:  1.1810E+00  6.5566E-01  3.1609E+00  1.3526E+00  1.1080E+00  2.5781E+00  1.0000E-02  1.0000E-02  2.5374E+00  2.6284E+00
             3.1724E+00
 PARAMETER:  2.6633E-01 -3.2211E-01  1.2509E+00  4.0201E-01  2.0255E-01  1.0471E+00 -8.3086E+00 -6.7995E+00  1.0312E+00  1.0664E+00
             1.2545E+00
 GRADIENT:   8.8500E+00 -4.7248E-01 -1.7135E+01 -2.2804E+00 -2.1570E+01 -3.5749E+01  0.0000E+00  0.0000E+00  1.1607E+01  2.5089E+01
             1.0426E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1464.52991634338        NO. OF FUNC. EVALS.: 202
 CUMULATIVE NO. OF FUNC. EVALS.:     1631             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1940E+00  6.5725E-01  3.1673E+00  1.3528E+00  1.1161E+00  2.6119E+00  1.0000E-02  1.0000E-02  2.5302E+00  2.6061E+00
             3.1680E+00
 PARAMETER:  2.7727E-01 -3.1969E-01  1.2529E+00  4.0220E-01  2.0982E-01  1.0601E+00 -8.3086E+00 -6.7995E+00  1.0283E+00  1.0579E+00
             1.2531E+00
 GRADIENT:   1.3314E+02  8.3476E+00 -1.4018E+01  5.1048E+01 -1.5257E+01  2.2160E+02  0.0000E+00  0.0000E+00  6.9796E+01  5.7712E+01
             2.3582E+01

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1464.52991634338        NO. OF FUNC. EVALS.:  69
 CUMULATIVE NO. OF FUNC. EVALS.:     1700
 NPARAMETR:  1.1906E+00  6.5933E-01  3.1695E+00  1.3582E+00  1.1138E+00  2.5872E+00  1.0000E-02  1.0000E-02  2.5318E+00  2.6327E+00
             3.1662E+00
 PARAMETER:  2.7727E-01 -3.1969E-01  1.2529E+00  4.0220E-01  2.0982E-01  1.0601E+00 -8.3086E+00 -6.7995E+00  1.0283E+00  1.0579E+00
             1.2531E+00
 GRADIENT:   8.3772E+02 -1.4323E+03 -3.8640E+01 -5.7384E+02  2.1543E+03  4.0548E+02  0.0000E+00  0.0000E+00 -5.5047E+01 -7.9906E+02
             3.3954E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1700
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.2462E-02 -3.7775E-04 -8.7828E-06 -4.2890E-03  1.1953E-03
 SE:             3.2538E-02  8.7085E-05  3.8188E-05  2.8368E-02  2.4398E-02
 N:                     100         100         100         100         100

 P VAL.:         7.0171E-01  1.4411E-05  8.1810E-01  8.7983E-01  9.6092E-01

 ETASHRINKSD(%)  1.0000E-10  9.9708E+01  9.9872E+01  4.9634E+00  1.8264E+01
 ETASHRINKVR(%)  1.0000E-10  9.9999E+01  1.0000E+02  9.6804E+00  3.3192E+01
 EBVSHRINKSD(%)  4.1646E-01  9.9786E+01  9.9824E+01  2.6738E+00  8.0341E+00
 EBVSHRINKVR(%)  8.3118E-01  1.0000E+02  1.0000E+02  5.2762E+00  1.5423E+01
 RELATIVEINF(%)  9.9144E+01  1.6497E-04  1.3873E-04  3.6284E+01  3.8965E+01
 EPSSHRINKSD(%)  2.4585E+01
 EPSSHRINKVR(%)  4.3126E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1464.5299163433824     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -361.80367649777531     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    33.46
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.16
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1464.530       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.19E+00  6.57E-01  3.17E+00  1.35E+00  1.12E+00  2.61E+00  1.00E-02  1.00E-02  2.53E+00  2.61E+00  3.17E+00
 


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
+        1.05E+05
 
 TH 2
+        1.46E+02  2.59E+05
 
 TH 3
+        2.91E+02 -4.48E+02  1.39E+03
 
 TH 4
+       -6.78E+02 -1.54E+01 -2.03E+02  7.41E+04
 
 TH 5
+       -1.25E+03  1.94E+03  3.91E+02  7.59E+02  2.07E+05
 
 TH 6
+        7.18E+02 -1.14E+03 -5.63E+01  6.46E+03  1.02E+03  2.98E+03
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.30E+03 -4.73E+01  9.75E+02  7.59E+02  1.53E+02 -8.60E+01  0.00E+00  0.00E+00  3.39E+03
 
 TH10
+        2.14E+02 -3.60E+01 -2.48E+01 -7.63E+03  1.55E+02 -9.44E+01  0.00E+00  0.00E+00  1.73E+02  2.86E+03
 
 TH11
+        1.14E+01 -3.84E+00 -6.59E+02 -3.71E+01 -9.90E+01  5.70E+01  0.00E+00  0.00E+00 -1.44E+02  1.05E+03  7.47E+02
 
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
 #CPUT: Total CPU Time in Seconds,       42.693
Stop Time:
Thu Sep 30 08:45:11 CDT 2021
