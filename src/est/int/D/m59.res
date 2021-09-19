Sat Sep 18 07:12:41 CDT 2021
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
$DATA ../../../../data/int/D/dat59.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m59.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   22444.0388734648        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.5122E+02  2.8460E+02 -3.2624E+01  1.2663E+02  3.6850E+02 -3.2619E+03 -1.1871E+03 -9.8674E+01 -2.0461E+03 -1.0895E+03
            -4.5329E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1143.75382074277        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  3.8155E+00  2.1396E+00  9.7632E-01  1.9557E+00  1.3272E+00  6.1608E+00  8.8451E+00  1.0411E+00  5.5732E+00  3.2057E+00
             8.8649E+00
 PARAMETER:  1.4391E+00  8.6062E-01  7.6038E-02  7.7076E-01  3.8304E-01  1.9182E+00  2.2799E+00  1.4032E-01  1.8180E+00  1.2649E+00
             2.2821E+00
 GRADIENT:   6.7842E+01  4.2055E+00 -4.6753E+01  5.6714E+01 -2.0831E+01  1.0669E+02  1.3871E+02  2.5079E+00  8.3735E+01  7.4044E+01
             3.2887E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1236.87678174841        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  2.7713E+00  6.8284E+00  5.4628E+01  4.0277E-01  2.7573E+00  6.3384E+00  7.7207E+00  9.8570E-01  5.6110E+01  2.6031E+00
             7.5228E+00
 PARAMETER:  1.1193E+00  2.0211E+00  4.1005E+00 -8.0939E-01  1.1143E+00  1.9466E+00  2.1439E+00  8.5592E-02  4.1273E+00  1.0567E+00
             2.1179E+00
 GRADIENT:   4.9895E+01  3.7362E+01  6.7385E+00  7.4957E+00 -1.1749E+01  1.3673E+02  1.4351E+02  1.7906E-01  5.4703E+01  1.0174E+02
             1.5896E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1339.47488820094        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  1.4491E+00  6.3583E-01  1.6201E+02  3.1711E+00  1.8540E+00  2.5767E+00  5.0394E+00  3.2536E+00  5.2004E+00  2.6372E+00
             9.5137E+00
 PARAMETER:  4.7094E-01 -3.5282E-01  5.1876E+00  1.2541E+00  7.1733E-01  1.0465E+00  1.7173E+00  1.2798E+00  1.7487E+00  1.0697E+00
             2.3527E+00
 GRADIENT:   1.3593E+01 -5.3882E+00 -4.1586E-01  3.7497E+01 -5.4467E+01 -3.8963E+01  1.2454E+01  2.2993E+00  7.0281E+01  8.0316E+01
             3.9560E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1457.55180515947        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.3076E+00  1.6668E-01  3.0763E+02  2.0197E+00  2.1204E+00  4.0775E+00  2.9012E+00  1.0000E-02  3.4435E+00  1.3936E+00
             6.9384E+00
 PARAMETER:  3.6816E-01 -1.6917E+00  5.8289E+00  8.0296E-01  8.5159E-01  1.5055E+00  1.1651E+00 -6.6772E+00  1.3365E+00  4.3190E-01
             2.0371E+00
 GRADIENT:  -2.6737E+00 -4.0148E+00 -1.4042E+00 -1.4775E+01 -2.4353E+01  6.2910E+01  8.8843E-01  0.0000E+00 -1.6784E+01  2.3910E+01
             1.9585E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1467.52623061595        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.3555E+00  1.3856E-01  6.2450E+02  2.2369E+00  2.3004E+00  3.4504E+00  3.0497E+00  1.0000E-02  3.5200E+00  9.9705E-01
             6.9207E+00
 PARAMETER:  4.0420E-01 -1.8764E+00  6.5369E+00  9.0509E-01  9.3309E-01  1.3385E+00  1.2150E+00 -6.3143E+00  1.3585E+00  9.7045E-02
             2.0345E+00
 GRADIENT:   8.4330E-01 -2.9325E+00 -6.3988E-01  1.7552E+00  8.1063E+00  1.0595E+00  7.1835E-01  0.0000E+00 -1.7106E+00  5.6613E+00
            -1.1379E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1473.49416041297        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      449
 NPARAMETR:  1.3390E+00  4.6019E-01  2.4867E+03  2.0050E+00  2.4599E+00  3.5197E+00  2.7309E+00  9.1118E-02  3.7253E+00  1.2744E-01
             7.1192E+00
 PARAMETER:  3.9195E-01 -6.7612E-01  7.9187E+00  7.9565E-01  1.0001E+00  1.3584E+00  1.1046E+00 -2.2956E+00  1.4152E+00 -1.9601E+00
             2.0628E+00
 GRADIENT:  -2.3511E+00 -5.0541E+00 -4.6102E-02  3.3604E+00  1.9463E+01  9.4746E+00  3.2371E+00  9.7059E-05 -1.1896E+01 -1.2187E-01
             4.3655E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1492.49328939056        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      520
 NPARAMETR:  1.3365E+00  1.6099E+00  1.6318E+03  1.0132E+00  2.4394E+00  3.1539E+00  2.0789E+00  1.3419E-01  5.5623E+00  6.0999E-02
             7.2390E+00
 PARAMETER:  3.9002E-01  5.7618E-01  7.4974E+00  1.1313E-01  9.9174E-01  1.2486E+00  8.3186E-01 -1.9085E+00  1.8160E+00 -2.6969E+00
             2.0795E+00
 GRADIENT:   3.9552E+00 -1.1846E+00  1.1975E+00  3.9413E+00  5.4791E+00  2.0816E+00  5.2073E+00  7.8499E-05  5.5200E+00 -3.5754E-02
            -2.0974E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1497.13074279423        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      614
 NPARAMETR:  1.2498E+00  1.8570E+00  3.9231E+02  6.8511E-01  2.3696E+00  2.9354E+00  1.7618E+00  1.9833E-02  6.3436E+00  1.4898E-01
             7.3326E+00
 PARAMETER:  3.2297E-01  7.1898E-01  6.0720E+00 -2.7817E-01  9.6273E-01  1.1768E+00  6.6634E-01 -3.8204E+00  1.9474E+00 -1.8039E+00
             2.0923E+00
 GRADIENT:  -9.4678E+00 -5.8427E-01  7.6315E-01 -1.2720E+00 -6.5625E+00 -7.4261E+00 -6.5182E+00  1.6691E-04 -1.5489E+01 -2.6024E-01
            -2.6735E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1499.75301970652        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      791
 NPARAMETR:  1.2712E+00  2.2643E+00  2.8172E+02  4.6477E-01  2.3758E+00  2.9436E+00  2.0810E+00  1.0000E-02  7.7935E+00  2.3861E-01
             7.3279E+00
 PARAMETER:  3.3993E-01  9.1726E-01  5.7409E+00 -6.6622E-01  9.6531E-01  1.1796E+00  8.3287E-01 -5.2647E+00  2.1533E+00 -1.3329E+00
             2.0917E+00
 GRADIENT:  -5.2226E+00 -4.7358E-01 -8.0523E-01  2.3332E+00 -1.3935E-01 -6.3650E+00  3.2770E-01  0.0000E+00  2.2578E+00 -5.3690E-01
             1.2091E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1500.19836915803        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      967
 NPARAMETR:  1.2861E+00  2.4132E+00  2.9490E+02  3.5791E-01  2.3808E+00  2.9697E+00  2.1598E+00  1.0000E-02  8.1924E+00  2.6304E-01
             7.3177E+00
 PARAMETER:  3.5161E-01  9.8096E-01  5.7866E+00 -9.2748E-01  9.6745E-01  1.1885E+00  8.7003E-01 -5.7119E+00  2.2032E+00 -1.2355E+00
             2.0903E+00
 GRADIENT:  -2.3064E+00  2.3407E+00 -5.0238E-01  4.2364E-01  8.0269E-01 -3.1036E+00  1.3813E+00  0.0000E+00 -1.7984E+00 -6.1691E-01
            -2.5499E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1500.96179615320        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:     1092
 NPARAMETR:  1.2858E+00  2.4112E+00  2.9317E+02  3.5756E-01  2.3334E+00  2.9665E+00  2.1577E+00  1.0000E-02  8.2139E+00  6.5973E-01
             7.3326E+00
 PARAMETER:  3.5137E-01  9.8012E-01  5.7807E+00 -9.2844E-01  9.4733E-01  1.1874E+00  8.6905E-01 -5.7172E+00  2.2058E+00 -3.1592E-01
             2.0923E+00
 GRADIENT:   1.1843E+00  8.4516E+00 -3.8156E-01  1.7478E+00  1.4542E+00  5.5593E+00  4.1659E+00  0.0000E+00  3.0377E+01 -1.1079E-01
             2.2757E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1500.96581621334        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1162
 NPARAMETR:  1.2858E+00  2.4111E+00  2.9324E+02  3.5755E-01  2.3186E+00  2.9665E+00  2.1577E+00  1.0000E-02  8.2047E+00  6.9675E-01
             7.3294E+00
 PARAMETER:  3.5138E-01  9.8010E-01  5.7810E+00 -9.2848E-01  9.4095E-01  1.1874E+00  8.6906E-01 -5.7172E+00  2.2047E+00 -2.6133E-01
             2.0919E+00
 GRADIENT:   1.2156E+00  8.7631E+00 -3.1818E-01  1.7382E+00 -1.1346E+00  5.5903E+00  4.1745E+00  0.0000E+00  3.0132E+01  1.7297E-01
             2.3702E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1500.96647076247        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1233
 NPARAMETR:  1.2859E+00  2.4110E+00  2.9354E+02  3.5749E-01  2.3080E+00  2.9668E+00  2.1579E+00  1.0000E-02  8.1731E+00  7.2416E-01
             7.3181E+00
 PARAMETER:  3.5143E-01  9.8005E-01  5.7820E+00 -9.2864E-01  9.3640E-01  1.1875E+00  8.6911E-01 -5.7172E+00  2.2008E+00 -2.2275E-01
             2.0904E+00
 GRADIENT:   1.3116E+00  9.3472E+00 -2.7130E-01  1.6290E+00 -2.9402E+00  5.5962E+00  4.0540E+00  0.0000E+00  2.9323E+01  3.9401E-01
             2.1971E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1500.96696881176        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:     1353
 NPARAMETR:  1.2858E+00  2.4104E+00  2.9320E+02  3.5756E-01  2.3071E+00  2.9660E+00  2.1574E+00  1.0000E-02  8.1692E+00  7.2794E-01
             7.3189E+00
 PARAMETER:  3.5144E-01  9.8003E-01  5.7823E+00 -9.2868E-01  9.3575E-01  1.1875E+00  8.6913E-01 -5.7172E+00  2.1998E+00 -2.1731E-01
             2.0899E+00
 GRADIENT:   4.9200E+03  8.8556E+02  2.9836E+02 -1.8642E+03 -9.2810E+02  1.4405E+03  9.9846E+02  0.0000E+00 -3.9571E+02  4.1972E-01
            -8.0120E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1353
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.8424E-03 -8.9770E-02  5.1982E-06  5.9616E-02 -6.8626E-03
 SE:             2.9556E-02  1.7360E-02  1.1623E-05  2.1719E-02  1.1737E-02
 N:                     100         100         100         100         100

 P VAL.:         7.6481E-01  2.3305E-07  6.5471E-01  6.0532E-03  5.5874E-01

 ETASHRINKSD(%)  9.8359E-01  4.1843E+01  9.9961E+01  2.7239E+01  6.0681E+01
 ETASHRINKVR(%)  1.9575E+00  6.6178E+01  1.0000E+02  4.7059E+01  8.4540E+01
 EBVSHRINKSD(%)  2.0325E+00  3.1954E+01  9.9933E+01  2.8938E+01  6.0679E+01
 EBVSHRINKVR(%)  4.0236E+00  5.3698E+01  1.0000E+02  4.9502E+01  8.4539E+01
 RELATIVEINF(%)  9.5893E+01  2.4081E+01  4.2895E-05  2.6432E+01  1.4781E+01
 EPSSHRINKSD(%)  9.6858E+00
 EPSSHRINKVR(%)  1.8433E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1500.9669688117640     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       153.12239095664677     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    36.56
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    18.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1500.967       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.29E+00  2.41E+00  2.94E+02  3.57E-01  2.31E+00  2.97E+00  2.16E+00  1.00E-02  8.16E+00  7.28E-01  7.32E+00
 


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
+        2.12E+06
 
 TH 2
+       -1.25E+01  7.76E+04
 
 TH 3
+       -3.06E+00 -2.42E-03  1.49E-01
 
 TH 4
+        4.91E+03  3.18E+01  4.16E+00  3.93E+06
 
 TH 5
+        6.63E+02  1.22E+02  1.69E-01 -9.06E+02  9.28E+04
 
 TH 6
+       -4.61E+02 -1.27E+00 -3.92E-01  2.05E+01  8.50E+01  3.42E+04
 
 TH 7
+       -8.69E+02 -1.07E+01 -7.36E-01  3.50E+01  1.61E+02 -4.36E-01  1.23E+05
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.08E+01 -1.93E+00  7.70E-02  2.43E+00 -1.65E+01  1.36E+02  5.28E-01  0.00E+00  1.34E+03
 
 TH10
+        1.53E+02  2.51E+01  3.67E-02 -1.84E+02 -2.26E+01  2.01E+01  3.58E+01  0.00E+00 -3.40E+00  1.05E+00
 
 TH11
+        1.03E+02 -2.29E+00  9.07E-02 -7.15E+00 -2.08E+01  1.61E+02  2.15E+00  0.00E+00 -3.67E+01 -1.15E+00  1.83E+03
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       55.015
Stop Time:
Sat Sep 18 07:13:37 CDT 2021
