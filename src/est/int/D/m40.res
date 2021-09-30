Wed Sep 29 08:44:06 CDT 2021
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
$DATA ../../../../data/int/D/dat40.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25071.7836339058        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2167E+02  4.1036E+02 -2.7624E+01  3.8815E+02  2.1909E+02 -1.9320E+03 -9.5626E+02 -5.1115E+01 -1.3701E+03 -5.4204E+02
            -5.2483E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -952.508659081232        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1585E+00  1.9111E+00  9.2890E-01  2.1095E+00  8.5340E-01  3.7354E+00  3.7740E+00  9.8413E-01  2.3110E+00  1.5877E+00
             1.2813E+01
 PARAMETER:  2.4710E-01  7.4767E-01  2.6250E-02  8.4644E-01 -5.8527E-02  1.4179E+00  1.4281E+00  8.4001E-02  9.3769E-01  5.6231E-01
             2.6505E+00
 GRADIENT:  -2.2364E+01  4.7200E+01 -2.8923E+01  1.5624E+02 -3.3501E+01  1.4354E+02 -2.6979E+01  4.6454E+00  9.8461E+00  3.6095E+01
             5.2958E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1054.52537277760        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0844E+00  1.5011E+00  4.0802E+01  3.3591E+00  2.1277E+00  2.8567E+00  9.4819E+00  6.8162E-01  2.7941E+00  2.6419E+00
             1.2190E+01
 PARAMETER:  1.8100E-01  5.0621E-01  3.8087E+00  1.3117E+00  8.5505E-01  1.1497E+00  2.3494E+00 -2.8329E-01  1.1275E+00  1.0715E+00
             2.6006E+00
 GRADIENT:  -4.0343E+01  3.1876E+01 -9.1354E-01  1.5505E+02 -1.9722E+01  1.1048E+02  7.8999E+01  3.6548E-02  2.6046E+01  8.3932E+01
             5.2397E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1268.82621742940        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.1521E+00  1.2612E+00  4.3188E+00  9.4579E-01  1.9385E+00  1.8242E+00  3.9939E+00  9.0941E-01  1.3954E+00  6.3406E-01
             9.5406E+00
 PARAMETER:  2.4159E-01  3.3208E-01  1.5630E+00  4.4269E-02  7.6192E-01  7.0116E-01  1.4848E+00  5.0407E-03  4.3316E-01 -3.5561E-01
             2.3556E+00
 GRADIENT:  -3.8326E+00 -3.3016E+01 -1.1056E+01 -2.1069E+01  3.7887E+01 -1.6986E+01  6.8134E+00  5.9595E-01  1.7803E+01  7.8122E+00
             3.8501E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1303.34695946703        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0504E+00  1.6575E+00  3.2158E+00  6.8737E-01  1.8910E+00  1.8887E+00  3.4410E+00  7.5386E-01  1.3051E+00  4.3758E-01
             7.5898E+00
 PARAMETER:  1.4913E-01  6.0529E-01  1.2681E+00 -2.7489E-01  7.3711E-01  7.3590E-01  1.3358E+00 -1.8255E-01  3.6625E-01 -7.2651E-01
             2.1268E+00
 GRADIENT:  -2.0963E+01 -6.9086E+00 -1.1555E+00 -9.9601E+00  9.3754E+00 -1.1656E+01 -6.7833E-03  1.8100E-01  8.7019E+00  2.9505E+00
             9.4972E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1304.92739297271        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:      453
 NPARAMETR:  1.0524E+00  1.6683E+00  3.2345E+00  6.8995E-01  1.8981E+00  1.8851E+00  3.4908E+00  1.7521E-01  1.2869E+00  4.3679E-01
             7.7045E+00
 PARAMETER:  1.5109E-01  6.1178E-01  1.2739E+00 -2.7114E-01  7.4085E-01  7.3396E-01  1.3501E+00 -1.6418E+00  3.5225E-01 -7.2830E-01
             2.1418E+00
 GRADIENT:  -2.2183E+01 -5.8655E+00 -1.2460E+00 -1.1628E+01  1.0487E+01 -1.2475E+01  5.6136E+00  1.0154E-02  8.9630E+00  2.9848E+00
             4.0222E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1311.39816542215        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      634
 NPARAMETR:  1.0524E+00  1.6975E+00  3.2674E+00  6.9221E-01  1.8994E+00  1.9129E+00  4.1782E+00  1.0000E-02  1.2825E+00  4.3772E-01
             7.6747E+00
 PARAMETER:  1.5108E-01  6.2918E-01  1.2840E+00 -2.6786E-01  7.4153E-01  7.4860E-01  1.5299E+00 -6.4898E+00  3.4885E-01 -7.2619E-01
             2.1379E+00
 GRADIENT:  -3.4141E+01 -1.1287E+01 -3.1080E+00 -2.9755E+01  8.9441E+00 -3.5129E+01  5.1885E+00  0.0000E+00  9.9299E+00  2.7819E+00
            -2.6753E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1320.11082578182        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      791
 NPARAMETR:  1.1506E+00  1.6316E+00  3.7143E+00  8.2892E-01  1.8620E+00  2.2923E+00  4.5192E+00  1.0000E-02  5.6150E-01  2.0747E-01
             7.7397E+00
 PARAMETER:  2.4026E-01  5.8954E-01  1.4122E+00 -8.7635E-02  7.2166E-01  9.2958E-01  1.6083E+00 -6.4898E+00 -4.7714E-01 -1.4728E+00
             2.1464E+00
 GRADIENT:   3.9897E-01  4.4628E+00 -2.1991E+00 -2.8920E+00  6.8920E+00  3.2380E+00 -4.4103E+00  0.0000E+00 -2.5589E-01  6.1259E-01
            -5.1874E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1323.26075814306        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      967
 NPARAMETR:  1.1486E+00  1.1069E+00  5.9454E+00  1.0611E+00  1.8584E+00  2.2791E+00  5.3407E+00  1.0000E-02  9.2524E-01  6.3445E-02
             7.7502E+00
 PARAMETER:  2.3857E-01  2.0156E-01  1.8826E+00  1.5929E-01  7.1969E-01  9.2378E-01  1.7754E+00 -6.4898E+00  2.2297E-02 -2.6576E+00
             2.1477E+00
 GRADIENT:  -3.9329E-01 -1.1863E+00 -2.7690E-02  1.6419E-01 -2.3777E+00  1.4180E+00 -3.9440E-01  0.0000E+00 -1.1519E-01  4.4443E-02
            -1.9808E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1323.38521679139        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1135             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1508E+00  1.1091E+00  6.0410E+00  1.0635E+00  1.8664E+00  2.2789E+00  5.4600E+00  1.0000E-02  9.0541E-01  1.0000E-02
             7.7571E+00
 PARAMETER:  2.4043E-01  2.0352E-01  1.8986E+00  1.6161E-01  7.2402E-01  9.2367E-01  1.7974E+00 -6.4898E+00  6.3162E-04 -4.6208E+00
             2.1486E+00
 GRADIENT:   2.3067E+01  3.4866E+00  2.4354E-01  4.5015E+00  3.2898E+00  3.5300E+01  1.1029E+02  0.0000E+00 -2.9094E-01  0.0000E+00
             3.7682E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1323.39778545417        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1313
 NPARAMETR:  1.1537E+00  1.1076E+00  6.3631E+00  1.0753E+00  1.8812E+00  2.2775E+00  5.4449E+00  1.0000E-02  9.4544E-01  1.0000E-02
             7.7724E+00
 PARAMETER:  2.4301E-01  2.0224E-01  1.9505E+00  1.7256E-01  7.3189E-01  9.2310E-01  1.7947E+00 -6.4898E+00  4.3900E-02 -4.6382E+00
             2.1506E+00
 GRADIENT:   1.0223E+00  2.0577E-01 -1.0246E-01 -7.7359E-03 -7.9452E-01  1.1936E+00  3.3586E+00  0.0000E+00  2.4615E-01  0.0000E+00
             2.7398E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1325.93258220303        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     1455
 NPARAMETR:  1.1437E+00  1.1347E+00  6.5347E+00  1.0706E+00  1.8969E+00  2.2299E+00  5.2780E+00  1.0000E-02  9.5631E-01  1.0000E-02
             7.9881E+00
 PARAMETER:  2.3425E-01  2.2639E-01  1.9771E+00  1.6824E-01  7.4025E-01  9.0195E-01  1.7636E+00 -6.4898E+00  5.5324E-02 -4.6634E+00
             2.1779E+00
 GRADIENT:   2.6053E+01  2.7766E+00  4.2394E-01  2.9318E+00 -2.3524E+00  5.4715E+01  1.0219E+02  0.0000E+00  1.7070E+00  0.0000E+00
             7.8640E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1327.17851805022        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1592
 NPARAMETR:  1.1276E+00  1.1540E+00  6.6573E+00  1.0517E+00  1.9213E+00  2.0458E+00  5.2693E+00  1.0000E-02  9.1438E-01  1.0000E-02
             7.8936E+00
 PARAMETER:  2.2011E-01  2.4319E-01  1.9957E+00  1.5038E-01  7.5298E-01  8.1581E-01  1.7619E+00 -6.4898E+00  1.0495E-02 -4.6634E+00
             2.1660E+00
 GRADIENT:  -1.1711E-01 -1.0163E+00 -2.5519E-01 -3.3396E+00 -2.3683E+00 -4.1746E+00 -2.1179E-01  0.0000E+00  8.7878E-01  0.0000E+00
             1.3847E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1327.38798155337        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1769
 NPARAMETR:  1.1260E+00  1.1566E+00  7.7405E+00  1.0639E+00  1.9760E+00  2.0748E+00  5.3126E+00  1.0000E-02  9.1595E-01  1.0000E-02
             7.8435E+00
 PARAMETER:  2.1870E-01  2.4549E-01  2.1465E+00  1.6193E-01  7.8107E-01  8.2985E-01  1.7701E+00 -6.4898E+00  1.2208E-02 -4.6634E+00
             2.1597E+00
 GRADIENT:  -2.1856E-02  2.5173E-01 -5.7562E-02 -1.3572E-01  5.9694E-01  2.7770E-01  1.0347E-01  0.0000E+00 -1.0304E-01  0.0000E+00
            -7.3571E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1327.63468926014        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1939
 NPARAMETR:  1.1283E+00  1.0807E+00  8.4951E+00  1.0954E+00  1.9836E+00  2.0769E+00  5.5899E+00  1.0000E-02  9.5853E-01  1.0000E-02
             7.8436E+00
 PARAMETER:  2.2072E-01  1.7757E-01  2.2395E+00  1.9115E-01  7.8491E-01  8.3086E-01  1.8210E+00 -6.4898E+00  5.7648E-02 -4.6634E+00
             2.1597E+00
 GRADIENT:   1.1241E+00 -1.3026E-01  1.1562E-01 -5.3344E+00 -7.2195E-02  4.9770E-01  7.9569E+00  0.0000E+00  7.4865E-01  0.0000E+00
             3.9600E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1327.67236923394        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2121
 NPARAMETR:  1.1300E+00  1.0669E+00  8.5382E+00  1.1134E+00  1.9820E+00  2.0816E+00  5.5880E+00  1.0000E-02  9.6553E-01  1.0000E-02
             7.8569E+00
 PARAMETER:  2.2223E-01  1.6472E-01  2.2446E+00  2.0745E-01  7.8410E-01  8.3316E-01  1.8206E+00 -6.4898E+00  6.4923E-02 -4.6634E+00
             2.1614E+00
 GRADIENT:   1.4739E+00  2.7850E-01 -5.9005E-02 -1.0016E-01 -2.1010E-01  1.5091E+00  5.2876E+00  0.0000E+00 -1.3228E-01  0.0000E+00
             1.4574E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1327.70162309772        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2303
 NPARAMETR:  1.1264E+00  1.0560E+00  8.6647E+00  1.1178E+00  1.9842E+00  2.0777E+00  5.6951E+00  1.0000E-02  9.7661E-01  1.0000E-02
             7.8503E+00
 PARAMETER:  2.1899E-01  1.5453E-01  2.2593E+00  2.1134E-01  7.8521E-01  8.3126E-01  1.8396E+00 -6.4898E+00  7.6327E-02 -4.6634E+00
             2.1606E+00
 GRADIENT:   1.2953E-01  5.7782E-01 -1.6171E-01 -2.8574E+00  5.3840E-01  7.1481E-01  9.4787E+00  0.0000E+00  4.9338E-01  0.0000E+00
             9.6305E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1327.71480138552        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2478
 NPARAMETR:  1.1298E+00  1.0386E+00  8.8112E+00  1.1279E+00  1.9857E+00  2.0837E+00  5.6694E+00  1.0000E-02  9.7700E-01  1.0000E-02
             7.8607E+00
 PARAMETER:  2.2205E-01  1.3791E-01  2.2760E+00  2.2031E-01  7.8595E-01  8.3414E-01  1.8351E+00 -6.4898E+00  7.6733E-02 -4.6634E+00
             2.1619E+00
 GRADIENT:   1.3458E+00  1.3096E-01 -9.2125E-02 -1.8996E-02 -1.2517E-01  1.8815E+00  6.5555E+00  0.0000E+00 -2.4446E-01  0.0000E+00
             2.1388E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1327.73031084430        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2659
 NPARAMETR:  1.1285E+00  1.0326E+00  8.9046E+00  1.1315E+00  1.9867E+00  2.0811E+00  5.7042E+00  1.0000E-02  9.8490E-01  1.0000E-02
             7.8565E+00
 PARAMETER:  2.2089E-01  1.3205E-01  2.2866E+00  2.2354E-01  7.8647E-01  8.3290E-01  1.8412E+00 -6.4898E+00  8.4780E-02 -4.6634E+00
             2.1613E+00
 GRADIENT:   8.8283E-01  2.3962E-01 -9.0870E-02 -6.1171E-01 -1.3787E-01  1.4050E+00  7.6251E+00  0.0000E+00 -3.5475E-02  0.0000E+00
             1.3284E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1327.73985599185        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2840
 NPARAMETR:  1.1259E+00  1.0233E+00  9.0215E+00  1.1346E+00  1.9890E+00  2.0777E+00  5.7324E+00  1.0000E-02  9.9123E-01  1.0000E-02
             7.8516E+00
 PARAMETER:  2.1860E-01  1.2303E-01  2.2996E+00  2.2632E-01  7.8764E-01  8.3126E-01  1.8461E+00 -6.4898E+00  9.1192E-02 -4.6634E+00
             2.1607E+00
 GRADIENT:  -9.5295E-02  1.0200E-01 -1.0052E-01 -1.3352E+00  8.0923E-02  7.9650E-01  8.3381E+00  0.0000E+00  1.2943E-01  0.0000E+00
             4.9027E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1327.74847900094        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3029             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1266E+00  1.0089E+00  9.1672E+00  1.1417E+00  1.9904E+00  2.0779E+00  5.7696E+00  1.0000E-02  1.0040E+00  1.0000E-02
             7.8521E+00
 PARAMETER:  2.1922E-01  1.0888E-01  2.3156E+00  2.3255E-01  7.8833E-01  8.3137E-01  1.8526E+00 -6.4898E+00  1.0394E-01 -4.6634E+00
             2.1608E+00
 GRADIENT:   2.0380E+01  1.6660E+00  1.8585E-01  6.4377E+00  4.4294E+00  2.8869E+01  1.2827E+02  0.0000E+00  6.3304E-01  0.0000E+00
             3.9607E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1327.75275193117        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     3207
 NPARAMETR:  1.1271E+00  1.0064E+00  9.2292E+00  1.1465E+00  1.9913E+00  2.0792E+00  5.7677E+00  1.0000E-02  9.9935E-01  1.0000E-02
             7.8528E+00
 PARAMETER:  2.1968E-01  1.0641E-01  2.3224E+00  2.3668E-01  7.8878E-01  8.3197E-01  1.8523E+00 -6.4898E+00  9.9352E-02 -4.6634E+00
             2.1609E+00
 GRADIENT:   3.6380E-01  1.3986E-01 -8.3067E-02  1.6395E-01 -2.9834E-01  1.0950E+00  8.0571E+00  0.0000E+00 -2.0792E-01  0.0000E+00
            -1.3599E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1327.75566933980        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3392
 NPARAMETR:  1.1266E+00  1.0020E+00  9.3247E+00  1.1499E+00  1.9927E+00  2.0780E+00  5.7848E+00  1.0000E-02  1.0010E+00  1.0000E-02
             7.8508E+00
 PARAMETER:  2.1923E-01  1.0198E-01  2.3327E+00  2.3968E-01  7.8947E-01  8.3143E-01  1.8552E+00 -6.4898E+00  1.0102E-01 -4.6634E+00
             2.1606E+00
 GRADIENT:   1.8023E-01  2.3099E-01 -4.8855E-02  5.3654E-01 -5.0177E-01  9.0686E-01  8.2528E+00  0.0000E+00 -3.3364E-01  0.0000E+00
            -9.2990E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1327.75889053227        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3579             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1267E+00  9.9724E-01  9.4199E+00  1.1524E+00  1.9946E+00  2.0780E+00  5.7973E+00  1.0000E-02  1.0059E+00  1.0000E-02
             7.8513E+00
 PARAMETER:  2.1925E-01  9.7234E-02  2.3428E+00  2.4182E-01  7.9042E-01  8.3140E-01  1.8574E+00 -6.4898E+00  1.0586E-01 -4.6634E+00
             2.1607E+00
 GRADIENT:   2.0363E+01  1.7287E+00  2.5269E-01  8.9904E+00  3.9084E+00  2.8919E+01  1.2862E+02  0.0000E+00  4.3389E-02  0.0000E+00
             3.8192E+01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1327.76020199571        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     3757
 NPARAMETR:  1.1271E+00  9.9373E-01  9.4654E+00  1.1529E+00  1.9960E+00  2.0789E+00  5.8047E+00  1.0000E-02  1.0126E+00  1.0000E-02
             7.8541E+00
 PARAMETER:  2.1964E-01  9.3707E-02  2.3476E+00  2.4227E-01  7.9114E-01  8.3182E-01  1.8587E+00 -6.4898E+00  1.1250E-01 -4.6634E+00
             2.1610E+00
 GRADIENT:   3.4820E-01  1.6245E-03 -7.7809E-02 -6.9794E-01 -1.1246E-01  1.0263E+00  8.8090E+00  0.0000E+00  1.0195E-01  0.0000E+00
             4.8946E-01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1327.76113172628        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3942
 NPARAMETR:  1.1267E+00  9.8889E-01  9.5358E+00  1.1545E+00  1.9976E+00  2.0779E+00  5.8147E+00  1.0000E-02  1.0168E+00  1.0000E-02
             7.8532E+00
 PARAMETER:  2.1928E-01  8.8823E-02  2.3551E+00  2.4364E-01  7.9196E-01  8.3133E-01  1.8604E+00 -6.4898E+00  1.1667E-01 -4.6634E+00
             2.1609E+00
 GRADIENT:   2.0056E-01 -1.2443E-01 -9.0684E-02 -1.1273E+00  7.4436E-02  8.4505E-01  9.0055E+00  0.0000E+00  2.3079E-01  0.0000E+00
             4.6420E-01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1327.76289042152        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     4127
 NPARAMETR:  1.1266E+00  9.8680E-01  9.6015E+00  1.1565E+00  1.9991E+00  2.0780E+00  5.8207E+00  1.0000E-02  1.0174E+00  1.0000E-02
             7.8528E+00
 PARAMETER:  2.1921E-01  8.6709E-02  2.3619E+00  2.4540E-01  7.9268E-01  8.3139E-01  1.8614E+00 -6.4898E+00  1.1729E-01 -4.6634E+00
             2.1609E+00
 GRADIENT:   1.6723E-01 -7.6457E-02 -8.0552E-02 -7.7257E-01  4.4242E-02  8.7153E-01  8.9598E+00  0.0000E+00  1.2897E-01  0.0000E+00
             1.5836E-01

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1327.76390573725        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     4314
 NPARAMETR:  1.1266E+00  9.8603E-01  9.6503E+00  1.1580E+00  1.9998E+00  2.0780E+00  5.8241E+00  1.0000E-02  1.0176E+00  1.0000E-02
             7.8526E+00
 PARAMETER:  2.1924E-01  8.5933E-02  2.3670E+00  2.4668E-01  7.9307E-01  8.3138E-01  1.8620E+00 -6.4898E+00  1.1742E-01 -4.6634E+00
             2.1608E+00
 GRADIENT:   1.7518E-01  9.6475E-05 -5.3889E-02 -4.1330E-01 -1.1690E-01  8.7714E-01  8.8962E+00  0.0000E+00  2.7760E-02  0.0000E+00
            -1.1169E-01

0ITERATION NO.:  140    OBJECTIVE VALUE:  -1327.76483779467        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     4508             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1266E+00  9.8386E-01  9.7386E+00  1.1607E+00  2.0012E+00  2.0780E+00  5.8314E+00  1.0000E-02  1.0182E+00  1.0000E-02
             7.8521E+00
 PARAMETER:  2.1921E-01  8.3727E-02  2.3761E+00  2.4899E-01  7.9375E-01  8.3142E-01  1.8632E+00 -6.4898E+00  1.1804E-01 -4.6634E+00
             2.1608E+00
 GRADIENT:   2.0331E+01  1.6357E+00  2.7396E-01  9.2161E+00  3.9983E+00  2.8917E+01  1.3015E+02  0.0000E+00  1.8913E-01  0.0000E+00
             3.8333E+01

0ITERATION NO.:  145    OBJECTIVE VALUE:  -1327.76516358638        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     4693
 NPARAMETR:  1.1266E+00  9.8260E-01  9.7906E+00  1.1621E+00  2.0021E+00  2.0780E+00  5.8356E+00  1.0000E-02  1.0195E+00  1.0000E-02
             7.8521E+00
 PARAMETER:  2.1922E-01  8.2448E-02  2.3814E+00  2.5023E-01  7.9422E-01  8.3142E-01  1.8640E+00 -6.4898E+00  1.1930E-01 -4.6634E+00
             2.1608E+00
 GRADIENT:   1.6475E-01  1.4220E-01  4.8670E-03  3.0906E-01 -4.5539E-01  9.0458E-01  8.8046E+00  0.0000E+00 -1.6835E-01  0.0000E+00
            -6.8489E-01

0ITERATION NO.:  150    OBJECTIVE VALUE:  -1327.76576875007        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     4880
 NPARAMETR:  1.1266E+00  9.8079E-01  9.8253E+00  1.1627E+00  2.0031E+00  2.0778E+00  5.8395E+00  1.0000E-02  1.0215E+00  1.0000E-02
             7.8524E+00
 PARAMETER:  2.1921E-01  8.0601E-02  2.3850E+00  2.5077E-01  7.9472E-01  8.3132E-01  1.8646E+00 -6.4898E+00  1.2128E-01 -4.6634E+00
             2.1608E+00
 GRADIENT:   1.5729E-01  8.7913E-02 -8.1964E-03  9.0518E-02 -3.1048E-01  8.6781E-01  8.8960E+00  0.0000E+00 -9.5140E-02  0.0000E+00
            -5.0959E-01

0ITERATION NO.:  153    OBJECTIVE VALUE:  -1327.76592203965        NO. OF FUNC. EVALS.: 107
 CUMULATIVE NO. OF FUNC. EVALS.:     4987
 NPARAMETR:  1.1266E+00  9.7961E-01  9.8850E+00  1.1644E+00  2.0041E+00  2.0782E+00  5.8437E+00  1.0000E-02  1.0220E+00  1.0000E-02
             7.8524E+00
 PARAMETER:  2.1921E-01  7.9247E-02  2.3852E+00  2.5064E-01  7.9493E-01  8.3135E-01  1.8649E+00 -6.4898E+00  1.2298E-01 -4.6634E+00
             2.1609E+00
 GRADIENT:  -5.4003E-03 -2.8059E-03 -3.4866E-02 -2.4683E-01 -6.3406E-02 -1.6049E-02 -3.6308E-02  0.0000E+00  1.9923E-02  0.0000E+00
             5.1477E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4987
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.7602E-03  3.0728E-02 -1.8116E-05 -6.9119E-02  5.9430E-05
 SE:             2.8656E-02  2.4948E-02  1.8863E-05  1.3291E-02  1.6997E-04
 N:                     100         100         100         100         100

 P VAL.:         7.5983E-01  2.1807E-01  3.3685E-01  1.9918E-07  7.2660E-01

 ETASHRINKSD(%)  3.9999E+00  1.6420E+01  9.9937E+01  5.5474E+01  9.9431E+01
 ETASHRINKVR(%)  7.8398E+00  3.0144E+01  1.0000E+02  8.0174E+01  9.9997E+01
 EBVSHRINKSD(%)  4.0785E+00  1.1013E+01  9.9929E+01  5.9663E+01  9.9426E+01
 EBVSHRINKVR(%)  7.9907E+00  2.0813E+01  1.0000E+02  8.3729E+01  9.9997E+01
 RELATIVEINF(%)  9.1832E+01  4.0323E+01  9.8240E-06  8.2756E+00  6.5005E-04
 EPSSHRINKSD(%)  7.2631E+00
 EPSSHRINKVR(%)  1.3999E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1327.7659220396470     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       326.32343772876379     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   172.01
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    16.61
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1327.766       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.13E+00  9.79E-01  9.83E+00  1.16E+00  2.00E+00  2.08E+00  5.84E+00  1.00E-02  1.02E+00  1.00E-02  7.85E+00
 


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
+        1.67E+02
 
 TH 2
+       -5.26E+00  6.20E+00
 
 TH 3
+       -4.63E-02  5.97E-03  6.69E-02
 
 TH 4
+       -5.44E+01  3.15E+01 -6.88E-01  1.73E+02
 
 TH 5
+        8.52E+00 -4.08E+00 -2.28E+00  3.12E+00  8.04E+01
 
 TH 6
+        1.55E+01 -1.65E+00  6.88E-02 -1.16E+01 -9.96E-01  1.77E+00
 
 TH 7
+        4.90E+00 -2.01E+00 -2.93E-02 -1.05E+01  2.33E+00  7.87E-01  7.28E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.57E+01 -7.97E+00 -3.45E-02 -4.18E+01  6.36E+00  2.89E+00  2.78E+00  0.00E+00  1.08E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.83E-01 -1.65E+00  4.36E-02 -8.66E+00 -5.13E-01  4.07E-01  4.87E-01  0.00E+00  2.05E+00  0.00E+00  9.31E-01
 
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
+        1.83E+02
 
 TH 2
+       -9.08E-01  2.79E+01
 
 TH 3
+        1.19E-01  1.76E-01  8.83E-02
 
 TH 4
+       -5.21E+00  3.20E+01 -6.45E-01  1.68E+02
 
 TH 5
+       -3.23E+00 -8.13E+00 -2.35E+00  3.85E-01  8.35E+01
 
 TH 6
+        1.19E+00 -1.42E-01  1.10E-02  1.97E+00 -4.76E-01  3.87E+01
 
 TH 7
+        2.52E-01  2.94E+00 -4.50E-02 -1.02E+01  1.59E+00 -2.14E-01  3.43E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.60E-01 -2.79E+00 -2.24E-01 -3.78E+01  6.17E+00 -9.82E-01  2.30E+00  0.00E+00  2.26E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.57E+00 -2.37E+00 -1.60E-02 -1.14E+01  6.14E-01  1.62E+00  4.26E-01  0.00E+00  3.63E+00  0.00E+00  1.73E+01
 
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
+        1.87E+02
 
 TH 2
+        4.86E+01  2.75E+01
 
 TH 3
+        4.41E-01  2.63E-01  5.64E-02
 
 TH 4
+        8.61E+01  3.44E+01  1.45E-01  1.80E+02
 
 TH 5
+       -2.72E+01 -1.69E+00 -1.87E+00 -1.38E+01  8.40E+01
 
 TH 6
+        2.71E+01  6.57E+00  1.20E-01 -2.20E+01 -5.84E+00  3.79E+01
 
 TH 7
+       -2.32E+00  4.83E+00 -7.81E-02 -1.23E+01  1.15E+01  3.35E+00  6.41E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.13E+01 -3.38E+00 -5.54E-02 -4.11E+01  4.59E+00  7.36E+00  2.51E+00  0.00E+00  2.01E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.09E+02 -8.93E+00 -2.69E-01 -1.10E+02  2.49E+01  2.93E+01  1.73E+01  0.00E+00  2.06E+01  0.00E+00  6.61E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,      187.658
Stop Time:
Wed Sep 29 08:47:17 CDT 2021
