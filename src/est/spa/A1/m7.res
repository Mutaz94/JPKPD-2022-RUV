Sat Sep 25 07:51:07 CDT 2021
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
$DATA ../../../../data/spa/A1/dat7.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1177.69255780186        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.0429E+01  1.0577E+01  1.5818E+01  2.7651E+01  6.5779E+01 -4.5343E+01  6.8116E+00 -3.8072E+01  1.1956E+01 -1.7670E+01
            -8.6659E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1434.76597291697        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0764E+00  1.0450E+00  9.5130E-01  9.9917E-01  9.5797E-01  1.2419E+00  8.9247E-01  1.1703E+00  8.1334E-01  8.5553E-01
             1.8853E+00
 PARAMETER:  1.7361E-01  1.4400E-01  5.0077E-02  9.9172E-02  5.7063E-02  3.1666E-01 -1.3764E-02  2.5723E-01 -1.0660E-01 -5.6030E-02
             7.3406E-01
 GRADIENT:   1.3614E+02  2.7921E+01  2.6632E+00  4.2557E+01  9.0386E+00  3.2780E+01 -2.6395E+00 -1.1954E+01 -7.4939E+00  2.6361E+00
            -9.9673E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1445.17027566345        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.0287E+00  9.0368E-01  8.0373E-01  1.0553E+00  7.8600E-01  1.0850E+00  1.4033E+00  1.1689E+00  6.7271E-01  3.8045E-01
             2.0118E+00
 PARAMETER:  1.2827E-01 -1.2761E-03 -1.1850E-01  1.5378E-01 -1.4080E-01  1.8160E-01  4.3885E-01  2.5607E-01 -2.9643E-01 -8.6639E-01
             7.9902E-01
 GRADIENT:   7.9142E+01  3.5812E+01  2.4252E+01  9.4897E+00 -3.5155E+01 -8.7435E+00  1.9280E+01 -6.5438E+00  3.7342E+00 -8.5723E-01
            -6.9269E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1455.16829013038        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      237
 NPARAMETR:  9.8295E-01  8.8430E-01  1.7036E+00  1.1018E+00  1.1272E+00  1.1037E+00  1.0734E+00  2.1922E+00  6.1753E-01  3.8814E-01
             2.4437E+00
 PARAMETER:  8.2800E-02 -2.2963E-02  6.3271E-01  1.9699E-01  2.1971E-01  1.9865E-01  1.7086E-01  8.8491E-01 -3.8203E-01 -8.4638E-01
             9.9350E-01
 GRADIENT:  -9.0876E+00  7.0668E-01 -6.4182E+00  1.5008E+01  3.7663E+00  3.4057E+00 -1.9071E+00  2.3194E+00  1.6230E+00  6.1788E-01
             1.1264E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1456.78105369297        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      311
 NPARAMETR:  9.9753E-01  8.2396E-01  2.4667E+00  1.1328E+00  1.2273E+00  1.0814E+00  1.3280E+00  2.7041E+00  4.6040E-01  3.0016E-01
             2.4724E+00
 PARAMETER:  9.7530E-02 -9.3632E-02  1.0029E+00  2.2466E-01  3.0481E-01  1.7825E-01  3.8364E-01  1.0948E+00 -6.7566E-01 -1.1034E+00
             1.0052E+00
 GRADIENT:   1.9798E+01 -1.0531E+01 -1.7233E+00 -1.4429E+01  8.2300E+00 -3.0802E+00  1.2906E+00 -8.6033E-01  3.4541E+00  2.1986E-01
             4.1914E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1457.07320954114        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  9.9755E-01  8.2397E-01  2.4671E+00  1.1328E+00  1.2272E+00  1.0891E+00  1.3280E+00  2.7038E+00  4.1336E-01  2.1785E-01
             2.4721E+00
 PARAMETER:  9.7543E-02 -9.3619E-02  1.0030E+00  2.2469E-01  3.0477E-01  1.8537E-01  3.8369E-01  1.0946E+00 -7.8344E-01 -1.4240E+00
             1.0051E+00
 GRADIENT:   1.9413E+01 -9.2245E+00 -1.9872E+00 -5.5797E+00  7.0825E+00 -4.7814E-01 -2.5700E+00 -9.8612E-01  1.1602E+00  1.0212E-01
            -8.5004E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1457.10327654658        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:      570
 NPARAMETR:  9.9760E-01  8.2401E-01  2.5529E+00  1.1329E+00  1.2271E+00  1.0890E+00  1.3283E+00  2.7023E+00  4.1319E-01  2.1656E-01
             2.4708E+00
 PARAMETER:  9.7593E-02 -9.3569E-02  1.0372E+00  2.2480E-01  3.0462E-01  1.8529E-01  3.8388E-01  1.0941E+00 -7.8384E-01 -1.4299E+00
             1.0046E+00
 GRADIENT:   1.9756E+01 -7.9268E+00  5.5526E-01 -6.2914E+00  2.6989E+00 -4.0032E-01 -2.3262E+00 -3.1847E+00  1.2126E+00  1.0183E-01
            -1.4741E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1457.36201768102        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  9.8836E-01  8.3238E-01  2.5525E+00  1.1337E+00  1.2252E+00  1.0946E+00  1.3660E+00  2.7911E+00  3.8353E-01  1.3350E-01
             2.4815E+00
 PARAMETER:  8.8296E-02 -8.3472E-02  1.0371E+00  2.2547E-01  3.0311E-01  1.9036E-01  4.1189E-01  1.1264E+00 -8.5834E-01 -1.9137E+00
             1.0089E+00
 GRADIENT:   1.7328E+00  3.4804E+00 -1.1163E+00  1.5324E+01 -3.9563E+00  1.4110E+00 -3.5141E-01  5.7778E-01  7.9002E-01  3.6155E-02
            -2.7083E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1457.38352943795        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:      833
 NPARAMETR:  9.8823E-01  8.3242E-01  2.5548E+00  1.1334E+00  1.2254E+00  1.0937E+00  1.3672E+00  2.7855E+00  3.5713E-01  8.1782E-02
             2.4805E+00
 PARAMETER:  8.8159E-02 -8.3418E-02  1.0380E+00  2.2524E-01  3.0325E-01  1.8955E-01  4.1280E-01  1.1244E+00 -9.2966E-01 -2.4037E+00
             1.0084E+00
 GRADIENT:   1.3003E+00  4.4847E+00 -1.0967E+00  2.1035E+01 -4.7894E+00  1.0256E+00 -2.4122E+00  2.4507E-01 -1.7719E-01  1.2897E-02
            -2.7471E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1457.39328194786        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      990
 NPARAMETR:  9.8823E-01  8.3242E-01  2.5548E+00  1.1334E+00  1.2254E+00  1.0961E+00  1.3672E+00  2.7926E+00  3.6884E-01  4.0945E-02
             2.4805E+00
 PARAMETER:  8.8159E-02 -8.3418E-02  1.0380E+00  2.2524E-01  3.0325E-01  1.9177E-01  4.1280E-01  1.1270E+00 -8.9739E-01 -3.0955E+00
             1.0084E+00
 GRADIENT:  -5.8018E+00  3.6160E+00 -1.3383E+00  1.4455E+01 -4.6993E+00  1.8745E-01 -1.6566E+00  1.7928E-01  6.7667E-02  2.8421E-03
            -3.1427E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1457.39502629880        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1165
 NPARAMETR:  9.8823E-01  8.3242E-01  2.5548E+00  1.1334E+00  1.2254E+00  1.0955E+00  1.3673E+00  2.7887E+00  3.6706E-01  1.0000E-02
             2.4805E+00
 PARAMETER:  8.8159E-02 -8.3418E-02  1.0380E+00  2.2524E-01  3.0325E-01  1.9122E-01  4.1280E-01  1.1256E+00 -9.0223E-01 -4.8867E+00
             1.0085E+00
 GRADIENT:  -5.8192E+00  3.7154E+00 -1.2629E+00  1.4851E+01 -4.7036E+00 -1.9907E-02 -1.8057E+00  1.6836E-02 -4.1657E-04  0.0000E+00
            -3.3340E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1457.39786002686        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1352
 NPARAMETR:  9.8825E-01  8.3242E-01  2.5553E+00  1.1333E+00  1.2255E+00  1.0970E+00  1.3674E+00  2.7793E+00  3.6665E-01  1.0000E-02
             2.4816E+00
 PARAMETER:  8.8177E-02 -8.3421E-02  1.0382E+00  2.2513E-01  3.0331E-01  1.9257E-01  4.1290E-01  1.1222E+00 -9.0334E-01 -6.2205E+01
             1.0089E+00
 GRADIENT:  -5.7688E+00  3.5598E+00 -1.0421E+00  1.4283E+01 -4.5302E+00  4.8955E-01 -1.7844E+00 -4.0283E-01  6.5750E-03  0.0000E+00
            -2.9392E+00

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1457.39786529840        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:     1435
 NPARAMETR:  9.8826E-01  8.3242E-01  2.5550E+00  1.1333E+00  1.2254E+00  1.0970E+00  1.3674E+00  2.7825E+00  3.6668E-01  1.0000E-02
             2.4814E+00
 PARAMETER:  8.8177E-02 -8.3421E-02  1.0382E+00  2.2513E-01  3.0331E-01  1.9257E-01  4.1290E-01  1.1222E+00 -9.0335E-01 -6.2292E+01
             1.0089E+00
 GRADIENT:  -3.3079E+04 -6.6149E+04  6.3745E+03 -2.9366E+04  2.1806E+04  1.7177E+04 -8.0126E+03 -3.7073E-01 -7.3227E+03  0.0000E+00
             6.5685E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1435
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.2708E-03  1.3883E-02 -3.5921E-02 -1.9314E-02 -2.1588E-04
 SE:             2.9474E-02  2.2758E-02  1.5233E-02  1.2106E-02  1.6262E-04
 N:                     100         100         100         100         100

 P VAL.:         8.3152E-01  5.4185E-01  1.8372E-02  1.1062E-01  1.8433E-01

 ETASHRINKSD(%)  1.2571E+00  2.3757E+01  4.8966E+01  5.9444E+01  9.9455E+01
 ETASHRINKVR(%)  2.4984E+00  4.1870E+01  7.3956E+01  8.3552E+01  9.9997E+01
 EBVSHRINKSD(%)  1.6969E+00  2.4963E+01  5.8616E+01  5.9884E+01  9.9407E+01
 EBVSHRINKVR(%)  3.3650E+00  4.3695E+01  8.2874E+01  8.3907E+01  9.9996E+01
 RELATIVEINF(%)  9.5938E+01  4.2006E+00  4.6589E+00  1.2365E+00  9.8909E-04
 EPSSHRINKSD(%)  2.7266E+01
 EPSSHRINKVR(%)  4.7098E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1457.3978652983958     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -722.24703873465762     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.44
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.68
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1457.398       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.88E-01  8.32E-01  2.56E+00  1.13E+00  1.23E+00  1.10E+00  1.37E+00  2.78E+00  3.67E-01  1.00E-02  2.48E+00
 


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
+        1.69E+08
 
 TH 2
+       -6.22E+03  2.39E+08
 
 TH 3
+        1.95E+02 -1.03E+03  2.35E+05
 
 TH 4
+       -2.09E+03  1.14E+04 -3.02E+03  2.54E+07
 
 TH 5
+        1.38E+03 -7.54E+03  2.05E+03 -2.98E+03  1.20E+07
 
 TH 6
+       -6.10E+03 -7.25E+03  2.28E+02 -2.38E+03  1.62E+03  3.71E+07
 
 TH 7
+        5.70E+03  6.79E+03 -2.11E+02  2.17E+03 -1.51E+03 -1.07E+03  5.19E+06
 
 TH 8
+        4.82E+04  5.72E+04 -1.80E+03  1.87E+04 -1.28E+04 -2.52E+06  8.44E+03  7.57E+00
 
 TH 9
+        2.11E+03  2.48E+03 -7.73E+01  7.08E+02 -5.45E+02 -1.82E+03  4.00E+02  1.44E+04  1.51E+07
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.97E+02 -1.11E+03  3.06E+02 -1.20E+04  8.19E+03  2.43E+02 -2.16E+02 -1.90E+03 -7.11E+01  0.00E+00  2.65E+05
 
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
 #CPUT: Total CPU Time in Seconds,       25.178
Stop Time:
Sat Sep 25 07:51:34 CDT 2021
