Sat Sep 18 09:58:29 CDT 2021
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
$DATA ../../../../data/spa/A2/dat63.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m63.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1192.44380041210        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -4.2648E+01 -8.2172E+01 -1.8638E+01 -1.3977E+02  6.4422E+01 -2.4933E+01 -1.7600E+01  8.1110E+00 -7.1042E+01 -1.8999E+00
            -8.7611E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1463.93665185245        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0665E+00  1.0077E+00  1.2210E+00  1.0940E+00  1.0091E+00  1.0527E+00  9.8301E-01  8.7500E-01  1.2190E+00  7.3221E-01
             1.8473E+00
 PARAMETER:  1.6438E-01  1.0765E-01  2.9965E-01  1.8982E-01  1.0909E-01  1.5138E-01  8.2865E-02 -3.3535E-02  2.9802E-01 -2.1169E-01
             7.1372E-01
 GRADIENT:   6.1102E+01  8.0857E+00  2.9984E+01 -1.3384E+01 -2.4297E+01  2.1923E+00  2.2228E+00 -1.5585E+00  7.0986E+00 -5.1291E+00
            -1.2164E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1472.86320268078        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0666E+00  1.0713E+00  7.5776E-01  1.0499E+00  8.4190E-01  1.0495E+00  1.1624E+00  4.6333E-01  1.1109E+00  5.3448E-01
             1.9526E+00
 PARAMETER:  1.6449E-01  1.6891E-01 -1.7739E-01  1.4870E-01 -7.2098E-02  1.4836E-01  2.5045E-01 -6.6931E-01  2.0519E-01 -5.2646E-01
             7.6917E-01
 GRADIENT:   4.8565E+01  5.2885E+00 -8.5017E+00  9.1688E+00  1.4846E+01  7.7348E-01  3.1209E+00  1.3515E+00 -1.5782E+00  1.7596E+00
            -7.4467E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1480.24649785654        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0439E+00  7.8863E-01  8.8029E-01  1.2247E+00  7.9116E-01  1.0373E+00  1.3342E+00  3.0424E-01  9.9583E-01  4.6275E-01
             2.2671E+00
 PARAMETER:  1.4300E-01 -1.3746E-01 -2.7506E-02  3.0267E-01 -1.3426E-01  1.3661E-01  3.8834E-01 -1.0899E+00  9.5817E-02 -6.7057E-01
             9.1849E-01
 GRADIENT:  -1.0634E-01  3.1403E+00 -5.6465E-01  1.4396E+00 -3.0328E+00  8.8156E-01 -3.1430E-01  5.8297E-01 -1.1458E+00  1.6679E+00
             3.5381E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1481.12992413547        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.0415E+00  6.2045E-01  9.4322E-01  1.3285E+00  7.7133E-01  1.0300E+00  1.5781E+00  1.5320E-01  9.4476E-01  4.2201E-01
             2.3005E+00
 PARAMETER:  1.4063E-01 -3.7731E-01  4.1544E-02  3.8403E-01 -1.5964E-01  1.2959E-01  5.5621E-01 -1.7760E+00  4.3176E-02 -7.6273E-01
             9.3311E-01
 GRADIENT:  -5.8984E-01  1.9517E+00  1.1733E-01  5.5660E+00  1.3335E+00 -4.1867E-01 -2.2441E-01  9.0696E-02 -5.3644E-01 -3.4093E-01
             1.6507E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1481.63665524928        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.0398E+00  4.5914E-01  8.2892E-01  1.3982E+00  6.6234E-01  1.0323E+00  2.0026E+00  3.0823E-02  8.9872E-01  4.7944E-01
             2.2327E+00
 PARAMETER:  1.3908E-01 -6.7840E-01 -8.7628E-02  4.3516E-01 -3.1198E-01  1.3175E-01  7.9444E-01 -3.3795E+00 -6.7823E-03 -6.3513E-01
             9.0320E-01
 GRADIENT:  -5.1052E-01  1.6051E+00  1.3376E+00  3.3541E+00 -2.6885E+00  1.1431E-01  1.2045E-01  7.2531E-03 -3.0582E-01  1.9476E-02
            -4.1170E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1481.76289897961        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      490
 NPARAMETR:  1.0442E+00  4.0560E-01  8.1319E-01  1.4282E+00  6.4095E-01  1.0334E+00  2.1927E+00  1.4899E-02  8.8489E-01  4.8310E-01
             2.2326E+00
 PARAMETER:  1.4329E-01 -8.0238E-01 -1.0679E-01  4.5640E-01 -3.4480E-01  1.3285E-01  8.8513E-01 -4.1064E+00 -2.2288E-02 -6.2753E-01
             9.0315E-01
 GRADIENT:  -1.8984E+00  2.9073E-01  3.8864E-02 -1.0466E+00 -1.9367E+00 -2.3007E-01 -2.8437E-01  1.7913E-03 -1.6671E-01 -1.5432E-01
             1.7150E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1482.04996346065        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      667
 NPARAMETR:  1.0383E+00  3.2089E-01  9.3627E-01  1.5055E+00  6.8110E-01  1.0208E+00  2.2825E+00  1.0000E-02  8.7524E-01  5.3428E-01
             2.2815E+00
 PARAMETER:  1.3757E-01 -1.0367E+00  3.4147E-02  5.0910E-01 -2.8404E-01  1.2057E-01  9.2527E-01 -5.7893E+00 -3.3263E-02 -5.2683E-01
             9.2485E-01
 GRADIENT:  -8.6009E+00  3.2096E+00  7.8255E+00  1.2112E+01 -1.2457E+01 -3.3127E+00 -9.3278E-01  0.0000E+00  4.1373E-01 -3.0681E-01
             3.0031E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1482.73938050246        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      844
 NPARAMETR:  1.0410E+00  1.5674E-01  8.5247E-01  1.5776E+00  6.0666E-01  1.0301E+00  3.3490E+00  1.0000E-02  8.5689E-01  5.7625E-01
             2.2228E+00
 PARAMETER:  1.4015E-01 -1.7531E+00 -5.9613E-02  5.5593E-01 -3.9979E-01  1.2962E-01  1.3087E+00 -1.1315E+01 -5.4448E-02 -4.5122E-01
             8.9878E-01
 GRADIENT:   4.9958E+00  5.4157E+00 -2.4521E+00  3.9242E+00  7.6056E-02 -1.3190E-01  8.5424E+00  0.0000E+00 -1.3702E+00 -1.6314E+00
            -4.6028E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1485.87497859345        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1021
 NPARAMETR:  1.0333E+00  1.4670E-02  6.9075E-01  1.5893E+00  5.0400E-01  1.0340E+00  9.6893E+00  1.0000E-02  8.1054E-01  5.6411E-01
             2.1670E+00
 PARAMETER:  1.3280E-01 -4.1220E+00 -2.6998E-01  5.6330E-01 -5.8517E-01  1.3344E-01  2.3710E+00 -3.2223E+01 -1.1006E-01 -4.7250E-01
             8.7336E-01
 GRADIENT:  -2.3222E+00 -2.7131E+00  2.4978E+00  6.0033E+00 -1.6097E+00  1.9723E+00 -7.5771E+00  0.0000E+00 -5.9682E+00  4.1650E+00
             1.3620E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1486.66165909959        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1200
 NPARAMETR:  1.0353E+00  1.0000E-02  6.5919E-01  1.5824E+00  4.8704E-01  1.0315E+00  1.2147E+01  1.0000E-02  8.2340E-01  5.2200E-01
             2.1739E+00
 PARAMETER:  1.3472E-01 -4.6420E+00 -3.1674E-01  5.5894E-01 -6.1942E-01  1.3102E-01  2.5971E+00 -3.6868E+01 -9.4316E-02 -5.5009E-01
             8.7651E-01
 GRADIENT:   1.9853E+00  0.0000E+00  4.9811E-01  8.4899E+00 -6.0450E-01  9.7635E-01 -3.4461E+00  0.0000E+00 -2.0371E+00  7.1539E-01
            -4.3118E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1486.70336891394        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1378
 NPARAMETR:  1.0341E+00  1.0000E-02  6.4002E-01  1.5722E+00  4.7626E-01  1.0291E+00  1.2128E+01  1.0000E-02  8.3140E-01  5.1845E-01
             2.1782E+00
 PARAMETER:  1.3349E-01 -4.6548E+00 -3.4625E-01  5.5245E-01 -6.4179E-01  1.2873E-01  2.5955E+00 -3.7015E+01 -8.4642E-02 -5.5690E-01
             8.7849E-01
 GRADIENT:  -3.3220E-01  0.0000E+00  6.8785E-01 -5.3239E-01 -7.1342E-01 -5.1780E-02 -1.6958E-01  0.0000E+00  1.3389E-01  3.9017E-02
             2.7734E-01

0ITERATION NO.:   58    OBJECTIVE VALUE:  -1486.70371141347        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:     1501
 NPARAMETR:  1.0342E+00  1.0000E-02  6.3764E-01  1.5719E+00  4.7501E-01  1.0294E+00  1.2115E+01  1.0000E-02  8.3124E-01  5.1876E-01
             2.1772E+00
 PARAMETER:  1.3364E-01 -4.6559E+00 -3.5002E-01  5.5220E-01 -6.4430E-01  1.2895E-01  2.5949E+00 -3.7024E+01 -8.4736E-02 -5.5631E-01
             8.7789E-01
 GRADIENT:  -8.3107E-02  0.0000E+00 -8.4957E-01 -1.8352E+03  1.5732E+03  9.9006E-03  3.8775E+02  0.0000E+00  8.6822E-02 -8.4608E-03
            -1.1549E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1501
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.3982E-05  1.3053E-02  4.3010E-05 -1.0289E-02 -7.7950E-03
 SE:             2.9379E-02  7.1295E-03  2.3968E-04  2.7561E-02  1.7808E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9908E-01  6.7119E-02  8.5758E-01  7.0891E-01  6.6158E-01

 ETASHRINKSD(%)  1.5779E+00  7.6115E+01  9.9197E+01  7.6657E+00  4.0341E+01
 ETASHRINKVR(%)  3.1309E+00  9.4295E+01  9.9994E+01  1.4744E+01  6.4409E+01
 EBVSHRINKSD(%)  1.6474E+00  8.1571E+01  9.9137E+01  7.5491E+00  3.9967E+01
 EBVSHRINKVR(%)  3.2677E+00  9.6604E+01  9.9993E+01  1.4528E+01  6.3960E+01
 RELATIVEINF(%)  9.6451E+01  2.9841E+00  2.7378E-04  4.3425E+01  1.3263E+00
 EPSSHRINKSD(%)  3.2221E+01
 EPSSHRINKVR(%)  5.4060E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1486.7037114134662     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -751.55288484972800     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.11
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.59
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1486.704       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.00E-02  6.38E-01  1.57E+00  4.75E-01  1.03E+00  1.21E+01  1.00E-02  8.31E-01  5.19E-01  2.18E+00
 


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
+        1.01E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -8.23E+06  0.00E+00  5.09E+06
 
 TH 4
+        1.47E+02  0.00E+00 -2.26E+03  3.37E+05
 
 TH 5
+       -4.47E+02  0.00E+00  2.72E+03 -9.53E+05  2.71E+06
 
 TH 6
+       -3.31E+01  0.00E+00  8.57E+06 -8.19E+01  2.20E+02  1.82E+02
 
 TH 7
+       -4.10E+00  0.00E+00  5.38E+01  3.12E+01 -1.07E+02  1.75E+00  2.56E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.24E+02  0.00E+00  1.37E+07 -5.10E+00  2.48E+01  6.24E+01 -1.91E+00  0.00E+00  5.07E+02
 
 TH10
+       -4.29E+02  0.00E+00  3.94E+06 -8.53E+02  2.53E+03  1.83E+02  2.02E+01  0.00E+00  1.03E+03  3.05E+06
 
 TH11
+        6.92E+01  0.00E+00 -9.71E+02  1.53E+05 -4.33E+05 -2.96E+01  9.82E+00  0.00E+00  1.29E+01 -3.45E+02  6.94E+04
 
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
 #CPUT: Total CPU Time in Seconds,       24.770
Stop Time:
Sat Sep 18 09:58:55 CDT 2021
