Fri Sep 24 22:24:10 CDT 2021
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
$DATA ../../../../data/int/A3/dat49.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -627.255815771601        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.1365E+02  9.6178E+01  3.2755E+02 -9.8568E+01  1.9380E+02  1.0400E-01 -2.4896E+02 -4.0089E+02 -1.0413E+02 -1.7416E+02
            -5.4223E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2767.67092271825        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.3276E-01  9.8925E-01  8.3552E-01  1.0209E+00  9.1482E-01  9.3460E-01  1.6119E+00  8.9185E-01  1.1588E+00  1.0718E+00
             2.3475E+00
 PARAMETER:  3.0389E-02  8.9190E-02 -7.9701E-02  1.2071E-01  1.0971E-02  3.2359E-02  5.7741E-01 -1.4461E-02  2.4734E-01  1.6934E-01
             9.5333E-01
 GRADIENT:  -1.3302E+00 -1.4604E+00 -7.8404E+00 -2.6735E+01  1.8932E+01 -1.5073E+01  1.2051E+01  1.7243E+01  2.7727E+01  7.1475E+00
            -4.2522E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2790.03372782506        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.4212E-01  5.7618E-01  4.9425E-01  1.3318E+00  4.8604E-01  9.7456E-01  1.9178E+00  2.3949E-01  1.0475E+00  7.8979E-01
             2.2685E+00
 PARAMETER:  4.0383E-02 -4.5134E-01 -6.0472E-01  3.8655E-01 -6.2146E-01  7.4229E-02  7.5120E-01 -1.3292E+00  1.4640E-01 -1.3598E-01
             9.1912E-01
 GRADIENT:   2.0850E+01  5.2612E+01  3.8412E+01  1.6644E+02 -1.7767E+01 -1.1466E+00  3.0922E+01  2.1467E+00 -2.9397E+01 -1.7589E+00
            -3.7829E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2812.26112961066        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.2785E-01  3.5979E-01  2.5316E-01  1.4399E+00  2.6914E-01  9.9127E-01  1.5575E+00  2.0129E-02  1.3207E+00  8.0620E-01
             2.1972E+00
 PARAMETER:  2.5110E-02 -9.2224E-01 -1.2737E+00  4.6454E-01 -1.2125E+00  9.1231E-02  5.4305E-01 -3.8056E+00  3.7813E-01 -1.1543E-01
             8.8717E-01
 GRADIENT:  -1.2254E+01  1.1049E+02  4.5454E+01  2.5435E+02 -1.0892E+02  1.9412E+00 -2.5249E+00  7.6082E-03 -4.6700E+01 -1.0199E+01
             4.7851E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2813.06085959674        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      367
 NPARAMETR:  9.3478E-01  3.5789E-01  2.5179E-01  1.4383E+00  2.6795E-01  9.8701E-01  1.5731E+00  1.0000E-02  1.3221E+00  8.4560E-01
             2.1964E+00
 PARAMETER:  3.2554E-02 -9.2753E-01 -1.2792E+00  4.6346E-01 -1.2169E+00  8.6921E-02  5.5307E-01 -5.5083E+00  3.7919E-01 -6.7712E-02
             8.8683E-01
 GRADIENT:   3.8023E+00  1.0673E+02  3.9887E+01  2.5422E+02 -9.8541E+01  6.4295E-01  1.1600E+00  0.0000E+00 -4.6593E+01  8.1740E-01
             5.0796E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2813.09695493874        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  9.3980E-01  3.5787E-01  2.5178E-01  1.4382E+00  2.6802E-01  9.8317E-01  1.5619E+00  1.0000E-02  1.3221E+00  8.3977E-01
             2.1964E+00
 PARAMETER:  3.7907E-02 -9.2760E-01 -1.2792E+00  4.6339E-01 -1.2167E+00  8.3026E-02  5.4588E-01 -5.4748E+00  3.7920E-01 -7.4632E-02
             8.8680E-01
 GRADIENT:   8.9167E+00  1.0126E+02  3.2869E+01  2.4391E+02 -1.3240E+02 -1.6154E+00 -1.2426E+00  0.0000E+00 -4.9845E+01 -1.0978E+00
             4.8141E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2821.51189863539        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      728
 NPARAMETR:  9.3575E-01  3.5174E-01  2.4883E-01  1.4115E+00  2.8763E-01  9.8584E-01  1.5507E+00  1.0000E-02  1.3245E+00  8.3023E-01
             2.1782E+00
 PARAMETER:  3.3598E-02 -9.4486E-01 -1.2910E+00  4.4465E-01 -1.1461E+00  8.5743E-02  5.3870E-01 -1.8706E+01  3.8105E-01 -8.6057E-02
             8.7848E-01
 GRADIENT:   2.7561E-01  6.9376E+01 -8.8579E+01  2.1916E+02  6.4131E+01 -1.0488E+00  5.8116E-01  0.0000E+00 -4.1846E+01 -2.9585E-01
             3.5773E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2821.77910312050        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      910
 NPARAMETR:  9.2787E-01  3.5083E-01  2.4883E-01  1.4072E+00  2.8945E-01  1.0099E+00  1.5301E+00  1.0000E-02  1.3249E+00  8.2675E-01
             2.1754E+00
 PARAMETER:  2.5133E-02 -9.4745E-01 -1.2910E+00  4.4159E-01 -1.1398E+00  1.0981E-01  5.2533E-01 -1.8927E+01  3.8133E-01 -9.0256E-02
             8.7723E-01
 GRADIENT:  -1.7325E+01  6.6333E+01 -9.7899E+01  2.1489E+02  8.0725E+01  7.6153E+00 -2.0492E+00  0.0000E+00 -4.0488E+01 -2.2462E+00
             3.3237E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2841.41559622222        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1090
 NPARAMETR:  9.3482E-01  3.2176E-01  2.6999E-01  1.2563E+00  2.9320E-01  9.8162E-01  1.5440E+00  1.0000E-02  1.3375E+00  8.4216E-01
             2.0862E+00
 PARAMETER:  3.2599E-02 -1.0339E+00 -1.2094E+00  3.2820E-01 -1.1269E+00  8.1447E-02  5.3434E-01 -6.1891E+01  3.9082E-01 -7.1791E-02
             8.3536E-01
 GRADIENT:   9.9974E-01 -1.4664E+01  4.8914E+01  7.0784E+01  1.4650E+01 -5.3782E-01 -1.9451E-01  0.0000E+00  5.1571E+00  3.5744E-01
            -6.0271E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2841.43003407282        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1270
 NPARAMETR:  9.3539E-01  3.2191E-01  2.6994E-01  1.2569E+00  2.9288E-01  9.7998E-01  1.5435E+00  1.0000E-02  1.3375E+00  8.4513E-01
             2.0869E+00
 PARAMETER:  3.3208E-02 -1.0335E+00 -1.2096E+00  3.2863E-01 -1.1280E+00  7.9780E-02  5.3408E-01 -6.1688E+01  3.9077E-01 -6.8261E-02
             8.3567E-01
 GRADIENT:   2.3190E+00 -1.3956E+01  5.0176E+01  7.1430E+01  1.2080E+01 -1.1530E+00 -8.9788E-02  0.0000E+00  5.0722E+00  1.1923E+00
            -5.9520E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2847.74463251334        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1452
 NPARAMETR:  9.3707E-01  3.1695E-01  2.4859E-01  1.1608E+00  2.7520E-01  9.7721E-01  1.5108E+00  1.0000E-02  1.3420E+00  8.4755E-01
             2.1876E+00
 PARAMETER:  3.5004E-02 -1.0490E+00 -1.2920E+00  2.4912E-01 -1.1903E+00  7.6951E-02  5.1267E-01 -6.8017E+01  3.9416E-01 -6.5403E-02
             8.8283E-01
 GRADIENT:   3.1767E+00  1.5143E+01  3.1250E+01 -8.3913E+00  4.5749E+00 -1.5675E+00  2.7029E-01  0.0000E+00  4.7381E+00  2.0377E+00
             5.6659E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2850.26090843240        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1587
 NPARAMETR:  9.3365E-01  2.9821E-01  2.3346E-01  1.1700E+00  2.6651E-01  9.8139E-01  1.5173E+00  1.0000E-02  1.3090E+00  8.4480E-01
             2.0636E+00
 PARAMETER:  3.1351E-02 -1.1100E+00 -1.3548E+00  2.5700E-01 -1.2223E+00  8.1214E-02  5.1695E-01 -6.8252E+01  3.6930E-01 -6.8658E-02
             8.2448E-01
 GRADIENT:   5.2892E+00  7.6349E+00  2.1252E+01  2.1903E+01  7.9195E+01 -5.9724E-02 -7.3993E-01  0.0000E+00 -8.3930E+00 -1.8457E+00
            -5.4961E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2851.41029501000        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:     1690
 NPARAMETR:  9.3260E-01  2.7802E-01  2.1381E-01  1.1599E+00  2.5022E-01  9.8276E-01  1.5250E+00  1.0000E-02  1.3196E+00  8.5322E-01
             2.0025E+00
 PARAMETER:  3.0221E-02 -1.1801E+00 -1.4427E+00  2.4831E-01 -1.2854E+00  8.2613E-02  5.2201E-01 -6.8252E+01  3.7734E-01 -5.8732E-02
             7.9442E-01
 GRADIENT:  -3.8077E+00 -7.9876E+00  1.0466E+01  2.6772E+01  4.1321E+01 -1.0701E+00 -1.6891E+00  0.0000E+00 -1.7359E+01 -2.0019E+00
            -1.0888E+02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2851.77167760803        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1858
 NPARAMETR:  9.3425E-01  2.8095E-01  2.1389E-01  1.1342E+00  2.5030E-01  9.8475E-01  1.5318E+00  1.0000E-02  1.3195E+00  8.6165E-01
             2.0030E+00
 PARAMETER:  3.1991E-02 -1.1696E+00 -1.4423E+00  2.2589E-01 -1.2851E+00  8.4634E-02  5.2644E-01 -6.8252E+01  3.7725E-01 -4.8903E-02
             7.9462E-01
 GRADIENT:   1.4071E-01  1.5216E-01  1.2213E+01 -6.9314E-02  3.9285E+01 -5.9711E-02 -1.4000E-01  0.0000E+00 -1.6402E+01  9.2145E-02
            -1.0741E+02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2854.97758562279        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2033             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3357E-01  2.7983E-01  2.1193E-01  1.1330E+00  2.4633E-01  9.8422E-01  1.5283E+00  1.0000E-02  1.3588E+00  8.5942E-01
             2.0681E+00
 PARAMETER:  3.1265E-02 -1.1736E+00 -1.4515E+00  2.2483E-01 -1.3011E+00  8.4092E-02  5.2416E-01 -6.8252E+01  4.0657E-01 -5.1498E-02
             8.2665E-01
 GRADIENT:   4.6034E+00  1.2427E+01  3.1064E+01  5.3154E+00  6.1077E+01  8.1640E-01  2.0825E+00  0.0000E+00 -3.5411E+00  1.5665E+00
            -3.2696E+01

0ITERATION NO.:   72    OBJECTIVE VALUE:  -2854.97758562279        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     2104
 NPARAMETR:  9.3358E-01  2.7983E-01  2.1194E-01  1.1329E+00  2.4634E-01  9.8422E-01  1.5275E+00  1.0000E-02  1.3587E+00  8.5942E-01
             2.0682E+00
 PARAMETER:  3.1265E-02 -1.1736E+00 -1.4515E+00  2.2483E-01 -1.3011E+00  8.4092E-02  5.2416E-01 -6.8252E+01  4.0657E-01 -5.1498E-02
             8.2665E-01
 GRADIENT:  -4.2712E+05  3.8408E+00 -5.8805E+04  1.8997E+05 -6.5634E+04 -4.2712E+05  1.1743E+00  0.0000E+00  2.1009E+05 -8.5424E+05
            -5.1784E+04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2104
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1369E-03  1.1403E-03 -5.2563E-05 -4.8121E-03 -5.4682E-03
 SE:             2.9488E-02  2.5136E-02  2.5914E-04  2.8904E-02  2.6882E-02
 N:                     100         100         100         100         100

 P VAL.:         9.1528E-01  9.6382E-01  8.3926E-01  8.6777E-01  8.3881E-01

 ETASHRINKSD(%)  1.2100E+00  1.5792E+01  9.9132E+01  3.1685E+00  9.9422E+00
 ETASHRINKVR(%)  2.4054E+00  2.9090E+01  9.9992E+01  6.2366E+00  1.8896E+01
 EBVSHRINKSD(%)  1.2977E+00  1.3420E+01  9.9138E+01  3.5733E+00  1.0100E+01
 EBVSHRINKVR(%)  2.5785E+00  2.5039E+01  9.9993E+01  7.0189E+00  1.9180E+01
 RELATIVEINF(%)  9.7402E+01  1.8497E+01  6.1551E-04  6.2948E+01  7.0580E+00
 EPSSHRINKSD(%)  1.9324E+01
 EPSSHRINKVR(%)  3.4914E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2854.9775856227920     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1200.8882258543813     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    55.28
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2854.978       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.34E-01  2.80E-01  2.12E-01  1.13E+00  2.46E-01  9.84E-01  1.53E+00  1.00E-02  1.36E+00  8.59E-01  2.07E+00
 


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
+        2.45E+09
 
 TH 2
+       -6.97E+08  3.96E+08
 
 TH 3
+       -1.82E+04 -2.71E+03  2.25E+08
 
 TH 4
+        2.20E+04  2.55E+08  2.02E+04  3.29E+08
 
 TH 5
+       -1.75E+04 -2.46E+03 -2.17E+08  1.96E+04  2.08E+08
 
 TH 6
+       -4.73E+04 -5.99E+00 -1.43E+04  1.73E+04 -1.38E+04  2.20E+09
 
 TH 7
+       -3.65E+02 -8.12E+07 -1.13E+02  1.33E+02 -1.37E+02 -3.46E+02  4.68E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.02E+04  2.01E+02  1.63E+04 -1.14E+04  1.55E+04  7.99E+03  6.09E+01  0.00E+00  7.00E+07
 
 TH10
+       -2.93E+04 -4.71E+00 -8.79E+03  1.08E+04 -8.64E+03  2.52E+09 -3.91E+02  0.00E+00  4.97E+03  2.89E+09
 
 TH11
+       -3.30E+03 -8.06E+01 -4.06E+07  3.69E+03  3.90E+07 -2.58E+03 -1.18E+01  0.00E+00  2.91E+03 -1.59E+03  7.33E+06
 
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
 #CPUT: Total CPU Time in Seconds,       68.698
Stop Time:
Fri Sep 24 22:25:21 CDT 2021
