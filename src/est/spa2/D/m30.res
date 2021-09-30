Thu Sep 30 08:53:23 CDT 2021
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
$DATA ../../../../data/spa2/D/dat30.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   14633.4655003125        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.2022E+02  1.8972E+02  6.2417E+00 -3.1562E+01  1.1308E+02 -1.1283E+03 -6.5458E+02 -8.3714E+01 -1.2154E+03 -5.7995E+02
            -3.0338E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -795.374698021337        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2516E+00  1.1738E+00  9.9396E-01  1.5329E+00  1.1474E+00  2.5187E+00  1.8176E+00  9.6479E-01  1.8805E+00  1.1596E+00
             1.3526E+01
 PARAMETER:  3.2442E-01  2.6024E-01  9.3937E-02  5.2717E-01  2.3748E-01  1.0237E+00  6.9753E-01  6.4159E-02  7.3156E-01  2.4806E-01
             2.7046E+00
 GRADIENT:  -5.5271E+01 -4.7079E+01 -3.8632E+01 -1.6338E+01  4.7995E+01  5.2627E+01 -3.0992E+00  2.6915E-01 -1.2687E+01  1.5601E+01
             4.1769E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -902.888393307291        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      201
 NPARAMETR:  1.3578E+00  5.9643E+00  6.3788E+01  1.1055E+00  3.9176E+01  2.7478E+00  3.1863E+00  5.1146E-01  2.0719E+01  3.7852E+00
             1.0839E+01
 PARAMETER:  4.0585E-01  1.8858E+00  4.2556E+00  2.0033E-01  3.7681E+00  1.1108E+00  1.2589E+00 -5.7049E-01  3.1310E+00  1.4311E+00
             2.4831E+00
 GRADIENT:  -6.1330E+00  5.6488E+01 -1.1795E-01  1.6328E+01  9.1369E-02  6.3875E+01  4.1479E+01  7.1693E-05  2.9344E+01 -3.9969E-04
             2.6534E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -925.500151855242        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  1.3733E+00  4.9963E+00  1.1848E+02  1.2805E+00  4.5412E+01  2.3955E+00  2.9717E+00  2.9845E-01  1.7517E+01  3.0478E+01
             9.9505E+00
 PARAMETER:  4.1720E-01  1.7087E+00  4.8747E+00  3.4728E-01  3.9158E+00  9.7359E-01  1.1891E+00 -1.1092E+00  2.9632E+00  3.5170E+00
             2.3976E+00
 GRADIENT:   1.6273E+01  3.7417E+01 -6.7089E-02  2.0556E+01 -1.1115E+00  4.3423E+01  4.7430E+01  4.0334E-06  3.1891E+01  1.2633E+01
             1.9349E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -934.666307072643        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:      506
 NPARAMETR:  1.3649E+00  4.9743E+00  1.2862E+02  9.1908E-01  4.5242E+01  2.3938E+00  2.9623E+00  2.9528E-01  1.4605E+01  3.0405E+01
             9.9446E+00
 PARAMETER:  4.1110E-01  1.7043E+00  4.9569E+00  1.5615E-02  3.9120E+00  9.7286E-01  1.1860E+00 -1.1198E+00  2.7814E+00  3.5146E+00
             2.3970E+00
 GRADIENT:   3.7326E+01  1.3280E+02 -5.1045E-02  2.1013E+01 -1.7777E+00  6.8645E+01  6.5143E+01  3.6308E-06  1.2227E+02  1.5713E+01
             2.1041E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -958.841845732974        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      670            RESET HESSIAN, TYPE II
 NPARAMETR:  1.3640E+00  4.8568E+00  1.3754E+02  8.4464E-01  5.0253E+01  2.3798E+00  2.9070E+00  2.9693E-01  1.4063E+01  1.4789E+01
             9.8021E+00
 PARAMETER:  4.1039E-01  1.6804E+00  5.0239E+00 -6.8841E-02  4.0171E+00  9.6700E-01  1.1671E+00 -1.1143E+00  2.7436E+00  2.7939E+00
             2.3826E+00
 GRADIENT:   3.9376E+01  1.3655E+02 -4.4557E-02  2.0953E+01  1.1930E-01  6.9269E+01  6.4712E+01  3.1964E-06  1.2440E+02  6.0215E-04
             2.0924E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1046.16297720430        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      787
 NPARAMETR:  1.2572E+00  2.9060E+00  1.0815E+02  2.5146E-01  1.0927E+02  2.0852E+00  1.8748E+00  2.9683E-01  8.4507E+00  1.0000E-02
             7.6805E+00
 PARAMETER:  3.2890E-01  1.1668E+00  4.7836E+00 -1.2805E+00  4.7938E+00  8.3488E-01  7.2848E-01 -1.1146E+00  2.2343E+00 -8.3930E+00
             2.1387E+00
 GRADIENT:   7.0320E+00 -1.8095E+01  1.0844E-01  3.8905E+00 -2.9597E-02 -1.6817E+01  1.7776E+01 -1.9575E-05  2.1238E+01  0.0000E+00
            -1.3449E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1048.23467416949        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      964
 NPARAMETR:  1.2478E+00  2.9221E+00  2.1175E+01  2.3320E-01  1.2392E+02  2.1611E+00  1.6829E+00  2.9721E-01  7.8522E+00  1.0000E-02
             7.7593E+00
 PARAMETER:  3.2138E-01  1.1723E+00  3.1528E+00 -1.3558E+00  4.9196E+00  8.7064E-01  6.2049E-01 -1.1133E+00  2.1608E+00 -9.2339E+00
             2.1489E+00
 GRADIENT:   9.4281E-01 -3.2682E+00 -4.3856E-02  3.9969E-01 -2.6493E-02 -2.4416E+00  3.2792E+00 -4.8094E-04  2.3982E+00  0.0000E+00
            -1.9205E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1048.55928010664        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1146
 NPARAMETR:  1.2469E+00  2.9372E+00  2.1288E+01  2.1849E-01  4.0584E+03  2.1830E+00  1.6394E+00  4.1698E-01  8.0149E+00  1.0000E-02
             7.7822E+00
 PARAMETER:  3.2067E-01  1.1775E+00  3.1581E+00 -1.4210E+00  8.4085E+00  8.8070E-01  5.9435E-01 -7.7471E-01  2.1813E+00 -9.2520E+00
             2.1518E+00
 GRADIENT:   4.6357E-01 -2.8670E+00  2.9067E-05  2.6171E-01 -7.3992E-04  1.3797E+00 -7.3566E-01 -7.9865E-04  4.3232E+00  0.0000E+00
             1.1207E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1048.61935432847        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1332
 NPARAMETR:  1.2448E+00  2.9638E+00  2.2358E+01  2.1615E-01  1.6224E+05  2.1923E+00  1.6479E+00  9.1094E-01  8.0882E+00  1.0000E-02
             7.7221E+00
 PARAMETER:  3.1900E-01  1.1865E+00  3.2072E+00 -1.4318E+00  1.2097E+01  8.8495E-01  5.9948E-01  6.7205E-03  2.1904E+00 -9.2520E+00
             2.1441E+00
 GRADIENT:   1.9673E-01  1.1522E+00  1.4217E-01  9.7625E-01 -1.4523E-05  2.2552E+00 -1.5254E+00 -3.1262E-03  5.4131E+00  0.0000E+00
            -9.5848E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1048.67839490381        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1522
 NPARAMETR:  1.2450E+00  2.9622E+00  2.2724E+01  2.1305E-01  9.5251E+05  2.1915E+00  1.6510E+00  5.0019E+00  8.0662E+00  1.0000E-02
             7.7649E+00
 PARAMETER:  3.1913E-01  1.1859E+00  3.2234E+00 -1.4462E+00  1.3867E+01  8.8458E-01  6.0138E-01  1.7098E+00  2.1877E+00 -9.2520E+00
             2.1496E+00
 GRADIENT:  -9.7881E-02 -5.9706E-01  1.6052E-01  6.8890E-01 -2.6998E-06  2.6841E+00 -1.0386E+00 -6.1577E-03  4.6754E+00  0.0000E+00
            -2.0072E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1048.81511984599        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1690
 NPARAMETR:  1.2459E+00  2.9643E+00  2.2560E+01  2.0832E-01  1.7817E+06  2.1811E+00  1.6608E+00  4.8948E+00  8.1678E+00  1.0000E-02
             7.7744E+00
 PARAMETER:  3.1983E-01  1.1866E+00  3.2162E+00 -1.4687E+00  1.4493E+01  8.7983E-01  6.0729E-01  1.6882E+00  2.2002E+00 -9.2520E+00
             2.1508E+00
 GRADIENT:   2.9009E+01  6.0642E+01  3.1263E-01  6.2201E+00 -1.4869E-06  4.1189E+01  7.7622E+00  1.9888E-02  1.5875E+02  0.0000E+00
             2.5571E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1048.88405320640        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1878             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2447E+00  2.9739E+00  2.2329E+01  2.0614E-01  1.8255E+06  2.1825E+00  1.6660E+00  4.8939E+00  8.2460E+00  1.0000E-02
             7.7329E+00
 PARAMETER:  3.1893E-01  1.1899E+00  3.2059E+00 -1.4792E+00  1.4517E+01  8.8045E-01  6.1045E-01  1.6880E+00  2.2097E+00 -9.2520E+00
             2.1455E+00
 GRADIENT:   2.9206E+01  6.2646E+01  2.3458E-01  6.6398E+00 -1.4650E-06  4.1498E+01  7.7947E+00  3.9631E-02  1.6259E+02  0.0000E+00
             1.8635E+01

0ITERATION NO.:   63    OBJECTIVE VALUE:  -1048.92296525532        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:     1957
 NPARAMETR:  1.2446E+00  2.9772E+00  2.2290E+01  2.0122E-01  2.1176E+06  2.1869E+00  1.6687E+00  4.8920E+00  8.2998E+00  1.0000E-02
             7.6979E+00
 PARAMETER:  3.1873E-01  1.1894E+00  3.2055E+00 -1.5071E+00  1.4521E+01  8.8031E-01  6.1052E-01  1.6879E+00  2.2105E+00 -9.2520E+00
             2.1460E+00
 GRADIENT:  -4.1206E-02 -1.8543E+00  2.5410E-02 -4.4759E+02 -7.6428E-06 -7.6597E+02 -5.5682E+02  5.3551E-04 -1.5406E+02  0.0000E+00
             3.1527E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1957
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.6165E-03 -8.8297E-02 -3.4485E-03  6.6446E-02 -1.2568E-11
 SE:             2.9373E-02  1.7827E-02  1.4964E-03  2.1148E-02  1.2453E-11
 N:                     100         100         100         100         100

 P VAL.:         7.4337E-01  7.3138E-07  2.1188E-02  1.6781E-03  3.1285E-01

 ETASHRINKSD(%)  1.5961E+00  4.0279E+01  9.4987E+01  2.9152E+01  1.0000E+02
 ETASHRINKVR(%)  3.1668E+00  6.4334E+01  9.9749E+01  4.9806E+01  1.0000E+02
 EBVSHRINKSD(%)  3.7139E+00  2.5910E+01  9.7611E+01  4.0464E+01  1.0000E+02
 EBVSHRINKVR(%)  7.2900E+00  4.5107E+01  9.9943E+01  6.4555E+01  1.0000E+02
 RELATIVEINF(%)  9.0023E+01  3.1503E+01  5.2211E-02  1.9491E+01  0.0000E+00
 EPSSHRINKSD(%)  1.1636E+01
 EPSSHRINKVR(%)  2.1918E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1048.9229652553240     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       53.803274590283081     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    56.06
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.64
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1048.923       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.24E+00  2.97E+00  2.23E+01  2.00E-01  1.83E+06  2.18E+00  1.67E+00  4.89E+00  8.25E+00  1.00E-02  7.74E+00
 


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
+        1.33E+02
 
 TH 2
+       -3.01E+00  3.38E+01
 
 TH 3
+       -3.89E-03 -6.52E-02  1.18E-02
 
 TH 4
+        1.72E+01  2.95E+02  3.56E-01  1.84E+05
 
 TH 5
+        1.03E-10 -4.26E-12 -6.78E-13  4.46E-10 -2.56E-18
 
 TH 6
+       -1.15E+01  2.50E+03  5.49E-02  2.89E+04 -1.48E-11  4.55E+03
 
 TH 7
+       -4.21E+04  4.75E+03  7.78E-02  5.51E+04 -5.06E-11  1.80E+01  1.65E+04
 
 TH 8
+        6.87E-03  4.46E-02 -9.73E-03  3.63E-01 -9.01E-12  7.51E-02  1.14E-01  2.34E-02
 
 TH 9
+       -1.27E+00  6.45E+02  5.91E-03 -2.59E+01  4.88E-12  4.61E+02  2.76E+00  7.57E-03  6.54E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.83E+00 -8.10E+00 -8.85E-03 -9.75E+02  7.27E-13 -1.41E+02  3.85E-01 -4.93E-03 -7.15E+01  0.00E+00  9.07E+01
 
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
 #CPUT: Total CPU Time in Seconds,       69.780
Stop Time:
Thu Sep 30 08:54:35 CDT 2021
