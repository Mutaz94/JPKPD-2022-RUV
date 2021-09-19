Sat Sep 18 06:43:00 CDT 2021
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
$DATA ../../../../data/int/D/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   12660.5656904942        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.0377E+02  3.9278E+01 -8.9887E+01 -1.6085E+02  1.6624E+02 -1.4685E+03 -5.8932E+02 -3.4840E+01 -1.0019E+03 -4.3980E+02
            -2.9760E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1122.27962499006        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  2.7049E+00  1.5015E+00  1.0779E+00  2.6768E+00  1.3809E+00  7.6676E+00  5.9466E+00  1.0204E+00  7.4749E+00  2.3543E+00
             9.7949E+00
 PARAMETER:  1.0951E+00  5.0649E-01  1.7498E-01  1.0846E+00  4.2276E-01  2.1370E+00  1.8828E+00  1.2017E-01  2.1115E+00  9.5624E-01
             2.3819E+00
 GRADIENT:   3.2361E+01 -2.7979E+01 -4.2473E+01  4.4958E+01 -3.9523E+01  1.6955E+02  1.0503E+02  2.3503E+00  1.1528E+02  6.1910E+01
             6.7406E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1270.40872666726        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      206
 NPARAMETR:  1.4268E+00  1.1788E+01  2.2405E+01  4.8268E+00  2.0478E+01  2.8126E+00  3.5732E+00  9.7818E-01  6.1254E+01  7.4570E+00
             9.1712E+00
 PARAMETER:  4.5546E-01  2.5671E+00  3.2093E+00  1.6742E+00  3.1194E+00  1.1341E+00  1.3735E+00  7.7939E-02  4.2150E+00  2.1092E+00
             2.3161E+00
 GRADIENT:   3.5122E+01  1.4650E+02 -3.9932E+00  5.6632E+00  1.9011E+01  8.8229E+01 -1.0096E+02 -9.8824E-02 -5.7782E+00  1.3721E+00
             5.1055E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1330.70471291477        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      390
 NPARAMETR:  1.1953E+00  1.0685E+01  4.1459E+01  7.3003E+00  2.3963E+01  2.0632E+00  3.4964E+00  7.9500E-01  8.0118E+01  9.7353E+00
             8.8880E+00
 PARAMETER:  2.7837E-01  2.4689E+00  3.8247E+00  2.0879E+00  3.2765E+00  8.2426E-01  1.3517E+00 -1.2942E-01  4.4835E+00  2.3758E+00
             2.2847E+00
 GRADIENT:  -1.2554E+01  1.2592E+02 -5.9665E+00  4.4532E+00  1.2655E+00  1.6517E+01 -6.9919E+01 -9.4923E-02 -2.3605E+00  1.0582E+01
             4.6381E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1439.31550681334        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      571
 NPARAMETR:  1.1844E+00  6.0732E+00  6.9642E+01  6.5807E+00  2.4301E+01  2.0874E+00  3.0115E+00  1.1670E+00  7.4835E+01  9.0616E+00
             6.6426E+00
 PARAMETER:  2.6922E-01  1.9039E+00  4.3434E+00  1.9841E+00  3.2905E+00  8.3593E-01  1.2025E+00  2.5444E-01  4.4153E+00  2.3040E+00
             1.9935E+00
 GRADIENT:   1.3891E+01  4.5094E+00 -1.2487E+00  5.9581E+00  8.6335E+00  1.8211E+01  1.3515E+01 -1.5409E-02  1.6627E+01  1.0216E+01
             5.7280E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1564.88686282223        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      757
 NPARAMETR:  1.2490E+00  1.9981E+00  7.5675E+01  1.0724E+00  1.3824E+01  2.0558E+00  1.6331E+00  1.1930E+01  7.1182E+00  2.4517E+00
             6.4171E+00
 PARAMETER:  3.2231E-01  7.9222E-01  4.4264E+00  1.6986E-01  2.7264E+00  8.2066E-01  5.9045E-01  2.5790E+00  2.0626E+00  9.9680E-01
             1.9590E+00
 GRADIENT:   4.3801E+01 -4.2626E+01 -9.3569E-01  1.7664E+01  1.3598E+01  7.2966E+00  2.9326E+01  1.6192E+00  6.1806E+01  2.1051E+00
            -8.6406E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1589.74788107700        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      932
 NPARAMETR:  1.1156E+00  1.8102E+00  6.5570E+01  9.5980E-01  1.0905E+01  1.8877E+00  1.4734E+00  1.0363E+00  4.3732E+00  1.4703E+00
             6.6575E+00
 PARAMETER:  2.0936E-01  6.9345E-01  4.2831E+00  5.8973E-02  2.4893E+00  7.3535E-01  4.8758E-01  1.3562E-01  1.5755E+00  4.8549E-01
             1.9957E+00
 GRADIENT:  -1.8360E+01 -1.3970E+01 -3.2328E+00 -6.6384E+00  3.3059E+01 -1.8105E+01  1.5193E+01 -6.4745E-02 -2.0082E+01 -2.9487E-01
             1.8630E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1778.32264235822        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1117
 NPARAMETR:  1.2263E+00  1.1525E+00  1.4117E+02  8.9970E-01  2.2138E+00  2.5611E+00  2.0771E-01  2.6546E+00  4.1933E+00  3.3253E-02
             5.8730E+00
 PARAMETER:  3.0402E-01  2.4191E-01  5.0500E+00 -5.6966E-03  8.9470E-01  1.0404E+00 -1.4716E+00  1.0763E+00  1.5335E+00 -3.3036E+00
             1.8704E+00
 GRADIENT:   1.9073E+01 -2.7907E+01 -3.2848E+00 -3.3537E+01 -8.8131E+01  4.7437E+01 -1.9107E+00 -2.9346E+00 -1.0132E+01 -4.0407E-02
             9.7118E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1806.28170644083        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1295
 NPARAMETR:  1.1319E+00  1.1991E+00  1.8291E+02  1.1498E+00  2.4219E+00  2.1860E+00  2.3243E-01  8.3891E+00  4.2799E+00  3.9531E-02
             5.6038E+00
 PARAMETER:  2.2386E-01  2.8154E-01  5.3090E+00  2.3962E-01  9.8454E-01  8.8207E-01 -1.3591E+00  2.2269E+00  1.5539E+00 -3.1307E+00
             1.8235E+00
 GRADIENT:  -4.2261E+00  6.9741E+00 -2.0979E+00 -1.7846E+00  5.5315E+00 -1.0016E+00 -3.2285E+00 -8.9871E-01 -1.4810E-01 -1.8220E-02
             8.5593E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1807.40791674948        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1475
 NPARAMETR:  1.1447E+00  1.1096E+00  3.6744E+02  1.2826E+00  2.4152E+00  2.1891E+00  2.1236E-01  1.1515E+01  4.0609E+00  3.9051E-02
             5.5934E+00
 PARAMETER:  2.3518E-01  2.0396E-01  6.0066E+00  3.4886E-01  9.8177E-01  8.8348E-01 -1.4495E+00  2.5436E+00  1.5014E+00 -3.1429E+00
             1.8216E+00
 GRADIENT:   1.2975E+00 -4.0793E+00 -1.7438E+00  1.8719E+00 -2.5494E-01 -1.4894E-01 -1.1382E+00  1.6555E+00  8.3520E-01 -1.0616E-02
             8.0526E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1807.92531640005        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1652
 NPARAMETR:  1.1339E+00  1.1660E+00  1.1243E+03  1.1960E+00  2.4224E+00  2.1871E+00  2.3751E-01  1.8629E+01  4.1737E+00  3.9701E-02
             5.5971E+00
 PARAMETER:  2.2565E-01  2.5359E-01  7.1249E+00  2.7899E-01  9.8476E-01  8.8257E-01 -1.3375E+00  3.0247E+00  1.5288E+00 -3.1264E+00
             1.8222E+00
 GRADIENT:  -3.0335E+00  2.0336E+00 -2.2248E-01 -9.0828E-01 -3.2723E-01 -5.0399E-01 -2.3687E+00  1.2190E+00 -1.0347E+00 -1.0991E-02
             2.8301E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1814.57150425823        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1832
 NPARAMETR:  1.1256E+00  1.4702E+00  3.3477E+04  9.5049E-01  2.4571E+00  2.1808E+00  7.4297E-01  5.3655E+01  4.7435E+00  4.1680E-02
             5.6278E+00
 PARAMETER:  2.1828E-01  4.8542E-01  1.0519E+01  4.9227E-02  9.9896E-01  8.7969E-01 -1.9710E-01  4.0826E+00  1.6568E+00 -3.0777E+00
             1.8277E+00
 GRADIENT:  -8.4614E+00 -8.5101E+00 -1.8288E-01  1.2264E+00  1.0209E+01 -1.6978E+00  1.2126E-01  1.6311E+00 -3.9349E-01 -2.7284E-03
             1.5551E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1814.65739965410        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:     1983
 NPARAMETR:  1.1412E+00  1.4725E+00  3.3385E+04  9.4555E-01  2.4562E+00  2.1813E+00  7.4496E-01  5.3670E+01  4.7509E+00  4.1736E-02
             5.6246E+00
 PARAMETER:  2.3211E-01  4.8697E-01  1.0516E+01  4.4016E-02  9.9862E-01  8.7990E-01 -1.9443E-01  4.0829E+00  1.6583E+00 -3.0764E+00
             1.8272E+00
 GRADIENT:   2.4206E+00 -6.9157E+00 -1.2338E-01  1.0651E+00  1.0509E+01  4.7749E+00  1.9933E-01  3.6228E+00  1.4205E+01 -2.6933E-03
             1.6446E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1814.66544530190        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2161
 NPARAMETR:  1.1436E+00  1.4725E+00  3.3390E+04  9.3085E-01  2.4561E+00  2.1813E+00  7.4238E-01  5.3664E+01  4.7509E+00  4.2829E-02
             5.6231E+00
 PARAMETER:  2.3416E-01  4.8699E-01  1.0516E+01  2.8341E-02  9.9856E-01  8.7991E-01 -1.9790E-01  4.0827E+00  1.6583E+00 -3.0505E+00
             1.8269E+00
 GRADIENT:  -1.1669E-01 -1.0493E+01 -1.8607E-01 -4.9415E-01  1.0114E+01 -1.5001E+00  1.5813E-01  1.6364E+00 -9.1584E-01 -3.2782E-03
             1.2866E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1814.69093485325        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2341
 NPARAMETR:  1.1429E+00  1.4727E+00  3.3434E+04  9.3199E-01  2.4547E+00  2.1813E+00  7.2060E-01  5.3604E+01  4.7510E+00  5.5260E-02
             5.6086E+00
 PARAMETER:  2.3359E-01  4.8709E-01  1.0517E+01  2.9565E-02  9.9801E-01  8.7992E-01 -2.2767E-01  4.0816E+00  1.6583E+00 -2.7957E+00
             1.8243E+00
 GRADIENT:   5.7244E-02 -7.8725E+00 -1.9035E-01 -1.3707E-01  9.0378E+00 -1.4861E+00 -7.7430E-01  1.6428E+00 -1.7847E+00 -7.9645E-03
             6.1957E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1814.74210358968        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2521
 NPARAMETR:  1.1429E+00  1.4732E+00  3.3569E+04  9.3114E-01  2.4503E+00  2.1815E+00  7.4475E-01  5.3414E+01  4.7518E+00  1.7584E-01
             5.5670E+00
 PARAMETER:  2.3359E-01  4.8741E-01  1.0521E+01  2.8657E-02  9.9622E-01  8.8000E-01 -1.9471E-01  4.0781E+00  1.6585E+00 -1.6382E+00
             1.8169E+00
 GRADIENT:  -7.5148E-01 -9.0640E+00 -1.9724E-01  2.2662E-02  9.7188E+00 -1.1178E+00  8.3471E-02  1.6518E+00 -1.3301E+00 -3.7582E-02
            -9.0924E+00

0ITERATION NO.:   79    OBJECTIVE VALUE:  -1814.74386335008        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     2655
 NPARAMETR:  1.1430E+00  1.4730E+00  3.3478E+04  9.3339E-01  2.4498E+00  2.1820E+00  7.4617E-01  5.3473E+01  4.7538E+00  1.7915E-01
             5.5658E+00
 PARAMETER:  2.3360E-01  4.8740E-01  1.0521E+01  3.1102E-02  9.9626E-01  8.8000E-01 -1.9261E-01  4.0782E+00  1.6585E+00 -1.6191E+00
             1.8171E+00
 GRADIENT:  -1.6611E+04  7.9527E+03  3.6860E+02  1.9053E-01  1.9517E+03 -2.2058E+03  4.0638E-01 -9.5055E+02 -1.1705E+03  2.3964E+03
             2.1140E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2655
 NO. OF SIG. DIGITS UNREPORTABLE

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.9334E-03 -5.1347E-02 -1.5516E-03  1.9191E-02 -1.1726E-03
 SE:             2.9486E-02  1.0839E-02  1.7030E-03  2.6582E-02  3.8430E-03
 N:                     100         100         100         100         100

 P VAL.:         7.8789E-01  2.1708E-06  3.6222E-01  4.7033E-01  7.6027E-01

 ETASHRINKSD(%)  1.2170E+00  6.3687E+01  9.4295E+01  1.0947E+01  8.7126E+01
 ETASHRINKVR(%)  2.4191E+00  8.6813E+01  9.9675E+01  2.0695E+01  9.8343E+01
 EBVSHRINKSD(%)  1.9850E+00  7.1135E+01  8.9205E+01  6.3217E+00  8.7382E+01
 EBVSHRINKVR(%)  3.9306E+00  9.1668E+01  9.8835E+01  1.2244E+01  9.8408E+01
 RELATIVEINF(%)  9.6000E+01  4.0156E+00  1.1612E+00  4.2469E+01  1.5570E+00
 EPSSHRINKSD(%)  8.8652E+00
 EPSSHRINKVR(%)  1.6944E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1814.7438633500824     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -160.65450358167163     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   105.68
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1814.744       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.14E+00  1.47E+00  3.36E+04  9.33E-01  2.45E+00  2.18E+00  7.46E-01  5.34E+01  4.75E+00  1.79E-01  5.57E+00
 


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
+        1.36E+07
 
 TH 2
+        2.94E+02  1.88E+06
 
 TH 3
+        6.88E-04  6.20E-04  7.78E-06
 
 TH 4
+       -3.89E+07  1.84E+02  2.47E-04  1.11E+08
 
 TH 5
+        9.95E+01  8.10E+01 -6.15E-03  1.89E+01  1.62E+05
 
 TH 6
+       -1.30E+02  3.61E+02  7.36E-04 -4.34E+01  1.07E+02  2.63E+05
 
 TH 7
+       -3.01E+02 -7.16E+01 -4.46E-05 -8.70E+02 -3.80E+00  3.64E+01  4.70E+07
 
 TH 8
+       -1.01E+00 -9.64E-01 -1.80E-05 -7.26E-01  9.99E+00 -1.19E+00 -5.41E-02  2.05E+01
 
 TH 9
+       -2.66E+01  9.48E+01  2.14E-04 -7.22E+00  3.23E+01 -3.27E+01  1.60E+01 -3.50E-01  1.56E+04
 
 TH10
+        7.36E+02 -9.99E+02 -2.12E-03  1.92E+02 -2.70E+02  3.78E+02 -1.65E+01  3.24E+00  9.21E+01  1.15E+07
 
 TH11
+        1.31E+01  1.11E+01 -3.30E-03  7.01E+00 -2.16E+02  2.64E+01  7.55E+00  5.35E+00  8.63E+00 -6.75E+01  9.40E+03
 
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
 #CPUT: Total CPU Time in Seconds,      121.724
Stop Time:
Sat Sep 18 06:45:03 CDT 2021
