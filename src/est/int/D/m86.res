Sat Sep 18 07:37:48 CDT 2021
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
$DATA ../../../../data/int/D/dat86.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   42462.7171029778        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8770E+02  7.5123E+02 -1.2109E+01  6.0611E+02  1.1277E+02 -3.3877E+03 -1.7459E+03 -3.2793E+01 -2.5137E+03 -8.7704E+02
            -8.3391E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -582.242059766629        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1244E+00  1.7647E+00  9.2426E-01  2.3742E+00  9.7870E-01  4.9912E+00  2.9672E+00  9.5553E-01  2.2672E+00  1.1424E+00
             1.2774E+01
 PARAMETER:  2.1725E-01  6.6799E-01  2.1235E-02  9.6465E-01  7.8468E-02  1.7077E+00  1.1876E+00  5.4513E-02  9.1855E-01  2.3314E-01
             2.6474E+00
 GRADIENT:  -1.4797E+01  4.8964E+01 -4.0264E+01  2.0633E+02 -6.2731E+00  1.2918E+02 -9.9089E+01  4.9696E+00 -6.9587E+01  2.5517E+01
             1.3156E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -650.918895129145        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.1381E+00  1.8508E+00  2.6822E+01  5.1962E+00  3.0763E+00  3.8082E+00  4.1615E+00  5.4983E-01  6.5606E+00  9.8358E-01
             1.2417E+01
 PARAMETER:  2.2937E-01  7.1562E-01  3.3892E+00  1.7479E+00  1.2237E+00  1.4372E+00  1.5259E+00 -4.9814E-01  1.9811E+00  8.3442E-02
             2.6190E+00
 GRADIENT:  -1.2135E+01  1.1161E+01 -1.6900E+01  7.1796E+01  5.4365E+01  1.0489E+02  2.7792E+01  1.6770E-01  2.0340E+01  1.0755E+01
             1.9476E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -727.677041299854        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0341E+00  8.4846E-01  3.1840E+01  1.9090E+00  1.9962E+00  4.5101E+00  2.2677E+00  4.3056E+00  3.5701E+00  1.4393E+00
             1.1682E+01
 PARAMETER:  1.3349E-01 -6.4338E-02  3.5607E+00  7.4656E-01  7.9126E-01  1.6063E+00  9.1878E-01  1.5599E+00  1.3726E+00  4.6418E-01
             2.5581E+00
 GRADIENT:  -2.3166E+01 -3.3253E+00 -4.2013E+00  2.1678E+01 -2.7422E+01  1.1138E+02  1.1808E+01 -7.3287E-01  9.7107E+00  3.4581E+01
             1.9048E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -820.616188283853        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.1418E+00  1.6792E+00  2.7652E+01  8.7431E-01  2.3421E+00  2.3549E+00  9.2389E-01  4.9135E+00  4.8322E+00  3.5416E-01
             1.0693E+01
 PARAMETER:  2.3262E-01  6.1832E-01  3.4197E+00 -3.4322E-02  9.5106E-01  9.5651E-01  2.0833E-02  1.6920E+00  1.6753E+00 -9.3802E-01
             2.4696E+00
 GRADIENT:  -7.2325E+00  2.7466E+01 -8.5383E+00  6.2851E+00  1.2452E+01 -1.5729E+01  8.0850E-02  6.3157E-01 -2.5169E+01  1.9989E+00
             8.7537E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -829.373794692640        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  1.1477E+00  1.6042E+00  8.7266E+01  7.4255E-01  2.3166E+00  2.4561E+00  7.1028E-01  9.0963E+00  5.4164E+00  1.5668E-01
             1.0382E+01
 PARAMETER:  2.3772E-01  5.7261E-01  4.5690E+00 -1.9766E-01  9.4009E-01  9.9857E-01 -2.4209E-01  2.3079E+00  1.7894E+00 -1.7536E+00
             2.4401E+00
 GRADIENT:   1.2236E-01  6.0268E+00 -1.3629E+00 -1.2285E+00  8.9589E+00 -3.2447E+00 -4.3958E-01  2.0948E+00 -7.8007E-01  4.7095E-01
            -2.0575E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -831.376935334310        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  1.1496E+00  1.5720E+00  1.5318E+02  7.7247E-01  2.2757E+00  2.4711E+00  7.0410E-01  9.1493E+00  5.3545E+00  9.9766E-02
             1.0527E+01
 PARAMETER:  2.3942E-01  5.5235E-01  5.1316E+00 -1.5816E-01  9.2230E-01  1.0047E+00 -2.5083E-01  2.3137E+00  1.7779E+00 -2.2049E+00
             2.4539E+00
 GRADIENT:  -5.7426E-01 -8.0482E-01 -1.7569E-01 -1.1952E+00 -1.2958E+00 -3.0204E-01  9.7915E-01  2.2137E-01  1.5095E+00  1.8649E-01
             7.5950E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -831.503493459880        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      568
 NPARAMETR:  1.1530E+00  1.5714E+00  1.5422E+02  7.8692E-01  2.2888E+00  2.4830E+00  6.7871E-01  9.1724E+00  5.3592E+00  5.8577E-02
             1.0512E+01
 PARAMETER:  2.4234E-01  5.5194E-01  5.1384E+00 -1.3962E-01  9.2804E-01  1.0095E+00 -2.8756E-01  2.3162E+00  1.7788E+00 -2.7374E+00
             2.4525E+00
 GRADIENT:   5.0295E-01  3.1846E+00 -1.8384E-01  6.6104E-02 -5.5969E-03  1.1627E+00  2.3061E-01  2.0395E-01  1.5978E+00  6.3943E-02
             2.1314E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -831.512149592550        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      640
 NPARAMETR:  1.1507E+00  1.5713E+00  1.5421E+02  7.8496E-01  2.2903E+00  2.4722E+00  6.6829E-01  9.1728E+00  5.3592E+00  3.1439E-02
             1.0510E+01
 PARAMETER:  2.4036E-01  5.5192E-01  5.1383E+00 -1.4212E-01  9.2868E-01  1.0051E+00 -3.0303E-01  2.3162E+00  1.7788E+00 -3.3597E+00
             2.4523E+00
 GRADIENT:  -1.6679E-01  4.1299E+00 -1.8483E-01 -2.1956E-02  1.6019E-01 -3.0136E-01 -5.0592E-02  1.9514E-01  1.2134E+00  1.8419E-02
             1.2821E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -831.522866681782        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  1.1510E+00  1.5713E+00  1.5419E+02  7.8510E-01  2.2888E+00  2.4750E+00  6.7080E-01  9.1734E+00  5.3588E+00  1.0000E-02
             1.0506E+01
 PARAMETER:  2.4067E-01  5.5188E-01  5.1382E+00 -1.4194E-01  9.2803E-01  1.0062E+00 -2.9928E-01  2.3163E+00  1.7787E+00 -4.8518E+00
             2.4520E+00
 GRADIENT:  -1.4451E+00  2.8606E+00 -2.3286E-01 -1.6430E-01 -4.1573E-01 -3.0782E+00 -3.1312E-02  1.3281E-02 -6.7534E+00  0.0000E+00
            -3.7005E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -831.544622146027        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      913
 NPARAMETR:  1.1559E+00  1.5713E+00  1.5422E+02  7.8849E-01  2.2914E+00  2.4972E+00  6.7266E-01  9.1729E+00  5.3604E+00  1.0504E-02
             1.0511E+01
 PARAMETER:  2.4491E-01  5.5190E-01  5.1384E+00 -1.3763E-01  9.2918E-01  1.0152E+00 -2.9651E-01  2.3163E+00  1.7790E+00 -4.4560E+00
             2.4525E+00
 GRADIENT:  -9.7524E-02  2.7239E+00 -2.3698E-01  2.8565E-02 -9.8340E-02 -1.3472E-01  4.6881E-02  1.4797E-02 -6.4355E+00  2.0403E-03
            -2.9071E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -831.570534190352        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1091
 NPARAMETR:  1.1584E+00  1.5711E+00  1.5433E+02  7.8960E-01  2.2919E+00  2.4946E+00  6.5603E-01  9.1734E+00  5.3916E+00  1.0000E-02
             1.0566E+01
 PARAMETER:  2.4707E-01  5.5177E-01  5.1391E+00 -1.3623E-01  9.2938E-01  1.0141E+00 -3.2155E-01  2.3163E+00  1.7849E+00 -5.2120E+00
             2.4576E+00
 GRADIENT:  -5.5707E-02  1.9122E+00 -2.3157E-01  4.5474E-02  2.7447E-02 -2.0632E-01 -3.1250E-02  2.5031E-02 -4.9635E+00  0.0000E+00
             6.0131E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -831.648373606288        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1269
 NPARAMETR:  1.1558E+00  1.5686E+00  1.5451E+02  7.7543E-01  2.2882E+00  2.4990E+00  6.4166E-01  9.1901E+00  5.5271E+00  1.0000E-02
             1.0516E+01
 PARAMETER:  2.4482E-01  5.5016E-01  5.1403E+00 -1.5434E-01  9.2775E-01  1.0159E+00 -3.4370E-01  2.3181E+00  1.8097E+00 -4.3510E+01
             2.4529E+00
 GRADIENT:  -1.3412E-03 -1.3123E+00 -2.3720E-01  4.2484E-03  9.8879E-04  1.5892E-02 -4.4678E-03  1.0309E-01 -6.0212E-02  0.0000E+00
            -4.1645E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -831.652167511401        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1433
 NPARAMETR:  1.1559E+00  1.5692E+00  1.5686E+02  7.7546E-01  2.2881E+00  2.4988E+00  6.4182E-01  9.1817E+00  5.5265E+00  1.0000E-02
             1.0517E+01
 PARAMETER:  2.4485E-01  5.5060E-01  5.1554E+00 -1.5430E-01  9.2773E-01  1.0158E+00 -3.4345E-01  2.3172E+00  1.8096E+00 -4.3296E+01
             2.4530E+00
 GRADIENT:  -5.6108E-03 -1.0843E+00 -1.9159E-01  3.7063E-02 -7.5140E-02 -3.2208E-03 -3.2633E-02  5.9189E-03 -1.4907E-01  0.0000E+00
            -6.1397E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -831.657075356476        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1617
 NPARAMETR:  1.1561E+00  1.5764E+00  1.5646E+02  7.7151E-01  2.2898E+00  2.4982E+00  6.5029E-01  9.1922E+00  5.5423E+00  1.0000E-02
             1.0521E+01
 PARAMETER:  2.4503E-01  5.5513E-01  5.1528E+00 -1.5941E-01  9.2849E-01  1.0156E+00 -3.3033E-01  2.3184E+00  1.8124E+00 -4.3296E+01
             2.4534E+00
 GRADIENT:  -2.4234E-02 -7.8724E-01 -2.0373E-01  1.1710E-01  1.5369E-01 -3.5509E-02 -3.6991E-02  2.3538E-02 -8.6187E-02  0.0000E+00
             3.9211E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -831.657246526668        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     1762
 NPARAMETR:  1.1561E+00  1.5763E+00  1.5640E+02  7.6952E-01  2.2888E+00  2.4984E+00  6.5158E-01  9.1938E+00  5.5430E+00  1.0000E-02
             1.0519E+01
 PARAMETER:  2.4505E-01  5.5509E-01  5.1524E+00 -1.6199E-01  9.2801E-01  1.0157E+00 -3.2836E-01  2.3185E+00  1.8125E+00 -4.3296E+01
             2.4532E+00
 GRADIENT:   2.0494E-02 -1.0931E+00 -2.0476E-01 -2.2713E-02  3.9310E-02 -1.3581E-02 -1.1555E-03  3.2011E-02 -1.1824E-01  0.0000E+00
             2.9096E-01

0ITERATION NO.:   77    OBJECTIVE VALUE:  -831.657251779299        NO. OF FUNC. EVALS.:  87
 CUMULATIVE NO. OF FUNC. EVALS.:     1849
 NPARAMETR:  1.1560E+00  1.5759E+00  1.5600E+02  7.6976E-01  2.2886E+00  2.4987E+00  6.5149E-01  9.2044E+00  5.5462E+00  1.0000E-02
             1.0506E+01
 PARAMETER:  2.4505E-01  5.5509E-01  5.1524E+00 -1.6158E-01  9.2792E-01  1.0157E+00 -3.2817E-01  2.3185E+00  1.8125E+00 -4.3296E+01
             2.4532E+00
 GRADIENT:   1.8032E-02  1.6057E+03  1.7215E+02  3.2868E-03 -4.6896E-03 -3.0701E-02  6.1641E-03 -3.8248E+02 -1.3845E-01  0.0000E+00
             3.5964E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1849
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1663E-02 -5.3182E-02 -6.2783E-03  2.4239E-02  1.6106E-06
 SE:             2.8488E-02  1.0268E-02  6.5034E-03  2.5530E-02  1.0444E-04
 N:                     100         100         100         100         100

 P VAL.:         6.8225E-01  2.2281E-07  3.3435E-01  3.4240E-01  9.8770E-01

 ETASHRINKSD(%)  4.5612E+00  6.5602E+01  7.8213E+01  1.4471E+01  9.9650E+01
 ETASHRINKVR(%)  8.9143E+00  8.8168E+01  9.5253E+01  2.6848E+01  9.9999E+01
 EBVSHRINKSD(%)  4.0206E+00  7.2171E+01  7.7844E+01  9.5594E+00  9.9547E+01
 EBVSHRINKVR(%)  7.8796E+00  9.2255E+01  9.5091E+01  1.8205E+01  9.9998E+01
 RELATIVEINF(%)  9.1906E+01  4.0631E+00  4.6077E+00  4.4000E+01  1.8198E-03
 EPSSHRINKSD(%)  4.7110E+00
 EPSSHRINKVR(%)  9.2001E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -831.65725177929869     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       822.43210798911207     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    54.88
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.64
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -831.657       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.16E+00  1.58E+00  1.56E+02  7.70E-01  2.29E+00  2.50E+00  6.52E-01  9.19E+00  5.54E+00  1.00E-02  1.05E+01
 


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
+        1.18E+02
 
 TH 2
+        6.19E+02  2.91E+05
 
 TH 3
+        2.83E-01 -7.25E-02  3.41E-01
 
 TH 4
+       -5.78E+00  3.69E+02  2.13E-01  1.97E+01
 
 TH 5
+       -7.51E-01 -1.02E+03 -1.11E+00 -1.01E+01  4.12E+01
 
 TH 6
+        6.95E-01 -1.00E+05  1.79E-01 -3.16E+00 -1.09E+00  2.58E+01
 
 TH 7
+        6.73E+00 -7.65E+01 -4.19E-02  3.27E+00  4.41E+00 -1.91E+00  4.71E+00
 
 TH 8
+       -1.05E+01  3.07E+00  1.35E-01 -7.76E+00  4.15E+01 -6.78E+00  1.44E+00  4.86E+02
 
 TH 9
+        3.35E-01 -2.53E+04  5.61E-02  3.18E+00  1.02E+00 -1.48E-01  2.45E+00 -2.12E+00  3.93E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        4.36E+00  9.75E+03 -2.42E-01  3.31E+00 -3.30E+01  7.07E+00  2.50E+00  9.19E+00  2.65E+00  0.00E+00  3.36E+02
 
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
 #CPUT: Total CPU Time in Seconds,       71.642
Stop Time:
Sat Sep 18 07:39:02 CDT 2021
