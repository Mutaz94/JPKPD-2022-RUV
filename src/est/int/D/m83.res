Sat Sep 18 07:34:25 CDT 2021
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
$DATA ../../../../data/int/D/dat83.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m83.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   52685.3715682745        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.1676E+02  7.2834E+02  6.7597E+01  5.9965E+02 -1.7033E+02 -2.4956E+03 -1.6655E+03 -1.7488E+02 -2.2466E+03 -3.1864E+02
            -1.0510E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -615.456671473831        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3922E+00  1.6359E+00  8.3689E-01  1.4676E+00  1.1358E+00  3.0348E+00  3.6784E+00  1.0146E+00  1.6804E+00  1.1907E+00
             1.3895E+01
 PARAMETER:  4.3091E-01  5.9220E-01 -7.8062E-02  4.8361E-01  2.2730E-01  1.2101E+00  1.4025E+00  1.1453E-01  6.1903E-01  2.7452E-01
             2.7315E+00
 GRADIENT:  -2.4859E+00 -1.0546E+00 -4.0393E+01  1.0730E+02  2.2673E+01  7.4265E+01  2.3730E+00  3.0189E+00  1.0905E+01  2.4814E+01
             1.4752E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -675.633105361596        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.3115E+00  1.1824E+00  2.2917E+01  2.0936E+00  2.1316E+00  2.7576E+00  9.9305E+00  9.8483E-01  1.8887E+00  7.7381E-01
             1.3496E+01
 PARAMETER:  3.7119E-01  2.6757E-01  3.2319E+00  8.3890E-01  8.5685E-01  1.1144E+00  2.3956E+00  8.4713E-02  7.3588E-01 -1.5643E-01
             2.7024E+00
 GRADIENT:   6.0746E+00  2.1452E+01  5.7637E-02  2.0881E+01 -2.7580E+01  7.8501E+01  4.8856E+01  1.0320E-01  1.1582E+01  7.5234E+00
             1.2702E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -724.917551394505        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      224
 NPARAMETR:  1.1865E+00  1.0192E+00  2.5025E+01  1.4269E+00  2.4457E+00  1.8334E+00  6.1794E+00  2.8118E+00  9.0097E-01  2.5119E-01
             1.3168E+01
 PARAMETER:  2.7099E-01  1.1900E-01  3.3199E+00  4.5551E-01  9.9434E-01  7.0616E-01  1.9212E+00  1.1338E+00 -4.2826E-03 -1.2815E+00
             2.6778E+00
 GRADIENT:  -2.3220E+01  1.0789E+00 -1.2554E-01  1.5998E+01 -7.7195E+00 -2.8877E+00  1.2532E+01  5.5729E-02 -2.6407E+00  4.9366E-01
             1.0977E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -731.091525222748        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.2356E+00  1.1324E+00  1.8947E+01  1.2066E+00  2.4780E+00  1.8721E+00  5.3788E+00  3.7673E+00  8.9345E-01  1.1886E-01
             1.2162E+01
 PARAMETER:  3.1154E-01  2.2432E-01  3.0417E+00  2.8778E-01  1.0075E+00  7.2705E-01  1.7825E+00  1.4264E+00 -1.2666E-02 -2.0298E+00
             2.5984E+00
 GRADIENT:   1.6955E+01 -6.6803E+00 -6.2059E-01 -1.2873E+01  5.9792E+00 -4.8187E+00  8.5341E+00  6.8616E-02  4.7357E+00  1.4855E-01
             8.7806E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -733.311799764014        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      366
 NPARAMETR:  1.1798E+00  1.6124E+00  9.1955E+00  9.7449E-01  2.3516E+00  1.9108E+00  4.4565E+00  2.9371E+00  4.5830E-01  1.0671E-01
             1.2065E+01
 PARAMETER:  2.6535E-01  5.7775E-01  2.3187E+00  7.4161E-02  9.5508E-01  7.4754E-01  1.5944E+00  1.1774E+00 -6.8023E-01 -2.1376E+00
             2.5903E+00
 GRADIENT:  -6.7049E+00 -1.9252E+00  2.0456E-01 -3.5319E+00  1.6376E+00  1.4107E+00 -1.6456E+00  1.4827E-03  1.3431E+00  1.4214E-01
             5.3131E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -733.477733106342        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  1.1854E+00  1.6409E+00  8.6163E+00  9.6194E-01  2.3363E+00  1.9085E+00  4.4299E+00  2.8513E+00  4.0831E-01  1.0320E-01
             1.2069E+01
 PARAMETER:  2.7009E-01  5.9522E-01  2.2537E+00  6.1194E-02  9.4856E-01  7.4634E-01  1.5884E+00  1.1478E+00 -7.9573E-01 -2.1711E+00
             2.5906E+00
 GRADIENT:  -4.3656E+00 -1.4178E+00  2.5411E-01 -3.4149E+00  1.2150E+00  1.0810E+00 -1.1902E+00  1.1221E-02  1.0801E+00  1.3484E-01
             4.6557E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -733.547200957019        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      568
 NPARAMETR:  1.1852E+00  1.6413E+00  8.2800E+00  9.6199E-01  2.3374E+00  1.9078E+00  4.4335E+00  2.7748E+00  4.0815E-01  2.1293E-02
             1.2084E+01
 PARAMETER:  2.6995E-01  5.9552E-01  2.2138E+00  6.1245E-02  9.4903E-01  7.4596E-01  1.5892E+00  1.1206E+00 -7.9613E-01 -3.7494E+00
             2.5919E+00
 GRADIENT:  -4.7424E+00 -1.6019E+00 -2.7882E-02 -3.6541E+00  2.8822E+00  1.0285E+00 -7.5533E-01  6.7950E-03  1.1378E+00  5.7952E-03
             7.0075E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -733.550354330099        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      640
 NPARAMETR:  1.1853E+00  1.6414E+00  8.2170E+00  9.6199E-01  2.3373E+00  1.9078E+00  4.4338E+00  2.7858E+00  4.0813E-01  1.0000E-02
             1.2075E+01
 PARAMETER:  2.6996E-01  5.9554E-01  2.2062E+00  6.1246E-02  9.4900E-01  7.4594E-01  1.5893E+00  1.1245E+00 -7.9616E-01 -7.5373E+00
             2.5911E+00
 GRADIENT:  -4.5966E+00 -1.5490E+00 -1.2217E-01 -3.4667E+00  3.2272E+00  9.8609E-01 -7.0939E-01  6.1675E-02  1.1389E+00  0.0000E+00
             5.7152E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -733.550684401422        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      711
 NPARAMETR:  1.1853E+00  1.6414E+00  8.1975E+00  9.6199E-01  2.3372E+00  1.9077E+00  4.4343E+00  2.7873E+00  4.0812E-01  1.0000E-02
             1.2061E+01
 PARAMETER:  2.6996E-01  5.9557E-01  2.2038E+00  6.1249E-02  9.4894E-01  7.4590E-01  1.5894E+00  1.1251E+00 -7.9620E-01 -1.2960E+01
             2.5899E+00
 GRADIENT:  -4.3496E+00 -1.4289E+00 -1.4211E-01 -3.1686E+00  3.2898E+00  9.1673E-01 -7.2309E-01  7.1170E-02  1.1279E+00  0.0000E+00
             3.6090E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -733.550689909546        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      786
 NPARAMETR:  1.1853E+00  1.6414E+00  8.2032E+00  9.6199E-01  2.3372E+00  1.9077E+00  4.4344E+00  2.7869E+00  4.0811E-01  1.0000E-02
             1.2058E+01
 PARAMETER:  2.6996E-01  5.9558E-01  2.2045E+00  6.1250E-02  9.4894E-01  7.4590E-01  1.5894E+00  1.1249E+00 -7.9621E-01 -1.3815E+01
             2.5898E+00
 GRADIENT:  -4.3084E+00 -1.4060E+00 -1.3241E-01 -3.1225E+00  3.2495E+00  9.0557E-01 -7.3345E-01  6.6300E-02  1.1248E+00  0.0000E+00
             3.2691E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -733.552493422634        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      971
 NPARAMETR:  1.1853E+00  1.6415E+00  8.3062E+00  9.6199E-01  2.3371E+00  1.9077E+00  4.4353E+00  2.7464E+00  4.0811E-01  1.0000E-02
             1.2066E+01
 PARAMETER:  2.6997E-01  5.9559E-01  2.2170E+00  6.1249E-02  9.4893E-01  7.4591E-01  1.5896E+00  1.1103E+00 -7.9621E-01 -1.1154E+01
             2.5904E+00
 GRADIENT:  -5.5698E+00 -2.0400E+00  4.7593E-02 -3.4144E+00  2.3356E+00  5.6292E-02 -3.2706E+00 -7.4472E-02  1.1023E+00  0.0000E+00
            -3.6667E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -733.582100355213        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1150
 NPARAMETR:  1.1857E+00  1.6437E+00  8.2945E+00  9.6205E-01  2.3325E+00  1.9074E+00  4.5002E+00  2.7660E+00  4.0757E-01  1.0000E-02
             1.2064E+01
 PARAMETER:  2.7033E-01  5.9692E-01  2.2156E+00  6.1316E-02  9.4692E-01  7.4575E-01  1.6041E+00  1.1174E+00 -7.9753E-01 -2.7797E+01
             2.5902E+00
 GRADIENT:  -5.1670E+00 -1.4358E+00 -1.1713E-03 -5.7917E+00  2.0742E+00 -9.6634E-02  6.7307E-01  2.4575E-03  1.1957E+00  0.0000E+00
            -1.6440E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -733.936367892433        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  1.1964E+00  1.6353E+00  7.4461E+00  9.5490E-01  2.2946E+00  1.8950E+00  4.4831E+00  2.7187E+00  2.2423E-01  1.0000E-02
             1.2092E+01
 PARAMETER:  2.7930E-01  5.9185E-01  2.1077E+00  5.3852E-02  9.3056E-01  7.3924E-01  1.6003E+00  1.1001E+00 -1.3951E+00 -2.7797E+01
             2.5926E+00
 GRADIENT:   8.4047E-02 -1.5923E+00  1.6663E-01 -1.3343E+00  1.1961E+00 -7.7004E-01  1.3098E+00  1.6689E-01  2.0203E-01  0.0000E+00
             3.9393E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -734.002208605936        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1400
 NPARAMETR:  1.1960E+00  1.6368E+00  6.8933E+00  9.5478E-01  2.2963E+00  1.8937E+00  4.4897E+00  2.5333E+00  6.8261E-02  1.0000E-02
             1.2123E+01
 PARAMETER:  2.7900E-01  5.9271E-01  2.0305E+00  5.3725E-02  9.3129E-01  7.3852E-01  1.6018E+00  1.0295E+00 -2.5844E+00 -2.7797E+01
             2.5951E+00
 GRADIENT:  -7.9755E-01 -1.6996E+00 -2.4717E-01  2.5620E+00  3.9922E+00 -7.7226E-01  8.9991E-01  1.3340E-01  6.6728E-03  0.0000E+00
             6.0418E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -734.006183542605        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1475
 NPARAMETR:  1.1960E+00  1.6367E+00  6.9729E+00  9.5477E-01  2.2961E+00  1.8937E+00  4.4892E+00  2.5338E+00  2.4650E-02  1.0000E-02
             1.2115E+01
 PARAMETER:  2.7902E-01  5.9271E-01  2.0420E+00  5.3720E-02  9.3121E-01  7.3855E-01  1.6017E+00  1.0297E+00 -3.6030E+00 -2.7797E+01
             2.5944E+00
 GRADIENT:  -6.2770E-01 -1.5545E+00 -3.0909E-02  3.2931E+00  3.2914E+00 -7.8701E-01  5.9480E-01  4.1851E-02  4.7735E-04  0.0000E+00
             4.3485E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -734.006341739691        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:     1614
 NPARAMETR:  1.1960E+00  1.6367E+00  6.9821E+00  9.5477E-01  2.2961E+00  1.8937E+00  4.4892E+00  2.5192E+00  2.3977E-02  1.0000E-02
             1.2115E+01
 PARAMETER:  2.7902E-01  5.9271E-01  2.0434E+00  5.3720E-02  9.3121E-01  7.3856E-01  1.6017E+00  1.0239E+00 -3.6307E+00 -2.7797E+01
             2.5944E+00
 GRADIENT:  -1.8141E+00 -2.1705E+00  3.7143E-03  3.2349E+00  2.9119E+00 -1.6406E+00 -1.9674E+00 -3.7044E-04  2.0753E-04  0.0000E+00
            -4.1297E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -734.022559603014        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1797
 NPARAMETR:  1.1962E+00  1.6388E+00  6.8084E+00  9.5472E-01  2.2912E+00  1.8951E+00  4.5245E+00  2.4349E+00  1.0000E-02  1.0000E-02
             1.2118E+01
 PARAMETER:  2.7913E-01  5.9399E-01  2.0182E+00  5.3668E-02  9.2907E-01  7.3928E-01  1.6095E+00  9.8991E-01 -9.2703E+00 -2.7797E+01
             2.5947E+00
 GRADIENT:  -1.7528E+00 -1.8746E+00 -1.1951E-01  2.0334E+00  3.2728E+00 -1.4749E+00  2.4385E-01  5.9181E-03  0.0000E+00  0.0000E+00
             4.4184E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -734.075291363430        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:     1952
 NPARAMETR:  1.2003E+00  1.6740E+00  6.8285E+00  9.4969E-01  2.2651E+00  1.9066E+00  4.5049E+00  2.4797E+00  1.0000E-02  1.0000E-02
             1.2111E+01
 PARAMETER:  2.8258E-01  6.1522E-01  2.0211E+00  4.8385E-02  9.1761E-01  7.4533E-01  1.6052E+00  1.0081E+00 -1.2956E+01 -2.7797E+01
             2.5941E+00
 GRADIENT:   1.3295E+00  6.3383E-01  5.7993E-01  2.6899E+00 -8.1218E-01  1.1194E+00  2.8769E+00 -7.4375E-02  0.0000E+00  0.0000E+00
             2.6281E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -734.093207310208        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     2086             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2002E+00  1.6789E+00  6.5268E+00  9.4728E-01  2.2627E+00  1.9064E+00  4.4929E+00  2.5017E+00  1.0000E-02  1.0000E-02
             1.2117E+01
 PARAMETER:  2.8250E-01  6.1812E-01  1.9759E+00  4.5842E-02  9.1654E-01  7.4524E-01  1.6025E+00  1.0170E+00 -1.2956E+01 -2.7797E+01
             2.5946E+00
 GRADIENT:   1.1035E+00  3.9850E-01 -7.9522E-02  1.9976E+00  1.3125E+00  1.1067E+00  3.1704E+00  2.4697E-01  0.0000E+00  0.0000E+00
             4.2942E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -734.102965686684        NO. OF FUNC. EVALS.:  91
 CUMULATIVE NO. OF FUNC. EVALS.:     2177
 NPARAMETR:  1.2002E+00  1.6789E+00  6.2773E+00  9.4728E-01  2.2627E+00  1.9064E+00  4.4929E+00  2.2224E+00  1.0000E-02  1.0000E-02
             1.2117E+01
 PARAMETER:  2.8249E-01  6.1812E-01  1.9369E+00  4.5843E-02  9.1654E-01  7.4524E-01  1.6025E+00  8.9860E-01 -1.2956E+01 -2.7797E+01
             2.5946E+00
 GRADIENT:  -1.4438E-01 -3.0519E-01 -3.0959E-02  2.7736E+00  1.8408E+00  2.0909E-01  4.5830E-01  3.9412E-03  0.0000E+00  0.0000E+00
            -8.9553E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -734.121313817595        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     2312
 NPARAMETR:  1.1993E+00  1.6829E+00  6.2325E+00  9.4163E-01  2.2493E+00  1.9019E+00  4.4552E+00  2.2426E+00  1.0000E-02  1.0000E-02
             1.2133E+01
 PARAMETER:  2.8171E-01  6.2049E-01  1.9298E+00  3.9855E-02  9.1060E-01  7.4286E-01  1.5941E+00  9.0763E-01 -1.2956E+01 -2.7797E+01
             2.5959E+00
 GRADIENT:   4.4864E-01 -1.5311E-01  2.3803E-01  1.4000E+00  5.5390E-01  4.7395E-01  1.9183E+00  1.4101E-02  0.0000E+00  0.0000E+00
             7.2697E+00

0ITERATION NO.:  110    OBJECTIVE VALUE:  -734.126531881199        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     2439
 NPARAMETR:  1.1992E+00  1.6830E+00  5.8310E+00  9.4164E-01  2.2494E+00  1.9018E+00  4.4558E+00  1.8917E+00  1.0000E-02  1.0000E-02
             1.2135E+01
 PARAMETER:  2.8168E-01  6.2055E-01  1.8632E+00  3.9864E-02  9.1068E-01  7.4279E-01  1.5942E+00  7.3748E-01 -1.2956E+01 -2.7797E+01
             2.5961E+00
 GRADIENT:   3.0542E-01 -3.0916E-01 -1.0330E-01  2.4440E+00  2.7512E+00  4.2100E-01  1.9201E+00  3.5856E-02  0.0000E+00  0.0000E+00
             7.1466E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -734.186206508903        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:     2520
 NPARAMETR:  1.1992E+00  1.6830E+00  5.2477E+00  9.4164E-01  2.2495E+00  1.9017E+00  4.4555E+00  4.8872E-01  1.0000E-02  1.0000E-02
             1.2126E+01
 PARAMETER:  2.8167E-01  6.2059E-01  1.7578E+00  3.9870E-02  9.1070E-01  7.4275E-01  1.5941E+00 -6.1596E-01 -1.2956E+01 -2.7797E+01
             2.5954E+00
 GRADIENT:   3.3915E-01 -3.6255E-01 -2.9705E-01  5.4851E+00  4.9857E+00  2.8419E-01  1.3649E+00  2.8540E-02  0.0000E+00  0.0000E+00
             4.2570E+00

0ITERATION NO.:  119    OBJECTIVE VALUE:  -734.186292207206        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:     2635
 NPARAMETR:  1.1991E+00  1.6833E+00  5.2456E+00  9.4167E-01  2.2500E+00  1.9013E+00  4.4573E+00  4.8718E-01  1.0000E-02  1.0000E-02
             1.2134E+01
 PARAMETER:  2.8167E-01  6.2059E-01  1.7578E+00  3.9870E-02  9.1070E-01  7.4275E-01  1.5941E+00 -6.1851E-01 -1.2956E+01 -2.7797E+01
             2.5954E+00
 GRADIENT:   8.1132E+03 -1.8429E+03  1.3005E+03 -2.2861E+04 -2.5051E+03  3.0719E+03 -7.1664E+02  2.8688E-02  0.0000E+00  0.0000E+00
            -8.7587E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2635
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6520E-02  9.6889E-03 -4.1434E-04 -7.7733E-04 -1.0963E-05
 SE:             2.7952E-02  2.7168E-02  1.0897E-03  1.6529E-04  1.0347E-04
 N:                     100         100         100         100         100

 P VAL.:         5.5453E-01  7.2137E-01  7.0378E-01  2.5706E-06  9.1562E-01

 ETASHRINKSD(%)  6.3565E+00  8.9825E+00  9.6349E+01  9.9446E+01  9.9653E+01
 ETASHRINKVR(%)  1.2309E+01  1.7158E+01  9.9867E+01  9.9997E+01  9.9999E+01
 EBVSHRINKSD(%)  6.8501E+00  5.9967E+00  9.6335E+01  9.9625E+01  9.9601E+01
 EBVSHRINKVR(%)  1.3231E+01  1.1634E+01  9.9866E+01  9.9999E+01  9.9998E+01
 RELATIVEINF(%)  8.6377E+01  4.9062E+01  3.0870E-02  7.9095E-04  3.6601E-04
 EPSSHRINKSD(%)  3.5474E+00
 EPSSHRINKVR(%)  6.9689E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -734.18629220720618     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       919.90306756120458     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    73.87
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    16.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -734.186       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.20E+00  1.68E+00  5.25E+00  9.42E-01  2.25E+00  1.90E+00  4.46E+00  4.87E-01  1.00E-02  1.00E-02  1.21E+01
 


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
+        1.11E+04
 
 TH 2
+       -3.74E+02  1.26E+01
 
 TH 3
+        1.79E-01 -6.04E-03  2.89E-06
 
 TH 4
+       -5.11E+05  1.73E+04 -8.26E+00  2.36E+07
 
 TH 5
+       -4.82E+01  1.63E+00 -7.80E-04  2.23E+03  2.10E-01
 
 TH 6
+        1.48E+02 -5.01E+00  2.39E-03 -6.84E+03 -6.46E-01  1.98E+00
 
 TH 7
+       -3.92E-01  1.33E-02 -6.35E-06  1.81E+01  1.71E-03 -5.26E-03  1.39E-05
 
 TH 8
+        1.23E-02 -4.14E-04  1.98E-07 -5.66E-01 -5.35E-05  1.64E-04 -4.35E-07  1.36E-08
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        2.26E-01 -7.64E-03  3.65E-06 -1.04E+01 -9.85E-04  3.02E-03 -8.02E-06  2.50E-07  0.00E+00  0.00E+00  4.61E-06
 
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
+        5.01E+06
 
 TH 2
+        6.83E+02  5.24E+05
 
 TH 3
+       -1.17E+01  4.08E+00  6.72E+03
 
 TH 4
+       -2.78E+02  1.19E+02  4.22E+01  6.45E+07
 
 TH 5
+        6.50E+02 -2.16E+02 -7.10E-01  3.98E+01  1.36E+05
 
 TH 6
+       -1.23E+03  1.64E+02 -2.82E+00 -5.94E+01  1.55E+02  2.85E+05
 
 TH 7
+        2.46E+02 -3.12E+01  3.83E-01 -1.12E+01 -2.97E+01  2.47E+02  1.13E+04
 
 TH 8
+        3.45E-01 -3.54E-01  6.54E-03 -1.42E+00  1.30E-01  6.63E-02 -1.92E-02 -4.32E-01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        4.80E+01 -9.90E+00  5.75E-02 -1.21E+01 -6.17E+00  5.27E+01 -1.04E+01  1.51E-02  0.00E+00  0.00E+00  5.76E+02
 
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
+        1.14E+07
 
 TH 2
+       -3.71E+06  1.20E+06
 
 TH 3
+        4.19E+05 -1.36E+05  1.54E+04
 
 TH 4
+       -4.11E+07  1.33E+07 -1.50E+06  1.47E+08
 
 TH 5
+       -1.89E+06  6.11E+05 -6.92E+04  6.78E+06  3.11E+05
 
 TH 6
+        2.73E+06 -8.85E+05  1.00E+05 -9.81E+06 -4.51E+05  6.53E+05
 
 TH 7
+       -5.46E+05  1.77E+05 -2.00E+04  1.96E+06  9.00E+04 -1.30E+05  2.60E+04
 
 TH 8
+       -1.21E-01 -9.98E-02 -2.50E-02 -1.12E+00  2.37E-01  1.99E-01  7.41E-02  6.31E-03
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.30E+05  4.21E+04 -4.76E+03  4.67E+05  2.15E+04 -3.10E+04  6.21E+03 -5.88E-02  0.00E+00  0.00E+00  1.60E+03
 
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
 #CPUT: Total CPU Time in Seconds,       90.102
Stop Time:
Sat Sep 18 07:35:56 CDT 2021
