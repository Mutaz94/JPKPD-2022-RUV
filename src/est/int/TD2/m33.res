Wed Sep 29 07:16:52 CDT 2021
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
$DATA ../../../../data/int/TD2/dat33.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m33.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3066.13392763739        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3087E+02  1.7731E+01  1.0280E+02  1.5674E+02  9.6998E+01  4.2250E+01 -1.5594E+01 -1.3048E+02 -8.0539E+01  6.5483E+00
            -1.2797E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3321.39475562590        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7019E-01  1.1136E+00  7.8439E-01  9.2501E-01  9.6129E-01  9.2772E-01  1.0006E+00  1.2203E+00  1.1705E+00  9.4983E-01
             1.6011E+00
 PARAMETER:  6.9741E-02  2.0759E-01 -1.4285E-01  2.2045E-02  6.0520E-02  2.4973E-02  1.0061E-01  2.9909E-01  2.5740E-01  4.8525E-02
             5.7067E-01
 GRADIENT:   2.1281E+02  8.3662E+01 -4.2797E+01  2.1362E+01 -1.3205E+01  1.7222E+00  8.2460E+00 -4.8231E+00  9.5838E+00 -3.8449E+00
             1.2615E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3331.23073569578        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.5507E-01  1.2238E+00  8.8849E-01  8.7458E-01  1.0841E+00  9.7945E-01  9.7102E-01  2.1928E+00  1.0861E+00  1.1062E+00
             1.5402E+00
 PARAMETER:  5.4027E-02  3.0193E-01 -1.8226E-02 -3.4009E-02  1.8075E-01  7.9233E-02  7.0592E-02  8.8517E-01  1.8259E-01  2.0089E-01
             5.3189E-01
 GRADIENT:   1.8586E+02  1.2290E+02 -2.8251E+01  2.9103E+01 -1.4963E+01  2.5510E+01  1.6725E+01  2.7387E+01 -5.2680E+00 -4.7273E+00
             9.8145E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3353.62327441378        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      337
 NPARAMETR:  9.4436E-01  1.4220E+00  1.0978E+00  7.6516E-01  1.3226E+00  9.5793E-01  7.1740E-01  2.2999E+00  1.2171E+00  1.2870E+00
             1.4714E+00
 PARAMETER:  4.2755E-02  4.5203E-01  1.9327E-01 -1.6767E-01  3.7959E-01  5.7018E-02 -2.3213E-01  9.3286E-01  2.9643E-01  3.5228E-01
             4.8623E-01
 GRADIENT:  -1.1527E+01 -8.6309E+00 -8.8789E+00  1.1406E+01  9.5591E-03  1.2202E+00 -9.5069E+00 -1.0547E+01 -7.4418E+00 -8.6951E-01
             2.2512E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3356.11888710096        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      501
 NPARAMETR:  9.3637E-01  1.4890E+00  1.1841E+00  7.2788E-01  1.4099E+00  9.4560E-01  8.6828E-01  2.6652E+00  1.2097E+00  1.3442E+00
             1.4804E+00
 PARAMETER:  3.4254E-02  4.9810E-01  2.6901E-01 -2.1762E-01  4.4351E-01  4.4068E-02 -4.1237E-02  1.0803E+00  2.9035E-01  3.9581E-01
             4.9232E-01
 GRADIENT:  -3.2797E+01 -1.1029E+01 -3.7043E+00  1.2924E+01  5.8009E+00 -4.3717E+00  6.3876E+00 -3.8311E+00  8.4619E+00  3.0264E+00
             2.1681E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3356.76514499023        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      661
 NPARAMETR:  9.4862E-01  1.4928E+00  1.1858E+00  7.2707E-01  1.4067E+00  9.5108E-01  8.5356E-01  2.6507E+00  1.1582E+00  1.3335E+00
             1.4669E+00
 PARAMETER:  4.7258E-02  5.0063E-01  2.7038E-01 -2.1873E-01  4.4128E-01  4.9845E-02 -5.8345E-02  1.0748E+00  2.4690E-01  3.8784E-01
             4.8318E-01
 GRADIENT:  -6.7488E-01 -4.6700E+00 -2.9019E+00  1.1709E+01  2.7351E+00 -1.6148E+00  9.6545E-01 -5.9256E+00  1.0863E+00  7.4378E-02
             2.0737E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3357.14084213603        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      837
 NPARAMETR:  9.4351E-01  1.4979E+00  1.2234E+00  7.1510E-01  1.4217E+00  9.5175E-01  8.5835E-01  2.8182E+00  1.1394E+00  1.3453E+00
             1.4628E+00
 PARAMETER:  4.1847E-02  5.0405E-01  3.0164E-01 -2.3533E-01  4.5183E-01  5.0548E-02 -5.2744E-02  1.1361E+00  2.3050E-01  3.9662E-01
             4.8036E-01
 GRADIENT:  -1.3584E+01 -1.8338E+01 -6.8367E-01 -3.7531E-01  1.5017E-01 -1.4117E+00 -4.2154E-01 -1.4326E+00  9.6819E-01 -7.5432E-01
            -2.8214E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3357.23725187006        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  9.4406E-01  1.5009E+00  1.2514E+00  7.1306E-01  1.4241E+00  9.5237E-01  8.8489E-01  2.8978E+00  1.0919E+00  1.3628E+00
             1.4629E+00
 PARAMETER:  4.2438E-02  5.0606E-01  3.2428E-01 -2.3819E-01  4.5357E-01  5.1197E-02 -2.2296E-02  1.1640E+00  1.8793E-01  4.0952E-01
             4.8044E-01
 GRADIENT:  -1.2123E+01 -1.8686E+01  1.1480E+00 -2.3602E+00 -5.6541E+00 -1.1122E+00 -1.9553E-01 -3.8216E-01 -1.8688E-01  2.7032E-01
            -1.5168E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3359.47600531500        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  9.4965E-01  1.5336E+00  1.3042E+00  6.9421E-01  1.4591E+00  9.5476E-01  9.2227E-01  3.0518E+00  1.0261E+00  1.3971E+00
             1.4654E+00
 PARAMETER:  4.8336E-02  5.2760E-01  3.6560E-01 -2.6498E-01  4.7780E-01  5.3700E-02  1.9087E-02  1.2157E+00  1.2575E-01  4.3438E-01
             4.8213E-01
 GRADIENT:   2.2331E+00 -1.7640E+01 -4.1237E+00 -1.7336E-01 -1.3483E+01 -3.7497E-03 -9.6469E-01 -1.0542E+00 -2.6340E-01  4.8508E+00
            -1.5888E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3359.50377852707        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1380
 NPARAMETR:  9.4965E-01  1.5361E+00  1.3045E+00  6.9204E-01  1.4599E+00  9.5704E-01  9.3513E-01  3.0522E+00  1.0044E+00  1.3968E+00
             1.4652E+00
 PARAMETER:  4.8334E-02  5.2922E-01  3.6581E-01 -2.6811E-01  4.7837E-01  5.6094E-02  3.2934E-02  1.2159E+00  1.0442E-01  4.3421E-01
             4.8199E-01
 GRADIENT:   2.2373E+00 -1.7324E+01 -3.7827E+00 -1.3167E+00 -1.4899E+01  9.2592E-01  1.1759E-02 -1.7716E+00 -9.4537E-01  3.6904E+00
            -2.2413E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3359.69230866079        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1540
 NPARAMETR:  9.5014E-01  1.5537E+00  1.3070E+00  6.9108E-01  1.4636E+00  9.5417E-01  9.2698E-01  3.0700E+00  1.0265E+00  1.3717E+00
             1.4688E+00
 PARAMETER:  4.8851E-02  5.4065E-01  3.6771E-01 -2.6950E-01  4.8087E-01  5.3084E-02  2.4177E-02  1.2217E+00  1.2611E-01  4.1604E-01
             4.8448E-01
 GRADIENT:   3.2791E+00  6.0212E-01 -3.1564E+00  7.6896E+00 -1.5552E+01 -2.9487E-01  6.9639E-01 -7.5093E-01 -1.3429E-01 -7.4039E-01
             1.7845E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3359.85681867252        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1721            RESET HESSIAN, TYPE II
 NPARAMETR:  9.4907E-01  1.5573E+00  1.3094E+00  6.8914E-01  1.4878E+00  9.5499E-01  9.2215E-01  3.0532E+00  1.0304E+00  1.3728E+00
             1.4685E+00
 PARAMETER:  4.7732E-02  5.4297E-01  3.6954E-01 -2.7232E-01  4.9731E-01  5.3947E-02  1.8948E-02  1.2162E+00  1.2996E-01  4.1688E-01
             4.8424E-01
 GRADIENT:   1.9155E+02  3.4462E+02  1.4145E-01  6.6863E+01  8.6108E+01  1.7803E+01  5.0457E+00  7.2297E+00  2.7878E+00  1.3125E+01
             9.1416E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3359.90171668174        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1899
 NPARAMETR:  9.4884E-01  1.5599E+00  1.3147E+00  6.8801E-01  1.4864E+00  9.5452E-01  9.2294E-01  3.0625E+00  1.0255E+00  1.3829E+00
             1.4685E+00
 PARAMETER:  4.7486E-02  5.4464E-01  3.7359E-01 -2.7395E-01  4.9632E-01  5.3452E-02  1.9808E-02  1.2192E+00  1.2521E-01  4.2418E-01
             4.8424E-01
 GRADIENT:  -1.2291E-02 -1.0548E+00 -3.1426E+00  1.1425E+01 -4.5767E+00 -1.8465E-01 -3.2967E-01 -2.4339E+00 -1.1025E-01  1.4121E-01
             6.5121E-01

0ITERATION NO.:   62    OBJECTIVE VALUE:  -3359.90219118792        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1956
 NPARAMETR:  9.4884E-01  1.5606E+00  1.3147E+00  6.8801E-01  1.4864E+00  9.5464E-01  9.2486E-01  3.0625E+00  1.0245E+00  1.3822E+00
             1.4685E+00
 PARAMETER:  4.7490E-02  5.4509E-01  3.7359E-01 -2.7395E-01  4.9634E-01  5.3574E-02  2.1892E-02  1.2192E+00  1.2419E-01  4.2371E-01
             4.8424E-01
 GRADIENT:  -4.7530E-02 -1.9766E-01 -3.1201E+00  1.1814E+01 -4.6920E+00 -1.3832E-01  3.9447E-03 -2.4783E+00 -5.4997E-02 -5.9378E-02
             6.1419E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1956
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1827E-03 -2.9804E-02 -2.9991E-02  2.7029E-02 -3.4615E-02
 SE:             2.9806E-02  2.3872E-02  2.0097E-02  2.1956E-02  2.4746E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6835E-01  2.1185E-01  1.3561E-01  2.1830E-01  1.6186E-01

 ETASHRINKSD(%)  1.4693E-01  2.0027E+01  3.2674E+01  2.6444E+01  1.7099E+01
 ETASHRINKVR(%)  2.9365E-01  3.6043E+01  5.4672E+01  4.5895E+01  3.1274E+01
 EBVSHRINKSD(%)  6.0271E-01  1.9318E+01  3.5390E+01  3.2440E+01  1.3227E+01
 EBVSHRINKVR(%)  1.2018E+00  3.4904E+01  5.8256E+01  5.4357E+01  2.4705E+01
 RELATIVEINF(%)  9.8793E+01  1.4877E+01  3.1045E+01  9.9747E+00  3.7196E+01
 EPSSHRINKSD(%)  2.0262E+01
 EPSSHRINKVR(%)  3.6419E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3359.9021911879208     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1705.8128314195101     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    61.17
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.68
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3359.902       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.49E-01  1.56E+00  1.31E+00  6.88E-01  1.49E+00  9.55E-01  9.25E-01  3.06E+00  1.02E+00  1.38E+00  1.47E+00
 


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
+        5.73E+06
 
 TH 2
+        6.38E+05  7.14E+04
 
 TH 3
+        2.21E+06  1.37E+01  2.14E+05
 
 TH 4
+       -5.09E+03  3.23E+05 -1.03E+03  1.45E+06
 
 TH 5
+        7.37E+05 -6.05E+01  3.67E+02 -5.27E+02  1.90E+05
 
 TH 6
+       -1.65E+02  1.71E+01  2.15E-01 -8.68E+01 -2.23E+01  2.14E+02
 
 TH 7
+       -3.68E+02  6.54E+05  1.14E+06 -1.88E+02 -4.97E+01 -3.10E-01  1.10E+02
 
 TH 8
+       -3.85E+00  4.26E+01 -5.62E+04  1.35E+02 -1.91E+04 -4.30E+00 -1.01E+01  7.41E+03
 
 TH 9
+        4.23E+01  9.53E+05  8.26E+05  1.98E+01  1.16E+01 -1.86E-01  4.38E+06  5.24E+00  3.18E+06
 
 TH10
+       -3.39E+01 -1.03E+05 -1.79E+05  1.22E+01 -1.69E+01  2.74E+01 -1.90E+06  2.31E+04 -1.38E+06  1.50E+05
 
 TH11
+        1.11E+03 -8.53E+04 -2.18E+02  3.85E+05 -1.49E+02  2.45E+01  5.55E+01  7.16E+00 -5.69E+05  1.24E+05  1.03E+05
 
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
 #CPUT: Total CPU Time in Seconds,       76.943
Stop Time:
Wed Sep 29 07:18:11 CDT 2021
