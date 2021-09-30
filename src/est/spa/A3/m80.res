Wed Sep 29 13:49:44 CDT 2021
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
$DATA ../../../../data/spa/A3/dat80.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m80.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   362.903948397772        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3340E+02  1.4322E+02  1.1355E+02  9.2012E+01  2.1095E+02  7.5071E+01 -5.7518E+01 -6.7684E+01 -1.1788E+02 -1.6225E+02
            -3.6353E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -949.357487788707        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0927E+00  8.0107E-01  5.9651E-01  1.0947E+00  6.7554E-01  8.5261E-01  9.8828E-01  9.3015E-01  1.0837E+00  8.8077E-01
             1.8413E+00
 PARAMETER:  1.8865E-01 -1.2181E-01 -4.1666E-01  1.9044E-01 -2.9224E-01 -5.9448E-02  8.8207E-02  2.7586E-02  1.8038E-01 -2.6956E-02
             7.1045E-01
 GRADIENT:   4.5991E+02  2.4918E+01 -1.3875E+01  8.3155E+01  1.5509E+02 -1.4726E+01 -1.2630E+01  1.6905E+00 -4.3356E+01 -3.1600E+01
            -9.2415E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1114.89729385432        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.1080E+00  8.1074E-01  1.4858E-01  1.1241E+00  3.3317E-01  1.0209E+00  7.7604E-01  8.6655E-01  1.8706E+00  1.4447E-01
             1.8786E+00
 PARAMETER:  2.0259E-01 -1.0980E-01 -1.8066E+00  2.1699E-01 -9.9911E-01  1.2069E-01 -1.5355E-01 -4.3240E-02  7.2628E-01 -1.8347E+00
             7.3052E-01
 GRADIENT:   4.5919E+02  1.8543E+02  1.0519E+02  1.7075E+02  1.0458E+02  2.8361E+01 -1.3790E+01 -6.8131E+01  2.1416E+01 -4.7473E+00
            -4.5457E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1183.23461869162        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.9079E-01  4.8704E-01  1.1075E-01  1.1441E+00  2.2097E-01  8.4809E-01  9.1280E-01  7.1381E-01  1.7336E+00  6.1551E-02
             2.1158E+00
 PARAMETER:  9.0750E-02 -6.1940E-01 -2.1005E+00  2.3461E-01 -1.4097E+00 -6.4764E-02  8.7610E-03 -2.3713E-01  6.5022E-01 -2.6879E+00
             8.4943E-01
 GRADIENT:   7.5771E+01  1.2605E+02  6.0790E+01  1.6472E+02  5.0456E+01 -1.4080E+01 -1.3579E+01 -5.5430E+01 -6.3549E+01 -1.0461E+00
            -3.2360E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1253.11242100545        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.6650E-01  3.1668E-01  9.6352E-02  9.4243E-01  1.8249E-01  8.6961E-01  8.5581E-01  1.5811E+00  1.7779E+00  1.8861E-02
             2.8017E+00
 PARAMETER:  6.5924E-02 -1.0499E+00 -2.2398E+00  4.0708E-02 -1.6011E+00 -3.9705E-02 -5.5709E-02  5.5813E-01  6.7546E-01 -3.8707E+00
             1.1302E+00
 GRADIENT:  -9.7210E+01 -1.9671E+01 -6.4005E+00  5.3779E+01  4.2734E+00 -1.1570E+01 -2.3554E+00  1.1038E+01 -2.9746E+01 -1.6993E-01
            -9.2841E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1254.75193118592        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0013E+00  3.0528E-01  8.7205E-02  8.9931E-01  1.7560E-01  8.5622E-01  8.7508E-01  1.5879E+00  1.9342E+00  1.2530E-02
             3.0493E+00
 PARAMETER:  1.0132E-01 -1.0865E+00 -2.3395E+00 -6.1281E-03 -1.6395E+00 -5.5229E-02 -3.3435E-02  5.6242E-01  7.5967E-01 -4.2796E+00
             1.2149E+00
 GRADIENT:  -1.0119E+01 -3.2016E+01 -2.5226E+01  3.5929E+01  3.9585E+01 -1.2746E+01 -7.0395E-01  1.2739E+01 -1.3314E+01 -7.3051E-02
            -5.4144E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1255.09383933014        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  1.0108E+00  3.1180E-01  8.1757E-02  8.5807E-01  1.7205E-01  8.6536E-01  8.7185E-01  1.5393E+00  2.0598E+00  1.0000E-02
             3.2271E+00
 PARAMETER:  1.1075E-01 -1.0654E+00 -2.4040E+00 -5.3071E-02 -1.6600E+00 -4.4611E-02 -3.7136E-02  5.3133E-01  8.2261E-01 -4.5366E+00
             1.2716E+00
 GRADIENT:   1.7115E+01 -2.2016E+01 -2.4486E+01  2.1080E+01  3.4800E+01 -1.0244E+01  7.7853E-01  1.1001E+01 -3.3930E+00 -1.4397E-02
            -2.8738E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1260.84063170479        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      516
 NPARAMETR:  9.9747E-01  3.4318E-01  8.5475E-02  7.9012E-01  1.7476E-01  9.1509E-01  6.3371E-01  1.1795E+00  2.1377E+00  1.0000E-02
             3.5071E+00
 PARAMETER:  9.7471E-02 -9.6951E-01 -2.3595E+00 -1.3557E-01 -1.6443E+00  1.1266E-02 -3.5616E-01  2.6507E-01  8.5974E-01 -4.5811E+00
             1.3548E+00
 GRADIENT:  -2.0825E+01  2.5490E+01  2.2656E+01 -1.5825E+01 -4.5360E+01  5.2327E+00 -2.7552E+00 -6.9537E+00  4.3562E+00  0.0000E+00
             1.4047E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1261.04536986834        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      588
 NPARAMETR:  1.0028E+00  3.3514E-01  8.4119E-02  8.0602E-01  1.7371E-01  9.0343E-01  7.0322E-01  1.2879E+00  2.1168E+00  1.0000E-02
             3.4382E+00
 PARAMETER:  1.0275E-01 -9.9320E-01 -2.3755E+00 -1.1564E-01 -1.6504E+00 -1.5548E-03 -2.5209E-01  3.5300E-01  8.4991E-01 -4.5893E+00
             1.3349E+00
 GRADIENT:  -6.5356E+00  1.0356E+01  8.9085E+00 -7.7400E+00 -1.8534E+01  2.3487E+00 -1.5272E+00 -3.1040E+00  2.2763E+00  0.0000E+00
             6.7307E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1261.84867019849        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      662
 NPARAMETR:  1.0044E+00  3.3311E-01  8.4788E-02  8.1194E-01  1.7435E-01  9.0099E-01  7.2825E-01  1.3045E+00  2.0801E+00  1.0000E-02
             3.4173E+00
 PARAMETER:  1.0435E-01 -9.9927E-01 -2.3676E+00 -1.0833E-01 -1.6467E+00 -4.2643E-03 -2.1711E-01  3.6581E-01  8.3240E-01 -4.5550E+00
             1.3288E+00
 GRADIENT:  -3.1869E+00  3.9079E+00  4.8050E+00 -6.2073E+00 -6.0095E+00  1.7553E+00 -1.0563E+00 -2.4390E+00  7.8461E-01  0.0000E+00
             5.6609E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1290.56809925379        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      841
 NPARAMETR:  1.0154E+00  4.8323E-01  1.3193E-01  9.2710E-01  2.3803E-01  8.6825E-01  6.8415E-01  1.3987E+00  1.6597E+00  2.9520E-02
             3.0858E+00
 PARAMETER:  1.1529E-01 -6.2727E-01 -1.9255E+00  2.4302E-02 -1.3354E+00 -4.1279E-02 -2.7958E-01  4.3556E-01  6.0661E-01 -3.4227E+00
             1.2268E+00
 GRADIENT:   1.1376E+01  4.2023E+00  1.2516E+01 -2.6638E+00 -2.9179E+01 -3.0180E+00  8.5556E-01 -2.0748E+00  1.2113E+01 -2.4499E-02
            -5.8605E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1291.93864126686        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1016
 NPARAMETR:  1.0097E+00  5.3908E-01  1.3584E-01  9.2451E-01  2.5764E-01  8.7343E-01  6.3180E-01  1.4215E+00  1.5133E+00  3.8895E-02
             3.1463E+00
 PARAMETER:  1.0967E-01 -5.1789E-01 -1.8962E+00  2.1507E-02 -1.2562E+00 -3.5326E-02 -3.5919E-01  4.5175E-01  5.1426E-01 -3.1469E+00
             1.2462E+00
 GRADIENT:   2.1599E-01 -8.3997E-01  3.8973E-01 -4.8656E-01  1.2496E+00  3.2038E-01  3.1876E-01 -4.5639E-01 -1.6634E-01 -1.4842E-02
             1.9922E-01

0ITERATION NO.:   59    OBJECTIVE VALUE:  -1291.94401428777        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1144
 NPARAMETR:  1.0096E+00  5.3966E-01  1.3560E-01  9.2479E-01  2.5751E-01  8.7243E-01  6.2139E-01  1.4281E+00  1.5167E+00  3.9315E-02
             3.1447E+00
 PARAMETER:  1.0956E-01 -5.1682E-01 -1.8980E+00  2.1811E-02 -1.2567E+00 -3.6470E-02 -3.7580E-01  4.5635E-01  5.1653E-01 -3.1362E+00
             1.2457E+00
 GRADIENT:   2.6696E-01  1.1820E+00  6.7855E-01  2.5358E-01 -2.1267E+00 -4.2276E-02 -8.9096E-02 -1.3316E-01 -2.9296E-01 -1.8741E-02
            -3.2399E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1144
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.9153E-03 -1.5053E-02  1.0510E-02 -6.9200E-03  1.9776E-03
 SE:             2.8452E-02  1.3461E-02  1.9337E-02  2.5993E-02  1.5475E-03
 N:                     100         100         100         100         100

 P VAL.:         9.1839E-01  2.6347E-01  5.8677E-01  7.9007E-01  2.0127E-01

 ETASHRINKSD(%)  4.6831E+00  5.4904E+01  3.5219E+01  1.2919E+01  9.4816E+01
 ETASHRINKVR(%)  9.1469E+00  7.9663E+01  5.8034E+01  2.4169E+01  9.9731E+01
 EBVSHRINKSD(%)  4.6770E+00  5.4197E+01  3.5603E+01  1.1207E+01  9.5193E+01
 EBVSHRINKVR(%)  9.1353E+00  7.9021E+01  5.8531E+01  2.1159E+01  9.9769E+01
 RELATIVEINF(%)  8.5137E+01  9.3284E-01  1.1383E+01  5.8714E+01  7.8433E-03
 EPSSHRINKSD(%)  3.2088E+01
 EPSSHRINKVR(%)  5.3880E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1291.9440142877665     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -556.79318772402837     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.02
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1291.944       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  5.40E-01  1.36E-01  9.25E-01  2.58E-01  8.72E-01  6.21E-01  1.43E+00  1.52E+00  3.93E-02  3.14E+00
 


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
+        1.31E+03
 
 TH 2
+       -3.81E+01  2.21E+03
 
 TH 3
+       -5.18E+02  3.07E+03  1.27E+04
 
 TH 4
+       -4.13E+01  1.87E+02 -4.72E+02  4.44E+02
 
 TH 5
+        4.48E+02 -7.31E+03 -1.27E+04 -3.53E+02  2.70E+04
 
 TH 6
+       -1.27E+00 -2.04E+01  2.41E+01 -1.09E+01  1.16E+02  2.16E+02
 
 TH 7
+       -3.16E+00 -4.32E+01 -8.75E+01 -4.19E+00  2.07E+02  7.22E-01  3.24E+01
 
 TH 8
+        2.27E+00 -3.28E-01 -7.38E+01 -1.34E+00 -1.04E+01  5.09E+00  4.38E+00  2.42E+01
 
 TH 9
+        1.68E+01 -4.76E+01  9.54E+01 -1.25E+01  2.06E+02  5.39E+00  8.18E+00 -1.56E+00  4.43E+01
 
 TH10
+       -1.32E-01 -1.05E+01 -1.33E+01 -6.03E-01  5.65E+01  7.17E-01  4.94E+00  2.26E+00  1.12E+00 -1.14E+01
 
 TH11
+       -2.30E+01 -2.05E+01 -1.22E+00  3.92E-01  5.26E+01  6.28E-01  8.09E+00  5.55E+00  7.11E+00  1.73E+00  2.91E+01
 
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
 #CPUT: Total CPU Time in Seconds,       18.720
Stop Time:
Wed Sep 29 13:50:04 CDT 2021
