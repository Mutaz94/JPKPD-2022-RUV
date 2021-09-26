Sat Sep 25 02:01:11 CDT 2021
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
$DATA ../../../../data/int/SL3/dat20.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      982
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

 TOT. NO. OF OBS RECS:      882
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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   664.977177904467        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1493E+02  6.5050E+01  6.9757E+01  6.0649E+01  1.8326E+02 -2.5757E+00 -1.1212E+02 -4.2716E+02 -1.8779E+02 -7.2963E+01
            -8.0801E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2289.04198875700        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0720E+00  1.1302E+00  1.7634E+00  1.0415E+00  1.1647E+00  9.5014E-01  9.5637E-01  1.0751E+00  1.2754E+00  9.5984E-01
             5.1959E+00
 PARAMETER:  1.6953E-01  2.2239E-01  6.6723E-01  1.4064E-01  2.5245E-01  4.8854E-02  5.5389E-02  1.7238E-01  3.4329E-01  5.9013E-02
             1.7479E+00
 GRADIENT:   8.3237E+01 -4.3087E-01  8.7788E+00  8.3251E+00 -1.3507E+01 -1.0443E+01  1.3540E+01 -9.1193E-01  4.4515E+01 -4.7975E-01
             7.7661E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2370.48725175002        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.1894E-01  1.3416E+00  1.6697E+02  9.0420E-01  3.6145E+00  1.0815E+00  7.8042E-01  6.7484E+00  6.8758E-01  3.3034E+00
             4.1329E+00
 PARAMETER:  1.5466E-02  3.9386E-01  5.2178E+00 -7.0174E-04  1.3849E+00  1.7834E-01 -1.4792E-01  2.0093E+00 -2.7457E-01  1.2950E+00
             1.5190E+00
 GRADIENT:  -1.7441E+02 -3.3550E+01 -4.1208E-01 -2.5810E+01  9.7915E+01  8.8213E+00 -2.9577E+01 -5.5713E-03 -2.1142E+01  4.6234E+01
             5.9274E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2463.62936209639        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.3911E-01  1.2191E+00  9.8003E+01  9.4898E-01  2.9490E+00  9.9194E-01  1.0141E+00  1.8593E+01  8.7592E-01  2.8164E+00
             3.4314E+00
 PARAMETER:  3.7181E-02  2.9815E-01  4.6850E+00  4.7628E-02  1.1815E+00  9.1904E-02  1.1404E-01  3.0228E+00 -3.2485E-02  1.1355E+00
             1.3330E+00
 GRADIENT:  -1.3132E+02 -5.2375E+01 -7.4645E-01 -1.9622E+01  9.8967E+01 -1.1701E+01 -7.2006E+00  2.9853E+01  5.5641E+00  3.8310E+01
             3.6736E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2564.67149799116        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.9879E-01  9.4472E-01  1.5189E+01  1.1168E+00  1.6961E+00  1.0019E+00  6.1114E-01  4.3462E+00  1.0186E+00  1.6471E+00
             2.8990E+00
 PARAMETER:  9.8792E-02  4.3130E-02  2.8206E+00  2.1051E-01  6.2835E-01  1.0187E-01 -3.9243E-01  1.5693E+00  1.1841E-01  5.9904E-01
             1.1644E+00
 GRADIENT:   2.6141E+01 -6.1458E+00  7.4174E+00 -1.2890E+01 -2.1556E+01 -7.1406E-02  1.7350E+00  1.4507E+01  3.7554E+00 -8.2323E+00
             1.9277E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2580.87165525498        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  9.9165E-01  9.2374E-01  3.8327E+00  1.1128E+00  1.3659E+00  1.0126E+00  3.9338E-01  2.7104E+00  1.0905E+00  1.3599E+00
             2.8487E+00
 PARAMETER:  9.1611E-02  2.0670E-02  1.4436E+00  2.0689E-01  4.1182E-01  1.1250E-01 -8.3298E-01  1.0971E+00  1.8662E-01  4.0742E-01
             1.1468E+00
 GRADIENT:   1.2267E+01 -1.2743E+01 -1.4339E-01  6.9700E+00  1.0389E+01  3.5430E+00  1.1839E+00 -7.7122E+00  1.1061E+01  6.3911E+00
             1.2413E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2582.29870980233        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  9.8575E-01  1.0192E+00  3.5695E+00  1.0394E+00  1.3511E+00  1.0046E+00  3.3446E-01  3.0380E+00  1.1255E+00  1.3000E+00
             2.8328E+00
 PARAMETER:  8.5645E-02  1.1900E-01  1.3724E+00  1.3869E-01  4.0089E-01  1.0458E-01 -9.9524E-01  1.2112E+00  2.1819E-01  3.6240E-01
             1.1413E+00
 GRADIENT:  -4.7032E-01 -8.9756E-01  9.0217E-01 -1.9325E+00 -2.7694E+00  2.1392E-01  5.4514E-01 -6.4481E-02  9.9229E-01 -1.3118E+00
             1.6837E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2582.40007918745        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      578
 NPARAMETR:  9.8912E-01  1.0520E+00  3.6366E+00  1.0205E+00  1.3833E+00  1.0077E+00  3.0444E-01  3.0759E+00  1.1523E+00  1.3411E+00
             2.8341E+00
 PARAMETER:  8.9057E-02  1.5069E-01  1.3911E+00  1.2027E-01  4.2445E-01  1.0767E-01 -1.0893E+00  1.2236E+00  2.4174E-01  3.9353E-01
             1.1417E+00
 GRADIENT:   3.2580E-01  7.0024E-01 -7.8719E-02  6.0582E-01  3.8243E-01 -2.0116E-02  2.2764E-01 -1.4644E-01  2.3710E-01  1.7546E-01
             6.7517E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2582.46794651207        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      756
 NPARAMETR:  9.8875E-01  1.0192E+00  3.6854E+00  1.0397E+00  1.3685E+00  1.0073E+00  1.6255E-01  3.0992E+00  1.1506E+00  1.3221E+00
             2.8321E+00
 PARAMETER:  8.8685E-02  1.1901E-01  1.4044E+00  1.3894E-01  4.1372E-01  1.0727E-01 -1.7167E+00  1.2311E+00  2.4026E-01  3.7919E-01
             1.1410E+00
 GRADIENT:  -2.8778E-01 -3.3437E-01  7.2357E-02 -2.1579E-01 -2.1781E-01 -1.8366E-01  3.6005E-02  6.2298E-02  1.7421E-01 -4.1914E-02
            -3.8908E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2582.47770752281        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      931
 NPARAMETR:  9.8896E-01  1.0141E+00  3.6890E+00  1.0426E+00  1.3669E+00  1.0079E+00  5.9055E-02  3.0989E+00  1.1543E+00  1.3201E+00
             2.8329E+00
 PARAMETER:  8.8898E-02  1.1405E-01  1.4054E+00  1.4173E-01  4.1252E-01  1.0791E-01 -2.7293E+00  1.2310E+00  2.4353E-01  3.7770E-01
             1.1413E+00
 GRADIENT:   1.2145E-01  1.4215E-01 -7.2697E-03  3.2533E-01  7.6823E-02  3.1371E-02  1.7435E-03 -2.3261E-02  8.3394E-03  3.3173E-02
             1.9937E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2582.47855480026        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1106
 NPARAMETR:  9.8889E-01  1.0135E+00  3.6884E+00  1.0427E+00  1.3663E+00  1.0079E+00  1.0000E-02  3.1006E+00  1.1552E+00  1.3193E+00
             2.8326E+00
 PARAMETER:  8.8832E-02  1.1340E-01  1.4052E+00  1.4177E-01  4.1207E-01  1.0784E-01 -4.6277E+00  1.2316E+00  2.4423E-01  3.7710E-01
             1.1412E+00
 GRADIENT:   2.3807E-03 -1.4438E-02 -1.1845E-03 -6.6246E-03  1.6456E-03  3.4892E-03  0.0000E+00  7.0540E-04  8.2826E-03  5.0820E-04
             3.4891E-03

0ITERATION NO.:   52    OBJECTIVE VALUE:  -2582.47855509990        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  9.8889E-01  1.0135E+00  3.6884E+00  1.0426E+00  1.3663E+00  1.0079E+00  1.0000E-02  3.1006E+00  1.1552E+00  1.3193E+00
             2.8326E+00
 PARAMETER:  8.8832E-02  1.1344E-01  1.4052E+00  1.4175E-01  4.1208E-01  1.0783E-01 -4.6217E+00  1.2316E+00  2.4424E-01  3.7711E-01
             1.1412E+00
 GRADIENT:   1.7367E-03 -1.7137E-03 -2.0215E-04 -1.2303E-03  3.3031E-04  7.7591E-04  0.0000E+00  7.9336E-05  3.9464E-04  1.7141E-05
             2.6210E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1163
 NO. OF SIG. DIGITS IN FINAL EST.:  3.9
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0946E-03 -8.0641E-04 -3.4036E-02 -3.4538E-04 -2.5095E-02
 SE:             2.9382E-02  2.2968E-04  1.8177E-02  2.8125E-02  2.3144E-02
 N:                     100         100         100         100         100

 P VAL.:         9.1612E-01  4.4643E-04  6.1146E-02  9.9020E-01  2.7823E-01

 ETASHRINKSD(%)  1.5672E+00  9.9231E+01  3.9103E+01  5.7778E+00  2.2464E+01
 ETASHRINKVR(%)  3.1098E+00  9.9994E+01  6.2916E+01  1.1222E+01  3.9881E+01
 EBVSHRINKSD(%)  1.8016E+00  9.9345E+01  4.1701E+01  5.8865E+00  1.9998E+01
 EBVSHRINKVR(%)  3.5706E+00  9.9996E+01  6.6012E+01  1.1426E+01  3.5996E+01
 RELATIVEINF(%)  9.6324E+01  6.9311E-04  1.6253E+01  1.6216E+01  3.2642E+01
 EPSSHRINKSD(%)  1.6529E+01
 EPSSHRINKVR(%)  3.0326E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          882
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1621.0075725730426     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2582.4785550999050     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -961.47098252686237     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.89
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.22
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2582.479       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.89E-01  1.01E+00  3.69E+00  1.04E+00  1.37E+00  1.01E+00  1.00E-02  3.10E+00  1.16E+00  1.32E+00  2.83E+00
 


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
+        1.08E+03
 
 TH 2
+       -2.17E+01  5.07E+02
 
 TH 3
+        9.78E-02  1.05E+01  3.24E+00
 
 TH 4
+       -2.56E+01  4.76E+02 -2.31E+00  6.80E+02
 
 TH 5
+       -3.51E+00 -1.21E+02 -1.57E+01 -3.17E+01  2.17E+02
 
 TH 6
+        3.60E+00 -7.43E+00 -1.37E-01 -4.32E+00 -5.15E-01  1.81E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -1.13E-01 -5.77E+00 -1.63E+00 -4.33E-01  3.75E+00 -3.78E-02  0.00E+00  4.57E+00
 
 TH 9
+        3.86E+00 -8.63E+01  3.06E-01  9.14E+00  1.69E+00  2.14E+00  0.00E+00 -7.58E-02  1.15E+02
 
 TH10
+        1.23E+00 -1.66E+01 -2.32E+00  4.67E+00 -2.21E+01  1.14E+00  0.00E+00  2.20E+00 -5.97E-02  5.02E+01
 
 TH11
+       -1.58E+01 -2.42E+01  8.75E-02 -1.87E+01 -1.14E+00  1.54E+00  0.00E+00  2.56E+00  7.00E+00  5.99E+00  1.40E+02
 
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
 #CPUT: Total CPU Time in Seconds,       39.214
Stop Time:
Sat Sep 25 02:01:51 CDT 2021
