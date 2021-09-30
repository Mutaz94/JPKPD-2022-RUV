Thu Sep 30 00:10:20 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat33.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   170.449789046849        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3778E+02  1.9053E+02  2.5857E+02  1.3586E+02  3.2972E+02  4.5398E+01 -1.0843E+02 -2.2388E+02 -4.2668E+01 -2.7168E+02
            -3.7658E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1465.65832788557        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.6235E-01  9.4299E-01  9.1039E-01  1.0492E+00  9.0749E-01  8.5250E-01  1.0041E+00  9.9432E-01  9.4406E-01  1.0616E+00
             3.8072E+00
 PARAMETER:  6.1625E-02  4.1300E-02  6.1231E-03  1.4803E-01  2.9289E-03 -5.9583E-02  1.0407E-01  9.4302E-02  4.2434E-02  1.5981E-01
             1.4369E+00
 GRADIENT:   7.4496E-01 -1.4474E+01 -1.5957E+01 -7.0901E+00  7.1062E+00 -1.4978E+01  5.3020E+00  7.3300E+00  1.3179E+01  2.2337E+01
             9.3493E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1476.29196627110        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.6245E-01  6.5934E-01  6.6065E-01  1.2509E+00  6.2703E-01  8.6316E-01  1.3977E+00  3.0078E-01  7.7470E-01  5.6435E-01
             3.7773E+00
 PARAMETER:  6.1730E-02 -3.1651E-01 -3.1454E-01  3.2387E-01 -3.6677E-01 -4.7161E-02  4.3482E-01 -1.1014E+00 -1.5528E-01 -4.7209E-01
             1.4290E+00
 GRADIENT:  -9.7856E+00  3.7720E+01 -7.8477E+00  1.1999E+02 -8.7325E-01 -1.4393E+01  5.5340E-01  9.7340E-01 -2.0056E+00  6.1624E+00
             7.7515E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1481.61296382023        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      273
 NPARAMETR:  9.5675E-01  6.5444E-01  3.8477E-01  1.1869E+00  4.5366E-01  8.9781E-01  1.0659E+00  8.3666E-02  1.0061E+00  5.5045E-01
             3.5119E+00
 PARAMETER:  5.5789E-02 -3.2398E-01 -8.5510E-01  2.7133E-01 -6.9041E-01 -7.8008E-03  1.6378E-01 -2.3809E+00  1.0603E-01 -4.9702E-01
             1.3562E+00
 GRADIENT:  -4.6235E+01  2.4765E+01 -4.4112E+01  1.2080E+02  6.3504E+01 -7.3629E+00 -2.1283E+00  1.6003E-02  2.3509E+00 -3.3175E+00
             5.2686E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1517.29775056655        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      453
 NPARAMETR:  9.7204E-01  3.8882E-01  2.6384E-01  1.1225E+00  2.8541E-01  9.3786E-01  1.0428E+00  1.0000E-02  1.1654E+00  7.5929E-01
             2.7079E+00
 PARAMETER:  7.1644E-02 -8.4463E-01 -1.2324E+00  2.1554E-01 -1.1538E+00  3.5842E-02  1.4187E-01 -7.0331E+00  2.5305E-01 -1.7537E-01
             1.0962E+00
 GRADIENT:   2.6147E+01  3.7611E+01  4.3959E+01  4.5887E+01 -4.8990E+01  7.4262E+00 -1.7025E-01  0.0000E+00  6.4299E+00 -2.7819E+00
            -7.2484E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1528.03276339107        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      629
 NPARAMETR:  9.6424E-01  2.7559E-01  1.9063E-01  1.0170E+00  2.2392E-01  9.1387E-01  5.3149E-01  1.0000E-02  1.2166E+00  7.8212E-01
             2.8574E+00
 PARAMETER:  6.3583E-02 -1.1888E+00 -1.5574E+00  1.1689E-01 -1.3965E+00  9.9287E-03 -5.3208E-01 -1.0630E+01  2.9610E-01 -1.4575E-01
             1.1499E+00
 GRADIENT:   2.4274E-01 -2.0905E-01  6.3791E-01  2.3474E+00 -5.2426E-02  4.8741E-01  2.1710E-01  0.0000E+00  3.0493E+00  3.3021E+00
             3.2351E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1528.11354234220        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      806
 NPARAMETR:  9.6436E-01  2.8193E-01  1.8633E-01  1.0065E+00  2.2191E-01  9.1341E-01  4.0622E-01  1.0000E-02  1.2148E+00  7.7407E-01
             2.8516E+00
 PARAMETER:  6.3705E-02 -1.1661E+00 -1.5802E+00  1.0650E-01 -1.4055E+00  9.4295E-03 -8.0086E-01 -1.0651E+01  2.9454E-01 -1.5609E-01
             1.1479E+00
 GRADIENT:  -2.1332E-01  1.4118E+00  1.0073E+00 -3.3944E-01 -2.6722E+00 -1.5395E-02  1.3116E-01  0.0000E+00  9.3579E-02  4.6380E-01
             1.0777E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1528.17255183033        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      981
 NPARAMETR:  9.6438E-01  2.7952E-01  1.8566E-01  1.0059E+00  2.2130E-01  9.1309E-01  6.2794E-02  1.0000E-02  1.2156E+00  7.7754E-01
             2.8502E+00
 PARAMETER:  6.3728E-02 -1.1747E+00 -1.5838E+00  1.0587E-01 -1.4082E+00  9.0799E-03 -2.6679E+00 -1.0889E+01  2.9522E-01 -1.5162E-01
             1.1474E+00
 GRADIENT:  -1.2417E-01  1.7595E-01  2.8661E-01  2.3540E-01 -5.9239E-01 -9.2703E-02  2.4892E-03  0.0000E+00 -1.5227E-03  2.3568E-01
             1.1647E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1528.17402206984        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1159
 NPARAMETR:  9.6437E-01  2.7928E-01  1.8528E-01  1.0051E+00  2.2107E-01  9.1335E-01  1.0000E-02  1.0000E-02  1.2164E+00  7.7692E-01
             2.8496E+00
 PARAMETER:  6.3719E-02 -1.1755E+00 -1.5859E+00  1.0509E-01 -1.4093E+00  9.3640E-03 -9.2080E+00 -1.1519E+01  2.9586E-01 -1.5242E-01
             1.1472E+00
 GRADIENT:  -1.2302E-01 -1.0741E-04 -5.8857E-02  1.0016E-02 -3.0613E-02  4.5385E-03  0.0000E+00  0.0000E+00  4.0931E-03  4.0963E-02
            -5.7479E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1528.17402206984        NO. OF FUNC. EVALS.:  31
 CUMULATIVE NO. OF FUNC. EVALS.:     1190
 NPARAMETR:  9.6497E-01  2.7999E-01  1.8540E-01  1.0049E+00  2.2097E-01  9.1334E-01  1.0000E-02  1.0000E-02  1.2171E+00  7.7574E-01
             2.8497E+00
 PARAMETER:  6.3719E-02 -1.1755E+00 -1.5859E+00  1.0509E-01 -1.4093E+00  9.3640E-03 -9.2080E+00 -1.1519E+01  2.9586E-01 -1.5242E-01
             1.1472E+00
 GRADIENT:  -2.0877E-01 -5.7457E-02 -8.0504E-02  3.1103E-02  1.5190E-01  5.1918E-04  0.0000E+00  0.0000E+00 -2.1026E-02  4.0089E-02
            -3.8404E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1190
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.9754E-04 -5.4050E-05  1.5976E-04 -9.1782E-03  1.7152E-03
 SE:             2.8961E-02  8.1328E-05  2.0599E-04  2.6453E-02  2.5092E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7803E-01  5.0631E-01  4.3800E-01  7.2862E-01  9.4550E-01

 ETASHRINKSD(%)  2.9779E+00  9.9728E+01  9.9310E+01  1.1380E+01  1.5938E+01
 ETASHRINKVR(%)  5.8672E+00  9.9999E+01  9.9995E+01  2.1464E+01  2.9335E+01
 EBVSHRINKSD(%)  2.8565E+00  9.9727E+01  9.9350E+01  9.2472E+00  1.6317E+01
 EBVSHRINKVR(%)  5.6313E+00  9.9999E+01  9.9996E+01  1.7639E+01  2.9971E+01
 RELATIVEINF(%)  9.4068E+01  2.0682E-04  3.1482E-04  2.7276E+01  3.8505E+00
 EPSSHRINKSD(%)  2.6771E+01
 EPSSHRINKVR(%)  4.6376E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1528.1740220698430     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -609.23548886517028     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.87
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1528.174       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.64E-01  2.79E-01  1.85E-01  1.01E+00  2.21E-01  9.13E-01  1.00E-02  1.00E-02  1.22E+00  7.77E-01  2.85E+00
 


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
+        1.35E+03
 
 TH 2
+       -5.85E+01  1.09E+03
 
 TH 3
+       -4.48E+00  2.50E+03  1.37E+04
 
 TH 4
+       -2.16E+01  1.27E+02 -8.25E+02  6.17E+02
 
 TH 5
+        1.78E+02 -4.33E+03 -1.67E+04 -4.06E+02  2.68E+04
 
 TH 6
+        1.89E+00 -1.74E+01  3.03E+01 -1.14E+01  3.24E+00  2.11E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.05E+01 -3.54E+01  1.85E+02 -9.32E+00  1.91E+01  1.89E+00  0.00E+00  0.00E+00  8.81E+01
 
 TH10
+       -4.31E+00 -6.91E+00 -2.63E+01  1.44E+01  2.37E+01  3.76E+00  0.00E+00  0.00E+00  5.80E+00  1.62E+02
 
 TH11
+       -2.01E+01 -4.67E+00 -4.67E+01 -8.46E+00  5.27E+01  2.02E+00  0.00E+00  0.00E+00  4.19E+00  1.40E+01  5.83E+01
 
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
 #CPUT: Total CPU Time in Seconds,       27.170
Stop Time:
Thu Sep 30 00:10:59 CDT 2021
