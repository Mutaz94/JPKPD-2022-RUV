Wed Sep 29 07:03:01 CDT 2021
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
$DATA ../../../../data/int/TD2/dat14.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m14.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3594.48284983908        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3499E+02 -2.6014E+01  4.8851E+00  1.2861E+02  1.7355E+02  3.0354E+01 -2.1408E+01 -4.1776E+01 -1.4255E+01  1.1266E+01
            -5.3726E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3657.10621919552        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      129
 NPARAMETR:  9.6173E-01  1.0483E+00  9.4678E-01  9.9031E-01  9.2216E-01  1.0260E+00  1.0290E+00  1.0401E+00  1.0127E+00  9.2840E-01
             1.2592E+00
 PARAMETER:  6.0979E-02  1.4716E-01  4.5315E-02  9.0261E-02  1.8959E-02  1.2566E-01  1.2862E-01  1.3935E-01  1.1263E-01  2.5710E-02
             3.3048E-01
 GRADIENT:  -9.6392E+01 -1.2355E+01 -9.5561E+00  6.3100E+00 -1.4463E+01 -1.3117E+01 -5.1415E+00 -3.0154E+00 -1.3450E+00  7.9201E+00
             2.3022E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3660.97206428506        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.7919E-01  1.1797E+00  1.2753E+00  9.2071E-01  1.0906E+00  1.0079E+00  1.0233E+00  1.9006E+00  1.1453E+00  9.8853E-01
             1.2369E+00
 PARAMETER:  7.8971E-02  2.6524E-01  3.4318E-01  1.7389E-02  1.8671E-01  1.0790E-01  1.2305E-01  7.4217E-01  2.3564E-01  8.8468E-02
             3.1262E-01
 GRADIENT:  -6.0693E+01 -5.3051E+00  8.5714E+00 -1.2152E+01 -6.5012E+00 -1.8166E+01  5.4921E+00  9.0427E-01  2.9631E+01  2.8089E+00
             3.0522E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3663.08745366446        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      481
 NPARAMETR:  9.9283E-01  1.1372E+00  1.5684E+00  9.7099E-01  1.1189E+00  1.0597E+00  9.0368E-01  2.9689E+00  9.5742E-01  9.6667E-01
             1.2411E+00
 PARAMETER:  9.2804E-02  2.2861E-01  5.5003E-01  7.0565E-02  2.1239E-01  1.5799E-01 -1.2840E-03  1.1882E+00  5.6482E-02  6.6099E-02
             3.1599E-01
 GRADIENT:  -2.7772E+01 -4.3332E+00 -6.3793E-01 -8.8560E+00 -1.1768E+01  3.6143E+00 -1.0316E+01  2.9466E+01  1.3216E+01 -1.2366E+00
             5.6374E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3663.61669465946        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      665
 NPARAMETR:  9.9415E-01  1.1334E+00  1.5721E+00  9.7081E-01  1.1196E+00  1.0444E+00  9.0720E-01  2.9120E+00  9.5480E-01  9.6566E-01
             1.2434E+00
 PARAMETER:  9.4136E-02  2.2522E-01  5.5240E-01  7.0375E-02  2.1294E-01  1.4346E-01  2.6081E-03  1.1689E+00  5.3748E-02  6.5059E-02
             3.1788E-01
 GRADIENT:  -2.5983E+01 -9.4968E+00 -5.3017E-02 -1.1619E+01 -9.0771E+00 -2.0273E+00 -1.0294E+01  2.6600E+01  1.2778E+01 -5.1069E-01
             5.8905E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3667.36833384182        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      847            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0005E+00  1.1369E+00  1.5534E+00  9.7192E-01  1.1222E+00  1.0506E+00  9.4826E-01  2.4867E+00  9.3365E-01  9.6825E-01
             1.2239E+00
 PARAMETER:  1.0053E-01  2.2827E-01  5.4042E-01  7.1516E-02  2.1526E-01  1.4939E-01  4.6876E-02  1.0110E+00  3.1348E-02  6.7733E-02
             3.0204E-01
 GRADIENT:   2.8954E+02  1.2341E+02  1.8400E+01  6.0480E+01  6.1342E+01  5.1915E+01 -3.1981E+00  1.9184E+01  1.2144E+01  8.1825E+00
             1.8246E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3667.75429241802        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      990
 NPARAMETR:  1.0007E+00  1.1356E+00  1.5498E+00  9.7149E-01  1.1226E+00  1.0386E+00  1.0607E+00  2.4724E+00  8.6248E-01  9.2888E-01
             1.2236E+00
 PARAMETER:  1.0066E-01  2.2716E-01  5.3815E-01  7.1076E-02  2.1567E-01  1.3790E-01  1.5894E-01  1.0052E+00 -4.7949E-02  2.6221E-02
             3.0182E-01
 GRADIENT:  -1.2748E+01  5.1242E+00  2.9370E+00 -1.1249E+01  8.2152E+00 -4.0274E+00  3.3832E+00 -1.7310E+00 -4.7731E+00 -5.0500E+00
             2.2353E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3668.01959623861        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1167
 NPARAMETR:  1.0008E+00  1.1356E+00  1.5491E+00  9.7152E-01  1.1225E+00  1.0461E+00  1.0189E+00  2.4747E+00  8.8807E-01  9.5993E-01
             1.2226E+00
 PARAMETER:  1.0077E-01  2.2718E-01  5.3768E-01  7.1108E-02  2.1558E-01  1.4505E-01  1.1869E-01  1.0061E+00 -1.8708E-02  5.9103E-02
             3.0098E-01
 GRADIENT:  -1.2251E+01  4.0727E-02  3.8177E+00 -6.9046E+00  5.4373E+00 -1.1448E+00 -3.9641E-01 -9.1895E-01  1.0089E-01 -7.2810E-01
             2.1273E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3668.13927858797        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1350
 NPARAMETR:  1.0072E+00  1.1352E+00  1.5504E+00  9.7486E-01  1.1229E+00  1.0492E+00  1.0222E+00  2.4760E+00  8.8750E-01  9.6129E-01
             1.2184E+00
 PARAMETER:  1.0713E-01  2.2677E-01  5.3849E-01  7.4539E-02  2.1593E-01  1.4801E-01  1.2197E-01  1.0066E+00 -1.9352E-02  6.0517E-02
             2.9751E-01
 GRADIENT:   7.9328E-01  2.6840E+00  3.4241E+00 -9.9786E-01  7.1962E+00  1.0954E-01 -1.2075E-01 -1.3819E+00  3.7882E-01 -4.3900E-01
             1.4534E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3668.19444136665        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1535
 NPARAMETR:  1.0087E+00  1.1275E+00  1.5506E+00  9.8011E-01  1.1219E+00  1.0506E+00  1.0310E+00  2.4843E+00  8.7826E-01  9.6125E-01
             1.2146E+00
 PARAMETER:  1.0864E-01  2.2003E-01  5.3865E-01  7.9912E-02  2.1505E-01  1.4936E-01  1.3057E-01  1.0100E+00 -2.9818E-02  6.0484E-02
             2.9438E-01
 GRADIENT:   3.9006E+00  5.3879E-01  2.0359E+00  2.0722E+00  1.0877E+01  6.3471E-01  1.8608E-01 -1.7642E+00 -5.9546E-01  1.5803E-01
             9.4015E+00

0ITERATION NO.:   46    OBJECTIVE VALUE:  -3668.19444136665        NO. OF FUNC. EVALS.:  26
 CUMULATIVE NO. OF FUNC. EVALS.:     1561
 NPARAMETR:  1.0076E+00  1.1282E+00  1.5529E+00  9.8037E-01  1.1224E+00  1.0498E+00  1.0299E+00  2.4792E+00  8.7802E-01  9.6099E-01
             1.2152E+00
 PARAMETER:  1.0864E-01  2.2003E-01  5.3865E-01  7.9912E-02  2.1505E-01  1.4936E-01  1.3057E-01  1.0100E+00 -2.9818E-02  6.0484E-02
             2.9438E-01
 GRADIENT:   1.9918E+00 -3.4162E+01 -3.9766E+02 -7.3904E+01 -2.8642E+01  2.7766E-01  1.3343E-01  7.5842E+00  7.3926E+01  4.0376E+01
            -2.0622E+01
 NUMSIGDIG:         1.7         2.3         2.3         2.3         2.5         2.0         1.8         2.4         2.3         2.3
                    2.5

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1561
 NO. OF SIG. DIGITS IN FINAL EST.:  1.7
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -6.9053E-04 -2.5915E-02 -3.0254E-02  1.5829E-02 -4.6363E-02
 SE:             2.9840E-02  2.2775E-02  2.3387E-02  2.5237E-02  2.2568E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8154E-01  2.5518E-01  1.9578E-01  5.3052E-01  3.9935E-02

 ETASHRINKSD(%)  3.1524E-02  2.3700E+01  2.1652E+01  1.5452E+01  2.4396E+01
 ETASHRINKVR(%)  6.3038E-02  4.1784E+01  3.8616E+01  2.8517E+01  4.2840E+01
 EBVSHRINKSD(%)  3.4217E-01  2.3991E+01  2.3491E+01  1.7536E+01  2.2938E+01
 EBVSHRINKVR(%)  6.8317E-01  4.2226E+01  4.1464E+01  3.1997E+01  4.0615E+01
 RELATIVEINF(%)  9.9314E+01  2.8287E+01  4.7545E+01  3.6946E+01  3.5106E+01
 EPSSHRINKSD(%)  2.1844E+01
 EPSSHRINKVR(%)  3.8916E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3668.1944413666502     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2014.1050815982394     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    50.47
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.19
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3668.194       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.13E+00  1.55E+00  9.80E-01  1.12E+00  1.05E+00  1.03E+00  2.48E+00  8.78E-01  9.61E-01  1.21E+00
 


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
+        9.83E+02
 
 TH 2
+       -4.79E+01  3.66E+03
 
 TH 3
+        4.34E+00  8.55E+01  7.67E+03
 
 TH 4
+       -4.68E+01  7.17E+02 -2.33E+02  2.09E+04
 
 TH 5
+       -3.50E+00 -1.05E+02 -2.30E+02  3.95E+01  4.38E+03
 
 TH 6
+       -1.11E+00  1.92E+01  3.00E+00  1.78E+01 -3.17E-02  1.78E+02
 
 TH 7
+        2.68E-01  1.55E+02  2.45E+01  1.71E+02  4.63E+01 -1.14E-01  1.04E+04
 
 TH 8
+        2.98E-01 -5.76E+00 -7.10E+01 -1.67E+00 -1.36E+02 -1.64E-01 -3.77E+00  5.59E+01
 
 TH 9
+        1.95E+04 -8.85E+03  9.56E+01 -1.57E+02  2.50E+01 -1.65E+02 -1.59E+04  3.94E+00  2.46E+04
 
 TH10
+        1.55E+01 -3.00E+02  1.09E+02  9.77E+01 -3.33E+02 -6.73E+00 -1.18E+02  3.42E+01  4.75E+01  2.25E+04
 
 TH11
+       -8.58E+00  2.85E+01  9.57E+01 -7.48E+01  5.35E+02  2.03E+00  4.68E+01 -7.44E+01  2.07E+01 -1.80E+02  2.53E+03
 
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
 #CPUT: Total CPU Time in Seconds,       66.776
Stop Time:
Wed Sep 29 07:04:09 CDT 2021
