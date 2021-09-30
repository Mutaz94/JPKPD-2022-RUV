Wed Sep 29 12:51:01 CDT 2021
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
$DATA ../../../../data/spa/A2/dat48.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m48.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -870.563700748527        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5566E+02  4.6059E+01  7.8070E+01  4.6792E+00  6.3583E+01  1.8347E+01 -2.4135E+01 -3.5744E+01 -2.9937E+01 -4.9905E+01
            -1.3968E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1354.88460992979        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0867E+00  9.9268E-01  9.3373E-01  1.0510E+00  9.2398E-01  1.1228E+00  1.0183E+00  1.0203E+00  1.0455E+00  9.7203E-01
             2.2196E+00
 PARAMETER:  1.8317E-01  9.2652E-02  3.1436E-02  1.4975E-01  2.0938E-02  2.1578E-01  1.1811E-01  1.2012E-01  1.4450E-01  7.1627E-02
             8.9735E-01
 GRADIENT:   2.9949E+02  1.6627E+01  2.1478E+01  8.4616E+00 -1.3651E+01  3.1899E+01  3.2667E-01 -2.4031E+00  2.4927E+00  6.2301E+00
            -1.1651E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1361.58446838007        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      201
 NPARAMETR:  1.0842E+00  9.4550E-01  5.8457E-01  1.0619E+00  7.0499E-01  1.0833E+00  1.1216E+00  6.9155E-01  9.1382E-01  4.1106E-01
             2.2652E+00
 PARAMETER:  1.8085E-01  4.3956E-02 -4.3688E-01  1.6007E-01 -2.4957E-01  1.8005E-01  2.1475E-01 -2.6881E-01  9.8789E-03 -7.8903E-01
             9.1767E-01
 GRADIENT:   1.3791E+02  1.5426E+01  1.2076E+00  2.3305E+01  2.1292E+01  4.8911E+00 -2.8788E+00 -4.0421E+00 -1.5367E+01 -4.3176E+00
            -1.3374E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1379.68500288436        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      379
 NPARAMETR:  1.0089E+00  8.3408E-01  6.7942E-01  1.1251E+00  7.0628E-01  1.0327E+00  1.1921E+00  7.9238E-01  9.0121E-01  3.1123E-01
             2.8115E+00
 PARAMETER:  1.0882E-01 -8.1427E-02 -2.8651E-01  2.1783E-01 -2.4774E-01  1.3216E-01  2.7575E-01 -1.3271E-01 -4.0168E-03 -1.0672E+00
             1.1337E+00
 GRADIENT:  -3.0622E+00  1.3224E+01  1.5827E+01 -6.5011E+00 -2.6748E+01  2.7350E+00  2.6137E+00  2.4091E-02 -7.2969E-01  8.5970E-01
             9.6491E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1380.33204622829        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      554
 NPARAMETR:  1.0077E+00  6.5468E-01  7.8958E-01  1.2481E+00  7.1803E-01  1.0230E+00  1.3640E+00  8.9834E-01  8.4070E-01  2.1482E-01
             2.8128E+00
 PARAMETER:  1.0771E-01 -3.2361E-01 -1.3626E-01  3.2161E-01 -2.3124E-01  1.2270E-01  4.1044E-01 -7.2104E-03 -7.3520E-02 -1.4380E+00
             1.1342E+00
 GRADIENT:  -6.3324E-02  4.0173E+00  2.2798E+00  5.2880E+00 -3.1350E+00  1.9529E-01  3.1177E-01  1.4339E-01 -3.1246E+00  2.3304E-01
             2.9561E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1380.54057516063        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      729
 NPARAMETR:  1.0070E+00  5.6055E-01  7.1872E-01  1.2849E+00  6.4521E-01  1.0234E+00  1.5189E+00  8.6348E-01  8.4541E-01  1.2589E-01
             2.7688E+00
 PARAMETER:  1.0702E-01 -4.7883E-01 -2.3028E-01  3.5071E-01 -3.3818E-01  1.2309E-01  5.1799E-01 -4.6779E-02 -6.7931E-02 -1.9724E+00
             1.1184E+00
 GRADIENT:   1.5625E-01  4.8560E-01  6.8299E-01  4.9768E-01 -1.2853E+00 -1.8294E-02 -1.3953E-01 -5.4935E-02  1.9644E-02 -1.6818E-03
             2.8632E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1380.54170172457        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      907             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0068E+00  5.4888E-01  7.2402E-01  1.2915E+00  6.4547E-01  1.0232E+00  1.5459E+00  8.7061E-01  8.4209E-01  1.3126E-01
             2.7673E+00
 PARAMETER:  1.0674E-01 -4.9988E-01 -2.2294E-01  3.5582E-01 -3.3778E-01  1.2292E-01  5.3560E-01 -3.8565E-02 -7.1870E-02 -1.9306E+00
             1.1179E+00
 GRADIENT:   5.6466E+01  6.9563E+00  1.3189E+00  5.0973E+01  6.2482E+00  6.2388E+00  1.5840E+00  2.2108E-01  1.4897E+00  5.6968E-02
             9.1666E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1380.54208238848        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1083
 NPARAMETR:  1.0066E+00  5.4841E-01  7.2430E-01  1.2922E+00  6.4543E-01  1.0232E+00  1.5441E+00  8.6516E-01  8.4247E-01  1.4185E-01
             2.7673E+00
 PARAMETER:  1.0659E-01 -5.0073E-01 -2.2254E-01  3.5634E-01 -3.3784E-01  1.2292E-01  5.3445E-01 -4.4841E-02 -7.1415E-02 -1.8530E+00
             1.1179E+00
 GRADIENT:  -9.4158E-02  8.2776E-02 -6.7728E-02 -6.2799E-02  4.6140E-02 -5.6975E-03 -5.7222E-03  4.6685E-03  6.8573E-03  5.7265E-04
             8.2553E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1380.54210015942        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     1142
 NPARAMETR:  1.0066E+00  5.4788E-01  7.2447E-01  1.2925E+00  6.4540E-01  1.0232E+00  1.5440E+00  8.6378E-01  8.4254E-01  1.4733E-01
             2.7669E+00
 PARAMETER:  1.0661E-01 -5.0169E-01 -2.2231E-01  3.5657E-01 -3.3788E-01  1.2295E-01  5.3439E-01 -4.6432E-02 -7.1336E-02 -1.8151E+00
             1.1177E+00
 GRADIENT:  -2.3498E-01  4.3754E-02 -8.3460E-02  3.7578E-01  1.3598E-01 -1.7637E-02 -2.5966E-02 -9.8534E-03  8.0734E-03 -4.3076E-04
             8.9624E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1142
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1475E-03  9.1057E-03 -1.0384E-02 -1.2323E-02 -1.7659E-03
 SE:             2.9090E-02  1.5644E-02  1.5131E-02  2.3350E-02  3.7895E-03
 N:                     100         100         100         100         100

 P VAL.:         9.6853E-01  5.6054E-01  4.9257E-01  5.9768E-01  6.4122E-01

 ETASHRINKSD(%)  2.5458E+00  4.7589E+01  4.9308E+01  2.1774E+01  8.7305E+01
 ETASHRINKVR(%)  5.0267E+00  7.2531E+01  7.4303E+01  3.8806E+01  9.8388E+01
 EBVSHRINKSD(%)  2.6142E+00  5.0475E+01  4.8930E+01  2.0994E+01  8.7312E+01
 EBVSHRINKVR(%)  5.1601E+00  7.5473E+01  7.3919E+01  3.7581E+01  9.8390E+01
 RELATIVEINF(%)  9.2384E+01  1.2683E+00  1.2863E+00  5.3043E+00  8.1744E-02
 EPSSHRINKSD(%)  2.9421E+01
 EPSSHRINKVR(%)  5.0186E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1380.5421001594243     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -645.39127359568613     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.31
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1380.542       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  5.48E-01  7.24E-01  1.29E+00  6.45E-01  1.02E+00  1.54E+00  8.64E-01  8.43E-01  1.47E-01  2.77E+00
 


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
+        9.94E+02
 
 TH 2
+       -5.08E+01  3.83E+02
 
 TH 3
+        1.28E+01  3.06E+02  7.68E+02
 
 TH 4
+       -5.15E+01  3.27E+02 -8.64E+01  5.83E+02
 
 TH 5
+        2.86E+01 -5.75E+02 -1.16E+03 -3.05E+01  1.86E+03
 
 TH 6
+       -2.38E-02 -9.00E+00  6.49E+00 -1.43E+01  1.53E+00  1.69E+02
 
 TH 7
+        2.81E-01  1.30E+01  1.58E+00 -2.29E-01  1.77E+00  5.12E-01  8.29E+00
 
 TH 8
+       -2.72E+00 -1.26E+01 -5.00E+01 -1.37E+00  4.65E+01 -7.23E-01  1.50E+00  2.63E+01
 
 TH 9
+        6.37E+00 -1.26E+01 -1.36E+00 -1.12E+01  1.50E+01  2.14E+00  1.09E+01  1.53E+00  1.15E+02
 
 TH10
+       -1.11E+00 -1.31E+00 -1.05E+01 -2.70E+00  1.17E+01 -5.25E-01  1.10E+00  6.58E+00 -1.96E-01  2.51E+00
 
 TH11
+       -1.24E+01  1.16E+00 -8.02E+00 -1.05E+01 -1.86E+01  3.21E+00  3.70E+00  1.15E+01  9.14E+00  5.07E+00  4.15E+01
 
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
 #CPUT: Total CPU Time in Seconds,       22.163
Stop Time:
Wed Sep 29 12:51:25 CDT 2021
