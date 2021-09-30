Thu Sep 30 08:13:36 CDT 2021
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
$DATA ../../../../data/spa2/TD2/dat60.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1942.13974762993        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0807E+02 -8.4993E+01  7.3064E+01  5.5638E+00  9.2906E+01  4.7030E+01 -2.4552E+01 -2.3190E+02 -2.5882E+01  2.1603E+01
            -8.0771E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2368.05063505299        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:       94
 NPARAMETR:  8.1877E-01  1.0613E+00  9.4790E-01  9.7871E-01  9.4037E-01  1.0207E+00  1.0278E+00  1.1744E+00  1.0056E+00  9.6870E-01
             1.5704E+00
 PARAMETER: -9.9957E-02  1.5950E-01  4.6496E-02  7.8476E-02  3.8514E-02  1.2049E-01  1.2747E-01  2.6072E-01  1.0558E-01  6.8204E-02
             5.5132E-01
 GRADIENT:  -2.5283E+02 -1.1435E+01 -1.8257E+00  4.0779E+00 -2.6375E+01 -1.4900E+01 -1.1009E+01  1.7614E+01  2.4200E+01  2.7704E+01
             3.5462E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2399.08682599065        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  8.3769E-01  9.5529E-01  6.2168E-01  1.0004E+00  7.4587E-01  9.9972E-01  1.4389E+00  3.9595E-01  1.0118E+00  5.0891E-01
             1.4790E+00
 PARAMETER: -7.7104E-02  5.4255E-02 -3.7533E-01  1.0035E-01 -1.9320E-01  9.9718E-02  4.6385E-01 -8.2647E-01  1.1172E-01 -5.7548E-01
             4.9138E-01
 GRADIENT:  -2.0725E+02  1.4712E+01 -2.2467E+01  1.6226E+01  9.9621E-02 -2.3891E+01  2.4600E+01  3.8117E+00  4.1075E+01  2.5586E+00
             3.1157E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2500.93117144302        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      289
 NPARAMETR:  1.0349E+00  1.0851E+00  7.5780E-01  9.7272E-01  8.9750E-01  9.2493E-01  1.2526E+00  1.9334E-01  9.2771E-01  7.6031E-01
             1.1308E+00
 PARAMETER:  1.3430E-01  1.8171E-01 -1.7734E-01  7.2341E-02 -8.1454E-03  2.1957E-02  3.2523E-01 -1.5433E+00  2.4962E-02 -1.7402E-01
             2.2293E-01
 GRADIENT:   7.6536E+01 -1.0187E+01 -2.0466E+01  1.7733E+01  2.3017E+01 -4.7708E+01 -2.9520E+00 -4.8664E-01  7.8320E+00 -8.7859E-01
             1.0134E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2512.27742998454        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      465
 NPARAMETR:  9.9147E-01  1.1019E+00  8.1636E-01  9.6645E-01  9.1879E-01  1.0283E+00  1.2871E+00  5.6767E-01  9.1178E-01  7.8448E-01
             9.9952E-01
 PARAMETER:  9.1430E-02  1.9707E-01 -1.0290E-01  6.5879E-02  1.5301E-02  1.2795E-01  3.5238E-01 -4.6621E-01  7.6478E-03 -1.4274E-01
             9.9518E-02
 GRADIENT:  -2.5457E+01  4.8409E+00  2.5132E+00  2.2051E+00 -7.1281E+00  5.2568E-01  3.4565E+00 -2.7898E+00  2.5035E+00  1.3840E+00
            -1.0217E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2512.60129179049        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  1.0038E+00  1.1029E+00  8.1781E-01  9.6297E-01  9.2432E-01  1.0269E+00  1.2657E+00  5.8513E-01  9.0428E-01  7.8334E-01
             9.9936E-01
 PARAMETER:  1.0376E-01  1.9798E-01 -1.0113E-01  6.2267E-02  2.1300E-02  1.2658E-01  3.3565E-01 -4.3592E-01 -6.1507E-04 -1.4419E-01
             9.9358E-02
 GRADIENT:   4.2093E-01 -7.9101E-01 -2.2413E-01 -2.4945E-01 -9.2835E-01  2.9219E-01  2.7960E-01 -2.6171E+00  1.8903E-01  2.3107E-01
            -1.0275E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2512.69268304911        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      823
 NPARAMETR:  1.0043E+00  1.1024E+00  8.1803E-01  9.6273E-01  9.2458E-01  1.0284E+00  1.2667E+00  6.0375E-01  9.0398E-01  7.8002E-01
             1.0009E+00
 PARAMETER:  1.0428E-01  1.9753E-01 -1.0085E-01  6.2017E-02  2.1582E-02  1.2799E-01  3.3643E-01 -4.0459E-01 -9.5097E-04 -1.4843E-01
             1.0088E-01
 GRADIENT:   1.5077E+00 -1.5803E+00 -1.5135E+00 -6.7640E-01 -1.4351E-01  8.7539E-01  3.4265E-01 -2.2337E+00  1.4802E-01  6.5602E-02
            -7.7284E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2512.72110055934        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1003
 NPARAMETR:  1.0037E+00  1.1055E+00  8.2111E-01  9.6243E-01  9.2669E-01  1.0260E+00  1.2637E+00  6.0253E-01  9.0352E-01  7.8120E-01
             1.0076E+00
 PARAMETER:  1.0368E-01  2.0026E-01 -9.7103E-02  6.1704E-02  2.3869E-02  1.2568E-01  3.3401E-01 -4.0662E-01 -1.4523E-03 -1.4693E-01
             1.0757E-01
 GRADIENT:   1.2454E-01 -2.7139E-01 -3.2735E-02  5.5610E-01 -1.3030E+00 -6.2481E-02  9.1996E-02 -2.2379E+00  8.5407E-02 -1.0191E-01
            -6.0169E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2512.78300233877        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1184
 NPARAMETR:  1.0033E+00  1.1168E+00  8.2519E-01  9.5698E-01  9.3560E-01  1.0259E+00  1.2554E+00  6.1735E-01  9.0374E-01  7.9017E-01
             1.0088E+00
 PARAMETER:  1.0334E-01  2.1043E-01 -9.2141E-02  5.6026E-02  3.3431E-02  1.2558E-01  3.2742E-01 -3.8232E-01 -1.2131E-03 -1.3550E-01
             1.0871E-01
 GRADIENT:  -6.2983E-01  2.9793E-01 -1.2768E+00  1.0188E+00  5.0906E-01 -1.0098E-01  1.3326E-01 -2.0772E+00 -2.7488E-01  2.7554E-01
             7.9332E-01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -2512.78423022817        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1285
 NPARAMETR:  1.0046E+00  1.1162E+00  8.2630E-01  9.5672E-01  9.3570E-01  1.0285E+00  1.2558E+00  6.1591E-01  9.0516E-01  7.8846E-01
             1.0094E+00
 PARAMETER:  1.0360E-01  2.1024E-01 -9.0544E-02  5.5624E-02  3.3181E-02  1.2793E-01  3.2739E-01 -3.8328E-01 -3.3649E-05 -1.3736E-01
             1.0898E-01
 GRADIENT:  -1.4961E+00  4.8324E-01  2.3781E-01 -3.3020E-01 -1.7063E+05 -1.1920E-01 -9.6266E-02  8.8862E+04 -5.4420E-02  5.3308E-02
            -3.1355E+05
 NUMSIGDIG:         1.9         2.6         2.5         2.7         2.3         2.6         2.8         2.3         2.3         2.5
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1285
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1928E-03 -8.2134E-03 -1.6227E-02  5.9793E-03 -1.5649E-02
 SE:             2.9898E-02  2.4727E-02  1.1324E-02  2.4918E-02  2.1924E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6818E-01  7.3977E-01  1.5188E-01  8.1036E-01  4.7538E-01

 ETASHRINKSD(%)  1.0000E-10  1.7162E+01  6.2063E+01  1.6523E+01  2.6552E+01
 ETASHRINKVR(%)  1.0000E-10  3.1378E+01  8.5608E+01  3.0315E+01  4.6053E+01
 EBVSHRINKSD(%)  3.0376E-01  1.7087E+01  6.5301E+01  1.7452E+01  2.6313E+01
 EBVSHRINKVR(%)  6.0659E-01  3.1254E+01  8.7960E+01  3.1858E+01  4.5702E+01
 RELATIVEINF(%)  9.9385E+01  2.0372E+01  5.9012E+00  2.2581E+01  1.1037E+01
 EPSSHRINKSD(%)  2.9549E+01
 EPSSHRINKVR(%)  5.0367E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2512.7842302281747     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1410.0579903825676     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.95
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.79
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2512.784       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.12E+00  8.27E-01  9.57E-01  9.35E-01  1.03E+00  1.26E+00  6.17E-01  9.05E-01  7.89E-01  1.01E+00
 


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
+        1.04E+03
 
 TH 2
+        6.57E+02  3.95E+02
 
 TH 3
+       -1.22E+07 -4.40E+07  5.06E+02
 
 TH 4
+       -7.07E+00  2.50E+02 -2.37E+02  1.00E+03
 
 TH 5
+       -3.58E+00  3.89E+07 -5.17E+02  4.49E+02  1.10E+03
 
 TH 6
+        2.49E+00 -6.41E-01  1.58E+00 -1.52E+00 -1.03E+00  1.87E+02
 
 TH 7
+        5.56E-01  2.16E+01 -2.51E+07  2.17E+07  1.11E+01 -3.90E-02  6.10E+01
 
 TH 8
+        6.49E+02  1.53E+07  4.36E+07 -9.11E+03  9.37E+02  1.05E+03  2.44E+02  1.52E+07
 
 TH 9
+       -9.07E+07 -1.24E+01 -1.14E+08  1.11E+08  3.51E+00  7.17E+07 -2.29E+07  3.98E+07  1.04E+08
 
 TH10
+       -1.75E-01 -1.56E+01 -2.27E+01 -2.03E+01 -5.65E+01 -5.99E+07  2.01E+01 -2.09E+01 -7.04E+01  9.95E+01
 
 TH11
+        7.47E+07 -1.64E+01 -4.05E+04  1.97E+04 -8.31E+07 -2.28E+03 -1.89E+07 -3.28E+07 -8.58E+07  7.17E+07  7.07E+07
 
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
 #CPUT: Total CPU Time in Seconds,       32.814
Stop Time:
Thu Sep 30 08:14:10 CDT 2021
