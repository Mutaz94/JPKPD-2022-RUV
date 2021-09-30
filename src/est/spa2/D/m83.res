Thu Sep 30 09:51:47 CDT 2021
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
$DATA ../../../../data/spa2/D/dat83.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:   36206.8176997487        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.9705E+02  5.9581E+02  1.0159E+02  5.9744E+02 -1.2501E+02 -2.1022E+03 -1.3762E+03 -1.4047E+02 -1.7776E+03 -2.5268E+02
            -7.0822E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -453.218871954753        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2928E+00  1.3247E+00  9.3933E-01  1.4618E+00  1.1050E+00  1.8438E+00  1.4353E+00  9.9786E-01  1.1620E+00  9.4678E-01
             1.4682E+01
 PARAMETER:  3.5680E-01  3.8115E-01  3.7409E-02  4.7968E-01  1.9986E-01  7.1182E-01  4.6134E-01  9.7861E-02  2.5015E-01  4.5310E-02
             2.7866E+00
 GRADIENT:  -1.8976E+01  2.1749E+01 -6.6314E+00  8.1910E+01  6.5891E+00  3.1408E+01 -3.3061E+01  1.9586E+00 -1.7626E+01  7.1749E+00
             3.5665E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -501.204328197298        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.2666E+00  1.7125E+00  2.1168E+00  1.2983E+00  1.9401E+00  2.1120E+00  4.6486E+00  7.9759E-01  1.1712E+00  3.0157E-01
             1.4116E+01
 PARAMETER:  3.3635E-01  6.3793E-01  8.4993E-01  3.6105E-01  7.6274E-01  8.4762E-01  1.6366E+00 -1.2616E-01  2.5804E-01 -1.0987E+00
             2.7473E+00
 GRADIENT:  -1.4011E+01  1.3684E+01 -1.3918E+01  1.7319E+01  1.3986E+01  3.3628E+01  2.9852E+01  6.4277E-01  1.4356E+01  5.6756E-01
             1.1592E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -507.314897748914        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      254
 NPARAMETR:  1.2554E+00  1.3087E+00  2.4368E+00  1.3905E+00  1.6032E+00  1.9433E+00  5.0297E+00  7.7104E-01  9.6713E-01  3.7636E-01
             1.4098E+01
 PARAMETER:  3.2749E-01  3.6902E-01  9.9069E-01  4.2965E-01  5.7200E-01  7.6441E-01  1.7154E+00 -1.6001E-01  6.6575E-02 -8.7721E-01
             2.7460E+00
 GRADIENT:  -2.7709E+01  7.1064E+00 -4.1840E+00  1.2355E+01 -4.8129E+00  9.9456E+00  1.1010E+01  6.1343E-01  8.1604E+00  1.0359E+00
             8.3006E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -511.769480940263        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      373
 NPARAMETR:  1.2543E+00  1.3085E+00  2.4700E+00  1.3923E+00  1.6133E+00  1.9308E+00  5.0445E+00  1.4763E-01  8.0858E-01  7.5716E-02
             1.3157E+01
 PARAMETER:  3.2661E-01  3.6887E-01  1.0042E+00  4.3092E-01  5.7828E-01  7.5795E-01  1.7183E+00 -1.8131E+00 -1.1247E-01 -2.4808E+00
             2.6770E+00
 GRADIENT:  -1.6627E+00  1.5188E+01 -1.1329E+00  5.0819E+01 -9.7410E+00  1.5049E+01  1.6376E+01  2.3299E-02  1.0866E-01  3.8604E-02
             2.9019E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -512.743069349013        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      556
 NPARAMETR:  1.2577E+00  1.3042E+00  2.4709E+00  1.3672E+00  1.6239E+00  1.8891E+00  4.9799E+00  1.0000E-02  6.9161E-01  1.0000E-02
             1.3333E+01
 PARAMETER:  3.2925E-01  3.6561E-01  1.0046E+00  4.1278E-01  5.8484E-01  7.3612E-01  1.7054E+00 -9.1410E+00 -2.6873E-01 -1.0093E+01
             2.6902E+00
 GRADIENT:  -1.5181E+01  1.1703E+01 -3.4177E-01  4.4850E+01 -1.0105E+01  2.9391E+00  4.5878E-01  0.0000E+00 -1.1188E+00  0.0000E+00
             9.2621E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -512.883393686769        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      735
 NPARAMETR:  1.2570E+00  1.3049E+00  2.4648E+00  1.3656E+00  1.6228E+00  1.8808E+00  4.9873E+00  1.0000E-02  7.5845E-01  1.0000E-02
             1.3155E+01
 PARAMETER:  3.2871E-01  3.6612E-01  1.0021E+00  4.1161E-01  5.8417E-01  7.3168E-01  1.7069E+00 -9.9315E+00 -1.7648E-01 -1.0912E+01
             2.6768E+00
 GRADIENT:   5.7906E-02  1.3785E+01 -7.9290E-01  4.6635E+01 -8.9071E+00  8.3731E+00  1.6070E+01  0.0000E+00  2.1860E-01  0.0000E+00
             3.0093E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -513.216261782814        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      917
 NPARAMETR:  1.2578E+00  1.3036E+00  2.4578E+00  1.3532E+00  1.6271E+00  1.8723E+00  4.9806E+00  1.0000E-02  6.8422E-01  1.0000E-02
             1.3325E+01
 PARAMETER:  3.2940E-01  3.6516E-01  9.9928E-01  4.0248E-01  5.8679E-01  7.2714E-01  1.7055E+00 -9.9315E+00 -2.7947E-01 -1.0912E+01
             2.6896E+00
 GRADIENT:  -1.4639E+01  1.0904E+01 -4.9580E-01  3.9924E+01 -9.1768E+00  5.4187E-01  2.0347E+00  0.0000E+00 -4.4094E-01  0.0000E+00
             1.1535E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -515.297811931520        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:     1020
 NPARAMETR:  1.2552E+00  1.3066E+00  5.4959E+00  1.3561E+00  2.1511E+00  1.8631E+00  5.0203E+00  1.0000E-02  7.8155E-01  1.0000E-02
             1.2954E+01
 PARAMETER:  3.2728E-01  3.6741E-01  1.8040E+00  4.0462E-01  8.6598E-01  7.2226E-01  1.7135E+00 -9.9315E+00 -1.4647E-01 -1.0912E+01
             2.6614E+00
 GRADIENT:   3.8815E+00  1.1358E+01 -2.9246E-01  2.9533E+01 -2.0951E+00  4.7016E+00  1.9454E+01  0.0000E+00  2.5843E-01  0.0000E+00
             2.5338E+01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -515.314128028685        NO. OF FUNC. EVALS.:  82
 CUMULATIVE NO. OF FUNC. EVALS.:     1102
 NPARAMETR:  1.2550E+00  1.3068E+00  5.6674E+00  1.3564E+00  2.1610E+00  1.8624E+00  5.0241E+00  1.0000E-02  7.8045E-01  1.0000E-02
             1.2932E+01
 PARAMETER:  3.2727E-01  3.6742E-01  1.8220E+00  4.0463E-01  8.7103E-01  7.2223E-01  1.7134E+00 -9.9315E+00 -1.4642E-01 -1.0912E+01
             2.6609E+00
 GRADIENT:   1.3743E+03 -1.2242E+03 -3.2299E-01 -1.0981E+03  5.1875E+02  6.2211E+02 -2.5438E+02  0.0000E+00  1.4482E-01  0.0000E+00
             1.6816E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1102
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4588E-02  3.3224E-03  9.6971E-06 -6.0936E-02  9.3597E-06
 SE:             2.8179E-02  2.3490E-02  6.6601E-06  1.0131E-02  4.3054E-05
 N:                     100         100         100         100         100

 P VAL.:         3.8291E-01  8.8752E-01  1.4539E-01  1.8058E-09  8.2790E-01

 ETASHRINKSD(%)  5.5953E+00  2.1306E+01  9.9978E+01  6.6061E+01  9.9856E+01
 ETASHRINKVR(%)  1.0878E+01  3.8072E+01  1.0000E+02  8.8481E+01  1.0000E+02
 EBVSHRINKSD(%)  7.4183E+00  1.5354E+01  9.9962E+01  7.0494E+01  9.9792E+01
 EBVSHRINKVR(%)  1.4286E+01  2.8350E+01  1.0000E+02  9.1294E+01  1.0000E+02
 RELATIVEINF(%)  8.4104E+01  3.2434E+01  1.4211E-06  3.7148E+00  4.2233E-05
 EPSSHRINKSD(%)  4.8809E+00
 EPSSHRINKVR(%)  9.5235E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -515.31412802868522     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       587.41211181692188     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -515.314       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.26E+00  1.31E+00  5.60E+00  1.36E+00  2.16E+00  1.86E+00  5.02E+00  1.00E-02  7.82E-01  1.00E-02  1.29E+01
 


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
+        6.71E+04
 
 TH 2
+        9.93E+01  4.92E+04
 
 TH 3
+       -4.82E-01  5.60E-01  1.09E+02
 
 TH 4
+       -5.54E+01  6.44E+01  2.95E-01  3.78E+04
 
 TH 5
+        1.79E+01 -1.93E+01 -8.52E-01 -2.27E+01  3.21E+03
 
 TH 6
+       -2.14E+02 -1.75E+04 -1.25E-01 -9.04E+00  5.19E+00  6.24E+03
 
 TH 7
+        1.64E+01 -3.81E+00 -6.09E-03 -8.53E+00 -2.42E-01  8.64E+00  1.50E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.44E-01 -2.03E+00  9.72E+03 -3.60E+01  2.98E+00 -1.73E+00  3.10E+00  0.00E+00  1.62E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.07E+01 -8.83E-01 -6.65E-03 -9.41E+00  8.70E-01 -1.53E+00 -1.47E+00  0.00E+00  2.25E+00  0.00E+00  1.37E+01
 
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
 #CPUT: Total CPU Time in Seconds,       37.938
Stop Time:
Thu Sep 30 09:52:33 CDT 2021
