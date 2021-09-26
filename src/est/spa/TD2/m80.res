Sat Sep 25 13:48:25 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat80.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1717.06900539931        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.3124E+01 -2.7312E+01  9.6319E+00 -4.8811E+01 -3.8320E+01  3.7023E+01  1.2463E+01  6.6540E+00  3.4355E+01  1.9376E+01
            -3.2604E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1725.56990650489        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0365E+00  1.0671E+00  1.0411E+00  9.9931E-01  1.0499E+00  8.5119E-01  9.1948E-01  9.5444E-01  8.0794E-01  8.9340E-01
             1.1367E+00
 PARAMETER:  1.3585E-01  1.6491E-01  1.4028E-01  9.9312E-02  1.4865E-01 -6.1121E-02  1.6057E-02  5.3370E-02 -1.1327E-01 -1.2726E-02
             2.2812E-01
 GRADIENT:   7.2456E+01  1.8228E+01  2.0208E+01 -4.1272E+00 -2.1591E+01 -2.3899E+01 -5.5660E-01 -1.8550E+00 -7.6615E+00  7.0099E-01
             1.1522E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1726.50724721505        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0294E+00  1.0680E+00  7.3964E-01  9.7597E-01  8.8014E-01  8.8291E-01  1.0308E+00  7.3166E-01  7.4620E-01  6.5415E-01
             1.1049E+00
 PARAMETER:  1.2897E-01  1.6575E-01 -2.0159E-01  7.5675E-02 -2.7676E-02 -2.4530E-02  1.3030E-01 -2.1244E-01 -1.9276E-01 -3.2442E-01
             1.9978E-01
 GRADIENT:   4.3551E+01  1.7455E+01  1.4495E+01 -5.8633E+00 -2.0620E+01 -9.3024E+00  2.1335E+00  1.0834E+00 -1.2842E+01 -4.1695E+00
             2.4417E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1728.04773590061        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0146E+00  1.0986E+00  6.6866E-01  9.4694E-01  8.6692E-01  9.0478E-01  9.5669E-01  4.5064E-01  8.2809E-01  7.1881E-01
             1.0869E+00
 PARAMETER:  1.1448E-01  1.9406E-01 -3.0249E-01  4.5478E-02 -4.2814E-02 -6.7655E-05  5.5719E-02 -6.9709E-01 -8.8639E-02 -2.3016E-01
             1.8331E-01
 GRADIENT:  -6.5290E-01 -1.2382E+00 -1.8142E+00 -6.7422E-02  7.6656E-01 -3.0725E-01 -3.2591E-01  9.4904E-01  5.8776E-02  3.1883E-01
             7.8364E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1728.12708098254        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  1.0138E+00  1.0772E+00  6.4348E-01  9.5566E-01  8.4249E-01  9.0738E-01  9.8879E-01  2.5443E-01  8.2030E-01  7.2019E-01
             1.0847E+00
 PARAMETER:  1.1371E-01  1.7433E-01 -3.4086E-01  5.4650E-02 -7.1392E-02  2.8112E-03  8.8727E-02 -1.2687E+00 -9.8082E-02 -2.2825E-01
             1.8135E-01
 GRADIENT:  -3.1993E+00 -3.9026E+00 -3.9642E+00  9.2163E-01  4.9472E+00  7.0100E-01  1.4021E+00  3.0498E-01  1.2471E+00  1.3200E+00
             8.1182E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1728.26662240952        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  1.0156E+00  1.1068E+00  5.9463E-01  9.3195E-01  8.2054E-01  9.0583E-01  9.6479E-01  8.7750E-02  8.1973E-01  6.7808E-01
             1.0837E+00
 PARAMETER:  1.1544E-01  2.0146E-01 -4.1981E-01  2.9521E-02 -9.7792E-02  1.0999E-03  6.4154E-02 -2.3333E+00 -9.8785E-02 -2.8849E-01
             1.8040E-01
 GRADIENT:   7.7528E-01  1.0260E+00  1.4465E+00 -4.6273E-01 -2.1913E+00 -1.4799E-01 -6.2293E-01  3.5223E-02 -4.5679E-01 -3.9454E-01
            -2.6456E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1728.71460008957        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.0246E+00  1.1422E+00  6.0666E-01  9.1562E-01  8.4988E-01  9.1512E-01  9.3872E-01  1.1818E-02  8.4324E-01  7.0737E-01
             1.0872E+00
 PARAMETER:  1.2431E-01  2.3294E-01 -3.9979E-01  1.1841E-02 -6.2656E-02  1.1296E-02  3.6765E-02 -4.3381E+00 -7.0507E-02 -2.4620E-01
             1.8365E-01
 GRADIENT:  -1.7220E+01 -5.6585E+00 -1.4850E+00 -1.2506E+00  1.7278E+00  5.7262E-01 -4.4434E-01  4.2439E-04  3.8848E-01 -2.7073E-01
             5.4818E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1729.00469960551        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      674
 NPARAMETR:  1.0308E+00  1.3269E+00  5.5852E-01  8.0490E-01  9.2019E-01  9.1351E-01  8.3421E-01  1.0000E-02  9.2987E-01  7.3168E-01
             1.0892E+00
 PARAMETER:  1.3031E-01  3.8286E-01 -4.8246E-01 -1.1704E-01  1.6827E-02  9.5378E-03 -8.1267E-02 -1.0696E+01  2.7294E-02 -2.1241E-01
             1.8541E-01
 GRADIENT:  -2.3426E+00  1.8337E+00  5.5558E-01  1.1804E+00 -9.4955E-01 -3.2210E-01  6.2824E-02  0.0000E+00  5.1990E-02 -5.1519E-02
            -6.5859E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1729.01167292624        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  1.0317E+00  1.3571E+00  5.4377E-01  7.8450E-01  9.2867E-01  9.1436E-01  8.1902E-01  1.0000E-02  9.4541E-01  7.3059E-01
             1.0895E+00
 PARAMETER:  1.3124E-01  4.0532E-01 -5.0923E-01 -1.4271E-01  2.6000E-02  1.0467E-02 -9.9641E-02 -1.1963E+01  4.3868E-02 -2.1391E-01
             1.8575E-01
 GRADIENT:   1.3275E-02 -1.8463E-02  3.5302E-03 -1.5689E-02  3.3007E-03 -2.0457E-03 -4.0053E-03  0.0000E+00  5.5344E-04 -1.2995E-03
            -3.6797E-03

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1729.01167292624        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:      871
 NPARAMETR:  1.0317E+00  1.3571E+00  5.4377E-01  7.8450E-01  9.2867E-01  9.1436E-01  8.1902E-01  1.0000E-02  9.4541E-01  7.3059E-01
             1.0895E+00
 PARAMETER:  1.3124E-01  4.0532E-01 -5.0923E-01 -1.4271E-01  2.6000E-02  1.0467E-02 -9.9641E-02 -1.1963E+01  4.3868E-02 -2.1391E-01
             1.8575E-01
 GRADIENT:   1.3275E-02 -1.8463E-02  3.5302E-03 -1.5689E-02  3.3007E-03 -2.0457E-03 -4.0053E-03  0.0000E+00  5.5344E-04 -1.2995E-03
            -3.6797E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      871
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1269E-05 -1.3936E-02 -3.3946E-04  8.5909E-03 -2.1839E-02
 SE:             2.9809E-02  2.3249E-02  1.5516E-04  2.3770E-02  2.0848E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9916E-01  5.4889E-01  2.8685E-02  7.1778E-01  2.9486E-01

 ETASHRINKSD(%)  1.3755E-01  2.2112E+01  9.9480E+01  2.0369E+01  3.0155E+01
 ETASHRINKVR(%)  2.7491E-01  3.9334E+01  9.9997E+01  3.6589E+01  5.1217E+01
 EBVSHRINKSD(%)  5.7365E-01  2.2002E+01  9.9516E+01  2.0844E+01  3.0114E+01
 EBVSHRINKVR(%)  1.1440E+00  3.9163E+01  9.9998E+01  3.7343E+01  5.1160E+01
 RELATIVEINF(%)  9.8786E+01  2.6548E+00  1.5215E-04  2.9612E+00  4.2774E+00
 EPSSHRINKSD(%)  4.2732E+01
 EPSSHRINKVR(%)  6.7204E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1729.0116729262390     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -993.86084636250087     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.46
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1729.012       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.36E+00  5.44E-01  7.84E-01  9.29E-01  9.14E-01  8.19E-01  1.00E-02  9.45E-01  7.31E-01  1.09E+00
 


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
+        1.23E+03
 
 TH 2
+       -1.04E+01  5.63E+02
 
 TH 3
+        1.33E+01  3.34E+02  9.32E+02
 
 TH 4
+       -2.48E+01  3.98E+02 -5.54E+02  1.28E+03
 
 TH 5
+       -5.26E+00 -4.45E+02 -9.13E+02  4.73E+02  1.18E+03
 
 TH 6
+       -1.43E+00 -1.78E+00  3.72E+00 -5.39E+00  1.90E-01  2.34E+02
 
 TH 7
+        8.04E-01  2.45E+01 -2.81E+01 -1.16E+01 -3.36E+00 -1.43E+00  1.15E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.47E+00 -2.45E+01 -4.58E+01  5.65E+01 -5.21E+00  1.05E+00  2.13E+01  0.00E+00  9.68E+01
 
 TH10
+       -2.20E-01 -1.26E+01 -5.75E+01 -2.25E+01 -6.86E+01 -1.50E-01  2.99E+01  0.00E+00  1.32E+01  9.46E+01
 
 TH11
+       -7.67E+00 -1.95E+01 -3.43E+01 -3.14E+00 -7.19E+00  2.52E+00  9.14E+00  0.00E+00  1.34E+01  2.55E+01  1.82E+02
 
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
 #CPUT: Total CPU Time in Seconds,       14.052
Stop Time:
Sat Sep 25 13:48:42 CDT 2021
