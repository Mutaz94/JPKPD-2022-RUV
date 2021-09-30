Wed Sep 29 17:38:05 CDT 2021
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
$DATA ../../../../data/spa/S2/dat69.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m69.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1659.48055427081        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4532E+02  5.2532E+01 -7.6554E+01  2.0161E+02  1.2924E+02  5.2600E+01 -5.7133E+00  8.0873E+00  6.1227E+00 -3.4941E+00
            -4.8016E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1671.19175112204        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  1.0546E+00  1.0293E+00  1.1086E+00  9.0828E-01  9.3808E-01  1.0073E+00  1.0304E+00  9.7489E-01  1.0036E+00  9.8108E-01
             1.1314E+00
 PARAMETER:  1.5317E-01  1.2887E-01  2.0314E-01  3.7982E-03  3.6078E-02  1.0724E-01  1.2993E-01  7.4574E-02  1.0363E-01  8.0894E-02
             2.2345E-01
 GRADIENT:   6.3047E+01 -8.7397E+00  1.4550E+01 -3.3076E+01 -5.6136E+01  3.8227E+00 -4.2360E+00  2.6422E+00  2.0544E-02  7.3620E+00
             9.3744E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1672.63168345488        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  1.0389E+00  1.2748E+00  1.2176E+00  7.8234E-01  1.0781E+00  9.8519E-01  9.4074E-01  1.1144E+00  1.0741E+00  1.0815E+00
             1.1125E+00
 PARAMETER:  1.3813E-01  3.4280E-01  2.9687E-01 -1.4547E-01  1.7523E-01  8.5082E-02  3.8906E-02  2.0828E-01  1.7151E-01  1.7835E-01
             2.0661E-01
 GRADIENT:   2.9545E+01  2.2604E+01  1.3807E+01  9.1191E-01 -2.7099E+01 -3.7575E+00 -1.5580E+00 -3.5044E+00 -9.1420E+00  6.0843E+00
            -2.9111E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1673.87134908812        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  1.0261E+00  1.2506E+00  9.6215E-01  7.7570E-01  1.0126E+00  1.0156E+00  9.3727E-01  8.7911E-01  1.1240E+00  9.7445E-01
             1.1094E+00
 PARAMETER:  1.2575E-01  3.2366E-01  6.1415E-02 -1.5399E-01  1.1248E-01  1.1552E-01  3.5217E-02 -2.8847E-02  2.1689E-01  7.4122E-02
             2.0384E-01
 GRADIENT:  -3.2509E-01 -9.1652E+00  6.7394E-02 -6.1446E+00 -6.1618E+00  8.0460E+00 -3.0950E+00  1.3409E+00  9.3943E-01  6.3167E-01
             6.0272E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1674.40580441437        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  1.0267E+00  1.4462E+00  7.3626E-01  6.5795E-01  1.0257E+00  9.8907E-01  8.9203E-01  4.9013E-01  1.2267E+00  9.5388E-01
             1.1076E+00
 PARAMETER:  1.2632E-01  4.6893E-01 -2.0618E-01 -3.1863E-01  1.2533E-01  8.9009E-02 -1.4250E-02 -6.1308E-01  3.0436E-01  5.2780E-02
             2.0224E-01
 GRADIENT:  -2.4260E+00  9.3120E+00 -3.0746E+00  1.1927E+01  4.2474E+00 -2.6753E+00 -6.7975E-01  4.9844E-01 -1.0281E+00 -8.8311E-01
            -6.3401E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1674.55678588990        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  1.0253E+00  1.6064E+00  6.3493E-01  5.5794E-01  1.0725E+00  9.9257E-01  8.2775E-01  2.4130E-01  1.3838E+00  9.8659E-01
             1.1086E+00
 PARAMETER:  1.2498E-01  5.7398E-01 -3.5425E-01 -4.8351E-01  1.6999E-01  9.2544E-02 -8.9040E-02 -1.3217E+00  4.2487E-01  8.6504E-02
             2.0307E-01
 GRADIENT:  -6.4176E+00  1.4250E+01 -4.8517E+00  1.7028E+01  8.1783E+00 -1.3658E+00 -1.3638E+00  1.3930E-01 -6.0070E-01  1.3707E-01
            -2.2341E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1674.65796409332        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  1.0254E+00  1.7417E+00  5.6131E-01  4.7069E-01  1.1193E+00  9.9650E-01  7.8855E-01  1.0074E-01  1.5554E+00  1.0204E+00
             1.1102E+00
 PARAMETER:  1.2513E-01  6.5487E-01 -4.7748E-01 -6.5355E-01  2.1270E-01  9.6494E-02 -1.3757E-01 -2.1953E+00  5.4171E-01  1.2019E-01
             2.0453E-01
 GRADIENT:  -6.2949E+00  1.4847E+01 -5.0684E+00  1.6340E+01  6.8802E+00  2.1324E-01 -9.0012E-01  2.6289E-02  2.2447E-01  1.0984E+00
             1.7349E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1674.93024962668        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1234             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0294E+00  1.8210E+00  5.3169E-01  3.9202E-01  1.1523E+00  9.9642E-01  7.6873E-01  2.7538E-02  1.7026E+00  1.0405E+00
             1.1114E+00
 PARAMETER:  1.2895E-01  6.9937E-01 -5.3169E-01 -8.3643E-01  2.4177E-01  9.6411E-02 -1.6302E-01 -3.4922E+00  6.3214E-01  1.3966E-01
             2.0565E-01
 GRADIENT:   5.4703E+02  8.7788E+02  4.5641E+00  9.8894E+01  3.1462E+00  4.1625E+01  1.3266E+01  1.8035E-03  1.6800E+01  1.0273E+00
             1.4176E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1675.06332162495        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1416             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0297E+00  1.8307E+00  5.2330E-01  3.9288E-01  1.1618E+00  9.9574E-01  7.6515E-01  1.0000E-02  1.7546E+00  1.0451E+00
             1.1118E+00
 PARAMETER:  1.2929E-01  7.0470E-01 -5.4760E-01 -8.3426E-01  2.5000E-01  9.5728E-02 -1.6768E-01 -4.6975E+00  6.6223E-01  1.4410E-01
             2.0602E-01
 GRADIENT:   5.4814E+02  8.9628E+02  2.4798E-01  1.0583E+02  1.5231E+01  4.1206E+01  1.3738E+01  0.0000E+00  2.2758E+01  1.3196E+00
             2.0080E+00

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1675.06977822629        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:     1516
 NPARAMETR:  1.0300E+00  1.8299E+00  5.2721E-01  3.9078E-01  1.1636E+00  9.9582E-01  7.6348E-01  1.0000E-02  1.7473E+00  1.0457E+00
             1.1116E+00
 PARAMETER:  1.2828E-01  7.0684E-01 -5.4416E-01 -8.3454E-01  2.4944E-01  9.5635E-02 -1.6936E-01 -4.6975E+00  6.5725E-01  1.4325E-01
             2.0579E-01
 GRADIENT:  -1.3244E+00  3.1886E+00 -2.0557E-01  7.9410E-01 -1.2682E+00 -3.2341E-02  4.7115E-02  0.0000E+00 -5.4694E-02 -1.2029E-01
            -9.8425E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1516
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.6504E-04 -3.3759E-02 -1.8877E-04  3.0951E-02 -3.6393E-02
 SE:             2.9807E-02  2.3706E-02  6.7906E-05  2.1159E-02  2.3075E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8755E-01  1.5443E-01  5.4378E-03  1.4352E-01  1.1476E-01

 ETASHRINKSD(%)  1.4118E-01  2.0582E+01  9.9773E+01  2.9115E+01  2.2696E+01
 ETASHRINKVR(%)  2.8215E-01  3.6928E+01  9.9999E+01  4.9754E+01  4.0241E+01
 EBVSHRINKSD(%)  5.3756E-01  1.8995E+01  9.9802E+01  3.3326E+01  2.0198E+01
 EBVSHRINKVR(%)  1.0722E+00  3.4382E+01  1.0000E+02  5.5546E+01  3.6316E+01
 RELATIVEINF(%)  9.8855E+01  5.4958E+00  6.9902E-05  3.5544E+00  1.9582E+01
 EPSSHRINKSD(%)  4.3388E+01
 EPSSHRINKVR(%)  6.7951E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1675.0697782262932     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -939.91895166255506     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.07
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.47
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1675.070       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.83E+00  5.25E-01  3.93E-01  1.16E+00  9.96E-01  7.64E-01  1.00E-02  1.75E+00  1.04E+00  1.11E+00
 


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
+        1.05E+03
 
 TH 2
+       -8.12E+00  3.93E+02
 
 TH 3
+        6.67E+00  8.58E+01  2.00E+02
 
 TH 4
+       -1.98E+01  3.97E+02 -2.15E+02  1.09E+03
 
 TH 5
+       -5.98E+00 -1.29E+02 -1.83E+02  2.28E+02  4.84E+02
 
 TH 6
+        1.21E+00 -1.62E+00  9.65E-01 -4.39E+00 -1.43E+00  1.97E+02
 
 TH 7
+        1.64E+00  9.49E-01 -2.96E+00 -2.02E+01 -1.37E+01 -6.41E-01  1.68E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.53E-01 -1.45E+01 -2.33E+01  5.90E+01  1.56E+00 -1.72E-01  1.37E+01  0.00E+00  2.35E+01
 
 TH10
+       -6.37E-01 -1.54E+01 -2.40E+01  5.67E+00 -5.54E+01  3.30E-01  6.25E+00  0.00E+00  3.77E+00  8.26E+01
 
 TH11
+       -7.54E+00 -1.51E+01 -1.24E+01  4.32E+00  1.71E+00  2.34E+00  7.52E+00  0.00E+00  2.72E+00  1.64E+01  1.80E+02
 
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
 #CPUT: Total CPU Time in Seconds,       26.604
Stop Time:
Wed Sep 29 17:38:34 CDT 2021
