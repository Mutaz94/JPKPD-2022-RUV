Sat Sep 18 09:53:53 CDT 2021
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
$DATA ../../../../data/spa/A2/dat49.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -567.639257537473        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7517E+02 -2.2682E+01  8.2285E+01 -7.4049E+01  1.0022E+02  1.1269E+00 -7.2616E+01 -2.2830E+02 -1.3598E+02 -6.0245E+01
            -1.6107E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1362.44545621507        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.5324E-01  1.0351E+00  1.0643E+00  1.0458E+00  1.0420E+00  9.9268E-01  1.6693E+00  7.3743E-01  1.2785E+00  1.1799E+00
             1.9119E+00
 PARAMETER:  5.2116E-02  1.3453E-01  1.6228E-01  1.4477E-01  1.4117E-01  9.2658E-02  6.1241E-01 -2.0458E-01  3.4572E-01  2.6543E-01
             7.4808E-01
 GRADIENT:   2.2580E+01  1.5656E+01 -3.1787E+00  2.5386E+01  2.2123E+01  2.2406E+00  1.0868E+01  7.8561E-01  1.3757E+01 -4.8034E+00
            -6.8691E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1369.13895917808        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.4987E-01  8.1131E-01  9.7070E-01  1.1783E+00  8.9718E-01  9.7753E-01  1.6914E+00  3.3535E-01  1.1660E+00  1.0780E+00
             2.1459E+00
 PARAMETER:  4.8570E-02 -1.0911E-01  7.0261E-02  2.6403E-01 -8.5003E-03  7.7270E-02  6.2555E-01 -9.9259E-01  2.5360E-01  1.7512E-01
             8.6354E-01
 GRADIENT:   6.5137E+00  2.6179E+00 -1.7234E+01  3.4956E+01  2.6693E+01 -2.3390E+00 -2.8673E+00  4.4403E-01  8.9705E+00  6.7317E-01
            -1.8940E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1371.02216098533        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  9.4643E-01  5.9513E-01  8.1659E-01  1.2611E+00  7.0736E-01  9.8293E-01  2.1187E+00  2.6933E-01  1.0276E+00  8.8769E-01
             2.1467E+00
 PARAMETER:  4.4938E-02 -4.1898E-01 -1.0262E-01  3.3200E-01 -2.4621E-01  8.2787E-02  8.5082E-01 -1.2118E+00  1.2718E-01 -1.9133E-02
             8.6392E-01
 GRADIENT:  -1.3528E-01  3.1842E+00  3.8026E-01  1.2632E+01 -2.6548E+00 -8.1192E-02 -1.2451E+00  4.9438E-01 -1.0028E+00 -1.6634E+00
             2.1725E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1371.24337235993        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      308
 NPARAMETR:  9.4506E-01  5.6490E-01  8.7139E-01  1.2734E+00  7.3177E-01  9.8177E-01  2.1915E+00  2.4249E-01  1.0347E+00  9.4854E-01
             2.1445E+00
 PARAMETER:  4.3497E-02 -4.7111E-01 -3.7660E-02  3.4170E-01 -2.1228E-01  8.1606E-02  8.8458E-01 -1.3168E+00  1.3407E-01  4.7166E-02
             8.6289E-01
 GRADIENT:  -5.1051E-01 -6.3203E-01 -8.6636E-01 -8.6931E-01  9.6678E-01 -5.2656E-03  4.9691E-01  3.4632E-01  5.8242E-01  4.6831E-01
             2.3462E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1371.24681533883        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      403
 NPARAMETR:  9.4533E-01  5.6583E-01  8.6547E-01  1.2727E+00  7.2770E-01  9.8198E-01  2.1746E+00  2.3125E-01  1.0321E+00  9.3793E-01
             2.1451E+00
 PARAMETER:  4.3777E-02 -4.6946E-01 -4.4479E-02  3.4111E-01 -2.1787E-01  8.1815E-02  8.7684E-01 -1.3643E+00  1.3162E-01  3.5925E-02
             8.6321E-01
 GRADIENT:  -8.7539E+00 -1.8214E+00  2.3861E-01 -9.1806E+00 -9.6398E-01 -1.1091E+00 -1.3433E+00  2.9962E-01 -5.9510E-01 -3.4337E-01
             1.2450E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1372.31115511133        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      582
 NPARAMETR:  9.4985E-01  3.4019E-01  1.0233E+00  1.4470E+00  7.4129E-01  9.8260E-01  2.6653E+00  6.5641E-02  1.0114E+00  1.0326E+00
             2.1729E+00
 PARAMETER:  4.8544E-02 -9.7824E-01  1.2304E-01  4.6952E-01 -1.9937E-01  8.2448E-02  1.0803E+00 -2.6236E+00  1.1135E-01  1.3210E-01
             8.7605E-01
 GRADIENT:   1.1123E+01  4.4891E+00 -1.2526E+00  1.6483E+01  1.6297E+00  9.0729E-01 -3.8960E-01  7.3218E-03 -2.2807E+00 -1.9853E+00
             3.1316E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1373.82118968227        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      757
 NPARAMETR:  9.3951E-01  1.2065E-01  9.2731E-01  1.5462E+00  6.3793E-01  9.7768E-01  4.5171E+00  1.0000E-02  9.7475E-01  1.0010E+00
             2.1495E+00
 PARAMETER:  3.7600E-02 -2.0148E+00  2.4535E-02  5.3577E-01 -3.4953E-01  7.7427E-02  1.6079E+00 -4.9902E+00  7.4429E-02  1.0097E-01
             8.6523E-01
 GRADIENT:  -3.5383E+00  1.3663E-01  7.8847E+00  1.1907E+01 -1.2875E+01 -4.6660E-03 -1.9854E+00  0.0000E+00  8.7308E-01 -3.9540E-01
             1.3688E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1373.94336670051        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      932
 NPARAMETR:  9.4025E-01  1.0212E-01  9.2826E-01  1.5482E+00  6.4122E-01  9.7747E-01  4.9196E+00  1.0000E-02  9.7137E-01  1.0089E+00
             2.1462E+00
 PARAMETER:  3.8390E-02 -2.1816E+00  2.5552E-02  5.3708E-01 -3.4438E-01  7.7214E-02  1.6932E+00 -5.3875E+00  7.0951E-02  1.0884E-01
             8.6368E-01
 GRADIENT:  -9.3240E-02 -4.7245E-02  8.7068E-02  3.4175E-01 -2.2635E-02  3.7718E-02 -1.3877E-01  0.0000E+00  1.9282E-01  6.7168E-02
             4.1852E-02

0ITERATION NO.:   44    OBJECTIVE VALUE:  -1373.94381596574        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1059
 NPARAMETR:  9.4027E-01  1.0044E-01  9.2384E-01  1.5481E+00  6.3871E-01  9.7739E-01  4.9541E+00  1.0000E-02  9.7077E-01  1.0067E+00
             2.1457E+00
 PARAMETER:  3.8407E-02 -2.1982E+00  2.0779E-02  5.3705E-01 -3.4830E-01  7.7135E-02  1.7002E+00 -5.4478E+00  7.0331E-02  1.0669E-01
             8.6346E-01
 GRADIENT:   1.9456E-02 -1.6479E-03 -5.3631E-03 -4.1209E-04  1.0392E-02  3.4869E-04 -2.9056E-03  0.0000E+00  7.0829E-05  1.4451E-03
             3.1097E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1059
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.2781E-04  2.2704E-02 -8.4998E-05 -2.0967E-02 -1.6299E-02
 SE:             2.9285E-02  1.2062E-02  1.4232E-04  2.7007E-02  2.1081E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8834E-01  5.9806E-02  5.5034E-01  4.3754E-01  4.3944E-01

 ETASHRINKSD(%)  1.8912E+00  5.9590E+01  9.9523E+01  9.5216E+00  2.9374E+01
 ETASHRINKVR(%)  3.7466E+00  8.3670E+01  9.9998E+01  1.8137E+01  5.0120E+01
 EBVSHRINKSD(%)  1.8776E+00  7.0589E+01  9.9469E+01  8.2846E+00  2.5035E+01
 EBVSHRINKVR(%)  3.7199E+00  9.1350E+01  9.9997E+01  1.5883E+01  4.3802E+01
 RELATIVEINF(%)  9.5731E+01  3.3483E+00  2.2826E-04  3.7100E+01  4.4165E+00
 EPSSHRINKSD(%)  3.6198E+01
 EPSSHRINKVR(%)  5.9292E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1373.9438159657402     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -638.79298940200204     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.05
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1373.944       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.40E-01  1.00E-01  9.24E-01  1.55E+00  6.39E-01  9.77E-01  4.95E+00  1.00E-02  9.71E-01  1.01E+00  2.15E+00
 


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
+        1.27E+03
 
 TH 2
+       -1.21E+01  2.63E+03
 
 TH 3
+        6.09E+00 -1.05E+02  3.93E+02
 
 TH 4
+       -2.82E+01 -1.02E+02 -1.52E+01  4.69E+02
 
 TH 5
+        2.70E+01  1.52E+02 -7.59E+02 -1.05E+02  1.62E+03
 
 TH 6
+        4.98E+00 -1.52E+01  9.60E+00 -4.81E+00 -4.30E+00  1.91E+02
 
 TH 7
+        1.83E+00  1.18E+02 -1.05E+01 -1.58E+01  2.28E+01 -3.13E-01  5.76E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.13E+01 -2.19E+02  1.35E+01  2.52E+01 -4.19E+01 -5.20E+00 -8.90E+00  0.00E+00  1.78E+02
 
 TH10
+       -3.79E+00 -1.66E+02  2.16E+01  1.84E+01 -8.23E+01 -6.86E+00 -8.13E+00  0.00E+00  1.23E+01  7.43E+01
 
 TH11
+       -1.41E+01 -4.11E+01 -9.70E+00 -3.16E+00 -6.97E+00  2.93E+00 -1.72E+00  0.00E+00  9.31E+00  2.05E+01  5.87E+01
 
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
 #CPUT: Total CPU Time in Seconds,       20.565
Stop Time:
Sat Sep 18 09:54:15 CDT 2021
