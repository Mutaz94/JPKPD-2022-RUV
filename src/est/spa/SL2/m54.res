Sat Sep 18 12:19:51 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat54.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1691.10092649705        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -4.2182E+01  4.7230E+01 -3.6105E+01  1.0557E+02 -1.7926E+01 -2.6408E+01 -1.3410E+01  1.3853E+01 -1.3729E+01  2.1310E+01
            -2.1457E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1705.46785879272        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0452E+00  1.0889E+00  1.2778E+00  8.9655E-01  1.1482E+00  1.0839E+00  1.1039E+00  8.8570E-01  1.0838E+00  8.6576E-01
             1.1038E+00
 PARAMETER:  1.4420E-01  1.8513E-01  3.4511E-01 -9.2008E-03  2.3820E-01  1.8061E-01  1.9884E-01 -2.1377E-02  1.8052E-01 -4.4146E-02
             1.9878E-01
 GRADIENT:   6.7952E+01  1.2468E+01  5.0891E+00  6.5413E+00  1.4333E+00  1.5337E+01  3.5338E+00  2.3019E+00 -8.2965E-01 -1.1693E+01
             7.4387E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1706.86041911567        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0415E+00  1.0645E+00  1.4075E+00  9.1059E-01  1.1941E+00  1.0758E+00  8.7341E-01  5.4888E-01  1.2158E+00  1.0401E+00
             1.0717E+00
 PARAMETER:  1.4070E-01  1.6248E-01  4.4183E-01  6.3417E-03  2.7743E-01  1.7309E-01 -3.5352E-02 -4.9987E-01  2.9538E-01  1.3936E-01
             1.6928E-01
 GRADIENT:   6.5225E+01  1.3192E+00  7.3704E+00  9.0029E+00  3.9088E+00  1.2850E+01  2.1641E+00 -9.9080E-01  1.2278E+01 -5.8971E+00
            -1.1772E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1709.04104414632        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      318
 NPARAMETR:  1.0454E+00  1.3434E+00  1.0596E+00  7.3583E-01  1.2053E+00  1.0767E+00  8.5119E-01  3.0437E-01  1.3091E+00  1.0492E+00
             1.0633E+00
 PARAMETER:  1.4442E-01  3.9521E-01  1.5789E-01 -2.0675E-01  2.8669E-01  1.7390E-01 -6.1121E-02 -1.0895E+00  3.6934E-01  1.4803E-01
             1.6135E-01
 GRADIENT:   7.0301E-01  1.0665E+01 -2.6480E+00  1.1929E+01 -7.9013E+00  2.7032E-01 -6.8120E-01  5.4422E-01 -3.2681E+00  4.0096E+00
             6.9843E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1711.86225563074        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      495
 NPARAMETR:  1.0451E+00  1.7564E+00  9.2336E-01  4.5654E-01  1.3798E+00  1.0759E+00  6.6449E-01  8.0168E-02  1.9997E+00  1.1087E+00
             1.0738E+00
 PARAMETER:  1.4413E-01  6.6325E-01  2.0266E-02 -6.8407E-01  4.2191E-01  1.7316E-01 -3.0874E-01 -2.4236E+00  7.9298E-01  2.0319E-01
             1.7117E-01
 GRADIENT:  -2.3252E+00  1.2050E+01 -3.0623E+00  1.1181E+01  8.3841E+00 -3.0967E-01  5.7731E-01  3.1214E-02  1.8753E+00 -1.4013E+00
            -1.0116E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1712.51914556358        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      670
 NPARAMETR:  1.0470E+00  1.9164E+00  8.5914E-01  3.3197E-01  1.4247E+00  1.0765E+00  6.2544E-01  2.4235E-02  2.4532E+00  1.1323E+00
             1.0820E+00
 PARAMETER:  1.4596E-01  7.5046E-01 -5.1820E-02 -1.0027E+00  4.5394E-01  1.7371E-01 -3.6931E-01 -3.6200E+00  9.9740E-01  2.2429E-01
             1.7877E-01
 GRADIENT:   5.3556E-01  2.1175E+00  4.1071E-01  3.7967E-01 -1.0829E+00 -1.1369E-02  3.5793E-02  2.3919E-03  4.9181E-02  1.1796E-02
            -3.2702E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1712.52126060483        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      848             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0468E+00  1.9198E+00  8.5275E-01  3.2861E-01  1.4265E+00  1.0766E+00  6.2570E-01  1.5976E-02  2.4644E+00  1.1320E+00
             1.0826E+00
 PARAMETER:  1.4573E-01  7.5222E-01 -5.9293E-02 -1.0129E+00  4.5519E-01  1.7379E-01 -3.6888E-01 -4.0367E+00  1.0020E+00  2.2398E-01
             1.7939E-01
 GRADIENT:   7.0204E+01  1.2888E+02 -1.1734E-02  1.0812E+01  1.8910E+00  1.2510E+01  3.0467E+00  1.0645E-03  5.8241E+00  1.2218E-01
             1.5716E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1712.52161436174        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  1.0467E+00  1.9200E+00  8.5284E-01  3.2873E-01  1.4265E+00  1.0765E+00  6.2578E-01  1.0000E-02  2.4637E+00  1.1320E+00
             1.0826E+00
 PARAMETER:  1.4569E-01  7.5234E-01 -5.9179E-02 -1.0125E+00  4.5521E-01  1.7374E-01 -3.6876E-01 -4.5983E+00  1.0017E+00  2.2400E-01
             1.7941E-01
 GRADIENT:  -1.4744E-03 -6.5634E-02 -2.0637E-03 -2.0637E-02  2.4352E-02 -2.0843E-03  8.8273E-03  0.0000E+00 -7.2176E-03  3.0033E-03
             7.0742E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1712.52168759940        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1215
 NPARAMETR:  1.0467E+00  1.9192E+00  8.5319E-01  3.2928E-01  1.4262E+00  1.0765E+00  6.2527E-01  1.0000E-02  2.4612E+00  1.1321E+00
             1.0827E+00
 PARAMETER:  1.4565E-01  7.5193E-01 -5.8772E-02 -1.0109E+00  4.5502E-01  1.7369E-01 -3.6957E-01 -4.5983E+00  1.0006E+00  2.2405E-01
             1.7942E-01
 GRADIENT:  -8.3066E-02 -9.7149E-03 -3.3512E-02 -1.6634E-03  3.6992E-02 -1.9684E-02 -8.9515E-02  0.0000E+00 -4.6897E-02  1.0203E-02
             2.8454E-02

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1712.52173260232        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1272
 NPARAMETR:  1.0467E+00  1.9192E+00  8.5338E-01  3.2928E-01  1.4262E+00  1.0765E+00  6.2569E-01  1.0000E-02  2.4617E+00  1.1320E+00
             1.0826E+00
 PARAMETER:  1.4568E-01  7.5193E-01 -5.8556E-02 -1.0108E+00  4.5499E-01  1.7373E-01 -3.6890E-01 -4.5983E+00  1.0009E+00  2.2399E-01
             1.7937E-01
 GRADIENT:  -1.0368E-02 -4.9676E-02 -1.0980E-03  5.2375E-03 -3.8761E-02 -5.3088E-03 -1.6589E-04  0.0000E+00  1.6073E-02  1.0904E-02
             3.2265E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1272
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.6561E-04 -4.9164E-02 -9.8048E-05  3.4360E-02 -3.9697E-02
 SE:             2.9848E-02  1.9728E-02  5.4648E-05  2.2430E-02  2.3382E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9557E-01  1.2697E-02  7.2786E-02  1.2556E-01  8.9544E-02

 ETASHRINKSD(%)  5.8671E-03  3.3910E+01  9.9817E+01  2.4855E+01  2.1669E+01
 ETASHRINKVR(%)  1.1734E-02  5.6321E+01  1.0000E+02  4.3533E+01  3.8643E+01
 EBVSHRINKSD(%)  4.4062E-01  3.0073E+01  9.9803E+01  2.7464E+01  2.0289E+01
 EBVSHRINKVR(%)  8.7930E-01  5.1102E+01  1.0000E+02  4.7385E+01  3.6462E+01
 RELATIVEINF(%)  9.9064E+01  7.4593E+00  2.1007E-04  8.8363E+00  3.9466E+01
 EPSSHRINKSD(%)  4.1392E+01
 EPSSHRINKVR(%)  6.5652E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1712.5217326023183     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -977.37090603858007     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.32
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1712.522       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.92E+00  8.53E-01  3.29E-01  1.43E+00  1.08E+00  6.26E-01  1.00E-02  2.46E+00  1.13E+00  1.08E+00
 


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
+        8.67E+02
 
 TH 2
+       -9.18E+00  4.11E+02
 
 TH 3
+        3.23E-01  4.64E+01  5.27E+01
 
 TH 4
+       -1.15E+01  4.40E+02 -3.92E+01  9.29E+02
 
 TH 5
+       -4.27E+00 -8.43E+01 -4.67E+01  4.31E+01  2.77E+02
 
 TH 6
+        1.30E+00 -7.89E-01  5.41E+00 -3.08E+00 -1.20E+00  1.68E+02
 
 TH 7
+        3.20E+00 -4.40E+01  2.09E+01 -2.98E+01 -2.44E+01  1.09E-01  1.43E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.80E-01 -1.49E+01 -6.36E+00  4.82E+01  2.52E+00 -3.19E-01  1.50E+01  0.00E+00  1.41E+01
 
 TH10
+        4.07E-01 -3.93E+00 -8.31E+00 -6.19E+00 -4.25E+01  2.83E-01  1.19E+01  0.00E+00  9.29E-01  6.20E+01
 
 TH11
+       -5.68E+00 -2.45E+01 -1.16E+01 -2.10E+00  4.22E+00  1.49E+00  5.91E+00  0.00E+00  1.70E+00  2.68E+01  1.94E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.402
Stop Time:
Sat Sep 18 12:20:18 CDT 2021
