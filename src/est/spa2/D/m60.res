Thu Sep 30 09:23:40 CDT 2021
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
$DATA ../../../../data/spa2/D/dat60.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   18808.8687315684        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.5599E+02  2.8822E+02  5.9090E+01  5.4318E+01  3.9714E+02 -2.1026E+03 -9.4938E+02 -3.2017E+02 -1.6983E+03 -8.8870E+02
            -3.6271E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -671.916635661794        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2239E+00  1.2670E+00  9.1330E-01  2.1863E+00  1.0731E+00  3.5493E+00  2.8485E+00  9.5516E-01  3.5198E+00  1.4044E+00
             1.1297E+01
 PARAMETER:  3.0206E-01  3.3665E-01  9.3111E-03  8.8221E-01  1.7056E-01  1.3668E+00  1.1468E+00  5.4127E-02  1.3584E+00  4.3962E-01
             2.5245E+00
 GRADIENT:  -1.1775E+01 -4.2263E+01 -4.9962E+01  7.8035E+01  2.7707E+01  1.1813E+02 -5.5121E+01  3.6751E+00  5.2510E+01  2.1024E+01
             2.3269E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -748.816021070820        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      177
 NPARAMETR:  1.4460E+00  2.2409E+00  1.0760E+01  2.8607E+00  7.2774E+00  2.8551E+00  6.8283E+00  5.7763E-01  6.0892E+00  5.1405E+00
             1.0077E+01
 PARAMETER:  4.6880E-01  9.0688E-01  2.4759E+00  1.1511E+00  2.0848E+00  1.1491E+00  2.0211E+00 -4.4882E-01  1.9065E+00  1.7371E+00
             2.4103E+00
 GRADIENT:   1.8441E+01  1.5763E+01 -7.1973E+00  3.1077E+01 -2.7832E+01  1.7845E+00  3.3823E+01  2.9433E-02  7.0558E+01  4.2681E+01
             1.2429E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -814.953126144800        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  1.0864E+00  1.2356E+00  3.1900E+01  1.1764E+00  7.5750E+00  3.2542E+00  4.7156E+00  2.3008E-01  1.9577E+00  2.4834E+00
             9.3013E+00
 PARAMETER:  1.8290E-01  3.1157E-01  3.5626E+00  2.6245E-01  2.1249E+00  1.2799E+00  1.6509E+00 -1.3693E+00  7.7178E-01  1.0096E+00
             2.3302E+00
 GRADIENT:  -4.1999E+01 -2.6770E+01  4.2229E-02 -2.1553E+01  2.8029E+00  3.1149E+01 -3.6267E+01  2.2981E-06  2.0351E+01  3.4460E-01
             1.3707E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -842.682157916718        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  1.3259E+00  1.3803E+00  8.6186E+01  1.0958E+00  1.9412E+01  2.7736E+00  5.5607E+00  5.0170E+00  8.4382E-01  3.7637E+00
             8.4223E+00
 PARAMETER:  3.8212E-01  4.2229E-01  4.5565E+00  1.9148E-01  3.0659E+00  1.1201E+00  1.8157E+00  1.7128E+00 -6.9819E-02  1.4254E+00
             2.2309E+00
 GRADIENT:   1.6379E+01 -3.6083E+00  4.4934E-04 -7.5012E+00  1.0720E-01 -5.9071E-01 -5.2983E-01 -3.0951E-05  5.2050E+00 -8.8259E-03
            -2.3532E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -845.427500653121        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  1.2503E+00  1.5988E+00  1.2870E+02  9.9386E-01  2.9135E+01  2.7667E+00  5.3078E+00  3.1409E+01  2.6085E-01  5.5177E+00
             8.5971E+00
 PARAMETER:  3.2337E-01  5.6927E-01  4.9575E+00  9.3840E-02  3.4720E+00  1.1176E+00  1.7692E+00  3.5471E+00 -1.2438E+00  1.8080E+00
             2.2514E+00
 GRADIENT:  -9.3100E-01 -1.2336E+00 -4.1652E-04 -1.5954E+00  9.2535E-02 -2.4937E-01 -1.3910E+00 -4.3362E-03  3.9139E-01 -2.2251E-02
            -3.5891E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -845.591363771037        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      890             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2585E+00  1.6289E+00  1.8340E+02  9.9092E-01  3.4045E+01  2.7705E+00  5.3376E+00  1.1809E+02  1.1754E-01  7.0670E+00
             8.6557E+00
 PARAMETER:  3.2995E-01  5.8789E-01  5.3117E+00  9.0883E-02  3.6277E+00  1.1190E+00  1.7748E+00  4.8714E+00 -2.0409E+00  2.0554E+00
             2.2582E+00
 GRADIENT:   2.2967E+01  7.0727E+00  1.9924E-02  2.3694E+00  6.2705E-02  4.4124E+01  6.3201E+01 -1.0882E-02  1.1008E-01  7.6348E-02
             2.8855E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -846.320258076526        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      968
 NPARAMETR:  1.2455E+00  1.6126E+00  1.9426E+02  9.8524E-01  4.5500E+00  2.7560E+00  5.2438E+00  1.1202E+02  1.6939E-02  2.5945E+00
             8.6114E+00
 PARAMETER:  3.1957E-01  5.7783E-01  5.3692E+00  8.5130E-02  1.6151E+00  1.1138E+00  1.7570E+00  4.8186E+00 -3.9782E+00  1.0534E+00
             2.2531E+00
 GRADIENT:   1.6981E+01  4.7423E+00 -5.1122E-02 -1.5850E+00 -9.8061E+00  4.2766E+01  6.1160E+01  3.7974E-02  4.2683E-03  3.2037E+01
             2.0332E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -847.245801558643        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1124
 NPARAMETR:  1.2474E+00  1.6171E+00  1.9467E+02  9.8573E-01  4.4163E+00  2.7564E+00  5.2390E+00  1.1180E+02  1.6075E-02  2.4955E+00
             8.6156E+00
 PARAMETER:  3.2107E-01  5.8066E-01  5.3713E+00  8.5625E-02  1.5853E+00  1.1139E+00  1.7561E+00  4.8167E+00 -4.0305E+00  1.0145E+00
             2.2536E+00
 GRADIENT:  -5.0219E+00 -2.3271E+00 -2.0743E-01 -3.1006E+00 -7.6360E+00 -1.7200E+00 -2.9753E+00  2.8609E-01  1.9450E-03  2.8047E+01
            -3.6754E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -847.308134258156        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1301
 NPARAMETR:  1.2718E+00  1.6187E+00  1.9635E+02  9.8680E-01  4.4320E+00  2.7518E+00  5.2263E+00  1.1095E+02  1.0000E-02  2.4972E+00
             8.6311E+00
 PARAMETER:  3.4064E-01  5.8120E-01  5.3761E+00  8.5714E-02  1.5877E+00  1.1130E+00  1.7549E+00  4.8125E+00 -4.8220E+00  1.0145E+00
             2.2539E+00
 GRADIENT:   8.7300E+02 -2.5808E+02 -3.5884E+01 -2.1945E+00 -1.2776E+02  3.9270E+02  1.6550E+02  3.1431E+01  0.0000E+00 -2.6485E+02
            -1.7613E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1301
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1279E-02  1.4345E-02  3.7773E-03 -8.4495E-04 -4.7467E-02
 SE:             2.9122E-02  2.6868E-02  3.9252E-04  1.7525E-04  5.2765E-03
 N:                     100         100         100         100         100

 P VAL.:         6.9854E-01  5.9341E-01  6.5015E-22  1.4269E-06  2.3781E-19

 ETASHRINKSD(%)  2.4387E+00  9.9879E+00  9.8685E+01  9.9413E+01  8.2323E+01
 ETASHRINKVR(%)  4.8179E+00  1.8978E+01  9.9983E+01  9.9997E+01  9.6875E+01
 EBVSHRINKSD(%)  2.9646E+00  6.7403E+00  9.7824E+01  9.9604E+01  6.8006E+01
 EBVSHRINKVR(%)  5.8413E+00  1.3026E+01  9.9953E+01  9.9998E+01  8.9764E+01
 RELATIVEINF(%)  9.0853E+01  5.0751E+01  1.0912E-02  7.3266E-04  2.5086E+00
 EPSSHRINKSD(%)  9.5301E+00
 EPSSHRINKVR(%)  1.8152E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -847.30813425815575     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       255.41810558745135     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    32.76
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.32
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -847.308       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.27E+00  1.62E+00  1.96E+02  9.86E-01  4.43E+00  2.75E+00  5.23E+00  1.11E+02  1.00E-02  2.50E+00  8.62E+00
 


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
+        3.97E+04
 
 TH 2
+        2.48E-01  8.42E+03
 
 TH 3
+       -7.60E-03  1.91E-02  6.72E-03
 
 TH 4
+        7.99E+01 -1.07E+01 -1.97E-02  3.14E+02
 
 TH 5
+       -8.70E-01  2.79E+00 -1.01E+00 -3.58E+00  1.50E+02
 
 TH 6
+        2.71E+03  9.82E-01  6.40E-03  1.22E+01 -1.01E+02  1.17E+03
 
 TH 7
+        7.86E+00 -2.01E+00 -1.18E-02 -1.83E+01 -1.65E+00 -1.26E+02  9.36E+01
 
 TH 8
+       -9.42E+00  4.32E+00  2.76E-06  6.50E-02  5.86E-01  1.32E-03  4.60E-01  2.67E-02
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        6.31E+00 -3.13E+00 -7.99E-04 -1.57E+01  1.19E+02  1.81E+00 -9.85E+00  4.46E-03  0.00E+00  1.17E+03
 
 TH11
+       -2.66E+02  1.19E+02  1.09E-01 -1.80E+01  1.66E+01 -5.43E+01  4.20E+01  2.07E-01  0.00E+00  4.28E+01  3.52E+01
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       46.160
Stop Time:
Thu Sep 30 09:24:28 CDT 2021
