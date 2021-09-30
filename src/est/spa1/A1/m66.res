Wed Sep 29 22:42:57 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat66.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1667.89662886758        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9087E+02 -7.8287E+00  1.1269E+01 -1.2788E+01  6.9044E+01  3.9881E+01 -2.3526E+01  1.1680E+01 -1.1526E+01 -6.3279E+01
            -7.1571E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1798.70775650612        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0941E+00  1.0065E+00  1.0141E+00  1.1405E+00  9.6751E-01  1.0844E+00  1.3553E+00  7.0796E-01  1.0146E+00  1.6120E+00
             1.9304E+00
 PARAMETER:  1.8998E-01  1.0643E-01  1.1398E-01  2.3149E-01  6.6971E-02  1.8099E-01  4.0400E-01 -2.4537E-01  1.1449E-01  5.7746E-01
             7.5772E-01
 GRADIENT:   4.4296E+02  5.5838E+01 -2.5378E+01  1.3176E+02 -1.9652E+01  2.3448E+01  1.1675E+01  6.9590E+00  1.4488E+01  6.2388E+01
             1.3060E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1826.88607623128        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.0264E+00  9.1270E-01  8.4288E-01  1.1255E+00  8.9131E-01  1.0497E+00  1.3786E+00  3.8760E-02  9.2926E-01  1.0910E+00
             1.8506E+00
 PARAMETER:  1.2605E-01  8.6544E-03 -7.0935E-02  2.1824E-01 -1.5059E-02  1.4850E-01  4.2107E-01 -3.1504E+00  2.6628E-02  1.8712E-01
             7.1552E-01
 GRADIENT:   2.7856E+02  6.9447E+00 -4.3899E+01  8.7152E+01  6.3499E+01  3.5251E+01  6.3785E+00  3.0963E-02  6.6291E+00 -2.7423E+00
             1.0233E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1843.40967133047        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      292
 NPARAMETR:  9.5196E-01  8.6050E-01  1.0241E+00  1.1513E+00  9.3844E-01  9.6894E-01  1.4480E+00  7.1919E-02  9.1962E-01  1.2822E+00
             1.6189E+00
 PARAMETER:  5.0769E-02 -5.0240E-02  1.2382E-01  2.4091E-01  3.6463E-02  6.8449E-02  4.7018E-01 -2.5322E+00  1.6205E-02  3.4858E-01
             5.8175E-01
 GRADIENT:   3.4014E+00  1.2490E+01  6.6091E-01  2.0253E+00 -2.0825E-01  3.5517E+00  1.6667E+00  4.3631E-02 -6.6868E-01  2.1693E+00
            -6.1408E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1847.86754576049        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      471
 NPARAMETR:  9.4361E-01  3.9948E-01  9.5061E-01  1.4241E+00  7.2346E-01  9.4185E-01  2.0849E+00  4.7334E-02  8.1628E-01  1.1240E+00
             1.6333E+00
 PARAMETER:  4.1958E-02 -8.1759E-01  4.9347E-02  4.5356E-01 -2.2371E-01  4.0088E-02  8.3474E-01 -2.9505E+00 -1.0300E-01  2.1692E-01
             5.9061E-01
 GRADIENT:  -5.7809E+00  8.5498E+00 -1.4855E+00  3.6637E+01 -8.3469E+00 -5.6468E+00 -4.1542E+00  3.2558E-02 -7.4880E+00  1.2611E+00
             1.2563E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1849.96541799767        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  9.3720E-01  1.4485E-01  1.0480E+00  1.5698E+00  7.1049E-01  9.5124E-01  3.5911E+00  4.9280E-02  7.9666E-01  1.1657E+00
             1.6429E+00
 PARAMETER:  3.5146E-02 -1.8321E+00  1.4687E-01  5.5094E-01 -2.4180E-01  5.0011E-02  1.3784E+00 -2.9102E+00 -1.2733E-01  2.5332E-01
             5.9645E-01
 GRADIENT:  -7.7983E+00  1.9871E+00  5.0151E+00  2.6163E+01 -9.9627E+00 -1.2136E-01 -2.1619E+00  1.2736E-02  8.0821E-01  1.0559E+00
             2.7097E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1851.42712780661        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      828
 NPARAMETR:  9.3827E-01  4.1649E-02  1.0608E+00  1.6131E+00  6.9906E-01  9.4946E-01  6.2268E+00  2.8932E-02  7.7536E-01  1.1608E+00
             1.6411E+00
 PARAMETER:  3.6284E-02 -3.0785E+00  1.5899E-01  5.7814E-01 -2.5802E-01  4.8137E-02  1.9289E+00 -3.4428E+00 -1.5443E-01  2.4914E-01
             5.9538E-01
 GRADIENT:   1.4017E+00 -1.0586E-01  1.1717E+00  8.8578E-01 -1.4835E+00 -9.8076E-02 -6.7750E-01  2.7736E-03 -5.0307E-01  3.9064E-01
            -4.8183E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1852.13868146727        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      992
 NPARAMETR:  9.3300E-01  1.5494E-02  1.0308E+00  1.6281E+00  6.8073E-01  9.4921E-01  9.8379E+00  1.2900E-02  7.7392E-01  1.1307E+00
             1.6527E+00
 PARAMETER:  3.0646E-02 -4.0673E+00  1.3032E-01  5.8739E-01 -2.8459E-01  4.7879E-02  2.3862E+00 -4.2505E+00 -1.5629E-01  2.2281E-01
             6.0240E-01
 GRADIENT:  -1.0874E+01 -5.6772E-01 -2.6372E+00  1.5387E+01  3.5704E-01 -9.2708E-02  4.4475E+00  7.0948E-04  6.1107E-01  1.5785E-02
             5.4636E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1852.17340398668        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1056
 NPARAMETR:  9.3460E-01  1.5732E-02  1.0314E+00  1.6275E+00  6.8084E-01  9.4935E-01  9.8100E+00  1.1511E-02  7.7382E-01  1.1311E+00
             1.6500E+00
 PARAMETER:  3.2270E-02 -4.0517E+00  1.3084E-01  5.8701E-01 -2.8443E-01  4.8020E-02  2.3836E+00 -4.3212E+00 -1.5639E-01  2.2317E-01
             6.0088E-01
 GRADIENT:  -6.6977E+00  1.8684E+01 -2.0813E+00 -1.1503E+02  1.3233E-01 -3.9341E-02  3.0242E+01  5.7502E-04  3.5263E-01 -1.5943E-01
             2.6430E+00
 NUMSIGDIG:         1.3         2.3         1.4         2.4         3.4         3.0         2.3         0.2         1.9         2.5
                    2.0

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1056
 NO. OF SIG. DIGITS IN FINAL EST.:  0.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.1986E-03  6.0111E-03 -1.3514E-04 -1.3566E-02 -2.1388E-02
 SE:             2.9628E-02  5.5195E-03  1.9716E-04  2.8495E-02  2.4235E-02
 N:                     100         100         100         100         100

 P VAL.:         8.8731E-01  2.7613E-01  4.9307E-01  6.3402E-01  3.7748E-01

 ETASHRINKSD(%)  7.4082E-01  8.1509E+01  9.9339E+01  4.5389E+00  1.8811E+01
 ETASHRINKVR(%)  1.4761E+00  9.6581E+01  9.9996E+01  8.8717E+00  3.4083E+01
 EBVSHRINKSD(%)  9.6690E-01  8.7308E+01  9.9324E+01  4.0464E+00  1.6446E+01
 EBVSHRINKVR(%)  1.9245E+00  9.8389E+01  9.9995E+01  7.9291E+00  3.0187E+01
 RELATIVEINF(%)  9.7842E+01  1.0136E+00  5.6081E-04  4.9199E+01  8.8045E+00
 EPSSHRINKSD(%)  2.9897E+01
 EPSSHRINKVR(%)  5.0856E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1852.1734039866831     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -933.23487078201038     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1852.173       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.35E-01  1.57E-02  1.03E+00  1.63E+00  6.81E-01  9.49E-01  9.81E+00  1.20E-02  7.74E-01  1.13E+00  1.65E+00
 


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
+        1.38E+03
 
 TH 2
+       -1.22E+02  8.98E+05
 
 TH 3
+       -1.13E+01  4.67E+02  3.40E+02
 
 TH 4
+       -6.66E+00  1.43E+03 -1.12E+02  2.77E+03
 
 TH 5
+        1.90E+01 -1.38E+03 -5.99E+02  4.48E+01  1.34E+03
 
 TH 6
+        9.81E+00 -3.65E+00  2.04E-01 -4.42E+00 -7.37E-01  2.13E+02
 
 TH 7
+       -2.03E-01  2.46E+03  1.08E+00  2.96E+00 -3.05E+00  4.95E-03  6.75E+00
 
 TH 8
+        3.23E-01  2.31E+00 -2.44E-01 -2.68E-01  1.18E-01 -2.19E-01  8.17E-03  4.57E+00
 
 TH 9
+        5.08E+00 -3.40E+02  1.53E+01  1.92E+00 -6.38E+00  2.81E-01 -9.10E-01  4.06E-01  2.83E+02
 
 TH10
+        6.16E-01  6.49E+02 -1.22E+01 -4.70E+01 -7.23E+01  3.95E-01  1.78E+00  2.51E-01  2.95E+00  9.21E+01
 
 TH11
+       -1.36E+01  2.66E+02 -7.46E+00  2.06E+03  3.91E+00  1.68E+00  7.39E-01  2.20E-01  1.21E+01  1.87E+01  1.66E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.196
Stop Time:
Wed Sep 29 22:43:24 CDT 2021
