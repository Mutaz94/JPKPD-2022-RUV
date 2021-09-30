Wed Sep 29 20:47:48 CDT 2021
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
$DATA ../../../../data/spa1/B/dat4.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m4.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2085.44741580083        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9013E+02  6.5400E+00 -2.2323E+01  7.7969E+01  8.3684E+01  2.5033E+01  1.5783E+00 -2.9104E+00  3.2050E+00 -1.8467E+01
             2.7445E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2089.15779662563        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      178
 NPARAMETR:  9.8157E-01  1.0275E+00  1.0017E+00  1.0097E+00  9.2277E-01  1.0886E+00  1.0039E+00  1.0149E+00  1.0176E+00  1.0481E+00
             9.4743E-01
 PARAMETER:  8.1399E-02  1.2710E-01  1.0166E-01  1.0962E-01  1.9621E-02  1.8490E-01  1.0387E-01  1.1481E-01  1.1741E-01  1.4701E-01
             4.5998E-02
 GRADIENT:   9.1477E+00  3.2196E+01  2.4963E+01  1.5727E+01 -4.7667E+01  3.5861E+00  1.0969E+00 -3.5464E+00 -1.8488E+00 -5.4917E-01
            -1.3262E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2089.50152626590        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      359
 NPARAMETR:  9.8217E-01  1.1105E+00  1.0033E+00  9.5759E-01  9.6210E-01  1.0839E+00  9.2340E-01  1.0528E+00  1.0792E+00  1.1127E+00
             9.5195E-01
 PARAMETER:  8.2007E-02  2.0484E-01  1.0329E-01  5.6665E-02  6.1364E-02  1.8061E-01  2.0311E-02  1.5146E-01  1.7622E-01  2.0675E-01
             5.0762E-02
 GRADIENT:   9.4688E+00  3.3647E+01  2.2993E+01  1.3549E+01 -4.4900E+01  1.6401E+00  1.7834E+00 -3.1916E+00  4.4137E-01  3.8524E+00
            -8.5838E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2090.68388185408        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      540
 NPARAMETR:  9.7738E-01  1.1084E+00  1.0698E+00  9.4538E-01  1.0207E+00  1.0804E+00  7.3612E-01  1.2340E+00  1.1353E+00  1.1235E+00
             9.6167E-01
 PARAMETER:  7.7115E-02  2.0288E-01  1.6748E-01  4.3834E-02  1.2051E-01  1.7737E-01 -2.0636E-01  3.1026E-01  2.2688E-01  2.1642E-01
             6.0914E-02
 GRADIENT:   5.2625E-01  3.9090E+00  3.8327E+00  1.1025E+00 -4.8449E+00  4.0441E-01 -7.0041E-01 -1.0725E+00  1.7286E-01 -1.8165E+00
            -6.8444E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2090.71013933080        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  9.7738E-01  1.1083E+00  1.0698E+00  9.4514E-01  1.0207E+00  1.0791E+00  7.6790E-01  1.2468E+00  1.1260E+00  1.1303E+00
             9.6184E-01
 PARAMETER:  7.7115E-02  2.0282E-01  1.6748E-01  4.3578E-02  1.2051E-01  1.7616E-01 -1.6409E-01  3.2055E-01  2.1867E-01  2.2247E-01
             6.1093E-02
 GRADIENT:   5.6327E-01  2.4258E+00  3.0584E+00 -4.4353E-01 -5.7041E+00 -5.6020E-02  1.6130E-02  3.4240E-02  1.3967E-01  4.8568E-02
            -2.4477E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2090.73146383210        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      900
 NPARAMETR:  9.7758E-01  1.1072E+00  1.0622E+00  9.4566E-01  1.0230E+00  1.0798E+00  7.7245E-01  1.2310E+00  1.1243E+00  1.1326E+00
             9.6176E-01
 PARAMETER:  7.7326E-02  2.0180E-01  1.6034E-01  4.4130E-02  1.2272E-01  1.7678E-01 -1.5819E-01  3.0780E-01  2.1712E-01  2.2454E-01
             6.1006E-02
 GRADIENT:   9.8077E-01 -1.4159E-01  4.6336E-01  2.5081E-01  7.0542E-02  2.0402E-01 -3.4535E-02  3.3385E-02  1.4042E-01  9.0933E-02
             2.6199E-04

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2090.73618581103        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1080
 NPARAMETR:  9.7556E-01  1.1071E+00  1.0363E+00  9.4476E-01  1.0128E+00  1.0771E+00  8.0524E-01  1.1863E+00  1.1131E+00  1.1220E+00
             9.6163E-01
 PARAMETER:  7.5254E-02  2.0178E-01  1.3570E-01  4.3175E-02  1.1275E-01  1.7431E-01 -1.1662E-01  2.7085E-01  2.0718E-01  2.1511E-01
             6.0875E-02
 GRADIENT:  -3.1399E+00 -8.8370E-01  3.0122E-01 -2.5417E-01  4.1139E-01 -8.3784E-01  6.5562E-03 -2.2988E-02 -1.6299E-01  1.7985E-01
            -2.4496E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2090.73979923507        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  9.7602E-01  1.1060E+00  1.0347E+00  9.4420E-01  1.0124E+00  1.0811E+00  8.0710E-01  1.1811E+00  1.1135E+00  1.1220E+00
             9.6165E-01
 PARAMETER:  7.5724E-02  2.0077E-01  1.3412E-01  4.2588E-02  1.1236E-01  1.7799E-01 -1.1431E-01  2.6641E-01  2.0748E-01  2.1515E-01
             6.0896E-02
 GRADIENT:  -2.1941E+00 -2.4995E+00  1.5441E-01 -1.8401E+00  1.0223E+00  6.6298E-01  7.3385E-02 -7.6662E-02  6.4723E-02  2.1861E-01
             3.3218E-02

0ITERATION NO.:   39    OBJECTIVE VALUE:  -2090.74667651630        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1379
 NPARAMETR:  9.7805E-01  1.1077E+00  1.0349E+00  9.4434E-01  1.0119E+00  1.0814E+00  8.0625E-01  1.1813E+00  1.1140E+00  1.1221E+00
             9.6170E-01
 PARAMETER:  7.7637E-02  2.0261E-01  1.3390E-01  4.3574E-02  1.1173E-01  1.7801E-01 -1.1653E-01  2.6747E-01  2.0762E-01  2.1455E-01
             6.0872E-02
 GRADIENT:  -5.1932E-01  5.6711E-01 -8.1735E+04  9.2884E-01 -1.9833E-01 -1.3115E-01 -2.5466E-02  4.0773E+04 -1.3025E-01 -5.1032E+04
            -9.3392E-02
 NUMSIGDIG:         2.6         2.5         2.3         1.9         2.8         2.7         1.8         2.3         2.6         2.3
                    2.9

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1379
 NO. OF SIG. DIGITS IN FINAL EST.:  1.8
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.7308E-04 -2.9551E-02 -3.3552E-02  1.0985E-02 -3.6402E-02
 SE:             2.9920E-02  1.5526E-02  1.4632E-02  2.6449E-02  2.2828E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7939E-01  5.6994E-02  2.1841E-02  6.7789E-01  1.1080E-01

 ETASHRINKSD(%)  1.0000E-10  4.7986E+01  5.0982E+01  1.1394E+01  2.3524E+01
 ETASHRINKVR(%)  1.0000E-10  7.2945E+01  7.5973E+01  2.1489E+01  4.1515E+01
 EBVSHRINKSD(%)  2.7643E-01  4.8150E+01  5.4825E+01  1.1793E+01  2.0280E+01
 EBVSHRINKVR(%)  5.5209E-01  7.3116E+01  7.9592E+01  2.2194E+01  3.6446E+01
 RELATIVEINF(%)  9.9000E+01  1.3733E+00  3.5364E+00  4.7549E+00  1.6590E+01
 EPSSHRINKSD(%)  3.4844E+01
 EPSSHRINKVR(%)  5.7547E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2090.7466765163040     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1171.8081433116313     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.02
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2090.747       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.78E-01  1.11E+00  1.03E+00  9.45E-01  1.01E+00  1.08E+00  8.05E-01  1.18E+00  1.11E+00  1.12E+00  9.62E-01
 


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
+        9.89E+02
 
 TH 2
+       -8.53E+00  4.00E+02
 
 TH 3
+       -1.47E+02  3.94E+03  1.43E+07
 
 TH 4
+       -2.26E+00  1.29E+07 -2.09E+07  3.06E+07
 
 TH 5
+        3.03E+00 -2.02E+02 -4.93E+02  6.18E+01  5.61E+02
 
 TH 6
+        9.29E-01 -2.15E+00 -1.51E+02 -9.24E-01  6.02E-01  1.69E+02
 
 TH 7
+        7.85E-01 -1.36E+01  5.48E+02 -8.51E+00 -5.50E+00  2.64E-01  2.12E+01
 
 TH 8
+        6.49E+01 -1.70E+03 -6.22E+06 -3.63E+03  1.23E+02  6.60E+01 -2.34E+02  2.72E+06
 
 TH 9
+        7.89E-01 -2.84E+01  1.10E+04  2.89E+01 -3.69E+00 -2.76E-01  2.92E+01 -4.77E+03  9.99E+01
 
 TH10
+       -8.49E+01  2.21E+03  8.21E+06 -1.20E+07 -2.17E+02 -8.68E+01  3.21E+02 -3.58E+06  6.33E+03  4.73E+06
 
 TH11
+       -6.72E+00 -1.34E+01  7.10E+03 -8.95E+00 -2.61E-01  9.76E-01  5.15E+00 -3.10E+03  3.77E+00  4.11E+03  4.28E+02
 
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
 
 Elapsed finaloutput time in seconds:     2.86
 #CPUT: Total CPU Time in Seconds,       29.385
Stop Time:
Wed Sep 29 20:48:25 CDT 2021
