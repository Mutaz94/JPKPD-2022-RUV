Sat Sep 25 01:24:35 CDT 2021
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
$DATA ../../../../data/int/SL2/dat60.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      998
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      898
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -558.197103083377        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3883E+01  3.2317E+01  1.6320E+02 -7.7603E+01  1.1944E+02  3.2504E+00 -8.1058E+01 -3.9478E+02 -1.1784E+02 -1.8892E+01
            -6.0457E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2843.22652666931        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0082E+00  1.3494E+00  1.1087E+00  8.7472E-01  1.2273E+00  1.0048E+00  1.2001E+00  9.0305E-01  9.4060E-01  8.1141E-01
             2.2694E+00
 PARAMETER:  1.0818E-01  3.9963E-01  2.0316E-01 -3.3857E-02  3.0478E-01  1.0479E-01  2.8242E-01 -1.9739E-03  3.8763E-02 -1.0898E-01
             9.1950E-01
 GRADIENT:  -9.3149E+00  4.5728E+01 -7.0902E+00  3.0713E+01  3.5642E+01  2.9657E-01  9.5826E+00 -3.5998E+00 -1.6380E+01 -2.7769E+01
            -2.5566E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2850.28477777855        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0271E+00  1.4070E+00  8.9553E-01  8.3618E-01  1.1938E+00  9.5471E-01  9.8176E-01  2.4411E-01  9.7801E-01  1.0878E+00
             2.3698E+00
 PARAMETER:  1.2675E-01  4.4149E-01 -1.0337E-02 -7.8907E-02  2.7713E-01  5.3653E-02  8.1593E-02 -1.3102E+00  7.7761E-02  1.8414E-01
             9.6281E-01
 GRADIENT:   2.9651E+01  4.1288E+01 -2.6040E+01  3.6421E+01 -7.1680E+00 -1.9618E+01 -1.1983E+01 -2.2603E-01 -1.0710E+01  9.0231E+00
            -1.4996E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2858.59643963600        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0149E+00  1.2423E+00  1.0547E+00  9.1161E-01  1.1390E+00  9.9643E-01  1.1448E+00  4.3939E-01  9.8825E-01  9.3032E-01
             2.5036E+00
 PARAMETER:  1.1482E-01  3.1695E-01  1.5324E-01  7.4569E-03  2.3018E-01  9.6419E-02  2.3524E-01 -7.2237E-01  8.8185E-02  2.7776E-02
             1.0177E+00
 GRADIENT:  -1.9578E+00  2.1349E+00  3.1425E+00  6.9406E-01  8.4577E-01 -1.6687E+00  1.4135E+00 -9.1510E-01  2.6435E-01 -1.1736E+00
            -2.8979E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2859.20288584125        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0160E+00  1.2557E+00  1.0951E+00  9.0495E-01  1.1601E+00  9.9881E-01  1.1183E+00  8.5651E-01  9.9897E-01  9.3865E-01
             2.4978E+00
 PARAMETER:  1.1591E-01  3.2772E-01  1.9089E-01  1.2706E-04  2.4851E-01  9.8813E-02  2.1178E-01 -5.4889E-02  9.8972E-02  3.6691E-02
             1.0154E+00
 GRADIENT:   7.1891E-01 -7.2701E-01 -1.0936E+00  9.8374E-01  4.9708E-01 -6.9908E-01 -1.7678E-01  6.2119E-02  9.6297E-01  3.9645E-01
             2.7849E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2859.30075184789        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      378
 NPARAMETR:  1.0156E+00  1.2707E+00  1.1151E+00  8.9639E-01  1.1767E+00  1.0008E+00  1.1141E+00  9.0730E-01  9.9506E-01  9.4536E-01
             2.4972E+00
 PARAMETER:  1.1551E-01  3.3956E-01  2.0896E-01 -9.3745E-03  2.6274E-01  1.0078E-01  2.0802E-01  2.7156E-03  9.5048E-02  4.3815E-02
             1.0152E+00
 GRADIENT:  -9.5641E-02 -1.5740E-01 -1.4378E-01 -1.2465E-01 -1.2607E-02  5.5888E-02  3.4075E-02  8.4838E-02 -4.0655E-02  1.9910E-03
             4.2147E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2859.83788263179        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  1.0220E+00  1.4783E+00  1.2274E+00  7.9199E-01  1.3626E+00  1.0054E+00  1.0126E+00  1.2973E+00  1.0581E+00  1.0672E+00
             2.4961E+00
 PARAMETER:  1.2178E-01  4.9088E-01  3.0492E-01 -1.3320E-01  4.0937E-01  1.0535E-01  1.1249E-01  3.6030E-01  1.5651E-01  1.6500E-01
             1.0147E+00
 GRADIENT:   3.8178E+00  1.3064E+01 -1.6539E+00  1.7481E+01  4.8443E+00 -9.8728E-02 -2.4219E-01 -3.0236E-02 -2.3576E-01  7.6453E-01
            -1.7856E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2860.07198823610        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      693
 NPARAMETR:  1.0199E+00  1.5695E+00  1.2777E+00  7.2237E-01  1.4485E+00  1.0055E+00  9.7209E-01  1.5525E+00  1.1059E+00  1.1119E+00
             2.4926E+00
 PARAMETER:  1.1972E-01  5.5074E-01  3.4506E-01 -2.2522E-01  4.7053E-01  1.0548E-01  7.1693E-02  5.3988E-01  2.0064E-01  2.0605E-01
             1.0133E+00
 GRADIENT:   1.3614E-02  2.1669E-01 -4.4712E-02  2.6940E-01  2.7970E-01  1.0452E-03 -6.3592E-03 -3.6058E-02 -2.7897E-02 -4.0424E-03
            -2.5138E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2860.07218074725        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      868
 NPARAMETR:  1.0199E+00  1.5707E+00  1.2848E+00  7.2156E-01  1.4508E+00  1.0055E+00  9.7169E-01  1.5714E+00  1.1064E+00  1.1125E+00
             2.4924E+00
 PARAMETER:  1.1970E-01  5.5150E-01  3.5058E-01 -2.2634E-01  4.7214E-01  1.0547E-01  7.1278E-02  5.5196E-01  2.0107E-01  2.0664E-01
             1.0132E+00
 GRADIENT:  -6.4554E-03  4.9459E-02  1.0509E-02  2.5145E-02 -1.1781E-02 -1.2891E-03 -1.5876E-03 -6.6962E-03  1.1517E-02 -2.1496E-03
            -5.4779E-03

0ITERATION NO.:   42    OBJECTIVE VALUE:  -2860.07218406320        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      925
 NPARAMETR:  1.0199E+00  1.5708E+00  1.2857E+00  7.2147E-01  1.4512E+00  1.0055E+00  9.7168E-01  1.5742E+00  1.1062E+00  1.1126E+00
             2.4924E+00
 PARAMETER:  1.1971E-01  5.5159E-01  3.5134E-01 -2.2647E-01  4.7239E-01  1.0547E-01  7.1275E-02  5.5374E-01  2.0094E-01  2.0673E-01
             1.0132E+00
 GRADIENT:   1.1097E-03  8.1914E-03  2.9130E-03  4.2405E-03 -7.8632E-03  4.4522E-05 -6.9638E-04 -4.3971E-04  5.1980E-03  2.1951E-04
            -3.0005E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      925
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4822E-03 -1.7490E-02 -1.6218E-02  1.3907E-02 -2.1469E-02
 SE:             2.9475E-02  2.3836E-02  1.1274E-02  2.1109E-02  2.3411E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5989E-01  4.6308E-01  1.5028E-01  5.1001E-01  3.5911E-01

 ETASHRINKSD(%)  1.2558E+00  2.0148E+01  6.2231E+01  2.9282E+01  2.1570E+01
 ETASHRINKVR(%)  2.4959E+00  3.6237E+01  8.5735E+01  4.9989E+01  3.8487E+01
 EBVSHRINKSD(%)  1.4204E+00  2.0536E+01  6.4923E+01  3.2081E+01  1.9737E+01
 EBVSHRINKVR(%)  2.8206E+00  3.6855E+01  8.7696E+01  5.3870E+01  3.5579E+01
 RELATIVEINF(%)  9.7138E+01  9.6627E+00  5.9724E+00  7.0098E+00  1.5505E+01
 EPSSHRINKSD(%)  1.6803E+01
 EPSSHRINKVR(%)  3.0782E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          898
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1650.4136056355921     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2860.0721840631991     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1209.6585784276069     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.55
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2860.072       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.57E+00  1.29E+00  7.21E-01  1.45E+00  1.01E+00  9.72E-01  1.57E+00  1.11E+00  1.11E+00  2.49E+00
 


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
+        1.03E+03
 
 TH 2
+       -7.89E+00  3.09E+02
 
 TH 3
+        1.99E+00  2.39E+01  3.47E+01
 
 TH 4
+       -1.79E+01  3.24E+02 -6.08E+01  7.86E+02
 
 TH 5
+       -3.96E+00 -1.02E+02 -6.20E+01  1.63E+02  2.74E+02
 
 TH 6
+        3.01E+00 -2.68E+00  3.92E-01 -5.99E+00 -1.61E+00  1.85E+02
 
 TH 7
+        8.06E-01  1.22E+01 -3.21E+00 -9.37E+00  7.15E-01  3.10E+00  8.78E+01
 
 TH 8
+       -1.56E-01 -4.36E+00 -8.34E+00  5.33E+00 -1.84E+00 -9.63E-02  1.53E+00  4.32E+00
 
 TH 9
+        1.24E+00 -6.82E+00  7.64E+07  2.36E+01  6.87E+00 -1.09E+00  2.11E+01  1.93E+00  3.81E+01
 
 TH10
+       -1.08E-01 -1.06E+01  4.45E-01  1.26E+01 -2.07E+01  1.54E+00  1.14E+01  1.59E+00  4.87E+00  6.83E+01
 
 TH11
+       -1.37E+01 -1.23E+01 -4.29E-01 -9.71E+00 -1.91E+00  2.25E+00  2.95E+00  2.62E+00  7.06E+00  5.95E+00  1.86E+02
 
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
 #CPUT: Total CPU Time in Seconds,       32.316
Stop Time:
Sat Sep 25 01:25:09 CDT 2021
