Fri Sep 24 21:39:36 CDT 2021
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
$DATA ../../../../data/int/A2/dat67.csv ignore=@
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
Current Date:       24 SEP 2021
Days until program expires : 205
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3000.15367982591        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.0287E+01  3.9634E+01  1.4144E+02 -1.4708E+02  1.2689E+02 -9.4991E+00 -1.2896E+02 -9.3471E+01 -2.6691E+01 -8.3098E+01
            -1.4273E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3333.46042911647        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0080E+00  7.0020E-01  6.4763E-01  1.2291E+00  6.4627E-01  1.1253E+00  1.2822E+00  6.4865E-01  9.5644E-01  6.4655E-01
             1.3864E+00
 PARAMETER:  1.0798E-01 -2.5638E-01 -3.3443E-01  3.0631E-01 -3.3654E-01  2.1808E-01  3.4855E-01 -3.3286E-01  5.5468E-02 -3.3611E-01
             4.2669E-01
 GRADIENT:   4.7207E+01  3.9583E+01 -2.1598E+01  3.4701E+01  1.0303E+02  3.5545E+01 -1.1641E+01 -4.0545E-02 -2.2092E+01 -3.3729E+01
            -2.9178E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3338.89591557967        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      212
 NPARAMETR:  1.0197E+00  5.5621E-01  4.7754E-01  1.2949E+00  4.7812E-01  1.0691E+00  1.4950E+00  4.3575E-01  1.0387E+00  5.2391E-01
             1.3424E+00
 PARAMETER:  1.1952E-01 -4.8661E-01 -6.3910E-01  3.5843E-01 -6.3789E-01  1.6677E-01  5.0210E-01 -7.3070E-01  1.3800E-01 -5.4644E-01
             3.9444E-01
 GRADIENT:   5.2158E+01  7.4731E+01 -4.1827E+01  6.2394E+01  5.7322E+01  1.1224E+01  7.1673E+00 -1.3195E+01 -1.2384E+00 -4.6181E+01
            -3.2406E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3375.42616955613        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      390
 NPARAMETR:  1.0111E+00  3.8360E-01  3.4873E-01  1.2875E+00  3.4375E-01  1.0463E+00  1.4298E+00  1.6324E-01  1.0447E+00  6.4100E-01
             1.5081E+00
 PARAMETER:  1.1101E-01 -8.5817E-01 -9.5345E-01  3.5268E-01 -9.6784E-01  1.4525E-01  4.5754E-01 -1.7125E+00  1.4373E-01 -3.4473E-01
             5.1084E-01
 GRADIENT:   3.4938E+01  1.3490E+01  6.3177E+01  6.1308E+01 -9.0320E+01  4.7841E+00  8.7472E+00 -6.3145E-01  1.1654E-01 -2.2371E+01
            -8.3019E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3379.11151048811        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      565
 NPARAMETR:  9.9360E-01  3.6904E-01  3.3220E-01  1.2402E+00  3.3799E-01  1.0314E+00  1.3765E+00  9.7270E-02  1.0434E+00  7.1037E-01
             1.5055E+00
 PARAMETER:  9.3582E-02 -8.9685E-01 -1.0020E+00  3.1530E-01 -9.8475E-01  1.3094E-01  4.1951E-01 -2.2303E+00  1.4249E-01 -2.4197E-01
             5.0914E-01
 GRADIENT:  -3.0538E-02 -7.9035E-02 -2.3924E-01 -6.0349E-02  3.4603E-01  7.2753E-02 -4.4814E-02  1.9493E-02 -8.0508E-03 -2.6191E-02
            -8.7016E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3379.11792415287        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      684
 NPARAMETR:  9.9366E-01  3.6897E-01  3.3208E-01  1.2402E+00  3.3781E-01  1.0312E+00  1.3773E+00  1.0000E-02  1.0427E+00  7.1354E-01
             1.5062E+00
 PARAMETER:  9.3644E-02 -8.9705E-01 -1.0024E+00  3.1528E-01 -9.8527E-01  1.3070E-01  4.2012E-01 -4.5529E+00  1.4179E-01 -2.3752E-01
             5.0960E-01
 GRADIENT:   1.5837E+01  1.2768E+01  1.0914E+01  1.6706E+01  4.7260E+01  2.2342E+00  1.1698E+00  0.0000E+00  1.0689E+00  4.7693E-01
             8.2035E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3379.11812842002        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      761
 NPARAMETR:  9.9351E-01  3.6825E-01  3.3133E-01  1.2397E+00  3.3722E-01  1.0310E+00  1.3765E+00  1.0000E-02  1.0425E+00  7.1416E-01
             1.5061E+00
 PARAMETER:  9.3489E-02 -8.9899E-01 -1.0046E+00  3.1488E-01 -9.8703E-01  1.3056E-01  4.1955E-01 -1.3538E+01  1.4167E-01 -2.3665E-01
             5.0951E-01
 GRADIENT:   1.5511E+01  1.2574E+01  1.0601E+01  1.6425E+01  4.7782E+01  2.1796E+00  1.0067E+00  0.0000E+00  9.5050E-01  7.0650E-01
             8.9643E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3379.11812971077        NO. OF FUNC. EVALS.:  85
 CUMULATIVE NO. OF FUNC. EVALS.:      846
 NPARAMETR:  9.9350E-01  3.6815E-01  3.3121E-01  1.2397E+00  3.3712E-01  1.0310E+00  1.3765E+00  1.0000E-02  1.0426E+00  7.1417E-01
             1.5061E+00
 PARAMETER:  9.3483E-02 -8.9928E-01 -1.0050E+00  3.1484E-01 -9.8732E-01  1.3056E-01  4.1956E-01 -1.4947E+01  1.4169E-01 -2.3663E-01
             5.0950E-01
 GRADIENT:   1.5499E+01  1.2567E+01  1.0595E+01  1.6414E+01  4.7790E+01  2.1777E+00  1.0020E+00  0.0000E+00  9.4674E-01  7.1442E-01
             8.9973E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3379.11813896276        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:      973
 NPARAMETR:  9.9351E-01  3.6822E-01  3.3131E-01  1.2397E+00  3.3720E-01  1.0310E+00  1.3765E+00  1.0000E-02  1.0426E+00  7.1415E-01
             1.5061E+00
 PARAMETER:  9.3488E-02 -8.9907E-01 -1.0047E+00  3.1488E-01 -9.8707E-01  1.3057E-01  4.1956E-01 -1.3771E+01  1.4170E-01 -2.3667E-01
             5.0952E-01
 GRADIENT:  -2.2999E-01 -3.8463E-03 -2.0261E-01 -2.2189E-01  2.7304E-01 -5.6145E-02 -1.1297E-01  0.0000E+00 -1.1668E-01  1.9773E-01
            -2.5826E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3379.11832212909        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1135
 NPARAMETR:  9.9362E-01  3.6815E-01  3.3125E-01  1.2398E+00  3.3712E-01  1.0312E+00  1.3773E+00  1.0000E-02  1.0429E+00  7.1353E-01
             1.5061E+00
 PARAMETER:  9.3597E-02 -8.9928E-01 -1.0049E+00  3.1498E-01 -9.8732E-01  1.3072E-01  4.2010E-01 -1.4500E+01  1.4204E-01 -2.3753E-01
             5.0952E-01
 GRADIENT:  -1.8498E-03 -1.0197E-03  1.3800E-03 -5.2963E-04 -4.3244E-04 -2.9371E-05  5.1469E-04  0.0000E+00  3.4413E-04 -1.5850E-03
             1.7484E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1135
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.0166E-04  2.1801E-04 -1.3545E-04 -2.1261E-03 -1.9893E-03
 SE:             2.9745E-02  2.5647E-02  3.3133E-04  2.9198E-02  2.6488E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7850E-01  9.9322E-01  6.8268E-01  9.4195E-01  9.4013E-01

 ETASHRINKSD(%)  3.5183E-01  1.4081E+01  9.8890E+01  2.1822E+00  1.1262E+01
 ETASHRINKVR(%)  7.0242E-01  2.6179E+01  9.9988E+01  4.3168E+00  2.1256E+01
 EBVSHRINKSD(%)  6.0014E-01  1.2184E+01  9.8910E+01  2.1169E+00  1.1759E+01
 EBVSHRINKVR(%)  1.1967E+00  2.2883E+01  9.9988E+01  4.1890E+00  2.2135E+01
 RELATIVEINF(%)  9.8796E+01  2.2743E+01  1.1196E-03  8.6898E+01  7.3932E+00
 EPSSHRINKSD(%)  2.1125E+01
 EPSSHRINKVR(%)  3.7787E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3379.1183221290858     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1725.0289623606750     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.88
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3379.118       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.94E-01  3.68E-01  3.31E-01  1.24E+00  3.37E-01  1.03E+00  1.38E+00  1.00E-02  1.04E+00  7.14E-01  1.51E+00
 


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
+       -3.48E-01  3.36E+03
 
 TH 3
+        4.43E+00 -1.51E+03  1.10E+04
 
 TH 4
+       -3.00E+00  1.35E+01 -2.81E+02  6.37E+02
 
 TH 5
+       -2.00E+00 -1.70E+03 -1.10E+04 -1.44E+01  1.50E+04
 
 TH 6
+        4.39E+00 -2.78E+00  8.69E+00 -2.64E+00 -5.18E+00  1.81E+02
 
 TH 7
+       -8.20E-01  2.64E+01  2.26E+01  1.90E+00 -3.96E+01 -2.70E-01  6.03E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.41E+00 -1.61E+01  8.56E+01 -1.41E+00 -2.81E+01 -7.41E-01 -1.77E+00  0.00E+00  1.65E+02
 
 TH10
+       -5.46E-01 -1.55E+01  2.36E+00  9.30E+00 -2.49E+01  1.25E+00  5.09E+00  0.00E+00  5.59E-01  2.41E+02
 
 TH11
+       -9.48E+00 -8.95E+00 -1.65E+02 -1.09E+01  9.02E+01  2.11E+00  1.35E+01  0.00E+00  7.22E+00  2.43E+01  4.59E+02
 
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
 #CPUT: Total CPU Time in Seconds,       36.034
Stop Time:
Fri Sep 24 21:40:14 CDT 2021
