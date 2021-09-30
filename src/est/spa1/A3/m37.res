Thu Sep 30 00:12:26 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat37.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m37.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -600.732113066380        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2222E+02  4.5952E+01  1.4809E+02 -2.7310E+01  5.6885E+01  3.6049E+01 -1.8533E+01 -6.2078E+01  4.9437E+00 -6.8069E+01
            -2.8338E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1625.98528249541        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0129E+00  1.1169E+00  9.3301E-01  1.0548E+00  1.1154E+00  8.4794E-01  1.0214E+00  9.6974E-01  7.6209E-01  8.8583E-01
             3.2202E+00
 PARAMETER:  1.1283E-01  2.1054E-01  3.0663E-02  1.5335E-01  2.0924E-01 -6.4946E-02  1.2120E-01  6.9270E-02 -1.7169E-01 -2.1229E-02
             1.2694E+00
 GRADIENT:   3.0253E+01  1.9185E+01 -2.6876E+01  7.2427E+01  2.2222E+01 -3.7499E+01 -1.2183E+01  4.1437E+00 -9.3127E-01  1.1056E+01
             1.0560E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1634.24030884152        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0081E+00  1.0064E+00  4.6103E-01  1.0340E+00  6.4270E-01  9.8885E-01  1.3829E+00  7.4490E-01  4.7227E-01  2.9664E-01
             3.1238E+00
 PARAMETER:  1.0810E-01  1.0640E-01 -6.7429E-01  1.3344E-01 -3.4207E-01  8.8783E-02  4.2420E-01 -1.9451E-01 -6.5020E-01 -1.1152E+00
             1.2390E+00
 GRADIENT:   1.3477E+01  4.1240E+01  1.3205E+01  3.6014E+01 -3.4686E+01  1.5098E+01  6.1086E-01 -2.3328E-01 -1.5081E+01  2.2129E+00
             9.6914E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1644.91077966955        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9354E-01  1.0251E+00  3.4904E-01  9.7157E-01  5.6891E-01  9.2163E-01  1.2227E+00  8.7136E-01  8.0934E-01  1.5434E-01
             2.6810E+00
 PARAMETER:  9.3523E-02  1.2482E-01 -9.5257E-01  7.1157E-02 -4.6404E-01  1.8393E-02  3.0102E-01 -3.7696E-02 -1.1154E-01 -1.7686E+00
             1.0862E+00
 GRADIENT:  -5.3246E+00  3.1084E+01  1.9431E+01 -1.1752E+01 -2.7783E+01 -1.7248E+01  3.6849E+00 -2.5868E-01  3.9496E+00  6.4928E-01
             1.0634E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1647.38785857214        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      367
 NPARAMETR:  1.0126E+00  8.8292E-01  3.5523E-01  1.0526E+00  5.2653E-01  9.7016E-01  1.3139E+00  8.2149E-01  7.7495E-01  9.8708E-02
             2.6678E+00
 PARAMETER:  1.1249E-01 -2.4517E-02 -9.3499E-01  1.5122E-01 -5.4145E-01  6.9702E-02  3.7304E-01 -9.6641E-02 -1.5496E-01 -2.2156E+00
             1.0813E+00
 GRADIENT:  -7.2156E+00  6.8623E+00 -6.3859E+00  1.6380E+01  1.7474E+00 -1.5962E+00 -2.7756E+00 -1.0880E+00  1.7606E+00  1.5308E-01
            -2.9070E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1648.06904289055        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      542
 NPARAMETR:  1.0146E+00  7.4080E-01  3.7148E-01  1.1121E+00  4.8676E-01  9.6928E-01  1.5182E+00  8.6779E-01  7.2479E-01  4.8666E-02
             2.6640E+00
 PARAMETER:  1.1454E-01 -2.0003E-01 -8.9027E-01  2.0624E-01 -6.1999E-01  6.8794E-02  5.1752E-01 -4.1808E-02 -2.2187E-01 -2.9228E+00
             1.0798E+00
 GRADIENT:  -2.0939E-01  2.4078E-01 -1.0071E-01 -7.2948E-01  3.6506E-01 -1.1362E+00  5.1794E-02  3.5050E-01 -1.3904E+00  2.9295E-02
             1.3834E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1648.11828111520        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      717
 NPARAMETR:  1.0145E+00  7.0052E-01  3.5866E-01  1.1243E+00  4.6349E-01  9.7352E-01  1.5558E+00  8.6422E-01  7.3871E-01  3.7700E-02
             2.6445E+00
 PARAMETER:  1.1439E-01 -2.5594E-01 -9.2538E-01  2.1716E-01 -6.6896E-01  7.3160E-02  5.4199E-01 -4.5927E-02 -2.0285E-01 -3.1781E+00
             1.0725E+00
 GRADIENT:  -8.2692E-03  1.1414E-02  9.8890E-02 -1.3776E-01 -1.7509E-01  1.5098E-02  8.4692E-03  8.5248E-03  4.0978E-02  1.1129E-02
             1.0915E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1648.12318639499        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      837
 NPARAMETR:  1.0145E+00  7.0060E-01  3.5828E-01  1.1240E+00  4.6323E-01  9.7342E-01  1.5562E+00  8.6489E-01  7.3837E-01  1.0000E-02
             2.6432E+00
 PARAMETER:  1.1435E-01 -2.5582E-01 -9.2643E-01  2.1687E-01 -6.6954E-01  7.3059E-02  5.4226E-01 -4.5152E-02 -2.0331E-01 -6.9125E+00
             1.0720E+00
 GRADIENT:  -5.7860E-02  9.7412E-02  1.7895E-01 -3.5365E-01 -3.8891E-01 -3.9664E-02  5.9865E-02 -5.2944E-02 -8.0256E-02  0.0000E+00
            -3.4293E-01

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1648.12328126013        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:      970
 NPARAMETR:  1.0147E+00  7.0074E-01  3.5843E-01  1.1240E+00  4.6332E-01  9.7360E-01  1.5559E+00  8.6631E-01  7.3883E-01  1.0000E-02
             2.6442E+00
 PARAMETER:  1.1439E-01 -2.5583E-01 -9.2611E-01  2.1704E-01 -6.6921E-01  7.3128E-02  5.4191E-01 -4.4511E-02 -2.0306E-01 -5.8414E+00
             1.0723E+00
 GRADIENT:  -1.8804E-01 -4.9647E-02 -2.9630E-02  1.1071E-01  1.0069E-01 -1.8923E-02 -1.1904E-02 -1.4906E-02 -2.3639E-02  0.0000E+00
            -3.9556E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      970
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0620E-03  1.5869E-02 -1.9987E-02 -1.7625E-02  1.1692E-04
 SE:             2.9269E-02  2.2849E-02  1.5150E-02  2.3197E-02  2.7052E-04
 N:                     100         100         100         100         100

 P VAL.:         9.4384E-01  4.8735E-01  1.8707E-01  4.4736E-01  6.6560E-01

 ETASHRINKSD(%)  1.9438E+00  2.3454E+01  4.9247E+01  2.2287E+01  9.9094E+01
 ETASHRINKVR(%)  3.8499E+00  4.1407E+01  7.4241E+01  3.9607E+01  9.9992E+01
 EBVSHRINKSD(%)  2.1857E+00  2.2937E+01  5.0056E+01  2.2501E+01  9.9086E+01
 EBVSHRINKVR(%)  4.3237E+00  4.0614E+01  7.5056E+01  3.9939E+01  9.9992E+01
 RELATIVEINF(%)  9.5465E+01  7.5039E+00  3.4701E+00  1.9651E+01  5.2501E-04
 EPSSHRINKSD(%)  2.6536E+01
 EPSSHRINKVR(%)  4.6031E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1648.1232812601281     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -729.18474805545543     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.99
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.67
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1648.123       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  7.01E-01  3.58E-01  1.12E+00  4.63E-01  9.73E-01  1.56E+00  8.65E-01  7.39E-01  1.00E-02  2.64E+00
 


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
+        1.09E+03
 
 TH 2
+       -3.18E+01  5.67E+02
 
 TH 3
+        8.76E+00  6.78E+02  2.99E+03
 
 TH 4
+       -3.58E+01  2.50E+02 -6.79E+02  9.72E+02
 
 TH 5
+        4.14E+01 -1.18E+03 -3.28E+03  4.07E+02  4.23E+03
 
 TH 6
+        3.22E+00 -7.25E+00  1.44E+01 -1.33E+01  2.88E+00  1.90E+02
 
 TH 7
+        2.63E+00  3.91E+01 -4.93E+01 -1.15E+01  1.74E+01  2.33E-01  3.40E+01
 
 TH 8
+       -1.69E+00 -1.45E+01 -5.39E+01  4.62E-01  6.57E+01 -3.11E-01  5.66E+00  2.30E+01
 
 TH 9
+        4.89E+00 -2.26E+01  8.63E+00 -1.33E+01  8.17E+01 -2.03E+00  1.46E+01  3.59E+00  1.35E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.50E+01 -8.48E+00 -2.37E+01 -1.50E+01 -1.37E+01  4.22E+00  2.14E+00  1.06E+01  1.40E+01  0.00E+00  6.58E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.727
Stop Time:
Thu Sep 30 00:12:50 CDT 2021
