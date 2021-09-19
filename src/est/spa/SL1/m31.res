Sat Sep 18 11:37:24 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat31.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m31.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1648.64877345839        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0233E+02 -4.4773E+01 -4.7947E+01 -2.3814E+00  8.7043E+01  1.4018E+00  6.6248E+00  8.8911E+00  4.8309E+00  2.8011E+00
            -1.3394E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1654.13295274233        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  9.6483E-01  9.7871E-01  1.0445E+00  1.0175E+00  9.3761E-01  9.9733E-01  8.8460E-01  9.0351E-01  9.9290E-01  9.0505E-01
             1.0256E+00
 PARAMETER:  6.4192E-02  7.8484E-02  1.4355E-01  1.1731E-01  3.5580E-02  9.7325E-02 -2.2617E-02 -1.4707E-03  9.2873E-02  2.3386E-04
             1.2527E-01
 GRADIENT:   2.1700E+01 -3.6332E+00  7.4155E+00 -1.1457E+01  3.2221E-01  2.4272E+00  1.3201E+00  1.7788E+00 -1.9354E+00 -5.8309E+00
            -6.0212E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1655.05662981502        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  9.6956E-01  9.9969E-01  9.7102E-01  1.0129E+00  9.2191E-01  9.9220E-01  7.8252E-01  6.4913E-01  1.0181E+00  9.9941E-01
             1.0339E+00
 PARAMETER:  6.9088E-02  9.9688E-02  7.0591E-02  1.1277E-01  1.8690E-02  9.2170E-02 -1.4523E-01 -3.3212E-01  1.1797E-01  9.9411E-02
             1.3333E-01
 GRADIENT:   2.9961E+01  6.4972E+00  4.9674E-01  9.6167E+00 -2.5958E+00  1.0573E-01 -4.6490E-01  1.1918E+00 -2.8361E+00  5.8378E+00
            -8.5064E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1655.52187370113        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.6340E-01  1.0632E+00  7.9912E-01  9.6490E-01  8.6482E-01  9.9355E-01  9.0155E-01  4.4484E-01  1.0198E+00  8.8880E-01
             1.0322E+00
 PARAMETER:  6.2717E-02  1.6129E-01 -1.2424E-01  6.4270E-02 -4.5239E-02  9.3526E-02 -3.6368E-03 -7.1003E-01  1.1962E-01 -1.7878E-02
             1.3170E-01
 GRADIENT:   1.2322E+01  7.7799E+00 -4.1699E+00  1.1299E+01 -1.1581E+00 -1.4732E-01  3.6314E-01  1.3274E+00  1.1702E+00  2.6255E+00
             1.0480E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1655.56444084728        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.5869E-01  1.0345E+00  7.9869E-01  9.7627E-01  8.5409E-01  9.9391E-01  9.3092E-01  3.4300E-01  1.0020E+00  8.8129E-01
             1.0333E+00
 PARAMETER:  5.7814E-02  1.3388E-01 -1.2478E-01  7.5988E-02 -5.7719E-02  9.3892E-02  2.8416E-02 -9.7002E-01  1.0201E-01 -2.6373E-02
             1.3274E-01
 GRADIENT:   1.4882E+00  8.1150E-01 -1.9692E+00  2.2536E+00  3.2557E-01 -2.0904E-01  2.6882E-01  6.6100E-01  4.9179E-01  1.0688E+00
             3.0773E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1655.63503599411        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.5692E-01  1.0163E+00  7.8018E-01  9.8263E-01  8.3798E-01  9.9465E-01  9.5862E-01  1.6325E-01  9.8803E-01  8.7026E-01
             1.0333E+00
 PARAMETER:  5.5969E-02  1.1617E-01 -1.4823E-01  8.2480E-02 -7.6765E-02  9.4633E-02  5.7741E-02 -1.7125E+00  8.7955E-02 -3.8962E-02
             1.3275E-01
 GRADIENT:  -2.8891E+00 -2.7450E+00 -8.2666E-01 -2.2067E+00  1.0280E+00 -9.6325E-02  1.8299E-01  1.3904E-01  9.1797E-02  2.5236E-01
             3.2288E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1655.81934855501        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      472
 NPARAMETR:  9.5934E-01  1.0381E+00  7.8529E-01  9.7205E-01  8.4888E-01  9.9516E-01  9.3963E-01  3.0778E-02  9.9938E-01  8.7478E-01
             1.0329E+00
 PARAMETER:  5.8492E-02  1.3739E-01 -1.4170E-01  7.1649E-02 -6.3840E-02  9.5152E-02  3.7726E-02 -3.3810E+00  9.9377E-02 -3.3788E-02
             1.3237E-01
 GRADIENT:  -3.5716E+01 -1.8359E+00  4.1124E+00 -5.7452E+00 -3.7453E+00 -3.8023E+00 -4.1304E-01  2.2785E-03 -9.7326E-01 -1.9049E+00
            -1.0170E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1656.24479243108        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      649
 NPARAMETR:  9.7427E-01  1.1737E+00  7.4753E-01  8.9011E-01  8.9367E-01  1.0030E+00  8.5156E-01  1.0000E-02  1.0794E+00  8.9728E-01
             1.0354E+00
 PARAMETER:  7.3930E-02  2.6017E-01 -1.9098E-01 -1.6414E-02 -1.2419E-02  1.0301E-01 -6.0684E-02 -6.9758E+00  1.7636E-01 -8.3843E-03
             1.3478E-01
 GRADIENT:  -2.3193E+00  6.9615E-01  3.8909E-01  4.1581E-01 -4.5913E-01 -2.5633E-01 -4.2429E-02  0.0000E+00 -2.0297E-01 -2.1099E-01
             1.8615E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1656.24786455094        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      811
 NPARAMETR:  9.7534E-01  1.1912E+00  7.4239E-01  8.7859E-01  8.9991E-01  1.0037E+00  8.4106E-01  1.0000E-02  1.0918E+00  9.0079E-01
             1.0353E+00
 PARAMETER:  7.5030E-02  2.7496E-01 -1.9788E-01 -2.9432E-02 -5.4601E-03  1.0369E-01 -7.3091E-02 -7.1856E+00  1.8783E-01 -4.4879E-03
             1.3471E-01
 GRADIENT:   2.5004E-02  1.5000E-02  5.2973E-03  1.2704E-02 -8.5158E-03  3.8236E-03 -1.0710E-03  0.0000E+00  5.7721E-04 -4.1904E-04
            -2.0568E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      811
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8992E-05 -2.0098E-02 -3.3597E-04  8.5949E-03 -2.3044E-02
 SE:             2.9816E-02  1.9247E-02  1.5355E-04  2.5263E-02  2.3760E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9949E-01  2.9639E-01  2.8663E-02  7.3369E-01  3.3213E-01

 ETASHRINKSD(%)  1.1237E-01  3.5521E+01  9.9486E+01  1.5366E+01  2.0400E+01
 ETASHRINKVR(%)  2.2461E-01  5.8424E+01  9.9997E+01  2.8371E+01  3.6638E+01
 EBVSHRINKSD(%)  4.6488E-01  3.5094E+01  9.9525E+01  1.5618E+01  1.9539E+01
 EBVSHRINKVR(%)  9.2761E-01  5.7873E+01  9.9998E+01  2.8796E+01  3.5260E+01
 RELATIVEINF(%)  9.8896E+01  1.6227E+00  2.4894E-04  3.7171E+00  6.7057E+00
 EPSSHRINKSD(%)  4.3642E+01
 EPSSHRINKVR(%)  6.8238E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1656.2478645509432     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -921.09703798720500     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.61
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1656.248       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  1.19E+00  7.42E-01  8.79E-01  9.00E-01  1.00E+00  8.41E-01  1.00E-02  1.09E+00  9.01E-01  1.04E+00
 


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
+        1.15E+03
 
 TH 2
+       -6.70E+00  4.90E+02
 
 TH 3
+        1.24E+01  2.46E+02  4.91E+02
 
 TH 4
+       -1.64E+01  3.91E+02 -1.99E+02  8.63E+02
 
 TH 5
+        5.50E-01 -4.09E+02 -6.06E+02  2.19E+02  1.09E+03
 
 TH 6
+        3.24E+00 -1.16E+00  3.05E+00 -5.54E+00 -9.99E-01  1.96E+02
 
 TH 7
+        4.61E-01  1.17E+01  4.72E+00 -6.22E+00 -1.81E+01 -3.74E-01  5.29E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.25E+00 -2.83E+01 -2.32E+01  3.73E+01  6.48E-01 -9.08E-01  2.55E+01  0.00E+00  8.90E+01
 
 TH10
+        6.41E-01 -8.42E+00 -5.48E+01 -1.19E+01 -5.74E+01  1.46E-01  2.30E+01  0.00E+00  4.65E+00  1.06E+02
 
 TH11
+       -9.37E+00 -1.93E+01 -3.33E+01 -2.95E+00  7.89E+00  1.48E+00  7.97E+00  0.00E+00  7.81E+00  2.32E+01  1.99E+02
 
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
 #CPUT: Total CPU Time in Seconds,       13.196
Stop Time:
Sat Sep 18 11:37:39 CDT 2021
