Sat Sep 25 14:58:26 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	One-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/All/dat25.csv ignore=@
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
$PK

ET1 = EXP(ETA(1)*THETA(4))
ET2 = EXP(ETA(2)*THETA(5))
ET3 = EXP(ETA(3)*THETA(6))


CL = 5.0 * THETA(1) * ET1
V = 85  * THETA(2) * ET2
KA = 0.7 * THETA(3) * ET3

SC = V
$ERROR
CVERR 	= 0.05
W  	= THETA(7)*F*CVERR
Y  	= F + W * ERR(1)
$THETA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvKA
(0,1) ; tvCL
(0,1) ; tvV
(0,1) ; tvK
(0,1) ; RUV
$OMEGA
0.9 FIX ;     IIV CL
0.9 FIX  ;     IIV V
0.9 FIX ;      IIV KA
$SIGMA  1  FIX;        [P]
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
0LENGTH OF THETA:   7
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   3
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
0INITIAL ESTIMATE OF OMEGA:
 0.9000E+00
 0.0000E+00   0.9000E+00
 0.0000E+00   0.0000E+00   0.9000E+00
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

 ONE COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN2)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   3
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   ELIMINATION RATE (K) IS BASIC PK PARAMETER NO.:  1
   ABSORPTION RATE (KA) IS BASIC PK PARAMETER NO.:  3

 TRANSLATOR WILL CONVERT PARAMETERS
 CLEARANCE (CL) AND VOLUME (V) TO K (TRANS2)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            4           *           *           *           *
    3            *           -           -           -           -
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   5610.68765026675        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   4.7473E+03 -2.0840E+03  1.3141E+03  5.9019E+03 -1.2321E+02 -1.9924E+03 -1.3593E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -654.589629058330        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:       60
 NPARAMETR:  1.8417E+00  1.4519E+00  2.8288E+00  4.9484E-01  3.5456E-01  6.1870E-01  1.5063E+01
 PARAMETER:  7.1068E-01  4.7287E-01  1.1399E+00 -6.0353E-01 -9.3688E-01 -3.8013E-01  2.8122E+00
 GRADIENT:   2.1712E+02 -1.0301E+02 -1.2025E+00 -7.2769E+00  1.3642E+01  8.9162E-01  2.1896E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -767.736416924683        NO. OF FUNC. EVALS.:  53
 CUMULATIVE NO. OF FUNC. EVALS.:      113
 NPARAMETR:  1.0883E+00  1.1643E+00  3.3841E+01  4.2815E-01  2.5111E-01  7.3206E+00  9.5540E+00
 PARAMETER:  1.8461E-01  2.5212E-01  3.6217E+00 -7.4827E-01 -1.2819E+00  2.0907E+00  2.3570E+00
 GRADIENT:  -2.7736E+01 -5.2109E+01 -8.9294E-01  2.1682E+01  1.7905E-01  1.4151E+01  6.7264E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -787.495926398094        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      164
 NPARAMETR:  1.1422E+00  1.2857E+00  4.9518E+01  4.8345E-01  4.5892E-01  1.3211E+00  7.8603E+00
 PARAMETER:  2.3293E-01  3.5128E-01  4.0023E+00 -6.2680E-01 -6.7887E-01  3.7844E-01  2.1618E+00
 GRADIENT:  -3.5343E-01 -2.3286E+00 -4.6753E-02  2.1271E+00 -1.0419E+00 -3.4354E-03 -5.0346E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -787.537398829627        NO. OF FUNC. EVALS.:  53
 CUMULATIVE NO. OF FUNC. EVALS.:      217
 NPARAMETR:  1.1428E+00  1.2900E+00  5.1006E+01  4.7800E-01  4.5634E-01  1.3773E+00  7.9377E+00
 PARAMETER:  2.3345E-01  3.5467E-01  4.0319E+00 -6.3814E-01 -6.8451E-01  4.2013E-01  2.1716E+00
 GRADIENT:  -1.1190E+00  6.4703E-01 -4.4535E-02  4.7585E-02  1.0402E-01 -3.5687E-03  3.9212E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -787.574174744100        NO. OF FUNC. EVALS.:  54
 CUMULATIVE NO. OF FUNC. EVALS.:      271
 NPARAMETR:  1.1427E+00  1.2856E+00  1.2708E+03  4.7685E-01  4.5539E-01  2.0626E+01  7.9413E+00
 PARAMETER:  2.3343E-01  3.5126E-01  7.2474E+00 -6.4055E-01 -6.8659E-01  3.1266E+00  2.1721E+00
 GRADIENT:   4.1267E-01 -6.3123E-01 -6.0424E-04 -1.2645E-01  1.8191E-01 -1.2437E-03  1.2850E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -787.593183298257        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  1.1464E+00  1.2927E+00  1.7551E+03  4.8149E-01  4.5981E-01  2.6721E+01  7.9403E+00
 PARAMETER:  2.3663E-01  3.5671E-01  7.5703E+00 -6.3087E-01 -6.7694E-01  3.3855E+00  2.1720E+00
 GRADIENT:  -1.8156E-01 -2.4786E-01 -1.2410E-04 -1.8382E-01 -2.8851E-02 -1.0671E-03 -2.0892E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -787.593897840463        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      469
 NPARAMETR:  1.1469E+00  1.2936E+00  3.4436E+03  4.8212E-01  4.6037E-01  4.7026E+01  7.9408E+00
 PARAMETER:  2.3704E-01  3.5745E-01  8.2443E+00 -6.2957E-01 -6.7573E-01  3.9507E+00  2.1720E+00
 GRADIENT:  -9.5892E-03  1.0735E-02 -9.6026E-04  6.4805E-03  1.9609E-03  1.5346E-03  1.6054E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -787.593904654702        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      585            RESET HESSIAN, TYPE II
 NPARAMETR:  1.1469E+00  1.2936E+00  3.5263E+03  4.8210E-01  4.6036E-01  4.7494E+01  7.9406E+00
 PARAMETER:  2.3705E-01  3.5743E-01  8.2680E+00 -6.2960E-01 -6.7575E-01  3.9606E+00  2.1720E+00
 GRADIENT:   1.9741E+00  1.4234E+00 -4.8173E-04  1.3790E+00  4.6000E-01  7.0135E-04  1.7309E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -787.594025095152        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  1.1468E+00  1.2935E+00  1.6625E+04  4.8186E-01  4.6009E-01  1.0013E+02  7.9422E+00
 PARAMETER:  2.3699E-01  3.5732E-01  9.8187E+00 -6.3011E-01 -6.7634E-01  4.7064E+00  2.1722E+00
 GRADIENT:  -1.5177E-02  1.6738E-02  2.0223E-05 -6.9412E-02 -1.3183E-02 -1.2374E-04  2.3453E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -787.594088740467        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      820
 NPARAMETR:  1.1469E+00  1.2936E+00  3.5642E+05  4.8209E-01  4.6034E-01  4.4970E+02  7.9407E+00
 PARAMETER:  2.3702E-01  3.5740E-01  1.2884E+01 -6.2963E-01 -6.7579E-01  6.2086E+00  2.1720E+00
 GRADIENT:  -2.2429E-03  8.4572E-04  1.1317E-06 -9.4157E-04 -1.1593E-03 -6.5175E-06 -1.1343E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -787.594090990800        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      936
 NPARAMETR:  1.1469E+00  1.2936E+00  8.7889E+07  4.8209E-01  4.6035E-01  6.6603E+03  7.9407E+00
 PARAMETER:  2.3702E-01  3.5740E-01  1.8392E+01 -6.2962E-01 -6.7577E-01  8.9039E+00  2.1720E+00
 GRADIENT:  -1.1968E-03 -3.2173E-04  2.8311E-09  7.1292E-04  4.5666E-04 -2.5664E-08  1.4848E-03

0ITERATION NO.:   56    OBJECTIVE VALUE:  -787.594090990800        NO. OF FUNC. EVALS.:  17
 CUMULATIVE NO. OF FUNC. EVALS.:      953
 NPARAMETR:  1.1469E+00  1.2936E+00  8.7889E+07  4.8209E-01  4.6035E-01  6.6603E+03  7.9407E+00
 PARAMETER:  2.3702E-01  3.5740E-01  1.8392E+01 -6.2962E-01 -6.7577E-01  8.9039E+00  2.1720E+00
 GRADIENT:  -2.5244E-03  1.1919E-03 -7.2152E-07 -1.9536E-03  1.3343E-05 -3.8653E-06  3.1641E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      953
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.2025E-03 -7.3204E-02 -2.2020E-06
 SE:             8.6929E-02  7.9131E-02  2.0001E-06
 N:                     100         100         100

 P VAL.:         9.8896E-01  3.5491E-01  2.7094E-01

 ETASHRINKSD(%)  7.9076E+00  1.6169E+01  9.9998E+01
 ETASHRINKVR(%)  1.5190E+01  2.9723E+01  1.0000E+02
 EBVSHRINKSD(%)  6.4890E+00  1.4686E+01  9.9998E+01
 EBVSHRINKVR(%)  1.2557E+01  2.7214E+01  1.0000E+02
 RELATIVEINF(%)  1.7175E+01  1.5354E+01  3.5796E-09
 EPSSHRINKSD(%)  1.5474E+01
 EPSSHRINKVR(%)  2.8553E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -787.59409099080005     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -52.443264427061877     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     8.79
 Elapsed covariance  time in seconds:     2.15
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -787.594       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.15E+00  1.29E+00  8.79E+07  4.82E-01  4.60E-01  6.66E+03  7.94E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.00E-01
 
 ETA2
+        0.00E+00  9.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-01
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.49E-01
 
 ETA2
+        0.00E+00  9.49E-01
 
 ETA3
+        0.00E+00  0.00E+00  9.49E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         8.06E-02  1.30E-01  4.50E+08  9.45E-02  1.81E-01  1.53E+04  1.54E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        6.49E-03
 
 TH 2
+        6.41E-03  1.69E-02
 
 TH 3
+       -1.51E+06  4.37E+07  2.03E+17
 
 TH 4
+        1.20E-03  9.35E-03  3.81E+07  8.92E-03
 
 TH 5
+        3.42E-04  1.72E-02  7.43E+07  1.60E-02  3.28E-02
 
 TH 6
+        4.06E+02 -8.61E+02 -5.89E+12 -1.24E+03 -2.36E+03  2.33E+08
 
 TH 7
+        4.15E-02 -8.17E-02 -5.25E+08 -1.10E-01 -2.37E-01  2.19E+04  2.39E+00
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        8.06E-02
 
 TH 2
+        6.12E-01  1.30E-01
 
 TH 3
+       -4.17E-02  7.47E-01  4.50E+08
 
 TH 4
+        1.57E-01  7.61E-01  8.97E-01  9.45E-02
 
 TH 5
+        2.34E-02  7.31E-01  9.11E-01  9.38E-01  1.81E-01
 
 TH 6
+        3.30E-01 -4.33E-01 -8.58E-01 -8.59E-01 -8.55E-01  1.53E+04
 
 TH 7
+        3.34E-01 -4.06E-01 -7.55E-01 -7.54E-01 -8.48E-01  9.28E-01  1.54E+00
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      OM11      OM12      OM13      OM22      OM23  
             OM33      SG11  
 
 TH 1
+        2.60E+07
 
 TH 2
+       -1.58E+07  1.03E+07
 
 TH 3
+        3.58E-03 -2.46E-03  6.14E-13
 
 TH 4
+       -2.61E+07  1.41E+07 -2.80E-03  3.17E+07
 
 TH 5
+        6.93E+06 -3.49E+06  6.34E-04 -9.12E+06  2.70E+06
 
 TH 6
+       -1.02E+02  4.61E+01 -7.04E-09  1.49E+02 -4.57E+01  8.03E-04
 
 TH 7
+        2.10E+05 -3.29E+04 -1.27E-05 -4.90E+05  1.67E+05 -3.24E+00  1.61E+04
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       10.977
Stop Time:
Sat Sep 25 14:58:47 CDT 2021
