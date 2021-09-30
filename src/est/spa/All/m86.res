Wed Sep 29 20:42:54 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	One-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/28/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa/All/dat86.csv ignore=@
$SUBR ADVAN2 TRANS2
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
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
$OMEGA (0.09 FIX)x3
$SIGMA  1  FIX;        [P]
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
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   20064.7948181704        NO. OF FUNC. EVALS.:   9
 CUMULATIVE NO. OF FUNC. EVALS.:        9
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
 GRADIENT:   9.2035E+02  6.1311E+02 -5.6654E+02 -1.5498E+03 -1.8513E+03 -5.4899E+02 -3.8280E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -486.549940169127        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:       59
 NPARAMETR:  1.0612E+00  1.2201E+00  2.1637E+00  1.8413E+00  6.7204E-01  9.2642E-01  1.5386E+01
 PARAMETER:  1.5937E-01  2.9891E-01  8.7183E-01  7.1048E-01 -2.9744E-01  2.3567E-02  2.8335E+00
 GRADIENT:  -8.0923E+01  2.6562E+01 -4.7237E+00  5.8633E+01  4.1224E+00  7.3480E-01  1.1202E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -505.086203390282        NO. OF FUNC. EVALS.:  52
 CUMULATIVE NO. OF FUNC. EVALS.:      111
 NPARAMETR:  1.1848E+00  1.1264E+00  4.6977E+00  1.3845E+00  5.0631E-01  5.7715E+00  1.4896E+01
 PARAMETER:  2.6961E-01  2.1899E-01  1.6471E+00  4.2535E-01 -5.8060E-01  1.8529E+00  2.8011E+00
 GRADIENT:   2.4384E+01 -2.5479E+01 -2.4487E+00  9.6200E+00  4.2643E+00  1.0346E-01  6.2447E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -507.055585779193        NO. OF FUNC. EVALS.:  51
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  1.1337E+00  1.1059E+00  8.3310E+00  1.3665E+00  5.1706E-01  1.1102E+01  1.4001E+01
 PARAMETER:  2.2550E-01  2.0063E-01  2.2200E+00  4.1226E-01 -5.5959E-01  2.5071E+00  2.7392E+00
 GRADIENT:   7.2972E+00 -3.3041E+00 -1.3531E+00  2.8940E+00  3.5063E+00  1.0481E+00  1.9052E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -507.094944076223        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      212
 NPARAMETR:  1.1164E+00  1.0919E+00  9.9806E+00  1.3503E+00  4.7944E-01  1.2046E+01  1.3834E+01
 PARAMETER:  2.1011E-01  1.8789E-01  2.4006E+00  4.0031E-01 -6.3514E-01  2.5887E+00  2.7271E+00
 GRADIENT:   1.6495E+00 -2.0537E-01 -9.2914E-01  5.5817E-01  2.6533E+00  9.7511E-01  9.8890E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -507.122203934771        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      262
 NPARAMETR:  1.1066E+00  1.0802E+00  1.1506E+01  1.3392E+00  4.3249E-01  1.2625E+01  1.3765E+01
 PARAMETER:  2.0128E-01  1.7715E-01  2.5428E+00  3.9205E-01 -7.3819E-01  2.6357E+00  2.7221E+00
 GRADIENT:   1.0983E-01  2.3685E-01 -6.1187E-01 -1.8373E-01  1.7313E+00  6.7348E-01  4.9803E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -507.128771237966        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      312
 NPARAMETR:  1.1003E+00  1.0717E+00  1.3001E+01  1.3316E+00  3.9039E-01  1.3043E+01  1.3739E+01
 PARAMETER:  1.9561E-01  1.6926E-01  2.6650E+00  3.8640E-01 -8.4061E-01  2.6683E+00  2.7202E+00
 GRADIENT:  -3.5987E-01  3.3800E-01 -4.2963E-01 -4.1417E-01  1.0580E+00  4.5194E-01  2.1588E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -507.130263150604        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      362
 NPARAMETR:  1.0973E+00  1.0672E+00  1.4100E+01  1.3279E+00  3.6441E-01  1.3308E+01  1.3734E+01
 PARAMETER:  1.9281E-01  1.6500E-01  2.7462E+00  3.8361E-01 -9.0948E-01  2.6883E+00  2.7199E+00
 GRADIENT:  -3.9157E-01  3.1710E-01 -3.3322E-01 -4.4282E-01  7.1874E-01  3.1426E-01  1.0196E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -507.131228874765        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      412
 NPARAMETR:  1.0946E+00  1.0630E+00  1.5418E+01  1.3249E+00  3.3792E-01  1.3653E+01  1.3736E+01
 PARAMETER:  1.9040E-01  1.6111E-01  2.8355E+00  3.8135E-01 -9.8495E-01  2.7140E+00  2.7200E+00
 GRADIENT:  -3.3361E-01  1.9962E-01 -2.2132E-01 -3.6515E-01  4.3974E-01  1.9001E-01  3.5840E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -507.131581231054        NO. OF FUNC. EVALS.:  50
 CUMULATIVE NO. OF FUNC. EVALS.:      462
 NPARAMETR:  1.0929E+00  1.0602E+00  1.6562E+01  1.3230E+00  3.1775E-01  1.3965E+01  1.3744E+01
 PARAMETER:  1.8884E-01  1.5848E-01  2.9071E+00  3.7987E-01 -1.0465E+00  2.7365E+00  2.7206E+00
 GRADIENT:  -2.5868E-01  9.1334E-02 -1.3038E-01 -2.5205E-01  2.7063E-01  1.1411E-01  1.6400E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -507.847238453118        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      581            RESET HESSIAN, TYPE II
 NPARAMETR:  1.1313E+00  1.0910E+00  1.6704E+01  1.3453E+00  3.5829E-01  1.3863E+01  1.4333E+01
 PARAMETER:  2.2339E-01  1.8706E-01  2.9156E+00  3.9660E-01 -9.2641E-01  2.7292E+00  2.7626E+00
 GRADIENT:   1.1600E+01 -3.2138E+00 -3.1703E-02  4.1570E+00  7.9297E-01 -1.0879E-01  2.8486E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -507.922647501188        NO. OF FUNC. EVALS.: 103
 CUMULATIVE NO. OF FUNC. EVALS.:      684
 NPARAMETR:  1.1349E+00  1.0994E+00  2.0171E+01  1.3534E+00  3.5696E-01  1.5167E+01  1.4293E+01
 PARAMETER:  2.2654E-01  1.9480E-01  3.1042E+00  4.0260E-01 -9.3013E-01  2.8191E+00  2.7598E+00
 GRADIENT:   8.5299E-02  2.3233E-01 -1.3671E-02 -7.8047E-03  8.8875E-03 -2.0356E-02 -2.9180E-01

0ITERATION NO.:   56    OBJECTIVE VALUE:  -507.922647501188        NO. OF FUNC. EVALS.:  14
 CUMULATIVE NO. OF FUNC. EVALS.:      698
 NPARAMETR:  1.1349E+00  1.0994E+00  2.0171E+01  1.3534E+00  3.5696E-01  1.5167E+01  1.4293E+01
 PARAMETER:  2.2654E-01  1.9480E-01  3.1042E+00  4.0260E-01 -9.3013E-01  2.8191E+00  2.7598E+00
 GRADIENT:   8.5299E-02  2.3233E-01 -1.3671E-02 -7.8047E-03  8.8875E-03 -2.0356E-02 -2.9180E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      698
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1540E-02 -1.3550E-02 -7.6160E-03
 SE:             2.7273E-02  8.8862E-03  3.4467E-03
 N:                     100         100         100

 P VAL.:         6.7220E-01  1.2729E-01  2.7129E-02

 ETASHRINKSD(%)  8.6324E+00  7.0230E+01  8.8453E+01
 ETASHRINKVR(%)  1.6520E+01  9.1138E+01  9.8667E+01
 EBVSHRINKSD(%)  8.0126E+00  7.1355E+01  9.2416E+01
 EBVSHRINKVR(%)  1.5383E+01  9.1795E+01  9.9425E+01
 RELATIVEINF(%)  3.1914E+01  2.7512E+00  2.0012E-01
 EPSSHRINKSD(%)  3.7949E+00
 EPSSHRINKVR(%)  7.4458E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -507.92264750118812     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       227.22817906255005     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           300
  
 #TERE:
 Elapsed estimation  time in seconds:     5.43
 Elapsed covariance  time in seconds:     2.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -507.923       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7     
 
         1.13E+00  1.10E+00  2.02E+01  1.35E+00  3.57E-01  1.52E+01  1.43E+01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 


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
 
         1.02E-01  1.01E-01  9.67E+00  1.44E-01  3.94E-01  3.33E+00  1.76E+00
 


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
+        1.03E-02
 
 TH 2
+        8.76E-03  1.02E-02
 
 TH 3
+       -3.04E-01 -2.87E-01  9.35E+01
 
 TH 4
+        5.16E-04  4.52E-03  2.97E-01  2.07E-02
 
 TH 5
+        6.30E-03  1.91E-02 -1.25E+00  3.80E-02  1.55E-01
 
 TH 6
+       -4.35E-02 -6.52E-02  2.88E+01 -1.28E-02 -7.57E-01  1.11E+01
 
 TH 7
+        9.03E-02  5.15E-02  4.61E-01 -7.06E-02 -3.76E-01  2.13E+00  3.09E+00
 
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
+        1.02E-01
 
 TH 2
+        8.56E-01  1.01E-01
 
 TH 3
+       -3.10E-01 -2.95E-01  9.67E+00
 
 TH 4
+        3.53E-02  3.12E-01  2.13E-01  1.44E-01
 
 TH 5
+        1.57E-01  4.82E-01 -3.27E-01  6.68E-01  3.94E-01
 
 TH 6
+       -1.29E-01 -1.94E-01  8.93E-01 -2.67E-02 -5.77E-01  3.33E+00
 
 TH 7
+        5.05E-01  2.91E-01  2.71E-02 -2.79E-01 -5.43E-01  3.64E-01  1.76E+00
 
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
+        6.39E+02
 
 TH 2
+       -6.53E+02  1.02E+03
 
 TH 3
+       -1.85E+00  4.17E+00  1.60E-01
 
 TH 4
+        2.01E+01  1.88E+01 -1.65E+00  1.65E+02
 
 TH 5
+        6.58E+01 -1.75E+02 -6.65E-01 -6.31E+01  7.05E+01
 
 TH 6
+        8.98E+00 -1.79E+01 -4.79E-01  1.41E+00  4.74E+00  1.67E+00
 
 TH 7
+       -5.24E+00 -7.11E+00  1.72E-01 -5.55E+00  4.96E+00 -4.36E-01  1.35E+00
 
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
 #CPUT: Total CPU Time in Seconds,        7.609
Stop Time:
Wed Sep 29 20:43:03 CDT 2021
