Wed Sep 29 13:51:30 CDT 2021
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
$DATA ../../../../data/spa/A3/dat84.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m84.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   118.116034339739        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3437E+02  1.1476E+02  4.9615E+01  9.5445E+01  2.3958E+02  9.0949E+01 -8.7471E+01 -1.5596E+01 -1.4644E+02 -1.6839E+02
            -3.1604E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1267.08866697445        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0566E+00  9.6383E-01  1.0368E+00  1.0516E+00  1.0647E+00  6.3516E-01  9.8135E-01  9.4274E-01  1.0964E+00  7.8484E-01
             3.7835E+00
 PARAMETER:  1.5506E-01  6.3156E-02  1.3616E-01  1.5034E-01  1.6271E-01 -3.5387E-01  8.1173E-02  4.1039E-02  1.9208E-01 -1.4228E-01
             1.4307E+00
 GRADIENT:   1.4536E+02 -4.5608E+01 -4.1031E+01 -2.2628E+01  6.1781E+01 -4.8511E+01  3.3724E+00  5.1108E+00  3.5093E+00  1.1504E+01
            -5.5042E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1281.96039075914        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0543E+00  1.0784E+00  1.0647E+00  9.8314E-01  9.0257E-01  6.9312E-01  8.1517E-01  3.9922E-01  1.1228E+00  2.7651E-01
             4.2429E+00
 PARAMETER:  1.5291E-01  1.7545E-01  1.6273E-01  8.2994E-02 -2.5082E-03 -2.6655E-01 -1.0435E-01 -8.1825E-01  2.1580E-01 -1.1855E+00
             1.5452E+00
 GRADIENT:   9.0983E+01  2.6946E+01  2.2694E+01 -2.8565E+00 -5.7928E+01 -1.3725E+01  1.6220E+00  6.6006E-01  2.2404E+00  1.2940E+00
             4.1068E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1286.93387052547        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.0223E+00  8.5352E-01  9.9338E-01  1.1097E+00  8.6513E-01  7.4238E-01  7.8266E-01  9.5275E-02  1.0993E+00  2.7766E-01
             3.9619E+00
 PARAMETER:  1.2209E-01 -5.8387E-02  9.3362E-02  2.0412E-01 -4.4871E-02 -1.9790E-01 -1.4505E-01 -2.2510E+00  1.9466E-01 -1.1814E+00
             1.4767E+00
 GRADIENT:  -8.3604E+00 -9.5838E-01  7.5701E-01 -2.5824E+00 -9.6037E-01  2.7925E+00  6.2922E-01  7.8897E-02  2.3593E+00  1.4349E+00
             4.2450E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1288.05697907010        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0230E+00  5.7424E-01  8.3005E-01  1.2594E+00  6.8795E-01  7.3150E-01  1.0817E+00  3.3577E-02  9.9276E-01  4.6554E-02
             3.9211E+00
 PARAMETER:  1.2271E-01 -4.5470E-01 -8.6275E-02  3.3065E-01 -2.7404E-01 -2.1266E-01  1.7858E-01 -3.2939E+00  9.2737E-02 -2.9672E+00
             1.4664E+00
 GRADIENT:  -7.3550E-01  5.4328E+00  1.1462E+00  1.1752E+01 -3.3116E+00 -2.0734E+00 -1.6252E-01  1.5886E-02 -4.7266E-01  3.2466E-02
            -2.4732E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1288.81639959869        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      460
 NPARAMETR:  1.0378E+00  5.1909E-01  8.7689E-01  1.3156E+00  6.9831E-01  7.4795E-01  1.1388E+00  2.5666E-02  9.5391E-01  2.5473E-02
             4.0100E+00
 PARAMETER:  1.3707E-01 -5.5567E-01 -3.1378E-02  3.7430E-01 -2.5909E-01 -1.9042E-01  2.2994E-01 -3.5626E+00  5.2812E-02 -3.5701E+00
             1.4888E+00
 GRADIENT:   9.6292E+00  2.9603E+00  1.3657E+00  1.9898E+00 -4.8432E+00  1.7653E+00 -3.3641E-01  8.6587E-03 -6.2346E-01  9.0207E-03
             2.4959E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1289.43366394179        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      636
 NPARAMETR:  1.0289E+00  2.9889E-01  8.9122E-01  1.4356E+00  6.5158E-01  7.4268E-01  1.8492E+00  1.0000E-02  8.8641E-01  1.0000E-02
             3.9773E+00
 PARAMETER:  1.2851E-01 -1.1077E+00 -1.5164E-02  4.6158E-01 -3.2836E-01 -1.9749E-01  7.1476E-01 -5.5667E+00 -2.0581E-02 -7.4820E+00
             1.4806E+00
 GRADIENT:  -1.4526E+00  5.6242E-01  6.4686E-01 -4.1952E-01  1.5196E+00  5.9322E-01  2.0418E-01  0.0000E+00  5.9948E-01  0.0000E+00
            -3.4845E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1292.61568107925        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      817
 NPARAMETR:  1.0247E+00  1.1874E-01  4.1488E-01  1.4009E+00  3.5828E-01  7.4194E-01  3.6455E+00  1.0000E-02  9.9247E-01  1.0000E-02
             3.9037E+00
 PARAMETER:  1.2443E-01 -2.0308E+00 -7.7976E-01  4.3712E-01 -9.2645E-01 -1.9848E-01  1.3935E+00 -9.7357E+00  9.2438E-02 -1.5494E+01
             1.4619E+00
 GRADIENT:  -1.7291E+01 -3.6097E-01 -1.9802E+01  6.2934E+01  2.4681E+01 -4.0938E+00 -6.9568E+00  0.0000E+00 -7.5021E+00  0.0000E+00
             2.1102E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1296.14747523597        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      997
 NPARAMETR:  1.0224E+00  9.0516E-02  3.2788E-01  1.3379E+00  2.9647E-01  7.4826E-01  4.1687E+00  1.0000E-02  1.0816E+00  1.0000E-02
             3.8072E+00
 PARAMETER:  1.2215E-01 -2.3022E+00 -1.0151E+00  3.9108E-01 -1.1158E+00 -1.9001E-01  1.5276E+00 -1.1192E+01  1.7841E-01 -1.8142E+01
             1.4369E+00
 GRADIENT:  -2.2689E+01  5.1465E-01  7.9073E+00  6.4844E+01 -2.9905E+01 -5.4090E+00 -7.1871E+00  0.0000E+00 -3.4925E+00  0.0000E+00
             1.4654E+01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1296.35154571740        NO. OF FUNC. EVALS.:  68
 CUMULATIVE NO. OF FUNC. EVALS.:     1065
 NPARAMETR:  1.0236E+00  8.9434E-02  3.2670E-01  1.3341E+00  2.9537E-01  7.4928E-01  4.1916E+00  1.0000E-02  1.0844E+00  1.0000E-02
             3.8034E+00
 PARAMETER:  1.2213E-01 -2.3128E+00 -1.0192E+00  3.8814E-01 -1.1193E+00 -1.8877E-01  1.5340E+00 -1.1259E+01  1.8089E-01 -1.8262E+01
             1.4352E+00
 GRADIENT:  -1.9511E+01  4.8185E+01 -4.9789E+01 -9.5640E+01  2.3123E+01 -1.6827E+02  6.5090E+01  0.0000E+00 -1.7335E+02  0.0000E+00
            -6.5307E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1065
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.8720E-03  5.0937E-03  8.4111E-05 -4.6175E-02  2.5254E-04
 SE:             2.7824E-02  9.7250E-03  2.4190E-04  2.5986E-02  3.9373E-04
 N:                     100         100         100         100         100

 P VAL.:         8.3286E-01  6.0044E-01  7.2806E-01  7.5583E-02  5.2127E-01

 ETASHRINKSD(%)  6.7860E+00  6.7420E+01  9.9190E+01  1.2943E+01  9.8681E+01
 ETASHRINKVR(%)  1.3112E+01  8.9385E+01  9.9993E+01  2.4211E+01  9.9983E+01
 EBVSHRINKSD(%)  7.5327E+00  7.8313E+01  9.9190E+01  1.2117E+01  9.8842E+01
 EBVSHRINKVR(%)  1.4498E+01  9.5297E+01  9.9993E+01  2.2765E+01  9.9987E+01
 RELATIVEINF(%)  7.9830E+01  1.9460E+00  2.3773E-04  2.5095E+01  4.5227E-04
 EPSSHRINKSD(%)  2.2849E+01
 EPSSHRINKVR(%)  4.0478E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1296.3515457174024     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -561.20071915366418     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.37
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1296.352       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  8.96E-02  3.27E-01  1.33E+00  2.95E-01  7.49E-01  4.20E+00  1.00E-02  1.08E+00  1.00E-02  3.80E+00
 


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
+        1.00E+05
 
 TH 2
+       -6.22E+02  6.32E+04
 
 TH 3
+       -8.53E+01  5.63E+03  2.61E+04
 
 TH 4
+       -5.17E+01  8.60E+02  4.18E+03  9.16E+03
 
 TH 5
+        3.89E+02 -7.06E+03 -9.16E+03 -9.14E+03  3.09E+04
 
 TH 6
+       -8.69E+04  3.70E+02 -3.07E+01 -6.71E+01  1.57E+02  7.70E+04
 
 TH 7
+       -1.42E+01  8.84E+02  1.62E+02  2.29E+01 -1.84E+02  1.30E+01  6.53E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -6.27E+04  2.75E+02  1.02E+00 -6.90E+01  2.51E+02  5.53E+04  1.07E+01  0.00E+00  4.00E+04
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.32E+01  2.90E+02  1.14E+02  1.65E+02  5.29E+02 -1.02E+01  9.85E+00  0.00E+00  5.07E+02  0.00E+00  1.19E+02
 
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
 #CPUT: Total CPU Time in Seconds,       20.325
Stop Time:
Wed Sep 29 13:51:52 CDT 2021
