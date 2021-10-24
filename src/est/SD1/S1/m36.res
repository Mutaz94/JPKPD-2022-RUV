Sat Oct 23 14:32:36 CDT 2021
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
$DATA ../../../../data/SD1/S1/dat36.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
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
 (2E4.0,E21.0,E4.0,2E2.0)

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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m36.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3382.43933111327        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3112E+02  1.4471E+01  1.3729E+02  5.1661E+01  1.3301E+02  3.3858E+01 -3.7579E+01 -4.6320E+02 -8.8794E+01 -2.6680E+01
            -3.7251E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3596.20121263671        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.5973E-01  9.3447E-01  7.8932E-01  1.0862E+00  8.0593E-01  1.0066E+00  1.0844E+00  2.7411E+00  9.2986E-01  1.0865E+00
             1.2329E+00
 PARAMETER:  5.8898E-02  3.2228E-02 -1.3658E-01  1.8267E-01 -1.1576E-01  1.0660E-01  1.8102E-01  1.1084E+00  2.7279E-02  1.8294E-01
             3.0941E-01
 GRADIENT:   1.8194E+02  7.5570E+01 -6.2523E+00  1.7680E+02 -6.5757E+01  1.9822E+01 -1.2315E+00  5.2589E+01  1.3494E+01 -4.5040E+00
             1.9756E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3604.78779086890        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      266
 NPARAMETR:  9.6437E-01  1.1226E+00  9.5542E-01  1.0312E+00  9.8691E-01  1.0055E+00  1.0920E+00  2.9283E+00  9.4686E-01  1.2698E+00
             1.2812E+00
 PARAMETER:  6.3722E-02  2.1567E-01  5.4398E-02  1.3077E-01  8.6824E-02  1.0551E-01  1.8801E-01  1.1744E+00  4.5392E-02  3.3886E-01
             3.4779E-01
 GRADIENT:  -1.0720E+02  5.0220E+01 -4.9845E+00  1.1216E+02 -5.0673E+01 -3.8730E+01 -6.4787E+00  1.1479E+01  1.5285E+01  9.0712E-01
             2.3437E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3626.68675798591        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  9.8240E-01  1.1515E+00  9.1596E-01  9.8522E-01  1.0333E+00  1.0152E+00  1.1536E+00  2.9071E+00  8.4395E-01  1.5028E+00
             1.1647E+00
 PARAMETER:  8.2244E-02  2.4104E-01  1.2221E-02  8.5109E-02  1.3274E-01  1.1508E-01  2.4287E-01  1.1671E+00 -6.9663E-02  5.0734E-01
             2.5250E-01
 GRADIENT:  -6.3470E+01  2.8152E+01 -1.2682E+01  8.3655E+01 -2.5776E+01 -3.1314E+01 -5.8382E+00 -6.7832E-01  8.9763E+00  1.9980E+01
             8.0166E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3627.84469393859        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:      639
 NPARAMETR:  9.8258E-01  1.1515E+00  9.1613E-01  9.8533E-01  1.0332E+00  1.0686E+00  1.1765E+00  2.9010E+00  8.4410E-01  1.5023E+00
             1.1644E+00
 PARAMETER:  8.2425E-02  2.4104E-01  1.2401E-02  8.5225E-02  1.3266E-01  1.6634E-01  2.6253E-01  1.1650E+00 -6.9484E-02  5.0702E-01
             2.5220E-01
 GRADIENT:  -5.6697E+01  2.8340E+01 -1.2721E+01  8.3227E+01 -2.5699E+01 -8.6879E+00 -2.8301E+00 -1.3248E+00  1.0208E+01  1.9861E+01
             8.0720E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3628.81504641224        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      802
 NPARAMETR:  9.8432E-01  1.1515E+00  9.1612E-01  9.8022E-01  1.0332E+00  1.0827E+00  1.1764E+00  2.9014E+00  8.4409E-01  1.5022E+00
             1.1587E+00
 PARAMETER:  8.4196E-02  2.4104E-01  1.2390E-02  8.0017E-02  1.3267E-01  1.7949E-01  2.6249E-01  1.1652E+00 -6.9495E-02  5.0692E-01
             2.4726E-01
 GRADIENT:   2.9118E+02  2.0735E+02 -7.8164E+00  1.3936E+02  2.2873E+01  1.0511E+02  3.1869E+01  3.9010E+01  1.4826E+01  5.7311E+01
             7.5706E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3628.87894371884        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:      997
 NPARAMETR:  9.8420E-01  1.1518E+00  9.1612E-01  9.8009E-01  1.0332E+00  1.0902E+00  1.1867E+00  2.8972E+00  8.4420E-01  1.5012E+00
             1.1583E+00
 PARAMETER:  8.4071E-02  2.4134E-01  1.2390E-02  7.9891E-02  1.3267E-01  1.8636E-01  2.7121E-01  1.1638E+00 -6.9370E-02  5.0628E-01
             2.4695E-01
 GRADIENT:  -5.1189E+01  2.5163E+01 -1.2426E+01  7.2821E+01 -2.7261E+01 -2.1762E-01 -1.6217E+00 -1.8689E+00  1.1399E+01  1.9049E+01
             7.1223E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3630.03016077304        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1188             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8407E-01  1.1515E+00  9.2417E-01  9.7996E-01  1.0334E+00  1.0915E+00  1.1925E+00  2.9014E+00  8.4409E-01  1.5012E+00
             1.1364E+00
 PARAMETER:  8.3947E-02  2.4104E-01  2.1137E-02  7.9758E-02  1.3284E-01  1.8753E-01  2.7602E-01  1.1652E+00 -6.9495E-02  5.0628E-01
             2.2789E-01
 GRADIENT:   3.0548E+02  2.1644E+02 -5.8434E+00  1.3990E+02  2.4702E+01  1.1787E+02  3.6665E+01  3.9448E+01  1.5652E+01  5.7854E+01
             3.9883E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3630.15918394878        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1359             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8498E-01  1.1512E+00  9.2409E-01  9.7988E-01  1.0333E+00  1.0818E+00  1.1876E+00  2.9043E+00  8.4402E-01  1.5012E+00
             1.1329E+00
 PARAMETER:  8.4870E-02  2.4083E-01  2.1052E-02  7.9672E-02  1.3278E-01  1.7859E-01  2.7189E-01  1.1662E+00 -6.9581E-02  5.0628E-01
             2.2480E-01
 GRADIENT:   3.0853E+02  2.1750E+02 -5.7938E+00  1.3989E+02  2.4989E+01  1.0913E+02  3.5433E+01  3.9873E+01  1.5425E+01  5.8008E+01
             3.3759E+01

0ITERATION NO.:   44    OBJECTIVE VALUE:  -3630.17250260404        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:     1502
 NPARAMETR:  9.8486E-01  1.1509E+00  9.2398E-01  9.7976E-01  1.0332E+00  1.0864E+00  1.1872E+00  2.9084E+00  8.4392E-01  1.5012E+00
             1.1330E+00
 PARAMETER:  8.4808E-02  2.4068E-01  2.0990E-02  7.9610E-02  1.3278E-01  1.8106E-01  2.7174E-01  1.1669E+00 -6.9643E-02  5.0628E-01
             2.2484E-01
 GRADIENT:   2.5270E+04  1.0475E+04  2.5153E+04  2.5226E+04  1.8913E+04 -3.0205E+00  2.8875E+01 -2.1905E+03  1.2680E+04  2.3499E+01
            -1.4571E+01
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.3         1.1         2.3         2.3         2.3         4.9
                    2.8

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1502
 NO. OF SIG. DIGITS IN FINAL EST.:  1.1
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5486E-02 -5.3445E-02 -2.2392E-02  2.4879E-03 -2.9002E-02
 SE:             3.0023E-02  2.3661E-02  2.4995E-02  2.4105E-02  2.3533E-02
 N:                     100         100         100         100         100

 P VAL.:         3.9596E-01  2.3898E-02  3.7032E-01  9.1780E-01  2.1780E-01

 ETASHRINKSD(%)  1.0000E-10  2.0732E+01  1.6265E+01  1.9246E+01  2.1163E+01
 ETASHRINKVR(%)  1.0000E-10  3.7166E+01  2.9884E+01  3.4788E+01  3.7847E+01
 EBVSHRINKSD(%)  2.7784E-01  2.1677E+01  1.5398E+01  1.9114E+01  1.4946E+01
 EBVSHRINKVR(%)  5.5491E-01  3.8655E+01  2.8425E+01  3.4575E+01  2.7659E+01
 RELATIVEINF(%)  9.9443E+01  3.0506E+01  6.6710E+01  3.5204E+01  4.3529E+01
 EPSSHRINKSD(%)  2.4496E+01
 EPSSHRINKVR(%)  4.2991E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3630.1725026040385     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1976.0831428356278     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3630.173       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  1.15E+00  9.24E-01  9.80E-01  1.03E+00  1.08E+00  1.19E+00  2.91E+00  8.44E-01  1.50E+00  1.13E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,      125.594
Stop Time:
Sat Oct 23 14:32:56 CDT 2021
