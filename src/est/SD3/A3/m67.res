Sat Oct 23 22:43:03 CDT 2021
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
$DATA ../../../../data/SD3/A3/dat67.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   809.001506040216        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2668E+02 -6.5730E+01  2.8512E+02 -2.1089E+02  2.3123E+02  5.1009E+01 -7.8444E+01 -3.2835E+02 -1.3103E+02 -1.3162E+02
            -5.0841E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1412.03400213441        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.8818E-01  1.1014E+00  9.4001E-01  1.2252E+00  9.9133E-01  7.4131E-01  9.8074E-01  9.9019E-01  9.4963E-01  9.3745E-01
             5.2971E+00
 PARAMETER:  8.8113E-02  1.9663E-01  3.8138E-02  3.0310E-01  9.1294E-02 -1.9934E-01  8.0552E-02  9.0144E-02  4.8317E-02  3.5407E-02
             1.7672E+00
 GRADIENT:  -9.2023E+01  9.8504E+00 -1.8418E+01  3.6234E+01 -6.8254E+00 -3.8826E+01  8.9450E+00  7.8029E+00  1.7904E+01  1.7955E+01
             2.9169E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1457.15985069004        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.8533E-01  8.2108E-01  2.4989E-01  1.2885E+00  3.5814E-01  8.9149E-01  9.2111E-01  5.1708E-02  1.1722E+00  2.2650E-01
             4.3968E+00
 PARAMETER:  8.5219E-02 -9.7133E-02 -1.2867E+00  3.5348E-01 -9.2683E-01 -1.4858E-02  1.7826E-02 -2.8622E+00  2.5890E-01 -1.3850E+00
             1.5809E+00
 GRADIENT:  -5.6477E+01  1.2211E+02  3.2446E+01  9.4504E+01 -1.0455E+02 -3.9357E+00  9.1242E-01  8.1133E-03 -8.1361E-01  5.4500E-01
             1.8327E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1484.56540559072        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.7830E-01  6.2425E-01  2.0081E-01  1.2862E+00  2.9238E-01  9.3489E-01  9.7042E-01  1.9106E-02  1.2644E+00  1.7917E-01
             3.5960E+00
 PARAMETER:  7.8063E-02 -3.7120E-01 -1.5054E+00  3.5168E-01 -1.1297E+00  3.2676E-02  6.9974E-02 -3.8578E+00  3.3459E-01 -1.6194E+00
             1.3798E+00
 GRADIENT:  -1.1830E+01  6.1031E+01 -1.4858E+01  8.4967E+01  5.5729E+00  7.1462E+00 -2.7443E+00 -9.8255E-03 -1.4799E+00 -2.0327E+00
             2.4822E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1494.41568999011        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      389
 NPARAMETR:  9.9473E-01  5.2709E-01  2.4786E-01  1.2742E+00  3.0302E-01  8.9545E-01  1.0249E+00  1.0000E-02  1.1092E+00  3.5865E-01
             3.3968E+00
 PARAMETER:  9.4718E-02 -5.4038E-01 -1.2949E+00  3.4231E-01 -1.0940E+00 -1.0427E-02  1.2455E-01 -4.6598E+00  2.0360E-01 -9.2539E-01
             1.3228E+00
 GRADIENT:   1.0552E+01  1.6239E+01  8.1831E+00  1.7900E+01 -1.2350E+01 -4.6585E+00  2.0155E-01  0.0000E+00 -1.2226E+01 -3.5699E+00
            -3.8923E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1500.96809419056        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      568
 NPARAMETR:  9.8449E-01  3.5303E-01  1.5720E-01  1.2055E+00  2.0486E-01  9.3038E-01  4.3156E-01  1.0000E-02  1.5695E+00  6.0399E-01
             3.4271E+00
 PARAMETER:  8.4373E-02 -9.4120E-01 -1.7503E+00  2.8693E-01 -1.4854E+00  2.7841E-02 -7.4035E-01 -7.4870E+00  5.5074E-01 -4.0419E-01
             1.3317E+00
 GRADIENT:  -1.1425E+01  8.5089E+00  3.2324E+01  1.5591E+01 -4.4335E+01  5.7849E+00  1.1921E+00  0.0000E+00  1.3426E+01  8.6896E-01
             3.4424E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1503.75763807922        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      744
 NPARAMETR:  9.8318E-01  3.1748E-01  1.2559E-01  1.1213E+00  1.8070E-01  9.1306E-01  3.2110E-01  1.0000E-02  1.6850E+00  6.5593E-01
             3.2524E+00
 PARAMETER:  8.3040E-02 -1.0473E+00 -1.9747E+00  2.1447E-01 -1.6109E+00  9.0409E-03 -1.0360E+00 -8.7353E+00  6.2175E-01 -3.2171E-01
             1.2794E+00
 GRADIENT:   8.0525E-01 -1.7148E+00 -5.2443E-01 -2.0106E+00 -1.3936E+00  7.4085E-02  1.1626E+00  0.0000E+00  3.0112E-01 -1.3442E+00
             2.9129E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1504.36342473507        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      921
 NPARAMETR:  9.8322E-01  3.2778E-01  1.2742E-01  1.1324E+00  1.8377E-01  9.1149E-01  7.7424E-02  1.0000E-02  1.6793E+00  6.6499E-01
             3.2423E+00
 PARAMETER:  8.3079E-02 -1.0154E+00 -1.9603E+00  2.2435E-01 -1.5941E+00  7.3201E-03 -2.4585E+00 -7.5010E+00  6.1837E-01 -3.0798E-01
             1.2763E+00
 GRADIENT:   1.5801E+00  5.8255E-01 -1.5283E+00  1.2358E+00  1.6702E+00 -8.1558E-01  6.7658E-02  0.0000E+00  1.3671E-01  3.6079E-01
             4.1317E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1504.40103347880        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1096
 NPARAMETR:  9.8263E-01  3.2598E-01  1.2733E-01  1.1299E+00  1.8333E-01  9.1348E-01  1.0868E-02  1.0000E-02  1.6777E+00  6.6426E-01
             3.2409E+00
 PARAMETER:  8.2479E-02 -1.0209E+00 -1.9610E+00  2.2209E-01 -1.5965E+00  9.5112E-03 -4.4220E+00 -6.0742E+00  6.1745E-01 -3.0908E-01
             1.2758E+00
 GRADIENT:  -4.3218E-02  6.0591E-02 -3.0074E-02  1.6547E-01 -1.1430E-02  1.8996E-03  1.3173E-03  0.0000E+00  4.6744E-03  2.4315E-02
            -6.7281E-02

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1504.40106030842        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1153
 NPARAMETR:  9.8261E-01  3.2588E-01  1.2730E-01  1.1297E+00  1.8330E-01  9.1353E-01  1.0000E-02  1.0000E-02  1.6779E+00  6.6423E-01
             3.2409E+00
 PARAMETER:  8.2460E-02 -1.0212E+00 -1.9612E+00  2.2195E-01 -1.5966E+00  9.5640E-03 -4.5765E+00 -5.9620E+00  6.1755E-01 -3.0913E-01
             1.2759E+00
 GRADIENT:  -8.6764E-02  3.3057E-02 -1.9351E-02  1.1376E-01 -2.2713E-02  2.3209E-02  0.0000E+00  0.0000E+00  1.1310E-02  7.2090E-03
            -3.9866E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1153
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.7560E-03 -6.7341E-05  3.6447E-04 -2.0570E-02  4.8354E-03
 SE:             2.8500E-02  1.1398E-04  2.1014E-04  2.5022E-02  2.4077E-02
 N:                     100         100         100         100         100

 P VAL.:         8.9515E-01  5.5463E-01  8.2840E-02  4.1103E-01  8.4083E-01

 ETASHRINKSD(%)  4.5204E+00  9.9618E+01  9.9296E+01  1.6175E+01  1.9338E+01
 ETASHRINKVR(%)  8.8365E+00  9.9999E+01  9.9995E+01  2.9733E+01  3.4936E+01
 EBVSHRINKSD(%)  4.1898E+00  9.9553E+01  9.9347E+01  1.1239E+01  2.0503E+01
 EBVSHRINKVR(%)  8.2040E+00  9.9998E+01  9.9996E+01  2.1215E+01  3.6802E+01
 RELATIVEINF(%)  9.1142E+01  4.4826E-04  5.2524E-04  5.1381E+01  4.1709E+00
 EPSSHRINKSD(%)  2.4518E+01
 EPSSHRINKVR(%)  4.3025E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1504.4010603084225     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -585.46252710374984     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1504.401       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  3.26E-01  1.27E-01  1.13E+00  1.83E-01  9.14E-01  1.00E-02  1.00E-02  1.68E+00  6.64E-01  3.24E+00
 


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
 #CPUT: Total CPU Time in Seconds,       95.333
Stop Time:
Sat Oct 23 22:43:18 CDT 2021
