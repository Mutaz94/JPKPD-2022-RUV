Sun Oct 24 00:04:05 CDT 2021
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
$DATA ../../../../data/SD3/SL2/dat85.csv ignore=@
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
Current Date:       24 OCT 2021
Days until program expires : 175
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
 RAW OUTPUT FILE (FILE): m85.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2038.69573981675        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5709E+02 -2.6776E+01 -1.7042E+01  9.1981E-01  2.8125E+00  6.7277E+01 -1.7265E+01  1.2738E+01 -5.2805E+00  8.1344E+00
            -1.2671E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2048.79050500368        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      114
 NPARAMETR:  9.6952E-01  1.0186E+00  1.0212E+00  9.9841E-01  1.0168E+00  9.1820E-01  1.0259E+00  9.8062E-01  1.0059E+00  9.8745E-01
             1.0822E+00
 PARAMETER:  6.9047E-02  1.1847E-01  1.2093E-01  9.8405E-02  1.1669E-01  1.4664E-02  1.2558E-01  8.0431E-02  1.0593E-01  8.7366E-02
             1.7904E-01
 GRADIENT:  -2.5610E+01 -5.3814E+01 -1.5809E+01 -6.1104E+01 -6.4031E+00 -1.0833E+01 -1.8706E+01  1.1034E+01 -7.6344E+00  5.4085E+00
            -4.9563E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2053.60284345071        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.7214E-01  1.0812E+00  1.0487E+00  1.0020E+00  1.0613E+00  9.4793E-01  1.3033E+00  8.6079E-01  8.6934E-01  1.0914E+00
             1.0956E+00
 PARAMETER:  7.1741E-02  1.7808E-01  1.4754E-01  1.0198E-01  1.5947E-01  4.6524E-02  3.6487E-01 -4.9899E-02 -4.0016E-02  1.8747E-01
             1.9127E-01
 GRADIENT:  -1.7482E+01  8.5115E-01 -5.4032E+00 -8.1731E+00 -1.6468E+01  2.0740E+00 -3.0104E+00  6.3676E+00 -1.4371E+01  1.5500E+01
            -3.8565E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2057.42716144471        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      479
 NPARAMETR:  9.7887E-01  9.0521E-01  1.0763E+00  1.1168E+00  9.9573E-01  9.4469E-01  1.4614E+00  6.0219E-01  9.0056E-01  9.6453E-01
             1.1389E+00
 PARAMETER:  7.8639E-02  4.0733E-04  1.7349E-01  2.1048E-01  9.5720E-02  4.3097E-02  4.7938E-01 -4.0718E-01 -4.7427E-03  6.3882E-02
             2.3011E-01
 GRADIENT:   2.0664E+00  3.4345E+00 -9.3804E-01  7.8570E+00 -1.1742E+00  1.4640E+00 -9.9080E-01  6.3039E-01 -1.3029E+00  1.6761E+00
            -6.6621E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2057.68963183114        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      657
 NPARAMETR:  9.7953E-01  1.0360E+00  8.8852E-01  1.0219E+00  9.5597E-01  9.4197E-01  1.3558E+00  3.4403E-01  9.3534E-01  8.7898E-01
             1.1475E+00
 PARAMETER:  7.9322E-02  1.3535E-01 -1.8196E-02  1.2165E-01  5.4974E-02  4.0223E-02  4.0438E-01 -9.6703E-01  3.3150E-02 -2.8995E-02
             2.3761E-01
 GRADIENT:  -1.6057E+00  2.8654E+00 -4.7695E-01  3.9581E+00 -1.8593E+00 -5.2243E-01 -2.4401E-01  2.4190E-01 -7.7186E-01  3.7535E-01
            -1.9961E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2057.73414168878        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      832
 NPARAMETR:  9.8092E-01  1.0997E+00  8.4884E-01  9.7737E-01  9.6788E-01  9.4425E-01  1.2943E+00  2.4194E-01  9.6644E-01  8.7492E-01
             1.1476E+00
 PARAMETER:  8.0733E-02  1.9505E-01 -6.3885E-02  7.7107E-02  6.7348E-02  4.2633E-02  3.5796E-01 -1.3191E+00  6.5865E-02 -3.3621E-02
             2.3771E-01
 GRADIENT:   3.3521E-01 -5.1052E-01 -4.7704E-01 -3.7026E-01  4.8584E-01  1.4800E-01 -1.1756E-01  7.5777E-02  2.3334E-01  4.4476E-02
             3.1523E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2057.74919620311        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1008
 NPARAMETR:  9.8093E-01  1.0957E+00  8.2715E-01  9.7789E-01  9.5283E-01  9.4417E-01  1.3031E+00  1.4006E-01  9.6073E-01  8.6007E-01
             1.1478E+00
 PARAMETER:  8.0741E-02  1.9138E-01 -8.9772E-02  7.7640E-02  5.1677E-02  4.2556E-02  3.6477E-01 -1.8657E+00  5.9941E-02 -5.0737E-02
             2.3785E-01
 GRADIENT:   1.9568E-02 -2.1798E-01 -9.3236E-02 -2.0507E-01 -1.4716E-01  6.0983E-02 -4.7839E-02  1.1894E-02  2.0841E-03  6.6217E-02
             1.5647E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2057.75222702428        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1183
 NPARAMETR:  9.8106E-01  1.1114E+00  8.1365E-01  9.6776E-01  9.5290E-01  9.4423E-01  1.2906E+00  4.7284E-02  9.6631E-01  8.5653E-01
             1.1477E+00
 PARAMETER:  8.0874E-02  2.0558E-01 -1.0622E-01  6.7225E-02  5.1752E-02  4.2615E-02  3.5508E-01 -2.9516E+00  6.5727E-02 -5.4870E-02
             2.3773E-01
 GRADIENT:  -9.5887E-02  4.1220E-02 -2.2822E-01  1.3051E-01  5.9671E-02  3.7544E-03  3.7184E-02  1.4366E-03  2.0158E-02  7.2342E-02
             7.7186E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2057.75274594705        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1359
 NPARAMETR:  9.8109E-01  1.1105E+00  8.1385E-01  9.6817E-01  9.5262E-01  9.4419E-01  1.2912E+00  1.1371E-02  9.6601E-01  8.5629E-01
             1.1476E+00
 PARAMETER:  8.0914E-02  2.0481E-01 -1.0598E-01  6.7650E-02  5.1462E-02  4.2572E-02  3.5556E-01 -4.3766E+00  6.5422E-02 -5.5147E-02
             2.3769E-01
 GRADIENT:   1.9824E-02 -1.4455E-02 -1.3663E-02 -5.7115E-03  2.2230E-02 -1.0187E-02  1.0981E-03  6.6975E-05 -1.8874E-04 -1.3964E-04
            -1.1890E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -2057.75274821317        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1453
 NPARAMETR:  9.8109E-01  1.1105E+00  8.1383E-01  9.6817E-01  9.5260E-01  9.4420E-01  1.2912E+00  1.0000E-02  9.6601E-01  8.5628E-01
             1.1476E+00
 PARAMETER:  8.0911E-02  2.0481E-01 -1.0601E-01  6.7653E-02  5.1442E-02  4.2580E-02  3.5557E-01 -4.5294E+00  6.5418E-02 -5.5162E-02
             2.3770E-01
 GRADIENT:   1.2337E-02 -1.0654E-02 -1.2117E-02 -3.5288E-03  1.6540E-02 -7.1597E-03  1.0350E-03  7.2812E-05  6.4976E-05  7.9699E-04
            -6.7463E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1453
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.2452E-04 -1.2757E-03 -3.9776E-04 -1.4135E-03 -1.5270E-02
 SE:             2.9843E-02  2.2967E-02  1.7149E-04  2.3786E-02  2.1628E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7529E-01  9.5570E-01  2.0367E-02  9.5261E-01  4.8017E-01

 ETASHRINKSD(%)  2.2006E-02  2.3057E+01  9.9426E+01  2.0314E+01  2.7542E+01
 ETASHRINKVR(%)  4.4007E-02  4.0798E+01  9.9997E+01  3.6501E+01  4.7499E+01
 EBVSHRINKSD(%)  4.7922E-01  2.2703E+01  9.9489E+01  2.0739E+01  2.6800E+01
 EBVSHRINKVR(%)  9.5614E-01  4.0252E+01  9.9997E+01  3.7176E+01  4.6418E+01
 RELATIVEINF(%)  9.8692E+01  5.1420E+00  3.5574E-04  5.6879E+00  6.9897E+00
 EPSSHRINKSD(%)  3.2250E+01
 EPSSHRINKVR(%)  5.4099E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2057.7527482131654     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1138.8142150084927     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2057.753       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.81E-01  1.11E+00  8.14E-01  9.68E-01  9.53E-01  9.44E-01  1.29E+00  1.00E-02  9.66E-01  8.56E-01  1.15E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       49.023
Stop Time:
Sun Oct 24 00:04:16 CDT 2021
