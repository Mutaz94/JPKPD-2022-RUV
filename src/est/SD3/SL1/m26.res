Sat Oct 23 23:24:07 CDT 2021
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
$DATA ../../../../data/SD3/SL1/dat26.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1599.65393021195        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2321E+02 -1.2980E+01  5.9745E+01  9.8821E+00  1.2058E+01  4.4239E+01 -6.5309E+00 -2.0126E+02 -1.2350E+01 -1.6923E+01
            -7.5879E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2089.57865625624        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:      145
 NPARAMETR:  9.9331E-01  1.0298E+00  1.0339E+00  9.9661E-01  1.0194E+00  1.0405E+00  1.0960E+00  5.8515E-01  8.6711E-01  1.1081E+00
             1.3453E+00
 PARAMETER:  9.3286E-02  1.2935E-01  1.3331E-01  9.6601E-02  1.1925E-01  1.3970E-01  1.9168E-01 -4.3589E-01 -4.2591E-02  2.0268E-01
             3.9661E-01
 GRADIENT:   2.3673E+02  2.0938E+01  8.8666E+00  2.0195E+01 -2.2495E+01  6.7853E+01  1.1349E+00  1.8608E+00 -2.7451E+00 -3.4431E+00
             1.9601E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2092.05739789784        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.9651E-01  1.0383E+00  1.0648E+00  1.0040E+00  1.0664E+00  9.5250E-01  1.1064E+00  4.5697E-01  8.9769E-01  1.1969E+00
             1.3453E+00
 PARAMETER:  9.6506E-02  1.3761E-01  1.6282E-01  1.0400E-01  1.6432E-01  5.1338E-02  2.0110E-01 -6.8313E-01 -7.9288E-03  2.7972E-01
             3.9665E-01
 GRADIENT:   6.2838E+00 -2.6435E+00 -4.5302E+00  4.9977E-01 -3.2568E+00 -2.9642E+00  6.8288E-01  6.5062E-01  3.3415E-01  1.5428E-01
             1.9481E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2092.42927327061        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      477
 NPARAMETR:  9.9265E-01  1.0917E+00  1.1437E+00  9.7561E-01  1.1321E+00  9.7061E-01  1.0163E+00  4.7063E-01  9.3044E-01  1.2747E+00
             1.3451E+00
 PARAMETER:  9.2620E-02  1.8775E-01  2.3424E-01  7.5312E-02  2.2410E-01  7.0168E-02  1.1621E-01 -6.5368E-01  2.7907E-02  3.4268E-01
             3.9649E-01
 GRADIENT:  -2.7444E+00  4.1721E-01  2.2978E+00 -1.5795E-01 -4.9269E+00  4.3018E+00 -1.2738E+00 -2.6075E-01 -1.4555E+00  6.6159E-01
             1.9355E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2092.67035595640        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      654
 NPARAMETR:  9.9403E-01  1.2635E+00  1.1428E+00  8.7017E-01  1.2273E+00  9.6448E-01  9.2289E-01  6.5223E-01  1.0171E+00  1.3320E+00
             1.3442E+00
 PARAMETER:  9.4011E-02  3.3392E-01  2.3344E-01 -3.9070E-02  3.0486E-01  6.3832E-02  1.9758E-02 -3.2735E-01  1.1697E-01  3.8665E-01
             3.9579E-01
 GRADIENT:  -1.6488E+00  3.2656E+00  1.6262E+00  2.0460E+00 -1.5845E+00  1.5123E+00 -3.6588E-01  2.8019E-01 -5.6075E-01  5.6362E-01
             1.9271E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2093.20727055275        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      830
 NPARAMETR:  9.9765E-01  1.4761E+00  8.4518E-01  7.2641E-01  1.2070E+00  9.6299E-01  8.6952E-01  2.6001E-01  1.1149E+00  1.2604E+00
             1.3418E+00
 PARAMETER:  9.7643E-02  4.8941E-01 -6.8203E-02 -2.1964E-01  2.8815E-01  6.2290E-02 -3.9814E-02 -1.2470E+00  2.0876E-01  3.3147E-01
             3.9403E-01
 GRADIENT:   2.1071E+00  2.3497E+00 -2.1926E+00  3.0828E+00 -2.9766E+00  1.2143E-01 -1.0675E-02  2.3776E-01  2.7951E-01  5.0108E-01
             1.9336E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2093.54195475844        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1006
 NPARAMETR:  9.9731E-01  1.7301E+00  7.3279E-01  5.6405E-01  1.3381E+00  9.6535E-01  7.8292E-01  4.8303E-02  1.3254E+00  1.3506E+00
             1.3383E+00
 PARAMETER:  9.7303E-02  6.4820E-01 -2.1090E-01 -4.7261E-01  3.9128E-01  6.4738E-02 -1.4473E-01 -2.9303E+00  3.8172E-01  4.0057E-01
             3.9142E-01
 GRADIENT:  -3.5582E-01  5.2629E+00 -1.6368E+00  4.3724E+00  2.9891E+00  7.4170E-01 -2.4412E-01  6.8368E-03 -3.9437E-01  1.9287E-01
             1.8996E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2124.27042868925        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1189
 NPARAMETR:  1.0005E+00  1.7862E+00  6.8157E-01  5.3745E-01  1.3611E+00  9.7026E-01  7.8132E-01  1.9665E-02  1.4014E+00  1.3881E+00
             9.5696E-01
 PARAMETER:  1.0045E-01  6.8010E-01 -2.8335E-01 -5.2093E-01  4.0829E-01  6.9804E-02 -1.4677E-01 -3.8289E+00  4.3748E-01  4.2793E-01
             5.6004E-02
 GRADIENT:   1.5122E+01  4.1971E+01  6.3853E+00  1.3917E+01  5.4078E+00  7.0761E-01 -3.0406E+00  1.6196E-03  2.2980E+00 -1.1073E+01
            -2.6562E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2125.26174493528        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1364
 NPARAMETR:  9.9575E-01  1.8014E+00  6.3014E-01  5.0937E-01  1.3565E+00  9.7334E-01  7.8936E-01  1.0000E-02  1.3741E+00  1.4514E+00
             9.8135E-01
 PARAMETER:  9.5737E-02  6.8856E-01 -3.6181E-01 -5.7459E-01  4.0493E-01  7.2982E-02 -1.3653E-01 -4.7082E+00  4.1778E-01  4.7255E-01
             8.1174E-02
 GRADIENT:   3.3497E+00  4.5465E+00  1.6816E-01 -4.9939E-01 -8.0427E-01  2.1615E+00 -3.6603E-01  0.0000E+00  5.7170E-01  5.5289E-01
            -1.4625E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2125.32542131194        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1528
 NPARAMETR:  9.9464E-01  1.7946E+00  6.3130E-01  5.1237E-01  1.3547E+00  9.6820E-01  7.9239E-01  1.0000E-02  1.3627E+00  1.4461E+00
             9.8271E-01
 PARAMETER:  9.4622E-02  6.8479E-01 -3.5998E-01 -5.6870E-01  4.0358E-01  6.7681E-02 -1.3270E-01 -4.7199E+00  4.0944E-01  4.6887E-01
             8.2557E-02
 GRADIENT:   7.1761E-01  1.6396E+00 -2.3866E-01 -1.3930E+00  5.7104E-01  1.0388E-01 -1.9315E-01  0.0000E+00  2.1340E-01  1.4766E-01
            -1.1875E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1528
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0979E-03 -3.4834E-02 -2.9717E-04  3.1609E-02 -4.0876E-02
 SE:             2.9889E-02  2.3669E-02  8.5023E-05  2.1850E-02  2.3304E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7070E-01  1.4109E-01  4.7387E-04  1.4800E-01  7.9418E-02

 ETASHRINKSD(%)  1.0000E-10  2.0707E+01  9.9715E+01  2.6799E+01  2.1929E+01
 ETASHRINKVR(%)  1.0000E-10  3.7126E+01  9.9999E+01  4.6417E+01  3.9050E+01
 EBVSHRINKSD(%)  3.4518E-01  1.9201E+01  9.9742E+01  3.1650E+01  1.6414E+01
 EBVSHRINKVR(%)  6.8918E-01  3.4716E+01  9.9999E+01  5.3283E+01  3.0134E+01
 RELATIVEINF(%)  9.9148E+01  4.1276E+00  9.7515E-05  2.7010E+00  2.4008E+01
 EPSSHRINKSD(%)  3.3293E+01
 EPSSHRINKVR(%)  5.5502E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2125.3254213119385     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1206.3868881072658     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.80
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2125.325       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.95E-01  1.79E+00  6.31E-01  5.12E-01  1.35E+00  9.68E-01  7.92E-01  1.00E-02  1.36E+00  1.45E+00  9.83E-01
 


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
 #CPUT: Total CPU Time in Seconds,      122.818
Stop Time:
Sat Oct 23 23:24:26 CDT 2021
