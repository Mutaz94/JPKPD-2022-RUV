Sun Oct 24 03:28:28 CDT 2021
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
$DATA ../../../../data/SD4/SL3/dat23.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1637.48358137561        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1303E+02  2.6321E+01 -5.2455E+01  1.2510E+02  3.5414E+01  6.5720E+01  6.4043E+00  1.5207E+01  2.7233E+01  8.7174E+00
            -1.3873E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1661.07840700069        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  1.0406E+00  1.0568E+00  1.1396E+00  9.5602E-01  1.0680E+00  9.7597E-01  9.8237E-01  9.4278E-01  9.2708E-01  9.5658E-01
             1.3581E+00
 PARAMETER:  1.3976E-01  1.5528E-01  2.3066E-01  5.5027E-02  1.6579E-01  7.5677E-02  8.2217E-02  4.1078E-02  2.4287E-02  5.5607E-02
             4.0610E-01
 GRADIENT:   3.6677E+01 -1.6746E+01 -2.7627E+01  5.7011E+00  1.1317E+01  2.2956E+00  1.9239E+00  8.6959E+00  6.8983E+00  2.7796E+00
             2.4279E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1668.16596252414        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  1.0203E+00  1.0601E+00  1.9211E+00  9.6266E-01  1.3074E+00  9.4291E-01  7.7536E-01  6.7491E-01  1.0089E+00  1.2371E+00
             1.2947E+00
 PARAMETER:  1.2006E-01  1.5839E-01  7.5290E-01  6.1941E-02  3.6807E-01  4.1214E-02 -1.5443E-01 -2.9318E-01  1.0888E-01  3.1275E-01
             3.5830E-01
 GRADIENT:  -5.3497E+00 -1.4883E+01 -5.4988E+00 -1.4257E+01  2.0842E+01 -9.7530E+00  5.2661E+00  7.5199E-01  1.7556E+00  6.3951E+00
             5.2047E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1669.90803029108        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  1.0217E+00  1.0490E+00  1.7727E+00  9.6989E-01  1.2410E+00  9.7306E-01  5.4153E-01  4.4880E-01  1.0636E+00  1.1470E+00
             1.2888E+00
 PARAMETER:  1.2146E-01  1.4781E-01  6.7250E-01  6.9429E-02  3.1591E-01  7.2694E-02 -5.1335E-01 -7.0118E-01  1.6163E-01  2.3718E-01
             3.5371E-01
 GRADIENT:  -1.4331E+00 -2.1278E+00 -3.0142E+00  1.5492E+00  4.6961E+00  2.5237E+00  6.1392E-01  3.9901E-01  6.8212E-01 -1.8377E+00
             2.8276E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1670.22483733071        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      701
 NPARAMETR:  1.0237E+00  1.1878E+00  1.8034E+00  8.7849E-01  1.2840E+00  9.6570E-01  4.5663E-01  2.6288E-01  1.1668E+00  1.1802E+00
             1.2920E+00
 PARAMETER:  1.2343E-01  2.7209E-01  6.8965E-01 -2.9552E-02  3.4996E-01  6.5094E-02 -6.8388E-01 -1.2361E+00  2.5425E-01  2.6572E-01
             3.5621E-01
 GRADIENT:   1.1813E+00  5.7447E+00  1.2263E+00  4.7577E+00 -3.6737E+00 -6.4336E-01 -4.7472E-01  8.7220E-02 -1.7049E+00  1.4068E-01
             1.5965E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1670.38047576793        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      878
 NPARAMETR:  1.0241E+00  1.3550E+00  1.5962E+00  7.6695E-01  1.3087E+00  9.7103E-01  4.8263E-01  6.6334E-02  1.3078E+00  1.1829E+00
             1.2882E+00
 PARAMETER:  1.2382E-01  4.0379E-01  5.6765E-01 -1.6533E-01  3.6906E-01  7.0597E-02 -6.2850E-01 -2.6131E+00  3.6837E-01  2.6796E-01
             3.5327E-01
 GRADIENT:  -3.0021E-02  5.4360E+00 -4.7134E-01  6.7952E+00 -7.2272E-01  9.1238E-01 -4.3962E-02  7.5995E-03  8.9397E-01  7.0621E-01
            -5.5624E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1670.40659999022        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  1.0247E+00  1.4491E+00  1.5130E+00  7.0279E-01  1.3282E+00  9.6969E-01  4.8419E-01  2.3864E-02  1.3918E+00  1.1835E+00
             1.2899E+00
 PARAMETER:  1.2436E-01  4.7091E-01  5.1408E-01 -2.5270E-01  3.8383E-01  6.9222E-02 -6.2528E-01 -3.6354E+00  4.3063E-01  2.6851E-01
             3.5455E-01
 GRADIENT:   9.2334E-02  5.3684E+00 -8.3804E-01  5.5907E+00  5.4640E-01  1.7169E-01 -3.7958E-01  1.0452E-03 -2.6401E-01  4.9650E-01
            -4.1890E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1670.40741773626        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1232
 NPARAMETR:  1.0248E+00  1.4781E+00  1.4903E+00  6.8305E-01  1.3343E+00  9.6960E-01  4.8648E-01  1.6545E-02  1.4203E+00  1.1834E+00
             1.2905E+00
 PARAMETER:  1.2448E-01  4.9073E-01  4.9895E-01 -2.8119E-01  3.8838E-01  6.9131E-02 -6.2056E-01 -4.0016E+00  4.5089E-01  2.6841E-01
             3.5507E-01
 GRADIENT:   6.8592E-02  5.2859E+00 -7.4595E-01  5.1994E+00  4.9933E-01  8.1167E-02 -3.7292E-01  5.0480E-04 -3.8840E-01  4.0080E-01
            -3.8573E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1670.43107211868        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1415
 NPARAMETR:  1.0254E+00  1.4758E+00  1.5036E+00  6.8004E-01  1.3345E+00  9.6960E-01  4.8873E-01  1.0000E-02  1.4218E+00  1.1800E+00
             1.2916E+00
 PARAMETER:  1.2512E-01  4.8920E-01  5.0786E-01 -2.8561E-01  3.8853E-01  6.9124E-02 -6.1595E-01 -4.5739E+00  4.5190E-01  2.6550E-01
             3.5591E-01
 GRADIENT:   1.8544E+00 -2.6429E+00 -1.2079E-01 -2.9212E-01  1.1757E-01  1.3354E-01  7.5186E-02  0.0000E+00  4.1807E-04 -3.8464E-02
            -4.5108E-02

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1670.43107211868        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1437
 NPARAMETR:  1.0254E+00  1.4758E+00  1.5036E+00  6.8004E-01  1.3345E+00  9.6960E-01  4.8873E-01  1.0000E-02  1.4218E+00  1.1800E+00
             1.2916E+00
 PARAMETER:  1.2512E-01  4.8920E-01  5.0786E-01 -2.8561E-01  3.8853E-01  6.9124E-02 -6.1595E-01 -4.5739E+00  4.5190E-01  2.6550E-01
             3.5591E-01
 GRADIENT:   1.8544E+00 -2.6429E+00 -1.2079E-01 -2.9212E-01  1.1757E-01  1.3354E-01  7.5186E-02  0.0000E+00  4.1807E-04 -3.8464E-02
            -4.5108E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1437
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.5824E-04 -3.1891E-02 -3.6287E-05  8.0879E-03 -3.4845E-02
 SE:             2.9750E-02  1.3606E-02  4.3531E-05  2.5325E-02  2.3133E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8771E-01  1.9084E-02  4.0452E-01  7.4945E-01  1.3200E-01

 ETASHRINKSD(%)  3.3247E-01  5.4418E+01  9.9854E+01  1.5157E+01  2.2501E+01
 ETASHRINKVR(%)  6.6383E-01  7.9222E+01  1.0000E+02  2.8016E+01  3.9939E+01
 EBVSHRINKSD(%)  6.9270E-01  5.5012E+01  9.9825E+01  1.4665E+01  2.0824E+01
 EBVSHRINKVR(%)  1.3806E+00  7.9761E+01  1.0000E+02  2.7179E+01  3.7312E+01
 RELATIVEINF(%)  9.8437E+01  1.4842E+00  1.2139E-04  5.7959E+00  2.2601E+01
 EPSSHRINKSD(%)  3.8655E+01
 EPSSHRINKVR(%)  6.2368E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1670.4310721186755     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -935.28024555493732     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1670.431       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.48E+00  1.50E+00  6.80E-01  1.33E+00  9.70E-01  4.89E-01  1.00E-02  1.42E+00  1.18E+00  1.29E+00
 


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
 #CPUT: Total CPU Time in Seconds,       40.311
Stop Time:
Sun Oct 24 03:28:37 CDT 2021
