Sun Oct 24 02:21:55 CDT 2021
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
$DATA ../../../../data/SD4/A2/dat94.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m94.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1010.23553535338        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2756E+02  1.1329E+01  4.1981E+01 -4.3477E+01  1.1672E+02  4.5726E+01 -3.7635E+01 -1.8160E+00 -7.1464E+01 -8.4438E+01
            -1.0911E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1367.50740572498        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1684E+00  9.4652E-01  9.6347E-01  1.1242E+00  8.5877E-01  1.1692E+00  1.0419E+00  9.2839E-01  1.2037E+00  1.0687E+00
             2.5907E+00
 PARAMETER:  2.5565E-01  4.5039E-02  6.2781E-02  2.1704E-01 -5.2257E-02  2.5630E-01  1.4108E-01  2.5695E-02  2.8544E-01  1.6642E-01
             1.0519E+00
 GRADIENT:   3.8895E+02  2.1340E+01  1.7273E+01  1.5974E+01 -2.9272E+01  4.3184E+01  5.7706E+00  8.7201E+00  1.2475E+01  1.2263E+01
             3.6545E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1386.46846709331        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      181
 NPARAMETR:  1.1397E+00  5.5847E-01  4.9131E-01  1.4019E+00  4.8117E-01  1.1148E+00  8.3485E-01  2.0088E-01  1.1434E+00  6.2478E-01
             2.4715E+00
 PARAMETER:  2.3078E-01 -4.8255E-01 -6.1068E-01  4.3784E-01 -6.3154E-01  2.0865E-01 -8.0503E-02 -1.5050E+00  2.3400E-01 -3.7035E-01
             1.0048E+00
 GRADIENT:   1.9995E+02  2.7897E+01 -4.0770E+01  1.4639E+02  4.2543E+01  2.3922E+01 -2.9376E+00  5.6324E-01  4.6846E+00 -1.9162E+00
             8.3948E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1414.35177551723        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      359
 NPARAMETR:  1.0000E+00  4.6601E-01  4.1500E-01  1.3095E+00  4.0290E-01  9.6549E-01  1.2758E+00  1.2627E-02  9.9372E-01  5.9657E-01
             2.2781E+00
 PARAMETER:  1.0002E-01 -6.6354E-01 -7.7949E-01  3.6964E-01 -8.0907E-01  6.4875E-02  3.4357E-01 -4.2719E+00  9.3700E-02 -4.1657E-01
             9.2335E-01
 GRADIENT:  -2.9206E+01  2.0218E+01 -2.3644E+00  6.1851E+01 -4.5930E+00 -1.3386E+00 -3.6642E+00  1.5108E-03 -1.4057E+01 -8.1491E-01
            -7.9488E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1418.39233330084        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  1.0143E+00  2.9767E-01  3.9954E-01  1.3233E+00  3.5621E-01  9.5176E-01  1.6834E+00  1.0000E-02  1.0165E+00  6.6396E-01
             2.2573E+00
 PARAMETER:  1.1422E-01 -1.1118E+00 -8.1744E-01  3.8012E-01 -9.3223E-01  5.0554E-02  6.2083E-01 -6.6212E+00  1.1635E-01 -3.0953E-01
             9.1417E-01
 GRADIENT:   1.1544E+01  8.5481E+00  2.3195E+01 -1.0040E+01 -3.4570E+01 -4.5566E+00  1.4851E-01  0.0000E+00  1.4525E+00  5.6628E+00
            -2.6980E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1420.57088329619        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  9.9933E-01  1.5359E-01  4.5175E-01  1.4191E+00  3.7557E-01  9.5434E-01  2.5188E+00  1.0000E-02  9.7765E-01  6.8674E-01
             2.2996E+00
 PARAMETER:  9.9326E-02 -1.7735E+00 -6.9463E-01  4.5004E-01 -8.7932E-01  5.3262E-02  1.0238E+00 -8.1243E+00  7.7395E-02 -2.7579E-01
             9.3274E-01
 GRADIENT:  -5.5419E+00  7.0322E-01 -4.1844E+00 -3.1142E+00  5.9704E+00 -1.0754E+00 -8.4958E-03  0.0000E+00  2.4475E+00  3.5314E+00
             4.5372E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1421.02113273820        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      885
 NPARAMETR:  9.9936E-01  1.1386E-01  4.4435E-01  1.4370E+00  3.6525E-01  9.5681E-01  3.1181E+00  1.0000E-02  9.5976E-01  6.5257E-01
             2.2932E+00
 PARAMETER:  9.9362E-02 -2.0728E+00 -7.1115E-01  4.6253E-01 -9.0718E-01  5.5850E-02  1.2372E+00 -9.1118E+00  5.8930E-02 -3.2684E-01
             9.2996E-01
 GRADIENT:  -7.4142E-01  9.4733E+00 -4.2694E+00  1.0929E+00 -4.1975E-01 -7.9204E-01  1.2019E+01  0.0000E+00 -5.7979E+00 -5.9815E+00
            -3.2392E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1421.14583334567        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1053
 NPARAMETR:  9.9774E-01  9.9453E-02  4.4347E-01  1.4465E+00  3.6011E-01  9.5782E-01  3.3082E+00  1.0000E-02  9.6082E-01  6.5079E-01
             2.3103E+00
 PARAMETER:  9.8685E-02 -2.1869E+00 -7.1870E-01  4.6546E-01 -9.1607E-01  5.5944E-02  1.3077E+00 -9.5196E+00  5.9051E-02 -3.3289E-01
             9.2893E-01
 GRADIENT:   2.1901E+02  9.4396E+00 -5.2564E+01 -8.1625E+01  3.5438E+01 -2.2144E+02  2.8627E+01  0.0000E+00 -2.2537E+02 -6.8592E+01
            -4.6656E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1053
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.7595E-04  4.5701E-03 -2.9336E-05 -1.4849E-02 -1.3172E-03
 SE:             2.9176E-02  7.8110E-03  2.4931E-04  2.7592E-02  2.1510E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7605E-01  5.5849E-01  9.0633E-01  5.9046E-01  9.5117E-01

 ETASHRINKSD(%)  2.2558E+00  7.3832E+01  9.9165E+01  7.5632E+00  2.7938E+01
 ETASHRINKVR(%)  4.4607E+00  9.3152E+01  9.9993E+01  1.4554E+01  4.8070E+01
 EBVSHRINKSD(%)  2.1586E+00  8.1710E+01  9.9144E+01  7.3081E+00  2.6967E+01
 EBVSHRINKVR(%)  4.2705E+00  9.6655E+01  9.9993E+01  1.4082E+01  4.6662E+01
 RELATIVEINF(%)  9.2463E+01  9.3186E-01  2.8808E-04  2.9312E+01  1.9360E+00
 EPSSHRINKSD(%)  3.4622E+01
 EPSSHRINKVR(%)  5.7257E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1421.1458333456742     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -685.99500678193601     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1421.146       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  1.02E-01  4.41E-01  1.44E+00  3.62E-01  9.57E-01  3.35E+00  1.00E-02  9.60E-01  6.49E-01  2.29E+00
 


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
 #CPUT: Total CPU Time in Seconds,       31.760
Stop Time:
Sun Oct 24 02:22:03 CDT 2021
