Sat Oct 23 21:05:29 CDT 2021
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
$DATA ../../../../data/SD3/B/dat7.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1301.49329637210        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6924E+02  7.8311E+00  1.4106E+02  1.3149E+01  8.3579E+01  6.9567E+00  1.9250E+01 -6.4764E+02 -9.6770E+01 -6.0342E+00
            -8.0573E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1835.61476703589        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  8.2891E-01  1.1197E+00  8.1155E-01  8.7581E-01  9.8204E-01  1.1081E+00  8.9902E-01  4.1051E+00  9.9264E-01  8.9925E-01
             1.8794E+00
 PARAMETER: -8.7644E-02  2.1310E-01 -1.0881E-01 -3.2606E-02  8.1881E-02  2.0265E-01 -6.4523E-03  1.5122E+00  9.2611E-02 -6.1946E-03
             7.3094E-01
 GRADIENT:  -1.8904E+02 -1.2178E+02 -2.0084E+01 -6.1952E+01 -4.8972E+01  9.0741E+00  2.4802E+00  3.4301E+01  2.0647E+01  2.3842E+01
             2.9374E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1922.81369345213        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  8.3323E-01  1.2348E+00  2.5299E+00  9.2951E-01  1.3690E+00  9.8003E-01  1.8523E-01  5.4814E+00  1.4506E+00  9.3722E-01
             1.7406E+00
 PARAMETER: -8.2451E-02  3.1091E-01  1.0282E+00  2.6898E-02  4.1408E-01  7.9824E-02 -1.5861E+00  1.8014E+00  4.7200E-01  3.5168E-02
             6.5422E-01
 GRADIENT:  -2.4002E+02  6.2147E+00 -7.0675E+00  7.9349E+01  4.5471E+01 -8.3710E+01 -2.5785E+00  8.4006E+01  6.3491E+01 -1.0646E+00
             2.4785E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1928.83283260498        NO. OF FUNC. EVALS.: 102
 CUMULATIVE NO. OF FUNC. EVALS.:      259
 NPARAMETR:  8.3514E-01  1.2394E+00  2.4261E+00  9.2611E-01  1.3434E+00  9.8262E-01  1.7924E-01  5.4077E+00  1.4438E+00  8.9134E-01
             1.7130E+00
 PARAMETER: -8.0153E-02  3.1463E-01  9.8627E-01  2.3235E-02  3.9522E-01  8.2468E-02 -1.6190E+00  1.7878E+00  4.6728E-01 -1.5030E-02
             6.3827E-01
 GRADIENT:  -3.8936E+02 -7.4685E+01 -1.0344E+01  6.2904E+01  3.4436E+01 -1.2134E+02 -4.3676E+00  1.1874E+01  4.9318E+01 -2.7586E+00
             2.3628E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1982.90030796738        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  8.6398E-01  1.2119E+00  2.5235E+00  8.9382E-01  1.3892E+00  1.0196E+00  5.9633E-01  5.5561E+00  1.2103E+00  1.2667E+00
             1.3078E+00
 PARAMETER: -4.6207E-02  2.9217E-01  1.0257E+00 -1.2249E-02  4.2872E-01  1.1945E-01 -4.1696E-01  1.8149E+00  2.9088E-01  3.3645E-01
             3.6832E-01
 GRADIENT:  -2.8480E+02 -9.4468E+01 -5.8306E+00 -7.1278E+00  2.2624E+01 -7.7956E+01  1.3120E+00 -2.8903E+01  3.2872E+01  1.5394E+01
             1.4395E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2040.76340993391        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      627             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2083E-01  1.3952E+00  2.5725E+00  7.7889E-01  1.3756E+00  1.0223E+00  4.6196E-01  5.8339E+00  1.3189E+00  1.1429E+00
             1.1694E+00
 PARAMETER:  1.7522E-02  4.3305E-01  1.0449E+00 -1.4989E-01  4.1888E-01  1.2205E-01 -6.7227E-01  1.8637E+00  3.7683E-01  2.3355E-01
             2.5647E-01
 GRADIENT:   2.0010E+02  2.8379E+02  4.2820E-01  7.0305E+01  4.2534E+01  3.3677E+01  1.2784E+01  6.4609E+01  3.1061E+01  5.9200E+00
             7.7513E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2041.22544227147        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      772
 NPARAMETR:  9.2056E-01  1.3952E+00  2.5725E+00  7.7888E-01  1.3756E+00  1.0299E+00  4.5281E-01  5.8345E+00  1.3111E+00  1.1429E+00
             1.1694E+00
 PARAMETER:  1.7229E-02  4.3302E-01  1.0449E+00 -1.4990E-01  4.1887E-01  1.2945E-01 -6.9229E-01  1.8638E+00  3.7086E-01  2.3355E-01
             2.5647E-01
 GRADIENT:  -1.4946E+02 -3.2579E+01 -3.8609E+00  2.0698E+01  2.7268E+01 -4.8538E+01  2.2838E+00  5.2917E+00  1.2369E+01  4.5925E+00
             7.5010E+01

0ITERATION NO.:   32    OBJECTIVE VALUE:  -2041.30803050083        NO. OF FUNC. EVALS.:  62
 CUMULATIVE NO. OF FUNC. EVALS.:      834
 NPARAMETR:  9.2057E-01  1.3952E+00  2.5725E+00  7.7888E-01  1.3756E+00  1.0316E+00  4.5279E-01  5.8343E+00  1.3111E+00  1.1429E+00
             1.1694E+00
 PARAMETER:  1.7234E-02  4.3302E-01  1.0449E+00 -1.4990E-01  4.1887E-01  1.3112E-01 -6.9232E-01  1.8638E+00  3.7084E-01  2.3355E-01
             2.5647E-01
 GRADIENT:   9.3652E+03 -1.1228E+02  9.0930E+02  1.0671E+01  2.9036E+01  7.2092E+03 -1.3735E+03 -6.1908E+02  1.2957E+03  4.0953E+00
             7.3969E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:      834
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.0658E-02 -3.3685E-02 -5.3450E-02 -1.9807E-03 -8.5193E-02
 SE:             3.2581E-02  1.2131E-02  1.7260E-02  2.5530E-02  1.8908E-02
 N:                     100         100         100         100         100

 P VAL.:         3.0107E-02  5.4883E-03  1.9568E-03  9.3816E-01  6.6249E-06

 ETASHRINKSD(%)  1.0000E-10  5.9361E+01  4.2176E+01  1.4471E+01  3.6655E+01
 ETASHRINKVR(%)  1.0000E-10  8.3485E+01  6.6564E+01  2.6849E+01  5.9874E+01
 EBVSHRINKSD(%)  4.2594E-01  5.9678E+01  4.0485E+01  1.0813E+01  2.8327E+01
 EBVSHRINKVR(%)  8.5006E-01  8.3741E+01  6.4579E+01  2.0457E+01  4.8630E+01
 RELATIVEINF(%)  9.9079E+01  2.2196E+00  2.3082E+01  1.0996E+01  3.6583E+01
 EPSSHRINKSD(%)  4.0591E+01
 EPSSHRINKVR(%)  6.4706E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2041.3080305008264     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1122.3694972961537     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2041.308       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.21E-01  1.40E+00  2.57E+00  7.79E-01  1.38E+00  1.03E+00  4.53E-01  5.83E+00  1.31E+00  1.14E+00  1.17E+00
 


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
 #CPUT: Total CPU Time in Seconds,       68.965
Stop Time:
Sat Oct 23 21:05:41 CDT 2021
