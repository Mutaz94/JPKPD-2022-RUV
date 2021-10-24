Sun Oct 24 00:06:28 CDT 2021
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
$DATA ../../../../data/SD3/SL2/dat97.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2060.03563928812        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.5051E+02 -4.1776E+01 -3.5010E+01  2.3220E+01  4.7172E+01  3.4464E+01 -7.3895E+00  9.0310E+00  2.9472E+01  9.2596E+00
            -2.9599E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2068.74494895278        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.3706E-01  1.1769E+00  1.1183E+00  9.7179E-01  1.0636E+00  1.0945E+00  1.0826E+00  9.4269E-01  8.5016E-01  9.3089E-01
             1.1263E+00
 PARAMETER:  3.4993E-02  2.6288E-01  2.1179E-01  7.1381E-02  1.6170E-01  1.9030E-01  1.7938E-01  4.0987E-02 -6.2333E-02  2.8389E-02
             2.1891E-01
 GRADIENT:   5.7166E+00  4.7726E+01  2.0502E+01  3.7339E+01 -2.2327E+01  2.1515E+01 -6.3098E+00 -1.8265E+00 -4.7240E+00 -1.1740E+01
             5.4639E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2070.24601558418        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.3776E-01  1.0934E+00  1.1569E+00  1.0129E+00  1.0472E+00  1.0832E+00  1.2409E+00  8.5803E-01  7.2103E-01  1.0015E+00
             1.0975E+00
 PARAMETER:  3.5739E-02  1.8933E-01  2.4575E-01  1.1281E-01  1.4615E-01  1.7988E-01  3.1587E-01 -5.3117E-02 -2.2708E-01  1.0147E-01
             1.9301E-01
 GRADIENT:   9.1210E+00  3.8134E+01  2.5309E+01  1.6530E+01 -3.2837E+01  1.7955E+01 -2.1721E+00 -2.9781E+00 -1.4513E+01 -4.3038E+00
             3.5795E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2074.33932504773        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.3333E-01  1.0760E+00  9.1661E-01  9.8910E-01  9.6451E-01  1.0316E+00  1.2039E+00  5.4138E-01  8.2957E-01  9.3428E-01
             1.0387E+00
 PARAMETER:  3.1003E-02  1.7325E-01  1.2925E-02  8.9039E-02  6.3870E-02  1.3107E-01  2.8556E-01 -5.1364E-01 -8.6848E-02  3.2025E-02
             1.3796E-01
 GRADIENT:  -6.0347E-01 -1.9614E+00 -2.3284E+00  7.5503E-01  2.7501E+00 -1.0230E+00  1.2976E-01  3.5745E-01 -2.4079E-01  1.8706E-01
            -1.8760E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2074.43974380866        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.3461E-01  1.2084E+00  8.3389E-01  9.0512E-01  9.8856E-01  1.0355E+00  1.0940E+00  3.3626E-01  8.8908E-01  9.5154E-01
             1.0415E+00
 PARAMETER:  3.2377E-02  2.8932E-01 -8.1656E-02  3.0950E-04  8.8495E-02  1.3490E-01  1.8980E-01 -9.8987E-01 -1.7567E-02  5.0322E-02
             1.4070E-01
 GRADIENT:   9.6267E-02 -8.4415E-02 -7.6913E-01  8.3854E-01  1.1345E+00  1.4001E-01 -7.4940E-03  8.8575E-02  1.0592E-01 -2.1440E-02
             3.7056E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2074.44238156535        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      883
 NPARAMETR:  9.3470E-01  1.2349E+00  8.1709E-01  8.8745E-01  9.9374E-01  1.0352E+00  1.0749E+00  2.7002E-01  9.0208E-01  9.5724E-01
             1.0411E+00
 PARAMETER:  3.2467E-02  3.1101E-01 -1.0201E-01 -1.9405E-02  9.3724E-02  1.3462E-01  1.7219E-01 -1.2093E+00 -3.0570E-03  5.6299E-02
             1.4026E-01
 GRADIENT:  -7.5056E-02 -3.9751E-01 -2.0108E-01 -5.3393E-01  4.0721E-02 -3.6746E-02  3.4011E-02  5.1456E-02  1.3732E-01  2.1465E-01
             3.4775E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2074.45462156484        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1061
 NPARAMETR:  9.3478E-01  1.2337E+00  8.0565E-01  8.8764E-01  9.8736E-01  1.0354E+00  1.0766E+00  1.0266E-01  9.0244E-01  9.5428E-01
             1.0415E+00
 PARAMETER:  3.2557E-02  3.1004E-01 -1.1610E-01 -1.9189E-02  8.7279E-02  1.3480E-01  1.7377E-01 -2.1763E+00 -2.6556E-03  5.3202E-02
             1.4064E-01
 GRADIENT:  -6.1657E-04  5.8672E-02 -5.0124E-02  1.8580E-01  5.7111E-02  1.2114E-02 -2.5290E-03  1.8738E-03  3.0554E-03 -1.4510E-03
             1.8086E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2074.45569163856        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  9.3542E-01  1.2332E+00  8.0537E-01  8.8751E-01  9.8736E-01  1.0367E+00  1.0767E+00  6.8210E-02  9.0267E-01  9.5493E-01
             1.0415E+00
 PARAMETER:  3.3237E-02  3.0962E-01 -1.1645E-01 -1.9334E-02  8.7281E-02  1.3609E-01  1.7394E-01 -2.5852E+00 -2.4000E-03  5.3879E-02
             1.4070E-01
 GRADIENT:   1.4090E+00 -5.4137E-01 -5.6985E-02 -3.2867E-01  2.5970E-01  5.3187E-01 -7.6196E-04  6.3299E-04  2.3311E-02  1.5424E-02
             3.6075E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2074.45577155298        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1426
 NPARAMETR:  9.3542E-01  1.2332E+00  8.0518E-01  8.8766E-01  9.8721E-01  1.0367E+00  1.0769E+00  5.6623E-02  9.0267E-01  9.5505E-01
             1.0415E+00
 PARAMETER:  3.3235E-02  3.0957E-01 -1.1669E-01 -1.9166E-02  8.7126E-02  1.3608E-01  1.7409E-01 -2.7713E+00 -2.4033E-03  5.4006E-02
             1.4070E-01
 GRADIENT:   1.4040E+00 -3.9160E-01 -8.6433E-02 -1.1862E-01  2.6684E-01  5.3036E-01  1.1218E-02  4.4606E-04  3.3248E-02  3.7548E-02
             2.9464E-02

0ITERATION NO.:   42    OBJECTIVE VALUE:  -2074.45578594209        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     1485
 NPARAMETR:  9.3535E-01  1.2332E+00  8.0520E-01  8.8767E-01  9.8719E-01  1.0366E+00  1.0769E+00  4.2753E-02  9.0265E-01  9.5503E-01
             1.0415E+00
 PARAMETER:  3.3167E-02  3.0961E-01 -1.1667E-01 -1.9158E-02  8.7109E-02  1.3595E-01  1.7409E-01 -3.0523E+00 -2.4258E-03  5.3983E-02
             1.4070E-01
 GRADIENT:  -1.4069E-01  2.7952E-02  8.1108E-03 -4.7201E-02  1.9751E-01 -5.3064E-02 -1.3128E-02  2.5242E-06  5.2211E-03 -2.0229E-03
             3.5129E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1485
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.4612E-04 -7.4520E-03 -1.5577E-03  3.7430E-03 -1.8783E-02
 SE:             2.9861E-02  2.3014E-02  6.6131E-04  2.2876E-02  2.3149E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8808E-01  7.4609E-01  1.8497E-02  8.7003E-01  4.1713E-01

 ETASHRINKSD(%)  1.0000E-10  2.2899E+01  9.7785E+01  2.3363E+01  2.2448E+01
 ETASHRINKVR(%)  1.0000E-10  4.0554E+01  9.9951E+01  4.1268E+01  3.9857E+01
 EBVSHRINKSD(%)  3.4909E-01  2.2673E+01  9.8068E+01  2.4387E+01  2.0847E+01
 EBVSHRINKVR(%)  6.9696E-01  4.0205E+01  9.9963E+01  4.2827E+01  3.7348E+01
 RELATIVEINF(%)  9.9028E+01  3.5384E+00  3.9105E-03  3.2693E+00  9.5036E+00
 EPSSHRINKSD(%)  3.2827E+01
 EPSSHRINKVR(%)  5.4878E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2074.4557859420943     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1155.5172527374216     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2074.456       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.35E-01  1.23E+00  8.05E-01  8.88E-01  9.87E-01  1.04E+00  1.08E+00  4.28E-02  9.03E-01  9.55E-01  1.04E+00
 


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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       49.770
Stop Time:
Sun Oct 24 00:06:39 CDT 2021
