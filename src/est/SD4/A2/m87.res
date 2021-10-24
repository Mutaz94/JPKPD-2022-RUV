Sun Oct 24 02:20:57 CDT 2021
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
$DATA ../../../../data/SD4/A2/dat87.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m87.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -456.143988337372        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.0140E+02  1.0228E+01  2.3630E+01 -1.9795E+01  1.4776E+02  3.5384E+01 -3.4664E+01 -5.7845E+00 -7.7419E+01 -7.6014E+01
            -2.2284E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1320.18442683200        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0208E+00  9.6701E-01  1.0385E+00  1.0769E+00  9.4027E-01  8.8437E-01  9.8531E-01  9.2632E-01  1.1769E+00  7.8658E-01
             2.5677E+00
 PARAMETER:  1.2055E-01  6.6452E-02  1.3779E-01  1.7406E-01  3.8413E-02 -2.2883E-02  8.5202E-02  2.3460E-02  2.6286E-01 -1.4006E-01
             1.0430E+00
 GRADIENT:   5.3239E+01 -7.0341E+00  7.8560E-01 -5.2177E+00  4.0369E+01 -3.8260E+01  1.8376E+00  3.1909E+00 -6.3439E-02 -2.2672E+00
            -1.6620E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1337.68500239489        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0087E+00  7.9640E-01  6.2193E-01  1.1418E+00  5.8944E-01  1.0171E+00  8.3249E-01  5.2581E-01  1.1320E+00  2.3404E-01
             2.7149E+00
 PARAMETER:  1.0867E-01 -1.2766E-01 -3.7493E-01  2.3260E-01 -4.2858E-01  1.1696E-01 -8.3332E-02 -5.4282E-01  2.2398E-01 -1.3523E+00
             1.0988E+00
 GRADIENT:   3.7600E+00  7.3126E+01  1.0253E+02  1.1189E+01 -1.4340E+02  1.1527E+01 -8.9246E+00 -3.3258E+00 -4.0553E+00 -2.4224E+00
            -1.2755E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1343.73622310168        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      255
 NPARAMETR:  1.0069E+00  6.8535E-01  3.4050E-01  1.1880E+00  4.1109E-01  1.0166E+00  1.0242E+00  1.2954E-01  1.2351E+00  1.7744E-01
             2.7662E+00
 PARAMETER:  1.0690E-01 -2.7783E-01 -9.7735E-01  2.7226E-01 -7.8893E-01  1.1642E-01  1.2389E-01 -1.9437E+00  3.1112E-01 -1.6291E+00
             1.1175E+00
 GRADIENT:  -7.0472E+01  3.7560E+01  1.8373E+01  7.2014E+01 -5.5379E+01 -3.1577E+00 -4.8131E+00 -2.7466E-01  2.1826E+01 -2.0968E+00
            -7.7643E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1357.59206550410        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      431
 NPARAMETR:  1.0409E+00  5.2977E-01  6.5158E-01  1.3307E+00  5.5933E-01  9.7456E-01  1.1100E+00  4.1483E-02  1.0304E+00  2.1972E-01
             3.2313E+00
 PARAMETER:  1.4012E-01 -5.3532E-01 -3.2835E-01  3.8570E-01 -4.8102E-01  7.4233E-02  2.0436E-01 -3.0825E+00  1.2994E-01 -1.4154E+00
             1.2729E+00
 GRADIENT:  -8.9793E-01  1.3877E+01  2.0304E+01  1.6493E+01 -3.4277E+01 -2.0384E+00 -8.2489E-01  1.6964E-02  7.4938E-01  1.3589E-01
            -8.4321E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1358.36212050490        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  1.0373E+00  4.4180E-01  6.9142E-01  1.3826E+00  5.7093E-01  9.7174E-01  1.4263E+00  1.1431E-02  9.8914E-01  1.7737E-01
             3.2635E+00
 PARAMETER:  1.3662E-01 -7.1690E-01 -2.6901E-01  4.2399E-01 -4.6048E-01  7.1335E-02  4.5505E-01 -4.3714E+00  8.9079E-02 -1.6295E+00
             1.2828E+00
 GRADIENT:  -2.1874E+00  3.6692E+00 -1.2556E+00  1.0547E+01  2.2824E+00 -9.2783E-01  1.2684E+00  1.6526E-03  1.7723E+00  2.5361E-01
            -3.7037E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1358.79612448952        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  1.0353E+00  3.2006E-01  6.6722E-01  1.4313E+00  5.3040E-01  9.7115E-01  1.6943E+00  1.0000E-02  9.4845E-01  7.8222E-02
             3.2894E+00
 PARAMETER:  1.3473E-01 -1.0392E+00 -3.0464E-01  4.5858E-01 -5.3413E-01  7.0727E-02  6.2725E-01 -5.3307E+00  4.7069E-02 -2.4482E+00
             1.2907E+00
 GRADIENT:  -2.5928E-01  1.2478E+00  3.3745E+00  3.2735E+00 -5.4045E+00 -7.0854E-02 -1.8683E-01  0.0000E+00 -4.8963E-01 -8.9152E-03
            -1.9367E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1358.81993505122        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      962             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0342E+00  2.7797E-01  6.6412E-01  1.4483E+00  5.2237E-01  9.7002E-01  1.9034E+00  1.0000E-02  9.3843E-01  5.9767E-02
             3.2982E+00
 PARAMETER:  1.3366E-01 -1.1802E+00 -3.0929E-01  4.7037E-01 -5.4937E-01  6.9560E-02  7.4365E-01 -5.9833E+00  3.6457E-02 -2.7173E+00
             1.2934E+00
 GRADIENT:   4.9773E+01  3.3942E+00  1.4783E+00  5.9994E+01  9.9649E+00  3.3890E+00  5.6245E-01  0.0000E+00  1.8490E+00  8.9812E-03
             1.0495E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1358.82352098653        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1139
 NPARAMETR:  1.0342E+00  2.7632E-01  6.6424E-01  1.4499E+00  5.2210E-01  9.7032E-01  1.9034E+00  1.0000E-02  9.3829E-01  1.1213E-01
             3.2943E+00
 PARAMETER:  1.3359E-01 -1.1862E+00 -3.0911E-01  4.7152E-01 -5.4990E-01  6.9867E-02  7.4363E-01 -5.9833E+00  3.6304E-02 -2.0881E+00
             1.2922E+00
 GRADIENT:   1.4909E-01 -2.3868E-01 -1.2786E+00 -6.0480E-01  2.3497E+00 -5.1069E-03 -2.7311E-02  0.0000E+00  9.0017E-02  1.9709E-03
             7.9876E-01

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1358.82547705912        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1231
 NPARAMETR:  1.0341E+00  2.7614E-01  6.6455E-01  1.4504E+00  5.2165E-01  9.7039E-01  1.9055E+00  1.0000E-02  9.3789E-01  1.2828E-01
             3.2894E+00
 PARAMETER:  1.3356E-01 -1.1868E+00 -3.0864E-01  4.7186E-01 -5.5077E-01  6.9945E-02  7.4474E-01 -5.9833E+00  3.5876E-02 -1.9535E+00
             1.2907E+00
 GRADIENT:   2.9083E-02  3.1604E-02  1.2492E-01 -1.0382E-02  4.8190E-01 -3.9251E-02 -1.5105E-02  0.0000E+00 -4.2496E-02  2.2912E-03
             1.2983E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1231
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.1578E-04 -1.8121E-03  1.1358E-04 -1.2929E-02 -8.8313E-04
 SE:             2.8623E-02  8.2308E-03  2.0141E-04  2.5399E-02  3.8246E-03
 N:                     100         100         100         100         100

 P VAL.:         9.8841E-01  8.2574E-01  5.7281E-01  6.1071E-01  8.1739E-01

 ETASHRINKSD(%)  4.1087E+00  7.2426E+01  9.9325E+01  1.4912E+01  8.7187E+01
 ETASHRINKVR(%)  8.0487E+00  9.2397E+01  9.9995E+01  2.7600E+01  9.8358E+01
 EBVSHRINKSD(%)  3.8567E+00  7.3644E+01  9.9271E+01  1.4345E+01  8.7209E+01
 EBVSHRINKVR(%)  7.5647E+00  9.3054E+01  9.9995E+01  2.6631E+01  9.8364E+01
 RELATIVEINF(%)  8.3372E+01  2.5228E-01  2.5900E-04  5.6997E+00  5.4607E-02
 EPSSHRINKSD(%)  2.3030E+01
 EPSSHRINKVR(%)  4.0756E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1358.8254770591177     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -623.67465049537952     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     5.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1358.825       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.76E-01  6.65E-01  1.45E+00  5.22E-01  9.70E-01  1.91E+00  1.00E-02  9.38E-01  1.28E-01  3.29E+00
 


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
 #CPUT: Total CPU Time in Seconds,       36.079
Stop Time:
Sun Oct 24 02:21:06 CDT 2021
