Sun Oct 24 00:33:25 CDT 2021
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
$DATA ../../../../data/SD3/TD1/dat77.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2162.19539609672        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9218E+02  3.3812E+01  8.9416E+00  7.8301E+01 -3.7944E+01  4.0013E+01  5.5524E+00  7.7722E+00  4.8749E+01  1.5726E+00
             2.8696E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2168.39876439871        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  1.0122E+00  1.1251E+00  1.0564E+00  9.4083E-01  1.1273E+00  1.0925E+00  1.0091E+00  9.3100E-01  7.0376E-01  1.0469E+00
             9.4504E-01
 PARAMETER:  1.1209E-01  2.1789E-01  1.5485E-01  3.9009E-02  2.1985E-01  1.8848E-01  1.0903E-01  2.8505E-02 -2.5131E-01  1.4584E-01
             4.3467E-02
 GRADIENT:  -2.9869E+01  2.2246E+01  1.6756E+01  9.7273E+00 -2.7011E+00  1.9096E+01 -1.3572E+01 -2.6160E+00 -1.8340E+01 -1.4599E+01
            -2.3718E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2169.81073813808        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0108E+00  1.0150E+00  1.0750E+00  1.0085E+00  1.0965E+00  1.0775E+00  1.1526E+00  7.5090E-01  6.6122E-01  1.1477E+00
             9.4407E-01
 PARAMETER:  1.1077E-01  1.1490E-01  1.7233E-01  1.0842E-01  1.9208E-01  1.7463E-01  2.4198E-01 -1.8649E-01 -3.1366E-01  2.3778E-01
             4.2440E-02
 GRADIENT:  -3.1034E+01  1.9179E+01  1.0117E+01  1.6408E+01 -4.4298E+00  1.4206E+01 -6.2964E+00 -2.3129E+00 -1.6526E+01 -1.0680E-02
            -2.4463E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2172.28935779063        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  1.0289E+00  1.0470E+00  8.9410E-01  9.6916E-01  1.0149E+00  1.0384E+00  1.1127E+00  4.9313E-01  7.6590E-01  1.0402E+00
             9.7128E-01
 PARAMETER:  1.2849E-01  1.4592E-01 -1.1938E-02  6.8673E-02  1.1479E-01  1.3765E-01  2.0678E-01 -6.0697E-01 -1.6670E-01  1.3943E-01
             7.0856E-02
 GRADIENT:  -7.9798E-02 -4.5531E-01 -1.0080E+00  1.1349E+00  3.6194E-01 -7.9863E-02 -9.1882E-02  2.2095E-01  3.6108E-01  3.0924E-01
             3.7293E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2172.33498466614        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  1.0300E+00  1.1254E+00  8.2568E-01  9.1660E-01  1.0182E+00  1.0391E+00  1.0547E+00  3.4738E-01  7.9458E-01  1.0354E+00
             9.7099E-01
 PARAMETER:  1.2957E-01  2.1811E-01 -9.1547E-02  1.2911E-02  1.1803E-01  1.3838E-01  1.5327E-01 -9.5733E-01 -1.2994E-01  1.3475E-01
             7.0564E-02
 GRADIENT:   2.0242E-01 -1.8965E+00 -8.0388E-01 -1.7449E+00  7.6707E-01 -1.3296E-01  2.6757E-01  1.1731E-01  2.4777E-01  2.9275E-01
             3.9625E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2172.36030119802        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      897
 NPARAMETR:  1.0303E+00  1.1622E+00  7.8232E-01  8.9263E-01  1.0117E+00  1.0404E+00  1.0277E+00  1.2331E-01  8.0951E-01  1.0259E+00
             9.7082E-01
 PARAMETER:  1.2989E-01  2.5035E-01 -1.4549E-01 -1.3587E-02  1.1158E-01  1.3960E-01  1.2730E-01 -1.9930E+00 -1.1133E-01  1.2560E-01
             7.0389E-02
 GRADIENT:  -2.0586E-01  1.9672E-02  1.8732E-01  2.1590E-01  1.1055E-01  1.5400E-01 -1.6424E-01  3.6048E-03  3.5766E-02 -1.7274E-01
            -1.7216E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2172.36283404430        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1072
 NPARAMETR:  1.0306E+00  1.1731E+00  7.7318E-01  8.8536E-01  1.0123E+00  1.0401E+00  1.0216E+00  3.7751E-02  8.1327E-01  1.0259E+00
             9.7092E-01
 PARAMETER:  1.3012E-01  2.5967E-01 -1.5725E-01 -2.1759E-02  1.1219E-01  1.3928E-01  1.2134E-01 -3.1767E+00 -1.0669E-01  1.2553E-01
             7.0490E-02
 GRADIENT:   2.6220E-02 -1.2604E-01 -3.9580E-02 -1.4976E-01  5.0711E-02 -1.6167E-02  1.8413E-02  8.3263E-04  5.0314E-03  1.8661E-02
             3.2308E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2172.36553071069        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1233             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0321E+00  1.1730E+00  7.7310E-01  8.8553E-01  1.0122E+00  1.0409E+00  1.0217E+00  1.0000E-02  8.1335E-01  1.0259E+00
             9.7092E-01
 PARAMETER:  1.3161E-01  2.5952E-01 -1.5734E-01 -2.1564E-02  1.1208E-01  1.4008E-01  1.2145E-01 -4.8254E+00 -1.0660E-01  1.2558E-01
             7.0489E-02
 GRADIENT:   6.3605E+02  1.4695E+02  2.8457E+00  7.2001E+01  8.9974E+00  7.8481E+01  6.0110E+00  0.0000E+00  4.4122E+00  1.1438E+00
             1.0009E+00

0ITERATION NO.:   37    OBJECTIVE VALUE:  -2172.36553071069        NO. OF FUNC. EVALS.:  55
 CUMULATIVE NO. OF FUNC. EVALS.:     1288
 NPARAMETR:  1.0321E+00  1.1730E+00  7.7310E-01  8.8553E-01  1.0122E+00  1.0409E+00  1.0217E+00  1.0000E-02  8.1335E-01  1.0259E+00
             9.7092E-01
 PARAMETER:  1.3161E-01  2.5952E-01 -1.5734E-01 -2.1564E-02  1.1208E-01  1.4008E-01  1.2145E-01 -4.8254E+00 -1.0660E-01  1.2558E-01
             7.0489E-02
 GRADIENT:   3.1455E+00 -6.7811E-02 -2.7253E-02  1.8730E-02  1.1437E-01  3.0671E-01  1.8391E-02  0.0000E+00  2.6878E-02  1.9950E-02
             9.2776E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1288
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3911E-04 -7.9504E-03 -4.2293E-04  4.1771E-03 -1.9322E-02
 SE:             2.9897E-02  2.2626E-02  1.7221E-04  2.2914E-02  2.3739E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9362E-01  7.2531E-01  1.4051E-02  8.5535E-01  4.1569E-01

 ETASHRINKSD(%)  1.0000E-10  2.4198E+01  9.9423E+01  2.3234E+01  2.0471E+01
 ETASHRINKVR(%)  1.0000E-10  4.2541E+01  9.9997E+01  4.1070E+01  3.6751E+01
 EBVSHRINKSD(%)  2.7624E-01  2.3999E+01  9.9507E+01  2.4432E+01  1.8215E+01
 EBVSHRINKVR(%)  5.5171E-01  4.2239E+01  9.9998E+01  4.2895E+01  3.3112E+01
 RELATIVEINF(%)  9.9039E+01  3.0731E+00  2.6935E-04  2.9137E+00  1.0518E+01
 EPSSHRINKSD(%)  3.3092E+01
 EPSSHRINKVR(%)  5.5233E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2172.3655307106883     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1253.4269975060156     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2172.366       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.17E+00  7.73E-01  8.86E-01  1.01E+00  1.04E+00  1.02E+00  1.00E-02  8.13E-01  1.03E+00  9.71E-01
 


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
 #CPUT: Total CPU Time in Seconds,       40.921
Stop Time:
Sun Oct 24 00:33:34 CDT 2021
