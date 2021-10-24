Sun Oct 24 02:21:06 CDT 2021
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
$DATA ../../../../data/SD4/A2/dat88.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m88.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -805.122141656283        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3700E+02  6.1920E+00 -2.2962E+01  6.0756E+01  2.1298E+02  5.2917E+01 -4.9250E+01 -1.7380E+01 -6.2996E+01 -6.8661E+01
            -1.5055E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1376.36962532643        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0020E+00  1.1017E+00  1.1135E+00  1.0526E+00  9.5794E-01  8.1260E-01  9.6141E-01  9.6315E-01  9.2083E-01  8.0057E-01
             2.8604E+00
 PARAMETER:  1.0201E-01  1.9686E-01  2.0752E-01  1.5126E-01  5.7026E-02 -1.0752E-01  6.0645E-02  6.2453E-02  1.7515E-02 -1.2243E-01
             1.1509E+00
 GRADIENT:  -2.0854E+00  6.6513E+01 -1.2170E+00  8.8620E+01 -2.0889E+01 -5.2178E+01 -1.1115E+00  2.0329E+00 -4.2576E+00  1.1959E+01
            -2.4548E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1386.48373291303        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  9.8991E-01  1.0322E+00  2.4416E+00  1.1056E+00  1.2874E+00  9.4456E-01  8.0439E-01  7.2695E-01  9.2230E-01  4.5512E-01
             3.0491E+00
 PARAMETER:  8.9858E-02  1.3171E-01  9.9265E-01  2.0040E-01  3.5262E-01  4.2967E-02 -1.1767E-01 -2.1890E-01  1.9111E-02 -6.8720E-01
             1.2148E+00
 GRADIENT:  -2.0901E+01  2.3035E+01 -4.9330E+00  5.3171E+01  2.8269E+01  6.8774E+00  1.6812E+00  3.0742E-01  3.4001E+00  1.5587E+00
            -1.0467E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1389.28476500270        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.9598E-01  1.0301E+00  1.4849E+00  1.0356E+00  1.0629E+00  9.1859E-01  7.0957E-01  3.0303E-01  9.8704E-01  3.3661E-01
             3.0776E+00
 PARAMETER:  9.5976E-02  1.2964E-01  4.9536E-01  1.3497E-01  1.6097E-01  1.5082E-02 -2.4310E-01 -1.0939E+00  8.6953E-02 -9.8883E-01
             1.2242E+00
 GRADIENT:   1.7624E+00 -2.6279E+00 -1.0261E+00 -2.2995E+00  7.2860E-01  7.7779E-01  4.1603E-01  2.3726E-01  7.2607E-01  8.8703E-01
             2.6440E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1389.79441640867        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:      321
 NPARAMETR:  9.9497E-01  1.0406E+00  1.6731E+00  1.0337E+00  1.1087E+00  9.1728E-01  6.1012E-01  4.5974E-02  1.0057E+00  9.2523E-02
             3.0975E+00
 PARAMETER:  9.4957E-02  1.3984E-01  6.1471E-01  1.3316E-01  2.0319E-01  1.3656E-02 -3.9411E-01 -2.9797E+00  1.0570E-01 -2.2803E+00
             1.2306E+00
 GRADIENT:  -4.5616E+01 -4.4914E+00 -3.3875E-01 -7.6555E+00 -6.5896E-01 -3.7890E+00 -3.9081E-01  4.3461E-03 -4.7474E-01  2.8106E-02
            -9.5295E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1390.43522351558        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      498
 NPARAMETR:  1.0150E+00  1.0396E+00  1.8364E+00  1.0500E+00  1.1409E+00  9.2555E-01  6.6369E-01  2.4048E-02  9.6923E-01  5.9587E-02
             3.1434E+00
 PARAMETER:  1.1492E-01  1.3887E-01  7.0781E-01  1.4881E-01  2.3184E-01  2.2633E-02 -3.0995E-01 -3.6277E+00  6.8750E-02 -2.7203E+00
             1.2453E+00
 GRADIENT:  -3.7609E-01  2.5434E-01  2.2030E-01  5.3247E-02 -4.4815E-01 -2.1288E-01 -7.8566E-02  9.0288E-04 -9.2186E-02  1.6802E-02
            -8.1539E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1390.44690791186        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      673
 NPARAMETR:  1.0151E+00  1.0248E+00  1.8068E+00  1.0596E+00  1.1325E+00  9.2584E-01  6.8860E-01  1.0000E-02  9.5859E-01  1.0000E-02
             3.1435E+00
 PARAMETER:  1.1497E-01  1.2449E-01  6.9156E-01  1.5792E-01  2.2440E-01  2.2945E-02 -2.7309E-01 -6.2088E+00  5.7704E-02 -4.5188E+00
             1.2453E+00
 GRADIENT:  -6.1009E-02  1.1883E-01  4.2798E-02  1.1155E-01 -8.7655E-02 -7.0193E-02 -6.5853E-03  0.0000E+00  1.4590E-02  2.2023E-04
             1.9216E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1390.44713055768        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:      844
 NPARAMETR:  1.0153E+00  1.0084E+00  1.8051E+00  1.0700E+00  1.1281E+00  9.2611E-01  6.9609E-01  1.0000E-02  9.4999E-01  1.0000E-02
             3.1426E+00
 PARAMETER:  1.1489E-01  1.0860E-01  6.9172E-01  1.6778E-01  2.2027E-01  2.3014E-02 -2.5968E-01 -8.3360E+00  4.9238E-02 -6.0059E+00
             1.2452E+00
 GRADIENT:  -1.9992E-01  6.0671E-02  2.0389E-02  7.5821E-02 -7.3079E-02 -2.5696E-02  8.2881E-03  0.0000E+00  2.1779E-02  0.0000E+00
             4.9460E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      844
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.7060E-04 -1.3106E-02  5.7972E-05 -3.4385E-03 -3.4204E-05
 SE:             2.8800E-02  1.1441E-02  4.3191E-05  2.2444E-02  1.7041E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7865E-01  2.5198E-01  1.7952E-01  8.7824E-01  8.4092E-01

 ETASHRINKSD(%)  3.5148E+00  6.1671E+01  9.9855E+01  2.4811E+01  9.9429E+01
 ETASHRINKVR(%)  6.9060E+00  8.5309E+01  1.0000E+02  4.3465E+01  9.9997E+01
 EBVSHRINKSD(%)  3.5372E+00  6.1784E+01  9.9830E+01  2.4426E+01  9.9409E+01
 EBVSHRINKVR(%)  6.9493E+00  8.5396E+01  1.0000E+02  4.2886E+01  9.9997E+01
 RELATIVEINF(%)  9.0435E+01  2.0293E-01  3.5560E-05  9.9427E-01  2.6724E-04
 EPSSHRINKSD(%)  2.1648E+01
 EPSSHRINKVR(%)  3.8610E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1390.4471305576833     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -655.29630399394512     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     3.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1390.447       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.01E+00  1.81E+00  1.07E+00  1.13E+00  9.26E-01  6.98E-01  1.00E-02  9.51E-01  1.00E-02  3.14E+00
 


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
 #CPUT: Total CPU Time in Seconds,       23.235
Stop Time:
Sun Oct 24 02:21:13 CDT 2021
