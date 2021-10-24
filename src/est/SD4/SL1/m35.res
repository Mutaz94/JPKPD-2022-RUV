Sun Oct 24 02:58:51 CDT 2021
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
$DATA ../../../../data/SD4/SL1/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1687.64029348646        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4536E+02 -1.5262E+00 -5.2667E+01  7.8756E+01  1.0026E+02  6.8749E+01  1.7579E+00  1.1148E+01  2.1456E+00 -1.4693E+01
            -1.5738E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1696.28307702005        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0500E+00  1.0243E+00  1.1286E+00  1.0296E+00  9.9062E-01  9.3233E-01  9.7906E-01  8.8293E-01  1.0378E+00  1.0697E+00
             1.0561E+00
 PARAMETER:  1.4881E-01  1.2398E-01  2.2100E-01  1.2921E-01  9.0574E-02  2.9933E-02  7.8833E-02 -2.4506E-02  1.3712E-01  1.6735E-01
             1.5460E-01
 GRADIENT:  -1.2058E+01  3.0923E+01  3.1711E-01  4.0898E+01 -7.3722E+00 -9.5248E+00  2.3598E+00  4.0586E+00  1.4010E+00 -5.8810E+00
             5.7468E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1697.72859500523        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0514E+00  9.4810E-01  1.1947E+00  1.0560E+00  9.9642E-01  9.3986E-01  8.3640E-01  6.4781E-01  1.0463E+00  1.1828E+00
             1.0293E+00
 PARAMETER:  1.5010E-01  4.6700E-02  2.7788E-01  1.5448E-01  9.6416E-02  3.7976E-02 -7.8654E-02 -3.3416E-01  1.4527E-01  2.6786E-01
             1.2887E-01
 GRADIENT:  -5.3366E+00  1.1961E+01  9.1990E+00  1.1757E+01 -1.1557E+01 -5.6990E+00 -7.5778E-01 -1.9579E+00  1.8475E+00  4.5366E-01
            -4.2466E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1698.11139446880        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      547
 NPARAMETR:  1.0529E+00  8.7769E-01  1.1380E+00  1.0949E+00  9.4998E-01  9.5240E-01  1.0027E+00  6.2142E-01  9.8121E-01  1.1112E+00
             1.0404E+00
 PARAMETER:  1.5155E-01 -3.0459E-02  2.2931E-01  1.9069E-01  4.8688E-02  5.1227E-02  1.0274E-01 -3.7575E-01  8.1031E-02  2.0546E-01
             1.3958E-01
 GRADIENT:  -1.7163E+00  5.2070E+00  3.1906E+00  5.5093E+00 -3.6642E+00 -5.0816E-01  9.1486E-02 -5.4019E-01  4.9464E-01 -1.2855E+00
             8.9066E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1698.20443296406        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      723
 NPARAMETR:  1.0519E+00  7.1881E-01  1.2764E+00  1.1995E+00  9.4813E-01  9.5315E-01  1.0221E+00  7.7405E-01  9.2429E-01  1.1409E+00
             1.0372E+00
 PARAMETER:  1.5064E-01 -2.3016E-01  3.4407E-01  2.8187E-01  4.6738E-02  5.2019E-02  1.2190E-01 -1.5611E-01  2.1274E-02  2.3181E-01
             1.3649E-01
 GRADIENT:   9.4312E-01  4.2251E+00  8.5660E-01  6.5256E+00 -1.5501E+00  6.8977E-01 -3.1417E-01  3.1656E-01 -9.2440E-01 -8.0105E-01
            -1.2481E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1698.23112689224        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  1.0512E+00  6.1163E-01  1.3177E+00  1.2687E+00  9.2816E-01  9.5245E-01  1.0544E+00  7.8561E-01  8.8558E-01  1.1511E+00
             1.0339E+00
 PARAMETER:  1.4992E-01 -3.9163E-01  3.7592E-01  3.3797E-01  2.5444E-02  5.1278E-02  1.5294E-01 -1.4130E-01 -2.1510E-02  2.4069E-01
             1.3337E-01
 GRADIENT:   2.1972E+00  4.8943E+00  7.5395E-01  9.6005E+00 -2.7142E+00  8.9057E-01 -6.3051E-01 -3.6785E-02 -2.3085E+00  4.3065E-01
            -1.4767E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1698.30356595381        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1080
 NPARAMETR:  1.0498E+00  5.7889E-01  1.3247E+00  1.2839E+00  9.2300E-01  9.4966E-01  1.1460E+00  7.9338E-01  8.7663E-01  1.1440E+00
             1.0373E+00
 PARAMETER:  1.4860E-01 -4.4665E-01  3.8120E-01  3.4992E-01  1.9874E-02  4.8354E-02  2.3628E-01 -1.3145E-01 -3.1674E-02  2.3457E-01
             1.3661E-01
 GRADIENT:  -7.4976E-02  8.2506E-01 -5.6617E-01 -7.8729E-01  1.1940E+00 -1.1893E-01 -3.1044E-02  1.0633E-01  2.2813E-01  1.4767E-01
             2.4018E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1698.30602155084        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1256
 NPARAMETR:  1.0484E+00  5.5673E-01  1.3260E+00  1.2987E+00  9.1468E-01  9.4964E-01  1.1993E+00  7.9318E-01  8.6303E-01  1.1350E+00
             1.0369E+00
 PARAMETER:  1.4730E-01 -4.8567E-01  3.8219E-01  3.6133E-01  1.0818E-02  4.8330E-02  2.8177E-01 -1.3171E-01 -4.7307E-02  2.2664E-01
             1.3622E-01
 GRADIENT:  -2.7294E+00  1.9284E+00  6.2242E-01  9.0590E-01 -5.0629E-01 -6.1863E-02 -1.1400E-01 -3.9710E-02 -1.0041E+00 -3.8500E-01
            -1.9036E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1698.30645336660        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1436
 NPARAMETR:  1.0479E+00  5.4187E-01  1.3266E+00  1.3082E+00  9.0973E-01  9.4955E-01  1.2293E+00  7.9356E-01  8.5593E-01  1.1309E+00
             1.0368E+00
 PARAMETER:  1.4683E-01 -5.1274E-01  3.8265E-01  3.6868E-01  5.3897E-03  4.8230E-02  3.0642E-01 -1.3123E-01 -5.5565E-02  2.2300E-01
             1.3615E-01
 GRADIENT:  -3.4728E+00  2.2771E+00  8.6856E-01  1.6379E+00 -1.0054E+00 -4.9091E-02 -1.4847E-01 -7.6014E-02 -1.4292E+00 -5.4469E-01
            -3.0768E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1698.33208256525        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1619             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0513E+00  5.3434E-01  1.3260E+00  1.3092E+00  9.0863E-01  9.4987E-01  1.2296E+00  7.9223E-01  8.5785E-01  1.1338E+00
             1.0369E+00
 PARAMETER:  1.4999E-01 -5.2672E-01  3.8217E-01  3.6940E-01  4.1796E-03  4.8566E-02  3.0671E-01 -1.3290E-01 -5.3330E-02  2.2556E-01
             1.3628E-01
 GRADIENT:   6.4776E+02  5.8999E+01  6.3249E+00  4.5396E+02  7.7241E+00  4.5307E+01  3.7259E+00  1.4877E-01  6.9631E+00  2.1172E+00
             1.0473E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1698.33246676215        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1786
 NPARAMETR:  1.0505E+00  5.3426E-01  1.3256E+00  1.3101E+00  9.0832E-01  9.4978E-01  1.2335E+00  7.9175E-01  8.5850E-01  1.1342E+00
             1.0370E+00
 PARAMETER:  1.4928E-01 -5.2688E-01  3.8188E-01  3.7009E-01  3.8451E-03  4.8473E-02  3.0986E-01 -1.3350E-01 -5.2572E-02  2.2594E-01
             1.3635E-01
 GRADIENT:  -7.4708E-01  4.0495E-01  8.8277E-02  8.1145E-01  1.4500E-01 -2.8563E-02  1.1942E-02  1.8495E-02  9.6374E-02  3.7569E-02
             2.6046E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1786
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3482E-04 -1.1486E-02 -2.1635E-02 -2.0835E-03 -2.9174E-02
 SE:             2.9846E-02  1.0798E-02  1.1547E-02  2.7499E-02  2.3775E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9640E-01  2.8745E-01  6.0982E-02  9.3960E-01  2.1979E-01

 ETASHRINKSD(%)  1.0978E-02  6.3827E+01  6.1316E+01  7.8762E+00  2.0352E+01
 ETASHRINKVR(%)  2.1954E-02  8.6915E+01  8.5036E+01  1.5132E+01  3.6562E+01
 EBVSHRINKSD(%)  4.8608E-01  6.4018E+01  6.4713E+01  8.0132E+00  1.6917E+01
 EBVSHRINKVR(%)  9.6980E-01  8.7053E+01  8.7548E+01  1.5384E+01  3.0972E+01
 RELATIVEINF(%)  9.6659E+01  3.3647E-01  2.7004E+00  2.6420E+00  7.5337E+00
 EPSSHRINKSD(%)  4.3781E+01
 EPSSHRINKVR(%)  6.8394E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1698.3324667621528     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -963.18164019841458     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.19
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1698.332       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  5.34E-01  1.33E+00  1.31E+00  9.08E-01  9.50E-01  1.23E+00  7.92E-01  8.58E-01  1.13E+00  1.04E+00
 


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
 #CPUT: Total CPU Time in Seconds,       52.846
Stop Time:
Sun Oct 24 02:59:02 CDT 2021
