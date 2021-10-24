Sun Oct 24 00:28:55 CDT 2021
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
$DATA ../../../../data/SD3/TD1/dat46.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m46.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2133.83399020605        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6427E+02 -7.3282E+01 -5.0122E+01 -2.1748E+01  4.7316E+01  4.8396E+01  3.2485E+00  1.0238E+01  3.0990E+01  3.1180E+01
             7.9277E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2148.09155613971        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.6903E-01  1.1258E+00  1.1399E+00  1.0092E+00  1.0592E+00  9.5368E-01  9.7930E-01  9.6315E-01  8.8214E-01  8.6198E-01
             9.9089E-01
 PARAMETER:  6.8543E-02  2.1847E-01  2.3091E-01  1.0917E-01  1.5754E-01  5.2569E-02  7.9084E-02  6.2458E-02 -2.5407E-02 -4.8524E-02
             9.0847E-02
 GRADIENT:   3.2531E+00  9.1119E+00  3.5659E+00  9.9412E+00 -2.3846E-01 -5.5072E+00 -3.0117E+00 -5.0587E+00 -8.3025E-01 -9.5662E-01
            -8.3442E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2148.49285663344        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.6734E-01  1.0127E+00  1.2545E+00  1.0814E+00  1.0493E+00  9.6360E-01  1.0405E+00  1.1454E+00  8.1820E-01  8.4102E-01
             9.9861E-01
 PARAMETER:  6.6792E-02  1.1261E-01  3.2671E-01  1.7822E-01  1.4811E-01  6.2916E-02  1.3972E-01  2.3577E-01 -1.0065E-01 -7.3141E-02
             9.8608E-02
 GRADIENT:   1.1270E+00  3.9526E+00 -8.0456E-01  4.5430E+00 -4.4911E+00 -9.1555E-01 -4.8153E+00  1.4119E+00 -5.7077E+00 -2.0369E+00
            -2.2879E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2149.09629464605        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      527
 NPARAMETR:  9.6377E-01  8.4696E-01  1.6966E+00  1.2044E+00  1.1281E+00  9.6459E-01  1.1579E+00  1.3679E+00  8.0195E-01  9.5164E-01
             1.0010E+00
 PARAMETER:  6.3102E-02 -6.6100E-02  6.2861E-01  2.8595E-01  2.2052E-01  6.3946E-02  2.4660E-01  4.1327E-01 -1.2071E-01  5.0430E-02
             1.0104E-01
 GRADIENT:  -2.8101E+00  8.9753E+00  4.0767E+00  9.0461E+00 -5.5938E+00  3.1437E-01  7.0409E-01 -1.2070E+00  9.3610E-01 -6.8052E-01
            -1.4562E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2149.38384083182        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.6574E-01  6.5741E-01  1.8647E+00  1.3339E+00  1.1076E+00  9.6141E-01  1.3668E+00  1.3983E+00  7.5097E-01  9.5071E-01
             1.0023E+00
 PARAMETER:  6.5139E-02 -3.1945E-01  7.2309E-01  3.8814E-01  2.0221E-01  6.0646E-02  4.1247E-01  4.3529E-01 -1.8640E-01  4.9459E-02
             1.0230E-01
 GRADIENT:   6.2346E+00  1.0600E+01  4.5423E+00  2.0924E+01 -5.6772E+00 -2.4606E-01  8.2292E-01 -3.0587E+00  8.5855E-01 -8.7361E-01
             9.2056E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2149.83455229322        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  9.6076E-01  4.3623E-01  2.1497E+00  1.4779E+00  1.1117E+00  9.5913E-01  1.6227E+00  1.6046E+00  7.1075E-01  9.7390E-01
             9.9923E-01
 PARAMETER:  5.9968E-02 -7.2958E-01  8.6532E-01  4.9060E-01  2.0592E-01  5.8270E-02  5.8411E-01  5.7285E-01 -2.4144E-01  7.3553E-02
             9.9232E-02
 GRADIENT:  -1.0459E+00  4.2032E+00 -1.8062E+00  1.3205E+01 -3.2753E-01 -3.4014E-01 -1.1257E-01  1.0088E+00 -5.3421E-01  4.4158E-01
            -8.3675E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2149.90117594118        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1059             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6086E-01  3.4856E-01  2.2907E+00  1.5238E+00  1.1170E+00  9.5933E-01  1.8066E+00  1.6735E+00  7.0069E-01  9.7856E-01
             1.0000E+00
 PARAMETER:  6.0077E-02 -9.5394E-01  9.2888E-01  5.2119E-01  2.1067E-01  5.8483E-02  6.9146E-01  6.1493E-01 -2.5569E-01  7.8328E-02
             1.0000E-01
 GRADIENT:   3.9426E+02  4.1861E+01  8.1844E+00  7.8201E+02  1.3334E+01  4.1284E+01  6.9699E+00  3.1328E+00  2.2563E+01  4.3511E-01
             1.1016E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2149.92546807414        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1244             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6097E-01  3.5181E-01  2.2911E+00  1.5274E+00  1.1155E+00  9.5937E-01  1.8364E+00  1.6738E+00  6.9604E-01  9.7741E-01
             1.0000E+00
 PARAMETER:  6.0192E-02 -9.4465E-01  9.2901E-01  5.2355E-01  2.0933E-01  5.8518E-02  7.0780E-01  6.1512E-01 -2.6235E-01  7.7150E-02
             1.0004E-01
 GRADIENT:   3.9428E+02  4.4302E+01  8.6824E+00  7.9930E+02  1.0545E+01  4.1262E+01  7.4180E+00  3.1157E+00  2.1934E+01  6.0309E-01
             8.9015E-01

0ITERATION NO.:   37    OBJECTIVE VALUE:  -2149.92546807414        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1301
 NPARAMETR:  9.6097E-01  3.5181E-01  2.2911E+00  1.5274E+00  1.1155E+00  9.5937E-01  1.8364E+00  1.6738E+00  6.9604E-01  9.7741E-01
             1.0000E+00
 PARAMETER:  6.0192E-02 -9.4465E-01  9.2901E-01  5.2355E-01  2.0933E-01  5.8518E-02  7.0780E-01  6.1512E-01 -2.6235E-01  7.7150E-02
             1.0004E-01
 GRADIENT:  -6.9385E-02  4.2050E-01  1.3788E-01  2.9389E+00 -4.0956E-02 -1.3767E-02 -2.6330E-02  3.0115E-02  1.1250E-01  3.1769E-02
            -1.0982E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1301
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.8626E-04 -3.3184E-04 -3.8650E-02 -7.0762E-03 -4.4395E-02
 SE:             2.9859E-02  1.1351E-02  1.7924E-02  2.7310E-02  2.0217E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9235E-01  9.7668E-01  3.1058E-02  7.9555E-01  2.8093E-02

 ETASHRINKSD(%)  1.0000E-10  6.1972E+01  3.9952E+01  8.5076E+00  3.2272E+01
 ETASHRINKVR(%)  1.0000E-10  8.5539E+01  6.3942E+01  1.6291E+01  5.4129E+01
 EBVSHRINKSD(%)  3.7486E-01  6.2983E+01  4.4922E+01  8.5137E+00  2.8366E+01
 EBVSHRINKVR(%)  7.4831E-01  8.6297E+01  6.9664E+01  1.6303E+01  4.8686E+01
 RELATIVEINF(%)  9.6610E+01  2.5522E-01  7.4833E+00  1.6677E+00  1.0448E+01
 EPSSHRINKSD(%)  3.3620E+01
 EPSSHRINKVR(%)  5.5937E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2149.9254680741351     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1230.9869348694624     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     6.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2149.925       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.61E-01  3.52E-01  2.29E+00  1.53E+00  1.12E+00  9.59E-01  1.84E+00  1.67E+00  6.96E-01  9.77E-01  1.00E+00
 


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
 #CPUT: Total CPU Time in Seconds,       45.626
Stop Time:
Sun Oct 24 00:29:04 CDT 2021
