Thu Sep 30 22:57:42 CDT 2021
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
$DATA ../../../../data/SD1/SL2/dat1.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      998
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      898
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2189.49043776422        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3504E+02  1.3360E+02  1.4460E+02  9.9526E+01  6.7148E+01  5.3414E+01 -7.5162E+01 -1.6939E+02 -5.9846E+01 -1.3401E+01
            -3.0053E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3051.39770075123        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0356E+00  1.0771E+00  9.0909E-01  9.9229E-01  1.0236E+00  9.5123E-01  1.0662E+00  1.0227E+00  9.7307E-01  1.0105E+00
             1.9893E+00
 PARAMETER:  1.3500E-01  1.7423E-01  4.6898E-03  9.2255E-02  1.2330E-01  5.0006E-02  1.6412E-01  1.2241E-01  7.2702E-02  1.1043E-01
             7.8776E-01
 GRADIENT:   2.0397E+02  5.0138E+01 -1.2575E+01  4.8020E+01  2.7033E+00 -1.3513E+00  5.1410E+00  1.7700E+00 -3.2448E+00 -3.3011E+00
            -3.7445E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3054.86150195612        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      263
 NPARAMETR:  1.0423E+00  1.3722E+00  1.2357E+00  8.6736E-01  1.3318E+00  1.0035E+00  8.2419E-01  9.6758E-01  1.1195E+00  1.2802E+00
             2.0376E+00
 PARAMETER:  1.4140E-01  4.1640E-01  3.1160E-01 -4.2306E-02  3.8650E-01  1.0346E-01 -9.3359E-02  6.7043E-02  2.1286E-01  3.4705E-01
             8.1176E-01
 GRADIENT:   5.6877E+01  6.1615E+01  2.1709E+01  7.4843E+01 -2.6159E+01  7.6357E+00 -2.5379E+00 -4.1541E+00  1.6804E+00  4.6813E+00
            -2.5421E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3061.28274248418        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  1.0112E+00  1.4954E+00  1.2646E+00  7.6067E-01  1.4945E+00  9.8338E-01  7.5743E-01  1.5515E+00  1.1884E+00  1.3302E+00
             2.0446E+00
 PARAMETER:  1.1118E-01  5.0238E-01  3.3472E-01 -1.7355E-01  5.0178E-01  8.3241E-02 -1.7783E-01  5.3923E-01  2.7262E-01  3.8534E-01
             8.1522E-01
 GRADIENT:  -9.4813E+00  1.7811E+01  6.4648E+00  3.6971E+01  6.5244E+00  1.3766E+00 -1.9063E+00 -2.6656E+00  1.0486E+00  2.3635E+00
             4.2411E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3065.34864417065        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      616
 NPARAMETR:  1.0216E+00  1.8759E+00  9.9911E-01  5.0093E-01  1.7615E+00  9.8117E-01  6.7674E-01  2.1956E+00  1.5300E+00  1.4843E+00
             2.0340E+00
 PARAMETER:  1.2133E-01  7.2910E-01  9.9107E-02 -5.9129E-01  6.6614E-01  8.0992E-02 -2.9046E-01  8.8644E-01  5.2529E-01  4.9495E-01
             8.1001E-01
 GRADIENT:   1.3583E+01  2.1295E+01  4.7756E-01  1.1175E+01 -2.6426E+00  2.6597E-01 -1.0863E+00  3.7192E-01  7.9949E-01  2.7917E+00
             1.8110E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3065.79913154715        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      793
 NPARAMETR:  1.0108E+00  2.0922E+00  6.6056E-01  3.5303E-01  1.8799E+00  9.8437E-01  6.5584E-01  1.8845E+00  1.8816E+00  1.5521E+00
             2.0367E+00
 PARAMETER:  1.1079E-01  8.3820E-01 -3.1467E-01 -9.4120E-01  7.3120E-01  8.4245E-02 -3.2183E-01  7.3368E-01  7.3210E-01  5.3958E-01
             8.1133E-01
 GRADIENT:  -1.0347E+01  2.4597E+01  7.6661E-01  5.7455E+00 -6.9491E+00  1.3199E+00 -4.2670E-01  1.6516E+00  5.0234E-03  9.9492E-01
            -5.0192E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3065.85577824533        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      968
 NPARAMETR:  1.0130E+00  2.1864E+00  5.4971E-01  3.0298E-01  1.9468E+00  9.8213E-01  6.4925E-01  1.6711E+00  2.0761E+00  1.5948E+00
             2.0442E+00
 PARAMETER:  1.1292E-01  8.8224E-01 -4.9837E-01 -1.0941E+00  7.6620E-01  8.1971E-02 -3.3193E-01  6.1351E-01  8.3051E-01  5.6678E-01
             8.1499E-01
 GRADIENT:  -5.9835E+00  5.1578E+01 -3.0700E+00  1.4442E+01 -3.7562E+00  3.6084E-01  7.7257E-01  3.4791E+00 -1.6251E+00  3.0347E+00
            -2.8511E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3066.44461702709        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1147
 NPARAMETR:  1.0185E+00  2.2657E+00  3.7398E-01  2.3760E-01  2.0265E+00  9.7650E-01  6.4018E-01  1.1715E+00  2.3777E+00  1.6255E+00
             2.0484E+00
 PARAMETER:  1.1835E-01  9.1789E-01 -8.8355E-01 -1.3372E+00  8.0631E-01  7.6224E-02 -3.4601E-01  2.5827E-01  9.6612E-01  5.8579E-01
             8.1706E-01
 GRADIENT:   6.3018E+00  1.9051E+01 -6.9847E+00  1.1491E+01  7.8368E+00 -1.9137E+00 -7.3898E-01  1.4827E+00 -1.8847E+00  2.2240E+00
             3.8574E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3067.56678348974        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1323
 NPARAMETR:  1.0159E+00  2.3312E+00  3.4691E-01  1.7533E-01  2.0862E+00  9.8128E-01  6.2978E-01  8.6267E-01  2.9291E+00  1.6366E+00
             2.0445E+00
 PARAMETER:  1.1575E-01  9.4640E-01 -9.5869E-01 -1.6411E+00  8.3536E-01  8.1101E-02 -3.6238E-01 -4.7728E-02  1.1747E+00  5.9262E-01
             8.1514E-01
 GRADIENT:   1.1638E+00 -5.4830E+00 -5.3112E-01  4.4787E-02  9.1533E-01  7.8800E-02  6.2972E-01  6.1414E-01  9.1182E-01  3.6128E-01
             6.5596E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3067.85390002513        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1505
 NPARAMETR:  1.0155E+00  2.3262E+00  3.6276E-01  1.7786E-01  2.0830E+00  9.8106E-01  6.2704E-01  2.0914E-01  2.9123E+00  1.6318E+00
             2.0475E+00
 PARAMETER:  1.1535E-01  9.4424E-01 -9.1402E-01 -1.6268E+00  8.3379E-01  8.0879E-02 -3.6674E-01 -1.4648E+00  1.1689E+00  5.8970E-01
             8.1662E-01
 GRADIENT:   3.0126E-01 -3.4830E+00 -1.3226E-01 -6.1193E-01  5.9027E-01  4.4553E-02  3.0257E-02  2.9941E-02  1.9046E-01  4.8461E-02
             1.1365E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3067.93360954712        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1681
 NPARAMETR:  1.0153E+00  2.3137E+00  3.7714E-01  1.8780E-01  2.0694E+00  9.8090E-01  6.2749E-01  3.3479E-02  2.8277E+00  1.6249E+00
             2.0466E+00
 PARAMETER:  1.1516E-01  9.3885E-01 -8.7515E-01 -1.5724E+00  8.2725E-01  8.0717E-02 -3.6603E-01 -3.2968E+00  1.1395E+00  5.8543E-01
             8.1618E-01
 GRADIENT:  -1.6582E-01  5.1443E-01 -3.2899E-02  9.5331E-03  7.9553E-02 -1.9925E-02 -3.0481E-01  8.8734E-04 -1.5618E-01 -1.6447E-01
            -3.8428E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3067.95553738117        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1866
 NPARAMETR:  1.0155E+00  2.3072E+00  3.8234E-01  1.9089E-01  2.0629E+00  9.8095E-01  6.2918E-01  1.0000E-02  2.8214E+00  1.6249E+00
             2.0464E+00
 PARAMETER:  1.1542E-01  9.3603E-01 -8.6146E-01 -1.5560E+00  8.2412E-01  8.0771E-02 -3.6334E-01 -5.2391E+00  1.1372E+00  5.8547E-01
             8.1610E-01
 GRADIENT:  -3.8087E-01 -2.3094E+00 -1.0064E-02  1.8681E-01 -1.1873E+00 -8.9714E-02  1.4038E-01  0.0000E+00  7.6765E-01  1.9160E-01
            -3.6417E-01

0ITERATION NO.:   56    OBJECTIVE VALUE:  -3067.95553738117        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1888
 NPARAMETR:  1.0155E+00  2.3072E+00  3.8234E-01  1.9089E-01  2.0629E+00  9.8095E-01  6.2918E-01  1.0000E-02  2.8214E+00  1.6249E+00
             2.0464E+00
 PARAMETER:  1.1542E-01  9.3603E-01 -8.6146E-01 -1.5560E+00  8.2412E-01  8.0771E-02 -3.6334E-01 -5.2391E+00  1.1372E+00  5.8547E-01
             8.1610E-01
 GRADIENT:  -3.8087E-01 -2.3094E+00 -1.0064E-02  1.8681E-01 -1.1873E+00 -8.9714E-02  1.4038E-01  0.0000E+00  7.6765E-01  1.9160E-01
            -3.6417E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1888
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY
 THIS MUST BE ADDRESSED BEFORE THE COVARIANCE STEP CAN BE IMPLEMENTED

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1174E-03 -3.8758E-02 -3.4892E-05  4.5785E-02 -2.4508E-02
 SE:             2.9642E-02  2.4869E-02  5.0455E-05  1.9463E-02  2.6712E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6993E-01  1.1912E-01  4.8923E-01  1.8654E-02  3.5889E-01

 ETASHRINKSD(%)  6.9527E-01  1.6686E+01  9.9831E+01  3.4795E+01  1.0511E+01
 ETASHRINKVR(%)  1.3857E+00  3.0587E+01  1.0000E+02  5.7484E+01  1.9918E+01
 EBVSHRINKSD(%)  1.0199E+00  1.5152E+01  9.9847E+01  4.1056E+01  7.7084E+00
 EBVSHRINKVR(%)  2.0295E+00  2.8008E+01  1.0000E+02  6.5256E+01  1.4823E+01
 RELATIVEINF(%)  9.7943E+01  2.0384E+01  1.4769E-04  9.4803E+00  5.5620E+01
 EPSSHRINKSD(%)  1.7428E+01
 EPSSHRINKVR(%)  3.1819E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          898
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1650.4136056355921     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3067.9555373811731     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1417.5419317455810     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.74
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3067.956       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.31E+00  3.82E-01  1.91E-01  2.06E+00  9.81E-01  6.29E-01  1.00E-02  2.82E+00  1.62E+00  2.05E+00
 


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
 #CPUT: Total CPU Time in Seconds,       31.769
Stop Time:
Thu Sep 30 22:58:14 CDT 2021
