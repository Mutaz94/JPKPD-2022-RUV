Sat Sep 18 11:01:48 CDT 2021
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
$DATA ../../../../data/spa/S1/dat33.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m33.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1615.10648349155        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6393E+02 -1.3231E+01  5.5751E+00  1.2686E+00 -4.3962E+00  1.5268E+01 -8.3979E+00 -2.8266E+00  2.7530E+01 -8.6216E+00
            -3.3356E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1621.96287617846        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.6766E-01  1.1889E+00  9.6004E-01  8.9083E-01  1.1081E+00  9.3536E-01  1.2440E+00  1.0241E+00  6.6489E-01  1.0858E+00
             1.0964E+00
 PARAMETER:  6.7123E-02  2.7305E-01  5.9221E-02 -1.5603E-02  2.0268E-01  3.3180E-02  3.1829E-01  1.2378E-01 -3.0813E-01  1.8230E-01
             1.9208E-01
 GRADIENT:   8.8463E+01  2.4114E+01  1.5492E+00  5.7123E+00  1.7481E+01 -4.4292E+00  1.1750E+01 -3.3635E+00 -1.1919E+01 -2.8611E+00
             1.9849E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1622.77935172565        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.6048E-01  1.0692E+00  1.0202E+00  9.7514E-01  1.0572E+00  9.3428E-01  1.3833E+00  1.1111E+00  6.5055E-01  1.0529E+00
             1.0936E+00
 PARAMETER:  5.9676E-02  1.6694E-01  1.2002E-01  7.4823E-02  1.5558E-01  3.2024E-02  4.2446E-01  2.0532E-01 -3.2994E-01  1.5154E-01
             1.8947E-01
 GRADIENT:   7.1751E+01  2.8454E+01 -1.4493E+00  2.0712E+01  5.7188E+00 -3.8863E+00  1.5501E+01 -2.3601E-02 -9.3955E+00  6.2881E-02
             4.3913E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1625.15282398111        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.3254E-01  1.0526E+00  8.3544E-01  9.5930E-01  9.5449E-01  9.3851E-01  1.2170E+00  7.2825E-01  7.7696E-01  9.7983E-01
             1.0787E+00
 PARAMETER:  3.0160E-02  1.5128E-01 -7.9799E-02  5.8453E-02  5.3418E-02  3.6533E-02  2.9643E-01 -2.1712E-01 -1.5236E-01  7.9625E-02
             1.7579E-01
 GRADIENT:  -4.6988E+00  8.6414E-01 -1.7512E+00  2.8599E+00  2.6655E+00 -2.3344E+00  2.6788E-01  4.9489E-01  8.1482E-01  5.7726E-01
             2.3831E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1625.26632969757        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  9.3593E-01  1.0338E+00  7.1211E-01  9.5726E-01  8.6891E-01  9.4698E-01  1.2358E+00  4.4261E-01  7.6295E-01  9.0865E-01
             1.0691E+00
 PARAMETER:  3.3785E-02  1.3329E-01 -2.3953E-01  5.6321E-02 -4.0513E-02  4.5523E-02  3.1173E-01 -7.1506E-01 -1.7056E-01  4.2032E-03
             1.6680E-01
 GRADIENT:   2.4303E+00 -1.9876E+00 -1.4072E+00 -1.2025E+00 -2.7212E-01  9.3779E-01  1.0654E-01  5.1712E-01 -5.7631E-02  2.9679E-01
            -5.2955E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1625.55630835312        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      460
 NPARAMETR:  9.4809E-01  1.0468E+00  7.2955E-01  9.5589E-01  8.8695E-01  9.5115E-01  1.2355E+00  4.1181E-01  7.7134E-01  9.3338E-01
             1.0733E+00
 PARAMETER:  4.6698E-02  1.4576E-01 -2.1533E-01  5.4890E-02 -1.9969E-02  4.9916E-02  3.1146E-01 -7.8718E-01 -1.5962E-01  3.1058E-02
             1.7075E-01
 GRADIENT:   6.7362E-01 -2.7209E-01 -2.4265E-01  2.7754E-01 -5.1835E-02 -1.7222E-01  1.5318E-01  4.6083E-02  7.0499E-02  2.1461E-01
             2.2968E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1625.57699350352        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      635
 NPARAMETR:  9.4834E-01  1.1110E+00  6.7931E-01  9.1316E-01  8.8887E-01  9.5225E-01  1.1764E+00  2.7082E-01  7.9307E-01  9.2645E-01
             1.0734E+00
 PARAMETER:  4.6953E-02  2.0522E-01 -2.8668E-01  9.1563E-03 -1.7809E-02  5.1076E-02  2.6246E-01 -1.2063E+00 -1.3184E-01  2.3601E-02
             1.7085E-01
 GRADIENT:  -1.2473E-01 -4.9330E-01 -3.2324E-01 -3.0221E-01  4.1664E-01  4.4376E-03  4.9600E-02  2.7141E-02  4.2212E-02  9.8114E-02
             1.6701E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1625.58349830878        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      813
 NPARAMETR:  9.4857E-01  1.1165E+00  6.6039E-01  9.0822E-01  8.7913E-01  9.5257E-01  1.1713E+00  1.3157E-01  7.9475E-01  9.1851E-01
             1.0731E+00
 PARAMETER:  4.7199E-02  2.1019E-01 -3.1493E-01  3.7335E-03 -2.8819E-02  5.1412E-02  2.5810E-01 -1.9282E+00 -1.2972E-01  1.5000E-02
             1.7058E-01
 GRADIENT:   3.3406E-02  7.1952E-02  6.0079E-02  7.0436E-02 -8.6901E-02  4.8743E-02 -1.2749E-02  2.6256E-03 -2.7588E-02  1.0619E-02
             4.8366E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1625.58512871160        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      988
 NPARAMETR:  9.4861E-01  1.1245E+00  6.5298E-01  9.0269E-01  8.7865E-01  9.5265E-01  1.1643E+00  3.7962E-02  7.9775E-01  9.1705E-01
             1.0731E+00
 PARAMETER:  4.7242E-02  2.1734E-01 -3.2621E-01 -2.3799E-03 -2.9366E-02  5.1489E-02  2.5212E-01 -3.1712E+00 -1.2595E-01  1.3407E-02
             1.7051E-01
 GRADIENT:  -3.5581E-02 -2.9445E-02 -4.0603E-03 -9.9948E-03  2.4417E-02  4.2765E-02 -1.4579E-02  2.5573E-04 -1.0914E-02  8.5267E-03
             1.9701E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1625.58525919415        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1163
 NPARAMETR:  9.4863E-01  1.1250E+00  6.5225E-01  9.0234E-01  8.7840E-01  9.5256E-01  1.1640E+00  1.0000E-02  7.9796E-01  9.1670E-01
             1.0730E+00
 PARAMETER:  4.7262E-02  2.1776E-01 -3.2732E-01 -2.7630E-03 -2.9651E-02  5.1396E-02  2.5183E-01 -4.8332E+00 -1.2569E-01  1.3021E-02
             1.7046E-01
 GRADIENT:  -3.1582E-03 -5.2260E-03 -8.5164E-03  1.3014E-03  3.7814E-03  3.1412E-03 -2.1329E-04  0.0000E+00  2.2409E-04  5.2498E-03
             4.0419E-03

0ITERATION NO.:   47    OBJECTIVE VALUE:  -1625.58525943577        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1220
 NPARAMETR:  9.4863E-01  1.1249E+00  6.5233E-01  9.0240E-01  8.7841E-01  9.5255E-01  1.1640E+00  1.0000E-02  7.9794E-01  9.1669E-01
             1.0730E+00
 PARAMETER:  4.7262E-02  2.1769E-01 -3.2720E-01 -2.7024E-03 -2.9646E-02  5.1388E-02  2.5190E-01 -4.8217E+00 -1.2572E-01  1.3015E-02
             1.7045E-01
 GRADIENT:   5.9815E-05  2.6004E-03  6.9330E-03 -2.3775E-03 -2.5670E-03 -5.1359E-05 -8.7055E-04  0.0000E+00 -1.1429E-03 -2.7730E-03
            -1.6500E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1220
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.0972E-04 -3.5862E-03 -4.1207E-04 -2.2338E-03 -1.5431E-02
 SE:             2.9804E-02  2.3695E-02  1.6922E-04  2.2017E-02  2.2685E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9706E-01  8.7970E-01  1.4887E-02  9.1919E-01  4.9637E-01

 ETASHRINKSD(%)  1.5426E-01  2.0618E+01  9.9433E+01  2.6242E+01  2.4002E+01
 ETASHRINKVR(%)  3.0829E-01  3.6985E+01  9.9997E+01  4.5597E+01  4.2243E+01
 EBVSHRINKSD(%)  5.2093E-01  2.0255E+01  9.9502E+01  2.7053E+01  2.2831E+01
 EBVSHRINKVR(%)  1.0391E+00  3.6407E+01  9.9998E+01  4.6787E+01  4.0450E+01
 RELATIVEINF(%)  9.8721E+01  3.9418E+00  2.4174E-04  3.0974E+00  5.9162E+00
 EPSSHRINKSD(%)  4.3112E+01
 EPSSHRINKVR(%)  6.7638E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1625.5852594357741     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -890.43443287203593     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.29
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1625.585       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.49E-01  1.12E+00  6.52E-01  9.02E-01  8.78E-01  9.53E-01  1.16E+00  1.00E-02  7.98E-01  9.17E-01  1.07E+00
 


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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.35E+03
 
 TH 2
+       -7.24E+00  4.24E+02
 
 TH 3
+        2.00E+01  1.86E+02  7.27E+02
 
 TH 4
+       -1.93E+01  3.90E+02 -4.54E+02  1.16E+03
 
 TH 5
+       -5.73E+00 -2.81E+02 -7.34E+02  4.47E+02  1.02E+03
 
 TH 6
+        5.22E-01 -2.47E+00  3.12E+00 -5.40E+00 -1.04E+00  2.20E+02
 
 TH 7
+        1.38E+00  2.72E+01 -2.16E+01 -1.61E+01  4.27E+00  6.63E-02  6.46E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.63E-01 -1.85E+01 -3.72E+01  3.28E+01  8.82E+00 -2.93E-01  2.00E+01  0.00E+00  1.12E+02
 
 TH10
+       -5.06E+00 -1.12E+01 -5.50E+01 -1.58E+01 -5.77E+01 -9.00E-01  1.07E+01  0.00E+00  1.20E+01  9.27E+01
 
 TH11
+       -7.10E+00 -1.54E+01 -3.76E+01 -3.44E-02  3.45E+00  3.38E+00  6.99E+00  0.00E+00  1.41E+01  1.93E+01  1.89E+02
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       19.072
Stop Time:
Sat Sep 18 11:02:09 CDT 2021
