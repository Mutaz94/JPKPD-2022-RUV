Thu Sep 30 01:39:41 CDT 2021
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
$DATA ../../../../data/spa1/TD1/dat81.csv ignore=@
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
$COVARIANCE UNCOND

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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1623.10764295573        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9773E+02 -3.4772E+01  4.1441E+01 -7.2416E-01  7.1622E+01  2.9298E+01 -9.1386E-01 -2.0300E+02 -7.4250E+00  8.1922E+00
            -7.7878E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2099.19016490671        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  9.1599E-01  1.0398E+00  9.8586E-01  1.0121E+00  9.5079E-01  1.0858E+00  9.9367E-01  1.1537E+00  9.7710E-01  9.6373E-01
             1.3489E+00
 PARAMETER:  1.2253E-02  1.3904E-01  8.5755E-02  1.1200E-01  4.9537E-02  1.8229E-01  9.3646E-02  2.4300E-01  7.6835E-02  6.3052E-02
             3.9930E-01
 GRADIENT:   4.5390E+01  2.4647E+01 -9.3207E+00  4.6122E+01 -1.1046E+01  9.1318E+01  5.1441E+00  7.8355E+00  1.9295E+01  1.9269E+01
             1.8560E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2102.57164481586        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      210
 NPARAMETR:  9.1875E-01  9.6030E-01  9.9428E-01  1.0633E+00  9.1728E-01  1.0754E+00  1.1344E+00  9.7483E-01  9.2421E-01  8.1609E-01
             1.3446E+00
 PARAMETER:  1.5255E-02  5.9494E-02  9.4261E-02  1.6139E-01  1.3661E-02  1.7273E-01  2.2614E-01  7.4504E-02  2.1181E-02 -1.0323E-01
             3.9609E-01
 GRADIENT:  -2.1279E+02 -1.0307E+00  1.6087E+00 -6.6645E+00 -8.2784E+00 -3.0436E+01  2.9906E+00 -6.7293E-01  1.5399E+01  3.1213E+00
             1.7567E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2132.85947834827        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      387
 NPARAMETR:  9.6239E-01  9.3486E-01  9.3975E-01  1.0747E+00  8.9487E-01  1.0234E+00  1.2728E+00  7.5188E-01  8.7504E-01  9.8131E-01
             1.0689E+00
 PARAMETER:  6.1667E-02  3.2643E-02  3.7854E-02  1.7200E-01 -1.1082E-02  1.2310E-01  3.4118E-01 -1.8518E-01 -3.3483E-02  8.1137E-02
             1.6665E-01
 GRADIENT:  -1.3211E+02  6.1980E+00  3.6940E-02 -5.2731E+00 -6.9689E+00 -3.7438E+01  9.8336E+00 -3.6089E+00  1.0240E+01  1.5376E+01
             4.0616E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2136.12932172778        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  9.7218E-01  9.3767E-01  9.5441E-01  1.0746E+00  9.0180E-01  1.0500E+00  1.2518E+00  8.2243E-01  8.5843E-01  9.6063E-01
             1.0590E+00
 PARAMETER:  7.1788E-02  3.5640E-02  5.3339E-02  1.7194E-01 -3.3656E-03  1.4874E-01  3.2458E-01 -9.5496E-02 -5.2654E-02  5.9836E-02
             1.5737E-01
 GRADIENT:  -1.0507E+02  5.9016E+00 -2.8073E-01 -4.9957E+00 -5.5608E+00 -2.3089E+01  7.5103E+00 -2.7208E+00  5.9431E+00  1.2805E+01
             3.3555E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2137.27228290551        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      677
 NPARAMETR:  9.7262E-01  9.3809E-01  9.5484E-01  1.0754E+00  9.0139E-01  1.1084E+00  1.2536E+00  8.2206E-01  8.5881E-01  8.1488E-01
             1.0598E+00
 PARAMETER:  7.2239E-02  3.6091E-02  5.3790E-02  1.7271E-01 -3.8164E-03  2.0290E-01  3.2604E-01 -9.5946E-02 -5.2203E-02 -1.0472E-01
             1.5808E-01
 GRADIENT:   3.2342E+02  3.8395E+01  7.7587E+00  1.6475E+02  2.0216E+01  1.7722E+02  1.8241E+01 -7.1067E+00  1.0686E+01 -6.3060E+00
             3.0320E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2137.35230909378        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      833
 NPARAMETR:  9.7262E-01  9.3809E-01  9.5484E-01  1.0754E+00  9.0139E-01  1.1050E+00  1.2536E+00  8.2206E-01  8.5881E-01  8.2659E-01
             1.0598E+00
 PARAMETER:  7.2239E-02  3.6091E-02  5.3790E-02  1.7271E-01 -3.8168E-03  1.9980E-01  3.2604E-01 -9.5947E-02 -5.2203E-02 -9.0447E-02
             1.5808E-01
 GRADIENT:  -9.4258E+01  6.1168E+00  4.9962E+00  2.8469E-01  1.0684E+01 -9.4870E-01  5.2870E+00 -6.8271E+00  4.9996E+00 -5.2487E+00
             2.9331E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2138.11090729692        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1025
 NPARAMETR:  9.7739E-01  9.3732E-01  9.5454E-01  1.0754E+00  9.0056E-01  1.1036E+00  1.2460E+00  8.3951E-01  8.5881E-01  8.2887E-01
             1.0553E+00
 PARAMETER:  7.7126E-02  3.5267E-02  5.3476E-02  1.7270E-01 -4.7357E-03  1.9857E-01  3.1994E-01 -7.4943E-02 -5.2213E-02 -8.7695E-02
             1.5384E-01
 GRADIENT:  -8.5573E+01  5.1742E+00  3.8017E+00 -4.9826E-01  9.9378E+00 -5.8665E-01  4.8323E+00 -6.0953E+00  4.8589E+00 -4.4446E+00
             2.6566E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2138.65237972487        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1223
 NPARAMETR:  9.7763E-01  9.3406E-01  9.5478E-01  1.0756E+00  8.9790E-01  1.1054E+00  1.1958E+00  8.3930E-01  8.4758E-01  8.2866E-01
             1.0369E+00
 PARAMETER:  7.7376E-02  3.1789E-02  5.3726E-02  1.7291E-01 -7.6910E-03  2.0019E-01  2.7883E-01 -7.5193E-02 -6.5372E-02 -8.7945E-02
             1.3621E-01
 GRADIENT:  -8.4619E+01  3.2702E+00  6.1562E+00 -3.8826E+00  7.4406E+00  1.7226E-01  3.3999E-02 -6.7380E+00  1.9412E-01 -5.7372E+00
             1.2931E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2138.69516150993        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1400
 NPARAMETR:  9.7774E-01  9.3417E-01  9.5489E-01  1.0756E+00  8.9779E-01  1.1060E+00  1.1951E+00  8.4255E-01  8.4768E-01  8.2856E-01
             1.0363E+00
 PARAMETER:  7.7493E-02  3.1906E-02  5.3843E-02  1.7291E-01 -7.8161E-03  2.0072E-01  2.7827E-01 -7.1320E-02 -6.5255E-02 -8.8062E-02
             1.3567E-01
 GRADIENT:  -8.4312E+01  3.3312E+00  6.0669E+00 -3.8444E+00  7.1690E+00  4.1094E-01  1.8982E-02 -6.6256E+00  1.9565E-01 -5.6433E+00
             1.2563E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2138.75157199328        NO. OF FUNC. EVALS.: 205
 CUMULATIVE NO. OF FUNC. EVALS.:     1605             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7798E-01  9.3417E-01  9.5512E-01  1.0756E+00  8.9779E-01  1.1092E+00  1.1951E+00  8.4607E-01  8.4768E-01  8.2836E-01
             1.0356E+00
 PARAMETER:  7.7733E-02  3.1906E-02  5.4083E-02  1.7291E-01 -7.8161E-03  2.0367E-01  2.7824E-01 -6.7148E-02 -6.5255E-02 -8.8302E-02
             1.3494E-01
 GRADIENT:   3.5362E+02  3.7073E+01  8.4130E+00  1.6758E+02  1.5163E+01  1.8507E+02  1.0986E+01 -6.3300E+00  6.2946E+00 -5.0027E+00
             1.3380E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2138.77552043569        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1773
 NPARAMETR:  9.7780E-01  9.3417E-01  9.5494E-01  1.0756E+00  8.9779E-01  1.1088E+00  1.1947E+00  8.4623E-01  8.4768E-01  8.2852E-01
             1.0337E+00
 PARAMETER:  7.7765E-02  3.1906E-02  5.4115E-02  1.7291E-01 -7.8161E-03  2.0123E-01  2.7761E-01 -6.7180E-02 -6.5255E-02 -8.8333E-02
             1.3286E-01
 GRADIENT:   2.6115E+05 -1.0885E+01  1.3065E+05 -4.3563E+01  1.7525E+01 -9.2140E-01 -5.7142E-02 -1.3069E+05 -3.6856E+01 -1.3062E+05
            -1.9700E+05
 NUMSIGDIG:         2.3         7.0         2.3         6.1         6.8         1.6         2.7         2.3         6.5         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1773
 NO. OF SIG. DIGITS IN FINAL EST.:  1.6
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.2953E-02 -4.7514E-03 -2.8634E-02 -4.9552E-05 -2.6161E-02
 SE:             2.9674E-02  2.0022E-02  1.4628E-02  2.4532E-02  2.1795E-02
 N:                     100         100         100         100         100

 P VAL.:         1.4776E-01  8.1242E-01  5.0284E-02  9.9839E-01  2.3003E-01

 ETASHRINKSD(%)  5.8866E-01  3.2925E+01  5.0996E+01  1.7815E+01  2.6983E+01
 ETASHRINKVR(%)  1.1739E+00  5.5009E+01  7.5986E+01  3.2456E+01  4.6685E+01
 EBVSHRINKSD(%)  2.9816E-01  3.3394E+01  5.6777E+01  1.8064E+01  2.7427E+01
 EBVSHRINKVR(%)  5.9542E-01  5.5637E+01  8.1318E+01  3.2866E+01  4.7331E+01
 RELATIVEINF(%)  9.8760E+01  2.2055E+00  2.0092E+00  3.7432E+00  8.5667E+00
 EPSSHRINKSD(%)  3.4642E+01
 EPSSHRINKVR(%)  5.7284E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2138.7755204356854     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1219.8369872310127     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2138.776       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.78E-01  9.34E-01  9.55E-01  1.08E+00  8.98E-01  1.11E+00  1.19E+00  8.46E-01  8.48E-01  8.28E-01  1.03E+00
 


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
+        1.37E+08
 
 TH 2
+       -8.62E+00  1.50E+08
 
 TH 3
+        1.40E+08  7.32E+07  1.43E+08
 
 TH 4
+       -3.53E-01  4.12E+02 -1.89E+02  3.78E+07
 
 TH 5
+        1.07E+00 -1.56E+08 -4.83E+02  2.13E+02  1.62E+08
 
 TH 6
+        8.51E+02 -2.45E+00  7.94E+02  6.48E-01 -7.80E-01  1.62E+02
 
 TH 7
+       -5.37E+02  2.96E+01 -2.06E+07 -3.31E+00 -1.26E+01  1.81E-01  3.22E+01
 
 TH 8
+       -7.90E+07 -1.88E+01 -8.08E+07  4.15E+07 -4.12E+00 -8.96E+02  6.24E+02  1.82E+08
 
 TH 9
+        7.88E+07  1.65E+08  8.07E+07  1.60E+01 -1.72E+08  4.45E+00  2.48E+01  1.17E+02  1.82E+08
 
 TH10
+       -1.61E+08 -8.44E+07 -1.65E+08 -1.72E+01 -9.30E+01 -9.13E+02  6.46E+02  9.32E+07 -9.30E+07  1.90E+08
 
 TH11
+       -4.87E+07 -8.77E+00 -4.99E+07  2.56E+07  8.00E+00 -5.50E+02  3.86E+02  5.35E+04  1.73E+01  5.75E+07  3.48E+07
 
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
 #CPUT: Total CPU Time in Seconds,       37.238
Stop Time:
Thu Sep 30 01:40:20 CDT 2021
