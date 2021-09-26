Sat Sep 25 09:37:13 CDT 2021
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
$DATA ../../../../data/spa/A3/dat96.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -29.4639560207362        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.6370E+02  1.6515E+01  6.0661E+01 -5.7784E+01  1.0024E+02  6.0743E+00 -1.4926E+01 -3.5255E+01 -7.2653E+01 -5.8936E+01
            -2.9715E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1205.94820626877        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.6809E-01  1.0102E+00  9.5118E-01  1.1347E+00  9.7300E-01  7.8674E-01  8.9806E-01  1.0068E+00  9.3839E-01  9.1068E-01
             5.3750E+00
 PARAMETER:  6.7568E-02  1.1018E-01  4.9952E-02  2.2636E-01  7.2627E-02 -1.3986E-01 -7.5216E-03  1.0676E-01  3.6409E-02  6.4400E-03
             1.7818E+00
 GRADIENT:  -1.0469E+02  8.3361E+00 -1.6244E+01  3.3994E+01 -4.8783E+00 -2.4309E+01  1.3964E+01  6.4900E+00  2.8018E+01  1.8023E+01
             2.3650E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1240.99255132454        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      164
 NPARAMETR:  9.6229E-01  8.3372E-01  3.1358E-01  1.1158E+00  4.5779E-01  8.9617E-01  2.9548E-01  5.0686E-01  7.9047E-01  3.9789E-01
             4.7595E+00
 PARAMETER:  6.1556E-02 -8.1863E-02 -1.0597E+00  2.0955E-01 -6.8134E-01 -9.6283E-03 -1.1192E+00 -5.7951E-01 -1.3513E-01 -8.2158E-01
             1.6601E+00
 GRADIENT:  -9.3454E+01  4.5688E+01 -6.2730E+00  9.6629E+01 -2.6137E+01 -4.9322E+00  1.1587E+00  4.0883E+00  9.4358E+00  9.2335E+00
             1.6508E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1273.84372724878        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  9.7481E-01  6.1065E-01  3.8215E-01  1.1801E+00  4.4154E-01  9.1197E-01  6.1015E-01  2.3034E-01  8.0276E-01  8.7422E-02
             3.6598E+00
 PARAMETER:  7.4491E-02 -3.9324E-01 -8.6194E-01  2.6563E-01 -7.1748E-01  7.8506E-03 -3.9404E-01 -1.3682E+00 -1.1970E-01 -2.3370E+00
             1.3974E+00
 GRADIENT:  -1.6627E+00  5.4730E+00  1.7264E-01  1.0474E+01 -6.1139E-01  8.8414E-03 -5.0886E-01  2.8794E-01  1.0989E+00  6.5610E-02
            -1.0754E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1274.11781974125        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  9.7607E-01  5.0953E-01  2.7989E-01  1.1624E+00  3.4192E-01  9.0452E-01  7.8089E-01  1.3740E-01  7.6991E-01  6.1691E-02
             3.6289E+00
 PARAMETER:  7.5780E-02 -5.7427E-01 -1.1733E+00  2.5047E-01 -9.7317E-01 -3.5473E-04 -1.4732E-01 -1.8849E+00 -1.6148E-01 -2.6856E+00
             1.3889E+00
 GRADIENT:  -2.1184E+00  5.9237E+00 -2.7868E+00  2.5660E+01 -1.1904E+00 -7.6231E+00  3.4584E-01 -7.4335E-02 -5.2254E+00 -2.5532E-02
            -6.0632E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1274.26514945450        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      376
 NPARAMETR:  9.7447E-01  4.8053E-01  2.2870E-01  1.1151E+00  2.9998E-01  9.1856E-01  7.6285E-01  9.2554E-02  8.1535E-01  5.2850E-02
             3.5865E+00
 PARAMETER:  7.4136E-02 -6.3286E-01 -1.3754E+00  2.0898E-01 -1.1040E+00  1.5049E-02 -1.7069E-01 -2.2800E+00 -1.0414E-01 -2.8403E+00
             1.3772E+00
 GRADIENT:  -1.1301E+00  3.0399E+00  1.8595E+00  1.3718E+01 -7.8011E+00 -5.6691E+00  1.6788E-02 -9.5076E-02 -4.6199E+00 -7.5358E-02
            -8.8644E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1274.38203155392        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  9.7062E-01  4.9468E-01  1.9501E-01  1.0610E+00  2.8189E-01  9.3874E-01  7.1961E-01  7.7895E-02  8.8150E-01  5.1288E-02
             3.5401E+00
 PARAMETER:  7.0182E-02 -6.0384E-01 -1.5347E+00  1.5920E-01 -1.1662E+00  3.6780E-02 -2.2905E-01 -2.4524E+00 -2.6126E-02 -2.8703E+00
             1.3642E+00
 GRADIENT:   3.1758E-01 -8.3032E-01  3.1166E+00 -1.6331E+00 -5.3483E+00 -6.1034E-01 -1.2994E-01 -8.4172E-02 -9.9762E-01 -8.7895E-02
             5.1594E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1274.82867584837        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      521
 NPARAMETR:  9.6539E-01  5.7249E-01  1.8645E-01  1.0277E+00  2.9850E-01  9.4258E-01  6.6215E-01  1.3794E-01  9.1259E-01  6.7129E-02
             3.5000E+00
 PARAMETER:  6.4778E-02 -4.5776E-01 -1.5796E+00  1.2729E-01 -1.1090E+00  4.0871E-02 -3.1226E-01 -1.8810E+00  8.5260E-03 -2.6011E+00
             1.3528E+00
 GRADIENT:  -9.1088E-02 -7.9159E-01 -2.3883E+00  1.6753E+00  4.4604E+00  9.8655E-01  4.8593E-02 -2.0026E-01  2.1317E-01 -5.1501E-02
             9.5603E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1275.61845007007        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      600
 NPARAMETR:  9.6487E-01  5.5559E-01  1.7572E-01  1.0203E+00  2.8228E-01  9.4629E-01  6.8008E-01  6.5349E-01  9.3456E-01  7.3006E-02
             3.3671E+00
 PARAMETER:  6.4235E-02 -4.8773E-01 -1.6389E+00  1.2010E-01 -1.1649E+00  4.4789E-02 -2.8554E-01 -3.2543E-01  3.2321E-02 -2.5172E+00
             1.3141E+00
 GRADIENT:  -1.4938E+00  1.7216E+01  1.0566E+01  4.1601E+00 -2.5248E+01  3.6222E-01 -2.6844E-01 -1.0082E+00 -2.3055E+00  6.4896E-02
            -2.9309E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1275.66167684468        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      670
 NPARAMETR:  9.6336E-01  5.3219E-01  1.6370E-01  1.0101E+00  2.6830E-01  9.5249E-01  6.7039E-01  8.2990E-01  9.8047E-01  6.4851E-02
             3.3015E+00
 PARAMETER:  6.2672E-02 -5.3075E-01 -1.7097E+00  1.1004E-01 -1.2156E+00  5.1327E-02 -2.9989E-01 -8.6451E-02  8.0275E-02 -2.6357E+00
             1.2944E+00
 GRADIENT:  -8.0186E-01  1.1497E+01  8.2532E+00  2.8334E+00 -1.7944E+01  1.4258E+00 -8.2149E-02 -1.0021E-01 -1.3125E+00  7.1386E-02
            -1.6662E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1275.82807521934        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  9.6422E-01  5.6209E-01  1.7099E-01  1.0103E+00  2.8360E-01  9.4463E-01  6.3882E-01  7.5421E-01  9.6224E-01  6.7155E-02
             3.3596E+00
 PARAMETER:  6.3563E-02 -4.7609E-01 -1.6662E+00  1.1021E-01 -1.1602E+00  4.3040E-02 -3.4814E-01 -1.8209E-01  6.1505E-02 -2.6008E+00
             1.3118E+00
 GRADIENT:   2.8305E-02 -1.4341E+00 -9.7113E-01 -5.3872E-01  2.3387E+00 -1.8595E-01  1.0484E-01 -4.2694E-02  1.4749E-01  1.0171E-01
             6.8406E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1275.87905299307        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  9.6365E-01  5.7122E-01  1.7055E-01  1.0071E+00  2.8574E-01  9.4503E-01  6.3739E-01  7.8724E-01  9.6429E-01  1.3981E-02
             3.3479E+00
 PARAMETER:  6.2971E-02 -4.5998E-01 -1.6687E+00  1.0711E-01 -1.1527E+00  4.3461E-02 -3.5037E-01 -1.3922E-01  6.3635E-02 -4.1700E+00
             1.3083E+00
 GRADIENT:   9.5330E-02  1.5293E-01  5.4967E-02 -4.6524E-02  2.5570E-02  8.7190E-02 -2.1102E-02  2.0105E-02 -1.0372E-01  4.6854E-03
            -2.3573E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1275.88091413815        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1188
 NPARAMETR:  9.6379E-01  5.6768E-01  1.7042E-01  1.0081E+00  2.8464E-01  9.4499E-01  6.4008E-01  7.8351E-01  9.6472E-01  1.0000E-02
             3.3491E+00
 PARAMETER:  6.3113E-02 -4.6619E-01 -1.6695E+00  1.0807E-01 -1.1565E+00  4.3423E-02 -3.4617E-01 -1.4397E-01  6.4079E-02 -4.6107E+00
             1.3087E+00
 GRADIENT:  -2.2045E-03  4.9369E-03 -6.9797E-03  1.1245E-02 -4.3060E-03 -8.5992E-03  8.3594E-04  1.3861E-04 -8.8759E-03  0.0000E+00
            -2.1306E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1188
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.7333E-04 -5.3503E-03  3.8045E-03 -1.3264E-02  2.2751E-04
 SE:             2.8706E-02  1.4190E-02  1.2210E-02  2.3628E-02  3.6220E-04
 N:                     100         100         100         100         100

 P VAL.:         9.7295E-01  7.0614E-01  7.5536E-01  5.7454E-01  5.2991E-01

 ETASHRINKSD(%)  3.8308E+00  5.2462E+01  5.9094E+01  2.0843E+01  9.8787E+01
 ETASHRINKVR(%)  7.5148E+00  7.7401E+01  8.3267E+01  3.7342E+01  9.9985E+01
 EBVSHRINKSD(%)  3.6907E+00  5.2495E+01  5.8487E+01  1.9870E+01  9.8761E+01
 EBVSHRINKVR(%)  7.2451E+00  7.7433E+01  8.2766E+01  3.5792E+01  9.9985E+01
 RELATIVEINF(%)  8.4174E+01  4.9562E-01  1.5738E+00  2.1118E+01  2.3279E-04
 EPSSHRINKSD(%)  2.6134E+01
 EPSSHRINKVR(%)  4.5438E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1275.8809141381485     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -540.73008757441028     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.79
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1275.881       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.64E-01  5.68E-01  1.70E-01  1.01E+00  2.85E-01  9.45E-01  6.40E-01  7.84E-01  9.65E-01  1.00E-02  3.35E+00
 


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
+        1.25E+03
 
 TH 2
+       -9.24E+01  1.85E+03
 
 TH 3
+       -4.18E+02  3.27E+03  1.06E+04
 
 TH 4
+       -7.39E+01  2.74E+02 -9.53E+02  7.52E+02
 
 TH 5
+        4.69E+02 -6.03E+03 -1.28E+04 -1.62E+02  2.11E+04
 
 TH 6
+        8.45E-01 -1.12E+01  3.18E+01 -2.10E+01  5.40E+01  1.93E+02
 
 TH 7
+       -2.75E+00 -4.89E+01 -1.40E+02 -3.60E+00  2.15E+02  1.39E+00  2.69E+01
 
 TH 8
+       -9.61E-01 -3.27E+01 -8.80E+00 -3.00E+00  9.81E+01  7.50E-01  5.07E+00  1.20E+01
 
 TH 9
+        6.69E+00 -3.24E+01  9.73E+01 -2.31E+01  1.07E+02  2.02E+00  1.05E+01 -1.76E+00  7.98E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.69E+01 -9.86E+00 -2.46E+01 -5.93E+00  3.00E+01  3.36E+00  1.11E+01  8.61E+00  1.29E+01  0.00E+00  2.96E+01
 
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
 #CPUT: Total CPU Time in Seconds,       17.956
Stop Time:
Sat Sep 25 09:37:32 CDT 2021
