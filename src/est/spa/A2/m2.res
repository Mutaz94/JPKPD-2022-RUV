Sat Sep 18 09:35:27 CDT 2021
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
$DATA ../../../../data/spa/A2/dat2.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -972.953451100235        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0411E+02  2.3874E+01  1.3802E+01  1.1088E+01  7.9555E+01  1.4479E+01 -3.1176E+01  3.5303E+00 -3.7752E+01 -8.8595E+01
            -1.2127E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1359.31969238265        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0175E+00  8.4792E-01  1.1836E+00  1.1061E+00  1.0082E+00  9.7883E-01  1.1037E+00  7.3491E-01  1.1301E+00  1.3753E+00
             1.8248E+00
 PARAMETER:  1.1736E-01 -6.4970E-02  2.6859E-01  2.0086E-01  1.0819E-01  7.8601E-02  1.9870E-01 -2.0801E-01  2.2229E-01  4.1867E-01
             7.0145E-01
 GRADIENT:   1.0716E+02  1.1231E+00 -5.9592E-01  2.2945E+00  5.1507E+00  5.8345E+00 -1.3002E+00  3.3662E+00  1.5629E+01  6.0937E+00
            -2.0162E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1373.77767629061        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.9650E-01  7.0198E-01  6.0094E-01  1.1819E+00  6.4556E-01  8.9451E-01  1.1278E+00  6.6579E-02  9.9489E-01  7.9994E-01
             1.9858E+00
 PARAMETER:  9.6490E-02 -2.5385E-01 -4.0925E-01  2.6715E-01 -3.3763E-01 -1.1485E-02  2.2031E-01 -2.6094E+00  9.4877E-02 -1.2321E-01
             7.8603E-01
 GRADIENT:   4.0314E+01 -1.0816E+01 -8.1674E+01  7.6736E+01  1.2706E+02 -2.8518E+01 -6.9170E+00  1.0019E-01  1.9228E+00 -1.6423E+01
            -1.3714E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1399.67460160499        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9667E-01  4.6280E-01  4.5198E-01  1.2372E+00  4.3698E-01  9.5703E-01  1.5375E+00  1.3129E-02  8.7508E-01  5.6817E-01
             2.5135E+00
 PARAMETER:  9.6666E-02 -6.7047E-01 -6.9413E-01  3.1288E-01 -7.2786E-01  5.6083E-02  5.3017E-01 -4.2330E+00 -3.3444E-02 -4.6533E-01
             1.0217E+00
 GRADIENT:   1.3612E+01  1.5231E+01 -1.0113E+00  2.4159E+01  4.4009E+00 -3.0228E-01  9.8463E-01  3.7416E-03 -1.1738E+00 -8.7914E-01
             1.3387E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1401.66233257972        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  9.8443E-01  3.0743E-01  3.7064E-01  1.2367E+00  3.5129E-01  9.6961E-01  2.0347E+00  1.0000E-02  8.8932E-01  5.0499E-01
             2.4410E+00
 PARAMETER:  8.4311E-02 -1.0795E+00 -8.9251E-01  3.1244E-01 -9.4615E-01  6.9136E-02  8.1035E-01 -5.2195E+00 -1.7293E-02 -5.8322E-01
             9.9240E-01
 GRADIENT:  -9.9034E+00  6.3088E+00  1.4789E+01 -1.2141E+01 -1.8740E+01  3.8202E+00  3.9497E+00  0.0000E+00  2.4337E+00 -3.0335E+00
             1.6184E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1402.93594110121        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  9.8629E-01  2.3789E-01  3.9602E-01  1.2877E+00  3.5782E-01  9.5668E-01  2.2571E+00  1.0000E-02  8.7061E-01  5.7355E-01
             2.4285E+00
 PARAMETER:  8.6191E-02 -1.3359E+00 -8.2629E-01  3.5287E-01 -9.2774E-01  5.5712E-02  9.1409E-01 -5.3382E+00 -3.8564E-02 -4.5591E-01
             9.8729E-01
 GRADIENT:   8.2389E-01  1.5689E+00 -1.3901E+00  5.7422E-01 -1.1940E-01 -1.9599E-01 -1.9813E-02  0.0000E+00 -6.8511E-01  4.4543E-01
            -7.9606E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1404.42142640562        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  9.7527E-01  8.7173E-02  4.4562E-01  1.3868E+00  3.6682E-01  9.5208E-01  4.4551E+00  1.0000E-02  8.4406E-01  5.9076E-01
             2.4327E+00
 PARAMETER:  7.4959E-02 -2.3399E+00 -7.0830E-01  4.2703E-01 -9.0290E-01  5.0898E-02  1.5940E+00 -6.5034E+00 -6.9527E-02 -4.2635E-01
             9.8899E-01
 GRADIENT:  -5.5812E+00  3.7418E+00  2.6449E+01  2.2407E+01 -3.9750E+01 -1.4689E+00  3.8482E+00  0.0000E+00  3.4588E+00 -1.7837E+00
            -2.5862E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1404.71862906506        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      522
 NPARAMETR:  9.7749E-01  9.1241E-02  4.3450E-01  1.3672E+00  3.6513E-01  9.5527E-01  4.1909E+00  1.0000E-02  8.3217E-01  5.8522E-01
             2.4409E+00
 PARAMETER:  7.7236E-02 -2.2942E+00 -7.3356E-01  4.1273E-01 -9.0750E-01  5.4241E-02  1.5329E+00 -6.5269E+00 -8.3713E-02 -4.3578E-01
             9.9237E-01
 GRADIENT:  -2.1561E-01  1.1133E+00  3.6368E+00 -2.0503E-01 -4.2166E+00 -1.6571E-01  9.0797E-01  0.0000E+00 -1.7083E-01 -5.2771E-01
            -1.7897E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1405.07855406882        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      680
 NPARAMETR:  9.7845E-01  7.5246E-02  4.6669E-01  1.3967E+00  3.8504E-01  9.5808E-01  4.7837E+00  1.0000E-02  8.2223E-01  6.0325E-01
             2.4424E+00
 PARAMETER:  7.8215E-02 -2.4870E+00 -6.6209E-01  4.3408E-01 -8.5440E-01  5.7178E-02  1.6652E+00 -6.6188E+00 -9.5732E-02 -4.0543E-01
             9.9297E-01
 GRADIENT:  -3.2798E-01  9.7288E-01 -4.4659E+00  4.5929E-02  4.5314E+00  8.2725E-01  1.2480E+00  0.0000E+00 -1.6493E+00  4.7876E-01
             6.4830E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1405.16871453112        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      857
 NPARAMETR:  9.7710E-01  4.8716E-02  4.6613E-01  1.4057E+00  3.8115E-01  9.5736E-01  5.8272E+00  1.0000E-02  8.1778E-01  5.9613E-01
             2.4422E+00
 PARAMETER:  7.6832E-02 -2.9218E+00 -6.6329E-01  4.4056E-01 -8.6457E-01  5.6424E-02  1.8625E+00 -7.3025E+00 -1.0116E-01 -4.1729E-01
             9.9289E-01
 GRADIENT:  -1.0427E+00 -1.1291E+00  4.6737E+00  5.5728E+00 -8.0031E+00  8.9832E-01 -3.3633E+00  0.0000E+00 -6.3657E-01  1.5970E+00
             1.4800E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1406.41190190645        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1034
 NPARAMETR:  9.7517E-01  2.7036E-02  4.6052E-01  1.4134E+00  3.7489E-01  9.5802E-01  7.4871E+00  1.0000E-02  8.1850E-01  5.8516E-01
             2.4447E+00
 PARAMETER:  7.4856E-02 -3.5106E+00 -6.7539E-01  4.4601E-01 -8.8113E-01  5.7110E-02  2.1132E+00 -8.2544E+00 -1.0028E-01 -4.3587E-01
             9.9393E-01
 GRADIENT:  -3.9025E+00 -2.9009E+00  8.6061E+00  1.6140E+01 -1.6659E+01  1.9678E+00 -8.0539E+00  0.0000E+00 -1.9739E+00  2.9626E+00
             2.9179E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1406.51784175639        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     1226            RESET HESSIAN, TYPE II
 NPARAMETR:  9.7540E-01  2.7046E-02  4.6118E-01  1.4121E+00  3.7572E-01  9.5786E-01  7.5516E+00  1.0000E-02  8.1870E-01  5.8420E-01
             2.4401E+00
 PARAMETER:  7.5097E-02 -3.5102E+00 -6.7397E-01  4.4509E-01 -8.7890E-01  5.6947E-02  2.1218E+00 -8.2658E+00 -1.0004E-01 -4.3751E-01
             9.9203E-01
 GRADIENT:   3.2382E+00 -1.3174E+00  6.6189E+00  2.0524E+01 -6.4004E+00  2.1993E+00 -2.9940E+00  0.0000E+00 -1.7152E+00  1.9178E+00
             1.2491E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1406.55697893279        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1297
 NPARAMETR:  9.7595E-01  2.6466E-02  4.6084E-01  1.4075E+00  3.7644E-01  9.5197E-01  7.5844E+00  1.0000E-02  8.2614E-01  5.7778E-01
             2.4487E+00
 PARAMETER:  7.5660E-02 -3.5319E+00 -6.7471E-01  4.4183E-01 -8.7700E-01  5.0783E-02  2.1261E+00 -8.2658E+00 -9.0987E-02 -4.4856E-01
             9.9555E-01
 GRADIENT:   5.8675E+00  1.0357E+00  3.2610E-01  5.2755E+00  6.0543E+00 -2.4788E-01  2.3737E+00  0.0000E+00  8.4932E-01 -3.6570E-01
             2.7372E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1406.55945705010        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1367
 NPARAMETR:  9.7457E-01  2.5574E-02  4.5956E-01  1.4047E+00  3.7562E-01  9.5142E-01  7.6693E+00  1.0000E-02  8.2487E-01  5.7249E-01
             2.4485E+00
 PARAMETER:  7.4246E-02 -3.5662E+00 -6.7749E-01  4.3980E-01 -8.7917E-01  5.0202E-02  2.1372E+00 -8.2658E+00 -9.2535E-02 -4.5777E-01
             9.9547E-01
 GRADIENT:   2.8755E+00  9.4687E-01  1.1993E+00  1.1656E+00  6.0278E+00 -3.9458E-01  2.2447E+00  0.0000E+00  3.6330E-01 -8.3656E-01
            -2.2893E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1406.58750067666        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1555             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7610E-01  2.4760E-02  4.5692E-01  1.4049E+00  3.7372E-01  9.5348E-01  7.7944E+00  1.0000E-02  8.2147E-01  5.7356E-01
             2.4459E+00
 PARAMETER:  7.5811E-02 -3.5985E+00 -6.8325E-01  4.4000E-01 -8.8425E-01  5.2359E-02  2.1534E+00 -8.2658E+00 -9.6657E-02 -4.5589E-01
             9.9443E-01
 GRADIENT:   6.1919E+00  6.0189E-01  2.2648E+00  5.0342E+00  2.9253E+00  3.4491E-01  1.5343E+00  0.0000E+00 -9.0771E-01 -2.6779E-01
            -3.9280E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1406.58780105682        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1716            RESET HESSIAN, TYPE II
 NPARAMETR:  9.7595E-01  2.4759E-02  4.5691E-01  1.4049E+00  3.7371E-01  9.5407E-01  7.7942E+00  1.0000E-02  8.2176E-01  5.7357E-01
             2.4459E+00
 PARAMETER:  7.5658E-02 -3.5986E+00 -6.8327E-01  4.3999E-01 -8.8426E-01  5.2980E-02  2.1534E+00 -8.2658E+00 -9.6306E-02 -4.5587E-01
             9.9441E-01
 GRADIENT:   5.7533E+00  2.4883E-01  2.7062E+00  5.7832E+00  2.2801E+00  6.3955E-01  7.7407E-01  0.0000E+00 -7.9764E-01 -2.1687E-02
             1.2843E-02

0ITERATION NO.:   79    OBJECTIVE VALUE:  -1406.58816968629        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:     1831
 NPARAMETR:  9.7594E-01  2.4719E-02  4.5720E-01  1.4043E+00  3.7357E-01  9.5410E-01  7.7868E+00  1.0000E-02  8.2194E-01  5.7332E-01
             2.4435E+00
 PARAMETER:  7.5542E-02 -3.5986E+00 -6.8332E-01  4.3999E-01 -8.8426E-01  5.2914E-02  2.1534E+00 -8.2658E+00 -9.5991E-02 -4.5586E-01
             9.9441E-01
 GRADIENT:  -1.8223E+03  1.6288E+02 -5.3558E+02  8.1690E+02  6.6975E+02 -3.6447E+03  2.7081E+02  0.0000E+00  3.6430E+03  7.9697E+02
             3.6007E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1831
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.1867E-04  1.7298E-02  1.7061E-05 -1.1859E-02 -2.4048E-03
 SE:             2.9070E-02  8.8704E-03  2.3330E-04  2.6814E-02  1.8958E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8851E-01  5.1166E-02  9.4171E-01  6.5830E-01  8.9906E-01

 ETASHRINKSD(%)  2.6104E+00  7.0283E+01  9.9218E+01  1.0171E+01  3.6489E+01
 ETASHRINKVR(%)  5.1527E+00  9.1169E+01  9.9994E+01  1.9308E+01  5.9663E+01
 EBVSHRINKSD(%)  2.3608E+00  7.9860E+01  9.9160E+01  9.7020E+00  3.5004E+01
 EBVSHRINKVR(%)  4.6659E+00  9.5944E+01  9.9993E+01  1.8463E+01  5.7755E+01
 RELATIVEINF(%)  9.4555E+01  2.9715E+00  2.6458E-04  2.8939E+01  1.6055E+00
 EPSSHRINKSD(%)  3.1905E+01
 EPSSHRINKVR(%)  5.3631E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1406.5881696862866     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -671.43734312254844     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.13
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.50
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1406.588       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  2.48E-02  4.57E-01  1.40E+00  3.74E-01  9.54E-01  7.79E+00  1.00E-02  8.22E-01  5.74E-01  2.45E+00
 


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
+        9.56E+06
 
 TH 2
+       -6.18E+03  4.11E+07
 
 TH 3
+       -2.65E+03  5.74E+04  3.40E+06
 
 TH 4
+        3.50E+02  2.97E+04  1.23E+06  8.59E+05
 
 TH 5
+       -2.04E+03 -1.23E+05 -2.30E+06 -1.15E+06  3.00E+06
 
 TH 6
+       -6.95E+03 -1.30E+03 -8.80E+02  7.48E+02 -4.93E+02  1.00E+07
 
 TH 7
+       -3.13E+01  1.57E+05  2.97E+02  1.51E+02 -6.34E+02 -7.23E+00  1.16E+03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.16E+04  1.29E+03  2.05E+03 -3.72E+02  1.10E+03  6.95E+03  6.20E+00  0.00E+00  1.35E+07
 
 TH10
+        2.44E+03  2.43E+04  9.24E+03 -1.14E+03 -2.75E+06  1.82E+03  1.26E+02  0.00E+00 -2.59E+03  1.33E+06
 
 TH11
+        3.50E+01  1.25E+04  3.12E+05  1.58E+05 -2.92E+05  2.03E+02  6.48E+01  0.00E+00 -4.27E+01 -4.50E+01  5.53E+04
 
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
 #CPUT: Total CPU Time in Seconds,       31.694
Stop Time:
Sat Sep 18 09:36:01 CDT 2021
