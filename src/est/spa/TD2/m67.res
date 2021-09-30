Wed Sep 29 19:11:55 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat67.csv ignore=@
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
Current Date:       29 SEP 2021
Days until program expires : 200
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1652.03888666310        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1364E+02 -7.3924E+01  2.7459E+00 -1.0702E+02 -1.9598E+00  4.1315E+01 -1.4715E+01  5.6182E+00 -1.2678E+01  4.7553E+00
             6.4213E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1663.47433515611        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.9553E-01  1.0737E+00  1.0415E+00  1.0729E+00  1.0468E+00  1.0526E+00  1.0876E+00  9.7083E-01  1.0201E+00  9.9288E-01
             9.7961E-01
 PARAMETER:  9.5518E-02  1.7112E-01  1.4067E-01  1.7037E-01  1.4570E-01  1.5130E-01  1.8396E-01  7.0399E-02  1.1990E-01  9.2857E-02
             7.9396E-02
 GRADIENT:   4.2186E-01 -2.7630E+00  1.9363E+00 -2.4063E+00 -5.1738E-01  9.7997E+00 -7.4696E+00  7.3676E-01  3.1154E+00 -5.4958E+00
            -5.3279E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1664.40006314616        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.9968E-01  9.9916E-01  1.1058E+00  1.1263E+00  1.0633E+00  1.0289E+00  1.2836E+00  8.4765E-01  9.2473E-01  1.0779E+00
             9.9755E-01
 PARAMETER:  9.9676E-02  9.9160E-02  2.0056E-01  2.1891E-01  1.6137E-01  1.2845E-01  3.4966E-01 -6.5285E-02  2.1744E-02  1.7499E-01
             9.7548E-02
 GRADIENT:   1.1006E+01  2.8889E+00  2.6147E+00  2.0298E+00  5.0864E+00  1.3054E+00 -1.6854E+00 -1.5887E+00 -3.3861E+00 -5.2704E-01
            -9.6784E-02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1664.73785398090        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  9.9238E-01  8.6161E-01  1.1138E+00  1.2112E+00  9.9465E-01  1.0209E+00  1.4107E+00  8.6434E-01  9.0491E-01  1.0320E+00
             9.8945E-01
 PARAMETER:  9.2354E-02 -4.8949E-02  2.0779E-01  2.9164E-01  9.4636E-02  1.2067E-01  4.4407E-01 -4.5791E-02  8.1608E-05  1.3148E-01
             8.9393E-02
 GRADIENT:  -2.2698E+00  3.0378E+00  1.0587E+00  2.0982E+00 -2.8758E+00 -1.3075E+00 -2.8731E-02  2.1766E-02  6.3326E-01  4.6756E-01
            -1.0215E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1664.78158404015        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      711
 NPARAMETR:  9.9132E-01  7.2733E-01  1.2097E+00  1.3023E+00  9.8641E-01  1.0211E+00  1.5588E+00  9.3648E-01  8.7136E-01  1.0427E+00
             9.9135E-01
 PARAMETER:  9.1285E-02 -2.1838E-01  2.9034E-01  3.6416E-01  8.6313E-02  1.2086E-01  5.4390E-01  3.4373E-02 -3.7699E-02  1.4183E-01
             9.1310E-02
 GRADIENT:  -6.4071E-01  5.1525E+00  2.3355E+00  6.7320E+00 -4.1567E+00 -4.6213E-01  2.9125E-01 -1.5575E-01 -3.5325E-01  1.5936E-02
            -5.5063E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1664.79278814213        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  9.9053E-01  6.1644E-01  1.3022E+00  1.3773E+00  9.8849E-01  1.0223E+00  1.6896E+00  1.0112E+00  8.5039E-01  1.0590E+00
             9.9406E-01
 PARAMETER:  9.0488E-02 -3.8380E-01  3.6407E-01  4.2012E-01  8.8425E-02  1.2203E-01  6.2451E-01  1.1117E-01 -6.2063E-02  1.5730E-01
             9.4045E-02
 GRADIENT:   1.3591E+00  5.8375E+00  3.0756E+00  9.2356E+00 -4.3263E+00  6.8603E-01  4.7652E-01 -4.0786E-01 -1.1889E+00 -3.5514E-01
             9.3355E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1664.80711769076        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1061
 NPARAMETR:  9.9007E-01  5.1334E-01  1.4070E+00  1.4491E+00  9.9824E-01  1.0248E+00  1.8225E+00  1.1067E+00  8.3318E-01  1.0779E+00
             9.9772E-01
 PARAMETER:  9.0025E-02 -5.6682E-01  4.4144E-01  4.7095E-01  9.8236E-02  1.2445E-01  7.0022E-01  2.0141E-01 -8.2501E-02  1.7498E-01
             9.7716E-02
 GRADIENT:   4.0181E+00  6.9538E+00  3.8341E+00  1.3543E+01 -4.6153E+00  2.2856E+00  7.4606E-01 -5.4880E-01 -2.2786E+00 -6.4802E-01
             1.1173E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1664.81461395658        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1237
 NPARAMETR:  9.9020E-01  4.3680E-01  1.4909E+00  1.5025E+00  1.0077E+00  1.0283E+00  1.9514E+00  1.1839E+00  8.1808E-01  1.0903E+00
             1.0016E+00
 PARAMETER:  9.0151E-02 -7.2828E-01  4.9937E-01  5.0714E-01  1.0767E-01  1.2789E-01  7.6857E-01  2.6886E-01 -1.0079E-01  1.8641E-01
             1.0162E-01
 GRADIENT:   6.9127E+00  7.2497E+00  4.3070E+00  1.7351E+01 -4.1478E+00  4.0886E+00  7.0052E-01 -7.3808E-01 -3.4397E+00 -1.0050E+00
             2.2728E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1665.14570074725        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1418
 NPARAMETR:  9.8620E-01  3.9457E-01  1.5013E+00  1.5162E+00  1.0079E+00  1.0170E+00  1.9403E+00  1.2079E+00  8.2549E-01  1.1026E+00
             9.9372E-01
 PARAMETER:  8.6101E-02 -8.2997E-01  5.0635E-01  5.1617E-01  1.0782E-01  1.1681E-01  7.6283E-01  2.8892E-01 -9.1780E-02  1.9770E-01
             9.3697E-02
 GRADIENT:   8.1199E-02  1.6610E+00 -9.4000E-01 -5.3999E+00  3.3537E+00  4.4463E-03 -1.4292E-02  1.6636E-01 -1.6923E-01 -9.5764E-02
             9.6549E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1665.19597455312        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1594
 NPARAMETR:  9.8387E-01  3.2241E-01  1.5025E+00  1.5691E+00  9.8237E-01  1.0150E+00  1.9600E+00  1.2026E+00  8.2158E-01  1.0951E+00
             9.9042E-01
 PARAMETER:  8.3743E-02 -1.0319E+00  5.0713E-01  5.5047E-01  8.2214E-02  1.1493E-01  7.7293E-01  2.8446E-01 -9.6523E-02  1.9086E-01
             9.0375E-02
 GRADIENT:  -2.6359E+00  3.6446E+00  1.6316E-01  1.0533E+01 -9.4230E-01 -3.6555E-01 -6.5198E-01 -6.6198E-01  7.6155E-01  9.5817E-02
            -1.1941E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1665.20084546902        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1770
 NPARAMETR:  9.8307E-01  2.8529E-01  1.5064E+00  1.5946E+00  9.7156E-01  1.0143E+00  2.0086E+00  1.2092E+00  8.1629E-01  1.0911E+00
             9.8977E-01
 PARAMETER:  8.2928E-02 -1.1543E+00  5.0975E-01  5.6659E-01  7.1151E-02  1.1422E-01  7.9746E-01  2.8994E-01 -1.0299E-01  1.8715E-01
             8.9722E-02
 GRADIENT:  -3.1238E+00  3.9570E+00  2.6918E-01  1.5553E+01 -2.5622E+00 -4.3777E-01 -7.0452E-01 -7.3426E-01  7.4770E-01  1.8802E-01
            -1.3228E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1665.20264878962        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1948
 NPARAMETR:  9.8263E-01  2.6021E-01  1.5109E+00  1.6115E+00  9.6519E-01  1.0139E+00  2.0481E+00  1.2158E+00  8.1240E-01  1.0885E+00
             9.8947E-01
 PARAMETER:  8.2478E-02 -1.2462E+00  5.1270E-01  5.7717E-01  6.4565E-02  1.1380E-01  8.1692E-01  2.9536E-01 -1.0777E-01  1.8484E-01
             8.9414E-02
 GRADIENT:  -3.2127E+00  4.0066E+00  3.4894E-01  1.8341E+01 -3.5365E+00 -4.5522E-01 -7.0502E-01 -7.6652E-01  7.0712E-01  2.1816E-01
            -1.3653E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1665.20408235372        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2126
 NPARAMETR:  9.8221E-01  2.3672E-01  1.5163E+00  1.6273E+00  9.5985E-01  1.0135E+00  2.0898E+00  1.2234E+00  8.0881E-01  1.0864E+00
             9.8930E-01
 PARAMETER:  8.2052E-02 -1.3409E+00  5.1629E-01  5.8690E-01  5.9018E-02  1.1339E-01  8.3709E-01  3.0162E-01 -1.1219E-01  1.8287E-01
             8.9244E-02
 GRADIENT:  -3.2926E+00  3.9544E+00  4.2229E-01  2.0613E+01 -4.3746E+00 -4.7768E-01 -6.8105E-01 -7.7955E-01  7.4958E-01  2.3190E-01
            -1.3571E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1665.20688884690        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2302
 NPARAMETR:  9.8175E-01  2.0939E-01  1.5242E+00  1.6455E+00  9.5435E-01  1.0130E+00  2.1447E+00  1.2339E+00  8.0443E-01  1.0842E+00
             9.8922E-01
 PARAMETER:  8.1580E-02 -1.4635E+00  5.2144E-01  5.9803E-01  5.3271E-02  1.1293E-01  8.6301E-01  3.1017E-01 -1.1762E-01  1.8084E-01
             8.9163E-02
 GRADIENT:  -3.3141E+00  3.7951E+00  5.1497E-01  2.2909E+01 -5.2960E+00 -4.9912E-01 -6.3526E-01 -7.7417E-01  7.8528E-01  2.4407E-01
            -1.3053E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1665.21036518870        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2478
 NPARAMETR:  9.8128E-01  1.7551E-01  1.5361E+00  1.6681E+00  9.4849E-01  1.0126E+00  2.2252E+00  1.2492E+00  7.9837E-01  1.0820E+00
             9.8911E-01
 PARAMETER:  8.1101E-02 -1.6400E+00  5.2922E-01  6.1167E-01  4.7113E-02  1.1257E-01  8.9986E-01  3.2252E-01 -1.2518E-01  1.7880E-01
             8.9047E-02
 GRADIENT:  -3.0978E+00  3.4932E+00  6.5064E-01  2.5580E+01 -6.4216E+00 -4.3550E-01 -5.6098E-01 -7.4007E-01  6.6363E-01  2.7276E-01
            -1.2631E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1665.22685202053        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2653
 NPARAMETR:  9.8077E-01  1.3839E-01  1.5522E+00  1.6925E+00  9.4348E-01  1.0121E+00  2.3377E+00  1.2690E+00  7.9150E-01  1.0800E+00
             9.8928E-01
 PARAMETER:  8.0587E-02 -1.8777E+00  5.3969E-01  6.2623E-01  4.1824E-02  1.1200E-01  9.4917E-01  3.3824E-01 -1.3382E-01  1.7701E-01
             8.9220E-02
 GRADIENT:  -2.8013E+00  2.9836E+00  7.8878E-01  2.7546E+01 -7.4468E+00 -4.2864E-01 -4.4935E-01 -6.6942E-01  6.1323E-01  2.8157E-01
            -1.1093E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1665.30047924185        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2828
 NPARAMETR:  9.8015E-01  8.5012E-02  1.5857E+00  1.7273E+00  9.4054E-01  1.0117E+00  2.5920E+00  1.3069E+00  7.8034E-01  1.0793E+00
             9.8987E-01
 PARAMETER:  7.9951E-02 -2.3650E+00  5.6103E-01  6.4653E-01  3.8694E-02  1.1162E-01  1.0524E+00  3.6762E-01 -1.4803E-01  1.7627E-01
             8.9822E-02
 GRADIENT:  -2.1194E+00  1.9857E+00  9.3826E-01  2.8073E+01 -8.2873E+00 -2.4606E-01 -2.5910E-01 -5.0679E-01  4.4197E-01  2.6703E-01
            -7.9923E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1665.67038325742        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     3007
 NPARAMETR:  9.7990E-01  2.2054E-02  1.6816E+00  1.7644E+00  9.6001E-01  1.0090E+00  3.5277E+00  1.4025E+00  7.6285E-01  1.0918E+00
             9.9142E-01
 PARAMETER:  7.9698E-02 -3.7142E+00  6.1977E-01  6.6779E-01  5.9191E-02  1.0898E-01  1.3606E+00  4.3825E-01 -1.7070E-01  1.8781E-01
             9.1383E-02
 GRADIENT:  -2.1771E-01  4.4360E-01  2.8549E-01  1.3191E+01 -4.5375E+00 -8.9575E-01 -4.1893E-02 -1.4560E-01 -1.0839E-01  1.1945E-01
            -2.9320E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1665.85908618050        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3189
 NPARAMETR:  9.8082E-01  1.0000E-02  1.7454E+00  1.7442E+00  9.8384E-01  1.0119E+00  4.4005E+00  1.4628E+00  7.5891E-01  1.1043E+00
             9.9266E-01
 PARAMETER:  8.0632E-02 -4.5395E+00  6.5699E-01  6.5627E-01  8.3704E-02  1.1184E-01  1.5817E+00  4.8037E-01 -1.7587E-01  1.9924E-01
             9.2634E-02
 GRADIENT:   2.6336E+00  6.2857E-03 -1.9068E-01 -5.2016E+01  4.0329E+00  3.9962E-01 -1.1152E-02  2.1735E-01  3.3360E-01 -4.7948E-01
             5.4262E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1665.94409818072        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3376
 NPARAMETR:  9.8039E-01  1.0000E-02  1.7404E+00  1.7551E+00  9.8052E-01  1.0115E+00  4.7179E+00  1.4582E+00  7.5822E-01  1.1041E+00
             9.9176E-01
 PARAMETER:  8.0197E-02 -4.5633E+00  6.5412E-01  6.6253E-01  8.0328E-02  1.1140E-01  1.6514E+00  4.7720E-01 -1.7678E-01  1.9907E-01
             9.1731E-02
 GRADIENT:   1.5194E+00  0.0000E+00 -5.7795E-02 -2.8611E+01  1.3481E+00  1.8862E-01 -1.3401E-02  1.5107E-01  1.2002E-01 -1.4225E-01
             6.4031E-03

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1665.94814463541        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     3571             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8038E-01  1.0000E-02  1.7365E+00  1.7548E+00  9.7745E-01  1.0115E+00  5.1308E+00  1.4532E+00  7.5830E-01  1.1049E+00
             9.9173E-01
 PARAMETER:  8.0188E-02 -4.5633E+00  6.5185E-01  6.6238E-01  7.7196E-02  1.1139E-01  1.7353E+00  4.7380E-01 -1.7667E-01  1.9979E-01
             9.1697E-02
 GRADIENT:   4.0411E+02  0.0000E+00  1.0797E+01  1.2356E+03  4.9220E+00  5.6217E+01  1.1657E-01  1.7616E+00  2.8337E+01  2.0828E+00
             8.4049E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1665.95164607868        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     3767             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8038E-01  1.0000E-02  1.7295E+00  1.7545E+00  9.7727E-01  1.0114E+00  5.6590E+00  1.4485E+00  7.5832E-01  1.1019E+00
             9.9171E-01
 PARAMETER:  8.0181E-02 -4.5633E+00  6.4783E-01  6.6220E-01  7.7005E-02  1.1138E-01  1.8333E+00  4.7052E-01 -1.7665E-01  1.9701E-01
             9.1677E-02
 GRADIENT:   4.0398E+02  0.0000E+00  9.5972E+00  1.2347E+03  7.0961E+00  5.6190E+01  1.4636E-01  1.7232E+00  2.8316E+01  1.5854E+00
             8.2753E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1665.95489662403        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     3957
 NPARAMETR:  9.8037E-01  1.0000E-02  1.7253E+00  1.7543E+00  9.7604E-01  1.0114E+00  6.2068E+00  1.4439E+00  7.5834E-01  1.1011E+00
             9.9158E-01
 PARAMETER:  8.0176E-02 -4.5633E+00  6.4539E-01  6.6206E-01  7.5744E-02  1.1137E-01  1.9257E+00  4.6736E-01 -1.7662E-01  1.9632E-01
             9.1541E-02
 GRADIENT:   1.5233E+00  0.0000E+00 -4.7470E-01 -2.8740E+01  1.7388E+00  2.0775E-01 -1.7805E-02  9.6288E-02  1.2335E-01 -2.3503E-01
            -2.8581E-02

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1665.95998000873        NO. OF FUNC. EVALS.: 202
 CUMULATIVE NO. OF FUNC. EVALS.:     4159             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8036E-01  1.0000E-02  1.7225E+00  1.7541E+00  9.7314E-01  1.0114E+00  6.9792E+00  1.4398E+00  7.5835E-01  1.1020E+00
             9.9161E-01
 PARAMETER:  8.0166E-02 -4.5633E+00  6.4380E-01  6.6194E-01  7.2772E-02  1.1137E-01  2.0429E+00  4.6452E-01 -1.7661E-01  1.9708E-01
             9.1574E-02
 GRADIENT:   4.0394E+02  0.0000E+00  1.0554E+01  1.2343E+03  5.1454E+00  5.6156E+01  2.3261E-01  1.6511E+00  2.8354E+01  1.9704E+00
             8.2220E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1665.96243325688        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     4348
 NPARAMETR:  9.8036E-01  1.0000E-02  1.7204E+00  1.7539E+00  9.7203E-01  1.0114E+00  7.5796E+00  1.4374E+00  7.5834E-01  1.1014E+00
             9.9160E-01
 PARAMETER:  8.0162E-02 -4.5633E+00  6.4255E-01  6.6186E-01  7.1629E-02  1.1137E-01  2.1255E+00  4.6284E-01 -1.7663E-01  1.9660E-01
             9.1566E-02
 GRADIENT:   1.5801E+00  0.0000E+00  7.2831E-01 -2.8588E+01 -7.5516E-01  2.0149E-01 -1.9390E-02  8.0085E-02  1.7511E-01  2.2322E-01
             2.4267E-02

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1665.96683908074        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     4542             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8035E-01  1.0000E-02  1.7141E+00  1.7537E+00  9.7223E-01  1.0114E+00  8.8219E+00  1.4332E+00  7.5822E-01  1.0987E+00
             9.9149E-01
 PARAMETER:  8.0158E-02 -4.5633E+00  6.3890E-01  6.6171E-01  7.1833E-02  1.1135E-01  2.2772E+00  4.5989E-01 -1.7678E-01  1.9410E-01
             9.1456E-02
 GRADIENT:   4.0388E+02  0.0000E+00  9.4587E+00  1.2336E+03  7.1247E+00  5.6122E+01  3.8166E-01  1.5861E+00  2.8379E+01  1.5028E+00
             7.6924E-01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1665.97400281920        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     4684
 NPARAMETR:  9.8031E-01  1.0000E-02  1.7148E+00  1.7535E+00  9.7134E-01  1.0114E+00  1.1934E+01  1.4329E+00  7.5820E-01  1.0997E+00
             9.9152E-01
 PARAMETER:  8.0114E-02 -4.5633E+00  6.3929E-01  6.6164E-01  7.0922E-02  1.1132E-01  2.5794E+00  4.5970E-01 -1.7681E-01  1.9502E-01
             9.1479E-02
 GRADIENT:   1.4456E+00  0.0000E+00  6.8858E-02 -2.9024E+01  4.1352E-01  2.0387E-01  1.4791E-02  8.2758E-02  3.9782E-01  9.6198E-03
             4.0361E-03

0ITERATION NO.:  133    OBJECTIVE VALUE:  -1665.97429409133        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     4785
 NPARAMETR:  9.8031E-01  1.0000E-02  1.7096E+00  1.7533E+00  9.7080E-01  1.0114E+00  1.1959E+01  1.4308E+00  7.5838E-01  1.0985E+00
             9.9178E-01
 PARAMETER:  8.0192E-02 -4.5633E+00  6.3867E-01  6.6177E-01  6.9922E-02  1.1139E-01  2.5971E+00  4.5779E-01 -1.7835E-01  1.9494E-01
             9.1455E-02
 GRADIENT:   8.5633E-02  0.0000E+00  4.6504E-01  5.1901E-01 -2.0305E-01  1.4587E-02  4.1595E-03 -1.8968E-02 -3.1026E-01  6.6452E-02
            -5.7352E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4785
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9523E-04  4.4792E-04 -3.5269E-02 -7.3871E-03 -4.7788E-02
 SE:             2.9843E-02  1.8321E-03  1.8727E-02  2.9250E-02  1.9970E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9478E-01  8.0686E-01  5.9658E-02  8.0062E-01  1.6709E-02

 ETASHRINKSD(%)  2.2462E-02  9.3862E+01  3.7262E+01  2.0089E+00  3.3099E+01
 ETASHRINKVR(%)  4.4920E-02  9.9623E+01  6.0640E+01  3.9775E+00  5.5243E+01
 EBVSHRINKSD(%)  3.9983E-01  9.4078E+01  4.0706E+01  2.3721E+00  2.9584E+01
 EBVSHRINKVR(%)  7.9806E-01  9.9649E+01  6.4843E+01  4.6880E+00  5.0416E+01
 RELATIVEINF(%)  9.3617E+01  9.2626E-03  9.4004E+00  2.9694E+00  7.0881E+00
 EPSSHRINKSD(%)  4.5212E+01
 EPSSHRINKVR(%)  6.9983E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1665.9742940913336     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -930.82346752759543     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    66.69
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0INVERSE COVARIANCE MATRIX SET TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S
 Elapsed covariance  time in seconds:     8.41
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1665.974       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.00E-02  1.71E+00  1.75E+00  9.70E-01  1.01E+00  1.21E+01  1.43E+00  7.57E-01  1.10E+00  9.91E-01
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         2.98E-02  0.00E+00  2.92E-01  4.67E-02  9.86E-02  6.57E-02  3.17E-03  3.23E-01  6.15E-02  1.95E-01  7.08E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        8.88E-04
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        1.76E-04  0.00E+00  8.51E-02
 
 TH 4
+        4.41E-06  0.00E+00  6.88E-03  2.18E-03
 
 TH 5
+        8.92E-05  0.00E+00  2.54E-02  2.34E-03  9.72E-03
 
 TH 6
+        1.06E-04  0.00E+00  2.74E-03  1.50E-05  9.71E-04  4.31E-03
 
 TH 7
+        2.31E-06  0.00E+00  2.48E-04  1.38E-05  8.66E-05 -1.87E-05  1.01E-05
 
 TH 8
+        1.00E-04  0.00E+00  6.70E-02  4.43E-03  1.95E-02  3.58E-03  2.96E-04  1.04E-01
 
 TH 9
+       -4.33E-07  0.00E+00 -1.53E-03 -8.10E-04 -7.68E-04 -2.38E-04  9.83E-05 -2.29E-04  3.78E-03
 
 TH10
+       -2.32E-04  0.00E+00  8.71E-03  8.85E-04  5.44E-03 -1.20E-03 -1.20E-05 -1.07E-02 -1.45E-03  3.81E-02
 
 TH11
+       -2.35E-04  0.00E+00  2.42E-03  4.34E-04  6.10E-04 -2.61E-04  5.45E-06 -1.86E-03 -3.52E-04  1.56E-03  5.02E-03
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.98E-02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        2.03E-02  0.00E+00  2.92E-01
 
 TH 4
+        3.17E-03  0.00E+00  5.05E-01  4.67E-02
 
 TH 5
+        3.04E-02  0.00E+00  8.82E-01  5.10E-01  9.86E-02
 
 TH 6
+        5.40E-02  0.00E+00  1.43E-01  4.90E-03  1.50E-01  6.57E-02
 
 TH 7
+        2.44E-02  0.00E+00  2.68E-01  9.35E-02  2.77E-01 -8.97E-02  3.17E-03
 
 TH 8
+        1.04E-02  0.00E+00  7.12E-01  2.95E-01  6.12E-01  1.69E-01  2.90E-01  3.23E-01
 
 TH 9
+       -2.36E-04  0.00E+00 -8.53E-02 -2.82E-01 -1.27E-01 -5.90E-02  5.04E-01 -1.15E-02  6.15E-02
 
 TH10
+       -3.99E-02  0.00E+00  1.53E-01  9.71E-02  2.83E-01 -9.38E-02 -1.93E-02 -1.71E-01 -1.21E-01  1.95E-01
 
 TH11
+       -1.11E-01  0.00E+00  1.17E-01  1.31E-01  8.74E-02 -5.60E-02  2.42E-02 -8.16E-02 -8.08E-02  1.13E-01  7.08E-02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.15E+03
 
 TH 2
+        1.03E-13  1.09E-28
 
 TH 3
+       -9.92E-01  2.33E-14  4.57E+01
 
 TH 4
+        3.41E+00  1.65E-15 -1.39E-01  6.42E+02
 
 TH 5
+       -2.62E+01 -1.14E-13 -1.47E+02 -8.20E+01  5.06E+02
 
 TH 6
+       -2.01E+01  4.28E-14  1.51E+01  1.80E+01 -3.93E+01  2.29E+02
 
 TH 7
+        1.99E-03  4.67E-18  5.50E-04  8.82E-03 -1.68E-04  2.07E-03  4.90E-06
 
 TH 8
+        5.53E+00  9.61E-15 -1.02E+00 -6.29E+00 -1.38E+00 -2.86E+00 -1.47E-04  1.40E+00
 
 TH 9
+        3.29E+00  3.01E-14  1.41E+00  1.06E+02  3.10E+00  8.92E+00  3.71E-02 -1.65E+00  2.83E+02
 
 TH10
+        1.05E+01  2.16E-14  1.79E+01 -1.92E+00 -6.55E+01 -3.60E+00 -7.35E-04  1.53E+00 -6.80E+00  1.01E+01
 
 TH11
+        5.76E+01  1.21E-13 -2.74E+01 -1.48E+01  2.83E+01  1.60E+01  3.01E-03  1.63E+01  1.90E+01  9.71E+00  2.15E+02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.12E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -3.71E+00  0.00E+00  6.79E+01
 
 TH 4
+       -2.18E+01  0.00E+00 -3.02E+00  6.02E+02
 
 TH 5
+        4.54E+01  0.00E+00 -1.43E+02 -4.56E+00  5.14E+02
 
 TH 6
+        2.14E+01  0.00E+00 -2.89E+00 -2.07E+01  2.91E+01  1.58E+02
 
 TH 7
+        4.04E-03  0.00E+00  2.16E-03 -2.29E-02  1.60E-03 -2.78E-03  1.96E-05
 
 TH 8
+       -7.88E+00  0.00E+00 -1.85E+01 -2.16E+01 -7.73E+00  3.13E-02  2.40E-03  2.00E+01
 
 TH 9
+        5.35E+00  0.00E+00  2.01E+01 -1.28E+02 -1.91E+01 -1.20E+01  7.64E-02  3.37E+00  3.80E+02
 
 TH10
+       -1.75E+01  0.00E+00  2.29E+00 -1.72E+01 -1.02E+02 -1.62E+01 -4.43E-03  1.41E+01 -1.36E+01  8.84E+01
 
 TH11
+       -7.54E+01  0.00E+00  1.04E+00  2.95E+00 -4.77E+01 -8.98E+00  1.68E-03  9.61E+00  4.23E+00  2.31E+01  2.07E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.04
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       73.755
Stop Time:
Wed Sep 29 19:13:16 CDT 2021
