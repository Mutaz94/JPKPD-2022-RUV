Wed Sep 29 18:54:14 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1618.37031595902        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6771E+02 -2.8144E+01 -4.8985E+01  2.4239E+01  6.7718E+01  2.3213E+01 -2.0543E+01  7.3189E+00 -2.4112E+01  1.6884E+01
            -3.4125E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1629.41239278992        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  9.8641E-01  1.0428E+00  1.0977E+00  1.0126E+00  1.0144E+00  1.0930E+00  1.0965E+00  9.7274E-01  1.1329E+00  9.1435E-01
             1.1104E+00
 PARAMETER:  8.6313E-02  1.4188E-01  1.9318E-01  1.1252E-01  1.1426E-01  1.8895E-01  1.9216E-01  7.2364E-02  2.2481E-01  1.0457E-02
             2.0474E-01
 GRADIENT:  -1.2358E+00 -1.2934E+01 -1.3880E+01 -1.7587E+00  1.5381E+01  8.6851E+00 -6.8827E+00  2.2662E+00  1.1142E-01  4.6066E+00
             1.0472E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1630.19050662998        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      354
 NPARAMETR:  9.8390E-01  9.3493E-01  1.0749E+00  1.0809E+00  9.4482E-01  1.0600E+00  1.3418E+00  9.1721E-01  1.0777E+00  8.1723E-01
             1.0755E+00
 PARAMETER:  8.3773E-02  3.2720E-02  1.7220E-01  1.7782E-01  4.3239E-02  1.5823E-01  3.9399E-01  1.3577E-02  1.7485E-01 -1.0184E-01
             1.7276E-01
 GRADIENT:  -4.7606E+00 -1.3104E+00 -4.9024E+00 -3.5646E-02  8.9089E-03 -3.4994E+00  2.3605E+00  1.6251E+00  5.5833E+00 -6.2867E-01
            -3.0258E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1630.57920448935        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      534
 NPARAMETR:  9.8555E-01  8.4385E-01  1.2409E+00  1.1547E+00  9.7297E-01  1.0681E+00  1.3625E+00  9.8947E-01  1.0296E+00  8.7829E-01
             1.0866E+00
 PARAMETER:  8.5443E-02 -6.9784E-02  3.1583E-01  2.4384E-01  7.2599E-02  1.6584E-01  4.0929E-01  8.9417E-02  1.2914E-01 -2.9778E-02
             1.8309E-01
 GRADIENT:   9.7175E-01  6.3958E+00  2.1616E+00  7.1517E+00 -4.7792E+00  1.1299E-01 -7.6290E-02 -1.7550E-01 -3.3693E-01  1.1621E-01
            -1.8454E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1630.78900678484        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      713
 NPARAMETR:  9.8268E-01  5.8937E-01  1.4511E+00  1.3254E+00  9.5830E-01  1.0638E+00  1.5912E+00  1.1053E+00  9.5441E-01  8.9211E-01
             1.0883E+00
 PARAMETER:  8.2526E-02 -4.2870E-01  4.7232E-01  3.8168E-01  5.7410E-02  1.6183E-01  5.6450E-01  2.0008E-01  5.3335E-02 -1.4166E-02
             1.8459E-01
 GRADIENT:   1.0698E+00  6.4829E+00  4.9756E+00  1.1065E+01 -6.1345E+00 -4.8435E-01 -2.4036E-01 -1.2407E+00 -1.1651E+00 -1.3195E+00
            -4.8285E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1630.97807151409        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      891
 NPARAMETR:  9.7940E-01  4.1985E-01  1.5517E+00  1.4354E+00  9.4193E-01  1.0632E+00  1.8425E+00  1.1498E+00  9.1124E-01  9.1228E-01
             1.0880E+00
 PARAMETER:  7.9183E-02 -7.6785E-01  5.3936E-01  4.6142E-01  4.0174E-02  1.6124E-01  7.1114E-01  2.3961E-01  7.0515E-03  8.1906E-03
             1.8437E-01
 GRADIENT:  -1.2869E+00  4.8026E+00  4.0334E+00  1.1292E+01 -5.7205E+00 -9.3712E-02 -3.3073E-01 -1.7236E+00 -1.1013E+00 -1.3286E-02
            -5.8555E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1631.00158276497        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1066
 NPARAMETR:  9.7814E-01  3.4006E-01  1.6147E+00  1.4920E+00  9.3682E-01  1.0627E+00  2.0174E+00  1.2137E+00  8.9180E-01  9.1845E-01
             1.0876E+00
 PARAMETER:  7.7898E-02 -9.7862E-01  5.7915E-01  5.0011E-01  3.4734E-02  1.6085E-01  8.0179E-01  2.9371E-01 -1.4510E-02  1.4937E-02
             1.8400E-01
 GRADIENT:  -1.7703E+00  5.2179E+00  3.6327E+00  1.7581E+01 -8.2215E+00  2.2899E-02 -3.4751E-01 -9.6204E-01 -1.1082E+00  9.7204E-01
            -3.7202E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1631.06222938002        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1241
 NPARAMETR:  9.7700E-01  2.4903E-01  1.6859E+00  1.5535E+00  9.3217E-01  1.0623E+00  2.3215E+00  1.2922E+00  8.7061E-01  9.2132E-01
             1.0873E+00
 PARAMETER:  7.6727E-02 -1.2902E+00  6.2229E-01  5.4051E-01  2.9758E-02  1.6042E-01  9.4222E-01  3.5632E-01 -3.8559E-02  1.8049E-02
             1.8374E-01
 GRADIENT:  -1.6961E+00  4.3011E+00  2.1118E+00  1.9784E+01 -8.7161E+00  1.3940E-01 -3.0382E-01  2.0056E-01 -9.6833E-01  1.7234E+00
            -2.6448E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1631.18320737983        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1416
 NPARAMETR:  9.7595E-01  1.5324E-01  1.7685E+00  1.6153E+00  9.3067E-01  1.0616E+00  2.9148E+00  1.3732E+00  8.4994E-01  9.2290E-01
             1.0874E+00
 PARAMETER:  7.5651E-02 -1.7757E+00  6.7013E-01  5.7955E-01  2.8148E-02  1.5980E-01  1.1698E+00  4.1717E-01 -6.2592E-02  1.9760E-02
             1.8377E-01
 GRADIENT:  -1.2635E+00  2.5867E+00  3.6861E-01  1.5983E+01 -6.8423E+00  1.8012E-01 -1.8010E-01  1.0856E+00 -3.5491E-01  1.9421E+00
             2.8089E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1631.42606454342        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1595             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7628E-01  7.8459E-02  1.8277E+00  1.6371E+00  9.3487E-01  1.0615E+00  3.9552E+00  1.4010E+00  8.3723E-01  9.1078E-01
             1.0874E+00
 PARAMETER:  7.5997E-02 -2.4452E+00  7.0305E-01  5.9294E-01  3.2654E-02  1.5964E-01  1.4750E+00  4.3722E-01 -7.7656E-02  6.5478E-03
             1.8377E-01
 GRADIENT:   3.6881E+02  1.0817E+01  5.7730E+00  9.6684E+02  1.4090E+01  7.0775E+01  3.1710E+00  7.9170E-01  1.7152E+01 -6.7105E-01
             1.2853E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1631.49886141361        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1777
 NPARAMETR:  9.7601E-01  7.3242E-02  1.8312E+00  1.6475E+00  9.3015E-01  1.0611E+00  4.3305E+00  1.4029E+00  8.3203E-01  9.1545E-01
             1.0877E+00
 PARAMETER:  7.5713E-02 -2.5140E+00  7.0495E-01  5.9926E-01  2.7596E-02  1.5932E-01  1.5657E+00  4.3855E-01 -8.3888E-02  1.1658E-02
             1.8404E-01
 GRADIENT:   1.2679E+00  2.5351E-01  1.2256E-02 -2.2129E+01  1.0305E+00  2.9682E-01 -2.7769E-03  2.2229E-02  9.4677E-02  1.1107E-01
             2.2709E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1631.51526431014        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1958
 NPARAMETR:  9.7581E-01  5.7002E-02  1.8341E+00  1.6568E+00  9.2720E-01  1.0610E+00  4.8427E+00  1.4038E+00  8.2942E-01  9.1402E-01
             1.0877E+00
 PARAMETER:  7.5516E-02 -2.7647E+00  7.0654E-01  6.0488E-01  2.4419E-02  1.5922E-01  1.6775E+00  4.3919E-01 -8.7032E-02  1.0102E-02
             1.8404E-01
 GRADIENT:   1.3636E+00  1.6625E-01 -2.1816E-01 -2.4363E+01  1.4253E+00  3.1421E-01 -2.3313E-02 -3.5769E-02  3.0888E-01  2.0248E-02
             2.4144E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1631.52247910075        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2139
 NPARAMETR:  9.7552E-01  5.0037E-02  1.8354E+00  1.6642E+00  9.2565E-01  1.0607E+00  5.2202E+00  1.4047E+00  8.2819E-01  9.1365E-01
             1.0877E+00
 PARAMETER:  7.5212E-02 -2.8950E+00  7.0725E-01  6.0935E-01  2.2744E-02  1.5891E-01  1.7525E+00  4.3985E-01 -8.8517E-02  9.6882E-03
             1.8403E-01
 GRADIENT:   9.3551E-01  2.7647E-01 -3.0911E-01 -1.9126E+01  8.0689E-01  1.9802E-01 -3.1072E-03 -3.9486E-02  4.2405E-01  7.6987E-02
            -4.0346E-03

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1631.52571922225        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2318
 NPARAMETR:  9.7512E-01  3.5253E-02  1.8404E+00  1.6762E+00  9.2358E-01  1.0603E+00  6.0628E+00  1.4071E+00  8.2553E-01  9.1336E-01
             1.0877E+00
 PARAMETER:  7.4809E-02 -3.2452E+00  7.0996E-01  6.1652E-01  2.0497E-02  1.5855E-01  1.9022E+00  4.4156E-01 -9.1729E-02  9.3764E-03
             1.8410E-01
 GRADIENT:   5.1318E-01  2.6575E-01 -6.4551E-01 -1.5002E+01  7.7553E-01  9.1133E-02 -2.4199E-02 -9.7544E-02  4.7353E-01  6.5318E-02
            -1.3983E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1631.54724410248        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2500
 NPARAMETR:  9.7525E-01  2.4692E-02  1.8447E+00  1.6786E+00  9.2211E-01  1.0606E+00  7.1379E+00  1.4103E+00  8.2305E-01  9.1236E-01
             1.0877E+00
 PARAMETER:  7.4943E-02 -3.6013E+00  7.1234E-01  6.1795E-01  1.8911E-02  1.5883E-01  2.0654E+00  4.4378E-01 -9.4736E-02  8.2833E-03
             1.8407E-01
 GRADIENT:   1.1380E+00  1.2650E-01 -5.3075E-01 -2.3242E+01  1.2122E+00  2.5620E-01 -2.1623E-02 -9.1412E-02  3.2992E-01  8.6560E-03
             1.8225E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1631.55424235683        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2682
 NPARAMETR:  9.7515E-01  1.8217E-02  1.8474E+00  1.6831E+00  9.2115E-01  1.0605E+00  8.5410E+00  1.4123E+00  8.2176E-01  9.1185E-01
             1.0877E+00
 PARAMETER:  7.4834E-02 -3.9054E+00  7.1379E-01  6.2063E-01  1.7868E-02  1.5876E-01  2.2449E+00  4.4520E-01 -9.6306E-02  7.7225E-03
             1.8406E-01
 GRADIENT:   1.1336E+00  1.1973E-01 -5.4748E-01 -2.2846E+01  1.0780E+00  2.4613E-01 -9.5573E-04 -9.2421E-02  3.9680E-01  8.0867E-03
             9.2794E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1631.56330090814        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2869             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7518E-01  1.0000E-02  1.8539E+00  1.6875E+00  9.1970E-01  1.0607E+00  1.2120E+01  1.4183E+00  8.1925E-01  9.1102E-01
             1.0876E+00
 PARAMETER:  7.4866E-02 -5.8704E+00  7.1730E-01  6.2326E-01  1.6291E-02  1.5893E-01  2.5949E+00  4.4943E-01 -9.9360E-02  6.8080E-03
             1.8402E-01
 GRADIENT:   3.6750E+02  0.0000E+00  8.3351E+00  1.0874E+03  5.3528E+00  7.0155E+01  8.7592E-01  9.8301E-01  1.6566E+01  6.5349E-01
             1.4896E+00

0ITERATION NO.:   84    OBJECTIVE VALUE:  -1631.56341292621        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     3006
 NPARAMETR:  9.7518E-01  1.0000E-02  1.8551E+00  1.6875E+00  9.1936E-01  1.0607E+00  1.1887E+01  1.4210E+00  8.1926E-01  9.1062E-01
             1.0877E+00
 PARAMETER:  7.4868E-02 -5.8704E+00  7.1706E-01  6.2326E-01  1.6469E-02  1.5894E-01  2.5815E+00  4.5010E-01 -9.9444E-02  5.3751E-03
             1.8397E-01
 GRADIENT:   2.5240E-03  0.0000E+00 -1.5443E-01 -1.1640E-02  3.3280E-01 -2.6176E-03  1.8621E-03 -4.2414E-02 -1.4372E-02 -5.3826E-02
            -1.4594E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3006
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.1426E-06 -9.5312E-05 -2.7074E-02 -5.6603E-03 -3.8609E-02
 SE:             2.9809E-02  1.7940E-03  1.7834E-02  2.9268E-02  1.9939E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9997E-01  9.5763E-01  1.2898E-01  8.4665E-01  5.2818E-02

 ETASHRINKSD(%)  1.3726E-01  9.3990E+01  4.0255E+01  1.9483E+00  3.3203E+01
 ETASHRINKVR(%)  2.7433E-01  9.9639E+01  6.4305E+01  3.8587E+00  5.5382E+01
 EBVSHRINKSD(%)  4.4926E-01  9.4161E+01  4.3465E+01  2.2928E+00  3.0821E+01
 EBVSHRINKVR(%)  8.9650E-01  9.9659E+01  6.8038E+01  4.5330E+00  5.2142E+01
 RELATIVEINF(%)  9.5185E+01  7.6342E-03  6.7049E+00  2.5965E+00  5.8319E+00
 EPSSHRINKSD(%)  4.3969E+01
 EPSSHRINKVR(%)  6.8605E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1631.5634129262103     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -896.41258636247210     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    42.16
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1631.563       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  1.00E-02  1.85E+00  1.69E+00  9.20E-01  1.06E+00  1.20E+01  1.42E+00  8.19E-01  9.10E-01  1.09E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.10E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -8.40E+00  0.00E+00  3.98E+01
 
 TH 4
+       -7.26E+00  0.00E+00  4.04E-01  6.01E+02
 
 TH 5
+        5.97E+01  0.00E+00 -1.74E+02 -7.68E+01  8.14E+02
 
 TH 6
+        5.83E+01  0.00E+00 -1.32E+01 -1.13E+01  6.76E+01  1.40E+02
 
 TH 7
+       -8.75E-03  0.00E+00  1.91E-03  1.59E-04 -5.70E-03  7.46E-04  3.57E-06
 
 TH 8
+       -8.62E+00  0.00E+00 -8.16E-01 -9.53E+00 -4.14E+00 -1.77E+00  1.05E-04  2.11E+00
 
 TH 9
+       -9.34E+01  0.00E+00  1.73E+01  1.05E+02 -5.66E+01  5.99E+00  3.42E-02 -1.99E+00  3.46E+02
 
 TH10
+       -1.88E+01  0.00E+00  1.68E+01 -1.09E+01 -9.05E+01 -6.70E+00  5.86E-04  3.72E+00  4.96E-01  1.52E+01
 
 TH11
+       -1.09E+02  0.00E+00 -1.26E+01 -5.70E+01 -3.43E+01 -3.96E+01  4.44E-03  2.30E+01  1.78E+01  3.84E+01  2.62E+02
 
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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.03E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.43E+00  0.00E+00  5.42E+01
 
 TH 4
+       -7.60E+00  0.00E+00 -1.19E+01  5.56E+02
 
 TH 5
+       -1.55E+00  0.00E+00 -1.64E+02 -5.46E+01  7.80E+02
 
 TH 6
+        8.73E-01  0.00E+00  1.11E-01 -2.34E+00 -1.89E+00  1.73E+02
 
 TH 7
+        1.94E-03  0.00E+00  1.02E-03 -9.58E-03 -8.76E-05 -2.77E-05  1.46E+00
 
 TH 8
+        1.88E-01  0.00E+00 -1.43E+01 -2.44E+00 -5.49E+00 -5.84E-02  7.80E-04  1.81E+01
 
 TH 9
+        2.73E+00  0.00E+00  5.12E+00 -6.16E-01 -1.33E+00 -6.17E-01  2.83E-02  9.20E-01  2.73E+02
 
 TH10
+       -4.71E-02  0.00E+00  6.82E-01 -6.04E-01 -8.67E+01  9.33E-01  2.75E-03  1.52E+01  1.47E+00  6.85E+01
 
 TH11
+       -8.07E+00  0.00E+00 -7.29E+00 -6.77E+00 -9.55E+00  2.44E+00  1.50E-03  9.93E+00  5.70E+00  1.36E+01  1.81E+02
 
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
+        1.03E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        2.09E+01  0.00E+00  5.28E+01
 
 TH 4
+       -2.78E+01  0.00E+00 -7.69E+00  5.55E+02
 
 TH 5
+       -1.57E+01  0.00E+00 -1.70E+02 -5.20E+01  7.90E+02
 
 TH 6
+       -4.77E+01  0.00E+00  6.25E+00  2.66E+01 -9.04E+01  2.31E+02
 
 TH 7
+       -6.95E-03  0.00E+00  1.89E-03 -2.61E-02  8.10E-03 -2.07E-03  1.20E-05
 
 TH 8
+       -3.52E+01  0.00E+00 -1.20E+01 -1.71E+01 -8.69E+00  1.30E+01  7.73E-04  1.72E+01
 
 TH 9
+        7.23E+01  0.00E+00  6.89E+00 -9.11E+01  4.60E+01 -1.56E+01  4.17E-02 -8.97E+00  2.32E+02
 
 TH10
+       -2.95E+01  0.00E+00  1.12E+00  1.82E+00 -7.89E+01  3.38E+01  9.56E-04  1.40E+01 -1.72E+01  5.46E+01
 
 TH11
+        6.30E+01  0.00E+00 -1.08E+01  3.80E+01 -1.21E+01  3.16E+01 -4.92E-04  9.77E+00  1.52E+00  1.18E+01  1.33E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       48.728
Stop Time:
Wed Sep 29 18:55:04 CDT 2021
