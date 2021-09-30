Wed Sep 29 15:41:44 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat25.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m25.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1679.62430576556        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2556E+02 -2.8141E+01  5.1946E+00 -4.8991E+01 -2.4673E+01  6.4060E+01  5.4554E+00  8.6821E+00  3.8027E+00  8.6201E+00
             1.5278E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1685.42873385927        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  9.8778E-01  1.0415E+00  1.0620E+00  1.0473E+00  1.0833E+00  8.5608E-01  9.7892E-01  9.5151E-01  9.8031E-01  9.7591E-01
             1.0021E+00
 PARAMETER:  8.7706E-02  1.4063E-01  1.6014E-01  1.4620E-01  1.7998E-01 -5.5395E-02  7.8692E-02  5.0294E-02  8.0116E-02  7.5620E-02
             1.0206E-01
 GRADIENT:  -2.3758E+00 -9.7309E+00 -1.1431E+01 -2.5846E+00  2.4460E+01 -2.6642E+01  7.7210E-01  1.3508E+00 -3.3252E+00 -1.1164E+01
            -3.1700E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1686.13390885562        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  9.8668E-01  9.3831E-01  1.1890E+00  1.1238E+00  1.0927E+00  8.7298E-01  7.9964E-01  9.5313E-01  9.9501E-01  1.0406E+00
             1.0023E+00
 PARAMETER:  8.6588E-02  3.6324E-02  2.7312E-01  2.1675E-01  1.8865E-01 -3.5847E-02 -1.2360E-01  5.2000E-02  9.4996E-02  1.3978E-01
             1.0230E-01
 GRADIENT:  -1.5136E+00  2.4929E+00 -2.1794E+00  9.2986E+00  1.3331E+01 -1.7643E+01 -3.2454E+00 -3.7601E+00 -2.2399E+00 -8.7030E+00
            -4.3745E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1687.12586911517        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.8639E-01  9.1741E-01  1.2977E+00  1.1383E+00  1.1247E+00  9.1220E-01  8.7977E-01  1.1276E+00  9.7737E-01  1.1130E+00
             1.0058E+00
 PARAMETER:  8.6296E-02  1.3800E-02  3.6059E-01  2.2950E-01  2.1753E-01  8.0986E-03 -2.8100E-02  2.2013E-01  7.7115E-02  2.0705E-01
             1.0581E-01
 GRADIENT:  -1.4543E-01  2.4399E+00  8.2028E-01  1.7082E+00 -2.2751E+00  6.4269E-01 -1.8205E-01 -4.5057E-02 -1.2063E-01  3.3776E-01
             7.8931E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1687.15990734469        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  9.8523E-01  8.0352E-01  1.4069E+00  1.2171E+00  1.1195E+00  9.0963E-01  9.2750E-01  1.1903E+00  9.3294E-01  1.1191E+00
             1.0044E+00
 PARAMETER:  8.5116E-02 -1.1875E-01  4.4136E-01  2.9648E-01  2.1293E-01  5.2805E-03  2.4738E-02  2.7418E-01  3.0590E-02  2.1254E-01
             1.0443E-01
 GRADIENT:   2.8586E-01  5.8245E+00  3.3719E+00  7.4178E+00 -4.6188E+00  8.0994E-02  2.5912E-01 -5.2741E-01  4.0173E-01 -3.7377E-01
            -5.7064E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1687.21077297468        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      877
 NPARAMETR:  9.8378E-01  6.9147E-01  1.5046E+00  1.2915E+00  1.1156E+00  9.0600E-01  9.7479E-01  1.2540E+00  8.9333E-01  1.1264E+00
             1.0042E+00
 PARAMETER:  8.3651E-02 -2.6894E-01  5.0853E-01  3.5583E-01  2.0936E-01  1.2895E-03  7.4463E-02  3.2633E-01 -1.2797E-02  2.1906E-01
             1.0423E-01
 GRADIENT:   7.1646E-02  5.5440E+00  3.1584E+00  9.1925E+00 -3.5552E+00 -9.0073E-01  4.2450E-01 -5.7507E-01  6.7873E-01 -7.0944E-01
            -5.7501E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1687.32749901610        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1052
 NPARAMETR:  9.8096E-01  4.9153E-01  1.6801E+00  1.4234E+00  1.1106E+00  9.0129E-01  1.0695E+00  1.3929E+00  8.2797E-01  1.1413E+00
             1.0044E+00
 PARAMETER:  8.0780E-02 -6.1024E-01  6.1886E-01  4.5306E-01  2.0487E-01 -3.9265E-03  1.6716E-01  4.3137E-01 -8.8781E-02  2.3214E-01
             1.0435E-01
 GRADIENT:  -5.0719E-01  4.6059E+00  1.7642E+00  1.1196E+01 -1.5122E+00 -1.7878E+00  3.3224E-01 -4.6213E-01  3.6481E-01 -6.6999E-01
            -3.3394E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1687.53374646335        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1228
 NPARAMETR:  9.7718E-01  2.6040E-01  1.9643E+00  1.5810E+00  1.1296E+00  8.9447E-01  1.2736E+00  1.6647E+00  7.5691E-01  1.1705E+00
             1.0050E+00
 PARAMETER:  7.6918E-02 -1.2456E+00  7.7514E-01  5.5804E-01  2.2187E-01 -1.1518E-02  3.4182E-01  6.0964E-01 -1.7851E-01  2.5746E-01
             1.0497E-01
 GRADIENT:  -2.2581E+00  4.0011E+00 -1.2237E+00  1.7685E+01  1.3713E+00 -3.2885E+00  9.2473E-02  1.9014E-01 -5.1335E-01 -4.0928E-01
             8.9203E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1687.63223464604        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1403
 NPARAMETR:  9.7580E-01  1.5744E-01  2.1516E+00  1.6577E+00  1.1463E+00  8.9610E-01  1.4516E+00  1.8301E+00  7.2586E-01  1.1897E+00
             1.0043E+00
 PARAMETER:  7.5499E-02 -1.7487E+00  8.6619E-01  6.0542E-01  2.3656E-01 -9.6981E-03  4.7268E-01  7.0435E-01 -2.2040E-01  2.7368E-01
             1.0430E-01
 GRADIENT:  -2.4987E+00  3.7807E+00 -7.6518E-01  3.1877E+01 -1.2273E+00 -1.9012E+00 -1.1663E-02 -3.3962E-02 -1.4899E+00  1.0503E-01
            -2.8842E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1687.97009716285        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1580
 NPARAMETR:  9.7529E-01  7.4285E-02  2.3401E+00  1.7115E+00  1.1647E+00  8.9989E-01  1.7191E+00  1.9813E+00  6.9948E-01  1.2051E+00
             1.0040E+00
 PARAMETER:  7.4984E-02 -2.4998E+00  9.5017E-01  6.3735E-01  2.5247E-01 -5.4777E-03  6.4179E-01  7.8375E-01 -2.5742E-01  2.8660E-01
             1.0397E-01
 GRADIENT:  -1.0094E+00  1.6971E+00  4.0736E-01  2.1904E+01 -2.2814E+00  3.5671E-01 -1.9797E-02 -6.1221E-01 -2.6922E+00  2.2908E-01
            -4.1838E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1688.24560869286        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1758
 NPARAMETR:  9.7528E-01  1.9870E-02  2.5918E+00  1.7457E+00  1.2015E+00  8.9848E-01  2.1323E+00  2.1859E+00  6.8392E-01  1.2275E+00
             1.0043E+00
 PARAMETER:  7.4974E-02 -3.8185E+00  1.0523E+00  6.5717E-01  2.8355E-01 -7.0543E-03  8.5722E-01  8.8201E-01 -2.7991E-01  3.0495E-01
             1.0428E-01
 GRADIENT:   3.1785E-01  3.4129E-01  1.0757E+00  3.4963E+00 -1.2290E+00  2.2815E-01 -2.4722E-03 -3.8825E-01 -1.3592E+00 -2.1596E-01
            -1.7034E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1688.38930923742        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1946
 NPARAMETR:  9.7539E-01  1.0000E-02  2.5655E+00  1.7411E+00  1.1975E+00  8.9783E-01  2.9032E+00  2.1820E+00  6.8233E-01  1.2290E+00
             1.0030E+00
 PARAMETER:  7.5078E-02 -4.7984E+00  1.0421E+00  6.5453E-01  2.8022E-01 -7.7764E-03  1.1658E+00  8.8025E-01 -2.8225E-01  3.0624E-01
             1.0295E-01
 GRADIENT:   1.4912E+00  0.0000E+00  1.9970E-02 -2.3596E+01  6.9038E-01  5.6140E-02 -2.8933E-04  5.6228E-01 -8.4554E-01  2.7024E-01
            -2.4816E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1688.39771147040        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2137
 NPARAMETR:  9.7531E-01  1.0000E-02  2.5544E+00  1.7405E+00  1.1932E+00  8.9787E-01  3.3018E+00  2.1639E+00  6.8461E-01  1.2259E+00
             1.0034E+00
 PARAMETER:  7.5000E-02 -4.7984E+00  1.0378E+00  6.5418E-01  2.7664E-01 -7.7279E-03  1.2945E+00  8.7191E-01 -2.7891E-01  3.0364E-01
             1.0340E-01
 GRADIENT:   1.3661E+00  0.0000E+00  8.1913E-01 -2.4122E+01 -1.6782E-01  7.4930E-02  4.0973E-04  1.0716E-01  2.5644E-01  2.7751E-01
            -4.0518E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1688.40386974725        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2335             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7534E-01  1.0000E-02  2.5321E+00  1.7401E+00  1.1912E+00  8.9789E-01  3.5067E+00  2.1539E+00  6.8480E-01  1.2212E+00
             1.0035E+00
 PARAMETER:  7.5029E-02 -4.7984E+00  1.0290E+00  6.5392E-01  2.7494E-01 -7.7056E-03  1.3547E+00  8.6728E-01 -2.7862E-01  2.9984E-01
             1.0346E-01
 GRADIENT:   4.0000E+02  0.0000E+00  7.4654E+00  1.1829E+03  1.3579E+01  3.2002E+01  5.2427E-02  4.3531E+00  3.3801E+01  2.6686E+00
             8.1561E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1688.40571769812        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     2508
 NPARAMETR:  9.7528E-01  1.0000E-02  2.5319E+00  1.7400E+00  1.1887E+00  8.9788E-01  3.5538E+00  2.1481E+00  6.8488E-01  1.2226E+00
             1.0035E+00
 PARAMETER:  7.4969E-02 -4.7984E+00  1.0290E+00  6.5387E-01  2.7282E-01 -7.7171E-03  1.3680E+00  8.6458E-01 -2.7851E-01  3.0097E-01
             1.0349E-01
 GRADIENT:   1.1884E+00  0.0000E+00  6.6217E-01 -2.3908E+01 -3.1155E-01  6.0744E-02  6.9703E-04  2.1536E-01  2.2190E-01  2.8768E-01
             1.7284E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1688.41018019448        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     2706             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7530E-01  1.0000E-02  2.5096E+00  1.7395E+00  1.1870E+00  8.9792E-01  3.5397E+00  2.1368E+00  6.8514E-01  1.2164E+00
             1.0034E+00
 PARAMETER:  7.4992E-02 -4.7984E+00  1.0201E+00  6.5359E-01  2.7146E-01 -7.6748E-03  1.3640E+00  8.5929E-01 -2.7813E-01  2.9592E-01
             1.0342E-01
 GRADIENT:   3.9984E+02  0.0000E+00  7.2828E+00  1.1816E+03  1.3857E+01  3.2002E+01  5.3708E-02  4.2650E+00  3.3765E+01  2.3817E+00
             7.8508E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1688.41086697489        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     2877
 NPARAMETR:  9.7526E-01  1.0000E-02  2.5112E+00  1.7404E+00  1.1838E+00  8.9794E-01  2.6449E+00  2.1309E+00  6.8520E-01  1.2203E+00
             1.0035E+00
 PARAMETER:  7.4951E-02 -4.7984E+00  1.0208E+00  6.5411E-01  2.6870E-01 -7.6557E-03  1.0726E+00  8.5655E-01 -2.7804E-01  2.9912E-01
             1.0351E-01
 GRADIENT:   1.4134E+00  0.0000E+00  7.1239E-01 -2.1451E+01 -9.5974E-01  5.9569E-02 -1.3494E-04  1.9311E-01  1.8805E-01  4.5564E-01
             2.2854E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1688.41504147878        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     3050
 NPARAMETR:  9.7518E-01  1.0000E-02  2.5021E+00  1.7389E+00  1.1833E+00  8.9789E-01  3.5763E+00  2.1266E+00  6.8536E-01  1.2175E+00
             1.0035E+00
 PARAMETER:  7.4871E-02 -4.7984E+00  1.0171E+00  6.5323E-01  2.6834E-01 -7.7048E-03  1.3743E+00  8.5451E-01 -2.7780E-01  2.9683E-01
             1.0349E-01
 GRADIENT:   1.1764E+00  0.0000E+00  3.3760E-01 -2.4931E+01  1.7289E-01  6.2358E-02  7.7380E-04  2.8719E-01  2.3506E-01  1.1754E-01
             4.0018E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1688.41820203904        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     3216
 NPARAMETR:  9.7522E-01  1.0000E-02  2.4912E+00  1.7387E+00  1.1810E+00  8.9794E-01  3.6972E+00  2.1162E+00  6.8556E-01  1.2166E+00
             1.0035E+00
 PARAMETER:  7.4906E-02 -4.7984E+00  1.0128E+00  6.5316E-01  2.6640E-01 -7.6522E-03  1.4076E+00  8.4964E-01 -2.7752E-01  2.9604E-01
             1.0348E-01
 GRADIENT:   3.9952E+02  0.0000E+00  7.9599E+00  1.1797E+03  1.2021E+01  3.1983E+01  6.0055E-02  3.9778E+00  3.3738E+01  2.9624E+00
             8.6258E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1688.41943317078        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     3384
 NPARAMETR:  9.7516E-01  1.0000E-02  2.4787E+00  1.7389E+00  1.1805E+00  8.9791E-01  3.4840E+00  2.1109E+00  6.8555E-01  1.2133E+00
             1.0033E+00
 PARAMETER:  7.4845E-02 -4.7984E+00  1.0077E+00  6.5327E-01  2.6595E-01 -7.6858E-03  1.3482E+00  8.4713E-01 -2.7753E-01  2.9335E-01
             1.0330E-01
 GRADIENT:   1.2736E+00  0.0000E+00 -3.1316E-01 -2.3469E+01  1.4137E+00  5.8425E-02  5.9644E-04  3.3236E-01  1.5828E-01 -2.1746E-01
            -6.9469E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1688.42099578306        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     3551
 NPARAMETR:  9.7521E-01  1.0000E-02  2.4792E+00  1.7389E+00  1.1775E+00  8.9798E-01  3.4388E+00  2.1048E+00  6.8578E-01  1.2147E+00
             1.0035E+00
 PARAMETER:  7.4897E-02 -4.7984E+00  1.0079E+00  6.5328E-01  2.6343E-01 -7.6092E-03  1.3351E+00  8.4420E-01 -2.7720E-01  2.9452E-01
             1.0347E-01
 GRADIENT:   1.2931E+00  0.0000E+00  5.9581E-01 -2.2938E+01 -6.1889E-01  7.6573E-02  5.4016E-04  1.1209E-01  2.3295E-01  2.9510E-01
             1.8668E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1688.42278717130        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3740
 NPARAMETR:  9.7466E-01  1.0000E-02  2.4704E+00  1.7387E+00  1.1762E+00  8.9771E-01  3.1239E+00  2.0985E+00  6.8561E-01  1.2127E+00
             1.0034E+00
 PARAMETER:  7.4335E-02 -4.7984E+00  1.0044E+00  6.5314E-01  2.6226E-01 -7.9047E-03  1.2391E+00  8.4123E-01 -2.7745E-01  2.9286E-01
             1.0342E-01
 GRADIENT:  -9.4808E-02  0.0000E+00  4.3638E-01 -2.3003E+01 -2.7818E-01 -5.2051E-02  1.6574E-04  1.2207E-01  6.3327E-02  1.5573E-01
            -6.5059E-03

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1688.42491334505        NO. OF FUNC. EVALS.: 201
 CUMULATIVE NO. OF FUNC. EVALS.:     3941             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7521E-01  1.0000E-02  2.4511E+00  1.7379E+00  1.1751E+00  8.9800E-01  3.4326E+00  2.0908E+00  6.8610E-01  1.2093E+00
             1.0034E+00
 PARAMETER:  7.4895E-02 -4.7984E+00  9.9655E-01  6.5269E-01  2.6139E-01 -7.5840E-03  1.3333E+00  8.3753E-01 -2.7673E-01  2.9002E-01
             1.0339E-01
 GRADIENT:   3.9942E+02  0.0000E+00  7.0265E+00  1.1780E+03  1.3351E+01  3.1957E+01  4.9563E-02  4.0981E+00  3.3651E+01  2.4613E+00
             8.2853E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1688.42544353669        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     4110
 NPARAMETR:  9.7516E-01  1.0000E-02  2.4568E+00  1.7383E+00  1.1721E+00  8.9800E-01  3.4921E+00  2.0848E+00  6.8618E-01  1.2113E+00
             1.0034E+00
 PARAMETER:  7.4851E-02 -4.7984E+00  9.9887E-01  6.5292E-01  2.5878E-01 -7.5861E-03  1.3505E+00  8.3466E-01 -2.7661E-01  2.9171E-01
             1.0338E-01
 GRADIENT:   1.4123E+00  0.0000E+00  7.7659E-01 -2.2923E+01 -1.1327E+00  7.9501E-02  6.0012E-04  1.3065E-02  2.5417E-01  3.7169E-01
            -1.4609E-02

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1688.42739890781        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     4285
 NPARAMETR:  9.7513E-01  1.0000E-02  2.4501E+00  1.7375E+00  1.1714E+00  8.9799E-01  3.5278E+00  2.0821E+00  6.8623E-01  1.2107E+00
             1.0034E+00
 PARAMETER:  7.4815E-02 -4.7984E+00  9.9612E-01  6.5246E-01  2.5823E-01 -7.6005E-03  1.3607E+00  8.3340E-01 -2.7654E-01  2.9116E-01
             1.0344E-01
 GRADIENT:   1.3239E+00  0.0000E+00  4.9833E-01 -2.4581E+01 -7.0014E-01  6.9663E-02  6.8809E-04  1.5482E-01  2.3558E-01  3.3799E-01
             4.4607E-02

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1688.42832199104        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     4477
 NPARAMETR:  9.7518E-01  1.0000E-02  2.4452E+00  1.7376E+00  1.1705E+00  8.9802E-01  3.4924E+00  2.0777E+00  6.8637E-01  1.2100E+00
             1.0035E+00
 PARAMETER:  7.4863E-02 -4.7984E+00  9.9415E-01  6.5249E-01  2.5739E-01 -7.5612E-03  1.3506E+00  8.3128E-01 -2.7635E-01  2.9059E-01
             1.0344E-01
 GRADIENT:   1.5140E+00  0.0000E+00  4.9072E-01 -2.4158E+01 -6.8898E-01  8.6062E-02  6.5725E-04  1.1709E-01  2.6318E-01  3.3416E-01
             4.0269E-02

0ITERATION NO.:  129    OBJECTIVE VALUE:  -1688.42884842775        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     4622
 NPARAMETR:  9.7517E-01  1.0000E-02  2.4416E+00  1.7373E+00  1.1689E+00  8.9804E-01  3.6116E+00  2.0727E+00  6.8646E-01  1.2100E+00
             1.0035E+00
 PARAMETER:  7.4854E-02 -4.7984E+00  9.9150E-01  6.5255E-01  2.5818E-01 -7.5703E-03  1.3844E+00  8.3074E-01 -2.7634E-01  2.8775E-01
             1.0335E-01
 GRADIENT:  -1.1827E-02  0.0000E+00 -1.5195E-01  5.9864E-01  5.3520E-01 -7.4761E-03  2.4510E-06  1.5845E-01 -3.1991E-02 -1.3724E-01
            -4.0053E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4622
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.7077E-05 -1.1103E-03 -4.0707E-02 -9.4246E-03 -6.0210E-02
 SE:             2.9825E-02  6.2953E-04  1.8795E-02  2.9195E-02  1.9328E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9740E-01  7.7779E-02  3.0324E-02  7.4683E-01  1.8385E-03

 ETASHRINKSD(%)  8.0953E-02  9.7891E+01  3.7034E+01  2.1933E+00  3.5249E+01
 ETASHRINKVR(%)  1.6184E-01  9.9956E+01  6.0352E+01  4.3385E+00  5.8073E+01
 EBVSHRINKSD(%)  5.2694E-01  9.8030E+01  4.0436E+01  2.7366E+00  3.0867E+01
 EBVSHRINKVR(%)  1.0511E+00  9.9961E+01  6.4521E+01  5.3983E+00  5.2206E+01
 RELATIVEINF(%)  9.6289E+01  1.8518E-03  1.4878E+01  4.9438E+00  1.3878E+01
 EPSSHRINKSD(%)  4.4918E+01
 EPSSHRINKVR(%)  6.9659E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1688.4288484277452     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -953.27802186400697     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    64.50
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.11
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1688.429       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.75E-01  1.00E-02  2.44E+00  1.74E+00  1.17E+00  8.98E-01  3.61E+00  2.08E+00  6.86E-01  1.21E+00  1.00E+00
 


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
+        1.43E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -9.73E+00  0.00E+00  5.31E+00
 
 TH 4
+        8.24E+01  0.00E+00 -2.00E+00  7.69E+02
 
 TH 5
+        3.42E+01  0.00E+00 -3.57E+01 -8.74E+01  2.58E+02
 
 TH 6
+       -9.50E+00  0.00E+00 -8.34E-02 -2.55E+01 -1.01E+01  2.90E+02
 
 TH 7
+        4.65E-03  0.00E+00 -6.30E-05 -4.19E-03  4.01E-03  3.45E-03  3.21E-06
 
 TH 8
+       -1.47E+00  0.00E+00  1.47E+00 -1.27E+00 -9.95E+00  2.89E+00  3.43E-05  4.38E-01
 
 TH 9
+        1.88E+01  0.00E+00  2.67E+00  7.07E+01  6.96E+00  3.50E+00  3.30E-02  8.55E-01  3.64E+02
 
 TH10
+       -4.30E+00  0.00E+00  7.47E+00  9.11E+00 -5.32E+01  6.22E+00 -9.47E-04  2.13E+00 -4.49E+00  1.12E+01
 
 TH11
+        4.10E+00  0.00E+00 -6.35E+00 -5.37E+01  4.39E+01  1.00E+02  2.05E-04 -7.30E-01 -2.03E+01 -7.14E+00  4.57E+01
 
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
+        1.43E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.37E+00  0.00E+00  1.74E+01
 
 TH 4
+       -1.36E+01  0.00E+00 -9.68E+00  7.40E+02
 
 TH 5
+        4.99E-01  0.00E+00 -4.11E+01 -3.81E+01  2.91E+02
 
 TH 6
+        2.98E-01  0.00E+00  2.68E-02 -3.59E+00  7.34E-01  2.43E+02
 
 TH 7
+        4.22E-03  0.00E+00  5.68E-04 -9.94E-03  3.43E-03  2.01E-03  3.21E-04
 
 TH 8
+        1.15E-02  0.00E+00 -9.58E+00 -3.06E+00 -1.33E+01  1.73E-01 -3.90E-04  1.51E+01
 
 TH 9
+        4.38E+00  0.00E+00  4.18E+00 -4.14E-01  3.78E+00 -6.38E-01  3.57E-02  1.93E+00  3.82E+02
 
 TH10
+        2.06E+00  0.00E+00  1.23E+00 -5.70E-01 -5.34E+01  7.29E-01  7.47E-04  2.97E+00 -1.01E+00  5.09E+01
 
 TH11
+       -1.28E+01  0.00E+00 -2.45E+00 -1.31E+01 -4.88E+00  1.87E+00  2.39E-03  3.29E+00  1.21E+01  7.80E+00  2.15E+02
 
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
+        1.43E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -8.50E+00  0.00E+00  1.63E+01
 
 TH 4
+       -1.05E+02  0.00E+00 -2.05E+01  7.40E+02
 
 TH 5
+       -6.80E+01  0.00E+00 -4.26E+01  3.80E+01  2.95E+02
 
 TH 6
+       -3.66E+01  0.00E+00 -1.87E+00  2.87E+01 -2.02E+01  1.78E+02
 
 TH 7
+        4.51E-03  0.00E+00  7.17E-04 -2.02E-02 -1.49E-03  1.81E-04  7.12E-06
 
 TH 8
+        1.08E+01  0.00E+00 -7.36E+00 -1.07E+01 -1.19E+01  4.90E+00  2.48E-04  1.09E+01
 
 TH 9
+        2.13E+01  0.00E+00  9.99E-01 -8.33E+01  1.08E+01  1.62E+01  4.58E-02  2.30E+00  3.93E+02
 
 TH10
+        6.30E+01  0.00E+00  2.95E+00 -6.20E-01 -7.05E+01  1.23E+00 -4.63E-04  5.66E+00 -1.37E+01  5.73E+01
 
 TH11
+        8.14E+01  0.00E+00 -3.47E+00 -2.99E+01  3.16E+01  3.75E+00 -1.38E-03 -3.20E+00 -2.04E+01 -2.75E+00  1.56E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       70.657
Stop Time:
Wed Sep 29 15:42:56 CDT 2021
