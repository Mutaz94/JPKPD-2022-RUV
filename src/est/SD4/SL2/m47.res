Sun Oct 24 03:16:35 CDT 2021
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
$DATA ../../../../data/SD4/SL2/dat47.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m47.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1624.99010657629        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2978E+02 -2.5955E+01  1.1791E+01 -5.4603E+01  4.5241E+01  6.2572E+01 -3.0973E+00 -7.6154E+00 -3.4710E+01 -1.4696E+01
             1.4154E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1634.37755918635        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      174
 NPARAMETR:  9.8839E-01  9.7812E-01  9.5289E-01  1.0792E+00  9.3337E-01  9.2790E-01  1.0115E+00  1.0300E+00  1.1214E+00  1.0433E+00
             9.9125E-01
 PARAMETER:  8.8326E-02  7.7877E-02  5.1742E-02  1.7623E-01  3.1047E-02  2.5165E-02  1.1141E-01  1.2960E-01  2.1458E-01  1.4239E-01
             9.1210E-02
 GRADIENT:   5.8739E+00  1.3437E+00  2.6350E+00 -8.8359E+00 -9.8144E-01 -2.5710E+00  3.3168E+00 -7.1446E-01  3.9692E-01  3.0148E+00
             3.1714E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1635.01057012437        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.8770E-01  8.4345E-01  9.8364E-01  1.1849E+00  8.8227E-01  9.3614E-01  9.2412E-01  1.0301E+00  1.0901E+00  1.0075E+00
             9.8025E-01
 PARAMETER:  8.7620E-02 -7.0249E-02  8.3506E-02  2.6969E-01 -2.5262E-02  3.4014E-02  2.1091E-02  1.2970E-01  1.8625E-01  1.0747E-01
             8.0052E-02
 GRADIENT:   6.1686E+00  1.6542E+01  5.3140E+00  2.1347E+01 -1.1032E+01  1.1592E+00  2.0665E-01 -1.1727E+00  3.1405E+00  1.6081E-01
            -8.8561E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1635.77666531323        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      531
 NPARAMETR:  9.8173E-01  6.6040E-01  1.1187E+00  1.3000E+00  8.8127E-01  9.2761E-01  8.3328E-01  1.1166E+00  1.0029E+00  1.0248E+00
             9.8366E-01
 PARAMETER:  8.1558E-02 -3.1491E-01  2.1219E-01  3.6239E-01 -2.6396E-02  2.4858E-02 -8.2381E-02  2.1026E-01  1.0285E-01  1.2452E-01
             8.3530E-02
 GRADIENT:  -2.7393E+00  1.3747E+01  9.4147E+00  1.8718E+01 -1.1528E+01 -1.6503E+00 -6.7187E-01 -1.6643E+00 -4.5741E+00 -1.3286E+00
             4.1109E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1636.75192703096        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      707
 NPARAMETR:  9.7966E-01  3.6597E-01  1.0893E+00  1.4687E+00  7.8708E-01  9.2733E-01  7.1140E-01  1.0634E+00  9.1504E-01  9.9775E-01
             9.7462E-01
 PARAMETER:  7.9453E-02 -9.0521E-01  1.8551E-01  4.8437E-01 -1.3943E-01  2.4554E-02 -2.4052E-01  1.6146E-01  1.1210E-02  9.7751E-02
             7.4295E-02
 GRADIENT:   2.3869E+00  3.5089E+00 -1.0222E+00  1.0638E+01 -1.8950E+00 -1.1275E+00 -3.7941E-01  1.8767E-01  5.6159E-01  1.1339E+00
            -5.8806E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1636.78829408891        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      882
 NPARAMETR:  9.7744E-01  2.9140E-01  1.1044E+00  1.5160E+00  7.7131E-01  9.2718E-01  7.0770E-01  1.0806E+00  8.8257E-01  9.8438E-01
             9.7556E-01
 PARAMETER:  7.7187E-02 -1.1331E+00  1.9927E-01  5.1607E-01 -1.5967E-01  2.4393E-02 -2.4573E-01  1.7753E-01 -2.4916E-02  8.4255E-02
             7.5254E-02
 GRADIENT:   2.9930E-01  3.9150E+00  2.9141E+00  1.4439E+01 -5.8024E+00 -1.0176E+00 -2.8203E-01 -6.9711E-01 -1.9710E+00  3.6184E-01
            -2.9099E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1636.82032688664        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1057
 NPARAMETR:  9.7572E-01  2.3058E-01  1.1222E+00  1.5523E+00  7.6321E-01  9.2784E-01  7.1839E-01  1.1104E+00  8.5999E-01  9.7520E-01
             9.7610E-01
 PARAMETER:  7.5419E-02 -1.3672E+00  2.1525E-01  5.3975E-01 -1.7022E-01  2.5108E-02 -2.3075E-01  2.0475E-01 -5.0838E-02  7.4883E-02
             7.5813E-02
 GRADIENT:  -8.3484E-01  3.0810E+00  3.8108E+00  1.2208E+01 -5.8096E+00 -6.0198E-01 -1.7742E-01 -8.5505E-01 -2.6278E+00 -1.5691E-01
            -6.8773E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1636.85004857878        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1232
 NPARAMETR:  9.7409E-01  1.6446E-01  1.1397E+00  1.5911E+00  7.5483E-01  9.2835E-01  7.4809E-01  1.1438E+00  8.3823E-01  9.6733E-01
             9.7627E-01
 PARAMETER:  7.3753E-02 -1.7051E+00  2.3081E-01  5.6445E-01 -1.8127E-01  2.5650E-02 -1.9024E-01  2.3432E-01 -7.6463E-02  6.6782E-02
             7.5989E-02
 GRADIENT:  -1.2374E+00  2.0573E+00  3.6085E+00  8.3524E+00 -4.7042E+00 -2.7082E-01 -8.7938E-02 -8.0122E-01 -2.4507E+00 -3.7125E-01
             6.0914E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1636.98897143247        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1373
 NPARAMETR:  9.7297E-01  1.5638E-01  1.1371E+00  1.5849E+00  7.5529E-01  9.2797E-01  1.2323E+00  1.1469E+00  8.3904E-01  9.6818E-01
             9.7610E-01
 PARAMETER:  7.2596E-02 -1.7555E+00  2.2852E-01  5.6052E-01 -1.8065E-01  2.5249E-02  3.0892E-01  2.3709E-01 -7.5494E-02  6.7666E-02
             7.5806E-02
 GRADIENT:  -3.4658E+00  2.4207E-01  2.2894E-01 -1.3455E+01  2.0573E+00 -4.1192E-01 -4.0299E-02 -7.9340E-02  1.6357E-01 -1.4688E-01
             3.3190E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1637.00180028369        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1556             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7486E-01  1.5321E-01  1.1353E+00  1.5855E+00  7.5325E-01  9.2918E-01  1.3952E+00  1.1470E+00  8.3652E-01  9.6649E-01
             9.7522E-01
 PARAMETER:  7.4542E-02 -1.7760E+00  2.2686E-01  5.6089E-01 -1.8336E-01  2.6547E-02  4.3305E-01  2.3718E-01 -7.8504E-02  6.5912E-02
             7.4910E-02
 GRADIENT:   4.0396E+02  1.8470E+01  8.1154E+00  8.7377E+02  1.8388E+01  3.8478E+01  6.0937E-01  1.1663E+00  1.6939E+01  1.0449E+00
             7.7983E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1637.00715434469        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1738             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7517E-01  1.5185E-01  1.1299E+00  1.5851E+00  7.5118E-01  9.2941E-01  1.4288E+00  1.1431E+00  8.3683E-01  9.6361E-01
             9.7510E-01
 PARAMETER:  7.4859E-02 -1.7849E+00  2.2216E-01  5.6063E-01 -1.8611E-01  2.6794E-02  4.5685E-01  2.3372E-01 -7.8140E-02  6.2926E-02
             7.4780E-02
 GRADIENT:   4.0474E+02  1.8153E+01  7.2512E+00  8.7286E+02  1.9878E+01  3.8520E+01  6.6699E-01  1.2076E+00  1.7238E+01  8.5038E-01
             7.7231E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1637.01265753374        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1921             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7496E-01  1.5182E-01  1.1276E+00  1.5862E+00  7.4813E-01  9.2930E-01  1.4759E+00  1.1403E+00  8.3614E-01  9.6466E-01
             9.7512E-01
 PARAMETER:  7.4638E-02 -1.7851E+00  2.2005E-01  5.6135E-01 -1.9018E-01  2.6675E-02  4.8924E-01  2.3130E-01 -7.8954E-02  6.4019E-02
             7.4807E-02
 GRADIENT:   4.0398E+02  1.8541E+01  9.1250E+00  8.7686E+02  1.6554E+01  3.8475E+01  7.6081E-01  1.2805E+00  1.7077E+01  1.4620E+00
             8.3918E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1637.01452228680        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     2078
 NPARAMETR:  9.7291E-01  1.4391E-01  1.1185E+00  1.5917E+00  7.4224E-01  9.2819E-01  1.6072E+00  1.1325E+00  8.3329E-01  9.6202E-01
             9.7490E-01
 PARAMETER:  7.2531E-02 -1.8385E+00  2.1195E-01  5.6479E-01 -1.9808E-01  2.5485E-02  5.7450E-01  2.2442E-01 -8.2379E-02  6.1276E-02
             7.4576E-02
 GRADIENT:  -3.0370E+00  5.8800E-01  1.3348E+00 -1.1119E+01 -1.4632E+00 -3.2187E-01  3.0425E-02  2.5262E-01 -2.3370E-01  6.1214E-01
             5.2715E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1637.01517549633        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2257
 NPARAMETR:  9.7207E-01  1.3325E-01  1.1077E+00  1.5976E+00  7.3482E-01  9.2781E-01  1.7723E+00  1.1232E+00  8.2981E-01  9.5753E-01
             9.7466E-01
 PARAMETER:  7.1668E-02 -1.9155E+00  2.0224E-01  5.6853E-01 -2.0813E-01  2.5075E-02  6.7228E-01  2.1614E-01 -8.6556E-02  5.6605E-02
             7.4334E-02
 GRADIENT:  -4.7057E+00  6.8437E-01  1.1823E+00 -9.5633E+00 -1.8972E+00 -4.8190E-01  5.1888E-02  2.3365E-01 -3.3275E-01  7.0274E-01
             5.1761E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1637.01541830954        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2437
 NPARAMETR:  9.7171E-01  1.2588E-01  1.1014E+00  1.6015E+00  7.3029E-01  9.2767E-01  1.8684E+00  1.1181E+00  8.2763E-01  9.5458E-01
             9.7452E-01
 PARAMETER:  7.1301E-02 -1.9724E+00  1.9661E-01  5.7096E-01 -2.1432E-01  2.4924E-02  7.2508E-01  2.1167E-01 -8.9194E-02  5.3512E-02
             7.4193E-02
 GRADIENT:  -5.2518E+00  6.9954E-01  1.0785E+00 -9.0729E+00 -2.0642E+00 -5.3557E-01  5.8277E-02  2.1482E-01 -3.6827E-01  7.1106E-01
             5.0110E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1637.03963131398        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2621
 NPARAMETR:  9.7383E-01  1.1799E-01  1.0990E+00  1.6021E+00  7.2915E-01  9.2890E-01  1.7362E+00  1.1155E+00  8.2754E-01  9.4978E-01
             9.7444E-01
 PARAMETER:  7.3479E-02 -2.0372E+00  1.9441E-01  5.7128E-01 -2.1587E-01  2.6243E-02  6.5169E-01  2.0927E-01 -8.9295E-02  4.8480E-02
             7.4106E-02
 GRADIENT:   8.5082E-01  1.5003E-01 -3.8644E-01 -1.6437E+01  1.7990E+00  3.0931E-02  1.1026E-02 -1.3186E-01  2.3554E-01 -1.9218E-01
            -3.8228E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1637.04129461495        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2796
 NPARAMETR:  9.7406E-01  1.1809E-01  1.0985E+00  1.6014E+00  7.2843E-01  9.2905E-01  1.7059E+00  1.1166E+00  8.2734E-01  9.4973E-01
             9.7445E-01
 PARAMETER:  7.3721E-02 -2.0363E+00  1.9396E-01  5.7088E-01 -2.1686E-01  2.6411E-02  6.3408E-01  2.1033E-01 -8.9541E-02  4.8418E-02
             7.4120E-02
 GRADIENT:   2.3648E-02 -5.3287E-02  7.5032E-02 -2.0500E-01  9.1353E-01  2.8174E-03 -2.4583E-04  2.0218E-02  1.2470E-02 -3.2109E-02
             4.7433E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2796
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4522E-04 -3.5029E-03 -2.8879E-02 -4.8369E-03 -3.3318E-02
 SE:             2.9841E-02  3.1327E-03  1.8844E-02  2.9219E-02  2.0995E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9344E-01  2.6349E-01  1.2538E-01  8.6852E-01  1.1252E-01

 ETASHRINKSD(%)  2.9170E-02  8.9505E+01  3.6872E+01  2.1116E+00  2.9665E+01
 ETASHRINKVR(%)  5.8331E-02  9.8899E+01  6.0148E+01  4.1786E+00  5.0530E+01
 EBVSHRINKSD(%)  4.4195E-01  8.9935E+01  3.9128E+01  2.4394E+00  2.7382E+01
 EBVSHRINKVR(%)  8.8194E-01  9.8987E+01  6.2946E+01  4.8192E+00  4.7267E+01
 RELATIVEINF(%)  9.4248E+01  5.4535E-02  6.5145E+00  7.1798E+00  4.8742E+00
 EPSSHRINKSD(%)  4.6149E+01
 EPSSHRINKVR(%)  7.1001E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1637.0412946149540     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -901.89046805121586     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.21
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1637.041       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.74E-01  1.18E-01  1.10E+00  1.60E+00  7.28E-01  9.29E-01  1.71E+00  1.12E+00  8.27E-01  9.50E-01  9.74E-01
 


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
 #CPUT: Total CPU Time in Seconds,       78.946
Stop Time:
Sun Oct 24 03:16:50 CDT 2021
