Sat Sep 25 02:44:58 CDT 2021
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
$DATA ../../../../data/int/SL3/dat84.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      985
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

 TOT. NO. OF OBS RECS:      885
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
 RAW OUTPUT FILE (FILE): m84.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   221.443156395659        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -3.6490E+01 -3.4990E+01  1.3389E+02  8.4440E+01  1.4623E+02  4.2315E+01 -1.3092E+02 -2.6161E+02 -1.0609E+02 -6.6747E+01
            -7.6152E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2379.74602294190        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.1142E+00  1.3276E+00  9.8608E-01  9.2148E-01  1.0735E+00  7.4110E-01  1.2005E+00  1.0096E+00  9.4122E-01  9.7845E-01
             5.1942E+00
 PARAMETER:  2.0814E-01  3.8339E-01  8.5983E-02  1.8230E-02  1.7092E-01 -1.9962E-01  2.8273E-01  1.0954E-01  3.9423E-02  7.8214E-02
             1.7475E+00
 GRADIENT:   2.9464E+01  7.8632E+00 -7.5093E+00  5.2089E+00 -1.8517E+01 -2.9886E+01  1.6730E+01  5.7490E+00  1.8484E+01  1.4368E+01
             7.8462E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2558.57625016905        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0177E+00  1.5693E+00  2.1064E+01  8.1642E-01  3.0378E+00  6.9865E-01  1.1722E+00  2.4568E+00  1.4133E+00  2.0001E+00
             3.4683E+00
 PARAMETER:  1.1753E-01  5.5065E-01  3.1476E+00 -1.0282E-01  1.2111E+00 -2.5861E-01  2.5886E-01  9.9888E-01  4.4593E-01  7.9318E-01
             1.3437E+00
 GRADIENT:  -2.1756E+02  2.1861E+01 -5.5135E+00  2.6337E+01  1.8946E+02 -9.2761E+01  3.7288E+01 -1.1132E-01  3.1329E+01 -4.3127E+01
             3.5128E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2658.35330328026        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0698E+00  1.4820E+00  4.2572E+00  8.0708E-01  1.7634E+00  8.9033E-01  1.0659E+00  3.7484E+00  8.8936E-01  1.5255E+00
             2.9000E+00
 PARAMETER:  1.6743E-01  4.9341E-01  1.5486E+00 -1.1433E-01  6.6722E-01 -1.6164E-02  1.6384E-01  1.4213E+00 -1.7250E-02  5.2232E-01
             1.1647E+00
 GRADIENT:   3.6318E+01  1.2900E+01 -6.9627E+00  2.3099E+01  3.5400E+01  1.2606E+01  3.0602E+00  7.7111E+00 -6.0471E-01  4.8207E+00
             2.0955E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2660.47055151439        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0549E+00  1.4256E+00  3.7919E+00  8.2756E-01  1.5875E+00  8.5842E-01  1.0512E+00  3.1188E+00  9.3385E-01  1.3920E+00
             2.8676E+00
 PARAMETER:  1.5345E-01  4.5461E-01  1.4329E+00 -8.9275E-02  5.6213E-01 -5.2660E-02  1.4997E-01  1.2374E+00  3.1565E-02  4.3077E-01
             1.1535E+00
 GRADIENT:   4.9549E-01  8.5958E+00  3.6222E+00 -7.9031E-02 -9.9318E+00 -1.6840E-01  3.4402E-01 -2.7161E+00 -8.0359E-01  2.9234E-01
            -2.3104E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2660.70305524563        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      417
 NPARAMETR:  1.0585E+00  1.4624E+00  3.6935E+00  8.0525E-01  1.6229E+00  8.5700E-01  1.0377E+00  3.1993E+00  9.7225E-01  1.4077E+00
             2.8716E+00
 PARAMETER:  1.5684E-01  4.8011E-01  1.4066E+00 -1.1660E-01  5.8420E-01 -5.4313E-02  1.3705E-01  1.2629E+00  7.1853E-02  4.4197E-01
             1.1549E+00
 GRADIENT:   7.0900E-01  1.4853E+00  6.8173E-01  2.9036E+00 -4.9133E-03 -1.2451E+00  2.2486E+00 -9.3823E-01  1.2462E+00 -1.9670E-01
            -1.2367E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2661.14966434691        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      593
 NPARAMETR:  1.0569E+00  1.7850E+00  2.7591E+00  5.9294E-01  1.7294E+00  8.5550E-01  8.7101E-01  3.1660E+00  1.2077E+00  1.5380E+00
             2.8798E+00
 PARAMETER:  1.5530E-01  6.7940E-01  1.1149E+00 -4.2267E-01  6.4778E-01 -5.6064E-02 -3.8106E-02  1.2525E+00  2.8868E-01  5.3050E-01
             1.1577E+00
 GRADIENT:  -5.5536E+00  1.1092E+01  4.9713E-02  6.7901E+00  5.8601E-01 -1.8965E+00 -1.8499E+00  4.6107E-03 -3.4302E-01  6.8308E-01
             1.1613E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2661.34157255540        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      768
 NPARAMETR:  1.0598E+00  1.9659E+00  2.0754E+00  4.6652E-01  1.7750E+00  8.6212E-01  8.1607E-01  2.9828E+00  1.4365E+00  1.5922E+00
             2.8792E+00
 PARAMETER:  1.5812E-01  7.7596E-01  8.3018E-01 -6.6246E-01  6.7377E-01 -4.8358E-02 -1.0325E-01  1.1929E+00  4.6221E-01  5.6513E-01
             1.1575E+00
 GRADIENT:   2.5076E+00  9.0423E+00 -5.7713E-02  3.9244E+00 -9.5442E-01  8.7001E-01 -8.5563E-01  2.7753E-01  1.1240E-03  8.3612E-01
            -7.7842E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2661.59903254980        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      948
 NPARAMETR:  1.0580E+00  2.1356E+00  1.3733E+00  3.4470E-01  1.8205E+00  8.5877E-01  7.8185E-01  2.5685E+00  1.7352E+00  1.6248E+00
             2.8816E+00
 PARAMETER:  1.5635E-01  8.5874E-01  4.1720E-01 -9.6508E-01  6.9914E-01 -5.2254E-02 -1.4609E-01  1.0433E+00  6.5112E-01  5.8537E-01
             1.1583E+00
 GRADIENT:  -2.8545E+00  1.7143E+00 -3.5491E-01  1.1186E+00  1.2861E+00 -6.1631E-01  1.0607E+00  7.1261E-01 -4.3718E-01 -2.6686E-01
             3.9866E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2662.62719943983        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1125
 NPARAMETR:  1.0588E+00  2.3747E+00  5.0756E-01  1.7278E-01  1.8768E+00  8.6112E-01  7.2530E-01  1.0887E+00  2.6791E+00  1.6820E+00
             2.8733E+00
 PARAMETER:  1.5718E-01  9.6487E-01 -5.7814E-01 -1.6557E+00  7.2956E-01 -4.9526E-02 -2.2117E-01  1.8495E-01  1.0855E+00  6.1996E-01
             1.1555E+00
 GRADIENT:  -5.7607E-01 -9.0724E+00 -5.0131E-01 -8.2352E-01 -7.1064E-02  1.2172E-02  4.5889E-01  7.3617E-01  4.1036E-01 -4.1328E-02
            -2.3838E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2663.00784594334        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1304
 NPARAMETR:  1.0593E+00  2.3613E+00  5.3577E-01  1.8360E-01  1.8667E+00  8.6177E-01  7.2606E-01  2.9010E-01  2.6388E+00  1.6731E+00
             2.8748E+00
 PARAMETER:  1.5763E-01  9.5922E-01 -5.2406E-01 -1.5950E+00  7.2416E-01 -4.8761E-02 -2.2012E-01 -1.1375E+00  1.0703E+00  6.1470E-01
             1.1560E+00
 GRADIENT:   4.8304E-01 -3.1484E+00  2.2167E-01 -7.2167E-01 -3.7679E-01  2.0795E-01  5.2842E-01  5.6587E-02  7.8141E-02 -6.1972E-02
            -6.6969E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2663.04569379788        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1479
 NPARAMETR:  1.0592E+00  2.3447E+00  5.4785E-01  1.9645E-01  1.8542E+00  8.6133E-01  7.2672E-01  5.6053E-02  2.5572E+00  1.6627E+00
             2.8760E+00
 PARAMETER:  1.5755E-01  9.5216E-01 -5.0175E-01 -1.5274E+00  7.1747E-01 -4.9280E-02 -2.1922E-01 -2.7815E+00  1.0389E+00  6.0842E-01
             1.1564E+00
 GRADIENT:   2.7163E-02  6.9524E-02 -3.7430E-02  7.3513E-02 -1.0110E-02 -2.0048E-03 -1.0978E-02  2.1864E-03  1.5468E-02  8.6238E-04
             1.1260E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2663.04686210304        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1657
 NPARAMETR:  1.0592E+00  2.3476E+00  5.4552E-01  1.9444E-01  1.8564E+00  8.6132E-01  7.2641E-01  1.5019E-02  2.5693E+00  1.6643E+00
             2.8760E+00
 PARAMETER:  1.5753E-01  9.5340E-01 -5.0602E-01 -1.5376E+00  7.1861E-01 -4.9285E-02 -2.1964E-01 -4.0985E+00  1.0436E+00  6.0942E-01
             1.1564E+00
 GRADIENT:  -1.5682E-03 -2.1453E-03  9.2513E-04 -4.7675E-03  2.7387E-05 -1.6672E-03  6.7012E-03  1.5633E-04 -3.0228E-03 -5.1032E-04
             1.4953E-02

0ITERATION NO.:   64    OBJECTIVE VALUE:  -2663.04691179089        NO. OF FUNC. EVALS.: 133
 CUMULATIVE NO. OF FUNC. EVALS.:     1790
 NPARAMETR:  1.0593E+00  2.3471E+00  5.4594E-01  1.9462E-01  1.8561E+00  8.6140E-01  7.2657E-01  1.0000E-02  2.5676E+00  1.6641E+00
             2.8760E+00
 PARAMETER:  1.5753E-01  9.5325E-01 -5.0532E-01 -1.5365E+00  7.1849E-01 -4.9287E-02 -2.1964E-01 -4.8235E+00  1.0432E+00  6.0933E-01
             1.1564E+00
 GRADIENT:  -1.6954E-02  3.0174E-02 -1.3812E-04  1.8757E-03 -4.7986E-04 -3.4621E-03 -5.4936E-03  0.0000E+00  1.7332E-03  1.2593E-03
            -5.6007E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1790
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4420E-03 -2.6184E-02 -3.4472E-05  3.0414E-02 -2.2161E-02
 SE:             2.9195E-02  2.5199E-02  1.9041E-05  1.6620E-02  2.5934E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6061E-01  2.9875E-01  7.0234E-02  6.7259E-02  3.9283E-01

 ETASHRINKSD(%)  2.1947E+00  1.5581E+01  9.9936E+01  4.4320E+01  1.3117E+01
 ETASHRINKVR(%)  4.3412E+00  2.8734E+01  1.0000E+02  6.8997E+01  2.4514E+01
 EBVSHRINKSD(%)  2.3322E+00  1.3929E+01  9.9924E+01  5.3545E+01  1.0310E+01
 EBVSHRINKVR(%)  4.6100E+00  2.5917E+01  1.0000E+02  7.8419E+01  1.9557E+01
 RELATIVEINF(%)  9.5273E+01  1.2855E+01  3.4037E-05  3.4420E+00  4.4974E+01
 EPSSHRINKSD(%)  1.6026E+01
 EPSSHRINKVR(%)  2.9483E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          885
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1626.5212037722706     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2663.0469117908929     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1036.5257080186223     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    43.50
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2663.047       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  2.35E+00  5.46E-01  1.95E-01  1.86E+00  8.61E-01  7.26E-01  1.00E-02  2.57E+00  1.66E+00  2.88E+00
 


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
+        1.27E+03
 
 TH 2
+       -1.76E+01  2.92E+02
 
 TH 3
+        1.40E+00  1.00E+01  2.98E+01
 
 TH 4
+       -3.30E+01  3.75E+02 -6.97E+01  1.06E+03
 
 TH 5
+       -4.47E+00 -2.57E+01 -1.66E+01  7.78E+01  9.33E+01
 
 TH 6
+        8.40E+00 -4.32E+00  2.90E+00 -9.53E+00 -1.84E+00  2.37E+02
 
 TH 7
+        5.38E+00 -6.91E+00 -2.95E+00 -1.61E+01 -3.38E+00 -1.40E+00  2.29E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.01E-02 -2.51E+00 -1.90E+00  3.70E+01  5.03E-01 -5.08E-01  5.39E+00  0.00E+00  6.11E+00
 
 TH10
+       -3.65E-01 -4.27E+00 -4.51E+00  2.80E+01 -6.43E+00 -2.05E-01  4.04E+00  0.00E+00  1.69E-01  4.26E+01
 
 TH11
+       -2.00E+01 -1.40E+01 -5.09E-01 -8.97E+00  1.30E-01  4.89E+00  3.64E+00  0.00E+00  1.63E+00  5.03E+00  1.40E+02
 
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
 #CPUT: Total CPU Time in Seconds,       56.700
Stop Time:
Sat Sep 25 02:45:57 CDT 2021
