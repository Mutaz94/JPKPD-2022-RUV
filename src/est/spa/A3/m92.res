Sat Sep 25 09:35:48 CDT 2021
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
$DATA ../../../../data/spa/A3/dat92.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   153.555992874165        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.2078E+02 -1.2933E+01  7.3808E+01 -1.3170E+02  1.4688E+02  1.7054E+01 -4.2555E+01 -4.4974E+01 -1.2785E+02 -1.1811E+02
            -3.2401E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1190.92187642214        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0692E+00  1.0214E+00  9.4694E-01  1.1783E+00  9.3239E-01  8.7522E-01  9.9984E-01  1.0169E+00  1.1055E+00  1.0324E+00
             5.3660E+00
 PARAMETER:  1.6691E-01  1.2122E-01  4.5482E-02  2.6407E-01  3.0000E-02 -3.3277E-02  9.9845E-02  1.1680E-01  2.0027E-01  1.3192E-01
             1.7801E+00
 GRADIENT:   1.0159E+02 -7.5451E+00 -1.0773E+01  4.0534E-02 -1.2545E+01  7.6852E+00  1.0555E+01  6.5551E+00  2.2372E+01  2.1428E+01
             1.8694E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1221.42170242021        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0047E+00  8.5236E-01  4.0314E-01  1.1643E+00  5.0965E-01  8.8643E-01  2.5952E-01  3.1495E-01  1.1598E+00  5.3750E-01
             4.5372E+00
 PARAMETER:  1.0470E-01 -5.9747E-02 -8.0848E-01  2.5213E-01 -5.7402E-01 -2.0548E-02 -1.2489E+00 -1.0553E+00  2.4822E-01 -5.2083E-01
             1.6123E+00
 GRADIENT:  -2.3958E+01  3.2834E+01  1.5088E+00  4.2341E+01 -3.2912E+01  1.0604E+00 -6.4302E-02  1.1735E+00  1.4273E+01  8.7421E+00
             1.0269E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1241.85139288212        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.9311E-01  6.7954E-01  1.0457E+00  1.3296E+00  8.6658E-01  8.7144E-01  6.8237E-01  1.8635E-01  9.7628E-01  3.9454E-01
             3.9720E+00
 PARAMETER:  9.3090E-02 -2.8634E-01  1.4469E-01  3.8491E-01 -4.3204E-02 -3.7610E-02 -2.8219E-01 -1.5801E+00  7.5989E-02 -8.3003E-01
             1.4793E+00
 GRADIENT:   6.4984E+00  4.1955E+00 -2.5970E-01  1.1710E+01 -1.7353E+00  2.7972E+00  1.0544E+00  3.0050E-01  1.6772E+00  2.3190E+00
             8.4072E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1243.41422179179        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  9.8541E-01  5.1686E-01  1.0779E+00  1.4092E+00  8.3402E-01  8.6030E-01  2.4484E-01  4.0438E-02  9.2319E-01  1.1838E-01
             3.9847E+00
 PARAMETER:  8.5302E-02 -5.5999E-01  1.7502E-01  4.4304E-01 -8.1499E-02 -5.0473E-02 -1.3071E+00 -3.1080E+00  2.0084E-02 -2.0338E+00
             1.4825E+00
 GRADIENT:  -2.1431E+00 -8.7326E-01  2.5086E-02 -1.6833E+00  9.9217E-01  7.5682E-01  7.2478E-02  1.4828E-02  4.3907E-01  1.4638E-01
            -1.9438E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1243.50989038577        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.8519E-01  4.4114E-01  1.0374E+00  1.4520E+00  7.9006E-01  8.5872E-01  1.2276E-01  1.4804E-02  8.9652E-01  5.1497E-02
             3.9865E+00
 PARAMETER:  8.5075E-02 -7.1838E-01  1.3667E-01  4.7297E-01 -1.3565E-01 -5.2308E-02 -1.9975E+00 -4.1129E+00 -9.2293E-03 -2.8662E+00
             1.4829E+00
 GRADIENT:  -9.6630E-01  2.8543E-01  5.5827E-01  1.8326E+00 -9.6093E-01  6.3225E-02  1.1214E-02  2.1248E-03 -4.8472E-02  2.5908E-02
            -3.4847E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1243.51023748802        NO. OF FUNC. EVALS.:  91
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  9.8536E-01  4.4192E-01  1.0357E+00  1.4506E+00  7.8994E-01  8.5856E-01  1.2089E-01  1.4637E-02  8.9694E-01  5.0593E-02
             3.9872E+00
 PARAMETER:  8.5251E-02 -7.1662E-01  1.3504E-01  4.7196E-01 -1.3580E-01 -5.2503E-02 -2.0129E+00 -4.1242E+00 -8.7703E-03 -2.8839E+00
             1.4831E+00
 GRADIENT:  -2.6791E+00 -3.1087E-01  1.4672E-01 -2.3938E+00 -3.5446E-01 -1.4913E-01  1.0009E-02  2.0687E-03 -7.0695E-02  2.4833E-02
            -1.3120E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1243.52910442729        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      638
 NPARAMETR:  9.8563E-01  3.7760E-01  1.0272E+00  1.4894E+00  7.6819E-01  8.5890E-01  4.8359E-02  1.0000E-02  8.7259E-01  1.7654E-02
             3.9939E+00
 PARAMETER:  8.5530E-02 -8.7392E-01  1.2685E-01  4.9835E-01 -1.6372E-01 -5.2104E-02 -2.9291E+00 -5.3302E+00 -3.6293E-02 -3.9368E+00
             1.4848E+00
 GRADIENT:  -1.0674E-02 -8.1859E-02 -5.1129E-02 -2.9449E-01  9.2154E-02 -1.1811E-02  1.1547E-03  0.0000E+00  9.9424E-03  3.0394E-03
            -2.5684E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1243.52987678946        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      773
 NPARAMETR:  9.8548E-01  3.7286E-01  1.0264E+00  1.4922E+00  7.6624E-01  8.5889E-01  4.4027E-02  1.0000E-02  8.7080E-01  1.0000E-02
             3.9934E+00
 PARAMETER:  8.5376E-02 -8.8655E-01  1.2602E-01  5.0024E-01 -1.6626E-01 -5.2116E-02 -3.0230E+00 -5.4459E+00 -3.8345E-02 -4.8800E+00
             1.4846E+00
 GRADIENT:   2.0290E+00  2.3651E-01  1.3882E-01  3.3791E+00 -9.7247E-02  1.5973E-01  1.0715E-03  0.0000E+00  2.0949E-02  0.0000E+00
             1.0701E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1243.53047727075        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      956
 NPARAMETR:  9.8577E-01  3.8638E-01  1.0236E+00  1.4840E+00  7.6862E-01  8.5905E-01  3.8301E-02  1.0000E-02  8.7665E-01  1.0000E-02
             3.9936E+00
 PARAMETER:  8.5669E-02 -8.5093E-01  1.2332E-01  4.9474E-01 -1.6316E-01 -5.1933E-02 -3.1623E+00 -5.4459E+00 -3.1652E-02 -7.9176E+00
             1.4847E+00
 GRADIENT:  -1.0817E-01  1.8340E-03  1.0328E-02 -9.9857E-02 -2.0530E-01 -4.7900E-04  7.4264E-04  0.0000E+00  2.3150E-02  0.0000E+00
            -8.2388E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1243.53411598988        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1134
 NPARAMETR:  9.8656E-01  4.5136E-01  1.0533E+00  1.4498E+00  8.0136E-01  8.5882E-01  1.7078E-02  1.0000E-02  8.9842E-01  1.0000E-02
             3.9956E+00
 PARAMETER:  8.6467E-02 -6.9549E-01  1.5196E-01  4.7146E-01 -1.2144E-01 -5.2200E-02 -3.9700E+00 -5.4459E+00 -7.1151E-03 -2.0151E+01
             1.4852E+00
 GRADIENT:  -5.1457E-01  3.5093E-01  6.5762E-01  9.9095E-01 -1.0757E+00 -3.1763E-02  2.1091E-04  0.0000E+00 -4.7668E-02  0.0000E+00
            -2.2375E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1243.53620156367        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1310
 NPARAMETR:  9.8729E-01  4.8640E-01  1.0595E+00  1.4299E+00  8.1632E-01  8.5895E-01  1.1537E-02  1.0000E-02  9.1152E-01  1.0000E-02
             3.9977E+00
 PARAMETER:  8.7210E-02 -6.2072E-01  1.5783E-01  4.5761E-01 -1.0294E-01 -5.2046E-02 -4.3622E+00 -5.4459E+00  7.3546E-03 -2.8810E+01
             1.4857E+00
 GRADIENT:   1.3757E-01  6.7451E-02 -1.5746E-01  4.4794E-01  1.4666E-01 -6.0266E-03  1.1471E-04  0.0000E+00 -6.8916E-02  0.0000E+00
             5.4926E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1243.53704437972        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1486
 NPARAMETR:  9.8775E-01  5.3255E-01  1.0725E+00  1.4032E+00  8.3693E-01  8.5886E-01  1.0000E-02  1.0000E-02  9.3003E-01  1.0000E-02
             3.9983E+00
 PARAMETER:  8.7676E-02 -5.3009E-01  1.6996E-01  4.3877E-01 -7.8018E-02 -5.2145E-02 -4.8748E+00 -5.4459E+00  2.7464E-02 -4.1584E+01
             1.4859E+00
 GRADIENT:   2.8602E-03 -1.0025E-02 -9.9718E-03 -2.6050E-02  1.6235E-02 -8.3619E-04  0.0000E+00  0.0000E+00  1.9489E-03  0.0000E+00
             8.8794E-03

0ITERATION NO.:   63    OBJECTIVE VALUE:  -1243.53704694121        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1578
 NPARAMETR:  9.8777E-01  5.3401E-01  1.0729E+00  1.4024E+00  8.3760E-01  8.5887E-01  1.0000E-02  1.0000E-02  9.3061E-01  1.0000E-02
             3.9984E+00
 PARAMETER:  8.7694E-02 -5.2733E-01  1.7037E-01  4.3819E-01 -7.7215E-02 -5.2140E-02 -4.8907E+00 -5.4459E+00  2.8087E-02 -4.1982E+01
             1.4859E+00
 GRADIENT:   3.1260E-05  5.1574E-04  9.4050E-04  1.5668E-03 -1.4180E-03  4.3424E-05  0.0000E+00  0.0000E+00  2.5295E-04  0.0000E+00
            -7.8398E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1578
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6902E-03 -1.7983E-04  9.7214E-05 -1.2805E-02 -8.5037E-06
 SE:             2.7993E-02  7.9865E-05  9.6847E-05  2.3817E-02  1.6334E-04
 N:                     100         100         100         100         100

 P VAL.:         9.5185E-01  2.4341E-02  3.1548E-01  5.9082E-01  9.5848E-01

 ETASHRINKSD(%)  6.2208E+00  9.9732E+01  9.9676E+01  2.0211E+01  9.9453E+01
 ETASHRINKVR(%)  1.2055E+01  9.9999E+01  9.9999E+01  3.6337E+01  9.9997E+01
 EBVSHRINKSD(%)  6.0514E+00  9.9740E+01  9.9610E+01  1.9523E+01  9.9412E+01
 EBVSHRINKVR(%)  1.1737E+01  9.9999E+01  9.9998E+01  3.5235E+01  9.9997E+01
 RELATIVEINF(%)  8.3400E+01  2.2533E-05  1.1200E-04  3.6914E+00  1.3240E-04
 EPSSHRINKSD(%)  1.8991E+01
 EPSSHRINKVR(%)  3.4376E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1243.5370469412057     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -508.38622037746757     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.87
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1243.537       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.88E-01  5.34E-01  1.07E+00  1.40E+00  8.38E-01  8.59E-01  1.00E-02  1.00E-02  9.31E-01  1.00E-02  4.00E+00
 


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
+        1.36E+03
 
 TH 2
+       -9.54E+01  2.79E+02
 
 TH 3
+        9.28E+00  9.00E+01  1.26E+02
 
 TH 4
+       -1.11E+02  2.92E+02  2.42E+01  4.15E+02
 
 TH 5
+        2.02E+01 -2.53E+02 -2.47E+02 -1.23E+02  5.29E+02
 
 TH 6
+       -7.51E+00 -1.64E+01  9.37E+00 -2.69E+01 -1.72E+00  2.05E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.20E+01 -4.95E+01  3.88E+00 -9.64E+00  1.63E+01  1.10E+00  0.00E+00  0.00E+00  8.76E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.15E+01 -1.31E+01 -2.06E+00 -9.50E+00  6.90E+00  5.64E+00  0.00E+00  0.00E+00  1.16E+01  0.00E+00  2.79E+01
 
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
 #CPUT: Total CPU Time in Seconds,       24.024
Stop Time:
Sat Sep 25 09:36:14 CDT 2021
