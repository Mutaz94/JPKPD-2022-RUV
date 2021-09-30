Wed Sep 29 06:17:06 CDT 2021
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
$DATA ../../../../data/int/TD1/dat38.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
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

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3456.46535351202        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6332E+02  6.0604E+01  1.4177E+02  6.5434E+01  7.5788E+01  1.9580E+01 -4.9796E+01 -4.4593E+02 -1.2899E+02 -1.0860E+01
            -1.1568E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3611.42455995123        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:      148
 NPARAMETR:  9.6890E-01  1.0896E+00  1.0446E+00  9.2608E-01  1.0725E+00  6.8976E-01  1.1218E+00  1.8726E+00  1.0130E+00  9.4880E-01
             1.1215E+00
 PARAMETER:  6.8405E-02  1.8583E-01  1.4362E-01  2.3200E-02  1.7003E-01 -2.7141E-01  2.1491E-01  7.2735E-01  1.1288E-01  4.7443E-02
             2.1469E-01
 GRADIENT:  -1.1749E+02 -2.9337E+01  1.5368E+01 -6.2511E+01 -5.6942E+01 -2.7205E+02 -1.9502E+01 -1.3579E+02 -2.8060E+01 -9.8784E-01
             1.5034E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3629.28626204903        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      326
 NPARAMETR:  9.4942E-01  1.1212E+00  1.2801E+00  9.4027E-01  1.2359E+00  6.9701E-01  1.2575E+00  2.2406E+00  8.4482E-01  1.1314E+00
             1.1228E+00
 PARAMETER:  4.8098E-02  2.1442E-01  3.4691E-01  3.8410E-02  3.1179E-01 -2.6096E-01  3.2909E-01  9.0674E-01 -6.8634E-02  2.2345E-01
             2.1583E-01
 GRADIENT:  -2.0786E+02 -4.4138E+01  1.2028E+01  2.0906E+01  2.2602E+01 -2.6957E+02 -9.9011E+00 -1.1266E+02 -3.0294E+01  1.7114E+01
             1.5341E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3648.12189287081        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:      476
 NPARAMETR:  9.2515E-01  1.1648E+00  1.1999E+00  8.9746E-01  1.2567E+00  7.7782E-01  1.1671E+00  2.3260E+00  7.6512E-01  1.1892E+00
             1.1227E+00
 PARAMETER:  2.2200E-02  2.5254E-01  2.8223E-01 -8.1858E-03  3.2851E-01 -1.5125E-01  2.5451E-01  9.4416E-01 -1.6772E-01  2.7326E-01
             2.1576E-01
 GRADIENT:   1.0953E+02  1.5365E+02  1.1042E+01  5.0828E+01  1.2210E+02 -1.0538E+02 -4.3584E+00 -6.0807E+01 -3.8771E+01  2.8148E+01
             1.5657E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3677.54198421728        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      595
 NPARAMETR:  9.4042E-01  1.1642E+00  1.2005E+00  8.8650E-01  1.2574E+00  9.7048E-01  1.2088E+00  2.3302E+00  7.6491E-01  1.1886E+00
             1.1231E+00
 PARAMETER:  3.8572E-02  2.5206E-01  2.8272E-01 -2.0471E-02  3.2902E-01  7.0036E-02  2.8964E-01  9.4597E-01 -1.6800E-01  2.7278E-01
             2.1609E-01
 GRADIENT:  -1.2795E+02 -6.5359E+01  2.5112E+00 -3.7492E+01  2.0871E+00 -4.8839E+01 -2.8036E+01 -9.9273E+01 -4.1279E+01  1.4886E+01
             1.5501E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3680.00078238726        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:      786
 NPARAMETR:  9.5031E-01  1.1643E+00  1.2004E+00  8.9091E-01  1.2573E+00  9.9490E-01  1.2108E+00  2.3303E+00  7.6493E-01  1.1887E+00
             1.1231E+00
 PARAMETER:  4.9031E-02  2.5211E-01  2.8266E-01 -1.5509E-02  3.2896E-01  9.4890E-02  2.9127E-01  9.4601E-01 -1.6797E-01  2.7283E-01
             2.1605E-01
 GRADIENT:  -9.8337E+01 -6.1126E+01  1.8281E+00 -2.8528E+01  4.3285E+00 -3.4644E+01 -2.7394E+01 -9.9259E+01 -4.1186E+01  1.5294E+01
             1.5522E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3681.30190211505        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      904
 NPARAMETR:  9.4979E-01  1.1643E+00  1.1986E+00  8.9043E-01  1.2550E+00  1.0788E+00  1.2089E+00  2.3307E+00  7.6493E-01  1.1887E+00
             1.1230E+00
 PARAMETER:  4.8487E-02  2.5211E-01  2.8112E-01 -1.6053E-02  3.2716E-01  1.7590E-01  2.8968E-01  9.4615E-01 -1.6797E-01  2.7283E-01
             2.1604E-01
 GRADIENT:  -8.4309E+01 -6.1552E+01  1.9018E+00 -3.0516E+01  2.2258E+00  2.6125E-01 -2.7845E+01 -9.8997E+01 -4.1464E+01  1.5402E+01
             1.5565E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3687.11353165005        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  9.7254E-01  1.1638E+00  1.1990E+00  8.9867E-01  1.2557E+00  1.0783E+00  1.2180E+00  2.3773E+00  7.8648E-01  1.1813E+00
             1.1165E+00
 PARAMETER:  7.2153E-02  2.5171E-01  2.8148E-01 -6.8362E-03  3.2768E-01  1.7534E-01  2.9721E-01  9.6595E-01 -1.4018E-01  2.6658E-01
             2.1022E-01
 GRADIENT:  -3.9327E+01 -5.3344E+01 -8.8178E-01 -1.5167E+01  7.5646E+00  2.9393E+00 -2.2958E+01 -9.1249E+01 -3.7426E+01  1.5658E+01
             1.4706E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3688.83727986172        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1258
 NPARAMETR:  9.9418E-01  1.1646E+00  1.1986E+00  9.0816E-01  1.2550E+00  1.1414E+00  1.2182E+00  2.4273E+00  7.8669E-01  1.1813E+00
             1.1153E+00
 PARAMETER:  9.4162E-02  2.5241E-01  2.8116E-01  3.6693E-03  3.2712E-01  2.3227E-01  2.9738E-01  9.8680E-01 -1.3992E-01  2.6660E-01
             2.0912E-01
 GRADIENT:   2.4922E+00 -4.4515E+01 -3.4468E+00  3.7365E+00  1.0019E+01  2.4984E+01 -2.2750E+01 -8.4188E+01 -3.6961E+01  1.6409E+01
             1.4641E+02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3690.19408191136        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1439
 NPARAMETR:  9.9456E-01  1.1656E+00  1.1985E+00  9.0771E-01  1.2548E+00  1.0594E+00  1.2174E+00  2.4391E+00  7.8703E-01  1.1820E+00
             1.1142E+00
 PARAMETER:  9.4546E-02  2.5328E-01  2.8110E-01  3.1736E-03  3.2699E-01  1.5770E-01  2.9673E-01  9.9162E-01 -1.3948E-01  2.6723E-01
             2.0816E-01
 GRADIENT:   3.2533E+00 -4.3891E+01 -3.6555E+00  3.5778E+00  9.0905E+00 -3.2836E+00 -2.2795E+01 -8.2524E+01 -3.6752E+01  1.6417E+01
             1.4448E+02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3695.41824425812        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1576
 NPARAMETR:  9.9549E-01  1.1657E+00  1.1985E+00  9.0505E-01  1.2548E+00  1.0635E+00  1.2546E+00  2.4524E+00  1.0068E+00  1.1804E+00
             1.1154E+00
 PARAMETER:  9.5479E-02  2.5328E-01  2.8110E-01  2.3431E-04  3.2699E-01  1.6157E-01  3.2685E-01  9.9705E-01  1.0680E-01  2.6581E-01
             2.0925E-01
 GRADIENT:   5.3326E+00 -4.5386E+01 -3.2885E+00 -1.6104E+00  1.6061E+01 -1.8330E+00  4.3137E+00 -7.1143E+01  1.4033E-01  2.2217E+01
             1.5107E+02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3697.00213613807        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1756
 NPARAMETR:  9.8899E-01  1.1665E+00  1.1986E+00  9.0990E-01  1.2542E+00  1.0694E+00  1.0930E+00  2.5364E+00  1.0283E+00  1.1790E+00
             1.1124E+00
 PARAMETER:  8.8926E-02  2.5399E-01  2.8116E-01  5.5763E-03  3.2654E-01  1.6707E-01  1.8897E-01  1.0307E+00  1.2789E-01  2.6463E-01
             2.0649E-01
 GRADIENT:  -7.3396E+00 -5.6001E+01 -5.2487E+00  2.1544E+00  1.9737E+01  3.6449E-01 -1.9401E+01 -5.9382E+01 -4.8383E+00  2.2353E+01
             1.4515E+02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3698.72772767778        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1955
 NPARAMETR:  9.8892E-01  1.1669E+00  1.2003E+00  9.1032E-01  1.2539E+00  1.0693E+00  1.0846E+00  2.5875E+00  1.0311E+00  1.1714E+00
             1.1077E+00
 PARAMETER:  8.8860E-02  2.5439E-01  2.8261E-01  6.0444E-03  3.2628E-01  1.6697E-01  1.8121E-01  1.0507E+00  1.3062E-01  2.5821E-01
             2.0230E-01
 GRADIENT:  -7.3887E+00 -5.5384E+01 -6.0151E+00  2.1827E+00  1.8666E+01  2.8965E-01 -2.0625E+01 -5.3519E+01 -4.3658E+00  2.1192E+01
             1.3771E+02

0ITERATION NO.:   64    OBJECTIVE VALUE:  -3698.85998585425        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:     2077
 NPARAMETR:  9.8904E-01  1.1666E+00  1.2016E+00  9.1022E-01  1.2541E+00  1.0667E+00  1.0923E+00  2.5844E+00  1.0315E+00  1.1718E+00
             1.1080E+00
 PARAMETER:  8.8875E-02  2.5435E-01  2.8334E-01  6.0292E-03  3.2628E-01  1.6292E-01  1.8844E-01  1.0505E+00  1.3111E-01  2.5825E-01
             2.0233E-01
 GRADIENT:  -1.2135E+04  9.4749E+03 -4.2817E+03  1.2125E+04 -1.0524E+04 -1.5956E+00  1.2852E+04  2.2129E+03  1.8495E+04 -9.3736E+03
            -1.1850E+04
 NUMSIGDIG:         2.3         2.3         2.3         2.3         2.7         1.3         2.3         2.3         2.3         2.3
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2077
 NO. OF SIG. DIGITS IN FINAL EST.:  1.3
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.6858E-03 -1.1857E-02 -3.0802E-02  2.6046E-02 -6.4604E-02
 SE:             3.0034E-02  2.4840E-02  2.7916E-02  2.6155E-02  2.1330E-02
 N:                     100         100         100         100         100

 P VAL.:         8.7602E-01  6.3313E-01  2.6986E-01  3.1933E-01  2.4555E-03

 ETASHRINKSD(%)  1.0000E-10  1.6784E+01  6.4778E+00  1.2379E+01  2.8542E+01
 ETASHRINKVR(%)  1.0000E-10  3.0750E+01  1.2536E+01  2.3226E+01  4.8937E+01
 EBVSHRINKSD(%)  2.8502E-01  2.2620E+01  2.2308E+01  1.6288E+01  1.9740E+01
 EBVSHRINKVR(%)  5.6923E-01  4.0123E+01  3.9639E+01  2.9923E+01  3.5583E+01
 RELATIVEINF(%)  9.9429E+01  2.9841E+01  5.2961E+01  3.9037E+01  3.6513E+01
 EPSSHRINKSD(%)  2.7513E+01
 EPSSHRINKVR(%)  4.7457E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3698.8599858542475     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2044.7706260858367     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    66.02
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.95
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3698.860       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.89E-01  1.17E+00  1.20E+00  9.10E-01  1.25E+00  1.06E+00  1.09E+00  2.59E+00  1.03E+00  1.17E+00  1.11E+00
 


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
+        6.20E+06
 
 TH 2
+        5.84E+01  1.66E+06
 
 TH 3
+       -3.43E+03  1.16E+03  5.23E+05
 
 TH 4
+        1.95E+02 -1.27E+03  3.69E+03  7.31E+06
 
 TH 5
+       -2.12E+06 -5.00E+05 -4.78E+01 -1.63E+06  1.24E+06
 
 TH 6
+       -7.38E+01  2.70E+01 -2.39E+01  8.92E+01 -3.65E+01  1.75E+02
 
 TH 7
+       -4.26E+02  1.86E+02 -1.24E+02  4.81E+02  4.64E+02  3.96E+01  1.43E+06
 
 TH 8
+        6.68E+00 -1.03E+02 -9.13E+04 -1.73E+02 -6.44E+01  2.87E+00  1.50E+01  7.93E+03
 
 TH 9
+        3.17E+01 -2.14E+06  1.86E+06 -1.80E+01  1.10E+01  6.06E+01  1.91E+01  5.64E+00  3.31E+06
 
 TH10
+        1.39E+03 -9.57E+05  8.33E+05 -1.50E+03 -2.63E+01 -2.70E+01 -1.36E+02 -5.00E+01  1.23E+01  6.63E+05
 
 TH11
+       -8.69E+01 -1.93E+02 -1.52E+03  2.11E+03  5.99E+02 -3.44E+01 -1.79E+02 -1.32E+01 -2.83E+06 -1.27E+06  1.21E+06
 
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
 #CPUT: Total CPU Time in Seconds,       81.077
Stop Time:
Wed Sep 29 06:18:29 CDT 2021
