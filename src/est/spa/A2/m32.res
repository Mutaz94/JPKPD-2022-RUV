Sat Sep 25 08:35:38 CDT 2021
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
$DATA ../../../../data/spa/A2/dat32.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m32.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -842.388162129039        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -2.5106E+01  1.2703E+01  3.2841E+01 -1.1087E+01  1.3693E+02  4.1841E+01 -2.6183E+01 -2.1726E+01 -5.0483E+01 -5.7085E+01
            -1.5390E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1408.27097196210        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0550E+00  9.8762E-01  9.8325E-01  1.0973E+00  8.9717E-01  7.5441E-01  9.4369E-01  9.6759E-01  9.7884E-01  8.4149E-01
             3.3657E+00
 PARAMETER:  1.5356E-01  8.7546E-02  8.3106E-02  1.9286E-01 -8.5074E-03 -1.8182E-01  4.2042E-02  6.7051E-02  7.8618E-02 -7.2584E-02
             1.3136E+00
 GRADIENT:  -1.2788E+01  3.9082E+01 -1.2827E+00  5.7331E+01 -2.4319E+01 -3.0957E+01  7.9523E+00  6.5273E+00  9.8271E+00  2.0768E+01
             7.9782E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1422.71144166642        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0606E+00  7.7177E-01  5.7304E-01  1.1976E+00  6.1265E-01  8.3889E-01  7.0575E-01  3.2057E-01  1.0411E+00  4.2535E-01
             3.1723E+00
 PARAMETER:  1.5884E-01 -1.5907E-01 -4.5680E-01  2.8034E-01 -3.8996E-01 -7.5672E-02 -2.4849E-01 -1.0377E+00  1.4025E-01 -7.5483E-01
             1.2544E+00
 GRADIENT:  -4.3851E-01  2.6804E+01 -2.3787E+01  9.3821E+01  1.4992E+01 -3.7113E+00  2.0582E+00  1.3353E+00  1.4146E+01  6.3367E+00
             5.4999E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1431.13852053084        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0495E+00  6.8002E-01  4.8808E-01  1.1741E+00  5.2020E-01  8.6312E-01  1.3226E+00  3.3410E-01  8.2312E-01  2.5002E-01
             2.7842E+00
 PARAMETER:  1.4835E-01 -2.8564E-01 -6.1727E-01  2.6047E-01 -5.5355E-01 -4.7197E-02  3.7960E-01 -9.9631E-01 -9.4648E-02 -1.2862E+00
             1.1240E+00
 GRADIENT:  -1.2319E+01  2.6541E+01 -5.1150E+00  5.8905E+01  6.7745E+00  2.7964E+00  3.4668E+00  1.1908E-01 -6.8617E+00  6.0903E-01
            -6.8300E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1435.36447000490        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0472E+00  5.4635E-01  2.5158E-01  1.0549E+00  3.3146E-01  8.7242E-01  1.1987E+00  5.5142E-01  8.8515E-01  6.9273E-02
             2.5792E+00
 PARAMETER:  1.4613E-01 -5.0450E-01 -1.2800E+00  1.5349E-01 -1.0042E+00 -3.6488E-02  2.8123E-01 -4.9526E-01 -2.1997E-02 -2.5697E+00
             1.0475E+00
 GRADIENT:  -7.1554E+00 -4.9717E+00 -1.0336E+01  6.8816E-01  2.1888E+01 -2.6295E+00 -1.7527E+00  5.6226E-01 -2.2967E+00  2.2258E-02
            -1.2007E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1436.22357660302        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      478
 NPARAMETR:  1.0543E+00  4.7441E-01  2.9880E-01  1.1278E+00  3.4162E-01  8.7359E-01  1.4475E+00  5.6527E-01  8.5435E-01  7.2225E-02
             2.6122E+00
 PARAMETER:  1.5285E-01 -6.4569E-01 -1.1080E+00  2.2029E-01 -9.7407E-01 -3.5149E-02  4.6984E-01 -4.7046E-01 -5.7417E-02 -2.5280E+00
             1.0602E+00
 GRADIENT:  -5.1704E+00  4.4880E-01 -1.4903E+00  1.6521E+00  8.3562E-01 -3.2931E-01 -9.7633E-01  5.6418E-01 -1.8192E-01 -5.1902E-02
            -3.8725E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1436.34237456792        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      653
 NPARAMETR:  1.0564E+00  4.2989E-01  3.2553E-01  1.1661E+00  3.4853E-01  8.7065E-01  1.6383E+00  5.6238E-01  8.3357E-01  7.0035E-02
             2.6349E+00
 PARAMETER:  1.5487E-01 -7.4423E-01 -1.0223E+00  2.5364E-01 -9.5402E-01 -3.8519E-02  5.9365E-01 -4.7558E-01 -8.2040E-02 -2.5588E+00
             1.0688E+00
 GRADIENT:   6.5011E-02  5.6440E-01  7.2060E-01  5.9632E-01 -1.1760E+00 -1.9237E-02 -3.4774E-02 -1.9776E-02 -9.7240E-02 -9.4357E-02
            -1.3267E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1436.72783674800        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      836
 NPARAMETR:  1.0549E+00  3.9506E-01  3.3274E-01  1.1833E+00  3.4444E-01  8.6990E-01  1.7400E+00  5.0594E-01  8.2457E-01  2.6540E-01
             2.6280E+00
 PARAMETER:  1.5343E-01 -8.2873E-01 -1.0004E+00  2.6833E-01 -9.6584E-01 -3.9382E-02  6.5388E-01 -5.8134E-01 -9.2893E-02 -1.2265E+00
             1.0662E+00
 GRADIENT:  -3.6080E+00  2.1142E+00  3.5485E+00 -6.0428E+00 -4.8280E+00  2.4181E-01  3.9757E-01  1.4346E+00 -2.0277E+00  2.0311E-01
             3.3875E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1437.16079497188        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1014
 NPARAMETR:  1.0541E+00  3.0447E-01  3.5666E-01  1.2417E+00  3.4263E-01  8.6415E-01  2.1073E+00  3.1455E-01  8.2315E-01  4.0024E-01
             2.6262E+00
 PARAMETER:  1.5264E-01 -1.0892E+00 -9.3096E-01  3.1650E-01 -9.7112E-01 -4.6004E-02  8.4542E-01 -1.0566E+00 -9.4617E-02 -8.1569E-01
             1.0655E+00
 GRADIENT:   3.2358E-01  1.1769E+00  6.1000E+00  2.7463E+00 -9.1942E+00  3.4750E-02  2.0498E-01  8.3307E-02 -1.1608E-02 -7.6826E-01
            -1.0959E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1437.24861249761        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1190
 NPARAMETR:  1.0537E+00  2.9505E-01  3.5540E-01  1.2445E+00  3.4126E-01  8.6333E-01  2.1237E+00  1.7748E-01  8.2626E-01  4.4927E-01
             2.6307E+00
 PARAMETER:  1.5229E-01 -1.1206E+00 -9.3450E-01  3.1870E-01 -9.7512E-01 -4.6962E-02  8.5317E-01 -1.6289E+00 -9.0851E-02 -7.0013E-01
             1.0672E+00
 GRADIENT:   7.3858E-02 -2.5283E-01 -1.4624E+00  1.6254E+00  1.6951E+00  4.1963E-02  1.8592E-01  6.6739E-02  1.0846E-01  2.1515E-01
             6.8967E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1437.27229239168        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1366
 NPARAMETR:  1.0542E+00  3.0865E-01  3.5113E-01  1.2355E+00  3.4064E-01  8.6374E-01  2.0406E+00  7.9109E-02  8.2974E-01  4.5235E-01
             2.6293E+00
 PARAMETER:  1.5276E-01 -1.0755E+00 -9.4660E-01  3.1149E-01 -9.7694E-01 -4.6480E-02  8.1324E-01 -2.4369E+00 -8.6645E-02 -6.9330E-01
             1.0667E+00
 GRADIENT:   5.5758E-02  1.7768E-01  8.3644E-02  2.2781E-01 -2.4421E-01  2.3756E-03  3.8442E-02  5.0761E-03  1.8924E-02  2.8551E-02
             7.5215E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1437.27513411683        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1541
 NPARAMETR:  1.0540E+00  3.0450E-01  3.5197E-01  1.2376E+00  3.4050E-01  8.6347E-01  2.0581E+00  1.4935E-02  8.2936E-01  4.5734E-01
             2.6298E+00
 PARAMETER:  1.5259E-01 -1.0891E+00 -9.4421E-01  3.1318E-01 -9.7733E-01 -4.6799E-02  8.2179E-01 -4.1040E+00 -8.7101E-02 -6.8234E-01
             1.0669E+00
 GRADIENT:  -3.6642E-03  8.6615E-03  2.8430E-02 -2.2937E-02 -3.7062E-02 -1.8713E-03  4.1424E-03  1.6510E-04  3.7056E-03  6.4364E-03
             7.9601E-03

0ITERATION NO.:   58    OBJECTIVE VALUE:  -1437.27517950232        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1633
 NPARAMETR:  1.0540E+00  3.0446E-01  3.5197E-01  1.2376E+00  3.4050E-01  8.6347E-01  2.0582E+00  1.0000E-02  8.2935E-01  4.5733E-01
             2.6298E+00
 PARAMETER:  1.5259E-01 -1.0892E+00 -9.4421E-01  3.1320E-01 -9.7734E-01 -4.6795E-02  8.2185E-01 -4.6360E+00 -8.7117E-02 -6.8235E-01
             1.0669E+00
 GRADIENT:  -1.0499E-03 -1.2808E-03 -7.2918E-03  4.8700E-03  9.2403E-03  4.8710E-04  5.2009E-04  0.0000E+00  7.4220E-04 -1.5609E-04
            -2.8888E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1633
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.5088E-04  3.7324E-02 -2.8137E-04 -2.1914E-02  1.2365E-02
 SE:             2.8917E-02  1.7297E-02  2.3134E-04  2.5141E-02  1.5186E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7377E-01  3.0937E-02  2.2388E-01  3.8341E-01  4.1552E-01

 ETASHRINKSD(%)  3.1235E+00  4.2054E+01  9.9225E+01  1.5775E+01  4.9124E+01
 ETASHRINKVR(%)  6.1495E+00  6.6423E+01  9.9994E+01  2.9062E+01  7.4116E+01
 EBVSHRINKSD(%)  3.1660E+00  4.9598E+01  9.9178E+01  1.4397E+01  4.6126E+01
 EBVSHRINKVR(%)  6.2318E+00  7.4597E+01  9.9993E+01  2.6722E+01  7.0976E+01
 RELATIVEINF(%)  9.2185E+01  5.8087E+00  1.7870E-04  1.8498E+01  7.1789E-01
 EPSSHRINKSD(%)  3.2073E+01
 EPSSHRINKVR(%)  5.3859E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1437.2751795023173     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -702.12435293857914     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.73
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1437.275       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  3.04E-01  3.52E-01  1.24E+00  3.41E-01  8.63E-01  2.06E+00  1.00E-02  8.29E-01  4.57E-01  2.63E+00
 


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
+        1.26E+03
 
 TH 2
+       -8.17E+01  7.52E+02
 
 TH 3
+        5.39E+00  1.21E+03  6.50E+03
 
 TH 4
+       -6.34E+01  2.69E+02 -6.40E+02  7.80E+02
 
 TH 5
+        1.59E+02 -2.31E+03 -8.90E+03  5.92E+01  1.37E+04
 
 TH 6
+       -2.22E+00 -1.11E+01  3.06E+01 -1.66E+01  2.31E+00  2.40E+02
 
 TH 7
+        1.37E+00  4.30E+01  3.68E-01 -1.33E+00 -2.53E+01  4.71E-01  8.28E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.77E+00 -6.04E+00  3.03E+00 -2.21E+01  8.59E+01  6.52E-01  7.01E+00  0.00E+00  1.63E+02
 
 TH10
+       -8.35E+00  2.47E+01 -2.04E+02 -1.68E+01  2.64E+02 -3.55E-01  5.15E+00  0.00E+00 -7.74E+00  1.03E+02
 
 TH11
+       -1.55E+01  5.60E+00 -6.37E+01 -1.15E+01  5.68E+01  4.83E+00  2.69E+00  0.00E+00  9.28E+00  2.36E+01  4.23E+01
 
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
 #CPUT: Total CPU Time in Seconds,       26.319
Stop Time:
Sat Sep 25 08:36:06 CDT 2021
