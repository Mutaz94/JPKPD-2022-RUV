Sat Sep 18 13:26:09 CDT 2021
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
$DATA ../../../../data/spa/S2/dat40.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1702.41300558551        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3065E+01 -1.2049E+01 -3.7330E+01  3.1032E+01  4.3306E+01  5.4589E+01 -6.4263E+00  8.8568E+00 -8.3368E+00  1.2839E+01
             1.2359E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1708.05078392510        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0247E+00  1.0283E+00  1.0586E+00  9.6364E-01  1.0215E+00  8.3209E-01  1.0652E+00  9.2766E-01  1.0754E+00  8.9139E-01
             9.8918E-01
 PARAMETER:  1.2438E-01  1.2792E-01  1.5696E-01  6.2965E-02  1.2126E-01 -8.3815E-02  1.6319E-01  2.4907E-02  1.7266E-01 -1.4972E-02
             8.9117E-02
 GRADIENT:   1.2000E+02 -1.5715E+01 -6.5674E+00  3.1937E-01  3.5842E+01 -1.1310E+01 -9.1518E-01  1.4781E+00  6.8740E+00 -6.6919E+00
             2.3630E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1709.99026154360        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0204E+00  1.0361E+00  8.1161E-01  9.5853E-01  8.9086E-01  8.4051E-01  1.2233E+00  5.8324E-01  1.0412E+00  7.4940E-01
             9.8252E-01
 PARAMETER:  1.2015E-01  1.3550E-01 -1.0873E-01  5.7650E-02 -1.5568E-02 -7.3746E-02  3.0153E-01 -4.3915E-01  1.4034E-01 -1.8848E-01
             8.2367E-02
 GRADIENT:   9.8156E+01  7.7854E+00 -4.2281E+00  2.2289E+01  9.7973E+00 -7.1357E+00  4.9325E+00  1.6771E+00  1.2006E+01 -5.4061E+00
             2.2866E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1711.36206599263        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9797E-01  1.0791E+00  7.6047E-01  9.1777E-01  8.8623E-01  8.5149E-01  1.1578E+00  3.6287E-01  1.0126E+00  8.0987E-01
             9.7350E-01
 PARAMETER:  9.7963E-02  1.7613E-01 -1.7382E-01  1.4191E-02 -2.0774E-02 -6.0764E-02  2.4650E-01 -9.1372E-01  1.1251E-01 -1.1088E-01
             7.3138E-02
 GRADIENT:   1.9869E+01  4.6662E+00 -5.4698E-01  2.2123E+00 -3.5426E+00 -1.4047E+00  1.0399E+00  8.6555E-01  2.3553E+00  1.8197E+00
             8.5058E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1711.36509716914        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  9.9553E-01  1.0653E+00  7.5619E-01  9.2470E-01  8.7851E-01  8.5297E-01  1.1714E+00  3.1950E-01  9.9997E-01  8.0459E-01
             9.7276E-01
 PARAMETER:  9.5524E-02  1.6325E-01 -1.7946E-01  2.1710E-02 -2.9533E-02 -5.9036E-02  2.5816E-01 -1.0410E+00  9.9973E-02 -1.1742E-01
             7.2381E-02
 GRADIENT:   1.2561E+01  2.8346E+00 -7.3571E-01  1.4369E+00 -2.1728E+00 -8.5450E-01  7.2846E-01  6.6414E-01  1.4622E+00  1.4251E+00
             5.9369E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1711.36606077391        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  9.9382E-01  1.0552E+00  7.5175E-01  9.2960E-01  8.7220E-01  8.5405E-01  1.1818E+00  2.7392E-01  9.9083E-01  8.0020E-01
             9.7223E-01
 PARAMETER:  9.3801E-02  1.5371E-01 -1.8536E-01  2.6994E-02 -3.6741E-02 -5.7767E-02  2.6701E-01 -1.1949E+00  9.0784E-02 -1.2290E-01
             7.1832E-02
 GRADIENT:   7.3966E+00  1.5496E+00 -7.4166E-01  8.5938E-01 -1.2044E+00 -4.7702E-01  4.8647E-01  4.7878E-01  8.6167E-01  1.0458E+00
             3.9520E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1711.36709798094        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      441
 NPARAMETR:  9.9284E-01  1.0495E+00  7.4837E-01  9.3225E-01  8.6819E-01  8.5469E-01  1.1879E+00  2.3786E-01  9.8560E-01  7.9739E-01
             9.7190E-01
 PARAMETER:  9.2819E-02  1.4827E-01 -1.8985E-01  2.9850E-02 -4.1344E-02 -5.7021E-02  2.7216E-01 -1.3361E+00  8.5497E-02 -1.2641E-01
             7.1503E-02
 GRADIENT:   4.4434E+00  8.2651E-01 -6.8866E-01  5.1774E-01 -6.6309E-01 -2.6580E-01  3.3714E-01  3.5520E-01  5.1177E-01  7.8997E-01
             2.7325E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1711.36848588373        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      511
 NPARAMETR:  9.9202E-01  1.0448E+00  7.4490E-01  9.3430E-01  8.6457E-01  8.5524E-01  1.1930E+00  1.9766E-01  9.8130E-01  7.9485E-01
             9.7162E-01
 PARAMETER:  9.1993E-02  1.4383E-01 -1.9451E-01  3.2047E-02 -4.5520E-02 -5.6376E-02  2.7644E-01 -1.5212E+00  8.1123E-02 -1.2960E-01
             7.1214E-02
 GRADIENT:   1.9510E+00  2.2439E-01 -6.1367E-01  2.2337E-01 -2.1358E-01 -8.9396E-02  2.0553E-01  2.4133E-01  2.1944E-01  5.5306E-01
             1.6644E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1711.37222542073        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      581
 NPARAMETR:  9.9133E-01  1.0411E+00  7.4113E-01  9.3580E-01  8.6118E-01  8.5572E-01  1.1972E+00  1.4707E-01  9.7774E-01  7.9244E-01
             9.7137E-01
 PARAMETER:  9.1290E-02  1.4027E-01 -1.9959E-01  3.3649E-02 -4.9451E-02 -5.5812E-02  2.7999E-01 -1.8169E+00  7.7488E-02 -1.3264E-01
             7.0952E-02
 GRADIENT:  -1.8725E-01 -2.8130E-01 -5.1891E-01 -3.1856E-02  1.6588E-01  5.8498E-02  8.6027E-02  1.3173E-01 -2.6852E-02  3.2772E-01
             7.1169E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1711.38548142906        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      654
 NPARAMETR:  9.9087E-01  1.0390E+00  7.3739E-01  9.3642E-01  8.5845E-01  8.5605E-01  1.1999E+00  6.9306E-02  9.7551E-01  7.9031E-01
             9.7118E-01
 PARAMETER:  9.0825E-02  1.3829E-01 -2.0464E-01  3.4304E-02 -5.2622E-02 -5.5422E-02  2.8221E-01 -2.5692E+00  7.5208E-02 -1.3533E-01
             7.0753E-02
 GRADIENT:  -1.6286E+00 -5.4302E-01 -2.8770E-01 -2.2208E-01  3.3653E-01  1.5427E-01 -1.9908E-02  2.9066E-02 -2.0027E-01  8.2997E-02
            -1.7729E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1711.42902018819        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      729
 NPARAMETR:  9.9161E-01  1.0450E+00  7.3241E-01  9.3260E-01  8.5865E-01  8.5545E-01  1.1947E+00  1.1181E-02  9.7978E-01  7.9087E-01
             9.7133E-01
 PARAMETER:  9.1577E-02  1.4398E-01 -2.1142E-01  3.0219E-02 -5.2399E-02 -5.6125E-02  2.7787E-01 -4.3935E+00  7.9572E-02 -1.3462E-01
             7.0908E-02
 GRADIENT:   3.7007E-01 -5.0191E-01 -1.1820E+00  2.7966E-01  8.3418E-01 -8.9108E-02  1.2188E-01  8.5611E-04  2.3595E-01  5.2896E-01
             2.2291E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1711.86211564246        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      909
 NPARAMETR:  1.0061E+00  1.0446E+00  7.5082E-01  9.4003E-01  8.6824E-01  8.6316E-01  1.2086E+00  1.0000E-02  9.8012E-01  8.0539E-01
             9.7135E-01
 PARAMETER:  1.0604E-01  1.4359E-01 -1.8659E-01  3.8152E-02 -4.1286E-02 -4.7153E-02  2.8944E-01 -8.3191E+00  7.9915E-02 -1.1643E-01
             7.0936E-02
 GRADIENT:  -2.9451E+00 -1.3461E+00 -4.6656E-01 -6.5064E-01  4.0565E-01  5.0325E-02  2.4145E-01  0.0000E+00 -1.4895E-01  1.2793E-01
            -3.6241E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1711.90554026820        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1084
 NPARAMETR:  1.0075E+00  1.1438E+00  7.2887E-01  8.7988E-01  9.0330E-01  8.6332E-01  1.1153E+00  1.0000E-02  1.0383E+00  8.1692E-01
             9.7319E-01
 PARAMETER:  1.0748E-01  2.3438E-01 -2.1625E-01 -2.7974E-02 -1.6955E-03 -4.6964E-02  2.0915E-01 -1.2944E+01  1.3755E-01 -1.0221E-01
             7.2820E-02
 GRADIENT:   3.5501E-01 -1.8298E-02 -1.9983E-02  3.1165E-02  6.4315E-02  5.1098E-02  9.1290E-03  0.0000E+00  5.7402E-03  3.5473E-03
             2.0618E-02

0ITERATION NO.:   63    OBJECTIVE VALUE:  -1711.90556892093        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1176
 NPARAMETR:  1.0074E+00  1.1448E+00  7.2843E-01  8.7926E-01  9.0349E-01  8.6322E-01  1.1145E+00  1.0000E-02  1.0388E+00  8.1683E-01
             9.7314E-01
 PARAMETER:  1.0736E-01  2.3522E-01 -2.1687E-01 -2.8679E-02 -1.4909E-03 -4.7091E-02  2.0840E-01 -1.2953E+01  1.3805E-01 -1.0233E-01
             7.2775E-02
 GRADIENT:  -1.7268E-03  6.6939E-03  1.2734E-03  6.5676E-03 -2.6097E-03  6.0043E-04 -8.6401E-04  0.0000E+00 -1.3038E-03 -3.0426E-04
             3.9640E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1176
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.3738E-04 -9.0137E-03 -3.9797E-04  4.4044E-03 -1.8527E-02
 SE:             2.9810E-02  2.2521E-02  1.6020E-04  2.4221E-02  2.2303E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9632E-01  6.8898E-01  1.2986E-02  8.5570E-01  4.0615E-01

 ETASHRINKSD(%)  1.3154E-01  2.4553E+01  9.9463E+01  1.8858E+01  2.5282E+01
 ETASHRINKVR(%)  2.6291E-01  4.3077E+01  9.9997E+01  3.4160E+01  4.4171E+01
 EBVSHRINKSD(%)  5.3853E-01  2.4112E+01  9.9507E+01  1.9286E+01  2.4686E+01
 EBVSHRINKVR(%)  1.0742E+00  4.2410E+01  9.9998E+01  3.4853E+01  4.3278E+01
 RELATIVEINF(%)  9.8624E+01  3.6135E+00  2.6657E-04  4.6987E+00  6.1303E+00
 EPSSHRINKSD(%)  4.3990E+01
 EPSSHRINKVR(%)  6.8629E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1711.9055689209308     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -976.75474235719264     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.52
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1711.906       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.14E+00  7.28E-01  8.79E-01  9.03E-01  8.63E-01  1.11E+00  1.00E-02  1.04E+00  8.17E-01  9.73E-01
 


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
+        1.45E+03
 
 TH 2
+       -6.90E+00  4.11E+02
 
 TH 3
+        2.09E+01  2.14E+02  5.52E+02
 
 TH 4
+       -1.60E+01  3.05E+02 -3.02E+02  8.82E+02
 
 TH 5
+       -5.87E+00 -3.50E+02 -6.71E+02  3.57E+02  1.15E+03
 
 TH 6
+        1.38E-01 -2.00E+00  3.50E+00 -5.66E+00 -1.03E+00  2.66E+02
 
 TH 7
+        1.19E+00  2.77E+01 -8.18E+00 -1.31E+01 -3.19E+00 -1.44E-01  5.84E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.09E+00 -2.45E+01 -3.34E+01  4.10E+01 -1.22E+00  1.16E+00  1.57E+01  0.00E+00  8.83E+01
 
 TH10
+       -2.47E+00 -1.16E+01 -5.91E+01 -1.95E+01 -5.98E+01  7.83E-01  1.67E+01  0.00E+00  1.08E+01  1.02E+02
 
 TH11
+       -9.68E+00 -1.51E+01 -3.38E+01 -7.10E-01  5.01E+00  6.24E-01  6.14E+00  0.00E+00  8.61E+00  2.43E+01  2.27E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.332
Stop Time:
Sat Sep 18 13:26:27 CDT 2021
