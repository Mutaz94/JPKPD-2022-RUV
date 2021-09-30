Wed Sep 29 09:55:16 CDT 2021
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
$DATA ../../../../data/int/D/dat86.csv ignore=@
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
 (2E4.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   42455.4255233879        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.6667E+02  8.0924E+02 -5.1658E+00  6.9075E+02  1.5973E+02 -3.1180E+03 -1.7250E+03 -3.2087E+01 -2.4632E+03 -8.1884E+02
            -8.3310E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -606.434532081803        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9469E-01  1.9514E+00  9.4787E-01  2.3240E+00  8.7855E-01  4.4076E+00  3.4796E+00  9.4838E-01  2.4666E+00  1.1231E+00
             1.2839E+01
 PARAMETER:  9.4680E-02  7.6855E-01  4.6457E-02  9.4328E-01 -2.9481E-02  1.5833E+00  1.3469E+00  4.7000E-02  1.0028E+00  2.1607E-01
             2.6525E+00
 GRADIENT:  -3.0872E+01  7.4654E+01 -2.8911E+01  1.7752E+02 -4.9952E+01  1.7280E+02 -4.0791E+01  5.3284E+00 -2.1787E+01  2.3013E+01
             1.9559E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -668.571539115274        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0338E+00  1.8971E+00  1.9968E+01  4.2545E+00  2.5084E+00  3.6689E+00  5.6136E+00  4.9411E-01  4.6747E+00  7.5178E-01
             1.2668E+01
 PARAMETER:  1.3325E-01  7.4033E-01  3.0941E+00  1.5480E+00  1.0196E+00  1.3999E+00  1.8252E+00 -6.0500E-01  1.6422E+00 -1.8531E-01
             2.6391E+00
 GRADIENT:  -2.5019E+01  3.7009E+01 -1.2851E+01  1.2925E+02  2.5862E+01  1.3848E+02  2.8980E+01  7.4050E-03  4.0297E+01  7.4977E+00
             2.5145E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -759.784350569307        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0529E+00  1.6050E+00  6.2890E+00  1.0510E+00  1.4419E+00  1.9480E+00  4.0589E+00  1.4779E+00  9.9365E-01  1.0537E+00
             1.2680E+01
 PARAMETER:  1.5152E-01  5.7313E-01  1.9388E+00  1.4976E-01  4.6598E-01  7.6680E-01  1.5009E+00  4.9065E-01  9.3632E-02  1.5233E-01
             2.6401E+00
 GRADIENT:  -6.6014E+01  2.1012E+01  1.7452E+01 -1.0900E+01 -1.2463E+02 -4.0183E+01 -1.1089E+00  6.2835E-01  4.8883E+00  1.4041E+01
             2.7538E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -815.901921296337        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.1482E+00  1.4211E+00  8.1955E+00  9.5553E-01  2.2915E+00  2.3406E+00  3.9545E+00  3.1215E-01  9.2141E-01  5.2730E-01
             1.1343E+01
 PARAMETER:  2.3820E-01  4.5144E-01  2.2036E+00  5.4507E-02  9.2921E-01  9.5042E-01  1.4748E+00 -1.0643E+00  1.8146E-02 -5.3999E-01
             2.5286E+00
 GRADIENT:  -2.3540E+00 -1.2013E+01 -3.3218E+00 -1.6070E+01  2.8454E+01  1.1836E+01  8.8852E+00  2.2202E-02  7.1347E+00  4.3238E+00
             1.6881E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -822.503456523302        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:      423
 NPARAMETR:  1.1476E+00  1.4199E+00  1.9866E+01  9.5642E-01  2.2928E+00  2.3310E+00  3.9971E+00  7.1133E-02  2.4928E-01  5.2474E-01
             1.0733E+01
 PARAMETER:  2.3770E-01  4.5060E-01  3.0890E+00  5.5438E-02  9.2977E-01  9.4629E-01  1.4856E+00 -2.5432E+00 -1.2892E+00 -5.4486E-01
             2.4733E+00
 GRADIENT:   7.2052E+00  1.6688E+00  3.9830E-01  3.8617E+01 -2.7711E+00  9.5891E+00 -8.8331E+00  2.6762E-04 -1.1035E+00  3.8721E+00
             4.3171E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -823.104823264856        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.1470E+00  1.4212E+00  1.8962E+01  9.5654E-01  2.2972E+00  2.3244E+00  4.0331E+00  4.1651E-02  3.7202E-01  5.2527E-01
             1.0480E+01
 PARAMETER:  2.3715E-01  4.5154E-01  3.0424E+00  5.5566E-02  9.3171E-01  9.4345E-01  1.4945E+00 -3.0784E+00 -8.8882E-01 -5.4384E-01
             2.4494E+00
 GRADIENT:   1.0740E+01  3.9581E+00  2.6943E-01  3.4144E+01 -5.2810E-01  7.9780E+00 -2.8995E+00  1.0065E-04 -1.3605E+00  3.8974E+00
             3.9844E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -826.065073214346        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      661
 NPARAMETR:  1.1440E+00  1.4290E+00  1.9598E+01  9.5759E-01  2.3212E+00  2.3092E+00  4.2776E+00  1.0746E-02  3.7278E-01  5.2840E-01
             1.0530E+01
 PARAMETER:  2.3454E-01  4.5697E-01  3.0754E+00  5.6660E-02  9.4207E-01  9.3689E-01  1.5534E+00 -4.4332E+00 -8.8676E-01 -5.3790E-01
             2.4543E+00
 GRADIENT:  -4.9230E+00 -7.4927E-01 -2.1453E-02  1.9024E+01  1.1075E+00 -2.3752E+01 -3.8283E+01  6.0654E-06 -7.6647E-01  3.8844E+00
            -2.4560E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -829.917063079306        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      799
 NPARAMETR:  1.1475E+00  1.4163E+00  2.1255E+01  9.6246E-01  2.3200E+00  2.3770E+00  4.8453E+00  1.0000E-02  3.5781E-01  4.5166E-01
             1.0573E+01
 PARAMETER:  2.3758E-01  4.4808E-01  3.1566E+00  6.1737E-02  9.4155E-01  9.6583E-01  1.6780E+00 -5.9449E+00 -9.2775E-01 -6.9482E-01
             2.4583E+00
 GRADIENT:   1.0833E+01  8.9821E+00 -8.1452E-02 -2.8352E+00  4.5438E+00  1.5505E+01  6.5521E+01  0.0000E+00  3.0447E-01  2.7566E+00
             4.0925E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -830.868900177715        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      987            RESET HESSIAN, TYPE II
 NPARAMETR:  1.1732E+00  1.3693E+00  2.4933E+01  9.8541E-01  2.3220E+00  2.3522E+00  5.0122E+00  1.0000E-02  3.3286E-01  1.7644E-01
             1.0787E+01
 PARAMETER:  2.5970E-01  4.1430E-01  3.3162E+00  8.5305E-02  9.4244E-01  9.5534E-01  1.7119E+00 -5.9449E+00 -1.0000E+00 -1.6348E+00
             2.4783E+00
 GRADIENT:   1.7225E+01  7.3597E+00  2.9980E-02 -5.0358E+00  2.4324E+00  1.1439E+01  7.2244E+01  0.0000E+00 -6.7581E-02  4.0782E-01
             7.1365E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -831.385212987225        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1125
 NPARAMETR:  1.1595E+00  1.3053E+00  2.4602E+01  1.0148E+00  2.3211E+00  2.3526E+00  4.9782E+00  1.0000E-02  4.9551E-01  1.2279E-01
             1.0655E+01
 PARAMETER:  2.4796E-01  3.6647E-01  3.3028E+00  1.1470E-01  9.4203E-01  9.5551E-01  1.7051E+00 -5.9449E+00 -6.0216E-01 -1.9972E+00
             2.4660E+00
 GRADIENT:  -5.4006E-01  7.1307E-01 -3.2283E-02 -4.6705E+00  7.5828E-01 -1.6435E+01 -1.9222E+00  0.0000E+00  4.7724E-01  1.9365E-01
             6.6370E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -832.182782590810        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1301
 NPARAMETR:  1.1627E+00  1.0672E+00  2.9943E+01  1.1464E+00  2.2883E+00  2.3897E+00  5.5106E+00  1.0000E-02  7.6579E-01  5.3076E-02
             1.0651E+01
 PARAMETER:  2.5076E-01  1.6507E-01  3.4993E+00  2.3661E-01  9.2780E-01  9.7117E-01  1.8067E+00 -5.9449E+00 -1.6684E-01 -2.8360E+00
             2.4657E+00
 GRADIENT:   6.6566E-01  2.7076E-01  2.6511E-01 -5.7232E+00 -3.8180E+00 -1.0262E+01  4.3845E+00  0.0000E+00  1.4066E+00  3.7634E-02
             6.0253E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -833.346850369061        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1478
 NPARAMETR:  1.1592E+00  7.1251E-01  8.8100E+00  1.4669E+00  2.0604E+00  2.4945E+00  6.1346E+00  1.0000E-02  1.3080E+00  1.0000E-02
             1.0635E+01
 PARAMETER:  2.4774E-01 -2.3895E-01  2.2759E+00  4.8315E-01  8.2290E-01  1.0141E+00  1.9139E+00 -5.9449E+00  3.6853E-01 -5.2588E+00
             2.4641E+00
 GRADIENT:  -1.7871E+00  1.7654E+00 -5.8498E+00  6.5474E+00  2.6918E+01  3.8902E+00  3.7829E-01  0.0000E+00 -2.7060E+00  0.0000E+00
             9.0551E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -835.034112537566        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1653
 NPARAMETR:  1.1592E+00  6.0656E-01  1.0225E+01  1.5341E+00  1.9655E+00  2.5226E+00  6.3903E+00  1.0000E-02  1.4519E+00  1.0000E-02
             1.0558E+01
 PARAMETER:  2.4769E-01 -3.9995E-01  2.4249E+00  5.2796E-01  7.7576E-01  1.0253E+00  1.9548E+00 -5.9449E+00  4.7288E-01 -5.7824E+00
             2.4569E+00
 GRADIENT:  -3.8254E-01  2.2278E+00  4.7972E-01  4.9320E+00  2.3416E+00  7.6390E+00 -2.0868E+00  0.0000E+00 -4.4851E-01  0.0000E+00
             4.5111E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -835.355638566703        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1837             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1605E+00  5.7800E-01  9.5300E+00  1.5425E+00  1.9235E+00  2.5286E+00  6.6040E+00  1.0000E-02  1.4525E+00  1.0000E-02
             1.0599E+01
 PARAMETER:  2.4881E-01 -4.4818E-01  2.3544E+00  5.3342E-01  7.5417E-01  1.0276E+00  1.9877E+00 -5.9449E+00  4.7332E-01 -5.8546E+00
             2.4608E+00
 GRADIENT:   1.4135E+01  5.6635E+00  5.5048E-01  1.9544E+01  4.0556E+00  4.1439E+01  1.1944E+02  0.0000E+00  2.2866E+00  0.0000E+00
             5.3443E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -835.574770864964        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2013
 NPARAMETR:  1.1605E+00  5.3421E-01  9.4861E+00  1.5329E+00  1.9185E+00  2.5100E+00  6.7501E+00  1.0000E-02  1.4586E+00  1.0000E-02
             1.0590E+01
 PARAMETER:  2.4889E-01 -5.2698E-01  2.3498E+00  5.2717E-01  7.5157E-01  1.0203E+00  2.0096E+00 -5.9449E+00  4.7751E-01 -5.8546E+00
             2.4599E+00
 GRADIENT:  -1.8812E-01  5.0684E-01 -3.3987E-01 -7.5831E+00  3.8952E+00  6.1998E+00  4.8093E+00  0.0000E+00  1.4512E+00  0.0000E+00
             1.3243E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -835.789546719648        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2196             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1592E+00  5.0768E-01  9.4245E+00  1.5585E+00  1.8916E+00  2.4834E+00  6.9033E+00  1.0000E-02  1.4514E+00  1.0000E-02
             1.0604E+01
 PARAMETER:  2.4773E-01 -5.7791E-01  2.3433E+00  5.4370E-01  7.3742E-01  1.0096E+00  2.0320E+00 -5.9449E+00  4.7250E-01 -5.8546E+00
             2.4613E+00
 GRADIENT:   1.3797E+01  5.1640E+00  7.8603E-01  1.5305E+01  2.0140E+00  3.4971E+01  1.3269E+02  0.0000E+00  2.5307E+00  0.0000E+00
             5.6591E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -835.986523486525        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2372
 NPARAMETR:  1.1605E+00  4.6694E-01  9.5019E+00  1.5791E+00  1.8854E+00  2.4752E+00  7.0728E+00  1.0000E-02  1.4577E+00  1.0000E-02
             1.0518E+01
 PARAMETER:  2.4884E-01 -6.6156E-01  2.3515E+00  5.5688E-01  7.3415E-01  1.0063E+00  2.0563E+00 -5.9449E+00  4.7683E-01 -5.8546E+00
             2.4531E+00
 GRADIENT:   6.0737E-01  6.8410E-01  2.8407E-01 -7.9650E-01 -1.3566E+00  1.4155E+00  5.6507E+00  0.0000E+00 -5.2727E-01  0.0000E+00
            -2.4169E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -836.067218540470        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2554
 NPARAMETR:  1.1647E+00  4.2936E-01  9.6080E+00  1.6060E+00  1.8879E+00  2.4757E+00  7.2289E+00  1.0000E-02  1.4657E+00  1.0000E-02
             1.0555E+01
 PARAMETER:  2.5248E-01 -7.4545E-01  2.3626E+00  5.7373E-01  7.3544E-01  1.0065E+00  2.0781E+00 -5.9449E+00  4.8232E-01 -5.8546E+00
             2.4566E+00
 GRADIENT:   1.2664E+00 -2.8883E-02 -2.0895E-01  8.9123E-01  4.1038E-02  1.6596E+00  5.5427E+00  0.0000E+00 -3.3901E-01  0.0000E+00
             2.9687E+00

0ITERATION NO.:   92    OBJECTIVE VALUE:  -836.067542446875        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     2619
 NPARAMETR:  1.1600E+00  4.2970E-01  9.6855E+00  1.6035E+00  1.8904E+00  2.4632E+00  7.2855E+00  1.0000E-02  1.4731E+00  1.0000E-02
             1.0588E+01
 PARAMETER:  2.5096E-01 -7.4221E-01  2.3694E+00  5.7497E-01  7.3597E-01  1.0052E+00  2.0771E+00 -5.9449E+00  4.8714E-01 -5.8546E+00
             2.4558E+00
 GRADIENT:   5.2370E-01  8.0800E-02 -4.9272E-02  9.6157E-01 -4.9540E-01  5.1710E+01 -1.7164E+01  0.0000E+00 -3.4130E-02  0.0000E+00
            -1.7076E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2619
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.5132E-03  5.7844E-02 -3.0676E-05 -8.2655E-02  7.2464E-05
 SE:             2.8645E-02  2.1730E-02  2.4360E-05  1.7373E-02  1.1588E-04
 N:                     100         100         100         100         100

 P VAL.:         7.6631E-01  7.7708E-03  2.0793E-01  1.9599E-06  5.3177E-01

 ETASHRINKSD(%)  4.0363E+00  2.7200E+01  9.9918E+01  4.1800E+01  9.9612E+01
 ETASHRINKVR(%)  7.9098E+00  4.7002E+01  1.0000E+02  6.6127E+01  9.9998E+01
 EBVSHRINKSD(%)  3.9623E+00  2.7410E+01  9.9906E+01  3.5918E+01  9.9524E+01
 EBVSHRINKVR(%)  7.7676E+00  4.7307E+01  1.0000E+02  5.8935E+01  9.9998E+01
 RELATIVEINF(%)  9.2032E+01  2.6691E+01  1.7107E-05  2.1582E+01  4.4793E-04
 EPSSHRINKSD(%)  5.6313E+00
 EPSSHRINKVR(%)  1.0945E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -836.06754244687477     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       818.02181732153599     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    81.97
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.85
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -836.068       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.16E+00  4.31E-01  9.67E+00  1.61E+00  1.89E+00  2.47E+00  7.22E+00  1.00E-02  1.47E+00  1.00E-02  1.05E+01
 


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
+        1.21E+02
 
 TH 2
+       -1.09E+00  6.97E+01
 
 TH 3
+        1.40E-01  5.03E-01  1.71E-01
 
 TH 4
+       -3.32E+00  2.97E+01 -1.80E-01  1.05E+02
 
 TH 5
+       -2.00E+00 -1.30E+01 -2.81E+00 -1.09E+01  6.75E+01
 
 TH 6
+       -6.54E+01  1.54E+00 -1.85E-02  6.76E-01  2.92E+01  1.99E+02
 
 TH 7
+        1.12E+01  6.38E+00 -4.27E-02 -4.39E+00 -3.80E+00  1.36E+01  8.67E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.74E-01 -4.07E+00 -1.87E-01 -1.67E+01  4.25E+00  7.36E+02  5.53E-01  0.00E+00  2.52E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        1.10E+00 -3.59E+00 -6.69E-02 -8.29E+00 -1.16E+00  9.02E+00  1.63E-01  0.00E+00  2.73E+00  0.00E+00  9.43E+00
 
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
 #CPUT: Total CPU Time in Seconds,       98.921
Stop Time:
Wed Sep 29 09:56:56 CDT 2021
