Sun Oct 24 02:53:51 CDT 2021
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
$DATA ../../../../data/SD4/SL1/dat1.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m1.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1666.05957186984        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3207E+02  2.5353E+01  3.4963E+01  2.1758E+01 -4.2823E+01  4.9700E+01  2.6839E+00 -1.0614E+00  1.2661E+01 -2.2872E+01
             1.0445E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1669.71976665788        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  9.9876E-01  1.0244E+00  9.5673E-01  9.9758E-01  1.0598E+00  9.8200E-01  1.0067E+00  9.8668E-01  9.2908E-01  1.2624E+00
             9.2762E-01
 PARAMETER:  9.8755E-02  1.2409E-01  5.5762E-02  9.7581E-02  1.5809E-01  8.1840E-02  1.0665E-01  8.6590E-02  2.6435E-02  3.3303E-01
             2.4871E-02
 GRADIENT:  -1.5522E+00 -2.4844E+01 -1.7238E+01 -2.5750E+01 -3.6643E-01 -4.3930E+00  4.5754E-01  6.0267E+00 -8.7977E+00  7.3963E+00
            -1.6281E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1671.48081140529        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      352
 NPARAMETR:  9.9823E-01  1.2141E+00  1.0342E+00  8.9556E-01  1.2196E+00  9.9337E-01  8.0220E-01  1.0336E+00  1.1080E+00  1.4894E+00
             9.5252E-01
 PARAMETER:  9.8233E-02  2.9396E-01  1.3366E-01 -1.0305E-02  2.9852E-01  9.3351E-02 -1.2039E-01  1.3308E-01  2.0259E-01  4.9835E-01
             5.1352E-02
 GRADIENT:  -2.9353E+00 -4.5710E+00  1.3961E+00 -2.5284E+00 -4.2479E+00  5.1354E-01  1.7176E+00 -3.5320E+00 -1.2593E-01  1.5152E+01
            -4.1511E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1672.75971542675        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.9980E-01  1.2876E+00  1.1007E+00  8.5882E-01  1.2751E+00  9.9252E-01  7.2847E-01  1.3722E+00  1.1743E+00  1.3928E+00
             9.6197E-01
 PARAMETER:  9.9797E-02  3.5281E-01  1.9594E-01 -5.2194E-02  3.4303E-01  9.2490E-02 -2.1681E-01  4.1638E-01  2.6063E-01  4.3130E-01
             6.1231E-02
 GRADIENT:   1.6665E-01  7.1365E+00  1.6864E+00  1.0377E+01 -1.5990E+00  6.8082E-02  6.3790E-02 -4.3817E-01  7.7626E-01 -5.2008E-01
            -3.3255E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1673.34693607881        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  1.0012E+00  1.5949E+00  8.7725E-01  6.5445E-01  1.3395E+00  9.9612E-01  6.8180E-01  1.3932E+00  1.4128E+00  1.4074E+00
             9.6186E-01
 PARAMETER:  1.0120E-01  5.6683E-01 -3.0960E-02 -3.2396E-01  3.9227E-01  9.6113E-02 -2.8301E-01  4.3162E-01  4.4560E-01  4.4172E-01
             6.1112E-02
 GRADIENT:   4.4625E-02  2.1497E+01  9.5366E+00  5.5529E+00 -1.2335E+01  8.1599E-01  2.0106E-01 -2.1438E+00 -1.6122E+00 -6.5876E-01
            -8.0509E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1673.80439379860        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  1.0018E+00  1.8153E+00  6.4966E-01  5.1007E-01  1.3844E+00  9.9541E-01  6.6465E-01  1.2445E+00  1.6703E+00  1.4346E+00
             9.5872E-01
 PARAMETER:  1.0180E-01  6.9625E-01 -3.3131E-01 -5.7322E-01  4.2528E-01  9.5403E-02 -3.0849E-01  3.1876E-01  6.1301E-01  4.6089E-01
             5.7841E-02
 GRADIENT:  -5.4861E-01  2.9953E+01  7.1797E+00  8.5746E+00 -1.2526E+01  4.6277E-02  9.4570E-02 -1.2605E+00 -1.1714E+00  2.4518E+00
            -2.0366E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1674.10026621173        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1065             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0025E+00  1.8317E+00  5.7511E-01  4.6695E-01  1.3927E+00  9.9627E-01  6.6335E-01  1.1496E+00  1.7653E+00  1.4300E+00
             9.5839E-01
 PARAMETER:  1.0252E-01  7.0522E-01 -4.5319E-01 -6.6153E-01  4.3124E-01  9.6266E-02 -3.1046E-01  2.3945E-01  6.6832E-01  4.5767E-01
             5.7501E-02
 GRADIENT:   4.9592E+02  8.8636E+02  1.3440E+00  1.0807E+02  2.1621E+01  5.4023E+01  1.9342E+01  8.5316E-01  3.1617E+01  1.0485E+01
             2.1421E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1674.33858735170        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1246
 NPARAMETR:  1.0014E+00  1.8567E+00  5.6856E-01  4.6629E-01  1.3997E+00  9.9573E-01  6.5840E-01  1.1335E+00  1.7536E+00  1.4081E+00
             9.6172E-01
 PARAMETER:  1.0143E-01  7.1878E-01 -4.6464E-01 -6.6296E-01  4.3623E-01  9.5722E-02 -3.1795E-01  2.2529E-01  6.6168E-01  4.4226E-01
             6.0971E-02
 GRADIENT:  -1.8145E+00 -9.9152E+00  1.1563E+00  1.0720E+00  1.6074E-01  2.5052E-02 -3.8596E-01 -5.9817E-02 -2.9217E-01  8.1999E-02
            -3.1683E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1674.35211409771        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1425
 NPARAMETR:  1.0009E+00  1.8591E+00  5.6057E-01  4.6356E-01  1.3975E+00  9.9504E-01  6.6016E-01  1.1120E+00  1.7508E+00  1.4052E+00
             9.6221E-01
 PARAMETER:  1.0094E-01  7.2009E-01 -4.7881E-01 -6.6881E-01  4.3468E-01  9.5030E-02 -3.1527E-01  2.0616E-01  6.6006E-01  4.4020E-01
             6.1472E-02
 GRADIENT:  -2.9827E+00 -1.2072E+01  8.5892E-01  2.9426E-01  2.8447E-01 -2.6699E-01 -2.4959E-01  1.1477E-02 -7.8533E-01  1.6743E-01
            -1.0543E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1674.35971724320        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1602
 NPARAMETR:  9.9829E-01  1.8707E+00  5.3556E-01  4.5488E-01  1.3899E+00  9.9371E-01  6.6302E-01  1.0174E+00  1.7511E+00  1.3940E+00
             9.6270E-01
 PARAMETER:  9.8288E-02  7.2631E-01 -5.2444E-01 -6.8772E-01  4.2920E-01  9.3687E-02 -3.1095E-01  1.1722E-01  6.6026E-01  4.3216E-01
             6.1985E-02
 GRADIENT:  -9.2146E+00 -1.1804E+01  9.8964E-01 -4.9504E-01 -4.2465E-01 -9.0040E-01 -3.8837E-01 -6.1701E-02 -2.0816E+00  1.1160E-02
            -1.1262E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1674.36140864578        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1784
 NPARAMETR:  9.9769E-01  1.8791E+00  5.1654E-01  4.4862E-01  1.3841E+00  9.9343E-01  6.6505E-01  9.3984E-01  1.7515E+00  1.3853E+00
             9.6305E-01
 PARAMETER:  9.7691E-02  7.3079E-01 -5.6061E-01 -7.0158E-01  4.2505E-01  9.3410E-02 -3.0789E-01  3.7955E-02  6.6048E-01  4.2595E-01
             6.2353E-02
 GRADIENT:  -1.0751E+01 -1.2041E+01  8.3599E-01 -8.5095E-01 -5.8616E-01 -1.0769E+00 -4.8220E-01 -6.9674E-02 -2.9077E+00 -8.7870E-02
            -1.3699E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1674.41988561446        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     1977            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0034E+00  1.8775E+00  5.1303E-01  4.4949E-01  1.3827E+00  9.9661E-01  6.6770E-01  9.3196E-01  1.7791E+00  1.3857E+00
             9.6312E-01
 PARAMETER:  1.0336E-01  7.2997E-01 -5.6742E-01 -6.9963E-01  4.2401E-01  9.6603E-02 -3.0391E-01  2.9540E-02  6.7610E-01  4.2623E-01
             6.2420E-02
 GRADIENT:   4.9542E+02  9.6806E+02  2.1725E+00  1.1323E+02  2.2973E+01  5.3347E+01  1.8977E+01  2.3953E-01  2.7647E+01  6.3169E+00
             9.8477E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1674.42530794293        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2158
 NPARAMETR:  1.0025E+00  1.8787E+00  5.1287E-01  4.4882E-01  1.3825E+00  9.9609E-01  6.6612E-01  9.1015E-01  1.7770E+00  1.3837E+00
             9.6300E-01
 PARAMETER:  1.0252E-01  7.3056E-01 -5.6773E-01 -7.0113E-01  4.2388E-01  9.6087E-02 -3.0628E-01  5.8550E-03  6.7493E-01  4.2475E-01
             6.2295E-02
 GRADIENT:   1.3799E-01 -1.4873E+01 -7.0414E-02  8.0857E-01 -4.3789E-02  7.4017E-03  1.8694E-01  1.4584E-02 -1.7794E-01  1.6070E-01
            -5.4447E-02

0ITERATION NO.:   62    OBJECTIVE VALUE:  -1674.42531128090        NO. OF FUNC. EVALS.:  60
 CUMULATIVE NO. OF FUNC. EVALS.:     2218
 NPARAMETR:  1.0025E+00  1.8787E+00  5.1287E-01  4.4882E-01  1.3825E+00  9.9610E-01  6.6612E-01  9.0924E-01  1.7770E+00  1.3837E+00
             9.6304E-01
 PARAMETER:  1.0251E-01  7.3056E-01 -5.6773E-01 -7.0113E-01  4.2390E-01  9.6091E-02 -3.0628E-01  4.8550E-03  6.7493E-01  4.2475E-01
             6.2341E-02
 GRADIENT:   7.8301E+04  1.0981E+04  1.4169E+04 -1.1463E+04  1.0203E-01 -1.5847E-01  2.6207E+04 -8.0343E+04 -1.1972E+04  1.8899E+04
             8.0258E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2218
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1485E-03 -4.8131E-02 -2.1634E-02  3.6947E-02 -5.2437E-02
 SE:             2.9901E-02  2.1839E-02  7.2421E-03  2.3025E-02  2.2661E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6936E-01  2.7529E-02  2.8149E-03  1.0857E-01  2.0671E-02

 ETASHRINKSD(%)  1.0000E-10  2.6838E+01  7.5738E+01  2.2865E+01  2.4081E+01
 ETASHRINKVR(%)  1.0000E-10  4.6473E+01  9.4114E+01  4.0501E+01  4.2364E+01
 EBVSHRINKSD(%)  4.0599E-01  2.4977E+01  7.8733E+01  2.5049E+01  2.0144E+01
 EBVSHRINKVR(%)  8.1033E-01  4.3715E+01  9.5477E+01  4.3824E+01  3.6230E+01
 RELATIVEINF(%)  9.9144E+01  4.6504E+00  7.2768E-01  5.0104E+00  2.4780E+01
 EPSSHRINKSD(%)  4.5212E+01
 EPSSHRINKVR(%)  6.9983E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1674.4253112809008     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -939.27448471716264     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1674.425       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.88E+00  5.13E-01  4.49E-01  1.38E+00  9.96E-01  6.66E-01  9.09E-01  1.78E+00  1.38E+00  9.63E-01
 


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
 #CPUT: Total CPU Time in Seconds,       68.428
Stop Time:
Sun Oct 24 02:54:04 CDT 2021
