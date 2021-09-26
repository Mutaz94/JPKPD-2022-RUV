Sat Sep 25 09:43:39 CDT 2021
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
$DATA ../../../../data/spa/S1/dat15.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m15.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1698.41460764132        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.5646E+01 -1.0829E+02 -4.8862E+01 -9.3956E+01  1.1220E+02  3.7846E+01 -1.0191E+01 -2.5064E+00 -3.0629E+01  9.5771E-01
             3.0774E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1715.61816019410        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0423E+00  1.0605E+00  9.9141E-01  1.0319E+00  9.2928E-01  8.8772E-01  1.0105E+00  1.0784E+00  1.1416E+00  8.8924E-01
             9.2520E-01
 PARAMETER:  1.4138E-01  1.5873E-01  9.1371E-02  1.3136E-01  2.6657E-02 -1.9104E-02  1.1040E-01  1.7550E-01  2.3247E-01 -1.7387E-02
             2.2255E-02
 GRADIENT:   6.1181E+01 -6.8803E+00 -4.5650E+00  8.1697E+00  1.1563E+01 -1.0232E+00  9.1452E-01 -3.1015E+00  7.6881E+00  4.5170E-01
             1.0220E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1716.56092500869        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0366E+00  1.2133E+00  1.1251E+00  9.5045E-01  1.0310E+00  8.8960E-01  9.1786E-01  1.5436E+00  1.1956E+00  9.1815E-01
             9.1385E-01
 PARAMETER:  1.3593E-01  2.9337E-01  2.1788E-01  4.9180E-02  1.3057E-01 -1.6986E-02  1.4286E-02  5.3414E-01  2.7862E-01  1.4609E-02
             9.9121E-03
 GRADIENT:   4.6355E+01  1.8184E+01  2.9874E-02  1.1277E+01  9.1102E+00 -2.6050E-02  3.9866E+00  9.6062E-01 -1.5256E-02 -2.0511E+00
            -4.6454E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1717.50624978487        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      320
 NPARAMETR:  1.0524E+00  1.3833E+00  9.6747E-01  8.3848E-01  1.0394E+00  9.0211E-01  7.7661E-01  1.4685E+00  1.3417E+00  9.4601E-01
             9.2820E-01
 PARAMETER:  1.5103E-01  4.2445E-01  6.6930E-02 -7.6163E-02  1.3868E-01 -3.0231E-03 -1.5281E-01  4.8426E-01  3.9397E-01  4.4497E-02
             2.5489E-02
 GRADIENT:   1.4692E+01  9.9849E+00  4.1115E+00  8.3304E+00 -1.0343E+01  1.5422E+00 -2.0158E+00 -1.4050E+00 -1.2594E-02  7.9299E-01
             1.7040E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1719.07399345790        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      496
 NPARAMETR:  1.0445E+00  1.7534E+00  6.3863E-01  5.7797E-01  1.0995E+00  8.9581E-01  7.4455E-01  1.4467E+00  1.6724E+00  9.4291E-01
             9.2162E-01
 PARAMETER:  1.4356E-01  6.6155E-01 -3.4843E-01 -4.4823E-01  1.9490E-01 -1.0022E-02 -1.9498E-01  4.6927E-01  6.1428E-01  4.1213E-02
             1.8373E-02
 GRADIENT:  -8.0772E+00  3.7552E+00  2.6991E+00  4.2637E-01 -7.4155E+00 -1.6126E+00 -1.3953E+00 -2.6800E-01 -1.0474E-01 -1.1059E-01
            -1.5282E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1719.36173217036        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      663
 NPARAMETR:  1.0493E+00  1.8519E+00  5.5026E-01  5.0491E-01  1.1244E+00  9.0420E-01  7.3960E-01  1.4028E+00  1.7949E+00  9.5021E-01
             9.2528E-01
 PARAMETER:  1.4811E-01  7.1623E-01 -4.9736E-01 -5.8338E-01  2.1727E-01 -7.1011E-04 -2.0164E-01  4.3844E-01  6.8495E-01  4.8924E-02
             2.2341E-02
 GRADIENT:   7.9806E+01  9.2923E+01  2.0225E+00  8.2960E+00 -6.2930E-01  5.6075E+00  1.4093E+00 -8.7843E-02  2.5491E+00 -1.6687E-01
            -9.5708E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1719.38305369496        NO. OF FUNC. EVALS.: 155
 CUMULATIVE NO. OF FUNC. EVALS.:      818
 NPARAMETR:  1.0477E+00  1.8519E+00  5.5026E-01  5.0817E-01  1.1258E+00  8.9950E-01  7.3948E-01  1.4028E+00  1.8003E+00  9.5261E-01
             9.2550E-01
 PARAMETER:  1.4656E-01  7.1623E-01 -4.9736E-01 -5.7694E-01  2.1846E-01 -5.9196E-03 -2.0181E-01  4.3844E-01  6.8797E-01  5.1453E-02
             2.2579E-02
 GRADIENT:  -6.3152E-02 -2.0069E+00  4.9827E-01  1.3116E-01  2.4589E-01 -4.0793E-02 -5.8507E-02 -2.0134E-02 -3.4886E-02 -6.0517E-02
            -3.6191E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1719.38524264604        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      956
 NPARAMETR:  1.0471E+00  1.8526E+00  5.4785E-01  5.0768E-01  1.1254E+00  8.9924E-01  7.3956E-01  1.4017E+00  1.8000E+00  9.5266E-01
             9.2552E-01
 PARAMETER:  1.4607E-01  7.1657E-01 -5.0175E-01 -5.7791E-01  2.1814E-01 -6.2089E-03 -2.0170E-01  4.3768E-01  6.8777E-01  5.1506E-02
             2.2597E-02
 GRADIENT:  -1.4166E+00 -2.4721E+00  2.0221E-01  2.8118E-01  5.7888E-01 -1.5917E-01 -3.2676E-02  5.7250E-02 -7.3019E-02  6.8198E-02
             7.6343E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1719.42753096946        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1138
 NPARAMETR:  1.0479E+00  1.8806E+00  5.1013E-01  4.8865E-01  1.1258E+00  9.0097E-01  7.3807E-01  1.3186E+00  1.8239E+00  9.4734E-01
             9.2503E-01
 PARAMETER:  1.4676E-01  7.3159E-01 -5.7310E-01 -6.1610E-01  2.1847E-01 -4.2878E-03 -2.0372E-01  3.7660E-01  7.0099E-01  4.5901E-02
             2.2071E-02
 GRADIENT:   3.3438E-01 -5.5432E-01 -3.8969E-01  8.2590E-01  2.0985E+00  5.2602E-01 -4.8469E-01 -6.1042E-02 -1.2683E+00 -3.8938E-01
            -4.4199E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1719.43449753201        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1304
 NPARAMETR:  1.0476E+00  1.8808E+00  5.0985E-01  4.8819E-01  1.1245E+00  8.9972E-01  7.3950E-01  1.3188E+00  1.8359E+00  9.4893E-01
             9.2588E-01
 PARAMETER:  1.4650E-01  7.3170E-01 -5.7364E-01 -6.1704E-01  2.1737E-01 -5.6752E-03 -2.0178E-01  3.7674E-01  7.0754E-01  4.7577E-02
             2.2987E-02
 GRADIENT:  -2.7091E-01 -1.0746E+00 -2.8309E-01  8.4723E-01 -7.5576E-02 -1.8743E-02  1.2805E-01  3.0533E-02  1.2418E-03  1.9018E-01
             7.1193E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1719.43613741848        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1478
 NPARAMETR:  1.0476E+00  1.8815E+00  5.1047E-01  4.8716E-01  1.1246E+00  8.9988E-01  7.3909E-01  1.3195E+00  1.8366E+00  9.4780E-01
             9.2582E-01
 PARAMETER:  1.4649E-01  7.3206E-01 -5.7242E-01 -6.1917E-01  2.1744E-01 -5.4979E-03 -2.0233E-01  3.7723E-01  7.0791E-01  4.6391E-02
             2.2921E-02
 GRADIENT:  -3.9418E-01 -9.5248E-01  1.6481E-01  1.8031E-01 -8.4713E-01  7.3145E-02  2.3924E-02 -3.9133E-02 -2.0552E-01 -2.2154E-02
            -1.6453E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1719.43659391419        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1661
 NPARAMETR:  1.0477E+00  1.8812E+00  5.1040E-01  4.8675E-01  1.1253E+00  8.9972E-01  7.3882E-01  1.3207E+00  1.8390E+00  9.4834E-01
             9.2582E-01
 PARAMETER:  1.4662E-01  7.3192E-01 -5.7256E-01 -6.2000E-01  2.1809E-01 -5.6690E-03 -2.0270E-01  3.7818E-01  7.0923E-01  4.6961E-02
             2.2922E-02
 GRADIENT:   1.4591E-02 -2.4491E+00 -2.9550E-02  3.7984E-02 -1.2331E-01  2.9573E-04  1.7447E-02 -2.2947E-02 -2.3870E-02 -2.4643E-02
             9.4241E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1719.43728407553        NO. OF FUNC. EVALS.: 152
 CUMULATIVE NO. OF FUNC. EVALS.:     1813
 NPARAMETR:  1.0477E+00  1.8812E+00  5.1065E-01  4.8644E-01  1.1258E+00  8.9971E-01  7.3864E-01  1.3217E+00  1.8389E+00  9.4874E-01
             9.2583E-01
 PARAMETER:  1.4651E-01  7.3220E-01 -5.7231E-01 -6.2038E-01  2.1829E-01 -5.7032E-03 -2.0297E-01  3.7910E-01  7.0949E-01  4.7341E-02
             2.2890E-02
 GRADIENT:  -1.9626E+05  7.8533E+04 -1.0046E+05  9.2685E+04 -1.8536E-01 -1.3294E-02 -3.8163E-03  1.5162E+05  8.1019E+04 -2.8753E+05
            -5.7508E+05
 NUMSIGDIG:         3.3         3.3         3.3         3.3         2.9         3.5         4.0         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1813
 NO. OF SIG. DIGITS IN FINAL EST.:  2.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.3721E-04 -4.0467E-02 -2.4454E-02  3.5105E-02 -4.8608E-02
 SE:             2.9860E-02  2.2971E-02  9.8189E-03  2.3220E-02  2.1482E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8297E-01  7.8129E-02  1.2756E-02  1.3058E-01  2.3656E-02

 ETASHRINKSD(%)  1.0000E-10  2.3044E+01  6.7106E+01  2.2209E+01  2.8031E+01
 ETASHRINKVR(%)  1.0000E-10  4.0778E+01  8.9180E+01  3.9486E+01  4.8205E+01
 EBVSHRINKSD(%)  4.7575E-01  2.2541E+01  6.9525E+01  2.2672E+01  2.5816E+01
 EBVSHRINKVR(%)  9.4923E-01  4.0001E+01  9.0713E+01  4.0204E+01  4.4968E+01
 RELATIVEINF(%)  9.9026E+01  5.0794E+00  1.1057E+00  5.3866E+00  1.6175E+01
 EPSSHRINKSD(%)  4.6271E+01
 EPSSHRINKVR(%)  7.1132E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1719.4372840755320     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -984.28645751179386     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.89
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1719.437       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.88E+00  5.11E-01  4.87E-01  1.13E+00  9.00E-01  7.39E-01  1.32E+00  1.84E+00  9.49E-01  9.26E-01
 


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
+        6.10E+08
 
 TH 2
+        5.72E+02  7.57E+06
 
 TH 3
+       -2.72E+03  4.09E+03  1.68E+08
 
 TH 4
+        2.63E+03 -8.31E+03  4.03E+04  1.58E+08
 
 TH 5
+       -3.81E+08  4.12E+02 -3.01E+03  2.88E+03  5.56E+02
 
 TH 6
+        1.04E+09  6.53E+02 -3.08E+03  2.98E+03  1.13E+00  2.47E+02
 
 TH 7
+        6.25E+08  4.17E+02 -1.97E+03  1.89E+03 -3.90E+08  1.07E+09  6.40E+08
 
 TH 8
+        1.59E+03  2.08E+07  8.97E+04 -2.36E+04  1.61E+03  1.80E+03  1.15E+03  5.72E+07
 
 TH 9
+        6.12E+02 -6.82E+03  3.21E+04 -3.10E+04  6.13E+02  6.91E+02  4.49E+02 -1.87E+04  8.43E+06
 
 TH10
+       -8.40E+03  3.34E+02 -1.65E+03  1.57E+03 -8.60E+03 -9.49E+03 -6.03E+03  9.53E+02  3.68E+02  1.60E+09
 
 TH11
+       -8.61E+03 -5.20E+03  2.44E+04 -2.37E+04 -8.74E+03 -9.72E+03 -6.19E+03 -1.42E+04  1.01E+05 -5.13E+03  1.68E+09
 
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
 #CPUT: Total CPU Time in Seconds,       28.980
Stop Time:
Sat Sep 25 09:44:10 CDT 2021
