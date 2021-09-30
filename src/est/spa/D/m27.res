Wed Sep 29 19:54:11 CDT 2021
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
$DATA ../../../../data/spa/D/dat27.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m27.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   6446.25362896934        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.5914E+02  4.4201E+00 -6.8270E+01 -5.0330E+01  3.3794E+02 -9.3448E+02 -4.2504E+02 -5.3294E+01 -7.2552E+02 -4.6689E+02
            -1.3448E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -714.347260424375        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.3277E+00  1.1780E+00  1.0438E+00  1.8978E+00  1.1745E+00  2.1168E+00  1.7505E+00  9.9596E-01  2.0986E+00  1.4574E+00
             1.2801E+01
 PARAMETER:  3.8343E-01  2.6382E-01  1.4291E-01  7.4072E-01  2.6087E-01  8.4992E-01  6.5990E-01  9.5956E-02  8.4127E-01  4.7663E-01
             2.6495E+00
 GRADIENT:  -2.1382E+01  1.9441E+01 -1.7091E+01  5.0029E+01 -1.3823E+01  5.3054E+01  4.2735E+00  5.1009E+00  2.9784E+01  7.6853E+00
             2.4212E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -749.348200013181        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.2638E+00  1.7938E+00  6.2428E+00  1.7741E+00  8.3811E+00  2.1655E+00  3.5262E+00  9.8148E-01  2.8153E+00  9.2732E+00
             1.0165E+01
 PARAMETER:  3.3409E-01  6.8432E-01  1.9314E+00  6.7328E-01  2.2260E+00  8.7265E-01  1.3602E+00  8.1309E-02  1.1351E+00  2.3271E+00
             2.4190E+00
 GRADIENT:  -1.0714E+01  3.4997E+01  1.3309E+00  4.7232E+01  2.3831E-01  4.1461E+01  2.7552E+01 -2.0047E-02  3.9263E+01  2.5985E-02
             2.0047E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -796.183907214570        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  1.1589E+00  1.3932E+00  3.5400E+00  1.0488E+00  4.4084E+00  1.5792E+00  1.6086E+00  3.5567E+00  1.6706E+00  6.3057E+00
             8.4865E+00
 PARAMETER:  2.4748E-01  4.3161E-01  1.3641E+00  1.4767E-01  1.5835E+00  5.5693E-01  5.7535E-01  1.3688E+00  6.1321E-01  1.9414E+00
             2.2385E+00
 GRADIENT:   5.5129E+00 -1.8369E+01  5.1788E+00 -3.3311E+01 -1.3058E+01 -4.7486E-01  3.9080E+00 -3.6400E-01  3.7107E+00  1.7213E+01
             5.6382E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -822.355388218606        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      309
 NPARAMETR:  1.1286E+00  9.6934E-01  2.1013E+00  1.5272E+00  4.6449E+00  1.8295E+00  3.0002E+00  2.8828E+00  1.8745E+00  4.2739E+00
             6.5058E+00
 PARAMETER:  2.2096E-01  6.8864E-02  8.4257E-01  5.2347E-01  1.6358E+00  7.0403E-01  1.1987E+00  1.1587E+00  7.2832E-01  1.5525E+00
             1.9727E+00
 GRADIENT:   1.4494E+00  3.4219E+01  7.7962E+00  2.4624E+01 -6.3790E+00 -1.3203E+00  1.6809E+00 -7.7012E+00 -3.8037E+00 -3.1687E+00
            -2.2250E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -876.817782653205        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:      390
 NPARAMETR:  8.1314E-01  8.7013E-02  4.3672E-01  1.2202E+00  4.2192E+01  1.0437E+00  5.6435E-01  3.3448E+00  6.4062E-01  4.5645E+00
             6.8716E+00
 PARAMETER: -1.0685E-01 -2.3417E+00 -7.2847E-01  2.9899E-01  3.8422E+00  1.4278E-01 -4.7209E-01  1.3074E+00 -3.4533E-01  1.6183E+00
             2.0274E+00
 GRADIENT:  -1.9864E+02  3.0391E+01  3.7757E+01  8.2608E+01 -6.9026E-01 -1.6062E+02  5.8604E-01 -1.0257E+00  4.0220E+00 -1.2357E-03
             5.4870E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -930.449168269387        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      465
 NPARAMETR:  6.0639E-01  2.9348E-02  6.6908E-02  6.0922E-01  1.5917E+02  1.2577E+00  1.9846E-01  1.5124E+00  4.1806E-01  5.5289E+00
             6.1487E+00
 PARAMETER: -4.0024E-01 -3.4285E+00 -2.6044E+00 -3.9558E-01  5.1700E+00  3.2926E-01 -1.5172E+00  5.1372E-01 -7.7214E-01  1.8100E+00
             1.9162E+00
 GRADIENT:  -4.6223E+00 -6.7708E-01 -1.3191E+01  1.5318E+02  1.7041E-02 -5.1312E+01  8.1052E-03 -6.2153E+00 -1.7458E+00  5.4231E-05
            -3.4395E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -938.964997598070        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  6.0777E-01  2.8550E-02  6.6082E-02  5.6996E-01  1.7416E+02  1.3308E+00  1.7527E-01  1.4321E+00  3.9994E-01  6.0310E+00
             6.4036E+00
 PARAMETER: -3.9795E-01 -3.4561E+00 -2.6169E+00 -4.6219E-01  5.2600E+00  3.8581E-01 -1.6414E+00  4.5911E-01 -8.1645E-01  1.8969E+00
             1.9569E+00
 GRADIENT:  -2.7442E+01 -1.0425E+00 -7.5025E+00  2.7252E+01  5.6084E-03 -3.1985E+01  5.6847E-03 -5.7690E+00 -3.5198E-01  6.7949E-05
            -1.1425E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -940.867478993659        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  6.2884E-01  3.4780E-02  6.9084E-02  5.8234E-01  1.4882E+02  1.4364E+00  2.1633E-01  1.4841E+00  4.2587E-01  6.1887E+00
             6.5106E+00
 PARAMETER: -3.6387E-01 -3.2587E+00 -2.5724E+00 -4.4070E-01  5.1028E+00  4.6216E-01 -1.4309E+00  4.9481E-01 -7.5362E-01  1.9227E+00
             1.9734E+00
 GRADIENT:  -5.4124E+00  1.8872E+00 -4.8953E-01 -1.0315E+00 -7.6062E-03 -2.7092E+00  5.1742E-02 -4.0437E-01  1.3331E+00  1.2865E-03
             1.4959E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -941.016076661410        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      957
 NPARAMETR:  6.5390E-01  3.3599E-02  7.7291E-02  6.2236E-01  1.4929E+02  1.4563E+00  2.4743E-01  1.5827E+00  3.8108E-01  6.1891E+00
             6.4882E+00
 PARAMETER: -3.2481E-01 -3.2933E+00 -2.4602E+00 -3.7424E-01  5.1059E+00  4.7588E-01 -1.2966E+00  5.5913E-01 -8.6474E-01  1.9228E+00
             1.9700E+00
 GRADIENT:   2.8342E+00 -1.1937E+00 -5.6007E-01  3.2992E-01  1.4296E-02  8.2740E-01  3.0594E-02  4.4844E-01  1.2548E+00  5.1805E-04
            -9.3590E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -941.317899107185        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1138             RESET HESSIAN, TYPE I
 NPARAMETR:  6.5525E-01  3.4924E-02  7.8409E-02  6.2650E-01  1.4972E+02  1.4547E+00  4.1373E-02  1.6636E+00  1.5887E-01  1.2073E+00
             6.5063E+00
 PARAMETER: -3.2274E-01 -3.2546E+00 -2.4458E+00 -3.6761E-01  5.1088E+00  4.7479E-01 -3.0851E+00  6.0900E-01 -1.7397E+00  2.8838E-01
             1.9728E+00
 GRADIENT:   4.0455E+01  8.4793E-01  3.6743E+01  2.0796E+01  1.2602E-02  1.0893E+01  1.9940E-03 -1.3596E-01  1.5208E-01  4.8472E-05
             1.3987E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -941.365986907675        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1325
 NPARAMETR:  6.5405E-01  3.4859E-02  7.8223E-02  6.2442E-01  3.5372E+01  1.4563E+00  1.0000E-02  1.7079E+00  1.1888E-01  2.0377E-01
             6.5075E+00
 PARAMETER: -3.2458E-01 -3.2565E+00 -2.4482E+00 -3.7093E-01  3.6659E+00  4.7590E-01 -4.9732E+00  6.3526E-01 -2.0296E+00 -1.4908E+00
             1.9730E+00
 GRADIENT:  -7.5073E-01  5.3432E-01 -4.4808E-01 -1.3003E+00  1.6108E-02  6.7887E-02  0.0000E+00  1.5197E+00  2.0668E-02  4.7441E-05
            -1.3716E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -941.402063097969        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1502
 NPARAMETR:  6.5478E-01  3.2544E-02  7.8041E-02  6.2467E-01  1.3892E+01  1.4576E+00  1.0000E-02  1.6910E+00  8.0875E-02  5.4151E-02
             6.5089E+00
 PARAMETER: -3.2346E-01 -3.3251E+00 -2.4505E+00 -3.7054E-01  2.7313E+00  4.7679E-01 -5.9307E+00  6.2530E-01 -2.4149E+00 -2.8160E+00
             1.9732E+00
 GRADIENT:  -7.8844E-01  1.1785E-01 -5.2075E-01  4.9033E-02 -1.4396E-02 -1.7378E-01  0.0000E+00 -2.3861E-01 -3.8903E-03  5.3102E-05
            -4.4763E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -941.403712466434        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     1644
 NPARAMETR:  6.5522E-01  3.2578E-02  7.7920E-02  6.2456E-01  1.4016E+01  1.4582E+00  1.0000E-02  1.6926E+00  1.0850E-01  4.1411E-02
             6.5099E+00
 PARAMETER: -3.2278E-01 -3.3241E+00 -2.4521E+00 -3.7071E-01  2.7402E+00  4.7723E-01 -5.8617E+00  6.2626E-01 -2.1210E+00 -3.0842E+00
             1.9733E+00
 GRADIENT:   7.4841E-02  4.4991E-02 -1.1220E+00  2.2321E-01 -2.5689E-03  5.9876E-02  0.0000E+00  5.1356E-01  5.6268E-03  2.9740E-05
            -1.4485E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -941.404382596665        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1817
 NPARAMETR:  6.5541E-01  3.2601E-02  7.7939E-02  6.2450E-01  1.4120E+01  1.4585E+00  1.0000E-02  1.6890E+00  1.1676E-01  4.0851E-02
             6.5119E+00
 PARAMETER: -3.2249E-01 -3.3234E+00 -2.4518E+00 -3.7080E-01  2.7476E+00  4.7740E-01 -5.8617E+00  6.2412E-01 -2.0476E+00 -3.0978E+00
             1.9736E+00
 GRADIENT:   3.4465E-01  2.8641E-02 -8.3013E-01 -2.5597E-01  2.3229E-03  1.4577E-01  0.0000E+00  3.2606E-01  6.1710E-03  2.7748E-05
             7.3561E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -941.404731225139        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2003
 NPARAMETR:  6.5541E-01  3.2630E-02  7.7930E-02  6.2450E-01  1.4110E+01  1.4584E+00  1.0000E-02  1.6873E+00  1.1720E-01  2.7827E-02
             6.5113E+00
 PARAMETER: -3.2250E-01 -3.3225E+00 -2.4519E+00 -3.7081E-01  2.7469E+00  4.7733E-01 -5.8617E+00  6.2315E-01 -2.0439E+00 -3.4817E+00
             1.9735E+00
 GRADIENT:   3.2312E-01  7.5550E-02 -8.7015E-01 -1.0170E-01 -2.0574E-03  1.1329E-01  0.0000E+00  1.4414E-01  3.3481E-03  1.2982E-05
            -8.7711E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -941.404899498964        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     2194             RESET HESSIAN, TYPE I
 NPARAMETR:  6.5532E-01  3.2653E-02  7.7888E-02  6.2438E-01  1.4167E+01  1.4583E+00  1.0000E-02  1.6852E+00  1.1835E-01  1.0000E-02
             6.5111E+00
 PARAMETER: -3.2263E-01 -3.3218E+00 -2.4525E+00 -3.7100E-01  2.7509E+00  4.7725E-01 -5.8617E+00  6.2186E-01 -2.0341E+00 -4.9205E+00
             1.9735E+00
 GRADIENT:   4.0399E+01  1.3294E+00  3.8493E+01  1.6095E+01  2.0949E-02  1.1430E+01  0.0000E+00  2.6306E+00  1.3772E-01  0.0000E+00
             1.5583E+01

0ITERATION NO.:   84    OBJECTIVE VALUE:  -941.404930898542        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     2329
 NPARAMETR:  6.5510E-01  3.2668E-02  7.7787E-02  6.2416E-01  1.4220E+01  1.4581E+00  1.0000E-02  1.6831E+00  1.1838E-01  1.0000E-02
             6.5100E+00
 PARAMETER: -3.2263E-01 -3.3220E+00 -2.4525E+00 -3.7111E-01  2.7516E+00  4.7728E-01 -5.8617E+00  6.2209E-01 -2.0137E+00 -4.9205E+00
             1.9736E+00
 GRADIENT:   4.9896E-02 -3.4967E-03  1.0319E-01  6.3024E-02 -4.2569E-04  8.0875E-03  0.0000E+00  3.8984E-02  2.3039E-04  0.0000E+00
             2.7746E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2329
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.9602E-03 -7.8301E-05  4.3241E-03 -4.9495E-03  6.2932E-06
 SE:             2.9070E-02  1.6718E-05  2.4435E-02  3.4267E-03  2.4229E-06
 N:                     100         100         100         100         100

 P VAL.:         8.1078E-01  2.8216E-06  8.5954E-01  1.4862E-01  9.3934E-03

 ETASHRINKSD(%)  2.6103E+00  9.9944E+01  1.8140E+01  8.8520E+01  9.9992E+01
 ETASHRINKVR(%)  5.1525E+00  1.0000E+02  3.2990E+01  9.8682E+01  1.0000E+02
 EBVSHRINKSD(%)  2.8791E+00  9.9921E+01  1.6136E+01  8.9054E+01  9.9990E+01
 EBVSHRINKVR(%)  5.6752E+00  1.0000E+02  2.9668E+01  9.8802E+01  1.0000E+02
 RELATIVEINF(%)  6.7337E+00  1.0164E-06  2.1822E+00  1.8562E-02  1.1709E-08
 EPSSHRINKSD(%)  1.6916E+01
 EPSSHRINKVR(%)  3.0971E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -941.40493089854237     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -206.25410433480420     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    36.62
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     8.01
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -941.405       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         6.55E-01  3.26E-02  7.79E-02  6.24E-01  1.42E+01  1.46E+00  1.00E-02  1.69E+00  1.21E-01  1.00E-02  6.51E+00
 


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
+        3.39E+00
 
 TH 2
+       -7.60E+00  1.70E+01
 
 TH 3
+       -3.80E+02  8.51E+02  4.25E+04
 
 TH 4
+        7.79E+01 -1.75E+02 -8.72E+03  1.79E+03
 
 TH 5
+        1.26E-02 -2.82E-02 -1.41E+00  2.89E-01  4.67E-05
 
 TH 6
+       -6.23E-01  1.40E+00  6.98E+01 -1.43E+01 -2.31E-03  1.15E-01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        4.81E-01 -1.08E+00 -5.39E+01  1.11E+01  1.79E-03 -8.85E-02  0.00E+00  6.84E-02
 
 TH 9
+       -3.48E-01  7.81E-01  3.90E+01 -8.00E+00 -1.29E-03  6.41E-02  0.00E+00 -4.95E-02  3.58E-02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.94E-01  1.56E+00  7.77E+01 -1.59E+01 -2.58E-03  1.28E-01  0.00E+00 -9.86E-02  7.13E-02  0.00E+00  1.42E-01
 
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
+        1.14E+03
 
 TH 2
+       -3.86E+02  1.83E+04
 
 TH 3
+       -4.79E+02  2.13E+02  4.40E+04
 
 TH 4
+       -4.78E+02 -6.67E+02 -8.91E+03  2.31E+03
 
 TH 5
+        1.79E-01 -4.86E+00 -1.27E+00  4.73E-01  2.36E-03
 
 TH 6
+        2.89E+00  1.58E+01  7.02E+01 -2.49E+01  7.44E-03  8.10E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        7.11E+00 -7.21E+01 -6.51E+01 -4.43E+01 -4.09E-03  1.70E+00  0.00E+00  3.23E+01
 
 TH 9
+        1.69E+00 -3.69E+01  3.93E+01 -1.89E+01  2.75E-03  1.21E+00  0.00E+00  7.00E+00  2.67E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.37E+01 -4.51E+01  8.20E+01 -1.43E+01  5.84E-03  7.50E-01  0.00E+00  3.06E+00  1.08E+00  0.00E+00  9.73E+00
 
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
+        1.15E+03
 
 TH 2
+       -4.82E+02  1.81E+03
 
 TH 3
+       -7.76E+02  1.41E+03  4.54E+04
 
 TH 4
+       -5.12E+02  5.88E+01 -9.25E+03  2.40E+03
 
 TH 5
+        1.86E-01 -4.43E-01 -1.92E+00  3.34E-01  1.93E-04
 
 TH 6
+        8.57E+01  6.06E+01 -2.82E+02  3.49E+00  1.03E-02  8.41E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        2.22E+01 -1.03E+02  1.71E+02 -9.81E+01 -5.36E-03  8.77E-02  0.00E+00  3.34E+01
 
 TH 9
+        4.37E+00 -2.57E+01  9.57E+01 -3.06E+01 -1.63E-03 -7.57E-01  0.00E+00  7.05E+00  2.20E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.89E+00  2.07E+01  7.98E+01 -4.27E+00 -7.24E-03 -6.77E+00  0.00E+00 -7.38E+00 -1.59E+00  0.00E+00  1.75E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       44.709
Stop Time:
Wed Sep 29 19:54:57 CDT 2021
