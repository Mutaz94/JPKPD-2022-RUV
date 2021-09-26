Sat Sep 25 01:01:07 CDT 2021
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
$DATA ../../../../data/int/SL2/dat21.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      996
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

 TOT. NO. OF OBS RECS:      896
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
 RAW OUTPUT FILE (FILE): m21.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2188.20092200282        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -5.2409E+01  5.0990E+01  1.5045E+02  6.8733E+00  4.0624E+01  3.0966E+00 -6.9386E+01 -2.1552E+02 -6.9811E+01  3.6988E+00
            -3.0618E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3086.17345833720        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0372E+00  1.2915E+00  1.0312E+00  8.3042E-01  1.2350E+00  9.4267E-01  9.9421E-01  9.6499E-01  9.6765E-01  8.0252E-01
             1.9662E+00
 PARAMETER:  1.3651E-01  3.5579E-01  1.3077E-01 -8.5828E-02  3.1104E-01  4.0962E-02  9.4195E-02  6.4357E-02  6.7112E-02 -1.1999E-01
             7.7611E-01
 GRADIENT:  -2.5285E+01 -7.0641E+00 -1.5987E+01 -3.7550E+01  6.0146E+00 -1.7091E+01  6.5836E+00  1.7475E+00 -2.7579E+01 -4.1148E+01
            -1.3721E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3098.14832443009        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0669E+00  1.6779E+00  1.0323E+00  6.8207E-01  1.5322E+00  9.4186E-01  7.9967E-01  5.2327E-01  1.2041E+00  9.9577E-01
             2.0302E+00
 PARAMETER:  1.6476E-01  6.1755E-01  1.3178E-01 -2.8262E-01  5.2674E-01  4.0105E-02 -1.2355E-01 -5.4765E-01  2.8573E-01  9.5758E-02
             8.0811E-01
 GRADIENT:   4.2039E+01  1.2380E+02 -1.0247E+01  7.4265E+01  2.0616E+01 -1.7567E+01 -2.6139E+00 -1.2667E+00 -1.0200E+01 -4.7973E+01
            -6.2111E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3109.81911546933        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0406E+00  1.6481E+00  1.3325E+00  6.2854E-01  1.6187E+00  9.9214E-01  7.4331E-01  1.2104E+00  1.3282E+00  1.2801E+00
             2.0348E+00
 PARAMETER:  1.3981E-01  5.9959E-01  3.8703E-01 -3.6435E-01  5.8161E-01  9.2111E-02 -1.9664E-01  2.9099E-01  3.8382E-01  3.4692E-01
             8.1039E-01
 GRADIENT:  -1.5060E+01 -1.2151E+01 -1.9939E+00  2.2154E+00 -5.1674E+00  3.7496E+00  1.3534E-01 -2.1452E+00  9.3545E-01  7.2135E+00
            -2.3036E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3118.04168491928        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0411E+00  1.8495E+00  1.8853E+00  4.8749E-01  1.8813E+00  9.8976E-01  7.1651E-01  3.3622E+00  1.4342E+00  1.4190E+00
             2.0302E+00
 PARAMETER:  1.4028E-01  7.1493E-01  7.3410E-01 -6.1848E-01  7.3196E-01  8.9710E-02 -2.3336E-01  1.3126E+00  4.6063E-01  4.4993E-01
             8.0811E-01
 GRADIENT:  -1.3616E+01 -3.4388E+01 -7.7985E+00 -1.6633E+01 -7.5057E+00  3.1418E+00 -1.4359E+00  2.8757E+00  4.3870E+00  5.5299E+00
             7.9802E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3119.97914293830        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.0389E+00  1.7124E+00  2.8739E+00  5.8951E-01  1.8608E+00  9.8535E-01  8.6294E-01  3.8178E+00  8.9228E-01  1.4571E+00
             1.9797E+00
 PARAMETER:  1.3815E-01  6.3789E-01  1.1557E+00 -4.2847E-01  7.2098E-01  8.5243E-02 -4.7410E-02  1.4397E+00 -1.3977E-02  4.7648E-01
             7.8293E-01
 GRADIENT:  -1.6643E+01 -1.5344E+01  2.2332E+00 -2.9127E+01 -1.8750E+01  1.3380E+00 -6.1344E-01 -5.2117E+00 -3.2202E+00  8.0398E+00
            -3.7133E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3121.01349879876        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:      535
 NPARAMETR:  1.0497E+00  1.8689E+00  2.5722E+00  5.0561E-01  1.9624E+00  9.9570E-01  8.0463E-01  4.2636E+00  1.0413E+00  1.5098E+00
             2.0017E+00
 PARAMETER:  1.4847E-01  7.2534E-01  1.0448E+00 -5.8199E-01  7.7418E-01  9.5691E-02 -1.1738E-01  1.5501E+00  1.4050E-01  5.1198E-01
             7.9400E-01
 GRADIENT:  -1.2396E+01 -1.0154E+01 -1.2414E-01 -8.4298E+00 -7.2687E+00  4.0646E+00 -7.4565E-01 -1.3928E-01 -6.7414E-01  7.6875E+00
            -1.7886E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3121.04502119405        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:      731             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0498E+00  1.8676E+00  2.5815E+00  5.0661E-01  1.9610E+00  9.9421E-01  8.0570E-01  4.2696E+00  1.0365E+00  1.5074E+00
             2.0029E+00
 PARAMETER:  1.4858E-01  7.2466E-01  1.0484E+00 -5.8002E-01  7.7347E-01  9.4195E-02 -1.1604E-01  1.5515E+00  1.3589E-01  5.1040E-01
             7.9457E-01
 GRADIENT:   6.3324E+00  2.5651E+01  1.7519E-01 -4.1895E+00 -2.2612E+00  4.9187E+00 -3.8014E-01  3.2831E-01 -6.3605E-01  8.1224E+00
            -1.5093E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3121.07634912441        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      925            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0499E+00  1.8673E+00  2.5822E+00  5.0667E-01  1.9614E+00  9.8506E-01  8.0577E-01  4.2680E+00  1.0471E+00  1.5070E+00
             2.0035E+00
 PARAMETER:  1.4870E-01  7.2447E-01  1.0487E+00 -5.7990E-01  7.7367E-01  8.4950E-02 -1.1596E-01  1.5512E+00  1.4605E-01  5.1009E-01
             7.9487E-01
 GRADIENT:   6.3054E+00  2.5206E+01  1.4326E-01 -4.1756E+00 -1.9361E+00  1.4185E+00  1.5974E-01  5.0496E-01 -2.8303E-01  8.1512E+00
            -1.4387E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3121.16600093850        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     1072
 NPARAMETR:  1.0551E+00  1.8669E+00  2.5829E+00  5.0662E-01  1.9618E+00  9.8490E-01  8.0566E-01  4.2661E+00  1.0646E+00  1.5068E+00
             2.0169E+00
 PARAMETER:  1.5360E-01  7.2429E-01  1.0489E+00 -5.7999E-01  7.7387E-01  8.4785E-02 -1.1609E-01  1.5507E+00  1.6264E-01  5.1000E-01
             8.0158E-01
 GRADIENT:  -1.5546E+00 -1.1856E+01 -1.5265E-01 -8.5017E+00 -7.2127E+00  1.2068E-01  8.4432E-01  6.0851E-01  2.7138E-01  7.7548E+00
            -6.5021E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3121.23018510913        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1250            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0558E+00  1.8680E+00  2.5840E+00  5.0780E-01  1.9641E+00  9.8454E-01  8.0321E-01  4.2606E+00  1.0609E+00  1.5008E+00
             2.0177E+00
 PARAMETER:  1.5431E-01  7.2487E-01  1.0493E+00 -5.7768E-01  7.7503E-01  8.4422E-02 -1.1914E-01  1.5494E+00  1.5909E-01  5.0598E-01
             8.0195E-01
 GRADIENT:   1.8941E+01  2.6561E+01 -3.3815E-02 -2.7172E+00 -6.6715E-01  1.3410E+00  3.6058E-01  9.8981E-01  8.1094E-02  7.5817E+00
             1.1648E+00

0ITERATION NO.:   54    OBJECTIVE VALUE:  -3121.26610459944        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     1349
 NPARAMETR:  1.0555E+00  1.8677E+00  2.5854E+00  5.0861E-01  1.9653E+00  9.8444E-01  8.0315E-01  4.2537E+00  1.0612E+00  1.4957E+00
             2.0177E+00
 PARAMETER:  1.5402E-01  7.2472E-01  1.0499E+00 -5.7607E-01  7.7563E-01  8.4309E-02 -1.1920E-01  1.5478E+00  1.5922E-01  5.0256E-01
             8.0194E-01
 GRADIENT:  -6.8109E-01  2.0420E+04 -7.0503E+03  2.5695E+04 -9.5506E+03 -6.6341E-02  3.5716E-02  4.7714E+03 -6.2450E-02 -2.9458E+04
            -4.4287E-01
 NUMSIGDIG:         2.4         3.3         3.3         3.3         3.3         2.4         2.5         3.3         1.2         3.3
                    3.6

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1349
 NO. OF SIG. DIGITS IN FINAL EST.:  1.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2299E-03 -1.3742E-02 -3.3185E-02  2.1952E-02 -2.3936E-02
 SE:             2.9647E-02  2.5427E-02  1.4790E-02  1.6674E-02  2.4491E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6691E-01  5.8889E-01  2.4848E-02  1.8799E-01  3.2840E-01

 ETASHRINKSD(%)  6.7971E-01  1.4817E+01  5.0452E+01  4.4140E+01  1.7953E+01
 ETASHRINKVR(%)  1.3548E+00  2.7439E+01  7.5450E+01  6.8796E+01  3.2684E+01
 EBVSHRINKSD(%)  9.5989E-01  1.4617E+01  5.5413E+01  5.0156E+01  1.1588E+01
 EBVSHRINKVR(%)  1.9106E+00  2.7097E+01  8.0120E+01  7.5156E+01  2.1834E+01
 RELATIVEINF(%)  9.8059E+01  5.8624E+00  7.1487E+00  1.8223E+00  4.4964E+01
 EPSSHRINKSD(%)  1.7813E+01
 EPSSHRINKVR(%)  3.2454E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          896
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1646.7378515027735     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3121.2661045994378     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1474.5282530966642     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    37.99
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3121.266       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  1.87E+00  2.59E+00  5.09E-01  1.97E+00  9.84E-01  8.03E-01  4.25E+00  1.06E+00  1.50E+00  2.02E+00
 


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
+        1.40E+08
 
 TH 2
+        5.77E+01  2.02E+06
 
 TH 3
+        8.39E+06  5.78E+02  5.02E+05
 
 TH 4
+        2.95E+02 -1.80E+03  2.64E+03  4.31E+07
 
 TH 5
+       -6.10E+01 -2.39E+01 -5.18E+02  2.10E+03  1.59E+06
 
 TH 6
+        1.14E+00  8.04E+01 -4.14E+01  3.77E+02 -7.47E+01  1.83E+02
 
 TH 7
+        3.68E-01  3.09E+02 -1.45E+02  1.33E+03 -2.61E+02 -8.02E+00  1.68E+02
 
 TH 8
+        1.36E+01 -2.19E+03  1.17E+02 -4.65E+02  1.95E+03  1.71E+01  5.98E+01  8.50E+04
 
 TH 9
+       -1.51E+00 -1.28E+02  6.19E+01 -5.68E+02  1.14E+02 -1.01E-02  2.96E+01 -2.38E+01  1.45E+01
 
 TH10
+       -1.20E+02  1.49E+01  1.81E+06  1.16E+02 -2.65E+01 -1.50E+02 -5.32E+02  3.33E+00  2.29E+02  6.55E+06
 
 TH11
+        1.41E+07 -4.48E+02  8.43E+05 -2.00E+03  3.82E+02  2.74E+00  9.33E+00 -8.65E+01  3.22E+00  7.87E+02  1.41E+06
 
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
 #CPUT: Total CPU Time in Seconds,       53.329
Stop Time:
Sat Sep 25 01:02:01 CDT 2021
