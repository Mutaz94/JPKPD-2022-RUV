Sat Sep 18 11:36:52 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat30.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1683.21058187785        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.2780E+01 -1.7617E+01 -3.3897E+01  3.5296E+01  3.2846E+01  1.9667E+00 -1.5676E+01  2.0205E+00  5.4391E+00 -1.0722E+01
             2.0021E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1689.48507806382        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9889E-01  1.1749E+00  1.1264E+00  8.9421E-01  1.1297E+00  9.9813E-01  1.2375E+00  1.0171E+00  9.1454E-01  1.1755E+00
             1.0082E+00
 PARAMETER:  9.8884E-02  2.6122E-01  2.1899E-01 -1.1817E-02  2.2193E-01  9.8129E-02  3.1312E-01  1.1692E-01  1.0667E-02  2.6166E-01
             1.0816E-01
 GRADIENT:   2.0662E+01  2.3305E+01  1.1036E+00  1.9000E+01  2.8826E+00  1.4540E+00  7.0630E+00 -3.9517E+00  2.5719E+00  3.0844E+00
             5.0166E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1690.37624923124        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      266
 NPARAMETR:  1.0126E+00  1.3390E+00  1.0729E+00  7.9651E-01  1.1753E+00  1.0246E+00  1.1440E+00  1.3973E+00  9.6785E-01  1.1433E+00
             9.8068E-01
 PARAMETER:  1.1253E-01  3.9194E-01  1.7035E-01 -1.2752E-01  2.6151E-01  1.2434E-01  2.3454E-01  4.3452E-01  6.7319E-02  2.3393E-01
             8.0491E-02
 GRADIENT:   1.0578E-01  1.8107E+01 -2.4400E+00  1.6951E+01 -2.4127E-01  4.7303E+00  7.5273E+00  7.8226E-01  2.2921E+00 -2.4287E+00
            -5.1630E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1690.89129113494        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  1.0125E+00  1.3811E+00  1.0855E+00  7.5565E-01  1.2072E+00  1.0050E+00  1.0709E+00  1.4501E+00  9.7653E-01  1.1796E+00
             9.9309E-01
 PARAMETER:  1.1238E-01  4.2291E-01  1.8204E-01 -1.8018E-01  2.8827E-01  1.0503E-01  1.6848E-01  4.7164E-01  7.6248E-02  2.6513E-01
             9.3063E-02
 GRADIENT:  -4.1738E-01  9.9034E-01  1.3557E+00 -5.6547E-01 -1.2738E+00 -2.8302E+00  4.7382E-01 -6.6627E-01 -1.2603E+00 -5.1838E-01
            -8.5987E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1690.92857442240        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      627
 NPARAMETR:  1.0126E+00  1.4005E+00  1.0751E+00  7.4376E-01  1.2121E+00  1.0080E+00  1.0553E+00  1.4723E+00  9.9362E-01  1.1821E+00
             9.9334E-01
 PARAMETER:  1.1249E-01  4.3686E-01  1.7242E-01 -1.9603E-01  2.9232E-01  1.0793E-01  1.5381E-01  4.8679E-01  9.3600E-02  2.6727E-01
             9.3315E-02
 GRADIENT:  -3.8209E-01  1.6109E+00  1.3563E+00  1.4201E-01 -1.7835E+00 -1.7309E+00  1.6359E-01 -5.0601E-01 -8.6145E-01 -3.1055E-01
            -6.7897E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1690.94765403039        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  1.0127E+00  1.4160E+00  1.0610E+00  7.3367E-01  1.2150E+00  1.0088E+00  1.0460E+00  1.4785E+00  1.0035E+00  1.1831E+00
             9.9341E-01
 PARAMETER:  1.1260E-01  4.4783E-01  1.5917E-01 -2.0970E-01  2.9477E-01  1.0877E-01  1.4493E-01  4.9101E-01  1.0347E-01  2.6817E-01
             9.3391E-02
 GRADIENT:  -3.3632E-01  1.6399E+00  1.2724E+00  2.3189E-01 -1.7135E+00 -1.4471E+00  1.2237E-01 -4.4334E-01 -7.4114E-01 -2.6541E-01
            -6.0312E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1690.95518272589        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      981
 NPARAMETR:  1.0129E+00  1.4157E+00  1.0609E+00  7.3374E-01  1.2149E+00  1.0117E+00  1.0414E+00  1.4781E+00  1.0160E+00  1.1830E+00
             9.9440E-01
 PARAMETER:  1.1282E-01  4.4761E-01  1.5909E-01 -2.0960E-01  2.9462E-01  1.1159E-01  1.4057E-01  4.9076E-01  1.1589E-01  2.6804E-01
             9.4384E-02
 GRADIENT:   1.2888E-01  7.9547E-01  9.5621E-01  4.4779E-01 -1.5657E+00 -3.2976E-01  9.7235E-02 -2.8006E-01  7.1310E-02 -1.1359E-01
            -9.0036E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1690.95910929800        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1159
 NPARAMETR:  1.0129E+00  1.4158E+00  1.0578E+00  7.3293E-01  1.2158E+00  1.0125E+00  1.0411E+00  1.4847E+00  1.0135E+00  1.1834E+00
             9.9425E-01
 PARAMETER:  1.1278E-01  4.4772E-01  1.5615E-01 -2.1071E-01  2.9542E-01  1.1245E-01  1.4026E-01  4.9522E-01  1.1343E-01  2.6841E-01
             9.4232E-02
 GRADIENT:   6.8126E-03 -6.3200E-01  2.6121E-01  2.9681E-01 -3.2395E-01 -7.3707E-04 -7.2802E-02 -7.0151E-02 -4.7503E-02 -5.0134E-02
            -4.3336E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1691.03230603271        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1338
 NPARAMETR:  1.0138E+00  1.5797E+00  9.1502E-01  6.2849E-01  1.2511E+00  1.0138E+00  9.6107E-01  1.5684E+00  1.1093E+00  1.1979E+00
             9.9600E-01
 PARAMETER:  1.1371E-01  5.5721E-01  1.1187E-02 -3.6443E-01  3.2402E-01  1.1369E-01  6.0297E-02  5.5005E-01  2.0373E-01  2.8055E-01
             9.5987E-02
 GRADIENT:   1.8525E-01  5.2968E+00  5.2108E-01  2.8700E+00 -8.7373E-02  3.0510E-02 -5.9111E-01  5.2693E-03 -2.1282E-01 -5.6509E-03
             6.3532E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1691.05744222903        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1513
 NPARAMETR:  1.0141E+00  1.6383E+00  8.1527E-01  5.8455E-01  1.2511E+00  1.0142E+00  9.3928E-01  1.4888E+00  1.1537E+00  1.1950E+00
             9.9478E-01
 PARAMETER:  1.1400E-01  5.9366E-01 -1.0423E-01 -4.3692E-01  3.2404E-01  1.1410E-01  3.7356E-02  4.9794E-01  2.4301E-01  2.7819E-01
             9.4771E-02
 GRADIENT:   3.5651E-02 -5.0354E-01 -1.6547E-01 -1.5738E-01  4.1642E-01  2.7016E-03  5.0174E-02  1.3195E-01 -4.6700E-02  1.0949E-01
             3.2391E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1691.06051139318        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1688
 NPARAMETR:  1.0140E+00  1.6224E+00  8.2453E-01  5.9515E-01  1.2448E+00  1.0141E+00  9.4535E-01  1.4546E+00  1.1443E+00  1.1914E+00
             9.9413E-01
 PARAMETER:  1.1393E-01  5.8393E-01 -9.2944E-02 -4.1894E-01  3.1897E-01  1.1400E-01  4.3796E-02  4.7474E-01  2.3479E-01  2.7510E-01
             9.4108E-02
 GRADIENT:  -1.9161E-03  3.5505E-02 -5.3034E-03  2.5831E-02  2.0113E-02 -5.1950E-04 -3.7242E-03  5.4391E-04 -6.8255E-03  2.6840E-03
             3.8849E-03

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1691.06051139318        NO. OF FUNC. EVALS.:  25
 CUMULATIVE NO. OF FUNC. EVALS.:     1713
 NPARAMETR:  1.0140E+00  1.6224E+00  8.2453E-01  5.9515E-01  1.2448E+00  1.0141E+00  9.4535E-01  1.4546E+00  1.1443E+00  1.1914E+00
             9.9413E-01
 PARAMETER:  1.1393E-01  5.8393E-01 -9.2944E-02 -4.1894E-01  3.1897E-01  1.1400E-01  4.3796E-02  4.7474E-01  2.3479E-01  2.7510E-01
             9.4108E-02
 GRADIENT:  -3.0402E-02  1.0540E-01 -4.7415E-03  3.2162E-02  2.1190E-02 -3.5132E-03 -3.9203E-03  4.1921E-04 -7.2026E-03  2.8345E-03
             3.0963E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1713
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.6606E-05 -1.8068E-02 -3.4415E-02  1.6228E-02 -4.0963E-02
 SE:             2.9861E-02  2.4554E-02  1.1357E-02  1.9452E-02  2.2218E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9875E-01  4.6184E-01  2.4440E-03  4.0414E-01  6.5231E-02

 ETASHRINKSD(%)  1.0000E-10  1.7740E+01  6.1952E+01  3.4832E+01  2.5565E+01
 ETASHRINKVR(%)  1.0000E-10  3.2333E+01  8.5524E+01  5.7532E+01  4.4595E+01
 EBVSHRINKSD(%)  4.1341E-01  1.7323E+01  6.7213E+01  3.8432E+01  2.1214E+01
 EBVSHRINKVR(%)  8.2510E-01  3.1645E+01  8.9250E+01  6.2094E+01  3.7927E+01
 RELATIVEINF(%)  9.8806E+01  2.3479E+00  8.1060E-01  1.1317E+00  1.9367E+01
 EPSSHRINKSD(%)  4.5132E+01
 EPSSHRINKVR(%)  6.9895E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1691.0605113931829     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -955.90968482944470     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.41
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.58
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1691.061       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.62E+00  8.25E-01  5.95E-01  1.24E+00  1.01E+00  9.45E-01  1.45E+00  1.14E+00  1.19E+00  9.94E-01
 


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
+        1.04E+03
 
 TH 2
+       -6.43E+00  3.37E+02
 
 TH 3
+        7.12E+00  2.08E+08  1.08E+02
 
 TH 4
+       -1.24E+01  3.93E+02 -1.50E+02  9.22E+02
 
 TH 5
+       -1.19E+00 -7.31E+01 -1.15E+02  1.60E+02  3.23E+02
 
 TH 6
+        6.85E-01 -1.33E+00  2.58E+00 -3.01E+00 -2.18E-01  1.93E+02
 
 TH 7
+        1.13E+00  1.55E+01  2.02E+00 -2.13E+01 -3.91E+00 -2.06E+00  1.17E+02
 
 TH 8
+       -3.93E-01 -2.49E+07 -1.54E+01  1.65E+01  1.03E+00  1.39E-01  7.39E-01  4.53E+00
 
 TH 9
+       -9.23E-02 -1.03E+01 -1.57E+01  2.85E+01  4.23E+00  2.18E-01  2.30E+01  3.86E+00  3.12E+01
 
 TH10
+        1.02E+00  5.24E+07 -9.60E+00  6.90E+00 -4.83E+01  1.06E+00  2.27E-01  4.04E+00  4.90E+00  6.33E+01
 
 TH11
+       -5.01E+00 -1.18E+01 -1.24E+01  6.95E+00 -5.18E+00  2.65E+00  7.19E+00  4.36E+00  5.68E+00  7.92E+00  2.13E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.044
Stop Time:
Sat Sep 18 11:37:23 CDT 2021
