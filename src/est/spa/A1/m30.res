Sat Sep 25 07:59:12 CDT 2021
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
$DATA ../../../../data/spa/A1/dat30.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1212.55803288398        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.3332E+01  2.1934E+01  2.2792E+01  1.6474E+01  6.4025E+01  7.4195E+00 -3.7044E+01 -3.8153E+00 -2.5136E+01 -6.6558E+01
            -8.2472E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1458.87974084454        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0224E+00  9.5713E-01  1.0992E+00  1.0324E+00  9.8986E-01  9.6779E-01  1.2402E+00  8.1568E-01  1.1079E+00  1.1944E+00
             1.8214E+00
 PARAMETER:  1.2214E-01  5.6184E-02  1.9455E-01  1.3191E-01  8.9808E-02  6.7260E-02  3.1528E-01 -1.0373E-01  2.0249E-01  2.7761E-01
             6.9963E-01
 GRADIENT:   3.1212E+01  8.3592E+00  9.6505E+00 -1.8713E+00 -1.1970E+01 -5.2089E+00  2.5285E+00  3.4109E+00  1.5810E+01 -3.6470E+00
            -9.2607E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1464.86830335466        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0215E+00  8.6638E-01  9.0735E-01  1.0816E+00  9.0298E-01  1.0161E+00  1.1594E+00  1.5827E-01  9.5220E-01  1.1035E+00
             1.9171E+00
 PARAMETER:  1.2123E-01 -4.3429E-02  2.7680E-03  1.7848E-01 -2.0533E-03  1.1600E-01  2.4789E-01 -1.7435E+00  5.1017E-02  1.9848E-01
             7.5081E-01
 GRADIENT:   1.9085E+01 -8.8937E+00 -2.1320E+01  9.0071E+00  3.9253E+01  1.1995E+01 -1.1510E+01  2.0480E-01 -9.8626E+00 -2.8258E+00
            -5.8443E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1471.45036249992        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0169E+00  6.1449E-01  8.3854E-01  1.2273E+00  7.2440E-01  9.8473E-01  1.7103E+00  8.7405E-02  8.6345E-01  9.2279E-01
             2.1386E+00
 PARAMETER:  1.1678E-01 -3.8696E-01 -7.6089E-02  3.0480E-01 -2.2242E-01  8.4607E-02  6.3668E-01 -2.3372E+00 -4.6817E-02  1.9648E-02
             8.6017E-01
 GRADIENT:   2.9188E+00  1.1991E+01  6.1466E+00  1.5143E+01 -1.0493E+01  8.0338E-01  1.6705E+00  8.7456E-02 -4.6823E+00 -2.7184E-01
             7.2830E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1472.81076653030        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0121E+00  4.0780E-01  8.4217E-01  1.3289E+00  6.6873E-01  9.7965E-01  2.0491E+00  4.1379E-02  8.6177E-01  9.3123E-01
             2.0983E+00
 PARAMETER:  1.1206E-01 -7.9699E-01 -7.1779E-02  3.8436E-01 -3.0238E-01  7.9438E-02  8.1740E-01 -3.0850E+00 -4.8771E-02  2.8753E-02
             8.4112E-01
 GRADIENT:   7.4353E-01  1.6631E+00  2.1847E+00  4.5117E-01 -3.6407E+00 -1.7268E-01  4.7294E-01  2.1149E-02 -5.9611E-01 -2.0035E-01
             2.2825E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1474.07146004482        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      382
 NPARAMETR:  1.0078E+00  2.1887E-01  8.8853E-01  1.4387E+00  6.4835E-01  9.7876E-01  2.4754E+00  1.0000E-02  8.4397E-01  9.5588E-01
             2.1003E+00
 PARAMETER:  1.0775E-01 -1.4193E+00 -1.8191E-02  4.6377E-01 -3.3333E-01  7.8536E-02  1.0064E+00 -4.9112E+00 -6.9635E-02  5.4873E-02
             8.4208E-01
 GRADIENT:   8.1076E-01  6.1236E-01 -6.6052E-01  5.5906E+00  1.1175E+00  4.2464E-01 -1.0383E+00  0.0000E+00 -1.9056E+00 -1.4098E+00
            -9.6196E-03

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1474.59609525819        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      459
 NPARAMETR:  1.0049E+00  1.0623E-01  8.6072E-01  1.4980E+00  6.0785E-01  9.7649E-01  3.7573E+00  1.0000E-02  8.3057E-01  9.6570E-01
             2.0923E+00
 PARAMETER:  1.0491E-01 -2.1421E+00 -4.9990E-02  5.0411E-01 -3.9782E-01  7.6206E-02  1.4237E+00 -7.1421E+00 -8.5645E-02  6.5101E-02
             8.3825E-01
 GRADIENT:  -1.0192E-01  1.2288E+00  1.8401E+00  1.7526E+01 -7.1932E+00 -1.2387E-01  4.5766E-01  0.0000E+00 -1.3147E+00  2.0765E+00
             9.8260E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1474.60307520350        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      529
 NPARAMETR:  1.0047E+00  1.0773E-01  8.6713E-01  1.4925E+00  6.1275E-01  9.7631E-01  3.6947E+00  1.0000E-02  8.3299E-01  9.6156E-01
             2.0910E+00
 PARAMETER:  1.0467E-01 -2.1281E+00 -4.2572E-02  5.0044E-01 -3.8979E-01  7.6029E-02  1.4069E+00 -7.1274E+00 -8.2733E-02  6.0802E-02
             8.3762E-01
 GRADIENT:  -1.2870E-01  5.7069E-01  8.5734E-01  6.6852E+00 -3.0215E+00 -6.4542E-02  2.2656E-01  0.0000E+00 -4.2454E-01  8.5391E-01
             3.4245E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1474.86420638266        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  1.0016E+00  3.7492E-02  8.7269E-01  1.5389E+00  5.9677E-01  9.7465E-01  5.9770E+00  1.0000E-02  8.2404E-01  9.5756E-01
             2.0865E+00
 PARAMETER:  1.0157E-01 -3.1836E+00 -3.6177E-02  5.3104E-01 -4.1622E-01  7.4322E-02  1.8879E+00 -1.0829E+01 -9.3541E-02  5.6634E-02
             8.3549E-01
 GRADIENT:  -3.4664E+00  4.5097E-01  7.7174E+00  2.7681E+01 -1.5178E+01 -4.9121E-01 -1.5310E-01  0.0000E+00 -2.3317E-01  2.0559E-01
            -1.4133E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1475.16052074756        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      721
 NPARAMETR:  1.0062E+00  3.9816E-02  9.3008E-01  1.5475E+00  6.2913E-01  9.7829E-01  5.7331E+00  1.0000E-02  8.2300E-01  9.9315E-01
             2.0960E+00
 PARAMETER:  1.0618E-01 -3.1235E+00  2.7514E-02  5.3664E-01 -3.6341E-01  7.8055E-02  1.8463E+00 -1.0601E+01 -9.4804E-02  9.3122E-02
             8.4005E-01
 GRADIENT:  -3.1570E+00  1.2461E-01  2.0170E+00 -1.5412E+00 -6.1905E+00  7.0972E-02 -2.3022E-01  0.0000E+00  4.0822E-02  1.3279E+00
             5.7568E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1475.37533097506        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  1.0062E+00  1.0000E-02  9.9084E-01  1.5725E+00  6.5336E-01  9.7671E-01  1.1853E+01  1.0000E-02  8.1588E-01  1.0093E+00
             2.0978E+00
 PARAMETER:  1.0616E-01 -4.5126E+00  9.0793E-02  5.5265E-01 -3.2562E-01  7.6431E-02  2.5726E+00 -1.5556E+01 -1.0348E-01  1.0923E-01
             8.4090E-01
 GRADIENT:  -9.5673E-01  0.0000E+00  7.3192E-01 -3.4971E+00 -4.5379E-01 -1.4634E-01 -2.7388E-02  0.0000E+00  2.5295E-01 -5.1888E-01
            -1.7774E-01

0ITERATION NO.:   54    OBJECTIVE VALUE:  -1475.37986964741        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1026
 NPARAMETR:  1.0067E+00  1.0000E-02  9.9662E-01  1.5754E+00  6.5650E-01  9.7714E-01  1.2095E+01  1.0000E-02  8.1495E-01  1.0161E+00
             2.0979E+00
 PARAMETER:  1.0664E-01 -4.5538E+00  9.6615E-02  5.5450E-01 -3.2083E-01  7.6875E-02  2.5928E+00 -1.5698E+01 -1.0463E-01  1.1598E-01
             8.4094E-01
 GRADIENT:   7.8874E-03  0.0000E+00  5.3145E-04  2.5615E-02 -1.2689E-02  2.1953E-03 -1.2130E-03  0.0000E+00 -3.4005E-03  7.8460E-03
             7.5461E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1026
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.8736E-04  7.2470E-04 -1.1478E-05 -1.0727E-02 -2.4501E-02
 SE:             2.9375E-02  1.9674E-03  1.3721E-04  2.7696E-02  2.2379E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9219E-01  7.1262E-01  9.3333E-01  6.9854E-01  2.7358E-01

 ETASHRINKSD(%)  1.5914E+00  9.3409E+01  9.9540E+01  7.2147E+00  2.5029E+01
 ETASHRINKVR(%)  3.1575E+00  9.9566E+01  9.9998E+01  1.3909E+01  4.3794E+01
 EBVSHRINKSD(%)  1.7163E+00  9.4037E+01  9.9488E+01  6.6789E+00  2.3665E+01
 EBVSHRINKVR(%)  3.4032E+00  9.9644E+01  9.9997E+01  1.2912E+01  4.1729E+01
 RELATIVEINF(%)  8.8444E+01  1.0814E-02  1.7752E-04  3.5034E+00  3.3289E+00
 EPSSHRINKSD(%)  3.5301E+01
 EPSSHRINKVR(%)  5.8140E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1475.3798696474137     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -740.22904308367549     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.08
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.25
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1475.380       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.00E-02  9.97E-01  1.58E+00  6.56E-01  9.77E-01  1.21E+01  1.00E-02  8.15E-01  1.02E+00  2.10E+00
 


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
+        1.12E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.37E+00  0.00E+00  3.25E+02
 
 TH 4
+       -3.16E+01  0.00E+00 -5.48E+01  5.90E+02
 
 TH 5
+        9.23E+00  0.00E+00 -6.06E+02 -8.02E+01  1.46E+03
 
 TH 6
+        4.26E+00  0.00E+00  1.30E+01 -6.50E+00 -2.15E+00  2.24E+02
 
 TH 7
+        1.34E-02  0.00E+00  7.72E-03 -3.10E-02  4.26E-02 -1.22E-02  6.69E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.21E+00  0.00E+00  4.00E+01 -1.24E+01 -1.48E+01  1.03E+01  2.23E-02  0.00E+00  2.57E+02
 
 TH10
+       -2.06E+01  0.00E+00  2.32E+00 -2.66E+00 -6.06E+01  8.04E+00 -2.87E-02  0.00E+00 -8.43E+00  7.59E+01
 
 TH11
+       -1.36E+01  0.00E+00 -1.09E+01 -9.55E+00  4.66E+00  9.24E-02 -2.28E-03  0.00E+00  6.97E+00  1.89E+01  6.25E+01
 
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
 #CPUT: Total CPU Time in Seconds,       16.398
Stop Time:
Sat Sep 25 07:59:30 CDT 2021
