Sat Sep 25 09:31:21 CDT 2021
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
$DATA ../../../../data/spa/A3/dat81.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m81.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -106.993916370224        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4986E+01  4.2449E+01  6.6867E+01 -3.9164E+01  2.1050E+02 -3.9049E+01 -6.5777E+01 -3.0185E+01 -1.3562E+02 -1.8732E+02
            -2.6668E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1206.66805993032        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0797E+00  9.4327E-01  8.9426E-01  1.1681E+00  8.5729E-01  1.0104E+00  1.0037E+00  1.0040E+00  1.1401E+00  1.0719E+00
             5.4528E+00
 PARAMETER:  1.7664E-01  4.1597E-02 -1.1761E-02  2.5542E-01 -5.3975E-02  1.1034E-01  1.0370E-01  1.0404E-01  2.3110E-01  1.6943E-01
             1.7961E+00
 GRADIENT:   1.3083E+01 -1.1273E+01 -1.6102E+01  2.1474E+00 -5.9512E+00  6.3189E+00  1.1578E+01  7.2161E+00  2.9649E+01  2.7565E+01
             2.0631E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1225.38660169972        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0696E+00  6.0173E-01  3.6485E-01  1.3227E+00  3.7658E-01  1.0396E+00  5.4277E-01  4.1914E-01  1.3857E+00  4.2023E-01
             4.8170E+00
 PARAMETER:  1.6725E-01 -4.0794E-01 -9.0828E-01  3.7969E-01 -8.7662E-01  1.3887E-01 -5.1106E-01 -7.6954E-01  4.2622E-01 -7.6695E-01
             1.6721E+00
 GRADIENT:  -1.3528E+01  6.5841E+01  4.6741E+01  9.9182E+01 -1.1247E+02 -4.5303E+00  9.1177E-01  2.2228E+00  3.7673E+01  4.3630E+00
             1.6050E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1271.63442452539        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.0232E+00  5.0866E-01  2.3672E-01  1.1503E+00  3.0597E-01  1.0759E+00  9.6250E-01  6.4326E-01  1.1988E+00  2.1077E-01
             3.2587E+00
 PARAMETER:  1.2294E-01 -5.7598E-01 -1.3409E+00  2.4004E-01 -1.0843E+00  1.7317E-01  6.1777E-02 -3.4120E-01  2.8131E-01 -1.4570E+00
             1.2813E+00
 GRADIENT:  -2.2553E+01 -2.2232E+00 -2.9412E+01  7.5762E+01  4.9624E+01 -5.5764E+00 -2.0546E+00 -7.3088E-01  8.4947E-02 -1.5101E+00
            -1.7451E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1275.20885900152        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0359E+00  4.7701E-01  1.9709E-01  1.0415E+00  2.6583E-01  1.0912E+00  1.0082E+00  8.2287E-01  1.2336E+00  2.0692E-01
             3.2120E+00
 PARAMETER:  1.3529E-01 -6.4021E-01 -1.5241E+00  1.4065E-01 -1.2249E+00  1.8730E-01  1.0819E-01 -9.4955E-02  3.0994E-01 -1.4754E+00
             1.2669E+00
 GRADIENT:   3.8130E+00  1.5227E+01  1.3423E+01  4.8448E+00 -2.6745E+01 -1.6277E+00 -1.0708E+00  4.3711E-01 -1.5700E+00 -1.0007E+00
            -3.1269E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1275.55355591998        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      457
 NPARAMETR:  1.0369E+00  4.7655E-01  2.1009E-01  1.0596E+00  2.7688E-01  1.0948E+00  1.0407E+00  7.0552E-01  1.2081E+00  2.2149E-01
             3.2669E+00
 PARAMETER:  1.3624E-01 -6.4118E-01 -1.4602E+00  1.5791E-01 -1.1842E+00  1.9062E-01  1.3992E-01 -2.4883E-01  2.8907E-01 -1.4074E+00
             1.2838E+00
 GRADIENT:  -3.3022E-01  7.7730E-01 -6.2679E-01  1.3726E+00  2.2496E+00  2.3814E-02 -3.2735E-02 -1.9029E-01  2.7464E-02 -1.0956E+00
             1.1508E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1276.94369511644        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      636
 NPARAMETR:  1.0426E+00  3.7021E-01  2.2806E-01  1.1122E+00  2.6297E-01  1.0908E+00  9.7328E-01  3.7748E-01  1.1915E+00  5.2383E-01
             3.2313E+00
 PARAMETER:  1.4174E-01 -8.9369E-01 -1.3781E+00  2.0630E-01 -1.2357E+00  1.8691E-01  7.2917E-02 -8.7424E-01  2.7521E-01 -5.4658E-01
             1.2729E+00
 GRADIENT:   5.7244E-01 -9.6269E-01 -3.7550E+00 -1.0263E+00  4.9390E+00 -9.4421E-03  1.0052E+00  7.0264E-01 -2.3435E-01  1.5772E+00
             9.1266E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1277.30894209437        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      811
 NPARAMETR:  1.0434E+00  3.3540E-01  2.4771E-01  1.1491E+00  2.6738E-01  1.0908E+00  5.3058E-01  1.2850E-01  1.1766E+00  5.8636E-01
             3.2463E+00
 PARAMETER:  1.4249E-01 -9.9243E-01 -1.2955E+00  2.3899E-01 -1.2191E+00  1.8694E-01 -5.3379E-01 -1.9518E+00  2.6264E-01 -4.3382E-01
             1.2775E+00
 GRADIENT:  -2.2664E-01  2.2676E+00  4.2830E+00  1.6894E+00 -7.6045E+00  1.1041E+00 -5.1688E-02  7.6030E-02 -2.6132E-01 -1.7471E-01
            -2.5302E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1277.38576632957        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      987
 NPARAMETR:  1.0432E+00  2.8384E-01  2.7074E-01  1.1908E+00  2.7305E-01  1.0806E+00  1.7369E-01  1.5348E-02  1.1431E+00  6.0915E-01
             3.2783E+00
 PARAMETER:  1.4226E-01 -1.1594E+00 -1.2066E+00  2.7464E-01 -1.1981E+00  1.7756E-01 -1.6505E+00 -4.0768E+00  2.3378E-01 -3.9569E-01
             1.2873E+00
 GRADIENT:   9.5346E-02  1.6694E+00  5.0724E+00 -4.7365E-01 -7.9993E+00 -3.8643E-01 -2.5698E-05  1.5295E-03 -6.2797E-02  7.4295E-01
             1.2680E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1277.61795830733        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1164
 NPARAMETR:  1.0367E+00  1.8458E-01  3.0533E-01  1.2619E+00  2.8399E-01  1.0726E+00  1.0000E-02  1.0000E-02  1.0876E+00  6.1623E-01
             3.3002E+00
 PARAMETER:  1.3607E-01 -1.5897E+00 -1.0864E+00  3.3263E-01 -1.1588E+00  1.7010E-01 -5.0828E+00 -1.0368E+01  1.8394E-01 -3.8413E-01
             1.2940E+00
 GRADIENT:  -6.7180E-01  3.4251E-01 -1.1405E+00 -1.1732E+00  2.0903E+00 -1.3015E-01  0.0000E+00  0.0000E+00  3.1795E-01  8.9871E-01
             1.7471E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1277.80259749143        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1339
 NPARAMETR:  1.0304E+00  9.7826E-02  3.1869E-01  1.3039E+00  2.8489E-01  1.0668E+00  1.0000E-02  1.0000E-02  1.0524E+00  6.1717E-01
             3.3019E+00
 PARAMETER:  1.2990E-01 -2.2246E+00 -1.0435E+00  3.6539E-01 -1.1557E+00  1.6471E-01 -1.0189E+01 -1.9606E+01  1.5106E-01 -3.8261E-01
             1.2945E+00
 GRADIENT:   2.7877E-02  3.2997E-03 -6.5539E-02 -1.6256E+00  7.5920E-01 -2.6175E-01  0.0000E+00  0.0000E+00 -2.0650E-01 -7.8622E-02
            -6.1743E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1277.80715316205        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1501
 NPARAMETR:  1.0295E+00  8.6187E-02  3.1956E-01  1.3099E+00  2.8442E-01  1.0671E+00  1.0000E-02  1.0000E-02  1.0499E+00  6.1946E-01
             3.3016E+00
 PARAMETER:  1.2908E-01 -2.3512E+00 -1.0408E+00  3.6997E-01 -1.1573E+00  1.6493E-01 -1.1218E+01 -2.1455E+01  1.4866E-01 -3.7891E-01
             1.2944E+00
 GRADIENT:  -7.5778E-03  3.3176E-04  6.7997E-03 -7.6350E-04 -8.2575E-03 -4.8405E-03  0.0000E+00  0.0000E+00 -4.9521E-03  1.4123E-03
             1.1610E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1501
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.0227E-04 -1.5928E-05  9.3407E-05 -1.2157E-02  5.6078E-04
 SE:             2.8819E-02  1.0825E-05  2.2004E-04  2.6249E-02  1.9070E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8886E-01  1.4119E-01  6.7120E-01  6.4326E-01  9.7654E-01

 ETASHRINKSD(%)  3.4525E+00  9.9964E+01  9.9263E+01  1.2062E+01  3.6112E+01
 ETASHRINKVR(%)  6.7857E+00  1.0000E+02  9.9995E+01  2.2669E+01  5.9183E+01
 EBVSHRINKSD(%)  3.2415E+00  9.9964E+01  9.9255E+01  1.0787E+01  3.5748E+01
 EBVSHRINKVR(%)  6.3779E+00  1.0000E+02  9.9994E+01  2.0410E+01  5.8716E+01
 RELATIVEINF(%)  7.8752E+01  1.5550E-06  2.0631E-04  1.7302E+01  1.1610E+00
 EPSSHRINKSD(%)  2.9023E+01
 EPSSHRINKVR(%)  4.9622E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1277.8071531620544     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -542.65632659831624     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.99
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.13
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1277.807       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  8.62E-02  3.20E-01  1.31E+00  2.84E-01  1.07E+00  1.00E-02  1.00E-02  1.05E+00  6.19E-01  3.30E+00
 


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
+        8.65E+02
 
 TH 2
+       -9.49E+01  2.02E+02
 
 TH 3
+       -4.38E+01  5.74E+02  5.97E+03
 
 TH 4
+       -3.56E+01  1.80E+02 -3.35E+02  4.69E+02
 
 TH 5
+        2.02E+02 -1.26E+03 -8.58E+03 -3.09E+02  1.46E+04
 
 TH 6
+       -7.08E-01 -1.08E+01  2.29E+01 -1.04E+01  7.17E-01  1.52E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.66E+00 -3.51E+01  7.26E+01 -1.24E+01  5.00E+01  9.80E-01  0.00E+00  0.00E+00  1.15E+02
 
 TH10
+       -9.72E+00 -3.62E+00 -1.18E+02  1.03E+00  1.71E+02  1.37E+00  0.00E+00  0.00E+00 -1.20E-01  8.51E+01
 
 TH11
+       -1.33E+01 -2.25E+00 -1.17E+01 -6.44E+00  4.79E+00  2.77E+00  0.00E+00  0.00E+00  7.37E+00  2.17E+01  2.90E+01
 
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
 #CPUT: Total CPU Time in Seconds,       24.191
Stop Time:
Sat Sep 25 09:31:46 CDT 2021
