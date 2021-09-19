Sat Sep 18 10:20:57 CDT 2021
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
$DATA ../../../../data/spa/A3/dat24.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m24.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -487.438615163345        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3318E+02  1.7212E+00  1.4099E+02 -1.9712E+02  8.2852E+01  1.3819E+01 -3.8059E+01 -6.1481E+01 -1.6106E+02 -1.0863E+02
            -1.9123E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1265.12775106563        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9399E-01  9.5016E-01  8.7934E-01  1.2267E+00  9.0490E-01  8.1575E-01  9.7507E-01  1.0045E+00  1.2130E+00  9.8410E-01
             3.7041E+00
 PARAMETER:  9.3976E-02  4.8875E-02 -2.8585E-02  3.0433E-01  6.8697E-05 -1.0365E-01  7.4758E-02  1.0447E-01  2.9312E-01  8.3968E-02
             1.4094E+00
 GRADIENT:   9.2577E+00  1.5220E+01 -5.3671E+00  3.5107E+01 -9.9020E+00 -2.8927E+01  7.4301E+00  7.3016E+00  1.5808E+01  1.6426E+01
             7.0852E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1274.19848395320        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  9.9446E-01  7.0582E-01  4.3483E-01  1.3399E+00  4.9015E-01  8.7440E-01  6.9669E-01  4.7771E-01  1.2131E+00  5.2487E-01
             3.4723E+00
 PARAMETER:  9.4449E-02 -2.4840E-01 -7.3279E-01  3.9256E-01 -6.1305E-01 -3.4219E-02 -2.6142E-01 -6.3874E-01  2.9320E-01 -5.4460E-01
             1.3448E+00
 GRADIENT:   3.0179E+00  4.5552E+01  9.5420E+00  1.0022E+02 -3.2368E+01 -1.5960E+01 -1.1659E+00  1.6159E+00  1.2549E+01  2.9390E+00
             4.8423E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1297.76608464564        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      234
 NPARAMETR:  9.6608E-01  5.9821E-01  1.8312E-01  1.1279E+00  2.9228E-01  1.0060E+00  1.2207E+00  1.3678E-01  1.1317E+00  2.8057E-01
             2.7032E+00
 PARAMETER:  6.5496E-02 -4.1382E-01 -1.5976E+00  2.2038E-01 -1.1300E+00  1.0594E-01  2.9944E-01 -1.8894E+00  2.2375E-01 -1.1709E+00
             1.0944E+00
 GRADIENT:  -1.3773E+00  4.5361E+01  1.2338E+01  5.8110E+01 -4.1898E+01  2.0442E+01 -2.1629E+00 -4.7743E-01 -2.5175E+01 -3.4130E+00
            -5.2948E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1304.72296113563        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  9.6014E-01  3.6680E-01  1.3479E-01  1.0826E+00  2.0162E-01  9.1341E-01  1.2018E+00  1.0845E-02  1.4168E+00  6.9408E-01
             2.4855E+00
 PARAMETER:  5.9329E-02 -9.0294E-01 -1.9040E+00  1.7933E-01 -1.5013E+00  9.4337E-03  2.8385E-01 -4.4241E+00  4.4838E-01 -2.6517E-01
             1.0105E+00
 GRADIENT:  -1.5455E+01  3.0078E+01  6.0698E+00  2.3931E+01 -2.2148E+01 -1.6439E+01  1.1073E+01 -1.3773E-03 -1.7567E+00  8.3287E+00
            -6.6697E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1305.64923113337        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  9.6770E-01  2.9100E-01  1.1123E-01  1.0238E+00  1.6959E-01  9.3172E-01  1.0583E+00  1.0000E-02  1.5661E+00  7.5848E-01
             2.5093E+00
 PARAMETER:  6.7163E-02 -1.1344E+00 -2.0962E+00  1.2348E-01 -1.6744E+00  2.9275E-02  1.5671E-01 -5.6729E+00  5.4856E-01 -1.7644E-01
             1.0200E+00
 GRADIENT:   1.7761E-01  2.3439E+01  1.4983E+01  2.1635E+01 -4.5292E+01 -6.1543E+00  6.3288E+00  0.0000E+00 -2.4769E+00  6.9968E+00
             1.0416E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1307.84193768042        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  9.7040E-01  3.0626E-01  1.2397E-01  1.0304E+00  1.8576E-01  9.5310E-01  9.6519E-01  1.0000E-02  1.4894E+00  7.1579E-01
             2.5140E+00
 PARAMETER:  6.9950E-02 -1.0833E+00 -1.9877E+00  1.2998E-01 -1.5833E+00  5.1966E-02  6.4565E-02 -5.2286E+00  4.9840E-01 -2.3437E-01
             1.0219E+00
 GRADIENT:   1.8151E+00 -3.9457E+00 -2.9936E+00 -1.9959E+00  9.3423E+00  1.4528E+00  2.9584E+00  0.0000E+00  1.7594E+00  8.1970E-02
             6.1835E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1309.00993049048        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      709
 NPARAMETR:  9.6749E-01  3.2702E-01  1.2203E-01  1.0248E+00  1.8767E-01  9.4895E-01  3.6425E-01  1.0612E-02  1.4949E+00  7.6021E-01
             2.5077E+00
 PARAMETER:  6.6950E-02 -1.0177E+00 -2.0035E+00  1.2453E-01 -1.5731E+00  4.7599E-02 -9.0992E-01 -4.4458E+00  5.0205E-01 -1.7416E-01
             1.0194E+00
 GRADIENT:  -1.8545E-01  9.3705E-01  2.1502E+00  5.3791E-01 -1.9170E+00  2.0497E-01  2.5708E-01 -9.4362E-04  3.6214E-01 -8.7074E-01
             9.6364E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1309.16876209934        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      884
 NPARAMETR:  9.6696E-01  3.2298E-01  1.1879E-01  1.0159E+00  1.8485E-01  9.4788E-01  5.5269E-02  3.9937E-02  1.5130E+00  7.7517E-01
             2.5048E+00
 PARAMETER:  6.6397E-02 -1.0302E+00 -2.0304E+00  1.1578E-01 -1.5882E+00  4.6469E-02 -2.7955E+00 -3.1205E+00  5.1410E-01 -1.5468E-01
             1.0182E+00
 GRADIENT:   8.5754E-02  4.0271E-01  9.7849E-02  5.1376E-01 -9.3790E-01  6.1623E-02  6.0698E-03 -1.2147E-02 -5.3788E-02  4.6041E-01
             2.2430E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1310.50085403622        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1067
 NPARAMETR:  9.6436E-01  3.2447E-01  1.1674E-01  1.0173E+00  1.8367E-01  9.3934E-01  1.0000E-02  9.6088E-01  1.5306E+00  7.1379E-01
             2.4343E+00
 PARAMETER:  6.3709E-02 -1.0256E+00 -2.0478E+00  1.1714E-01 -1.5946E+00  3.7422E-02 -6.4509E+00  6.0098E-02  5.2567E-01 -2.3717E-01
             9.8965E-01
 GRADIENT:   5.7746E-01  1.6350E+00  2.1221E+00  8.0169E-01 -2.2300E+00 -3.1450E-01  0.0000E+00  1.1481E+00  1.2875E+00  4.0377E+00
             4.9057E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1310.59074322266        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1242
 NPARAMETR:  9.6375E-01  3.1823E-01  1.1490E-01  1.0124E+00  1.8150E-01  9.4032E-01  1.0000E-02  9.7304E-01  1.5390E+00  6.9428E-01
             2.4114E+00
 PARAMETER:  6.3077E-02 -1.0450E+00 -2.0637E+00  1.1228E-01 -1.6065E+00  3.8461E-02 -6.4763E+00  7.2672E-02  5.3115E-01 -2.6488E-01
             9.8021E-01
 GRADIENT:   1.6663E-02  6.7867E-03  1.1922E-02 -8.0606E-03 -3.0139E-02  1.8987E-04  0.0000E+00 -8.9640E-04 -4.4845E-03  1.2251E-03
            -8.7350E-03

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1310.59074393596        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1299
 NPARAMETR:  9.6374E-01  3.1824E-01  1.1490E-01  1.0124E+00  1.8151E-01  9.4031E-01  1.0000E-02  9.7315E-01  1.5391E+00  6.9425E-01
             2.4114E+00
 PARAMETER:  6.3069E-02 -1.0449E+00 -2.0637E+00  1.1230E-01 -1.6065E+00  3.8459E-02 -6.4762E+00  7.2786E-02  5.3117E-01 -2.6493E-01
             9.8022E-01
 GRADIENT:   2.7279E-03  1.1917E-03  1.0686E-03 -5.1584E-04 -5.3515E-03 -1.9827E-04  0.0000E+00 -4.6120E-04 -8.2295E-04 -5.8594E-04
            -2.5235E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1299
 NO. OF SIG. DIGITS IN FINAL EST.:  3.8
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6380E-03 -1.3629E-04  2.9090E-02 -8.7700E-03  2.0520E-02
 SE:             2.8927E-02  1.3357E-04  1.4709E-02  2.6118E-02  2.3172E-02
 N:                     100         100         100         100         100

 P VAL.:         9.2734E-01  3.0753E-01  4.7955E-02  7.3704E-01  3.7587E-01

 ETASHRINKSD(%)  3.0904E+00  9.9553E+01  5.0724E+01  1.2501E+01  2.2370E+01
 ETASHRINKVR(%)  6.0853E+00  9.9998E+01  7.5719E+01  2.3439E+01  3.9736E+01
 EBVSHRINKSD(%)  2.9102E+00  9.9551E+01  5.0381E+01  9.4445E+00  2.2767E+01
 EBVSHRINKVR(%)  5.7356E+00  9.9998E+01  7.5380E+01  1.7997E+01  4.0350E+01
 RELATIVEINF(%)  9.0557E+01  2.6666E-04  4.1759E+00  4.8731E+01  3.5444E+00
 EPSSHRINKSD(%)  3.8217E+01
 EPSSHRINKVR(%)  6.1828E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1310.5907439359569     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -575.43991737221870     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.65
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1310.591       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.64E-01  3.18E-01  1.15E-01  1.01E+00  1.82E-01  9.40E-01  1.00E-02  9.73E-01  1.54E+00  6.94E-01  2.41E+00
 


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
+        1.28E+03
 
 TH 2
+        1.06E+02  2.37E+03
 
 TH 3
+       -4.61E+02  4.21E+03  2.21E+04
 
 TH 4
+       -4.27E-02  5.04E+01 -6.30E+02  3.77E+02
 
 TH 5
+        1.91E+02 -8.91E+03 -2.39E+04 -4.22E+02  4.25E+04
 
 TH 6
+        2.75E+00 -1.80E+01  1.73E+01 -6.59E+00  2.69E+01  2.01E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        2.34E+00 -7.22E+00  5.68E+01 -9.43E+00  1.81E+01  4.27E+00  0.00E+00  1.05E+01
 
 TH 9
+        1.35E+01 -2.97E+01  2.58E+02 -1.12E+01  1.53E+02 -1.00E+00  0.00E+00  7.31E-01  5.11E+01
 
 TH10
+       -1.95E-01 -5.67E+01  4.06E+01  6.73E+00  2.01E+02  3.22E+00  0.00E+00  2.75E+01  8.84E+00  1.45E+02
 
 TH11
+       -2.21E+01 -1.09E+01 -8.85E-01 -1.30E-01  4.50E+01  2.28E+00  0.00E+00  7.74E+00  5.40E+00  1.24E+01  4.08E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.065
Stop Time:
Sat Sep 18 10:21:20 CDT 2021
