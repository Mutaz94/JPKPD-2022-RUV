Sat Sep 25 01:46:30 CDT 2021
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
$DATA ../../../../data/int/SL2/dat97.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -529.870850633489        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.0420E+02 -6.7347E+01  2.6396E+02  6.6888E+01  1.8244E+02 -1.9844E+01 -8.2625E+01 -2.8336E+02 -1.0141E+02 -1.9603E+01
            -6.1283E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2745.32882737270        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.3812E-01  1.3283E+00  9.1283E-01  8.6647E-01  1.0674E+00  1.0226E+00  9.7564E-01  9.6167E-01  9.6264E-01  8.3673E-01
             2.6777E+00
 PARAMETER:  3.6118E-02  3.8390E-01  8.7934E-03 -4.3325E-02  1.6520E-01  1.2236E-01  7.5335E-02  6.0911E-02  6.1925E-02 -7.8254E-02
             1.0849E+00
 GRADIENT:  -1.0070E+01  2.6655E+01 -1.3086E+00  1.4357E+01 -3.4698E+01  1.9827E+00  5.8813E+00  7.7784E+00 -8.7210E+00 -3.2488E+00
             6.5872E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2750.28251750471        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.4344E-01  1.5576E+00  1.0859E+00  7.6939E-01  1.2609E+00  1.0656E+00  8.9245E-01  4.6101E-01  9.9458E-01  1.0874E+00
             2.7178E+00
 PARAMETER:  4.1776E-02  5.4317E-01  1.8238E-01 -1.6216E-01  3.3179E-01  1.6354E-01 -1.3783E-02 -6.7433E-01  9.4563E-02  1.8383E-01
             1.0998E+00
 GRADIENT:  -2.1473E-01  7.9682E+01  3.8330E+00  5.9022E+01 -3.2709E+01  1.6693E+01  6.0698E+00  9.7325E-01 -1.0810E+01  1.1721E+01
             4.7158E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2761.31259029282        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.4304E-01  1.7551E+00  1.2979E+00  5.8899E-01  1.5431E+00  1.0105E+00  7.0915E-01  2.6434E-01  1.4003E+00  1.1806E+00
             2.6179E+00
 PARAMETER:  4.1353E-02  6.6250E-01  3.6071E-01 -4.2935E-01  5.3381E-01  1.1041E-01 -2.4369E-01 -1.2305E+00  4.3671E-01  2.6600E-01
             1.0624E+00
 GRADIENT:   3.4069E+00 -1.6982E+01 -7.1527E+00  1.4506E+01  2.5080E+01 -2.6637E+00 -1.8328E+00  1.1581E-01 -1.0582E+00 -2.0522E-01
            -1.2039E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2767.05869603931        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      326
 NPARAMETR:  9.4273E-01  2.0458E+00  1.8931E+00  4.0149E-01  1.8059E+00  1.0153E+00  6.4139E-01  1.1341E-01  1.8509E+00  1.3796E+00
             2.6080E+00
 PARAMETER:  4.1024E-02  8.1577E-01  7.3821E-01 -8.1258E-01  6.9104E-01  1.1518E-01 -3.4412E-01 -2.0768E+00  7.1570E-01  4.2176E-01
             1.0586E+00
 GRADIENT:  -3.0115E+00 -8.6625E+00 -1.8383E+00  4.1126E+00  8.9613E+00 -2.1905E+00 -2.7761E+00 -4.1272E-03 -3.3079E-01  2.8022E+00
            -5.0604E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2767.97172055624        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      501
 NPARAMETR:  9.4438E-01  2.2470E+00  2.6049E+00  2.8061E-01  1.9339E+00  1.0210E+00  6.9893E-01  3.8669E-02  1.9250E+00  1.4568E+00
             2.6136E+00
 PARAMETER:  4.2772E-02  9.0961E-01  1.0574E+00 -1.1708E+00  7.5954E-01  1.2081E-01 -2.5821E-01 -3.1527E+00  7.5495E-01  4.7621E-01
             1.0607E+00
 GRADIENT:   2.1947E-01  8.5476E+00 -3.4517E-01  2.8368E+00  2.3282E+00 -2.7385E-01  2.7120E+00 -1.4241E-04  4.8284E-02 -3.3581E-01
            -2.5710E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2768.48602785203        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      676
 NPARAMETR:  9.4437E-01  2.4372E+00  4.4248E+00  1.5076E-01  2.0019E+00  1.0254E+00  6.7684E-01  1.0000E-02  2.5516E+00  1.5037E+00
             2.6126E+00
 PARAMETER:  4.2765E-02  9.9084E-01  1.5872E+00 -1.7921E+00  7.9408E-01  1.2511E-01 -2.9032E-01 -5.5229E+00  1.0367E+00  5.0795E-01
             1.0603E+00
 GRADIENT:  -7.0279E-02  6.8197E+00 -2.1706E-01  1.2753E+00  5.9759E+00  1.2126E+00 -5.1172E-01  0.0000E+00  3.9678E-01  2.3695E+00
             1.0190E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2768.68893695989        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      852
 NPARAMETR:  9.4439E-01  2.5218E+00  7.0861E+00  8.9593E-02  1.9986E+00  1.0228E+00  6.7011E-01  1.0000E-02  3.2537E+00  1.4951E+00
             2.6112E+00
 PARAMETER:  4.2785E-02  1.0250E+00  2.0581E+00 -2.3125E+00  7.9246E-01  1.2256E-01 -3.0031E-01 -7.6940E+00  1.2798E+00  5.0218E-01
             1.0598E+00
 GRADIENT:   6.5710E-02  2.7803E+00 -5.2545E-02  2.3551E-01  1.8379E+00  3.1235E-01  4.2335E-01  0.0000E+00 -2.4240E-02  5.2836E-01
             1.6952E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2768.73984197408        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1027
 NPARAMETR:  9.4435E-01  2.5865E+00  1.4300E+01  4.4476E-02  1.9965E+00  1.0221E+00  6.6333E-01  1.0000E-02  4.3492E+00  1.4926E+00
             2.6110E+00
 PARAMETER:  4.2744E-02  1.0503E+00  2.7603E+00 -3.0128E+00  7.9138E-01  1.2181E-01 -3.1048E-01 -1.0720E+01  1.5700E+00  5.0054E-01
             1.0597E+00
 GRADIENT:  -2.9933E-02  3.2433E+00 -7.9955E-03  9.0399E-02  2.3511E-01  3.9989E-02 -1.3354E-01  0.0000E+00 -1.8447E-02  7.5726E-02
            -8.7702E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2768.74777649714        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1202
 NPARAMETR:  9.4434E-01  2.6142E+00  2.7404E+01  2.4172E-02  1.9956E+00  1.0218E+00  6.6079E-01  1.0000E-02  5.4982E+00  1.4919E+00
             2.6110E+00
 PARAMETER:  4.2736E-02  1.0610E+00  3.4107E+00 -3.6226E+00  7.9094E-01  1.2157E-01 -3.1432E-01 -1.3408E+01  1.8044E+00  5.0008E-01
             1.0597E+00
 GRADIENT:   1.8538E-03  9.9400E-01 -1.3995E-03  1.7434E-02 -1.2085E-01 -3.2964E-02  1.9923E-02  0.0000E+00  7.0823E-04 -3.8784E-02
            -2.3712E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2768.74927531162        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1377
 NPARAMETR:  9.4438E-01  2.6306E+00  5.6743E+01  1.2441E-02  1.9957E+00  1.0217E+00  6.5876E-01  1.0000E-02  7.0412E+00  1.4920E+00
             2.6109E+00
 PARAMETER:  4.2774E-02  1.0672E+00  4.1385E+00 -4.2867E+00  7.9099E-01  1.2150E-01 -3.1740E-01 -1.6360E+01  2.0518E+00  5.0013E-01
             1.0597E+00
 GRADIENT:   7.4114E-02  3.2554E-01 -1.9835E-04  4.9868E-03 -8.0142E-02 -6.0952E-02 -1.2135E-02  0.0000E+00  3.5456E-03 -2.6883E-02
            -3.8034E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2768.74946857809        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1561             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4435E-01  2.6338E+00  7.2438E+01  1.0000E-02  1.9959E+00  1.0219E+00  6.5834E-01  1.0000E-02  7.5898E+00  1.4922E+00
             2.6110E+00
 PARAMETER:  4.2743E-02  1.0684E+00  4.3827E+00 -4.5148E+00  7.9108E-01  1.2162E-01 -3.1803E-01 -1.7316E+01  2.1268E+00  5.0025E-01
             1.0597E+00
 GRADIENT:   6.0945E+00  5.2299E+01 -1.0291E-04  0.0000E+00  3.4479E+00  1.2041E+00  6.6343E-01  0.0000E+00  1.8938E-02  4.4566E-01
             1.9337E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2768.75427168573        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1744
 NPARAMETR:  9.4447E-01  2.6347E+00  1.0243E+02  1.0000E-02  1.9962E+00  1.0238E+00  6.5997E-01  1.0000E-02  1.8316E+00  1.4923E+00
             2.6108E+00
 PARAMETER:  4.2872E-02  1.0688E+00  4.7291E+00 -4.5148E+00  7.9122E-01  1.2354E-01 -3.1556E-01 -1.7316E+01  7.0522E-01  5.0033E-01
             1.0597E+00
 GRADIENT:   2.8671E-01  1.6722E+00 -1.0491E-04  0.0000E+00  4.3883E-02  7.1273E-01  1.7484E-01  0.0000E+00  1.4775E-03  1.6792E-02
            -1.2809E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2768.75620854968        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1919
 NPARAMETR:  9.4440E-01  2.6337E+00  1.5503E+02  1.0000E-02  1.9959E+00  1.0219E+00  6.5975E-01  1.0000E-02  3.4996E-01  1.4921E+00
             2.6109E+00
 PARAMETER:  4.2790E-02  1.0684E+00  5.1436E+00 -4.5148E+00  7.9107E-01  1.2165E-01 -3.1589E-01 -1.7316E+01 -9.4995E-01  5.0021E-01
             1.0597E+00
 GRADIENT:   7.1229E-02  3.6379E-02 -7.0468E-05  0.0000E+00 -6.8064E-03 -2.5273E-03  2.7457E-02  0.0000E+00  6.8620E-05 -6.1023E-03
            -2.0849E-02

0ITERATION NO.:   69    OBJECTIVE VALUE:  -2768.75624628378        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     2054
 NPARAMETR:  9.4433E-01  2.6336E+00  1.7975E+02  1.0000E-02  1.9959E+00  1.0219E+00  6.5970E-01  1.0000E-02  1.9473E-01  1.4922E+00
             2.6110E+00
 PARAMETER:  4.2798E-02  1.0684E+00  5.2943E+00 -4.5148E+00  7.9108E-01  1.2165E-01 -3.1599E-01 -1.7316E+01 -1.5346E+00  5.0023E-01
             1.0597E+00
 GRADIENT:   1.0741E-01  4.2544E-02  5.4900E-04  0.0000E+00 -1.9877E-03 -1.5161E-02 -6.9589E-03  0.0000E+00  2.3628E-03 -4.2924E-03
            -1.2615E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2054
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2772E-03 -7.8268E-03  3.3723E-10 -9.9670E-05 -1.2176E-02
 SE:             2.9329E-02  2.8357E-02  2.2104E-10  4.9537E-05  2.7001E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6526E-01  7.8254E-01  1.2710E-01  4.4219E-02  6.5203E-01

 ETASHRINKSD(%)  1.7451E+00  5.0013E+00  1.0000E+02  9.9834E+01  9.5433E+00
 ETASHRINKVR(%)  3.4598E+00  9.7524E+00  1.0000E+02  1.0000E+02  1.8176E+01
 EBVSHRINKSD(%)  1.6577E+00  5.1313E+00  1.0000E+02  9.9842E+01  8.4911E+00
 EBVSHRINKVR(%)  3.2878E+00  9.9992E+00  1.0000E+02  1.0000E+02  1.6261E+01
 RELATIVEINF(%)  9.6662E+01  3.0703E+00  0.0000E+00  8.5598E-06  8.2401E+01
 EPSSHRINKSD(%)  1.5472E+01
 EPSSHRINKVR(%)  2.8551E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2768.7562462837832     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1114.6668865153724     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    47.77
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.26
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2768.756       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.44E-01  2.63E+00  1.80E+02  1.00E-02  2.00E+00  1.02E+00  6.60E-01  1.00E-02  1.95E-01  1.49E+00  2.61E+00
 


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
+        1.39E+03
 
 TH 2
+       -1.65E+01  3.32E+02
 
 TH 3
+       -3.62E-03  7.36E-04  7.92E-07
 
 TH 4
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 5
+        6.59E+00 -1.46E+01  1.16E-04  0.00E+00  1.05E+02
 
 TH 6
+        3.05E+01 -6.90E+00  1.25E-03  0.00E+00 -7.27E+00  2.71E+02
 
 TH 7
+       -7.21E+01 -1.30E+00 -1.31E-02  0.00E+00  2.34E-01  1.31E+01  3.78E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.63E+01 -7.11E-01  2.61E-03  0.00E+00  1.36E+00 -3.26E+01 -6.22E+00  0.00E+00  7.69E-01
 
 TH10
+       -9.05E+00 -2.18E-01 -5.66E-04  0.00E+00 -5.94E+00 -1.58E+01 -3.25E+00  0.00E+00 -3.00E+00  6.06E+01
 
 TH11
+       -1.33E+01 -1.21E+01 -5.70E-04  0.00E+00  1.82E+00 -1.46E-02  1.09E+01  0.00E+00 -2.51E+00  9.14E+00  1.76E+02
 
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
 #CPUT: Total CPU Time in Seconds,       59.124
Stop Time:
Sat Sep 25 01:47:31 CDT 2021
