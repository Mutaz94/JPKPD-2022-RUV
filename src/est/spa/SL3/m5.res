Sat Sep 18 12:37:50 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat5.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m5.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1620.70682876488        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5728E+02 -4.9845E+01 -3.1274E+01 -2.3586E+01  4.1255E+01  1.3316E+01 -2.3647E+00  9.9845E+00  1.8576E+01 -2.2721E+01
            -1.4436E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1632.11155657246        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.4609E-01  1.1216E+00  1.2111E+00  9.4772E-01  1.1745E+00  9.4534E-01  1.0561E+00  8.3694E-01  8.0691E-01  1.3475E+00
             1.0406E+00
 PARAMETER:  4.4586E-02  2.1473E-01  2.9156E-01  4.6300E-02  2.6081E-01  4.3785E-02  1.5461E-01 -7.8002E-02 -1.1454E-01  3.9826E-01
             1.3982E-01
 GRADIENT:   3.3585E+01 -1.4069E+01 -9.3394E+00 -1.1163E+01  1.2817E+01 -1.0922E+00 -8.4972E+00  2.8699E+00 -1.0563E+01  1.3468E+00
            -7.2658E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1633.47375684236        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  9.4896E-01  9.7899E-01  1.4626E+00  1.0688E+00  1.1954E+00  9.3235E-01  1.0535E+00  6.5975E-01  8.8840E-01  1.3987E+00
             1.0765E+00
 PARAMETER:  4.7608E-02  7.8771E-02  4.8020E-01  1.6653E-01  2.7852E-01  2.9953E-02  1.5210E-01 -3.1589E-01 -1.8336E-02  4.3557E-01
             1.7368E-01
 GRADIENT:   4.2239E+01  1.0754E+01 -5.9179E+00  2.8464E+01  1.1633E+01 -6.6587E+00 -1.1533E-01  9.7949E-01  1.8058E-01  4.2465E-01
             9.3163E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1633.77417310244        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.3761E-01  1.0074E+00  1.2788E+00  1.0338E+00  1.1284E+00  9.4318E-01  1.0708E+00  4.5104E-01  8.8712E-01  1.3289E+00
             1.0529E+00
 PARAMETER:  3.5579E-02  1.0737E-01  3.4591E-01  1.3319E-01  2.2080E-01  4.1507E-02  1.6838E-01 -6.9620E-01 -1.9773E-02  3.8438E-01
             1.5154E-01
 GRADIENT:   1.2421E+01  1.4960E+00 -3.4384E+00  7.2798E+00  4.0082E+00 -2.2292E+00 -4.5562E-01  4.8940E-01  1.4190E-01 -4.3927E-01
             2.3523E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1633.78122540868        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.3475E-01  1.0089E+00  1.2847E+00  1.0308E+00  1.1301E+00  9.4652E-01  1.0678E+00  3.9024E-01  8.8986E-01  1.3361E+00
             1.0510E+00
 PARAMETER:  3.2520E-02  1.0887E-01  3.5056E-01  1.3029E-01  2.2230E-01  4.5039E-02  1.6556E-01 -8.4100E-01 -1.6696E-02  3.8975E-01
             1.4973E-01
 GRADIENT:   5.3939E+00  3.6586E-01 -1.8008E+00  3.1339E+00  1.9123E+00 -9.3700E-01 -2.6404E-01  3.2180E-01  1.0660E-01 -2.4443E-01
             1.1130E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1633.78618151840        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  9.3293E-01  1.0105E+00  1.2888E+00  1.0285E+00  1.1317E+00  9.4869E-01  1.0637E+00  3.0243E-01  8.9283E-01  1.3431E+00
             1.0502E+00
 PARAMETER:  3.0575E-02  1.1040E-01  3.5368E-01  1.2810E-01  2.2374E-01  4.7330E-02  1.6177E-01 -1.0959E+00 -1.3359E-02  3.9498E-01
             1.4900E-01
 GRADIENT:   9.5444E-01 -2.4486E-01 -6.5136E-01  5.1053E-01  5.0710E-01 -1.3149E-01 -1.1802E-01  1.7073E-01  7.0603E-02 -4.8221E-02
             3.1708E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1633.80091574435        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  9.3186E-01  1.0119E+00  1.2903E+00  1.0268E+00  1.1328E+00  9.5000E-01  1.0592E+00  1.6116E-01  8.9585E-01  1.3491E+00
             1.0502E+00
 PARAMETER:  2.9431E-02  1.1185E-01  3.5487E-01  1.2643E-01  2.2470E-01  4.8712E-02  1.5748E-01 -1.7253E+00 -9.9868E-03  3.9941E-01
             1.4897E-01
 GRADIENT:  -1.6693E+00 -5.2408E-01  1.1654E-01 -1.0380E+00 -3.8345E-01  3.3666E-01 -1.5954E-02  4.3865E-02  4.3528E-02  1.0960E-01
            -1.6078E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1634.19245106877        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      586
 NPARAMETR:  9.4631E-01  9.5890E-01  1.4034E+00  1.0707E+00  1.1529E+00  9.5611E-01  1.0668E+00  1.0000E-02  8.9599E-01  1.3908E+00
             1.0557E+00
 PARAMETER:  4.4811E-02  5.8031E-02  4.3891E-01  1.6835E-01  2.4224E-01  5.5123E-02  1.6466E-01 -1.3996E+01 -9.8237E-03  4.2985E-01
             1.5425E-01
 GRADIENT:   2.6411E+00  3.4503E+00  7.7377E-01  3.2505E+00 -1.7217E+00  6.4640E-01  6.7244E-01  0.0000E+00  1.2301E+00  4.7700E-01
            -1.7146E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1634.35871079808        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      763
 NPARAMETR:  9.4301E-01  7.6158E-01  1.5538E+00  1.1993E+00  1.1270E+00  9.5254E-01  1.1212E+00  1.0000E-02  8.4728E-01  1.4013E+00
             1.0592E+00
 PARAMETER:  4.1319E-02 -1.7236E-01  5.4073E-01  2.8177E-01  2.1953E-01  5.1372E-02  2.1437E-01 -8.2257E+00 -6.5725E-02  4.3737E-01
             1.5750E-01
 GRADIENT:  -1.7236E+00  2.2558E+00  1.0435E-01  4.2970E+00 -4.1579E-01 -2.4481E-01 -2.8827E-01  0.0000E+00 -3.5261E-01 -2.5894E-01
            -3.1207E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1634.40524739745        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      938
 NPARAMETR:  9.4334E-01  6.4438E-01  1.6287E+00  1.2739E+00  1.1105E+00  9.5083E-01  1.1973E+00  1.3348E-01  8.1870E-01  1.4085E+00
             1.0614E+00
 PARAMETER:  4.1668E-02 -3.3946E-01  5.8776E-01  3.4206E-01  2.0478E-01  4.9581E-02  2.8004E-01 -1.9138E+00 -1.0003E-01  4.4251E-01
             1.5957E-01
 GRADIENT:   1.8083E+00  5.7717E-01 -4.9206E-01  1.9914E+00  6.3385E-01 -5.7428E-01 -1.1542E-01  4.7945E-03  3.7692E-03 -2.2498E-02
             3.7079E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1634.41588134743        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1116
 NPARAMETR:  9.4221E-01  5.9243E-01  1.6494E+00  1.3064E+00  1.0985E+00  9.5179E-01  1.2877E+00  8.8124E-02  8.0062E-01  1.4071E+00
             1.0608E+00
 PARAMETER:  4.0473E-02 -4.2352E-01  6.0039E-01  3.6726E-01  1.9394E-01  5.0587E-02  3.5286E-01 -2.3290E+00 -1.2237E-01  4.4154E-01
             1.5901E-01
 GRADIENT:   1.5974E-01  1.3975E-01  6.2643E-03  3.2059E-01 -6.2698E-04 -1.9239E-02  6.6778E-02  2.0424E-04  7.6017E-02 -1.3377E-02
             5.2883E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1634.41660854021        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1292
 NPARAMETR:  9.4191E-01  5.6132E-01  1.6821E+00  1.3267E+00  1.0982E+00  9.5155E-01  1.2851E+00  2.2248E-01  7.9630E-01  1.4121E+00
             1.0612E+00
 PARAMETER:  4.0154E-02 -4.7746E-01  6.2007E-01  3.8272E-01  1.9366E-01  5.0336E-02  3.5082E-01 -1.4029E+00 -1.2778E-01  4.4510E-01
             1.5937E-01
 GRADIENT:   1.8617E-01 -3.9933E-04 -1.6568E-01  1.2689E-01  1.3290E-01 -3.5449E-03 -3.5905E-02  6.8564E-04 -4.1733E-02  1.3263E-01
             1.3952E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1634.41713283398        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1470
 NPARAMETR:  9.4137E-01  5.3454E-01  1.7121E+00  1.3448E+00  1.0971E+00  9.5151E-01  1.3120E+00  3.4055E-01  7.8976E-01  1.4119E+00
             1.0604E+00
 PARAMETER:  3.9580E-02 -5.2636E-01  6.3770E-01  3.9625E-01  1.9266E-01  5.0299E-02  3.7158E-01 -9.7719E-01 -1.3602E-01  4.4493E-01
             1.5861E-01
 GRADIENT:  -4.9302E-01  3.0279E-01  2.0962E-01  7.8502E-01 -4.7014E-01  7.6004E-02 -2.4135E-02 -6.2303E-04 -3.9234E-02 -1.6453E-02
            -1.8296E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1634.41856694345        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1649
 NPARAMETR:  9.4087E-01  4.8633E-01  1.7608E+00  1.3769E+00  1.0952E+00  9.5089E-01  1.3846E+00  4.7425E-01  7.7673E-01  1.4120E+00
             1.0600E+00
 PARAMETER:  3.9051E-02 -6.2087E-01  6.6576E-01  4.1987E-01  1.9093E-01  4.9647E-02  4.2540E-01 -6.4603E-01 -1.5266E-01  4.4498E-01
             1.5829E-01
 GRADIENT:  -5.2762E-01  4.4569E-01  2.4408E-01  1.3882E+00 -4.1498E-01 -1.1636E-02  1.4574E-02  1.2668E-04 -5.0313E-02 -1.0383E-01
            -1.5510E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1634.42311494138        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1827
 NPARAMETR:  9.4038E-01  4.0747E-01  1.8234E+00  1.4280E+00  1.0880E+00  9.5037E-01  1.5053E+00  5.8925E-01  7.5883E-01  1.4142E+00
             1.0599E+00
 PARAMETER:  3.8529E-02 -7.9779E-01  7.0068E-01  4.5625E-01  1.8433E-01  4.9097E-02  5.0896E-01 -4.2891E-01 -1.7597E-01  4.4653E-01
             1.5819E-01
 GRADIENT:   2.6335E-01  1.1759E-01 -9.5318E-02  4.8586E-01 -1.5343E-01  4.2666E-02  2.2679E-02  3.9873E-03  3.9816E-02  2.0652E-01
             1.2453E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1634.42813499739        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2006
 NPARAMETR:  9.3959E-01  3.3092E-01  1.9243E+00  1.4794E+00  1.0911E+00  9.4961E-01  1.5749E+00  7.3350E-01  7.4510E-01  1.4188E+00
             1.0598E+00
 PARAMETER:  3.7691E-02 -1.0059E+00  7.5458E-01  4.9161E-01  1.8723E-01  4.8301E-02  5.5418E-01 -2.0992E-01 -1.9423E-01  4.4983E-01
             1.5804E-01
 GRADIENT:   1.9058E-01  2.1555E-01  6.2514E-03  1.1358E+00  1.2472E-03 -1.4522E-02  4.0332E-03  1.6048E-02 -3.4255E-02 -5.7092E-02
            -2.1785E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1634.43222686548        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2184
 NPARAMETR:  9.3907E-01  2.6805E-01  2.0015E+00  1.5218E+00  1.0913E+00  9.4911E-01  1.5681E+00  8.0748E-01  7.3576E-01  1.4246E+00
             1.0602E+00
 PARAMETER:  3.7136E-02 -1.2166E+00  7.9389E-01  5.1992E-01  1.8740E-01  4.7768E-02  5.4990E-01 -1.1384E-01 -2.0685E-01  4.5390E-01
             1.5846E-01
 GRADIENT:   3.8498E-01  4.5491E-01  1.8424E-01  3.3169E+00 -4.8530E-01 -2.6209E-02  1.5650E-02 -8.2969E-03 -4.2811E-02 -6.7657E-04
            -4.8275E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1634.45257254362        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2364
 NPARAMETR:  9.3793E-01  1.7260E-01  2.1067E+00  1.5852E+00  1.0901E+00  9.4833E-01  1.3358E+00  9.1059E-01  7.2137E-01  1.4291E+00
             1.0601E+00
 PARAMETER:  3.5915E-02 -1.6568E+00  8.4511E-01  5.6072E-01  1.8623E-01  4.6946E-02  3.8952E-01  6.3326E-03 -2.2660E-01  4.5703E-01
             1.5838E-01
 GRADIENT:  -2.0821E-01  3.7513E-01 -2.1823E-01  5.2244E+00  1.5110E-01 -6.2050E-02  5.2216E-02 -1.2427E-02  3.1633E-01 -8.7136E-02
            -1.0129E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1634.47113942735        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2541
 NPARAMETR:  9.3759E-01  1.1866E-01  2.1756E+00  1.6197E+00  1.0888E+00  9.4812E-01  9.9474E-01  9.7263E-01  7.1028E-01  1.4304E+00
             1.0601E+00
 PARAMETER:  3.5558E-02 -2.0315E+00  8.7732E-01  5.8223E-01  1.8508E-01  4.6722E-02  9.4724E-02  7.2251E-02 -2.4210E-01  4.5797E-01
             1.5837E-01
 GRADIENT:   2.7595E-01  9.5992E-02  5.8363E-01  1.6704E+00 -1.0014E+00  2.1456E-02  3.1175E-02 -1.9921E-02  5.5845E-01 -6.8553E-03
            -1.3464E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1634.47595038067        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2716
 NPARAMETR:  9.3716E-01  9.3725E-02  2.1761E+00  1.6351E+00  1.0836E+00  9.4783E-01  7.8440E-01  9.7736E-01  7.0310E-01  1.4283E+00
             1.0599E+00
 PARAMETER:  3.5103E-02 -2.2674E+00  8.7751E-01  5.9168E-01  1.8028E-01  4.6416E-02 -1.4284E-01  7.7103E-02 -2.5225E-01  4.5650E-01
             1.5820E-01
 GRADIENT:  -8.0884E-02  1.3370E-02 -7.5360E-02  6.4185E-01  2.9201E-02 -9.3145E-03  1.3641E-02 -6.4697E-03 -1.5177E-02  3.7557E-02
             3.4373E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1634.47605445481        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2897
 NPARAMETR:  9.3716E-01  9.0063E-02  2.1807E+00  1.6373E+00  1.0836E+00  9.4780E-01  7.4743E-01  9.8183E-01  7.0219E-01  1.4282E+00
             1.0599E+00
 PARAMETER:  3.5097E-02 -2.3073E+00  8.7965E-01  5.9307E-01  1.8028E-01  4.6393E-02 -1.9112E-01  8.1664E-02 -2.5355E-01  4.5644E-01
             1.5818E-01
 GRADIENT:  -4.6212E-03 -6.8734E-03 -3.8389E-02  1.4177E-01  3.2388E-02 -8.0269E-03  1.1882E-02 -5.5975E-03  1.0829E-02  2.2737E-02
             3.0260E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1634.47644411262        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3082
 NPARAMETR:  9.3717E-01  9.2111E-02  2.1814E+00  1.6361E+00  1.0840E+00  9.4782E-01  7.4026E-01  9.8363E-01  7.0299E-01  1.4289E+00
             1.0600E+00
 PARAMETER:  3.5108E-02 -2.2848E+00  8.7997E-01  5.9231E-01  1.8062E-01  4.6406E-02 -2.0076E-01  8.3494E-02 -2.5242E-01  4.5691E-01
             1.5827E-01
 GRADIENT:  -4.1633E-02 -9.6822E-04  5.9495E-02  2.4742E-01 -2.4438E-01 -9.4111E-03  1.2386E-02  8.9949E-03  1.3560E-01  1.2544E-01
             8.6801E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1634.48056059418        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     3260
 NPARAMETR:  9.3754E-01  1.2425E-01  2.1828E+00  1.6163E+00  1.0930E+00  9.4812E-01  6.2914E-01  9.9772E-01  7.1140E-01  1.4308E+00
             1.0599E+00
 PARAMETER:  3.5506E-02 -1.9854E+00  8.8061E-01  5.8012E-01  1.8894E-01  4.6724E-02 -3.6341E-01  9.7712E-02 -2.4052E-01  4.5823E-01
             1.5821E-01
 GRADIENT:  -4.7170E-02  1.0197E-01  8.2270E-02  1.7060E+00 -1.2662E-01 -1.8625E-03  1.5085E-02  9.8089E-02 -1.6556E-01 -6.6237E-02
            -9.6257E-02

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1634.48894988001        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     3440
 NPARAMETR:  9.3806E-01  1.7297E-01  2.1536E+00  1.5850E+00  1.1012E+00  9.4847E-01  4.6018E-01  9.2734E-01  7.2576E-01  1.4398E+00
             1.0612E+00
 PARAMETER:  3.6064E-02 -1.6546E+00  8.6713E-01  5.6056E-01  1.9638E-01  4.7090E-02 -6.7613E-01  2.4567E-02 -2.2053E-01  4.6449E-01
             1.5941E-01
 GRADIENT:  -7.9852E-02  2.8642E-01  2.5054E-01  3.5979E+00 -1.3709E-02 -1.7365E-02  1.4696E-02 -6.0050E-02 -3.3472E-01  1.4932E-02
            -1.6968E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1634.49708805363        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     3616
 NPARAMETR:  9.3826E-01  1.9840E-01  2.0853E+00  1.5658E+00  1.0921E+00  9.4864E-01  3.6669E-01  8.7120E-01  7.3467E-01  1.4334E+00
             1.0609E+00
 PARAMETER:  3.6270E-02 -1.5175E+00  8.3491E-01  5.4841E-01  1.8812E-01  4.7269E-02 -9.0323E-01 -3.7884E-02 -2.0834E-01  4.6007E-01
             1.5909E-01
 GRADIENT:  -5.6382E-02  1.8386E-01 -5.3067E-03  2.0163E+00 -3.7849E-01 -2.8615E-03  1.1812E-02 -4.4325E-02 -1.8806E-01  7.2872E-03
            -3.8104E-02

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1634.50009594116        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     3793
 NPARAMETR:  9.3838E-01  2.0878E-01  2.0889E+00  1.5585E+00  1.0962E+00  9.4872E-01  3.3527E-01  8.8631E-01  7.3832E-01  1.4344E+00
             1.0608E+00
 PARAMETER:  3.6403E-02 -1.4665E+00  8.3665E-01  5.4374E-01  1.9181E-01  4.7359E-02 -9.9283E-01 -2.0690E-02 -2.0338E-01  4.6075E-01
             1.5901E-01
 GRADIENT:  -1.1354E-02  7.8532E-03 -4.3800E-02  3.1650E-02 -8.6652E-02 -1.1968E-03  1.1627E-02  1.3535E-02 -2.8479E-03  3.9170E-04
             1.0663E-03

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1634.50335661975        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     3970
 NPARAMETR:  9.3873E-01  2.4191E-01  2.0686E+00  1.5368E+00  1.1012E+00  9.4897E-01  2.4839E-01  8.5215E-01  7.4920E-01  1.4381E+00
             1.0614E+00
 PARAMETER:  3.6769E-02 -1.3192E+00  8.2688E-01  5.2968E-01  1.9644E-01  4.7622E-02 -1.2928E+00 -5.9990E-02 -1.8875E-01  4.6333E-01
             1.5954E-01
 GRADIENT:  -2.1522E-02  5.7985E-02  1.6570E-03  3.5442E-01 -8.8250E-02 -4.4165E-03  8.2413E-03 -6.9012E-03 -3.4291E-02  1.0688E-02
            -6.9660E-03

0ITERATION NO.:  135    OBJECTIVE VALUE:  -1634.50376850605        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     4149
 NPARAMETR:  9.3878E-01  2.4699E-01  2.0667E+00  1.5333E+00  1.1024E+00  9.4902E-01  2.3289E-01  8.4919E-01  7.5097E-01  1.4387E+00
             1.0614E+00
 PARAMETER:  3.6831E-02 -1.2984E+00  8.2594E-01  5.2744E-01  1.9745E-01  4.7673E-02 -1.3572E+00 -6.3475E-02 -1.8640E-01  4.6373E-01
             1.5961E-01
 GRADIENT:  -7.6858E-03  2.6004E-02  1.0856E-03  5.9256E-02 -2.8692E-02 -6.2898E-04  7.5956E-03 -4.9198E-03 -8.0870E-03  1.9352E-03
            -3.3933E-03

0ITERATION NO.:  140    OBJECTIVE VALUE:  -1634.50725371595        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     4330
 NPARAMETR:  9.3866E-01  2.3441E-01  2.0793E+00  1.5417E+00  1.1017E+00  9.4891E-01  7.3860E-02  8.6536E-01  7.4736E-01  1.4383E+00
             1.0613E+00
 PARAMETER:  3.6703E-02 -1.3507E+00  8.3203E-01  5.3288E-01  1.9682E-01  4.7563E-02 -2.5056E+00 -4.4609E-02 -1.9121E-01  4.6347E-01
             1.5949E-01
 GRADIENT:   4.0551E-03  3.4411E-04  9.4665E-03 -5.6671E-03  2.2850E-04 -5.4661E-03  7.5647E-04  6.2057E-04  6.2108E-02  5.6186E-03
            -1.4319E-03

0ITERATION NO.:  145    OBJECTIVE VALUE:  -1634.50726108328        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     4505
 NPARAMETR:  9.3862E-01  2.2987E-01  2.0821E+00  1.5447E+00  1.1010E+00  9.4888E-01  6.1902E-02  8.6931E-01  7.4583E-01  1.4379E+00
             1.0612E+00
 PARAMETER:  3.6652E-02 -1.3702E+00  8.3339E-01  5.3483E-01  1.9620E-01  4.7531E-02 -2.6822E+00 -4.0053E-02 -1.9325E-01  4.6318E-01
             1.5944E-01
 GRADIENT:   1.2456E-03  2.6021E-03  5.5948E-03  1.8123E-02 -3.1831E-03 -2.8341E-03  5.1422E-04  7.8608E-04  3.6259E-02  3.8858E-03
            -1.0520E-03

0ITERATION NO.:  150    OBJECTIVE VALUE:  -1634.50743375469        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     4689
 NPARAMETR:  9.3858E-01  2.2968E-01  2.0817E+00  1.5448E+00  1.1009E+00  9.4888E-01  3.3541E-02  8.6898E-01  7.4572E-01  1.4377E+00
             1.0613E+00
 PARAMETER:  3.6614E-02 -1.3711E+00  8.3320E-01  5.3491E-01  1.9609E-01  4.7529E-02 -3.2950E+00 -4.0438E-02 -1.9341E-01  4.6305E-01
             1.5946E-01
 GRADIENT:  -8.8968E-02  9.9892E-03 -9.7586E-03  7.5501E-02  2.2342E-02 -3.0657E-03  1.5061E-04 -1.1204E-04 -1.4513E-03 -1.3540E-02
             9.6393E-03

0ITERATION NO.:  155    OBJECTIVE VALUE:  -1634.50748591591        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     4864
 NPARAMETR:  9.3861E-01  2.2906E-01  2.0813E+00  1.5452E+00  1.1006E+00  9.4888E-01  1.0500E-02  8.6804E-01  7.4552E-01  1.4377E+00
             1.0612E+00
 PARAMETER:  3.6645E-02 -1.3738E+00  8.3301E-01  5.3514E-01  1.9583E-01  4.7532E-02 -4.4563E+00 -4.1512E-02 -1.9367E-01  4.6305E-01
             1.5943E-01
 GRADIENT:   7.6191E-03  2.6130E-03 -3.0421E-03 -1.6797E-03 -1.5238E-03  2.2105E-04  1.4872E-05 -2.3567E-03 -5.6939E-03  9.4093E-04
            -2.6769E-03

0ITERATION NO.:  157    OBJECTIVE VALUE:  -1634.50748597046        NO. OF FUNC. EVALS.:  62
 CUMULATIVE NO. OF FUNC. EVALS.:     4926
 NPARAMETR:  9.3861E-01  2.2905E-01  2.0813E+00  1.5452E+00  1.1006E+00  9.4888E-01  1.0363E-02  8.6804E-01  7.4552E-01  1.4377E+00
             1.0612E+00
 PARAMETER:  3.6645E-02 -1.3738E+00  8.3301E-01  5.3514E-01  1.9582E-01  4.7532E-02 -4.4695E+00 -4.1519E-02 -1.9367E-01  4.6305E-01
             1.5943E-01
 GRADIENT:  -2.2218E-03  6.7381E-04 -2.7251E-03  8.6027E-02 -1.8972E-03 -3.8247E-03  1.4456E-04 -3.4613E-03 -7.3002E-03  3.5504E-04
            -3.7624E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     4926
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.6772E-05 -1.0085E-04 -1.1879E-02 -8.6923E-03 -3.6044E-02
 SE:             2.9766E-02  4.4664E-05  7.6629E-03  2.9020E-02  2.3741E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9901E-01  2.3944E-02  1.2110E-01  7.6454E-01  1.2897E-01

 ETASHRINKSD(%)  2.7917E-01  9.9850E+01  7.4328E+01  2.7806E+00  2.0463E+01
 ETASHRINKVR(%)  5.5756E-01  1.0000E+02  9.3410E+01  5.4838E+00  3.6739E+01
 EBVSHRINKSD(%)  5.0296E-01  9.9861E+01  7.4988E+01  2.8324E+00  1.6169E+01
 EBVSHRINKVR(%)  1.0034E+00  1.0000E+02  9.3744E+01  5.5845E+00  2.9724E+01
 RELATIVEINF(%)  9.8088E+01  9.2688E-06  1.2743E+00  4.8057E+00  1.2760E+01
 EPSSHRINKSD(%)  4.1090E+01
 EPSSHRINKVR(%)  6.5296E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1634.5074859704648     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -899.35665940672664     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    58.97
 Elapsed covariance  time in seconds:     5.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1634.507       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.39E-01  2.29E-01  2.08E+00  1.55E+00  1.10E+00  9.49E-01  1.04E-02  8.68E-01  7.46E-01  1.44E+00  1.06E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         2.80E-02  1.09E+00  3.91E-01  7.04E-01  3.51E-01  5.88E-02  1.12E-01  2.30E-01  3.50E-01  3.04E-01  7.10E-02
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        7.84E-04
 
 TH 2
+        9.16E-03  1.20E+00
 
 TH 3
+        1.61E-03  4.04E-02  1.53E-01
 
 TH 4
+       -5.80E-03 -7.69E-01 -1.68E-02  4.95E-01
 
 TH 5
+        3.01E-03  3.67E-01  4.68E-02 -2.33E-01  1.23E-01
 
 TH 6
+       -2.29E-04  1.09E-02  3.21E-03 -6.92E-03  3.54E-03  3.46E-03
 
 TH 7
+       -7.69E-04 -1.12E-01 -1.60E-03  7.19E-02 -3.48E-02  1.33E-03  1.26E-02
 
 TH 8
+       -1.75E-03 -2.36E-01  2.38E-03  1.52E-01 -7.07E-02  1.32E-03  2.52E-02  5.31E-02
 
 TH 9
+        2.84E-03  3.80E-01  1.43E-02 -2.44E-01  1.17E-01  3.71E-03 -3.53E-02 -7.42E-02  1.23E-01
 
 TH10
+        2.27E-03  2.66E-01  4.32E-02 -1.68E-01  9.18E-02  1.90E-03 -2.44E-02 -4.91E-02  8.38E-02  9.23E-02
 
 TH11
+        3.08E-04  2.11E-02  5.23E-03 -1.33E-02  7.32E-03 -6.29E-04 -2.49E-03 -7.20E-03  6.08E-03  4.70E-03  5.05E-03
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.80E-02
 
 TH 2
+        2.99E-01  1.09E+00
 
 TH 3
+        1.47E-01  9.45E-02  3.91E-01
 
 TH 4
+       -2.94E-01 -9.98E-01 -6.12E-02  7.04E-01
 
 TH 5
+        3.07E-01  9.54E-01  3.41E-01 -9.42E-01  3.51E-01
 
 TH 6
+       -1.39E-01  1.69E-01  1.39E-01 -1.67E-01  1.71E-01  5.88E-02
 
 TH 7
+       -2.44E-01 -9.12E-01 -3.64E-02  9.08E-01 -8.81E-01  2.00E-01  1.12E-01
 
 TH 8
+       -2.71E-01 -9.37E-01  2.64E-02  9.36E-01 -8.73E-01  9.76E-02  9.73E-01  2.30E-01
 
 TH 9
+        2.90E-01  9.91E-01  1.05E-01 -9.90E-01  9.48E-01  1.80E-01 -8.96E-01 -9.19E-01  3.50E-01
 
 TH10
+        2.66E-01  7.98E-01  3.63E-01 -7.85E-01  8.60E-01  1.06E-01 -7.15E-01 -7.01E-01  7.87E-01  3.04E-01
 
 TH11
+        1.55E-01  2.72E-01  1.88E-01 -2.66E-01  2.93E-01 -1.51E-01 -3.11E-01 -4.40E-01  2.44E-01  2.18E-01  7.10E-02
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.44E+07
 
 TH 2
+       -1.66E+06  8.08E+04
 
 TH 3
+       -1.98E+05  9.48E+03  1.18E+03
 
 TH 4
+       -4.12E+06  2.00E+05  2.36E+04  4.95E+05
 
 TH 5
+       -2.28E+07  1.10E+06  1.31E+05  2.73E+06  1.51E+07
 
 TH 6
+        5.12E+07 -2.47E+06 -2.94E+05 -6.13E+06 -3.40E+07  7.61E+07
 
 TH 7
+       -1.47E+08  7.10E+06  8.45E+05  1.76E+07  9.76E+07 -2.19E+08  6.29E+08
 
 TH 8
+        6.38E+07 -3.08E+06 -3.67E+05 -7.64E+06 -4.23E+07  9.49E+07 -2.73E+08  1.18E+08
 
 TH 9
+        7.93E+06 -3.83E+05 -4.55E+04 -9.49E+05 -5.26E+06  1.18E+07 -3.39E+07  1.47E+07  1.83E+06
 
 TH10
+        4.00E+06 -1.93E+05 -2.29E+04 -4.79E+05 -2.65E+06  5.94E+06 -1.71E+07  7.41E+06  9.21E+05  4.64E+05
 
 TH11
+        3.91E+07 -1.88E+06 -2.24E+05 -4.68E+06 -2.59E+07  5.81E+07 -1.67E+08  7.24E+07  9.00E+06  4.54E+06  4.43E+07
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       64.831
Stop Time:
Sat Sep 18 12:38:56 CDT 2021
