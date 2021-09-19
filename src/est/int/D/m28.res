Sat Sep 18 06:46:24 CDT 2021
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
$DATA ../../../../data/int/D/dat28.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   20098.9163471511        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.7827E+01  8.2325E+01 -2.0770E+01 -2.0926E+02  1.2391E+02 -1.1866E+03 -5.8258E+02 -1.1326E+02 -1.1002E+03 -2.8790E+02
            -4.4643E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1029.26167318387        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.5503E+00  1.7873E+00  9.1848E-01  2.5483E+00  9.5649E-01  5.8146E+00  4.0386E+00  1.0360E+00  2.8679E+00  1.6796E+00
             1.2615E+01
 PARAMETER:  5.3846E-01  6.8068E-01  1.4966E-02  1.0354E+00  5.5519E-02  1.8604E+00  1.4959E+00  1.3536E-01  1.1536E+00  6.1857E-01
             2.6348E+00
 GRADIENT:   5.0305E+00  9.8323E+00 -4.2936E+01  1.1202E+02  8.3894E+00  1.6165E+02  3.8741E+01  4.6216E+00  3.7449E+01  3.6454E+01
             6.0410E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1128.33794323868        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.6600E-01  1.2495E+00  1.6870E+02  7.1709E+00  2.3895E+00  3.1916E+00  1.2992E+01  1.4383E+00  3.7942E+00  1.5304E+00
             1.1696E+01
 PARAMETER:  6.5413E-02  3.2277E-01  5.2281E+00  2.0700E+00  9.7107E-01  1.2605E+00  2.6643E+00  4.6343E-01  1.4335E+00  5.2556E-01
             2.5592E+00
 GRADIENT:  -6.9935E+01  1.7321E+01 -1.3274E-01  1.5113E+02 -3.5976E+01  7.5691E+01  1.5191E+00  6.3966E-03 -3.2564E+01  2.8888E+01
             5.8162E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1461.10093088317        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  8.5998E-01  2.2253E+00  1.0755E+01  8.7813E-01  2.6433E+00  2.2357E+00  4.1288E+00  2.9821E+00  1.5309E+00  2.6028E-01
             7.7262E+00
 PARAMETER: -5.0841E-02  8.9991E-01  2.4754E+00 -2.9959E-02  1.0720E+00  9.0454E-01  1.5180E+00  1.1926E+00  5.2585E-01 -1.2460E+00
             2.1446E+00
 GRADIENT:  -1.2605E+02  1.2160E+01 -1.3290E+00 -2.0472E+01  5.3525E+01 -1.6528E+01  5.1302E+01 -1.2334E+00  1.4550E+01  3.0076E-01
             2.1414E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1504.12314530709        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.1873E+00  2.3074E+00  1.8691E+00  6.3417E-01  1.9209E+00  2.3304E+00  2.7127E+00  2.6511E+00  1.3845E+00  2.2805E-01
             6.7997E+00
 PARAMETER:  2.7166E-01  9.3613E-01  7.2544E-01 -3.5544E-01  7.5279E-01  9.4605E-01  1.0979E+00  1.0750E+00  4.2537E-01 -1.3782E+00
             2.0169E+00
 GRADIENT:  -1.2477E+00  7.8706E+00  1.8727E+00 -3.9053E+00 -1.5894E+00  9.6189E+00 -1.8964E+00  5.0148E-01  5.2342E+00  2.1669E-01
            -5.9565E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1506.93416742637        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  1.1795E+00  1.9266E+00  1.6543E+00  7.7899E-01  1.7046E+00  2.2281E+00  3.0866E+00  1.1801E+00  8.7473E-01  2.0207E-01
             6.8510E+00
 PARAMETER:  2.6505E-01  7.5577E-01  6.0336E-01 -1.4976E-01  6.3334E-01  9.0116E-01  1.2271E+00  2.6563E-01 -3.3844E-02 -1.4991E+00
             2.0244E+00
 GRADIENT:  -4.2387E+00 -1.4412E+00  1.1860E+00 -6.3643E+00 -3.7336E+00 -5.9119E+00  2.8883E+00 -1.6321E+00  4.6401E+00 -6.1294E-02
             6.9577E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1515.98387225873        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  1.1646E+00  1.8225E+00  2.1146E+00  8.2921E-01  1.7995E+00  2.1868E+00  3.4758E+00  4.4658E-01  2.0898E-01  1.5876E-01
             6.9296E+00
 PARAMETER:  2.5242E-01  7.0022E-01  8.4885E-01 -8.7277E-02  6.8752E-01  8.8245E-01  1.3458E+00 -7.0614E-01 -1.4655E+00 -1.7403E+00
             2.0358E+00
 GRADIENT:  -3.1190E+00  3.7509E-01  2.8232E-01 -3.4573E-01  2.4948E+00 -3.0730E+00  1.8892E+00 -3.7974E-01  4.2873E-01 -3.2514E-01
             6.9296E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1517.04801170439        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  1.1790E+00  1.9717E+00  1.9404E+00  7.8638E-01  1.7852E+00  2.2311E+00  3.2481E+00  1.4081E+00  4.2507E-02  3.2192E-01
             6.9016E+00
 PARAMETER:  2.6468E-01  7.7891E-01  7.6288E-01 -1.4032E-01  6.7953E-01  9.0248E-01  1.2781E+00  4.4221E-01 -3.0581E+00 -1.0335E+00
             2.0318E+00
 GRADIENT:   2.1244E+00  5.9550E+00 -6.1934E-01  2.7444E+00  3.6086E-01  4.8328E+00 -5.2497E+00 -5.3463E-01  7.7670E-03 -3.8818E-01
             6.3219E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1517.57122471106        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      590
 NPARAMETR:  1.1718E+00  1.7991E+00  2.2624E+00  8.3235E-01  1.7909E+00  2.2032E+00  3.4591E+00  1.5368E+00  1.0000E-02  4.8432E-01
             6.8538E+00
 PARAMETER:  2.5855E-01  6.8730E-01  9.1642E-01 -8.3497E-02  6.8274E-01  8.8991E-01  1.3410E+00  5.2973E-01 -4.6770E+00 -6.2500E-01
             2.0248E+00
 GRADIENT:  -8.4796E-02 -1.0590E+00 -4.2448E-01  5.4565E-01  9.9352E-01 -1.5046E-01 -7.7043E-02  2.7252E-01  0.0000E+00 -1.4266E-01
             7.6097E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1517.63781061216        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      660
 NPARAMETR:  1.1721E+00  1.8150E+00  2.2207E+00  8.2820E-01  1.7854E+00  2.2062E+00  3.4596E+00  1.1697E+00  1.0000E-02  6.1491E-01
             6.8402E+00
 PARAMETER:  2.5876E-01  6.9608E-01  8.9781E-01 -8.8503E-02  6.7963E-01  8.9125E-01  1.3412E+00  2.5674E-01 -5.8609E+00 -3.8628E-01
             2.0228E+00
 GRADIENT:   2.4062E-02  3.3243E-02  4.5872E-02 -4.1343E-01 -1.3004E-01 -2.6267E-02  2.8129E-01  4.4209E-02  0.0000E+00  5.6397E-02
             6.1323E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1517.70811716613        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      839
 NPARAMETR:  1.1790E+00  1.8205E+00  2.2797E+00  8.3646E-01  1.7941E+00  2.2284E+00  3.5316E+00  9.8734E-01  1.0000E-02  6.4861E-01
             6.8515E+00
 PARAMETER:  2.6470E-01  6.9911E-01  9.2403E-01 -7.8578E-02  6.8449E-01  9.0129E-01  1.3618E+00  8.7261E-02 -6.2352E+00 -3.3292E-01
             2.0245E+00
 GRADIENT:  -1.4146E-03  3.5630E-02 -1.0482E-01  1.0266E-02  1.3228E-01 -1.0798E-01  2.6374E-01  1.7380E-02  0.0000E+00  2.1815E-02
             1.0439E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1517.71818699808        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1017
 NPARAMETR:  1.1792E+00  1.8150E+00  2.3086E+00  8.3809E-01  1.8000E+00  2.2318E+00  3.5445E+00  4.5264E-01  1.0000E-02  7.2072E-01
             6.8443E+00
 PARAMETER:  2.6484E-01  6.9611E-01  9.3662E-01 -7.6628E-02  6.8778E-01  9.0280E-01  1.3654E+00 -6.9266E-01 -7.4031E+00 -2.2751E-01
             2.0234E+00
 GRADIENT:   7.3328E-02 -1.9530E-03 -5.9035E-02  2.7315E-02  2.9123E-01  1.7433E-01  4.7095E-02  1.8439E-02  0.0000E+00  9.9233E-02
             1.6088E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1517.72727927305        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1192
 NPARAMETR:  1.1790E+00  1.8127E+00  2.3140E+00  8.3888E-01  1.8002E+00  2.2310E+00  3.5490E+00  5.4021E-02  1.0000E-02  7.3110E-01
             6.8433E+00
 PARAMETER:  2.6466E-01  6.9480E-01  9.3896E-01 -7.5690E-02  6.8788E-01  9.0246E-01  1.3667E+00 -2.8184E+00 -9.5163E+00 -2.1320E-01
             2.0233E+00
 GRADIENT:   1.4423E-02  7.6699E-03 -1.3762E-03  1.9887E-02  2.3399E-02  7.2165E-03  7.0910E-04  3.0898E-04  0.0000E+00  9.5815E-03
             4.2460E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1517.72743184552        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1354
 NPARAMETR:  1.1789E+00  1.8127E+00  2.3127E+00  8.3881E-01  1.7999E+00  2.2310E+00  3.5488E+00  1.0000E-02  1.0000E-02  7.3071E-01
             6.8432E+00
 PARAMETER:  2.6462E-01  6.9480E-01  9.3841E-01 -7.5769E-02  6.8774E-01  9.0244E-01  1.3666E+00 -4.9135E+00 -1.1540E+01 -2.1374E-01
             2.0233E+00
 GRADIENT:  -9.4282E-04  7.0536E-04  1.3290E-03  2.4922E-03 -6.8078E-03 -8.8180E-04 -3.5660E-03  0.0000E+00  0.0000E+00 -1.3698E-04
            -2.9216E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1354
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.2910E-03  2.8175E-03 -2.4681E-05 -7.1080E-04 -4.8910E-03
 SE:             2.9165E-02  2.8248E-02  5.0011E-05  1.8123E-04  1.2955E-02
 N:                     100         100         100         100         100

 P VAL.:         8.0260E-01  9.2055E-01  6.2166E-01  8.7863E-05  7.0578E-01

 ETASHRINKSD(%)  2.2934E+00  5.3671E+00  9.9832E+01  9.9393E+01  5.6598E+01
 ETASHRINKVR(%)  4.5343E+00  1.0446E+01  1.0000E+02  9.9996E+01  8.1162E+01
 EBVSHRINKSD(%)  2.1991E+00  3.5221E+00  9.9836E+01  9.9551E+01  5.7572E+01
 EBVSHRINKVR(%)  4.3497E+00  6.9201E+00  1.0000E+02  9.9998E+01  8.1999E+01
 RELATIVEINF(%)  9.5545E+01  4.7958E+01  6.2642E-05  7.5804E-04  4.7267E+00
 EPSSHRINKSD(%)  8.5014E+00
 EPSSHRINKVR(%)  1.6280E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1517.7274318455206     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       136.36192792289012     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.41
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.25
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1517.727       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.18E+00  1.81E+00  2.31E+00  8.39E-01  1.80E+00  2.23E+00  3.55E+00  1.00E-02  1.00E-02  7.31E-01  6.84E+00
 


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
+        1.52E+02
 
 TH 2
+       -8.18E-01  2.48E+01
 
 TH 3
+        6.17E-01  1.43E+00  6.09E+00
 
 TH 4
+       -6.08E+00  4.29E+01 -2.14E+01  4.42E+02
 
 TH 5
+       -2.74E+00 -8.39E+00 -2.18E+01  7.92E+01  1.06E+02
 
 TH 6
+        5.05E-01 -1.76E-01  7.45E-02  8.96E-01 -1.16E+00  3.60E+01
 
 TH 7
+        2.05E-01  1.19E+00 -1.15E+00 -3.68E+01  2.50E+00 -5.07E-01  1.24E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -1.18E+00 -1.03E+00 -1.98E+00  7.31E-01  5.99E+00 -3.89E-01  5.72E-01  0.00E+00  0.00E+00  1.11E+01
 
 TH11
+       -5.62E+00 -3.22E+00  3.05E-01 -1.68E+01 -1.67E+00  5.50E-01  1.08E+00  0.00E+00  0.00E+00  4.29E+00  2.36E+01
 
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
 #CPUT: Total CPU Time in Seconds,       50.792
Stop Time:
Sat Sep 18 06:47:17 CDT 2021
