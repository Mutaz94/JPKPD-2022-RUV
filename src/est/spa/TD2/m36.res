Wed Sep 29 18:57:39 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat36.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m36.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1663.44218474532        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1095E+02 -3.5048E+01 -1.5233E+01  2.9902E+01  9.7508E+01  4.5670E+01 -5.0308E+00 -1.0541E+00  3.3643E+01 -5.9729E+01
             3.4043E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1679.56234657798        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      193
 NPARAMETR:  1.0251E+00  1.0494E+00  9.1614E-01  9.9919E-01  9.5307E-01  1.0986E+00  1.0364E+00  9.8851E-01  8.2674E-01  1.3587E+00
             8.7281E-01
 PARAMETER:  1.2478E-01  1.4822E-01  1.2415E-02  9.9191E-02  5.1938E-02  1.9400E-01  1.3578E-01  8.8443E-02 -9.0262E-02  4.0650E-01
            -3.6039E-02
 GRADIENT:   1.8443E+01 -2.1417E+01 -1.5002E+01 -1.6949E+01  3.7771E+00  9.8766E+00 -4.3669E+00  9.9094E+00 -1.9357E+00  9.5528E+00
            -1.1608E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1682.14132342620        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  1.0098E+00  8.8992E-01  8.7697E-01  1.1045E+00  8.7827E-01  1.0953E+00  1.2718E+00  5.4264E-01  7.3969E-01  1.3628E+00
             8.5730E-01
 PARAMETER:  1.0971E-01 -1.6619E-02 -3.1277E-02  1.9939E-01 -2.9797E-02  1.9105E-01  3.4043E-01 -5.1132E-01 -2.0153E-01  4.0952E-01
            -5.3966E-02
 GRADIENT:  -9.4991E+00 -3.0291E+00 -2.2737E+01  1.3288E+01  1.7652E+01  9.0846E+00 -2.1406E+00  3.9435E+00 -3.4748E+00  1.2187E+01
            -1.6550E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1684.32578707128        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      546
 NPARAMETR:  1.0160E+00  7.5036E-01  8.5028E-01  1.1816E+00  7.8190E-01  1.0677E+00  1.4646E+00  2.9099E-01  7.2771E-01  1.2289E+00
             9.0128E-01
 PARAMETER:  1.1586E-01 -1.8721E-01 -6.2193E-02  2.6689E-01 -1.4602E-01  1.6554E-01  4.8160E-01 -1.1345E+00 -2.1785E-01  3.0615E-01
            -3.9343E-03
 GRADIENT:   2.4213E+00  3.7695E+00  3.2138E+00  4.2139E+00 -6.9952E+00 -7.7053E-01 -4.1442E-01  6.6973E-01  1.6491E+00 -4.8727E-01
             5.5275E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1684.60054174035        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      726             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0163E+00  6.7711E-01  8.4386E-01  1.2183E+00  7.5813E-01  1.0731E+00  1.6080E+00  7.6561E-02  7.0171E-01  1.2316E+00
             8.8769E-01
 PARAMETER:  1.1612E-01 -2.8992E-01 -6.9772E-02  2.9746E-01 -1.7691E-01  1.7052E-01  5.7501E-01 -2.4697E+00 -2.5424E-01  3.0834E-01
            -1.9138E-02
 GRADIENT:   6.2933E+02  5.9024E+01  8.9552E+00  4.3688E+02  2.1494E+01  1.4575E+02  2.4736E+01  8.6369E-02  1.5926E+01  5.3013E+00
             1.7199E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1684.61917217036        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      862
 NPARAMETR:  1.0146E+00  6.8012E-01  8.3925E-01  1.2178E+00  7.5812E-01  1.0682E+00  1.6012E+00  5.9973E-02  7.0207E-01  1.2364E+00
             8.8834E-01
 PARAMETER:  1.1450E-01 -2.8549E-01 -7.5241E-02  2.9701E-01 -1.7692E-01  1.6597E-01  5.7075E-01 -2.7139E+00 -2.5372E-01  3.1218E-01
            -1.8399E-02
 GRADIENT:   8.3476E-01  5.7993E-01  2.2420E-03  3.0637E+00  5.7121E-02 -2.5300E-01  2.3436E-01  2.2619E-02 -1.5890E-01  2.8281E-01
             1.4364E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1684.69044660401        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1040
 NPARAMETR:  1.0151E+00  7.0592E-01  8.2599E-01  1.2007E+00  7.5844E-01  1.0697E+00  1.5532E+00  1.0000E-02  7.0326E-01  1.2214E+00
             8.8786E-01
 PARAMETER:  1.1499E-01 -2.4825E-01 -9.1172E-02  2.8289E-01 -1.7650E-01  1.6742E-01  5.4029E-01 -5.2402E+00 -2.5203E-01  3.0002E-01
            -1.8946E-02
 GRADIENT:   1.0042E+00  1.3884E-01  5.8820E-01  7.6353E-01 -1.1547E+00  1.2992E-01 -5.9057E-01  0.0000E+00 -1.1846E+00 -5.2424E-01
            -1.0942E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1684.72417223514        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1218
 NPARAMETR:  1.0149E+00  7.2474E-01  8.2568E-01  1.1901E+00  7.6642E-01  1.0697E+00  1.5288E+00  1.0000E-02  7.1340E-01  1.2275E+00
             8.8771E-01
 PARAMETER:  1.1481E-01 -2.2195E-01 -9.1548E-02  2.7400E-01 -1.6603E-01  1.6739E-01  5.2450E-01 -6.4982E+00 -2.3771E-01  3.0501E-01
            -1.9114E-02
 GRADIENT:   3.7319E-01  2.6607E-02 -8.8452E-01  7.0731E-01  8.6982E-01  1.3852E-02  2.9893E-01  0.0000E+00  5.8892E-01  1.4554E-01
            -8.0057E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1684.73081931808        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1406
 NPARAMETR:  1.0161E+00  7.2558E-01  8.2528E-01  1.1885E+00  7.6648E-01  1.0724E+00  1.5289E+00  1.0000E-02  7.1154E-01  1.2272E+00
             8.8819E-01
 PARAMETER:  1.1600E-01 -2.2078E-01 -9.2030E-02  2.7273E-01 -1.6595E-01  1.6989E-01  5.2457E-01 -6.4982E+00 -2.4033E-01  3.0476E-01
            -1.8574E-02
 GRADIENT:   2.7022E+00 -5.7457E-01 -4.0546E-01 -1.5479E+00  4.3700E-01  1.0310E+00  2.3089E-01  0.0000E+00  1.9209E-01  2.1256E-01
             1.0519E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1684.73309681491        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1600             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0162E+00  7.2869E-01  8.2539E-01  1.1868E+00  7.6617E-01  1.0725E+00  1.5258E+00  1.0000E-02  7.1128E-01  1.2256E+00
             8.8797E-01
 PARAMETER:  1.1605E-01 -2.1650E-01 -9.1896E-02  2.7128E-01 -1.6635E-01  1.6995E-01  5.2249E-01 -6.4982E+00 -2.4069E-01  3.0345E-01
            -1.8815E-02
 GRADIENT:   6.2828E+02  5.0563E+01  6.2440E+00  3.7719E+02  2.0945E+01  1.4578E+02  2.3197E+01  0.0000E+00  1.4628E+01  6.1882E+00
             6.3980E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1684.73506573239        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1780
 NPARAMETR:  1.0162E+00  7.3046E-01  8.2405E-01  1.1867E+00  7.6693E-01  1.0721E+00  1.5222E+00  1.0000E-02  7.1185E-01  1.2253E+00
             8.8807E-01
 PARAMETER:  1.1610E-01 -2.1408E-01 -9.3519E-02  2.7115E-01 -1.6535E-01  1.6960E-01  5.2015E-01 -6.4982E+00 -2.3989E-01  3.0321E-01
            -1.8704E-02
 GRADIENT:   2.7860E+00  4.1889E-01 -1.0564E-01  5.5015E-01 -5.6810E-02  8.7684E-01  1.9983E-01  0.0000E+00  7.0800E-03 -2.1241E-02
            -3.0009E-02

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1684.73506573239        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1802
 NPARAMETR:  1.0162E+00  7.3046E-01  8.2405E-01  1.1867E+00  7.6693E-01  1.0721E+00  1.5222E+00  1.0000E-02  7.1185E-01  1.2253E+00
             8.8807E-01
 PARAMETER:  1.1610E-01 -2.1408E-01 -9.3519E-02  2.7115E-01 -1.6535E-01  1.6960E-01  5.2015E-01 -6.4982E+00 -2.3989E-01  3.0321E-01
            -1.8704E-02
 GRADIENT:   2.7860E+00  4.1889E-01 -1.0564E-01  5.5015E-01 -5.6810E-02  8.7684E-01  1.9983E-01  0.0000E+00  7.0800E-03 -2.1241E-02
            -3.0009E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1802
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -4.8757E-05  6.7591E-03 -5.2665E-04 -1.2320E-02 -9.0959E-03
 SE:             2.9900E-02  1.9970E-02  1.9652E-04  2.3723E-02  2.5559E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9870E-01  7.3502E-01  7.3661E-03  6.0354E-01  7.2193E-01

 ETASHRINKSD(%)  1.0000E-10  3.3097E+01  9.9342E+01  2.0526E+01  1.4375E+01
 ETASHRINKVR(%)  1.0000E-10  5.5240E+01  9.9996E+01  3.6839E+01  2.6684E+01
 EBVSHRINKSD(%)  2.9817E-01  3.4244E+01  9.9435E+01  2.0292E+01  1.1019E+01
 EBVSHRINKVR(%)  5.9544E-01  5.6762E+01  9.9997E+01  3.6467E+01  2.0823E+01
 RELATIVEINF(%)  9.8790E+01  3.2324E+00  5.3346E-04  5.0527E+00  1.1160E+01
 EPSSHRINKSD(%)  4.4684E+01
 EPSSHRINKVR(%)  6.9402E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1684.7350657323866     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -949.58423916864842     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.10
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.49
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1684.735       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  7.30E-01  8.24E-01  1.19E+00  7.67E-01  1.07E+00  1.52E+00  1.00E-02  7.12E-01  1.23E+00  8.88E-01
 


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
+        9.31E+02
 
 TH 2
+       -8.35E+00  4.10E+02
 
 TH 3
+        1.77E+01  1.29E+02  6.16E+02
 
 TH 4
+       -6.68E+00  4.62E+02 -2.55E+02  9.88E+02
 
 TH 5
+       -1.85E+00 -2.53E+02 -6.27E+02  2.27E+02  9.90E+02
 
 TH 6
+       -1.68E-01 -2.44E+00  3.24E+00 -1.23E+00 -8.68E-01  1.72E+02
 
 TH 7
+        1.33E+00  2.64E+01  1.01E+00 -2.12E+00 -7.61E+00  1.80E-01  2.00E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.01E+00 -1.53E+01 -3.06E+01 -1.59E+01  2.81E+01 -5.69E-01  2.66E+01  0.00E+00  1.67E+02
 
 TH10
+       -8.09E-01 -2.19E+00 -6.68E+01 -1.70E+01 -3.10E+01  1.47E-01  6.31E+00  0.00E+00  4.00E+00  7.43E+01
 
 TH11
+       -6.05E+00 -1.59E+01 -4.98E+01 -9.34E-01  1.20E+01  1.28E+00  9.76E-01  0.00E+00  1.40E+01  1.70E+01  2.72E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.649
Stop Time:
Wed Sep 29 18:58:12 CDT 2021
