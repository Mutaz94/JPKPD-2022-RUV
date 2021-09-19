Sat Sep 18 11:28:42 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat7.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1158.72863401990        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.6769E+01 -1.0764E+01  4.7652E+01  4.2596E+00  7.5200E+01 -3.0677E+01  1.5175E+01 -2.0308E+02 -4.7517E+00 -1.1173E+01
            -7.6598E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1652.84177949727        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:       92
 NPARAMETR:  9.9362E-01  9.9609E-01  9.8620E-01  9.7285E-01  9.6251E-01  1.0411E+00  9.4640E-01  1.0743E+00  9.5300E-01  9.8511E-01
             1.1405E+00
 PARAMETER:  9.3602E-02  9.6082E-02  8.6099E-02  7.2474E-02  6.1791E-02  1.4031E-01  4.4905E-02  1.7167E-01  5.1856E-02  8.4995E-02
             2.3147E-01
 GRADIENT:   3.7713E+01 -2.2185E+01 -7.1833E+00 -2.0218E+01  1.3876E+01 -1.5366E+01  8.9231E+00  2.9326E+00  8.4290E+00  4.5127E+00
             3.3313E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1653.81787464882        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      275
 NPARAMETR:  9.9338E-01  1.0071E+00  9.6566E-01  9.6527E-01  9.5878E-01  1.0690E+00  7.2882E-01  9.4478E-01  1.0495E+00  1.0580E+00
             1.1286E+00
 PARAMETER:  9.3357E-02  1.0707E-01  6.5059E-02  6.4654E-02  5.7908E-02  1.6677E-01 -2.1633E-01  4.3201E-02  1.4828E-01  1.5640E-01
             2.2095E-01
 GRADIENT:   1.0714E+00 -2.9392E+01 -7.1447E+00 -2.0143E+01  7.3794E+00 -1.2232E+01  4.6874E+00  1.7043E-01  1.6953E+01  9.4075E+00
             3.0174E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1656.41540597662        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      463
 NPARAMETR:  9.9524E-01  1.1476E+00  8.6600E-01  8.9211E-01  9.6229E-01  1.0927E+00  6.8570E-01  8.8838E-01  1.0495E+00  1.0190E+00
             1.0602E+00
 PARAMETER:  9.5227E-02  2.3765E-01 -4.3871E-02 -1.4163E-02  6.1563E-02  1.8862E-01 -2.7732E-01 -1.8355E-02  1.4835E-01  1.1886E-01
             1.5848E-01
 GRADIENT:   2.4851E+00  7.5431E+00  3.6970E+00  7.8523E+00 -6.2252E+00 -4.3041E+00  2.2830E-01 -2.3274E+00  2.9753E-02  3.5773E+00
             2.9920E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1656.46160169080        NO. OF FUNC. EVALS.: 154
 CUMULATIVE NO. OF FUNC. EVALS.:      617
 NPARAMETR:  9.9519E-01  1.1474E+00  8.6600E-01  8.9211E-01  9.6231E-01  1.0996E+00  6.8567E-01  8.8843E-01  1.0496E+00  9.8488E-01
             1.0603E+00
 PARAMETER:  9.5176E-02  2.3753E-01 -4.3873E-02 -1.4161E-02  6.1583E-02  1.9498E-01 -2.7736E-01 -1.8297E-02  1.4836E-01  8.4760E-02
             1.5856E-01
 GRADIENT:   2.3463E+00  7.4298E+00  4.9827E+00  7.9966E+00 -2.3527E+00 -1.7976E+00 -5.1153E-01 -3.1615E+00  2.9182E-02 -2.0707E+00
             2.1228E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1656.50669727285        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      763
 NPARAMETR:  9.9502E-01  1.1471E+00  8.6557E-01  8.9189E-01  9.6241E-01  1.1043E+00  6.8583E-01  8.9516E-01  1.0495E+00  9.9727E-01
             1.0600E+00
 PARAMETER:  9.5003E-02  2.3726E-01 -4.4363E-02 -1.4415E-02  6.1690E-02  1.9921E-01 -2.7713E-01 -1.0757E-02  1.4834E-01  9.7270E-02
             1.5827E-01
 GRADIENT:   4.3671E+01  1.8289E+01  3.9611E+00  1.3086E+01 -2.3532E+00  1.3188E+01  4.1494E-01 -2.6424E+00  9.2701E-01  2.7321E-01
             2.6537E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1656.50861952932        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:      937
 NPARAMETR:  9.9492E-01  1.1468E+00  8.6549E-01  8.9180E-01  9.6232E-01  1.1046E+00  6.8602E-01  8.9507E-01  1.0494E+00  9.9513E-01
             1.0602E+00
 PARAMETER:  9.4903E-02  2.3702E-01 -4.4463E-02 -1.4515E-02  6.1590E-02  1.9950E-01 -2.7685E-01 -1.0857E-02  1.4819E-01  9.5120E-02
             1.5843E-01
 GRADIENT:   1.8409E+00  5.8497E+00  3.8520E+00  6.9347E+00 -2.7483E+00 -5.5077E-03 -2.0282E-01 -2.6941E+00  6.7488E-02 -1.3034E-01
             2.5512E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1656.54523868178        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1117             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9395E-01  1.1467E+00  8.5776E-01  8.8787E-01  9.6413E-01  1.1046E+00  6.8962E-01  8.9561E-01  1.0491E+00  9.9586E-01
             1.0602E+00
 PARAMETER:  9.3935E-02  2.3691E-01 -5.3434E-02 -1.8931E-02  6.3469E-02  1.9950E-01 -2.7161E-01 -1.0255E-02  1.4792E-01  9.5849E-02
             1.5842E-01
 GRADIENT:   4.1699E+01  8.5753E+00 -1.2486E-01  7.5310E+00  5.4118E+00  1.3344E+01  5.0196E-01 -2.1384E+00  9.7584E-01  3.3393E-01
             3.1163E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1656.59058785688        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1277
 NPARAMETR:  9.9396E-01  1.1467E+00  8.2877E-01  8.8590E-01  9.4448E-01  1.1052E+00  7.1529E-01  8.9561E-01  1.0433E+00  9.7060E-01
             1.0602E+00
 PARAMETER:  9.3940E-02  2.3691E-01 -8.7813E-02 -2.1146E-02  4.2878E-02  2.0005E-01 -2.3507E-01 -1.0255E-02  1.4242E-01  7.0155E-02
             1.5842E-01
 GRADIENT:  -3.1984E-01 -3.3109E+00  1.6905E-01  6.3224E-01 -1.4210E+00  8.9587E-02  2.6813E-01 -1.0237E+00  1.2596E+00  5.6174E-01
             3.8387E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1656.62377175201        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1453            RESET HESSIAN, TYPE II
 NPARAMETR:  9.9415E-01  1.1479E+00  8.2832E-01  8.8567E-01  9.4491E-01  1.1050E+00  7.1859E-01  9.0949E-01  1.0364E+00  9.6677E-01
             1.0565E+00
 PARAMETER:  9.4128E-02  2.3793E-01 -8.8358E-02 -2.1408E-02  4.3333E-02  1.9984E-01 -2.3046E-01  5.1274E-03  1.3575E-01  6.6207E-02
             1.5500E-01
 GRADIENT:   4.1823E+01  9.4134E+00 -1.7984E-01  7.0570E+00  2.2998E-01  1.3402E+01  5.9713E-01 -7.8598E-01  6.7879E-01  3.4967E-01
             2.5328E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1656.62990181301        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     1609
 NPARAMETR:  9.9412E-01  1.1481E+00  8.2920E-01  8.8564E-01  9.4528E-01  1.1050E+00  7.1669E-01  9.1296E-01  1.0369E+00  9.6457E-01
             1.0557E+00
 PARAMETER:  9.4100E-02  2.3810E-01 -8.7297E-02 -2.1449E-02  4.3726E-02  1.9988E-01 -2.3312E-01  8.9407E-03  1.3624E-01  6.3924E-02
             1.5419E-01
 GRADIENT:  -1.5572E-02 -2.4854E+00 -2.1084E-01  1.3596E+00 -3.6887E-01 -9.5185E-03 -3.1134E-02 -8.4891E-01 -1.5681E-01 -1.5955E-01
             2.0085E+00

0ITERATION NO.:   53    OBJECTIVE VALUE:  -1656.63009941588        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1707
 NPARAMETR:  9.9412E-01  1.1481E+00  8.2940E-01  8.8564E-01  9.4540E-01  1.1050E+00  7.1689E-01  9.1296E-01  1.0372E+00  9.6500E-01
             1.0557E+00
 PARAMETER:  9.4103E-02  2.3810E-01 -8.7047E-02 -2.1449E-02  4.3851E-02  1.9989E-01 -2.3283E-01  8.9407E-03  1.3652E-01  6.4371E-02
             1.5419E-01
 GRADIENT:   2.2953E+06  9.6400E+05  2.2954E+06  1.1477E+06 -1.1477E+06 -1.4779E-02 -9.8584E+05  2.2952E+06  1.6813E+06  2.2953E+06
            -1.4889E+06

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1707
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2627E-04 -2.5284E-02 -2.4503E-02  9.1951E-03 -3.1645E-02
 SE:             2.9875E-02  1.6511E-02  1.2360E-02  2.5638E-02  2.2825E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9129E-01  1.2567E-01  4.7437E-02  7.1986E-01  1.6563E-01

 ETASHRINKSD(%)  1.0000E-10  4.4688E+01  5.8592E+01  1.4108E+01  2.3532E+01
 ETASHRINKVR(%)  1.0000E-10  6.9405E+01  8.2854E+01  2.6226E+01  4.1527E+01
 EBVSHRINKSD(%)  3.8671E-01  4.4778E+01  6.2067E+01  1.4458E+01  2.1801E+01
 EBVSHRINKVR(%)  7.7193E-01  6.9506E+01  8.5611E+01  2.6826E+01  3.8849E+01
 RELATIVEINF(%)  9.8893E+01  9.1491E-01  1.3803E+00  2.8032E+00  8.7071E+00
 EPSSHRINKSD(%)  4.4897E+01
 EPSSHRINKVR(%)  6.9637E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1656.6300994158785     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -921.47927285214030     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.09
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1656.630       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.94E-01  1.15E+00  8.29E-01  8.86E-01  9.45E-01  1.11E+00  7.17E-01  9.13E-01  1.04E+00  9.65E-01  1.06E+00
 


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
+        5.81E+09
 
 TH 2
+        3.89E+03  7.68E+08
 
 TH 3
+        2.23E+05  8.13E+04  8.34E+09
 
 TH 4
+       -6.52E+09 -3.74E+04  2.50E+05  7.32E+09
 
 TH 5
+        4.49E+04  1.60E+04  7.32E+09  6.85E+09  6.42E+09
 
 TH 6
+        1.21E+04  4.39E+03  1.45E+04  1.36E+04 -1.27E+04  1.61E+02
 
 TH 7
+        9.03E+03  3.28E+03  1.08E+04  1.01E+04 -3.64E+09 -7.20E+03  2.06E+09
 
 TH 8
+        1.17E+04 -2.00E+05  2.43E+05 -1.14E+05  4.89E+04  1.32E+04  9.83E+03  6.88E+09
 
 TH 9
+       -8.39E+04 -3.05E+04 -4.89E+09 -4.58E+09 -4.29E+09  8.49E+03  2.43E+09 -9.13E+04  2.86E+09
 
 TH10
+       -1.16E+04 -4.22E+03 -1.40E+04 -6.71E+09 -6.29E+09  1.24E+04  9.31E+03 -1.26E+04  4.20E+09  6.16E+09
 
 TH11
+       -6.55E+03 -5.82E+05 -1.36E+05  6.37E+04 -2.74E+04 -7.38E+03 -5.51E+03  3.36E+05  5.12E+04  7.11E+03  2.17E+09
 
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
 #CPUT: Total CPU Time in Seconds,       28.033
Stop Time:
Sat Sep 18 11:29:12 CDT 2021
