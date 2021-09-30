Thu Sep 30 04:56:03 CDT 2021
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
$DATA ../../../../data/spa2/A1/dat26.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2158.48857517738        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1818E+02  5.8126E+01  5.9348E+01  3.6562E+01  3.0400E+01  5.6171E+01 -2.7058E+01 -1.0557E+02  3.6956E+01 -2.0691E+01
            -4.7757E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2200.75411455818        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2435E+00  1.0100E+00  8.8219E-01  1.0614E+00  9.3019E-01  1.3498E+00  1.0959E+00  1.3142E+00  8.2419E-01  9.5576E-01
             1.4686E+00
 PARAMETER:  3.1794E-01  1.0991E-01 -2.5346E-02  1.5958E-01  2.7630E-02  3.9998E-01  1.9154E-01  3.7324E-01 -9.3357E-02  5.4747E-02
             4.8430E-01
 GRADIENT:   8.6350E+02  7.6918E+01 -1.2260E+01  1.2964E+02 -1.2423E+01  1.0462E+02 -4.7618E+00  6.9774E-01  4.7132E-01  8.2403E+00
             4.0190E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2238.52398339791        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0095E+00  8.6542E-01  1.1677E+00  1.1064E+00  9.6572E-01  9.3767E-01  1.0048E+00  1.6997E+00  7.9311E-01  8.7814E-01
             1.0900E+00
 PARAMETER:  1.0947E-01 -4.4545E-02  2.5508E-01  2.0111E-01  6.5119E-02  3.5642E-02  1.0477E-01  6.3042E-01 -1.3180E-01 -2.9950E-02
             1.8616E-01
 GRADIENT:   4.1651E+02  7.3526E+00  4.8245E+01  1.1766E+02 -3.5462E+01  3.2288E+01 -3.6102E+01 -2.1430E+01 -1.2657E+01 -2.3513E+01
            -2.0640E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2280.00663296066        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      313
 NPARAMETR:  1.0123E+00  1.0358E+00  1.6247E+00  1.0490E+00  1.2589E+00  9.2086E-01  1.1569E+00  2.4504E+00  7.8872E-01  1.1451E+00
             1.2386E+00
 PARAMETER:  1.1224E-01  1.3521E-01  5.8533E-01  1.4787E-01  3.3025E-01  1.7555E-02  2.4574E-01  9.9624E-01 -1.3735E-01  2.3546E-01
             3.1395E-01
 GRADIENT:   5.3357E+01  5.3426E+00  9.3172E+00  3.6807E+01 -4.3018E+00  1.5495E+00 -1.0279E+01 -7.2903E+00 -6.2209E+00 -1.9449E-01
            -4.8001E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2287.81140863302        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      490
 NPARAMETR:  9.8406E-01  1.4535E+00  1.2743E+00  7.8198E-01  1.4482E+00  9.2272E-01  9.7742E-01  2.6674E+00  8.5656E-01  1.2961E+00
             1.3085E+00
 PARAMETER:  8.3930E-02  4.7400E-01  3.4236E-01 -1.4593E-01  4.7030E-01  1.9572E-02  7.7160E-02  1.0811E+00 -5.4837E-02  3.5935E-01
             3.6890E-01
 GRADIENT:  -2.5020E+01  2.8344E+01  1.8528E+00  3.3658E+01 -1.8865E+00  2.3623E+00 -5.5732E+00 -1.4207E+00 -2.8753E+00 -1.3041E+00
            -3.2241E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2289.64654068620        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      654
 NPARAMETR:  9.9277E-01  1.5750E+00  1.1480E+00  6.8687E-01  1.5268E+00  9.0950E-01  9.5538E-01  2.7742E+00  9.1854E-01  1.3589E+00
             1.3117E+00
 PARAMETER:  9.2739E-02  5.5426E-01  2.3798E-01 -2.7561E-01  5.2315E-01  5.1353E-03  5.4349E-02  1.1204E+00  1.5026E-02  4.0668E-01
             3.7135E-01
 GRADIENT:   2.3478E+02  3.1752E+02  1.4474E+00  7.0612E+01  4.9358E+01  1.4262E+01  7.4209E+00  4.8379E+01  3.9590E+00  6.5321E+00
            -2.5199E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2289.78358321113        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      831            RESET HESSIAN, TYPE II
 NPARAMETR:  9.9409E-01  1.5706E+00  1.1531E+00  6.7582E-01  1.5230E+00  9.1728E-01  9.5104E-01  2.7587E+00  8.9139E-01  1.3688E+00
             1.3142E+00
 PARAMETER:  9.4072E-02  5.5146E-01  2.4248E-01 -2.9182E-01  5.2069E-01  1.3654E-02  4.9806E-02  1.1148E+00 -1.4969E-02  4.1394E-01
             3.7322E-01
 GRADIENT:   2.3797E+02  2.9814E+02  4.8533E+00  5.5336E+01  4.2857E+01  1.7555E+01  5.1340E+00  4.4365E+01  1.4774E+00  7.5909E+00
            -2.2975E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2289.87811099587        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1008
 NPARAMETR:  9.9387E-01  1.5751E+00  1.1382E+00  6.7693E-01  1.5256E+00  9.1715E-01  9.4916E-01  2.7926E+00  8.9039E-01  1.3700E+00
             1.3168E+00
 PARAMETER:  9.3854E-02  5.5430E-01  2.2941E-01 -2.9019E-01  5.2236E-01  1.3513E-02  4.7818E-02  1.1270E+00 -1.6100E-02  4.1480E-01
             3.7522E-01
 GRADIENT:   4.4603E-01 -7.1328E+00 -7.1908E-01  3.6790E+00 -1.4477E+00  4.2986E-02  2.1311E-01  7.4377E-01 -1.0612E-02 -3.5082E-03
            -2.7530E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2289.92528458753        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1184
 NPARAMETR:  9.8885E-01  1.5950E+00  1.1107E+00  6.5767E-01  1.5349E+00  9.1510E-01  9.3522E-01  2.7534E+00  9.0406E-01  1.3724E+00
             1.3217E+00
 PARAMETER:  8.8785E-02  5.6687E-01  2.0500E-01 -3.1905E-01  5.2850E-01  1.1276E-02  3.3030E-02  1.1129E+00 -8.5526E-04  4.1660E-01
             3.7889E-01
 GRADIENT:  -1.2913E+01 -1.6079E+01 -2.4501E-01 -3.4994E+00 -2.2710E+00 -9.4682E-01 -3.1622E-01 -2.5533E+00 -5.8104E-02 -4.2393E-01
            -1.6239E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2290.00180191412        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1364
 NPARAMETR:  9.8807E-01  1.6557E+00  1.0376E+00  6.1731E-01  1.5657E+00  9.1349E-01  9.0163E-01  2.6962E+00  9.5823E-01  1.3998E+00
             1.3177E+00
 PARAMETER:  8.7999E-02  6.0423E-01  1.3693E-01 -3.8239E-01  5.4832E-01  9.5188E-03 -3.5484E-03  1.0918E+00  5.7337E-02  4.3631E-01
             3.7588E-01
 GRADIENT:  -1.4996E+01 -1.7149E+01 -4.4536E+00 -1.5252E+00 -3.5912E+00 -1.8868E+00 -9.5742E-01 -5.4414E+00  1.1301E+00  9.1222E-01
            -9.2860E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2290.22828567591        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1551
 NPARAMETR:  9.9403E-01  1.7374E+00  9.5604E-01  5.6602E-01  1.6034E+00  9.1816E-01  8.6464E-01  2.6882E+00  1.0207E+00  1.4222E+00
             1.3220E+00
 PARAMETER:  9.4014E-02  6.5237E-01  5.5040E-02 -4.6913E-01  5.7215E-01  1.4616E-02 -4.5447E-02  1.0889E+00  1.2046E-01  4.5221E-01
             3.7912E-01
 GRADIENT:   2.8676E-01 -1.1662E+01 -3.1050E+00  8.0311E-02 -3.1087E+00  6.9893E-02 -2.0299E+00 -7.8097E+00  1.8709E+00  1.3107E+00
            -1.2412E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2290.35511333166        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1720
 NPARAMETR:  9.9388E-01  1.7430E+00  9.7114E-01  5.6735E-01  1.6100E+00  9.1808E-01  8.7862E-01  2.7030E+00  9.8635E-01  1.4206E+00
             1.3195E+00
 PARAMETER:  9.3865E-02  6.5564E-01  7.0713E-02 -4.6678E-01  5.7624E-01  1.4529E-02 -2.9404E-02  1.0944E+00  8.6251E-02  4.5110E-01
             3.7723E-01
 GRADIENT:  -1.1477E+00 -1.0424E+01 -3.0419E+04 -6.5237E+03  5.2797E+03 -5.5095E-02  3.0430E+04 -1.1515E+02 -3.0429E+04  6.7458E+03
             8.0425E+03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1720
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1716E-03 -1.7505E-02 -4.3712E-02  2.1210E-02 -3.8031E-02
 SE:             2.9766E-02  2.5857E-02  1.7550E-02  1.7434E-02  2.3774E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6860E-01  4.9841E-01  1.2747E-02  2.2376E-01  1.0967E-01

 ETASHRINKSD(%)  2.7896E-01  1.3377E+01  4.1206E+01  4.1593E+01  2.0355E+01
 ETASHRINKVR(%)  5.5714E-01  2.4965E+01  6.5433E+01  6.5886E+01  3.6566E+01
 EBVSHRINKSD(%)  6.2077E-01  1.3817E+01  5.0013E+01  4.5356E+01  1.5392E+01
 EBVSHRINKVR(%)  1.2377E+00  2.5725E+01  7.5013E+01  7.0140E+01  2.8415E+01
 RELATIVEINF(%)  9.8722E+01  1.0455E+01  1.0183E+01  3.6322E+00  3.2447E+01
 EPSSHRINKSD(%)  2.8484E+01
 EPSSHRINKVR(%)  4.8855E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2290.3551133316619     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1187.6288734860548     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    32.96
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.60
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2290.355       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.94E-01  1.74E+00  9.71E-01  5.67E-01  1.61E+00  9.18E-01  8.79E-01  2.70E+00  9.86E-01  1.42E+00  1.32E+00
 


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
+        1.32E+03
 
 TH 2
+       -9.41E+00  5.88E+04
 
 TH 3
+       -2.00E+02  3.09E+01  8.06E+06
 
 TH 4
+       -2.89E+06 -2.48E+02 -1.88E+03  2.17E+06
 
 TH 5
+        1.36E+01 -3.01E+01  5.00E+02 -3.09E+05  8.85E+04
 
 TH 6
+        2.80E+00 -1.95E+00 -1.95E+02 -7.37E+01  1.94E+01  2.30E+02
 
 TH 7
+        2.17E+02  1.82E+01 -8.91E+06 -6.89E+02  1.87E+02  2.14E+02  9.86E+06
 
 TH 8
+        2.17E+00 -2.21E+04 -4.05E+01  9.89E+04 -1.93E+01  3.23E-01  1.65E+00  1.74E+04
 
 TH 9
+       -1.92E+02 -6.90E+00 -4.56E+03  2.91E+06 -8.31E+05 -1.92E+02 -1.76E+03  1.47E+00  7.82E+06
 
 TH10
+        3.05E+01 -2.62E+00  7.01E+02 -4.48E+05 -5.63E+01  2.98E+01  2.72E+02 -8.50E-01 -9.58E+01  1.85E+05
 
 TH11
+        2.37E+01 -7.69E+00  9.42E+02 -5.76E+05 -2.20E+01  4.07E+01  3.60E+02 -5.06E+04 -1.27E+02  6.92E+01  3.07E+05
 
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
 #CPUT: Total CPU Time in Seconds,       42.635
Stop Time:
Thu Sep 30 04:56:47 CDT 2021
