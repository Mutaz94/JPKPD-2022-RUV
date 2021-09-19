Sat Sep 18 12:23:49 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat66.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1619.89954616663        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4543E+02 -7.8639E+01 -8.9787E+00 -7.9203E+01  2.4222E+01  1.0789E+01 -9.6338E+00  3.7660E+00  2.2707E+01 -1.7090E+01
            -2.4638E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1629.42248799970        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.8126E-01  1.1800E+00  9.9382E-01  9.4522E-01  1.0967E+00  9.4085E-01  1.1473E+00  9.4315E-01  7.2846E-01  1.1972E+00
             1.0944E+00
 PARAMETER:  8.1084E-02  2.6548E-01  9.3798E-02  4.3659E-02  1.9235E-01  3.9024E-02  2.3741E-01  4.1466E-02 -2.1682E-01  2.8001E-01
             1.9018E-01
 GRADIENT:   1.0182E+02  2.1565E+00 -2.4579E+00 -8.7313E+00  1.0445E+01 -9.8116E+00 -3.3097E-01 -5.2769E-02 -7.2915E+00  2.7936E+00
             9.7487E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1630.85057876386        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.7638E-01  9.5907E-01  1.0018E+00  1.0895E+00  9.8163E-01  9.5399E-01  1.4732E+00  8.1311E-01  6.4159E-01  1.0867E+00
             1.1038E+00
 PARAMETER:  7.6100E-02  5.8210E-02  1.0176E-01  1.8573E-01  8.1461E-02  5.2897E-02  4.8741E-01 -1.0689E-01 -3.4380E-01  1.8316E-01
             1.9875E-01
 GRADIENT:   8.9015E+01  1.3663E+01 -5.0838E+00  1.3170E+01  1.0503E+01 -3.3293E+00  1.0030E+01  4.1924E-01 -4.4740E+00 -2.7552E-01
             1.3813E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1632.55386512343        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.3931E-01  9.5925E-01  8.2280E-01  1.0708E+00  8.7228E-01  9.5708E-01  1.3231E+00  5.1146E-01  7.3881E-01  9.8008E-01
             1.0616E+00
 PARAMETER:  3.7387E-02  5.8399E-02 -9.5037E-02  1.6838E-01 -3.6640E-02  5.6134E-02  3.7997E-01 -5.7049E-01 -2.0272E-01  7.9875E-02
             1.5974E-01
 GRADIENT:  -5.4611E+00  1.4596E+00 -4.2602E+00  6.0365E+00  3.7452E+00 -1.3933E+00 -3.1262E-01  8.6716E-01  8.0465E-01  1.1179E+00
             1.5266E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1633.44118390839        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.4007E-01  7.7510E-01  7.0651E-01  1.1618E+00  7.1624E-01  9.6842E-01  1.5965E+00  1.7660E-01  6.6974E-01  8.2674E-01
             1.0533E+00
 PARAMETER:  3.8199E-02 -1.5476E-01 -2.4742E-01  2.5000E-01 -2.3374E-01  6.7912E-02  5.6783E-01 -1.6338E+00 -3.0086E-01 -9.0260E-02
             1.5193E-01
 GRADIENT:  -3.6878E+00  8.3299E+00  7.2501E+00  1.0608E+01 -7.3292E+00  3.2326E+00  8.9460E-01  1.1316E-02 -3.5064E+00 -2.6816E+00
            -4.2422E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1633.45934174368        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  9.4174E-01  7.7551E-01  6.8661E-01  1.1552E+00  7.0573E-01  9.6403E-01  1.5784E+00  1.5867E-01  6.7969E-01  8.1519E-01
             1.0572E+00
 PARAMETER:  3.9970E-02 -1.5424E-01 -2.7599E-01  2.4432E-01 -2.4852E-01  6.3363E-02  5.5641E-01 -1.7410E+00 -2.8612E-01 -1.0434E-01
             1.5562E-01
 GRADIENT:  -7.8320E-01  3.5325E+00  2.3804E+00  4.7818E+00 -2.8809E+00  1.3898E+00  3.6326E-01  1.0267E-01 -1.4664E+00 -9.0498E-01
            -1.6482E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1633.46368208987        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  9.4231E-01  7.5866E-01  6.7026E-01  1.1599E+00  6.8968E-01  9.6217E-01  1.5978E+00  1.3115E-01  6.8018E-01  7.9851E-01
             1.0590E+00
 PARAMETER:  4.0583E-02 -1.7620E-01 -3.0009E-01  2.4837E-01 -2.7153E-01  6.1441E-02  5.6862E-01 -1.9314E+00 -2.8539E-01 -1.2501E-01
             1.5729E-01
 GRADIENT:   8.0360E-02  1.5217E+00  5.9467E-01  2.3012E+00 -1.0897E+00  5.9746E-01  1.5119E-01  1.0267E-01 -6.1550E-01 -2.4945E-01
            -6.2718E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1633.46590555605        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      514
 NPARAMETR:  9.4262E-01  7.3556E-01  6.5457E-01  1.1684E+00  6.7140E-01  9.6112E-01  1.6314E+00  1.0408E-01  6.7758E-01  7.7940E-01
             1.0601E+00
 PARAMETER:  4.0903E-02 -2.0713E-01 -3.2377E-01  2.5566E-01 -2.9839E-01  6.0341E-02  5.8943E-01 -2.1626E+00 -2.8923E-01 -1.4923E-01
             1.5835E-01
 GRADIENT:   4.5850E-01  5.0613E-01 -2.5610E-01  1.0313E+00 -1.7058E-01  1.4203E-01  2.9369E-02  8.1206E-02 -1.5559E-01  6.8473E-02
            -8.0905E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1633.46609456089        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      585
 NPARAMETR:  9.4265E-01  7.2444E-01  6.4831E-01  1.1728E+00  6.6339E-01  9.6093E-01  1.6493E+00  9.3425E-02  6.7576E-01  7.7097E-01
             1.0603E+00
 PARAMETER:  4.0943E-02 -2.2236E-01 -3.3339E-01  2.5941E-01 -3.1039E-01  6.0150E-02  6.0033E-01 -2.2706E+00 -2.9192E-01 -1.6011E-01
             1.5860E-01
 GRADIENT:   4.5747E-01  3.0403E-01 -3.3328E-01  7.3330E-01 -3.7939E-02  6.1975E-02  9.7230E-03  6.9938E-02 -7.5811E-02  1.0055E-01
             1.4434E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1634.05653601832        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.5278E-01  7.1803E-01  7.0840E-01  1.1901E+00  7.0127E-01  9.6237E-01  1.6770E+00  1.0972E-01  6.8396E-01  8.2590E-01
             1.0697E+00
 PARAMETER:  5.1627E-02 -2.3124E-01 -2.4475E-01  2.7400E-01 -2.5486E-01  6.1642E-02  6.1699E-01 -2.1098E+00 -2.7986E-01 -9.1285E-02
             1.6739E-01
 GRADIENT:  -3.4598E+00 -1.2523E+00 -1.6130E-01 -8.4335E+00  5.8869E-02 -1.8550E+00 -1.5371E+00  8.6862E-03  2.0009E+00 -1.0615E+00
             2.7641E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1634.48065359095        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      878
 NPARAMETR:  9.5217E-01  5.7776E-01  7.3067E-01  1.2710E+00  6.6987E-01  9.6555E-01  2.0224E+00  3.0305E-02  6.4161E-01  8.4170E-01
             1.0659E+00
 PARAMETER:  5.0990E-02 -4.4859E-01 -2.1379E-01  3.3981E-01 -3.0067E-01  6.4942E-02  8.0430E-01 -3.3965E+00 -3.4378E-01 -7.2329E-02
             1.6386E-01
 GRADIENT:  -7.3519E-01  7.0364E-01  1.5417E+00  1.5597E-01 -1.6063E+00  2.7865E-01  6.0412E-02 -4.1458E-03 -4.2716E-01 -4.5392E-01
            -4.7830E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1634.48311319950        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1040
 NPARAMETR:  9.5243E-01  5.7365E-01  7.3024E-01  1.2728E+00  6.6896E-01  9.6481E-01  2.0305E+00  2.8210E-02  6.4234E-01  8.4317E-01
             1.0667E+00
 PARAMETER:  5.1265E-02 -4.5573E-01 -2.1438E-01  3.4121E-01 -3.0202E-01  6.4176E-02  8.0830E-01 -3.4681E+00 -3.4263E-01 -7.0589E-02
             1.6457E-01
 GRADIENT:  -2.4989E-03  3.5212E-02  1.7660E-02  4.4546E-02 -2.0356E-02 -3.4467E-04 -4.7346E-03 -3.1470E-03  3.5921E-03  5.6592E-04
             9.4484E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1040
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.1096E-04  2.5831E-02 -1.4313E-03 -2.7315E-02  5.2233E-03
 SE:             2.9817E-02  2.1486E-02  6.6158E-04  2.3158E-02  2.3344E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8365E-01  2.2929E-01  3.0505E-02  2.3818E-01  8.2295E-01

 ETASHRINKSD(%)  1.0844E-01  2.8018E+01  9.7784E+01  2.2419E+01  2.1795E+01
 ETASHRINKVR(%)  2.1675E-01  4.8186E+01  9.9951E+01  3.9812E+01  3.8840E+01
 EBVSHRINKSD(%)  5.2435E-01  2.8654E+01  9.7975E+01  2.1468E+01  1.9819E+01
 EBVSHRINKVR(%)  1.0459E+00  4.9098E+01  9.9959E+01  3.8327E+01  3.5710E+01
 RELATIVEINF(%)  9.8372E+01  8.0726E+00  3.4966E-03  1.0769E+01  4.8190E+00
 EPSSHRINKSD(%)  4.2956E+01
 EPSSHRINKVR(%)  6.7459E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1634.4831131994963     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -899.33228663575812     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     9.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.84
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1634.483       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.52E-01  5.74E-01  7.30E-01  1.27E+00  6.69E-01  9.65E-01  2.03E+00  2.82E-02  6.42E-01  8.43E-01  1.07E+00
 


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
+        1.30E+03
 
 TH 2
+       -1.60E+01  4.44E+02
 
 TH 3
+        2.79E+01  2.19E+02  1.17E+03
 
 TH 4
+       -1.11E+01  3.74E+02 -4.27E+02  1.06E+03
 
 TH 5
+       -6.35E+00 -4.34E+02 -1.49E+03  4.36E+02  2.25E+03
 
 TH 6
+        1.27E+00 -3.35E+00  5.61E+00 -3.32E+00 -3.81E+00  2.09E+02
 
 TH 7
+        2.05E+00  3.81E+01 -6.32E+00 -1.54E+01  4.47E+00 -1.44E-01  1.87E+01
 
 TH 8
+       -8.36E-01 -2.50E+00 -5.64E+00  2.42E-01  1.77E+00 -3.39E-01  7.15E-02 -7.71E+00
 
 TH 9
+        3.16E+00 -2.38E+01 -3.92E+01 -1.54E+01  4.97E+01 -5.12E-01  1.33E+01  3.74E-01  1.95E+02
 
 TH10
+       -1.77E+00 -5.96E+00 -8.79E+01 -3.14E+01 -3.71E+01  4.04E+00  3.71E+00  4.27E+00  1.74E+01  1.09E+02
 
 TH11
+       -5.61E+00 -8.85E+00 -4.59E+01 -8.65E+00  1.64E+01  2.45E+00  1.15E+00  1.57E+00  1.90E+01  2.43E+01  1.91E+02
 
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
 #CPUT: Total CPU Time in Seconds,       15.365
Stop Time:
Sat Sep 18 12:24:06 CDT 2021
