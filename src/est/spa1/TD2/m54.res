Thu Sep 30 02:12:13 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat54.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1670.24012658992        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8058E+02  5.3497E+01  3.8301E+01  1.0553E+02  1.0805E+01  3.2785E+01  2.0866E+00 -1.8409E+02 -3.1451E+01  1.5534E+01
            -7.2763E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2090.44003653193        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:      145
 NPARAMETR:  9.9692E-01  9.7415E-01  1.1018E+00  9.3469E-01  1.0267E+00  1.0708E+00  1.0452E+00  1.1578E+00  1.0258E+00  6.6270E-01
             1.6832E+00
 PARAMETER:  9.6918E-02  7.3810E-02  1.9690E-01  3.2457E-02  1.2639E-01  1.6837E-01  1.4417E-01  2.4652E-01  1.2547E-01 -3.1143E-01
             6.2071E-01
 GRADIENT:   6.8529E+01 -4.8471E+01 -7.9921E+00 -5.2849E+01 -6.4399E+00  4.3211E+01 -2.8454E+00  1.3205E+01  5.1789E+00  2.3069E+00
             2.9338E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2095.46604746690        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      324
 NPARAMETR:  1.0632E+00  9.7411E-01  1.1958E+00  9.3469E-01  1.0623E+00  1.0823E+00  1.1548E+00  1.1579E+00  1.0258E+00  6.2557E-01
             1.6826E+00
 PARAMETER:  1.6124E-01  7.3770E-02  2.7886E-01  3.2458E-02  1.6039E-01  1.7905E-01  2.4392E-01  2.4662E-01  1.2552E-01 -3.6909E-01
             6.2034E-01
 GRADIENT:   2.1494E+00 -6.2632E+01  1.7667E+00 -9.4143E+01 -9.0616E+00 -2.2124E+00  4.7477E-01  9.6185E+00  5.3032E+00 -1.0600E+00
             2.8178E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2096.06437911701        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      500
 NPARAMETR:  1.0607E+00  9.7415E-01  1.4317E+00  9.3469E-01  1.1515E+00  1.0815E+00  9.5000E-01  1.1578E+00  1.0258E+00  8.2225E-01
             1.6808E+00
 PARAMETER:  1.5895E-01  7.3808E-02  4.5887E-01  3.2463E-02  2.4104E-01  1.7839E-01  4.8708E-02  2.4655E-01  1.2549E-01 -9.5709E-02
             6.1929E-01
 GRADIENT:  -8.3236E-01 -7.7038E+01  1.8603E-01 -1.2034E+02  1.3721E+00 -1.6434E+00 -1.1228E-01  5.6621E+00 -1.2091E+01  6.7977E-02
             2.8606E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2135.52798763984        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      687
 NPARAMETR:  1.0800E+00  9.8399E-01  2.7785E+00  9.3630E-01  1.2378E+00  1.2158E+00  1.0000E-02  1.1391E+00  1.0190E+00  1.7825E+00
             1.2447E+00
 PARAMETER:  1.7692E-01  8.3861E-02  1.1219E+00  3.4177E-02  3.1332E-01  2.9537E-01 -4.6416E+01  2.3023E-01  1.1879E-01  6.7803E-01
             3.1890E-01
 GRADIENT:   3.4194E+01 -3.3305E+01  4.1753E+01 -1.4038E+02 -9.5142E+01  4.0273E+01  0.0000E+00 -1.0442E+01 -7.6119E+01  4.7907E+01
             1.6098E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2170.81362375966        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:      880             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0604E+00  9.8807E-01  2.0270E+00  9.3695E-01  1.1529E+00  1.1280E+00  1.0000E-02  1.1314E+00  1.0174E+00  1.3575E+00
             1.1010E+00
 PARAMETER:  1.5863E-01  8.7997E-02  8.0655E-01  3.4876E-02  2.4229E-01  2.2047E-01 -6.4919E+01  2.2348E-01  1.1726E-01  4.0561E-01
             1.9623E-01
 GRADIENT:   7.3183E+02  5.9553E+01  5.9063E+01 -2.9004E+01 -1.1503E+02  1.6881E+02  0.0000E+00 -1.0584E+01 -7.0919E+01  1.7882E+01
             1.0399E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2173.13444418126        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:     1011
 NPARAMETR:  1.0604E+00  9.8807E-01  2.0270E+00  9.3697E-01  1.1530E+00  1.0807E+00  1.0000E-02  1.1314E+00  1.0470E+00  1.3575E+00
             1.1011E+00
 PARAMETER:  1.5865E-01  8.7997E-02  8.0655E-01  3.4893E-02  2.4233E-01  1.7757E-01 -6.4919E+01  2.2348E-01  1.4590E-01  4.0561E-01
             1.9626E-01
 GRADIENT:   7.3253E+02  5.4383E+01  5.9083E+01 -2.2665E+01 -1.1502E+02  1.2129E+02  0.0000E+00 -1.0638E+01 -5.3644E+01  1.7948E+01
             1.0462E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2177.81101463971        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1191
 NPARAMETR:  1.0604E+00  9.8807E-01  2.0270E+00  9.3696E-01  1.1530E+00  1.0901E+00  1.0000E-02  1.1314E+00  1.1686E+00  1.3575E+00
             1.1010E+00
 PARAMETER:  1.5864E-01  8.7997E-02  8.0654E-01  3.4890E-02  2.4233E-01  1.8624E-01 -6.4919E+01  2.2348E-01  2.5582E-01  4.0561E-01
             1.9626E-01
 GRADIENT:   9.6748E+00 -1.4534E+01  5.6135E+01 -8.1764E+01 -1.2581E+02  9.5010E-02  0.0000E+00 -1.1052E+01 -1.8416E+01  1.5992E+01
             1.0510E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2184.71925529323        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1351
 NPARAMETR:  1.0597E+00  9.8766E-01  2.0269E+00  9.3696E-01  1.4415E+00  1.0420E+00  1.0000E-02  2.0853E+00  1.1806E+00  1.3575E+00
             8.4277E-01
 PARAMETER:  1.5798E-01  8.7581E-02  8.0653E-01  3.4890E-02  4.6569E-01  1.4111E-01 -6.4919E+01  8.3491E-01  2.6601E-01  4.0561E-01
            -7.1063E-02
 GRADIENT:   1.2425E+03  1.4955E+01 -1.1349E+01  5.1623E+01  9.9004E+01  1.4888E+02  0.0000E+00  2.5304E+01  2.6439E+01  6.2171E+00
            -9.7661E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2193.32271165046        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1422
 NPARAMETR:  1.0596E+00  9.8762E-01  2.0269E+00  9.3696E-01  1.3601E+00  1.0041E+00  1.0000E-02  1.7936E+00  1.1806E+00  1.3575E+00
             8.9143E-01
 PARAMETER:  1.5792E-01  8.7545E-02  8.0653E-01  3.4890E-02  4.0759E-01  1.0412E-01 -6.4919E+01  6.8423E-01  2.6601E-01  4.0561E-01
            -1.4927E-02
 GRADIENT:   1.1108E+03  1.8834E+01  3.6448E+00  3.8633E+01  4.5197E+01  7.4177E+01  0.0000E+00  1.5290E+01  2.2588E+01  1.4710E+01
            -4.4301E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2195.38107615506        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1603
 NPARAMETR:  1.0587E+00  9.8904E-01  2.0186E+00  9.3648E-01  1.3008E+00  1.1538E+00  1.0000E-02  1.7999E+00  1.2807E+00  1.2196E+00
             8.9097E-01
 PARAMETER:  1.5702E-01  8.8982E-02  8.0240E-01  3.4378E-02  3.6297E-01  2.4304E-01 -6.4919E+01  6.8773E-01  3.4742E-01  2.9856E-01
            -1.5439E-02
 GRADIENT:   1.1015E+03  1.5170E+01  1.3625E+01  5.2284E+01  1.8432E+01  2.8346E+02  0.0000E+00  1.4476E+01  6.6144E+01  1.1403E+00
            -4.6854E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2196.14236884957        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1780
 NPARAMETR:  1.0587E+00  9.8903E-01  2.0184E+00  9.3647E-01  1.3108E+00  1.0933E+00  1.0000E-02  1.8001E+00  1.2368E+00  1.2394E+00
             8.9096E-01
 PARAMETER:  1.5700E-01  8.8970E-02  8.0229E-01  3.4365E-02  3.7062E-01  1.8922E-01 -6.4919E+01  6.8782E-01  3.1252E-01  3.1465E-01
            -1.5452E-02
 GRADIENT:   1.1958E+01 -5.7236E+01  5.7438E+00 -7.3609E+01  8.1245E-01  1.0489E+00  0.0000E+00  1.2521E+01  9.2608E-02 -2.6038E-01
            -4.7604E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2196.71344624062        NO. OF FUNC. EVALS.: 217
 CUMULATIVE NO. OF FUNC. EVALS.:     1997             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0585E+00  9.9284E-01  2.0154E+00  9.3931E-01  1.3107E+00  1.0934E+00  1.0000E-02  1.7943E+00  1.2368E+00  1.2395E+00
             8.9298E-01
 PARAMETER:  1.5688E-01  9.2812E-02  8.0083E-01  3.7394E-02  3.7059E-01  1.8932E-01 -6.4919E+01  6.8464E-01  3.1253E-01  3.1469E-01
            -1.3187E-02
 GRADIENT:   1.0956E+03  2.5696E+01  1.1954E+01  5.4076E+01  2.3460E+01  2.0335E+02  0.0000E+00  1.4112E+01  4.7262E+01  3.1416E+00
            -4.4582E+01

0ITERATION NO.:   62    OBJECTIVE VALUE:  -2196.71344624062        NO. OF FUNC. EVALS.:  69
 CUMULATIVE NO. OF FUNC. EVALS.:     2066
 NPARAMETR:  1.0602E+00  9.9185E-01  1.9994E+00  9.3838E-01  1.3101E+00  1.0955E+00  1.0000E-02  1.8067E+00  1.2329E+00  1.2421E+00
             8.9209E-01
 PARAMETER:  1.5688E-01  9.2812E-02  8.0083E-01  3.7394E-02  3.7059E-01  1.8932E-01 -6.4919E+01  6.8464E-01  3.1253E-01  3.1469E-01
            -1.3187E-02
 GRADIENT:  -1.7739E+05  2.7823E+05  1.7388E+04  2.7822E+05  2.5688E-01 -7.3501E+04  0.0000E+00 -4.0692E+04  8.9048E+04 -1.8466E-01
             2.7776E+05

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2066
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.1632E-03 -5.5004E-04 -4.6584E-02  3.3616E-02 -4.2333E-02
 SE:             3.0007E-02  1.9810E-04  1.4407E-02  2.9433E-02  2.2451E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6908E-01  5.4936E-03  1.2233E-03  2.5341E-01  5.9358E-02

 ETASHRINKSD(%)  1.0000E-10  9.9336E+01  5.1734E+01  1.3947E+00  2.4785E+01
 ETASHRINKVR(%)  1.0000E-10  9.9996E+01  7.6704E+01  2.7700E+00  4.3427E+01
 EBVSHRINKSD(%)  2.2009E-01  9.9376E+01  5.4926E+01  1.1493E+00  1.8958E+01
 EBVSHRINKVR(%)  4.3969E-01  9.9996E+01  7.9683E+01  2.2854E+00  3.4322E+01
 RELATIVEINF(%)  9.9408E+01  4.7956E-04  1.0146E+01  1.2796E+01  3.3228E+01
 EPSSHRINKSD(%)  3.1011E+01
 EPSSHRINKVR(%)  5.2405E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2196.7134462406166     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1277.7749130359439     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.01
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2196.713       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  9.93E-01  2.02E+00  9.39E-01  1.31E+00  1.09E+00  1.00E-02  1.79E+00  1.24E+00  1.24E+00  8.93E-01
 


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
+        2.52E+07
 
 TH 2
+        9.79E+02  7.06E+07
 
 TH 3
+        4.07E+01 -5.60E+02  2.02E+01
 
 TH 4
+        7.56E+02 -9.74E+03 -5.19E+02  7.88E+07
 
 TH 5
+       -1.30E+03  2.14E+03  9.79E+01  2.37E+03  2.95E+06
 
 TH 6
+       -8.19E+02  1.27E+03  5.21E+01  9.57E+02 -1.06E+03  1.62E+07
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+       -3.87E+01  6.62E+02 -1.54E+05  4.14E+02 -1.84E+02 -4.92E+01  0.00E+00  4.62E+05
 
 TH 9
+        2.47E+02  1.81E+07  5.11E+01  8.70E+02  5.75E+02  3.15E+02  0.00E+00 -8.56E+01  4.66E+06
 
 TH10
+        1.32E+03 -2.31E+03 -1.47E+02 -2.47E+03 -3.67E+06  1.09E+03  0.00E+00  1.95E+02 -5.95E+02  4.57E+06
 
 TH11
+        7.47E+02 -1.05E+04  4.82E+06 -9.38E+03  2.51E+03  9.50E+02  0.00E+00  2.22E+04  9.22E+02 -2.59E+03  8.69E+07
 
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
 #CPUT: Total CPU Time in Seconds,       44.325
Stop Time:
Thu Sep 30 02:12:59 CDT 2021
