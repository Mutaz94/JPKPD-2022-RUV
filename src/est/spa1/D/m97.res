Thu Sep 30 03:49:55 CDT 2021
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
$DATA ../../../../data/spa1/D/dat97.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:   31508.9449452360        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1114E+03  7.6303E+02 -1.6757E+01  7.5754E+02  1.9062E+02 -3.1502E+03 -1.3698E+03 -4.6244E+01 -2.0500E+03 -7.6183E+02
            -5.8992E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -356.764229588960        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.8562E-01  9.6950E-01  9.0932E-01  1.2941E+00  1.3912E+00  2.0472E+00  1.0950E+00  9.6095E-01  1.0788E+00  9.6685E-01
             1.4540E+01
 PARAMETER:  8.5520E-02  6.9025E-02  4.9379E-03  3.5778E-01  4.3019E-01  8.1647E-01  1.9075E-01  6.0164E-02  1.7589E-01  6.6288E-02
             2.7769E+00
 GRADIENT:  -5.4439E+01  1.7684E+01 -2.0382E+00  1.4692E+01 -7.4506E+00  4.3748E+01 -1.7432E+00  3.8669E+00 -5.2864E+00  1.7057E+00
             1.3170E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -376.060265612625        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0754E+00  8.0381E-01  9.3308E-01  1.4591E+00  6.3857E+00  1.6501E+00  1.9232E-01  2.1461E-01  1.6232E+00  1.4840E+00
             1.5411E+01
 PARAMETER:  1.7268E-01 -1.1839E-01  3.0740E-02  4.7784E-01  1.9541E+00  6.0085E-01 -1.5486E+00 -1.4389E+00  5.8440E-01  4.9472E-01
             2.8351E+00
 GRADIENT:  -3.8781E+01  1.7638E+01  9.9786E+00  2.3671E+01 -4.9333E+00 -4.6154E+01  1.1970E-01  5.4980E-02  1.1439E+01  1.5368E-01
             5.7897E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -429.082626225504        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.2437E+00  2.7118E-01  2.6909E-01  1.5862E+00  1.5360E+01  2.5886E+00  2.1306E-01  1.0000E-02  1.2379E+00  4.5412E+00
             1.3855E+01
 PARAMETER:  3.1808E-01 -1.2050E+00 -1.2127E+00  5.6137E-01  2.8318E+00  1.0511E+00 -1.4462E+00 -8.8292E+00  3.1343E-01  1.6132E+00
             2.7286E+00
 GRADIENT:   8.1347E+01  4.6583E+01 -5.0923E+01  1.4424E+02 -4.6928E+00 -7.5957E+00  5.4030E-01  0.0000E+00 -5.3970E+01  1.7633E+00
            -5.2204E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -464.396599458783        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.8610E-01  1.1673E-01  1.5004E-01  1.0878E+00  6.5267E+01  2.3241E+00  1.0000E-02  1.0000E-02  9.0095E-01  1.6315E+01
             1.4240E+01
 PARAMETER:  8.5999E-02 -2.0479E+00 -1.7968E+00  1.8418E-01  4.2785E+00  9.4331E-01 -4.7878E+00 -1.1569E+01 -4.3027E-03  2.8921E+00
             2.7561E+00
 GRADIENT:   6.0489E+01  3.4426E+00  3.2745E+00  7.3903E+01 -1.1169E+00 -1.2939E+01  0.0000E+00  0.0000E+00 -3.3726E+01  4.7539E+00
             8.6064E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -506.710376123534        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      382
 NPARAMETR:  5.2284E-01  2.3174E-02  3.1514E-02  4.1298E-01  1.0356E+03  1.9559E+00  1.0000E-02  1.0000E-02  3.1734E-01  5.9570E+01
             1.4506E+01
 PARAMETER: -5.4849E-01 -3.6647E+00 -3.3573E+00 -7.8436E-01  7.0428E+00  7.7085E-01 -9.6905E+00 -1.8479E+01 -1.0478E+00  4.1871E+00
             2.7746E+00
 GRADIENT:   5.0276E+01 -3.6606E+00 -2.0459E+01  1.1282E+02  1.4351E-03 -3.9924E+01  0.0000E+00  0.0000E+00  3.2912E+00  6.7528E-03
             5.1441E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -519.168646983912        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      565
 NPARAMETR:  3.9731E-01  1.3999E-02  2.1527E-02  2.9172E-01  2.4020E+03  2.2186E+00  1.0000E-02  1.0000E-02  1.9607E-01  9.0501E+01
             1.3103E+01
 PARAMETER: -8.2304E-01 -4.1687E+00 -3.7384E+00 -1.1320E+00  7.8840E+00  8.9687E-01 -9.1272E+00 -2.0997E+01 -1.5293E+00  4.6054E+00
             2.6728E+00
 GRADIENT:  -1.9817E+01  9.2977E-01 -1.4046E+00  1.8159E+01 -1.8898E-04  2.0119E+00  0.0000E+00  0.0000E+00  8.3839E-01  3.4111E-04
            -1.2551E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -522.619464901899        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      750
 NPARAMETR:  3.6056E-01  1.0000E-02  1.4617E-02  2.1212E-01  1.8763E+04  2.1511E+00  1.0000E-02  1.0000E-02  4.2718E-02  5.6952E+01
             1.3273E+01
 PARAMETER: -9.2011E-01 -4.5712E+00 -4.1256E+00 -1.4506E+00  9.9396E+00  8.6599E-01 -1.0328E+01 -2.3097E+01 -3.0531E+00  4.1422E+00
             2.6857E+00
 GRADIENT:   3.3786E+00  0.0000E+00 -4.4222E+00  1.0152E+00  6.6738E-06 -1.8392E+00  0.0000E+00  0.0000E+00  4.7358E-02  8.1655E-07
            -2.4541E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -522.649809861815        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      938             RESET HESSIAN, TYPE I
 NPARAMETR:  3.5924E-01  1.0000E-02  1.4627E-02  2.1273E-01  4.2043E+04  2.1656E+00  1.0000E-02  1.0000E-02  1.1290E-02  4.9351E+01
             1.3285E+01
 PARAMETER: -9.2378E-01 -4.5712E+00 -4.1249E+00 -1.4477E+00  1.0746E+01  8.7268E-01 -1.0328E+01 -2.3097E+01 -4.3838E+00  3.9990E+00
             2.6867E+00
 GRADIENT:   5.6799E+01  0.0000E+00  8.0690E+01  2.6475E+01  3.0097E-06  1.8957E+01  0.0000E+00  0.0000E+00  4.2230E-03  1.8243E-07
             2.5648E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -522.653723321821        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:     1110
 NPARAMETR:  3.5851E-01  1.0000E-02  1.4634E-02  2.1255E-01  2.9341E+03  2.1644E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.4013E+01
             1.3288E+01
 PARAMETER: -9.2580E-01 -4.5712E+00 -4.1244E+00 -1.4486E+00  8.0842E+00  8.7213E-01 -1.0328E+01 -2.3097E+01 -4.6701E+00  3.6267E+00
             2.6869E+00
 GRADIENT:   5.5859E+01  0.0000E+00  8.2848E+01  2.4195E+01  1.4449E-05  1.8835E+01  0.0000E+00  0.0000E+00  0.0000E+00  1.7693E-05
             2.6287E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -522.655393863262        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1291
 NPARAMETR:  3.5841E-01  1.0000E-02  1.4623E-02  2.1242E-01  1.0243E+03  2.1646E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.2093E+01
             1.3289E+01
 PARAMETER: -9.2608E-01 -4.5712E+00 -4.1251E+00 -1.4492E+00  7.0318E+00  8.7226E-01 -1.0328E+01 -2.3097E+01 -4.6701E+00  3.1952E+00
             2.6869E+00
 GRADIENT:   4.8180E-01  0.0000E+00 -3.9823E+00  1.7954E+00 -1.5339E-05  3.6911E-01  0.0000E+00  0.0000E+00  0.0000E+00  4.4124E-05
            -5.3691E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -522.664542398332        NO. OF FUNC. EVALS.: 218
 CUMULATIVE NO. OF FUNC. EVALS.:     1509             RESET HESSIAN, TYPE I
 NPARAMETR:  3.5809E-01  1.0000E-02  1.4567E-02  2.1142E-01  1.3305E+03  2.1663E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.4696E-01
             1.3297E+01
 PARAMETER: -9.2698E-01 -4.5712E+00 -4.1290E+00 -1.4539E+00  7.2933E+00  8.7304E-01 -1.0328E+01 -2.3097E+01 -4.6701E+00 -7.0528E-01
             2.6875E+00
 GRADIENT:   5.6152E+01  0.0000E+00  8.5332E+01  2.1121E+01 -2.8220E-05  1.9391E+01  0.0000E+00  0.0000E+00  0.0000E+00  1.1267E-08
             2.7019E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -522.670874262817        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1707             RESET HESSIAN, TYPE I
 NPARAMETR:  3.5813E-01  1.0000E-02  1.4501E-02  2.1109E-01  4.6048E+03  2.1676E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.3754E-01
             1.3291E+01
 PARAMETER: -9.2687E-01 -4.5712E+00 -4.1336E+00 -1.4555E+00  8.5349E+00  8.7360E-01 -1.0328E+01 -2.3097E+01 -4.6701E+00 -7.2658E-01
             2.6871E+00
 GRADIENT:   5.7154E+01  0.0000E+00  8.1901E+01  2.5201E+01  9.3032E-06  1.9491E+01  0.0000E+00  0.0000E+00  0.0000E+00  9.3881E-10
             2.6248E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -522.677737926500        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1884
 NPARAMETR:  3.5706E-01  1.0000E-02  1.4454E-02  2.1052E-01  4.0455E+06  2.1619E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.5044E-01
             1.3293E+01
 PARAMETER: -9.2986E-01 -4.5712E+00 -4.1368E+00 -1.4582E+00  1.5313E+01  8.7097E-01 -1.0328E+01 -2.3097E+01 -4.6701E+00 -6.9752E-01
             2.6873E+00
 GRADIENT:   8.4430E-01  0.0000E+00 -5.2986E+00  2.8924E+00 -3.5087E-09  1.1834E-02  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -2.5964E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1884
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6811E-03  6.2847E-06  5.7456E-05 -1.8307E-04 -3.8540E-11
 SE:             2.9374E-02  3.7020E-06  2.9304E-04  3.2269E-04  9.8345E-11
 N:                     100         100         100         100         100

 P VAL.:         9.5436E-01  8.9573E-02  8.4456E-01  5.7050E-01  6.9514E-01

 ETASHRINKSD(%)  1.5943E+00  9.9988E+01  9.9018E+01  9.8919E+01  1.0000E+02
 ETASHRINKVR(%)  3.1632E+00  1.0000E+02  9.9990E+01  9.9988E+01  1.0000E+02
 EBVSHRINKSD(%)  1.9727E+00  9.9971E+01  9.9051E+01  9.8881E+01  1.0000E+02
 EBVSHRINKVR(%)  3.9066E+00  1.0000E+02  9.9991E+01  9.9987E+01  1.0000E+02
 RELATIVEINF(%)  1.2803E+01  5.6301E-06  5.4860E-05  7.7905E-05  0.0000E+00
 EPSSHRINKSD(%)  4.6011E+00
 EPSSHRINKVR(%)  8.9905E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -522.67773792649996     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       396.26079527817274     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    32.13
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.52
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -522.678       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.57E-01  1.00E-02  1.45E-02  2.11E-01  4.05E+06  2.16E+00  1.00E-02  1.00E-02  1.00E-02  4.50E-01  1.33E+01
 


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
+        1.78E+03
 
 TH 2
+        0.00E+00  2.69E+04
 
 TH 3
+       -1.93E+04  0.00E+00  3.25E+06
 
 TH 4
+        1.67E+02  0.00E+00 -2.54E+05  2.13E+04
 
 TH 5
+        2.90E-12  0.00E+00 -6.62E-11  5.55E-12 -6.71E-21
 
 TH 6
+        1.14E+00  0.00E+00  9.77E+02 -9.82E+01 -4.55E-13  3.77E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -2.50E-04  0.00E+00 -1.14E-03  1.30E-04 -2.86E-12  5.69E-05  0.00E+00  0.00E+00  0.00E+00  3.08E-04
 
 TH11
+       -1.60E+01  0.00E+00  3.61E+02 -2.18E+01 -2.41E-14  5.35E-01  0.00E+00  0.00E+00  0.00E+00 -1.83E-06  2.40E+00
 
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
 #CPUT: Total CPU Time in Seconds,       39.715
Stop Time:
Thu Sep 30 03:50:36 CDT 2021
