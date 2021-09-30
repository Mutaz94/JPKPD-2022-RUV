Thu Sep 30 03:31:55 CDT 2021
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
$DATA ../../../../data/spa1/D/dat71.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m71.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   16012.1538332111        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.0770E+02  3.3687E+02 -9.6111E+01 -1.1722E+02  3.5917E+02 -2.4314E+03 -8.2847E+02 -1.1340E+02 -2.0330E+03 -6.7874E+02
            -2.9806E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -627.310903215426        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.8171E+00  9.0416E-01  9.7916E-01  2.6106E+00  1.3705E+00  3.8989E+00  1.4014E+00  9.7485E-01  3.7035E+00  1.2343E+00
             1.0918E+01
 PARAMETER:  6.9725E-01 -7.4531E-04  7.8942E-02  1.0596E+00  4.1516E-01  1.4607E+00  4.3750E-01  7.4526E-02  1.4093E+00  3.1051E-01
             2.4905E+00
 GRADIENT:   7.3242E+01  1.8610E+01 -6.4444E+01  7.4162E+01 -7.7034E+00  1.1502E+02  3.9345E+00  6.9628E+00  4.3855E+01  5.2917E+00
             1.0630E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -675.496642693993        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.6857E+00  9.7880E-01  3.0054E+00  2.8828E+00  6.6557E+00  3.2289E+00  5.6598E+00  7.3102E-01  4.0218E+00  8.7739E+00
             9.4076E+00
 PARAMETER:  6.2218E-01  7.8576E-02  1.2004E+00  1.1588E+00  1.9955E+00  1.2721E+00  1.8334E+00 -2.1332E-01  1.4917E+00  2.2718E+00
             2.3415E+00
 GRADIENT:   9.1250E+01  1.3196E+01  7.3014E-01  8.1583E+01 -6.5050E+00  5.8857E+01  1.4237E+01 -4.5035E-01  7.9933E+01  7.9118E+00
             6.7413E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -702.153132164446        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.2733E+00  1.0025E+00  2.3440E+00  1.7839E+00  5.5336E+00  2.8550E+00  4.4647E+00  1.3333E+00  2.8951E+00  1.0072E+01
             9.5317E+00
 PARAMETER:  3.4162E-01  1.0252E-01  9.5186E-01  6.7879E-01  1.8108E+00  1.1491E+00  1.5962E+00  3.8768E-01  1.1630E+00  2.4098E+00
             2.3546E+00
 GRADIENT:   2.0129E+01  1.9934E+01  1.8197E+01  9.7296E+00 -9.4765E+00  3.5586E+01  1.6376E+01 -4.4157E+00  1.2084E+01  1.8090E+01
             9.4600E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -742.680637178583        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0161E+00  1.1769E-01  9.2063E-01  1.8272E+00  7.7843E+00  2.2083E+00  5.1011E+00  3.0599E+00  1.2263E+00  1.4134E+01
             8.9908E+00
 PARAMETER:  1.1599E-01 -2.0397E+00  1.7299E-02  7.0278E-01  2.1521E+00  8.9224E-01  1.7295E+00  1.2184E+00  3.0399E-01  2.7486E+00
             2.2962E+00
 GRADIENT:  -1.7007E+01  9.1321E+00  2.1778E+01  1.0841E+02 -2.1710E+00 -9.6090E+00 -7.4511E-01 -3.4235E+01 -3.1165E+01  5.2254E+01
            -6.7646E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -854.081337096557        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      424
 NPARAMETR:  6.8004E-01  1.0000E-02  8.4381E-02  6.9859E-01  1.5335E+01  1.5086E+00  1.7632E+00  2.8708E+00  3.8931E-01  7.4921E+00
             8.3086E+00
 PARAMETER: -2.8560E-01 -4.8510E+00 -2.3724E+00 -2.5869E-01  2.8301E+00  5.1115E-01  6.6712E-01  1.1546E+00 -8.4337E-01  2.1139E+00
             2.2173E+00
 GRADIENT:   2.7922E+01  0.0000E+00 -1.0911E+01  1.9139E+01  1.7023E-01 -6.3119E+01  2.3067E-03  3.0682E+01  4.0074E+00 -1.3575E-02
            -4.7841E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -887.035010456892        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      603
 NPARAMETR:  6.1221E-01  1.0000E-02  4.3372E-02  4.7842E-01  9.6184E+00  2.2446E+00  1.2172E+00  1.6456E+00  2.1392E-01  1.9038E+00
             8.4177E+00
 PARAMETER: -3.9069E-01 -5.9775E+00 -3.0379E+00 -6.3726E-01  2.3637E+00  9.0851E-01  2.9652E-01  5.9810E-01 -1.4421E+00  7.4385E-01
             2.2303E+00
 GRADIENT:   3.7924E+01  0.0000E+00 -5.1273E+01  6.4978E+01  3.2668E-01  2.3140E+01  1.1549E-03  1.4847E+01  6.1506E-01  1.6534E-02
             2.3904E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -898.179698562362        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      779
 NPARAMETR:  4.5683E-01  1.0000E-02  2.5665E-02  3.0998E-01  1.0360E+01  2.0475E+00  8.4691E-01  1.2264E+00  1.4701E-01  7.1460E-01
             8.1997E+00
 PARAMETER: -6.8344E-01 -6.8242E+00 -3.5626E+00 -1.0713E+00  2.4380E+00  8.1660E-01 -6.6163E-02  3.0409E-01 -1.8173E+00 -2.3603E-01
             2.2041E+00
 GRADIENT:   1.4767E+00  0.0000E+00 -1.1111E+01  1.3837E+01 -2.9241E-02  4.6572E-01  5.1782E-03 -3.3861E+00 -7.3617E-02  2.9262E-03
            -4.3595E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -898.255963479269        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      967
 NPARAMETR:  4.5576E-01  1.0000E-02  2.5766E-02  3.0867E-01  1.0535E+01  2.0442E+00  6.6213E-01  1.2455E+00  1.4373E-01  5.6139E-01
             8.2146E+00
 PARAMETER: -6.8580E-01 -6.8242E+00 -3.5587E+00 -1.0755E+00  2.4547E+00  8.1499E-01 -3.1230E-01  3.1957E-01 -1.8398E+00 -4.7735E-01
             2.2059E+00
 GRADIENT:   6.3533E-01  0.0000E+00 -2.0652E+00  1.2813E+00 -1.8826E-02  4.3482E-01  3.0802E-03 -6.4303E-01  1.3720E-02  1.3670E-03
            -4.6548E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -898.287932643326        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1157             RESET HESSIAN, TYPE I
 NPARAMETR:  4.5583E-01  1.0000E-02  2.5567E-02  3.0774E-01  1.0473E+01  2.0458E+00  1.2763E-01  1.2630E+00  3.0051E-02  1.6287E-01
             8.2155E+00
 PARAMETER: -6.8564E-01 -6.8242E+00 -3.5665E+00 -1.0785E+00  2.4488E+00  8.1577E-01 -1.9586E+00  3.3351E-01 -3.4049E+00 -1.7148E+00
             2.2060E+00
 GRADIENT:   7.5430E+01  0.0000E+00  1.0531E+02  3.8582E+01  4.4813E-02  3.0221E+01  1.5423E-04  1.2674E+00  1.0466E-02  1.5827E-04
             2.1374E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -898.308384713412        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1319
 NPARAMETR:  4.5133E-01  1.0000E-02  2.5391E-02  3.0592E-01  1.0434E+01  2.0351E+00  8.1524E-02  1.2642E+00  1.0000E-02  5.7560E-02
             8.2089E+00
 PARAMETER: -6.9556E-01 -6.8242E+00 -3.5734E+00 -1.0844E+00  2.4451E+00  8.1053E-01 -2.4069E+00  3.3443E-01 -5.0922E+00 -2.7549E+00
             2.2052E+00
 GRADIENT:  -1.8365E+00  0.0000E+00 -5.1830E+00  5.9114E+00 -1.3137E-02 -9.4770E-01  6.2809E-05  7.8091E-02  0.0000E+00  2.0228E-05
            -6.0864E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -898.335381411606        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     1510
 NPARAMETR:  4.5307E-01  1.0000E-02  2.5274E-02  3.0495E-01  1.0478E+01  2.0431E+00  4.8833E-02  1.2600E+00  1.0000E-02  4.4168E-02
             8.2154E+00
 PARAMETER: -6.9172E-01 -6.8242E+00 -3.5780E+00 -1.0876E+00  2.4493E+00  8.1446E-01 -2.9194E+00  3.3114E-01 -5.0922E+00 -3.0197E+00
             2.2060E+00
 GRADIENT:   8.7019E-01  0.0000E+00 -5.8236E+00  5.3368E+00  5.8983E-03  4.3663E-01  2.2966E-05 -1.4621E-01  0.0000E+00  1.1197E-05
            -3.3613E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -898.385314258648        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     1705             RESET HESSIAN, TYPE I
 NPARAMETR:  4.5237E-01  1.0000E-02  2.4831E-02  3.0120E-01  1.0401E+01  2.0492E+00  1.9707E-02  1.2581E+00  1.0000E-02  2.2160E-02
             8.2152E+00
 PARAMETER: -6.9326E-01 -6.8242E+00 -3.5956E+00 -1.1000E+00  2.4419E+00  8.1743E-01 -3.8268E+00  3.2961E-01 -5.0922E+00 -3.7095E+00
             2.2060E+00
 GRADIENT:   7.8133E+01  0.0000E+00  1.0638E+02  3.8992E+01  3.0289E-02  3.1186E+01  5.5392E-06  1.8207E+00  0.0000E+00  4.2486E-06
             2.1650E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -898.397362890045        NO. OF FUNC. EVALS.: 171
 CUMULATIVE NO. OF FUNC. EVALS.:     1876
 NPARAMETR:  4.4863E-01  1.0000E-02  2.4886E-02  3.0053E-01  1.0667E+01  2.0385E+00  1.9079E-02  1.2537E+00  1.0000E-02  2.1974E-02
             8.2098E+00
 PARAMETER: -7.0155E-01 -6.8242E+00 -3.5935E+00 -1.1022E+00  2.4672E+00  8.1223E-01 -3.8592E+00  3.2610E-01 -5.0922E+00 -3.7179E+00
             2.2053E+00
 GRADIENT:  -6.4423E-01  0.0000E+00 -2.4786E+00  1.2320E+00  5.2171E-03 -1.1477E-01  3.9103E-06 -2.2428E-01  0.0000E+00  2.1961E-06
            -1.6688E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1876
 NO. OF SIG. DIGITS IN FINAL EST.:  2.4
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.5556E-03 -8.1814E-06  1.0543E-02 -3.7280E-04  5.0658E-06
 SE:             2.9319E-02  1.6233E-06  2.4603E-02  2.7494E-04  3.1985E-06
 N:                     100         100         100         100         100

 P VAL.:         7.9664E-01  4.6633E-07  6.6828E-01  1.7513E-01  1.1324E-01

 ETASHRINKSD(%)  1.7773E+00  9.9995E+01  1.7576E+01  9.9079E+01  9.9989E+01
 ETASHRINKVR(%)  3.5230E+00  1.0000E+02  3.2063E+01  9.9992E+01  1.0000E+02
 EBVSHRINKSD(%)  2.0505E+00  9.9986E+01  1.7319E+01  9.9106E+01  9.9986E+01
 EBVSHRINKVR(%)  4.0590E+00  1.0000E+02  3.1639E+01  9.9992E+01  1.0000E+02
 RELATIVEINF(%)  2.0261E+01  6.2707E-07  1.4584E+00  1.5167E-04  4.5322E-07
 EPSSHRINKSD(%)  1.1980E+01
 EPSSHRINKVR(%)  2.2525E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -898.39736289004463     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       20.541170314628062     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    38.39
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.18
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -898.397       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.49E-01  1.00E-02  2.49E-02  3.01E-01  1.07E+01  2.04E+00  1.91E-02  1.25E+00  1.00E-02  2.20E-02  8.21E+00
 


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
+        1.26E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.84E+03  0.00E+00  7.55E+05
 
 TH 4
+       -6.12E+02  0.00E+00 -7.99E+04  9.39E+03
 
 TH 5
+        2.04E-01  0.00E+00 -7.38E+00  7.64E-01  7.25E-03
 
 TH 6
+        4.91E+00  0.00E+00  1.52E+02 -3.40E+01  7.76E-03  4.04E+01
 
 TH 7
+       -3.77E-03  0.00E+00  1.64E-02  3.99E-03 -3.10E-04 -2.21E-03 -7.82E-02
 
 TH 8
+        5.17E+00  0.00E+00 -2.43E+01 -7.94E+01  2.29E-03  3.01E+00 -1.91E-04  5.39E+01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -1.01E-02  0.00E+00  2.48E-03  1.32E-02  5.66E-05  4.39E-04  2.47E-02  1.92E-02  0.00E+00 -5.95E-02
 
 TH11
+       -1.28E+01  0.00E+00  2.39E+02 -2.52E+01 -2.56E-03  7.54E-01  3.94E-04  3.07E+00  0.00E+00 -9.48E-05  7.42E+00
 
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
 #CPUT: Total CPU Time in Seconds,       48.638
Stop Time:
Thu Sep 30 03:32:45 CDT 2021
