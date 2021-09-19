Sat Sep 18 15:15:37 CDT 2021
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
$DATA ../../../../data/spa/D/dat30.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m30.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   10408.9658070723        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.0445E+02  4.3039E+01 -4.5132E+01 -1.0898E+02  1.6746E+02 -1.3022E+03 -5.5509E+02 -1.1732E+02 -1.0288E+03 -3.0154E+02
            -2.0950E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -720.140822760423        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.7007E+00  1.1170E+00  1.0285E+00  1.6027E+00  1.1870E+00  2.3510E+00  1.3239E+00  9.6869E-01  1.5200E+00  1.0417E+00
             1.3691E+01
 PARAMETER:  6.3102E-01  2.1065E-01  1.2813E-01  5.7170E-01  2.7142E-01  9.5484E-01  3.8062E-01  6.8188E-02  5.1871E-01  1.4087E-01
             2.7167E+00
 GRADIENT:   5.4921E+01 -8.3348E-01 -6.2184E+00 -1.5930E+01 -5.4710E+00  6.5613E+01 -1.2524E+00  4.6002E+00  6.5864E+00  3.1470E+00
             1.6606E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -737.222247032377        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.6598E+00  5.9783E-01  1.5239E+00  2.4594E+00  1.6012E+00  2.1960E+00  6.7320E+00  5.5250E-01  1.4705E+00  1.9868E+00
             1.2215E+01
 PARAMETER:  6.0670E-01 -4.1445E-01  5.2124E-01  9.9994E-01  5.7075E-01  8.8663E-01  2.0069E+00 -4.9330E-01  4.8558E-01  7.8653E-01
             2.6027E+00
 GRADIENT:   5.4546E+01  1.9865E+01 -4.4090E+00  7.4762E+01 -7.2201E+00  2.1635E+01  9.5169E+00  8.6340E-01  1.7571E+00  3.3818E+00
             9.7183E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -775.572020232282        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.1609E+00  4.9152E-01  1.4731E+00  1.5847E+00  2.1254E+00  1.7276E+00  3.8780E+00  1.7669E-01  1.0534E+00  2.2446E+00
             1.0203E+01
 PARAMETER:  2.4922E-01 -6.1024E-01  4.8740E-01  5.6042E-01  8.5394E-01  6.4675E-01  1.4553E+00 -1.6333E+00  1.5206E-01  9.0853E-01
             2.4227E+00
 GRADIENT:  -4.5514E+01  2.2517E+00  9.7962E+00 -5.4562E+01 -9.6436E+00  1.0152E+01  8.6083E+00  3.4816E-02  2.7683E+00  4.7725E-01
             6.3701E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -801.536932936702        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.1619E+00  3.9945E-01  5.3247E-01  1.4699E+00  1.9377E+01  1.5313E+00  2.7947E+00  1.0000E-02  1.0461E+00  2.2275E+00
             9.0725E+00
 PARAMETER:  2.5004E-01 -8.1767E-01 -5.3023E-01  4.8517E-01  3.0641E+00  5.2608E-01  1.1277E+00 -5.8173E+00  1.4507E-01  9.0089E-01
             2.3052E+00
 GRADIENT:   2.8829E+01  2.9982E+01 -5.9446E+00  4.1396E+01 -1.1903E+00 -4.7304E+01  3.6853E+00  0.0000E+00 -2.1883E+01 -1.3477E-02
            -4.1445E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -854.447945836170        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  7.2431E-01  5.1407E-02  9.3886E-02  7.7587E-01  1.1546E+03  2.1224E+00  3.9053E-01  1.0000E-02  1.0239E+00  2.4416E+00
             7.9510E+00
 PARAMETER: -2.2254E-01 -2.8680E+00 -2.2657E+00 -1.5377E-01  7.1515E+00  8.5254E-01 -8.4026E-01 -1.6236E+01  1.2364E-01  9.9267E-01
             2.1733E+00
 GRADIENT:  -2.1581E+01  2.3483E-01 -4.0343E+01  1.0369E+02  5.5794E-03  3.7097E+01  5.2822E-02  0.0000E+00 -2.5533E+01  6.5624E-07
            -2.5488E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -866.033379576739        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  5.9261E-01  2.5571E-02  4.1460E-02  4.2428E-01  1.1008E+04  1.7301E+00  7.9787E-02  1.0000E-02  1.0364E+00  4.0534E+00
             8.0051E+00
 PARAMETER: -4.2322E-01 -3.5663E+00 -3.0830E+00 -7.5735E-01  9.4064E+00  6.4820E-01 -2.4284E+00 -2.0653E+01  1.3573E-01  1.4996E+00
             2.1801E+00
 GRADIENT:   2.2986E+01  2.2275E+00 -1.8283E+01  2.4143E+01  5.9015E-05 -7.3612E+00  2.1707E-03  0.0000E+00 -4.1443E+00  3.5329E-08
            -7.2013E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -866.757770454362        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      598
 NPARAMETR:  5.8044E-01  2.4201E-02  4.2884E-02  4.2802E-01  1.1045E+04  1.7843E+00  8.2287E-02  1.0000E-02  1.0522E+00  4.1075E+00
             8.0305E+00
 PARAMETER: -4.4396E-01 -3.6214E+00 -3.0493E+00 -7.4860E-01  9.4097E+00  6.7904E-01 -2.3975E+00 -2.0634E+01  1.5088E-01  1.5128E+00
             2.1832E+00
 GRADIENT:  -3.1044E+00  1.4207E+00 -5.7750E-01  7.2644E-01 -2.4851E-06  1.4702E+00  7.8482E-04  0.0000E+00 -9.1107E-02  1.3001E-08
            -2.4172E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -866.896052593976        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      773
 NPARAMETR:  5.7612E-01  1.8112E-02  4.1864E-02  4.2021E-01  1.1582E+04  1.7761E+00  7.8205E-02  1.0000E-02  1.0505E+00  3.8088E+00
             8.0321E+00
 PARAMETER: -4.5144E-01 -3.9112E+00 -3.0733E+00 -7.6701E-01  9.4572E+00  6.7442E-01 -2.4484E+00 -2.0900E+01  1.4930E-01  1.4373E+00
             2.1834E+00
 GRADIENT:  -8.2652E-02  3.5177E-02 -1.5710E-03  1.8188E-02  1.0094E-04  2.8373E-03  2.0422E-05  0.0000E+00  8.3709E-02 -2.2938E-09
            -2.9183E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -866.939345100833        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      935
 NPARAMETR:  5.7508E-01  1.6473E-02  4.1706E-02  4.1946E-01  2.6970E+01  1.7759E+00  6.3474E-02  1.0000E-02  1.0503E+00  3.8737E+00
             8.0315E+00
 PARAMETER: -4.5324E-01 -4.0060E+00 -3.0771E+00 -7.6878E-01  3.3947E+00  6.7428E-01 -2.6571E+00 -2.0926E+01  1.4908E-01  1.4542E+00
             2.1834E+00
 GRADIENT:   3.2464E+00  2.6140E-03  6.0521E+00  2.8268E+00  4.1726E-02  1.4146E+00  9.8607E-06  0.0000E+00  3.3630E-01 -3.4516E-04
             1.7301E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -866.977069038140        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1008
 NPARAMETR:  5.7459E-01  1.5950E-02  4.1601E-02  4.1923E-01  1.2187E+01  1.7741E+00  6.3465E-02  1.0000E-02  1.0496E+00  3.8679E+00
             8.0258E+00
 PARAMETER: -4.5410E-01 -4.0383E+00 -3.0796E+00 -7.6933E-01  2.6004E+00  6.7328E-01 -2.6573E+00 -2.0926E+01  1.4843E-01  1.4527E+00
             2.1827E+00
 GRADIENT:   1.9151E+00  1.0042E-01  6.1808E+00  3.1605E+00 -1.2098E-02  9.2855E-01  1.2397E-05  0.0000E+00  3.7963E-01  3.4441E-03
             1.5308E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -866.986385183448        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1186
 NPARAMETR:  5.7658E-01  1.4485E-02  4.1616E-02  4.1926E-01  1.1703E+01  1.7786E+00  6.3386E-02  1.0000E-02  1.0475E+00  3.8219E+00
             8.0347E+00
 PARAMETER: -4.5064E-01 -4.1347E+00 -3.0793E+00 -7.6927E-01  2.5599E+00  6.7583E-01 -2.6585E+00 -2.0926E+01  1.4641E-01  1.4407E+00
             2.1838E+00
 GRADIENT:   1.0991E-02  5.5389E-03  6.9273E-02 -3.3523E-01 -2.0540E-03  4.9963E-02  6.6795E-06  0.0000E+00  1.7825E-02  1.1121E-03
             1.0928E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -866.987438458934        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1362
 NPARAMETR:  5.7965E-01  1.4564E-02  4.2232E-02  4.2386E-01  1.1625E+01  1.7801E+00  6.2546E-02  1.0000E-02  1.0476E+00  3.4646E+00
             8.0351E+00
 PARAMETER: -4.4533E-01 -4.1292E+00 -3.0646E+00 -7.5835E-01  2.5532E+00  6.7665E-01 -2.6719E+00 -2.0926E+01  1.4650E-01  1.3426E+00
             2.1838E+00
 GRADIENT:   2.3616E-02 -1.4828E-03 -3.7029E-02  5.4248E-02  6.5504E-04 -6.3880E-03  6.3441E-06  0.0000E+00 -3.4384E-03  7.0534E-04
            -2.5779E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -866.987475820072        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1455
 NPARAMETR:  5.7953E-01  1.4857E-02  4.2223E-02  4.2368E-01  1.1616E+01  1.7802E+00  6.3582E-02  1.0000E-02  1.0477E+00  2.8494E+00
             8.0362E+00
 PARAMETER: -4.4554E-01 -4.1093E+00 -3.0648E+00 -7.5877E-01  2.5524E+00  6.7673E-01 -2.6554E+00 -2.0926E+01  1.4662E-01  1.1471E+00
             2.1840E+00
 GRADIENT:   4.1023E+00  1.2002E-02  5.9624E+00  2.4607E+00 -5.9547E-04  1.6728E+00  7.9568E-06  0.0000E+00  1.0864E-01  8.8864E-04
             1.8865E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -866.987491405922        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1530
 NPARAMETR:  5.7942E-01  1.5028E-02  4.2200E-02  4.2347E-01  1.1599E+01  1.7802E+00  6.5331E-02  1.0000E-02  1.0478E+00  2.1039E+00
             8.0366E+00
 PARAMETER: -4.4572E-01 -4.0978E+00 -3.0653E+00 -7.5927E-01  2.5509E+00  6.7672E-01 -2.6283E+00 -2.0926E+01  1.4669E-01  8.4381E-01
             2.1840E+00
 GRADIENT:   4.0444E+00  1.6218E-02  6.0594E+00  2.3401E+00 -1.9509E-03  1.6910E+00  8.8063E-06  0.0000E+00  1.2079E-01  5.0883E-04
             1.9468E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -866.987495866411        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:     1608
 NPARAMETR:  5.7929E-01  1.5110E-02  4.2172E-02  4.2325E-01  1.1582E+01  1.7801E+00  6.7373E-02  1.0000E-02  1.0478E+00  1.5192E+00
             8.0367E+00
 PARAMETER: -4.4595E-01 -4.0924E+00 -3.0660E+00 -7.5978E-01  2.5495E+00  6.7669E-01 -2.5975E+00 -2.0926E+01  1.4671E-01  5.1819E-01
             2.1840E+00
 GRADIENT:   4.0176E+00  1.8279E-02  6.1023E+00  2.2848E+00 -2.8362E-03  1.6998E+00  9.5764E-06  0.0000E+00  1.2614E-01  2.6366E-04
             1.9740E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -866.987781938080        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     1755
 NPARAMETR:  5.7938E-01  1.4839E-02  4.2154E-02  4.2325E-01  1.1556E+01  1.7798E+00  6.8233E-02  1.0000E-02  1.0477E+00  1.2571E+00
             8.0351E+00
 PARAMETER: -4.4579E-01 -4.1105E+00 -3.0664E+00 -7.5978E-01  2.5472E+00  6.7651E-01 -2.5848E+00 -2.0926E+01  1.4658E-01  3.2880E-01
             2.1838E+00
 GRADIENT:   3.9811E-02  1.7304E-04  1.1745E-02 -5.6930E-02  1.5930E-03 -8.9733E-03  8.0622E-06  0.0000E+00  4.1826E-03  9.0850E-05
             3.1518E-03

0ITERATION NO.:   85    OBJECTIVE VALUE:  -866.987804803781        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1934
 NPARAMETR:  5.7921E-01  1.4814E-02  4.2128E-02  4.2308E-01  1.1527E+01  1.7798E+00  7.0333E-02  1.0000E-02  1.0477E+00  9.0292E-01
             8.0349E+00
 PARAMETER: -4.4610E-01 -4.1122E+00 -3.0670E+00 -7.6020E-01  2.5446E+00  6.7649E-01 -2.5545E+00 -2.0926E+01  1.4657E-01 -2.1260E-03
             2.1838E+00
 GRADIENT:  -1.0453E-02 -3.5028E-05 -1.3167E-02 -7.2275E-03  2.4638E-05 -4.5741E-03  8.5237E-06  0.0000E+00 -4.9096E-05  4.6053E-05
            -4.7620E-03

0ITERATION NO.:   90    OBJECTIVE VALUE:  -866.987810038916        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2111
 NPARAMETR:  5.7923E-01  1.4789E-02  4.2127E-02  4.2310E-01  1.1517E+01  1.7798E+00  6.8345E-02  1.0000E-02  1.0476E+00  7.0227E-01
             8.0349E+00
 PARAMETER: -4.4605E-01 -4.1139E+00 -3.0671E+00 -7.6015E-01  2.5439E+00  6.7653E-01 -2.5832E+00 -2.0926E+01  1.4654E-01 -2.5344E-01
             2.1838E+00
 GRADIENT:   2.8441E-02 -1.3177E-03 -7.7614E-02  6.3913E-02  5.0757E-04  3.3690E-03  8.0014E-06  0.0000E+00 -9.5963E-03  2.5610E-05
            -1.6530E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -866.987827143148        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2287
 NPARAMETR:  5.7953E-01  1.4843E-02  4.2186E-02  4.2352E-01  1.1513E+01  1.7800E+00  5.4472E-02  1.0000E-02  1.0477E+00  2.9681E-01
             8.0351E+00
 PARAMETER: -4.4554E-01 -4.1102E+00 -3.0657E+00 -7.5916E-01  2.5435E+00  6.7661E-01 -2.8101E+00 -2.0926E+01  1.4657E-01 -1.1147E+00
             2.1838E+00
 GRADIENT:   1.7626E-02 -6.4211E-04 -3.8702E-02  3.5996E-02  1.9110E-04  2.8853E-03  5.1348E-06  0.0000E+00 -5.1304E-03  4.8555E-06
            -7.3830E-03

0ITERATION NO.:  100    OBJECTIVE VALUE:  -866.987827396453        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     2454
 NPARAMETR:  5.7955E-01  1.4862E-02  4.2193E-02  4.2355E-01  1.1513E+01  1.7800E+00  5.4068E-02  1.0000E-02  1.0477E+00  2.7350E-01
             8.0352E+00
 PARAMETER: -4.4551E-01 -4.1090E+00 -3.0655E+00 -7.5909E-01  2.5435E+00  6.7661E-01 -2.8175E+00 -2.0926E+01  1.4660E-01 -1.1965E+00
             2.1838E+00
 GRADIENT:  -5.9281E-03  5.0499E-05  1.8903E-02 -3.1118E-02 -4.2818E-04  4.2751E-03  5.1032E-06  0.0000E+00  2.9613E-03  4.3462E-06
             6.1150E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2454
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.0913E-04  2.0151E-06  6.6901E-05 -2.1027E-02 -3.4487E-06
 SE:             2.9089E-02  5.6307E-06  2.0164E-04  2.5270E-02  3.0875E-05
 N:                     100         100         100         100         100

 P VAL.:         9.8329E-01  7.2044E-01  7.4005E-01  4.0537E-01  9.1106E-01

 ETASHRINKSD(%)  2.5493E+00  9.9981E+01  9.9324E+01  1.5342E+01  9.9897E+01
 ETASHRINKVR(%)  5.0336E+00  1.0000E+02  9.9995E+01  2.8331E+01  1.0000E+02
 EBVSHRINKSD(%)  2.2328E+00  9.9973E+01  9.9386E+01  1.5489E+01  9.9899E+01
 EBVSHRINKVR(%)  4.4157E+00  1.0000E+02  9.9996E+01  2.8579E+01  1.0000E+02
 RELATIVEINF(%)  6.0299E+00  8.4847E-07  4.8804E-05  8.2301E-01  1.0140E-05
 EPSSHRINKSD(%)  1.6286E+01
 EPSSHRINKVR(%)  2.9919E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -866.98782739645344     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -131.83700083271526     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -866.988       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.80E-01  1.49E-02  4.22E-02  4.24E-01  1.15E+01  1.78E+00  5.41E-02  1.00E-02  1.05E+00  2.73E-01  8.04E+00
 


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
+        9.90E+02
 
 TH 2
+       -2.99E+02  1.37E+03
 
 TH 3
+       -4.19E+03 -1.12E+02  2.19E+05
 
 TH 4
+       -8.12E+01  1.40E+02 -2.70E+04  3.79E+03
 
 TH 5
+        2.36E-01 -9.06E-01 -3.13E+00  3.35E-01  4.17E-03
 
 TH 6
+       -5.15E+00  1.82E+01  3.34E+01 -2.05E+01  4.02E-03  5.11E+01
 
 TH 7
+        1.38E+00 -1.81E+00 -4.31E+00 -2.28E-01 -8.73E-03  5.90E-01  1.59E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.74E-01 -1.51E+01  4.30E+02 -6.27E+01 -3.01E-02 -2.71E+00 -1.27E+00  0.00E+00  9.22E+01
 
 TH10
+       -1.58E-01  3.10E+00 -2.70E-01 -7.49E-01 -3.26E-03  7.09E-03  2.29E+00  0.00E+00 -1.27E+00  9.45E-02
 
 TH11
+       -1.30E+01  3.38E+00  1.16E+02 -9.67E+00 -4.43E-03  6.91E-01 -1.41E-02  0.00E+00  3.28E+00 -2.45E-03  5.78E+00
 
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
 #CPUT: Total CPU Time in Seconds,       44.341
Stop Time:
Sat Sep 18 15:16:23 CDT 2021
