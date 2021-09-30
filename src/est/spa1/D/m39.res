Thu Sep 30 03:01:30 CDT 2021
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
$DATA ../../../../data/spa1/D/dat39.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m39.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   16929.7528320537        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.2651E+02  3.5109E+02 -1.5177E+01  1.3548E+02 -4.2355E+01 -1.4381E+03 -6.9327E+02 -4.5732E+01 -1.2588E+03 -2.9686E+02
            -3.3953E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -639.392908896283        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2060E+00  9.8420E-01  9.6046E-01  1.5920E+00  1.2797E+00  2.3028E+00  1.2331E+00  9.5310E-01  1.4391E+00  1.0238E+00
             1.4059E+01
 PARAMETER:  2.8731E-01  8.4078E-02  5.9662E-02  5.6498E-01  3.4660E-01  9.3414E-01  3.0956E-01  5.1963E-02  4.6400E-01  1.2348E-01
             2.7432E+00
 GRADIENT:  -3.0028E+01  8.4125E+00 -9.0483E+00  3.8036E-01 -2.7320E+00  7.7839E+01  9.4463E-01  5.7796E+00  1.6202E+00  2.6333E+00
             2.5746E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -671.125479385471        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2416E+00  7.9365E-01  1.3855E+00  1.9866E+00  7.7156E+00  1.9848E+00  4.1629E+00  3.8326E-01  1.9736E+00  4.4701E+00
             1.2531E+01
 PARAMETER:  3.1639E-01 -1.3111E-01  4.2605E-01  7.8643E-01  2.1432E+00  7.8551E-01  1.5262E+00 -8.5903E-01  7.7983E-01  1.5974E+00
             2.6282E+00
 GRADIENT:  -3.7853E+00  2.3085E+01  7.7174E+00  6.1946E+01 -8.4968E-01  6.6072E+00  1.1845E+01 -1.9561E-01  2.6005E+01 -5.9378E-02
             2.0605E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -694.319960554882        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.1708E+00  8.4538E-01  1.0408E+00  1.5435E+00  9.2060E+00  1.8254E+00  1.5014E+00  7.2293E-02  1.7870E+00  8.4162E+00
             1.2121E+01
 PARAMETER:  2.5771E-01 -6.7969E-02  1.4000E-01  5.3405E-01  2.3199E+00  7.0180E-01  5.0637E-01 -2.5270E+00  6.8053E-01  2.2302E+00
             2.5950E+00
 GRADIENT:  -1.8256E+01  1.6140E+01  1.1075E+01  1.6037E+01 -3.6350E+00 -1.6772E+00  5.2289E+00  4.7713E-05  1.2354E+01  1.7448E+00
             2.0137E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -756.329162202013        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  8.8132E-01  1.7215E-01  2.0248E-01  1.1482E+00  1.7598E+01  1.6088E+00  6.2890E-01  1.0000E-02  1.0742E+00  1.0181E+01
             8.3226E+00
 PARAMETER: -2.6336E-02 -1.6594E+00 -1.4971E+00  2.3821E-01  2.9678E+00  5.7550E-01 -3.6378E-01 -7.8017E+00  1.7158E-01  2.4205E+00
             2.2190E+00
 GRADIENT:   3.5322E+01  3.8514E+01 -4.1638E+01  1.5290E+02 -2.0605E+01 -2.9112E+01  4.1115E+00  0.0000E+00 -3.4607E+01  3.4351E+01
            -1.2506E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -784.451994907843        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  6.7775E-01  6.4900E-02  9.1894E-02  7.2340E-01  7.3791E+01  1.5055E+00  3.9217E-02  1.0000E-02  9.2501E-01  1.3027E+01
             7.1993E+00
 PARAMETER: -2.8898E-01 -2.6349E+00 -2.2871E+00 -2.2379E-01  4.4012E+00  5.0912E-01 -3.1386E+00 -1.2452E+01  2.2044E-02  2.6670E+00
             2.0740E+00
 GRADIENT:   4.9815E+01  1.6832E+01 -4.6323E+00  9.9562E+01 -2.8195E-01 -2.5194E+01  1.0693E-02  0.0000E+00 -1.1768E+01  2.2432E-01
            -2.5669E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -784.529099496193        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  6.6354E-01  6.0225E-02  8.5843E-02  6.9712E-01  8.3207E+01  1.4961E+00  3.1275E-02  1.0000E-02  9.0518E-01  1.3408E+01
             7.1622E+00
 PARAMETER: -3.1016E-01 -2.7097E+00 -2.3552E+00 -2.6079E-01  4.5213E+00  5.0286E-01 -3.3649E+00 -1.2876E+01  3.7769E-04  2.6959E+00
             2.0688E+00
 GRADIENT:   5.3557E+01  1.4783E+01 -8.7286E+00  1.1208E+02 -1.9683E-01 -2.5437E+01  6.3882E-03  0.0000E+00 -1.5038E+01  1.7021E-01
            -2.6417E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -817.836923870206        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      632
 NPARAMETR:  6.5082E-01  2.9719E-02  6.7898E-02  5.9463E-01  8.4324E+02  1.5284E+00  1.0000E-02  1.0000E-02  6.1763E-01  2.2094E+01
             9.1479E+00
 PARAMETER: -3.2953E-01 -3.4160E+00 -2.5897E+00 -4.1982E-01  6.8373E+00  5.2422E-01 -6.6137E+00 -1.7016E+01 -3.8186E-01  3.1953E+00
             2.3135E+00
 GRADIENT:   3.5655E+01 -1.3989E-01 -3.3524E+01  6.6803E+01  7.1482E-03 -2.3237E+00  0.0000E+00  0.0000E+00 -2.3533E+01 -6.4112E-05
            -3.3129E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -827.857309778883        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      809
 NPARAMETR:  5.4693E-01  1.0000E-02  4.5986E-02  4.3872E-01  5.5828E+03  1.5510E+00  1.0000E-02  1.0000E-02  7.6023E-01  2.6866E+01
             9.0124E+00
 PARAMETER: -5.0343E-01 -4.9450E+00 -2.9794E+00 -7.2390E-01  8.7274E+00  5.3893E-01 -1.0731E+01 -2.0716E+01 -1.7413E-01  3.3909E+00
             2.2986E+00
 GRADIENT:   3.1847E+00  0.0000E+00  3.2083E+01 -3.8375E+01  1.2628E-04  1.1888E+01  0.0000E+00  0.0000E+00 -5.6269E+00  6.8837E-07
             7.6134E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -830.775723751731        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      998
 NPARAMETR:  4.7580E-01  1.0000E-02  3.1295E-02  3.3771E-01  1.7212E+04  1.4562E+00  1.0000E-02  1.0000E-02  7.9115E-01  3.3256E+01
             8.9151E+00
 PARAMETER: -6.4275E-01 -5.9211E+00 -3.3643E+00 -9.8558E-01  9.8534E+00  4.7585E-01 -1.3514E+01 -2.3667E+01 -1.3427E-01  3.6042E+00
             2.2877E+00
 GRADIENT:   1.3861E+00  0.0000E+00  6.6524E-02 -1.3707E+00  3.8868E-05  3.0966E-01  0.0000E+00  0.0000E+00 -9.6733E-02  1.8142E-07
             9.1021E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -830.783828499294        NO. OF FUNC. EVALS.: 198
 CUMULATIVE NO. OF FUNC. EVALS.:     1196             RESET HESSIAN, TYPE I
 NPARAMETR:  4.7527E-01  1.0000E-02  3.1157E-02  3.3744E-01  7.9870E+03  1.4555E+00  1.0000E-02  1.0000E-02  7.9239E-01  2.8734E+01
             8.9009E+00
 PARAMETER: -6.4386E-01 -5.9211E+00 -3.3687E+00 -9.8637E-01  9.0856E+00  4.7534E-01 -1.3514E+01 -2.3667E+01 -1.3270E-01  3.4581E+00
             2.2862E+00
 GRADIENT:   6.3568E+01  0.0000E+00  7.5937E+01  3.2792E+01  9.7922E-05  8.0966E+00  0.0000E+00  0.0000E+00  1.1470E-01  1.6430E-06
             2.1859E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -830.897707313204        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:     1347
 NPARAMETR:  4.6752E-01  1.0000E-02  3.0049E-02  3.2806E-01  2.1343E+01  1.4484E+00  1.0000E-02  1.0000E-02  7.9193E-01  2.5808E-01
             8.9146E+00
 PARAMETER: -6.6031E-01 -5.9211E+00 -3.4049E+00 -1.0146E+00  3.1607E+00  4.7046E-01 -1.3514E+01 -2.3667E+01 -1.3328E-01 -1.2545E+00
             2.2877E+00
 GRADIENT:   5.9608E+01  0.0000E+00  8.2624E+01  2.9555E+01  1.8261E-02  7.2641E+00  0.0000E+00  0.0000E+00  4.4692E-01  1.1396E-05
             2.5032E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -830.916584583539        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1509
 NPARAMETR:  4.7055E-01  1.0000E-02  3.0164E-02  3.2899E-01  1.1743E+01  1.4518E+00  1.0000E-02  1.0000E-02  7.9070E-01  1.1945E+00
             8.9069E+00
 PARAMETER: -6.5385E-01 -5.9211E+00 -3.4011E+00 -1.0117E+00  2.5633E+00  4.7283E-01 -1.3514E+01 -2.3667E+01 -1.3484E-01  2.7770E-01
             2.2868E+00
 GRADIENT:  -2.5720E-01  0.0000E+00  8.5088E-02 -1.6546E+00 -1.6234E-02 -5.4806E-02  0.0000E+00  0.0000E+00  1.4796E-02  7.6139E-04
             8.2554E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -830.922078978076        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     1706             RESET HESSIAN, TYPE I
 NPARAMETR:  4.7056E-01  1.0000E-02  3.0100E-02  3.2848E-01  1.1720E+01  1.4525E+00  1.0000E-02  1.0000E-02  7.9189E-01  1.0278E+00
             8.9036E+00
 PARAMETER: -6.5383E-01 -5.9211E+00 -3.4032E+00 -1.0133E+00  2.5613E+00  4.7330E-01 -1.3514E+01 -2.3667E+01 -1.3333E-01  1.2742E-01
             2.2865E+00
 GRADIENT:   6.3235E+01  0.0000E+00  8.2616E+01  2.8112E+01  6.4734E-04  8.1019E+00  0.0000E+00  0.0000E+00  4.1748E-01  5.9783E-04
             2.3292E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -830.924728729746        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1890
 NPARAMETR:  4.7043E-01  1.0000E-02  3.0067E-02  3.2855E-01  1.1922E+01  1.4522E+00  1.0000E-02  1.0000E-02  7.9120E-01  1.0123E+00
             8.8999E+00
 PARAMETER: -6.5411E-01 -5.9211E+00 -3.4043E+00 -1.0131E+00  2.5784E+00  4.7306E-01 -1.3514E+01 -2.3667E+01 -1.3421E-01  1.1218E-01
             2.2860E+00
 GRADIENT:   8.4093E-01  0.0000E+00 -2.1681E+00  7.1998E-01  1.7047E-03  4.8571E-02  0.0000E+00  0.0000E+00 -2.5643E-02  4.9960E-04
            -2.4069E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -830.957040743971        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2069             RESET HESSIAN, TYPE I
 NPARAMETR:  4.6755E-01  1.0000E-02  2.9600E-02  3.2428E-01  1.2027E+01  1.4507E+00  1.0000E-02  1.0000E-02  7.9190E-01  9.4134E-01
             8.8987E+00
 PARAMETER: -6.6024E-01 -5.9211E+00 -3.4200E+00 -1.0261E+00  2.5872E+00  4.7204E-01 -1.3514E+01 -2.3667E+01 -1.3332E-01  3.9548E-02
             2.2859E+00
 GRADIENT:   6.4082E+01  0.0000E+00  8.3687E+01  2.8138E+01  8.2414E-03  8.1156E+00  0.0000E+00  0.0000E+00  3.9520E-01  4.4414E-04
             2.2952E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -830.977119934512        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:     2200             RESET HESSIAN, TYPE I
 NPARAMETR:  4.6590E-01  1.0000E-02  2.9194E-02  3.2181E-01  1.2098E+01  1.4493E+00  1.0000E-02  1.0000E-02  7.9201E-01  9.1844E-01
             8.8874E+00
 PARAMETER: -6.6379E-01 -5.9211E+00 -3.4338E+00 -1.0338E+00  2.5930E+00  4.7111E-01 -1.3514E+01 -2.3667E+01 -1.3319E-01  1.4925E-02
             2.2846E+00
 GRADIENT:   6.6351E+01  0.0000E+00  7.9654E+01  3.4067E+01  2.4564E-02  8.0410E+00  0.0000E+00  0.0000E+00  1.8004E-01  4.1783E-04
             2.1078E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -830.982436137172        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2383
 NPARAMETR:  4.6577E-01  1.0000E-02  2.9180E-02  3.2133E-01  1.2309E+01  1.4490E+00  1.0000E-02  1.0000E-02  7.9109E-01  8.8967E-01
             8.8993E+00
 PARAMETER: -6.6406E-01 -5.9211E+00 -3.4343E+00 -1.0353E+00  2.6104E+00  4.7087E-01 -1.3514E+01 -2.3667E+01 -1.3434E-01 -1.6910E-02
             2.2860E+00
 GRADIENT:   2.0778E+00  0.0000E+00 -3.7082E+00  1.4610E+00  1.5256E-02  1.3975E-01  0.0000E+00  0.0000E+00 -6.1279E-02  3.4869E-04
            -4.0093E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -830.988663296811        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2580             RESET HESSIAN, TYPE I
 NPARAMETR:  4.6485E-01  1.0000E-02  2.9112E-02  3.2051E-01  1.2164E+01  1.4486E+00  1.0000E-02  1.0000E-02  7.9153E-01  6.6161E-01
             8.8972E+00
 PARAMETER: -6.6605E-01 -5.9211E+00 -3.4366E+00 -1.0378E+00  2.5985E+00  4.7058E-01 -1.3514E+01 -2.3667E+01 -1.3379E-01 -3.1308E-01
             2.2857E+00
 GRADIENT:   6.5203E+01  0.0000E+00  8.3026E+01  3.0222E+01  1.3621E-02  8.0284E+00  0.0000E+00  0.0000E+00  2.8704E-01  2.2808E-04
             2.2649E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -830.988743321066        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     2676
 NPARAMETR:  4.6500E-01  1.0000E-02  2.9107E-02  3.2038E-01  1.2303E+01  1.4488E+00  1.0000E-02  1.0000E-02  7.9156E-01  6.4414E-01
             8.8985E+00
 PARAMETER: -6.6573E-01 -5.9211E+00 -3.4368E+00 -1.0383E+00  2.6098E+00  4.7071E-01 -1.3514E+01 -2.3667E+01 -1.3375E-01 -3.3983E-01
             2.2859E+00
 GRADIENT:   1.4165E+00  0.0000E+00 -1.8681E+00 -7.0534E-01  8.3662E-03  1.7745E-01  0.0000E+00  0.0000E+00  6.2927E-02  1.8471E-04
            -1.0267E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -830.992447795259        NO. OF FUNC. EVALS.: 205
 CUMULATIVE NO. OF FUNC. EVALS.:     2881             RESET HESSIAN, TYPE I
 NPARAMETR:  4.6452E-01  1.0000E-02  2.9035E-02  3.2009E-01  1.2103E+01  1.4481E+00  1.0000E-02  1.0000E-02  7.9115E-01  6.1545E-01
             8.8970E+00
 PARAMETER: -6.6675E-01 -5.9211E+00 -3.4392E+00 -1.0392E+00  2.5934E+00  4.7028E-01 -1.3514E+01 -2.3667E+01 -1.3426E-01 -3.8541E-01
             2.2857E+00
 GRADIENT:   6.5517E+01  0.0000E+00  8.2079E+01  3.1585E+01  1.4522E-02  7.9747E+00  0.0000E+00  0.0000E+00  1.9582E-01  2.0724E-04
             2.2458E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -830.994516596604        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     3076
 NPARAMETR:  4.6431E-01  1.0000E-02  2.8997E-02  3.1979E-01  1.2124E+01  1.4480E+00  1.0000E-02  1.0000E-02  7.9111E-01  6.0594E-01
             8.8966E+00
 PARAMETER: -6.6719E-01 -5.9211E+00 -3.4406E+00 -1.0401E+00  2.5952E+00  4.7018E-01 -1.3514E+01 -2.3667E+01 -1.3431E-01 -4.0098E-01
             2.2857E+00
 GRADIENT:   1.2980E+00  0.0000E+00 -3.4780E+00  1.4003E+00  4.6797E-03  5.8989E-02  0.0000E+00  0.0000E+00 -6.5650E-02  1.7618E-04
            -4.6677E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -830.998788303179        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     3251             RESET HESSIAN, TYPE I
 NPARAMETR:  4.6378E-01  1.0000E-02  2.8943E-02  3.1909E-01  1.2195E+01  1.4480E+00  1.0000E-02  1.0000E-02  7.9166E-01  4.1823E-01
             8.8955E+00
 PARAMETER: -6.6835E-01 -5.9211E+00 -3.4424E+00 -1.0423E+00  2.6010E+00  4.7015E-01 -1.3514E+01 -2.3667E+01 -1.3363E-01 -7.7172E-01
             2.2855E+00
 GRADIENT:   6.5410E+01  0.0000E+00  8.3399E+01  3.0285E+01  1.2757E-02  8.0307E+00  0.0000E+00  0.0000E+00  2.9834E-01  1.0419E-04
             2.2554E+01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -831.010669620515        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:     3407             RESET HESSIAN, TYPE I
 NPARAMETR:  4.6151E-01  1.0000E-02  2.8661E-02  3.1620E-01  1.1665E+01  1.4471E+00  1.0000E-02  1.0000E-02  7.9307E-01  1.7468E-01
             8.8944E+00
 PARAMETER: -6.7326E-01 -5.9211E+00 -3.4522E+00 -1.0514E+00  2.5566E+00  4.6958E-01 -1.3514E+01 -2.3667E+01 -1.3184E-01 -1.6448E+00
             2.2854E+00
 GRADIENT:   6.4657E+01  0.0000E+00  8.6979E+01  2.7060E+01 -2.3337E-02  8.1176E+00  0.0000E+00  0.0000E+00  6.1125E-01  2.8679E-05
             2.3210E+01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -831.016635083665        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3589
 NPARAMETR:  4.6175E-01  1.0000E-02  2.8621E-02  3.1652E-01  1.2259E+01  1.4464E+00  1.0000E-02  1.0000E-02  7.9130E-01  1.6959E-01
             8.8944E+00
 PARAMETER: -6.7272E-01 -5.9211E+00 -3.4536E+00 -1.0504E+00  2.6063E+00  4.6909E-01 -1.3514E+01 -2.3667E+01 -1.3407E-01 -1.6744E+00
             2.2854E+00
 GRADIENT:   7.3785E-01  0.0000E+00 -2.8925E+00  6.1315E-01  2.7465E-03  5.2884E-02  0.0000E+00  0.0000E+00 -2.7209E-02  1.3547E-05
            -3.9383E-01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -831.020053997512        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     3785
 NPARAMETR:  4.6132E-01  1.0000E-02  2.8562E-02  3.1583E-01  1.2136E+01  1.4459E+00  1.0000E-02  1.0000E-02  7.9126E-01  1.3068E-02
             8.8978E+00
 PARAMETER: -6.7366E-01 -5.9211E+00 -3.4557E+00 -1.0526E+00  2.5961E+00  4.6875E-01 -1.3514E+01 -2.3667E+01 -1.3412E-01 -4.2376E+00
             2.2858E+00
 GRADIENT:   4.8027E-01  0.0000E+00 -1.8349E+00 -7.0196E-01 -5.5655E-03  2.4747E-02  0.0000E+00  0.0000E+00  2.7751E-02  8.9498E-08
             1.1888E-01

0ITERATION NO.:  126    OBJECTIVE VALUE:  -831.020053997512        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     3807
 NPARAMETR:  4.6132E-01  1.0000E-02  2.8562E-02  3.1583E-01  1.2136E+01  1.4459E+00  1.0000E-02  1.0000E-02  7.9126E-01  1.3068E-02
             8.8978E+00
 PARAMETER: -6.7366E-01 -5.9211E+00 -3.4557E+00 -1.0526E+00  2.5961E+00  4.6875E-01 -1.3514E+01 -2.3667E+01 -1.3412E-01 -4.2376E+00
             2.2858E+00
 GRADIENT:   4.8027E-01  0.0000E+00 -1.8349E+00 -7.0196E-01 -5.5655E-03  2.4747E-02  0.0000E+00  0.0000E+00  2.7751E-02  8.9498E-08
             1.1888E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3807
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4404E-03  9.7924E-07  9.9209E-05 -1.9439E-02 -2.9398E-07
 SE:             2.8946E-02  1.2287E-06  2.5613E-04  2.3654E-02  1.2379E-06
 N:                     100         100         100         100         100

 P VAL.:         9.3281E-01  4.2545E-01  6.9851E-01  4.1119E-01  8.1228E-01

 ETASHRINKSD(%)  3.0288E+00  9.9996E+01  9.9142E+01  2.0756E+01  9.9996E+01
 ETASHRINKVR(%)  5.9659E+00  1.0000E+02  9.9993E+01  3.7203E+01  1.0000E+02
 EBVSHRINKSD(%)  2.8984E+00  9.9995E+01  9.9209E+01  2.1734E+01  9.9995E+01
 EBVSHRINKVR(%)  5.7128E+00  1.0000E+02  9.9994E+01  3.8744E+01  1.0000E+02
 RELATIVEINF(%)  3.2878E+00  1.6300E-08  5.0998E-05  4.5096E-01  1.1821E-08
 EPSSHRINKSD(%)  1.1284E+01
 EPSSHRINKVR(%)  2.1295E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -831.02005399751192     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       87.918479207160772     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    78.52
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     9.76
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -831.020       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.61E-01  1.00E-02  2.86E-02  3.16E-01  1.21E+01  1.45E+00  1.00E-02  1.00E-02  7.91E-01  1.31E-02  8.90E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.66E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.03E+04  0.00E+00  6.41E+05
 
 TH 4
+        1.10E+03  0.00E+00 -6.84E+04  7.30E+03
 
 TH 5
+        7.24E-02  0.00E+00 -4.50E+00  4.80E-01  3.16E-05
 
 TH 6
+        9.32E-01  0.00E+00 -5.80E+01  6.19E+00  4.07E-04  5.24E-03
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.04E+01  0.00E+00  6.49E+02 -6.93E+01 -4.56E-03 -5.87E-02  0.00E+00  0.00E+00  6.58E-01
 
 TH10
+        6.80E-04  0.00E+00 -4.23E-02  4.51E-03  2.97E-07  3.83E-06  0.00E+00  0.00E+00 -4.28E-05  2.79E-09
 
 TH11
+       -4.72E+00  0.00E+00  2.93E+02 -3.13E+01 -2.06E-03 -2.65E-02  0.00E+00  0.00E+00  2.97E-01 -1.94E-05  1.34E-01
 
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
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.33E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.29E+04  0.00E+00  7.97E+05
 
 TH 4
+       -1.24E+02  0.00E+00 -8.49E+04  1.02E+04
 
 TH 5
+        3.68E-01  0.00E+00 -5.61E+00  4.47E-01  1.68E-03
 
 TH 6
+        4.13E+00  0.00E+00 -7.58E+01 -3.36E+01  7.81E-03  8.10E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -7.81E+00  0.00E+00  8.07E+02 -9.02E+01 -1.34E-02 -3.50E+00  0.00E+00  0.00E+00  1.03E+02
 
 TH10
+       -3.82E-03  0.00E+00 -5.36E-02 -4.79E-03  1.37E-04 -4.01E-03  0.00E+00  0.00E+00 -1.06E-01 -4.47E-02
 
 TH11
+       -2.53E+01  0.00E+00  3.66E+02 -2.80E+01 -5.58E-03  1.64E+00  0.00E+00  0.00E+00  5.40E+00 -7.42E-04  6.24E+00
 
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
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        2.36E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.83E+04  0.00E+00  9.92E+05
 
 TH 4
+        2.88E+02  0.00E+00 -9.82E+04  1.10E+04
 
 TH 5
+        3.37E-01  0.00E+00 -6.44E+00  4.92E-01  7.23E-05
 
 TH 6
+        6.99E+01  0.00E+00 -3.31E+03  2.64E+02  2.13E-02  8.68E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.46E+02  0.00E+00  9.34E+02 -2.34E+02 -1.14E-03  1.37E+01  0.00E+00  0.00E+00  1.28E+02
 
 TH10
+        2.92E-04  0.00E+00 -7.48E-03  6.35E-04  7.39E-08 -1.74E-06  0.00E+00  0.00E+00 -5.96E-05  4.50E-09
 
 TH11
+       -8.13E+01  0.00E+00  5.22E+02 -2.22E+01 -1.34E-02  1.19E+01  0.00E+00  0.00E+00 -1.51E+01 -1.15E-04  1.11E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       88.328
Stop Time:
Thu Sep 30 03:03:00 CDT 2021
