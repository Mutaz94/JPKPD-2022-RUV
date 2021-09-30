Wed Sep 29 23:37:55 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat72.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m72.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1340.86721219362        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.4254E+02  2.8766E+01  1.8437E+02 -3.5903E+01  1.0131E+01  6.1672E+01 -5.4726E+00 -2.1142E+02  3.0490E+01 -3.3805E+01
            -1.3591E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1808.66471033069        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0170E+00  1.0533E+00  9.3360E-01  1.0460E+00  1.0473E+00  8.9261E-01  9.6307E-01  9.4644E-01  8.2645E-01  9.1547E-01
             2.0639E+00
 PARAMETER:  1.1682E-01  1.5191E-01  3.1294E-02  1.4495E-01  1.4624E-01 -1.3611E-02  6.2366E-02  4.4947E-02 -9.0622E-02  1.1686E-02
             8.2459E-01
 GRADIENT:   4.0446E+01  1.3408E+01  4.2387E+00  2.4261E+01  1.5141E+01 -4.2322E+00 -3.6449E+00 -4.7882E-01 -4.5083E-01  5.5515E-01
            -7.0567E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1820.33323362049        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      180
 NPARAMETR:  9.9388E-01  9.0286E-01  4.2655E-01  1.0841E+00  6.1375E-01  8.8540E-01  1.3307E+00  3.2176E-01  7.4736E-01  6.6704E-01
             2.2124E+00
 PARAMETER:  9.3860E-02 -2.1928E-03 -7.5202E-01  1.8079E-01 -3.8816E-01 -2.1718E-02  3.8571E-01 -1.0339E+00 -1.9121E-01 -3.0490E-01
             8.9408E-01
 GRADIENT:  -1.6292E+02  1.3103E+01 -3.1828E+01  5.5614E+01  5.9278E+01 -2.7739E+01  2.1898E+01  5.2360E-01 -2.6890E+00  1.5317E+01
             4.9450E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1824.83364614813        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      359
 NPARAMETR:  1.0527E+00  7.1470E-01  4.8017E-01  1.1874E+00  5.8288E-01  8.9718E-01  1.5802E+00  3.5746E-01  6.8094E-01  6.1231E-01
             2.2981E+00
 PARAMETER:  1.5139E-01 -2.3589E-01 -6.3361E-01  2.7178E-01 -4.3978E-01 -8.4988E-03  5.5754E-01 -9.2874E-01 -2.8428E-01 -3.9051E-01
             9.3210E-01
 GRADIENT:  -9.9860E-01  4.3668E+00 -3.7273E+01  4.9916E+01  6.6495E+01 -1.0656E+01  1.6090E+01  8.5436E-01 -6.4503E+00  1.0680E+01
             6.8850E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1832.25833036495        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      535
 NPARAMETR:  1.0201E+00  3.8894E-01  4.2008E-01  1.2850E+00  4.3659E-01  9.0816E-01  2.3505E+00  5.1407E-01  7.1001E-01  5.4118E-01
             2.1784E+00
 PARAMETER:  1.1987E-01 -8.4434E-01 -7.6730E-01  3.5074E-01 -7.2876E-01  3.6639E-03  9.5464E-01 -5.6541E-01 -2.4248E-01 -5.1400E-01
             8.7861E-01
 GRADIENT:  -6.6446E+01 -1.2904E+00 -1.1163E+01 -1.4060E+00  3.3841E+01 -6.9197E+00  1.2775E+01  3.4069E+00  5.9587E+00  1.0077E+01
             5.2322E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1836.09883698141        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      711
 NPARAMETR:  1.0310E+00  2.6507E-01  3.2963E-01  1.2657E+00  3.5150E-01  9.1622E-01  2.8418E+00  5.1768E-01  7.2706E-01  4.7077E-01
             2.0695E+00
 PARAMETER:  1.3054E-01 -1.2278E+00 -1.0098E+00  3.3560E-01 -9.4555E-01  1.2498E-02  1.1444E+00 -5.5840E-01 -2.1875E-01 -6.5340E-01
             8.2731E-01
 GRADIENT:  -3.0464E+01 -1.1392E+01 -1.9193E+01 -6.0840E+00  4.2190E+01 -3.1788E+00  9.3627E+00  1.1110E+00  1.7236E+00  5.4973E+00
             2.9584E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1844.38218463346        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      888
 NPARAMETR:  1.0512E+00  4.7167E-01  2.7037E-01  1.1591E+00  3.4156E-01  9.3753E-01  1.4611E+00  1.1243E-01  8.0889E-01  5.0618E-01
             1.9840E+00
 PARAMETER:  1.4989E-01 -6.5147E-01 -1.2080E+00  2.4765E-01 -9.7424E-01  3.5490E-02  4.7921E-01 -2.0854E+00 -1.1209E-01 -5.8087E-01
             7.8512E-01
 GRADIENT:  -7.7213E-01 -7.0202E+00 -6.5848E-01 -5.2130E+00  1.0544E+01 -2.2931E-01 -1.7018E+00 -3.8371E-01  1.1328E+00 -2.0743E-01
             1.7293E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1848.26287063609        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1024
 NPARAMETR:  1.0276E+00  4.7120E-01  2.6282E-01  1.1453E+00  3.3206E-01  9.3188E-01  1.3920E+00  9.3994E-01  8.0627E-01  4.6182E-01
             1.9376E+00
 PARAMETER:  1.2724E-01 -6.5247E-01 -1.2363E+00  2.3568E-01 -1.0024E+00  2.9452E-02  4.3073E-01  3.8063E-02 -1.1534E-01 -6.7259E-01
             7.6144E-01
 GRADIENT:   5.4324E+01  1.0249E+01  4.8016E+01  3.7639E+01  1.1791E+02  2.1264E+00  7.7349E+00  3.8770E+00 -3.2229E+00  1.4471E+01
             6.0042E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1851.53885842438        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:     1121
 NPARAMETR:  9.5204E-01  4.1930E-01  2.3047E-01  1.0982E+00  2.9464E-01  9.3096E-01  1.3154E+00  1.4497E+00  8.5195E-01  3.3866E-01
             1.6872E+00
 PARAMETER:  5.0850E-02 -7.6916E-01 -1.3676E+00  1.9371E-01 -1.1220E+00  2.8465E-02  3.7414E-01  4.7135E-01 -6.0233E-02 -9.8275E-01
             6.2306E-01
 GRADIENT:  -2.4675E+02 -1.0773E+01  2.8912E+01 -8.0689E+01  2.2327E+01 -2.9696E+01 -1.0879E+00  1.9489E+01 -2.1972E+01  5.8497E+00
             6.2638E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1872.01525866951        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1299
 NPARAMETR:  1.0394E+00  4.2809E-01  2.1567E-01  1.1083E+00  2.8652E-01  9.1694E-01  1.3443E+00  1.4218E+00  9.4172E-01  1.6963E-01
             1.4416E+00
 PARAMETER:  1.3860E-01 -7.4842E-01 -1.4340E+00  2.0286E-01 -1.1499E+00  1.3283E-02  3.9586E-01  4.5192E-01  3.9958E-02 -1.6742E+00
             4.6575E-01
 GRADIENT:  -2.0481E+01  1.7646E+01  2.3258E+01 -1.7314E+01  1.6256E+01 -1.3777E+01 -2.7805E+00  9.1308E+00 -3.7722E+00  8.2619E-01
            -3.1641E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1879.29922216211        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1475
 NPARAMETR:  1.0459E+00  3.4361E-01  1.7572E-01  1.0825E+00  2.4577E-01  9.4766E-01  1.4922E+00  1.3805E+00  1.0603E+00  3.7991E-02
             1.4400E+00
 PARAMETER:  1.4492E-01 -9.6825E-01 -1.6389E+00  1.7923E-01 -1.3034E+00  4.6242E-02  5.0028E-01  4.2246E-01  1.5857E-01 -3.1704E+00
             4.6462E-01
 GRADIENT:  -5.6781E-01  4.1060E-01 -1.6974E-01  4.2554E-01 -1.3855E-01 -2.4659E-02  1.4663E-01  1.4511E-01  7.8374E-02 -3.2092E-03
            -1.3537E-01

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1879.29922216211        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1497
 NPARAMETR:  1.0459E+00  3.4361E-01  1.7572E-01  1.0825E+00  2.4577E-01  9.4766E-01  1.4922E+00  1.3805E+00  1.0603E+00  3.7991E-02
             1.4400E+00
 PARAMETER:  1.4492E-01 -9.6825E-01 -1.6389E+00  1.7923E-01 -1.3034E+00  4.6242E-02  5.0028E-01  4.2246E-01  1.5857E-01 -3.1704E+00
             4.6462E-01
 GRADIENT:  -5.6781E-01  4.1060E-01 -1.6974E-01  4.2554E-01 -1.3855E-01 -2.4659E-02  1.4663E-01  1.4511E-01  7.8374E-02 -3.2092E-03
            -1.3537E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1497
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.9604E-04  2.9905E-02 -1.8608E-02 -9.4023E-03  2.0207E-03
 SE:             2.9785E-02  2.1336E-02  2.6335E-02  2.8446E-02  1.6220E-03
 N:                     100         100         100         100         100

 P VAL.:         9.7332E-01  1.6103E-01  4.7982E-01  7.4100E-01  2.1284E-01

 ETASHRINKSD(%)  2.1678E-01  2.8522E+01  1.1773E+01  4.7016E+00  9.4566E+01
 ETASHRINKVR(%)  4.3310E-01  4.8909E+01  2.2161E+01  9.1821E+00  9.9705E+01
 EBVSHRINKSD(%)  7.6123E-01  2.9626E+01  1.1054E+01  5.0379E+00  9.5145E+01
 EBVSHRINKVR(%)  1.5167E+00  5.0474E+01  2.0886E+01  9.8220E+00  9.9764E+01
 RELATIVEINF(%)  9.8433E+01  2.1061E+01  2.2179E+01  7.0689E+01  3.9235E-02
 EPSSHRINKSD(%)  3.9117E+01
 EPSSHRINKVR(%)  6.2932E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1879.2992221621098     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -960.36068895743711     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.49
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.02
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1879.299       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  3.44E-01  1.76E-01  1.08E+00  2.46E-01  9.48E-01  1.49E+00  1.38E+00  1.06E+00  3.80E-02  1.44E+00
 


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
+        1.11E+03
 
 TH 2
+       -2.17E+01  2.18E+03
 
 TH 3
+       -2.29E+01  1.72E+03  1.49E+04
 
 TH 4
+       -5.41E+00  1.88E+02 -2.42E+02  7.62E+02
 
 TH 5
+        4.18E+01 -5.52E+03 -1.63E+04 -1.08E+03  3.28E+04
 
 TH 6
+        2.19E+00 -8.90E+00  6.13E+00 -2.22E+00  2.32E+00  2.15E+02
 
 TH 7
+        1.62E+00  9.71E+01 -3.05E+01 -5.41E+00 -5.61E+01  1.98E-01  3.11E+01
 
 TH 8
+        6.43E-01  3.34E+01 -4.56E+01 -2.65E+00  1.82E+01  1.32E+00  4.18E+00  6.92E+01
 
 TH 9
+        3.86E+00 -2.34E+01  5.39E+01 -8.07E+00  5.48E+02  2.62E-01  7.45E+00 -6.66E+00  1.42E+02
 
 TH10
+       -1.45E-02 -7.45E-01 -1.38E+01 -3.59E+00  8.56E+01  1.13E-02  2.44E+00  3.74E+00  6.65E+00 -3.23E+00
 
 TH11
+       -1.31E+01 -2.84E+01 -1.05E+02 -3.00E+00  1.08E+01  1.17E+00  5.72E+00  1.09E+01  1.12E+01 -5.34E-02  1.52E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.306
Stop Time:
Wed Sep 29 23:38:29 CDT 2021
