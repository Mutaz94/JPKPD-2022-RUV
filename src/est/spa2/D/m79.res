Thu Sep 30 09:47:53 CDT 2021
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
$DATA ../../../../data/spa2/D/dat79.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   31411.7170809462        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.4861E+02  5.8223E+02 -1.8212E+01  3.6113E+02  1.4568E+02 -2.3622E+03 -1.2238E+03 -6.7740E+01 -1.8896E+03 -7.4678E+02
            -6.0853E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -537.450317974481        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.5387E+00  1.2835E+00  1.0110E+00  1.4572E+00  9.7465E-01  1.7012E+00  1.2488E+00  9.8514E-01  1.1718E+00  1.0273E+00
             1.4736E+01
 PARAMETER:  5.3095E-01  3.4958E-01  1.1096E-01  4.7650E-01  7.4319E-02  6.3133E-01  3.2219E-01  8.5026E-02  2.5852E-01  1.2694E-01
             2.7903E+00
 GRADIENT:   7.1066E+01  1.6350E+01  4.7351E+00  2.7399E+01 -1.5479E+01  1.7963E+01 -2.6846E+01  3.1163E+00 -1.4961E+01  8.1083E+00
             2.7463E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -627.361857461263        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.4463E+00  1.2530E+00  2.0129E+00  1.4758E+00  1.4165E+00  2.3265E+00  5.3504E+00  3.5101E-01  1.5543E+00  1.7331E-01
             1.4186E+01
 PARAMETER:  4.6903E-01  3.2554E-01  7.9959E-01  4.8923E-01  4.4816E-01  9.4435E-01  1.7772E+00 -9.4695E-01  5.4103E-01 -1.6527E+00
             2.7523E+00
 GRADIENT:   2.0616E+01 -2.7536E+00 -9.6663E+00 -8.3165E+00 -3.4749E+00  5.5178E+01  4.4557E+01  1.6819E-01  3.1565E+01  3.0251E-01
             1.9460E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -642.845481049319        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.3557E+00  1.0352E+00  3.3206E+00  1.6116E+00  1.5550E+00  1.9161E+00  5.1420E+00  1.3399E-01  1.4432E+00  6.5965E-01
             1.3468E+01
 PARAMETER:  4.0433E-01  1.3462E-01  1.3002E+00  5.7723E-01  5.4148E-01  7.5030E-01  1.7374E+00 -1.9100E+00  4.6684E-01 -3.1605E-01
             2.7003E+00
 GRADIENT:   9.8121E+00 -1.2837E+00  3.0548E-01  8.8292E+00 -1.9695E+01  1.2884E+01  1.4937E+01  1.5030E-02  2.0917E+01  3.6802E+00
             1.4574E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -651.579098131453        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  1.3456E+00  1.0417E+00  3.3529E+00  1.6035E+00  1.5802E+00  1.9153E+00  5.1440E+00  4.1663E-02  8.7511E-01  2.2640E-01
             1.3057E+01
 PARAMETER:  3.9686E-01  1.4085E-01  1.3098E+00  5.7221E-01  5.5754E-01  7.4985E-01  1.7378E+00 -3.0781E+00 -3.3410E-02 -1.3854E+00
             2.6693E+00
 GRADIENT:   1.1327E+01  4.3838E+00  5.1128E+00  6.1378E+01 -2.7154E+01  1.6157E+01 -1.6045E-01  1.5777E-03 -2.6829E+00  3.7099E-01
             9.0089E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -653.949156276547        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  1.3450E+00  1.0418E+00  3.3474E+00  1.6023E+00  1.5803E+00  1.9154E+00  5.1500E+00  1.1662E-02  1.1467E+00  7.9753E-02
             1.2095E+01
 PARAMETER:  3.9638E-01  1.4096E-01  1.3082E+00  5.7143E-01  5.5762E-01  7.4994E-01  1.7390E+00 -4.3515E+00  2.3685E-01 -2.4288E+00
             2.5928E+00
 GRADIENT:   1.2347E+01  6.3600E+00  3.1980E+00  4.4532E+01 -2.8248E+01  9.7253E-02 -2.8402E+01  1.1258E-04  2.9085E+00  4.8879E-02
            -8.8681E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -673.068044533636        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      622
 NPARAMETR:  1.3424E+00  1.0416E+00  3.3204E+00  1.5734E+00  1.5962E+00  1.9137E+00  6.2669E+00  1.0000E-02  9.3201E-01  1.0000E-02
             1.2112E+01
 PARAMETER:  3.9444E-01  1.4072E-01  1.3001E+00  5.5323E-01  5.6762E-01  7.4904E-01  1.9353E+00 -6.7705E+00  2.9591E-02 -5.4849E+00
             2.5942E+00
 GRADIENT:   1.9007E+01  7.5642E+00  2.6686E+00  3.8172E+01 -2.6430E+01  4.3594E+00 -1.0121E+01  0.0000E+00  2.1872E+00  0.0000E+00
            -1.2083E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -674.289896417635        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      797
 NPARAMETR:  1.3412E+00  1.0414E+00  3.3085E+00  1.5599E+00  1.6039E+00  1.9127E+00  6.8409E+00  1.0000E-02  7.8352E-01  1.0000E-02
             1.2252E+01
 PARAMETER:  3.9355E-01  1.4058E-01  1.2965E+00  5.4465E-01  5.7243E-01  7.4851E-01  2.0229E+00 -7.7062E+00 -1.4396E-01 -6.7324E+00
             2.6057E+00
 GRADIENT:   1.6843E+01  8.7995E+00  2.8232E+00  3.7067E+01 -2.5797E+01  5.1687E+00  5.0885E+00  0.0000E+00 -2.2866E-02  0.0000E+00
            -8.9884E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -674.823389129254        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      980
 NPARAMETR:  1.3361E+00  1.0409E+00  3.2852E+00  1.5341E+00  1.6219E+00  1.9051E+00  6.4939E+00  1.0000E-02  8.8865E-01  1.0000E-02
             1.2112E+01
 PARAMETER:  3.8975E-01  1.4012E-01  1.2894E+00  5.2797E-01  5.8361E-01  7.4451E-01  1.9709E+00 -7.8753E+00 -1.8055E-02 -6.8705E+00
             2.5942E+00
 GRADIENT:   1.7932E+01  6.6950E+00  1.5444E+00  2.3935E+01 -2.1681E+01  2.7823E+00 -2.6629E-01  0.0000E+00  4.4669E+00  0.0000E+00
            -2.9868E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -675.262265121754        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1174
 NPARAMETR:  1.3333E+00  1.0413E+00  3.2675E+00  1.5333E+00  1.6217E+00  1.9086E+00  6.5567E+00  1.0000E-02  7.3416E-01  1.0000E-02
             1.2224E+01
 PARAMETER:  3.8763E-01  1.4051E-01  1.2840E+00  5.2740E-01  5.8345E-01  7.4637E-01  1.9805E+00 -7.8895E+00 -2.0903E-01 -6.8875E+00
             2.6034E+00
 GRADIENT:   1.4348E+01  7.1351E+00  2.6738E+00  3.9652E+01 -2.4093E+01  5.0033E+00 -2.0015E+00  0.0000E+00 -7.2499E-01  0.0000E+00
            -2.0017E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -675.576851920563        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1349             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3324E+00  1.0334E+00  3.2479E+00  1.5288E+00  1.6307E+00  1.9095E+00  6.5779E+00  1.0000E-02  7.3846E-01  1.0000E-02
             1.2274E+01
 PARAMETER:  3.8696E-01  1.3288E-01  1.2780E+00  5.2447E-01  5.8902E-01  7.4685E-01  1.9837E+00 -7.8895E+00 -2.0318E-01 -6.8875E+00
             2.6075E+00
 GRADIENT:   2.8963E+01  6.9943E+00  2.2236E+00  4.2706E+01 -2.1070E+01  1.3799E+01  4.7935E+01  0.0000E+00  3.7797E-01  0.0000E+00
             3.6965E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -675.583279502132        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1483
 NPARAMETR:  1.3320E+00  1.0335E+00  3.2461E+00  1.5284E+00  1.6309E+00  1.9095E+00  6.5398E+00  1.0000E-02  7.3989E-01  1.0000E-02
             1.2217E+01
 PARAMETER:  3.8666E-01  1.3291E-01  1.2775E+00  5.2422E-01  5.8914E-01  7.4684E-01  1.9779E+00 -7.8895E+00 -2.0125E-01 -6.8875E+00
             2.6029E+00
 GRADIENT:   1.4093E+01  6.4606E+00  2.1085E+00  3.6482E+01 -2.2243E+01  5.0172E+00 -1.9242E+00  0.0000E+00 -8.1337E-02  0.0000E+00
             5.6228E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -676.626462071246        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1664
 NPARAMETR:  1.3077E+00  1.0346E+00  3.1278E+00  1.4868E+00  1.6630E+00  1.9076E+00  7.0843E+00  1.0000E-02  6.7205E-01  1.0000E-02
             1.1864E+01
 PARAMETER:  3.6826E-01  1.3400E-01  1.2403E+00  4.9663E-01  6.0863E-01  7.4586E-01  2.0579E+00 -7.8895E+00 -2.9743E-01 -6.8875E+00
             2.5735E+00
 GRADIENT:   1.2075E+01  7.5681E+00  2.2650E-01  2.7255E+01 -1.7103E+01  2.8376E+00  1.4548E+01  0.0000E+00  6.9992E-01  0.0000E+00
            -2.1731E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -680.148721802627        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1847
 NPARAMETR:  1.1996E+00  1.0431E+00  2.6217E+00  1.3514E+00  1.7648E+00  1.9406E+00  6.7478E+00  1.0000E-02  1.9566E-01  1.0000E-02
             1.1972E+01
 PARAMETER:  2.8202E-01  1.4215E-01  1.0638E+00  4.0111E-01  6.6805E-01  7.6299E-01  2.0092E+00 -7.8895E+00 -1.5314E+00 -6.8875E+00
             2.5826E+00
 GRADIENT:  -2.3047E+01  6.4704E-01 -4.2036E+00  1.9177E+01 -4.1560E+00  1.8864E+01  1.2234E+01  0.0000E+00 -1.1353E-01  0.0000E+00
             1.0437E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -682.427455769613        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:     1970
 NPARAMETR:  1.1985E+00  1.1297E+00  4.8885E+00  1.2773E+00  2.1239E+00  1.9150E+00  6.6117E+00  1.0000E-02  1.2512E-01  1.0000E-02
             1.1799E+01
 PARAMETER:  2.8108E-01  2.2197E-01  1.6869E+00  3.4471E-01  8.5323E-01  7.4974E-01  1.9888E+00 -7.8895E+00 -1.9785E+00 -6.8875E+00
             2.5680E+00
 GRADIENT:  -7.0767E+00  6.1887E-01 -3.2630E-01 -1.6070E+01 -9.3974E-01  2.3194E+01  7.4630E+01  0.0000E+00  2.3235E-01  0.0000E+00
             4.7481E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -683.361122439261        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:     2088
 NPARAMETR:  1.1943E+00  1.1679E+00  9.6424E+00  1.2923E+00  2.4785E+00  1.8764E+00  6.3891E+00  1.0000E-02  1.2609E-01  1.0000E-02
             1.1606E+01
 PARAMETER:  2.7757E-01  2.5524E-01  2.3662E+00  3.5644E-01  1.0076E+00  7.2935E-01  1.9546E+00 -7.8895E+00 -1.9708E+00 -6.8875E+00
             2.5515E+00
 GRADIENT:  -1.8501E+01  6.5736E-02 -1.9449E-01 -6.0224E+00  1.4402E+00  9.3607E+00  9.1197E+00  0.0000E+00  1.1556E-01  0.0000E+00
            -1.0039E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -684.363419880327        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2267
 NPARAMETR:  1.1898E+00  1.2868E+00  2.9250E+01  1.2779E+00  2.6340E+00  1.7750E+00  5.9117E+00  1.0000E-02  3.8409E-02  1.0000E-02
             1.1854E+01
 PARAMETER:  2.7382E-01  3.5216E-01  3.4759E+00  3.4524E-01  1.0685E+00  6.7380E-01  1.8769E+00 -7.8895E+00 -3.1595E+00 -6.8875E+00
             2.5726E+00
 GRADIENT:  -1.5998E+01 -2.8100E-01  3.5344E-02 -5.7824E+00 -7.4011E-01  1.4598E+01 -4.0075E+00  0.0000E+00  1.5104E-02  0.0000E+00
             2.6668E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -686.008404536138        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2444
 NPARAMETR:  1.1996E+00  1.2180E+00  3.3978E+01  1.3284E+00  2.7446E+00  1.5866E+00  6.1728E+00  1.0000E-02  1.0000E-02  1.0000E-02
             1.2018E+01
 PARAMETER:  2.8196E-01  2.9722E-01  3.6257E+00  3.8397E-01  1.1096E+00  5.6157E-01  1.9201E+00 -7.8895E+00 -8.5100E+00 -6.8875E+00
             2.5864E+00
 GRADIENT:  -3.2501E+00  2.1084E-01 -2.8582E-02 -1.2697E+00  9.4317E-02 -2.5021E+00  1.1382E-02  0.0000E+00  0.0000E+00  0.0000E+00
            -5.5716E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -686.093008787290        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     2618
 NPARAMETR:  1.2077E+00  1.1990E+00  4.5694E+01  1.3380E+00  2.7501E+00  1.6045E+00  6.2855E+00  1.0000E-02  1.0000E-02  1.0000E-02
             1.2040E+01
 PARAMETER:  2.8872E-01  2.8152E-01  3.9220E+00  3.9117E-01  1.1116E+00  5.7283E-01  1.9382E+00 -7.8895E+00 -8.2978E+00 -6.8875E+00
             2.5883E+00
 GRADIENT:   5.3121E-02  2.7953E-01 -9.7912E-03 -6.8348E-01 -4.3368E-01  5.1969E-01  2.6378E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -5.1278E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -686.109767204588        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2800
 NPARAMETR:  1.2092E+00  1.1876E+00  9.5147E+01  1.3444E+00  2.7911E+00  1.6035E+00  6.3173E+00  1.0000E-02  1.0000E-02  1.0000E-02
             1.2059E+01
 PARAMETER:  2.8999E-01  2.7194E-01  4.6554E+00  3.9592E-01  1.1264E+00  5.7218E-01  1.9433E+00 -7.8895E+00 -8.2978E+00 -6.8875E+00
             2.5898E+00
 GRADIENT:   3.9124E-01  9.4370E-02 -3.7177E-03 -3.7489E-01 -1.3883E-01  4.2382E-01  3.0304E+00  0.0000E+00  0.0000E+00  0.0000E+00
             1.2053E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -686.118270616007        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2976
 NPARAMETR:  1.2074E+00  1.1738E+00  2.6930E+02  1.3464E+00  2.8085E+00  1.6021E+00  6.3533E+00  1.0000E-02  1.0000E-02  1.0000E-02
             1.2039E+01
 PARAMETER:  2.8849E-01  2.6024E-01  5.6958E+00  3.9740E-01  1.1326E+00  5.7129E-01  1.9490E+00 -7.8895E+00 -8.2978E+00 -6.8875E+00
             2.5882E+00
 GRADIENT:  -3.0516E-03 -1.1942E-01 -6.1311E-04 -7.3266E-01 -5.8835E-02  1.8737E-01  3.7384E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -1.1142E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -686.119689114881        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     3167
 NPARAMETR:  1.2080E+00  1.1754E+00  8.5794E+02  1.3486E+00  2.8201E+00  1.6009E+00  6.3655E+00  1.0000E-02  1.0000E-02  1.0000E-02
             1.2038E+01
 PARAMETER:  2.8901E-01  2.6159E-01  6.8545E+00  3.9909E-01  1.1368E+00  5.7054E-01  1.9509E+00 -7.8895E+00 -8.2978E+00 -6.8875E+00
             2.5881E+00
 GRADIENT:   2.6557E-01  8.1620E-02 -2.3417E-04  1.7477E-01 -3.8658E-02 -3.7273E-02  3.9579E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -9.7987E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -686.119950585222        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     3347
 NPARAMETR:  1.2079E+00  1.1743E+00  6.8896E+03  1.3491E+00  2.8208E+00  1.6015E+00  6.3677E+00  1.0000E-02  1.0000E-02  1.0000E-02
             1.2040E+01
 PARAMETER:  2.8888E-01  2.6064E-01  8.9378E+00  3.9942E-01  1.1370E+00  5.7092E-01  1.9512E+00 -7.8895E+00 -8.2978E+00 -6.8875E+00
             2.5883E+00
 GRADIENT:   9.5080E-02  5.9470E-02 -1.8339E-05  2.3781E-01 -8.7320E-02  1.0291E-01  3.9600E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -6.8944E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     3347
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.0246E-02  2.3162E-02  1.0750E-08 -1.0161E-03  3.0265E-05
 SE:             2.7616E-02  2.4688E-02  4.0035E-09  1.9093E-04  3.5457E-05
 N:                     100         100         100         100         100

 P VAL.:         2.7341E-01  3.4815E-01  7.2490E-03  1.0295E-07  3.9334E-01

 ETASHRINKSD(%)  7.4846E+00  1.7291E+01  1.0000E+02  9.9360E+01  9.9881E+01
 ETASHRINKVR(%)  1.4409E+01  3.1593E+01  1.0000E+02  9.9996E+01  1.0000E+02
 EBVSHRINKSD(%)  1.1108E+01  1.1808E+01  1.0000E+02  9.9523E+01  9.9817E+01
 EBVSHRINKVR(%)  2.0983E+01  2.2222E+01  1.0000E+02  9.9998E+01  1.0000E+02
 RELATIVEINF(%)  7.5364E+01  4.6731E+01  0.0000E+00  1.1031E-03  6.5284E-05
 EPSSHRINKSD(%)  4.0650E+00
 EPSSHRINKVR(%)  7.9648E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -686.11995058522166     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       416.60628926038544     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    84.72
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    12.66
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -686.120       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.21E+00  1.17E+00  6.89E+03  1.35E+00  2.82E+00  1.60E+00  6.37E+00  1.00E-02  1.00E-02  1.00E-02  1.20E+01
 


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
+        2.20E+02
 
 TH 2
+       -1.31E+01  3.52E+00
 
 TH 3
+       -1.44E-06  1.25E-07  1.08E-14
 
 TH 4
+       -1.37E+02  3.53E+01  1.37E-06  3.63E+02
 
 TH 5
+        6.52E+00 -1.78E+00 -6.86E-08 -1.85E+01  9.51E-01
 
 TH 6
+       -4.70E+01  2.31E+00  4.94E-07  4.58E+01 -2.75E+00  6.04E+01
 
 TH 7
+        7.27E+00 -1.65E+00 -6.97E-08 -1.71E+01  8.74E-01 -2.58E+00  8.13E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -4.97E-01 -1.82E+00 -2.70E-08 -1.77E+01  9.06E-01  4.58E-01  7.87E-01  0.00E+00  0.00E+00  0.00E+00  1.47E+00
 
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
+        2.06E+02
 
 TH 2
+       -4.20E+00  1.57E+01
 
 TH 3
+       -1.30E-06  1.89E-06 -9.96E-12
 
 TH 4
+       -1.97E+01  2.51E+01  3.72E-07  2.62E+02
 
 TH 5
+        2.70E-01 -1.64E+00 -5.20E-07 -1.33E+01  4.87E+00
 
 TH 6
+       -1.26E+01 -1.28E+00  2.21E-07  1.15E+01 -1.03E+00  5.40E+01
 
 TH 7
+        1.39E+00  1.77E+00 -3.11E-08 -1.20E+01  4.26E-01 -9.18E-01  2.87E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.24E+00 -1.89E+00 -1.11E-08 -1.68E+01  7.66E-01  8.52E-01  5.88E-01  0.00E+00  0.00E+00  0.00E+00  4.98E+00
 
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
+        2.42E+02
 
 TH 2
+        3.43E+01  1.57E+01
 
 TH 3
+       -1.59E-07 -2.44E-08  3.65E-15
 
 TH 4
+        7.34E+01  2.60E+01 -1.78E-07  2.27E+02
 
 TH 5
+       -5.70E+00 -7.62E-01 -4.24E-08 -1.11E+01  1.70E+00
 
 TH 6
+        1.50E+01  2.49E+00 -4.32E-09 -3.32E+00 -4.52E-01  5.41E+01
 
 TH 7
+       -4.37E+00  2.77E+00  7.42E-09 -1.30E+01  1.27E+00  5.61E-02  4.07E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.02E+01 -5.17E+00  1.55E-08 -2.28E+01  1.90E+00 -5.10E+00  3.09E+00  0.00E+00  0.00E+00  0.00E+00  5.27E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       97.496
Stop Time:
Thu Sep 30 09:49:32 CDT 2021
