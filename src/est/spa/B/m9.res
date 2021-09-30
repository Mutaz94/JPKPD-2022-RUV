Wed Sep 29 10:56:53 CDT 2021
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
$DATA ../../../../data/spa/B/dat9.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m9.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1632.87432452960        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0997E+02 -2.8251E+01 -2.5027E+01  9.3370E+00  7.1826E+01  2.1530E+01  3.1522E+00  1.9823E+00  9.9784E+00 -5.0099E+00
             5.2332E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1640.64911797467        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      112
 NPARAMETR:  9.2696E-01  9.9263E-01  9.9042E-01  9.9854E-01  9.8027E-01  9.3477E-01  9.9196E-01  9.9862E-01  9.8894E-01  9.9656E-01
             9.2411E-01
 PARAMETER:  2.4154E-02  9.2605E-02  9.0378E-02  9.8541E-02  8.0074E-02  3.2542E-02  9.1927E-02  9.8621E-02  8.8874E-02  9.6550E-02
             2.1077E-02
 GRADIENT:  -2.0832E+01 -5.5710E+01 -2.0319E+01 -5.5503E+01  5.0100E+01 -2.3259E+01 -6.7496E-01  1.7092E+00  1.4912E+00 -4.7685E+00
             2.6462E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1644.61426258438        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      289
 NPARAMETR:  9.3238E-01  9.1231E-01  9.5710E-01  1.0928E+00  8.8437E-01  1.0121E+00  1.0215E+00  9.0476E-01  9.4753E-01  9.4377E-01
             8.6359E-01
 PARAMETER:  2.9989E-02  8.2210E-03  5.6148E-02  1.8878E-01 -2.2879E-02  1.1204E-01  1.2128E-01 -8.0215E-05  4.6102E-02  4.2126E-02
            -4.6652E-02
 GRADIENT:  -4.5338E+00  6.2469E+00  1.6327E+00  1.0167E+01 -2.5910E+00  8.8122E+00 -2.2692E+00 -1.7924E-01  2.3258E+00 -2.9864E+00
             1.0087E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1645.07931597212        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      471
 NPARAMETR:  9.3448E-01  8.6405E-01  8.7549E-01  1.1118E+00  8.2613E-01  9.8844E-01  1.2184E+00  7.5566E-01  8.7599E-01  9.0157E-01
             8.6326E-01
 PARAMETER:  3.2234E-02 -4.6120E-02 -3.2969E-02  2.0595E-01 -9.1000E-02  8.8377E-02  2.9756E-01 -1.8016E-01 -3.2401E-02 -3.6227E-03
            -4.7036E-02
 GRADIENT:  -8.1898E-01  3.7032E+00  2.2461E+00  2.5496E+00 -3.1368E+00 -4.4126E-01 -2.5918E-01 -3.8160E-01 -3.6202E-01 -5.6465E-01
             1.6206E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1645.18808268558        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      647
 NPARAMETR:  9.3502E-01  6.8543E-01  1.0258E+00  1.2344E+00  8.3034E-01  9.8896E-01  1.3700E+00  8.9276E-01  8.3360E-01  9.4565E-01
             8.5623E-01
 PARAMETER:  3.2815E-02 -2.7772E-01  1.2544E-01  3.1056E-01 -8.5915E-02  8.8896E-02  4.1478E-01 -1.3437E-02 -8.2003E-02  4.4119E-02
            -5.5213E-02
 GRADIENT:   6.4068E+00  5.7924E+00  3.0328E+00  7.4207E+00 -7.7002E+00  7.0585E-01  3.2697E-02  2.4791E-01 -3.5447E-01  4.9314E-01
            -1.7169E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1645.24328573078        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      822
 NPARAMETR:  9.3366E-01  5.7599E-01  1.1431E+00  1.3108E+00  8.4681E-01  9.8930E-01  1.4811E+00  9.9416E-01  8.1517E-01  9.9162E-01
             8.5209E-01
 PARAMETER:  3.1362E-02 -4.5166E-01  2.3375E-01  3.7061E-01 -6.6275E-02  8.9240E-02  4.9276E-01  9.4145E-02 -1.0436E-01  9.1586E-02
            -6.0066E-02
 GRADIENT:   7.6819E+00  6.7394E+00  2.6025E+00  1.1234E+01 -8.0478E+00  1.5210E+00  9.7227E-01  8.3659E-01  5.5253E-01  2.1858E+00
            -3.7128E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1645.34376186434        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      997
 NPARAMETR:  9.3057E-01  4.4029E-01  1.2473E+00  1.4013E+00  8.5062E-01  9.8754E-01  1.7275E+00  1.0736E+00  7.8390E-01  1.0225E+00
             8.5087E-01
 PARAMETER:  2.8045E-02 -7.2033E-01  3.2096E-01  4.3737E-01 -6.1784E-02  8.7460E-02  6.4670E-01  1.7104E-01 -1.4348E-01  1.2222E-01
            -6.1500E-02
 GRADIENT:   5.3761E+00  6.5494E+00  1.7606E+00  1.5338E+01 -5.8579E+00  1.5872E+00  1.5407E+00  6.9168E-01  7.3079E-01  2.9108E+00
            -4.3921E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1645.52570327487        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1172
 NPARAMETR:  9.2599E-01  2.8477E-01  1.3532E+00  1.5020E+00  8.4927E-01  9.8360E-01  2.2026E+00  1.1739E+00  7.4966E-01  1.0360E+00
             8.5485E-01
 PARAMETER:  2.3109E-02 -1.1561E+00  4.0247E-01  5.0681E-01 -6.3376E-02  8.3462E-02  8.8966E-01  2.6032E-01 -1.8813E-01  1.3541E-01
            -5.6829E-02
 GRADIENT:   8.8794E-02  4.8487E+00  3.6981E-01  1.6717E+01 -2.6739E+00  8.7128E-01  1.5018E+00  4.7718E-01  3.0883E-01  2.4708E+00
            -2.5962E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1645.59236596524        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1347
 NPARAMETR:  9.2347E-01  1.9028E-01  1.4133E+00  1.5637E+00  8.4403E-01  9.8005E-01  2.7046E+00  1.2440E+00  7.2927E-01  1.0273E+00
             8.5982E-01
 PARAMETER:  2.0380E-02 -1.5592E+00  4.4594E-01  5.4708E-01 -6.9570E-02  7.9851E-02  1.0950E+00  3.1836E-01 -2.1572E-01  1.2692E-01
            -5.1037E-02
 GRADIENT:  -2.5901E+00  3.6977E+00  6.4100E-01  2.0089E+01 -3.1010E+00 -2.9768E-02  9.9933E-01  1.9617E-01 -9.7853E-01  1.1894E+00
            -4.5036E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1645.67638594525        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1522
 NPARAMETR:  9.2206E-01  1.1538E-01  1.4742E+00  1.6102E+00  8.4411E-01  9.7794E-01  3.3284E+00  1.3201E+00  7.1781E-01  1.0205E+00
             8.6337E-01
 PARAMETER:  1.8859E-02 -2.0595E+00  4.8810E-01  5.7633E-01 -6.9469E-02  7.7696E-02  1.3025E+00  3.7772E-01 -2.3155E-01  1.2024E-01
            -4.6908E-02
 GRADIENT:  -3.0917E+00  2.0689E+00  6.7355E-01  1.4563E+01 -2.8049E+00 -4.5475E-01  4.3356E-01  2.8900E-01 -5.9510E-01  2.7058E-01
             1.0995E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1645.85564542424        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1704
 NPARAMETR:  9.2295E-01  8.3818E-02  1.4865E+00  1.6186E+00  8.4344E-01  9.7865E-01  3.1690E+00  1.3331E+00  7.1661E-01  1.0205E+00
             8.6092E-01
 PARAMETER:  1.9823E-02 -2.3791E+00  4.9644E-01  5.8155E-01 -7.0264E-02  7.8414E-02  1.2534E+00  3.8751E-01 -2.3322E-01  1.2025E-01
            -4.9749E-02
 GRADIENT:   6.3290E-01  3.9093E-01 -6.4994E-01 -1.3761E+01  2.5798E+00  4.8989E-02  1.0803E-01 -2.0836E-02 -1.3304E-02 -3.6516E-01
             7.8789E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1645.95132737687        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1886
 NPARAMETR:  9.2299E-01  7.3793E-02  1.4879E+00  1.6235E+00  8.3984E-01  9.7853E-01  3.9063E-01  1.3353E+00  7.1802E-01  1.0243E+00
             8.6072E-01
 PARAMETER:  1.9861E-02 -2.5065E+00  4.9737E-01  5.8456E-01 -7.4540E-02  7.8292E-02 -8.4000E-01  3.8914E-01 -2.3126E-01  1.2402E-01
            -4.9984E-02
 GRADIENT:   1.1074E+00  2.1654E-01  4.6942E-01 -1.4579E+01 -6.2078E-01  8.4933E-02  5.2011E-03 -4.5875E-02 -1.3331E+00  1.1220E-01
             4.2534E-03

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1645.95869150773        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:     2034
 NPARAMETR:  9.2269E-01  7.4923E-02  1.4867E+00  1.6218E+00  8.3951E-01  9.7822E-01  3.4322E-02  1.3354E+00  7.2047E-01  1.0237E+00
             8.6049E-01
 PARAMETER:  1.9533E-02 -2.4913E+00  4.9655E-01  5.8354E-01 -7.4939E-02  7.7976E-02 -3.2720E+00  3.8924E-01 -2.2785E-01  1.2338E-01
            -5.0248E-02
 GRADIENT:   4.8015E+02  6.0557E+00  1.2565E+01  1.1218E+03  8.7117E+00  3.9525E+01  9.0038E-04  1.9880E+00  2.2574E+01  1.3930E+00
             7.6334E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1645.96050807887        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2217
 NPARAMETR:  9.2323E-01  7.6041E-02  1.4852E+00  1.6213E+00  8.3953E-01  9.7855E-01  2.2984E-02  1.3344E+00  7.2127E-01  1.0234E+00
             8.6055E-01
 PARAMETER:  2.0123E-02 -2.4765E+00  4.9553E-01  5.8324E-01 -7.4908E-02  7.8315E-02 -3.6730E+00  3.8847E-01 -2.2675E-01  1.2310E-01
            -5.0179E-02
 GRADIENT:   1.6425E+00  1.3279E-01  4.1040E-01 -1.5994E+01 -5.0586E-01  8.3299E-02  3.5777E-05  6.0356E-02 -7.5738E-02  6.7367E-02
            -7.7309E-03

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1645.96205945675        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2406             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2325E-01  7.7205E-02  1.4825E+00  1.6197E+00  8.4009E-01  9.7858E-01  1.0000E-02  1.3326E+00  7.2199E-01  1.0228E+00
             8.6058E-01
 PARAMETER:  2.0140E-02 -2.4613E+00  4.9371E-01  5.8227E-01 -7.4250E-02  7.8347E-02 -5.0127E+00  3.8713E-01 -2.2575E-01  1.2252E-01
            -5.0147E-02
 GRADIENT:   4.8135E+02  6.2495E+00  1.1486E+01  1.1153E+03  1.0694E+01  3.9589E+01  0.0000E+00  2.0288E+00  2.2734E+01  1.1468E+00
             8.2382E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1645.96323785283        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2593
 NPARAMETR:  9.2326E-01  7.8339E-02  1.4815E+00  1.6190E+00  8.4009E-01  9.7859E-01  1.0000E-02  1.3316E+00  7.2232E-01  1.0227E+00
             8.6057E-01
 PARAMETER:  2.0160E-02 -2.4467E+00  4.9308E-01  5.8182E-01 -7.4250E-02  7.8360E-02 -5.0127E+00  3.8636E-01 -2.2529E-01  1.2241E-01
            -5.0157E-02
 GRADIENT:   1.6397E+00  3.6157E-02 -4.1801E-01 -1.7842E+01  1.2008E+00  8.6822E-02  0.0000E+00  8.8308E-02  9.1596E-02 -1.5014E-01
             1.6483E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1645.96500837343        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     2783             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2328E-01  8.0989E-02  1.4814E+00  1.6183E+00  8.3925E-01  9.7860E-01  1.0000E-02  1.3299E+00  7.2265E-01  1.0238E+00
             8.6055E-01
 PARAMETER:  2.0174E-02 -2.4134E+00  4.9300E-01  5.8140E-01 -7.5250E-02  7.8371E-02 -5.0127E+00  3.8512E-01 -2.2483E-01  1.2351E-01
            -5.0184E-02
 GRADIENT:   4.8132E+02  6.8855E+00  1.2558E+01  1.1126E+03  8.4371E+00  3.9596E+01  0.0000E+00  1.9687E+00  2.2440E+01  1.4247E+00
             7.8804E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1645.96617204253        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2966
 NPARAMETR:  9.2330E-01  8.2071E-02  1.4805E+00  1.6176E+00  8.3925E-01  9.7862E-01  1.0000E-02  1.3290E+00  7.2298E-01  1.0238E+00
             8.6055E-01
 PARAMETER:  2.0195E-02 -2.4002E+00  4.9238E-01  5.8096E-01 -7.5250E-02  7.8384E-02 -5.0127E+00  3.8439E-01 -2.2437E-01  1.2354E-01
            -5.0184E-02
 GRADIENT:   1.5485E+00  1.7798E-01  6.4251E-01 -1.5074E+01 -1.0784E+00  7.9257E-02  0.0000E+00  4.6735E-02 -1.3146E-01  1.4116E-01
            -1.1250E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1645.96773269512        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     3156             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2334E-01  8.2689E-02  1.4779E+00  1.6161E+00  8.4009E-01  9.7864E-01  1.0000E-02  1.3276E+00  7.2366E-01  1.0227E+00
             8.6057E-01
 PARAMETER:  2.0241E-02 -2.3927E+00  4.9063E-01  5.8003E-01 -7.4250E-02  7.8413E-02 -5.0127E+00  3.8337E-01 -2.2344E-01  1.2247E-01
            -5.0163E-02
 GRADIENT:   4.8139E+02  6.9636E+00  1.1406E+01  1.1054E+03  1.0785E+01  3.9597E+01  0.0000E+00  1.9876E+00  2.2485E+01  1.1080E+00
             8.1438E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1645.96887978694        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3343
 NPARAMETR:  9.2336E-01  8.3784E-02  1.4771E+00  1.6154E+00  8.4009E-01  9.7866E-01  1.0000E-02  1.3267E+00  7.2399E-01  1.0228E+00
             8.6057E-01
 PARAMETER:  2.0261E-02 -2.3795E+00  4.9006E-01  5.7960E-01 -7.4250E-02  7.8426E-02 -5.0127E+00  3.8266E-01 -2.2298E-01  1.2253E-01
            -5.0164E-02
 GRADIENT:   1.6391E+00  3.5421E-02 -4.5517E-01 -1.7653E+01  1.2313E+00  8.6026E-02  0.0000E+00  7.9529E-02  9.6384E-02 -1.6498E-01
             1.3940E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1645.97056328783        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     3533             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2337E-01  8.6444E-02  1.4771E+00  1.6148E+00  8.3925E-01  9.7867E-01  1.0000E-02  1.3251E+00  7.2430E-01  1.0240E+00
             8.6055E-01
 PARAMETER:  2.0274E-02 -2.3483E+00  4.9010E-01  5.7920E-01 -7.5250E-02  7.8437E-02 -5.0127E+00  3.8151E-01 -2.2255E-01  1.2376E-01
            -5.0188E-02
 GRADIENT:   4.8133E+02  7.6274E+00  1.2544E+01  1.1030E+03  8.3693E+00  3.9602E+01  0.0000E+00  1.9384E+00  2.2189E+01  1.4375E+00
             7.8725E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1645.97155956818        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3716
 NPARAMETR:  9.2339E-01  8.7534E-02  1.4764E+00  1.6141E+00  8.3925E-01  9.7868E-01  1.0000E-02  1.3242E+00  7.2463E-01  1.0242E+00
             8.6055E-01
 PARAMETER:  2.0293E-02 -2.3357E+00  4.8961E-01  5.7877E-01 -7.5250E-02  7.8450E-02 -5.0127E+00  3.8083E-01 -2.2210E-01  1.2387E-01
            -5.0188E-02
 GRADIENT:   1.5443E+00  1.9657E-01  7.0829E-01 -1.4744E+01 -1.2485E+00  7.8836E-02  0.0000E+00  3.4781E-02 -1.3794E-01  1.6099E-01
            -1.2937E-02

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1645.97355688028        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     3906             RESET HESSIAN, TYPE I
 NPARAMETR:  9.2342E-01  8.8197E-02  1.4743E+00  1.6128E+00  8.4009E-01  9.7871E-01  1.0000E-02  1.3232E+00  7.2522E-01  1.0232E+00
             8.6056E-01
 PARAMETER:  2.0333E-02 -2.3282E+00  4.8816E-01  5.7798E-01 -7.4250E-02  7.8475E-02 -5.0127E+00  3.8004E-01 -2.2129E-01  1.2291E-01
            -5.0170E-02
 GRADIENT:   4.8139E+02  7.7400E+00  1.1550E+01  1.0966E+03  1.0388E+01  3.9604E+01  0.0000E+00  1.9501E+00  2.2204E+01  1.1584E+00
             8.0597E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1645.97467740594        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     4093
 NPARAMETR:  9.2344E-01  8.9310E-02  1.4736E+00  1.6121E+00  8.4012E-01  9.7872E-01  1.0000E-02  1.3223E+00  7.2554E-01  1.0233E+00
             8.6056E-01
 PARAMETER:  2.0352E-02 -2.3156E+00  4.8772E-01  5.7755E-01 -7.4207E-02  7.8488E-02 -5.0127E+00  3.7941E-01 -2.2084E-01  1.2303E-01
            -5.0170E-02
 GRADIENT:   1.6124E+00  7.3602E-02 -2.5783E-01 -1.6732E+01  7.7570E-01  8.6469E-02  0.0000E+00  6.0098E-02  3.6623E-02 -1.1433E-01
             4.6464E-03

0ITERATION NO.:  116    OBJECTIVE VALUE:  -1645.97467740594        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     4121
 NPARAMETR:  9.2344E-01  9.0807E-02  1.4744E+00  1.6119E+00  8.3934E-01  9.7872E-01  1.0000E-02  1.3216E+00  7.2565E-01  1.0246E+00
             8.6055E-01
 PARAMETER:  2.0352E-02 -2.3156E+00  4.8772E-01  5.7755E-01 -7.4207E-02  7.8488E-02 -5.0127E+00  3.7941E-01 -2.2084E-01  1.2303E-01
            -5.0170E-02
 GRADIENT:  -3.4466E-03 -4.3361E-02 -2.1784E-01  4.2158E-01  4.8553E-01 -1.5011E-03  0.0000E+00  5.0914E-02 -4.9153E-02 -7.3538E-02
             5.0711E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4121
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7966E-04 -2.9206E-05 -3.4054E-02 -6.4476E-03 -3.8953E-02
 SE:             2.9864E-02  1.5948E-05  1.8665E-02  2.9371E-02  2.1144E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9520E-01  6.7053E-02  6.8078E-02  8.2624E-01  6.5440E-02

 ETASHRINKSD(%)  1.0000E-10  9.9947E+01  3.7470E+01  1.6032E+00  2.9164E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.0899E+01  3.1807E+00  4.9822E+01
 EBVSHRINKSD(%)  3.2407E-01  9.9949E+01  4.0872E+01  2.0808E+00  2.5501E+01
 EBVSHRINKVR(%)  6.4709E-01  1.0000E+02  6.5039E+01  4.1184E+00  4.4499E+01
 RELATIVEINF(%)  9.7251E+01  1.4577E-06  7.6921E+00  6.4391E+00  7.9470E+00
 EPSSHRINKSD(%)  4.6402E+01
 EPSSHRINKVR(%)  7.1272E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1645.9746774059354     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -910.82385084219720     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    55.28
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.60
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1645.975       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.23E-01  8.93E-02  1.47E+00  1.61E+00  8.40E-01  9.79E-01  1.00E-02  1.32E+00  7.26E-01  1.02E+00  8.61E-01
 


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
+        1.35E+03
 
 TH 2
+       -2.90E+01  3.80E+02
 
 TH 3
+       -4.07E-01  4.40E+01  1.09E+02
 
 TH 4
+       -7.89E+00  4.85E+02 -2.51E+01  7.79E+02
 
 TH 5
+        1.38E+00 -2.42E+02 -2.42E+02 -6.43E+01  8.60E+02
 
 TH 6
+        5.34E-01 -3.57E+00  1.22E-01 -2.30E+00 -1.41E+00  2.05E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.56E-01  7.21E-01 -2.69E+01 -4.11E+00 -6.40E+00  2.16E-02  0.00E+00  2.90E+01
 
 TH 9
+        2.61E+00 -1.03E+02  5.92E+00 -8.94E-01  8.74E-01 -9.51E-01  0.00E+00 -5.52E-01  3.50E+02
 
 TH10
+        6.44E-01  9.61E+00 -1.64E+00 -1.09E+00 -8.34E+01  1.07E-01  0.00E+00  1.66E+01  1.92E+00  6.65E+01
 
 TH11
+       -8.66E+00 -9.66E+00 -7.48E+00 -8.59E+00 -4.97E+00  2.51E+00  0.00E+00  6.48E+00  1.19E+01  1.03E+01  2.75E+02
 
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
 #CPUT: Total CPU Time in Seconds,       60.896
Stop Time:
Wed Sep 29 10:57:55 CDT 2021
