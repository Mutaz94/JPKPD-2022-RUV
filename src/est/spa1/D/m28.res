Thu Sep 30 02:52:36 CDT 2021
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
$DATA ../../../../data/spa1/D/dat28.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   11863.9354677993        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6741E+02 -5.7694E+01 -3.6878E+01 -1.3812E+02  1.8995E+02 -9.4301E+02 -5.9568E+02 -5.8107E+01 -8.6165E+02 -2.1088E+02
            -2.5039E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -736.884116617083        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.5337E+00  1.0697E+00  1.0240E+00  1.7478E+00  1.0843E+00  2.7138E+00  1.3402E+00  9.6346E-01  1.9492E+00  1.0693E+00
             1.3390E+01
 PARAMETER:  5.2771E-01  1.6740E-01  1.2371E-01  6.5837E-01  1.8093E-01  1.0983E+00  3.9283E-01  6.2772E-02  7.6742E-01  1.6704E-01
             2.6945E+00
 GRADIENT:   4.1539E+01  4.6740E-01 -1.5846E+01 -1.3091E+01 -1.8719E+00  9.7150E+01  2.5926E+00  6.7304E+00  1.0441E+01  5.1927E+00
             3.2785E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -771.529370100425        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.5192E+00  5.8708E-01  2.6711E+00  2.6412E+00  1.7690E+00  2.6159E+00  5.9908E+00  6.7995E-01  2.0294E+00  2.8988E+00
             1.1955E+01
 PARAMETER:  5.1816E-01 -4.3259E-01  1.0825E+00  1.0712E+00  6.7041E-01  1.0616E+00  1.8902E+00 -2.8574E-01  8.0774E-01  1.1643E+00
             2.5811E+00
 GRADIENT:   4.7066E+01  1.2191E+01  1.1431E+00  7.1719E+01 -1.6030E+01  7.6695E+01  6.0151E+00  5.1207E-01  3.2234E+01  1.0264E+01
             2.7274E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -857.800031388852        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0106E+00  6.3051E-01  2.8581E+00  1.5939E+00  3.7694E+00  1.7972E+00  3.2536E+00  7.0495E-01  1.6622E+00  4.7997E+00
             8.0131E+00
 PARAMETER:  1.1057E-01 -3.6123E-01  1.1502E+00  5.6616E-01  1.4269E+00  6.8624E-01  1.2798E+00 -2.4963E-01  6.0813E-01  1.6685E+00
             2.1811E+00
 GRADIENT:  -8.5097E+01 -1.8152E+00  1.8815E+01 -5.1723E+01 -7.7468E+00 -2.6144E+01  8.2479E+00 -4.3473E-01 -1.5910E+01 -5.3117E-01
             9.6248E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -911.889813881361        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.0315E+00  3.5906E-01  4.2689E-01  1.4165E+00  1.7998E+01  1.7842E+00  8.6076E-01  1.0000E-02  1.3667E+00  9.6848E+00
             7.5058E+00
 PARAMETER:  1.3101E-01 -9.2427E-01 -7.5122E-01  4.4821E-01  2.9902E+00  6.7896E-01 -4.9939E-02 -6.2271E+00  4.1237E-01  2.3706E+00
             2.1157E+00
 GRADIENT:   2.7236E+01  5.0701E+01 -3.7937E+01  1.4224E+01 -6.6524E+00  7.1063E+00  3.2292E+00  0.0000E+00  3.8282E+00  6.8736E+00
            -6.3744E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -986.697968127806        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  6.0965E-01  1.4702E-02  7.3252E-02  7.0288E-01  1.0262E+03  1.4558E+00  1.6847E-02  1.0000E-02  7.2354E-01  5.9690E+01
             7.3110E+00
 PARAMETER: -3.9487E-01 -4.1197E+00 -2.5138E+00 -2.5257E-01  7.0336E+00  4.7555E-01 -3.9836E+00 -1.5468E+01 -2.2360E-01  4.1892E+00
             2.0894E+00
 GRADIENT:  -3.2158E+00  7.6016E-01 -5.0591E+01  2.2816E+02  6.2283E-03 -1.4771E+01  6.4836E-07  0.0000E+00 -2.8774E+01 -3.2426E-04
            -1.8813E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -986.803241301202        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      447
 NPARAMETR:  5.9310E-01  1.2387E-02  6.5777E-02  6.5760E-01  1.3566E+03  1.4509E+00  1.2894E-02  1.0000E-02  6.8562E-01  6.8304E+01
             7.2326E+00
 PARAMETER: -4.2239E-01 -4.2911E+00 -2.6215E+00 -3.1916E-01  7.3128E+00  4.7221E-01 -4.2510E+00 -1.6072E+01 -2.7743E-01  4.3240E+00
             2.0786E+00
 GRADIENT:   9.6350E+00  6.3154E-01 -6.6394E+01  2.5308E+02  5.2749E-03 -1.1686E+01  2.9725E-07  0.0000E+00 -3.5597E+01 -2.3571E-04
            -3.7318E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -986.806956053141        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:      553
 NPARAMETR:  5.8960E-01  1.1956E-02  6.4413E-02  6.4879E-01  1.4354E+03  1.4496E+00  1.2201E-02  1.0000E-02  6.7824E-01  7.0204E+01
             7.2223E+00
 PARAMETER: -4.2832E-01 -4.3265E+00 -2.6424E+00 -3.3265E-01  7.3692E+00  4.7129E-01 -4.3062E+00 -1.6195E+01 -2.8826E-01  4.3514E+00
             2.0772E+00
 GRADIENT:  -3.2696E+01  5.8537E-01 -1.2873E+02  2.4647E+02  4.7573E-03 -2.1122E+01  9.6613E-08  0.0000E+00 -3.7156E+01 -4.9342E-04
            -6.2249E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -998.577626645250        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  6.1354E-01  1.0000E-02  5.6544E-02  5.6066E-01  3.3622E+03  1.6001E+00  1.0000E-02  1.0000E-02  7.8639E-01  1.0050E+02
             6.9769E+00
 PARAMETER: -3.8851E-01 -4.8355E+00 -2.7727E+00 -4.7864E-01  8.2204E+00  5.7008E-01 -5.0971E+00 -1.7159E+01 -1.4030E-01  4.7101E+00
             2.0426E+00
 GRADIENT:   4.5376E+01  0.0000E+00 -5.2103E+01  7.1564E+01  1.6995E-03  1.7093E+01  0.0000E+00  0.0000E+00 -2.6623E+00 -4.9569E-05
            -7.5926E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1008.01873177106        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      909            RESET HESSIAN, TYPE II
 NPARAMETR:  4.9112E-01  1.0000E-02  3.3697E-02  3.7958E-01  2.3069E+04  1.4279E+00  1.0000E-02  1.0000E-02  7.3679E-01  2.8509E+02
             7.3725E+00
 PARAMETER: -6.1106E-01 -7.2518E+00 -3.2903E+00 -8.6870E-01  1.0146E+01  4.5622E-01 -7.7023E+00 -2.1208E+01 -2.0546E-01  5.7528E+00
             2.0978E+00
 GRADIENT:   7.6834E+01  0.0000E+00  9.0362E+01  3.9931E+01  7.6383E-05  9.3679E+00  0.0000E+00  0.0000E+00 -1.4002E-01  2.4686E-05
             1.8903E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1008.05548745595        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1090
 NPARAMETR:  4.9043E-01  1.0000E-02  3.3770E-02  3.7903E-01  3.0851E+01  1.4278E+00  1.0000E-02  1.0000E-02  7.3766E-01  8.9051E+00
             7.3768E+00
 PARAMETER: -6.1248E-01 -7.2518E+00 -3.2882E+00 -8.7013E-01  3.5292E+00  4.5612E-01 -7.7023E+00 -2.1208E+01 -2.0428E-01  2.2866E+00
             2.0983E+00
 GRADIENT:  -4.6926E-01  0.0000E+00 -9.6141E-01  1.4268E+00  1.8450E-02 -6.5786E-02  0.0000E+00  0.0000E+00 -1.7479E-02 -2.4991E-04
            -1.4628E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1008.07865859641        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1273
 NPARAMETR:  4.9097E-01  1.0000E-02  3.3769E-02  3.7897E-01  1.2985E+01  1.4282E+00  1.0000E-02  1.0000E-02  7.3740E-01  5.5682E+00
             7.3763E+00
 PARAMETER: -6.1138E-01 -7.2518E+00 -3.2882E+00 -8.7029E-01  2.6638E+00  4.5639E-01 -7.7023E+00 -2.1208E+01 -2.0463E-01  1.8171E+00
             2.0983E+00
 GRADIENT:  -1.4890E+00  0.0000E+00  1.0786E+00 -9.8534E-01  5.8801E-03 -9.0851E-02  0.0000E+00  0.0000E+00  1.1101E-01  1.1973E-03
             1.6693E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1008.11098290759        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:     1421             RESET HESSIAN, TYPE I
 NPARAMETR:  4.9027E-01  1.0000E-02  3.3477E-02  3.7603E-01  1.2165E+01  1.4281E+00  1.0000E-02  1.0000E-02  7.3888E-01  3.7064E+00
             7.3769E+00
 PARAMETER: -6.1279E-01 -7.2518E+00 -3.2969E+00 -8.7808E-01  2.5986E+00  4.5634E-01 -7.7023E+00 -2.1208E+01 -2.0262E-01  1.4101E+00
             2.0983E+00
 GRADIENT:   7.3713E+01  0.0000E+00  1.0318E+02  2.5554E+01  1.7291E-03  9.5263E+00  0.0000E+00  0.0000E+00  9.1088E-01  4.3741E-03
             2.1415E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1008.12034445246        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1559
 NPARAMETR:  4.8774E-01  1.0000E-02  3.3343E-02  3.7558E-01  1.2930E+01  1.4262E+00  1.0000E-02  1.0000E-02  7.3792E-01  6.9668E+00
             7.3688E+00
 PARAMETER: -6.1798E-01 -7.2518E+00 -3.3009E+00 -8.7928E-01  2.6595E+00  4.5500E-01 -7.7023E+00 -2.1208E+01 -2.0391E-01  2.0412E+00
             2.0973E+00
 GRADIENT:  -3.3860E+00  0.0000E+00  7.4537E-01  3.0096E-02 -1.8965E-02 -1.7508E-01  0.0000E+00  0.0000E+00  1.9694E-01  7.2556E-03
            -4.9678E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1008.14380435374        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1732
 NPARAMETR:  4.8764E-01  1.0000E-02  3.3120E-02  3.7469E-01  1.1816E+01  1.4251E+00  1.0000E-02  1.0000E-02  7.3689E-01  3.7159E+00
             7.3707E+00
 PARAMETER: -6.1817E-01 -7.2518E+00 -3.3076E+00 -8.8165E-01  2.5694E+00  4.5425E-01 -7.7023E+00 -2.1208E+01 -2.0532E-01  1.4126E+00
             2.0975E+00
 GRADIENT:   7.2854E+01  0.0000E+00  9.6968E+01  3.6109E+01 -8.8587E-03  8.9851E+00  0.0000E+00  0.0000E+00  3.0415E-01  5.3359E-03
             2.0018E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1008.14634418222        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:     1815
 NPARAMETR:  4.8442E-01  1.0000E-02  3.2850E-02  3.7292E-01  1.2151E+01  1.4226E+00  1.0000E-02  1.0000E-02  7.3740E-01  2.8814E+00
             7.3607E+00
 PARAMETER: -6.2480E-01 -7.2518E+00 -3.3158E+00 -8.8639E-01  2.5974E+00  4.5245E-01 -7.7023E+00 -2.1208E+01 -2.0462E-01  1.1583E+00
             2.0962E+00
 GRADIENT:   7.0374E+01  0.0000E+00  9.6339E+01  3.9735E+01 -3.6749E-04  8.5636E+00  0.0000E+00  0.0000E+00  2.6449E-01  2.0734E-03
             1.9057E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1008.19362420665        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:     1947
 NPARAMETR:  4.8352E-01  1.0000E-02  3.2613E-02  3.6997E-01  1.3027E+01  1.4219E+00  1.0000E-02  1.0000E-02  7.3818E-01  2.3449E+00
             7.3601E+00
 PARAMETER: -6.2666E-01 -7.2518E+00 -3.3230E+00 -8.9434E-01  2.6671E+00  4.5202E-01 -7.7023E+00 -2.1208E+01 -2.0356E-01  9.5224E-01
             2.0961E+00
 GRADIENT:   7.1957E+01  0.0000E+00  1.0061E+02  3.4122E+01  1.6211E-02  8.7185E+00  0.0000E+00  0.0000E+00  5.6155E-01  7.9327E-04
             1.9280E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1008.21193594601        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     2143             RESET HESSIAN, TYPE I
 NPARAMETR:  4.8593E-01  1.0000E-02  3.2593E-02  3.6947E-01  1.1931E+01  1.4247E+00  1.0000E-02  1.0000E-02  7.3673E-01  2.2504E+00
             7.3775E+00
 PARAMETER: -6.2170E-01 -7.2518E+00 -3.3237E+00 -8.9569E-01  2.5792E+00  4.5400E-01 -7.7023E+00 -2.1208E+01 -2.0554E-01  9.1109E-01
             2.0984E+00
 GRADIENT:   7.5805E+01  0.0000E+00  1.0078E+02  3.1276E+01 -6.0415E-03  9.4131E+00  0.0000E+00  0.0000E+00  5.0432E-01  1.2278E-03
             2.1022E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1008.21538006982        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     2219
 NPARAMETR:  4.8185E-01  1.0000E-02  3.2294E-02  3.6801E-01  1.2087E+01  1.4209E+00  1.0000E-02  1.0000E-02  7.3697E-01  9.8665E-01
             7.3632E+00
 PARAMETER: -6.3012E-01 -7.2518E+00 -3.3329E+00 -8.9966E-01  2.5921E+00  4.5131E-01 -7.7023E+00 -2.1208E+01 -2.0520E-01  8.6555E-02
             2.0965E+00
 GRADIENT:   7.2113E+01  0.0000E+00  9.8097E+01  3.8646E+01 -8.3024E-03  8.6783E+00  0.0000E+00  0.0000E+00  3.0368E-01  4.7208E-05
             1.9426E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1008.21563489780        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:     2338
 NPARAMETR:  4.8110E-01  1.0000E-02  3.2223E-02  3.6751E-01  1.2157E+01  1.4204E+00  1.0000E-02  1.0000E-02  7.3708E-01  8.3926E-01
             7.3611E+00
 PARAMETER: -6.3167E-01 -7.2518E+00 -3.3351E+00 -9.0100E-01  2.5979E+00  4.5092E-01 -7.7023E+00 -2.1208E+01 -2.0507E-01 -7.5238E-02
             2.0962E+00
 GRADIENT:  -4.7955E+00  0.0000E+00 -5.8097E+00  8.6409E+00 -2.1593E-02 -6.8781E-01  0.0000E+00  0.0000E+00 -9.8761E-02  1.6616E-05
            -1.7694E+00

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1008.24976483624        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2535
 NPARAMETR:  4.8398E-01  1.0000E-02  3.2180E-02  3.6652E-01  1.1890E+01  1.4232E+00  1.0000E-02  1.0000E-02  7.3535E-01  8.5880E-01
             7.3761E+00
 PARAMETER: -6.2572E-01 -7.2518E+00 -3.3364E+00 -9.0370E-01  2.5757E+00  4.5291E-01 -7.7023E+00 -2.1208E+01 -2.0741E-01 -5.2214E-02
             2.0983E+00
 GRADIENT:   1.4381E+00  0.0000E+00 -4.7814E+00  4.0585E+00 -2.0393E-02  6.5174E-02  0.0000E+00  0.0000E+00 -2.0226E-01  4.4671E-05
            -5.3513E-01

0ITERATION NO.:  102    OBJECTIVE VALUE:  -1008.25116259749        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     2592
 NPARAMETR:  4.8382E-01  1.0000E-02  3.2205E-02  3.6637E-01  1.2223E+01  1.4231E+00  1.0000E-02  1.0000E-02  7.3583E-01  8.5858E-01
             7.3770E+00
 PARAMETER: -6.2605E-01 -7.2518E+00 -3.3356E+00 -9.0411E-01  2.6033E+00  4.5287E-01 -7.7023E+00 -2.1208E+01 -2.0675E-01 -5.2474E-02
             2.0984E+00
 GRADIENT:   9.2027E-01  0.0000E+00 -2.8843E+00  1.7954E+00 -6.5463E-03  6.4567E-02  0.0000E+00  0.0000E+00 -5.7280E-02  1.9704E-05
            -1.4005E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2592
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2788E-03  7.2177E-07  9.5766E-05 -1.9239E-02 -2.0800E-05
 SE:             2.9182E-02  1.0982E-06  2.6620E-04  2.3853E-02  8.8348E-05
 N:                     100         100         100         100         100

 P VAL.:         9.3776E-01  5.1102E-01  7.1903E-01  4.1990E-01  8.1388E-01

 ETASHRINKSD(%)  2.2371E+00  9.9996E+01  9.9108E+01  2.0091E+01  9.9704E+01
 ETASHRINKVR(%)  4.4242E+00  1.0000E+02  9.9992E+01  3.6145E+01  9.9999E+01
 EBVSHRINKSD(%)  2.3631E+00  9.9996E+01  9.9158E+01  2.1005E+01  9.9709E+01
 EBVSHRINKVR(%)  4.6703E+00  1.0000E+02  9.9993E+01  3.7598E+01  9.9999E+01
 RELATIVEINF(%)  3.0959E+00  1.1767E-08  5.7597E-05  4.3616E-01  4.5014E-05
 EPSSHRINKSD(%)  1.2810E+01
 EPSSHRINKVR(%)  2.3980E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1008.2511625974910     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -89.312629392818280     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    49.54
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.93
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1008.251       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.84E-01  1.00E-02  3.22E-02  3.66E-01  1.22E+01  1.42E+00  1.00E-02  1.00E-02  7.36E-01  8.59E-01  7.38E+00
 


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
+        2.20E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.16E+04  0.00E+00  7.04E+05
 
 TH 4
+       -5.94E+01  0.00E+00 -7.49E+04  8.93E+03
 
 TH 5
+        3.42E-01  0.00E+00 -5.36E+00  4.47E-01  2.13E-03
 
 TH 6
+       -5.77E+00  0.00E+00 -6.54E+01 -2.71E+01  6.46E-03  8.26E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.33E+00  0.00E+00  7.31E+02 -9.32E+01 -1.69E-02 -1.37E+00  0.00E+00  0.00E+00  1.35E+02
 
 TH10
+       -3.11E-02  0.00E+00  5.52E-03 -8.10E-03 -2.54E-05  5.16E-03  0.00E+00  0.00E+00  9.41E-03 -1.45E-02
 
 TH11
+       -2.38E+01  0.00E+00  3.52E+02 -2.80E+01 -5.45E-03  3.65E-01  0.00E+00  0.00E+00  7.00E+00  6.76E-05  9.87E+00
 
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
 #CPUT: Total CPU Time in Seconds,       59.529
Stop Time:
Thu Sep 30 02:53:37 CDT 2021
