Sat Sep 25 13:17:55 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat7.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m7.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1643.36852641838        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.5997E+01  4.6832E+00 -3.5049E+01  5.7728E+01  5.2906E+01 -3.9153E+01  1.0480E+01  6.8545E+00  1.2798E+01  1.4411E+00
            -3.1387E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1649.56456627412        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:       89
 NPARAMETR:  9.7223E-01  1.0159E+00  1.0419E+00  9.7889E-01  9.9327E-01  1.0789E+00  9.2963E-01  9.6249E-01  9.5010E-01  9.6893E-01
             1.0779E+00
 PARAMETER:  7.1839E-02  1.1582E-01  1.4104E-01  7.8662E-02  9.3248E-02  1.7598E-01  2.7033E-02  6.1764E-02  4.8814E-02  6.8441E-02
             1.7501E-01
 GRADIENT:   2.0203E+01  1.3227E+01 -6.4298E+00  2.7667E+01  1.2563E+01  4.2341E-02  5.0275E+00  2.7889E+00 -1.4541E+00 -2.0046E+00
            -1.7995E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1650.89534385433        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.8468E-01  1.0822E+00  1.0465E+00  9.2705E-01  1.0271E+00  1.0830E+00  7.5400E-01  8.5452E-01  1.0532E+00  1.0492E+00
             1.0715E+00
 PARAMETER:  8.4557E-02  1.7898E-01  1.4546E-01  2.4252E-02  1.2674E-01  1.7974E-01 -1.8236E-01 -5.7219E-02  1.5183E-01  1.4804E-01
             1.6904E-01
 GRADIENT:   4.4975E+01  8.2105E+00  2.3286E-01  1.2324E+01  5.2846E+00  1.7907E+00  4.1102E+00  1.0732E-01  4.2243E+00  7.3540E-01
            -3.5037E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1651.25959867238        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.7094E-01  1.1135E+00  1.0668E+00  9.0482E-01  1.0440E+00  1.0815E+00  6.0339E-01  9.8250E-01  1.0968E+00  1.0700E+00
             1.0748E+00
 PARAMETER:  7.0514E-02  2.0753E-01  1.6466E-01 -1.4064E-05  1.4310E-01  1.7836E-01 -4.0519E-01  8.2347E-02  1.9236E-01  1.6766E-01
             1.7217E-01
 GRADIENT:   1.8201E+01  1.1794E+01 -9.3953E-01  1.1300E+01 -6.0127E-02  1.2207E+00  2.1211E+00  7.6133E-01  1.3287E+00  1.0833E+00
            -1.6842E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1651.27017775807        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      306
 NPARAMETR:  9.6660E-01  1.0887E+00  1.1209E+00  9.1722E-01  1.0544E+00  1.0803E+00  5.1444E-01  1.0399E+00  1.1069E+00  1.0912E+00
             1.0771E+00
 PARAMETER:  6.6030E-02  1.8496E-01  2.1414E-01  1.3592E-02  1.5296E-01  1.7726E-01 -5.6468E-01  1.3913E-01  2.0152E-01  1.8728E-01
             1.7431E-01
 GRADIENT:   1.0248E+01  7.0477E+00 -6.1482E-01  6.7741E+00 -6.5964E-02  7.2529E-01  1.7845E+00  5.9460E-01  1.5491E+00  1.2586E+00
            -8.6010E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1651.83008600463        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  9.7789E-01  1.1360E+00  1.1025E+00  8.7959E-01  1.0727E+00  1.1178E+00  4.3881E-01  9.9285E-01  1.1607E+00  1.1122E+00
             1.0888E+00
 PARAMETER:  7.7642E-02  2.2752E-01  1.9756E-01 -2.8297E-02  1.7014E-01  2.1140E-01 -7.2368E-01  9.2828E-02  2.4898E-01  2.0631E-01
             1.8509E-01
 GRADIENT:  -9.0034E+00 -9.1935E+00 -5.5529E-01 -6.5839E+00  1.9515E+00  3.0385E+00  6.4428E-01 -3.8380E-01  6.4870E-01  4.7926E-01
             2.0538E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1652.10300513505        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      621
 NPARAMETR:  9.8380E-01  1.2278E+00  1.1097E+00  8.2154E-01  1.1117E+00  1.1095E+00  3.4275E-01  1.2025E+00  1.2576E+00  1.1213E+00
             1.0842E+00
 PARAMETER:  8.3667E-02  3.0524E-01  2.0413E-01 -9.6569E-02  2.0590E-01  2.0395E-01 -9.7075E-01  2.8443E-01  3.2919E-01  2.1453E-01
             1.8088E-01
 GRADIENT:   9.4177E-01  4.2761E-01 -1.9777E-01  9.5067E-01  2.9396E-01 -1.9043E-01 -6.4878E-02  5.8662E-02 -1.7007E-01 -9.8389E-02
            -1.4477E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1652.10394689811        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      807            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8366E-01  1.2284E+00  1.1091E+00  8.2103E-01  1.1115E+00  1.1102E+00  3.4458E-01  1.1999E+00  1.2581E+00  1.1214E+00
             1.0844E+00
 PARAMETER:  8.3524E-02  3.0570E-01  2.0352E-01 -9.7198E-02  2.0570E-01  2.0452E-01 -9.6543E-01  2.8224E-01  3.2957E-01  2.1459E-01
             1.8101E-01
 GRADIENT:   4.1798E+01  2.0144E+01  1.1984E-01  5.0445E+00  1.1065E+00  1.2964E+01  6.7351E-01  5.9829E-02  2.1253E+00  7.7973E-02
             4.2768E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1652.10433822751        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:      974
 NPARAMETR:  9.8368E-01  1.2285E+00  1.1091E+00  8.2105E-01  1.1114E+00  1.1101E+00  3.4907E-01  1.1955E+00  1.2574E+00  1.1216E+00
             1.0846E+00
 PARAMETER:  8.3549E-02  3.0577E-01  2.0357E-01 -9.7173E-02  2.0565E-01  2.0446E-01 -9.5250E-01  2.7858E-01  3.2909E-01  2.1476E-01
             1.8124E-01
 GRADIENT:   7.3756E-01  1.0880E-01  1.7804E-01  5.9468E-01 -1.9932E-01  2.3563E-02 -8.5256E-03 -1.5606E-02  1.3674E-01 -4.7155E-02
            -4.2251E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1652.10469883152        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1134
 NPARAMETR:  9.8337E-01  1.2284E+00  1.1086E+00  8.2089E-01  1.1115E+00  1.1099E+00  3.5188E-01  1.1951E+00  1.2558E+00  1.1216E+00
             1.0847E+00
 PARAMETER:  8.3234E-02  3.0571E-01  2.0312E-01 -9.7364E-02  2.0574E-01  2.0424E-01 -9.4447E-01  2.7824E-01  3.2781E-01  2.1477E-01
             1.8130E-01
 GRADIENT:   1.7611E-01 -4.4108E-01  5.5330E-02  1.3220E-01  1.0219E-01 -5.8904E-02 -2.5543E-05  1.8953E-02 -2.5341E-02 -7.0484E-04
             1.9292E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1652.10802874906        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1317
 NPARAMETR:  9.8310E-01  1.2364E+00  1.0963E+00  8.1610E-01  1.1107E+00  1.1121E+00  3.6276E-01  1.1774E+00  1.2599E+00  1.1216E+00
             1.0846E+00
 PARAMETER:  8.2955E-02  3.1220E-01  1.9192E-01 -1.0321E-01  2.0501E-01  2.0628E-01 -9.1402E-01  2.6335E-01  3.3105E-01  2.1474E-01
             1.8119E-01
 GRADIENT:  -4.4153E-01  2.2443E-01  1.1320E-01  7.6926E-01 -3.3986E-01  7.2256E-01 -1.8631E-02 -5.0564E-02  1.2382E-01  7.7750E-02
            -1.0306E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1652.17123668347        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1499
 NPARAMETR:  9.8487E-01  1.3576E+00  9.3460E-01  7.4283E-01  1.1060E+00  1.1111E+00  5.4048E-01  9.9491E-01  1.2993E+00  1.1005E+00
             1.0860E+00
 PARAMETER:  8.4757E-02  4.0569E-01  3.2368E-02 -1.9729E-01  2.0079E-01  2.0539E-01 -5.1530E-01  9.4898E-02  3.6182E-01  1.9578E-01
             1.8250E-01
 GRADIENT:   1.1964E+00  4.5620E+00  1.2929E+00  5.0445E+00 -4.4394E+00  1.2174E-01  8.8843E-01  9.2424E-02  2.0742E+00  1.8500E+00
             3.4763E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1652.47534622694        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1677
 NPARAMETR:  9.8386E-01  1.5225E+00  7.2449E-01  6.2996E-01  1.1099E+00  1.1113E+00  6.1198E-01  7.5964E-01  1.3663E+00  1.0408E+00
             1.0888E+00
 PARAMETER:  8.3730E-02  5.2038E-01 -2.2229E-01 -3.6210E-01  2.0429E-01  2.0553E-01 -3.9106E-01 -1.7491E-01  4.1210E-01  1.4003E-01
             1.8512E-01
 GRADIENT:  -2.4492E+00 -6.7640E-01 -3.5168E-01  3.2345E-01  2.8532E+00 -2.6000E-01 -6.0765E-01  1.2089E-01 -2.8273E+00 -9.8459E-01
             4.3089E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1652.69801341304        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1852
 NPARAMETR:  9.8614E-01  1.6806E+00  5.5741E-01  5.2216E-01  1.1194E+00  1.1122E+00  6.1567E-01  5.3189E-01  1.5389E+00  1.0167E+00
             1.0874E+00
 PARAMETER:  8.6040E-02  6.1916E-01 -4.8445E-01 -5.4979E-01  2.1283E-01  2.0638E-01 -3.8504E-01 -5.3131E-01  5.3110E-01  1.1661E-01
             1.8376E-01
 GRADIENT:   7.8167E-01 -7.1877E-03 -7.1425E-01  1.1406E+00 -2.0711E+00 -2.5480E-01 -2.5955E-01  1.5998E-01  5.5812E-01  8.0500E-02
             4.3096E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1652.74168556540        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2027
 NPARAMETR:  9.8570E-01  1.7465E+00  5.1777E-01  4.7672E-01  1.1457E+00  1.1126E+00  6.0954E-01  3.9559E-01  1.6284E+00  1.0286E+00
             1.0905E+00
 PARAMETER:  8.5597E-02  6.5760E-01 -5.5822E-01 -6.4083E-01  2.3600E-01  2.0668E-01 -3.9506E-01 -8.2738E-01  5.8763E-01  1.2820E-01
             1.8663E-01
 GRADIENT:  -3.5745E-02 -1.4858E+00 -1.2763E-01 -6.0145E-01  6.1716E-01 -1.3792E-01  9.2121E-03  3.7189E-02 -9.6555E-02 -4.3209E-03
             2.6721E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1652.74340138611        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2204
 NPARAMETR:  9.8570E-01  1.7399E+00  5.2450E-01  4.8186E-01  1.1444E+00  1.1127E+00  6.0978E-01  4.0107E-01  1.6198E+00  1.0298E+00
             1.0906E+00
 PARAMETER:  8.5598E-02  6.5382E-01 -5.4531E-01 -6.3011E-01  2.3489E-01  2.0683E-01 -3.9466E-01 -8.1363E-01  5.8232E-01  1.2934E-01
             1.8677E-01
 GRADIENT:  -2.6428E-02  7.5915E-02  1.9525E-02 -7.1232E-02  1.9751E-01 -7.2801E-02  4.3719E-03  3.0936E-02 -3.7648E-02  1.3793E-02
             1.8211E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1652.75093375422        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2383
 NPARAMETR:  9.8579E-01  1.7419E+00  5.1053E-01  4.8008E-01  1.1381E+00  1.1130E+00  6.1427E-01  2.1574E-01  1.6147E+00  1.0234E+00
             1.0908E+00
 PARAMETER:  8.5691E-02  6.5498E-01 -5.7230E-01 -6.3381E-01  2.2939E-01  2.0703E-01 -3.8732E-01 -1.4337E+00  5.7915E-01  1.2318E-01
             1.8692E-01
 GRADIENT:   7.9085E-02 -3.1435E-02 -9.0352E-03 -1.7715E-02 -1.2573E-01 -3.7148E-03 -7.9866E-03  2.8704E-03  4.1538E-03  3.5139E-02
            -2.4571E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1652.75124310156        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2558
 NPARAMETR:  9.8576E-01  1.7458E+00  5.0553E-01  4.7742E-01  1.1382E+00  1.1130E+00  6.1487E-01  1.2573E-01  1.6185E+00  1.0225E+00
             1.0911E+00
 PARAMETER:  8.5654E-02  6.5719E-01 -5.8215E-01 -6.3936E-01  2.2944E-01  2.0704E-01 -3.8635E-01 -1.9736E+00  5.8149E-01  1.2227E-01
             1.8723E-01
 GRADIENT:  -3.5502E-04 -2.1214E-02 -1.4909E-02  2.4841E-03  2.3569E-02 -3.6372E-03  1.4466E-03  2.9312E-04 -6.1658E-03  9.2621E-03
             8.1678E-03

0ITERATION NO.:   89    OBJECTIVE VALUE:  -1652.75126077241        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     2696
 NPARAMETR:  9.8576E-01  1.7460E+00  5.0495E-01  4.7724E-01  1.1380E+00  1.1130E+00  6.1501E-01  1.0071E-01  1.6186E+00  1.0223E+00
             1.0912E+00
 PARAMETER:  8.5654E-02  6.5732E-01 -5.8329E-01 -6.3974E-01  2.2929E-01  2.0705E-01 -3.8612E-01 -2.1955E+00  5.8156E-01  1.2206E-01
             1.8726E-01
 GRADIENT:  -1.9534E-02  1.1325E-01  2.5036E-03 -4.3476E-03 -2.9464E-03 -9.2251E-03 -7.0209E-04  5.4888E-05 -4.9775E-03 -1.7979E-03
            -2.6152E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2696
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.1643E-04 -4.0520E-02 -2.2977E-03  2.8486E-02 -3.7439E-02
 SE:             2.9861E-02  2.1267E-02  9.3494E-04  2.3425E-02  2.3271E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9155E-01  5.6736E-02  1.3985E-02  2.2395E-01  1.0765E-01

 ETASHRINKSD(%)  1.0000E-10  2.8754E+01  9.6868E+01  2.1524E+01  2.2040E+01
 ETASHRINKVR(%)  1.0000E-10  4.9240E+01  9.9902E+01  3.8415E+01  3.9223E+01
 EBVSHRINKSD(%)  4.1227E-01  2.8071E+01  9.7380E+01  2.2816E+01  2.0652E+01
 EBVSHRINKVR(%)  8.2285E-01  4.8263E+01  9.9931E+01  4.0427E+01  3.7038E+01
 RELATIVEINF(%)  9.9137E+01  3.5816E+00  9.4868E-03  4.6881E+00  1.5664E+01
 EPSSHRINKSD(%)  4.3597E+01
 EPSSHRINKVR(%)  6.8188E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1652.7512607724091     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -917.60043420867089     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    33.67
 Elapsed covariance  time in seconds:     6.39
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1652.751       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.75E+00  5.05E-01  4.77E-01  1.14E+00  1.11E+00  6.15E-01  1.01E-01  1.62E+00  1.02E+00  1.09E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.30E-02  2.75E-01  3.22E-01  1.96E-01  1.23E-01  7.25E-02  1.06E-01  1.47E-01  3.81E-01  1.80E-01  1.01E-01
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.09E-03
 
 TH 2
+        1.23E-03  7.56E-02
 
 TH 3
+       -1.03E-03 -7.10E-02  1.03E-01
 
 TH 4
+       -7.92E-04 -5.31E-02  5.14E-02  3.84E-02
 
 TH 5
+        5.00E-04  9.46E-03  1.20E-02 -6.19E-03  1.52E-02
 
 TH 6
+       -3.49E-04  2.40E-03 -6.14E-05 -1.98E-03  2.41E-03  5.26E-03
 
 TH 7
+        4.69E-04 -2.63E-05 -1.30E-02  1.73E-04 -7.29E-03 -6.92E-04  1.11E-02
 
 TH 8
+       -1.00E-03 -2.45E-02  1.96E-02  1.77E-02 -4.48E-03  2.31E-03  2.24E-03  2.15E-02
 
 TH 9
+        1.24E-03  8.04E-02 -4.15E-02 -5.67E-02  2.93E-02  3.99E-03 -1.55E-02 -3.16E-02  1.45E-01
 
 TH10
+       -5.19E-04 -1.70E-02  4.08E-02  1.31E-02  1.28E-02  2.05E-03 -9.76E-03  8.73E-03  9.00E-04  3.25E-02
 
 TH11
+       -2.95E-04  1.56E-03  7.83E-03 -1.13E-03  6.12E-03  1.04E-03 -2.17E-03  2.22E-03  1.01E-02  5.46E-03  1.02E-02
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        3.30E-02
 
 TH 2
+        1.36E-01  2.75E-01
 
 TH 3
+       -9.70E-02 -8.03E-01  3.22E-01
 
 TH 4
+       -1.23E-01 -9.86E-01  8.16E-01  1.96E-01
 
 TH 5
+        1.23E-01  2.79E-01  3.01E-01 -2.56E-01  1.23E-01
 
 TH 6
+       -1.46E-01  1.21E-01 -2.63E-03 -1.40E-01  2.69E-01  7.25E-02
 
 TH 7
+        1.35E-01 -9.06E-04 -3.84E-01  8.36E-03 -5.59E-01 -9.04E-02  1.06E-01
 
 TH 8
+       -2.07E-01 -6.07E-01  4.15E-01  6.14E-01 -2.47E-01  2.17E-01  1.45E-01  1.47E-01
 
 TH 9
+        9.85E-02  7.67E-01 -3.38E-01 -7.60E-01  6.22E-01  1.44E-01 -3.84E-01 -5.65E-01  3.81E-01
 
 TH10
+       -8.74E-02 -3.43E-01  7.04E-01  3.71E-01  5.75E-01  1.57E-01 -5.13E-01  3.30E-01  1.31E-02  1.80E-01
 
 TH11
+       -8.89E-02  5.61E-02  2.41E-01 -5.72E-02  4.92E-01  1.42E-01 -2.04E-01  1.50E-01  2.64E-01  3.00E-01  1.01E-01
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.12E+03
 
 TH 2
+        2.61E+01  6.35E+02
 
 TH 3
+        4.47E+01  1.86E+02  3.82E+02
 
 TH 4
+       -4.05E+01  5.40E+02 -4.04E+02  1.53E+03
 
 TH 5
+       -2.13E+02 -2.05E+02 -3.00E+02  2.19E+02  5.88E+02
 
 TH 6
+        1.03E+02  5.58E+01 -2.89E+01  1.69E+02 -9.52E+01  2.87E+02
 
 TH 7
+       -9.08E+01 -4.24E+01  5.97E+01 -1.25E+02 -1.61E+01 -1.29E+01  1.94E+02
 
 TH 8
+        1.12E+01  1.86E+01  9.15E+01 -1.63E+02  1.67E+01 -9.45E+01 -7.70E+00  1.55E+02
 
 TH 9
+        8.57E-01 -4.85E+01 -5.78E+01  7.54E+01  4.96E+00  1.16E+01  1.71E+01 -9.80E+00  4.44E+01
 
 TH10
+        2.71E+01 -5.85E+01 -9.09E+01  5.38E+01 -4.82E+01  3.05E+01  2.41E+01 -7.04E+01  1.95E+01  1.29E+02
 
 TH11
+        7.03E+01  4.82E+00 -8.35E+01  1.53E+02 -4.96E+01  5.34E+01 -3.20E+01 -7.94E+01  7.14E+00  4.36E+01  1.85E+02
 
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
 #CPUT: Total CPU Time in Seconds,       39.631
Stop Time:
Sat Sep 25 13:18:37 CDT 2021
