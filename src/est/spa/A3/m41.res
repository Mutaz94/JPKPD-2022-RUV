Wed Sep 29 13:31:38 CDT 2021
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
$DATA ../../../../data/spa/A3/dat41.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m41.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -22.3898243517026        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9804E+02  3.8887E+01  5.9998E+01 -2.5756E+01  2.0357E+02  5.0120E+01 -5.7563E+01 -2.2466E+01 -1.2049E+02 -1.4551E+02
            -2.9118E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1086.52597736044        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0174E+00  6.0058E-01  5.0814E-01  1.2386E+00  4.7130E-01  9.9386E-01  9.9462E-01  6.6043E-01  1.1154E+00  8.6596E-01
             1.8035E+00
 PARAMETER:  1.1720E-01 -4.0986E-01 -5.7699E-01  3.1396E-01 -6.5227E-01  9.3843E-02  9.4604E-02 -3.1487E-01  2.0923E-01 -4.3918E-02
             6.8971E-01
 GRADIENT:   1.3582E+02  7.8098E+01  8.6310E+01  1.0699E+02  2.0501E+01  2.7302E+01 -8.9098E+00  7.1376E+00 -1.7337E+01  1.5409E+01
            -6.8796E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1166.35410832498        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.0113E+00  6.1170E-01  1.3846E-01  1.2205E+00  2.3977E-01  1.4560E+00  7.2442E-01  8.5272E-01  1.7437E+00  1.6317E-01
             1.8333E+00
 PARAMETER:  1.1120E-01 -3.9151E-01 -1.8772E+00  2.9928E-01 -1.3281E+00  4.7566E-01 -2.2239E-01 -5.9325E-02  6.5599E-01 -1.7130E+00
             7.0612E-01
 GRADIENT:   1.1924E+02  4.0551E+02  2.2389E+02  1.6677E+02 -2.4989E+02  1.3744E+02 -6.0569E+01 -7.8364E+01 -2.0330E+01 -9.5094E+00
            -3.0777E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1204.86570687561        NO. OF FUNC. EVALS.: 116
 CUMULATIVE NO. OF FUNC. EVALS.:      277
 NPARAMETR:  1.0124E+00  6.3250E-01  1.5039E-01  1.2819E+00  2.6125E-01  1.4862E+00  8.4594E-01  1.2406E+00  2.1369E+00  7.0175E-02
             1.9635E+00
 PARAMETER:  1.1237E-01 -3.5807E-01 -1.7945E+00  3.4833E-01 -1.2423E+00  4.9623E-01 -6.7302E-02  3.1559E-01  8.5936E-01 -2.5568E+00
             7.7471E-01
 GRADIENT:   3.3734E+01  1.7621E+02  7.9678E+01  8.8965E+01 -7.5926E+01  1.1495E+02 -2.4958E+00 -5.1654E+01  4.1182E+00 -1.3862E-01
            -2.0347E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1309.63190389989        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      452
 NPARAMETR:  1.0049E+00  4.5335E-01  1.2914E-01  9.8904E-01  2.2013E-01  9.2743E-01  5.8246E-01  1.7253E+00  1.6649E+00  5.9366E-02
             2.5764E+00
 PARAMETER:  1.0490E-01 -6.9110E-01 -1.9469E+00  8.8976E-02 -1.4136E+00  2.4665E-02 -4.4049E-01  6.4538E-01  6.0976E-01 -2.7240E+00
             1.0464E+00
 GRADIENT:   3.6384E+00  7.7851E+00  5.4175E+00 -9.9911E+00 -3.2179E+01 -4.1316E+00 -1.0492E+00  3.9805E+00 -6.8642E+00 -1.9908E-01
            -2.6996E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1310.66123702196        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      628            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0039E+00  4.6670E-01  1.2846E-01  1.0043E+00  2.2384E-01  9.3710E-01  5.9010E-01  1.6782E+00  1.6870E+00  5.8450E-02
             2.7051E+00
 PARAMETER:  1.0394E-01 -6.6206E-01 -1.9521E+00  1.0427E-01 -1.3968E+00  3.5033E-02 -4.2746E-01  6.1773E-01  6.2297E-01 -2.7396E+00
             1.0951E+00
 GRADIENT:   3.8416E+01  1.3116E+01  3.1104E+01  5.8504E+00  1.3631E+02  2.7925E+00  5.2797E-01  4.1905E+00  1.1969E+01  2.8073E-02
             8.3739E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1310.72440908820        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      765
 NPARAMETR:  1.0061E+00  4.6579E-01  1.2823E-01  1.0055E+00  2.2301E-01  9.3401E-01  5.6850E-01  1.6835E+00  1.6943E+00  1.1524E-01
             2.7007E+00
 PARAMETER:  1.0608E-01 -6.6401E-01 -1.9540E+00  1.0548E-01 -1.4005E+00  3.1731E-02 -4.6476E-01  6.2086E-01  6.2726E-01 -2.0608E+00
             1.0935E+00
 GRADIENT:   5.1993E+00  3.6011E+00  1.7639E+00  1.0114E+00 -3.4131E+00 -1.0421E+00  3.0784E-01  1.2418E+00  7.2969E-01 -1.1317E-01
            -8.9244E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1310.88971733293        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      942
 NPARAMETR:  1.0026E+00  4.5727E-01  1.2534E-01  9.9967E-01  2.1979E-01  9.3949E-01  4.9056E-01  1.6510E+00  1.7019E+00  2.0551E-01
             2.7135E+00
 PARAMETER:  1.0263E-01 -6.8249E-01 -1.9768E+00  9.9670E-02 -1.4151E+00  3.7583E-02 -6.1220E-01  6.0139E-01  6.3174E-01 -1.4823E+00
             1.0982E+00
 GRADIENT:  -2.8995E+00 -2.0607E+00 -9.5752E-01 -7.6273E-01  9.8378E+00  8.3973E-01  2.6752E-01 -3.6485E-01  5.8045E-01  5.4544E-02
             2.2110E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1311.17737895115        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1122
 NPARAMETR:  1.0038E+00  4.2099E-01  1.1514E-01  9.8299E-01  2.0430E-01  9.4041E-01  2.4339E-01  1.5781E+00  1.7663E+00  3.5982E-01
             2.6943E+00
 PARAMETER:  1.0378E-01 -7.6515E-01 -2.0616E+00  8.2843E-02 -1.4882E+00  3.8563E-02 -1.3131E+00  5.5623E-01  6.6890E-01 -9.2215E-01
             1.0911E+00
 GRADIENT:   3.6286E-01 -3.8984E-01 -4.6853E-01 -3.0337E-01  3.0838E-01 -8.6967E-02  1.7720E-01  3.3487E-01  4.2710E-01  1.2381E-01
            -1.6285E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1311.20549826903        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1300
 NPARAMETR:  1.0046E+00  4.1598E-01  1.1516E-01  9.8382E-01  2.0307E-01  9.4056E-01  1.4366E-01  1.5739E+00  1.7616E+00  3.7393E-01
             2.6994E+00
 PARAMETER:  1.0454E-01 -7.7711E-01 -2.0614E+00  8.3689E-02 -1.4942E+00  3.8720E-02 -1.8403E+00  5.5355E-01  6.6624E-01 -8.8368E-01
             1.0930E+00
 GRADIENT:   6.7665E-01 -1.9513E+00  1.0658E+00  4.9111E-01  5.1776E-01  1.9083E-03  5.1198E-02  1.2528E-01 -2.4903E-01  2.0521E-01
             8.4722E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1313.82456465006        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1460
 NPARAMETR:  1.0010E+00  4.1878E-01  1.1365E-01  9.7723E-01  2.0647E-01  9.4010E-01  4.2296E-02  1.5349E+00  1.7468E+00  3.8336E-01
             2.6937E+00
 PARAMETER:  1.0101E-01 -7.7042E-01 -2.0747E+00  7.6968E-02 -1.4776E+00  3.8226E-02 -3.0631E+00  5.2849E-01  6.5778E-01 -8.5878E-01
             1.0909E+00
 GRADIENT:   2.8522E+01 -1.1719E+01  6.7784E+00  1.1494E+01  2.0333E+02  1.0067E+00  9.0587E-03  1.7412E+01  1.5795E+01 -1.3102E+00
             1.0489E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1319.11279087786        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1597
 NPARAMETR:  1.0249E+00  4.0926E-01  1.1537E-01  9.6895E-01  2.0203E-01  9.5489E-01  1.0000E-02  1.1234E+00  1.7726E+00  4.6926E-01
             2.7341E+00
 PARAMETER:  1.2461E-01 -7.9340E-01 -2.0596E+00  6.8461E-02 -1.4993E+00  5.3839E-02 -5.0266E+00  2.1637E-01  6.7242E-01 -6.5659E-01
             1.1058E+00
 GRADIENT:   3.8695E+01 -1.5463E+01 -1.1295E+01  7.8762E+00  3.0069E+01 -8.5175E-01  0.0000E+00  6.2867E+00  1.3866E+01 -3.5427E+00
            -1.1718E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1320.91714200740        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1774
 NPARAMETR:  1.0129E+00  3.9575E-01  1.1081E-01  9.4320E-01  1.9435E-01  9.5894E-01  1.0000E-02  9.0928E-01  1.6839E+00  5.4191E-01
             2.7979E+00
 PARAMETER:  1.1283E-01 -8.2697E-01 -2.0999E+00  4.1519E-02 -1.5381E+00  5.8075E-02 -6.3859E+00  4.8950E-03  6.2112E-01 -5.1266E-01
             1.1289E+00
 GRADIENT:   4.7328E+00 -4.6895E+00 -9.3080E-01  2.2552E+00  2.8156E+00 -1.9802E-01  0.0000E+00  3.6512E-01  1.4713E+00  9.9386E-01
            -4.1907E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1321.00558710225        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1949
 NPARAMETR:  1.0104E+00  4.0845E-01  1.1278E-01  9.4377E-01  1.9833E-01  9.5894E-01  1.0000E-02  9.1712E-01  1.6587E+00  5.2986E-01
             2.8076E+00
 PARAMETER:  1.1032E-01 -7.9539E-01 -2.0824E+00  4.2128E-02 -1.5178E+00  5.8070E-02 -5.6765E+00  1.3485E-02  6.0602E-01 -5.3515E-01
             1.1323E+00
 GRADIENT:  -1.3031E-01 -4.9146E-02  1.8327E-02 -8.7022E-02  9.5540E-02  8.6977E-03  0.0000E+00 -2.5756E-02 -1.4950E-02  1.3944E-02
            -7.4193E-04

0ITERATION NO.:   69    OBJECTIVE VALUE:  -1321.00573664584        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     2079
 NPARAMETR:  1.0105E+00  4.0851E-01  1.1275E-01  9.4392E-01  1.9837E-01  9.5890E-01  1.0000E-02  9.1810E-01  1.6601E+00  5.2958E-01
             2.8071E+00
 PARAMETER:  1.1040E-01 -7.9524E-01 -2.0826E+00  4.2288E-02 -1.5176E+00  5.8031E-02 -5.6645E+00  1.4554E-02  6.0685E-01 -5.3567E-01
             1.1322E+00
 GRADIENT:   1.8586E-01 -4.0777E-01 -2.8834E-01 -1.4930E-02  9.5145E-01 -3.3319E-04  0.0000E+00 -7.8577E-03  1.3909E-01  8.1922E-03
            -3.9673E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2079
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4394E-03 -1.6572E-04  2.3322E-02 -1.0320E-02  1.9712E-02
 SE:             2.8735E-02  1.7673E-04  1.3915E-02  2.5945E-02  2.0275E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3235E-01  3.4841E-01  9.3740E-02  6.9079E-01  3.3095E-01

 ETASHRINKSD(%)  3.7341E+00  9.9408E+01  5.3382E+01  1.3082E+01  3.2075E+01
 ETASHRINKVR(%)  7.3287E+00  9.9996E+01  7.8267E+01  2.4453E+01  5.3862E+01
 EBVSHRINKSD(%)  3.4633E+00  9.9362E+01  5.3308E+01  9.4432E+00  3.2761E+01
 EBVSHRINKVR(%)  6.8067E+00  9.9996E+01  7.8198E+01  1.7995E+01  5.4789E+01
 RELATIVEINF(%)  8.8878E+01  3.3565E-04  4.4028E+00  5.7602E+01  2.1118E+00
 EPSSHRINKSD(%)  3.4112E+01
 EPSSHRINKVR(%)  5.6588E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1321.0057366458368     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -585.85491008209863     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.73
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.64
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1321.006       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  4.09E-01  1.13E-01  9.44E-01  1.98E-01  9.59E-01  1.00E-02  9.18E-01  1.66E+00  5.30E-01  2.81E+00
 


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
+        4.57E+01  2.90E+03
 
 TH 3
+       -5.91E+02  4.83E+03  2.39E+04
 
 TH 4
+       -7.98E+00  8.58E+01 -5.98E+02  3.70E+02
 
 TH 5
+        3.64E+02 -1.06E+04 -2.46E+04 -3.91E+02  4.68E+04
 
 TH 6
+        1.75E+00 -1.12E+01  5.60E+01 -8.36E+00  4.56E+01  1.87E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.85E+00 -3.13E+00 -8.55E+01 -2.27E+00 -1.38E+01  4.04E+00  0.00E+00  1.99E+01
 
 TH 9
+        1.44E+01 -3.83E+01  2.33E+02 -1.24E+01  1.74E+02  9.79E-01  0.00E+00 -1.88E+00  4.01E+01
 
 TH10
+       -4.54E+00 -7.79E+01  4.65E+01  4.96E+00  4.72E+02  4.19E+00  0.00E+00  1.79E+01  1.06E+01  1.62E+02
 
 TH11
+       -2.07E+01 -1.52E+01 -8.28E+00  4.09E-01  4.29E+01  1.43E+00  0.00E+00  6.66E+00  6.08E+00  1.68E+01  3.41E+01
 
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
 #CPUT: Total CPU Time in Seconds,       34.416
Stop Time:
Wed Sep 29 13:32:14 CDT 2021
