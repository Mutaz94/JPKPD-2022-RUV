Sat Sep 25 12:40:57 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat9.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1632.99472768407        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.8240E+02 -5.6232E+01 -2.7674E+01 -4.0236E+01  6.4978E+01 -2.4672E+00  6.7274E-01  2.8963E+00  5.5433E+00 -6.5356E+00
             5.1872E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1638.14018479751        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9881E-01  9.7069E-01  9.7485E-01  1.0131E+00  9.5370E-01  9.9751E-01  9.9708E-01  9.8783E-01  9.8102E-01  1.0328E+00
             8.6504E-01
 PARAMETER:  9.8813E-02  7.0249E-02  7.4531E-02  1.1303E-01  5.2591E-02  9.7504E-02  9.7074E-02  8.7752E-02  8.0840E-02  1.3230E-01
            -4.4983E-02
 GRADIENT:   1.9444E+02 -4.7253E+01 -1.7458E+01 -4.8758E+01  3.2052E+01 -3.1377E+00 -3.6459E-01  4.4926E+00  1.5920E+00  3.4398E+00
             2.6347E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1640.08599309325        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      267
 NPARAMETR:  9.8957E-01  7.2833E-01  1.0221E+00  1.2297E+00  8.4403E-01  1.0465E+00  1.0602E+00  8.0498E-01  9.1639E-01  1.0250E+00
             8.6223E-01
 PARAMETER:  8.9514E-02 -2.1700E-01  1.2190E-01  3.0676E-01 -6.9564E-02  1.4545E-01  1.5850E-01 -1.1693E-01  1.2685E-02  1.2467E-01
            -4.8235E-02
 GRADIENT:   1.1712E+02  1.7510E+01  2.3557E+00  4.6822E+01 -1.3498E+01  1.4355E+01 -1.0979E+00 -1.9517E+00  7.0616E+00  4.3938E+00
             6.2642E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1645.13997233038        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      444
 NPARAMETR:  9.3672E-01  8.0029E-01  9.9842E-01  1.1586E+00  8.6213E-01  9.8336E-01  1.2080E+00  8.8517E-01  8.9056E-01  9.6524E-01
             8.5213E-01
 PARAMETER:  3.4634E-02 -1.2278E-01  9.8423E-02  2.4725E-01 -4.8354E-02  8.3225E-02  2.8894E-01 -2.1980E-02 -1.5902E-02  6.4619E-02
            -6.0021E-02
 GRADIENT:   6.0108E+00  1.6848E+00  8.5282E-01  6.3456E-01 -4.7894E+00 -2.4948E+00  1.2172E-01  1.3054E+00  2.8400E+00  6.2748E-01
            -3.4522E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1645.30772269808        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  9.3449E-01  6.8893E-01  1.0409E+00  1.2313E+00  8.4046E-01  9.8719E-01  1.3747E+00  8.8426E-01  8.3474E-01  9.6764E-01
             8.5988E-01
 PARAMETER:  3.2241E-02 -2.7262E-01  1.4006E-01  3.0809E-01 -7.3806E-02  8.7103E-02  4.1823E-01 -2.3007E-02 -8.0637E-02  6.7108E-02
            -5.0960E-02
 GRADIENT:   3.1618E+00  4.3286E+00  2.6081E+00  4.6904E+00 -7.3737E+00 -4.2479E-01  2.5879E-01  7.4428E-01 -6.2010E-01  1.0278E+00
             2.5252E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1645.61219387995        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      799
 NPARAMETR:  9.2848E-01  4.6504E-01  1.2535E+00  1.3806E+00  8.6474E-01  9.8639E-01  1.4824E+00  1.0541E+00  8.0727E-01  1.0317E+00
             8.6030E-01
 PARAMETER:  2.5791E-02 -6.6564E-01  3.2596E-01  4.2249E-01 -4.5324E-02  8.6301E-02  4.9365E-01  1.5268E-01 -1.1410E-01  1.3116E-01
            -5.0478E-02
 GRADIENT:  -2.7076E+00  2.6937E+00  1.8613E-01  8.1762E+00 -5.9989E-01  6.3180E-01  5.4018E-01  3.9685E-01  1.4313E+00  5.1792E-01
            -1.5305E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1645.91492481044        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      975
 NPARAMETR:  9.2611E-01  2.8433E-01  1.3462E+00  1.4872E+00  8.5031E-01  9.8067E-01  1.4863E+00  1.1425E+00  7.7227E-01  1.0408E+00
             8.6004E-01
 PARAMETER:  2.3239E-02 -1.1576E+00  3.9726E-01  4.9691E-01 -6.2159E-02  8.0478E-02  4.9627E-01  2.3326E-01 -1.5842E-01  1.3998E-01
            -5.0774E-02
 GRADIENT:  -1.7898E+00 -3.4901E-01 -2.9665E-01 -4.8909E+00  1.3639E+00 -6.2822E-01  2.0143E-01 -2.0862E-02  2.7745E-01 -1.5704E-01
            -1.7199E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1646.02330362778        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1150
 NPARAMETR:  9.2503E-01  1.6409E-01  1.4261E+00  1.5667E+00  8.4422E-01  9.8022E-01  1.3207E+00  1.2354E+00  7.4650E-01  1.0461E+00
             8.6059E-01
 PARAMETER:  2.2071E-02 -1.7074E+00  4.5497E-01  5.4896E-01 -6.9341E-02  8.0019E-02  3.7819E-01  3.1141E-01 -1.9236E-01  1.4505E-01
            -5.0139E-02
 GRADIENT:   3.2559E-01  3.8952E-01  8.8228E-01  3.6330E+00 -2.1302E+00 -1.5656E-01  1.4751E-01  7.1175E-02  4.4796E-01  5.1463E-01
             1.5099E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1646.05906616973        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1325
 NPARAMETR:  9.2398E-01  9.9239E-02  1.4632E+00  1.6060E+00  8.4090E-01  9.7971E-01  1.0515E+00  1.2794E+00  7.2967E-01  1.0418E+00
             8.6036E-01
 PARAMETER:  2.0937E-02 -2.2102E+00  4.8060E-01  5.7373E-01 -7.3283E-02  7.9504E-02  1.5026E-01  3.4635E-01 -2.1516E-01  1.4096E-01
            -5.0407E-02
 GRADIENT:   4.1018E-01 -2.6108E-02  1.0778E-01 -2.6726E-01 -5.7776E-02 -1.5522E-02  5.4533E-02  4.2556E-03  2.3203E-01  7.8242E-02
             7.8963E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1646.05959970040        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1503
 NPARAMETR:  9.2362E-01  8.9074E-02  1.4647E+00  1.6122E+00  8.3885E-01  9.7962E-01  9.8976E-01  1.2826E+00  7.2656E-01  1.0396E+00
             8.6014E-01
 PARAMETER:  2.0543E-02 -2.3183E+00  4.8165E-01  5.7760E-01 -7.5719E-02  7.9413E-02  8.9709E-02  3.4887E-01 -2.1944E-01  1.3887E-01
            -5.0666E-02
 GRADIENT:  -8.2701E-02 -2.4733E-02 -1.4917E-02 -8.8522E-02  6.8516E-02  6.3490E-04  4.0262E-02 -1.4439E-02  3.5570E-03 -2.9260E-02
            -2.1453E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1646.05961922504        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1679
 NPARAMETR:  9.2359E-01  8.6983E-02  1.4661E+00  1.6135E+00  8.3883E-01  9.7959E-01  9.7187E-01  1.2844E+00  7.2609E-01  1.0399E+00
             8.6016E-01
 PARAMETER:  2.0516E-02 -2.3420E+00  4.8263E-01  5.7840E-01 -7.5747E-02  7.9382E-02  7.1470E-02  3.5029E-01 -2.2009E-01  1.3908E-01
            -5.0632E-02
 GRADIENT:  -4.9311E-02 -3.3276E-02 -5.6935E-02 -2.1230E-01  9.8529E-02 -7.9040E-05  3.7693E-02  1.4377E-02  5.4322E-02  1.1894E-02
             2.7690E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1646.05979621424        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1858
 NPARAMETR:  9.2354E-01  8.1456E-02  1.4699E+00  1.6169E+00  8.3877E-01  9.7951E-01  9.0280E-01  1.2892E+00  7.2481E-01  1.0404E+00
             8.6024E-01
 PARAMETER:  2.0460E-02 -2.4077E+00  4.8522E-01  5.8052E-01 -7.5822E-02  7.9299E-02 -2.2592E-03  3.5401E-01 -2.2185E-01  1.3962E-01
            -5.0541E-02
 GRADIENT:   5.9978E-02 -5.3862E-02 -1.6364E-01 -5.2905E-01  1.7794E-01 -3.9810E-03  3.0004E-02  8.5761E-02  1.6289E-01  1.1426E-01
             6.5296E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1646.09211483146        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2040
 NPARAMETR:  9.2440E-01  1.0598E-01  1.4648E+00  1.6031E+00  8.4268E-01  9.8059E-01  1.2212E-01  1.2806E+00  7.3139E-01  1.0428E+00
             8.6031E-01
 PARAMETER:  2.1392E-02 -2.1445E+00  4.8173E-01  5.7192E-01 -7.1167E-02  8.0395E-02 -2.0028E+00  3.4736E-01 -2.1280E-01  1.4187E-01
            -5.0463E-02
 GRADIENT:   1.1418E+00  9.5737E-02  3.9009E-01  2.9378E+00 -8.5088E-01  2.9782E-01  1.0447E-03 -3.6463E-03 -4.5095E-01 -5.2174E-03
            -1.9989E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1646.14025030648        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2219
 NPARAMETR:  9.2468E-01  1.8136E-01  1.4259E+00  1.5527E+00  8.5114E-01  9.8046E-01  1.5303E-02  1.2334E+00  7.5681E-01  1.0530E+00
             8.6057E-01
 PARAMETER:  2.1694E-02 -1.6073E+00  4.5480E-01  5.4001E-01 -6.1179E-02  8.0267E-02 -4.0797E+00  3.0975E-01 -1.7864E-01  1.5165E-01
            -5.0164E-02
 GRADIENT:  -1.2157E+00 -3.8483E-01 -5.0389E-01 -3.1474E+00  7.5774E-01 -1.3264E-01  5.1355E-05  1.2356E-01  3.9957E-01  2.4865E-01
             1.5754E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1646.16937609612        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2397
 NPARAMETR:  9.2658E-01  2.4461E-01  1.4072E+00  1.5146E+00  8.6123E-01  9.8188E-01  1.0000E-02  1.2059E+00  7.7636E-01  1.0579E+00
             8.6061E-01
 PARAMETER:  2.3745E-02 -1.3081E+00  4.4158E-01  5.1512E-01 -4.9391E-02  8.1710E-02 -5.8046E+00  2.8724E-01 -1.5314E-01  1.5633E-01
            -5.0119E-02
 GRADIENT:   9.7640E-01 -4.8873E-03  2.3342E-01  2.4570E-01 -2.9291E-02  1.0672E-01  0.0000E+00 -1.1295E-01 -2.7403E-01 -2.1025E-01
            -6.8416E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1646.18244339729        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2574
 NPARAMETR:  9.2702E-01  3.0030E-01  1.3868E+00  1.4794E+00  8.6995E-01  9.8230E-01  1.0000E-02  1.1811E+00  7.9644E-01  1.0662E+00
             8.6082E-01
 PARAMETER:  2.4218E-02 -1.1030E+00  4.2697E-01  4.9167E-01 -3.9321E-02  8.2138E-02 -7.3426E+00  2.6643E-01 -1.2761E-01  1.6411E-01
            -4.9869E-02
 GRADIENT:  -5.3619E-02  3.0915E-02  1.2482E-02  3.6946E-01 -8.2352E-02 -5.3901E-03  0.0000E+00  5.4088E-03 -1.0007E-02  2.7310E-02
             3.8756E-03

0ITERATION NO.:   78    OBJECTIVE VALUE:  -1646.18436390703        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     2666
 NPARAMETR:  9.2737E-01  3.2320E-01  1.3797E+00  1.4648E+00  8.7407E-01  9.8259E-01  1.0000E-02  1.1725E+00  8.0478E-01  1.0691E+00
             8.6088E-01
 PARAMETER:  2.4595E-02 -1.0295E+00  4.2188E-01  4.8174E-01 -3.4591E-02  8.2440E-02 -8.0102E+00  2.5911E-01 -1.1718E-01  1.6679E-01
            -4.9801E-02
 GRADIENT:  -2.6308E-02 -1.5497E-02 -1.1396E-02 -9.2186E-02  1.3728E-02  4.5857E-05  0.0000E+00  8.8518E-04  8.0682E-03  7.4889E-03
            -5.8293E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2666
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.0966E-04 -1.2798E-04 -3.5119E-02 -5.1503E-03 -3.5841E-02
 SE:             2.9868E-02  6.0112E-05  1.7395E-02  2.9393E-02  2.2150E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8906E-01  3.3252E-02  4.3491E-02  8.6091E-01  1.0563E-01

 ETASHRINKSD(%)  1.0000E-10  9.9799E+01  4.1726E+01  1.5288E+00  2.5795E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.6041E+01  3.0342E+00  4.4936E+01
 EBVSHRINKSD(%)  3.2418E-01  9.9810E+01  4.6117E+01  1.9697E+00  2.1883E+01
 EBVSHRINKVR(%)  6.4732E-01  1.0000E+02  7.0966E+01  3.9005E+00  3.8977E+01
 RELATIVEINF(%)  9.8062E+01  2.2180E-05  6.9059E+00  7.2061E+00  9.4650E+00
 EPSSHRINKSD(%)  4.6342E+01
 EPSSHRINKVR(%)  7.1208E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1646.1843639070275     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -911.03353734328937     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    31.91
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1646.184       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.27E-01  3.23E-01  1.38E+00  1.46E+00  8.74E-01  9.83E-01  1.00E-02  1.17E+00  8.05E-01  1.07E+00  8.61E-01
 


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
+        1.33E+03
 
 TH 2
+       -2.64E+01  4.04E+02
 
 TH 3
+        2.64E+00  5.89E+01  1.32E+02
 
 TH 4
+       -8.40E+00  4.85E+02 -2.29E+01  7.69E+02
 
 TH 5
+        4.58E+00 -2.40E+02 -2.54E+02 -4.69E+01  7.93E+02
 
 TH 6
+        3.64E+00 -2.39E+00  1.06E+00 -1.82E+00 -3.92E+00  2.02E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        6.02E-01 -3.65E+00 -3.22E+01 -3.90E+00 -1.25E+00  1.31E+00  0.00E+00  3.01E+01
 
 TH 9
+        4.79E-01 -1.04E+02  3.61E+00  6.36E-02 -6.61E-01 -2.98E+00  0.00E+00 -1.73E-01  2.82E+02
 
 TH10
+        1.50E+00  8.50E+00 -5.74E+00 -8.55E-01 -7.85E+01 -1.54E+00  0.00E+00  1.77E+01  4.18E+00  6.99E+01
 
 TH11
+       -7.02E+00 -1.26E+01 -1.09E+01 -8.46E+00  2.23E+00  1.88E-02  0.00E+00  7.95E+00  1.37E+01  1.19E+01  2.75E+02
 
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
 #CPUT: Total CPU Time in Seconds,       37.443
Stop Time:
Sat Sep 25 12:41:36 CDT 2021
