Wed Sep 29 08:34:50 CDT 2021
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
$DATA ../../../../data/int/D/dat34.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m34.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   28819.3335580630        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2261E+02  2.9867E+02 -8.1042E+01  1.6327E+02  1.5275E+02 -1.7066E+03 -7.8537E+02 -2.6916E+01 -1.1935E+03 -7.5840E+02
            -6.0496E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1035.92674255273        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2516E+00  1.5241E+00  8.9665E-01  1.8438E+00  9.1750E-01  4.6096E+00  3.1207E+00  9.6694E-01  2.4619E+00  1.4912E+00
             1.3065E+01
 PARAMETER:  3.2439E-01  5.2141E-01 -9.0907E-03  7.1181E-01  1.3899E-02  1.6281E+00  1.2381E+00  6.6386E-02  1.0009E+00  4.9961E-01
             2.6700E+00
 GRADIENT:  -5.7369E+00 -3.2555E+01 -3.6968E+01  7.2457E+01  1.6217E+01  1.7921E+02 -3.4703E+01  4.2173E+00  6.1679E+01  3.8398E+01
             5.6872E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1106.91845066724        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0645E+00  1.2265E+00  2.4209E+01  2.7781E+00  2.0985E+00  3.3497E+00  4.7181E+00  5.1470E-01  5.9086E+00  1.4203E+00
             1.2394E+01
 PARAMETER:  1.6254E-01  3.0417E-01  3.2867E+00  1.1218E+00  8.4124E-01  1.3089E+00  1.6514E+00 -5.6417E-01  1.8764E+00  4.5087E-01
             2.6172E+00
 GRADIENT:  -2.6810E+01 -1.5668E+01 -1.6583E+01  4.2679E+01 -5.7196E-01  1.3015E+02  3.1690E+01  2.2107E-01  1.4596E+02  3.3400E+01
             5.3349E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1278.30710996586        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  1.1299E+00  2.1451E-01  2.7413E+03  2.1325E+00  2.7221E+00  1.7421E+00  7.4680E+00  2.3163E+00  2.4739E+00  1.1605E+00
             9.0946E+00
 PARAMETER:  2.2215E-01 -1.4394E+00  8.0162E+00  8.5727E-01  1.1014E+00  6.5511E-01  2.1106E+00  9.3997E-01  1.0058E+00  2.4886E-01
             2.3077E+00
 GRADIENT:  -7.0781E+00 -1.4398E+00 -1.7738E-01 -2.8277E+01  7.0501E+01 -6.1041E+01  6.3856E+00  2.2863E-03 -3.1281E+01  2.1741E+01
             2.2558E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1307.30105852551        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.1315E+00  1.3571E-01  3.4202E+03  2.0880E+00  2.2434E+00  2.1342E+00  7.8078E+00  3.2388E-01  2.5432E+00  5.1537E-01
             7.9498E+00
 PARAMETER:  2.2351E-01 -1.8972E+00  8.2375E+00  8.3620E-01  9.0800E-01  8.5807E-01  2.1551E+00 -1.0274E+00  1.0334E+00 -5.6287E-01
             2.1731E+00
 GRADIENT:   9.1188E+00 -1.0486E+00 -1.2682E-01 -3.0684E+00  1.2640E+01 -9.2109E+00  3.7791E+00  7.6637E-05  3.2176E+00  3.9711E+00
             3.3751E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1315.91795082935        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  1.0981E+00  5.5153E-01  1.1762E+02  1.8217E+00  2.3004E+00  2.1361E+00  1.8918E+00  1.0000E-02  2.8286E+00  3.1010E-02
             7.8510E+00
 PARAMETER:  1.9356E-01 -4.9506E-01  4.8675E+00  6.9979E-01  9.3308E-01  8.5897E-01  7.3754E-01 -6.4694E+00  1.1398E+00 -3.3734E+00
             2.1606E+00
 GRADIENT:  -5.0994E+00 -1.5674E+01 -6.9949E+00 -5.2913E+00  3.9038E+01 -1.3798E+01  5.6039E+00  0.0000E+00  9.8428E+00  1.3183E-02
            -3.1727E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1358.61866102364        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      449
 NPARAMETR:  1.2063E+00  1.9309E+00  3.5074E+02  1.1054E+00  2.1754E+00  1.9868E+00  8.0223E-01  1.4392E-02  4.2028E+00  3.6094E-02
             8.0516E+00
 PARAMETER:  2.8757E-01  7.5798E-01  5.9601E+00  2.0024E-01  8.7723E-01  7.8651E-01 -1.2036E-01 -4.1411E+00  1.5357E+00 -3.2216E+00
             2.1859E+00
 GRADIENT:   6.6606E+01  9.6274E+01 -6.4374E-03  2.8982E+01 -4.8897E+01  1.1799E+01 -9.8147E+00 -1.1412E-05  1.0982E+01  3.5429E-03
            -7.3327E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1374.90704167423        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  1.1041E+00  2.1629E+00  1.7159E+03  7.3760E-01  2.3178E+00  2.0031E+00  7.7904E-01  1.6691E-01  5.8836E+00  1.0032E-01
             8.3262E+00
 PARAMETER:  1.9902E-01  8.7144E-01  7.5477E+00 -2.0435E-01  9.4063E-01  7.9470E-01 -1.4969E-01 -1.6903E+00  1.8722E+00 -2.1994E+00
             2.2194E+00
 GRADIENT:  -5.1430E-01  1.7430E+01  1.9922E-01  1.0218E+01 -4.1449E+00 -1.7864E+00 -7.2211E+00 -1.4408E-04  7.5838E+00  9.5968E-02
            -1.3281E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1376.60421728537        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      781
 NPARAMETR:  1.1048E+00  2.2766E+00  2.6770E+03  5.5906E-01  2.3139E+00  2.0133E+00  9.5086E-01  7.6432E-01  5.9418E+00  1.9432E-01
             8.3310E+00
 PARAMETER:  1.9963E-01  9.2267E-01  7.9924E+00 -4.8150E-01  9.3895E-01  7.9976E-01  4.9611E-02 -1.6877E-01  1.8820E+00 -1.5382E+00
             2.2200E+00
 GRADIENT:  -3.9448E-01 -7.0710E-01  1.4022E-01  8.7838E-01 -5.1243E-02  2.3167E-01 -2.8095E-01 -1.2864E-03 -8.7421E-02  3.9974E-01
            -1.6392E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1377.52510915133        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      959
 NPARAMETR:  1.1076E+00  2.3265E+00  5.1321E+02  5.0977E-01  2.2979E+00  2.0127E+00  9.8818E-01  3.7245E-02  6.1519E+00  5.8686E-02
             8.3294E+00
 PARAMETER:  2.0223E-01  9.4435E-01  6.3407E+00 -5.7380E-01  9.3202E-01  7.9946E-01  8.8111E-02 -3.1902E+00  1.9168E+00 -2.7356E+00
             2.2198E+00
 GRADIENT:   1.1553E+00 -2.5893E+00  1.8992E-01  1.0508E-01 -8.7768E-01  2.3803E-01  4.9316E-01 -4.9124E-05  1.4446E+00  3.6509E-02
            -2.1361E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1377.82318451806        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1135            RESET HESSIAN, TYPE II
 NPARAMETR:  1.1056E+00  2.3442E+00  3.8962E+02  5.0044E-01  2.2980E+00  2.0143E+00  9.9390E-01  2.2345E-02  6.4497E+00  3.1407E-02
             8.3400E+00
 PARAMETER:  2.0039E-01  9.5196E-01  6.0652E+00 -5.9227E-01  9.3204E-01  8.0025E-01  9.3881E-02 -3.7011E+00  1.9640E+00 -3.3607E+00
             2.2211E+00
 GRADIENT:   1.3297E+01  3.1508E+01  8.3184E-02  6.7390E+00  5.0885E+00  1.9293E+01  1.6115E+00 -2.2247E-05  1.4426E+02  1.2447E-02
             4.0795E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1377.94798728042        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1316             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1058E+00  2.3824E+00  3.8226E+02  4.6471E-01  2.2996E+00  2.0136E+00  1.0210E+00  2.0968E-02  6.7508E+00  1.0000E-02
             8.3447E+00
 PARAMETER:  2.0057E-01  9.6809E-01  6.0461E+00 -6.6635E-01  9.3274E-01  7.9991E-01  1.2074E-01 -3.7648E+00  2.0097E+00 -5.5390E+00
             2.2216E+00
 GRADIENT:   1.3502E+01  2.8405E+01  8.8334E-02  7.2485E+00  6.3697E+00  1.9170E+01  2.9684E+00 -2.0378E-05  1.5845E+02  0.0000E+00
             4.2853E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1378.06374341106        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     1477
 NPARAMETR:  1.1050E+00  2.4097E+00  3.7467E+02  4.4383E-01  2.2951E+00  2.0143E+00  9.8555E-01  1.9891E-02  6.6040E+00  1.0000E-02
             8.3371E+00
 PARAMETER:  1.9983E-01  9.7952E-01  6.0260E+00 -7.1232E-01  9.3077E-01  8.0026E-01  8.5446E-02 -3.8175E+00  1.9877E+00 -5.5422E+00
             2.2207E+00
 GRADIENT:   7.3389E-02  2.1996E+00  5.8174E-03  4.4341E-01 -1.4730E+00  7.0134E-01 -2.4490E+00 -1.9581E-05  7.2916E+00  0.0000E+00
            -3.0330E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1378.16693156237        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1657
 NPARAMETR:  1.1058E+00  2.4294E+00  3.7640E+02  4.3450E-01  2.2988E+00  2.0135E+00  1.0282E+00  2.0463E-02  6.8288E+00  1.0000E-02
             8.3441E+00
 PARAMETER:  2.0053E-01  9.8763E-01  6.0307E+00 -7.3356E-01  9.3238E-01  7.9986E-01  1.2776E-01 -3.7891E+00  2.0211E+00 -5.5422E+00
             2.2216E+00
 GRADIENT:   4.2891E-01 -5.2269E+00  2.4892E-02  1.3729E+00  6.6658E-01  4.6755E-01  7.5106E-01 -2.0903E-05  1.3935E+01  0.0000E+00
             4.8733E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1378.24625338763        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1838             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1053E+00  2.4600E+00  3.7076E+02  4.0563E-01  2.3006E+00  2.0128E+00  1.0567E+00  2.0164E-02  6.9293E+00  1.0000E-02
             8.3494E+00
 PARAMETER:  2.0012E-01  1.0002E+00  6.0155E+00 -8.0232E-01  9.3317E-01  7.9955E-01  1.5511E-01 -3.8038E+00  2.0358E+00 -5.5422E+00
             2.2222E+00
 GRADIENT:   1.3031E+01  3.5914E+01  6.5907E-02  6.6406E+00  5.9480E+00  1.9218E+01  2.2749E+00 -1.9211E-05  1.6618E+02  0.0000E+00
             4.2602E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1378.30516048713        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2017
 NPARAMETR:  1.1048E+00  2.4831E+00  3.6765E+02  3.8892E-01  2.2997E+00  2.0129E+00  1.0574E+00  2.1336E-02  6.9538E+00  1.0000E-02
             8.3487E+00
 PARAMETER:  1.9966E-01  1.0095E+00  6.0071E+00 -8.4439E-01  9.3276E-01  7.9956E-01  1.5578E-01 -3.7474E+00  2.0393E+00 -5.5422E+00
             2.2221E+00
 GRADIENT:  -1.2970E-01 -3.9862E+00  3.0601E-03 -1.4637E-02  8.6190E-01  4.9431E-01  4.8849E-01 -2.3039E-05  1.2171E+01  0.0000E+00
             1.4518E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1378.34144556139        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2200
 NPARAMETR:  1.1068E+00  2.5192E+00  3.6268E+02  3.7574E-01  2.2995E+00  2.0166E+00  1.0663E+00  2.9088E-02  7.0717E+00  1.0000E-02
             8.3558E+00
 PARAMETER:  2.0143E-01  1.0240E+00  5.9935E+00 -8.7885E-01  9.3269E-01  8.0141E-01  1.6421E-01 -3.4374E+00  2.0561E+00 -5.5422E+00
             2.2230E+00
 GRADIENT:   5.1366E-01  7.4441E-01 -1.8141E-03  5.3581E-01 -2.4154E-01  1.1817E+00 -5.5203E-01 -4.3418E-05  1.2999E+01  0.0000E+00
             9.8201E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1378.35511086261        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2381
 NPARAMETR:  1.1054E+00  2.5251E+00  3.6321E+02  3.6777E-01  2.2980E+00  2.0123E+00  1.0729E+00  3.0269E-02  7.1188E+00  1.0000E-02
             8.3456E+00
 PARAMETER:  2.0016E-01  1.0263E+00  5.9950E+00 -9.0030E-01  9.3204E-01  7.9929E-01  1.7040E-01 -3.3976E+00  2.0627E+00 -5.5422E+00
             2.2217E+00
 GRADIENT:   1.0676E-01 -3.8073E-01  5.7788E-03  3.4961E-01 -3.1178E-01  4.5249E-01 -2.5143E-01 -4.7784E-05  1.3600E+01  0.0000E+00
            -1.0254E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1378.36120025888        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2567
 NPARAMETR:  1.1053E+00  2.5306E+00  3.6211E+02  3.6065E-01  2.2989E+00  2.0116E+00  1.0801E+00  5.3654E-02  7.1604E+00  1.0000E-02
             8.3464E+00
 PARAMETER:  2.0011E-01  1.0284E+00  5.9920E+00 -9.1984E-01  9.3242E-01  7.9894E-01  1.7708E-01 -2.8252E+00  2.0686E+00 -5.5422E+00
             2.2218E+00
 GRADIENT:   8.1306E-02 -1.9530E+00  3.0496E-03  1.3215E-01  2.4276E-01  3.3566E-01  2.0632E-01 -1.5050E-04  1.4133E+01  0.0000E+00
            -3.0215E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1378.36374484934        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2748
 NPARAMETR:  1.1055E+00  2.5360E+00  3.6141E+02  3.5680E-01  2.2992E+00  2.0127E+00  1.0815E+00  7.3193E-02  7.1578E+00  1.0000E-02
             8.3491E+00
 PARAMETER:  2.0034E-01  1.0306E+00  5.9900E+00 -9.3058E-01  9.3257E-01  7.9949E-01  1.7837E-01 -2.5147E+00  2.0682E+00 -5.5422E+00
             2.2222E+00
 GRADIENT:   1.2635E-01 -1.3025E+00 -6.2860E-04 -6.0055E-02  2.1459E-01  5.4546E-01  7.0722E-03 -2.7954E-04  1.3484E+01  0.0000E+00
             1.6278E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1378.36658807250        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2940             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1054E+00  2.5433E+00  3.6073E+02  3.5540E-01  2.2980E+00  2.0127E+00  1.0822E+00  9.3152E-02  7.1921E+00  1.0000E-02
             8.3471E+00
 PARAMETER:  2.0023E-01  1.0335E+00  5.9881E+00 -9.3452E-01  9.3203E-01  7.9946E-01  1.7903E-01 -2.2735E+00  2.0730E+00 -5.5422E+00
             2.2219E+00
 GRADIENT:   1.2984E+01  4.6863E+01  6.0754E-02  6.9209E+00  3.9015E+00  1.9352E+01  4.2085E-01 -4.4165E-04  1.7708E+02  0.0000E+00
             3.9106E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1378.52901675746        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     3108
 NPARAMETR:  1.1054E+00  2.5611E+00  3.5943E+02  3.4948E-01  2.2968E+00  2.0128E+00  1.0802E+00  2.7154E+00  7.2290E+00  1.0000E-02
             8.3454E+00
 PARAMETER:  2.0025E-01  1.0404E+00  5.9845E+00 -9.5131E-01  9.3152E-01  7.9951E-01  1.7717E-01  1.0989E+00  2.0781E+00 -5.5422E+00
             2.2217E+00
 GRADIENT:   1.2971E+01  5.1583E+01 -1.9197E-02  7.1363E+00  4.2110E+00  1.9401E+01 -8.0034E-01  6.9352E-03  1.7772E+02  0.0000E+00
             3.8307E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1378.54781059790        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     3270
 NPARAMETR:  1.1055E+00  2.5502E+00  3.7369E+02  3.4905E-01  2.2969E+00  2.0128E+00  1.0945E+00  2.7800E+00  7.2232E+00  1.0000E-02
             8.3483E+00
 PARAMETER:  2.0034E-01  1.0362E+00  6.0234E+00 -9.5254E-01  9.3157E-01  7.9952E-01  1.9028E-01  1.1224E+00  2.0773E+00 -5.5422E+00
             2.2221E+00
 GRADIENT:   1.6443E-01 -1.4022E+00 -8.9685E-02  2.0263E-02  9.2632E-01  5.6078E-01  4.1836E-01 -2.3962E-03  1.4162E+01  0.0000E+00
             1.1023E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1378.55222699881        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3453
 NPARAMETR:  1.1055E+00  2.5518E+00  3.9413E+02  3.4763E-01  2.2952E+00  2.0129E+00  1.0928E+00  2.8587E+00  7.2351E+00  1.0000E-02
             8.3460E+00
 PARAMETER:  2.0027E-01  1.0368E+00  6.0767E+00 -9.5660E-01  9.3081E-01  7.9957E-01  1.8876E-01  1.1504E+00  2.0789E+00 -5.5422E+00
             2.2218E+00
 GRADIENT:   1.6846E-01 -9.5044E-01 -9.6411E-03  3.6946E-02  1.6815E-01  5.9228E-01  1.8473E-01  1.2113E-03  1.4308E+01  0.0000E+00
             3.8133E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1378.55282813369        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     3629
 NPARAMETR:  1.1051E+00  2.5565E+00  3.9593E+02  3.4547E-01  2.2941E+00  2.0126E+00  1.0898E+00  2.8669E+00  7.2513E+00  1.0000E-02
             8.3427E+00
 PARAMETER:  2.0008E-01  1.0372E+00  6.0800E+00 -9.5857E-01  9.3068E-01  7.9939E-01  1.8787E-01  1.1513E+00  2.0801E+00 -5.5422E+00
             2.2216E+00
 GRADIENT:   1.4372E-02 -3.1063E-01 -5.8213E-04  3.3128E-02  5.7791E-02 -3.4014E-03  5.0552E-02 -9.3014E-04 -5.8023E-02  0.0000E+00
             1.1037E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3629
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4112E-02 -8.6240E-02 -8.2162E-04  5.5578E-02  4.3634E-06
 SE:             2.8536E-02  1.5788E-02  2.9363E-03  2.1958E-02  1.4462E-04
 N:                     100         100         100         100         100

 P VAL.:         6.2094E-01  4.7071E-08  7.7962E-01  1.1371E-02  9.7593E-01

 ETASHRINKSD(%)  4.4019E+00  4.7109E+01  9.0163E+01  2.6438E+01  9.9515E+01
 ETASHRINKVR(%)  8.6100E+00  7.2026E+01  9.9032E+01  4.5886E+01  9.9998E+01
 EBVSHRINKSD(%)  5.4246E+00  4.2917E+01  9.1914E+01  2.3608E+01  9.9464E+01
 EBVSHRINKVR(%)  1.0555E+01  6.7415E+01  9.9346E+01  4.1642E+01  9.9997E+01
 RELATIVEINF(%)  8.9311E+01  1.6937E+01  6.2896E-01  3.0982E+01  2.6412E-03
 EPSSHRINKSD(%)  6.6766E+00
 EPSSHRINKVR(%)  1.2907E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1378.5528281336860     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       275.53653163472472     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   128.97
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    18.10
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1378.553       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.11E+00  2.55E+00  3.95E+02  3.47E-01  2.29E+00  2.01E+00  1.09E+00  2.86E+00  7.24E+00  1.00E-02  8.34E+00
 


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
+        2.61E+02
 
 TH 2
+       -5.34E+01  5.32E+01
 
 TH 3
+       -1.41E-03  1.05E-03  2.15E-08
 
 TH 4
+       -1.36E+02  6.88E+01  1.55E-03  1.80E+02
 
 TH 5
+        2.70E+01 -2.72E+01 -5.38E-04 -3.45E+01  1.40E+01
 
 TH 6
+        1.95E+01 -1.50E+01 -2.91E-04 -8.80E+00  7.82E+00  6.48E+00
 
 TH 7
+        7.53E+00 -2.58E+01 -4.60E-04 -9.59E+00  1.35E+01  1.01E+01  1.89E+01
 
 TH 8
+        7.94E-02 -6.29E-02 -1.29E-06 -9.63E-02  3.20E-02  1.64E-02  2.66E-02  7.86E-05
 
 TH 9
+       -7.69E+00  1.17E+00  4.19E-05  1.04E+01 -5.17E-01  7.21E-01  1.77E+00 -2.84E-03  8.94E-01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -7.81E-01 -4.93E+00 -9.05E-05 -5.07E+00  2.49E+00  1.14E+00  2.73E+00  6.26E-03  1.70E-02  0.00E+00  1.04E+00
 
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
+        1.93E+02
 
 TH 2
+       -2.78E+00  5.90E+01
 
 TH 3
+       -9.27E-05  6.18E-05  5.44E-06
 
 TH 4
+       -1.00E+00  3.63E+01  2.79E-04  1.14E+02
 
 TH 5
+       -1.60E+00 -1.37E+01 -3.16E-03 -1.02E+01  5.64E+01
 
 TH 6
+        1.44E+00  1.78E-01  7.28E-05 -9.64E-01 -1.52E+00  4.14E+01
 
 TH 7
+        3.39E-01 -2.98E+01  2.92E-05 -2.74E+00  6.80E+00 -6.61E-01  3.94E+01
 
 TH 8
+       -4.86E-03  3.20E-04 -3.14E-04 -2.82E-02  1.70E-01 -3.11E-03  1.37E-02  1.00E-01
 
 TH 9
+        3.55E-01 -3.38E+00  1.09E-04  7.59E+00  5.38E-01 -4.70E-02  1.57E+00 -5.50E-03  1.86E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.67E+00 -6.45E+00 -1.14E-04 -4.86E+00  7.66E-01  8.53E-01  3.61E+00  3.11E-02  1.18E-01  0.00E+00  1.52E+01
 
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
+        2.03E+02
 
 TH 2
+        6.20E+01  6.45E+01
 
 TH 3
+       -5.66E-04  2.31E-05  6.14E-06
 
 TH 4
+        8.73E+01  3.66E+01  4.45E-04  1.13E+02
 
 TH 5
+       -2.43E+01 -1.93E+01 -3.88E-03 -1.62E+01  4.93E+01
 
 TH 6
+        3.88E+01 -6.46E+00 -3.64E-04  2.28E+01  6.19E-01  4.83E+01
 
 TH 7
+       -1.27E+01 -2.65E+01  4.39E-06 -3.28E+00  4.05E+00  7.36E+00  3.75E+01
 
 TH 8
+       -1.09E-02 -8.85E-04  8.75E-06  1.22E-02 -1.65E-02 -5.86E-03 -3.48E-05  3.78E-03
 
 TH 9
+       -6.11E-01 -3.96E+00  2.00E-04  8.41E+00  1.07E+00  2.17E+00  1.35E+00  5.36E-03  2.15E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.69E+01 -1.91E+01 -1.31E-03  1.09E+01  1.07E+01 -9.23E+00 -3.16E+00 -1.54E-02  4.42E+00  0.00E+00  4.87E+02
 
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
 #CPUT: Total CPU Time in Seconds,      147.153
Stop Time:
Wed Sep 29 08:37:19 CDT 2021
