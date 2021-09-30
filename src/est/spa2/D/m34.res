Thu Sep 30 08:57:07 CDT 2021
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
$DATA ../../../../data/spa2/D/dat34.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:   19901.9917370745        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8636E+02  1.8471E+02 -4.9021E+01  1.0459E+02  1.6071E+02 -1.6030E+03 -6.8149E+02 -2.1906E+01 -9.4901E+02 -6.3955E+02
            -4.0530E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -715.381212369923        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3039E+00  1.4714E+00  1.0308E+00  1.8013E+00  9.6529E-01  2.6472E+00  1.4693E+00  9.7192E-01  1.2959E+00  1.0704E+00
             1.3954E+01
 PARAMETER:  3.6533E-01  4.8622E-01  1.3029E-01  6.8849E-01  6.4672E-02  1.0735E+00  4.8481E-01  7.1514E-02  3.5923E-01  1.6803E-01
             2.7358E+00
 GRADIENT:  -1.3566E+01  2.3457E+01 -1.0218E+01  7.0046E+01  8.8658E+00  9.3609E+01 -1.3562E+01  4.7093E+00 -1.2890E+01  1.5684E+01
             2.6109E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -788.783423496276        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      153
 NPARAMETR:  1.2354E+00  1.5567E+00  6.8239E+00  1.8006E+00  2.1595E+00  2.4105E+00  4.8676E+00  2.9959E-01  3.3849E+00  1.2278E+00
             1.2426E+01
 PARAMETER:  3.1143E-01  5.4259E-01  2.0204E+00  6.8813E-01  8.6988E-01  9.7982E-01  1.6826E+00 -1.1053E+00  1.3193E+00  3.0520E-01
             2.6198E+00
 GRADIENT:  -1.0337E+01 -3.3461E-01 -7.4699E+00  1.1254E+01  4.2153E+00  5.4176E+01  4.1359E+01  1.7724E-02  7.0124E+01  1.0853E+01
             3.1401E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -819.697933752697        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.2305E+00  7.4015E-01  7.5096E+00  2.1552E+00  1.3481E+00  1.7684E+00  5.1531E+00  6.1684E-02  1.7530E+00  1.0088E+00
             1.2429E+01
 PARAMETER:  3.0738E-01 -2.0091E-01  2.1162E+00  8.6789E-01  3.9868E-01  6.7006E-01  1.7396E+00 -2.6857E+00  6.6132E-01  1.0872E-01
             2.6200E+00
 GRADIENT:  -1.3970E+01  4.1553E+00  7.8477E+00  4.4458E+01 -5.4845E+01 -1.6504E+01 -2.0685E+00  1.7363E-03  3.6660E+01  1.0336E+01
             2.9723E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -820.781140510072        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      387
 NPARAMETR:  1.2304E+00  7.4016E-01  6.8960E+00  2.1554E+00  1.3480E+00  1.7682E+00  5.1542E+00  5.9503E-02  1.7531E+00  9.8418E-01
             1.2420E+01
 PARAMETER:  3.0735E-01 -2.0088E-01  2.0309E+00  8.6798E-01  3.9863E-01  6.6999E-01  1.7398E+00 -2.7217E+00  6.6139E-01  8.4049E-02
             2.6193E+00
 GRADIENT:  -2.5822E+01  3.4651E+00  7.2833E+00  2.6664E+01 -5.2071E+01 -2.5389E+01 -3.7748E+01  1.8523E-03  3.4352E+01  9.9348E+00
             2.7321E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -868.147871045394        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      532            RESET HESSIAN, TYPE II
 NPARAMETR:  1.1831E+00  5.7941E-01  2.5293E+00  2.0044E+00  1.5302E+00  1.8399E+00  5.6470E+00  1.0000E-02  1.6731E+00  7.7173E-01
             9.3612E+00
 PARAMETER:  2.6817E-01 -4.4575E-01  1.0279E+00  7.9535E-01  5.2542E-01  7.0971E-01  1.8311E+00 -1.1806E+01  6.1468E-01 -1.5912E-01
             2.3366E+00
 GRADIENT:   2.3025E+01 -3.2760E-01 -4.7859E+01  5.3224E+01  3.8451E+01 -1.0532E+01  5.3538E+01  0.0000E+00  4.8567E+01  6.6497E+00
             6.1924E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -900.094464026417        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      604
 NPARAMETR:  1.0888E+00  7.7801E-01  6.4206E+00  1.5567E+00  1.7190E+00  1.8607E+00  5.5518E+00  1.0000E-02  8.6365E-01  2.8115E-01
             8.5994E+00
 PARAMETER:  1.8505E-01 -1.5102E-01  1.9595E+00  5.4256E-01  6.4174E-01  7.2095E-01  1.8141E+00 -1.1806E+01 -4.6591E-02 -1.1689E+00
             2.2517E+00
 GRADIENT:  -8.2098E+00 -5.6614E-01  2.5057E-01 -1.6463E+00 -5.1925E-01 -4.9465E+00  6.4643E+01  0.0000E+00  1.8922E+00  1.0054E+00
             5.5708E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -905.632718909519        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      765
 NPARAMETR:  1.1039E+00  8.5343E-01  6.5052E+00  1.5400E+00  1.7335E+00  1.8986E+00  5.4995E+00  1.0000E-02  8.1985E-01  2.4292E-01
             8.6854E+00
 PARAMETER:  1.9882E-01 -5.8491E-02  1.9726E+00  5.3177E-01  6.5014E-01  7.4113E-01  1.8047E+00 -1.1806E+01 -9.8634E-02 -1.3150E+00
             2.2616E+00
 GRADIENT:  -8.4952E+00 -5.6985E-01  4.8881E-01 -9.2622E+00 -4.0898E+00 -1.0361E+01 -1.1963E+01  0.0000E+00  8.6991E-02  7.2143E-01
            -2.0336E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -906.402800817272        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      927
 NPARAMETR:  1.1028E+00  8.5386E-01  7.0464E+00  1.5458E+00  1.8096E+00  1.8936E+00  5.5968E+00  1.0000E-02  8.0522E-01  3.9603E-02
             8.8271E+00
 PARAMETER:  1.9788E-01 -5.7983E-02  2.0525E+00  5.3556E-01  6.9309E-01  7.3846E-01  1.8222E+00 -1.1806E+01 -1.1664E-01 -3.1288E+00
             2.2778E+00
 GRADIENT:  -1.0907E+01 -1.1627E+00 -4.3169E-01 -1.5468E+01  2.8820E+00 -1.0917E+01 -8.1115E+00  0.0000E+00  3.4932E-01  1.8232E-02
             2.2569E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -907.127954124991        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1103
 NPARAMETR:  1.1038E+00  8.5380E-01  7.0225E+00  1.5630E+00  1.7921E+00  1.9117E+00  5.9701E+00  1.0000E-02  8.0876E-01  1.0000E-02
             8.7904E+00
 PARAMETER:  1.9871E-01 -5.8063E-02  2.0491E+00  5.4660E-01  6.8339E-01  7.4801E-01  1.8868E+00 -1.1806E+01 -1.1226E-01 -7.9875E+00
             2.2737E+00
 GRADIENT:  -9.0506E+00  2.5952E+00 -3.0737E-01 -1.6609E+01  1.0811E+00 -7.9941E+00  4.4256E+00  0.0000E+00  4.7367E-01  0.0000E+00
            -5.1065E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -907.594099574596        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1281
 NPARAMETR:  1.1056E+00  8.5373E-01  7.5465E+00  1.6254E+00  1.8278E+00  1.9526E+00  5.9973E+00  1.0000E-02  8.6682E-01  1.0000E-02
             8.8906E+00
 PARAMETER:  2.0038E-01 -5.8140E-02  2.1211E+00  5.8573E-01  7.0312E-01  7.6916E-01  1.8913E+00 -1.1806E+01 -4.2923E-02 -1.2916E+01
             2.2850E+00
 GRADIENT:  -1.0283E+01  5.0001E+00 -4.5089E-01 -1.9585E+00  4.8150E-01 -7.6355E-01  8.3656E-01  0.0000E+00  5.8625E-01  0.0000E+00
             3.3845E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -908.390848482737        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1464
 NPARAMETR:  1.1242E+00  7.4395E-01  9.8258E+00  1.6674E+00  1.8628E+00  1.9588E+00  6.2137E+00  1.0000E-02  8.4168E-01  1.0000E-02
             8.8544E+00
 PARAMETER:  2.1711E-01 -1.9578E-01  2.3850E+00  6.1128E-01  7.2210E-01  7.7234E-01  1.9268E+00 -1.1806E+01 -7.2360E-02 -1.4080E+01
             2.2809E+00
 GRADIENT:  -1.4768E+00  2.9543E+00  4.3095E-01  2.8654E+00 -3.2680E+00  9.3528E-01 -5.1017E-01  0.0000E+00 -3.2107E+00  0.0000E+00
            -5.4272E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -909.013591364360        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1640
 NPARAMETR:  1.1176E+00  6.4263E-01  1.0429E+01  1.7441E+00  1.8927E+00  1.9435E+00  6.4262E+00  1.0000E-02  9.7324E-01  1.0000E-02
             8.8804E+00
 PARAMETER:  2.1120E-01 -3.4218E-01  2.4446E+00  6.5626E-01  7.3799E-01  7.6447E-01  1.9604E+00 -1.1806E+01  7.2879E-02 -1.4080E+01
             2.2838E+00
 GRADIENT:  -4.8814E+00  1.5033E+00 -1.8182E-01 -8.8818E-01  5.7571E-01 -1.2628E+00 -1.0876E+00  0.0000E+00  1.8131E+00  0.0000E+00
            -3.6247E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -909.704563909268        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1824
 NPARAMETR:  1.1272E+00  5.2598E-01  1.4672E+01  1.8064E+00  1.9114E+00  1.9504E+00  7.1533E+00  1.0000E-02  9.7483E-01  1.0000E-02
             8.9025E+00
 PARAMETER:  2.1975E-01 -5.4249E-01  2.7859E+00  6.9132E-01  7.4785E-01  7.6806E-01  2.0676E+00 -1.1806E+01  7.4507E-02 -1.4080E+01
             2.2863E+00
 GRADIENT:   1.1139E-01  1.6581E+00  4.9826E-01 -3.1903E+00 -3.9617E+00  9.2905E-01  9.2946E+00  0.0000E+00 -3.1737E-01  0.0000E+00
            -7.3171E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -909.867752903279        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2007             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1280E+00  4.8823E-01  1.3228E+01  1.8328E+00  1.9297E+00  1.9460E+00  7.3092E+00  1.0000E-02  1.0013E+00  1.0000E-02
             8.9089E+00
 PARAMETER:  2.2045E-01 -6.1697E-01  2.6824E+00  7.0583E-01  7.5738E-01  7.6579E-01  2.0891E+00 -1.1806E+01  1.0134E-01 -1.4080E+01
             2.2870E+00
 GRADIENT:   1.3147E+01  3.7478E+00  8.3703E-03  1.7110E+01  1.5697E+00  1.5668E+01  1.1107E+02  0.0000E+00  1.0515E+00  0.0000E+00
             2.7968E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -909.909674980014        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2182
 NPARAMETR:  1.1309E+00  4.6519E-01  1.3410E+01  1.8539E+00  1.9324E+00  1.9516E+00  7.2601E+00  1.0000E-02  1.0042E+00  1.0000E-02
             8.9316E+00
 PARAMETER:  2.2301E-01 -6.6532E-01  2.6960E+00  7.1727E-01  7.5877E-01  7.6864E-01  2.0824E+00 -1.1806E+01  1.0416E-01 -1.4080E+01
             2.2896E+00
 GRADIENT:   1.0848E+00  2.2260E-01 -1.2773E-02  7.4274E-01 -2.6119E-01  1.5519E+00  6.1931E+00  0.0000E+00 -5.7092E-01  0.0000E+00
             1.8201E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -909.951788338468        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     2363
 NPARAMETR:  1.1303E+00  4.5011E-01  1.3713E+01  1.8636E+00  1.9341E+00  1.9492E+00  7.3691E+00  1.0000E-02  1.0121E+00  1.0000E-02
             8.9280E+00
 PARAMETER:  2.2246E-01 -6.9827E-01  2.7183E+00  7.2249E-01  7.5964E-01  7.6742E-01  2.0973E+00 -1.1806E+01  1.1203E-01 -1.4080E+01
             2.2892E+00
 GRADIENT:   9.9651E-01  2.7124E-01 -4.7154E-03  1.5301E-01 -3.1849E-01  1.2380E+00  7.6514E+00  0.0000E+00 -3.9585E-01  0.0000E+00
             1.2126E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -909.981282114492        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2551
 NPARAMETR:  1.1269E+00  4.2887E-01  1.4115E+01  1.8710E+00  1.9350E+00  1.9432E+00  7.4809E+00  1.0000E-02  1.0162E+00  1.0000E-02
             8.9079E+00
 PARAMETER:  2.1948E-01 -7.4661E-01  2.7472E+00  7.2650E-01  7.6008E-01  7.6433E-01  2.1124E+00 -1.1806E+01  1.1603E-01 -1.4080E+01
             2.2869E+00
 GRADIENT:   3.1269E-02 -1.0298E-02  2.4024E-02 -9.6494E-01 -3.8229E-01  3.5362E-01  8.7054E+00  0.0000E+00 -4.6897E-01  0.0000E+00
            -1.1539E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -909.992068725455        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2738             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1272E+00  4.1824E-01  1.4378E+01  1.8793E+00  1.9395E+00  1.9427E+00  7.5313E+00  1.0000E-02  1.0240E+00  1.0000E-02
             8.9129E+00
 PARAMETER:  2.1970E-01 -7.7171E-01  2.7657E+00  7.3092E-01  7.6242E-01  7.6407E-01  2.1191E+00 -1.1806E+01  1.2373E-01 -1.4080E+01
             2.2875E+00
 GRADIENT:   1.2796E+01  2.9818E+00  5.4098E-02  2.2309E+01  9.5951E-01  1.5776E+01  1.1725E+02  0.0000E+00  6.5204E-02  0.0000E+00
             2.6882E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -909.995034529842        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2918
 NPARAMETR:  1.1278E+00  4.2126E-01  1.4386E+01  1.8824E+00  1.9426E+00  1.9441E+00  7.5551E+00  1.0000E-02  1.0303E+00  1.0000E-02
             8.9200E+00
 PARAMETER:  2.2026E-01 -7.6449E-01  2.7663E+00  7.3254E-01  7.6404E-01  7.6482E-01  2.1222E+00 -1.1806E+01  1.2981E-01 -1.4080E+01
             2.2883E+00
 GRADIENT:   2.2160E-01  1.8122E-01 -3.1527E-02 -1.1234E+00  5.1960E-02  5.2950E-01  9.7230E+00  0.0000E+00  1.4723E-01  0.0000E+00
             2.1403E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -909.998284253463        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3103
 NPARAMETR:  1.1278E+00  4.1685E-01  1.4545E+01  1.8861E+00  1.9445E+00  1.9441E+00  7.5789E+00  1.0000E-02  1.0324E+00  1.0000E-02
             8.9207E+00
 PARAMETER:  2.2025E-01 -7.7504E-01  2.7772E+00  7.3451E-01  7.6502E-01  7.6478E-01  2.1254E+00 -1.1806E+01  1.3189E-01 -1.4080E+01
             2.2884E+00
 GRADIENT:   2.1122E-01  1.5508E-01 -2.6741E-02 -8.8056E-01  2.2129E-02  5.5416E-01  9.8346E+00  0.0000E+00  1.0733E-01  0.0000E+00
             1.6630E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -910.000803657104        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     3298             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1274E+00  4.1216E-01  1.4750E+01  1.8890E+00  1.9456E+00  1.9430E+00  7.5980E+00  1.0000E-02  1.0324E+00  1.0000E-02
             8.9175E+00
 PARAMETER:  2.1989E-01 -7.8635E-01  2.7912E+00  7.3606E-01  7.6555E-01  7.6425E-01  2.1279E+00 -1.1806E+01  1.3188E-01 -1.4080E+01
             2.2880E+00
 GRADIENT:   1.2790E+01  3.2383E+00  4.4065E-02  2.3285E+01  9.8743E-01  1.5834E+01  1.1913E+02  0.0000E+00  2.5919E-01  0.0000E+00
             2.7025E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -910.001268202412        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     3484
 NPARAMETR:  1.1274E+00  4.0975E-01  1.4956E+01  1.8922E+00  1.9460E+00  1.9429E+00  7.6123E+00  1.0000E-02  1.0328E+00  1.0000E-02
             8.9167E+00
 PARAMETER:  2.1989E-01 -7.9220E-01  2.8051E+00  7.3773E-01  7.6578E-01  7.6421E-01  2.1298E+00 -1.1806E+01  1.3229E-01 -1.4080E+01
             2.2879E+00
 GRADIENT:   9.1462E-02  1.2536E-01  3.0040E-02  3.8573E-01 -4.7809E-01  4.5004E-01  9.7952E+00  0.0000E+00 -3.1101E-01  0.0000E+00
            -8.5819E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -910.002517920190        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3671
 NPARAMETR:  1.1274E+00  4.0805E-01  1.4997E+01  1.8932E+00  1.9483E+00  1.9428E+00  7.6214E+00  1.0000E-02  1.0350E+00  1.0000E-02
             8.9183E+00
 PARAMETER:  2.1996E-01 -7.9636E-01  2.8079E+00  7.3828E-01  7.6696E-01  7.6411E-01  2.1310E+00 -1.1806E+01  1.3441E-01 -1.4080E+01
             2.2881E+00
 GRADIENT:   1.1246E-01  8.6793E-02  8.7724E-03 -3.8541E-02 -2.3435E-01  4.2148E-01  9.8968E+00  0.0000E+00 -1.4520E-01  0.0000E+00
            -5.3243E-01

0ITERATION NO.:  118    OBJECTIVE VALUE:  -910.002781776680        NO. OF FUNC. EVALS.: 105
 CUMULATIVE NO. OF FUNC. EVALS.:     3776
 NPARAMETR:  1.1273E+00  4.0480E-01  1.5241E+01  1.8963E+00  1.9494E+00  1.9429E+00  7.6388E+00  1.0000E-02  1.0350E+00  1.0000E-02
             8.9173E+00
 PARAMETER:  2.1994E-01 -7.9813E-01  2.8061E+00  7.3815E-01  7.6755E-01  7.6411E-01  2.1314E+00 -1.1806E+01  1.3575E-01 -1.4080E+01
             2.2882E+00
 GRADIENT:   1.1383E-02  3.3560E-02 -1.1451E-02 -3.5009E-01  2.2839E-03 -5.2198E-03 -8.4446E-02  0.0000E+00  2.8752E-02  0.0000E+00
             6.3738E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     3776
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7436E-02  6.0022E-02  3.3622E-06 -8.3292E-02  2.0432E-05
 SE:             2.8473E-02  2.1169E-02  4.6026E-06  1.5172E-02  5.9919E-05
 N:                     100         100         100         100         100

 P VAL.:         5.4029E-01  4.5772E-03  4.6508E-01  4.0291E-08  7.3310E-01

 ETASHRINKSD(%)  4.6112E+00  2.9081E+01  9.9985E+01  4.9173E+01  9.9799E+01
 ETASHRINKVR(%)  9.0097E+00  4.9706E+01  1.0000E+02  7.4166E+01  1.0000E+02
 EBVSHRINKSD(%)  6.4964E+00  2.7751E+01  9.9972E+01  4.4657E+01  9.9695E+01
 EBVSHRINKVR(%)  1.2571E+01  4.7801E+01  1.0000E+02  6.9371E+01  9.9999E+01
 RELATIVEINF(%)  8.6613E+01  2.8872E+01  1.0268E-06  1.5339E+01  1.1354E-04
 EPSSHRINKSD(%)  9.1186E+00
 EPSSHRINKVR(%)  1.7406E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -910.00278177668008     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       192.72345806892702     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    99.61
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    12.91
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -910.003       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.13E+00  4.07E-01  1.50E+01  1.89E+00  1.95E+00  1.94E+00  7.62E+00  1.00E-02  1.04E+00  1.00E-02  8.92E+00
 


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
+        3.51E+02
 
 TH 2
+       -1.08E+02  1.04E+02
 
 TH 3
+       -1.39E-01  4.11E-02  2.18E-04
 
 TH 4
+       -1.64E+02  8.19E+01  1.81E-01  1.78E+02
 
 TH 5
+        2.41E+01 -1.26E+01 -3.11E-02 -2.95E+01  4.97E+00
 
 TH 6
+       -3.99E+01 -1.00E+01  5.10E-02  2.40E+01 -4.63E+00  4.44E+01
 
 TH 7
+       -2.21E+00  6.70E+00 -6.31E-03 -1.71E+00  3.93E-01 -2.17E+00  8.69E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        6.22E+01 -8.26E+00 -9.55E-02 -7.53E+01  1.29E+01 -2.43E+01  3.64E+00  0.00E+00  4.34E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        2.79E+00 -4.88E+00 -8.69E-03 -9.52E+00  1.61E+00  8.83E-01  5.28E-02  0.00E+00  3.38E+00  0.00E+00  1.11E+00
 
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
+        1.91E+02
 
 TH 2
+       -4.54E+00  6.72E+01
 
 TH 3
+        2.20E-03  3.90E-02  5.94E-03
 
 TH 4
+       -7.32E+00  3.59E+01  6.09E-02  1.14E+02
 
 TH 5
+       -8.22E-01 -7.50E+00 -3.16E-01 -1.66E+01  2.57E+01
 
 TH 6
+       -8.93E-01 -3.54E+00  1.75E-02  9.44E-01 -1.58E+00  4.17E+01
 
 TH 7
+        4.29E-01  5.80E+00 -5.49E-03 -3.76E+00  4.11E-01 -2.95E-01  1.59E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.63E-01 -1.65E+00 -7.58E-02 -3.92E+01  6.94E+00 -1.29E+00  1.24E+00  0.00E+00  4.06E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.63E+00 -3.27E+00  1.82E-03 -8.86E+00  8.17E-01  1.30E+00  6.12E-02  0.00E+00  3.30E+00  0.00E+00  8.28E+00
 
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
+        2.05E+02
 
 TH 2
+        6.21E+01  6.86E+01
 
 TH 3
+        1.54E-01  6.55E-02  1.58E-03
 
 TH 4
+        8.41E+01  3.50E+01  1.54E-01  1.23E+02
 
 TH 5
+       -2.03E+01 -4.94E+00 -1.22E-01 -2.26E+01  1.24E+01
 
 TH 6
+        3.98E+01  1.87E+01 -1.43E-02 -8.24E+00  5.79E-01  5.18E+01
 
 TH 7
+       -5.39E-01  6.72E+00 -7.66E-03 -4.75E+00  1.66E+00  1.90E+00  1.86E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.34E+01 -2.20E+00 -1.60E-02 -3.65E+01  5.56E+00  7.18E+00  1.02E+00  0.00E+00  3.76E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.42E+01  4.03E+00 -4.68E-02 -1.03E+01  5.52E+00 -7.81E+00  2.42E+00  0.00E+00 -4.37E+00  0.00E+00  1.52E+02
 
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
 #CPUT: Total CPU Time in Seconds,      112.540
Stop Time:
Thu Sep 30 08:59:02 CDT 2021
