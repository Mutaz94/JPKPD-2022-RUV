Wed Sep 29 21:22:52 CDT 2021
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
$DATA ../../../../data/spa1/B/dat78.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2126.56916069826        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5697E+02 -7.4139E+01 -5.2368E+01 -2.9342E+01  3.7489E+01  3.2153E+01 -1.5038E+01  1.7939E+01  1.1352E+01  4.1216E+00
             2.8383E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2142.80398646741        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  9.7669E-01  1.2170E+00  1.3642E+00  9.7084E-01  1.2534E+00  9.9802E-01  1.1593E+00  8.3990E-01  9.1184E-01  1.0058E+00
             9.3791E-01
 PARAMETER:  7.6416E-02  2.9642E-01  4.1053E-01  7.0402E-02  3.2583E-01  9.8018E-02  2.4784E-01 -7.4474E-02  7.7143E-03  1.0581E-01
             3.5897E-02
 GRADIENT:   6.2711E+00  1.0784E+01 -4.5313E+00  2.4516E+01  4.2564E+01 -4.4905E-01 -4.1370E+00  2.3165E-01 -1.4552E+00 -4.4393E+01
            -3.8358E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2147.73820372130        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  9.6994E-01  1.1555E+00  1.7450E+00  9.9985E-01  1.3489E+00  1.0046E+00  1.2251E+00  6.9629E-01  9.2766E-01  1.2944E+00
             9.7251E-01
 PARAMETER:  6.9482E-02  2.4451E-01  6.5675E-01  9.9850E-02  3.9928E-01  1.0454E-01  3.0306E-01 -2.6199E-01  2.4910E-02  3.5808E-01
             7.2126E-02
 GRADIENT:  -7.4239E+00 -1.6517E+00  7.3650E+00 -1.1866E+01  7.0162E+00  2.6931E+00  5.2834E+00 -2.0634E+00  1.8557E+00 -7.8808E+00
            -8.4254E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2148.50038757287        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      525
 NPARAMETR:  9.7419E-01  1.0493E+00  1.5937E+00  1.0684E+00  1.2520E+00  9.9722E-01  1.3082E+00  6.0478E-01  8.5512E-01  1.2662E+00
             9.7808E-01
 PARAMETER:  7.3847E-02  1.4809E-01  5.6605E-01  1.6616E-01  3.2475E-01  9.7212E-02  3.6865E-01 -4.0288E-01 -5.6517E-02  3.3601E-01
             7.7835E-02
 GRADIENT:   3.3098E+00  1.2241E+00  2.6118E+00 -2.7788E+00 -9.7043E-01  3.6985E-01  3.7461E-01 -9.5440E-01 -1.4981E+00 -7.9461E-01
             1.1047E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2148.91399154956        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  9.6927E-01  8.9485E-01  2.0202E+00  1.1872E+00  1.2969E+00  9.9329E-01  1.3813E+00  1.0373E+00  8.5913E-01  1.3162E+00
             9.7047E-01
 PARAMETER:  6.8787E-02 -1.1097E-02  8.0320E-01  2.7156E-01  3.6001E-01  9.3263E-02  4.2303E-01  1.3658E-01 -5.1833E-02  3.7475E-01
             7.0030E-02
 GRADIENT:  -5.2778E+00  7.5274E+00  1.1664E+00  9.5533E+00 -2.1393E+00 -9.9653E-01  9.0761E-01 -4.4424E-01  3.5604E-01 -1.3378E-01
            -4.1238E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2149.18080925640        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  9.7053E-01  6.3147E-01  2.4814E+00  1.3693E+00  1.3024E+00  9.9465E-01  1.5200E+00  1.2389E+00  8.3406E-01  1.3399E+00
             9.7505E-01
 PARAMETER:  7.0089E-02 -3.5970E-01  1.0088E+00  4.1430E-01  3.6418E-01  9.4631E-02  5.1870E-01  3.1420E-01 -8.1448E-02  3.9257E-01
             7.4735E-02
 GRADIENT:   1.3198E+00  6.0724E+00  1.7645E-01  1.4135E+01  1.5917E-01  3.2013E-02 -1.7802E-01 -3.3046E-01 -2.4351E-01 -2.1001E-01
            -1.0200E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2149.24026440298        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  9.7000E-01  5.0081E-01  2.6839E+00  1.4577E+00  1.2939E+00  9.9461E-01  1.6053E+00  1.3144E+00  8.1850E-01  1.3427E+00
             9.7685E-01
 PARAMETER:  6.9538E-02 -5.9154E-01  1.0873E+00  4.7683E-01  3.5765E-01  9.4595E-02  5.7334E-01  3.7338E-01 -1.0029E-01  3.9466E-01
             7.6583E-02
 GRADIENT:   1.9190E+00  5.1297E+00  8.9965E-01  1.3843E+01 -1.6565E+00  2.5318E-01 -2.9434E-01 -4.2201E-01 -4.1010E-01 -2.7678E-01
             2.2758E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2149.27070402934        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1229
 NPARAMETR:  9.6904E-01  3.9519E-01  2.8410E+00  1.5284E+00  1.2875E+00  9.9400E-01  1.6870E+00  1.3783E+00  8.0477E-01  1.3458E+00
             9.7721E-01
 PARAMETER:  6.8552E-02 -8.2838E-01  1.1441E+00  5.2420E-01  3.5267E-01  9.3985E-02  6.2297E-01  4.2087E-01 -1.1719E-01  3.9698E-01
             7.6943E-02
 GRADIENT:   1.2737E+00  4.1608E+00  8.6968E-01  1.2311E+01 -2.3112E+00  2.0189E-01 -2.1017E-01 -2.1779E-01 -3.2349E-01 -1.4704E-01
             4.3447E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2149.29119908097        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1405
 NPARAMETR:  9.6804E-01  2.8589E-01  2.9943E+00  1.6018E+00  1.2799E+00  9.9339E-01  1.7851E+00  1.4367E+00  7.8877E-01  1.3478E+00
             9.7763E-01
 PARAMETER:  6.7519E-02 -1.1521E+00  1.1967E+00  5.7111E-01  3.4675E-01  9.3370E-02  6.7946E-01  4.6234E-01 -1.3728E-01  3.9850E-01
             7.7380E-02
 GRADIENT:   6.2045E-01  3.3016E+00  7.3228E-01  1.1021E+01 -2.8098E+00  1.6121E-01 -9.0471E-02 -1.8409E-02 -2.1550E-01 -4.7915E-03
             7.2660E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2149.40362876519        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1588
 NPARAMETR:  9.6753E-01  2.2462E-01  3.0530E+00  1.6318E+00  1.2779E+00  9.9287E-01  1.8440E+00  1.4570E+00  7.7882E-01  1.3493E+00
             9.7650E-01
 PARAMETER:  6.6995E-02 -1.3934E+00  1.2161E+00  5.8968E-01  3.4525E-01  9.2845E-02  7.1194E-01  4.7638E-01 -1.4997E-01  3.9959E-01
             7.6221E-02
 GRADIENT:   6.6334E-01  6.2552E-01  6.9772E-02 -1.1959E+01  9.0838E-01  1.3237E-01  5.4882E-03  1.3130E-01 -1.7115E-01 -1.4871E-01
             1.7400E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2149.42298730244        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1774             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6753E-01  1.9930E-01  3.0560E+00  1.6446E+00  1.2727E+00  9.9257E-01  1.5398E+00  1.4494E+00  7.7923E-01  1.3495E+00
             9.7630E-01
 PARAMETER:  6.6989E-02 -1.5129E+00  1.2171E+00  5.9752E-01  3.4114E-01  9.2541E-02  5.3168E-01  4.7112E-01 -1.4945E-01  3.9970E-01
             7.6011E-02
 GRADIENT:   4.2693E+02  2.7253E+01  5.7165E+00  1.1039E+03  1.7594E+01  3.6877E+01  1.5677E+00  5.3831E-01  1.7411E+01  5.5146E+00
             1.1062E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2149.45565621295        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1953
 NPARAMETR:  9.6627E-01  1.8003E-01  3.0637E+00  1.6645E+00  1.2676E+00  9.9210E-01  3.8374E-01  1.4466E+00  7.8241E-01  1.3483E+00
             9.7624E-01
 PARAMETER:  6.5692E-02 -1.6146E+00  1.2196E+00  6.0951E-01  3.3714E-01  9.2072E-02 -8.5780E-01  4.6924E-01 -1.4537E-01  3.9888E-01
             7.5956E-02
 GRADIENT:  -1.4483E+00  1.4074E+00 -4.1518E-01 -9.0305E-01  2.4450E-01 -6.7343E-02  1.1582E-02 -3.7145E-02 -7.3482E-01  1.1207E-01
            -2.3564E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2149.45579011586        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2135
 NPARAMETR:  9.6623E-01  1.7676E-01  3.0655E+00  1.6669E+00  1.2670E+00  9.9208E-01  3.2877E-01  1.4467E+00  7.8174E-01  1.3481E+00
             9.7623E-01
 PARAMETER:  6.5644E-02 -1.6330E+00  1.2202E+00  6.1095E-01  3.3664E-01  9.2044E-02 -1.0124E+00  4.6930E-01 -1.4623E-01  3.9869E-01
             7.5945E-02
 GRADIENT:  -1.4826E+00  1.4106E+00 -4.2198E-01 -5.3108E-01  2.1040E-01 -6.9348E-02  8.8821E-03 -4.0677E-02 -6.5085E-01  1.1249E-01
            -2.3794E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2149.49640742906        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2320
 NPARAMETR:  9.6710E-01  1.6664E-01  3.0707E+00  1.6666E+00  1.2662E+00  9.9223E-01  1.6816E-01  1.4476E+00  7.8181E-01  1.3477E+00
             9.7634E-01
 PARAMETER:  6.6544E-02 -1.6919E+00  1.2219E+00  6.1078E-01  3.3602E-01  9.2197E-02 -1.6829E+00  4.6991E-01 -1.4614E-01  3.9842E-01
             7.6053E-02
 GRADIENT:   8.3228E-01  2.1931E-01 -1.9391E-01 -1.5656E+01  1.1634E+00  5.8735E-02  3.4457E-03 -1.5766E-02  8.1795E-01  4.9582E-02
             9.5630E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2149.49917895691        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     2507
 NPARAMETR:  9.6743E-01  1.6516E-01  3.0677E+00  1.6663E+00  1.2645E+00  9.9229E-01  3.9623E-02  1.4487E+00  7.8044E-01  1.3464E+00
             9.7631E-01
 PARAMETER:  6.6893E-02 -1.7008E+00  1.2209E+00  6.1059E-01  3.3470E-01  9.2259E-02 -3.1283E+00  4.7064E-01 -1.4789E-01  3.9744E-01
             7.6030E-02
 GRADIENT:   1.6659E+00  9.3086E-02 -6.3431E-02 -1.8231E+01  7.1527E-01  9.7669E-02  3.0207E-04  7.5272E-02  3.2581E-01  2.6217E-02
             9.1933E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2149.50029326939        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2692
 NPARAMETR:  9.6741E-01  1.6488E-01  3.0628E+00  1.6670E+00  1.2635E+00  9.9226E-01  2.4753E-02  1.4463E+00  7.7989E-01  1.3458E+00
             9.7628E-01
 PARAMETER:  6.6866E-02 -1.7025E+00  1.2193E+00  6.1101E-01  3.3385E-01  9.2229E-02 -3.5988E+00  4.6902E-01 -1.4861E-01  3.9702E-01
             7.5994E-02
 GRADIENT:   1.6279E+00  2.0146E-01 -7.1897E-02 -1.6926E+01  5.0539E-01  9.1470E-02  1.3602E-04  7.8052E-02  8.5666E-02  4.0184E-02
             5.0425E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2149.50092920982        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     2877             RESET HESSIAN, TYPE I
 NPARAMETR:  9.6738E-01  1.6500E-01  3.0582E+00  1.6671E+00  1.2618E+00  9.9223E-01  1.8526E-02  1.4414E+00  7.7979E-01  1.3454E+00
             9.7621E-01
 PARAMETER:  6.6835E-02 -1.7018E+00  1.2178E+00  6.1110E-01  3.3257E-01  9.2195E-02 -3.8886E+00  4.6561E-01 -1.4873E-01  3.9667E-01
             7.5920E-02
 GRADIENT:   4.2721E+02  1.9762E+01  6.0598E+00  1.1570E+03  1.5706E+01  3.6808E+01  1.4934E-03  4.8771E-01  1.8125E+01  5.5016E+00
             9.5612E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2149.50126763661        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     3057
 NPARAMETR:  9.6737E-01  1.6490E-01  3.0538E+00  1.6671E+00  1.2611E+00  9.9221E-01  1.5878E-02  1.4392E+00  7.7976E-01  1.3450E+00
             9.7619E-01
 PARAMETER:  6.6823E-02 -1.7024E+00  1.2164E+00  6.1111E-01  3.3202E-01  9.2180E-02 -4.0428E+00  4.6407E-01 -1.4877E-01  3.9636E-01
             7.5907E-02
 GRADIENT:   1.5615E+00  2.8508E-01  1.5757E-01 -1.6071E+01 -1.6090E-01  8.2303E-02  6.4565E-05  1.4803E-03  5.3255E-03  1.0405E-01
            -4.6114E-02

0ITERATION NO.:   88    OBJECTIVE VALUE:  -2149.50141588322        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     3156
 NPARAMETR:  9.6738E-01  1.6436E-01  3.0494E+00  1.6669E+00  1.2613E+00  9.9221E-01  1.4645E-02  1.4392E+00  7.7986E-01  1.3445E+00
             9.7625E-01
 PARAMETER:  6.6832E-02 -1.7057E+00  1.2150E+00  6.1095E-01  3.3213E-01  9.2183E-02 -4.1237E+00  4.6408E-01 -1.4863E-01  3.9603E-01
             7.5961E-02
 GRADIENT:   4.5778E-02 -4.7802E-02 -1.6648E-02 -6.2255E-01  5.1758E-01  7.3244E-03  1.1681E-05  6.2151E-02  6.4878E-02 -8.2519E-03
             5.7144E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     3156
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.2665E-04 -9.5138E-05 -2.2715E-02 -4.8317E-03 -4.0430E-02
 SE:             2.9823E-02  4.0731E-05  1.1762E-02  2.9461E-02  2.2999E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8324E-01  1.9502E-02  5.3452E-02  8.6973E-01  7.8768E-02

 ETASHRINKSD(%)  8.9478E-02  9.9864E+01  6.0596E+01  1.3028E+00  2.2950E+01
 ETASHRINKVR(%)  1.7888E-01  1.0000E+02  8.4473E+01  2.5886E+00  4.0634E+01
 EBVSHRINKSD(%)  3.5915E-01  9.9872E+01  6.4493E+01  1.4341E+00  1.8522E+01
 EBVSHRINKVR(%)  7.1702E-01  1.0000E+02  8.7393E+01  2.8476E+00  3.3613E+01
 RELATIVEINF(%)  9.8404E+01  8.9104E-06  3.7027E+00  5.4541E+00  1.8476E+01
 EPSSHRINKSD(%)  3.2248E+01
 EPSSHRINKVR(%)  5.4097E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2149.5014158832182     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1230.5628826785455     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    52.70
 Elapsed covariance  time in seconds:     7.47
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2149.501       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.67E-01  1.64E-01  3.05E+00  1.67E+00  1.26E+00  9.92E-01  1.46E-02  1.44E+00  7.80E-01  1.34E+00  9.76E-01
 


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
 
         2.94E-02  1.05E+00  1.25E+00  6.94E-01  3.41E-01  7.72E-02  6.49E-02  6.89E-01  3.29E-01  2.41E-01  4.63E-02
 


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
+        8.64E-04
 
 TH 2
+        7.28E-03  1.11E+00
 
 TH 3
+       -1.63E-03 -3.13E-02  1.56E+00
 
 TH 4
+       -4.92E-03 -7.30E-01  6.18E-02  4.82E-01
 
 TH 5
+        1.60E-03  2.91E-01  2.30E-01 -1.85E-01  1.16E-01
 
 TH 6
+       -5.54E-04  2.42E-02  2.58E-02 -1.53E-02  9.74E-03  5.97E-03
 
 TH 7
+       -3.59E-05  6.07E-02 -2.44E-02 -4.03E-02  1.19E-02  1.47E-03  4.22E-03
 
 TH 8
+        2.01E-04  2.68E-01  7.02E-01 -1.57E-01  1.80E-01  1.84E-02  5.20E-03  4.75E-01
 
 TH 9
+        2.07E-03  3.43E-01 -2.93E-02 -2.25E-01  8.67E-02  7.10E-03  1.94E-02  7.23E-02  1.08E-01
 
 TH10
+       -3.53E-04  1.38E-01  1.47E-01 -8.69E-02  6.12E-02  8.19E-03  3.91E-03  9.90E-02  4.03E-02  5.82E-02
 
 TH11
+        4.50E-04  1.98E-02 -8.34E-03 -1.33E-02  4.08E-03 -1.34E-05  7.24E-04 -2.39E-04  6.40E-03  5.29E-04  2.15E-03
 
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
+        2.94E-02
 
 TH 2
+        2.35E-01  1.05E+00
 
 TH 3
+       -4.43E-02 -2.37E-02  1.25E+00
 
 TH 4
+       -2.41E-01 -9.97E-01  7.12E-02  6.94E-01
 
 TH 5
+        1.60E-01  8.10E-01  5.39E-01 -7.80E-01  3.41E-01
 
 TH 6
+       -2.44E-01  2.97E-01  2.67E-01 -2.85E-01  3.70E-01  7.72E-02
 
 TH 7
+       -1.88E-02  8.85E-01 -3.01E-01 -8.94E-01  5.38E-01  2.93E-01  6.49E-02
 
 TH 8
+        9.91E-03  3.68E-01  8.14E-01 -3.28E-01  7.65E-01  3.46E-01  1.16E-01  6.89E-01
 
 TH 9
+        2.14E-01  9.87E-01 -7.13E-02 -9.87E-01  7.72E-01  2.79E-01  9.06E-01  3.19E-01  3.29E-01
 
 TH10
+       -4.97E-02  5.42E-01  4.86E-01 -5.19E-01  7.44E-01  4.39E-01  2.50E-01  5.96E-01  5.08E-01  2.41E-01
 
 TH11
+        3.30E-01  4.06E-01 -1.44E-01 -4.12E-01  2.58E-01 -3.73E-03  2.41E-01 -7.47E-03  4.20E-01  4.74E-02  4.63E-02
 
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
+        2.49E+05
 
 TH 2
+       -2.53E+04  2.95E+03
 
 TH 3
+       -1.86E+03  1.90E+02  2.29E+01
 
 TH 4
+       -1.64E+03  6.28E+02 -2.09E+00  6.89E+02
 
 TH 5
+        3.92E+04 -4.09E+03 -3.29E+02 -3.11E+02  6.47E+03
 
 TH 6
+       -2.75E+04  2.81E+03  1.95E+02  2.14E+02 -4.31E+03  3.33E+03
 
 TH 7
+        4.90E+05 -5.01E+04 -3.64E+03 -3.49E+03  7.74E+04 -5.47E+04  9.69E+05
 
 TH 8
+       -9.02E+02  9.02E+01  2.29E+00  7.09E+00 -1.51E+02  1.01E+02 -1.81E+03  1.52E+01
 
 TH 9
+       -6.69E+04  6.72E+03  4.99E+02  4.79E+02 -1.06E+04  7.51E+03 -1.33E+05  2.54E+02  1.86E+04
 
 TH10
+        4.03E+04 -4.10E+03 -3.00E+02 -2.66E+02  6.33E+03 -4.52E+03  7.97E+04 -1.46E+02 -1.09E+04  6.59E+03
 
 TH11
+        1.14E+05 -1.16E+04 -8.45E+02 -8.01E+02  1.80E+04 -1.27E+04  2.26E+05 -4.13E+02 -3.10E+04  1.86E+04  5.32E+04
 
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
 #CPUT: Total CPU Time in Seconds,       60.230
Stop Time:
Wed Sep 29 21:23:53 CDT 2021
