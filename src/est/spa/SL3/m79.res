Sat Sep 18 13:03:37 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat79.csv ignore=@
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
Current Date:       18 SEP 2021
Days until program expires : 211
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1592.67420153708        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4764E+00 -6.2612E+01 -2.0871E+01 -3.9870E+01  3.7104E+01  2.4234E+01 -1.6700E+01 -1.4889E-01  9.7195E+00 -6.9416E+00
            -1.7710E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1622.77898959308        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0167E+00  1.1327E+00  1.0520E+00  9.5026E-01  1.0702E+00  9.1168E-01  1.1325E+00  9.9443E-01  9.1124E-01  9.8182E-01
             1.3402E+00
 PARAMETER:  1.1653E-01  2.2465E-01  1.5074E-01  4.8976E-02  1.6784E-01  7.5385E-03  2.2439E-01  9.4414E-02  7.0562E-03  8.1657E-02
             3.9282E-01
 GRADIENT:   1.9065E+01 -1.3552E+01 -2.0566E+00 -8.0030E+00  1.9511E+01 -9.8219E+00 -4.0009E+00 -2.6972E+00  1.4334E+00 -4.3617E+00
            -6.1583E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1623.19333862014        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0169E+00  1.1381E+00  1.0241E+00  9.4911E-01  1.0629E+00  9.3126E-01  1.1924E+00  1.0458E+00  8.7509E-01  9.7329E-01
             1.3505E+00
 PARAMETER:  1.1671E-01  2.2939E-01  1.2384E-01  4.7775E-02  1.6103E-01  2.8785E-02  2.7593E-01  1.4481E-01 -3.3427E-02  7.2931E-02
             4.0044E-01
 GRADIENT:   1.8633E+01 -8.2133E+00 -7.7657E+00 -3.2204E+00  2.2183E+01 -1.4363E+00  2.0807E+00 -5.6150E-01 -4.9750E-01 -2.2802E+00
            -1.4045E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1623.67511451318        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      300
 NPARAMETR:  1.0183E+00  8.8348E-01  1.3856E+00  1.1398E+00  1.0670E+00  9.3774E-01  1.4408E+00  1.3776E+00  7.7317E-01  1.0209E+00
             1.3553E+00
 PARAMETER:  1.1811E-01 -2.3882E-02  4.2611E-01  2.3084E-01  1.6487E-01  3.5720E-02  4.6517E-01  4.2035E-01 -1.5725E-01  1.2069E-01
             4.0402E-01
 GRADIENT:   1.2522E+00  9.3491E+00 -1.1506E+00  1.7488E+01 -3.0698E+00  7.0835E-01 -7.6282E-01  5.4225E-01 -1.0342E+00  1.0938E+00
            -2.6376E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1623.89932077769        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      476
 NPARAMETR:  1.0175E+00  7.7646E-01  1.7212E+00  1.2090E+00  1.1233E+00  9.3405E-01  1.5930E+00  1.6420E+00  7.3520E-01  1.0609E+00
             1.3609E+00
 PARAMETER:  1.1736E-01 -1.5302E-01  6.4303E-01  2.8977E-01  2.1628E-01  3.1774E-02  5.6563E-01  5.9592E-01 -2.0761E-01  1.5909E-01
             4.0813E-01
 GRADIENT:   3.0854E+00  4.8703E+00  2.8886E+00  4.2594E+00 -2.5597E+00  2.6262E-02  2.6840E-01 -1.3197E+00 -8.8114E-01 -2.4774E-01
            -6.0729E-03

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1625.60394436329        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      655
 NPARAMETR:  1.0094E+00  1.9877E-01  2.9934E+00  1.6025E+00  1.1885E+00  9.2295E-01  3.7137E+00  2.4954E+00  6.4664E-01  1.1372E+00
             1.3557E+00
 PARAMETER:  1.0937E-01 -1.5156E+00  1.1964E+00  5.7159E-01  2.7265E-01  1.9820E-02  1.4120E+00  1.0145E+00 -3.3596E-01  2.2861E-01
             4.0434E-01
 GRADIENT:  -3.4249E+00  4.3018E+00 -1.4228E+00  2.8967E+01  5.6643E-01 -7.3683E-01  3.2227E+00 -6.1835E-01  1.5581E+00  5.5452E-01
            -1.9872E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1626.08044419794        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      831
 NPARAMETR:  1.0108E+00  1.2059E-01  3.7138E+00  1.6525E+00  1.2412E+00  9.2378E-01  4.3962E+00  2.9379E+00  6.3826E-01  1.1565E+00
             1.3646E+00
 PARAMETER:  1.1073E-01 -2.0154E+00  1.4121E+00  6.0232E-01  3.1612E-01  2.0719E-02  1.5807E+00  1.1777E+00 -3.4901E-01  2.4536E-01
             4.1089E-01
 GRADIENT:  -4.4272E-02  5.5737E-01 -1.1553E+00  8.5444E+00  1.1747E+00  3.8516E-02 -8.1314E-01  8.0194E-01  9.4671E-01 -2.7192E-01
            -6.4719E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1626.52698362830        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1008
 NPARAMETR:  1.0121E+00  6.1869E-02  4.9704E+00  1.6953E+00  1.3137E+00  9.2053E-01  5.8507E+00  3.4983E+00  6.1995E-01  1.2182E+00
             1.3730E+00
 PARAMETER:  1.1205E-01 -2.6827E+00  1.7035E+00  6.2786E-01  3.7286E-01  1.7199E-02  1.8666E+00  1.3523E+00 -3.7812E-01  2.9739E-01
             4.1701E-01
 GRADIENT:   8.2324E-02 -8.3152E-02 -4.3511E-01 -4.6955E+00  1.2933E+00 -2.1085E-01  2.0538E-02  2.3558E-01  9.7975E-01  2.9900E-02
             1.0839E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1626.53700041248        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1183
 NPARAMETR:  1.0123E+00  5.8961E-02  5.1581E+00  1.6996E+00  1.3212E+00  9.2079E-01  5.9770E+00  3.5704E+00  6.1650E-01  1.2240E+00
             1.3710E+00
 PARAMETER:  1.1223E-01 -2.7309E+00  1.7406E+00  6.3041E-01  3.7856E-01  1.7476E-02  1.8879E+00  1.3727E+00 -3.8370E-01  3.0216E-01
             4.1555E-01
 GRADIENT:   7.7845E-02  1.0202E-01 -3.0523E-01 -1.1545E+00  5.3250E-01 -2.2488E-02  4.4743E-01  2.0290E-01 -2.1174E-01  4.0893E-02
             1.0717E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1626.61324277859        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1367
 NPARAMETR:  1.0134E+00  5.3470E-02  6.0402E+00  1.7079E+00  1.3605E+00  9.2170E-01  6.0992E+00  3.9506E+00  6.0487E-01  1.2463E+00
             1.3757E+00
 PARAMETER:  1.1328E-01 -2.8286E+00  1.8984E+00  6.3525E-01  4.0786E-01  1.8467E-02  1.9082E+00  1.4739E+00 -4.0274E-01  3.2015E-01
             4.1896E-01
 GRADIENT:   2.7123E-01  9.2680E-01 -8.7101E-01  1.1808E+00  2.5299E+00  8.6661E-01 -4.5341E-01  1.0995E+00 -4.0068E+00  1.1651E-01
             6.1592E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1628.32372816412        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1546
 NPARAMETR:  1.0108E+00  1.0000E-02  4.9510E+00  1.7171E+00  1.3066E+00  9.1891E-01  1.2154E+01  3.4221E+00  5.9327E-01  1.2025E+00
             1.3731E+00
 PARAMETER:  1.1078E-01 -4.6293E+00  1.6996E+00  6.4066E-01  3.6747E-01  1.5428E-02  2.5976E+00  1.3303E+00 -4.2211E-01  2.8443E-01
             4.1709E-01
 GRADIENT:  -1.7669E+00  0.0000E+00  1.3267E+00  3.4903E+01  4.1969E-01 -7.0121E-01 -9.1604E+00 -1.9773E+00  1.0670E+01  2.1024E+00
             1.3421E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1628.81884263289        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1724
 NPARAMETR:  1.0107E+00  1.0000E-02  4.7615E+00  1.7136E+00  1.2975E+00  9.1911E-01  1.2691E+01  3.3425E+00  5.9031E-01  1.1959E+00
             1.3727E+00
 PARAMETER:  1.1059E-01 -4.7642E+00  1.6606E+00  6.3861E-01  3.6043E-01  1.5648E-02  2.6409E+00  1.3067E+00 -4.2710E-01  2.7891E-01
             4.1677E-01
 GRADIENT:  -1.6407E+00  0.0000E+00  1.2183E+00  2.7146E+01  1.2221E+00 -6.2587E-01 -5.5531E+00 -2.1895E+00  7.7747E+00  1.6310E+00
             1.3918E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1628.90267349221        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1901
 NPARAMETR:  1.0108E+00  1.0000E-02  4.4575E+00  1.7053E+00  1.2785E+00  9.2020E-01  1.2451E+01  3.2548E+00  5.8491E-01  1.1821E+00
             1.3688E+00
 PARAMETER:  1.1077E-01 -4.7661E+00  1.5946E+00  6.3372E-01  3.4568E-01  1.6838E-02  2.6218E+00  1.2801E+00 -4.3629E-01  2.6731E-01
             4.1396E-01
 GRADIENT:  -1.7799E-01  0.0000E+00 -6.1461E-01 -3.9481E+01  4.0906E+00  3.0777E-01  2.1430E+01 -1.5367E+00 -1.6792E+01 -3.2019E+00
             2.1030E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1628.90673075238        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2087
 NPARAMETR:  1.0106E+00  1.0000E-02  4.3414E+00  1.7039E+00  1.2696E+00  9.2028E-01  1.2425E+01  3.2187E+00  5.8568E-01  1.1741E+00
             1.3668E+00
 PARAMETER:  1.1056E-01 -4.7655E+00  1.5682E+00  6.3293E-01  3.3869E-01  1.6924E-02  2.6197E+00  1.2690E+00 -4.3498E-01  2.6054E-01
             4.1246E-01
 GRADIENT:  -1.7974E-02  0.0000E+00 -2.2421E-01 -1.0885E+01  8.8901E-01  5.8299E-02  5.5553E+00 -2.4196E-01 -4.4629E+00 -8.9174E-01
            -8.1943E-02

0ITERATION NO.:   68    OBJECTIVE VALUE:  -1628.90676143494        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     2179
 NPARAMETR:  1.0106E+00  1.0000E-02  4.3414E+00  1.7039E+00  1.2696E+00  9.2036E-01  1.2425E+01  3.2187E+00  5.8568E-01  1.1741E+00
             1.3668E+00
 PARAMETER:  1.1056E-01 -4.7655E+00  1.5682E+00  6.3293E-01  3.3869E-01  1.7010E-02  2.6197E+00  1.2690E+00 -4.3498E-01  2.6054E-01
             4.1246E-01
 GRADIENT:   8.9719E-02  0.0000E+00 -1.7366E-02 -1.6708E+00  1.5837E-01  3.2996E-02  4.4134E-01 -2.2152E-02 -5.0499E-01 -6.2558E-02
             2.7686E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2179
 NO. OF SIG. DIGITS IN FINAL EST.:  3.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.2586E-03  2.4816E-02 -4.4333E-02 -2.2665E-02 -5.6747E-02
 SE:             2.9633E-02  9.8484E-03  1.7818E-02  2.6524E-02  1.9102E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6612E-01  1.1743E-02  1.2845E-02  3.9282E-01  2.9712E-03

 ETASHRINKSD(%)  7.2683E-01  6.7007E+01  4.0306E+01  1.1141E+01  3.6005E+01
 ETASHRINKVR(%)  1.4484E+00  8.9114E+01  6.4366E+01  2.1040E+01  5.9047E+01
 EBVSHRINKSD(%)  9.2731E-01  7.4361E+01  4.8965E+01  9.2846E+00  3.1656E+01
 EBVSHRINKVR(%)  1.8460E+00  9.3426E+01  7.3954E+01  1.7707E+01  5.3291E+01
 RELATIVEINF(%)  9.7831E+01  5.9050E+00  1.2819E+01  6.5942E+01  2.3174E+01
 EPSSHRINKSD(%)  4.2143E+01
 EPSSHRINKVR(%)  6.6526E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1628.9067614349365     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -893.75593487119829     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    32.10
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1628.907       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.00E-02  4.34E+00  1.70E+00  1.27E+00  9.20E-01  1.24E+01  3.22E+00  5.86E-01  1.17E+00  1.37E+00
 


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
+        1.50E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -6.83E+00  0.00E+00  4.95E+03
 
 TH 4
+       -5.12E+01  0.00E+00 -1.66E+01  1.97E+05
 
 TH 5
+        4.50E+01  0.00E+00 -7.84E+04  7.52E+01  1.24E+06
 
 TH 6
+       -1.96E+02  0.00E+00  7.52E+00  1.47E+01 -8.50E+01  3.85E+02
 
 TH 7
+        5.61E-01  0.00E+00  3.46E-01  3.41E+01 -4.01E+00 -1.00E+00  2.16E+02
 
 TH 8
+       -7.48E+00  0.00E+00  8.25E+03 -1.28E+01 -1.31E+05  1.09E+01  3.69E-01  1.38E+04
 
 TH 9
+       -6.99E+01  0.00E+00 -5.56E+01 -9.64E+02  6.87E+02  9.07E+01  2.56E+01 -6.48E+01  3.53E+06
 
 TH10
+       -9.54E+01  0.00E+00  1.10E+05 -1.47E+02 -1.75E+06  1.39E+02  4.36E+00  1.84E+05 -7.40E+02  2.45E+06
 
 TH11
+        3.04E+01  0.00E+00  1.74E+00  7.89E+00 -5.14E+01  7.71E+00 -1.09E+00  5.37E+00  1.33E+02  5.76E+01  1.46E+02
 
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
 #CPUT: Total CPU Time in Seconds,       40.538
Stop Time:
Sat Sep 18 13:04:19 CDT 2021
