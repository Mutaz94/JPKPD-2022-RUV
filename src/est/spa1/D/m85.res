Thu Sep 30 03:41:31 CDT 2021
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
$DATA ../../../../data/spa1/D/dat85.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m85.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   25212.7272047493        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.5231E+02  5.6883E+02 -2.1458E+01  4.6692E+02  1.7973E+02 -2.6487E+03 -1.2509E+03 -4.4615E+01 -1.9556E+03 -6.0453E+02
            -4.7435E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -480.367161775138        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1607E+00  1.0495E+00  9.5240E-01  1.4522E+00  1.2676E+00  1.6649E+00  1.1533E+00  9.5711E-01  1.1621E+00  1.0041E+00
             1.4840E+01
 PARAMETER:  2.4898E-01  1.4831E-01  5.1229E-02  4.7305E-01  3.3716E-01  6.0975E-01  2.4262E-01  5.6161E-02  2.5023E-01  1.0411E-01
             2.7974E+00
 GRADIENT:  -4.2432E+01  2.0516E+01 -5.3630E+00  3.0333E+01 -7.2677E+00  2.0842E+01 -1.1519E+00  4.3670E+00  9.0335E-01  2.8119E+00
             7.9339E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -489.153019840277        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.1849E+00  8.9567E-01  1.2481E+00  1.6088E+00  3.1817E+00  1.5093E+00  2.6056E+00  3.8252E-01  1.1933E+00  3.3202E+00
             1.3921E+01
 PARAMETER:  2.6966E-01 -1.0186E-02  3.2159E-01  5.7547E-01  1.2574E+00  5.1162E-01  1.0577E+00 -8.6098E-01  2.7669E-01  1.3000E+00
             2.7334E+00
 GRADIENT:  -4.4838E+00  2.4932E+01 -2.1916E+00  5.1932E+01 -5.4024E+00 -2.0970E+01  4.5738E+00  1.6487E-01  4.9039E+00 -1.6550E-01
             1.9713E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -499.318874119245        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.1357E+00  5.9121E-01  9.0557E-01  1.4490E+00  4.1333E+00  1.7160E+00  1.5950E+00  2.4368E-02  9.9437E-01  6.8173E+00
             1.3373E+01
 PARAMETER:  2.2728E-01 -4.2559E-01  8.0793E-04  4.7087E-01  1.5191E+00  6.4000E-01  5.6689E-01 -3.6145E+00  9.4359E-02  2.0195E+00
             2.6932E+00
 GRADIENT:   2.1699E+01  7.0160E+00  2.2278E+00 -3.3109E+01 -1.0960E+01  3.5739E+01  1.9462E+00 -2.3330E-04 -6.1696E+00  3.9752E+00
             2.8387E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -532.614079543211        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      304
 NPARAMETR:  1.0022E+00  1.2386E-01  3.9783E-01  1.6214E+00  1.7379E+01  1.8084E+00  1.6956E-01  1.0000E-02  1.5724E+00  1.5278E+01
             1.2212E+01
 PARAMETER:  1.0221E-01 -1.9886E+00 -8.2172E-01  5.8328E-01  2.9553E+00  6.9244E-01 -1.6745E+00 -9.5898E+00  5.5258E-01  2.8264E+00
             2.6024E+00
 GRADIENT:   1.9932E+00  8.4863E-01  3.2522E+00  7.7245E+01 -5.2687E-01 -2.5378E+01  2.6421E-03  0.0000E+00  2.7175E+01  1.4384E+01
             2.1905E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -544.798634752451        NO. OF FUNC. EVALS.: 126
 CUMULATIVE NO. OF FUNC. EVALS.:      430
 NPARAMETR:  9.9655E-01  1.2441E-01  3.9232E-01  1.4045E+00  1.7401E+01  2.0276E+00  5.2247E-02  1.0000E-02  1.5545E+00  1.4824E+01
             1.2336E+01
 PARAMETER:  9.6541E-02 -1.9842E+00 -8.3567E-01  4.3970E-01  2.9565E+00  8.0685E-01 -2.8518E+00 -9.6137E+00  5.4114E-01  2.7962E+00
             2.6125E+00
 GRADIENT:   2.9730E+00 -1.9704E+00  4.0730E+01 -2.0801E+01 -3.6060E-01  2.2827E+01  3.0695E-04  0.0000E+00  3.0287E+01  6.9083E+00
             6.5624E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -545.341760431999        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  9.9196E-01  1.2298E-01  3.9266E-01  1.4046E+00  1.7736E+01  1.9335E+00  2.8400E-02  1.0000E-02  1.5575E+00  1.4381E+01
             1.2145E+01
 PARAMETER:  9.1930E-02 -1.9957E+00 -8.3481E-01  4.3977E-01  2.9756E+00  7.5934E-01 -3.4614E+00 -9.6137E+00  5.4311E-01  2.7659E+00
             2.5969E+00
 GRADIENT:   2.8625E+00 -1.9502E+00  4.1006E+01 -1.6578E+01 -8.8076E-02  5.1625E+00  9.5600E-05  0.0000E+00  3.1474E+01  6.6488E+00
             4.9782E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -547.619961448793        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      621
 NPARAMETR:  9.5640E-01  1.1928E-01  3.8974E-01  1.3967E+00  1.8755E+01  1.9356E+00  1.0000E-02  1.0000E-02  1.5621E+00  1.2534E+01
             1.1209E+01
 PARAMETER:  5.5426E-02 -2.0263E+00 -8.4229E-01  4.3413E-01  3.0314E+00  7.6041E-01 -5.2612E+00 -9.6137E+00  5.4606E-01  2.6284E+00
             2.5168E+00
 GRADIENT:  -4.2765E-01 -1.1118E+00  4.1957E+01  2.3130E-01  1.0256E-01  7.1850E-01  0.0000E+00  0.0000E+00  2.4637E+01  3.8227E+00
            -1.7354E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -548.787117761984        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      783
 NPARAMETR:  9.7307E-01  1.1916E-01  3.8979E-01  1.3967E+00  1.8777E+01  1.9721E+00  1.0000E-02  1.0000E-02  1.5628E+00  1.2724E+01
             1.1613E+01
 PARAMETER:  7.2701E-02 -2.0273E+00 -8.4214E-01  4.3412E-01  3.0326E+00  7.7911E-01 -4.9471E+00 -9.6137E+00  5.4647E-01  2.6435E+00
             2.5521E+00
 GRADIENT:   2.3568E+00 -1.4481E+00  4.2322E+01 -8.0134E+00  6.4380E-02  1.0118E+01  0.0000E+00  0.0000E+00  2.6660E+01  3.7219E+00
             1.4699E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -549.626666989819        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      856
 NPARAMETR:  9.3933E-01  1.1913E-01  3.6442E-01  1.3957E+00  1.9107E+01  1.9390E+00  1.0000E-02  1.0000E-02  1.5384E+00  1.1271E+01
             1.1025E+01
 PARAMETER:  3.7414E-02 -2.0276E+00 -9.0945E-01  4.3340E-01  3.0500E+00  7.6218E-01 -4.9471E+00 -9.6137E+00  5.3075E-01  2.5222E+00
             2.5002E+00
 GRADIENT:  -2.2609E+00 -3.3093E-01  3.8586E+01  1.9236E+01  1.2936E-01  2.0747E+00  0.0000E+00  0.0000E+00  2.2307E+01  6.8345E-03
            -3.5244E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -559.524260918141        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      930
 NPARAMETR:  8.7271E-01  1.2000E-01  2.1225E-01  1.3556E+00  2.1263E+01  1.7351E+00  1.0000E-02  1.0000E-02  1.3487E+00  4.4900E+00
             1.1820E+01
 PARAMETER: -3.6156E-02 -2.0203E+00 -1.4500E+00  4.0428E-01  3.1570E+00  6.5105E-01 -4.9471E+00 -9.6137E+00  3.9915E-01  1.6018E+00
             2.5698E+00
 GRADIENT:  -9.2174E+00 -2.4537E+00 -1.6318E+01  1.5557E+02  5.6108E-01 -2.3523E+01  0.0000E+00  0.0000E+00  4.1980E+00  5.1096E-02
            -1.4261E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -570.502649178516        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1003
 NPARAMETR:  7.6262E-01  1.2481E-01  1.6383E-01  1.1651E+00  2.1685E+01  1.8207E+00  1.0000E-02  1.0000E-02  1.2348E+00  2.7596E+00
             1.1691E+01
 PARAMETER: -1.7100E-01 -1.9810E+00 -1.7089E+00  2.5281E-01  3.1766E+00  6.9923E-01 -4.9471E+00 -9.6137E+00  3.1091E-01  1.1151E+00
             2.5588E+00
 GRADIENT:  -5.1183E+01  1.7529E+01 -2.0436E+01  1.2981E+02 -6.9920E-01 -1.5052E+00  0.0000E+00  0.0000E+00  5.4237E+00  1.8019E-01
             5.1144E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -594.494987728259        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     1079
 NPARAMETR:  7.9513E-01  6.4096E-02  8.6003E-02  7.4611E-01  2.7258E+01  1.7731E+00  1.0000E-02  1.0000E-02  9.6326E-01  7.3482E-01
             1.1591E+01
 PARAMETER: -1.2925E-01 -2.6474E+00 -2.3534E+00 -1.9288E-01  3.4054E+00  6.7275E-01 -4.9471E+00 -9.6137E+00  6.2566E-02 -2.0812E-01
             2.5502E+00
 GRADIENT:   9.5998E+01  9.1537E+00 -3.2121E+01  7.0916E+01 -1.9266E-02  1.0186E+01  0.0000E+00  0.0000E+00 -1.0666E+00  5.6902E-03
            -1.6678E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -602.180403349779        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1175
 NPARAMETR:  6.2291E-01  2.7493E-02  4.4205E-02  4.8062E-01  3.6108E+01  1.6188E+00  1.0000E-02  1.0000E-02  7.4688E-01  1.8625E-01
             1.1364E+01
 PARAMETER: -3.7335E-01 -3.4938E+00 -3.0189E+00 -6.3268E-01  3.6865E+00  5.8171E-01 -4.9471E+00 -9.6137E+00 -1.9185E-01 -1.5807E+00
             2.5305E+00
 GRADIENT:   1.0073E+02 -2.6652E+00 -1.1184E+02  1.2518E+02  2.9572E-01 -1.4194E+01  0.0000E+00  0.0000E+00 -2.1053E+01  1.2685E-05
            -6.6158E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -617.320539997834        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1353
 NPARAMETR:  3.6367E-01  1.0000E-02  1.4413E-02  1.9853E-01  5.9752E+01  1.4566E+00  1.0000E-02  1.0000E-02  4.7585E-01  1.6630E-02
             1.2036E+01
 PARAMETER: -9.1151E-01 -5.1547E+00 -4.1397E+00 -1.5168E+00  4.1902E+00  4.7613E-01 -4.9471E+00 -9.6137E+00 -6.4265E-01 -3.9966E+00
             2.5879E+00
 GRADIENT:   2.2149E+01  0.0000E+00 -4.8006E+01  4.3958E+01  2.5953E-02 -1.4235E+00  0.0000E+00  0.0000E+00  2.5462E+00  4.0194E-09
            -6.6270E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -617.464850680218        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1481
 NPARAMETR:  3.7869E-01  1.0000E-02  1.7174E-02  2.2409E-01  5.2705E+01  1.4757E+00  1.0000E-02  1.0000E-02  5.0374E-01  2.3986E-02
             1.1874E+01
 PARAMETER: -8.7105E-01 -4.9372E+00 -3.9643E+00 -1.3957E+00  4.0647E+00  4.8914E-01 -4.9471E+00 -9.6137E+00 -5.8569E-01 -3.6303E+00
             2.5744E+00
 GRADIENT:   5.5165E+01  0.0000E+00  7.8213E+01  3.6460E+01  6.8341E-03  4.6831E+00  0.0000E+00  0.0000E+00  8.5793E-01  8.2111E-09
             2.0636E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -617.563700339088        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1657
 NPARAMETR:  3.7474E-01  1.0000E-02  1.6802E-02  2.2079E-01  3.3254E+01  1.4720E+00  1.0000E-02  1.0000E-02  4.1733E-01  2.0420E-02
             1.1853E+01
 PARAMETER: -8.8153E-01 -4.9372E+00 -3.9862E+00 -1.4105E+00  3.6042E+00  4.8659E-01 -4.9471E+00 -9.6137E+00 -7.7388E-01 -3.7912E+00
             2.5725E+00
 GRADIENT:   5.9258E+01  0.0000E+00  8.5168E+01  2.9749E+01  1.1604E-02  6.0833E+00  0.0000E+00  0.0000E+00  3.1355E-01  1.8073E-08
             1.0499E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -618.028964469553        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1730
 NPARAMETR:  3.7586E-01  1.0000E-02  1.6123E-02  2.1644E-01  1.9300E+01  1.4441E+00  1.0000E-02  1.0000E-02  1.0431E-01  1.8153E-02
             1.2673E+01
 PARAMETER: -8.7855E-01 -4.9372E+00 -4.0275E+00 -1.4304E+00  3.0601E+00  4.6747E-01 -4.9471E+00 -9.6137E+00 -2.1604E+00 -3.9089E+00
             2.6395E+00
 GRADIENT:   6.8767E+01  0.0000E+00  9.0596E+01  1.4868E+01  2.8790E-02  8.5108E+00  0.0000E+00  0.0000E+00  1.4963E-01  7.3248E-08
             4.5902E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -618.056980420741        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1802
 NPARAMETR:  3.7175E-01  1.0000E-02  1.5371E-02  2.0999E-01  1.0801E+01  1.4287E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.1056E-02
             1.2798E+01
 PARAMETER: -8.8954E-01 -4.9372E+00 -4.0753E+00 -1.4607E+00  2.4797E+00  4.5677E-01 -4.9471E+00 -9.6137E+00 -4.7164E+00 -4.4047E+00
             2.6493E+00
 GRADIENT:   7.3750E+01  0.0000E+00  8.0296E+01  2.7178E+01  1.2853E-02  5.9583E+00  0.0000E+00  0.0000E+00  0.0000E+00  1.8431E-06
             5.0908E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -618.277665739058        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:     1965
 NPARAMETR:  3.7478E-01  1.0000E-02  1.6154E-02  2.1823E-01  1.1142E+01  1.4229E+00  1.0000E-02  1.0000E-02  8.4526E-02  1.5416E-02
             1.2544E+01
 PARAMETER: -8.8141E-01 -4.9372E+00 -4.0256E+00 -1.4222E+00  2.5107E+00  4.5269E-01 -4.9471E+00 -9.6137E+00 -2.3707E+00 -4.0723E+00
             2.6292E+00
 GRADIENT:   9.5764E+00  0.0000E+00 -6.1904E+00  3.8811E+00  3.7014E-02 -2.5576E+00  0.0000E+00  0.0000E+00  8.9856E-02  1.6623E-06
             1.0054E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -618.598742959537        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2140
 NPARAMETR:  3.6566E-01  1.0000E-02  1.5791E-02  2.1396E-01  1.3513E+01  1.4327E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.2352E+01
 PARAMETER: -9.0606E-01 -4.9372E+00 -4.0483E+00 -1.4420E+00  2.7037E+00  4.5953E-01 -4.9471E+00 -9.6137E+00 -4.5215E+00 -5.2679E+00
             2.6138E+00
 GRADIENT:   1.9426E-01  0.0000E+00 -8.5832E-01  1.7672E+00  6.3558E-03 -5.8922E-01  0.0000E+00  0.0000E+00  5.0254E-04  0.0000E+00
             1.0777E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -618.822207194962        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2319
 NPARAMETR:  3.5773E-01  1.0000E-02  1.4979E-02  2.0383E-01  1.6030E+01  1.4442E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.2313E+01
 PARAMETER: -9.2799E-01 -4.9372E+00 -4.1011E+00 -1.4905E+00  2.8745E+00  4.6759E-01 -4.9471E+00 -9.6137E+00 -7.8816E+00 -6.7540E+00
             2.6106E+00
 GRADIENT:   9.2406E-01  0.0000E+00  8.9429E+00 -1.1944E+01 -2.5231E-02  3.2799E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -4.9548E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -618.856790871198        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2507
 NPARAMETR:  3.5834E-01  1.0000E-02  1.4987E-02  2.0514E-01  1.8869E+01  1.4323E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.2320E+01
 PARAMETER: -9.2626E-01 -4.9372E+00 -4.1006E+00 -1.4840E+00  3.0375E+00  4.5931E-01 -4.9471E+00 -9.6137E+00 -7.6144E+00 -6.5552E+00
             2.6112E+00
 GRADIENT:   2.6237E+00  0.0000E+00 -3.6649E+00  3.0774E+00  9.5736E-03  9.5564E-02  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -1.9232E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -618.871853062253        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2684
 NPARAMETR:  3.5755E-01  1.0000E-02  1.4947E-02  2.0465E-01  1.6254E+01  1.4313E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.2345E+01
 PARAMETER: -9.2849E-01 -4.9372E+00 -4.1033E+00 -1.4865E+00  2.8883E+00  4.5860E-01 -4.9471E+00 -9.6137E+00 -7.6144E+00 -6.5552E+00
             2.6132E+00
 GRADIENT:   5.0205E-01  0.0000E+00 -1.8069E+00  1.6935E+00 -4.7562E-03  5.6985E-02  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
             5.3761E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -618.904261995408        NO. OF FUNC. EVALS.: 201
 CUMULATIVE NO. OF FUNC. EVALS.:     2885
 NPARAMETR:  3.5679E-01  1.0000E-02  1.4780E-02  2.0342E-01  1.6337E+01  1.4312E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.2337E+01
 PARAMETER: -9.3060E-01 -4.9372E+00 -4.1145E+00 -1.4925E+00  2.8934E+00  4.5850E-01 -4.9471E+00 -9.6137E+00 -7.6144E+00 -6.5552E+00
             2.6126E+00
 GRADIENT:   3.0933E+00  0.0000E+00 -8.6446E+00  8.6288E+00  9.9178E-03 -4.8247E-02  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -1.0303E+00

0ITERATION NO.:  125    OBJECTIVE VALUE:  -618.927940306003        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     3073
 NPARAMETR:  3.5604E-01  1.0000E-02  1.4767E-02  2.0266E-01  1.9226E+01  1.4306E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.2334E+01
 PARAMETER: -9.3272E-01 -4.9372E+00 -4.1154E+00 -1.4962E+00  3.0563E+00  4.5809E-01 -4.9471E+00 -9.6137E+00 -7.6144E+00 -6.5552E+00
             2.6124E+00
 GRADIENT:   1.8136E+00  0.0000E+00 -3.0487E+00  2.4008E+00  2.9044E-03  4.7076E-02  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -4.4371E-01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -618.946064443600        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     3203
 NPARAMETR:  3.5283E-01  1.0000E-02  1.4630E-02  2.0159E-01  1.5451E+01  1.4290E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.2303E+01
 PARAMETER: -9.4177E-01 -4.9372E+00 -4.1247E+00 -1.5015E+00  2.8377E+00  4.5696E-01 -4.9471E+00 -9.6137E+00 -7.6144E+00 -6.5552E+00
             2.6098E+00
 GRADIENT:  -3.2974E+00  0.0000E+00 -4.9415E+00  7.0529E+00 -1.9356E-02 -5.4756E-01  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -1.2888E+00

0ITERATION NO.:  135    OBJECTIVE VALUE:  -619.043866338939        NO. OF FUNC. EVALS.: 154
 CUMULATIVE NO. OF FUNC. EVALS.:     3357
 NPARAMETR:  3.4649E-01  1.0000E-02  1.4071E-02  1.9560E-01  1.2399E+02  1.4291E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.2190E+01
 PARAMETER: -9.5991E-01 -4.9372E+00 -4.1636E+00 -1.5317E+00  4.9202E+00  4.5708E-01 -4.9471E+00 -9.6137E+00 -7.6144E+00 -6.5552E+00
             2.6006E+00
 GRADIENT:  -4.4943E-01  0.0000E+00 -1.2124E+01  1.4115E+01  2.6922E-03 -5.3222E-01  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -1.0013E+01

0ITERATION NO.:  139    OBJECTIVE VALUE:  -619.130142057709        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     3497
 NPARAMETR:  3.4813E-01  1.0000E-02  1.4095E-02  1.9495E-01  1.1792E+02  1.4283E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.2303E+01
 PARAMETER: -9.5517E-01 -4.9372E+00 -4.1620E+00 -1.5350E+00  4.8700E+00  4.5648E-01 -4.9471E+00 -9.6137E+00 -7.6144E+00 -6.5552E+00
             2.6098E+00
 GRADIENT:   8.1608E-01  0.0000E+00 -3.0609E+00  2.2171E+00  1.3995E-03  2.3849E-01  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -1.6043E+00

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     3497
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.9232E-03  8.8419E-06  4.7145E-05 -1.6917E-04 -2.0848E-08
 SE:             2.8753E-02  4.1252E-06  3.2916E-04  3.6498E-04  1.0172E-07
 N:                     100         100         100         100         100

 P VAL.:         9.4667E-01  3.2080E-02  8.8611E-01  6.4300E-01  8.3760E-01

 ETASHRINKSD(%)  3.6749E+00  9.9986E+01  9.8897E+01  9.8777E+01  1.0000E+02
 ETASHRINKVR(%)  7.2147E+00  1.0000E+02  9.9988E+01  9.9985E+01  1.0000E+02
 EBVSHRINKSD(%)  3.8495E+00  9.9977E+01  9.8947E+01  9.8774E+01  1.0000E+02
 EBVSHRINKVR(%)  7.5509E+00  1.0000E+02  9.9989E+01  9.9985E+01  1.0000E+02
 RELATIVEINF(%)  8.0406E-01  2.0600E-07  2.5588E-05  2.8492E-05  0.0000E+00
 EPSSHRINKSD(%)  4.8614E+00
 EPSSHRINKVR(%)  9.4866E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -619.13014205770867     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       299.80839114696403     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    55.59
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     7.31
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -619.130       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.48E-01  1.00E-02  1.41E-02  1.95E-01  1.18E+02  1.43E+00  1.00E-02  1.00E-02  1.00E-02  1.00E-02  1.23E+01
 


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
+        2.59E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -2.64E+04  0.00E+00  2.69E+06
 
 TH 4
+        2.11E+03  0.00E+00 -2.15E+05  1.72E+04
 
 TH 5
+        7.66E-04  0.00E+00 -7.80E-02  6.25E-03  2.26E-09
 
 TH 6
+       -7.02E+00  0.00E+00  7.15E+02 -5.73E+01 -2.08E-05  1.90E-01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.87E+00  0.00E+00  3.94E+02 -3.16E+01 -1.14E-05  1.05E-01  0.00E+00  0.00E+00  0.00E+00  0.00E+00  5.78E-02
 
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
+        4.11E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -4.16E+04  0.00E+00  4.21E+06
 
 TH 4
+        6.73E+01  0.00E+00 -3.37E+05  3.01E+04
 
 TH 5
+        4.56E-03  0.00E+00 -1.22E-01  7.65E-03 -2.82E-08
 
 TH 6
+        9.66E+00  0.00E+00  1.12E+03 -1.46E+02  8.99E-06  8.07E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.38E+01  0.00E+00  6.19E+02 -3.01E+01 -4.92E-05  1.76E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.16E+00
 
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
+        4.15E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -6.81E+04  0.00E+00  6.61E+06
 
 TH 4
+        1.96E+03  0.00E+00 -5.07E+05  4.23E+04
 
 TH 5
+        4.95E-03  0.00E+00 -1.81E-01  1.13E-02  8.45E-09
 
 TH 6
+        4.28E+02  0.00E+00 -8.29E+03  2.36E+02  5.28E-04  1.38E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.12E+02  0.00E+00  2.18E+03 -9.89E+01 -1.54E-04  1.43E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.16E+01
 
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
 #CPUT: Total CPU Time in Seconds,       62.962
Stop Time:
Thu Sep 30 03:42:35 CDT 2021
