Thu Sep 30 08:41:29 CDT 2021
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
$DATA ../../../../data/spa2/D/dat19.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m19.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   925.602854208906        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.3440E+02 -3.5029E+01 -7.5141E+01 -2.3873E+02  2.2304E+02 -5.6509E+02 -2.9638E+02 -2.4073E+01 -7.6924E+02 -2.6891E+02
            -4.9262E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1426.00067176348        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      128
 NPARAMETR:  2.1665E+00  9.9157E-01  1.4395E+00  1.8212E+00  1.1939E+00  4.0160E+00  2.1804E+00  1.0024E+00  4.1321E+00  2.1788E+00
             2.9218E+00
 PARAMETER:  8.7314E-01  9.1532E-02  4.6429E-01  6.9950E-01  2.7722E-01  1.4903E+00  8.7952E-01  1.0236E-01  1.5188E+00  8.7878E-01
             1.1722E+00
 GRADIENT:   9.9405E+01 -3.9228E+01 -7.7549E+01  5.3086E+01  4.1155E+01  8.0285E+01  1.6504E+01  1.6269E+00  9.5557E+01  2.6094E+01
            -6.6786E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1485.18714420228        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  3.3725E+00  1.0886E+00  1.3372E+01  1.4987E+00  2.2422E+00  5.0560E+00  4.5076E+00  8.3043E-01  2.1817E+00  1.5019E+00
             3.0863E+00
 PARAMETER:  1.3157E+00  1.8493E-01  2.6931E+00  5.0461E-01  9.0746E-01  1.7206E+00  1.6058E+00 -8.5817E-02  8.8009E-01  5.0671E-01
             1.2270E+00
 GRADIENT:   1.2595E+02  1.7513E+01 -3.6721E+00  2.9061E+01  6.9255E+01  9.3863E+00  2.0623E+01 -3.2029E-02  2.9743E+01 -1.2671E+00
            -2.8371E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1587.10636567977        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      482
 NPARAMETR:  9.1806E-01  6.1766E-01  1.0793E+01  1.8880E+00  1.5975E+00  2.7247E+00  5.8751E+00  2.9449E-01  2.7799E+00  1.5605E+00
             3.1931E+00
 PARAMETER:  1.4512E-02 -3.8182E-01  2.4789E+00  7.3550E-01  5.6847E-01  1.1024E+00  1.8707E+00 -1.1225E+00  1.1224E+00  5.4500E-01
             1.2610E+00
 GRADIENT:  -7.5069E+01  1.0775E+01  4.2962E+00  2.7255E+01 -6.7363E+01  7.9714E+01  2.2714E+01  4.6791E-02  6.6354E+01  2.0818E+01
             3.7735E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1639.95418060626        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      658
 NPARAMETR:  1.1690E+00  5.1717E-01  1.8350E+01  1.5306E+00  2.0080E+00  2.0813E+00  5.3347E+00  1.6387E+00  1.5706E+00  1.3145E+00
             3.1381E+00
 PARAMETER:  2.5611E-01 -5.5939E-01  3.0096E+00  5.2564E-01  7.9715E-01  8.3299E-01  1.7742E+00  5.9387E-01  5.5149E-01  3.7346E-01
             1.2436E+00
 GRADIENT:  -5.8818E+00 -5.9154E+00  8.7821E-01 -1.2406E+01  2.6371E+01  2.2994E+00 -2.9463E+00 -1.7673E-01  1.9165E+01 -3.1068E+00
            -3.2122E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1651.28433618198        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      835
 NPARAMETR:  1.1882E+00  5.8893E-01  8.2808E+00  1.4187E+00  1.7325E+00  2.0917E+00  5.3534E+00  3.3868E+00  1.2103E+00  1.1284E+00
             3.1260E+00
 PARAMETER:  2.7244E-01 -4.2945E-01  2.2139E+00  4.4975E-01  6.4955E-01  8.3797E-01  1.7777E+00  1.3199E+00  2.9090E-01  2.2082E-01
             1.2398E+00
 GRADIENT:   3.5139E+00 -2.8234E+00  7.9977E+00 -1.3366E+00 -4.7079E+00  4.6625E+00 -2.3005E+00 -1.4321E+00  7.1849E-02  5.2947E+00
             1.4140E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1658.35678153481        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1010
 NPARAMETR:  1.1813E+00  9.4164E-01  2.4617E+00  1.1632E+00  1.4241E+00  2.0715E+00  4.3096E+00  1.5935E+00  1.0831E+00  8.8879E-01
             3.0855E+00
 PARAMETER:  2.6661E-01  3.9864E-02  1.0009E+00  2.5114E-01  4.5353E-01  8.2828E-01  1.5608E+00  5.6596E-01  1.7980E-01 -1.7896E-02
             1.2267E+00
 GRADIENT:   1.0059E+00 -1.2271E+00  1.3437E-01  3.0777E-01 -2.4721E+00  8.6037E-01  3.5742E-01 -3.2085E-01 -3.4109E-01 -1.6195E-01
            -4.8950E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1658.55625690016        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1187
 NPARAMETR:  1.1841E+00  1.0273E+00  2.4405E+00  1.1244E+00  1.4503E+00  2.1037E+00  4.1282E+00  1.6373E+00  1.0884E+00  8.9595E-01
             3.0874E+00
 PARAMETER:  2.6900E-01  1.2695E-01  9.9220E-01  2.1721E-01  4.7174E-01  8.4369E-01  1.5178E+00  5.9306E-01  1.8467E-01 -9.8667E-03
             1.2273E+00
 GRADIENT:   2.0374E+00  1.1713E-01 -4.2236E-02 -3.5340E-01  2.3779E-01  7.0997E+00  4.1225E-01  3.3741E-02  1.6382E-01 -5.6357E-02
            -1.4817E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1658.56370160832        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1370             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1827E+00  1.0319E+00  2.4371E+00  1.1226E+00  1.4503E+00  2.1015E+00  4.0873E+00  1.6327E+00  1.0825E+00  8.9497E-01
             3.0869E+00
 PARAMETER:  2.6783E-01  1.3142E-01  9.9082E-01  2.1569E-01  4.7175E-01  8.4263E-01  1.5079E+00  5.9021E-01  1.7926E-01 -1.0961E-02
             1.2272E+00
 GRADIENT:   1.3729E+02  8.7238E+00  2.5677E+00  2.9710E+01  6.4797E+00  1.7720E+02  1.7858E+02  1.9615E-01  1.9862E+00 -4.2292E-02
             1.4149E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1658.56749086258        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1550             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1827E+00  1.0354E+00  2.4317E+00  1.1208E+00  1.4507E+00  2.1014E+00  4.0857E+00  1.6316E+00  1.0877E+00  8.9562E-01
             3.0872E+00
 PARAMETER:  2.6784E-01  1.3481E-01  9.8861E-01  2.1405E-01  4.7206E-01  8.4259E-01  1.5075E+00  5.8955E-01  1.8404E-01 -1.0241E-02
             1.2273E+00
 GRADIENT:   1.3727E+02  9.0134E+00  2.4113E+00  2.8944E+01  6.7672E+00  1.7715E+02  1.7872E+02  2.5325E-01  2.4435E+00  2.0772E-02
             1.4428E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1658.56896616282        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1734
 NPARAMETR:  1.1828E+00  1.0371E+00  2.4273E+00  1.1195E+00  1.4509E+00  2.1014E+00  4.0802E+00  1.6320E+00  1.0885E+00  8.9631E-01
             3.0873E+00
 PARAMETER:  2.6785E-01  1.3647E-01  9.8679E-01  2.1288E-01  4.7218E-01  8.4259E-01  1.5061E+00  5.8978E-01  1.8480E-01 -9.4708E-03
             1.2273E+00
 GRADIENT:   1.4329E+00 -3.6934E-02  1.8334E-01  3.5918E-01 -1.4437E-01  6.7012E+00 -1.1881E+00 -3.3325E-02 -4.9260E-02 -3.9144E-02
            -1.7861E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1658.57049849173        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1920             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1828E+00  1.0397E+00  2.4188E+00  1.1173E+00  1.4512E+00  2.1013E+00  4.0718E+00  1.6331E+00  1.0896E+00  8.9781E-01
             3.0877E+00
 PARAMETER:  2.6786E-01  1.3898E-01  9.8328E-01  2.1095E-01  4.7241E-01  8.4258E-01  1.5041E+00  5.9048E-01  1.8583E-01 -7.7974E-03
             1.2274E+00
 GRADIENT:   1.3722E+02  9.1109E+00  2.0730E+00  2.8041E+01  7.1880E+00  1.7709E+02  1.7803E+02  4.2517E-01  2.6736E+00  1.7863E-01
             1.4854E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1658.57095305345        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2104
 NPARAMETR:  1.1828E+00  1.0418E+00  2.4139E+00  1.1162E+00  1.4513E+00  2.1013E+00  4.0673E+00  1.6328E+00  1.0905E+00  8.9851E-01
             3.0879E+00
 PARAMETER:  2.6786E-01  1.4091E-01  9.8125E-01  2.0992E-01  4.7249E-01  8.4257E-01  1.5030E+00  5.9029E-01  1.8665E-01 -7.0168E-03
             1.2275E+00
 GRADIENT:   1.4156E+00 -2.0676E-01 -1.6757E-01  1.2212E-01  3.1062E-01  6.6853E+00 -1.1359E+00  1.3463E-01  1.6795E-01  1.1656E-01
             2.5952E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1658.57216333913        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2290             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1828E+00  1.0450E+00  2.4121E+00  1.1151E+00  1.4509E+00  2.1014E+00  4.0618E+00  1.6270E+00  1.0884E+00  8.9587E-01
             3.0874E+00
 PARAMETER:  2.6789E-01  1.4400E-01  9.8050E-01  2.0891E-01  4.7215E-01  8.4258E-01  1.5016E+00  5.8671E-01  1.8469E-01 -9.9581E-03
             1.2273E+00
 GRADIENT:   1.3728E+02  9.6641E+00  2.3114E+00  2.7639E+01  6.8459E+00  1.7714E+02  1.7711E+02  3.0002E-01  2.5161E+00  7.2641E-02
             1.4535E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1658.57253585937        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2474
 NPARAMETR:  1.1828E+00  1.0464E+00  2.4093E+00  1.1142E+00  1.4507E+00  2.1014E+00  4.0585E+00  1.6252E+00  1.0883E+00  8.9564E-01
             3.0874E+00
 PARAMETER:  2.6790E-01  1.4535E-01  9.7932E-01  2.0816E-01  4.7204E-01  8.4258E-01  1.5008E+00  5.8564E-01  1.8459E-01 -1.0218E-02
             1.2273E+00
 GRADIENT:   1.4422E+00  4.0446E-03  2.3612E-01  3.4436E-02 -2.7206E-01  6.7041E+00 -1.2703E+00 -5.1962E-02 -5.0635E-02 -4.5245E-02
            -2.0991E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1658.57309852890        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     2660             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1828E+00  1.0475E+00  2.4039E+00  1.1134E+00  1.4508E+00  2.1013E+00  4.0554E+00  1.6256E+00  1.0889E+00  8.9635E-01
             3.0876E+00
 PARAMETER:  2.6789E-01  1.4641E-01  9.7709E-01  2.0745E-01  4.7210E-01  8.4257E-01  1.5000E+00  5.8585E-01  1.8512E-01 -9.4201E-03
             1.2274E+00
 GRADIENT:   1.3725E+02  9.8014E+00  2.1501E+00  2.7332E+01  7.0642E+00  1.7711E+02  1.7673E+02  3.5733E-01  2.5786E+00  1.1873E-01
             1.4679E+01

0ITERATION NO.:   77    OBJECTIVE VALUE:  -1658.57309852890        NO. OF FUNC. EVALS.:  60
 CUMULATIVE NO. OF FUNC. EVALS.:     2720
 NPARAMETR:  1.1828E+00  1.0475E+00  2.4039E+00  1.1134E+00  1.4508E+00  2.1013E+00  4.0554E+00  1.6256E+00  1.0889E+00  8.9635E-01
             3.0876E+00
 PARAMETER:  2.6789E-01  1.4641E-01  9.7709E-01  2.0745E-01  4.7210E-01  8.4257E-01  1.5000E+00  5.8585E-01  1.8512E-01 -9.4201E-03
             1.2274E+00
 GRADIENT:  -3.0849E-03 -4.6415E-02  4.5266E-02  1.8305E-01  8.0534E-02 -8.6402E-04  8.2900E-02  2.2315E-02  1.5145E-02  7.6373E-03
             1.7259E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2720
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.8488E-03  2.3714E-02 -2.8946E-02 -4.6076E-02 -1.0967E-02
 SE:             2.9744E-02  2.5296E-02  1.1542E-02  1.6526E-02  1.6567E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5044E-01  3.4853E-01  1.2145E-02  5.3032E-03  5.0797E-01

 ETASHRINKSD(%)  3.5356E-01  1.5254E+01  6.1333E+01  4.4634E+01  4.4499E+01
 ETASHRINKVR(%)  7.0588E-01  2.8181E+01  8.5049E+01  6.9347E+01  6.9197E+01
 EBVSHRINKSD(%)  6.1302E-01  1.1322E+01  6.5736E+01  4.7461E+01  4.4307E+01
 EBVSHRINKVR(%)  1.2223E+00  2.1362E+01  8.8260E+01  7.2397E+01  6.8983E+01
 RELATIVEINF(%)  9.8713E+01  3.0120E+01  3.0251E+00  8.8033E+00  9.4318E+00
 EPSSHRINKSD(%)  2.1925E+01
 EPSSHRINKVR(%)  3.9043E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1658.5730985288972     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -555.84685868329007     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    65.19
 Elapsed covariance  time in seconds:    10.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1658.573       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.18E+00  1.05E+00  2.40E+00  1.11E+00  1.45E+00  2.10E+00  4.06E+00  1.63E+00  1.09E+00  8.96E-01  3.09E+00
 


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
 
         7.95E-02  2.88E-01  6.88E-01  1.29E-01  1.29E-01  1.58E-01  6.12E-01  5.11E-01  1.67E-01  6.68E-01  1.40E+00
 


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
+        6.32E-03
 
 TH 2
+        1.27E-02  8.32E-02
 
 TH 3
+        4.98E-03 -6.10E-02  4.74E-01
 
 TH 4
+       -6.20E-04 -2.87E-02  4.94E-02  1.67E-02
 
 TH 5
+        2.89E-03  1.03E-02  5.63E-02 -1.02E-03  1.67E-02
 
 TH 6
+        2.53E-03 -2.04E-03  3.31E-02  2.07E-03  6.31E-03  2.51E-02
 
 TH 7
+       -9.13E-03 -1.39E-01  2.06E-01  5.73E-02  2.95E-03  4.74E-02  3.75E-01
 
 TH 8
+        1.16E-02  1.74E-02  2.37E-01  1.01E-02  3.54E-02  7.61E-03  1.09E-02  2.61E-01
 
 TH 9
+        2.11E-04 -1.75E-03 -4.19E-03  2.00E-03 -1.66E-03  8.74E-04  1.34E-03 -2.66E-02  2.79E-02
 
 TH10
+       -1.46E-02 -6.87E-02 -5.20E-02  9.19E-03 -1.19E-02  3.08E-02  1.79E-01 -1.95E-01  3.96E-02  4.46E-01
 
 TH11
+        3.21E-02  1.39E-01  1.75E-01 -1.33E-02  3.99E-02 -6.10E-02 -3.47E-01  4.23E-01 -9.00E-02 -8.86E-01  1.95E+00
 
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
+        7.95E-02
 
 TH 2
+        5.56E-01  2.88E-01
 
 TH 3
+        9.10E-02 -3.07E-01  6.88E-01
 
 TH 4
+       -6.03E-02 -7.70E-01  5.55E-01  1.29E-01
 
 TH 5
+        2.81E-01  2.75E-01  6.34E-01 -6.13E-02  1.29E-01
 
 TH 6
+        2.01E-01 -4.46E-02  3.04E-01  1.01E-01  3.09E-01  1.58E-01
 
 TH 7
+       -1.88E-01 -7.88E-01  4.89E-01  7.23E-01  3.73E-02  4.89E-01  6.12E-01
 
 TH 8
+        2.85E-01  1.18E-01  6.75E-01  1.54E-01  5.37E-01  9.41E-02  3.49E-02  5.11E-01
 
 TH 9
+        1.59E-02 -3.62E-02 -3.64E-02  9.27E-02 -7.69E-02  3.30E-02  1.31E-02 -3.11E-01  1.67E-01
 
 TH10
+       -2.74E-01 -3.57E-01 -1.13E-01  1.06E-01 -1.38E-01  2.91E-01  4.39E-01 -5.73E-01  3.54E-01  6.68E-01
 
 TH11
+        2.89E-01  3.44E-01  1.82E-01 -7.35E-02  2.21E-01 -2.76E-01 -4.06E-01  5.93E-01 -3.86E-01 -9.49E-01  1.40E+00
 
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
+        5.71E+02
 
 TH 2
+       -2.47E+02  1.82E+02
 
 TH 3
+        1.50E+01  8.03E-01  1.51E+01
 
 TH 4
+       -2.99E+02  1.63E+02 -3.43E+01  4.34E+02
 
 TH 5
+        1.55E+01 -4.98E+01 -3.57E+01  6.50E+01  2.09E+02
 
 TH 6
+        2.86E+00 -4.18E+01 -3.83E+00  2.52E+01  2.29E+00  9.34E+01
 
 TH 7
+       -4.07E+01  4.14E+01 -9.74E-01 -4.60E+00 -9.41E+00 -2.85E+01  2.38E+01
 
 TH 8
+       -1.10E+01 -1.39E+00 -9.16E+00  1.52E+01  4.26E+00  2.00E+00 -1.18E+00  1.58E+01
 
 TH 9
+       -5.88E+00 -5.62E+00 -5.13E+00 -2.70E+01  3.93E+00 -2.15E+00  7.24E+00  6.04E+00  5.22E+01
 
 TH10
+       -3.55E+00  8.94E-01 -1.89E+00  3.23E+00 -1.35E+01  2.52E+00 -2.03E+00  5.65E+00  2.42E+00  2.69E+01
 
 TH11
+       -2.21E+00 -3.48E-01 -5.61E-01 -5.06E+00 -5.85E+00  1.85E+00  9.98E-01  3.93E-01  4.11E+00  1.13E+01  6.17E+00
 
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
 #CPUT: Total CPU Time in Seconds,       76.154
Stop Time:
Thu Sep 30 08:42:47 CDT 2021
