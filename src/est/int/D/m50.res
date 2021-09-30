Wed Sep 29 09:04:10 CDT 2021
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
$DATA ../../../../data/int/D/dat50.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   24698.3666081154        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7303E+02  3.0521E+02  7.8859E+00  1.9754E+02  1.8032E+02 -1.4913E+03 -9.4380E+02 -6.2251E+01 -1.2809E+03 -5.3731E+02
            -5.2349E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -983.356302592071        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2706E+00  1.8099E+00  8.7072E-01  2.0759E+00  8.4970E-01  4.2884E+00  4.8725E+00  9.8547E-01  2.6701E+00  2.0058E+00
             1.2699E+01
 PARAMETER:  3.3947E-01  6.9327E-01 -3.8436E-02  8.3038E-01 -6.2868E-02  1.5559E+00  1.6836E+00  8.5365E-02  1.0821E+00  7.9604E-01
             2.6415E+00
 GRADIENT:  -5.7129E+00  2.6407E+01 -4.0395E+01  9.9527E+01 -2.1260E+01  1.6555E+02  7.3080E+01  4.5857E+00  6.0684E+01  5.2001E+01
             5.6813E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1075.20378274431        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.2574E+00  1.5235E+00  2.9149E+01  3.5226E+00  2.4643E+00  2.6489E+00  1.0679E+01  6.3848E-01  2.7846E+00  2.8649E+00
             1.2468E+01
 PARAMETER:  3.2906E-01  5.2101E-01  3.4724E+00  1.3592E+00  1.0019E+00  1.0742E+00  2.4683E+00 -3.4866E-01  1.1241E+00  1.1525E+00
             2.6232E+00
 GRADIENT:  -5.6918E+00  2.8922E+01 -5.2957E+00  1.3082E+02 -4.1711E+00  9.5765E+01  8.0629E+01  5.5241E-02  3.4858E+01  8.6826E+01
             5.4585E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1337.44272083658        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  1.0358E+00  1.3679E+00  1.0845E+01  7.7164E-01  2.6798E+00  1.6540E+00  3.4842E+00  4.6522E-01  1.0082E+00  6.4293E-01
             8.0428E+00
 PARAMETER:  1.3516E-01  4.1328E-01  2.4837E+00 -1.5924E-01  1.0857E+00  6.0318E-01  1.3482E+00 -6.6524E-01  1.0813E-01 -3.4171E-01
             2.1848E+00
 GRADIENT:  -3.8604E+01 -5.9637E+01 -3.2927E+00 -5.0165E+01  7.4872E+01 -2.2641E+01 -1.4389E-01  8.3911E-03  9.7605E+00  7.0713E+00
             1.4270E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1365.29124646550        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  1.1635E+00  1.9431E+00  3.3175E+00  7.1520E-01  2.1246E+00  2.0181E+00  3.2055E+00  2.5303E-01  9.6391E-01  4.6260E-01
             7.4590E+00
 PARAMETER:  2.5146E-01  7.6426E-01  1.2992E+00 -2.3520E-01  8.5358E-01  8.0215E-01  1.2649E+00 -1.2742E+00  6.3247E-02 -6.7089E-01
             2.1094E+00
 GRADIENT:   2.0522E+01  1.2842E+01 -2.3927E+00  1.9763E+00  3.1104E+01 -1.1858E+01  1.3408E+00 -4.5837E-03  7.3832E+00  4.1216E+00
             1.4402E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1366.24955671274        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:      437             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1591E+00  1.9457E+00  3.2575E+00  7.1601E-01  2.1211E+00  2.0087E+00  3.2658E+00  2.6261E-01  9.4080E-01  4.4848E-01
             7.6153E+00
 PARAMETER:  2.4766E-01  7.6564E-01  1.2810E+00 -2.3406E-01  8.5192E-01  7.9747E-01  1.2835E+00 -1.2371E+00  3.8980E-02 -7.0188E-01
             2.1302E+00
 GRADIENT:   1.4966E+01  1.1680E+01 -2.7879E+00 -1.5959E+00  3.1596E+01 -1.4644E+01  8.2498E+00 -1.3221E-03  7.5004E+00  3.9953E+00
             5.9232E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1376.27605306560        NO. OF FUNC. EVALS.: 112
 CUMULATIVE NO. OF FUNC. EVALS.:      549
 NPARAMETR:  1.1264E+00  1.7762E+00  3.1744E+00  7.6608E-01  1.9918E+00  2.1044E+00  3.4667E+00  2.6277E-01  2.7103E-01  2.0233E-01
             7.4321E+00
 PARAMETER:  2.1905E-01  6.7449E-01  1.2551E+00 -1.6647E-01  7.8904E-01  8.4404E-01  1.3432E+00 -1.2365E+00 -1.2055E+00 -1.4978E+00
             2.1058E+00
 GRADIENT:   4.8834E+00  1.1441E+01  1.5646E+00 -1.5127E+00  9.5936E+00  8.7075E+00  1.2055E+01  7.0663E-03  9.4722E-01  7.0641E-01
             6.9415E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1381.75024713834        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      688
 NPARAMETR:  1.1468E+00  1.6566E+00  2.9260E+00  8.4850E-01  1.9023E+00  2.1344E+00  4.0204E+00  2.5068E-01  2.1358E-01  1.5954E-01
             7.5789E+00
 PARAMETER:  2.3702E-01  6.0478E-01  1.1737E+00 -6.4288E-02  7.4305E-01  8.5817E-01  1.4914E+00 -1.2836E+00 -1.4438E+00 -1.7355E+00
             2.1254E+00
 GRADIENT:  -8.1835E+00 -1.1322E+00 -3.5049E+00  1.1367E+00  1.0004E+01 -1.5760E+01 -1.2254E+01  2.2612E-02  3.5685E-01  4.5124E-01
             4.8051E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1382.30900318609        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      829
 NPARAMETR:  1.1485E+00  1.6255E+00  3.2527E+00  8.5562E-01  1.9033E+00  2.1309E+00  4.1108E+00  8.1656E-02  1.1014E-01  5.2753E-02
             7.6177E+00
 PARAMETER:  2.3847E-01  5.8582E-01  1.2795E+00 -5.5927E-02  7.4358E-01  8.5653E-01  1.5136E+00 -2.4052E+00 -2.1060E+00 -2.8421E+00
             2.1305E+00
 GRADIENT:   1.2849E+01  1.2917E+01  2.1528E+00 -1.4284E+00 -1.0993E+00  1.3454E+01  5.5272E+01  1.5999E-03  5.7488E-02  4.7424E-02
             5.2382E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1382.64479889534        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1011             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1483E+00  1.6308E+00  3.1524E+00  8.5690E-01  1.9060E+00  2.1385E+00  4.2055E+00  4.7461E-02  8.5324E-02  2.3077E-02
             7.5785E+00
 PARAMETER:  2.3830E-01  5.8909E-01  1.2482E+00 -5.4438E-02  7.4499E-01  8.6011E-01  1.5364E+00 -2.9478E+00 -2.3613E+00 -3.6689E+00
             2.1253E+00
 GRADIENT:   1.3498E+01  1.4129E+01  1.8085E-02 -1.2140E+00  5.5937E+00  1.4957E+01  6.3928E+01  8.1635E-04  4.1679E-02  9.7201E-03
             4.1912E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1382.76585451430        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1190
 NPARAMETR:  1.1485E+00  1.5858E+00  3.2641E+00  8.6808E-01  1.9056E+00  2.1480E+00  4.2387E+00  1.0000E-02  6.1525E-02  1.0000E-02
             7.5658E+00
 PARAMETER:  2.3848E-01  5.6106E-01  1.2830E+00 -4.1471E-02  7.4482E-01  8.6456E-01  1.5443E+00 -4.6815E+00 -2.6883E+00 -5.2484E+00
             2.1236E+00
 GRADIENT:  -6.6229E+00 -1.1349E+00  5.2536E-01 -3.4807E-01 -1.6028E+00 -1.3265E+01 -4.6573E+00  0.0000E+00 -1.3863E-02  0.0000E+00
            -1.4433E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1382.80770865656        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     1377
 NPARAMETR:  1.1478E+00  1.6071E+00  3.2382E+00  8.6775E-01  1.9097E+00  2.1477E+00  4.2686E+00  1.0000E-02  7.1749E-02  1.0000E-02
             7.5817E+00
 PARAMETER:  2.3824E-01  5.7790E-01  1.2747E+00 -4.0846E-02  7.4581E-01  8.6549E-01  1.5487E+00 -4.6815E+00 -2.5503E+00 -5.2484E+00
             2.1237E+00
 GRADIENT:   2.8725E+02  6.3180E-01 -4.9937E-02  1.0506E+00 -9.6938E+01  7.0398E+01 -2.6159E+01  0.0000E+00 -1.6682E-02  0.0000E+00
            -3.5920E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1377
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7606E-02  5.2151E-03 -2.3507E-05 -5.7859E-03  2.4512E-05
 SE:             3.0012E-02  2.8237E-02  4.5207E-05  1.4242E-03  1.6158E-04
 N:                     100         100         100         100         100

 P VAL.:         5.5745E-01  8.5347E-01  6.0307E-01  4.8550E-05  8.7943E-01

 ETASHRINKSD(%)  1.0000E-10  5.4026E+00  9.9849E+01  9.5229E+01  9.9459E+01
 ETASHRINKVR(%)  1.0000E-10  1.0513E+01  1.0000E+02  9.9772E+01  9.9997E+01
 EBVSHRINKSD(%)  3.3469E+00  3.7494E+00  9.9852E+01  9.6790E+01  9.9416E+01
 EBVSHRINKVR(%)  6.5818E+00  7.3582E+00  1.0000E+02  9.9897E+01  9.9997E+01
 RELATIVEINF(%)  9.3249E+01  5.2501E+01  5.0244E-05  5.0361E-02  8.4208E-04
 EPSSHRINKSD(%)  6.6666E+00
 EPSSHRINKVR(%)  1.2889E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1382.8077086565568     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       271.28165111185399     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    41.34
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    16.07
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1382.808       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.15E+00  1.61E+00  3.24E+00  8.69E-01  1.91E+00  2.15E+00  4.26E+00  1.00E-02  7.06E-02  1.00E-02  7.57E+00
 


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
+        2.39E+04
 
 TH 2
+       -4.97E+01  2.18E+01
 
 TH 3
+        7.33E+00  8.73E-01  2.56E+00
 
 TH 4
+       -5.42E+01  3.93E+01 -8.42E+00  4.30E+02
 
 TH 5
+       -3.33E+01 -4.57E+00 -1.41E+01  4.44E+01  1.00E+03
 
 TH 6
+       -3.30E+01 -3.45E-01  9.02E-02  2.16E+00 -2.07E+00  5.61E+02
 
 TH 7
+       -1.40E+01  2.22E+00 -7.40E-01 -2.90E+01  1.15E+01 -1.04E+00  5.08E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.13E-02 -4.37E-01 -3.38E-01 -1.58E+01  1.57E+00 -9.53E-02  7.74E-01  0.00E+00 -1.83E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.84E-01 -2.88E+00 -7.22E-02 -2.15E+01 -1.02E-01  3.22E+00  1.17E+00  0.00E+00  8.59E-01  0.00E+00  2.63E+01
 
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
 #CPUT: Total CPU Time in Seconds,       57.542
Stop Time:
Wed Sep 29 09:05:09 CDT 2021
