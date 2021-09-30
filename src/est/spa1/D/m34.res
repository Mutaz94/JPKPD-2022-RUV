Thu Sep 30 02:56:24 CDT 2021
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
$DATA ../../../../data/spa1/D/dat34.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   16654.5809471497        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0066E+02  1.8186E+02 -5.2066E+01  4.0213E+01  1.1766E+02 -1.4457E+03 -5.8459E+02 -6.0072E+01 -9.5552E+02 -4.2053E+02
            -3.3737E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -685.832116264021        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.3376E+00  1.3555E+00  1.0462E+00  1.9142E+00  9.9998E-01  2.1322E+00  1.1588E+00  1.0017E+00  1.2698E+00  1.0803E+00
             1.4314E+01
 PARAMETER:  3.9089E-01  4.0416E-01  1.4521E-01  7.4929E-01  9.9982E-02  8.5717E-01  2.4739E-01  1.0175E-01  3.3885E-01  1.7725E-01
             2.7613E+00
 GRADIENT:  -1.6565E+01  5.7540E+01 -1.7142E+00  1.0153E+02 -1.4750E+01  4.0441E+01 -6.4476E+00  5.3758E+00 -1.2444E+01  6.7555E+00
             1.8833E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -734.557054324716        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.2904E+00  1.4030E+00  1.4679E+00  2.0647E+00  7.5162E+00  2.1250E+00  6.6352E+00  5.6637E-01  1.4497E+00  2.9274E+00
             1.2799E+01
 PARAMETER:  3.5496E-01  4.3858E-01  4.8382E-01  8.2500E-01  2.1171E+00  8.5379E-01  1.9924E+00 -4.6851E-01  4.7134E-01  1.1741E+00
             2.6494E+00
 GRADIENT:   3.3504E+01  3.0641E+01  9.7204E-01  5.4753E+01  1.4443E-01  4.0309E+01  3.1441E+01  1.8769E-01  2.3740E+01  7.3177E-02
             2.0295E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -745.716300880761        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:      262
 NPARAMETR:  1.2463E+00  1.3527E+00  1.4376E+00  1.8821E+00  7.9680E+00  2.0680E+00  5.8546E+00  5.9909E-01  1.3190E+00  1.1502E+01
             1.2176E+01
 PARAMETER:  3.2018E-01  4.0210E-01  4.6299E-01  7.3238E-01  2.1754E+00  8.2657E-01  1.8672E+00 -4.1234E-01  3.7685E-01  2.5425E+00
             2.5994E+00
 GRADIENT:   1.7409E+01  2.9059E+01 -3.8213E+00  2.5757E+01 -6.7670E+00  2.3928E+01 -1.4785E+00  1.1230E-01  1.7509E+01  1.5895E+01
             1.6707E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -768.880252220278        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  1.2417E+00  1.3569E+00  1.5933E+00  1.8859E+00  6.4571E+00  2.0424E+00  5.8639E+00  3.3771E-02  1.1658E+00  9.5593E+00
             9.7421E+00
 PARAMETER:  3.1649E-01  4.0523E-01  5.6580E-01  7.3443E-01  1.9652E+00  8.1414E-01  1.8688E+00 -3.2881E+00  2.5343E-01  2.3575E+00
             2.3765E+00
 GRADIENT:   6.0289E+01  4.1418E+01 -9.9578E+00  1.0105E+02 -7.2145E+00  2.1469E+01  2.3948E+01  2.0477E-04 -9.2234E-01  2.1386E+01
             1.3135E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -773.623837009195        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      463
 NPARAMETR:  1.2390E+00  1.3559E+00  4.8362E+00  1.8709E+00  6.6707E+00  2.0306E+00  5.6305E+00  1.0000E-02  1.1549E+00  6.9697E+00
             9.3036E+00
 PARAMETER:  3.1434E-01  4.0449E-01  1.6761E+00  7.2644E-01  1.9977E+00  8.0835E-01  1.8282E+00 -4.8610E+00  2.4397E-01  2.0416E+00
             2.3304E+00
 GRADIENT:   3.3042E+01  3.5926E+01  2.8409E-01  6.8271E+01  1.0539E+00  1.6973E+00 -6.9750E+00  0.0000E+00 -1.0443E+00 -2.3767E-01
            -2.4376E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -774.252616258695        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      640
 NPARAMETR:  1.2394E+00  1.3555E+00  3.5822E+00  1.8706E+00  5.5894E+00  2.0344E+00  5.7261E+00  1.3149E-02  1.1821E+00  7.5215E+00
             9.5710E+00
 PARAMETER:  3.1460E-01  4.0420E-01  1.3760E+00  7.2624E-01  1.8209E+00  8.1018E-01  1.8450E+00 -4.2314E+00  2.6729E-01  2.1178E+00
             2.3587E+00
 GRADIENT:   3.1754E+01  3.5325E+01  1.1329E+00  5.9332E+01 -3.1523E+00  3.3336E+00 -3.5003E+00 -2.5937E-06  1.2404E+00  5.5590E+00
             1.1193E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -776.295003206721        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      817
 NPARAMETR:  1.2340E+00  1.3508E+00  2.1806E+00  1.8303E+00  4.1194E+00  2.0248E+00  5.7994E+00  1.0341E-02  1.1542E+00  5.6807E+00
             9.5388E+00
 PARAMETER:  3.1025E-01  4.0072E-01  8.7960E-01  7.0446E-01  1.5157E+00  8.0547E-01  1.8578E+00 -4.4716E+00  2.4345E-01  1.8371E+00
             2.3554E+00
 GRADIENT:   3.3653E+01  3.5866E+01 -1.0037E+00  5.5933E+01 -3.7117E-01  4.4732E+00 -2.3293E+00  1.2581E-05  1.7529E+00  1.7710E-01
            -1.2042E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -783.779905236595        NO. OF FUNC. EVALS.: 104
 CUMULATIVE NO. OF FUNC. EVALS.:      921
 NPARAMETR:  1.2340E+00  1.0674E+00  2.1294E+00  1.8294E+00  4.1203E+00  2.0247E+00  5.7929E+00  1.0209E-02  1.1355E+00  5.6793E+00
             9.5276E+00
 PARAMETER:  3.1024E-01  1.6520E-01  8.5585E-01  7.0399E-01  1.5159E+00  8.0544E-01  1.8566E+00 -4.4845E+00  2.2704E-01  1.8368E+00
             2.3542E+00
 GRADIENT:   5.3289E+01  2.7949E+01  4.0608E-01  5.8346E+01 -8.1116E-01  2.2644E+01  3.8268E+01  1.3427E-05  3.2764E+00  3.3893E-01
             2.7626E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -783.969237658875        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:     1070
 NPARAMETR:  1.2358E+00  1.0665E+00  2.1366E+00  1.8219E+00  4.0891E+00  2.0327E+00  5.8405E+00  1.0000E-02  1.1029E+00  5.6271E+00
             9.4142E+00
 PARAMETER:  3.1175E-01  1.6436E-01  8.5921E-01  6.9989E-01  1.5083E+00  8.0938E-01  1.8648E+00 -4.5517E+00  1.9795E-01  1.8276E+00
             2.3422E+00
 GRADIENT:   3.8867E+01  2.7588E+01 -3.2404E-02  4.4644E+01 -8.9687E-01  7.2861E+00 -1.4416E+00  0.0000E+00  9.5815E-01  2.1065E-01
            -7.6447E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -785.007051132095        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1249
 NPARAMETR:  1.2330E+00  1.0660E+00  2.1370E+00  1.7693E+00  4.0932E+00  2.0271E+00  5.8958E+00  1.0000E-02  9.5761E-01  5.6202E+00
             9.6628E+00
 PARAMETER:  3.0943E-01  1.6391E-01  8.5941E-01  6.7061E-01  1.5093E+00  8.0659E-01  1.8742E+00 -4.5684E+00  5.6682E-02  1.8264E+00
             2.3683E+00
 GRADIENT:   3.5130E+01  2.5480E+01  9.6067E-01  3.3439E+01 -1.1048E+00  8.3018E+00  2.5663E+00  0.0000E+00 -1.9796E+00  1.1494E-01
             1.9760E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -785.387939256222        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1387
 NPARAMETR:  1.2346E+00  1.0668E+00  1.9507E+00  1.7579E+00  4.0626E+00  2.0345E+00  5.9470E+00  1.0000E-02  9.8487E-01  5.5674E+00
             9.5315E+00
 PARAMETER:  3.1072E-01  1.6469E-01  7.6820E-01  6.6409E-01  1.5018E+00  8.1027E-01  1.8829E+00 -4.5701E+00  8.4753E-02  1.8169E+00
             2.3546E+00
 GRADIENT:   5.5932E+01  2.6587E+01  2.3878E-01  4.6124E+01 -1.5303E+00  2.4404E+01  4.6008E+01  0.0000E+00 -5.7128E-01  4.7320E-01
             3.0936E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -785.514274105394        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:     1523
 NPARAMETR:  1.2342E+00  1.0668E+00  1.9930E+00  1.7532E+00  4.0626E+00  2.0331E+00  5.9016E+00  1.0000E-02  1.0141E+00  5.5643E+00
             9.3805E+00
 PARAMETER:  3.1044E-01  1.6467E-01  7.8963E-01  6.6147E-01  1.5018E+00  8.0955E-01  1.8752E+00 -4.5701E+00  1.1399E-01  1.8164E+00
             2.3386E+00
 GRADIENT:   4.0744E+01  2.5933E+01  1.9597E-02  3.0699E+01 -1.5513E+00  7.4102E+00  3.3588E+00  0.0000E+00  3.5387E-01  2.8690E-01
            -6.1388E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -786.392271265914        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1703
 NPARAMETR:  1.2274E+00  1.0659E+00  1.6804E+00  1.6795E+00  4.0785E+00  2.0153E+00  5.5992E+00  1.0000E-02  1.0561E+00  5.5384E+00
             9.2377E+00
 PARAMETER:  3.0486E-01  1.6382E-01  6.1902E-01  6.1847E-01  1.5057E+00  8.0078E-01  1.8226E+00 -4.5701E+00  1.5460E-01  1.8117E+00
             2.3233E+00
 GRADIENT:   4.4083E+01  2.4052E+01 -1.8636E+00  1.6829E+01 -3.4836E+00  4.9834E+00  8.7807E-02  0.0000E+00  4.3428E+00  8.6663E-01
            -1.6402E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -786.548905268745        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1869
 NPARAMETR:  1.2292E+00  1.0647E+00  1.9931E+00  1.6737E+00  4.0480E+00  2.0233E+00  5.6453E+00  1.0000E-02  9.4493E-01  5.4855E+00
             9.1357E+00
 PARAMETER:  3.0633E-01  1.6270E-01  7.8970E-01  6.1504E-01  1.4982E+00  8.0471E-01  1.8308E+00 -4.5701E+00  4.3359E-02  1.8021E+00
             2.3122E+00
 GRADIENT:   6.0547E+01  2.3785E+01  1.1966E-01  3.6717E+01 -1.4945E+00  2.1023E+01  4.4250E+01  0.0000E+00 -1.5172E-01  3.7335E-01
            -2.7202E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -786.553694400884        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1941
 NPARAMETR:  1.2290E+00  1.0647E+00  1.9532E+00  1.6725E+00  4.0481E+00  2.0226E+00  5.6215E+00  1.0000E-02  9.7287E-01  5.4848E+00
             9.1370E+00
 PARAMETER:  3.0620E-01  1.6269E-01  7.6948E-01  6.1432E-01  1.4983E+00  8.0440E-01  1.8266E+00 -4.5701E+00  7.2498E-02  1.8020E+00
             2.3123E+00
 GRADIENT:   6.0815E+01  2.3692E+01  2.1445E-02  3.4227E+01 -1.6781E+00  2.1044E+01  4.4014E+01  0.0000E+00  1.1930E+00  4.0572E-01
            -1.8168E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -786.555841895232        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:     2017
 NPARAMETR:  1.2287E+00  1.0647E+00  1.9343E+00  1.6702E+00  4.0485E+00  2.0215E+00  5.5765E+00  1.0000E-02  9.8973E-01  5.4836E+00
             9.1396E+00
 PARAMETER:  3.0595E-01  1.6267E-01  7.5974E-01  6.1293E-01  1.4983E+00  8.0382E-01  1.8186E+00 -4.5701E+00  8.9672E-02  1.8018E+00
             2.3126E+00
 GRADIENT:   6.0780E+01  2.3485E+01 -3.7552E-02  3.2946E+01 -1.7703E+00  2.0975E+01  4.2893E+01  0.0000E+00  1.9819E+00  4.1810E-01
            -1.0936E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -786.556500614081        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:     2098
 NPARAMETR:  1.2283E+00  1.0646E+00  1.9255E+00  1.6673E+00  4.0490E+00  2.0200E+00  5.5223E+00  1.0000E-02  1.0005E+00  5.4822E+00
             9.1427E+00
 PARAMETER:  3.0564E-01  1.6264E-01  7.5520E-01  6.1123E-01  1.4985E+00  8.0309E-01  1.8088E+00 -4.5701E+00  1.0054E-01  1.8015E+00
             2.3130E+00
 GRADIENT:   6.0598E+01  2.3215E+01 -7.5637E-02  3.2273E+01 -1.8109E+00  2.0865E+01  4.1336E+01  0.0000E+00  2.4686E+00  4.1782E-01
            -4.7498E-01

0ITERATION NO.:   89    OBJECTIVE VALUE:  -786.585901257848        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     2257
 NPARAMETR:  1.2283E+00  1.0647E+00  1.9262E+00  1.6674E+00  4.1097E+00  2.0200E+00  5.5230E+00  1.0000E-02  9.9931E-01  5.4757E+00
             9.1413E+00
 PARAMETER:  3.0564E-01  1.6265E-01  7.5528E-01  6.1124E-01  1.5134E+00  8.0310E-01  1.8088E+00 -4.5701E+00  1.0032E-01  1.8013E+00
             2.3129E+00
 GRADIENT:   1.0054E+03 -3.6091E+03 -1.3622E-01 -9.4799E+02  3.8825E+02 -3.5442E+02 -3.9812E+02  0.0000E+00  2.2378E+00  2.1620E-01
             6.0549E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2257
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.1749E-02 -9.2050E-03  7.0961E-05 -8.2429E-02 -3.9013E-03
 SE:             2.8025E-02  2.2110E-02  8.5542E-06  1.3181E-02  3.5808E-03
 N:                     100         100         100         100         100

 P VAL.:         4.3771E-01  6.7716E-01  1.0939E-16  4.0314E-10  2.7594E-01

 ETASHRINKSD(%)  6.1118E+00  2.5930E+01  9.9971E+01  5.5841E+01  8.8004E+01
 ETASHRINKVR(%)  1.1850E+01  4.5137E+01  1.0000E+02  8.0500E+01  9.8561E+01
 EBVSHRINKSD(%)  6.1702E+00  2.4885E+01  9.9954E+01  5.4343E+01  8.6020E+01
 EBVSHRINKVR(%)  1.1960E+01  4.3577E+01  1.0000E+02  7.9155E+01  9.8046E+01
 RELATIVEINF(%)  4.5404E+01  1.9415E+01  7.5645E-06  3.3895E+00  6.1987E-01
 EPSSHRINKSD(%)  9.2674E+00
 EPSSHRINKVR(%)  1.7676E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -786.58590125784838     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       132.35263194682432     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    48.37
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.51
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -786.586       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.23E+00  1.06E+00  1.93E+00  1.67E+00  4.11E+00  2.02E+00  5.52E+00  1.00E-02  1.00E+00  5.48E+00  9.14E+00
 


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
+        1.04E+05
 
 TH 2
+        5.15E+01  4.92E+05
 
 TH 3
+       -4.28E+00 -9.34E-01  1.80E+00
 
 TH 4
+        2.75E+02  5.01E+00 -4.76E+00  1.43E+04
 
 TH 5
+       -3.23E+00  9.25E+00 -1.63E+03  2.39E+00  3.81E+02
 
 TH 6
+       -1.02E+03  1.08E+01 -5.36E-01  4.46E+02 -8.02E-01  1.37E+04
 
 TH 7
+        6.01E+01 -2.10E+00  6.14E-01 -8.71E+00  1.32E-01  9.66E+01  3.58E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  4.90E-01
 
 TH 9
+       -3.13E-01 -2.01E+00  1.53E+00 -3.26E+01 -3.29E-01 -5.91E-01  2.40E+00  0.00E+00  2.60E+01
 
 TH10
+       -1.61E+00  3.64E+00 -1.03E+03  6.59E-01  2.40E+02 -2.54E-01  8.79E-02  0.00E+00 -1.08E-02  1.52E+02
 
 TH11
+       -3.36E+01 -3.54E-01  7.93E-01 -9.26E+00 -1.87E-01 -4.30E+01 -8.17E-01  0.00E+00  2.43E+00 -6.07E-02  8.86E+01
 
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
 #CPUT: Total CPU Time in Seconds,       59.949
Stop Time:
Thu Sep 30 02:57:26 CDT 2021
