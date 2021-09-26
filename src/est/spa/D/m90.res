Sat Sep 25 14:45:59 CDT 2021
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
$DATA ../../../../data/spa/D/dat90.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m90.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   27553.1583264922        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.0792E+02  5.5355E+02 -5.9192E+00  5.0541E+02 -1.5887E+01 -2.6561E+03 -1.2704E+03 -5.1257E+01 -1.8151E+03 -5.0739E+02
            -5.1701E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -425.916168014018        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2402E+00  9.9437E-01  9.0194E-01  1.2895E+00  1.3578E+00  1.7402E+00  9.8410E-01  9.6497E-01  8.5417E-01  9.1718E-01
             1.5061E+01
 PARAMETER:  3.1530E-01  9.4351E-02 -3.2076E-03  3.5427E-01  4.0587E-01  6.5399E-01  8.3976E-02  6.4339E-02 -5.7624E-02  1.3544E-02
             2.8121E+00
 GRADIENT:   2.0447E+01  1.0371E+00  2.0199E-01 -1.2758E+01 -3.7035E+00  3.6438E+01  2.3983E-01  3.2166E+00  3.0082E+00  1.5390E+00
            -5.6769E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -439.912213823682        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.2231E+00  5.3388E-01  8.2998E-01  1.6854E+00  1.0424E+00  1.4749E+00  5.1923E-01  3.6722E-01  2.4760E-01  6.4459E-01
             1.6985E+01
 PARAMETER:  3.0138E-01 -5.2759E-01 -8.6352E-02  6.2201E-01  1.4150E-01  4.8858E-01 -5.5541E-01 -9.0179E-01 -1.2959E+00 -3.3913E-01
             2.9323E+00
 GRADIENT:  -3.5666E+01  1.5258E+01  2.4246E+00  4.8242E+01 -1.5875E+00 -8.1814E+00  3.3173E-03  9.1490E-01  7.3263E-01  1.0030E+00
            -1.3896E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -469.005417314178        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      237
 NPARAMETR:  1.4024E+00  2.2873E-01  1.7028E-01  1.5681E+00  1.6522E-01  1.9382E+00  3.4018E+00  1.0000E-02  4.4695E-02  1.0000E-02
             1.8270E+01
 PARAMETER:  4.3821E-01 -1.3752E+00 -1.6703E+00  5.4986E-01 -1.7004E+00  7.6175E-01  1.3243E+00 -9.9657E+00 -3.0079E+00 -5.5712E+00
             3.0053E+00
 GRADIENT:  -3.6502E+01  2.4268E+01  9.1753E+01  5.3445E+01 -1.4602E+02 -3.9275E+01  1.5947E+00  0.0000E+00  3.2241E-02  0.0000E+00
            -1.4428E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -527.555979372219        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      308
 NPARAMETR:  1.5734E+00  6.1520E-02  3.3532E-02  9.0177E-01  1.0241E-01  2.6238E+00  3.0586E+00  1.0000E-02  6.2900E-01  1.0000E-02
             1.5865E+01
 PARAMETER:  5.5327E-01 -2.6884E+00 -3.2952E+00 -3.3989E-03 -2.1787E+00  1.0646E+00  1.2180E+00 -2.4891E+01 -3.6362E-01 -1.8395E+01
             2.8641E+00
 GRADIENT:   3.8408E+01 -1.0820E-01 -1.7193E+00  2.9277E+00  3.2245E+01  6.1172E+00  1.2264E+00  0.0000E+00  2.4224E-01  0.0000E+00
             9.9376E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -532.280378558601        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      381
 NPARAMETR:  1.3396E+00  5.3539E-02  2.2953E-02  7.4290E-01  9.2601E-02  2.5794E+00  2.4801E+00  1.0000E-02  1.1994E+00  1.0000E-02
             1.4913E+01
 PARAMETER:  3.9236E-01 -2.8273E+00 -3.6743E+00 -1.9720E-01 -2.2795E+00  1.0476E+00  1.0083E+00 -2.6882E+01  2.8180E-01 -1.9452E+01
             2.8022E+00
 GRADIENT:   3.7507E-01 -9.3406E-01  1.1716E+00 -1.4633E+00 -9.3502E-01  9.0710E-01  7.3389E-01  0.0000E+00  2.9448E-01  0.0000E+00
             3.4413E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -534.272781741835        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      456
 NPARAMETR:  1.3230E+00  1.0924E-01  1.8936E-02  7.6741E-01  9.1879E-02  2.5806E+00  1.0090E+00  1.0000E-02  7.4113E-01  1.0000E-02
             1.4617E+01
 PARAMETER:  3.7986E-01 -2.1142E+00 -3.8667E+00 -1.6473E-01 -2.2873E+00  1.0480E+00  1.0897E-01 -2.4032E+01 -1.9959E-01 -1.3666E+01
             2.7822E+00
 GRADIENT:   1.9001E+00  6.0225E-01 -3.1527E+00 -9.7387E-01  9.0119E+00 -3.2349E+00  1.5842E+00  0.0000E+00  3.5137E-02  0.0000E+00
            -7.5212E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -536.516140024976        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  1.2682E+00  1.2174E-01  2.3798E-02  1.4905E+00  9.3295E-02  2.7444E+00  2.4487E-01  1.0000E-02  1.9488E-01  1.4144E-02
             1.4843E+01
 PARAMETER:  3.3756E-01 -2.0059E+00 -3.6382E+00  4.9909E-01 -2.2720E+00  1.1096E+00 -1.3070E+00 -1.6572E+01 -1.5353E+00 -4.1585E+00
             2.7975E+00
 GRADIENT:   2.0156E+00  6.1722E-01 -2.9648E+00 -8.3176E-01  1.0605E+01  1.0302E+00  1.1491E-01  0.0000E+00  7.8276E-03  1.2566E-02
             2.0278E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -537.015140781194        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      600
 NPARAMETR:  1.2390E+00  1.1661E-01  2.6288E-02  2.1610E+00  9.0246E-02  2.8825E+00  1.3836E-01  1.0000E-02  2.2590E-01  6.4691E-02
             1.4874E+01
 PARAMETER:  3.1432E-01 -2.0490E+00 -3.5386E+00  8.7056E-01 -2.3052E+00  1.1586E+00 -1.8779E+00 -1.3888E+01 -1.3877E+00 -2.6381E+00
             2.7996E+00
 GRADIENT:   3.6054E+00  1.5247E+00  1.9643E+00 -2.6679E+00 -1.4937E+01  3.8104E+00  3.2966E-02  0.0000E+00  1.2933E-02  2.5829E-01
             7.0044E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -538.499061771238        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      671
 NPARAMETR:  1.0572E+00  1.1439E-01  3.7055E-02  5.4841E+00  9.0007E-02  3.2003E+00  4.1578E-02  1.0000E-02  1.3241E+00  1.8728E-02
             1.4830E+01
 PARAMETER:  1.5562E-01 -2.0681E+00 -3.1953E+00  1.8019E+00 -2.3079E+00  1.2632E+00 -3.0802E+00 -7.9204E+00  3.8070E-01 -3.8777E+00
             2.7966E+00
 GRADIENT:  -2.2108E-01 -8.6312E-02 -4.1784E-01  5.5643E-02  2.6905E-01 -7.9174E-01  2.8251E-03  0.0000E+00  3.8170E-02  2.2938E-02
             1.2768E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -538.552393031853        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      828
 NPARAMETR:  1.0323E+00  1.1571E-01  3.9601E-02  6.2535E+00  9.0342E-02  3.2919E+00  3.5789E-02  1.0000E-02  1.4959E+00  1.7589E-02
             1.4914E+01
 PARAMETER:  1.3178E-01 -2.0567E+00 -3.1289E+00  1.9331E+00 -2.3042E+00  1.2915E+00 -3.2301E+00 -6.9906E+00  5.0271E-01 -3.9405E+00
             2.8023E+00
 GRADIENT:  -2.3366E-02 -8.7702E-02  9.4237E-04 -1.2500E-01  2.0163E-01  4.6736E-02  2.0659E-03  0.0000E+00  9.2524E-02  2.0171E-02
             1.5427E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -538.562183941149        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1005
 NPARAMETR:  1.0295E+00  1.1601E-01  3.9631E-02  6.2395E+00  9.0414E-02  3.2980E+00  4.0472E-02  1.0000E-02  1.4338E+00  1.0000E-02
             1.4875E+01
 PARAMETER:  1.2907E-01 -2.0541E+00 -3.1281E+00  1.9309E+00 -2.3034E+00  1.2933E+00 -3.1071E+00 -7.0249E+00  4.6034E-01 -4.6904E+00
             2.7997E+00
 GRADIENT:   1.4975E-02  1.4357E-02 -2.5598E-02  5.6119E-03  2.2369E-01  1.2104E-01  2.6773E-03  0.0000E+00 -7.6596E-03  0.0000E+00
            -1.4133E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -538.562410905308        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1184
 NPARAMETR:  1.0291E+00  1.1577E-01  3.9649E-02  6.2843E+00  9.0324E-02  3.2980E+00  3.8935E-02  1.0000E-02  1.4592E+00  1.0000E-02
             1.4879E+01
 PARAMETER:  1.2865E-01 -2.0562E+00 -3.1277E+00  1.9381E+00 -2.3044E+00  1.2933E+00 -3.1459E+00 -6.9721E+00  4.7791E-01 -4.5127E+00
             2.8000E+00
 GRADIENT:   1.8984E-02 -4.0234E-03 -5.6235E-03  3.5938E-03 -4.0075E-03 -1.2847E-02  2.4629E-03  0.0000E+00  5.7977E-03  0.0000E+00
             1.5144E-02

0ITERATION NO.:   64    OBJECTIVE VALUE:  -538.562426577245        NO. OF FUNC. EVALS.: 131
 CUMULATIVE NO. OF FUNC. EVALS.:     1315
 NPARAMETR:  1.0289E+00  1.1581E-01  3.9655E-02  6.2787E+00  9.0338E-02  3.2983E+00  3.9001E-02  1.0000E-02  1.4532E+00  1.0000E-02
             1.4879E+01
 PARAMETER:  1.2850E-01 -2.0558E+00 -3.1275E+00  1.9372E+00 -2.3042E+00  1.2934E+00 -3.1442E+00 -6.9725E+00  4.7380E-01 -4.5079E+00
             2.7999E+00
 GRADIENT:  -3.6734E-03  1.0039E-03 -3.4270E-03  2.9599E-03  3.8833E-03  2.3421E-03  2.4737E-03  0.0000E+00  2.3927E-04  1.3248E-03
             7.9611E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1315
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.6105E-03 -7.6135E-05 -1.3544E-04 -1.3748E-02  1.4140E-04
 SE:             2.7572E-02  8.0840E-05  6.3271E-05  6.0218E-03  2.3617E-04
 N:                     100         100         100         100         100

 P VAL.:         7.5482E-01  3.4630E-01  3.2301E-02  2.2428E-02  5.4936E-01

 ETASHRINKSD(%)  7.6299E+00  9.9729E+01  9.9788E+01  7.9826E+01  9.9209E+01
 ETASHRINKVR(%)  1.4678E+01  9.9999E+01  1.0000E+02  9.5930E+01  9.9994E+01
 EBVSHRINKSD(%)  6.1397E+00  9.9559E+01  9.9747E+01  8.2370E+01  9.9049E+01
 EBVSHRINKVR(%)  1.1902E+01  9.9998E+01  9.9999E+01  9.6892E+01  9.9991E+01
 RELATIVEINF(%)  3.9858E+01  1.5407E-03  1.3830E-04  6.5060E-01  2.7812E-03
 EPSSHRINKSD(%)  3.5883E+00
 EPSSHRINKVR(%)  7.0478E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -538.56242657724533     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       196.58839998649285     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.81
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.80
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -538.562       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.16E-01  3.97E-02  6.28E+00  9.03E-02  3.30E+00  3.90E-02  1.00E-02  1.45E+00  1.00E-02  1.49E+01
 


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
+        8.33E+01
 
 TH 2
+       -3.10E+00  1.54E+03
 
 TH 3
+        3.70E+02  2.28E+03  4.76E+04
 
 TH 4
+        2.13E+00 -6.26E+00 -9.00E+01  4.79E-01
 
 TH 5
+        2.15E+02 -5.87E+03 -4.34E+04  1.32E+02  8.28E+04
 
 TH 6
+        1.53E+00 -3.66E+00 -1.71E+02 -8.38E-01  4.75E+01  1.29E+01
 
 TH 7
+        8.24E-02  1.10E+00 -8.85E-02  9.09E-03  1.46E+00 -2.25E-03  7.98E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -9.75E-01  1.79E+00  2.12E+01 -3.41E-01 -1.86E+01  9.00E-01  1.09E-01  0.00E+00  6.49E-01
 
 TH10
+        1.57E+00  2.59E+00  3.14E-01 -4.77E-02  1.81E+00  6.49E-02 -4.43E+00  0.00E+00  1.59E-01  2.56E+03
 
 TH11
+       -4.53E+00 -8.66E+00 -7.23E+01 -5.75E-02  4.30E+01  4.20E-01 -4.50E-03  0.00E+00  1.02E-01  1.06E-02  1.65E+00
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       21.674
Stop Time:
Sat Sep 25 14:46:25 CDT 2021
