Sat Sep 18 11:49:34 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat68.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m68.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1738.03879270502        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.8060E+00  1.5323E+01 -2.3407E+01  5.9158E+01  2.3841E+01  3.9013E+00  2.3713E+01  9.2152E+00  4.6765E+01 -1.0964E+00
             6.7998E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1753.05172632861        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9247E-01  9.9164E-01  1.0134E+00  9.9537E-01  1.0111E+00  9.8710E-01  8.5206E-01  9.4739E-01  7.9383E-01  1.0366E+00
             8.2127E-01
 PARAMETER:  9.2440E-02  9.1605E-02  1.1335E-01  9.5357E-02  1.1107E-01  8.7018E-02 -6.0098E-02  4.5956E-02 -1.3089E-01  1.3592E-01
            -9.6909E-02
 GRADIENT:   1.2208E+01  2.4783E+01 -1.9418E+01  5.5794E+01  3.1999E+01  4.5390E-02  5.4797E+00  4.9559E+00 -4.5281E-01 -4.9000E+00
            -4.9987E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1754.51832053269        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      266
 NPARAMETR:  1.0128E+00  9.8116E-01  1.0225E+00  1.0014E+00  1.0085E+00  1.0049E+00  7.9020E-01  8.2100E-01  8.3932E-01  1.0800E+00
             8.1209E-01
 PARAMETER:  1.1275E-01  8.0980E-02  1.2223E-01  1.0142E-01  1.0844E-01  1.0492E-01 -1.3547E-01 -9.7232E-02 -7.5165E-02  1.7697E-01
            -1.0815E-01
 GRADIENT:  -6.8305E+00  1.9291E+01 -1.0363E+01  4.2828E+01  1.9334E+01  1.1750E+00  3.9777E+00  1.0113E+00  6.1608E+00 -6.2774E-01
            -9.7903E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1755.66465877629        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      443
 NPARAMETR:  1.0158E+00  8.7244E-01  9.8458E-01  1.0495E+00  9.2422E-01  1.0006E+00  8.1334E-01  6.7555E-01  7.7599E-01  1.0253E+00
             8.3251E-01
 PARAMETER:  1.1568E-01 -3.6457E-02  8.4456E-02  1.4835E-01  2.1199E-02  1.0064E-01 -1.0661E-01 -2.9223E-01 -1.5362E-01  1.2495E-01
            -8.3305E-02
 GRADIENT:   6.6070E-02  6.4508E+00  5.2868E+00  1.8242E+00 -6.3098E+00 -1.9480E-01  6.8660E-01 -8.0208E-01 -9.6210E-01 -5.3960E-01
            -7.9835E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1755.99650782185        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      618
 NPARAMETR:  1.0128E+00  6.9282E-01  1.1932E+00  1.1700E+00  9.4844E-01  9.9715E-01  7.1452E-01  8.8525E-01  7.4787E-01  1.0787E+00
             8.3531E-01
 PARAMETER:  1.1273E-01 -2.6698E-01  2.7661E-01  2.5704E-01  4.7068E-02  9.7146E-02 -2.3614E-01 -2.1885E-02 -1.9053E-01  1.7574E-01
            -7.9957E-02
 GRADIENT:   2.2821E-01  6.2379E+00  3.0322E+00  9.6133E+00 -4.8684E+00 -1.4126E-01 -2.1899E-01 -2.4976E-01  8.1212E-01  1.8909E-01
             7.9720E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1756.43334595207        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      794
 NPARAMETR:  1.0087E+00  4.5696E-01  1.3375E+00  1.3252E+00  9.2371E-01  9.9432E-01  6.3389E-01  1.0021E+00  6.8477E-01  1.0818E+00
             8.3293E-01
 PARAMETER:  1.0864E-01 -6.8317E-01  3.9083E-01  3.8157E-01  2.0644E-02  9.4301E-02 -3.5588E-01  1.0211E-01 -2.7868E-01  1.7862E-01
            -8.2809E-02
 GRADIENT:  -9.7956E-01  9.1245E+00  2.6688E+00  3.3141E+01 -7.0016E+00  3.0620E-01 -1.4131E-01 -6.0084E-01  8.8881E-01 -1.9202E-01
            -8.8992E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1757.18366243643        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      972
 NPARAMETR:  1.0056E+00  2.3941E-01  1.5191E+00  1.4614E+00  9.1910E-01  9.9044E-01  5.7567E-01  1.1567E+00  6.1001E-01  1.0846E+00
             8.3282E-01
 PARAMETER:  1.0555E-01 -1.3296E+00  5.1809E-01  4.7942E-01  1.5642E-02  9.0397E-02 -4.5222E-01  2.4553E-01 -3.9428E-01  1.8125E-01
            -8.2942E-02
 GRADIENT:   6.3084E-01  4.6448E+00  5.5138E+00  2.5385E+01 -6.8381E+00  4.9401E-01 -7.9900E-02 -2.9478E+00 -6.2888E+00 -1.4610E+00
            -1.6177E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1757.89408445863        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1148
 NPARAMETR:  1.0031E+00  1.1654E-01  1.7482E+00  1.5429E+00  9.5391E-01  9.8650E-01  4.8487E-01  1.4299E+00  5.8355E-01  1.1130E+00
             8.3271E-01
 PARAMETER:  1.0305E-01 -2.0495E+00  6.5857E-01  5.3365E-01  5.2815E-02  8.6408E-02 -6.2387E-01  4.5764E-01 -4.3862E-01  2.0708E-01
            -8.3073E-02
 GRADIENT:  -2.1282E-01  1.1812E+00 -3.8142E-01  8.9131E+00 -2.0119E+00 -4.9460E-02  2.0820E-03  1.3954E+00 -1.3247E+00  3.3565E-01
            -7.5875E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1758.27689697993        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1323
 NPARAMETR:  1.0016E+00  4.4016E-02  1.7512E+00  1.5858E+00  9.3749E-01  9.8505E-01  3.9424E-01  1.4286E+00  5.6734E-01  1.1044E+00
             8.3267E-01
 PARAMETER:  1.0155E-01 -3.0232E+00  6.6028E-01  5.6106E-01  3.5454E-02  8.4939E-02 -8.3078E-01  4.5670E-01 -4.6680E-01  1.9931E-01
            -8.3115E-02
 GRADIENT:  -6.3864E-01  2.8530E-01 -1.5066E+00  4.2341E+00  1.0282E-01 -1.2334E-01  5.4165E-04  8.2864E-01  1.3503E-02  5.0421E-01
             4.0397E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1758.44902838728        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1498
 NPARAMETR:  1.0014E+00  1.0000E-02  1.8480E+00  1.6106E+00  9.5126E-01  9.8462E-01  2.8023E-01  1.4975E+00  5.5732E-01  1.1109E+00
             8.3257E-01
 PARAMETER:  1.0137E-01 -4.6152E+00  7.1408E-01  5.7661E-01  5.0027E-02  8.4496E-02 -1.1721E+00  5.0382E-01 -4.8461E-01  2.0515E-01
            -8.3240E-02
 GRADIENT:   2.8494E-02  0.0000E+00  8.1027E-01  3.6597E+00 -1.1696E+00 -5.0693E-02  1.4643E-05 -1.4930E-01 -3.2891E-01 -2.1667E-01
            -9.1316E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1758.45349232068        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1674            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0013E+00  1.0000E-02  1.8204E+00  1.6085E+00  9.4566E-01  9.8469E-01  2.8552E-01  1.4767E+00  5.5814E-01  1.1084E+00
             8.3265E-01
 PARAMETER:  1.0131E-01 -4.5472E+00  6.9904E-01  5.7527E-01  4.4131E-02  8.4571E-02 -1.1535E+00  4.8978E-01 -4.8314E-01  2.0293E-01
            -8.3145E-02
 GRADIENT:   6.4303E+01  0.0000E+00  1.2619E+00  1.6811E+02  9.6052E-01  6.3372E+00  6.4449E-05  2.2818E-01  3.9640E+00  3.0527E-01
             9.0143E-02

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1758.45349232068        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1770
 NPARAMETR:  1.0013E+00  1.0000E-02  1.8204E+00  1.6085E+00  9.4566E-01  9.8469E-01  2.8585E-01  1.4767E+00  5.5814E-01  1.1084E+00
             8.3265E-01
 PARAMETER:  1.0131E-01 -4.5472E+00  6.9904E-01  5.7527E-01  4.4131E-02  8.4571E-02 -1.1535E+00  4.8978E-01 -4.8314E-01  2.0293E-01
            -8.3145E-02
 GRADIENT:  -1.9086E-03  0.0000E+00  1.6440E-02 -2.8362E-03  8.4929E-03 -3.9630E-04 -4.1156E-04  3.4286E-03 -1.7549E-03 -2.8719E-03
             4.8613E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1770
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -9.0279E-05 -1.1586E-04 -3.6506E-02 -9.1739E-03 -4.2317E-02
 SE:             2.9908E-02  6.3686E-05  1.7992E-02  2.9040E-02  2.1453E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9759E-01  6.8869E-02  4.2459E-02  7.5207E-01  4.8544E-02

 ETASHRINKSD(%)  1.0000E-10  9.9787E+01  3.9723E+01  2.7129E+00  2.8131E+01
 ETASHRINKVR(%)  1.0000E-10  1.0000E+02  6.3667E+01  5.3522E+00  4.8348E+01
 EBVSHRINKSD(%)  2.8811E-01  9.9798E+01  4.4087E+01  3.0917E+00  2.3377E+01
 EBVSHRINKVR(%)  5.7540E-01  1.0000E+02  6.8737E+01  6.0878E+00  4.1290E+01
 RELATIVEINF(%)  9.7608E+01  1.7065E-05  7.9724E+00  4.3113E+00  1.1542E+01
 EPSSHRINKSD(%)  4.5830E+01
 EPSSHRINKVR(%)  7.0657E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1758.4534923206816     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1023.3026657569434     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.21
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1758.453       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.00E-02  1.82E+00  1.61E+00  9.46E-01  9.85E-01  2.86E-01  1.48E+00  5.58E-01  1.11E+00  8.33E-01
 


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
+        1.15E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        5.79E-01  0.00E+00  5.23E+01
 
 TH 4
+       -1.15E+01  0.00E+00 -3.57E+01  1.30E+03
 
 TH 5
+       -3.11E+01  0.00E+00 -1.33E+02 -8.55E+01  6.48E+02
 
 TH 6
+        3.60E+01  0.00E+00  3.57E-01 -1.03E+00  1.08E+01  2.11E+02
 
 TH 7
+        6.79E+00  0.00E+00 -8.34E-01 -2.44E-01  4.67E+00 -6.22E+00  2.54E+00
 
 TH 8
+       -1.08E+00  0.00E+00 -1.87E+01 -5.31E+00 -6.62E+00  2.68E+00 -6.94E-01  2.36E+01
 
 TH 9
+       -6.88E+00  0.00E+00  5.81E+00 -7.76E-01  1.11E+01  3.74E+00 -3.83E+00  1.46E+00  5.69E+02
 
 TH10
+        4.76E+00  0.00E+00 -2.15E+00 -2.29E+00 -7.83E+01 -4.58E+00  8.86E-01  9.89E+00  2.27E+00  6.61E+01
 
 TH11
+       -1.31E+01  0.00E+00 -2.83E+00 -1.60E+01  3.45E+01 -8.65E+00 -2.96E+00  3.92E+00  1.22E+01 -3.07E+00  3.28E+02
 
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
 #CPUT: Total CPU Time in Seconds,       26.924
Stop Time:
Sat Sep 18 11:50:02 CDT 2021
